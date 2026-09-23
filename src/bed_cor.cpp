#include <Rcpp.h>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct ByteStats {
  unsigned char n, x, y, xy;
};

struct BedData {
  std::size_t samples, variants, stride;
  std::vector<std::string> sample_ids;
  std::vector<unsigned char> data;
  std::vector<double> mean, norm;
};

std::string stem(const std::string& path) {
  if (path.size() < 4 || path.substr(path.size() - 4) != ".bed") {
    Rcpp::stop("BED path must end in .bed: " + path);
  }
  return path.substr(0, path.size() - 4);
}

std::vector<std::string> read_samples(const std::string& path) {
  std::ifstream input(path);
  if (!input) Rcpp::stop("Cannot open FAM file: " + path);
  std::vector<std::string> ids;
  std::string line, fid, iid;
  while (std::getline(input, line)) {
    std::istringstream row(line);
    if (!(row >> fid >> iid)) Rcpp::stop("Invalid FAM row in: " + path);
    ids.push_back(fid + "\t" + iid);
  }
  if (ids.size() < 2) Rcpp::stop("BED correlation requires at least two samples.");
  return ids;
}

std::size_t count_variants(const std::string& path) {
  std::ifstream input(path);
  if (!input) Rcpp::stop("Cannot open BIM file: " + path);
  std::size_t count = 0;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty()) Rcpp::stop("Empty BIM row in: " + path);
    ++count;
  }
  if (!count) Rcpp::stop("BIM file has no variants: " + path);
  return count;
}

int dosage(unsigned char code) {
  // PLINK allele 1 copies: 00 -> 2, 10 -> 1, 11 -> 0, 01 -> missing.
  return code == 0 ? 2 : (code == 2 ? 1 : (code == 3 ? 0 : -1));
}

const std::array<ByteStats, 65536>& byte_table() {
  static const std::array<ByteStats, 65536> table = [] {
    std::array<ByteStats, 65536> result{};
    for (unsigned int a = 0; a < 256; ++a) {
      for (unsigned int b = 0; b < 256; ++b) {
        ByteStats s{};
        for (unsigned int k = 0; k < 4; ++k) {
          const int x = dosage((a >> (2 * k)) & 3);
          const int y = dosage((b >> (2 * k)) & 3);
          if (x < 0 || y < 0) continue;
          ++s.n;
          s.x += x;
          s.y += y;
          s.xy += x * y;
        }
        result[a | (b << 8)] = s;
      }
    }
    return result;
  }();
  return table;
}

BedData read_bed(const std::string& path) {
  const std::string base = stem(path);
  BedData bed;
  bed.sample_ids = read_samples(base + ".fam");
  bed.samples = bed.sample_ids.size();
  bed.variants = count_variants(base + ".bim");
  bed.stride = (bed.samples + 3) / 4;
  if (bed.variants > (std::numeric_limits<std::size_t>::max() - 3) / bed.stride)
    Rcpp::stop("BED dimensions exceed supported file size.");
  const std::size_t expected = 3 + bed.variants * bed.stride;
  std::ifstream input(path, std::ios::binary | std::ios::ate);
  if (!input) Rcpp::stop("Cannot open BED file: " + path);
  const std::streamoff actual = input.tellg();
  if (actual < 0 || static_cast<std::uint64_t>(actual) != expected)
    Rcpp::stop("BED size does not match its BIM and FAM files: " + path);
  input.seekg(0);
  unsigned char header[3];
  input.read(reinterpret_cast<char*>(header), 3);
  if (header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01)
    Rcpp::stop("BED must be a PLINK SNP-major BED file: " + path);
  bed.data.resize(expected - 3);
  input.read(reinterpret_cast<char*>(bed.data.data()), bed.data.size());
  if (!input) Rcpp::stop("Could not read complete BED file: " + path);

  bed.mean.resize(bed.variants);
  bed.norm.resize(bed.variants);
  for (std::size_t j = 0; j < bed.variants; ++j) {
    const unsigned char* column = bed.data.data() + j * bed.stride;
    std::uint64_t n = 0, sum = 0, square = 0;
    for (std::size_t i = 0; i < bed.samples; ++i) {
      const int x = dosage((column[i / 4] >> (2 * (i % 4))) & 3);
      if (x < 0) continue;
      ++n;
      sum += x;
      square += x * x;
    }
    bed.mean[j] = n ? static_cast<double>(sum) / n : NA_REAL;
    const double ss = n ? static_cast<double>(square) -
      static_cast<double>(sum) * sum / n : 0.0;
    bed.norm[j] = n >= 2 && ss > 0 ? std::sqrt(ss) : NA_REAL;
    if ((j & 255) == 0) Rcpp::checkUserInterrupt();
  }
  return bed;
}

double correlation(const unsigned char* a, const unsigned char* b,
                   std::size_t samples, std::size_t stride,
                   double mean_a, double mean_b, double norm_a, double norm_b) {
  if (Rcpp::NumericVector::is_na(norm_a) || Rcpp::NumericVector::is_na(norm_b))
    return NA_REAL;
  const auto& table = byte_table();
  std::uint64_t n = 0, x = 0, y = 0, xy = 0;
  const std::size_t full = samples / 4;
  for (std::size_t k = 0; k < full; ++k) {
    const ByteStats& s = table[a[k] | (static_cast<unsigned int>(b[k]) << 8)];
    n += s.n;
    x += s.x;
    y += s.y;
    xy += s.xy;
  }
  if (samples % 4) {
    const unsigned char aa = a[stride - 1];
    const unsigned char bb = b[stride - 1];
    for (std::size_t k = 0; k < samples % 4; ++k) {
      const int gx = dosage((aa >> (2 * k)) & 3);
      const int gy = dosage((bb >> (2 * k)) & 3);
      if (gx < 0 || gy < 0) continue;
      ++n;
      x += gx;
      y += gy;
      xy += gx * gy;
    }
  }
  const double dot = static_cast<double>(xy) - mean_a * y -
    mean_b * x + static_cast<double>(n) * mean_a * mean_b;
  double r = dot / (norm_a * norm_b);
  if (r > 1.0) r = 1.0;
  if (r < -1.0) r = -1.0;
  return r;
}

} // namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix bed_cor_cpp(const std::string& A, const std::string& B,
                                bool self) {
  const std::size_t variants_a = count_variants(stem(A) + ".bim");
  const std::size_t variants_b = self ? variants_a : count_variants(stem(B) + ".bim");
  if (variants_a > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      variants_b > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      variants_a > 100000000ULL / variants_b)
    Rcpp::stop("Requested correlation matrix is too large; subset BED variants first.");
  BedData a = read_bed(A);
  BedData b;
  if (!self) {
    b = read_bed(B);
    if (a.sample_ids != b.sample_ids)
      Rcpp::stop("A and B must have the same FID/IID samples in the same order.");
  }
  const BedData& right = self ? a : b;
  Rcpp::NumericMatrix result(a.variants, right.variants);
  for (std::size_t j = 0; j < right.variants; ++j) {
    const unsigned char* column_b = right.data.data() + j * right.stride;
    const std::size_t end = self ? j + 1 : a.variants;
    for (std::size_t i = 0; i < end; ++i) {
      const unsigned char* column_a = a.data.data() + i * a.stride;
      const double r = correlation(column_a, column_b, a.samples, a.stride,
                                   a.mean[i], right.mean[j], a.norm[i], right.norm[j]);
      result(i, j) = r;
      if (self) result(j, i) = r;
    }
    if ((j & 15) == 0) Rcpp::checkUserInterrupt();
  }
  return result;
}
