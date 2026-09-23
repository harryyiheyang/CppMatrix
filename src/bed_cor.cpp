#include <Rcpp.h>
#include <array>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include "geno_kernel.h"


namespace {

constexpr std::size_t kReadGroup = 64;     // variants per file read

struct Bed {
  std::string path;
  std::size_t samples = 0, variants = 0, stride = 0;
  std::vector<std::string> sample_ids, variant_ids;
};

std::string stem(const std::string& path) {
  if (path.size() < 4 || path.substr(path.size() - 4) != ".bed")
    Rcpp::stop("BED path must end in .bed: " + path);
  return path.substr(0, path.size() - 4);
}

std::vector<std::string> read_column(const std::string& path, int column,
                                     const char* what) {
  std::ifstream input(path);
  if (!input) Rcpp::stop(std::string("Cannot open ") + what + " file: " + path);
  std::vector<std::string> out;
  std::string line, a, b;
  while (std::getline(input, line)) {
    std::istringstream row(line);
    if (!(row >> a >> b)) Rcpp::stop(std::string("Invalid ") + what + " row in: " + path);
    out.push_back(column == 2 ? a + "\t" + b : b);
  }
  return out;
}

Bed open_bed(const std::string& path) {
  Bed bed;
  bed.path = path;
  const std::string base = stem(path);
  bed.sample_ids = read_column(base + ".fam", 2, "FAM");
  bed.variant_ids = read_column(base + ".bim", 1, "BIM");
  bed.samples = bed.sample_ids.size();
  bed.variants = bed.variant_ids.size();
  if (bed.samples < 2) Rcpp::stop("BED correlation requires at least two samples.");
  if (!bed.variants) Rcpp::stop("BIM file has no variants: " + base + ".bim");
  bed.stride = (bed.samples + 3) / 4;
  std::ifstream input(path, std::ios::binary | std::ios::ate);
  if (!input) Rcpp::stop("Cannot open BED file: " + path);
  const std::streamoff actual = input.tellg();
  if (actual < 0 || static_cast<std::uint64_t>(actual) != 3 + bed.variants * bed.stride)
    Rcpp::stop("BED size does not match its BIM and FAM files: " + path);
  input.seekg(0);
  unsigned char header[3];
  input.read(reinterpret_cast<char*>(header), 3);
  if (header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01)
    Rcpp::stop("BED must be a PLINK SNP-major BED file: " + path);
  return bed;
}

const NibbleTable& bed_table() {
  static const NibbleTable table({2, -1, 1, 0});  // 00:2, 01:NA, 10:1, 11:0
  return table;
}

void load_planes(const Bed& bed, std::size_t first, std::size_t m, Planes& out,
                 int threads) {
  reset_planes(out, bed.samples, m);
  std::ifstream input(bed.path, std::ios::binary);
  if (!input) Rcpp::stop("Cannot open BED file: " + bed.path);
  const NibbleTable& table = bed_table();
  std::vector<unsigned char> buffer(kReadGroup * bed.stride);
  for (std::size_t g = 0; g < m; g += kReadGroup) {
    const std::size_t count = std::min(kReadGroup, m - g);
    input.seekg(static_cast<std::streamoff>(3 + (first + g) * bed.stride));
    input.read(reinterpret_cast<char*>(buffer.data()),
               static_cast<std::streamsize>(count * bed.stride));
    if (!input) Rcpp::stop("Could not read BED file: " + bed.path);
    #pragma omp parallel for num_threads(threads) schedule(static)
    for (std::ptrdiff_t k = 0; k < static_cast<std::ptrdiff_t>(count); ++k) {
      const std::size_t j = g + k;
      encode_packed(table, buffer.data() + k * bed.stride, bed.samples, bed.stride,
                    out.words, out.bits.data() + 2 * j * out.words,
                    out.bits.data() + (2 * j + 1) * out.words, out.s[j], out.v[j]);
    }
  }
}

} // namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix bed_cor_cpp(const std::string& A, const std::string& B,
                                bool self, int threads) {
#ifdef _OPENMP
  if (threads <= 0) threads = omp_get_max_threads();
#else
  threads = 1;
#endif
  const Bed a = open_bed(A);
  const Bed b = self ? a : open_bed(B);
  if (a.sample_ids != b.sample_ids)
    Rcpp::stop("A and B must have the same FID/IID samples in the same order.");
  if (a.variants > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      b.variants > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      a.variants > 100000000ULL / b.variants)
    Rcpp::stop("Requested correlation matrix is too large; subset BED variants first.");

  Rcpp::NumericMatrix result = cor_driver(
    a.variants, b.variants, a.samples, self, threads,
    [&](std::size_t first, std::size_t m, Planes& out) { load_planes(a, first, m, out, threads); },
    [&](std::size_t first, std::size_t m, Planes& out) { load_planes(b, first, m, out, threads); });
  result.attr("dimnames") = Rcpp::List::create(Rcpp::wrap(a.variant_ids),
                                               Rcpp::wrap(b.variant_ids));
  return result;
}
