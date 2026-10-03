#include "definitions.hpp"
#include <argparse/argparse.hpp>
#include <chrono>
#include <format>
#include <print>
#include <zipper/expression/nullary/Random.hpp>
#include <zipper/utils/format.hpp>
#include <zipper/utils/matrix_norm.hpp>

namespace {
void run(auto const &sizes) {
  for (auto [m, n, k] : sizes) {

    // 1e-9 for gigaflop, 1e3 because the clock will be in ms
    double gflops = 2.0 * m * n * k * 1e-9 * 1e3;

    // TODO: original document had the ability to force these parameters
    index_type a_rows = m;
    index_type b_rows = k;
    index_type c_rows = m;

    AMat A(a_rows, k);
    BMat B(b_rows, n);
    CMat C(c_rows, n);

    A = zipper::expression::nullary::uniform_random<scalar_type>();
    B = zipper::expression::nullary::uniform_random<scalar_type>();
    C = zipper::expression::nullary::uniform_random<scalar_type>();

    // std::println("A:\n{}", A);
    // std::println("B:\n{}", B);
    // Baseline build-in matrix product. It might change over time but it should
    // be a more consistent baseline
    zipper::Matrix C_ref =
        C + zipper::MatrixBase(zipper::expression::binary::MatrixProduct(
                A.expression(), B.expression()));
    zipper::Matrix C_new = C;
    zipper::Matrix C_old = C;

    std::optional<double> best_time_opt;

    for (int _ : std::views::iota(0, num_repetitions)) {

      // TODO: verify this is a memcpy
      C = C_old;

      auto start = std::chrono::steady_clock::now();

      MULT_NAME(A, B, C);

      auto end = std::chrono::steady_clock::now();

      std::chrono::duration<double, std::milli> elapsed = end - start;

      if (const double e = elapsed.count(); best_time_opt.has_value()) {
        double &min = best_time_opt.value();
        min = std::min(min, e);
      } else {
        best_time_opt = e;
      }
    }

    std::println("{}x{}x{} {} {}", m, n, k, gflops / (*best_time_opt),
                 zipper::utils::frobenius_norm(C - C_ref));
  }
}

index_type read_token(std::string_view tok) {

  index_type r;
  auto begin = tok.data();
  auto end = tok.data() + tok.size();
  auto [p, ec] = std::from_chars(begin, end, r);
  if (ec != std::errc{} || p != end || r == 0) {
    throw std::invalid_argument(std::format("Bad dimension: {}", tok));
  }
  return r;
}
} // namespace

auto main(int argc, char *argv[]) -> int {

  argparse::ArgumentParser args(std::string{mult_name});

  args.add_argument("sizes")
      .help("N or MxNxK")
      .nargs(argparse::nargs_pattern::any)
      .default_value(std::vector<std::string>{});

  try {

    args.parse_args(argc, argv);
  } catch (const std::exception &e) {
    std::println(stderr, "{}\n{}", e.what(), args.help().str());
    return 1;
  }

  auto sizes = args.get<std::vector<std::string>>("sizes");
  if (sizes.empty()) {

    run(test_row_sizes | std::views::transform([](index_type a) {
          // return std::make_tuple(a, a, a);
          return std::array{a, a, a};
        }));
  } else {
    std::vector<std::array<index_type, 3>> dims;
    try {

      dims = sizes | std::views::transform([](const std::string &size_str) {
               auto dims = std::views::split(size_str, 'x') |
                           std::views::transform([](auto const &tok) {
                             return read_token(std::string_view(tok));
                           }) |
                           std::ranges::to<std::vector>();
               std::array<index_type, 3> r;
               if (dims.size() == 1) {
                 std::ranges::fill(r, dims[0]);
               } else if (dims.size() == 3) {
                 std::ranges::copy(dims, r.begin());
               } else {
                 throw std::invalid_argument(
                     std::format("Invalid matrix shape: {}", dims));
               }
               return r;
             }) |
             std::ranges::to<std::vector>();
    } catch (const std::exception &e) {
      std::println(stderr, "Dimensions weren't parsable: {}", e.what());
      return 1;
    }
    run(dims);
  }

  return 0;
}
