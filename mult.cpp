#include "definitions.hpp"
#include <ranges>
#include <zipper/expression/nullary/Constant.hpp>

namespace {
auto mymul_4x4(zipper::concepts::Matrix auto const &A,
               zipper::concepts::Matrix auto const &B) {

  // NOTE: by using a matrix here apparently the compiler is auto-vectorizing
  // better than the original tutorial which used register variables that
  // couldn't be realized as vectors easily
  zipper::Matrix<scalar_type, 4, 4, false> C =
      zipper::expression::nullary::Constant<scalar_type>(0);
  const index_type k = B.rows();
  for (auto p : std::views::iota(index_type{0}, k)) {
    zipper::Vector a = A.col(p);
    C(0, 0) += a(0) * B(p, 0);
    C(0, 1) += a(0) * B(p, 1);
    C(0, 2) += a(0) * B(p, 2);
    C(0, 3) += a(0) * B(p, 3);
    C(1, 0) += a(1) * B(p, 0);
    C(1, 1) += a(1) * B(p, 1);
    C(1, 2) += a(1) * B(p, 2);
    C(1, 3) += a(1) * B(p, 3);
    C(2, 0) += a(2) * B(p, 0);
    C(2, 1) += a(2) * B(p, 1);
    C(2, 2) += a(2) * B(p, 2);
    C(2, 3) += a(2) * B(p, 3);
    C(3, 0) += a(3) * B(p, 0);
    C(3, 1) += a(3) * B(p, 1);
    C(3, 2) += a(3) * B(p, 2);
    C(3, 3) += a(3) * B(p, 3);
  }
  return C;
}

} // namespace

void MULT_NAME(AMat const &A, BMat const &B, CMat &C) {
  const index_type m = C.rows();
  const index_type n = C.cols();
  for (auto j : std::views::iota(index_type{0}, n) | std::views::stride(4)) {
    for (auto i : std::views::iota(index_type{0}, m) | std::views::stride(4)) {
      auto c = C(zipper::slice(i, std::integral_constant<index_type, 4>{}),
                 zipper::slice(j, std::integral_constant<index_type, 4>{}));
      auto a = A(zipper::slice(i, std::integral_constant<index_type, 4>{}),
                 zipper::full_extent_t{});
      auto b = B(zipper::full_extent_t{},
                 zipper::slice(j, std::integral_constant<index_type, 4>{}));
      c += mymul_4x4(a, b);
    }
  }
}
