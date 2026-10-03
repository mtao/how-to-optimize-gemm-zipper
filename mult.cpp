#include "definitions.hpp"
#include <ranges>

namespace {
void mydot(zipper::concepts::Vector auto const &x,
           zipper::concepts::Vector auto const &y, scalar_type &c) {

  const index_type k = x.size();
  for (auto p : std::views::iota(index_type{0}, k)) {
    c += x(p) * y(p);
  }
}

void mydot_1x4(zipper::concepts::Vector auto const &x,
               zipper::concepts::Matrix auto const &B,
               zipper::concepts::Vector auto &c) {

  mydot(x, B.col(0), c(0));
  mydot(x, B.col(1), c(1));
  mydot(x, B.col(2), c(2));
  mydot(x, B.col(3), c(3));
}

} // namespace

void MULT_NAME(AMat const &A, BMat const &B, CMat &C) {
  const index_type m = C.rows();
  const index_type n = C.cols();
  for (auto j : std::views::iota(index_type{0}, n) | std::views::stride(4)) {
    for (auto i : std::views::iota(index_type{0}, m)) {
      auto c = C(i, zipper::slice(j, std::integral_constant<index_type, 4>{}));
      auto b = B(zipper::full_extent_t{},
                 zipper::slice(j, std::integral_constant<index_type, 4>{}));
      mydot_1x4(A.row(i), b, c);
    }
  }
}
