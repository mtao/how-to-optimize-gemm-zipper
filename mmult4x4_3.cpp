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

void mymul_4x4(zipper::concepts::Matrix auto const &A,
               zipper::concepts::Matrix auto const &B,
               zipper::concepts::Matrix auto &C) {

  mydot(A.row(0), B.col(0), C(0, 0));
  mydot(A.row(0), B.col(1), C(0, 1));
  mydot(A.row(0), B.col(2), C(0, 2));
  mydot(A.row(0), B.col(3), C(0, 3));

  mydot(A.row(1), B.col(0), C(1, 0));
  mydot(A.row(1), B.col(1), C(1, 1));
  mydot(A.row(1), B.col(2), C(1, 2));
  mydot(A.row(1), B.col(3), C(1, 3));

  mydot(A.row(2), B.col(0), C(2, 0));
  mydot(A.row(2), B.col(1), C(2, 1));
  mydot(A.row(2), B.col(2), C(2, 2));
  mydot(A.row(2), B.col(3), C(2, 3));

  mydot(A.row(3), B.col(0), C(3, 0));
  mydot(A.row(3), B.col(1), C(3, 1));
  mydot(A.row(3), B.col(2), C(3, 2));
  mydot(A.row(3), B.col(3), C(3, 3));
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
      mymul_4x4(a, b, c);
    }
  }
}
