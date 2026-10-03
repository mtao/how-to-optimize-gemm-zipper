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

} // namespace

void MULT_NAME(AMat const &A, BMat const &B, CMat &C) {
  const index_type m = C.rows();
  const index_type n = C.cols();
  for (auto j : std::views::iota(index_type{0}, n)) {
    for (auto i : std::views::iota(index_type{0}, m)) {
      mydot(A.row(i), B.col(j), C(i, j));
    }
  }
}
