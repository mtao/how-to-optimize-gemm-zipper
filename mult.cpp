#include "definitions.hpp"
#include <ranges>

void MULT_NAME(AMat const &A, BMat const &B, CMat &C) {
  const index_type m = C.rows();
  const index_type k = A.cols();
  const index_type n = C.cols();
  for (auto j : std::views::iota(index_type{0}, n)) {
    for (auto i : std::views::iota(index_type{0}, m)) {
      for (auto p : std::views::iota(index_type{0}, k)) {
        C(i, j) = C(i, j) + A(i, p) * B(p, j);
      }
    }
  }
}
