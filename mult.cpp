#include "definitions.hpp"
#include <ranges>

namespace {
void mymul_4x4(zipper::concepts::Matrix auto const &A,
               zipper::concepts::Matrix auto const &B,
               zipper::concepts::Matrix auto &C) {

  const index_type k = B.rows();
  for (auto p : std::views::iota(index_type{0}, k)) {
    C(0, 0) += A(0, p) * B(p, 0);
  }
  for (auto p : std::views::iota(index_type{0}, k)) {
    C(0, 1) += A(0, p) * B(p, 1);
  }
  for (auto p : std::views::iota(index_type{0}, k)) {
    C(0, 2) += A(0, p) * B(p, 2);
  }
  for (auto p : std::views::iota(index_type{0}, k)) {
    C(0, 3) += A(0, p) * B(p, 3);
  }
  for (auto p : std::views::iota(index_type{0}, k)) {
    C(1, 0) += A(1, p) * B(p, 0);
  }
  for (auto p : std::views::iota(index_type{0}, k)) {
    C(1, 1) += A(1, p) * B(p, 1);
  }
  for (auto p : std::views::iota(index_type{0}, k)) {
    C(1, 2) += A(1, p) * B(p, 2);
  }
  for (auto p : std::views::iota(index_type{0}, k)) {
    C(1, 3) += A(1, p) * B(p, 3);
  }
  for (auto p : std::views::iota(index_type{0}, k)) {
    C(2, 0) += A(2, p) * B(p, 0);
  }
  for (auto p : std::views::iota(index_type{0}, k)) {
    C(2, 1) += A(2, p) * B(p, 1);
  }
  for (auto p : std::views::iota(index_type{0}, k)) {
    C(2, 2) += A(2, p) * B(p, 2);
  }
  for (auto p : std::views::iota(index_type{0}, k)) {
    C(2, 3) += A(2, p) * B(p, 3);
  }
  for (auto p : std::views::iota(index_type{0}, k)) {
    C(3, 0) += A(3, p) * B(p, 0);
  }
  for (auto p : std::views::iota(index_type{0}, k)) {
    C(3, 1) += A(3, p) * B(p, 1);
  }
  for (auto p : std::views::iota(index_type{0}, k)) {
    C(3, 2) += A(3, p) * B(p, 2);
  }
  for (auto p : std::views::iota(index_type{0}, k)) {
    C(3, 3) += A(3, p) * B(p, 3);
  }
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
