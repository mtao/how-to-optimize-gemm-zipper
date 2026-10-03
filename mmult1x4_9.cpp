#include "definitions.hpp"
#include <ranges>
#include <zipper/expression/nullary/Constant.hpp>

namespace {
auto mydot_1x4(zipper::concepts::Vector auto const &x,
               zipper::concepts::Matrix auto const &B) {

  zipper::Vector<scalar_type, 4> r =
      zipper::expression::nullary::Constant<scalar_type>(0);
  const index_type k = x.size();
  ZIPPER_ASSERT(k % 4 == 0);
  const auto indices =
      std::views::iota(index_type{0}, k) | std::views::stride(4);

  // in case we ever experiment, this version is really dependent on colmajor to
  // work right
  static_assert(std::is_same_v<BMat::layout_type, zipper::storage::col_major>);
  std::array<scalar_type const *, 4> bptrs{{
      &B(0, 0),
      &B(0, 1),
      &B(0, 2),
      &B(0, 3),
  }};
  for (auto p : indices) {
    const scalar_type v0 = x(p);
    r(0) += v0 * *(bptrs[0]);
    r(1) += v0 * *(bptrs[1]);
    r(2) += v0 * *(bptrs[2]);
    r(3) += v0 * *(bptrs[3]);
    const scalar_type v1 = x(p + 1);
    r(0) += v1 * *(bptrs[0] + 1);
    r(1) += v1 * *(bptrs[1] + 1);
    r(2) += v1 * *(bptrs[2] + 1);
    r(3) += v1 * *(bptrs[3] + 1);
    const scalar_type v2 = x(p + 2);
    r(0) += v2 * *(bptrs[0] + 2);
    r(1) += v2 * *(bptrs[1] + 2);
    r(2) += v2 * *(bptrs[2] + 2);
    r(3) += v2 * *(bptrs[3] + 2);
    const scalar_type v3 = x(p + 3);
    r(0) += v3 * *(bptrs[0] + 3);
    r(1) += v3 * *(bptrs[1] + 3);
    r(2) += v3 * *(bptrs[2] + 3);
    r(3) += v3 * *(bptrs[3] + 3);
    bptrs[0] += 4;
    bptrs[1] += 4;
    bptrs[2] += 4;
    bptrs[3] += 4;
  }
  return r;
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
      c += mydot_1x4(A.row(i), b);
    }
  }
}
