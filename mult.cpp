#include "definitions.hpp"
#include <ranges>
#include <zipper/expression/nullary/Constant.hpp>

namespace {
auto mydot_1x4(zipper::concepts::Vector auto const &x,
               zipper::concepts::Matrix auto const &B) {

  zipper::Vector<scalar_type, 4> r =
      zipper::expression::nullary::Constant<scalar_type>(0);
  const index_type k = x.size();
  const auto indices = std::views::iota(index_type{0}, k);

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
    const scalar_type v = x(p);
    r(0) += v * *(bptrs[0]++);
    r(1) += v * *(bptrs[1]++);
    r(2) += v * *(bptrs[2]++);
    r(3) += v * *(bptrs[3]++);
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
