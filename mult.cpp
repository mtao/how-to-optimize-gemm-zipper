#include "definitions.hpp"
#include <experimental/simd>
#include <ranges>
#include <zipper/expression/nullary/Constant.hpp>

namespace {
auto mymul_4x4(zipper::concepts::Matrix auto const &A,
               zipper::concepts::Matrix auto const &B) {

  // in case we ever experiment, this version is really dependent on colmajor to
  // work right
  static_assert(std::is_same_v<BMat::layout_type, zipper::storage::col_major>);

  // simd_width (definitions.hpp / meson option) = 2 matches the original's SSE.
  using vec = std::experimental::fixed_size_simd<scalar_type, simd_width>;
  constexpr index_type W = vec::size();
  static_assert(W <= 4 && 4 % W == 0);
  constexpr index_type R = 4 / W;
  static_assert(std::is_same_v<AMat::layout_type, zipper::storage::col_major>);
  const index_type k = B.rows();

  std::array<std::array<vec, R>, 4> cv{};

  for (auto p : std::views::iota(index_type{0}, k)) {
    // we now use a local outer product
    zipper::VectorBase a = A.col(p);
    zipper::VectorBase b = B.row(p);

    std::array<vec, R> av;
    for (index_type r = 0; r < R; ++r) {
      av[r].copy_from(&a(r * W), std::experimental::element_aligned);
    }

// Unroll explicitly: cv[t][r] must use constant indices to stay in
// registers, and -O2 (unlike -O3) won't fully unroll this early enough.
#pragma GCC unroll 4
    for (index_type t : std::views::iota(0, 4)) {
      const vec bv = b(t);
      for (index_type r : std::views::iota(size_t{0}, R)) {
        cv[t][r] += av[r] * bv;
      }
    }
  }
  // we know this is colmajor so safe to copy W at once
  zipper::Matrix<scalar_type, 4, 4, false> C;
// Unroll explicitly: cv[t][r] must use constant indices to stay in
// registers, and -O2 (unlike -O3) won't fully unroll this early enough.
#pragma GCC unroll 4
  for (index_type t : std::views::iota(0, 4)) {
    for (index_type r : std::views::iota(size_t{0}, R)) {
      cv[t][r].copy_to(&C(r * W, t), std::experimental::element_aligned);
    }
  }

  return C;
}

void InnerKernel(zipper::concepts::Matrix auto const &A,
                 zipper::concepts::Matrix auto const &B,
                 zipper::concepts::Matrix auto &C) {
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

} // namespace

void MULT_NAME(AMat const &A, BMat const &B, CMat &C) {
  constexpr static index_type mc = 256;
  constexpr static index_type kc = 128;

  const index_type m = C.rows();
  const index_type k = A.cols();

  constexpr auto blocks = [](index_type n, index_type b) {
    constexpr auto block = [](auto r) {
      return std::make_pair(r.front(),
                            static_cast<index_type>(std::ranges::size(r)));
    };
    return std::views::iota(index_type{0}, n) | std::views::chunk(b) |
           std::views::transform(block);
  };

  for (auto [p, pc] : blocks(k, kc)) {
    for (auto [i, ic] : blocks(m, mc)) {

      auto islice = zipper::slice(i, ic);
      auto kslice = zipper::slice(p, pc);

      auto a = A(islice, kslice);
      auto b = B(kslice, zipper::full_extent_t{});
      auto c = C(islice, zipper::full_extent_t{});
      InnerKernel(a, b, c);
    }
  }
}
