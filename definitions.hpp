#pragma once
#include <ranges>
#include <string_view>
#include <zipper/Matrix.hpp>

#define STRINGIFY_IMPL(X) #X
#define STRINGIFY(X) STRINGIFY_IMPL(X)
constexpr static std::string_view mult_name = STRINGIFY(MULT_NAME);
#undef STRINGIFY

using scalar_type = double;
using index_type = zipper::index_type;
using AMat = zipper::Matrix<scalar_type, zipper::dynamic_extent,
                            zipper::dynamic_extent, false>;
using BMat = zipper::Matrix<scalar_type, zipper::dynamic_extent,
                            zipper::dynamic_extent, false>;
using CMat = zipper::Matrix<scalar_type, zipper::dynamic_extent,
                            zipper::dynamic_extent, false>;

void MULT_NAME(AMat const &A, BMat const &B, CMat &C);

constexpr static auto test_rows_range =
    std::views::iota(zipper::index_type{40}, zipper::index_type{801});
constexpr static zipper::index_type test_stride = 40;
constexpr static auto test_row_sizes =
    test_rows_range | std::views::stride(test_stride);

constexpr static int num_repetitions = 2;

// Number of doubles per SIMD vector used by the explicit-SIMD kernels
// (4x4_10 onward). Set via the meson option `simd_width`; the default of 2
// matches the original tutorial, which targets SSE3 (one 128-bit register).
#if !defined(SIMD_WIDTH)
#define SIMD_WIDTH 2
#endif
constexpr static zipper::index_type simd_width = SIMD_WIDTH;
static_assert(simd_width >= 1 && simd_width <= 4 && 4 % simd_width == 0,
              "simd_width must divide the 4-row micro-tile");
