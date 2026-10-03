# how-to-optimize-gemm-zipper

An implementation of [how-to-optimize-gemm](https://github.com/flame/how-to-optimize-gemm/) using [zipper](https://github.com/mtao/zipper)

## Repository layout

The tutorial is a sequence of small edits to a single kernel, `mult.cpp`.

* **`main`** has the benchmark harness and the starting point: `mult.cpp` is
  the naive triple loop (the original's `MMult0`).
* **`steps-1x4`** makes one commit per step, each a minimal diff to
  `mult.cpp`: `1`, `2`, then the 1x4 track `1x4_3` … `1x4_9`.
* **`steps-4x4`** branches off step `2` (as the original does) for the 4x4
  track `4x4_3` … `4x4_15`, which adds SIMD, cache blocking and packing.
* **`solutions`** branches off `main` and has every step as its own file
  (`mmult0.cpp` … `mmult4x4_15.cpp`, plus `mmult_.cpp`, which is just zipper's
  `C += A * B`), all built side by side for comparison.

Every step is tagged `step/<name>` (`step/0`, `step/1x4_6`, `step/4x4_15`, …),
so you can check one out or see what a step changed:

```sh
git checkout step/1x4_6
git diff step/1x4_5 step/1x4_6     # the change made by step 1x4_6
```

## Building and running

```sh
meson setup build
meson compile -C build
./build/mult                     # default sweep: 40..800 step 40 (square)
./build/mult 160 320 640 800x800x800  # explicit sizes; also accepts MxNxK
```

Output is one line per size: `MxNxK GFLOP/s error`. On the `solutions` branch
each step is built as its own executable, e.g. `./build/mmult4x4_15`.

The default build settings match the original's makefile (`gcc -O2 -msse3`,
no assertions). These can be adjusted using meson options (see
`meson_options.txt`).

## Differences from the original

* The original used a leading dimension of 1000 to avoid some cache line
idiosyncracies, whereas this tutorial will just do square matrices.
* Because there's a full linear algebra library backing this experiment here it
  just uses the frobeneous norm rather than a max-abs difference to measure
error.

### Kernel-level differences

* **`4x4_6`–`4x4_9`: early vectorization.** The original tutorial specifies a
ton of register variables, whereas I wanted to just use arrays to store things.
This gave the compiler opportunities to vectorize earlier than the tutorial
wants.
* **`4x4_10`+.** Explicit SIMD uses `std::experimental::fixed_size_simd` with
  `simd_width` doubles (default 2, i.e. the same SSE registers as the
  original's intrinsics) instead of SSE3 intrinsics.
* **Explicit unrolling (`4x4_10`+).** The original is unrolled by hand (16 named
  scalars / 8 named `__m128d`); here the 4x4 accumulator is an array updated in
  a loop over the 4 columns, which only stays in registers if the loop is fully
  unrolled early. GCC does that at `-O3` but not at `-O2`, so the column loops
  carry `#pragma GCC unroll 4` (without it these steps run ~2.5x slower at the
  original's `-O2`). The same applies inside zipper for small static objects
  (`ZIPPER_UNROLL_SMALL` in `for_each_index`).
* **`4x4_12`–`4x4_15`: packing buffers.** The original uses a stack array for
  packed A and a `static` array for packed B. Allocating the buffers per call
  (and zero-filling them) costs 20–30% at small sizes (144–176); keeping them
  across calls (`static thread_local`) and constructing them with
  `zipper::uninitialized` avoids this.
