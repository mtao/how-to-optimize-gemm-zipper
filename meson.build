project(
  'how-to-optimize-gemm-zipper',
  'cpp',
  version: '0.1',
  # Defaults match the original tutorial's makefile (gcc -O2 -msse3, no
  # assertions): optimization=2, no debug info, NDEBUG (which also drops
  # meson's -D_GLIBCXX_ASSERTIONS), and arch_flags=-msse3 (meson_options.txt).
  default_options: [
    'warning_level=3',
    'cpp_std=c++26,vc++latest',
    'buildtype=release',
    'optimization=2',
    'debug=false',
    'b_ndebug=true',
  ],
)

add_project_arguments(get_option('arch_flags'), language: 'cpp')
add_project_arguments(
  '-DSIMD_WIDTH=@0@'.format(get_option('simd_width')),
  language: 'cpp',
)

zipper_dep = dependency('zipper')
argparse_dep = dependency('argparse')
deps = [zipper_dep, argparse_dep]

# The kernel being optimized; each step of the tutorial edits mult.cpp.
executable(
  'mult',
  ['mult.cpp', 'run.cpp'],
  dependencies: deps,
  cpp_args: ['-DMULT_NAME=mult'],
)

# Every step as its own executable, for side-by-side comparison.
solutions = [
  '_',
  '0',
  '1',
  '2',
  '1x4_3',
  '1x4_4',
  '1x4_5',
  '1x4_6',
  '1x4_7',
  '1x4_8',
  '1x4_9',
  '4x4_3',
  '4x4_4',
  '4x4_5',
  '4x4_6',
  '4x4_7',
  '4x4_8',
  '4x4_9',
  '4x4_10',
]
foreach index : solutions
  mult = 'mmult@0@'.format(index)
  executable(
    mult,
    ['@0@.cpp'.format(mult), 'run.cpp'],
    dependencies: deps,
    cpp_args: ['-DMULT_NAME=@0@'.format(mult)],
  )
endforeach
