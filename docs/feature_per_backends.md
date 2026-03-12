# Feature availability per backend

Qaws detects the active backend at compile time via `QAWS_BACKEND` (see `qaws_platform.h`).
This document lists every cross-backend feature and whether it is available on each target.

## Platform layer

| Feature | C | HLSL | GLSL | Halide |
|---|:---:|:---:|:---:|:---:|
| `qaws_platform.h` (scalar type, math macros, qualifiers) | yes | yes | yes | yes |
| `qaws_core_types.h` (vec2, vec3, eval\_2d, eval\_3d) | yes | yes | yes | yes |
| `QAWS_OUT` / `QAWS_INOUT` output qualifiers | (empty) | `out` / `inout` | `out` / `inout` | (empty) |
| `QAWS_SELECT` branchless select | ternary | ternary | ternary | `Halide::select` |
| Bootstrap header | — | `qaws_hlsl.h` | `qaws_glsl.h` | `qaws_halide.h` |

## Core eval headers (`src/qaws/core/`)

All core functions are available on every backend.

### Bezier — `qaws_decasteljau_core.h`

| Function | C | HLSL | GLSL | Halide |
|---|:---:|:---:|:---:|:---:|
| `qaws_bezier3_eval_2d` / `3d` (cubic, CPs by value) | yes | yes | yes | yes |
| `qaws_decasteljau_2d` / `3d` (arbitrary degree) | yes | yes | yes | yes |
| `qaws_bezier_deriv_points_2d` / `3d` | yes | yes | yes | yes |
| `qaws_decasteljau_eval_2d` / `3d` (position + D1 + D2) | yes | yes | yes | yes |

### Rational Bezier — `qaws_rational_bezier_core.h`

| Function | C | HLSL | GLSL | Halide |
|---|:---:|:---:|:---:|:---:|
| `qaws_rational_bezier_eval_2d` / `3d` | yes | yes | yes | yes |

### Cubic polynomial — `qaws_cubic_poly_core.h`

Shared by Hermite, Catmull-Rom, and Trajectory curves.

| Function | C | HLSL | GLSL | Halide |
|---|:---:|:---:|:---:|:---:|
| `qaws_cubic_eval_2d` / `3d` (position + D1 + D2) | yes | yes | yes | yes |

### Horner polynomial — `qaws_horner_core.h`

| Function | C | HLSL | GLSL | Halide |
|---|:---:|:---:|:---:|:---:|
| `qaws_horner_2d` / `3d` (position only) | yes | yes | yes | yes |
| `qaws_horner_eval_2d` / `3d` (position + D1 + D2) | yes | yes | yes | yes |

### Arc — `qaws_arc_core.h`

| Function | C | HLSL | GLSL | Halide |
|---|:---:|:---:|:---:|:---:|
| `qaws_arc_eval_2d` / `3d` (position only) | yes | yes | yes | yes |
| `qaws_arc_eval_full_2d` / `3d` (position + D1 + D2) | yes | yes | yes | yes |

### B-spline basis — `qaws_bspline_basis_core.h`

| Function | C | HLSL | GLSL | Halide |
|---|:---:|:---:|:---:|:---:|
| `qaws_find_span` (knot span search) | yes (binary search) | yes (binary search) | yes (binary search) | yes (linear scan via `select`) |
| `qaws_bspline_basis` | yes | yes | yes | yes |
| `qaws_bspline_basis_derivs` | yes | yes | yes | yes |

### B-spline eval — `qaws_bspline_eval_core.h`

| Function | C | HLSL | GLSL | Halide |
|---|:---:|:---:|:---:|:---:|
| `qaws_bspline_eval_2d` / `3d` (position + D1 + D2) | yes | yes | yes | yes |

### NURBS eval — `qaws_nurbs_eval_core.h`

| Function | C | HLSL | GLSL | Halide |
|---|:---:|:---:|:---:|:---:|
| `qaws_nurbs_eval_2d` / `3d` (position + D1 + D2) | yes | yes | yes | yes |

### Surface — `qaws_surface_core.h`

| Function | C | HLSL | GLSL | Halide |
|---|:---:|:---:|:---:|:---:|
| `qaws_surface_bezier_eval_3d` | yes | yes | yes | yes |
| `qaws_surface_bspline_eval_3d` | yes | yes | yes | yes |

## C-only layers

| Feature | C | HLSL | GLSL | Halide |
|---|:---:|:---:|:---:|:---:|
| Prepare functions (`qaws_prepare.h` / `.c`) | yes | — | — | — |
| Full C runtime (`qaws.h`, create / eval / sample / traverse / inspect / operations) | yes | — | — | — |
| SIMD batch eval (`batch/*.c` — AVX2, SSE2, NEON) | yes | — | — | — |

## Backend totals

| | C | HLSL | GLSL | Halide |
|---|:---:|:---:|:---:|:---:|
| Core functions available | 22 | 22 | 22 | 22 |
| Prepare functions | 12 | — | — | — |
| Full runtime API | yes | — | — | — |
| SIMD batch eval | yes | — | — | — |

## Halide backend notes

All core eval functions compile on the Halide backend. Control flow in the
B-spline basis functions is entirely on plain C++ `int` values (degree, span,
swap indices), so loops unroll at pipeline definition time and build a
symbolic expression graph.

`qaws_find_span` has a dedicated Halide implementation that builds a chain of
`Halide::select()` calls, returning a symbolic `Halide::Expr` integer that
evaluates the correct knot span per sample at runtime.

`Halide::select` evaluates both branches, so divisions by a zero denominator
in `qaws_bspline_basis` execute silently at runtime; the result is discarded
by the surrounding select. The IR tree for high-degree splines can be large
but this is a one-time AOT compilation cost.
