# Qaws Overview

**Qaws** (Arabic: قوس, "arc") is a dependency-free C11 library for curve evaluation and traversal in 2D and 3D.

## What Qaws does

Qaws lets you:

- **Define** curves from control points, knots, weights, or keyframes
- **Evaluate** position and derivatives at any parameter
- **Sample** curves into point buffers for rendering or export
- **Traverse** curves with constant-speed, acceleration, or time-based motion
- **Inspect** curve properties: degree, continuity, arc length, bounds, closest point

All curve families share one evaluation interface. A Bezier curve and a NURBS curve are evaluated the same way.

## Design principles

1. **C first** -- plain structs, enums, typedefs, explicit memory ownership, opaque handles
2. **Geometry and motion are separate** -- curves define shape; traversal objects define motion
3. **One evaluation model** -- position + D1/D2/D3 derivatives, uniform across all families
4. **Piecewise is normal** -- multi-span curves are first-class, not a special case
5. **Public API is stable** -- internal caches, solvers, and basis math are hidden

## Supported curve families

| Family | Header | Continuity | Interpolating | Closed | 3D |
|---|---|---|---|---|---|
| Bezier | `qaws_bezier.h` | C0 | endpoints only | no | yes |
| Hermite | `qaws_hermite.h` | C1 | yes | no | yes |
| Catmull-Rom | `qaws_catmull_rom.h` | C1 | yes | yes | yes |
| B-Spline | `qaws_bspline.h` | up to C(p-1) | no | yes | yes |
| NURBS | `qaws_nurbs.h` | up to C(p-1) | no | yes | yes |
| Trajectory | `qaws_trajectory.h` | C1 | yes | yes | yes |
| Yuksel | `qaws_yuksel.h` | C2 | yes | yes | yes* |

\* Yuksel 3D supports Bezier mode only. Circular/Elliptical/Hybrid modes are 2D-only.

## Configurable scalar type

Qaws supports `float` (f32) or `double` (f64) via the `QAWS_SCALAR_IS_FLOAT` compile option. Default is `float`.

```c
// CMake
cmake .. -DQAWS_SCALAR_IS_FLOAT=ON   // float (default)
cmake .. -DQAWS_SCALAR_IS_FLOAT=OFF  // double
```

## Quick start

```c
#include "qaws.h"

// Create a cubic Bezier
qaws_vec2 points[] = { {0,0}, {0,1}, {1,1}, {1,0} };
qaws_bezier_desc desc = {
    .dimension = QAWS_DIMENSION_2D,
    .degree = 3,
    .control_points = points,
    .control_point_count = 4
};

qaws_curve* curve = NULL;
qaws_curve_create_bezier(&desc, &curve);

// Evaluate at midpoint
qaws_eval_result_2d result;
qaws_curve_evaluate_2d(curve, 0.5f,
    QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &result);

// result.position = {0.5, 0.75}
// result.d1       = first derivative at t=0.5

qaws_curve_destroy(curve);
```

## Three workflows

**A. Geometry-first** -- create a curve, evaluate it directly.

**B. Motion-on-top** -- wrap a curve with a traversal for constant-speed or accelerated motion.

**C. Tooling and inspection** -- query degree, span count, arc length; sample for rendering.

## Building

### CMake

```bash
mkdir build && cd build
cmake ..
cmake --build .
```

### Sharpmake (Visual Studio)

```bash
cd buildsystem
bootstrap.bat
```

## Repository structure

```
src/qaws/           Public headers and implementations
src/qaws/internal/  Internal numerical machinery
tests/              Test suite with SVG visual output
docs/               This documentation
docs/api/           Per-header API reference
buildsystem/        Sharpmake project generation
```
