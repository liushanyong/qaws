# Qaws Documentation

**Qaws** (قوس) -- Dependency-free C11 curve evaluation and traversal library for 2D and 3D.

## Feature documentation

- [Overview](overview.md) -- What Qaws is, design principles, quick start
- [Curve Families](curve_families.md) -- All 7 curve types with examples
- [Evaluation](evaluation.md) -- Position and derivative evaluation
- [Sampling](sampling.md) -- Generating point arrays from curves
- [Traversal](traversal.md) -- Motion profiles, constant-speed, arc-length
- [Inspection](inspection.md) -- Querying curve properties and analysis
- [Error Handling](error_handling.md) -- Status codes and error patterns
- [Building](building.md) -- CMake and Sharpmake build instructions

## API reference

Per-header documentation with all function signatures and parameters:

- [qaws_types.h](api/types.md) -- Scalar, vector, enum, and struct definitions
- [qaws_status.h](api/status.md) -- Error codes
- [qaws_curve.h](api/curve.md) -- Curve lifecycle (destroy)
- [qaws_bezier.h](api/bezier.md) -- Bezier creation
- [qaws_hermite.h](api/hermite.md) -- Hermite creation
- [qaws_catmull_rom.h](api/catmull_rom.md) -- Catmull-Rom creation
- [qaws_bspline.h](api/bspline.md) -- B-Spline creation
- [qaws_nurbs.h](api/nurbs.md) -- NURBS creation
- [qaws_trajectory.h](api/trajectory.md) -- Trajectory creation
- [qaws_yuksel.h](api/yuksel.md) -- Yuksel C2 spline creation
- [qaws_eval.h](api/eval.md) -- Evaluation functions
- [qaws_sampling.h](api/sampling.md) -- Sampling functions
- [qaws_traversal.h](api/traversal.md) -- Traversal functions
- [qaws_inspect.h](api/inspect.md) -- Inspection and analysis functions
