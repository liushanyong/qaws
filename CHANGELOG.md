# Changelog

All notable changes to Qaws are documented here.

## [1.0.0] - 2026-03-12

First public release.

### Curve families (13)

- Bezier (linear through arbitrary degree, 2D/3D)
- Hermite (cubic, 2D/3D)
- Catmull-Rom (uniform, chordal, centripetal parameterization, open/closed, 2D/3D)
- B-Spline (arbitrary degree, uniform/custom knots, open/closed, 2D/3D)
- NURBS (arbitrary degree, weighted, open/closed, 2D/3D)
- Trajectory (time-keyed cubic, optional velocity/acceleration, open/closed, 2D/3D)
- Yuksel C2 (Bezier, Circular, Elliptical, Hybrid modes, open/closed, 2D/3D)
- Rational Bezier (weighted De Casteljau, 2D/3D)
- Arc (circular/elliptical, multi-segment, 2D/3D)
- Polynomial (monomial form, Horner evaluation, 2D/3D)
- Clothoid (Euler spiral, linear curvature, 2D)
- Subdivision (Chaikin, Lane-Riesenfeld 3/4, open/closed, 2D/3D)
- Composite (heterogeneous multi-segment, 2D/3D)

### Surface types (5)

- Bezier patch (tensor product)
- B-Spline surface
- NURBS surface
- Swept surface (profile along path)
- Ruled surface (linear blend of two curves)

### Evaluation

- Position, D1, D2, D3 derivatives
- Per-span evaluation
- Batch evaluation (SIMD-accelerated: AVX2, SSE2, NEON)
- Surface evaluation with partial derivatives and normals

### Sampling

- Uniform parameter sampling
- Adaptive (error-tolerance based)
- Curvature-weighted
- Feature-preserving
- Streaming (callback-based)

### Traversal

- Parameter, arc-length, and timed modes
- Motion profiles: constant speed, constant acceleration, trapezoidal, S-curve, custom
- 10 easing functions (quad, cubic, sine -- in/out/in-out)
- Wrap modes: clamp, loop, ping-pong
- Multi-curve traversal

### Inspection

- Arc length, bounding box, closest point
- Tangent, normal, curvature, torsion, speed
- Frenet frame (3D)
- Inflection points, extrema
- Winding number (2D)
- Curvature comb (2D/3D)
- Curve-curve intersection, self-intersection

### Operations

- Split at parameter
- Join two curves
- Offset / parallel curve (2D, with self-intersection cleanup)
- Reverse
- Arc-length reparameterization

### Conversion

- Hermite -> Bezier, Catmull-Rom -> Bezier
- Bezier -> B-Spline, B-Spline -> NURBS
- Degree elevation, degree reduction (Bezier)

### Import / export

- SVG path data export (exact for Bezier/Hermite/Catmull-Rom, approximate for others)
- Polyline export (uniform / curvature-weighted, 2D/3D)
- Polyline import as Catmull-Rom or Trajectory
- B-spline fitting (least-squares)

### Multi-backend

- C (full runtime, SIMD batch eval, float/double)
- HLSL (22 core eval functions, SM 5.0+)
- GLSL (22 core eval functions, Vulkan/OpenGL)
- Halide (22 core eval functions, symbolic IR)
- Automatic backend detection via compiler macros
- Prepare functions for CPU-to-GPU coefficient upload

### Infrastructure

- Zero dependencies (C11 standard library only)
- Immutable, thread-safe curves and surfaces
- Custom allocator support (`_ex` variants)
- Stack-allocated inline curves (no malloc, degree <= 7)
- CMake install rules, pkg-config, CMake config export
- Sharpmake build system support
- MIT license
