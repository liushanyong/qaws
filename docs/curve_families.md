# Curve Families

Qaws supports seven curve families. Each is created through a descriptor struct and a `qaws_curve_create_*` function. After creation, all families share the same evaluation, sampling, traversal, and inspection APIs.

---

## Bezier

**Header:** `qaws_bezier.h`

Polynomial curves defined by control points. The curve interpolates the first and last control points and is influenced by the intermediate ones.

- **Degree:** 1 (linear), 2 (quadratic), 3 (cubic), or higher
- **Spans:** 1
- **Continuity:** C0 (single span, no joints)
- **Closed:** no
- **Dimension:** 2D, 3D

**When to use:** Single-segment curves where you need direct control over shape via a control polygon.

```c
qaws_vec2 pts[] = { {0,0}, {0,1}, {1,1}, {1,0} };
qaws_bezier_desc desc = {
    .dimension = QAWS_DIMENSION_2D,
    .degree = 3,
    .control_points = pts,
    .control_point_count = 4
};
qaws_curve* c = NULL;
qaws_curve_create_bezier(&desc, &c);
```

---

## Hermite

**Header:** `qaws_hermite.h`

Interpolating curves defined by points and tangent vectors at each point. Cubic Hermite is the most common variant.

- **Degree:** 3 (cubic)
- **Spans:** point_count - 1
- **Continuity:** C1 (tangent continuous at joints)
- **Closed:** no
- **Dimension:** 2D, 3D

**When to use:** When you know both positions and tangent directions at each keypoint.

```c
qaws_vec2 points[] = { {0,0}, {1,0} };
qaws_vec2 tangents[] = { {1,1}, {1,-1} };
qaws_hermite_desc desc = {
    .dimension = QAWS_DIMENSION_2D,
    .degree = 3,
    .points = points,
    .derivatives = tangents,
    .point_count = 2,
    .derivative_count = 2
};
```

---

## Catmull-Rom

**Header:** `qaws_catmull_rom.h`

Cubic interpolating splines that pass through control points. The tangent at each point is derived automatically from neighboring points.

- **Degree:** 3 (cubic)
- **Spans:** control_point_count - 3 (open), or control_point_count (closed)
- **Continuity:** C1
- **Closed:** yes (optional)
- **Dimension:** 2D, 3D

**Parameterization choices:**

| Mode | Description | Recommendation |
|---|---|---|
| `UNIFORM` | Equal parameter spacing | Fast but can produce cusps |
| `CHORDAL` | Spacing proportional to chord length | Avoids cusps, can overshoot |
| `CENTRIPETAL` | Spacing proportional to sqrt(chord length) | Best general-purpose choice |

**When to use:** Automatic smooth interpolation through a set of points, no tangent specification needed.

```c
qaws_vec2 pts[] = { {0,0}, {1,2}, {2,0}, {3,2} };
qaws_catmull_rom_desc desc = {
    .dimension = QAWS_DIMENSION_2D,
    .control_points = pts,
    .control_point_count = 4,
    .parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL,
    .closed = 0
};
```

Note: For an open Catmull-Rom, the first and last control points are "ghost" points used to define tangents at the second and second-to-last points. The curve interpolates the interior points only.

---

## B-Spline

**Header:** `qaws_bspline.h`

General-purpose piecewise polynomial curves with local control. Moving one control point only affects a limited portion of the curve.

- **Degree:** arbitrary (commonly 3)
- **Spans:** determined by knot vector and degree
- **Continuity:** up to C(degree-1) at internal knots
- **Closed:** yes (optional)
- **Dimension:** 2D, 3D

**Knot vector:** Can be uniform (auto-generated) or custom. A clamped uniform knot vector makes the curve interpolate the first and last control points.

**When to use:** CAD applications, complex multi-segment curves, situations where local control matters.

```c
qaws_vec2 pts[] = { {0,0}, {1,2}, {3,2}, {4,0} };
qaws_bspline_desc desc = {
    .dimension = QAWS_DIMENSION_2D,
    .degree = 3,
    .control_points = pts,
    .control_point_count = 4,
    .knots = NULL,    // auto-generate clamped uniform
    .knot_count = 0,
    .is_uniform = 1,
    .is_closed = 0
};
```

---

## NURBS

**Header:** `qaws_nurbs.h`

Non-Uniform Rational B-Splines. Extends B-Splines with per-control-point weights, enabling exact representation of circles, ellipses, and conic sections.

- **Degree:** arbitrary (commonly 2 or 3)
- **Spans:** determined by knot vector and degree
- **Continuity:** up to C(degree-1) at internal knots
- **Closed:** yes (optional)
- **Rational:** yes
- **Dimension:** 2D, 3D

**When to use:** Exact circles and conics, CAD interchange, any situation where B-Splines are insufficient.

```c
// Quarter circle
qaws_vec2 pts[] = { {1,0}, {1,1}, {0,1} };
qaws_scalar knots[] = { 0,0,0, 1,1,1 };
qaws_scalar weights[] = { 1, 0.70710678118f, 1 };
qaws_nurbs_desc desc = {
    .dimension = QAWS_DIMENSION_2D,
    .degree = 2,
    .control_points = pts,
    .control_point_count = 3,
    .knots = knots,
    .knot_count = 6,
    .weights = weights,
    .weight_count = 3,
    .is_closed = 0
};
```

---

## Trajectory

**Header:** `qaws_trajectory.h`

Time-keyed interpolating curves with optional velocity and acceleration constraints. Natural for animation and motion control.

- **Degree:** 3 (cubic)
- **Spans:** key_count - 1
- **Continuity:** C1
- **Closed:** yes (optional)
- **Dimension:** 2D, 3D

**When to use:** Animation keyframes, motion paths with timing, robotic trajectories.

```c
qaws_vec2 positions[] = { {0,0}, {5,3}, {10,0} };
qaws_scalar times[] = { 0.0f, 1.0f, 2.0f };
qaws_trajectory_desc desc = {
    .dimension = QAWS_DIMENSION_2D,
    .degree = 3,
    .key_positions = positions,
    .key_count = 3,
    .key_times = times,
    .key_time_count = 3,
    .key_velocities = NULL,
    .key_velocity_count = 0,
    .key_accelerations = NULL,
    .key_acceleration_count = 0,
    .closed = 0
};
```

---

## Yuksel C2 Interpolating Splines

**Header:** `qaws_yuksel.h`

Based on Cem Yuksel's "A Class of C2 Interpolating Splines" (2020). These curves interpolate all control points with C2 continuity, without requiring tangent input.

- **Degree:** 2 (quadratic sub-curves, blended)
- **Spans:** control_point_count - 1 (open), or control_point_count (closed)
- **Continuity:** C2 (achieved through cos^2/sin^2 blending)
- **Closed:** yes (optional)
- **Dimension:** 2D and 3D (Bezier mode); 2D only (Circular/Elliptical/Hybrid)

### Interpolation modes

| Mode | Description | Dimension |
|---|---|---|
| `BEZIER` | Quadratic Bezier sub-curves through point triplets | 2D, 3D |
| `CIRCULAR` | Circular arcs through point triplets | 2D only |
| `ELLIPTICAL` | Elliptical arcs with iterative fitting | 2D only |
| `HYBRID` | Circular when arc angle is small, elliptical otherwise | 2D only |

### How it works

For each control point P[i], a sub-curve is computed from the triplet (P[i-1], P[i], P[i+1]). Adjacent sub-curves are blended using cos^2 and sin^2 weights:

```
P(t) = cos^2(pi*t/2) * subcurve_left(t) + sin^2(pi*t/2) * subcurve_right(t)
```

This blending gives C2 continuity everywhere.

**When to use:** Smooth interpolation with higher continuity than Catmull-Rom, no tangent specification needed.

```c
qaws_vec2 pts[] = { {0,0}, {1,2}, {2,0}, {3,2}, {4,0} };
qaws_yuksel_desc desc = {
    .dimension = QAWS_DIMENSION_2D,
    .control_points = pts,
    .control_point_count = 5,
    .mode = QAWS_YUKSEL_MODE_BEZIER,
    .closed = 0
};
```

### Reference

Cem Yuksel, "A Class of C2 Interpolating Splines", ACM Transactions on Graphics, 2020.
