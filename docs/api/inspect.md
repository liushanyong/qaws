# API Reference: qaws_inspect.h

Curve inspection, analysis, and family-specific queries.

---

## Generic inspection

### qaws_curve_get_kind

```c
qaws_curve_kind qaws_curve_get_kind(qaws_curve const* curve);
```

Returns the curve family. Returns `QAWS_CURVE_KIND_INVALID` if curve is NULL.

### qaws_curve_get_dimension

```c
qaws_dimension qaws_curve_get_dimension(qaws_curve const* curve);
```

Returns `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D`.

### qaws_curve_get_degree

```c
unsigned int qaws_curve_get_degree(qaws_curve const* curve);
```

Returns the polynomial degree.

### qaws_curve_get_span_count

```c
unsigned int qaws_curve_get_span_count(qaws_curve const* curve);
```

Returns the number of piecewise spans.

### qaws_curve_get_parameter_range

```c
qaws_range qaws_curve_get_parameter_range(qaws_curve const* curve);
```

Returns the valid parameter domain.

### qaws_curve_is_closed

```c
int qaws_curve_is_closed(qaws_curve const* curve);
```

Returns 1 if the curve forms a closed loop.

### qaws_curve_is_periodic

```c
int qaws_curve_is_periodic(qaws_curve const* curve);
```

Returns 1 if the curve is periodic (wraps around).

### qaws_curve_is_rational

```c
int qaws_curve_is_rational(qaws_curve const* curve);
```

Returns 1 for NURBS curves (weighted control points).

### qaws_curve_get_continuity

```c
qaws_continuity qaws_curve_get_continuity(qaws_curve const* curve);
```

Returns the overall continuity class (C0, C1, C2, C3).

---

## Analysis helpers

### qaws_curve_compute_arc_length

```c
qaws_status qaws_curve_compute_arc_length(
    qaws_curve const* curve,
    qaws_scalar parameter_min,
    qaws_scalar parameter_max,
    qaws_scalar* out_length);
```

Computes the arc length (distance) along the curve between two parameter values. Uses Gauss-Legendre quadrature internally.

### qaws_curve_compute_bounds_2d / 3d

```c
qaws_status qaws_curve_compute_bounds_2d(
    qaws_curve const* curve,
    qaws_vec2* out_min,
    qaws_vec2* out_max);

qaws_status qaws_curve_compute_bounds_3d(
    qaws_curve const* curve,
    qaws_vec3* out_min,
    qaws_vec3* out_max);
```

Computes an axis-aligned bounding box by dense sampling (64 samples per span).

### qaws_curve_find_closest_parameter_2d / 3d

```c
qaws_status qaws_curve_find_closest_parameter_2d(
    qaws_curve const* curve,
    qaws_vec2 point,
    qaws_scalar* out_parameter);

qaws_status qaws_curve_find_closest_parameter_3d(
    qaws_curve const* curve,
    qaws_vec3 point,
    qaws_scalar* out_parameter);
```

Finds the curve parameter closest to a given point. Uses coarse sampling (32 points) followed by Newton iteration refinement.

### qaws_curve_get_span_continuity

```c
qaws_status qaws_curve_get_span_continuity(
    qaws_curve const* curve,
    unsigned int boundary_index,
    qaws_continuity* out_continuity);
```

Queries continuity at a specific span boundary (index 0 to span_count-2).

---

## Derived geometric helpers

### qaws_curve_compute_tangent_2d / 3d

```c
qaws_status qaws_curve_compute_tangent_2d(
    qaws_curve const* curve,
    qaws_scalar parameter,
    qaws_vec2* out_tangent);

qaws_status qaws_curve_compute_tangent_3d(
    qaws_curve const* curve,
    qaws_scalar parameter,
    qaws_vec3* out_tangent);
```

Computes the unit tangent vector at the given parameter. Returns a zero vector if the derivative magnitude is near zero.

### qaws_curve_compute_curvature_2d

```c
qaws_status qaws_curve_compute_curvature_2d(
    qaws_curve const* curve,
    qaws_scalar parameter,
    qaws_scalar* out_curvature);
```

Computes signed curvature: `(d1.x * d2.y - d1.y * d2.x) / |d1|^3`. Positive for counter-clockwise bending, negative for clockwise.

### qaws_curve_compute_curvature_3d

```c
qaws_status qaws_curve_compute_curvature_3d(
    qaws_curve const* curve,
    qaws_scalar parameter,
    qaws_scalar* out_curvature);
```

Computes unsigned curvature: `|d1 x d2| / |d1|^3`.

### qaws_curve_compute_torsion_3d

```c
qaws_status qaws_curve_compute_torsion_3d(
    qaws_curve const* curve,
    qaws_scalar parameter,
    qaws_scalar* out_torsion);
```

Computes torsion (rate of change of the osculating plane): `(d1 x d2) . d3 / |d1 x d2|^2`. Only meaningful for 3D curves.

### qaws_curve_compute_speed

```c
qaws_status qaws_curve_compute_speed(
    qaws_curve const* curve,
    qaws_scalar parameter,
    qaws_scalar* out_speed);
```

Computes the parametric speed `|d1|` (magnitude of the first derivative). Works for both 2D and 3D curves.

---

## Family-specific inspection

### qaws_bezier_get_control_points

```c
qaws_status qaws_bezier_get_control_points(
    qaws_curve const* curve,
    void* out_control_points,
    unsigned int point_capacity,
    unsigned int* out_point_count);
```

Extracts Bezier control points into a caller-provided buffer.

- `out_control_points` -- Buffer of `qaws_vec2` or `qaws_vec3` (matching the curve's dimension).
- `point_capacity` -- Size of the buffer in number of points.
- `out_point_count` -- Receives the actual number of control points.

Returns `QAWS_STATUS_UNSUPPORTED_OPERATION` if the curve is not a Bezier.

### qaws_bspline_get_knots

```c
qaws_status qaws_bspline_get_knots(
    qaws_curve const* curve,
    qaws_scalar* out_knots,
    unsigned int knot_capacity,
    unsigned int* out_knot_count);
```

Extracts the knot vector. Works for both B-Spline and NURBS curves.

Returns `QAWS_STATUS_UNSUPPORTED_OPERATION` for other curve types.

### qaws_nurbs_get_weights

```c
qaws_status qaws_nurbs_get_weights(
    qaws_curve const* curve,
    qaws_scalar* out_weights,
    unsigned int weight_capacity,
    unsigned int* out_weight_count);
```

Extracts NURBS weights.

Returns `QAWS_STATUS_UNSUPPORTED_OPERATION` if not a NURBS curve.
