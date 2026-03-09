# Inspection and Analysis

Qaws provides functions to query curve properties and perform geometric analysis.

## Generic inspection

These work on any curve regardless of family:

```c
qaws_curve_kind kind = qaws_curve_get_kind(curve);
qaws_dimension  dim  = qaws_curve_get_dimension(curve);
unsigned int  degree  = qaws_curve_get_degree(curve);
unsigned int  spans   = qaws_curve_get_span_count(curve);
qaws_range    range   = qaws_curve_get_parameter_range(curve);
int           closed  = qaws_curve_is_closed(curve);
int         periodic  = qaws_curve_is_periodic(curve);
int         rational  = qaws_curve_is_rational(curve);
qaws_continuity cont  = qaws_curve_get_continuity(curve);
```

| Property | Description |
|---|---|
| `kind` | Curve family (BEZIER, HERMITE, CATMULL_ROM, BSPLINE, NURBS, TRAJECTORY, YUKSEL) |
| `dimension` | 2D or 3D |
| `degree` | Polynomial degree of the curve |
| `span_count` | Number of piecewise segments |
| `parameter_range` | Valid parameter domain [min, max] |
| `is_closed` | 1 if the curve forms a closed loop |
| `is_periodic` | 1 if the curve wraps periodically |
| `is_rational` | 1 for NURBS (weighted control points) |
| `continuity` | C0, C1, C2, or C3 at span boundaries |

## Arc length

Compute the distance along the curve between two parameters:

```c
qaws_scalar length;
qaws_curve_compute_arc_length(curve,
    range.min_value, range.max_value, &length);
```

Internally uses Gauss-Legendre quadrature for accuracy.

## Bounding box

Compute an axis-aligned bounding box by sampling:

```c
qaws_vec2 bbmin, bbmax;
qaws_curve_compute_bounds_2d(curve, &bbmin, &bbmax);

qaws_vec3 bbmin3, bbmax3;
qaws_curve_compute_bounds_3d(curve, &bbmin3, &bbmax3);
```

## Closest point

Find the curve parameter nearest to a given point:

```c
qaws_scalar closest_t;
qaws_vec2 target = { 1.5f, 0.5f };
qaws_curve_find_closest_parameter_2d(curve, target, &closest_t);

// Evaluate at closest_t to get the actual closest point on the curve
qaws_eval_result_2d result;
qaws_curve_evaluate_2d(curve, closest_t,
    QAWS_EVAL_FLAG_POSITION, &result);
```

Uses coarse sampling followed by Newton refinement.

## Span continuity

Query continuity at a specific span boundary:

```c
qaws_continuity cont;
qaws_curve_get_span_continuity(curve, 0, &cont);
// boundary_index is 0 to span_count-2
```

## Derived geometric helpers

Convenience functions that compute geometric quantities from derivatives:

### Tangent

```c
qaws_vec2 tangent;
qaws_curve_compute_tangent_2d(curve, t, &tangent);
// tangent is a unit vector in the direction of travel
```

### Curvature

```c
qaws_scalar kappa;
qaws_curve_compute_curvature_2d(curve, t, &kappa);
// 2D: signed curvature (positive = CCW)
// 3D: unsigned curvature |d1 x d2| / |d1|^3
```

### Torsion (3D only)

```c
qaws_scalar tau;
qaws_curve_compute_torsion_3d(curve, t, &tau);
// Rate of change of the osculating plane
```

### Speed

```c
qaws_scalar speed;
qaws_curve_compute_speed(curve, t, &speed);
// |d1| -- magnitude of the first derivative
```

## Family-specific inspection

### Bezier control points

```c
qaws_vec2 cps[10];
unsigned int count;
qaws_status s = qaws_bezier_get_control_points(curve,
    cps, 10, &count);
// Returns QAWS_STATUS_UNSUPPORTED_OPERATION if not a Bezier curve
```

### B-Spline / NURBS knot vector

```c
qaws_scalar knots[64];
unsigned int knot_count;
qaws_bspline_get_knots(curve, knots, 64, &knot_count);
// Works for both B-Spline and NURBS curves
```

### NURBS weights

```c
qaws_scalar weights[32];
unsigned int weight_count;
qaws_nurbs_get_weights(curve, weights, 32, &weight_count);
// Returns QAWS_STATUS_UNSUPPORTED_OPERATION if not NURBS
```
