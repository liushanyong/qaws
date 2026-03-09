# Traversal

Traversal wraps a curve with a motion profile, allowing you to evaluate by distance or time instead of raw parameter.

## Concept

A curve defines *shape*. A traversal defines *how to move along that shape*.

The same curve can be:
- evaluated directly by parameter
- traversed at constant speed
- traversed with acceleration profiles
- traversed by time

## Creating a traversal

```c
qaws_traversal_desc desc = {
    .traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH,
    .motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED,
    .speed = 1.0f,
    .acceleration = 0,
    .max_speed = 0,
    .start_time = 0,
    .end_time = 0,
    .clamp_to_domain = 1
};

qaws_traversal* trav = NULL;
qaws_traversal_create(curve, &desc, &trav);
```

The traversal does **not** own the curve. The curve must remain alive for the traversal's lifetime.

## Traversal modes

| Mode | Input meaning |
|---|---|
| `QAWS_TRAVERSAL_MODE_PARAMETER` | Raw curve parameter |
| `QAWS_TRAVERSAL_MODE_ARC_LENGTH` | Distance along curve from start |
| `QAWS_TRAVERSAL_MODE_TIME` | Elapsed time (uses motion profile) |

## Motion profiles

| Profile | Description |
|---|---|
| `NONE` | No motion mapping, pass-through |
| `CONSTANT_SPEED` | Uniform speed along the curve |
| `CONSTANT_ACCELERATION` | Linearly increasing speed |
| `TRAPEZOIDAL_SPEED` | Accelerate, cruise, decelerate |

## Evaluating a traversal

```c
qaws_eval_result_2d result;
qaws_traversal_evaluate_2d(trav, 0.5f,
    QAWS_EVAL_FLAG_POSITION, &result);
```

The `input_value` is interpreted according to the traversal mode:
- Arc-length mode: `0.5` means 0.5 units of distance along the curve
- Time mode: `0.5` means 0.5 seconds elapsed

## Parameter mapping

Convert between parameter spaces:

```c
qaws_scalar param;

// Time -> curve parameter
qaws_traversal_map_time_to_parameter(trav, 1.0f, &param);

// Distance -> curve parameter
qaws_traversal_map_distance_to_parameter(trav, 5.0f, &param);

// Curve parameter -> distance
qaws_scalar distance;
qaws_traversal_map_parameter_to_distance(trav, 0.5f, &distance);
```

## Lifecycle

```c
qaws_traversal_destroy(trav);
// curve is still valid -- destroy it separately
qaws_curve_destroy(curve);
```

## Example: constant-speed motion along a Bezier

```c
qaws_vec2 pts[] = { {0,0}, {0,2}, {2,2}, {2,0} };
qaws_bezier_desc bdesc = {
    .dimension = QAWS_DIMENSION_2D,
    .degree = 3,
    .control_points = pts,
    .control_point_count = 4
};
qaws_curve* curve = NULL;
qaws_curve_create_bezier(&bdesc, &curve);

qaws_traversal_desc tdesc = {
    .traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH,
    .motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED,
    .speed = 1.0f,
    .clamp_to_domain = 1
};
qaws_traversal* trav = NULL;
qaws_traversal_create(curve, &tdesc, &trav);

// Evaluate at 25%, 50%, 75% of total arc length
qaws_scalar total_len;
qaws_curve_compute_arc_length(curve, 0, 1, &total_len);

for (int i = 1; i <= 3; i++) {
    qaws_eval_result_2d r;
    qaws_traversal_evaluate_2d(trav,
        total_len * (float)i / 4.0f,
        QAWS_EVAL_FLAG_POSITION, &r);
    // r.position is evenly spaced along the curve
}

qaws_traversal_destroy(trav);
qaws_curve_destroy(curve);
```
