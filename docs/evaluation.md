# Evaluation

Qaws provides a unified evaluation model for all curve families. Every curve can be evaluated for position and up to three orders of derivatives.

## Evaluation flags

Request what you need via a bitmask:

```c
QAWS_EVAL_FLAG_POSITION  // position on the curve
QAWS_EVAL_FLAG_D1        // first derivative (tangent direction)
QAWS_EVAL_FLAG_D2        // second derivative (curvature-related)
QAWS_EVAL_FLAG_D3        // third derivative
```

Combine flags with bitwise OR:

```c
unsigned int flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1;
```

Only requested quantities are computed. The `valid_flags` field in the result tells you what was actually filled.

## Evaluating by global parameter

```c
qaws_eval_result_2d result;
qaws_status s = qaws_curve_evaluate_2d(curve, 0.5f,
    QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &result);

// result.position.x, result.position.y
// result.d1.x, result.d1.y
```

The parameter must be within the curve's domain. Query the domain with:

```c
qaws_range range = qaws_curve_get_parameter_range(curve);
// range.min_value ... range.max_value
```

## Evaluating by span and local parameter

For piecewise curves, you can evaluate a specific span directly:

```c
qaws_eval_result_2d result;
qaws_curve_evaluate_span_2d(curve,
    0,       // span index (0 to span_count-1)
    0.5f,    // local parameter within span [0, 1]
    QAWS_EVAL_FLAG_POSITION, &result);
```

This is useful when you already know which span you're in.

## 3D evaluation

All evaluation functions have 3D variants:

```c
qaws_eval_result_3d result;
qaws_curve_evaluate_3d(curve, 0.5f,
    QAWS_EVAL_FLAG_POSITION, &result);
// result.position.x, .y, .z
```

## Derivative semantics

| Flag | Meaning | Typical use |
|---|---|---|
| `D1` | First derivative dP/dt | Tangent direction, velocity |
| `D2` | Second derivative d2P/dt2 | Curvature computation, acceleration |
| `D3` | Third derivative d3P/dt3 | Jerk, torsion in 3D |

Note: derivatives are with respect to the curve parameter, not arc length. For speed-independent derivatives, use traversal with arc-length mode.

## Computing derived quantities from derivatives

**Tangent vector (unit):**
```c
// normalize(d1)
float len = sqrtf(d1.x*d1.x + d1.y*d1.y);
float tx = d1.x / len, ty = d1.y / len;
```

**Normal vector (2D):**
```c
float nx = -ty, ny = tx;
```

**Curvature (2D):**
```c
// kappa = (d1.x*d2.y - d1.y*d2.x) / (d1.x^2 + d1.y^2)^(3/2)
```

**Speed:**
```c
float speed = sqrtf(d1.x*d1.x + d1.y*d1.y);
```
