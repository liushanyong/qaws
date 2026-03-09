# API Reference: qaws_eval.h

Curve evaluation functions.

---

## qaws_curve_evaluate_2d

```c
qaws_status qaws_curve_evaluate_2d(
    qaws_curve const* curve,
    qaws_scalar parameter,
    unsigned int eval_flags,
    qaws_eval_result_2d* out_result);
```

Evaluates a 2D curve at a global parameter.

**Parameters:**
- `curve` -- The curve. Must have dimension 2D.
- `parameter` -- Parameter value within the curve's domain.
- `eval_flags` -- Bitmask of `QAWS_EVAL_FLAG_*` specifying what to compute.
- `out_result` -- Filled with position and/or derivatives. Zeroed before filling.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_evaluate_3d

```c
qaws_status qaws_curve_evaluate_3d(
    qaws_curve const* curve,
    qaws_scalar parameter,
    unsigned int eval_flags,
    qaws_eval_result_3d* out_result);
```

Same as `evaluate_2d` but for 3D curves.

---

## qaws_curve_evaluate_span_2d

```c
qaws_status qaws_curve_evaluate_span_2d(
    qaws_curve const* curve,
    unsigned int span_index,
    qaws_scalar local_parameter,
    unsigned int eval_flags,
    qaws_eval_result_2d* out_result);
```

Evaluates a 2D curve at a specific span with a local parameter.

**Parameters:**
- `span_index` -- 0-based span index (0 to span_count-1).
- `local_parameter` -- Parameter within [0, 1] for that span.
- Other parameters same as `evaluate_2d`.

---

## qaws_curve_evaluate_span_3d

```c
qaws_status qaws_curve_evaluate_span_3d(
    qaws_curve const* curve,
    unsigned int span_index,
    qaws_scalar local_parameter,
    unsigned int eval_flags,
    qaws_eval_result_3d* out_result);
```

Same as `evaluate_span_2d` but for 3D curves.
