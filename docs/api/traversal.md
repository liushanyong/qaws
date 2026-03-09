# API Reference: qaws_traversal.h

Traversal creation, evaluation, and parameter mapping.

---

## qaws_traversal_create

```c
qaws_status qaws_traversal_create(
    qaws_curve const* curve,
    qaws_traversal_desc const* desc,
    qaws_traversal** out_traversal);
```

Creates a traversal object that wraps a curve with a motion profile.

**Parameters:**
- `curve` -- The curve to traverse. Must remain alive for the traversal's lifetime.
- `desc` -- Traversal configuration.
- `out_traversal` -- Receives the created traversal.

**Returns:** `QAWS_STATUS_OK` on success.

**Notes:**
- The traversal does **not** own the curve.
- Internally builds an arc-length lookup table for distance-based operations.

---

## qaws_traversal_destroy

```c
void qaws_traversal_destroy(qaws_traversal* traversal);
```

Destroys a traversal object. May be NULL (no-op).

---

## qaws_traversal_evaluate_2d

```c
qaws_status qaws_traversal_evaluate_2d(
    qaws_traversal const* traversal,
    qaws_scalar input_value,
    unsigned int eval_flags,
    qaws_eval_result_2d* out_result);
```

Evaluates the traversal at the given input value.

**Parameters:**
- `input_value` -- Interpreted based on traversal mode: parameter, distance, or time.
- `eval_flags` -- Bitmask of `QAWS_EVAL_FLAG_*`.
- `out_result` -- Receives position and derivatives.

---

## qaws_traversal_evaluate_3d

```c
qaws_status qaws_traversal_evaluate_3d(
    qaws_traversal const* traversal,
    qaws_scalar input_value,
    unsigned int eval_flags,
    qaws_eval_result_3d* out_result);
```

Same as `evaluate_2d` but for 3D curves.

---

## qaws_traversal_map_time_to_parameter

```c
qaws_status qaws_traversal_map_time_to_parameter(
    qaws_traversal const* traversal,
    qaws_scalar time_value,
    qaws_scalar* out_parameter);
```

Converts a time value to a curve parameter using the traversal's motion profile.

---

## qaws_traversal_map_distance_to_parameter

```c
qaws_status qaws_traversal_map_distance_to_parameter(
    qaws_traversal const* traversal,
    qaws_scalar distance_value,
    qaws_scalar* out_parameter);
```

Converts a distance-along-curve value to a curve parameter.

---

## qaws_traversal_map_parameter_to_distance

```c
qaws_status qaws_traversal_map_parameter_to_distance(
    qaws_traversal const* traversal,
    qaws_scalar parameter,
    qaws_scalar* out_distance);
```

Converts a curve parameter to a distance-along-curve value.
