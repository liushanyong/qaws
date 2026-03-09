# API Reference: qaws_sampling.h

Curve sampling functions.

---

## qaws_curve_get_sample_count

```c
qaws_status qaws_curve_get_sample_count(
    qaws_curve const* curve,
    qaws_sampling_desc const* desc,
    unsigned int* out_sample_count);
```

Queries how many samples will be generated for a given sampling configuration. Use this to allocate a buffer before calling `sample_2d`/`sample_3d`.

---

## qaws_curve_sample_2d

```c
qaws_status qaws_curve_sample_2d(
    qaws_curve const* curve,
    qaws_sampling_desc const* desc,
    qaws_vec2* out_positions,
    unsigned int position_capacity,
    unsigned int* out_position_count);
```

Samples a 2D curve into a caller-provided buffer of positions.

**Parameters:**
- `curve` -- The curve. Must be 2D.
- `desc` -- Sampling configuration.
- `out_positions` -- Buffer to receive sampled points.
- `position_capacity` -- Size of the buffer.
- `out_position_count` -- Receives the number of points written.

**Returns:** `QAWS_STATUS_OK` on success. `BUFFER_TOO_SMALL` if capacity is insufficient.

---

## qaws_curve_sample_3d

```c
qaws_status qaws_curve_sample_3d(
    qaws_curve const* curve,
    qaws_sampling_desc const* desc,
    qaws_vec3* out_positions,
    unsigned int position_capacity,
    unsigned int* out_position_count);
```

Same as `sample_2d` but for 3D curves.

---

## qaws_curve_sample_eval_2d

```c
qaws_status qaws_curve_sample_eval_2d(
    qaws_curve const* curve,
    qaws_sampling_desc const* desc,
    qaws_eval_result_2d* out_samples,
    unsigned int sample_capacity,
    unsigned int* out_sample_count);
```

Samples a 2D curve into a buffer of full evaluation results (position + derivatives).

---

## qaws_curve_sample_eval_3d

```c
qaws_status qaws_curve_sample_eval_3d(
    qaws_curve const* curve,
    qaws_sampling_desc const* desc,
    qaws_eval_result_3d* out_samples,
    unsigned int sample_capacity,
    unsigned int* out_sample_count);
```

Same as `sample_eval_2d` but for 3D curves.
