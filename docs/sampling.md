# Sampling

Sampling generates arrays of points (or full evaluation results) along a curve. This is the primary way to convert curves into renderable polylines.

## Sampling descriptor

```c
qaws_sampling_desc desc;
desc.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;  // or ARC_LENGTH
desc.sample_count = 100;
desc.error_tolerance = 0;       // for adaptive sampling
desc.use_adaptive_sampling = 0; // 0 = uniform, 1 = adaptive
```

## Sampling modes

| Mode | Description |
|---|---|
| `QAWS_TRAVERSAL_MODE_PARAMETER` | Uniform spacing in parameter space |
| `QAWS_TRAVERSAL_MODE_ARC_LENGTH` | Uniform spacing in distance along curve |

Arc-length sampling produces visually even spacing regardless of parameterization speed.

## Sample positions only

```c
qaws_vec2 positions[256];
unsigned int count = 0;

qaws_sampling_desc desc = {
    .traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER,
    .sample_count = 256,
    .error_tolerance = 0,
    .use_adaptive_sampling = 0
};

qaws_curve_sample_2d(curve, &desc, positions, 256, &count);
// positions[0..count-1] now contain sampled points
```

## Sample with derivatives

When you need tangent or curvature information at each sample:

```c
qaws_eval_result_2d samples[256];
unsigned int count = 0;

qaws_curve_sample_eval_2d(curve, &desc, samples, 256, &count);
// samples[i].position, samples[i].d1, etc.
```

## Query sample count before allocation

```c
unsigned int count;
qaws_curve_get_sample_count(curve, &desc, &count);
// allocate buffer of size `count`, then call sample_2d/3d
```

## Adaptive sampling

When `use_adaptive_sampling` is set to 1, the sampling functions use recursive subdivision instead of uniform spacing. Segments where the chord-to-curve midpoint deviation exceeds `error_tolerance` are subdivided (up to 10 levels deep).

```c
qaws_sampling_desc desc = {
    .traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER,
    .sample_count = 16,           // initial uniform segments
    .error_tolerance = 0.01f,     // max chord-to-curve deviation
    .use_adaptive_sampling = 1
};

qaws_vec2 positions[4096];
unsigned int count = 0;

qaws_status s = qaws_curve_sample_2d(
    curve, &desc, positions, 4096, &count);
// count may be larger than sample_count due to subdivision
```

- `sample_count` sets the initial number of uniform segments to subdivide from.
- Tighter `error_tolerance` produces more points in high-curvature regions.
- Returns `QAWS_STATUS_BUFFER_TOO_SMALL` if the output buffer is insufficient.
- Works with both `sample_2d`/`sample_3d` (positions only). Adaptive mode is not used with `sample_eval_*` or arc-length traversal mode.

## 3D variants

All sampling functions have 3D equivalents:

- `qaws_curve_sample_3d`
- `qaws_curve_sample_eval_3d`
