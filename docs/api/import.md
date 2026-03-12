# API Reference: qaws_import.h

Deserialization and import from raw point data.

---

## qaws_polyline_import_catmull_rom

```c
qaws_status qaws_polyline_import_catmull_rom(
    qaws_dimension dimension,
    void const* points,
    unsigned int point_count,
    qaws_parameterization parameterization,
    int closed,
    qaws_curve** out_curve);
```

Creates a Catmull-Rom spline that interpolates through a raw point array.

**Parameters:**
- `dimension` -- `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D`.
- `points` -- Flat array of `qaws_vec2` or `qaws_vec3`.
- `point_count` -- Number of points (>= 2).
- `parameterization` -- Knot spacing. Pass 0 for centripetal (default).
- `closed` -- Non-zero for a closed loop.
- `out_curve` -- Receives the created Catmull-Rom curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_polyline_import_trajectory

```c
qaws_status qaws_polyline_import_trajectory(
    qaws_dimension dimension,
    void const* points,
    unsigned int point_count,
    int closed,
    qaws_curve** out_curve);
```

Creates a cubic trajectory that interpolates through a raw point array with uniform timing.

**Parameters:**
- `dimension` -- `QAWS_DIMENSION_2D` or `QAWS_DIMENSION_3D`.
- `points` -- Flat array of `qaws_vec2` or `qaws_vec3`.
- `point_count` -- Number of points (>= 2).
- `closed` -- Non-zero for a closed loop.
- `out_curve` -- Receives the created trajectory curve.

**Returns:** `QAWS_STATUS_OK` on success.
