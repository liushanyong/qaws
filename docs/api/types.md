# API Reference: qaws_types.h

Core types, enums, and structs used throughout the Qaws library.

---

## Scalar type

```c
#if QAWS_SCALAR_IS_FLOAT
typedef float qaws_scalar;
#else
typedef double qaws_scalar;
#endif
```

Controlled by the `QAWS_SCALAR_IS_FLOAT` compile definition.

---

## Vector types

### qaws_vec2

```c
typedef struct qaws_vec2 {
    qaws_scalar x;
    qaws_scalar y;
} qaws_vec2;
```

### qaws_vec3

```c
typedef struct qaws_vec3 {
    qaws_scalar x;
    qaws_scalar y;
    qaws_scalar z;
} qaws_vec3;
```

---

## Enums

### qaws_dimension

```c
typedef enum qaws_dimension {
    QAWS_DIMENSION_2D = 2,
    QAWS_DIMENSION_3D = 3
} qaws_dimension;
```

### qaws_curve_kind

```c
typedef enum qaws_curve_kind {
    QAWS_CURVE_KIND_INVALID = 0,
    QAWS_CURVE_KIND_BEZIER,
    QAWS_CURVE_KIND_HERMITE,
    QAWS_CURVE_KIND_CATMULL_ROM,
    QAWS_CURVE_KIND_BSPLINE,
    QAWS_CURVE_KIND_NURBS,
    QAWS_CURVE_KIND_TRAJECTORY,
    QAWS_CURVE_KIND_YUKSEL
} qaws_curve_kind;
```

### qaws_parameterization

```c
typedef enum qaws_parameterization {
    QAWS_PARAMETERIZATION_DEFAULT = 0,
    QAWS_PARAMETERIZATION_UNIFORM,
    QAWS_PARAMETERIZATION_CHORDAL,
    QAWS_PARAMETERIZATION_CENTRIPETAL
} qaws_parameterization;
```

Used by Catmull-Rom curves to select knot spacing strategy.

### qaws_continuity

```c
typedef enum qaws_continuity {
    QAWS_CONTINUITY_NONE = 0,
    QAWS_CONTINUITY_C0,
    QAWS_CONTINUITY_C1,
    QAWS_CONTINUITY_C2,
    QAWS_CONTINUITY_C3
} qaws_continuity;
```

### qaws_traversal_mode

```c
typedef enum qaws_traversal_mode {
    QAWS_TRAVERSAL_MODE_PARAMETER = 0,
    QAWS_TRAVERSAL_MODE_ARC_LENGTH,
    QAWS_TRAVERSAL_MODE_TIME
} qaws_traversal_mode;
```

### qaws_motion_profile

```c
typedef enum qaws_motion_profile {
    QAWS_MOTION_PROFILE_NONE = 0,
    QAWS_MOTION_PROFILE_CONSTANT_SPEED,
    QAWS_MOTION_PROFILE_CONSTANT_ACCELERATION,
    QAWS_MOTION_PROFILE_TRAPEZOIDAL_SPEED
} qaws_motion_profile;
```

### qaws_eval_flags

```c
typedef enum qaws_eval_flags {
    QAWS_EVAL_FLAG_NONE     = 0,
    QAWS_EVAL_FLAG_POSITION = 1 << 0,
    QAWS_EVAL_FLAG_D1       = 1 << 1,
    QAWS_EVAL_FLAG_D2       = 1 << 2,
    QAWS_EVAL_FLAG_D3       = 1 << 3
} qaws_eval_flags;
```

Bitmask. Combine with `|`.

---

## Structs

### qaws_range

```c
typedef struct qaws_range {
    qaws_scalar min_value;
    qaws_scalar max_value;
} qaws_range;
```

Represents a parameter interval.

### qaws_eval_result_2d

```c
typedef struct qaws_eval_result_2d {
    qaws_vec2 position;
    qaws_vec2 d1;
    qaws_vec2 d2;
    qaws_vec2 d3;
    unsigned int valid_flags;
} qaws_eval_result_2d;
```

| Field | Description |
|---|---|
| `position` | Curve position at parameter |
| `d1` | First derivative (tangent) |
| `d2` | Second derivative |
| `d3` | Third derivative |
| `valid_flags` | Bitmask of `QAWS_EVAL_FLAG_*` indicating which fields are valid |

### qaws_eval_result_3d

Same layout as 2D but with `qaws_vec3` fields.

### qaws_sampling_desc

```c
typedef struct qaws_sampling_desc {
    qaws_traversal_mode traversal_mode;
    unsigned int sample_count;
    qaws_scalar error_tolerance;
    int use_adaptive_sampling;
} qaws_sampling_desc;
```

| Field | Description |
|---|---|
| `traversal_mode` | How to distribute samples (parameter or arc-length) |
| `sample_count` | Number of samples to generate |
| `error_tolerance` | Tolerance for adaptive sampling |
| `use_adaptive_sampling` | 0 = uniform, 1 = adaptive |

### qaws_traversal_desc

```c
typedef struct qaws_traversal_desc {
    qaws_traversal_mode traversal_mode;
    qaws_motion_profile motion_profile;
    qaws_scalar speed;
    qaws_scalar acceleration;
    qaws_scalar max_speed;
    qaws_scalar start_time;
    qaws_scalar end_time;
    int clamp_to_domain;
} qaws_traversal_desc;
```

| Field | Description |
|---|---|
| `traversal_mode` | Input interpretation (parameter, arc-length, time) |
| `motion_profile` | How speed varies along the curve |
| `speed` | Speed for constant-speed profiles |
| `acceleration` | Acceleration for acceleration-based profiles |
| `max_speed` | Speed cap for trapezoidal profiles |
| `start_time` | Start time for time-based traversal |
| `end_time` | End time for time-based traversal |
| `clamp_to_domain` | 1 = clamp input to valid range, 0 = out-of-range returns error |

---

## Opaque handles

```c
typedef struct qaws_curve qaws_curve;
typedef struct qaws_traversal qaws_traversal;
```

Internal layout is hidden. Interact through the public API only.
