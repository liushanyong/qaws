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
    QAWS_CURVE_KIND_YUKSEL,
    QAWS_CURVE_KIND_RATIONAL_BEZIER,
    QAWS_CURVE_KIND_COMPOSITE,
    QAWS_CURVE_KIND_ARC,
    QAWS_CURVE_KIND_POLYNOMIAL,
    QAWS_CURVE_KIND_CLOTHOID,
    QAWS_CURVE_KIND_SUBDIVISION,
    QAWS_CURVE_KIND_REPARAMETERIZED
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
    QAWS_MOTION_PROFILE_TRAPEZOIDAL_SPEED,
    QAWS_MOTION_PROFILE_S_CURVE,
    QAWS_MOTION_PROFILE_CUSTOM
} qaws_motion_profile;
```

| Value | Description |
|---|---|
| `NONE` | No speed control |
| `CONSTANT_SPEED` | Uniform speed along the curve |
| `CONSTANT_ACCELERATION` | Linearly increasing speed |
| `TRAPEZOIDAL_SPEED` | Accelerate, cruise, decelerate |
| `S_CURVE` | Smooth jerk-limited speed ramp |
| `CUSTOM` | User-provided speed callback via `qaws_speed_fn` |

### qaws_easing

```c
typedef enum qaws_easing {
    QAWS_EASING_LINEAR = 0,
    QAWS_EASING_QUAD_IN,
    QAWS_EASING_QUAD_OUT,
    QAWS_EASING_QUAD_IN_OUT,
    QAWS_EASING_CUBIC_IN,
    QAWS_EASING_CUBIC_OUT,
    QAWS_EASING_CUBIC_IN_OUT,
    QAWS_EASING_SINE_IN,
    QAWS_EASING_SINE_OUT,
    QAWS_EASING_SINE_IN_OUT
} qaws_easing;
```

Easing functions applied to traversal input values.

### qaws_wrap_mode

```c
typedef enum qaws_wrap_mode {
    QAWS_WRAP_CLAMP = 0,
    QAWS_WRAP_LOOP,
    QAWS_WRAP_PING_PONG
} qaws_wrap_mode;
```

Controls behavior when traversal input exceeds the curve domain.

| Value | Description |
|---|---|
| `CLAMP` | Clamp to domain boundaries |
| `LOOP` | Wrap around to the start |
| `PING_PONG` | Reverse direction at boundaries |

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
    qaws_easing easing;
    qaws_wrap_mode wrap_mode;
    qaws_scalar jerk;
    qaws_speed_fn custom_speed_fn;
    void* custom_speed_fn_user_data;
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
| `easing` | Easing function applied to input values |
| `wrap_mode` | Behavior when input exceeds domain |
| `jerk` | Maximum jerk for S-curve motion profile |
| `custom_speed_fn` | User-provided speed callback for `QAWS_MOTION_PROFILE_CUSTOM` |
| `custom_speed_fn_user_data` | Opaque pointer passed to the custom speed callback |

### qaws_speed_fn

```c
typedef qaws_scalar (*qaws_speed_fn)(qaws_scalar distance, void* user_data);
```

Custom speed function callback. Returns the desired speed given the distance along the curve. Used with `QAWS_MOTION_PROFILE_CUSTOM`.

### qaws_allocator

```c
typedef struct qaws_allocator {
    void* (*alloc)(unsigned long size, void* user_data);
    void  (*dealloc)(void* ptr, void* user_data);
    void* user_data;
} qaws_allocator;
```

Custom allocator. Pass to `_ex` creation functions to override malloc/free. If `alloc` returns NULL, creation fails with `QAWS_STATUS_ALLOCATION_FAILURE`.

---

## Inline curve

### qaws_curve_inline

```c
#define QAWS_INLINE_STORAGE_BYTES 512

typedef struct qaws_curve_inline {
    char storage[QAWS_INLINE_STORAGE_BYTES];
} qaws_curve_inline;
```

Stack-allocated curve for single-span curves (up to degree 7, 3D). No heap allocation. Use `qaws_curve_init_*_inline()` creation functions.

### qaws_curve_inline_get_curve

```c
qaws_curve* qaws_curve_inline_get_curve(qaws_curve_inline* ic);
```

Returns the usable `qaws_curve` pointer from an initialized inline curve. Do **not** call `qaws_curve_destroy()` on inline curves.

---

## Intersection types

### qaws_intersection_2d

```c
typedef struct qaws_intersection_2d {
    qaws_scalar parameter_a;
    qaws_scalar parameter_b;
    qaws_vec2 position;
} qaws_intersection_2d;
```

Result of a 2D curve-curve or self-intersection query. `parameter_a` and `parameter_b` are the parameter values on the respective curves (or both on the same curve for self-intersections).

### qaws_intersection_3d

```c
typedef struct qaws_intersection_3d {
    qaws_scalar parameter_a;
    qaws_scalar parameter_b;
    qaws_vec3 position;
} qaws_intersection_3d;
```

Same layout as `qaws_intersection_2d` but with `qaws_vec3` position.

---

## Curvature comb types

### qaws_curvature_sample_2d

```c
typedef struct qaws_curvature_sample_2d {
    qaws_vec2 position;
    qaws_scalar curvature;
    qaws_vec2 normal;
} qaws_curvature_sample_2d;
```

A single sample of the curvature comb for a 2D curve. Contains the curve position, the scalar curvature at that point, and the unit normal.

### qaws_curvature_sample_3d

```c
typedef struct qaws_curvature_sample_3d {
    qaws_vec3 position;
    qaws_scalar curvature;
    qaws_vec3 normal;
} qaws_curvature_sample_3d;
```

Same layout as `qaws_curvature_sample_2d` but with `qaws_vec3` fields.

---

## Surface types

Defined in `qaws_surface_types.h`.

### qaws_surface_kind

```c
typedef enum qaws_surface_kind {
    QAWS_SURFACE_KIND_INVALID = 0,
    QAWS_SURFACE_KIND_BEZIER,
    QAWS_SURFACE_KIND_BSPLINE,
    QAWS_SURFACE_KIND_NURBS,
    QAWS_SURFACE_KIND_SWEPT,
    QAWS_SURFACE_KIND_RULED
} qaws_surface_kind;
```

### qaws_surface_eval_flags

```c
typedef enum qaws_surface_eval_flags {
    QAWS_SURFACE_EVAL_NONE     = 0,
    QAWS_SURFACE_EVAL_POSITION = 1 << 0,
    QAWS_SURFACE_EVAL_DU       = 1 << 1,
    QAWS_SURFACE_EVAL_DV       = 1 << 2,
    QAWS_SURFACE_EVAL_NORMAL   = 1 << 3,
    QAWS_SURFACE_EVAL_DUU      = 1 << 4,
    QAWS_SURFACE_EVAL_DUV      = 1 << 5,
    QAWS_SURFACE_EVAL_DVV      = 1 << 6
} qaws_surface_eval_flags;
```

Bitmask. Combine with `|`. Controls which surface derivatives to compute.

### qaws_surface_eval_result

```c
typedef struct qaws_surface_eval_result {
    qaws_vec3 position;
    qaws_vec3 du;
    qaws_vec3 dv;
    qaws_vec3 normal;
    qaws_vec3 duu;
    qaws_vec3 duv;
    qaws_vec3 dvv;
    unsigned int valid_flags;
} qaws_surface_eval_result;
```

| Field | Description |
|---|---|
| `position` | Surface position at (u, v) |
| `du` | Partial derivative dS/du |
| `dv` | Partial derivative dS/dv |
| `normal` | Unit surface normal (du x dv, normalized) |
| `duu` | Second partial d^2S/du^2 |
| `duv` | Mixed partial d^2S/dudv |
| `dvv` | Second partial d^2S/dv^2 |
| `valid_flags` | Bitmask of `QAWS_SURFACE_EVAL_*` indicating which fields are valid |

---

## Opaque handles

```c
typedef struct qaws_curve qaws_curve;
typedef struct qaws_traversal qaws_traversal;
typedef struct qaws_surface qaws_surface;
```

Internal layout is hidden. Interact through the public API only.
