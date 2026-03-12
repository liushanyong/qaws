# API Reference: qaws_operations.h

Curve operations (split, join, offset, reparameterization).

Thread-safe: operates on immutable input curves.

---

## qaws_curve_split

```c
qaws_status qaws_curve_split(
    qaws_curve const* curve,
    qaws_scalar parameter,
    qaws_curve** out_left,
    qaws_curve** out_right);
```

Splits a curve at the given parameter into two sub-curves.

**Parameters:**
- `curve` -- The curve to split.
- `parameter` -- Parameter value at which to split. Must be within the curve's domain.
- `out_left` -- Receives the left sub-curve (from domain start to `parameter`).
- `out_right` -- Receives the right sub-curve (from `parameter` to domain end).

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_join

```c
qaws_status qaws_curve_join(
    qaws_curve const* curve_a,
    qaws_curve const* curve_b,
    qaws_curve** out_joined);
```

Joins two curves end-to-end into a single curve.

**Parameters:**
- `curve_a` -- First curve.
- `curve_b` -- Second curve. Must have the same dimension as `curve_a`.
- `out_joined` -- Receives the joined curve.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_offset_2d

```c
qaws_status qaws_curve_offset_2d(
    qaws_curve const* curve,
    qaws_scalar distance,
    int trim,
    qaws_curve** out_curves,
    unsigned int curve_capacity,
    unsigned int* out_count);
```

Computes the offset (parallel) curve at a signed distance. 2D only.

Positive distance offsets to the left of the travel direction. At cusps where the curvature radius is less than the absolute distance, backward loops are trimmed at their self-intersection points to produce a clean boundary.

**Parameters:**
- `curve` -- The 2D source curve.
- `distance` -- Signed offset distance. Positive = left, negative = right.
- `trim` -- Trim level controlling cleanup:

| Trim | Description |
|---|---|
| 0 | Raw offset -- cusp loop removal only. |
| 1 | + Distance-based trim: removes points closer than `|distance|` to the source curve. |
| 2 | + Self-intersection cleanup: detects remaining self-intersections, extracts each loop as a separate closed curve, and discards tails. Produces clean, non-self-intersecting closed loops. |

- `out_curves` -- Caller-allocated array to receive output curve pointers.
- `curve_capacity` -- Maximum number of curves that `out_curves` can hold.
- `out_count` -- Receives the number of curves actually written.

**Returns:** `QAWS_STATUS_OK` on success.

---

## qaws_curve_reparameterize_arc_length

```c
qaws_status qaws_curve_reparameterize_arc_length(
    qaws_curve const* curve,
    unsigned int table_resolution,
    qaws_curve** out_curve);
```

Creates an arc-length reparameterized wrapper curve.

The resulting curve has parameter domain [0, total_arc_length] and maps uniformly to arc length along the source curve.

**Parameters:**
- `curve` -- The source curve. Must outlive the returned wrapper (non-owning reference).
- `table_resolution` -- Size of the arc-length lookup table. Pass 0 for the default (256).
- `out_curve` -- Receives the reparameterized wrapper curve.

**Returns:** `QAWS_STATUS_OK` on success.

**Notes:**
- The source curve must remain valid for the lifetime of the returned wrapper curve.
- The wrapper does not take ownership of the source curve.
