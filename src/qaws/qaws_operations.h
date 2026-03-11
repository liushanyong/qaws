#ifndef QAWS_OPERATIONS_H
#define QAWS_OPERATIONS_H

#include "qaws_types.h"
#include "qaws_status.h"

/* Thread-safe: operates on immutable input curves. */

qaws_status qaws_curve_split(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_curve** out_left,
	qaws_curve** out_right);

qaws_status qaws_curve_join(
	qaws_curve const* curve_a,
	qaws_curve const* curve_b,
	qaws_curve** out_joined);

/*
 * Compute the offset (parallel) curve at signed distance.
 * Positive distance offsets to the left of the travel direction.
 * 2D only.
 *
 * At cusps (where curvature radius < |distance|), the raw offset
 * reverses direction. The backward loops are trimmed at their
 * self-intersection points so the result is a clean boundary.
 *
 * trim levels:
 *   0  Raw offset — cusp loop removal only.
 *   1  + Distance-based trim: removes points closer than |distance|
 *      to the curve (inside the swept boundary).
 *   2  + Self-intersection cleanup: detects remaining polyline
 *      self-intersections, extracts each loop as a separate
 *      closed curve (closed at the crossing point), and
 *      discards the tails. Produces clean, non-self-intersecting
 *      closed loops.
 */
qaws_status qaws_curve_offset_2d(
	qaws_curve const* curve,
	qaws_scalar distance,
	int trim,
	qaws_curve** out_curves,
	unsigned int curve_capacity,
	unsigned int* out_count);

#endif /* QAWS_OPERATIONS_H */
