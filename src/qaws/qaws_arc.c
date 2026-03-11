#include "qaws_arc.h"
#include "qaws_curve.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_validation.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ---------------------------------------------------------------------------
 * Impl struct
 * ------------------------------------------------------------------------- */

typedef struct qaws_arc_impl
{
	qaws_arc_segment* segments;
	unsigned int segment_count;
} qaws_arc_impl;

/* ---------------------------------------------------------------------------
 * Helpers
 * ------------------------------------------------------------------------- */

static qaws_scalar arc_segment_length(qaws_arc_segment const *seg)
{
	qaws_scalar sweep = seg->angle_end - seg->angle_start;
	return seg->radius * (qaws_scalar)fabs((double)sweep);
}

/* ---------------------------------------------------------------------------
 * Vtable: eval_span_2d
 *
 * Arc-length parameterization. local_t in [0, span_length].
 *   theta(s) = angle_start + s * sign(sweep) / radius
 *   dtheta/ds = sign(sweep) / radius
 *
 *   pos   = center + radius * (cos(theta), sin(theta))
 *   D1    = (-sin(theta), cos(theta)) * sign(sweep)
 *   D2    = (-cos(theta), -sin(theta)) / radius
 *   D3    = (sin(theta), -cos(theta)) * sign(sweep) / radius^2
 * ------------------------------------------------------------------------- */

static qaws_status arc_eval_span_2d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_2d *out_result)
{
	qaws_arc_impl *impl = (qaws_arc_impl *)curve->impl;
	qaws_arc_segment const *seg = &impl->segments[span_index];
	qaws_scalar sweep = seg->angle_end - seg->angle_start;
	qaws_scalar theta;
	qaws_scalar ct, st;

	memset(out_result, 0, sizeof(*out_result));

	/* local_t is in [0,1]; map to angle */
	theta = seg->angle_start + local_t * sweep;
	ct = (qaws_scalar)cos((double)theta);
	st = (qaws_scalar)sin((double)theta);

	/* Position */
	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position.x = seg->center[0] + seg->radius * ct;
		out_result->position.y = seg->center[1] + seg->radius * st;
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	/* D1 = dP/d(local_t) = radius * (-sin, cos) * sweep */
	if (eval_flags & (QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
	{
		out_result->d1.x = seg->radius * (-st) * sweep;
		out_result->d1.y = seg->radius * ct * sweep;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	/* D2 = d2P/d(local_t)^2 = radius * (-cos, -sin) * sweep^2 */
	if (eval_flags & (QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
	{
		qaws_scalar sweep2 = sweep * sweep;
		out_result->d2.x = seg->radius * (-ct) * sweep2;
		out_result->d2.y = seg->radius * (-st) * sweep2;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	/* D3 = d3P/d(local_t)^3 = radius * (sin, -cos) * sweep^3 */
	if (eval_flags & QAWS_EVAL_FLAG_D3)
	{
		qaws_scalar sweep3 = sweep * sweep * sweep;
		out_result->d3.x = seg->radius * st * sweep3;
		out_result->d3.y = seg->radius * (-ct) * sweep3;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Vtable: eval_span_3d
 *
 * Arc-length parameterization in 3D. Arcs lie in the plane spanned by
 * axis_u and axis_v:
 *   pos   = center + radius * (cos(theta) * axis_u + sin(theta) * axis_v)
 *   D1    = (-sin(theta) * axis_u + cos(theta) * axis_v) * sign(sweep)
 *   D2    = (-cos(theta) * axis_u - sin(theta) * axis_v) / radius
 *   D3    = (sin(theta) * axis_u - cos(theta) * axis_v) * sign(sweep) / r^2
 * ------------------------------------------------------------------------- */

static qaws_status arc_eval_span_3d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_3d *out_result)
{
	qaws_arc_impl *impl = (qaws_arc_impl *)curve->impl;
	qaws_arc_segment const *seg = &impl->segments[span_index];
	qaws_scalar sweep = seg->angle_end - seg->angle_start;
	qaws_scalar theta;
	qaws_scalar ct, st;

	memset(out_result, 0, sizeof(*out_result));

	/* local_t is in [0,1]; map to angle */
	theta = seg->angle_start + local_t * sweep;
	ct = (qaws_scalar)cos((double)theta);
	st = (qaws_scalar)sin((double)theta);

	/* Position: center + r * (cos * u + sin * v) */
	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position.x = seg->center[0]
			+ seg->radius * (ct * seg->axis_u[0] + st * seg->axis_v[0]);
		out_result->position.y = seg->center[1]
			+ seg->radius * (ct * seg->axis_u[1] + st * seg->axis_v[1]);
		out_result->position.z = seg->center[2]
			+ seg->radius * (ct * seg->axis_u[2] + st * seg->axis_v[2]);
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	/* D1 = radius * (-sin*u + cos*v) * sweep */
	if (eval_flags & (QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
	{
		out_result->d1.x = seg->radius * (-st * seg->axis_u[0] + ct * seg->axis_v[0]) * sweep;
		out_result->d1.y = seg->radius * (-st * seg->axis_u[1] + ct * seg->axis_v[1]) * sweep;
		out_result->d1.z = seg->radius * (-st * seg->axis_u[2] + ct * seg->axis_v[2]) * sweep;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	/* D2 = radius * (-cos*u - sin*v) * sweep^2 */
	if (eval_flags & (QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3))
	{
		qaws_scalar sweep2 = sweep * sweep;
		out_result->d2.x = seg->radius * (-ct * seg->axis_u[0] - st * seg->axis_v[0]) * sweep2;
		out_result->d2.y = seg->radius * (-ct * seg->axis_u[1] - st * seg->axis_v[1]) * sweep2;
		out_result->d2.z = seg->radius * (-ct * seg->axis_u[2] - st * seg->axis_v[2]) * sweep2;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	/* D3 = radius * (sin*u - cos*v) * sweep^3 */
	if (eval_flags & QAWS_EVAL_FLAG_D3)
	{
		qaws_scalar sweep3 = sweep * sweep * sweep;
		out_result->d3.x = seg->radius * (st * seg->axis_u[0] - ct * seg->axis_v[0]) * sweep3;
		out_result->d3.y = seg->radius * (st * seg->axis_u[1] - ct * seg->axis_v[1]) * sweep3;
		out_result->d3.z = seg->radius * (st * seg->axis_u[2] - ct * seg->axis_v[2]) * sweep3;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Vtable: destroy
 * ------------------------------------------------------------------------- */

static void arc_destroy_impl(void *impl, qaws_allocator const *allocator)
{
	qaws_arc_impl *ai = (qaws_arc_impl *)impl;
	if (ai)
	{
		qaws_internal_dealloc(allocator, ai->segments);
		qaws_internal_dealloc(allocator, ai);
	}
}

/* ---------------------------------------------------------------------------
 * Vtable: property queries
 * ------------------------------------------------------------------------- */

static int arc_is_closed(qaws_curve const *curve)
{
	qaws_arc_impl *impl = (qaws_arc_impl *)curve->impl;
	qaws_arc_segment const *first;
	qaws_arc_segment const *last;
	qaws_scalar end_theta;
	qaws_scalar start_x, start_y, start_z;
	qaws_scalar end_x, end_y, end_z;
	qaws_scalar dx, dy, dz;

	if (impl->segment_count < 1)
		return 0;

	first = &impl->segments[0];
	last = &impl->segments[impl->segment_count - 1];

	/* Compute start point of first segment */
	if (curve->dimension == QAWS_DIMENSION_2D)
	{
		start_x = first->center[0]
			+ first->radius * (qaws_scalar)cos((double)first->angle_start);
		start_y = first->center[1]
			+ first->radius * (qaws_scalar)sin((double)first->angle_start);

		end_theta = last->angle_end;
		end_x = last->center[0]
			+ last->radius * (qaws_scalar)cos((double)end_theta);
		end_y = last->center[1]
			+ last->radius * (qaws_scalar)sin((double)end_theta);

		dx = end_x - start_x;
		dy = end_y - start_y;
		return (dx * dx + dy * dy) < (qaws_scalar)1e-10;
	}
	else
	{
		qaws_scalar ct, st;

		ct = (qaws_scalar)cos((double)first->angle_start);
		st = (qaws_scalar)sin((double)first->angle_start);
		start_x = first->center[0]
			+ first->radius * (ct * first->axis_u[0] + st * first->axis_v[0]);
		start_y = first->center[1]
			+ first->radius * (ct * first->axis_u[1] + st * first->axis_v[1]);
		start_z = first->center[2]
			+ first->radius * (ct * first->axis_u[2] + st * first->axis_v[2]);

		end_theta = last->angle_end;
		ct = (qaws_scalar)cos((double)end_theta);
		st = (qaws_scalar)sin((double)end_theta);
		end_x = last->center[0]
			+ last->radius * (ct * last->axis_u[0] + st * last->axis_v[0]);
		end_y = last->center[1]
			+ last->radius * (ct * last->axis_u[1] + st * last->axis_v[1]);
		end_z = last->center[2]
			+ last->radius * (ct * last->axis_u[2] + st * last->axis_v[2]);

		dx = end_x - start_x;
		dy = end_y - start_y;
		dz = end_z - start_z;
		return (dx * dx + dy * dy + dz * dz) < (qaws_scalar)1e-10;
	}
}

static int arc_is_periodic(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static int arc_is_rational(qaws_curve const *curve)
{
	(void)curve;
	return 1;
}

static qaws_continuity arc_get_continuity(qaws_curve const *curve)
{
	(void)curve;
	return QAWS_CONTINUITY_C0;
}

/* ---------------------------------------------------------------------------
 * Vtable definition
 * ------------------------------------------------------------------------- */

static qaws_curve_vtable const arc_vtable = {
	arc_eval_span_2d,
	arc_eval_span_3d,
	arc_destroy_impl,
	arc_is_closed,
	arc_is_periodic,
	arc_is_rational,
	arc_get_continuity
};

/* ---------------------------------------------------------------------------
 * Creation
 * ------------------------------------------------------------------------- */

qaws_status qaws_curve_create_arc(
	qaws_arc_desc const *desc,
	qaws_curve **out_curve)
{
	unsigned int i;
	qaws_scalar total_length;
	qaws_scalar cumulative;
	qaws_range parameter_range;
	qaws_curve *curve;
	qaws_arc_impl *impl;
	size_t seg_size;

	if (!desc)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;

	/* Validate dimension */
	if (qaws_internal_validate_dimension(desc->dimension) != QAWS_STATUS_OK)
		return QAWS_STATUS_INVALID_DIMENSION;

	/* Need at least one segment */
	if (desc->segment_count < 1)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->segments)
		return QAWS_STATUS_INVALID_ARGUMENT;

	/* Validate segments: radius must be positive, sweep must be non-zero */
	total_length = (qaws_scalar)0;
	for (i = 0; i < desc->segment_count; ++i)
	{
		qaws_arc_segment const *seg = &desc->segments[i];
		qaws_scalar sweep = seg->angle_end - seg->angle_start;
		qaws_scalar len;

		if (seg->radius <= (qaws_scalar)0)
			return QAWS_STATUS_INVALID_ARGUMENT;

		if ((qaws_scalar)fabs((double)sweep) < (qaws_scalar)1e-15)
			return QAWS_STATUS_DEGENERATE_CURVE;

		len = arc_segment_length(seg);
		if (len < (qaws_scalar)1e-15)
			return QAWS_STATUS_DEGENERATE_CURVE;

		total_length += len;
	}

	/* Parameter range: [0, total_arc_length] */
	parameter_range.min_value = (qaws_scalar)0;
	parameter_range.max_value = total_length;

	/* Allocate curve: degree=2 by convention for circular arcs */
	curve = qaws_internal_curve_alloc(
		QAWS_CURVE_KIND_ARC,
		desc->dimension,
		2,
		desc->segment_count,
		parameter_range,
		&arc_vtable);

	if (!curve)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Set span boundaries as cumulative arc lengths */
	cumulative = (qaws_scalar)0;
	curve->span_boundaries[0] = (qaws_scalar)0;
	for (i = 0; i < desc->segment_count; ++i)
	{
		cumulative += arc_segment_length(&desc->segments[i]);
		curve->span_boundaries[i + 1] = cumulative;
	}

	/* Allocate impl */
	impl = (qaws_arc_impl *)calloc(1, sizeof(qaws_arc_impl));
	if (!impl)
	{
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->segment_count = desc->segment_count;

	/* Copy segments */
	seg_size = (size_t)desc->segment_count * sizeof(qaws_arc_segment);
	impl->segments = (qaws_arc_segment *)malloc(seg_size);
	if (!impl->segments)
	{
		free(impl);
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->segments, desc->segments, seg_size);

	curve->impl = impl;
	*out_curve = curve;

	return QAWS_STATUS_OK;
}
