#include "qaws_inspect.h"
#include "qaws_eval.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_arc_length.h"
#include "internal/qaws_internal_solver.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ------------------------------------------------------------------ */
/*  Generic inspection functions                                       */
/* ------------------------------------------------------------------ */

qaws_curve_kind qaws_curve_get_kind(qaws_curve const *curve)
{
	if (!curve)
		return QAWS_CURVE_KIND_INVALID;
	return curve->kind;
}

qaws_dimension qaws_curve_get_dimension(qaws_curve const *curve)
{
	if (!curve)
		return (qaws_dimension)0;
	return curve->dimension;
}

unsigned int qaws_curve_get_degree(qaws_curve const *curve)
{
	if (!curve)
		return 0;
	return curve->degree;
}

unsigned int qaws_curve_get_span_count(qaws_curve const *curve)
{
	if (!curve)
		return 0;
	return curve->span_count;
}

qaws_range qaws_curve_get_parameter_range(qaws_curve const *curve)
{
	qaws_range range;
	if (!curve)
	{
		range.min_value = (qaws_scalar)0.0;
		range.max_value = (qaws_scalar)0.0;
		return range;
	}
	return curve->parameter_range;
}

int qaws_curve_is_closed(qaws_curve const *curve)
{
	if (!curve || !curve->vtable || !curve->vtable->is_closed)
		return 0;
	return curve->vtable->is_closed(curve);
}

int qaws_curve_is_periodic(qaws_curve const *curve)
{
	if (!curve || !curve->vtable || !curve->vtable->is_periodic)
		return 0;
	return curve->vtable->is_periodic(curve);
}

int qaws_curve_is_rational(qaws_curve const *curve)
{
	if (!curve || !curve->vtable || !curve->vtable->is_rational)
		return 0;
	return curve->vtable->is_rational(curve);
}

qaws_continuity qaws_curve_get_continuity(qaws_curve const *curve)
{
	if (!curve || !curve->vtable || !curve->vtable->get_continuity)
		return QAWS_CONTINUITY_C0;
	return curve->vtable->get_continuity(curve);
}

/* ------------------------------------------------------------------ */
/*  Analysis helpers                                                    */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_compute_arc_length(
	qaws_curve const *curve,
	qaws_scalar parameter_min,
	qaws_scalar parameter_max,
	qaws_scalar *out_length)
{
	if (!curve || !out_length)
		return QAWS_STATUS_INVALID_ARGUMENT;

	return qaws_internal_arc_length(
		curve, parameter_min, parameter_max, out_length);
}

qaws_status qaws_curve_compute_bounds_2d(
	qaws_curve const *curve,
	qaws_vec2 *out_min,
	qaws_vec2 *out_max)
{
	unsigned int sample_count;
	unsigned int i;
	qaws_scalar t;
	qaws_scalar range_min;
	qaws_scalar range_max;
	qaws_eval_result_2d result;
	qaws_status status;

	if (!curve || !out_min || !out_max)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	sample_count = 64 * curve->span_count;
	if (sample_count < 64)
		sample_count = 64;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	/* Evaluate first point to initialize bounds */
	status = qaws_curve_evaluate_2d(
		curve, range_min, QAWS_EVAL_FLAG_POSITION, &result);
	if (status != QAWS_STATUS_OK)
		return status;

	out_min->x = result.position.x;
	out_min->y = result.position.y;
	out_max->x = result.position.x;
	out_max->y = result.position.y;

	for (i = 1; i < sample_count; ++i)
	{
		t = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)(sample_count - 1);

		status = qaws_curve_evaluate_2d(
			curve, t, QAWS_EVAL_FLAG_POSITION, &result);
		if (status != QAWS_STATUS_OK)
			return status;

		if (result.position.x < out_min->x)
			out_min->x = result.position.x;
		if (result.position.y < out_min->y)
			out_min->y = result.position.y;
		if (result.position.x > out_max->x)
			out_max->x = result.position.x;
		if (result.position.y > out_max->y)
			out_max->y = result.position.y;
	}

	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_compute_bounds_3d(
	qaws_curve const *curve,
	qaws_vec3 *out_min,
	qaws_vec3 *out_max)
{
	unsigned int sample_count;
	unsigned int i;
	qaws_scalar t;
	qaws_scalar range_min;
	qaws_scalar range_max;
	qaws_eval_result_3d result;
	qaws_status status;

	if (!curve || !out_min || !out_max)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	sample_count = 64 * curve->span_count;
	if (sample_count < 64)
		sample_count = 64;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	/* Evaluate first point to initialize bounds */
	status = qaws_curve_evaluate_3d(
		curve, range_min, QAWS_EVAL_FLAG_POSITION, &result);
	if (status != QAWS_STATUS_OK)
		return status;

	out_min->x = result.position.x;
	out_min->y = result.position.y;
	out_min->z = result.position.z;
	out_max->x = result.position.x;
	out_max->y = result.position.y;
	out_max->z = result.position.z;

	for (i = 1; i < sample_count; ++i)
	{
		t = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)(sample_count - 1);

		status = qaws_curve_evaluate_3d(
			curve, t, QAWS_EVAL_FLAG_POSITION, &result);
		if (status != QAWS_STATUS_OK)
			return status;

		if (result.position.x < out_min->x)
			out_min->x = result.position.x;
		if (result.position.y < out_min->y)
			out_min->y = result.position.y;
		if (result.position.z < out_min->z)
			out_min->z = result.position.z;
		if (result.position.x > out_max->x)
			out_max->x = result.position.x;
		if (result.position.y > out_max->y)
			out_max->y = result.position.y;
		if (result.position.z > out_max->z)
			out_max->z = result.position.z;
	}

	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Closest-point helpers                                              */
/* ------------------------------------------------------------------ */

struct closest_ctx_2d
{
	qaws_curve const *curve;
};

static void closest_eval_2d(
	void const *ctx,
	qaws_scalar t,
	qaws_scalar *out_pos,
	qaws_scalar *out_d1,
	unsigned int dim_count)
{
	struct closest_ctx_2d const *c = (struct closest_ctx_2d const *)ctx;
	qaws_eval_result_2d result;

	(void)dim_count;

	qaws_curve_evaluate_2d(c->curve, t,
		QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &result);

	out_pos[0] = result.position.x;
	out_pos[1] = result.position.y;
	out_d1[0] = result.d1.x;
	out_d1[1] = result.d1.y;
}

qaws_status qaws_curve_find_closest_parameter_2d(
	qaws_curve const *curve,
	qaws_vec2 point,
	qaws_scalar *out_parameter)
{
	unsigned int i;
	qaws_scalar t;
	qaws_scalar best_t;
	qaws_scalar best_dist_sq;
	qaws_scalar dx;
	qaws_scalar dy;
	qaws_scalar dist_sq;
	qaws_scalar range_min;
	qaws_scalar range_max;
	qaws_eval_result_2d result;
	qaws_status status;
	struct closest_ctx_2d ctx;
	qaws_scalar target_pos[2];

	if (!curve || !out_parameter)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	/* Initial coarse search: sample at 32 points */
	best_t = range_min;
	best_dist_sq = (qaws_scalar)1.0e30;

	for (i = 0; i < 32; ++i)
	{
		t = range_min + (qaws_scalar)i * (range_max - range_min) / (qaws_scalar)31.0;

		status = qaws_curve_evaluate_2d(
			curve, t, QAWS_EVAL_FLAG_POSITION, &result);
		if (status != QAWS_STATUS_OK)
			return status;

		dx = result.position.x - point.x;
		dy = result.position.y - point.y;
		dist_sq = dx * dx + dy * dy;

		if (dist_sq < best_dist_sq)
		{
			best_dist_sq = dist_sq;
			best_t = t;
		}
	}

	/* Refine with Newton iteration */
	ctx.curve = curve;
	target_pos[0] = point.x;
	target_pos[1] = point.y;

	best_t = qaws_internal_closest_point_newton(
		closest_eval_2d, &ctx, target_pos, 2,
		range_min, range_max, best_t, 32);

	*out_parameter = best_t;
	return QAWS_STATUS_OK;
}

struct closest_ctx_3d
{
	qaws_curve const *curve;
};

static void closest_eval_3d(
	void const *ctx,
	qaws_scalar t,
	qaws_scalar *out_pos,
	qaws_scalar *out_d1,
	unsigned int dim_count)
{
	struct closest_ctx_3d const *c = (struct closest_ctx_3d const *)ctx;
	qaws_eval_result_3d result;

	(void)dim_count;

	qaws_curve_evaluate_3d(c->curve, t,
		QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &result);

	out_pos[0] = result.position.x;
	out_pos[1] = result.position.y;
	out_pos[2] = result.position.z;
	out_d1[0] = result.d1.x;
	out_d1[1] = result.d1.y;
	out_d1[2] = result.d1.z;
}

qaws_status qaws_curve_find_closest_parameter_3d(
	qaws_curve const *curve,
	qaws_vec3 point,
	qaws_scalar *out_parameter)
{
	unsigned int i;
	qaws_scalar t;
	qaws_scalar best_t;
	qaws_scalar best_dist_sq;
	qaws_scalar dx;
	qaws_scalar dy;
	qaws_scalar dz;
	qaws_scalar dist_sq;
	qaws_scalar range_min;
	qaws_scalar range_max;
	qaws_eval_result_3d result;
	qaws_status status;
	struct closest_ctx_3d ctx;
	qaws_scalar target_pos[3];

	if (!curve || !out_parameter)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	/* Initial coarse search: sample at 32 points */
	best_t = range_min;
	best_dist_sq = (qaws_scalar)1.0e30;

	for (i = 0; i < 32; ++i)
	{
		t = range_min + (qaws_scalar)i * (range_max - range_min) / (qaws_scalar)31.0;

		status = qaws_curve_evaluate_3d(
			curve, t, QAWS_EVAL_FLAG_POSITION, &result);
		if (status != QAWS_STATUS_OK)
			return status;

		dx = result.position.x - point.x;
		dy = result.position.y - point.y;
		dz = result.position.z - point.z;
		dist_sq = dx * dx + dy * dy + dz * dz;

		if (dist_sq < best_dist_sq)
		{
			best_dist_sq = dist_sq;
			best_t = t;
		}
	}

	/* Refine with Newton iteration */
	ctx.curve = curve;
	target_pos[0] = point.x;
	target_pos[1] = point.y;
	target_pos[2] = point.z;

	best_t = qaws_internal_closest_point_newton(
		closest_eval_3d, &ctx, target_pos, 3,
		range_min, range_max, best_t, 32);

	*out_parameter = best_t;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Span continuity                                                    */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_get_span_continuity(
	qaws_curve const *curve,
	unsigned int boundary_index,
	qaws_continuity *out_continuity)
{
	if (!curve || !out_continuity)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->span_count == 0 || boundary_index >= curve->span_count - 1)
		return QAWS_STATUS_INVALID_ARGUMENT;

	/*
	 * For v1, return the overall curve continuity.
	 * A more precise implementation would check actual derivatives
	 * at the span boundary.
	 */
	if (!curve->vtable || !curve->vtable->get_continuity)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_continuity = curve->vtable->get_continuity(curve);
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Family-specific inspection: Bezier                                 */
/* ------------------------------------------------------------------ */

qaws_status qaws_bezier_get_control_points(
	qaws_curve const *curve,
	void *out_control_points,
	unsigned int point_capacity,
	unsigned int *out_point_count)
{
	qaws_bezier_impl const *impl;
	unsigned int count;
	unsigned int dim_count;

	if (!curve || !out_control_points || !out_point_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->kind != QAWS_CURVE_KIND_BEZIER)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	impl = (qaws_bezier_impl const *)curve->impl;
	count = impl->control_point_count;
	dim_count = (unsigned int)curve->dimension;

	*out_point_count = count;

	if (point_capacity < count)
		return QAWS_STATUS_BUFFER_TOO_SMALL;

	memcpy(out_control_points, impl->control_points,
		(size_t)count * (size_t)dim_count * sizeof(qaws_scalar));

	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Family-specific inspection: B-Spline knots                         */
/* ------------------------------------------------------------------ */

qaws_status qaws_bspline_get_knots(
	qaws_curve const *curve,
	qaws_scalar *out_knots,
	unsigned int knot_capacity,
	unsigned int *out_knot_count)
{
	if (!curve || !out_knots || !out_knot_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->kind == QAWS_CURVE_KIND_BSPLINE)
	{
		qaws_bspline_impl const *impl =
			(qaws_bspline_impl const *)curve->impl;
		unsigned int knot_count = impl->knot_count;

		*out_knot_count = knot_count;

		if (knot_capacity < knot_count)
			return QAWS_STATUS_BUFFER_TOO_SMALL;

		memcpy(out_knots, impl->knots,
			(size_t)knot_count * sizeof(qaws_scalar));

		return QAWS_STATUS_OK;
	}
	else if (curve->kind == QAWS_CURVE_KIND_NURBS)
	{
		qaws_nurbs_impl const *impl =
			(qaws_nurbs_impl const *)curve->impl;
		unsigned int knot_count = impl->knot_count;

		*out_knot_count = knot_count;

		if (knot_capacity < knot_count)
			return QAWS_STATUS_BUFFER_TOO_SMALL;

		memcpy(out_knots, impl->knots,
			(size_t)knot_count * sizeof(qaws_scalar));

		return QAWS_STATUS_OK;
	}

	return QAWS_STATUS_UNSUPPORTED_OPERATION;
}

/* ------------------------------------------------------------------ */
/*  Family-specific inspection: NURBS weights                          */
/* ------------------------------------------------------------------ */

qaws_status qaws_nurbs_get_weights(
	qaws_curve const *curve,
	qaws_scalar *out_weights,
	unsigned int weight_capacity,
	unsigned int *out_weight_count)
{
	qaws_nurbs_impl const *impl;
	unsigned int weight_count;

	if (!curve || !out_weights || !out_weight_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->kind != QAWS_CURVE_KIND_NURBS)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	impl = (qaws_nurbs_impl const *)curve->impl;
	weight_count = impl->weight_count;

	*out_weight_count = weight_count;

	if (weight_capacity < weight_count)
		return QAWS_STATUS_BUFFER_TOO_SMALL;

	memcpy(out_weights, impl->weights,
		(size_t)weight_count * sizeof(qaws_scalar));

	return QAWS_STATUS_OK;
}
