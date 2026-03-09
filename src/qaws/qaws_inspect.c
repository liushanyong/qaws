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

/* ------------------------------------------------------------------ */
/*  Derived geometric helpers                                          */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_compute_tangent_2d(
	qaws_curve const *curve,
	qaws_scalar parameter,
	qaws_vec2 *out_tangent)
{
	qaws_eval_result_2d result;
	qaws_status status;
	qaws_scalar len;

	if (!curve || !out_tangent)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	status = qaws_curve_evaluate_2d(
		curve, parameter, QAWS_EVAL_FLAG_D1, &result);
	if (status != QAWS_STATUS_OK)
		return status;

	len = (qaws_scalar)sqrt(
		(double)(result.d1.x * result.d1.x + result.d1.y * result.d1.y));

	if (len < (qaws_scalar)1.0e-12)
	{
		out_tangent->x = (qaws_scalar)0.0;
		out_tangent->y = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	out_tangent->x = result.d1.x / len;
	out_tangent->y = result.d1.y / len;
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_compute_tangent_3d(
	qaws_curve const *curve,
	qaws_scalar parameter,
	qaws_vec3 *out_tangent)
{
	qaws_eval_result_3d result;
	qaws_status status;
	qaws_scalar len;

	if (!curve || !out_tangent)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	status = qaws_curve_evaluate_3d(
		curve, parameter, QAWS_EVAL_FLAG_D1, &result);
	if (status != QAWS_STATUS_OK)
		return status;

	len = (qaws_scalar)sqrt((double)(
		result.d1.x * result.d1.x +
		result.d1.y * result.d1.y +
		result.d1.z * result.d1.z));

	if (len < (qaws_scalar)1.0e-12)
	{
		out_tangent->x = (qaws_scalar)0.0;
		out_tangent->y = (qaws_scalar)0.0;
		out_tangent->z = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	out_tangent->x = result.d1.x / len;
	out_tangent->y = result.d1.y / len;
	out_tangent->z = result.d1.z / len;
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_compute_curvature_2d(
	qaws_curve const *curve,
	qaws_scalar parameter,
	qaws_scalar *out_curvature)
{
	qaws_eval_result_2d result;
	qaws_status status;
	qaws_scalar speed;
	qaws_scalar cross;

	if (!curve || !out_curvature)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	status = qaws_curve_evaluate_2d(
		curve, parameter, QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &result);
	if (status != QAWS_STATUS_OK)
		return status;

	speed = (qaws_scalar)sqrt((double)(
		result.d1.x * result.d1.x + result.d1.y * result.d1.y));

	if (speed < (qaws_scalar)1.0e-12)
	{
		*out_curvature = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	/* Signed curvature: (d1.x * d2.y - d1.y * d2.x) / |d1|^3 */
	cross = result.d1.x * result.d2.y - result.d1.y * result.d2.x;
	*out_curvature = cross / (speed * speed * speed);
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_compute_curvature_3d(
	qaws_curve const *curve,
	qaws_scalar parameter,
	qaws_scalar *out_curvature)
{
	qaws_eval_result_3d result;
	qaws_status status;
	qaws_scalar speed;
	qaws_scalar cx, cy, cz;
	qaws_scalar cross_len;

	if (!curve || !out_curvature)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	status = qaws_curve_evaluate_3d(
		curve, parameter, QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &result);
	if (status != QAWS_STATUS_OK)
		return status;

	speed = (qaws_scalar)sqrt((double)(
		result.d1.x * result.d1.x +
		result.d1.y * result.d1.y +
		result.d1.z * result.d1.z));

	if (speed < (qaws_scalar)1.0e-12)
	{
		*out_curvature = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	/* |d1 x d2| / |d1|^3 */
	cx = result.d1.y * result.d2.z - result.d1.z * result.d2.y;
	cy = result.d1.z * result.d2.x - result.d1.x * result.d2.z;
	cz = result.d1.x * result.d2.y - result.d1.y * result.d2.x;
	cross_len = (qaws_scalar)sqrt((double)(cx * cx + cy * cy + cz * cz));

	*out_curvature = cross_len / (speed * speed * speed);
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_compute_torsion_3d(
	qaws_curve const *curve,
	qaws_scalar parameter,
	qaws_scalar *out_torsion)
{
	qaws_eval_result_3d result;
	qaws_status status;
	qaws_scalar cx, cy, cz;
	qaws_scalar cross_len_sq;
	qaws_scalar dot;

	if (!curve || !out_torsion)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	status = qaws_curve_evaluate_3d(
		curve, parameter,
		QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3,
		&result);
	if (status != QAWS_STATUS_OK)
		return status;

	/* (d1 x d2) . d3 / |d1 x d2|^2 */
	cx = result.d1.y * result.d2.z - result.d1.z * result.d2.y;
	cy = result.d1.z * result.d2.x - result.d1.x * result.d2.z;
	cz = result.d1.x * result.d2.y - result.d1.y * result.d2.x;
	cross_len_sq = cx * cx + cy * cy + cz * cz;

	if (cross_len_sq < (qaws_scalar)1.0e-24)
	{
		*out_torsion = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	dot = cx * result.d3.x + cy * result.d3.y + cz * result.d3.z;
	*out_torsion = dot / cross_len_sq;
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_compute_speed(
	qaws_curve const *curve,
	qaws_scalar parameter,
	qaws_scalar *out_speed)
{
	qaws_status status;

	if (!curve || !out_speed)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension == QAWS_DIMENSION_2D)
	{
		qaws_eval_result_2d result;
		status = qaws_curve_evaluate_2d(
			curve, parameter, QAWS_EVAL_FLAG_D1, &result);
		if (status != QAWS_STATUS_OK)
			return status;

		*out_speed = (qaws_scalar)sqrt((double)(
			result.d1.x * result.d1.x +
			result.d1.y * result.d1.y));
	}
	else
	{
		qaws_eval_result_3d result;
		status = qaws_curve_evaluate_3d(
			curve, parameter, QAWS_EVAL_FLAG_D1, &result);
		if (status != QAWS_STATUS_OK)
			return status;

		*out_speed = (qaws_scalar)sqrt((double)(
			result.d1.x * result.d1.x +
			result.d1.y * result.d1.y +
			result.d1.z * result.d1.z));
	}

	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Frenet frame                                                       */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_compute_normal_2d(
	qaws_curve const *curve,
	qaws_scalar parameter,
	qaws_vec2 *out_normal)
{
	qaws_eval_result_2d result;
	qaws_status status;
	qaws_scalar len;

	if (!curve || !out_normal)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	status = qaws_curve_evaluate_2d(
		curve, parameter, QAWS_EVAL_FLAG_D1, &result);
	if (status != QAWS_STATUS_OK)
		return status;

	len = (qaws_scalar)sqrt(
		(double)(result.d1.x * result.d1.x + result.d1.y * result.d1.y));

	if (len < (qaws_scalar)1.0e-12)
	{
		out_normal->x = (qaws_scalar)0.0;
		out_normal->y = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	/* Rotate normalized tangent 90 degrees CCW: (-ty, tx) */
	out_normal->x = -result.d1.y / len;
	out_normal->y =  result.d1.x / len;
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_compute_frenet_frame_3d(
	qaws_curve const *curve,
	qaws_scalar parameter,
	qaws_vec3 *out_tangent,
	qaws_vec3 *out_normal,
	qaws_vec3 *out_binormal)
{
	qaws_eval_result_3d result;
	qaws_status status;
	qaws_scalar d1_len;
	qaws_scalar tx, ty, tz;
	qaws_scalar bx, by, bz;
	qaws_scalar b_len;
	qaws_scalar ax, ay, az;

	if (!curve || !out_tangent || !out_normal || !out_binormal)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	status = qaws_curve_evaluate_3d(
		curve, parameter,
		QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &result);
	if (status != QAWS_STATUS_OK)
		return status;

	/* T = normalize(D1) */
	d1_len = (qaws_scalar)sqrt((double)(
		result.d1.x * result.d1.x +
		result.d1.y * result.d1.y +
		result.d1.z * result.d1.z));

	if (d1_len < (qaws_scalar)1.0e-12)
	{
		out_tangent->x  = (qaws_scalar)0.0;
		out_tangent->y  = (qaws_scalar)0.0;
		out_tangent->z  = (qaws_scalar)0.0;
		out_normal->x   = (qaws_scalar)0.0;
		out_normal->y   = (qaws_scalar)0.0;
		out_normal->z   = (qaws_scalar)0.0;
		out_binormal->x = (qaws_scalar)0.0;
		out_binormal->y = (qaws_scalar)0.0;
		out_binormal->z = (qaws_scalar)0.0;
		return QAWS_STATUS_OK;
	}

	tx = result.d1.x / d1_len;
	ty = result.d1.y / d1_len;
	tz = result.d1.z / d1_len;

	/* B = normalize(D1 x D2) */
	bx = result.d1.y * result.d2.z - result.d1.z * result.d2.y;
	by = result.d1.z * result.d2.x - result.d1.x * result.d2.z;
	bz = result.d1.x * result.d2.y - result.d1.y * result.d2.x;
	b_len = (qaws_scalar)sqrt((double)(bx * bx + by * by + bz * bz));

	if (b_len < (qaws_scalar)1.0e-12)
	{
		/*
		 * D1 x D2 is near zero (straight line or inflection).
		 * Pick an arbitrary vector not parallel to T, then
		 * use cross products to build a perpendicular normal.
		 */
		ax = (qaws_scalar)0.0;
		ay = (qaws_scalar)0.0;
		az = (qaws_scalar)0.0;

		/* Choose the axis least aligned with T */
		if (tx * tx <= ty * ty && tx * tx <= tz * tz)
			ax = (qaws_scalar)1.0;
		else if (ty * ty <= tz * tz)
			ay = (qaws_scalar)1.0;
		else
			az = (qaws_scalar)1.0;

		/* B = normalize(T x arbitrary) */
		bx = ty * az - tz * ay;
		by = tz * ax - tx * az;
		bz = tx * ay - ty * ax;
		b_len = (qaws_scalar)sqrt((double)(bx * bx + by * by + bz * bz));

		bx /= b_len;
		by /= b_len;
		bz /= b_len;
	}
	else
	{
		bx /= b_len;
		by /= b_len;
		bz /= b_len;
	}

	/* N = B x T */
	out_tangent->x = tx;
	out_tangent->y = ty;
	out_tangent->z = tz;

	out_normal->x = by * tz - bz * ty;
	out_normal->y = bz * tx - bx * tz;
	out_normal->z = bx * ty - by * tx;

	out_binormal->x = bx;
	out_binormal->y = by;
	out_binormal->z = bz;

	return QAWS_STATUS_OK;
}
