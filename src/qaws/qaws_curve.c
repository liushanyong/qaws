#include "qaws_curve.h"
#include "qaws_bezier.h"
#include "qaws_hermite.h"
#include "qaws_catmull_rom.h"
#include "qaws_bspline.h"
#include "qaws_nurbs.h"
#include "qaws_trajectory.h"
#include "qaws_yuksel.h"
#include "internal/qaws_internal_types.h"
#include <stdlib.h>
#include <string.h>

void qaws_curve_destroy(qaws_curve *curve)
{
	if (curve == NULL) {
		return;
	}

	if (curve->vtable != NULL && curve->vtable->destroy_impl != NULL) {
		curve->vtable->destroy_impl(curve->impl);
	}

	free(curve->span_boundaries);
	free(curve);
}

/* -------------------------------------------------------------------------- */
/*  Reverse: Bezier                                                           */
/* -------------------------------------------------------------------------- */

static qaws_status reverse_bezier(
	qaws_curve const *curve,
	qaws_curve **out_reversed)
{
	qaws_bezier_impl const *impl =
		(qaws_bezier_impl const *)curve->impl;
	unsigned int dim = (unsigned int)curve->dimension;
	unsigned int n = impl->control_point_count;
	qaws_scalar *rev_cp;
	qaws_bezier_desc desc;
	qaws_status status;
	unsigned int i;

	rev_cp = (qaws_scalar *)malloc(sizeof(qaws_scalar) * (size_t)(n * dim));
	if (!rev_cp)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	for (i = 0; i < n; i++)
		memcpy(&rev_cp[i * dim],
		       &impl->control_points[(n - 1 - i) * dim],
		       sizeof(qaws_scalar) * dim);

	desc.dimension = curve->dimension;
	desc.degree = curve->degree;
	desc.control_points = rev_cp;
	desc.control_point_count = n;

	status = qaws_curve_create_bezier(&desc, out_reversed);
	free(rev_cp);
	return status;
}

/* -------------------------------------------------------------------------- */
/*  Reverse: Hermite                                                          */
/* -------------------------------------------------------------------------- */

static qaws_status reverse_hermite(
	qaws_curve const *curve,
	qaws_curve **out_reversed)
{
	qaws_hermite_impl const *impl =
		(qaws_hermite_impl const *)curve->impl;
	unsigned int dim = (unsigned int)curve->dimension;
	unsigned int n = impl->point_count;
	qaws_scalar *rev_pts;
	qaws_scalar *rev_der;
	qaws_hermite_desc desc;
	qaws_status status;
	unsigned int i, c;

	rev_pts = (qaws_scalar *)malloc(sizeof(qaws_scalar) * (size_t)(n * dim));
	if (!rev_pts)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	rev_der = (qaws_scalar *)malloc(sizeof(qaws_scalar) * (size_t)(n * dim));
	if (!rev_der) {
		free(rev_pts);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	for (i = 0; i < n; i++) {
		memcpy(&rev_pts[i * dim],
		       &impl->points[(n - 1 - i) * dim],
		       sizeof(qaws_scalar) * dim);
		for (c = 0; c < dim; c++)
			rev_der[i * dim + c] =
				-impl->tangents[(n - 1 - i) * dim + c];
	}

	desc.dimension = curve->dimension;
	desc.degree = curve->degree;
	desc.points = rev_pts;
	desc.derivatives = rev_der;
	desc.point_count = n;
	desc.derivative_count = n;

	status = qaws_curve_create_hermite(&desc, out_reversed);
	free(rev_pts);
	free(rev_der);
	return status;
}

/* -------------------------------------------------------------------------- */
/*  Reverse: Catmull-Rom                                                      */
/* -------------------------------------------------------------------------- */

static qaws_status reverse_catmull_rom(
	qaws_curve const *curve,
	qaws_curve **out_reversed)
{
	qaws_catmull_rom_impl const *impl =
		(qaws_catmull_rom_impl const *)curve->impl;
	unsigned int dim = (unsigned int)curve->dimension;
	unsigned int n = impl->control_point_count;
	qaws_scalar *rev_cp;
	qaws_catmull_rom_desc desc;
	qaws_status status;
	unsigned int i;

	rev_cp = (qaws_scalar *)malloc(sizeof(qaws_scalar) * (size_t)(n * dim));
	if (!rev_cp)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	for (i = 0; i < n; i++)
		memcpy(&rev_cp[i * dim],
		       &impl->control_points[(n - 1 - i) * dim],
		       sizeof(qaws_scalar) * dim);

	desc.dimension = curve->dimension;
	desc.control_points = rev_cp;
	desc.control_point_count = n;
	desc.parameterization = impl->parameterization;
	desc.closed = impl->closed;

	status = qaws_curve_create_catmull_rom(&desc, out_reversed);
	free(rev_cp);
	return status;
}

/* -------------------------------------------------------------------------- */
/*  Reverse: B-Spline                                                         */
/* -------------------------------------------------------------------------- */

static qaws_status reverse_bspline(
	qaws_curve const *curve,
	qaws_curve **out_reversed)
{
	qaws_bspline_impl const *impl =
		(qaws_bspline_impl const *)curve->impl;
	unsigned int dim = (unsigned int)curve->dimension;
	unsigned int n = impl->control_point_count;
	unsigned int kn = impl->knot_count;
	qaws_scalar *rev_cp;
	qaws_scalar *rev_knots;
	qaws_bspline_desc desc;
	qaws_status status;
	unsigned int i;
	qaws_scalar knot_min;
	qaws_scalar knot_max;

	rev_cp = (qaws_scalar *)malloc(sizeof(qaws_scalar) * (size_t)(n * dim));
	if (!rev_cp)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	rev_knots = (qaws_scalar *)malloc(sizeof(qaws_scalar) * (size_t)kn);
	if (!rev_knots) {
		free(rev_cp);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	for (i = 0; i < n; i++)
		memcpy(&rev_cp[i * dim],
		       &impl->control_points[(n - 1 - i) * dim],
		       sizeof(qaws_scalar) * dim);

	knot_min = impl->knots[0];
	knot_max = impl->knots[kn - 1];
	for (i = 0; i < kn; i++)
		rev_knots[i] = knot_min + knot_max - impl->knots[kn - 1 - i];

	desc.dimension = curve->dimension;
	desc.degree = curve->degree;
	desc.control_points = rev_cp;
	desc.control_point_count = n;
	desc.knots = rev_knots;
	desc.knot_count = kn;
	desc.is_uniform = 0;
	desc.is_closed = 0;

	status = qaws_curve_create_bspline(&desc, out_reversed);
	free(rev_cp);
	free(rev_knots);
	return status;
}

/* -------------------------------------------------------------------------- */
/*  Reverse: NURBS                                                            */
/* -------------------------------------------------------------------------- */

static qaws_status reverse_nurbs(
	qaws_curve const *curve,
	qaws_curve **out_reversed)
{
	qaws_nurbs_impl const *impl =
		(qaws_nurbs_impl const *)curve->impl;
	unsigned int dim = (unsigned int)curve->dimension;
	unsigned int n = impl->control_point_count;
	unsigned int kn = impl->knot_count;
	unsigned int wn = impl->weight_count;
	qaws_scalar *rev_cp;
	qaws_scalar *rev_knots;
	qaws_scalar *rev_weights;
	qaws_nurbs_desc desc;
	qaws_status status;
	unsigned int i;
	qaws_scalar knot_min;
	qaws_scalar knot_max;

	rev_cp = (qaws_scalar *)malloc(sizeof(qaws_scalar) * (size_t)(n * dim));
	if (!rev_cp)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	rev_knots = (qaws_scalar *)malloc(sizeof(qaws_scalar) * (size_t)kn);
	if (!rev_knots) {
		free(rev_cp);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	rev_weights = (qaws_scalar *)malloc(sizeof(qaws_scalar) * (size_t)wn);
	if (!rev_weights) {
		free(rev_knots);
		free(rev_cp);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	for (i = 0; i < n; i++)
		memcpy(&rev_cp[i * dim],
		       &impl->control_points[(n - 1 - i) * dim],
		       sizeof(qaws_scalar) * dim);

	knot_min = impl->knots[0];
	knot_max = impl->knots[kn - 1];
	for (i = 0; i < kn; i++)
		rev_knots[i] = knot_min + knot_max - impl->knots[kn - 1 - i];

	for (i = 0; i < wn; i++)
		rev_weights[i] = impl->weights[wn - 1 - i];

	desc.dimension = curve->dimension;
	desc.degree = curve->degree;
	desc.control_points = rev_cp;
	desc.control_point_count = n;
	desc.knots = rev_knots;
	desc.knot_count = kn;
	desc.weights = rev_weights;
	desc.weight_count = wn;
	desc.is_closed = 0;

	status = qaws_curve_create_nurbs(&desc, out_reversed);
	free(rev_cp);
	free(rev_knots);
	free(rev_weights);
	return status;
}

/* -------------------------------------------------------------------------- */
/*  Reverse: Trajectory                                                       */
/* -------------------------------------------------------------------------- */

static qaws_status reverse_trajectory(
	qaws_curve const *curve,
	qaws_curve **out_reversed)
{
	qaws_trajectory_impl const *impl =
		(qaws_trajectory_impl const *)curve->impl;
	unsigned int dim = (unsigned int)curve->dimension;
	unsigned int n = impl->key_count;
	unsigned int vel_n = impl->key_velocity_count;
	unsigned int acc_n = impl->key_acceleration_count;
	qaws_scalar *rev_pos = NULL;
	qaws_scalar *rev_times = NULL;
	qaws_scalar *rev_vel = NULL;
	qaws_scalar *rev_acc = NULL;
	qaws_trajectory_desc desc;
	qaws_status status;
	unsigned int i, c;
	qaws_scalar time_min;
	qaws_scalar time_max;

	rev_pos = (qaws_scalar *)malloc(sizeof(qaws_scalar) * (size_t)(n * dim));
	if (!rev_pos)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	rev_times = (qaws_scalar *)malloc(sizeof(qaws_scalar) * (size_t)n);
	if (!rev_times) {
		free(rev_pos);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Reverse key positions */
	for (i = 0; i < n; i++)
		memcpy(&rev_pos[i * dim],
		       &impl->key_positions[(n - 1 - i) * dim],
		       sizeof(qaws_scalar) * dim);

	/* Reflect key times */
	time_min = impl->key_times[0];
	time_max = impl->key_times[n - 1];
	for (i = 0; i < n; i++)
		rev_times[i] = time_min + time_max - impl->key_times[n - 1 - i];

	/* Reverse and negate velocities */
	if (vel_n > 0 && impl->key_velocities) {
		rev_vel = (qaws_scalar *)malloc(
			sizeof(qaws_scalar) * (size_t)(vel_n * dim));
		if (!rev_vel) {
			free(rev_times);
			free(rev_pos);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		for (i = 0; i < vel_n; i++)
			for (c = 0; c < dim; c++)
				rev_vel[i * dim + c] =
					-impl->key_velocities[(vel_n - 1 - i) * dim + c];
	}

	/* Reverse and negate accelerations */
	if (acc_n > 0 && impl->key_accelerations) {
		rev_acc = (qaws_scalar *)malloc(
			sizeof(qaws_scalar) * (size_t)(acc_n * dim));
		if (!rev_acc) {
			free(rev_vel);
			free(rev_times);
			free(rev_pos);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		for (i = 0; i < acc_n; i++)
			for (c = 0; c < dim; c++)
				rev_acc[i * dim + c] =
					-impl->key_accelerations[(acc_n - 1 - i) * dim + c];
	}

	desc.dimension = curve->dimension;
	desc.degree = curve->degree;
	desc.key_positions = rev_pos;
	desc.key_count = n;
	desc.key_times = rev_times;
	desc.key_time_count = n;
	desc.key_velocities = rev_vel;
	desc.key_velocity_count = vel_n;
	desc.key_accelerations = rev_acc;
	desc.key_acceleration_count = acc_n;
	desc.closed = impl->closed;

	status = qaws_curve_create_trajectory(&desc, out_reversed);
	free(rev_acc);
	free(rev_vel);
	free(rev_times);
	free(rev_pos);
	return status;
}

/* -------------------------------------------------------------------------- */
/*  Reverse: Yuksel                                                           */
/* -------------------------------------------------------------------------- */

static qaws_status reverse_yuksel(
	qaws_curve const *curve,
	qaws_curve **out_reversed)
{
	qaws_yuksel_impl const *impl =
		(qaws_yuksel_impl const *)curve->impl;
	unsigned int dim = (unsigned int)curve->dimension;
	unsigned int n = impl->control_point_count;
	qaws_scalar *rev_cp;
	qaws_yuksel_desc desc;
	qaws_status status;
	unsigned int i;

	rev_cp = (qaws_scalar *)malloc(sizeof(qaws_scalar) * (size_t)(n * dim));
	if (!rev_cp)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	for (i = 0; i < n; i++)
		memcpy(&rev_cp[i * dim],
		       &impl->control_points[(n - 1 - i) * dim],
		       sizeof(qaws_scalar) * dim);

	desc.dimension = curve->dimension;
	desc.control_points = rev_cp;
	desc.control_point_count = n;
	desc.mode = (qaws_yuksel_mode)impl->mode;
	desc.closed = impl->closed;

	status = qaws_curve_create_yuksel(&desc, out_reversed);
	free(rev_cp);
	return status;
}

/* -------------------------------------------------------------------------- */
/*  qaws_curve_reverse                                                        */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_reverse(
	qaws_curve const *curve,
	qaws_curve **out_reversed)
{
	if (!curve || !out_reversed)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_reversed = NULL;

	switch (curve->kind) {
	case QAWS_CURVE_KIND_BEZIER:
		return reverse_bezier(curve, out_reversed);
	case QAWS_CURVE_KIND_HERMITE:
		return reverse_hermite(curve, out_reversed);
	case QAWS_CURVE_KIND_CATMULL_ROM:
		return reverse_catmull_rom(curve, out_reversed);
	case QAWS_CURVE_KIND_BSPLINE:
		return reverse_bspline(curve, out_reversed);
	case QAWS_CURVE_KIND_NURBS:
		return reverse_nurbs(curve, out_reversed);
	case QAWS_CURVE_KIND_TRAJECTORY:
		return reverse_trajectory(curve, out_reversed);
	case QAWS_CURVE_KIND_YUKSEL:
		return reverse_yuksel(curve, out_reversed);
	case QAWS_CURVE_KIND_RATIONAL_BEZIER:
	case QAWS_CURVE_KIND_COMPOSITE:
	case QAWS_CURVE_KIND_ARC:
	case QAWS_CURVE_KIND_POLYNOMIAL:
	case QAWS_CURVE_KIND_CLOTHOID:
	case QAWS_CURVE_KIND_SUBDIVISION:
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	default:
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	}
}
