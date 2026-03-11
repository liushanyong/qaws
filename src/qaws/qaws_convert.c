#include "qaws_convert.h"
#include "qaws_bezier.h"
#include "qaws_hermite.h"
#include "qaws_catmull_rom.h"
#include "qaws_bspline.h"
#include "qaws_nurbs.h"
#include "internal/qaws_internal_types.h"
#include <stdlib.h>
#include <string.h>

/* -------------------------------------------------------------------------- */
/*  Hermite to Bezier (per span)                                              */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_convert_hermite_to_bezier(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_curve **out_bezier)
{
	qaws_hermite_impl const *impl;
	unsigned int dim_count;
	qaws_scalar const *P0;
	qaws_scalar const *P1;
	qaws_scalar const *M0;
	qaws_scalar const *M1;
	qaws_scalar *cp;
	size_t cp_size;
	qaws_bezier_desc desc;
	qaws_status status;
	unsigned int d;

	if (!curve || !out_bezier)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_bezier = NULL;

	if (curve->kind != QAWS_CURVE_KIND_HERMITE)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	if (span_index >= curve->span_count)
		return QAWS_STATUS_OUT_OF_RANGE;

	impl = (qaws_hermite_impl const *)curve->impl;
	dim_count = (unsigned int)curve->dimension;

	P0 = impl->points + span_index * dim_count;
	P1 = impl->points + (span_index + 1) * dim_count;
	M0 = impl->tangents + span_index * dim_count;
	M1 = impl->tangents + (span_index + 1) * dim_count;

	cp_size = (size_t)(4 * dim_count) * sizeof(qaws_scalar);
	cp = (qaws_scalar *)malloc(cp_size);
	if (!cp)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	for (d = 0; d < dim_count; d++) {
		cp[0 * dim_count + d] = P0[d];
		cp[1 * dim_count + d] = P0[d] + M0[d] / (qaws_scalar)3.0;
		cp[2 * dim_count + d] = P1[d] - M1[d] / (qaws_scalar)3.0;
		cp[3 * dim_count + d] = P1[d];
	}

	memset(&desc, 0, sizeof(desc));
	desc.dimension = curve->dimension;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	status = qaws_curve_create_bezier(&desc, out_bezier);
	free(cp);
	return status;
}

/* -------------------------------------------------------------------------- */
/*  Catmull-Rom to Bezier (per span)                                          */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_convert_catmull_rom_to_bezier(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_curve **out_bezier)
{
	qaws_catmull_rom_impl const *impl;
	unsigned int dim_count;
	qaws_scalar const *coeffs;
	qaws_scalar *cp;
	size_t cp_size;
	qaws_bezier_desc desc;
	qaws_status status;
	unsigned int d;

	if (!curve || !out_bezier)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_bezier = NULL;

	if (curve->kind != QAWS_CURVE_KIND_CATMULL_ROM)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	if (span_index >= curve->span_count)
		return QAWS_STATUS_OUT_OF_RANGE;

	impl = (qaws_catmull_rom_impl const *)curve->impl;
	dim_count = (unsigned int)curve->dimension;

	cp_size = (size_t)(4 * dim_count) * sizeof(qaws_scalar);
	cp = (qaws_scalar *)malloc(cp_size);
	if (!cp)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/*
	 * Coefficients are stored per span as dim_count groups of 4 scalars.
	 * For dimension d: coeffs[d*4+0]=a, [d*4+1]=b, [d*4+2]=c, [d*4+3]=d_coeff
	 *
	 * Polynomial: P(t) = a*t^3 + b*t^2 + c*t + d_coeff  over [0,1]
	 *
	 * Hermite-to-Bezier mapping:
	 *   B0 = d_coeff
	 *   B1 = d_coeff + c/3
	 *   B2 = d_coeff + 2c/3 + b/3
	 *   B3 = a + b + c + d_coeff
	 */
	coeffs = impl->segment_coeffs + span_index * dim_count * 4;

	for (d = 0; d < dim_count; d++) {
		qaws_scalar a = coeffs[d * 4 + 0];
		qaws_scalar b = coeffs[d * 4 + 1];
		qaws_scalar c = coeffs[d * 4 + 2];
		qaws_scalar d_coeff = coeffs[d * 4 + 3];

		cp[0 * dim_count + d] = d_coeff;
		cp[1 * dim_count + d] = d_coeff + c / (qaws_scalar)3.0;
		cp[2 * dim_count + d] = d_coeff + (qaws_scalar)2.0 * c / (qaws_scalar)3.0
		                      + b / (qaws_scalar)3.0;
		cp[3 * dim_count + d] = a + b + c + d_coeff;
	}

	memset(&desc, 0, sizeof(desc));
	desc.dimension = curve->dimension;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	status = qaws_curve_create_bezier(&desc, out_bezier);
	free(cp);
	return status;
}

/* -------------------------------------------------------------------------- */
/*  Bezier to B-Spline                                                        */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_convert_bezier_to_bspline(
	qaws_curve const *curve,
	qaws_curve **out_bspline)
{
	qaws_bezier_impl const *impl;
	unsigned int degree;
	unsigned int dim_count;
	unsigned int knot_count;
	qaws_scalar *knots;
	qaws_bspline_desc desc;
	qaws_status status;
	unsigned int i;

	if (!curve || !out_bspline)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_bspline = NULL;

	if (curve->kind != QAWS_CURVE_KIND_BEZIER)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	impl = (qaws_bezier_impl const *)curve->impl;
	degree = curve->degree;
	dim_count = (unsigned int)curve->dimension;

	/* Clamped knot vector: [0,...,0, 1,...,1] with (degree+1) repeats each */
	knot_count = 2 * (degree + 1);
	knots = (qaws_scalar *)malloc((size_t)knot_count * sizeof(qaws_scalar));
	if (!knots)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	for (i = 0; i <= degree; i++)
		knots[i] = (qaws_scalar)0.0;
	for (i = degree + 1; i < knot_count; i++)
		knots[i] = (qaws_scalar)1.0;

	memset(&desc, 0, sizeof(desc));
	desc.dimension = curve->dimension;
	desc.degree = degree;
	desc.control_points = impl->control_points;
	desc.control_point_count = impl->control_point_count;
	desc.knots = knots;
	desc.knot_count = knot_count;
	desc.is_uniform = 0;
	desc.is_closed = 0;

	status = qaws_curve_create_bspline(&desc, out_bspline);
	free(knots);
	return status;
}

/* -------------------------------------------------------------------------- */
/*  B-Spline to NURBS                                                         */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_convert_bspline_to_nurbs(
	qaws_curve const *curve,
	qaws_curve **out_nurbs)
{
	qaws_bspline_impl const *impl;
	unsigned int cp_count;
	qaws_scalar *weights;
	qaws_nurbs_desc desc;
	qaws_status status;
	unsigned int i;

	if (!curve || !out_nurbs)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_nurbs = NULL;

	if (curve->kind != QAWS_CURVE_KIND_BSPLINE)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	impl = (qaws_bspline_impl const *)curve->impl;
	cp_count = impl->control_point_count;

	/* All weights = 1.0 for a non-rational B-Spline */
	weights = (qaws_scalar *)malloc((size_t)cp_count * sizeof(qaws_scalar));
	if (!weights)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	for (i = 0; i < cp_count; i++)
		weights[i] = (qaws_scalar)1.0;

	memset(&desc, 0, sizeof(desc));
	desc.dimension = curve->dimension;
	desc.degree = curve->degree;
	desc.control_points = impl->control_points;
	desc.control_point_count = cp_count;
	desc.knots = impl->knots;
	desc.knot_count = impl->knot_count;
	desc.weights = weights;
	desc.weight_count = cp_count;
	desc.is_closed = 0;

	status = qaws_curve_create_nurbs(&desc, out_nurbs);
	free(weights);
	return status;
}

/* -------------------------------------------------------------------------- */
/*  Degree elevation (Bezier only)                                            */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_elevate_degree(
	qaws_curve const *curve,
	qaws_curve **out_elevated)
{
	qaws_bezier_impl const *impl;
	unsigned int n;
	unsigned int dim_count;
	unsigned int new_count;
	qaws_scalar *new_cp;
	size_t new_cp_size;
	qaws_bezier_desc desc;
	qaws_status status;
	unsigned int i, d;

	if (!curve || !out_elevated)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_elevated = NULL;

	if (curve->kind != QAWS_CURVE_KIND_BEZIER)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	impl = (qaws_bezier_impl const *)curve->impl;
	n = curve->degree;
	dim_count = (unsigned int)curve->dimension;
	new_count = n + 2;

	new_cp_size = (size_t)(new_count * dim_count) * sizeof(qaws_scalar);
	new_cp = (qaws_scalar *)malloc(new_cp_size);
	if (!new_cp)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Q[0] = P[0] */
	for (d = 0; d < dim_count; d++)
		new_cp[0 * dim_count + d] = impl->control_points[0 * dim_count + d];

	/* Q[i] = (i/(n+1)) * P[i-1] + (1 - i/(n+1)) * P[i]  for i = 1..n */
	for (i = 1; i <= n; i++) {
		qaws_scalar alpha = (qaws_scalar)i / (qaws_scalar)(n + 1);
		qaws_scalar one_minus_alpha = (qaws_scalar)1.0 - alpha;
		for (d = 0; d < dim_count; d++) {
			new_cp[i * dim_count + d] =
				alpha * impl->control_points[(i - 1) * dim_count + d]
				+ one_minus_alpha * impl->control_points[i * dim_count + d];
		}
	}

	/* Q[n+1] = P[n] */
	for (d = 0; d < dim_count; d++)
		new_cp[(n + 1) * dim_count + d] = impl->control_points[n * dim_count + d];

	memset(&desc, 0, sizeof(desc));
	desc.dimension = curve->dimension;
	desc.degree = n + 1;
	desc.control_points = new_cp;
	desc.control_point_count = new_count;

	status = qaws_curve_create_bezier(&desc, out_elevated);
	free(new_cp);
	return status;
}
