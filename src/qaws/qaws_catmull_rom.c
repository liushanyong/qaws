#include "qaws_catmull_rom.h"
#include "qaws_curve.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_basis.h"
#include "internal/qaws_internal_validation.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* -------------------------------------------------------------------------- */
/*  Helpers                                                                   */
/* -------------------------------------------------------------------------- */

static qaws_scalar compute_alpha(qaws_parameterization param)
{
	switch (param) {
	case QAWS_PARAMETERIZATION_UNIFORM:
		return (qaws_scalar)0.0;
	case QAWS_PARAMETERIZATION_CHORDAL:
		return (qaws_scalar)1.0;
	case QAWS_PARAMETERIZATION_CENTRIPETAL:
		return (qaws_scalar)0.5;
	case QAWS_PARAMETERIZATION_DEFAULT:
	default:
		return (qaws_scalar)0.5; /* default is centripetal */
	}
}

static void compute_knot_params(
	qaws_scalar const *control_points,
	unsigned int control_point_count,
	unsigned int dim_count,
	qaws_scalar alpha,
	qaws_scalar *knot_params)
{
	unsigned int i, d;

	knot_params[0] = (qaws_scalar)0.0;

	for (i = 1; i < control_point_count; i++) {
		qaws_scalar dist_sq = (qaws_scalar)0.0;
		for (d = 0; d < dim_count; d++) {
			qaws_scalar diff = control_points[i * dim_count + d]
			                 - control_points[(i - 1) * dim_count + d];
			dist_sq += diff * diff;
		}

		if (alpha == (qaws_scalar)0.0) {
			/* Uniform: each interval = 1 */
			knot_params[i] = knot_params[i - 1] + (qaws_scalar)1.0;
		} else {
			qaws_scalar dist = (qaws_scalar)sqrt((double)dist_sq);
			qaws_scalar interval = (qaws_scalar)pow((double)dist, (double)alpha);
			knot_params[i] = knot_params[i - 1] + interval;
		}
	}
}

static void compute_segment_coefficients(
	qaws_scalar const *P0,
	qaws_scalar const *P1,
	qaws_scalar const *P2,
	qaws_scalar const *P3,
	qaws_scalar t0,
	qaws_scalar t1,
	qaws_scalar t2,
	qaws_scalar t3,
	unsigned int dim_count,
	qaws_scalar *coeffs_out)  /* dim_count * 4 scalars */
{
	unsigned int d;
	qaws_scalar dt10 = t1 - t0;
	qaws_scalar dt21 = t2 - t1;
	qaws_scalar dt20 = t2 - t0;
	qaws_scalar dt32 = t3 - t2;
	qaws_scalar dt31 = t3 - t1;

	/* Clamp tiny intervals to avoid division by zero */
	if (dt10 < (qaws_scalar)1e-10 && dt10 > (qaws_scalar)-1e-10)
		dt10 = (qaws_scalar)0.0;
	if (dt21 < (qaws_scalar)1e-10 && dt21 > (qaws_scalar)-1e-10)
		dt21 = (qaws_scalar)0.0;
	if (dt20 < (qaws_scalar)1e-10 && dt20 > (qaws_scalar)-1e-10)
		dt20 = (qaws_scalar)0.0;
	if (dt32 < (qaws_scalar)1e-10 && dt32 > (qaws_scalar)-1e-10)
		dt32 = (qaws_scalar)0.0;
	if (dt31 < (qaws_scalar)1e-10 && dt31 > (qaws_scalar)-1e-10)
		dt31 = (qaws_scalar)0.0;

	for (d = 0; d < dim_count; d++) {
		qaws_scalar p0 = P0[d];
		qaws_scalar p1 = P1[d];
		qaws_scalar p2 = P2[d];
		qaws_scalar p3 = P3[d];
		qaws_scalar m1, m2, M1, M2;
		qaws_scalar a_coeff, b_coeff, c_coeff, d_coeff;

		/* Compute tangent at P1 */
		m1 = (qaws_scalar)0.0;
		if (dt21 != (qaws_scalar)0.0 && dt20 != (qaws_scalar)0.0 && dt10 != (qaws_scalar)0.0) {
			m1 = (p2 - p1) / dt21
			   - (p2 - p0) / dt20
			   + (p1 - p0) / dt10;
		}

		/* Compute tangent at P2 */
		m2 = (qaws_scalar)0.0;
		if (dt32 != (qaws_scalar)0.0 && dt31 != (qaws_scalar)0.0 && dt21 != (qaws_scalar)0.0) {
			m2 = (p3 - p2) / dt32
			   - (p3 - p1) / dt31
			   + (p2 - p1) / dt21;
		}

		/* Scale tangents to unit interval */
		M1 = m1 * dt21;
		M2 = m2 * dt21;

		/* Hermite to power basis: P(u) = a*u^3 + b*u^2 + c*u + d */
		d_coeff = p1;
		c_coeff = M1;
		b_coeff = (qaws_scalar)3.0 * (p2 - p1) - (qaws_scalar)2.0 * M1 - M2;
		a_coeff = (qaws_scalar)2.0 * (p1 - p2) + M1 + M2;

		coeffs_out[d * 4 + 0] = a_coeff;
		coeffs_out[d * 4 + 1] = b_coeff;
		coeffs_out[d * 4 + 2] = c_coeff;
		coeffs_out[d * 4 + 3] = d_coeff;
	}
}

/* -------------------------------------------------------------------------- */
/*  Vtable: eval_span (2D)                                                    */
/* -------------------------------------------------------------------------- */

static qaws_status catmull_rom_eval_span_2d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_2d *out_result)
{
	qaws_catmull_rom_impl const *impl =
		(qaws_catmull_rom_impl const *)curve->impl;
	unsigned int const dim_count = 2;
	qaws_scalar t = local_t;
	qaws_scalar t2 = t * t;
	qaws_scalar t3 = t2 * t;
	qaws_scalar const *cx;
	qaws_scalar const *cy;
	qaws_scalar ax, bx, cx_coeff, dx_val;
	qaws_scalar ay, by, cy_coeff, dy_val;

	(void)curve;

	memset(out_result, 0, sizeof(*out_result));

	cx = &impl->segment_coeffs[span_index * dim_count * 4 + 0 * 4];
	cy = &impl->segment_coeffs[span_index * dim_count * 4 + 1 * 4];

	ax = cx[0]; bx = cx[1]; cx_coeff = cx[2]; dx_val = cx[3];
	ay = cy[0]; by = cy[1]; cy_coeff = cy[2]; dy_val = cy[3];

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position.x = ax * t3 + bx * t2 + cx_coeff * t + dx_val;
		out_result->position.y = ay * t3 + by * t2 + cy_coeff * t + dy_val;
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D1) {
		out_result->d1.x = (qaws_scalar)3.0 * ax * t2 + (qaws_scalar)2.0 * bx * t + cx_coeff;
		out_result->d1.y = (qaws_scalar)3.0 * ay * t2 + (qaws_scalar)2.0 * by * t + cy_coeff;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D2) {
		out_result->d2.x = (qaws_scalar)6.0 * ax * t + (qaws_scalar)2.0 * bx;
		out_result->d2.y = (qaws_scalar)6.0 * ay * t + (qaws_scalar)2.0 * by;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D3) {
		out_result->d3.x = (qaws_scalar)6.0 * ax;
		out_result->d3.y = (qaws_scalar)6.0 * ay;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* -------------------------------------------------------------------------- */
/*  Vtable: eval_span (3D)                                                    */
/* -------------------------------------------------------------------------- */

static qaws_status catmull_rom_eval_span_3d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_3d *out_result)
{
	qaws_catmull_rom_impl const *impl =
		(qaws_catmull_rom_impl const *)curve->impl;
	unsigned int const dim_count = 3;
	qaws_scalar t = local_t;
	qaws_scalar t2 = t * t;
	qaws_scalar t3 = t2 * t;
	qaws_scalar const *cx;
	qaws_scalar const *cy;
	qaws_scalar const *cz;
	qaws_scalar ax, bx, cx_coeff, dx_val;
	qaws_scalar ay, by, cy_coeff, dy_val;
	qaws_scalar az, bz, cz_coeff, dz_val;

	(void)curve;

	memset(out_result, 0, sizeof(*out_result));

	cx = &impl->segment_coeffs[span_index * dim_count * 4 + 0 * 4];
	cy = &impl->segment_coeffs[span_index * dim_count * 4 + 1 * 4];
	cz = &impl->segment_coeffs[span_index * dim_count * 4 + 2 * 4];

	ax = cx[0]; bx = cx[1]; cx_coeff = cx[2]; dx_val = cx[3];
	ay = cy[0]; by = cy[1]; cy_coeff = cy[2]; dy_val = cy[3];
	az = cz[0]; bz = cz[1]; cz_coeff = cz[2]; dz_val = cz[3];

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position.x = ax * t3 + bx * t2 + cx_coeff * t + dx_val;
		out_result->position.y = ay * t3 + by * t2 + cy_coeff * t + dy_val;
		out_result->position.z = az * t3 + bz * t2 + cz_coeff * t + dz_val;
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D1) {
		out_result->d1.x = (qaws_scalar)3.0 * ax * t2 + (qaws_scalar)2.0 * bx * t + cx_coeff;
		out_result->d1.y = (qaws_scalar)3.0 * ay * t2 + (qaws_scalar)2.0 * by * t + cy_coeff;
		out_result->d1.z = (qaws_scalar)3.0 * az * t2 + (qaws_scalar)2.0 * bz * t + cz_coeff;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D2) {
		out_result->d2.x = (qaws_scalar)6.0 * ax * t + (qaws_scalar)2.0 * bx;
		out_result->d2.y = (qaws_scalar)6.0 * ay * t + (qaws_scalar)2.0 * by;
		out_result->d2.z = (qaws_scalar)6.0 * az * t + (qaws_scalar)2.0 * bz;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D3) {
		out_result->d3.x = (qaws_scalar)6.0 * ax;
		out_result->d3.y = (qaws_scalar)6.0 * ay;
		out_result->d3.z = (qaws_scalar)6.0 * az;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* -------------------------------------------------------------------------- */
/*  Vtable: lifecycle and query functions                                     */
/* -------------------------------------------------------------------------- */

static void catmull_rom_destroy_impl(void *impl)
{
	qaws_catmull_rom_impl *cr = (qaws_catmull_rom_impl *)impl;
	if (cr) {
		free(cr->control_points);
		free(cr->knot_params);
		free(cr->segment_coeffs);
		free(cr);
	}
}

static int catmull_rom_is_closed(qaws_curve const *curve)
{
	qaws_catmull_rom_impl const *impl =
		(qaws_catmull_rom_impl const *)curve->impl;
	return impl->closed;
}

static int catmull_rom_is_periodic(qaws_curve const *curve)
{
	qaws_catmull_rom_impl const *impl =
		(qaws_catmull_rom_impl const *)curve->impl;
	return impl->closed;
}

static int catmull_rom_is_rational(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static qaws_continuity catmull_rom_get_continuity(qaws_curve const *curve)
{
	(void)curve;
	return QAWS_CONTINUITY_C1;
}

/* -------------------------------------------------------------------------- */
/*  Vtable                                                                    */
/* -------------------------------------------------------------------------- */

static qaws_curve_vtable const catmull_rom_vtable = {
	catmull_rom_eval_span_2d,
	catmull_rom_eval_span_3d,
	catmull_rom_destroy_impl,
	catmull_rom_is_closed,
	catmull_rom_is_periodic,
	catmull_rom_is_rational,
	catmull_rom_get_continuity,
};

/* -------------------------------------------------------------------------- */
/*  Creation                                                                  */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_create_catmull_rom(
	qaws_catmull_rom_desc const *desc,
	qaws_curve **out_curve)
{
	qaws_catmull_rom_impl *impl = NULL;
	qaws_curve *curve = NULL;
	unsigned int dim_count;
	unsigned int span_count;
	unsigned int N;
	unsigned int i, s;
	qaws_scalar alpha;
	qaws_range range;

	/* --- Validation ---------------------------------------------------- */

	if (!desc || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;

	if (desc->dimension != QAWS_DIMENSION_2D &&
	    desc->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	if (!desc->control_points)
		return QAWS_STATUS_INVALID_ARGUMENT;

	dim_count = (desc->dimension == QAWS_DIMENSION_2D) ? 2u : 3u;
	N = desc->control_point_count;

	if (desc->closed) {
		if (N < 3)
			return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
		span_count = N;
	} else {
		if (N < 4)
			return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
		span_count = N - 3;
	}

	/* --- Allocate implementation --------------------------------------- */

	impl = (qaws_catmull_rom_impl *)calloc(1, sizeof(qaws_catmull_rom_impl));
	if (!impl)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	impl->control_point_count = N;
	impl->parameterization = desc->parameterization;
	impl->closed = desc->closed ? 1 : 0;

	/* Copy control points */
	impl->control_points = (qaws_scalar *)malloc(
		sizeof(qaws_scalar) * (size_t)(N * dim_count));
	if (!impl->control_points) {
		free(impl);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->control_points, desc->control_points,
	       sizeof(qaws_scalar) * (size_t)(N * dim_count));

	/* Compute knot parameters */
	impl->knot_params = (qaws_scalar *)malloc(
		sizeof(qaws_scalar) * (size_t)N);
	if (!impl->knot_params) {
		free(impl->control_points);
		free(impl);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	alpha = compute_alpha(impl->parameterization);
	compute_knot_params(impl->control_points, N, dim_count,
	                    alpha, impl->knot_params);

	/* Allocate segment coefficients */
	impl->segment_coeffs = (qaws_scalar *)calloc(
		(size_t)(span_count * dim_count * 4), sizeof(qaws_scalar));
	if (!impl->segment_coeffs) {
		free(impl->knot_params);
		free(impl->control_points);
		free(impl);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Compute segment coefficients */
	if (desc->closed) {
		for (s = 0; s < span_count; s++) {
			int i0 = (((int)s - 1) % (int)N + (int)N) % (int)N;
			int i1 = (int)s;
			int i2 = (int)((s + 1) % N);
			int i3 = (int)((s + 2) % N);

			qaws_scalar const *P0 = &impl->control_points[i0 * (int)dim_count];
			qaws_scalar const *P1 = &impl->control_points[i1 * (int)dim_count];
			qaws_scalar const *P2 = &impl->control_points[i2 * (int)dim_count];
			qaws_scalar const *P3 = &impl->control_points[i3 * (int)dim_count];

			/*
			 * Build local knot values for these four points.
			 * Start from knot of i1 = 0, then compute relative
			 * distances.
			 */
			qaws_scalar lt0, lt1, lt2, lt3;
			lt1 = (qaws_scalar)0.0;

			/* lt0: distance from i0 to i1 */
			{
				qaws_scalar dist_sq = (qaws_scalar)0.0;
				unsigned int d;
				for (d = 0; d < dim_count; d++) {
					qaws_scalar diff = P1[d] - P0[d];
					dist_sq += diff * diff;
				}
				if (alpha == (qaws_scalar)0.0)
					lt0 = -(qaws_scalar)1.0;
				else
					lt0 = -(qaws_scalar)pow((double)sqrt((double)dist_sq),
					                        (double)alpha);
			}

			/* lt2: distance from i1 to i2 */
			{
				qaws_scalar dist_sq = (qaws_scalar)0.0;
				unsigned int d;
				for (d = 0; d < dim_count; d++) {
					qaws_scalar diff = P2[d] - P1[d];
					dist_sq += diff * diff;
				}
				if (alpha == (qaws_scalar)0.0)
					lt2 = (qaws_scalar)1.0;
				else
					lt2 = (qaws_scalar)pow((double)sqrt((double)dist_sq),
					                       (double)alpha);
			}

			/* lt3: distance from i2 to i3 */
			{
				qaws_scalar dist_sq = (qaws_scalar)0.0;
				unsigned int d;
				for (d = 0; d < dim_count; d++) {
					qaws_scalar diff = P3[d] - P2[d];
					dist_sq += diff * diff;
				}
				if (alpha == (qaws_scalar)0.0)
					lt3 = lt2 + (qaws_scalar)1.0;
				else
					lt3 = lt2 + (qaws_scalar)pow(
						(double)sqrt((double)dist_sq),
						(double)alpha);
			}

			compute_segment_coefficients(
				P0, P1, P2, P3,
				lt0, lt1, lt2, lt3,
				dim_count,
				&impl->segment_coeffs[s * dim_count * 4]);
		}
	} else {
		/* Open curve: segment s uses points [s, s+1, s+2, s+3] */
		for (s = 0; s < span_count; s++) {
			qaws_scalar const *P0 = &impl->control_points[(s + 0) * dim_count];
			qaws_scalar const *P1 = &impl->control_points[(s + 1) * dim_count];
			qaws_scalar const *P2 = &impl->control_points[(s + 2) * dim_count];
			qaws_scalar const *P3 = &impl->control_points[(s + 3) * dim_count];

			qaws_scalar t0 = impl->knot_params[s + 0];
			qaws_scalar t1 = impl->knot_params[s + 1];
			qaws_scalar t2 = impl->knot_params[s + 2];
			qaws_scalar t3 = impl->knot_params[s + 3];

			compute_segment_coefficients(
				P0, P1, P2, P3,
				t0, t1, t2, t3,
				dim_count,
				&impl->segment_coeffs[s * dim_count * 4]);
		}
	}

	/* --- Allocate curve ------------------------------------------------ */

	range.min_value = (qaws_scalar)0.0;
	range.max_value = (qaws_scalar)span_count;

	curve = qaws_internal_curve_alloc(
		QAWS_CURVE_KIND_CATMULL_ROM,
		desc->dimension,
		3,            /* degree */
		span_count,
		range,
		&catmull_rom_vtable);
	if (!curve) {
		free(impl->segment_coeffs);
		free(impl->knot_params);
		free(impl->control_points);
		free(impl);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Set span boundaries: uniform [0..span_count] */
	for (i = 0; i <= span_count; i++) {
		curve->span_boundaries[i] = (qaws_scalar)i;
	}

	curve->impl = impl;
	*out_curve = curve;

	return QAWS_STATUS_OK;
}
