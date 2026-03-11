#include "qaws_hermite.h"
#include "qaws_curve.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_basis.h"
#include "internal/qaws_internal_validation.h"
#include <stdlib.h>
#include <string.h>

/* -------------------------------------------------------------------------- */
/*  Vtable functions                                                          */
/* -------------------------------------------------------------------------- */

static qaws_status hermite_eval_span_2d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_2d *out_result)
{
	qaws_hermite_impl const *impl = (qaws_hermite_impl const *)curve->impl;
	unsigned int const dim_count = 2;
	qaws_scalar t = local_t;
	qaws_scalar t2 = t * t;
	qaws_scalar t3 = t2 * t;
	qaws_scalar const *cx;
	qaws_scalar const *cy;
	qaws_scalar ax, bx, cx_coeff, dx_val;
	qaws_scalar ay, by, cy_coeff, dy_val;

	memset(out_result, 0, sizeof(*out_result));

	cx = &impl->span_coeffs[span_index * dim_count * 4 + 0 * 4];
	cy = &impl->span_coeffs[span_index * dim_count * 4 + 1 * 4];

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

static qaws_status hermite_eval_span_3d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_3d *out_result)
{
	qaws_hermite_impl const *impl = (qaws_hermite_impl const *)curve->impl;
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

	memset(out_result, 0, sizeof(*out_result));

	cx = &impl->span_coeffs[span_index * dim_count * 4 + 0 * 4];
	cy = &impl->span_coeffs[span_index * dim_count * 4 + 1 * 4];
	cz = &impl->span_coeffs[span_index * dim_count * 4 + 2 * 4];

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

static void hermite_destroy_impl(void *impl, qaws_allocator const* allocator)
{
	qaws_hermite_impl *hi = (qaws_hermite_impl *)impl;
	if (hi) {
		qaws_internal_dealloc(allocator, hi->span_coeffs);
		qaws_internal_dealloc(allocator, hi->points);
		qaws_internal_dealloc(allocator, hi->tangents);
		qaws_internal_dealloc(allocator, hi);
	}
}

static int hermite_is_closed(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static int hermite_is_periodic(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static int hermite_is_rational(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static qaws_continuity hermite_get_continuity(qaws_curve const *curve)
{
	(void)curve;
	return QAWS_CONTINUITY_C1;
}

/* -------------------------------------------------------------------------- */
/*  Vtable                                                                    */
/* -------------------------------------------------------------------------- */

static qaws_curve_vtable const hermite_vtable = {
	hermite_eval_span_2d,
	hermite_eval_span_3d,
	hermite_destroy_impl,
	hermite_is_closed,
	hermite_is_periodic,
	hermite_is_rational,
	hermite_get_continuity,
};

/* -------------------------------------------------------------------------- */
/*  Creation function                                                         */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_create_hermite_ex(
	qaws_hermite_desc const *desc,
	qaws_allocator const *allocator,
	qaws_curve **out_curve)
{
	qaws_status status;
	qaws_curve *curve;
	qaws_hermite_impl *impl;
	unsigned int dim_count;
	unsigned int span_count;
	size_t pts_size;
	qaws_range range;
	unsigned int i;

	if (!desc || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;

	status = qaws_internal_validate_dimension(desc->dimension);
	if (status != QAWS_STATUS_OK)
		return status;

	if (desc->degree != 3)
		return QAWS_STATUS_INVALID_DEGREE;

	if (desc->point_count < 2)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;

	if (desc->derivative_count != desc->point_count)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;

	if (!desc->points || !desc->derivatives)
		return QAWS_STATUS_INVALID_ARGUMENT;

	dim_count = (desc->dimension == QAWS_DIMENSION_2D) ? 2u : 3u;
	span_count = desc->point_count - 1;

	range.min_value = (qaws_scalar)0;
	range.max_value = (qaws_scalar)span_count;

	curve = qaws_internal_curve_alloc_ex(
		QAWS_CURVE_KIND_HERMITE,
		desc->dimension,
		desc->degree,
		span_count,
		range,
		&hermite_vtable,
		allocator);
	if (!curve)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	for (i = 0; i <= span_count; i++) {
		curve->span_boundaries[i] = (qaws_scalar)i;
	}

	impl = (qaws_hermite_impl *)qaws_internal_alloc(allocator, (unsigned long)sizeof(qaws_hermite_impl));
	if (!impl) {
		qaws_internal_curve_free(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	pts_size = (size_t)desc->point_count * (size_t)dim_count * sizeof(qaws_scalar);

	impl->points = (qaws_scalar *)qaws_internal_alloc(allocator, (unsigned long)pts_size);
	if (!impl->points) {
		qaws_internal_dealloc(allocator, impl);
		qaws_internal_curve_free(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->points, desc->points, pts_size);

	impl->tangents = (qaws_scalar *)qaws_internal_alloc(allocator, (unsigned long)pts_size);
	if (!impl->tangents) {
		qaws_internal_dealloc(allocator, impl->points);
		qaws_internal_dealloc(allocator, impl);
		qaws_internal_curve_free(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->tangents, desc->derivatives, pts_size);

	impl->point_count = desc->point_count;

	/* Precompute per-span polynomial coefficients: a*t^3 + b*t^2 + c*t + d */
	{
		unsigned int s, d;
		impl->span_coeffs = (qaws_scalar *)qaws_internal_alloc(allocator,
			(unsigned long)(sizeof(qaws_scalar) * (size_t)span_count * (size_t)dim_count * 4));
		if (!impl->span_coeffs) {
			qaws_internal_dealloc(allocator, impl->tangents);
			qaws_internal_dealloc(allocator, impl->points);
			qaws_internal_dealloc(allocator, impl);
			qaws_internal_curve_free(curve);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		for (s = 0; s < span_count; s++) {
			qaws_scalar const *p0 = impl->points + s * dim_count;
			qaws_scalar const *p1 = impl->points + (s + 1) * dim_count;
			qaws_scalar const *m0 = impl->tangents + s * dim_count;
			qaws_scalar const *m1 = impl->tangents + (s + 1) * dim_count;
			for (d = 0; d < dim_count; d++) {
				qaws_scalar a_coeff = (qaws_scalar)2.0 * p0[d] + m0[d] - (qaws_scalar)2.0 * p1[d] + m1[d];
				qaws_scalar b_coeff = (qaws_scalar)-3.0 * p0[d] - (qaws_scalar)2.0 * m0[d] + (qaws_scalar)3.0 * p1[d] - m1[d];
				qaws_scalar c_coeff = m0[d];
				qaws_scalar d_coeff = p0[d];
				impl->span_coeffs[(s * dim_count + d) * 4 + 0] = a_coeff;
				impl->span_coeffs[(s * dim_count + d) * 4 + 1] = b_coeff;
				impl->span_coeffs[(s * dim_count + d) * 4 + 2] = c_coeff;
				impl->span_coeffs[(s * dim_count + d) * 4 + 3] = d_coeff;
			}
		}
	}

	curve->impl = impl;
	*out_curve = curve;

	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_create_hermite(
	qaws_hermite_desc const *desc,
	qaws_curve **out_curve)
{
	return qaws_curve_create_hermite_ex(desc, NULL, out_curve);
}
