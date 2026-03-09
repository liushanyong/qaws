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
	qaws_scalar const *p0 = &impl->points[span_index * dim_count];
	qaws_scalar const *p1 = &impl->points[(span_index + 1) * dim_count];
	qaws_scalar const *m0 = &impl->tangents[span_index * dim_count];
	qaws_scalar const *m1 = &impl->tangents[(span_index + 1) * dim_count];
	qaws_scalar h[4];
	qaws_scalar dh[4];
	qaws_scalar d2h[4];
	qaws_scalar d3h[4];

	memset(out_result, 0, sizeof(*out_result));

	qaws_internal_hermite_basis(
		local_t,
		h,
		(eval_flags & (QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3)) ? dh : NULL,
		(eval_flags & (QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3)) ? d2h : NULL,
		(eval_flags & QAWS_EVAL_FLAG_D3) ? d3h : NULL);

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position.x = h[0] * p0[0] + h[1] * m0[0] + h[2] * p1[0] + h[3] * m1[0];
		out_result->position.y = h[0] * p0[1] + h[1] * m0[1] + h[2] * p1[1] + h[3] * m1[1];
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D1) {
		out_result->d1.x = dh[0] * p0[0] + dh[1] * m0[0] + dh[2] * p1[0] + dh[3] * m1[0];
		out_result->d1.y = dh[0] * p0[1] + dh[1] * m0[1] + dh[2] * p1[1] + dh[3] * m1[1];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D2) {
		out_result->d2.x = d2h[0] * p0[0] + d2h[1] * m0[0] + d2h[2] * p1[0] + d2h[3] * m1[0];
		out_result->d2.y = d2h[0] * p0[1] + d2h[1] * m0[1] + d2h[2] * p1[1] + d2h[3] * m1[1];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D3) {
		out_result->d3.x = d3h[0] * p0[0] + d3h[1] * m0[0] + d3h[2] * p1[0] + d3h[3] * m1[0];
		out_result->d3.y = d3h[0] * p0[1] + d3h[1] * m0[1] + d3h[2] * p1[1] + d3h[3] * m1[1];
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
	qaws_scalar const *p0 = &impl->points[span_index * dim_count];
	qaws_scalar const *p1 = &impl->points[(span_index + 1) * dim_count];
	qaws_scalar const *m0 = &impl->tangents[span_index * dim_count];
	qaws_scalar const *m1 = &impl->tangents[(span_index + 1) * dim_count];
	qaws_scalar h[4];
	qaws_scalar dh[4];
	qaws_scalar d2h[4];
	qaws_scalar d3h[4];

	memset(out_result, 0, sizeof(*out_result));

	qaws_internal_hermite_basis(
		local_t,
		h,
		(eval_flags & (QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3)) ? dh : NULL,
		(eval_flags & (QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3)) ? d2h : NULL,
		(eval_flags & QAWS_EVAL_FLAG_D3) ? d3h : NULL);

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		out_result->position.x = h[0] * p0[0] + h[1] * m0[0] + h[2] * p1[0] + h[3] * m1[0];
		out_result->position.y = h[0] * p0[1] + h[1] * m0[1] + h[2] * p1[1] + h[3] * m1[1];
		out_result->position.z = h[0] * p0[2] + h[1] * m0[2] + h[2] * p1[2] + h[3] * m1[2];
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D1) {
		out_result->d1.x = dh[0] * p0[0] + dh[1] * m0[0] + dh[2] * p1[0] + dh[3] * m1[0];
		out_result->d1.y = dh[0] * p0[1] + dh[1] * m0[1] + dh[2] * p1[1] + dh[3] * m1[1];
		out_result->d1.z = dh[0] * p0[2] + dh[1] * m0[2] + dh[2] * p1[2] + dh[3] * m1[2];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D2) {
		out_result->d2.x = d2h[0] * p0[0] + d2h[1] * m0[0] + d2h[2] * p1[0] + d2h[3] * m1[0];
		out_result->d2.y = d2h[0] * p0[1] + d2h[1] * m0[1] + d2h[2] * p1[1] + d2h[3] * m1[1];
		out_result->d2.z = d2h[0] * p0[2] + d2h[1] * m0[2] + d2h[2] * p1[2] + d2h[3] * m1[2];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	if (eval_flags & QAWS_EVAL_FLAG_D3) {
		out_result->d3.x = d3h[0] * p0[0] + d3h[1] * m0[0] + d3h[2] * p1[0] + d3h[3] * m1[0];
		out_result->d3.y = d3h[0] * p0[1] + d3h[1] * m0[1] + d3h[2] * p1[1] + d3h[3] * m1[1];
		out_result->d3.z = d3h[0] * p0[2] + d3h[1] * m0[2] + d3h[2] * p1[2] + d3h[3] * m1[2];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

static void hermite_destroy_impl(void *impl)
{
	qaws_hermite_impl *hi = (qaws_hermite_impl *)impl;
	if (hi) {
		free(hi->points);
		free(hi->tangents);
		free(hi);
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

qaws_status qaws_curve_create_hermite(
	qaws_hermite_desc const *desc,
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

	curve = qaws_internal_curve_alloc(
		QAWS_CURVE_KIND_HERMITE,
		desc->dimension,
		desc->degree,
		span_count,
		range,
		&hermite_vtable);
	if (!curve)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	for (i = 0; i <= span_count; i++) {
		curve->span_boundaries[i] = (qaws_scalar)i;
	}

	impl = (qaws_hermite_impl *)malloc(sizeof(qaws_hermite_impl));
	if (!impl) {
		qaws_internal_curve_free(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	pts_size = (size_t)desc->point_count * (size_t)dim_count * sizeof(qaws_scalar);

	impl->points = (qaws_scalar *)malloc(pts_size);
	if (!impl->points) {
		free(impl);
		qaws_internal_curve_free(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->points, desc->points, pts_size);

	impl->tangents = (qaws_scalar *)malloc(pts_size);
	if (!impl->tangents) {
		free(impl->points);
		free(impl);
		qaws_internal_curve_free(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->tangents, desc->derivatives, pts_size);

	impl->point_count = desc->point_count;

	curve->impl = impl;
	*out_curve = curve;

	return QAWS_STATUS_OK;
}
