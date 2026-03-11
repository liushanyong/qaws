#include "qaws_bezier.h"
#include "qaws_curve.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_basis.h"
#include "internal/qaws_internal_validation.h"
#include <stdlib.h>
#include <string.h>

static qaws_status bezier_eval_span_2d(
	qaws_curve const* curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int eval_flags, qaws_eval_result_2d* out_result)
{
	qaws_bezier_impl const* impl = (qaws_bezier_impl const*)curve->impl;
	unsigned int degree = curve->degree;
	qaws_scalar const* cp = impl->control_points;
	qaws_scalar buf[2], d1p[128], d2p[128], d3p[128];
	(void)span_index;
	memset(out_result, 0, sizeof(*out_result));

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		qaws_internal_decasteljau(cp, degree, 2, local_t, buf);
		out_result->position.x = buf[0]; out_result->position.y = buf[1];
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D1) && degree >= 1) {
		qaws_internal_bezier_derivative_points(cp, degree, 2, d1p);
		qaws_internal_decasteljau(d1p, degree - 1, 2, local_t, buf);
		out_result->d1.x = buf[0]; out_result->d1.y = buf[1];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D2) && degree >= 2) {
		if (!(out_result->valid_flags & QAWS_EVAL_FLAG_D1))
			qaws_internal_bezier_derivative_points(cp, degree, 2, d1p);
		qaws_internal_bezier_derivative_points(d1p, degree - 1, 2, d2p);
		qaws_internal_decasteljau(d2p, degree - 2, 2, local_t, buf);
		out_result->d2.x = buf[0]; out_result->d2.y = buf[1];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D3) && degree >= 3) {
		if (!(out_result->valid_flags & QAWS_EVAL_FLAG_D1))
			qaws_internal_bezier_derivative_points(cp, degree, 2, d1p);
		if (!(out_result->valid_flags & QAWS_EVAL_FLAG_D2))
			qaws_internal_bezier_derivative_points(d1p, degree - 1, 2, d2p);
		qaws_internal_bezier_derivative_points(d2p, degree - 2, 2, d3p);
		qaws_internal_decasteljau(d3p, degree - 3, 2, local_t, buf);
		out_result->d3.x = buf[0]; out_result->d3.y = buf[1];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}
	return QAWS_STATUS_OK;
}

static qaws_status bezier_eval_span_3d(
	qaws_curve const* curve, unsigned int span_index, qaws_scalar local_t,
	unsigned int eval_flags, qaws_eval_result_3d* out_result)
{
	qaws_bezier_impl const* impl = (qaws_bezier_impl const*)curve->impl;
	unsigned int degree = curve->degree;
	qaws_scalar const* cp = impl->control_points;
	qaws_scalar buf[3], d1p[192], d2p[192], d3p[192];
	(void)span_index;
	memset(out_result, 0, sizeof(*out_result));

	if (eval_flags & QAWS_EVAL_FLAG_POSITION) {
		qaws_internal_decasteljau(cp, degree, 3, local_t, buf);
		out_result->position.x = buf[0]; out_result->position.y = buf[1]; out_result->position.z = buf[2];
		out_result->valid_flags |= QAWS_EVAL_FLAG_POSITION;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D1) && degree >= 1) {
		qaws_internal_bezier_derivative_points(cp, degree, 3, d1p);
		qaws_internal_decasteljau(d1p, degree - 1, 3, local_t, buf);
		out_result->d1.x = buf[0]; out_result->d1.y = buf[1]; out_result->d1.z = buf[2];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D2) && degree >= 2) {
		if (!(out_result->valid_flags & QAWS_EVAL_FLAG_D1))
			qaws_internal_bezier_derivative_points(cp, degree, 3, d1p);
		qaws_internal_bezier_derivative_points(d1p, degree - 1, 3, d2p);
		qaws_internal_decasteljau(d2p, degree - 2, 3, local_t, buf);
		out_result->d2.x = buf[0]; out_result->d2.y = buf[1]; out_result->d2.z = buf[2];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}
	if ((eval_flags & QAWS_EVAL_FLAG_D3) && degree >= 3) {
		if (!(out_result->valid_flags & QAWS_EVAL_FLAG_D1))
			qaws_internal_bezier_derivative_points(cp, degree, 3, d1p);
		if (!(out_result->valid_flags & QAWS_EVAL_FLAG_D2))
			qaws_internal_bezier_derivative_points(d1p, degree - 1, 3, d2p);
		qaws_internal_bezier_derivative_points(d2p, degree - 2, 3, d3p);
		qaws_internal_decasteljau(d3p, degree - 3, 3, local_t, buf);
		out_result->d3.x = buf[0]; out_result->d3.y = buf[1]; out_result->d3.z = buf[2];
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}
	return QAWS_STATUS_OK;
}

static void bezier_destroy_impl(void* impl, qaws_allocator const* allocator)
{
	qaws_bezier_impl* bi = (qaws_bezier_impl*)impl;
	if (bi) { qaws_internal_dealloc(allocator, bi->control_points); qaws_internal_dealloc(allocator, bi); }
}

static int bezier_is_closed(qaws_curve const* c)   { (void)c; return 0; }
static int bezier_is_periodic(qaws_curve const* c)  { (void)c; return 0; }
static int bezier_is_rational(qaws_curve const* c)  { (void)c; return 0; }
static qaws_continuity bezier_get_continuity(qaws_curve const* c) { (void)c; return QAWS_CONTINUITY_C3; }

static qaws_curve_vtable const bezier_vtable = {
	bezier_eval_span_2d, bezier_eval_span_3d, bezier_destroy_impl,
	bezier_is_closed, bezier_is_periodic, bezier_is_rational, bezier_get_continuity
};

qaws_status qaws_curve_create_bezier_ex(
	qaws_bezier_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve)
{
	unsigned int dim_count;
	size_t cp_size;
	qaws_range range;
	qaws_curve* curve;
	qaws_bezier_impl* impl;

	if (!desc || !out_curve) return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_internal_validate_dimension(desc->dimension) != QAWS_STATUS_OK)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (desc->degree < 1) return QAWS_STATUS_INVALID_DEGREE;
	if (desc->control_point_count != desc->degree + 1) return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;
	if (!desc->control_points) return QAWS_STATUS_INVALID_ARGUMENT;

	dim_count = (unsigned int)desc->dimension;
	range.min_value = 0; range.max_value = 1;

	curve = qaws_internal_curve_alloc_ex(QAWS_CURVE_KIND_BEZIER, desc->dimension, desc->degree, 1, range, &bezier_vtable, allocator);
	if (!curve) return QAWS_STATUS_ALLOCATION_FAILURE;
	curve->span_boundaries[0] = 0; curve->span_boundaries[1] = 1;

	impl = (qaws_bezier_impl*)qaws_internal_alloc(allocator, sizeof(qaws_bezier_impl));
	if (!impl) { qaws_curve_destroy(curve); return QAWS_STATUS_ALLOCATION_FAILURE; }

	cp_size = (size_t)dim_count * desc->control_point_count * sizeof(qaws_scalar);
	impl->control_points = (qaws_scalar*)qaws_internal_alloc(allocator, (unsigned long)cp_size);
	if (!impl->control_points) { qaws_internal_dealloc(allocator, impl); qaws_curve_destroy(curve); return QAWS_STATUS_ALLOCATION_FAILURE; }
	memcpy(impl->control_points, desc->control_points, cp_size);
	impl->control_point_count = desc->control_point_count;

	curve->impl = impl;
	*out_curve = curve;
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_create_bezier(
	qaws_bezier_desc const* desc, qaws_curve** out_curve)
{
	return qaws_curve_create_bezier_ex(desc, NULL, out_curve);
}
