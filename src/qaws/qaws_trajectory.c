#include "qaws_trajectory.h"
#include "qaws_curve.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_basis.h"
#include "internal/qaws_internal_validation.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ---------------------------------------------------------------------------
 * Helpers
 * ------------------------------------------------------------------------- */

static void trajectory_compute_velocities(
	qaws_scalar const *positions,
	qaws_scalar const *times,
	unsigned int key_count,
	unsigned int dim_count,
	qaws_scalar *out_velocities)
{
	unsigned int i, d;
	qaws_scalar dt, inv_dt;
	unsigned int last;

	if (key_count < 2)
		return;

	/* Endpoint 0: forward difference */
	dt = times[1] - times[0];
	if (dt < (qaws_scalar)1e-15) dt = (qaws_scalar)1e-15;
	inv_dt = (qaws_scalar)1.0 / dt;
	for (d = 0; d < dim_count; ++d)
	{
		out_velocities[0 * dim_count + d] =
			(positions[1 * dim_count + d] - positions[0 * dim_count + d]) * inv_dt;
	}

	/* Interior points: Catmull-Rom finite differences */
	for (i = 1; i < key_count - 1; ++i)
	{
		dt = times[i + 1] - times[i - 1];
		if (dt < (qaws_scalar)1e-15) dt = (qaws_scalar)1e-15;
		inv_dt = (qaws_scalar)1.0 / dt;
		for (d = 0; d < dim_count; ++d)
		{
			out_velocities[i * dim_count + d] =
				(positions[(i + 1) * dim_count + d] - positions[(i - 1) * dim_count + d]) * inv_dt;
		}
	}

	/* Endpoint N-1: backward difference */
	last = key_count - 1;
	dt = times[last] - times[last - 1];
	if (dt < (qaws_scalar)1e-15) dt = (qaws_scalar)1e-15;
	inv_dt = (qaws_scalar)1.0 / dt;
	for (d = 0; d < dim_count; ++d)
	{
		out_velocities[last * dim_count + d] =
			(positions[last * dim_count + d] - positions[(last - 1) * dim_count + d]) * inv_dt;
	}
}

/* ---------------------------------------------------------------------------
 * Vtable: eval_span_2d
 * ------------------------------------------------------------------------- */

static qaws_status trajectory_eval_span_2d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_2d *out_result)
{
	unsigned int const dim_count = 2;
	qaws_trajectory_impl *impl = (qaws_trajectory_impl *)curve->impl;
	unsigned int s = span_index;
	qaws_scalar dt;
	qaws_scalar const *p0;
	qaws_scalar const *p1;
	qaws_scalar const *v0;
	qaws_scalar const *v1;
	qaws_scalar m0[2], m1[2];
	qaws_scalar h[4], dh[4], d2h[4], d3h[4];
	int need_d1, need_d2, need_d3;
	unsigned int d;
	qaws_scalar inv_dt, inv_dt2, inv_dt3;

	memset(out_result, 0, sizeof(*out_result));

	dt = impl->key_times[s + 1] - impl->key_times[s];
	if (dt < (qaws_scalar)1e-10)
		return QAWS_STATUS_DEGENERATE_CURVE;

	/* Gather positions and velocities for the span */
	p0 = &impl->key_positions[s * dim_count];
	p1 = &impl->key_positions[(s + 1) * dim_count];
	v0 = &impl->key_velocities[s * dim_count];
	v1 = &impl->key_velocities[(s + 1) * dim_count];

	/* Scaled tangents for Hermite interpolation */
	for (d = 0; d < dim_count; ++d)
	{
		m0[d] = v0[d] * dt;
		m1[d] = v1[d] * dt;
	}

	need_d1 = (eval_flags & QAWS_EVAL_FLAG_D1) != 0;
	need_d2 = (eval_flags & QAWS_EVAL_FLAG_D2) != 0;
	need_d3 = (eval_flags & QAWS_EVAL_FLAG_D3) != 0;

	/* Compute Hermite basis values */
	qaws_internal_hermite_basis(
		local_t,
		h,
		(need_d1 || need_d2 || need_d3) ? dh : NULL,
		(need_d2 || need_d3) ? d2h : NULL,
		need_d3 ? d3h : NULL);

	/* Position */
	out_result->valid_flags = QAWS_EVAL_FLAG_POSITION;
	out_result->position.x = h[0] * p0[0] + h[1] * m0[0] + h[2] * p1[0] + h[3] * m1[0];
	out_result->position.y = h[0] * p0[1] + h[1] * m0[1] + h[2] * p1[1] + h[3] * m1[1];

	/* First derivative: dP/du / dt */
	if (need_d1 || need_d2 || need_d3)
	{
		inv_dt = (qaws_scalar)1.0 / dt;
		out_result->d1.x = (dh[0] * p0[0] + dh[1] * m0[0] + dh[2] * p1[0] + dh[3] * m1[0]) * inv_dt;
		out_result->d1.y = (dh[0] * p0[1] + dh[1] * m0[1] + dh[2] * p1[1] + dh[3] * m1[1]) * inv_dt;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	/* Second derivative: d2P/du2 / dt^2 */
	if (need_d2 || need_d3)
	{
		inv_dt2 = (qaws_scalar)1.0 / (dt * dt);
		out_result->d2.x = (d2h[0] * p0[0] + d2h[1] * m0[0] + d2h[2] * p1[0] + d2h[3] * m1[0]) * inv_dt2;
		out_result->d2.y = (d2h[0] * p0[1] + d2h[1] * m0[1] + d2h[2] * p1[1] + d2h[3] * m1[1]) * inv_dt2;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	/* Third derivative: d3P/du3 / dt^3 */
	if (need_d3)
	{
		inv_dt3 = (qaws_scalar)1.0 / (dt * dt * dt);
		out_result->d3.x = (d3h[0] * p0[0] + d3h[1] * m0[0] + d3h[2] * p1[0] + d3h[3] * m1[0]) * inv_dt3;
		out_result->d3.y = (d3h[0] * p0[1] + d3h[1] * m0[1] + d3h[2] * p1[1] + d3h[3] * m1[1]) * inv_dt3;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Vtable: eval_span_3d
 * ------------------------------------------------------------------------- */

static qaws_status trajectory_eval_span_3d(
	qaws_curve const *curve,
	unsigned int span_index,
	qaws_scalar local_t,
	unsigned int eval_flags,
	qaws_eval_result_3d *out_result)
{
	unsigned int const dim_count = 3;
	qaws_trajectory_impl *impl = (qaws_trajectory_impl *)curve->impl;
	unsigned int s = span_index;
	qaws_scalar dt;
	qaws_scalar const *p0;
	qaws_scalar const *p1;
	qaws_scalar const *v0;
	qaws_scalar const *v1;
	qaws_scalar m0[3], m1[3];
	qaws_scalar h[4], dh[4], d2h[4], d3h[4];
	int need_d1, need_d2, need_d3;
	unsigned int d;
	qaws_scalar inv_dt, inv_dt2, inv_dt3;

	memset(out_result, 0, sizeof(*out_result));

	dt = impl->key_times[s + 1] - impl->key_times[s];
	if (dt < (qaws_scalar)1e-10)
		return QAWS_STATUS_DEGENERATE_CURVE;

	/* Gather positions and velocities for the span */
	p0 = &impl->key_positions[s * dim_count];
	p1 = &impl->key_positions[(s + 1) * dim_count];
	v0 = &impl->key_velocities[s * dim_count];
	v1 = &impl->key_velocities[(s + 1) * dim_count];

	/* Scaled tangents for Hermite interpolation */
	for (d = 0; d < dim_count; ++d)
	{
		m0[d] = v0[d] * dt;
		m1[d] = v1[d] * dt;
	}

	need_d1 = (eval_flags & QAWS_EVAL_FLAG_D1) != 0;
	need_d2 = (eval_flags & QAWS_EVAL_FLAG_D2) != 0;
	need_d3 = (eval_flags & QAWS_EVAL_FLAG_D3) != 0;

	/* Compute Hermite basis values */
	qaws_internal_hermite_basis(
		local_t,
		h,
		(need_d1 || need_d2 || need_d3) ? dh : NULL,
		(need_d2 || need_d3) ? d2h : NULL,
		need_d3 ? d3h : NULL);

	/* Position */
	out_result->valid_flags = QAWS_EVAL_FLAG_POSITION;
	out_result->position.x = h[0] * p0[0] + h[1] * m0[0] + h[2] * p1[0] + h[3] * m1[0];
	out_result->position.y = h[0] * p0[1] + h[1] * m0[1] + h[2] * p1[1] + h[3] * m1[1];
	out_result->position.z = h[0] * p0[2] + h[1] * m0[2] + h[2] * p1[2] + h[3] * m1[2];

	/* First derivative: dP/du / dt */
	if (need_d1 || need_d2 || need_d3)
	{
		inv_dt = (qaws_scalar)1.0 / dt;
		out_result->d1.x = (dh[0] * p0[0] + dh[1] * m0[0] + dh[2] * p1[0] + dh[3] * m1[0]) * inv_dt;
		out_result->d1.y = (dh[0] * p0[1] + dh[1] * m0[1] + dh[2] * p1[1] + dh[3] * m1[1]) * inv_dt;
		out_result->d1.z = (dh[0] * p0[2] + dh[1] * m0[2] + dh[2] * p1[2] + dh[3] * m1[2]) * inv_dt;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	/* Second derivative: d2P/du2 / dt^2 */
	if (need_d2 || need_d3)
	{
		inv_dt2 = (qaws_scalar)1.0 / (dt * dt);
		out_result->d2.x = (d2h[0] * p0[0] + d2h[1] * m0[0] + d2h[2] * p1[0] + d2h[3] * m1[0]) * inv_dt2;
		out_result->d2.y = (d2h[0] * p0[1] + d2h[1] * m0[1] + d2h[2] * p1[1] + d2h[3] * m1[1]) * inv_dt2;
		out_result->d2.z = (d2h[0] * p0[2] + d2h[1] * m0[2] + d2h[2] * p1[2] + d2h[3] * m1[2]) * inv_dt2;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	/* Third derivative: d3P/du3 / dt^3 */
	if (need_d3)
	{
		inv_dt3 = (qaws_scalar)1.0 / (dt * dt * dt);
		out_result->d3.x = (d3h[0] * p0[0] + d3h[1] * m0[0] + d3h[2] * p1[0] + d3h[3] * m1[0]) * inv_dt3;
		out_result->d3.y = (d3h[0] * p0[1] + d3h[1] * m0[1] + d3h[2] * p1[1] + d3h[3] * m1[1]) * inv_dt3;
		out_result->d3.z = (d3h[0] * p0[2] + d3h[1] * m0[2] + d3h[2] * p1[2] + d3h[3] * m1[2]) * inv_dt3;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D3;
	}

	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Vtable: destroy
 * ------------------------------------------------------------------------- */

static void trajectory_destroy_impl(void *impl)
{
	qaws_trajectory_impl *ti = (qaws_trajectory_impl *)impl;
	if (ti)
	{
		free(ti->key_positions);
		free(ti->key_times);
		free(ti->key_velocities);
		free(ti->key_accelerations);
		free(ti);
	}
}

/* ---------------------------------------------------------------------------
 * Vtable: property queries
 * ------------------------------------------------------------------------- */

static int trajectory_is_closed(qaws_curve const *curve)
{
	qaws_trajectory_impl *impl = (qaws_trajectory_impl *)curve->impl;
	return impl->closed;
}

static int trajectory_is_periodic(qaws_curve const *curve)
{
	qaws_trajectory_impl *impl = (qaws_trajectory_impl *)curve->impl;
	return impl->closed;
}

static int trajectory_is_rational(qaws_curve const *curve)
{
	(void)curve;
	return 0;
}

static qaws_continuity trajectory_get_continuity(qaws_curve const *curve)
{
	(void)curve;
	return QAWS_CONTINUITY_C1;
}

/* ---------------------------------------------------------------------------
 * Vtable definition
 * ------------------------------------------------------------------------- */

static qaws_curve_vtable const trajectory_vtable = {
	trajectory_eval_span_2d,
	trajectory_eval_span_3d,
	trajectory_destroy_impl,
	trajectory_is_closed,
	trajectory_is_periodic,
	trajectory_is_rational,
	trajectory_get_continuity
};

/* ---------------------------------------------------------------------------
 * Creation
 * ------------------------------------------------------------------------- */

qaws_status qaws_curve_create_trajectory(
	qaws_trajectory_desc const *desc,
	qaws_curve **out_curve)
{
	unsigned int dim_count;
	unsigned int key_count;
	unsigned int span_count;
	unsigned int i;
	qaws_range parameter_range;
	qaws_curve *curve;
	qaws_trajectory_impl *impl;
	size_t pos_size;
	size_t time_size;
	size_t vel_size;
	size_t acc_size;

	if (!desc)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;

	/* Validate dimension */
	if (qaws_internal_validate_dimension(desc->dimension) != QAWS_STATUS_OK)
		return QAWS_STATUS_INVALID_DIMENSION;

	dim_count = (unsigned int)desc->dimension;

	/* Only degree 3 (cubic Hermite) is supported */
	if (desc->degree != 3)
		return QAWS_STATUS_INVALID_DEGREE;

	key_count = desc->key_count;
	if (key_count < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (!desc->key_positions)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (!desc->key_times)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (desc->key_time_count != key_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	/* Validate key_times are strictly increasing */
	for (i = 1; i < key_count; ++i)
	{
		if (desc->key_times[i] <= desc->key_times[i - 1])
			return QAWS_STATUS_INVALID_PARAMETER_RANGE;
	}

	span_count = key_count - 1;

	/* Parameter range */
	parameter_range.min_value = desc->key_times[0];
	parameter_range.max_value = desc->key_times[key_count - 1];

	/* Allocate curve */
	curve = qaws_internal_curve_alloc(
		QAWS_CURVE_KIND_TRAJECTORY,
		desc->dimension,
		desc->degree,
		span_count,
		parameter_range,
		&trajectory_vtable);

	if (!curve)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Set span boundaries from key_times */
	for (i = 0; i < key_count; ++i)
		curve->span_boundaries[i] = desc->key_times[i];

	/* Allocate impl */
	impl = (qaws_trajectory_impl *)calloc(1, sizeof(qaws_trajectory_impl));
	if (!impl)
	{
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	impl->key_count = key_count;
	impl->closed = desc->closed;

	/* Copy key_positions */
	pos_size = (size_t)key_count * (size_t)dim_count * sizeof(qaws_scalar);
	impl->key_positions = (qaws_scalar *)malloc(pos_size);
	if (!impl->key_positions)
	{
		free(impl);
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->key_positions, desc->key_positions, pos_size);

	/* Copy key_times */
	time_size = (size_t)key_count * sizeof(qaws_scalar);
	impl->key_times = (qaws_scalar *)malloc(time_size);
	if (!impl->key_times)
	{
		free(impl->key_positions);
		free(impl);
		qaws_curve_destroy(curve);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(impl->key_times, desc->key_times, time_size);

	/* Handle velocities */
	vel_size = (size_t)key_count * (size_t)dim_count * sizeof(qaws_scalar);

	if (desc->key_velocities && desc->key_velocity_count == key_count)
	{
		/* Copy user-provided velocities */
		impl->key_velocities = (qaws_scalar *)malloc(vel_size);
		if (!impl->key_velocities)
		{
			free(impl->key_times);
			free(impl->key_positions);
			free(impl);
			qaws_curve_destroy(curve);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		memcpy(impl->key_velocities, desc->key_velocities, vel_size);
		impl->key_velocity_count = key_count;
	}
	else
	{
		/* Compute velocities using finite differences */
		impl->key_velocities = (qaws_scalar *)malloc(vel_size);
		if (!impl->key_velocities)
		{
			free(impl->key_times);
			free(impl->key_positions);
			free(impl);
			qaws_curve_destroy(curve);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		trajectory_compute_velocities(
			impl->key_positions, impl->key_times,
			key_count, dim_count, impl->key_velocities);
		impl->key_velocity_count = key_count;
	}

	/* Copy key_accelerations if provided */
	impl->key_accelerations = NULL;
	impl->key_acceleration_count = 0;
	if (desc->key_accelerations && desc->key_acceleration_count > 0)
	{
		acc_size = (size_t)desc->key_acceleration_count * (size_t)dim_count * sizeof(qaws_scalar);
		impl->key_accelerations = (qaws_scalar *)malloc(acc_size);
		if (!impl->key_accelerations)
		{
			free(impl->key_velocities);
			free(impl->key_times);
			free(impl->key_positions);
			free(impl);
			qaws_curve_destroy(curve);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		memcpy(impl->key_accelerations, desc->key_accelerations, acc_size);
		impl->key_acceleration_count = desc->key_acceleration_count;
	}

	curve->impl = impl;
	*out_curve = curve;

	return QAWS_STATUS_OK;
}
