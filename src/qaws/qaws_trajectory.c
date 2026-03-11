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
	qaws_trajectory_impl *impl = (qaws_trajectory_impl *)curve->impl;
	unsigned int const dim_count = 2;
	unsigned int s = span_index;
	qaws_scalar dt;
	qaws_scalar t = local_t;
	qaws_scalar t2 = t * t;
	qaws_scalar t3 = t2 * t;
	qaws_scalar const *cx;
	qaws_scalar const *cy;
	qaws_scalar ax, bx, cx_coeff, dx_val;
	qaws_scalar ay, by, cy_coeff, dy_val;
	qaws_scalar inv_dt, inv_dt2, inv_dt3;

	memset(out_result, 0, sizeof(*out_result));

	dt = impl->key_times[s + 1] - impl->key_times[s];
	if (dt < (qaws_scalar)1e-10)
		return QAWS_STATUS_DEGENERATE_CURVE;

	cx = &impl->span_coeffs[s * dim_count * 4 + 0 * 4];
	cy = &impl->span_coeffs[s * dim_count * 4 + 1 * 4];

	ax = cx[0]; bx = cx[1]; cx_coeff = cx[2]; dx_val = cx[3];
	ay = cy[0]; by = cy[1]; cy_coeff = cy[2]; dy_val = cy[3];

	/* Position: P(u) */
	out_result->valid_flags = QAWS_EVAL_FLAG_POSITION;
	out_result->position.x = ax * t3 + bx * t2 + cx_coeff * t + dx_val;
	out_result->position.y = ay * t3 + by * t2 + cy_coeff * t + dy_val;

	/* D1: P'(u) / dt */
	if ((eval_flags & QAWS_EVAL_FLAG_D1) || (eval_flags & QAWS_EVAL_FLAG_D2) || (eval_flags & QAWS_EVAL_FLAG_D3)) {
		inv_dt = (qaws_scalar)1.0 / dt;
		out_result->d1.x = ((qaws_scalar)3.0 * ax * t2 + (qaws_scalar)2.0 * bx * t + cx_coeff) * inv_dt;
		out_result->d1.y = ((qaws_scalar)3.0 * ay * t2 + (qaws_scalar)2.0 * by * t + cy_coeff) * inv_dt;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	/* D2: P''(u) / dt^2 */
	if ((eval_flags & QAWS_EVAL_FLAG_D2) || (eval_flags & QAWS_EVAL_FLAG_D3)) {
		inv_dt2 = (qaws_scalar)1.0 / (dt * dt);
		out_result->d2.x = ((qaws_scalar)6.0 * ax * t + (qaws_scalar)2.0 * bx) * inv_dt2;
		out_result->d2.y = ((qaws_scalar)6.0 * ay * t + (qaws_scalar)2.0 * by) * inv_dt2;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	/* D3: P'''(u) / dt^3 */
	if (eval_flags & QAWS_EVAL_FLAG_D3) {
		inv_dt3 = (qaws_scalar)1.0 / (dt * dt * dt);
		out_result->d3.x = ((qaws_scalar)6.0 * ax) * inv_dt3;
		out_result->d3.y = ((qaws_scalar)6.0 * ay) * inv_dt3;
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
	qaws_trajectory_impl *impl = (qaws_trajectory_impl *)curve->impl;
	unsigned int const dim_count = 3;
	unsigned int s = span_index;
	qaws_scalar dt;
	qaws_scalar t = local_t;
	qaws_scalar t2 = t * t;
	qaws_scalar t3 = t2 * t;
	qaws_scalar const *cx;
	qaws_scalar const *cy;
	qaws_scalar const *cz;
	qaws_scalar ax, bx, cx_coeff, dx_val;
	qaws_scalar ay, by, cy_coeff, dy_val;
	qaws_scalar az, bz, cz_coeff, dz_val;
	qaws_scalar inv_dt, inv_dt2, inv_dt3;

	memset(out_result, 0, sizeof(*out_result));

	dt = impl->key_times[s + 1] - impl->key_times[s];
	if (dt < (qaws_scalar)1e-10)
		return QAWS_STATUS_DEGENERATE_CURVE;

	cx = &impl->span_coeffs[s * dim_count * 4 + 0 * 4];
	cy = &impl->span_coeffs[s * dim_count * 4 + 1 * 4];
	cz = &impl->span_coeffs[s * dim_count * 4 + 2 * 4];

	ax = cx[0]; bx = cx[1]; cx_coeff = cx[2]; dx_val = cx[3];
	ay = cy[0]; by = cy[1]; cy_coeff = cy[2]; dy_val = cy[3];
	az = cz[0]; bz = cz[1]; cz_coeff = cz[2]; dz_val = cz[3];

	/* Position: P(u) */
	out_result->valid_flags = QAWS_EVAL_FLAG_POSITION;
	out_result->position.x = ax * t3 + bx * t2 + cx_coeff * t + dx_val;
	out_result->position.y = ay * t3 + by * t2 + cy_coeff * t + dy_val;
	out_result->position.z = az * t3 + bz * t2 + cz_coeff * t + dz_val;

	/* D1: P'(u) / dt */
	if ((eval_flags & QAWS_EVAL_FLAG_D1) || (eval_flags & QAWS_EVAL_FLAG_D2) || (eval_flags & QAWS_EVAL_FLAG_D3)) {
		inv_dt = (qaws_scalar)1.0 / dt;
		out_result->d1.x = ((qaws_scalar)3.0 * ax * t2 + (qaws_scalar)2.0 * bx * t + cx_coeff) * inv_dt;
		out_result->d1.y = ((qaws_scalar)3.0 * ay * t2 + (qaws_scalar)2.0 * by * t + cy_coeff) * inv_dt;
		out_result->d1.z = ((qaws_scalar)3.0 * az * t2 + (qaws_scalar)2.0 * bz * t + cz_coeff) * inv_dt;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D1;
	}

	/* D2: P''(u) / dt^2 */
	if ((eval_flags & QAWS_EVAL_FLAG_D2) || (eval_flags & QAWS_EVAL_FLAG_D3)) {
		inv_dt2 = (qaws_scalar)1.0 / (dt * dt);
		out_result->d2.x = ((qaws_scalar)6.0 * ax * t + (qaws_scalar)2.0 * bx) * inv_dt2;
		out_result->d2.y = ((qaws_scalar)6.0 * ay * t + (qaws_scalar)2.0 * by) * inv_dt2;
		out_result->d2.z = ((qaws_scalar)6.0 * az * t + (qaws_scalar)2.0 * bz) * inv_dt2;
		out_result->valid_flags |= QAWS_EVAL_FLAG_D2;
	}

	/* D3: P'''(u) / dt^3 */
	if (eval_flags & QAWS_EVAL_FLAG_D3) {
		inv_dt3 = (qaws_scalar)1.0 / (dt * dt * dt);
		out_result->d3.x = ((qaws_scalar)6.0 * ax) * inv_dt3;
		out_result->d3.y = ((qaws_scalar)6.0 * ay) * inv_dt3;
		out_result->d3.z = ((qaws_scalar)6.0 * az) * inv_dt3;
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
		free(ti->span_coeffs);
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

	/* Precompute per-span polynomial coefficients: a*u^3 + b*u^2 + c*u + d */
	{
		unsigned int s, d;
		impl->span_coeffs = (qaws_scalar *)malloc(
			sizeof(qaws_scalar) * (size_t)span_count * (size_t)dim_count * 4);
		if (!impl->span_coeffs) {
			free(impl->key_accelerations);
			free(impl->key_velocities);
			free(impl->key_times);
			free(impl->key_positions);
			free(impl);
			qaws_curve_destroy(curve);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		for (s = 0; s < span_count; s++) {
			qaws_scalar dt = impl->key_times[s + 1] - impl->key_times[s];
			qaws_scalar const *p0 = impl->key_positions + s * dim_count;
			qaws_scalar const *p1 = impl->key_positions + (s + 1) * dim_count;
			qaws_scalar const *v0 = impl->key_velocities + s * dim_count;
			qaws_scalar const *v1 = impl->key_velocities + (s + 1) * dim_count;
			for (d = 0; d < dim_count; d++) {
				qaws_scalar m0 = v0[d] * dt;
				qaws_scalar m1 = v1[d] * dt;
				qaws_scalar a_coeff = (qaws_scalar)2.0 * p0[d] + m0 - (qaws_scalar)2.0 * p1[d] + m1;
				qaws_scalar b_coeff = (qaws_scalar)-3.0 * p0[d] - (qaws_scalar)2.0 * m0 + (qaws_scalar)3.0 * p1[d] - m1;
				qaws_scalar c_coeff = m0;
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
