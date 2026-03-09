#include "qaws_traversal.h"
#include "qaws_eval.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_arc_length.h"
#include "internal/qaws_internal_cache.h"
#include <stdlib.h>
#include <math.h>

/* ---------------------------------------------------------------------------
 * Easing
 * ------------------------------------------------------------------------- */

#define QAWS_PI (qaws_scalar)3.14159265358979323846

static qaws_scalar apply_easing(qaws_easing easing, qaws_scalar t)
{
	qaws_scalar inv;

	switch (easing)
	{
	case QAWS_EASING_LINEAR:
		return t;

	case QAWS_EASING_QUAD_IN:
		return t * t;

	case QAWS_EASING_QUAD_OUT:
		inv = (qaws_scalar)1.0 - t;
		return (qaws_scalar)1.0 - inv * inv;

	case QAWS_EASING_QUAD_IN_OUT:
		if (t < (qaws_scalar)0.5)
			return (qaws_scalar)2.0 * t * t;
		inv = (qaws_scalar)1.0 - t;
		return (qaws_scalar)1.0 - (qaws_scalar)2.0 * inv * inv;

	case QAWS_EASING_CUBIC_IN:
		return t * t * t;

	case QAWS_EASING_CUBIC_OUT:
		inv = (qaws_scalar)1.0 - t;
		return (qaws_scalar)1.0 - inv * inv * inv;

	case QAWS_EASING_CUBIC_IN_OUT:
		if (t < (qaws_scalar)0.5)
			return (qaws_scalar)4.0 * t * t * t;
		inv = (qaws_scalar)1.0 - t;
		return (qaws_scalar)1.0 - (qaws_scalar)4.0 * inv * inv * inv;

	case QAWS_EASING_SINE_IN:
		return (qaws_scalar)1.0 - (qaws_scalar)cos((double)(t * QAWS_PI / (qaws_scalar)2.0));

	case QAWS_EASING_SINE_OUT:
		return (qaws_scalar)sin((double)(t * QAWS_PI / (qaws_scalar)2.0));

	case QAWS_EASING_SINE_IN_OUT:
		return ((qaws_scalar)1.0 - (qaws_scalar)cos((double)(t * QAWS_PI))) / (qaws_scalar)2.0;

	default:
		return t;
	}
}

/* ---------------------------------------------------------------------------
 * Wrap mode
 * ------------------------------------------------------------------------- */

static qaws_scalar apply_wrap(qaws_wrap_mode mode, qaws_scalar distance, qaws_scalar total_length)
{
	qaws_scalar wrapped;
	qaws_scalar cycle;

	if (total_length <= (qaws_scalar)0.0)
		return (qaws_scalar)0.0;

	switch (mode)
	{
	case QAWS_WRAP_CLAMP:
		if (distance < (qaws_scalar)0.0)
			return (qaws_scalar)0.0;
		if (distance > total_length)
			return total_length;
		return distance;

	case QAWS_WRAP_LOOP:
		wrapped = (qaws_scalar)fmod((double)distance, (double)total_length);
		if (wrapped < (qaws_scalar)0.0)
			wrapped += total_length;
		return wrapped;

	case QAWS_WRAP_PING_PONG:
		/* Fold into [0, 2*total_length) then mirror */
		cycle = (qaws_scalar)2.0 * total_length;
		wrapped = (qaws_scalar)fmod((double)distance, (double)cycle);
		if (wrapped < (qaws_scalar)0.0)
			wrapped += cycle;
		if (wrapped > total_length)
			wrapped = cycle - wrapped;
		return wrapped;

	default:
		return distance;
	}
}

/* ---------------------------------------------------------------------------
 * Motion profile: map time to distance
 * ------------------------------------------------------------------------- */

static qaws_scalar traversal_map_time_to_distance(
	qaws_traversal_desc const *desc,
	qaws_scalar time)
{
	qaws_scalar dt;
	qaws_scalar total_time;
	qaws_scalar accel;
	qaws_scalar t_ramp;
	qaws_scalar d_accel;
	qaws_scalar d_cruise;
	qaws_scalar d_total;
	qaws_scalar dt_decel;

	switch (desc->motion_profile)
	{
	case QAWS_MOTION_PROFILE_CONSTANT_SPEED:
		dt = time - desc->start_time;
		return desc->speed * dt;

	case QAWS_MOTION_PROFILE_CONSTANT_ACCELERATION:
		dt = time - desc->start_time;
		return desc->speed * dt + (qaws_scalar)0.5 * desc->acceleration * dt * dt;

	case QAWS_MOTION_PROFILE_TRAPEZOIDAL_SPEED:
		dt = time - desc->start_time;
		total_time = desc->end_time - desc->start_time;

		if (total_time <= (qaws_scalar)0.0)
			return (qaws_scalar)0.0;

		accel = desc->acceleration > (qaws_scalar)0.0 ? desc->acceleration : (qaws_scalar)1.0;
		t_ramp = desc->max_speed / accel;

		if ((qaws_scalar)2.0 * t_ramp > total_time)
			t_ramp = total_time / (qaws_scalar)2.0;

		if (dt <= t_ramp)
		{
			/* Acceleration phase */
			return (qaws_scalar)0.5 * accel * dt * dt;
		}
		else if (dt <= total_time - t_ramp)
		{
			/* Cruise phase */
			d_accel = (qaws_scalar)0.5 * accel * t_ramp * t_ramp;
			return d_accel + desc->max_speed * (dt - t_ramp);
		}
		else
		{
			/* Deceleration phase */
			dt_decel = total_time - dt;
			d_accel = (qaws_scalar)0.5 * accel * t_ramp * t_ramp;
			d_cruise = desc->max_speed * (total_time - (qaws_scalar)2.0 * t_ramp);
			d_total = d_accel + d_cruise + d_accel;
			return d_total - (qaws_scalar)0.5 * accel * dt_decel * dt_decel;
		}

	case QAWS_MOTION_PROFILE_NONE:
	default:
		return time;
	}
}

/* ---------------------------------------------------------------------------
 * Create / Destroy
 * ------------------------------------------------------------------------- */

qaws_status qaws_traversal_create(
	qaws_curve const *curve,
	qaws_traversal_desc const *desc,
	qaws_traversal **out_traversal)
{
	qaws_traversal *traversal;
	unsigned int table_size;
	qaws_status status;

	if (!curve || !desc || !out_traversal)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_traversal = NULL;

	traversal = (qaws_traversal *)malloc(sizeof(qaws_traversal));
	if (!traversal)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	traversal->curve = curve;
	traversal->desc = *desc;

	table_size = QAWS_ARC_LENGTH_TABLE_SIZE;
	traversal->table_size = table_size;

	traversal->table_params = (qaws_scalar *)malloc(
		table_size * sizeof(qaws_scalar));
	if (!traversal->table_params)
	{
		free(traversal);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	traversal->table_distances = (qaws_scalar *)malloc(
		table_size * sizeof(qaws_scalar));
	if (!traversal->table_distances)
	{
		free(traversal->table_params);
		free(traversal);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	status = qaws_internal_build_arc_length_table(
		curve, table_size, traversal->table_params, traversal->table_distances);
	if (status != QAWS_STATUS_OK)
	{
		free(traversal->table_distances);
		free(traversal->table_params);
		free(traversal);
		return status;
	}

	traversal->total_arc_length = traversal->table_distances[table_size - 1];
	traversal->current_distance = (qaws_scalar)0.0;

	*out_traversal = traversal;
	return QAWS_STATUS_OK;
}

void qaws_traversal_destroy(qaws_traversal *traversal)
{
	if (!traversal)
		return;

	free(traversal->table_params);
	free(traversal->table_distances);
	free(traversal);
}

/* ---------------------------------------------------------------------------
 * Evaluate
 * ------------------------------------------------------------------------- */

qaws_status qaws_traversal_evaluate_2d(
	qaws_traversal const *traversal,
	qaws_scalar input_value,
	unsigned int eval_flags,
	qaws_eval_result_2d *out_result)
{
	qaws_scalar parameter;
	qaws_scalar distance;
	qaws_range range;
	qaws_scalar normalized_time;
	qaws_scalar eased_time;
	qaws_scalar time_span;

	if (!traversal || !out_result)
		return QAWS_STATUS_INVALID_ARGUMENT;

	switch (traversal->desc.traversal_mode)
	{
	case QAWS_TRAVERSAL_MODE_PARAMETER:
		parameter = input_value;
		break;

	case QAWS_TRAVERSAL_MODE_ARC_LENGTH:
		distance = input_value;
		distance = apply_wrap(traversal->desc.wrap_mode, distance,
			traversal->total_arc_length);
		parameter = qaws_internal_distance_to_parameter(
			traversal->table_params, traversal->table_distances,
			traversal->table_size, distance);
		break;

	case QAWS_TRAVERSAL_MODE_TIME:
		eased_time = input_value;
		if (traversal->desc.easing != QAWS_EASING_LINEAR)
		{
			time_span = traversal->desc.end_time - traversal->desc.start_time;
			if (time_span > (qaws_scalar)0.0)
			{
				normalized_time = (input_value - traversal->desc.start_time) / time_span;
				if (normalized_time < (qaws_scalar)0.0)
					normalized_time = (qaws_scalar)0.0;
				if (normalized_time > (qaws_scalar)1.0)
					normalized_time = (qaws_scalar)1.0;
				normalized_time = apply_easing(traversal->desc.easing, normalized_time);
				eased_time = traversal->desc.start_time + normalized_time * time_span;
			}
		}
		distance = traversal_map_time_to_distance(
			&traversal->desc, eased_time);
		distance = apply_wrap(traversal->desc.wrap_mode, distance,
			traversal->total_arc_length);
		parameter = qaws_internal_distance_to_parameter(
			traversal->table_params, traversal->table_distances,
			traversal->table_size, distance);
		break;

	default:
		parameter = input_value;
		break;
	}

	if (traversal->desc.clamp_to_domain)
	{
		range = traversal->curve->parameter_range;
		if (parameter < range.min_value)
			parameter = range.min_value;
		if (parameter > range.max_value)
			parameter = range.max_value;
	}

	return qaws_curve_evaluate_2d(
		traversal->curve, parameter, eval_flags, out_result);
}

qaws_status qaws_traversal_evaluate_3d(
	qaws_traversal const *traversal,
	qaws_scalar input_value,
	unsigned int eval_flags,
	qaws_eval_result_3d *out_result)
{
	qaws_scalar parameter;
	qaws_scalar distance;
	qaws_range range;
	qaws_scalar normalized_time;
	qaws_scalar eased_time;
	qaws_scalar time_span;

	if (!traversal || !out_result)
		return QAWS_STATUS_INVALID_ARGUMENT;

	switch (traversal->desc.traversal_mode)
	{
	case QAWS_TRAVERSAL_MODE_PARAMETER:
		parameter = input_value;
		break;

	case QAWS_TRAVERSAL_MODE_ARC_LENGTH:
		distance = input_value;
		distance = apply_wrap(traversal->desc.wrap_mode, distance,
			traversal->total_arc_length);
		parameter = qaws_internal_distance_to_parameter(
			traversal->table_params, traversal->table_distances,
			traversal->table_size, distance);
		break;

	case QAWS_TRAVERSAL_MODE_TIME:
		eased_time = input_value;
		if (traversal->desc.easing != QAWS_EASING_LINEAR)
		{
			time_span = traversal->desc.end_time - traversal->desc.start_time;
			if (time_span > (qaws_scalar)0.0)
			{
				normalized_time = (input_value - traversal->desc.start_time) / time_span;
				if (normalized_time < (qaws_scalar)0.0)
					normalized_time = (qaws_scalar)0.0;
				if (normalized_time > (qaws_scalar)1.0)
					normalized_time = (qaws_scalar)1.0;
				normalized_time = apply_easing(traversal->desc.easing, normalized_time);
				eased_time = traversal->desc.start_time + normalized_time * time_span;
			}
		}
		distance = traversal_map_time_to_distance(
			&traversal->desc, eased_time);
		distance = apply_wrap(traversal->desc.wrap_mode, distance,
			traversal->total_arc_length);
		parameter = qaws_internal_distance_to_parameter(
			traversal->table_params, traversal->table_distances,
			traversal->table_size, distance);
		break;

	default:
		parameter = input_value;
		break;
	}

	if (traversal->desc.clamp_to_domain)
	{
		range = traversal->curve->parameter_range;
		if (parameter < range.min_value)
			parameter = range.min_value;
		if (parameter > range.max_value)
			parameter = range.max_value;
	}

	return qaws_curve_evaluate_3d(
		traversal->curve, parameter, eval_flags, out_result);
}

/* ---------------------------------------------------------------------------
 * Mapping utilities
 * ------------------------------------------------------------------------- */

qaws_status qaws_traversal_map_time_to_parameter(
	qaws_traversal const *traversal,
	qaws_scalar time_value,
	qaws_scalar *out_parameter)
{
	qaws_scalar distance;
	qaws_scalar parameter;
	qaws_range range;

	if (!traversal || !out_parameter)
		return QAWS_STATUS_INVALID_ARGUMENT;

	distance = traversal_map_time_to_distance(&traversal->desc, time_value);
	parameter = qaws_internal_distance_to_parameter(
		traversal->table_params, traversal->table_distances,
		traversal->table_size, distance);

	if (traversal->desc.clamp_to_domain)
	{
		range = traversal->curve->parameter_range;
		if (parameter < range.min_value)
			parameter = range.min_value;
		if (parameter > range.max_value)
			parameter = range.max_value;
	}

	*out_parameter = parameter;
	return QAWS_STATUS_OK;
}

qaws_status qaws_traversal_map_distance_to_parameter(
	qaws_traversal const *traversal,
	qaws_scalar distance_value,
	qaws_scalar *out_parameter)
{
	qaws_scalar parameter;
	qaws_range range;

	if (!traversal || !out_parameter)
		return QAWS_STATUS_INVALID_ARGUMENT;

	parameter = qaws_internal_distance_to_parameter(
		traversal->table_params, traversal->table_distances,
		traversal->table_size, distance_value);

	if (traversal->desc.clamp_to_domain)
	{
		range = traversal->curve->parameter_range;
		if (parameter < range.min_value)
			parameter = range.min_value;
		if (parameter > range.max_value)
			parameter = range.max_value;
	}

	*out_parameter = parameter;
	return QAWS_STATUS_OK;
}

qaws_status qaws_traversal_map_parameter_to_distance(
	qaws_traversal const *traversal,
	qaws_scalar parameter,
	qaws_scalar *out_distance)
{
	if (!traversal || !out_distance)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_distance = qaws_internal_parameter_to_distance(
		traversal->table_params, traversal->table_distances,
		traversal->table_size, parameter);

	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Advance (iterator) functions
 * ------------------------------------------------------------------------- */

qaws_status qaws_traversal_advance_2d(
	qaws_traversal *traversal,
	qaws_scalar delta_time,
	unsigned int eval_flags,
	qaws_eval_result_2d *out_result)
{
	qaws_scalar distance;
	qaws_scalar parameter;
	qaws_range range;

	if (!traversal || !out_result)
		return QAWS_STATUS_INVALID_ARGUMENT;

	traversal->current_distance += delta_time * traversal->desc.speed;
	distance = apply_wrap(traversal->desc.wrap_mode,
		traversal->current_distance, traversal->total_arc_length);

	parameter = qaws_internal_distance_to_parameter(
		traversal->table_params, traversal->table_distances,
		traversal->table_size, distance);

	if (traversal->desc.clamp_to_domain)
	{
		range = traversal->curve->parameter_range;
		if (parameter < range.min_value)
			parameter = range.min_value;
		if (parameter > range.max_value)
			parameter = range.max_value;
	}

	return qaws_curve_evaluate_2d(
		traversal->curve, parameter, eval_flags, out_result);
}

qaws_status qaws_traversal_advance_3d(
	qaws_traversal *traversal,
	qaws_scalar delta_time,
	unsigned int eval_flags,
	qaws_eval_result_3d *out_result)
{
	qaws_scalar distance;
	qaws_scalar parameter;
	qaws_range range;

	if (!traversal || !out_result)
		return QAWS_STATUS_INVALID_ARGUMENT;

	traversal->current_distance += delta_time * traversal->desc.speed;
	distance = apply_wrap(traversal->desc.wrap_mode,
		traversal->current_distance, traversal->total_arc_length);

	parameter = qaws_internal_distance_to_parameter(
		traversal->table_params, traversal->table_distances,
		traversal->table_size, distance);

	if (traversal->desc.clamp_to_domain)
	{
		range = traversal->curve->parameter_range;
		if (parameter < range.min_value)
			parameter = range.min_value;
		if (parameter > range.max_value)
			parameter = range.max_value;
	}

	return qaws_curve_evaluate_3d(
		traversal->curve, parameter, eval_flags, out_result);
}

qaws_status qaws_traversal_reset(qaws_traversal *traversal)
{
	if (!traversal)
		return QAWS_STATUS_INVALID_ARGUMENT;

	traversal->current_distance = (qaws_scalar)0.0;
	return QAWS_STATUS_OK;
}
