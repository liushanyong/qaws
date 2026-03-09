#include "qaws_sampling.h"
#include "qaws_eval.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_arc_length.h"
#include "internal/qaws_internal_cache.h"
#include <stdlib.h>

/* ---------------------------------------------------------------------------
 * get_sample_count
 * ------------------------------------------------------------------------- */

qaws_status qaws_curve_get_sample_count(
	qaws_curve const *curve,
	qaws_sampling_desc const *desc,
	unsigned int *out_sample_count)
{
	if (!curve || !desc || !out_sample_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_sample_count = desc->sample_count;
	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * sample_2d
 * ------------------------------------------------------------------------- */

qaws_status qaws_curve_sample_2d(
	qaws_curve const *curve,
	qaws_sampling_desc const *desc,
	qaws_vec2 *out_positions,
	unsigned int position_capacity,
	unsigned int *out_position_count)
{
	unsigned int count;
	unsigned int i;
	qaws_scalar t;
	qaws_scalar range_min;
	qaws_scalar range_max;
	qaws_eval_result_2d result;
	qaws_status status;

	if (!curve || !desc || !out_positions || !out_position_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	count = desc->sample_count;
	if (count > position_capacity)
		count = position_capacity;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	if (desc->traversal_mode == QAWS_TRAVERSAL_MODE_ARC_LENGTH)
	{
		unsigned int table_size = QAWS_ARC_LENGTH_TABLE_SIZE;
		qaws_scalar *table_params;
		qaws_scalar *table_distances;
		qaws_scalar total_length;
		qaws_scalar target_dist;
		qaws_scalar param;

		table_params = (qaws_scalar *)malloc(table_size * sizeof(qaws_scalar));
		if (!table_params)
			return QAWS_STATUS_ALLOCATION_FAILURE;

		table_distances = (qaws_scalar *)malloc(table_size * sizeof(qaws_scalar));
		if (!table_distances)
		{
			free(table_params);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		status = qaws_internal_build_arc_length_table(
			curve, table_size, table_params, table_distances);
		if (status != QAWS_STATUS_OK)
		{
			free(table_params);
			free(table_distances);
			return status;
		}

		total_length = table_distances[table_size - 1];

		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				target_dist = (qaws_scalar)0.0;
			else
				target_dist = (qaws_scalar)i * total_length / (qaws_scalar)(count - 1);

			param = qaws_internal_distance_to_parameter(
				table_params, table_distances, table_size, target_dist);

			status = qaws_curve_evaluate_2d(
				curve, param, QAWS_EVAL_FLAG_POSITION, &result);
			if (status != QAWS_STATUS_OK)
			{
				free(table_params);
				free(table_distances);
				return status;
			}

			out_positions[i] = result.position;
		}

		free(table_params);
		free(table_distances);
	}
	else
	{
		/* PARAMETER mode (default) */
		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				t = range_min;
			else
				t = range_min + (qaws_scalar)i * (range_max - range_min) / (qaws_scalar)(count - 1);

			status = qaws_curve_evaluate_2d(
				curve, t, QAWS_EVAL_FLAG_POSITION, &result);
			if (status != QAWS_STATUS_OK)
				return status;

			out_positions[i] = result.position;
		}
	}

	*out_position_count = count;
	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * sample_3d
 * ------------------------------------------------------------------------- */

qaws_status qaws_curve_sample_3d(
	qaws_curve const *curve,
	qaws_sampling_desc const *desc,
	qaws_vec3 *out_positions,
	unsigned int position_capacity,
	unsigned int *out_position_count)
{
	unsigned int count;
	unsigned int i;
	qaws_scalar t;
	qaws_scalar range_min;
	qaws_scalar range_max;
	qaws_eval_result_3d result;
	qaws_status status;

	if (!curve || !desc || !out_positions || !out_position_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	count = desc->sample_count;
	if (count > position_capacity)
		count = position_capacity;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	if (desc->traversal_mode == QAWS_TRAVERSAL_MODE_ARC_LENGTH)
	{
		unsigned int table_size = QAWS_ARC_LENGTH_TABLE_SIZE;
		qaws_scalar *table_params;
		qaws_scalar *table_distances;
		qaws_scalar total_length;
		qaws_scalar target_dist;
		qaws_scalar param;

		table_params = (qaws_scalar *)malloc(table_size * sizeof(qaws_scalar));
		if (!table_params)
			return QAWS_STATUS_ALLOCATION_FAILURE;

		table_distances = (qaws_scalar *)malloc(table_size * sizeof(qaws_scalar));
		if (!table_distances)
		{
			free(table_params);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		status = qaws_internal_build_arc_length_table(
			curve, table_size, table_params, table_distances);
		if (status != QAWS_STATUS_OK)
		{
			free(table_params);
			free(table_distances);
			return status;
		}

		total_length = table_distances[table_size - 1];

		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				target_dist = (qaws_scalar)0.0;
			else
				target_dist = (qaws_scalar)i * total_length / (qaws_scalar)(count - 1);

			param = qaws_internal_distance_to_parameter(
				table_params, table_distances, table_size, target_dist);

			status = qaws_curve_evaluate_3d(
				curve, param, QAWS_EVAL_FLAG_POSITION, &result);
			if (status != QAWS_STATUS_OK)
			{
				free(table_params);
				free(table_distances);
				return status;
			}

			out_positions[i] = result.position;
		}

		free(table_params);
		free(table_distances);
	}
	else
	{
		/* PARAMETER mode (default) */
		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				t = range_min;
			else
				t = range_min + (qaws_scalar)i * (range_max - range_min) / (qaws_scalar)(count - 1);

			status = qaws_curve_evaluate_3d(
				curve, t, QAWS_EVAL_FLAG_POSITION, &result);
			if (status != QAWS_STATUS_OK)
				return status;

			out_positions[i] = result.position;
		}
	}

	*out_position_count = count;
	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * sample_eval_2d
 * ------------------------------------------------------------------------- */

qaws_status qaws_curve_sample_eval_2d(
	qaws_curve const *curve,
	qaws_sampling_desc const *desc,
	qaws_eval_result_2d *out_samples,
	unsigned int sample_capacity,
	unsigned int *out_sample_count)
{
	unsigned int count;
	unsigned int i;
	qaws_scalar t;
	qaws_scalar range_min;
	qaws_scalar range_max;
	unsigned int eval_flags;
	qaws_status status;

	if (!curve || !desc || !out_samples || !out_sample_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_2D)
		return QAWS_STATUS_INVALID_DIMENSION;

	count = desc->sample_count;
	if (count > sample_capacity)
		count = sample_capacity;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	eval_flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1
		| QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3;

	if (desc->traversal_mode == QAWS_TRAVERSAL_MODE_ARC_LENGTH)
	{
		unsigned int table_size = QAWS_ARC_LENGTH_TABLE_SIZE;
		qaws_scalar *table_params;
		qaws_scalar *table_distances;
		qaws_scalar total_length;
		qaws_scalar target_dist;
		qaws_scalar param;

		table_params = (qaws_scalar *)malloc(table_size * sizeof(qaws_scalar));
		if (!table_params)
			return QAWS_STATUS_ALLOCATION_FAILURE;

		table_distances = (qaws_scalar *)malloc(table_size * sizeof(qaws_scalar));
		if (!table_distances)
		{
			free(table_params);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		status = qaws_internal_build_arc_length_table(
			curve, table_size, table_params, table_distances);
		if (status != QAWS_STATUS_OK)
		{
			free(table_params);
			free(table_distances);
			return status;
		}

		total_length = table_distances[table_size - 1];

		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				target_dist = (qaws_scalar)0.0;
			else
				target_dist = (qaws_scalar)i * total_length / (qaws_scalar)(count - 1);

			param = qaws_internal_distance_to_parameter(
				table_params, table_distances, table_size, target_dist);

			status = qaws_curve_evaluate_2d(
				curve, param, eval_flags, &out_samples[i]);
			if (status != QAWS_STATUS_OK)
			{
				free(table_params);
				free(table_distances);
				return status;
			}
		}

		free(table_params);
		free(table_distances);
	}
	else
	{
		/* PARAMETER mode (default) */
		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				t = range_min;
			else
				t = range_min + (qaws_scalar)i * (range_max - range_min) / (qaws_scalar)(count - 1);

			status = qaws_curve_evaluate_2d(
				curve, t, eval_flags, &out_samples[i]);
			if (status != QAWS_STATUS_OK)
				return status;
		}
	}

	*out_sample_count = count;
	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * sample_eval_3d
 * ------------------------------------------------------------------------- */

qaws_status qaws_curve_sample_eval_3d(
	qaws_curve const *curve,
	qaws_sampling_desc const *desc,
	qaws_eval_result_3d *out_samples,
	unsigned int sample_capacity,
	unsigned int *out_sample_count)
{
	unsigned int count;
	unsigned int i;
	qaws_scalar t;
	qaws_scalar range_min;
	qaws_scalar range_max;
	unsigned int eval_flags;
	qaws_status status;

	if (!curve || !desc || !out_samples || !out_sample_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (curve->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	count = desc->sample_count;
	if (count > sample_capacity)
		count = sample_capacity;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	eval_flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1
		| QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3;

	if (desc->traversal_mode == QAWS_TRAVERSAL_MODE_ARC_LENGTH)
	{
		unsigned int table_size = QAWS_ARC_LENGTH_TABLE_SIZE;
		qaws_scalar *table_params;
		qaws_scalar *table_distances;
		qaws_scalar total_length;
		qaws_scalar target_dist;
		qaws_scalar param;

		table_params = (qaws_scalar *)malloc(table_size * sizeof(qaws_scalar));
		if (!table_params)
			return QAWS_STATUS_ALLOCATION_FAILURE;

		table_distances = (qaws_scalar *)malloc(table_size * sizeof(qaws_scalar));
		if (!table_distances)
		{
			free(table_params);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		status = qaws_internal_build_arc_length_table(
			curve, table_size, table_params, table_distances);
		if (status != QAWS_STATUS_OK)
		{
			free(table_params);
			free(table_distances);
			return status;
		}

		total_length = table_distances[table_size - 1];

		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				target_dist = (qaws_scalar)0.0;
			else
				target_dist = (qaws_scalar)i * total_length / (qaws_scalar)(count - 1);

			param = qaws_internal_distance_to_parameter(
				table_params, table_distances, table_size, target_dist);

			status = qaws_curve_evaluate_3d(
				curve, param, eval_flags, &out_samples[i]);
			if (status != QAWS_STATUS_OK)
			{
				free(table_params);
				free(table_distances);
				return status;
			}
		}

		free(table_params);
		free(table_distances);
	}
	else
	{
		/* PARAMETER mode (default) */
		for (i = 0; i < count; ++i)
		{
			if (count == 1)
				t = range_min;
			else
				t = range_min + (qaws_scalar)i * (range_max - range_min) / (qaws_scalar)(count - 1);

			status = qaws_curve_evaluate_3d(
				curve, t, eval_flags, &out_samples[i]);
			if (status != QAWS_STATUS_OK)
				return status;
		}
	}

	*out_sample_count = count;
	return QAWS_STATUS_OK;
}
