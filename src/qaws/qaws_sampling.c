#include "qaws_sampling.h"
#include "qaws_eval.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_arc_length.h"
#include "internal/qaws_internal_cache.h"
#include <stdlib.h>
#include <math.h>

/* ---------------------------------------------------------------------------
 * Adaptive sampling helpers
 * ------------------------------------------------------------------------- */

#define ADAPTIVE_MAX_DEPTH 10

struct adaptive_buf_2d
{
	qaws_vec2 *points;
	unsigned int count;
	unsigned int capacity;
};

static int adaptive_buf_2d_push(struct adaptive_buf_2d *buf, qaws_vec2 pt)
{
	if (buf->count >= buf->capacity)
	{
		unsigned int new_cap = buf->capacity * 2;
		qaws_vec2 *new_ptr = (qaws_vec2 *)realloc(
			buf->points, new_cap * sizeof(qaws_vec2));
		if (!new_ptr)
			return 0;
		buf->points = new_ptr;
		buf->capacity = new_cap;
	}
	buf->points[buf->count++] = pt;
	return 1;
}

static qaws_status adaptive_subdivide_2d(
	qaws_curve const *curve,
	qaws_scalar t0, qaws_vec2 p0,
	qaws_scalar t1, qaws_vec2 p1,
	qaws_scalar tolerance,
	int depth,
	struct adaptive_buf_2d *buf)
{
	qaws_scalar tmid;
	qaws_eval_result_2d result;
	qaws_status status;
	qaws_scalar mx, my, dx, dy, dist;

	if (depth >= ADAPTIVE_MAX_DEPTH)
	{
		if (!adaptive_buf_2d_push(buf, p1))
			return QAWS_STATUS_ALLOCATION_FAILURE;
		return QAWS_STATUS_OK;
	}

	tmid = (t0 + t1) * (qaws_scalar)0.5;
	status = qaws_curve_evaluate_2d(
		curve, tmid, QAWS_EVAL_FLAG_POSITION, &result);
	if (status != QAWS_STATUS_OK)
		return status;

	/* Distance from midpoint of chord to actual curve point */
	mx = (p0.x + p1.x) * (qaws_scalar)0.5;
	my = (p0.y + p1.y) * (qaws_scalar)0.5;
	dx = result.position.x - mx;
	dy = result.position.y - my;
	dist = (qaws_scalar)sqrt((double)(dx * dx + dy * dy));

	if (dist > tolerance)
	{
		status = adaptive_subdivide_2d(
			curve, t0, p0, tmid, result.position,
			tolerance, depth + 1, buf);
		if (status != QAWS_STATUS_OK)
			return status;

		status = adaptive_subdivide_2d(
			curve, tmid, result.position, t1, p1,
			tolerance, depth + 1, buf);
		if (status != QAWS_STATUS_OK)
			return status;
	}
	else
	{
		if (!adaptive_buf_2d_push(buf, p1))
			return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	return QAWS_STATUS_OK;
}

struct adaptive_buf_3d
{
	qaws_vec3 *points;
	unsigned int count;
	unsigned int capacity;
};

static int adaptive_buf_3d_push(struct adaptive_buf_3d *buf, qaws_vec3 pt)
{
	if (buf->count >= buf->capacity)
	{
		unsigned int new_cap = buf->capacity * 2;
		qaws_vec3 *new_ptr = (qaws_vec3 *)realloc(
			buf->points, new_cap * sizeof(qaws_vec3));
		if (!new_ptr)
			return 0;
		buf->points = new_ptr;
		buf->capacity = new_cap;
	}
	buf->points[buf->count++] = pt;
	return 1;
}

static qaws_status adaptive_subdivide_3d(
	qaws_curve const *curve,
	qaws_scalar t0, qaws_vec3 p0,
	qaws_scalar t1, qaws_vec3 p1,
	qaws_scalar tolerance,
	int depth,
	struct adaptive_buf_3d *buf)
{
	qaws_scalar tmid;
	qaws_eval_result_3d result;
	qaws_status status;
	qaws_scalar mx, my, mz, dx, dy, dz, dist;

	if (depth >= ADAPTIVE_MAX_DEPTH)
	{
		if (!adaptive_buf_3d_push(buf, p1))
			return QAWS_STATUS_ALLOCATION_FAILURE;
		return QAWS_STATUS_OK;
	}

	tmid = (t0 + t1) * (qaws_scalar)0.5;
	status = qaws_curve_evaluate_3d(
		curve, tmid, QAWS_EVAL_FLAG_POSITION, &result);
	if (status != QAWS_STATUS_OK)
		return status;

	mx = (p0.x + p1.x) * (qaws_scalar)0.5;
	my = (p0.y + p1.y) * (qaws_scalar)0.5;
	mz = (p0.z + p1.z) * (qaws_scalar)0.5;
	dx = result.position.x - mx;
	dy = result.position.y - my;
	dz = result.position.z - mz;
	dist = (qaws_scalar)sqrt((double)(dx * dx + dy * dy + dz * dz));

	if (dist > tolerance)
	{
		status = adaptive_subdivide_3d(
			curve, t0, p0, tmid, result.position,
			tolerance, depth + 1, buf);
		if (status != QAWS_STATUS_OK)
			return status;

		status = adaptive_subdivide_3d(
			curve, tmid, result.position, t1, p1,
			tolerance, depth + 1, buf);
		if (status != QAWS_STATUS_OK)
			return status;
	}
	else
	{
		if (!adaptive_buf_3d_push(buf, p1))
			return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	return QAWS_STATUS_OK;
}

static qaws_status adaptive_sample_2d(
	qaws_curve const *curve,
	qaws_sampling_desc const *desc,
	qaws_vec2 *out_positions,
	unsigned int position_capacity,
	unsigned int *out_position_count)
{
	unsigned int initial_count;
	unsigned int i;
	qaws_scalar range_min, range_max;
	qaws_scalar tolerance;
	qaws_eval_result_2d result;
	qaws_status status;
	struct adaptive_buf_2d buf;
	qaws_scalar t0, t1;
	qaws_vec2 p0, p1;
	unsigned int copy_count;

	initial_count = desc->sample_count;
	if (initial_count < 2)
		initial_count = 2;

	tolerance = desc->error_tolerance;
	if (tolerance <= (qaws_scalar)0.0)
		tolerance = (qaws_scalar)0.01;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	buf.capacity = initial_count * 4;
	buf.count = 0;
	buf.points = (qaws_vec2 *)malloc(buf.capacity * sizeof(qaws_vec2));
	if (!buf.points)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Evaluate first point */
	status = qaws_curve_evaluate_2d(
		curve, range_min, QAWS_EVAL_FLAG_POSITION, &result);
	if (status != QAWS_STATUS_OK)
	{
		free(buf.points);
		return status;
	}
	p0 = result.position;
	if (!adaptive_buf_2d_push(&buf, p0))
	{
		free(buf.points);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	for (i = 1; i < initial_count; ++i)
	{
		t0 = range_min + (qaws_scalar)(i - 1) * (range_max - range_min)
			/ (qaws_scalar)(initial_count - 1);
		t1 = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)(initial_count - 1);

		status = qaws_curve_evaluate_2d(
			curve, t1, QAWS_EVAL_FLAG_POSITION, &result);
		if (status != QAWS_STATUS_OK)
		{
			free(buf.points);
			return status;
		}
		p1 = result.position;

		status = adaptive_subdivide_2d(
			curve, t0, p0, t1, p1, tolerance, 0, &buf);
		if (status != QAWS_STATUS_OK)
		{
			free(buf.points);
			return status;
		}

		p0 = p1;
	}

	copy_count = buf.count;
	if (copy_count > position_capacity)
	{
		free(buf.points);
		*out_position_count = copy_count;
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	}

	for (i = 0; i < copy_count; ++i)
		out_positions[i] = buf.points[i];

	*out_position_count = copy_count;
	free(buf.points);
	return QAWS_STATUS_OK;
}

static qaws_status adaptive_sample_3d(
	qaws_curve const *curve,
	qaws_sampling_desc const *desc,
	qaws_vec3 *out_positions,
	unsigned int position_capacity,
	unsigned int *out_position_count)
{
	unsigned int initial_count;
	unsigned int i;
	qaws_scalar range_min, range_max;
	qaws_scalar tolerance;
	qaws_eval_result_3d result;
	qaws_status status;
	struct adaptive_buf_3d buf;
	qaws_scalar t0, t1;
	qaws_vec3 p0, p1;
	unsigned int copy_count;

	initial_count = desc->sample_count;
	if (initial_count < 2)
		initial_count = 2;

	tolerance = desc->error_tolerance;
	if (tolerance <= (qaws_scalar)0.0)
		tolerance = (qaws_scalar)0.01;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	buf.capacity = initial_count * 4;
	buf.count = 0;
	buf.points = (qaws_vec3 *)malloc(buf.capacity * sizeof(qaws_vec3));
	if (!buf.points)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	status = qaws_curve_evaluate_3d(
		curve, range_min, QAWS_EVAL_FLAG_POSITION, &result);
	if (status != QAWS_STATUS_OK)
	{
		free(buf.points);
		return status;
	}
	p0 = result.position;
	if (!adaptive_buf_3d_push(&buf, p0))
	{
		free(buf.points);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	for (i = 1; i < initial_count; ++i)
	{
		t0 = range_min + (qaws_scalar)(i - 1) * (range_max - range_min)
			/ (qaws_scalar)(initial_count - 1);
		t1 = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)(initial_count - 1);

		status = qaws_curve_evaluate_3d(
			curve, t1, QAWS_EVAL_FLAG_POSITION, &result);
		if (status != QAWS_STATUS_OK)
		{
			free(buf.points);
			return status;
		}
		p1 = result.position;

		status = adaptive_subdivide_3d(
			curve, t0, p0, t1, p1, tolerance, 0, &buf);
		if (status != QAWS_STATUS_OK)
		{
			free(buf.points);
			return status;
		}

		p0 = p1;
	}

	copy_count = buf.count;
	if (copy_count > position_capacity)
	{
		free(buf.points);
		*out_position_count = copy_count;
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	}

	for (i = 0; i < copy_count; ++i)
		out_positions[i] = buf.points[i];

	*out_position_count = copy_count;
	free(buf.points);
	return QAWS_STATUS_OK;
}

/* ---------------------------------------------------------------------------
 * Adaptive eval sampling helpers
 * ------------------------------------------------------------------------- */

struct adaptive_eval_buf_2d
{
	qaws_eval_result_2d *samples;
	unsigned int count;
	unsigned int capacity;
};

static int adaptive_eval_buf_2d_push(struct adaptive_eval_buf_2d *buf,
	qaws_eval_result_2d sample)
{
	if (buf->count >= buf->capacity)
	{
		unsigned int new_cap = buf->capacity * 2;
		qaws_eval_result_2d *new_ptr = (qaws_eval_result_2d *)realloc(
			buf->samples, new_cap * sizeof(qaws_eval_result_2d));
		if (!new_ptr)
			return 0;
		buf->samples = new_ptr;
		buf->capacity = new_cap;
	}
	buf->samples[buf->count++] = sample;
	return 1;
}

static qaws_status adaptive_eval_subdivide_2d(
	qaws_curve const *curve,
	qaws_scalar t0, qaws_vec2 p0,
	qaws_scalar t1, qaws_vec2 p1,
	qaws_scalar tolerance,
	int depth,
	struct adaptive_eval_buf_2d *buf)
{
	qaws_scalar tmid;
	qaws_eval_result_2d result;
	qaws_eval_result_2d end_result;
	qaws_status status;
	qaws_scalar mx, my, dx, dy, dist;
	unsigned int eval_flags;

	eval_flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1
		| QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3;

	if (depth >= ADAPTIVE_MAX_DEPTH)
	{
		status = qaws_curve_evaluate_2d(
			curve, t1, eval_flags, &end_result);
		if (status != QAWS_STATUS_OK)
			return status;
		if (!adaptive_eval_buf_2d_push(buf, end_result))
			return QAWS_STATUS_ALLOCATION_FAILURE;
		return QAWS_STATUS_OK;
	}

	tmid = (t0 + t1) * (qaws_scalar)0.5;
	status = qaws_curve_evaluate_2d(
		curve, tmid, eval_flags, &result);
	if (status != QAWS_STATUS_OK)
		return status;

	/* Distance from midpoint of chord to actual curve point */
	mx = (p0.x + p1.x) * (qaws_scalar)0.5;
	my = (p0.y + p1.y) * (qaws_scalar)0.5;
	dx = result.position.x - mx;
	dy = result.position.y - my;
	dist = (qaws_scalar)sqrt((double)(dx * dx + dy * dy));

	if (dist > tolerance)
	{
		status = adaptive_eval_subdivide_2d(
			curve, t0, p0, tmid, result.position,
			tolerance, depth + 1, buf);
		if (status != QAWS_STATUS_OK)
			return status;

		status = adaptive_eval_subdivide_2d(
			curve, tmid, result.position, t1, p1,
			tolerance, depth + 1, buf);
		if (status != QAWS_STATUS_OK)
			return status;
	}
	else
	{
		status = qaws_curve_evaluate_2d(
			curve, t1, eval_flags, &end_result);
		if (status != QAWS_STATUS_OK)
			return status;
		if (!adaptive_eval_buf_2d_push(buf, end_result))
			return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	return QAWS_STATUS_OK;
}

struct adaptive_eval_buf_3d
{
	qaws_eval_result_3d *samples;
	unsigned int count;
	unsigned int capacity;
};

static int adaptive_eval_buf_3d_push(struct adaptive_eval_buf_3d *buf,
	qaws_eval_result_3d sample)
{
	if (buf->count >= buf->capacity)
	{
		unsigned int new_cap = buf->capacity * 2;
		qaws_eval_result_3d *new_ptr = (qaws_eval_result_3d *)realloc(
			buf->samples, new_cap * sizeof(qaws_eval_result_3d));
		if (!new_ptr)
			return 0;
		buf->samples = new_ptr;
		buf->capacity = new_cap;
	}
	buf->samples[buf->count++] = sample;
	return 1;
}

static qaws_status adaptive_eval_subdivide_3d(
	qaws_curve const *curve,
	qaws_scalar t0, qaws_vec3 p0,
	qaws_scalar t1, qaws_vec3 p1,
	qaws_scalar tolerance,
	int depth,
	struct adaptive_eval_buf_3d *buf)
{
	qaws_scalar tmid;
	qaws_eval_result_3d result;
	qaws_eval_result_3d end_result;
	qaws_status status;
	qaws_scalar mx, my, mz, dx, dy, dz, dist;
	unsigned int eval_flags;

	eval_flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1
		| QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3;

	if (depth >= ADAPTIVE_MAX_DEPTH)
	{
		status = qaws_curve_evaluate_3d(
			curve, t1, eval_flags, &end_result);
		if (status != QAWS_STATUS_OK)
			return status;
		if (!adaptive_eval_buf_3d_push(buf, end_result))
			return QAWS_STATUS_ALLOCATION_FAILURE;
		return QAWS_STATUS_OK;
	}

	tmid = (t0 + t1) * (qaws_scalar)0.5;
	status = qaws_curve_evaluate_3d(
		curve, tmid, eval_flags, &result);
	if (status != QAWS_STATUS_OK)
		return status;

	mx = (p0.x + p1.x) * (qaws_scalar)0.5;
	my = (p0.y + p1.y) * (qaws_scalar)0.5;
	mz = (p0.z + p1.z) * (qaws_scalar)0.5;
	dx = result.position.x - mx;
	dy = result.position.y - my;
	dz = result.position.z - mz;
	dist = (qaws_scalar)sqrt((double)(dx * dx + dy * dy + dz * dz));

	if (dist > tolerance)
	{
		status = adaptive_eval_subdivide_3d(
			curve, t0, p0, tmid, result.position,
			tolerance, depth + 1, buf);
		if (status != QAWS_STATUS_OK)
			return status;

		status = adaptive_eval_subdivide_3d(
			curve, tmid, result.position, t1, p1,
			tolerance, depth + 1, buf);
		if (status != QAWS_STATUS_OK)
			return status;
	}
	else
	{
		status = qaws_curve_evaluate_3d(
			curve, t1, eval_flags, &end_result);
		if (status != QAWS_STATUS_OK)
			return status;
		if (!adaptive_eval_buf_3d_push(buf, end_result))
			return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	return QAWS_STATUS_OK;
}

static qaws_status adaptive_sample_eval_2d(
	qaws_curve const *curve,
	qaws_sampling_desc const *desc,
	qaws_eval_result_2d *out_samples,
	unsigned int sample_capacity,
	unsigned int *out_sample_count)
{
	unsigned int initial_count;
	unsigned int i;
	qaws_scalar range_min, range_max;
	qaws_scalar tolerance;
	qaws_eval_result_2d result;
	qaws_status status;
	struct adaptive_eval_buf_2d buf;
	qaws_scalar t0, t1;
	qaws_vec2 p0, p1;
	unsigned int copy_count;
	unsigned int eval_flags;

	eval_flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1
		| QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3;

	initial_count = desc->sample_count;
	if (initial_count < 2)
		initial_count = 2;

	tolerance = desc->error_tolerance;
	if (tolerance <= (qaws_scalar)0.0)
		tolerance = (qaws_scalar)0.01;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	buf.capacity = initial_count * 4;
	buf.count = 0;
	buf.samples = (qaws_eval_result_2d *)malloc(
		buf.capacity * sizeof(qaws_eval_result_2d));
	if (!buf.samples)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Evaluate first point */
	status = qaws_curve_evaluate_2d(
		curve, range_min, eval_flags, &result);
	if (status != QAWS_STATUS_OK)
	{
		free(buf.samples);
		return status;
	}
	p0 = result.position;
	if (!adaptive_eval_buf_2d_push(&buf, result))
	{
		free(buf.samples);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	for (i = 1; i < initial_count; ++i)
	{
		t0 = range_min + (qaws_scalar)(i - 1) * (range_max - range_min)
			/ (qaws_scalar)(initial_count - 1);
		t1 = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)(initial_count - 1);

		status = qaws_curve_evaluate_2d(
			curve, t1, eval_flags, &result);
		if (status != QAWS_STATUS_OK)
		{
			free(buf.samples);
			return status;
		}
		p1 = result.position;

		status = adaptive_eval_subdivide_2d(
			curve, t0, p0, t1, p1, tolerance, 0, &buf);
		if (status != QAWS_STATUS_OK)
		{
			free(buf.samples);
			return status;
		}

		p0 = p1;
	}

	copy_count = buf.count;
	if (copy_count > sample_capacity)
	{
		free(buf.samples);
		*out_sample_count = copy_count;
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	}

	for (i = 0; i < copy_count; ++i)
		out_samples[i] = buf.samples[i];

	*out_sample_count = copy_count;
	free(buf.samples);
	return QAWS_STATUS_OK;
}

static qaws_status adaptive_sample_eval_3d(
	qaws_curve const *curve,
	qaws_sampling_desc const *desc,
	qaws_eval_result_3d *out_samples,
	unsigned int sample_capacity,
	unsigned int *out_sample_count)
{
	unsigned int initial_count;
	unsigned int i;
	qaws_scalar range_min, range_max;
	qaws_scalar tolerance;
	qaws_eval_result_3d result;
	qaws_status status;
	struct adaptive_eval_buf_3d buf;
	qaws_scalar t0, t1;
	qaws_vec3 p0, p1;
	unsigned int copy_count;
	unsigned int eval_flags;

	eval_flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1
		| QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3;

	initial_count = desc->sample_count;
	if (initial_count < 2)
		initial_count = 2;

	tolerance = desc->error_tolerance;
	if (tolerance <= (qaws_scalar)0.0)
		tolerance = (qaws_scalar)0.01;

	range_min = curve->parameter_range.min_value;
	range_max = curve->parameter_range.max_value;

	buf.capacity = initial_count * 4;
	buf.count = 0;
	buf.samples = (qaws_eval_result_3d *)malloc(
		buf.capacity * sizeof(qaws_eval_result_3d));
	if (!buf.samples)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	/* Evaluate first point */
	status = qaws_curve_evaluate_3d(
		curve, range_min, eval_flags, &result);
	if (status != QAWS_STATUS_OK)
	{
		free(buf.samples);
		return status;
	}
	p0 = result.position;
	if (!adaptive_eval_buf_3d_push(&buf, result))
	{
		free(buf.samples);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	for (i = 1; i < initial_count; ++i)
	{
		t0 = range_min + (qaws_scalar)(i - 1) * (range_max - range_min)
			/ (qaws_scalar)(initial_count - 1);
		t1 = range_min + (qaws_scalar)i * (range_max - range_min)
			/ (qaws_scalar)(initial_count - 1);

		status = qaws_curve_evaluate_3d(
			curve, t1, eval_flags, &result);
		if (status != QAWS_STATUS_OK)
		{
			free(buf.samples);
			return status;
		}
		p1 = result.position;

		status = adaptive_eval_subdivide_3d(
			curve, t0, p0, t1, p1, tolerance, 0, &buf);
		if (status != QAWS_STATUS_OK)
		{
			free(buf.samples);
			return status;
		}

		p0 = p1;
	}

	copy_count = buf.count;
	if (copy_count > sample_capacity)
	{
		free(buf.samples);
		*out_sample_count = copy_count;
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	}

	for (i = 0; i < copy_count; ++i)
		out_samples[i] = buf.samples[i];

	*out_sample_count = copy_count;
	free(buf.samples);
	return QAWS_STATUS_OK;
}

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

	/* Adaptive sampling path */
	if (desc->use_adaptive_sampling)
		return adaptive_sample_2d(
			curve, desc, out_positions, position_capacity, out_position_count);

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

	/* Adaptive sampling path */
	if (desc->use_adaptive_sampling)
		return adaptive_sample_3d(
			curve, desc, out_positions, position_capacity, out_position_count);

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

	/* Adaptive sampling path */
	if (desc->use_adaptive_sampling)
		return adaptive_sample_eval_2d(
			curve, desc, out_samples, sample_capacity, out_sample_count);

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

	/* Adaptive sampling path */
	if (desc->use_adaptive_sampling)
		return adaptive_sample_eval_3d(
			curve, desc, out_samples, sample_capacity, out_sample_count);

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
