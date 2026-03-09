#ifndef QAWS_SAMPLING_H
#define QAWS_SAMPLING_H

#include "qaws_types.h"
#include "qaws_status.h"

/* All sampling functions are thread-safe on immutable curves. */

qaws_status qaws_curve_get_sample_count(
	qaws_curve const* curve,
	qaws_sampling_desc const* desc,
	unsigned int* out_sample_count);

qaws_status qaws_curve_sample_2d(
	qaws_curve const* curve,
	qaws_sampling_desc const* desc,
	qaws_vec2* out_positions,
	unsigned int position_capacity,
	unsigned int* out_position_count);

qaws_status qaws_curve_sample_3d(
	qaws_curve const* curve,
	qaws_sampling_desc const* desc,
	qaws_vec3* out_positions,
	unsigned int position_capacity,
	unsigned int* out_position_count);

qaws_status qaws_curve_sample_eval_2d(
	qaws_curve const* curve,
	qaws_sampling_desc const* desc,
	qaws_eval_result_2d* out_samples,
	unsigned int sample_capacity,
	unsigned int* out_sample_count);

qaws_status qaws_curve_sample_eval_3d(
	qaws_curve const* curve,
	qaws_sampling_desc const* desc,
	qaws_eval_result_3d* out_samples,
	unsigned int sample_capacity,
	unsigned int* out_sample_count);

#endif /* QAWS_SAMPLING_H */
