#ifndef QAWS_TRAVERSAL_H
#define QAWS_TRAVERSAL_H

#include "qaws_types.h"
#include "qaws_status.h"

qaws_status qaws_traversal_create(
	qaws_curve const* curve,
	qaws_traversal_desc const* desc,
	qaws_traversal** out_traversal);

void qaws_traversal_destroy(qaws_traversal* traversal);

qaws_status qaws_traversal_evaluate_2d(
	qaws_traversal const* traversal,
	qaws_scalar input_value,
	unsigned int eval_flags,
	qaws_eval_result_2d* out_result);

qaws_status qaws_traversal_evaluate_3d(
	qaws_traversal const* traversal,
	qaws_scalar input_value,
	unsigned int eval_flags,
	qaws_eval_result_3d* out_result);

qaws_status qaws_traversal_map_time_to_parameter(
	qaws_traversal const* traversal,
	qaws_scalar time_value,
	qaws_scalar* out_parameter);

qaws_status qaws_traversal_map_distance_to_parameter(
	qaws_traversal const* traversal,
	qaws_scalar distance_value,
	qaws_scalar* out_parameter);

qaws_status qaws_traversal_map_parameter_to_distance(
	qaws_traversal const* traversal,
	qaws_scalar parameter,
	qaws_scalar* out_distance);

#endif /* QAWS_TRAVERSAL_H */
