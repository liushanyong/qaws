#ifndef QAWS_EVAL_H
#define QAWS_EVAL_H

#include "qaws_types.h"
#include "qaws_status.h"

/* All evaluation functions are thread-safe on immutable curves. */

qaws_status qaws_curve_evaluate_2d(
	qaws_curve const* curve,
	qaws_scalar parameter,
	unsigned int eval_flags,
	qaws_eval_result_2d* out_result);

qaws_status qaws_curve_evaluate_3d(
	qaws_curve const* curve,
	qaws_scalar parameter,
	unsigned int eval_flags,
	qaws_eval_result_3d* out_result);

qaws_status qaws_curve_evaluate_span_2d(
	qaws_curve const* curve,
	unsigned int span_index,
	qaws_scalar local_parameter,
	unsigned int eval_flags,
	qaws_eval_result_2d* out_result);

qaws_status qaws_curve_evaluate_span_3d(
	qaws_curve const* curve,
	unsigned int span_index,
	qaws_scalar local_parameter,
	unsigned int eval_flags,
	qaws_eval_result_3d* out_result);

#endif /* QAWS_EVAL_H */
