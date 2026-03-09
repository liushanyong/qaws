#ifndef QAWS_INTERNAL_CURVE_H
#define QAWS_INTERNAL_CURVE_H

#include "qaws_internal_types.h"

qaws_curve* qaws_internal_curve_alloc(
	qaws_curve_kind kind,
	qaws_dimension dimension,
	unsigned int degree,
	unsigned int span_count,
	qaws_range parameter_range,
	qaws_curve_vtable const* vtable);

void qaws_internal_curve_free(qaws_curve* curve);

#endif /* QAWS_INTERNAL_CURVE_H */
