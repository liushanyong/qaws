#include "qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>

qaws_curve* qaws_internal_curve_alloc(
	qaws_curve_kind kind,
	qaws_dimension dimension,
	unsigned int degree,
	unsigned int span_count,
	qaws_range parameter_range,
	qaws_curve_vtable const* vtable)
{
	qaws_curve* curve = (qaws_curve*)malloc(sizeof(qaws_curve));
	if (!curve) return NULL;

	curve->span_boundaries = (qaws_scalar*)malloc((size_t)(span_count + 1) * sizeof(qaws_scalar));
	if (!curve->span_boundaries) {
		free(curve);
		return NULL;
	}

	memset(curve->span_boundaries, 0, (size_t)(span_count + 1) * sizeof(qaws_scalar));

	curve->kind = kind;
	curve->dimension = dimension;
	curve->degree = degree;
	curve->span_count = span_count;
	curve->parameter_range = parameter_range;
	curve->vtable = vtable;
	curve->impl = NULL;

	return curve;
}

void qaws_internal_curve_free(qaws_curve* curve)
{
	if (!curve) return;
	free(curve->span_boundaries);
	free(curve);
}
