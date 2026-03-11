#include "qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>

void* qaws_internal_alloc(qaws_allocator const* allocator, unsigned long size)
{
	if (allocator && allocator->alloc)
		return allocator->alloc(size, allocator->user_data);
	return malloc((size_t)size);
}

void qaws_internal_dealloc(qaws_allocator const* allocator, void* ptr)
{
	if (allocator && allocator->dealloc)
	{
		allocator->dealloc(ptr, allocator->user_data);
		return;
	}
	free(ptr);
}

qaws_curve* qaws_internal_curve_alloc_ex(
	qaws_curve_kind kind,
	qaws_dimension dimension,
	unsigned int degree,
	unsigned int span_count,
	qaws_range parameter_range,
	qaws_curve_vtable const* vtable,
	qaws_allocator const* allocator)
{
	qaws_curve* curve = (qaws_curve*)qaws_internal_alloc(allocator, sizeof(qaws_curve));
	if (!curve) return NULL;

	curve->span_boundaries = (qaws_scalar*)qaws_internal_alloc(
		allocator, (unsigned long)(span_count + 1) * sizeof(qaws_scalar));
	if (!curve->span_boundaries) {
		qaws_internal_dealloc(allocator, curve);
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
	curve->allocator = allocator;

	return curve;
}

qaws_curve* qaws_internal_curve_alloc(
	qaws_curve_kind kind,
	qaws_dimension dimension,
	unsigned int degree,
	unsigned int span_count,
	qaws_range parameter_range,
	qaws_curve_vtable const* vtable)
{
	return qaws_internal_curve_alloc_ex(
		kind, dimension, degree, span_count,
		parameter_range, vtable, NULL);
}

void qaws_internal_curve_free(qaws_curve* curve)
{
	if (!curve) return;
	qaws_internal_dealloc(curve->allocator, curve->span_boundaries);
	qaws_internal_dealloc(curve->allocator, curve);
}
