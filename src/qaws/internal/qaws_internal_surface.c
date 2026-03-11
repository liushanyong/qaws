#include "qaws_internal_surface.h"
#include <stdlib.h>
#include <string.h>

qaws_surface* qaws_internal_surface_alloc(
	qaws_surface_kind kind,
	unsigned int u_degree,
	unsigned int v_degree,
	qaws_range u_range,
	qaws_range v_range,
	qaws_surface_vtable const* vtable)
{
	qaws_surface* s = (qaws_surface*)malloc(sizeof(qaws_surface));
	if (!s) return NULL;
	s->kind = kind;
	s->u_degree = u_degree;
	s->v_degree = v_degree;
	s->u_range = u_range;
	s->v_range = v_range;
	s->vtable = vtable;
	s->impl = NULL;
	return s;
}

void qaws_internal_surface_free(qaws_surface* surface)
{
	if (surface)
	{
		if (surface->vtable && surface->vtable->destroy_impl && surface->impl)
			surface->vtable->destroy_impl(surface->impl);
		free(surface);
	}
}

unsigned int qaws_internal_surface_uniform_knots(
	unsigned int degree,
	unsigned int num_cp,
	qaws_scalar* out_knots,
	unsigned int capacity)
{
	unsigned int knot_count = num_cp + degree + 1;
	unsigned int i;
	unsigned int n_internal;
	if (knot_count > capacity) return 0;

	/* Clamped: first (degree+1) knots = 0, last (degree+1) knots = 1 */
	for (i = 0; i <= degree; i++)
		out_knots[i] = (qaws_scalar)0;
	n_internal = knot_count - 2 * (degree + 1);
	for (i = 0; i < n_internal; i++)
		out_knots[degree + 1 + i] = (qaws_scalar)(i + 1) / (qaws_scalar)(n_internal + 1);
	for (i = 0; i <= degree; i++)
		out_knots[knot_count - 1 - i] = (qaws_scalar)1;
	return knot_count;
}
