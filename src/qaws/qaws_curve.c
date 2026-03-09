#include "qaws_curve.h"
#include "internal/qaws_internal_types.h"
#include <stdlib.h>

void qaws_curve_destroy(qaws_curve *curve)
{
	if (curve == NULL) {
		return;
	}

	if (curve->vtable != NULL && curve->vtable->destroy_impl != NULL) {
		curve->vtable->destroy_impl(curve->impl);
	}

	free(curve->span_boundaries);
	free(curve);
}
