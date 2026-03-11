#ifndef QAWS_EXPORT_H
#define QAWS_EXPORT_H

#include "qaws_types.h"
#include "qaws_status.h"

typedef struct qaws_bspline_fit_desc
{
	qaws_dimension dimension;
	void const* data_points;           /* m * dim_count scalars (flat array of qaws_vec2 or qaws_vec3) */
	unsigned int data_point_count;     /* m >= degree + 1 */
	unsigned int degree;               /* desired B-spline degree (typically 3) */
	unsigned int control_point_count;  /* n: desired number of CPs (degree < n <= m) */
	qaws_scalar const* parameters;     /* optional: m parameter values in [0,1], NULL = chord-length */
} qaws_bspline_fit_desc;

qaws_status qaws_curve_fit_bspline(
	qaws_bspline_fit_desc const* desc,
	qaws_curve** out_curve);

#endif /* QAWS_EXPORT_H */
