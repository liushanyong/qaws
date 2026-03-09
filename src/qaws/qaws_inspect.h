#ifndef QAWS_INSPECT_H
#define QAWS_INSPECT_H

#include "qaws_types.h"
#include "qaws_status.h"

/* Generic inspection */
qaws_curve_kind qaws_curve_get_kind(qaws_curve const* curve);
qaws_dimension  qaws_curve_get_dimension(qaws_curve const* curve);
unsigned int    qaws_curve_get_degree(qaws_curve const* curve);
unsigned int    qaws_curve_get_span_count(qaws_curve const* curve);
qaws_range      qaws_curve_get_parameter_range(qaws_curve const* curve);
int             qaws_curve_is_closed(qaws_curve const* curve);
int             qaws_curve_is_periodic(qaws_curve const* curve);
int             qaws_curve_is_rational(qaws_curve const* curve);
qaws_continuity qaws_curve_get_continuity(qaws_curve const* curve);

/* Analysis helpers */
qaws_status qaws_curve_compute_arc_length(
	qaws_curve const* curve,
	qaws_scalar parameter_min,
	qaws_scalar parameter_max,
	qaws_scalar* out_length);

qaws_status qaws_curve_compute_bounds_2d(
	qaws_curve const* curve,
	qaws_vec2* out_min,
	qaws_vec2* out_max);

qaws_status qaws_curve_compute_bounds_3d(
	qaws_curve const* curve,
	qaws_vec3* out_min,
	qaws_vec3* out_max);

qaws_status qaws_curve_find_closest_parameter_2d(
	qaws_curve const* curve,
	qaws_vec2 point,
	qaws_scalar* out_parameter);

qaws_status qaws_curve_find_closest_parameter_3d(
	qaws_curve const* curve,
	qaws_vec3 point,
	qaws_scalar* out_parameter);

qaws_status qaws_curve_get_span_continuity(
	qaws_curve const* curve,
	unsigned int boundary_index,
	qaws_continuity* out_continuity);

/* Family-specific inspection */
qaws_status qaws_bezier_get_control_points(
	qaws_curve const* curve,
	void* out_control_points,
	unsigned int point_capacity,
	unsigned int* out_point_count);

qaws_status qaws_bspline_get_knots(
	qaws_curve const* curve,
	qaws_scalar* out_knots,
	unsigned int knot_capacity,
	unsigned int* out_knot_count);

qaws_status qaws_nurbs_get_weights(
	qaws_curve const* curve,
	qaws_scalar* out_weights,
	unsigned int weight_capacity,
	unsigned int* out_weight_count);

#endif /* QAWS_INSPECT_H */
