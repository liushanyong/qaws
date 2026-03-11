#ifndef QAWS_INSPECT_H
#define QAWS_INSPECT_H

#include "qaws_types.h"
#include "qaws_status.h"

/* All inspection functions are thread-safe on immutable curves. */

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

/* Derived geometric helpers */
qaws_status qaws_curve_compute_tangent_2d(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_vec2* out_tangent);

qaws_status qaws_curve_compute_tangent_3d(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_vec3* out_tangent);

qaws_status qaws_curve_compute_curvature_2d(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_scalar* out_curvature);

qaws_status qaws_curve_compute_curvature_3d(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_scalar* out_curvature);

qaws_status qaws_curve_compute_torsion_3d(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_scalar* out_torsion);

qaws_status qaws_curve_compute_speed(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_scalar* out_speed);

/* Frenet frame */
qaws_status qaws_curve_compute_normal_2d(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_vec2* out_normal);

qaws_status qaws_curve_compute_frenet_frame_3d(
	qaws_curve const* curve,
	qaws_scalar parameter,
	qaws_vec3* out_tangent,
	qaws_vec3* out_normal,
	qaws_vec3* out_binormal);

/* Inflection point detection */
qaws_status qaws_curve_find_inflection_points(
	qaws_curve const* curve,
	qaws_scalar* out_parameters,
	unsigned int parameter_capacity,
	unsigned int* out_count);

/* Extrema detection */
qaws_status qaws_curve_find_extrema(
	qaws_curve const* curve,
	unsigned int axis,
	qaws_scalar* out_parameters,
	unsigned int parameter_capacity,
	unsigned int* out_count);

/* Curvature comb data */
typedef struct qaws_curvature_sample_2d {
	qaws_vec2 position;
	qaws_scalar curvature;
	qaws_vec2 normal;
} qaws_curvature_sample_2d;

typedef struct qaws_curvature_sample_3d {
	qaws_vec3 position;
	qaws_scalar curvature;
	qaws_vec3 normal;
} qaws_curvature_sample_3d;

qaws_status qaws_curve_compute_curvature_comb_2d(
	qaws_curve const* curve,
	unsigned int sample_count,
	qaws_curvature_sample_2d* out_samples,
	unsigned int sample_capacity);

qaws_status qaws_curve_compute_curvature_comb_3d(
	qaws_curve const* curve,
	unsigned int sample_count,
	qaws_curvature_sample_3d* out_samples,
	unsigned int sample_capacity);

/* Winding number */
qaws_status qaws_curve_compute_winding_number_2d(
	qaws_curve const* curve,
	qaws_vec2 point,
	int* out_winding_number);

/* Intersection results */
typedef struct qaws_intersection_2d {
	qaws_scalar parameter_a;
	qaws_scalar parameter_b;
	qaws_vec2 position;
} qaws_intersection_2d;

typedef struct qaws_intersection_3d {
	qaws_scalar parameter_a;
	qaws_scalar parameter_b;
	qaws_vec3 position;
} qaws_intersection_3d;

/* Self-intersection detection */
qaws_status qaws_curve_find_self_intersections_2d(
	qaws_curve const* curve,
	qaws_intersection_2d* out_intersections,
	unsigned int intersection_capacity,
	unsigned int* out_count);

qaws_status qaws_curve_find_self_intersections_3d(
	qaws_curve const* curve,
	qaws_intersection_3d* out_intersections,
	unsigned int intersection_capacity,
	unsigned int* out_count);

/* Curve-curve intersection */
qaws_status qaws_curve_find_intersections_2d(
	qaws_curve const* curve_a,
	qaws_curve const* curve_b,
	qaws_intersection_2d* out_intersections,
	unsigned int intersection_capacity,
	unsigned int* out_count);

qaws_status qaws_curve_find_intersections_3d(
	qaws_curve const* curve_a,
	qaws_curve const* curve_b,
	qaws_intersection_3d* out_intersections,
	unsigned int intersection_capacity,
	unsigned int* out_count);

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
