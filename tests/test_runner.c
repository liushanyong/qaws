#include "qaws.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _WIN32
#include <direct.h>
#define MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#define MKDIR(path) mkdir(path, 0755)
#endif

#define SVG_OUTPUT_DIR "working_dir/tests_output"
#define SVG_SAMPLES 256

static int g_pass = 0;
static int g_fail = 0;

#define TEST_ASSERT(cond, msg) \
	do { \
		if (cond) { g_pass++; } \
		else { g_fail++; printf("  FAIL: %s (line %d)\n", msg, __LINE__); } \
	} while (0)

#define TEST_ASSERT_STATUS(status) \
	TEST_ASSERT((status) == QAWS_STATUS_OK, #status " == OK")

#if QAWS_SCALAR_IS_FLOAT
#define TOLERANCE 1e-4f
#else
#define TOLERANCE 1e-8
#endif

static int approx_eq(qaws_scalar a, qaws_scalar b)
{
	qaws_scalar diff = a - b;
	if (diff < 0) diff = -diff;
	return diff < TOLERANCE;
}

static void test_status(void)
{
	printf("test_status\n");
	TEST_ASSERT(qaws_status_to_string(QAWS_STATUS_OK) != NULL, "status_to_string OK");
	TEST_ASSERT(qaws_status_to_string(QAWS_STATUS_INTERNAL_ERROR) != NULL, "status_to_string INTERNAL_ERROR");
}

static void test_bezier_linear(void)
{
	printf("test_bezier_linear\n");

	qaws_vec2 points[] = { {0, 0}, {1, 1} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(curve != NULL, "curve allocated");

	TEST_ASSERT(qaws_curve_get_kind(curve) == QAWS_CURVE_KIND_BEZIER, "kind == BEZIER");
	TEST_ASSERT(qaws_curve_get_degree(curve) == 1, "degree == 1");
	TEST_ASSERT(qaws_curve_get_span_count(curve) == 1, "span_count == 1");

	qaws_eval_result_2d result;
	s = qaws_curve_evaluate_2d(curve, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &result);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(result.position.x, (qaws_scalar)0.5), "midpoint x");
	TEST_ASSERT(approx_eq(result.position.y, (qaws_scalar)0.5), "midpoint y");

	qaws_curve_destroy(curve);
}

static void test_bezier_cubic(void)
{
	printf("test_bezier_cubic\n");

	qaws_vec2 points[] = { {0, 0}, {0, 1}, {1, 1}, {1, 0} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = points;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	/* Evaluate at endpoints */
	qaws_eval_result_2d r0, r1;
	s = qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r0);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r0.position.x, (qaws_scalar)0.0), "start x");
	TEST_ASSERT(approx_eq(r0.position.y, (qaws_scalar)0.0), "start y");

	s = qaws_curve_evaluate_2d(curve, (qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION, &r1);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r1.position.x, (qaws_scalar)1.0), "end x");
	TEST_ASSERT(approx_eq(r1.position.y, (qaws_scalar)0.0), "end y");

	/* Evaluate at midpoint */
	qaws_eval_result_2d rm;
	s = qaws_curve_evaluate_2d(curve, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3, &rm);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(rm.position.x, (qaws_scalar)0.5), "mid x");
	TEST_ASSERT(approx_eq(rm.position.y, (qaws_scalar)0.75), "mid y");

	qaws_curve_destroy(curve);
}

static void test_bezier_3d(void)
{
	printf("test_bezier_3d\n");

	qaws_vec3 points[] = { {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {1, 1, 1} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_3D;
	desc.degree = 3;
	desc.control_points = points;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(qaws_curve_get_dimension(curve) == QAWS_DIMENSION_3D, "dim == 3D");

	qaws_eval_result_3d r;
	s = qaws_curve_evaluate_3d(curve, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0), "start x");

	s = qaws_curve_evaluate_3d(curve, (qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.z, (qaws_scalar)1.0), "end z");

	qaws_curve_destroy(curve);
}

static void test_hermite(void)
{
	printf("test_hermite\n");

	qaws_vec2 points[] = { {0, 0}, {1, 0} };
	qaws_vec2 tangents[] = { {1, 1}, {1, -1} };

	qaws_hermite_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.points = points;
	desc.derivatives = tangents;
	desc.point_count = 2;
	desc.derivative_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_hermite(&desc, &curve);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(qaws_curve_get_kind(curve) == QAWS_CURVE_KIND_HERMITE, "kind == HERMITE");

	qaws_eval_result_2d r;
	s = qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0), "start x");

	s = qaws_curve_evaluate_2d(curve, (qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)1.0), "end x");
	TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "end y");

	qaws_curve_destroy(curve);
}

static void test_catmull_rom(void)
{
	printf("test_catmull_rom\n");

	qaws_vec2 points[] = { {0, 0}, {1, 1}, {2, 0}, {3, 1} };

	qaws_catmull_rom_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = points;
	desc.control_point_count = 4;
	desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_catmull_rom(&desc, &curve);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(qaws_curve_get_kind(curve) == QAWS_CURVE_KIND_CATMULL_ROM, "kind == CATMULL_ROM");
	TEST_ASSERT(qaws_curve_get_span_count(curve) == 1, "span_count == 1");

	/* The curve interpolates the interior points (P1 and P2) */
	qaws_eval_result_2d r;
	s = qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)1.0), "interpolates P1 x");
	TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)1.0), "interpolates P1 y");

	s = qaws_curve_evaluate_2d(curve, (qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)2.0), "interpolates P2 x");
	TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "interpolates P2 y");

	qaws_curve_destroy(curve);
}

static void test_bspline(void)
{
	printf("test_bspline\n");

	qaws_vec2 points[] = { {0, 0}, {1, 2}, {3, 2}, {4, 0} };
	qaws_bspline_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = points;
	desc.control_point_count = 4;
	desc.knots = NULL;
	desc.knot_count = 0;
	desc.is_uniform = 1;
	desc.is_closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bspline(&desc, &curve);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(qaws_curve_get_kind(curve) == QAWS_CURVE_KIND_BSPLINE, "kind == BSPLINE");

	/* Evaluate at endpoints of clamped uniform B-spline */
	qaws_eval_result_2d r;
	qaws_range range = qaws_curve_get_parameter_range(curve);
	s = qaws_curve_evaluate_2d(curve, range.min_value, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0), "start x");
	TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "start y");

	s = qaws_curve_evaluate_2d(curve, range.max_value, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)4.0), "end x");
	TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "end y");

	qaws_curve_destroy(curve);
}

static void test_nurbs(void)
{
	printf("test_nurbs\n");

	/* Quarter-circle NURBS: degree 2, 3 control points */
	qaws_vec2 points[] = { {1, 0}, {1, 1}, {0, 1} };
	qaws_scalar knots[] = { 0, 0, 0, 1, 1, 1 };
	qaws_scalar weights[] = { 1, (qaws_scalar)0.70710678118, 1 };

	qaws_nurbs_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 2;
	desc.control_points = points;
	desc.control_point_count = 3;
	desc.knots = knots;
	desc.knot_count = 6;
	desc.weights = weights;
	desc.weight_count = 3;
	desc.is_closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_nurbs(&desc, &curve);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(qaws_curve_get_kind(curve) == QAWS_CURVE_KIND_NURBS, "kind == NURBS");
	TEST_ASSERT(qaws_curve_is_rational(curve) == 1, "is_rational");

	qaws_eval_result_2d r;
	s = qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)1.0), "start x");
	TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "start y");

	/* Midpoint of quarter circle should be near (cos(45), sin(45)) */
	s = qaws_curve_evaluate_2d(curve, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	qaws_scalar dist = (qaws_scalar)sqrt((double)(r.position.x * r.position.x + r.position.y * r.position.y));
	TEST_ASSERT(approx_eq(dist, (qaws_scalar)1.0), "midpoint on unit circle");

	qaws_curve_destroy(curve);
}

static void test_trajectory(void)
{
	printf("test_trajectory\n");

	qaws_vec2 positions[] = { {0, 0}, {5, 0}, {5, 5} };
	qaws_scalar times[] = { 0, 1, 2 };

	qaws_trajectory_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.key_positions = positions;
	desc.key_count = 3;
	desc.key_times = times;
	desc.key_time_count = 3;
	desc.key_velocities = NULL;
	desc.key_velocity_count = 0;
	desc.key_accelerations = NULL;
	desc.key_acceleration_count = 0;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_trajectory(&desc, &curve);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(qaws_curve_get_kind(curve) == QAWS_CURVE_KIND_TRAJECTORY, "kind == TRAJECTORY");

	qaws_eval_result_2d r;
	s = qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0), "start x");
	TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "start y");

	s = qaws_curve_evaluate_2d(curve, (qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)5.0), "key1 x");
	TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "key1 y");

	qaws_curve_destroy(curve);
}

static void test_arc_length(void)
{
	printf("test_arc_length\n");

	/* Straight line Bezier: arc length should equal distance */
	qaws_vec2 points[] = { {0, 0}, {3, 4} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_scalar length = 0;
	s = qaws_curve_compute_arc_length(curve, (qaws_scalar)0.0, (qaws_scalar)1.0, &length);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(length, (qaws_scalar)5.0), "straight line length == 5");

	qaws_curve_destroy(curve);
}

static void test_sampling(void)
{
	printf("test_sampling\n");

	qaws_vec2 points[] = { {0, 0}, {1, 1} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_sampling_desc sdesc;
	sdesc.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
	sdesc.sample_count = 11;
	sdesc.error_tolerance = 0;
	sdesc.use_adaptive_sampling = 0;

	qaws_vec2 samples[11];
	unsigned int count = 0;
	s = qaws_curve_sample_2d(curve, &sdesc, samples, 11, &count);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(count == 11, "sample count == 11");
	TEST_ASSERT(approx_eq(samples[0].x, (qaws_scalar)0.0), "first sample x");
	TEST_ASSERT(approx_eq(samples[10].x, (qaws_scalar)1.0), "last sample x");
	TEST_ASSERT(approx_eq(samples[5].x, (qaws_scalar)0.5), "mid sample x");

	qaws_curve_destroy(curve);
}

static void test_traversal(void)
{
	printf("test_traversal\n");

	qaws_vec2 points[] = { {0, 0}, {10, 0} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_traversal_desc tdesc;
	tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
	tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
	tdesc.speed = (qaws_scalar)5.0;
	tdesc.acceleration = 0;
	tdesc.max_speed = 0;
	tdesc.start_time = 0;
	tdesc.end_time = (qaws_scalar)2.0;
	tdesc.clamp_to_domain = 1;

	qaws_traversal* traversal = NULL;
	s = qaws_traversal_create(curve, &tdesc, &traversal);
	TEST_ASSERT_STATUS(s);

	/* At t=0, distance=0, should be at start */
	qaws_eval_result_2d r;
	s = qaws_traversal_evaluate_2d(traversal, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0), "traversal t=0 x");

	/* At t=1, distance=5, should be at x=5 */
	s = qaws_traversal_evaluate_2d(traversal, (qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)5.0), "traversal t=1 x");

	qaws_traversal_destroy(traversal);
	qaws_curve_destroy(curve);
}

static void test_inspect(void)
{
	printf("test_inspect\n");

	qaws_vec2 points[] = { {0, 0}, {1, 2}, {3, 1} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 2;
	desc.control_points = points;
	desc.control_point_count = 3;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	/* Bounds */
	qaws_vec2 bmin, bmax;
	s = qaws_curve_compute_bounds_2d(curve, &bmin, &bmax);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(bmin.x >= (qaws_scalar)-0.1, "bounds min x reasonable");
	TEST_ASSERT(bmax.x <= (qaws_scalar)3.1, "bounds max x reasonable");

	/* Get control points back */
	qaws_vec2 out_pts[3];
	unsigned int out_count = 0;
	s = qaws_bezier_get_control_points(curve, out_pts, 3, &out_count);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(out_count == 3, "got 3 control points");
	TEST_ASSERT(approx_eq(out_pts[0].x, (qaws_scalar)0.0), "cp0 x");

	qaws_curve_destroy(curve);
}

/* ------------------------------------------------------------------ */
/* Derivative value verification                                       */
/* ------------------------------------------------------------------ */

static void test_bezier_derivatives(void)
{
	printf("test_bezier_derivatives\n");

	/* Linear Bezier from (0,0) to (3,4): D1 should be (3,4) everywhere */
	qaws_vec2 pts_lin[] = { {0, 0}, {3, 4} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = pts_lin;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_eval_result_2d r;
	s = qaws_curve_evaluate_2d(curve, (qaws_scalar)0.5, QAWS_EVAL_FLAG_D1, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.d1.x, (qaws_scalar)3.0), "linear D1.x == 3");
	TEST_ASSERT(approx_eq(r.d1.y, (qaws_scalar)4.0), "linear D1.y == 4");

	qaws_curve_destroy(curve);

	/* Cubic Bezier: (0,0),(0,1),(1,1),(1,0) at t=0.5
	   D1 = 3*((0,1)-(0,0)) at t=0 => D1(0) = (0,3)
	   By symmetry, D1(0.5) should have x>0, y=0.
	   D1(t) = 3[(1-t)^2(P1-P0) + 2(1-t)t(P2-P1) + t^2(P3-P2)]
	   D1(0.5) = 3[0.25*(0,1) + 0.5*(1,0) + 0.25*(0,-1)] = 3*(0.5,0) = (1.5,0) */
	qaws_vec2 pts_cub[] = { {0, 0}, {0, 1}, {1, 1}, {1, 0} };
	desc.degree = 3;
	desc.control_points = pts_cub;
	desc.control_point_count = 4;

	s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	s = qaws_curve_evaluate_2d(curve, (qaws_scalar)0.5,
		QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.d1.x, (qaws_scalar)1.5), "cubic D1.x at 0.5");
	TEST_ASSERT(approx_eq(r.d1.y, (qaws_scalar)0.0), "cubic D1.y at 0.5");
	/* D2 by symmetry at midpoint: D2.x = 0, D2.y = -6 */
	TEST_ASSERT(approx_eq(r.d2.x, (qaws_scalar)0.0), "cubic D2.x at 0.5");
	TEST_ASSERT(approx_eq(r.d2.y, (qaws_scalar)-6.0), "cubic D2.y at 0.5");
	TEST_ASSERT(r.valid_flags & QAWS_EVAL_FLAG_D3, "D3 flag set");

	qaws_curve_destroy(curve);
}

/* ------------------------------------------------------------------ */
/* Span-level evaluation                                               */
/* ------------------------------------------------------------------ */

static void test_evaluate_span(void)
{
	printf("test_evaluate_span\n");

	qaws_vec2 points[] = { {0, 0}, {1, 1} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_eval_result_2d r;
	s = qaws_curve_evaluate_span_2d(curve, 0, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.5), "span eval midpoint x");
	TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.5), "span eval midpoint y");

	/* Out of range span index */
	s = qaws_curve_evaluate_span_2d(curve, 5, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT(s == QAWS_STATUS_OUT_OF_RANGE, "out-of-range span returns error");

	qaws_curve_destroy(curve);
}

/* ------------------------------------------------------------------ */
/* Inspection: is_closed, is_periodic, get_continuity                  */
/* ------------------------------------------------------------------ */

static void test_inspection_flags(void)
{
	printf("test_inspection_flags\n");

	qaws_vec2 points[] = { {0, 0}, {1, 1} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	TEST_ASSERT(qaws_curve_is_closed(curve) == 0, "bezier is_closed == 0");
	TEST_ASSERT(qaws_curve_is_periodic(curve) == 0, "bezier is_periodic == 0");
	TEST_ASSERT(qaws_curve_is_rational(curve) == 0, "bezier is_rational == 0");
	TEST_ASSERT(qaws_curve_get_continuity(curve) >= QAWS_CONTINUITY_C0, "bezier continuity >= C0");

	qaws_curve_destroy(curve);
}

/* ------------------------------------------------------------------ */
/* 3D bounds                                                           */
/* ------------------------------------------------------------------ */

static void test_bounds_3d(void)
{
	printf("test_bounds_3d\n");

	qaws_vec3 points[] = { {0, 0, 0}, {1, 2, 3}, {4, 1, 0} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_3D;
	desc.degree = 2;
	desc.control_points = points;
	desc.control_point_count = 3;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_vec3 bmin, bmax;
	s = qaws_curve_compute_bounds_3d(curve, &bmin, &bmax);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(bmin.x >= (qaws_scalar)-0.1, "3d bounds min x");
	TEST_ASSERT(bmax.x <= (qaws_scalar)4.1, "3d bounds max x");
	TEST_ASSERT(bmin.z >= (qaws_scalar)-0.1, "3d bounds min z");
	TEST_ASSERT(bmax.z <= (qaws_scalar)3.1, "3d bounds max z");

	qaws_curve_destroy(curve);
}

/* ------------------------------------------------------------------ */
/* Closest point                                                       */
/* ------------------------------------------------------------------ */

static void test_closest_point(void)
{
	printf("test_closest_point\n");

	/* Straight line from (0,0) to (10,0): closest to (5,3) should be t=0.5 */
	qaws_vec2 points[] = { {0, 0}, {10, 0} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_vec2 target = { 5, 3 };
	qaws_scalar param = 0;
	s = qaws_curve_find_closest_parameter_2d(curve, target, &param);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(param, (qaws_scalar)0.5), "closest param == 0.5");

	qaws_curve_destroy(curve);
}

/* ------------------------------------------------------------------ */
/* Span continuity                                                     */
/* ------------------------------------------------------------------ */

static void test_span_continuity(void)
{
	printf("test_span_continuity\n");

	/* Multi-span Hermite: 3 points => 2 spans => 1 interior boundary */
	qaws_vec2 pts[] = { {0, 0}, {1, 1}, {2, 0} };
	qaws_vec2 tans[] = { {1, 0}, {0, 0}, {1, 0} };

	qaws_hermite_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.points = pts;
	desc.derivatives = tans;
	desc.point_count = 3;
	desc.derivative_count = 3;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_hermite(&desc, &curve);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(qaws_curve_get_span_count(curve) == 2, "hermite 3pt span_count == 2");

	qaws_continuity cont;
	s = qaws_curve_get_span_continuity(curve, 0, &cont);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(cont >= QAWS_CONTINUITY_C0, "span continuity >= C0");

	qaws_curve_destroy(curve);
}

/* ------------------------------------------------------------------ */
/* Sampling: get_sample_count, sample_3d, sample_eval_2d               */
/* ------------------------------------------------------------------ */

static void test_sampling_extended(void)
{
	printf("test_sampling_extended\n");

	/* 2D: get_sample_count */
	qaws_vec2 points2d[] = { {0, 0}, {1, 1} };
	qaws_bezier_desc desc2;
	desc2.dimension = QAWS_DIMENSION_2D;
	desc2.degree = 1;
	desc2.control_points = points2d;
	desc2.control_point_count = 2;

	qaws_curve* curve2 = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc2, &curve2);
	TEST_ASSERT_STATUS(s);

	qaws_sampling_desc sdesc;
	sdesc.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
	sdesc.sample_count = 5;
	sdesc.error_tolerance = 0;
	sdesc.use_adaptive_sampling = 0;

	unsigned int count = 0;
	s = qaws_curve_get_sample_count(curve2, &sdesc, &count);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(count == 5, "get_sample_count == 5");

	/* sample_eval_2d */
	qaws_eval_result_2d evals[5];
	unsigned int eval_count = 0;
	s = qaws_curve_sample_eval_2d(curve2, &sdesc, evals, 5, &eval_count);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(eval_count == 5, "sample_eval count == 5");
	TEST_ASSERT(approx_eq(evals[0].position.x, (qaws_scalar)0.0), "eval sample 0 x");
	TEST_ASSERT(approx_eq(evals[4].position.x, (qaws_scalar)1.0), "eval sample 4 x");

	qaws_curve_destroy(curve2);

	/* 3D: sample_3d */
	qaws_vec3 points3d[] = { {0, 0, 0}, {1, 1, 1} };
	qaws_bezier_desc desc3;
	desc3.dimension = QAWS_DIMENSION_3D;
	desc3.degree = 1;
	desc3.control_points = points3d;
	desc3.control_point_count = 2;

	qaws_curve* curve3 = NULL;
	s = qaws_curve_create_bezier(&desc3, &curve3);
	TEST_ASSERT_STATUS(s);

	qaws_vec3 samples3[5];
	unsigned int count3 = 0;
	s = qaws_curve_sample_3d(curve3, &sdesc, samples3, 5, &count3);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(count3 == 5, "sample_3d count == 5");
	TEST_ASSERT(approx_eq(samples3[0].z, (qaws_scalar)0.0), "3d sample start z");
	TEST_ASSERT(approx_eq(samples3[4].z, (qaws_scalar)1.0), "3d sample end z");

	/* sample_eval_3d */
	qaws_eval_result_3d evals3[5];
	unsigned int eval3_count = 0;
	s = qaws_curve_sample_eval_3d(curve3, &sdesc, evals3, 5, &eval3_count);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(eval3_count == 5, "sample_eval_3d count == 5");

	qaws_curve_destroy(curve3);
}

/* ------------------------------------------------------------------ */
/* Traversal: mapping queries, arc-length mode, 3D eval                */
/* ------------------------------------------------------------------ */

static void test_traversal_mappings(void)
{
	printf("test_traversal_mappings\n");

	/* Straight line (0,0) to (10,0), length=10 */
	qaws_vec2 points[] = { {0, 0}, {10, 0} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_traversal_desc tdesc;
	tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
	tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
	tdesc.speed = (qaws_scalar)5.0;
	tdesc.acceleration = 0;
	tdesc.max_speed = 0;
	tdesc.start_time = 0;
	tdesc.end_time = (qaws_scalar)2.0;
	tdesc.clamp_to_domain = 1;

	qaws_traversal* traversal = NULL;
	s = qaws_traversal_create(curve, &tdesc, &traversal);
	TEST_ASSERT_STATUS(s);

	/* map_time_to_parameter: at t=1, distance=5, on length-10 line => param ~0.5 */
	qaws_scalar param = 0;
	s = qaws_traversal_map_time_to_parameter(traversal, (qaws_scalar)1.0, &param);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(param, (qaws_scalar)0.5), "time_to_param t=1 => 0.5");

	/* map_distance_to_parameter: distance 5 => param 0.5 */
	s = qaws_traversal_map_distance_to_parameter(traversal, (qaws_scalar)5.0, &param);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(param, (qaws_scalar)0.5), "dist_to_param d=5 => 0.5");

	/* map_parameter_to_distance: param 0.5 => distance 5 */
	qaws_scalar dist = 0;
	s = qaws_traversal_map_parameter_to_distance(traversal, (qaws_scalar)0.5, &dist);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(dist, (qaws_scalar)5.0), "param_to_dist p=0.5 => 5");

	qaws_traversal_destroy(traversal);
	qaws_curve_destroy(curve);
}

static void test_traversal_arc_length_mode(void)
{
	printf("test_traversal_arc_length_mode\n");

	/* Arc-length traversal: straight line, distance 5 should map to midpoint */
	qaws_vec2 points[] = { {0, 0}, {10, 0} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_traversal_desc tdesc;
	tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
	tdesc.motion_profile = QAWS_MOTION_PROFILE_NONE;
	tdesc.speed = 0;
	tdesc.acceleration = 0;
	tdesc.max_speed = 0;
	tdesc.start_time = 0;
	tdesc.end_time = 0;
	tdesc.clamp_to_domain = 1;

	qaws_traversal* traversal = NULL;
	s = qaws_traversal_create(curve, &tdesc, &traversal);
	TEST_ASSERT_STATUS(s);

	/* input_value is distance: 5.0 on a 10-length line => x=5 */
	qaws_eval_result_2d r;
	s = qaws_traversal_evaluate_2d(traversal, (qaws_scalar)5.0, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)5.0), "arc-length eval at d=5 => x=5");

	qaws_traversal_destroy(traversal);
	qaws_curve_destroy(curve);
}

static void test_traversal_3d(void)
{
	printf("test_traversal_3d\n");

	qaws_vec3 points[] = { {0, 0, 0}, {10, 0, 0} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_3D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_traversal_desc tdesc;
	tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
	tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
	tdesc.speed = (qaws_scalar)10.0;
	tdesc.acceleration = 0;
	tdesc.max_speed = 0;
	tdesc.start_time = 0;
	tdesc.end_time = (qaws_scalar)1.0;
	tdesc.clamp_to_domain = 1;

	qaws_traversal* traversal = NULL;
	s = qaws_traversal_create(curve, &tdesc, &traversal);
	TEST_ASSERT_STATUS(s);

	qaws_eval_result_3d r;
	s = qaws_traversal_evaluate_3d(traversal, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)5.0), "3d traversal t=0.5 x");

	qaws_traversal_destroy(traversal);
	qaws_curve_destroy(curve);
}

/* ------------------------------------------------------------------ */
/* Family-specific inspection: bspline_get_knots, nurbs_get_weights    */
/* ------------------------------------------------------------------ */

static void test_family_inspection(void)
{
	printf("test_family_inspection\n");

	/* bspline_get_knots */
	{
		qaws_vec2 points[] = { {0, 0}, {1, 2}, {3, 2}, {4, 0} };
		qaws_bspline_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = points;
		desc.control_point_count = 4;
		desc.knots = NULL;
		desc.knot_count = 0;
		desc.is_uniform = 1;
		desc.is_closed = 0;

		qaws_curve* curve = NULL;
		qaws_status s = qaws_curve_create_bspline(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		qaws_scalar knots[16];
		unsigned int knot_count = 0;
		s = qaws_bspline_get_knots(curve, knots, 16, &knot_count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(knot_count > 0, "bspline knot_count > 0");
		/* Clamped: first degree+1 knots should be equal */
		TEST_ASSERT(approx_eq(knots[0], knots[1]), "clamped knots start equal");

		/* Wrong family: bezier_get_control_points on bspline */
		qaws_vec2 dummy[4];
		unsigned int dummy_count = 0;
		s = qaws_bezier_get_control_points(curve, dummy, 4, &dummy_count);
		TEST_ASSERT(s == QAWS_STATUS_UNSUPPORTED_OPERATION, "bezier_get_cp on bspline => UNSUPPORTED");

		qaws_curve_destroy(curve);
	}

	/* nurbs_get_weights */
	{
		qaws_vec2 points[] = { {1, 0}, {1, 1}, {0, 1} };
		qaws_scalar knots[] = { 0, 0, 0, 1, 1, 1 };
		qaws_scalar weights[] = { 1, (qaws_scalar)0.707, 1 };

		qaws_nurbs_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = points;
		desc.control_point_count = 3;
		desc.knots = knots;
		desc.knot_count = 6;
		desc.weights = weights;
		desc.weight_count = 3;
		desc.is_closed = 0;

		qaws_curve* curve = NULL;
		qaws_status s = qaws_curve_create_nurbs(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		qaws_scalar out_w[3];
		unsigned int w_count = 0;
		s = qaws_nurbs_get_weights(curve, out_w, 3, &w_count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(w_count == 3, "nurbs weight_count == 3");
		TEST_ASSERT(approx_eq(out_w[0], (qaws_scalar)1.0), "weight 0 == 1");

		qaws_curve_destroy(curve);
	}
}

/* ------------------------------------------------------------------ */
/* Error handling / validation                                         */
/* ------------------------------------------------------------------ */

static void test_error_handling(void)
{
	printf("test_error_handling\n");

	qaws_status s;

	/* NULL arguments */
	s = qaws_curve_create_bezier(NULL, NULL);
	TEST_ASSERT(s != QAWS_STATUS_OK, "NULL desc rejected");

	qaws_curve* curve = NULL;
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = NULL;
	desc.control_point_count = 2;
	s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT(s != QAWS_STATUS_OK, "NULL control_points rejected");

	/* Invalid degree */
	qaws_vec2 points[] = { {0, 0} };
	desc.control_points = points;
	desc.control_point_count = 1;
	desc.degree = 0;
	s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT(s != QAWS_STATUS_OK, "degree 0 rejected");

	/* Wrong CP count */
	desc.degree = 1;
	desc.control_point_count = 5;
	s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT(s != QAWS_STATUS_OK, "wrong CP count rejected");

	/* Invalid dimension */
	desc.dimension = (qaws_dimension)99;
	desc.degree = 1;
	desc.control_point_count = 2;
	qaws_vec2 pts2[] = { {0,0}, {1,1} };
	desc.control_points = pts2;
	s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT(s != QAWS_STATUS_OK, "invalid dimension rejected");

	/* Evaluate NULL curve */
	qaws_eval_result_2d r;
	s = qaws_curve_evaluate_2d(NULL, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT(s != QAWS_STATUS_OK, "NULL curve evaluate rejected");

	/* Wrong dimension eval */
	desc.dimension = QAWS_DIMENSION_2D;
	s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_eval_result_3d r3;
	s = qaws_curve_evaluate_3d(curve, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &r3);
	TEST_ASSERT(s == QAWS_STATUS_INVALID_DIMENSION, "2D curve + evaluate_3d => INVALID_DIMENSION");

	qaws_curve_destroy(curve);
}

/* ------------------------------------------------------------------ */
/* Catmull-Rom parameterization variants                               */
/* ------------------------------------------------------------------ */

static void test_catmull_rom_variants(void)
{
	printf("test_catmull_rom_variants\n");

	qaws_vec2 points[] = { {0, 0}, {1, 1}, {2, 0}, {3, 1} };

	/* Uniform */
	{
		qaws_catmull_rom_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = points;
		desc.control_point_count = 4;
		desc.parameterization = QAWS_PARAMETERIZATION_UNIFORM;
		desc.closed = 0;

		qaws_curve* curve = NULL;
		qaws_status s = qaws_curve_create_catmull_rom(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		qaws_eval_result_2d r;
		s = qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)1.0), "uniform interpolates P1 x");

		qaws_curve_destroy(curve);
	}

	/* Chordal */
	{
		qaws_catmull_rom_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = points;
		desc.control_point_count = 4;
		desc.parameterization = QAWS_PARAMETERIZATION_CHORDAL;
		desc.closed = 0;

		qaws_curve* curve = NULL;
		qaws_status s = qaws_curve_create_catmull_rom(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		qaws_eval_result_2d r;
		s = qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)1.0), "chordal interpolates P1 x");

		qaws_curve_destroy(curve);
	}
}

/* ------------------------------------------------------------------ */
/* B-Spline with custom knots                                          */
/* ------------------------------------------------------------------ */

static void test_bspline_custom_knots(void)
{
	printf("test_bspline_custom_knots\n");

	qaws_vec2 points[] = { {0, 0}, {1, 2}, {3, 2}, {4, 0} };
	qaws_scalar knots[] = { 0, 0, 0, 0, 1, 1, 1, 1 }; /* clamped cubic, 4 CPs */

	qaws_bspline_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = points;
	desc.control_point_count = 4;
	desc.knots = knots;
	desc.knot_count = 8;
	desc.is_uniform = 0;
	desc.is_closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bspline(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_eval_result_2d r;
	qaws_range range = qaws_curve_get_parameter_range(curve);
	s = qaws_curve_evaluate_2d(curve, range.min_value, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0), "custom knots start x");

	s = qaws_curve_evaluate_2d(curve, range.max_value, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)4.0), "custom knots end x");

	qaws_curve_destroy(curve);
}

/* ------------------------------------------------------------------ */
/* 3D family tests                                                     */
/* ------------------------------------------------------------------ */

static void test_hermite_3d(void)
{
	printf("test_hermite_3d\n");

	qaws_vec3 points[] = { {0, 0, 0}, {1, 0, 0} };
	qaws_vec3 tangents[] = { {1, 0, 1}, {1, 0, -1} };

	qaws_hermite_desc desc;
	desc.dimension = QAWS_DIMENSION_3D;
	desc.degree = 3;
	desc.points = points;
	desc.derivatives = tangents;
	desc.point_count = 2;
	desc.derivative_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_hermite(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_eval_result_3d r;
	s = qaws_curve_evaluate_3d(curve, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0), "hermite 3d start x");
	TEST_ASSERT(approx_eq(r.position.z, (qaws_scalar)0.0), "hermite 3d start z");

	s = qaws_curve_evaluate_3d(curve, (qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)1.0), "hermite 3d end x");

	qaws_curve_destroy(curve);
}

/* ================================================================== */
/* Trivial SVG writer                                                  */
/* ================================================================== */

typedef struct svg_writer
{
	FILE* fp;
	qaws_scalar view_x, view_y, view_w, view_h;
	qaws_scalar svg_w, svg_h;
	qaws_scalar padding;
} svg_writer;

static int svg_open(svg_writer* w, char const* path,
	qaws_scalar view_x, qaws_scalar view_y,
	qaws_scalar view_w, qaws_scalar view_h,
	qaws_scalar svg_w, qaws_scalar svg_h)
{
	w->fp = fopen(path, "w");
	if (!w->fp) return 0;
	w->view_x = view_x; w->view_y = view_y;
	w->view_w = view_w; w->view_h = view_h;
	w->svg_w = svg_w;   w->svg_h = svg_h;
	w->padding = (qaws_scalar)20;
	fprintf(w->fp,
		"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
		"<svg xmlns=\"http://www.w3.org/2000/svg\" "
		"width=\"%.0f\" height=\"%.0f\" "
		"viewBox=\"0 0 %.0f %.0f\" "
		"style=\"background:#1a1a2e\">\n",
		(double)svg_w, (double)svg_h, (double)svg_w, (double)svg_h);
	/* Grid */
	fprintf(w->fp,
		"<defs><pattern id=\"grid\" width=\"%.4f\" height=\"%.4f\" patternUnits=\"userSpaceOnUse\">"
		"<path d=\"M %.4f 0 L 0 0 0 %.4f\" fill=\"none\" stroke=\"#2a2a4a\" stroke-width=\"0.5\"/>"
		"</pattern></defs>\n"
		"<rect width=\"100%%\" height=\"100%%\" fill=\"url(#grid)\"/>\n",
		(double)((svg_w - 2 * w->padding) / view_w),
		(double)((svg_h - 2 * w->padding) / view_h),
		(double)((svg_w - 2 * w->padding) / view_w),
		(double)((svg_h - 2 * w->padding) / view_h));
	/* Axes */
	{
		qaws_scalar scale_x = (svg_w - 2 * w->padding) / view_w;
		qaws_scalar scale_y = (svg_h - 2 * w->padding) / view_h;
		qaws_scalar ox = w->padding + (0 - view_x) * scale_x;
		qaws_scalar oy = svg_h - w->padding - (0 - view_y) * scale_y;
		fprintf(w->fp,
			"<line x1=\"%.2f\" y1=\"0\" x2=\"%.2f\" y2=\"%.0f\" stroke=\"#3a3a5a\" stroke-width=\"1\"/>\n"
			"<line x1=\"0\" y1=\"%.2f\" x2=\"%.0f\" y2=\"%.2f\" stroke=\"#3a3a5a\" stroke-width=\"1\"/>\n",
			(double)ox, (double)ox, (double)svg_h,
			(double)oy, (double)svg_w, (double)oy);
	}
	return 1;
}

static void svg_map(svg_writer const* w, qaws_scalar wx, qaws_scalar wy,
	qaws_scalar* sx, qaws_scalar* sy)
{
	qaws_scalar scale_x = (w->svg_w - 2 * w->padding) / w->view_w;
	qaws_scalar scale_y = (w->svg_h - 2 * w->padding) / w->view_h;
	*sx = w->padding + (wx - w->view_x) * scale_x;
	*sy = w->svg_h - w->padding - (wy - w->view_y) * scale_y;
}

static void svg_polyline(svg_writer* w, qaws_vec2 const* pts, unsigned int count,
	char const* stroke, qaws_scalar stroke_width)
{
	if (count < 2) return;
	fprintf(w->fp, "<polyline points=\"");
	for (unsigned int i = 0; i < count; i++) {
		qaws_scalar sx, sy;
		svg_map(w, pts[i].x, pts[i].y, &sx, &sy);
		fprintf(w->fp, "%.4f,%.4f ", (double)sx, (double)sy);
	}
	fprintf(w->fp, "\" fill=\"none\" stroke=\"%s\" stroke-width=\"%.1f\" "
		"stroke-linecap=\"round\" stroke-linejoin=\"round\"/>\n",
		stroke, (double)stroke_width);
}

static void svg_dots(svg_writer* w, qaws_vec2 const* pts, unsigned int count,
	char const* fill, qaws_scalar radius)
{
	for (unsigned int i = 0; i < count; i++) {
		qaws_scalar sx, sy;
		svg_map(w, pts[i].x, pts[i].y, &sx, &sy);
		fprintf(w->fp, "<circle cx=\"%.4f\" cy=\"%.4f\" r=\"%.1f\" fill=\"%s\"/>\n",
			(double)sx, (double)sy, (double)radius, fill);
	}
}

static void svg_label(svg_writer* w, qaws_scalar wx, qaws_scalar wy,
	char const* text, char const* fill)
{
	qaws_scalar sx, sy;
	svg_map(w, wx, wy, &sx, &sy);
	fprintf(w->fp,
		"<text x=\"%.2f\" y=\"%.2f\" font-family=\"monospace\" font-size=\"13\" "
		"fill=\"%s\">%s</text>\n",
		(double)sx, (double)(sy - 8), fill, text);
}

static void svg_close(svg_writer* w)
{
	fprintf(w->fp, "</svg>\n");
	fclose(w->fp);
	w->fp = NULL;
}

/* Helper: sample a 2D curve into a buffer */
static unsigned int svg_sample_curve(qaws_curve const* curve,
	qaws_vec2* buf, unsigned int capacity)
{
	qaws_sampling_desc sd;
	sd.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
	sd.sample_count = capacity;
	sd.error_tolerance = 0;
	sd.use_adaptive_sampling = 0;
	unsigned int count = 0;
	qaws_curve_sample_2d(curve, &sd, buf, capacity, &count);
	return count;
}

/* ================================================================== */
/* Visual SVG tests                                                    */
/* ================================================================== */

static void svg_ensure_output_dir(void)
{
	MKDIR("working_dir");
	MKDIR(SVG_OUTPUT_DIR);
}

static void test_svg_circle(void)
{
	printf("test_svg_circle\n");

	/* Full circle via 4 quarter-arc NURBS patches joined together.
	   We'll draw each quarter as a separate NURBS. */
	qaws_scalar w = (qaws_scalar)0.70710678118; /* 1/sqrt(2) */
	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/circle.svg",
		(qaws_scalar)-1.5, (qaws_scalar)-1.5, (qaws_scalar)3.0, (qaws_scalar)3.0,
		(qaws_scalar)500, (qaws_scalar)500))
		return;

	char const* colors[] = { "#e94560", "#f5a623", "#50fa7b", "#6272a4" };

	/* Quarter arcs: Q0 (1,0)->(0,1), Q1 (0,1)->(-1,0), Q2 (-1,0)->(0,-1), Q3 (0,-1)->(1,0) */
	qaws_vec2 q_cp[4][3] = {
		{ {1, 0}, {1, 1}, {0, 1} },
		{ {0, 1}, {-1, 1}, {-1, 0} },
		{ {-1, 0}, {-1, -1}, {0, -1} },
		{ {0, -1}, {1, -1}, {1, 0} }
	};
	qaws_scalar q_w[] = { 1, w, 1 };
	qaws_scalar knots[] = { 0, 0, 0, 1, 1, 1 };

	qaws_vec2 all_pts[SVG_SAMPLES * 4];
	unsigned int total = 0;

	for (int qi = 0; qi < 4; qi++) {
		qaws_nurbs_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = q_cp[qi];
		desc.control_point_count = 3;
		desc.knots = knots;
		desc.knot_count = 6;
		desc.weights = q_w;
		desc.weight_count = 3;
		desc.is_closed = 0;

		qaws_curve* curve = NULL;
		qaws_status s = qaws_curve_create_nurbs(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, colors[qi], (qaws_scalar)2.5);
		for (unsigned int i = 0; i < n; i++) all_pts[total++] = buf[i];
		svg_dots(&svg, q_cp[qi], 3, "#ffffff", (qaws_scalar)3);
		qaws_curve_destroy(curve);
	}

	/* Verify: check that sampled points are on the unit circle */
	int on_circle = 1;
	for (unsigned int i = 0; i < total; i++) {
		qaws_scalar r = (qaws_scalar)sqrt(
			(double)(all_pts[i].x * all_pts[i].x + all_pts[i].y * all_pts[i].y));
		if (r < (qaws_scalar)0.99 || r > (qaws_scalar)1.01) { on_circle = 0; break; }
	}
	TEST_ASSERT(on_circle, "circle: all samples on unit circle");

	svg_label(&svg, (qaws_scalar)-1.3, (qaws_scalar)1.3, "Circle (4x NURBS quarter-arcs)", "#8888aa");
	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/circle.svg\n");
}

static void test_svg_lemniscate(void)
{
	printf("test_svg_lemniscate\n");

	/* Lemniscate of Bernoulli (infinity symbol) approximated as a Bezier.
	   x(t) = cos(t) / (1 + sin^2(t))
	   y(t) = sin(t)*cos(t) / (1 + sin^2(t))
	   We sample the parametric form directly, then fit Bezier control points
	   for two lobes and draw them. Simpler: just build a Bezier curve that
	   approximates each half, using sampled points as guide. */

	/* Approach: we'll manually sample the lemniscate formula and create
	   a piecewise Bezier (degree-3) approximation by using Hermite
	   interpolation on several key points per lobe. */
	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/lemniscate.svg",
		(qaws_scalar)-1.5, (qaws_scalar)-0.8, (qaws_scalar)3.0, (qaws_scalar)1.6,
		(qaws_scalar)600, (qaws_scalar)320))
		return;

	/* Parametric lemniscate: sample N points over t=[0, 2*pi] */
	#define LEMN_N 512
	qaws_vec2 lemn_pts[LEMN_N];
	{
		double pi2 = 6.283185307179586;
		for (int i = 0; i < LEMN_N; i++) {
			double t = pi2 * (double)i / (double)(LEMN_N - 1);
			double s = sin(t);
			double c = cos(t);
			double denom = 1.0 + s * s;
			lemn_pts[i].x = (qaws_scalar)(c / denom);
			lemn_pts[i].y = (qaws_scalar)(s * c / denom);
		}
	}

	/* Draw the raw parametric lemniscate as a polyline */
	svg_polyline(&svg, lemn_pts, LEMN_N, "#50fa7b", (qaws_scalar)1.0);

	/* Now create a Bezier approximation of the right lobe using Hermite.
	   Right lobe: from (1,0) up through top, back to (0,0), down to (1,0).
	   We'll use two Hermite spans. */
	{
		/* Right lobe: 3 key points with tangents */
		qaws_vec2 pts[] = {
			{(qaws_scalar)0.0, (qaws_scalar)0.0},
			{(qaws_scalar)1.0, (qaws_scalar)0.0},  /* rightmost */
			{(qaws_scalar)0.0, (qaws_scalar)0.0}
		};
		qaws_vec2 tans[] = {
			{(qaws_scalar)1.0, (qaws_scalar)1.0},
			{(qaws_scalar)0.0, (qaws_scalar)-1.5},
			{(qaws_scalar)-1.0, (qaws_scalar)1.0}
		};

		qaws_hermite_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.points = pts;
		desc.derivatives = tans;
		desc.point_count = 3;
		desc.derivative_count = 3;

		qaws_curve* curve = NULL;
		qaws_status s = qaws_curve_create_hermite(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#e94560", (qaws_scalar)2.5);
		svg_dots(&svg, pts, 3, "#ffffff", (qaws_scalar)3);
		qaws_curve_destroy(curve);
	}

	/* Left lobe */
	{
		qaws_vec2 pts[] = {
			{(qaws_scalar)0.0, (qaws_scalar)0.0},
			{(qaws_scalar)-1.0, (qaws_scalar)0.0},
			{(qaws_scalar)0.0, (qaws_scalar)0.0}
		};
		qaws_vec2 tans[] = {
			{(qaws_scalar)-1.0, (qaws_scalar)-1.0},
			{(qaws_scalar)0.0, (qaws_scalar)1.5},
			{(qaws_scalar)1.0, (qaws_scalar)-1.0}
		};

		qaws_hermite_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.points = pts;
		desc.derivatives = tans;
		desc.point_count = 3;
		desc.derivative_count = 3;

		qaws_curve* curve = NULL;
		qaws_status s = qaws_curve_create_hermite(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#f5a623", (qaws_scalar)2.5);
		qaws_curve_destroy(curve);
	}

	svg_label(&svg, (qaws_scalar)-1.3, (qaws_scalar)0.65,
		"Lemniscate (Hermite lobes + parametric ref)", "#8888aa");
	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/lemniscate.svg\n");
}

static void test_svg_bezier_s(void)
{
	printf("test_svg_bezier_s\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/bezier_s_curve.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-0.5, (qaws_scalar)5.0, (qaws_scalar)4.0,
		(qaws_scalar)500, (qaws_scalar)400))
		return;

	/* S-shaped cubic Bezier */
	qaws_vec2 cp[] = { {0, 0}, {0, 3}, {4, 0}, {4, 3} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);

	/* Draw control polygon (dashed) */
	svg_polyline(&svg, cp, 4, "#444466", (qaws_scalar)1);
	/* Draw curve */
	svg_polyline(&svg, buf, n, "#ff6188", (qaws_scalar)2.5);
	/* Draw control points */
	svg_dots(&svg, cp, 4, "#a9dc76", (qaws_scalar)4);

	/* Verify endpoints */
	TEST_ASSERT(approx_eq(buf[0].x, (qaws_scalar)0.0), "s-curve start x");
	TEST_ASSERT(approx_eq(buf[n - 1].x, (qaws_scalar)4.0), "s-curve end x");

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.5, "Cubic Bezier S-curve", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/bezier_s_curve.svg\n");
}

static void test_svg_catmull_rom_wave(void)
{
	printf("test_svg_catmull_rom_wave\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/catmull_rom_wave.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-2.0, (qaws_scalar)8.0, (qaws_scalar)4.0,
		(qaws_scalar)600, (qaws_scalar)300))
		return;

	/* Wave-like Catmull-Rom through 6 control points */
	qaws_vec2 all_cp[] = {
		{0, 0}, {1, (qaws_scalar)1.5}, {2, (qaws_scalar)-1.5},
		{3, (qaws_scalar)1.5}, {4, (qaws_scalar)-1.5}, {5, (qaws_scalar)1.5},
		{6, 0}
	};

	/* Catmull-Rom needs at least 4 points for 1 interior span.
	   With 7 control points we get 5 interior spans (P1..P5). */
	qaws_catmull_rom_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = all_cp;
	desc.control_point_count = 7;
	desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_catmull_rom(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
	svg_polyline(&svg, buf, n, "#78dce8", (qaws_scalar)2.5);

	/* Show interpolated points (interior: indices 1..5) */
	svg_dots(&svg, all_cp + 1, 5, "#ff6188", (qaws_scalar)4);
	/* Show ghost points */
	svg_dots(&svg, all_cp, 1, "#666688", (qaws_scalar)3);
	svg_dots(&svg, all_cp + 6, 1, "#666688", (qaws_scalar)3);

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)1.7, "Catmull-Rom wave (centripetal)", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/catmull_rom_wave.svg\n");
}

static void test_svg_bspline(void)
{
	printf("test_svg_bspline\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/bspline.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-0.5, (qaws_scalar)5.0, (qaws_scalar)4.0,
		(qaws_scalar)500, (qaws_scalar)400))
		return;

	/* Cubic B-Spline with 6 control points */
	qaws_vec2 cp[] = {
		{0, 0}, {1, 3}, {2, 0}, {3, 3}, {4, 0}, {4, 3}
	};
	qaws_bspline_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 6;
	desc.knots = NULL;
	desc.knot_count = 0;
	desc.is_uniform = 1;
	desc.is_closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bspline(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);

	/* Control polygon */
	svg_polyline(&svg, cp, 6, "#444466", (qaws_scalar)1);
	/* Curve */
	svg_polyline(&svg, buf, n, "#ab9df2", (qaws_scalar)2.5);
	/* Control points */
	svg_dots(&svg, cp, 6, "#a9dc76", (qaws_scalar)4);

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.5, "Cubic B-Spline (6 CPs, uniform clamped)", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/bspline.svg\n");
}

static void test_svg_nurbs_ellipse(void)
{
	printf("test_svg_nurbs_ellipse\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/nurbs_ellipse.svg",
		(qaws_scalar)-3.0, (qaws_scalar)-2.0, (qaws_scalar)6.0, (qaws_scalar)4.0,
		(qaws_scalar)500, (qaws_scalar)340))
		return;

	/* Ellipse: scale unit circle x by 2 */
	qaws_scalar w = (qaws_scalar)0.70710678118;
	char const* colors[] = { "#e94560", "#f5a623", "#50fa7b", "#6272a4" };

	qaws_vec2 q_cp[4][3] = {
		{ {2, 0}, {2, 1}, {0, 1} },
		{ {0, 1}, {-2, 1}, {-2, 0} },
		{ {-2, 0}, {-2, -1}, {0, -1} },
		{ {0, -1}, {2, -1}, {2, 0} }
	};
	qaws_scalar q_w[] = { 1, w, 1 };
	qaws_scalar knots[] = { 0, 0, 0, 1, 1, 1 };

	for (int qi = 0; qi < 4; qi++) {
		qaws_nurbs_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = q_cp[qi];
		desc.control_point_count = 3;
		desc.knots = knots;
		desc.knot_count = 6;
		desc.weights = q_w;
		desc.weight_count = 3;
		desc.is_closed = 0;

		qaws_curve* curve = NULL;
		qaws_status s = qaws_curve_create_nurbs(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, colors[qi], (qaws_scalar)2.5);
		svg_dots(&svg, q_cp[qi], 3, "#ffffff80", (qaws_scalar)2);
		qaws_curve_destroy(curve);
	}

	svg_label(&svg, (qaws_scalar)-2.8, (qaws_scalar)1.7, "Ellipse (4x NURBS, a=2, b=1)", "#8888aa");
	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/nurbs_ellipse.svg\n");
}

static void test_svg_hermite_figure8(void)
{
	printf("test_svg_hermite_figure8\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/hermite_figure8.svg",
		(qaws_scalar)-1.5, (qaws_scalar)-1.5, (qaws_scalar)3.0, (qaws_scalar)3.0,
		(qaws_scalar)500, (qaws_scalar)500))
		return;

	/* Figure-8 using 5-point Hermite loop:
	   origin -> right lobe top -> origin -> left lobe bottom -> origin */
	qaws_vec2 pts[] = {
		{(qaws_scalar)0.0,  (qaws_scalar)0.0},
		{(qaws_scalar)1.0,  (qaws_scalar)0.5},
		{(qaws_scalar)0.0,  (qaws_scalar)0.0},
		{(qaws_scalar)-1.0, (qaws_scalar)-0.5},
		{(qaws_scalar)0.0,  (qaws_scalar)0.0}
	};
	qaws_vec2 tans[] = {
		{(qaws_scalar)1.0,  (qaws_scalar)1.0},
		{(qaws_scalar)0.0,  (qaws_scalar)-1.5},
		{(qaws_scalar)-1.0, (qaws_scalar)-1.0},
		{(qaws_scalar)0.0,  (qaws_scalar)1.5},
		{(qaws_scalar)1.0,  (qaws_scalar)1.0}
	};

	qaws_hermite_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.points = pts;
	desc.derivatives = tans;
	desc.point_count = 5;
	desc.derivative_count = 5;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_hermite(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
	svg_polyline(&svg, buf, n, "#ff6188", (qaws_scalar)2.5);
	svg_dots(&svg, pts, 5, "#a9dc76", (qaws_scalar)4);

	/* Verify: starts and ends at origin */
	TEST_ASSERT(approx_eq(buf[0].x, (qaws_scalar)0.0), "figure8 start x");
	TEST_ASSERT(approx_eq(buf[n - 1].x, (qaws_scalar)0.0), "figure8 end x");

	svg_label(&svg, (qaws_scalar)-1.3, (qaws_scalar)1.3, "Figure-8 (Hermite, 5 key points)", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/hermite_figure8.svg\n");
}

static void test_svg_all_families(void)
{
	printf("test_svg_all_families\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/all_families.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-1.5, (qaws_scalar)5.0, (qaws_scalar)5.0,
		(qaws_scalar)500, (qaws_scalar)500))
		return;

	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n;

	/* 1. Bezier cubic */
	{
		qaws_vec2 cp[] = { {0, 0}, {0, 3}, {4, 0}, {4, 3} };
		qaws_bezier_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = cp;
		desc.control_point_count = 4;

		qaws_curve* curve = NULL;
		qaws_curve_create_bezier(&desc, &curve);
		n = svg_sample_curve(curve, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#ff6188", (qaws_scalar)2);
		qaws_curve_destroy(curve);
	}

	/* 2. Hermite */
	{
		qaws_vec2 pts[] = { {0, (qaws_scalar)-0.5}, {4, (qaws_scalar)-0.5} };
		qaws_vec2 tans[] = { {2, 2}, {2, -2} };
		qaws_hermite_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.points = pts;
		desc.derivatives = tans;
		desc.point_count = 2;
		desc.derivative_count = 2;

		qaws_curve* curve = NULL;
		qaws_curve_create_hermite(&desc, &curve);
		n = svg_sample_curve(curve, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#f5a623", (qaws_scalar)2);
		qaws_curve_destroy(curve);
	}

	/* 3. B-Spline */
	{
		qaws_vec2 cp[] = { {0, (qaws_scalar)-1}, {1, 2}, {2, (qaws_scalar)-1}, {3, 2}, {4, (qaws_scalar)-1} };
		qaws_bspline_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = cp;
		desc.control_point_count = 5;
		desc.knots = NULL;
		desc.knot_count = 0;
		desc.is_uniform = 1;
		desc.is_closed = 0;

		qaws_curve* curve = NULL;
		qaws_curve_create_bspline(&desc, &curve);
		n = svg_sample_curve(curve, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#78dce8", (qaws_scalar)2);
		qaws_curve_destroy(curve);
	}

	/* Legend */
	svg_label(&svg, (qaws_scalar)0.0, (qaws_scalar)3.3, "Bezier", "#ff6188");
	svg_label(&svg, (qaws_scalar)0.0, (qaws_scalar)3.0, "Hermite", "#f5a623");
	svg_label(&svg, (qaws_scalar)0.0, (qaws_scalar)2.7, "B-Spline", "#78dce8");

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.7, "All families comparison", "#8888aa");
	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/all_families.svg\n");
}

/* ================================================================== */
/* Yuksel C2 Interpolating Splines tests                               */
/* ================================================================== */

static void test_yuksel_bezier_mode(void)
{
	printf("test_yuksel_bezier_mode\n");

	/* Triangle-ish open curve: 5 points */
	qaws_vec2 points[] = {
		{0, 0}, {1, 2}, {2, 0}, {3, 2}, {4, 0}
	};

	qaws_yuksel_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = points;
	desc.control_point_count = 5;
	desc.mode = QAWS_YUKSEL_MODE_BEZIER;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_yuksel(&desc, &curve);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(curve != NULL, "yuksel bezier curve allocated");
	TEST_ASSERT(qaws_curve_get_kind(curve) == QAWS_CURVE_KIND_YUKSEL, "kind == YUKSEL");
	TEST_ASSERT(qaws_curve_get_continuity(curve) == QAWS_CONTINUITY_C2, "continuity == C2");
	TEST_ASSERT(qaws_curve_get_span_count(curve) == 4, "span_count == 4");
	TEST_ASSERT(qaws_curve_is_closed(curve) == 0, "not closed");

	/* Evaluate at parameter 0 (first span start) */
	qaws_eval_result_2d r;
	s = qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
		QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &r);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(r.valid_flags & QAWS_EVAL_FLAG_POSITION, "position valid");
	TEST_ASSERT(r.valid_flags & QAWS_EVAL_FLAG_D1, "d1 valid");

	/* Evaluate at each integer parameter to check interpolation */
	{
		unsigned int pi;
		qaws_range range = qaws_curve_get_parameter_range(curve);
		/* Sample throughout to check no NaN or crash */
		for (pi = 0; pi < 40; pi++) {
			qaws_scalar t = range.min_value +
				(qaws_scalar)pi * (range.max_value - range.min_value) / (qaws_scalar)39;
			s = qaws_curve_evaluate_2d(curve, t, QAWS_EVAL_FLAG_POSITION, &r);
			TEST_ASSERT(s == QAWS_STATUS_OK, "eval ok across range");
		}
	}

	qaws_curve_destroy(curve);
}

static void test_yuksel_closed(void)
{
	printf("test_yuksel_closed\n");

	/* Pentagon */
	qaws_vec2 points[5];
	{
		unsigned int pi;
		for (pi = 0; pi < 5; pi++) {
			qaws_scalar angle = (qaws_scalar)(2.0 * 3.14159265358979323846) * (qaws_scalar)pi / (qaws_scalar)5;
			points[pi].x = (qaws_scalar)cos((double)angle);
			points[pi].y = (qaws_scalar)sin((double)angle);
		}
	}

	qaws_yuksel_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = points;
	desc.control_point_count = 5;
	desc.mode = QAWS_YUKSEL_MODE_BEZIER;
	desc.closed = 1;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_yuksel(&desc, &curve);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(qaws_curve_is_closed(curve) == 1, "closed");
	TEST_ASSERT(qaws_curve_get_span_count(curve) == 5, "5 spans for closed pentagon");

	/* Closed curve: first and last points should match */
	qaws_eval_result_2d r0, r1;
	qaws_range range = qaws_curve_get_parameter_range(curve);
	s = qaws_curve_evaluate_2d(curve, range.min_value, QAWS_EVAL_FLAG_POSITION, &r0);
	TEST_ASSERT_STATUS(s);
	s = qaws_curve_evaluate_2d(curve, range.max_value, QAWS_EVAL_FLAG_POSITION, &r1);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq(r0.position.x, r1.position.x), "closed: start.x == end.x");
	TEST_ASSERT(approx_eq(r0.position.y, r1.position.y), "closed: start.y == end.y");

	qaws_curve_destroy(curve);
}

static void test_yuksel_circular_mode(void)
{
	printf("test_yuksel_circular_mode\n");

	qaws_vec2 points[] = {
		{0, 0}, {1, 1}, {2, 0}, {3, 1}, {4, 0}
	};

	qaws_yuksel_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = points;
	desc.control_point_count = 5;
	desc.mode = QAWS_YUKSEL_MODE_CIRCULAR;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_yuksel(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	/* Sample throughout */
	qaws_eval_result_2d r;
	qaws_range range = qaws_curve_get_parameter_range(curve);
	unsigned int pi;
	for (pi = 0; pi < 20; pi++) {
		qaws_scalar t = range.min_value +
			(qaws_scalar)pi * (range.max_value - range.min_value) / (qaws_scalar)19;
		s = qaws_curve_evaluate_2d(curve, t, QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(s == QAWS_STATUS_OK, "circular eval ok");
	}

	qaws_curve_destroy(curve);
}

static void test_yuksel_elliptical_mode(void)
{
	printf("test_yuksel_elliptical_mode\n");

	qaws_vec2 points[] = {
		{0, 0}, {1, 2}, {2, 0}, {3, 2}, {4, 0}
	};

	qaws_yuksel_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = points;
	desc.control_point_count = 5;
	desc.mode = QAWS_YUKSEL_MODE_ELLIPTICAL;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_yuksel(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_eval_result_2d r;
	qaws_range range = qaws_curve_get_parameter_range(curve);
	unsigned int pi;
	for (pi = 0; pi < 20; pi++) {
		qaws_scalar t = range.min_value +
			(qaws_scalar)pi * (range.max_value - range.min_value) / (qaws_scalar)19;
		s = qaws_curve_evaluate_2d(curve, t, QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(s == QAWS_STATUS_OK, "elliptical eval ok");
	}

	qaws_curve_destroy(curve);
}

static void test_yuksel_hybrid_mode(void)
{
	printf("test_yuksel_hybrid_mode\n");

	qaws_vec2 points[] = {
		{0, 0}, {1, 2}, {2, 0}, {3, 2}, {4, 0}
	};

	qaws_yuksel_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = points;
	desc.control_point_count = 5;
	desc.mode = QAWS_YUKSEL_MODE_HYBRID;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_yuksel(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_eval_result_2d r;
	qaws_range range = qaws_curve_get_parameter_range(curve);
	unsigned int pi;
	for (pi = 0; pi < 20; pi++) {
		qaws_scalar t = range.min_value +
			(qaws_scalar)pi * (range.max_value - range.min_value) / (qaws_scalar)19;
		s = qaws_curve_evaluate_2d(curve, t, QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(s == QAWS_STATUS_OK, "hybrid eval ok");
	}

	qaws_curve_destroy(curve);
}

static void test_yuksel_3d(void)
{
	printf("test_yuksel_3d\n");

	qaws_vec3 points[] = {
		{0, 0, 0}, {1, 1, 1}, {2, 0, 2}, {3, 1, 1}, {4, 0, 0}
	};

	qaws_yuksel_desc desc;
	desc.dimension = QAWS_DIMENSION_3D;
	desc.control_points = points;
	desc.control_point_count = 5;
	desc.mode = QAWS_YUKSEL_MODE_BEZIER;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_yuksel(&desc, &curve);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(qaws_curve_get_dimension(curve) == QAWS_DIMENSION_3D, "dim == 3D");

	qaws_eval_result_3d r;
	qaws_range range = qaws_curve_get_parameter_range(curve);
	unsigned int pi;
	for (pi = 0; pi < 20; pi++) {
		qaws_scalar t = range.min_value +
			(qaws_scalar)pi * (range.max_value - range.min_value) / (qaws_scalar)19;
		s = qaws_curve_evaluate_3d(curve, t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &r);
		TEST_ASSERT(s == QAWS_STATUS_OK, "3d eval ok");
	}

	qaws_curve_destroy(curve);
}

static void test_yuksel_error_handling(void)
{
	printf("test_yuksel_error_handling\n");

	qaws_vec2 points[] = { {0, 0}, {1, 1} };
	qaws_yuksel_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = points;
	desc.control_point_count = 2; /* too few */
	desc.mode = QAWS_YUKSEL_MODE_BEZIER;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_yuksel(&desc, &curve);
	TEST_ASSERT(s == QAWS_STATUS_INVALID_CONTROL_POINT_COUNT, "too few points");
	TEST_ASSERT(curve == NULL, "curve is NULL on error");

	/* 3D with circular mode should fail */
	qaws_vec3 pts3[] = { {0,0,0}, {1,1,1}, {2,0,0} };
	desc.dimension = QAWS_DIMENSION_3D;
	desc.control_points = pts3;
	desc.control_point_count = 3;
	desc.mode = QAWS_YUKSEL_MODE_CIRCULAR;
	s = qaws_curve_create_yuksel(&desc, &curve);
	TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT, "circular 3d rejected");

	/* NULL desc */
	s = qaws_curve_create_yuksel(NULL, &curve);
	TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT, "null desc rejected");
}

/* SVG: Yuksel Bezier mode - wave pattern */
static void test_svg_yuksel_bezier(void)
{
	printf("test_svg_yuksel_bezier\n");

	qaws_vec2 points[] = {
		{0, 0}, {1, 2}, {2, -1}, {3, 2}, {4, -1}, {5, 2}, {6, 0}
	};

	qaws_yuksel_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = points;
	desc.control_point_count = 7;
	desc.mode = QAWS_YUKSEL_MODE_BEZIER;
	desc.closed = 0;

	qaws_curve* curve = NULL;
	qaws_curve_create_yuksel(&desc, &curve);

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/yuksel_bezier.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-2.0, (qaws_scalar)7.0, (qaws_scalar)5.0,
		(qaws_scalar)600, (qaws_scalar)400))
	{ qaws_curve_destroy(curve); return; }

	qaws_vec2 samples[SVG_SAMPLES];
	unsigned int count = svg_sample_curve(curve, samples, SVG_SAMPLES);
	svg_polyline(&svg, samples, count, "#e94560", (qaws_scalar)2);
	svg_dots(&svg, points, 7, "#50fa7b", (qaws_scalar)5);
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.8, "Yuksel Bezier", "#e94560");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/yuksel_bezier.svg\n");
	qaws_curve_destroy(curve);
}

/* SVG: Yuksel closed pentagon - all 4 modes comparison */
static void test_svg_yuksel_modes(void)
{
	printf("test_svg_yuksel_modes\n");

	qaws_vec2 points[6];
	{
		unsigned int pi;
		for (pi = 0; pi < 6; pi++) {
			qaws_scalar angle = (qaws_scalar)(2.0 * 3.14159265358979323846) * (qaws_scalar)pi / (qaws_scalar)6;
			points[pi].x = (qaws_scalar)(2.0 * cos((double)angle));
			points[pi].y = (qaws_scalar)(2.0 * sin((double)angle));
		}
	}

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/yuksel_modes.svg",
		(qaws_scalar)-3.0, (qaws_scalar)-3.0, (qaws_scalar)6.0, (qaws_scalar)6.0,
		(qaws_scalar)500, (qaws_scalar)500))
		return;

	char const* colors[] = { "#e94560", "#f5a623", "#50fa7b", "#78dce8" };
	char const* names[] = { "Bezier", "Circular", "Elliptical", "Hybrid" };
	int modes[] = {
		QAWS_YUKSEL_MODE_BEZIER,
		QAWS_YUKSEL_MODE_CIRCULAR,
		QAWS_YUKSEL_MODE_ELLIPTICAL,
		QAWS_YUKSEL_MODE_HYBRID
	};

	unsigned int m;
	for (m = 0; m < 4; m++) {
		qaws_yuksel_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = points;
		desc.control_point_count = 6;
		desc.mode = (qaws_yuksel_mode)modes[m];
		desc.closed = 1;

		qaws_curve* curve = NULL;
		if (qaws_curve_create_yuksel(&desc, &curve) == QAWS_STATUS_OK) {
			qaws_vec2 samples[SVG_SAMPLES];
			unsigned int count = svg_sample_curve(curve, samples, SVG_SAMPLES);
			svg_polyline(&svg, samples, count, colors[m], (qaws_scalar)1.5);
			qaws_curve_destroy(curve);
		}

		svg_label(&svg, (qaws_scalar)-2.8, (qaws_scalar)2.8 - (qaws_scalar)m * (qaws_scalar)0.35,
			names[m], colors[m]);
	}

	svg_dots(&svg, points, 6, "#ffffff", (qaws_scalar)4);
	svg_label(&svg, (qaws_scalar)-2.8, (qaws_scalar)3.2, "Yuksel modes (hexagon)", "#8888aa");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/yuksel_modes.svg\n");
}

/* SVG: Yuksel vs Catmull-Rom comparison */
static void test_svg_yuksel_vs_catmull(void)
{
	printf("test_svg_yuksel_vs_catmull\n");

	qaws_vec2 points[] = {
		{0, 0}, {1, 2}, {2, -0.5}, {3, 1.5}, {4, 0}
	};

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/yuksel_vs_catmull.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-1.5, (qaws_scalar)5.5, (qaws_scalar)4.5,
		(qaws_scalar)600, (qaws_scalar)400))
		return;

	/* Yuksel Bezier */
	{
		qaws_yuksel_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = points;
		desc.control_point_count = 5;
		desc.mode = QAWS_YUKSEL_MODE_BEZIER;
		desc.closed = 0;

		qaws_curve* curve = NULL;
		if (qaws_curve_create_yuksel(&desc, &curve) == QAWS_STATUS_OK) {
			qaws_vec2 samples[SVG_SAMPLES];
			unsigned int count = svg_sample_curve(curve, samples, SVG_SAMPLES);
			svg_polyline(&svg, samples, count, "#e94560", (qaws_scalar)2);
			qaws_curve_destroy(curve);
		}
	}

	/* Catmull-Rom centripetal (needs 2 extra ghost points) */
	{
		qaws_vec2 cr_pts[] = {
			{-1, -2}, {0, 0}, {1, 2}, {2, -0.5}, {3, 1.5}, {4, 0}, {5, -2}
		};
		qaws_catmull_rom_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = cr_pts;
		desc.control_point_count = 7;
		desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		desc.closed = 0;

		qaws_curve* curve = NULL;
		if (qaws_curve_create_catmull_rom(&desc, &curve) == QAWS_STATUS_OK) {
			qaws_vec2 samples[SVG_SAMPLES];
			unsigned int count = svg_sample_curve(curve, samples, SVG_SAMPLES);
			svg_polyline(&svg, samples, count, "#50fa7b", (qaws_scalar)2);
			qaws_curve_destroy(curve);
		}
	}

	svg_dots(&svg, points, 5, "#ffffff", (qaws_scalar)5);
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.8, "Yuksel (red) vs Catmull-Rom (green)", "#8888aa");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/yuksel_vs_catmull.svg\n");
}

int main(void)
{
	printf("=== Qaws Test Suite ===\n\n");

	test_status();
	test_bezier_linear();
	test_bezier_cubic();
	test_bezier_3d();
	test_bezier_derivatives();
	test_evaluate_span();
	test_hermite();
	test_hermite_3d();
	test_catmull_rom();
	test_catmull_rom_variants();
	test_bspline();
	test_bspline_custom_knots();
	test_nurbs();
	test_trajectory();
	test_arc_length();
	test_sampling();
	test_sampling_extended();
	test_traversal();
	test_traversal_mappings();
	test_traversal_arc_length_mode();
	test_traversal_3d();
	test_inspect();
	test_inspection_flags();
	test_bounds_3d();
	test_closest_point();
	test_span_continuity();
	test_family_inspection();
	test_error_handling();

	/* Yuksel C2 interpolating spline tests */
	test_yuksel_bezier_mode();
	test_yuksel_closed();
	test_yuksel_circular_mode();
	test_yuksel_elliptical_mode();
	test_yuksel_hybrid_mode();
	test_yuksel_3d();
	test_yuksel_error_handling();

	/* SVG visual tests */
	printf("\n--- SVG Visual Tests ---\n");
	svg_ensure_output_dir();
	test_svg_circle();
	test_svg_lemniscate();
	test_svg_bezier_s();
	test_svg_catmull_rom_wave();
	test_svg_bspline();
	test_svg_nurbs_ellipse();
	test_svg_hermite_figure8();
	test_svg_all_families();
	test_svg_yuksel_bezier();
	test_svg_yuksel_modes();
	test_svg_yuksel_vs_catmull();

	printf("\n=== Results: %d passed, %d failed ===\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
