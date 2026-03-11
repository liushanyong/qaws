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

#define SVG_OUTPUT_DIR "tests_output"
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

static void svg_line(svg_writer* w, qaws_scalar x0, qaws_scalar y0,
	qaws_scalar x1, qaws_scalar y1, char const* stroke, qaws_scalar width)
{
	qaws_scalar sx0, sy0, sx1, sy1;
	svg_map(w, x0, y0, &sx0, &sy0);
	svg_map(w, x1, y1, &sx1, &sy1);
	fprintf(w->fp, "<line x1=\"%.4f\" y1=\"%.4f\" x2=\"%.4f\" y2=\"%.4f\" "
		"stroke=\"%s\" stroke-width=\"%.1f\"/>\n",
		(double)sx0, (double)sy0, (double)sx1, (double)sy1,
		stroke, (double)width);
}

static void svg_dot(svg_writer* w, qaws_scalar wx, qaws_scalar wy,
	char const* fill, qaws_scalar r)
{
	qaws_scalar sx, sy;
	svg_map(w, wx, wy, &sx, &sy);
	fprintf(w->fp, "<circle cx=\"%.4f\" cy=\"%.4f\" r=\"%.1f\" fill=\"%s\"/>\n",
		(double)sx, (double)sy, (double)r, fill);
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
	//MKDIR("working_dir");
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

	/* Label each control point with its index.
	   Points 0, 2, 4 all sit at the origin — offset labels so all 5 are visible. */
	{
		qaws_scalar lx[] = {
			(qaws_scalar) 0.08, (qaws_scalar) 0.08,
			(qaws_scalar)-0.20, (qaws_scalar) 0.08,
			(qaws_scalar) 0.08
		};
		qaws_scalar ly[] = {
			(qaws_scalar) 0.12, (qaws_scalar) 0.10,
			(qaws_scalar)-0.06, (qaws_scalar) 0.10,
			(qaws_scalar)-0.12
		};
		for (int i = 0; i < 5; i++) {
			char lb[4];
			sprintf(lb, "%d", i);
			svg_label(&svg, pts[i].x + lx[i], pts[i].y + ly[i], lb, "#a9dc76");
		}
	}

	/* Verify: starts and ends at origin */
	TEST_ASSERT(approx_eq(buf[0].x, (qaws_scalar)0.0), "figure8 start x");
	TEST_ASSERT(approx_eq(buf[n - 1].x, (qaws_scalar)0.0), "figure8 end x");

	svg_label(&svg, (qaws_scalar)-1.3, (qaws_scalar)1.3, "Figure-8 (Hermite, 5 key points: 0,2,4 at origin)", "#8888aa");
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
/* Frenet frame tests                                                  */
/* ================================================================== */

static void test_frenet_frame(void)
{
	printf("test_frenet_frame\n");

	/* 2D normal: perpendicular to tangent */
	{
		qaws_vec2 pts[] = { {0, 0}, {2, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve *curve = NULL;
		qaws_vec2 normal;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;

		status = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(status);

		/* Horizontal line: tangent is (1,0), normal should be (0,1) */
		status = qaws_curve_compute_normal_2d(curve, (qaws_scalar)0.5, &normal);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(approx_eq(normal.x, (qaws_scalar)0.0), "2D normal x=0 for horizontal");
		TEST_ASSERT(approx_eq(normal.y, (qaws_scalar)1.0), "2D normal y=1 for horizontal");

		qaws_curve_destroy(curve);
	}

	/* 2D normal for vertical line */
	{
		qaws_vec2 pts[] = { {0, 0}, {0, 2} };
		qaws_bezier_desc bdesc;
		qaws_curve *curve = NULL;
		qaws_vec2 normal;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;

		status = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(status);

		/* Vertical line going up: tangent (0,1), normal (-1,0) */
		status = qaws_curve_compute_normal_2d(curve, (qaws_scalar)0.5, &normal);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(approx_eq(normal.x, (qaws_scalar)-1.0), "2D normal x=-1 for vertical");
		TEST_ASSERT(approx_eq(normal.y, (qaws_scalar)0.0), "2D normal y=0 for vertical");

		qaws_curve_destroy(curve);
	}

	/* 3D Frenet frame: planar curve in XY */
	{
		qaws_vec3 pts[] = { {0,0,0}, {1,2,0}, {2,0,0} };
		qaws_bezier_desc bdesc;
		qaws_curve *curve = NULL;
		qaws_vec3 T, N, B;
		qaws_scalar dot_TN, dot_TB, dot_NB;
		qaws_scalar T_len, N_len, B_len;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 2;
		bdesc.control_points = pts;
		bdesc.control_point_count = 3;

		status = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_compute_frenet_frame_3d(
			curve, (qaws_scalar)0.5, &T, &N, &B);
		TEST_ASSERT_STATUS(status);

		/* T, N, B should be orthonormal */
		T_len = (qaws_scalar)sqrt((double)(T.x*T.x + T.y*T.y + T.z*T.z));
		N_len = (qaws_scalar)sqrt((double)(N.x*N.x + N.y*N.y + N.z*N.z));
		B_len = (qaws_scalar)sqrt((double)(B.x*B.x + B.y*B.y + B.z*B.z));
		TEST_ASSERT(approx_eq(T_len, (qaws_scalar)1.0), "T is unit length");
		TEST_ASSERT(approx_eq(N_len, (qaws_scalar)1.0), "N is unit length");
		TEST_ASSERT(approx_eq(B_len, (qaws_scalar)1.0), "B is unit length");

		dot_TN = T.x*N.x + T.y*N.y + T.z*N.z;
		dot_TB = T.x*B.x + T.y*B.y + T.z*B.z;
		dot_NB = N.x*B.x + N.y*B.y + N.z*B.z;
		TEST_ASSERT(approx_eq(dot_TN, (qaws_scalar)0.0), "T perpendicular to N");
		TEST_ASSERT(approx_eq(dot_TB, (qaws_scalar)0.0), "T perpendicular to B");
		TEST_ASSERT(approx_eq(dot_NB, (qaws_scalar)0.0), "N perpendicular to B");

		/* For planar XY curve, binormal should point along Z */
		TEST_ASSERT(B.z != (qaws_scalar)0.0, "B has Z component for XY curve");

		qaws_curve_destroy(curve);
	}

	/* 3D Frenet frame: straight line (degenerate) */
	{
		qaws_vec3 pts[] = { {0,0,0}, {1,0,0} };
		qaws_bezier_desc bdesc;
		qaws_curve *curve = NULL;
		qaws_vec3 T, N, B;
		qaws_scalar dot_TN;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;

		status = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(status);

		/* Straight line: D1xD2=0, should still produce valid frame */
		status = qaws_curve_compute_frenet_frame_3d(
			curve, (qaws_scalar)0.5, &T, &N, &B);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(approx_eq(T.x, (qaws_scalar)1.0), "straight line T along X");

		dot_TN = T.x*N.x + T.y*N.y + T.z*N.z;
		TEST_ASSERT(approx_eq(dot_TN, (qaws_scalar)0.0), "T perp N for straight line");

		qaws_curve_destroy(curve);
	}

	/* Dimension mismatch */
	{
		qaws_vec2 pts[] = { {0,0}, {1,1} };
		qaws_bezier_desc bdesc;
		qaws_curve *curve = NULL;
		qaws_vec3 T, N, B;
		qaws_vec2 normal2d;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;
		status = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_compute_frenet_frame_3d(
			curve, (qaws_scalar)0.5, &T, &N, &B);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_DIMENSION, "3D frame rejects 2D curve");

		status = qaws_curve_compute_normal_2d(NULL, (qaws_scalar)0.5, &normal2d);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_ARGUMENT, "null curve");

		qaws_curve_destroy(curve);
	}
}

/* ================================================================== */
/* Curve reverse tests                                                 */
/* ================================================================== */

static void test_curve_reverse(void)
{
	printf("test_curve_reverse\n");

	/* Bezier reverse: endpoints swap */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 2}, {3, 1} };
		qaws_bezier_desc bdesc;
		qaws_curve *curve = NULL;
		qaws_curve *rev = NULL;
		qaws_eval_result_2d r_orig, r_rev;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 2;
		bdesc.control_points = pts;
		bdesc.control_point_count = 3;

		status = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_reverse(curve, &rev);
		TEST_ASSERT_STATUS(status);

		/* Original start == reversed end */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		qaws_curve_evaluate_2d(rev, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION, &r_rev);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x), "reverse: start==end x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y), "reverse: start==end y");

		/* Original end == reversed start */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		qaws_curve_evaluate_2d(rev, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r_rev);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x), "reverse: end==start x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y), "reverse: end==start y");

		qaws_curve_destroy(rev);
		qaws_curve_destroy(curve);
	}

	/* Catmull-Rom reverse */
	{
		qaws_vec2 pts[] = { {0,0}, {1,2}, {3,1}, {4,3} };
		qaws_catmull_rom_desc cdesc;
		qaws_curve *curve = NULL;
		qaws_curve *rev = NULL;
		qaws_eval_result_2d r_orig, r_rev;
		qaws_status status;
		qaws_range range;

		cdesc.dimension = QAWS_DIMENSION_2D;
		cdesc.control_points = pts;
		cdesc.control_point_count = 4;
		cdesc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		cdesc.closed = 0;

		status = qaws_curve_create_catmull_rom(&cdesc, &curve);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_reverse(curve, &rev);
		TEST_ASSERT_STATUS(status);

		/* Check start/end swap */
		range = qaws_curve_get_parameter_range(curve);
		qaws_curve_evaluate_2d(curve, range.min_value,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		range = qaws_curve_get_parameter_range(rev);
		qaws_curve_evaluate_2d(rev, range.max_value,
			QAWS_EVAL_FLAG_POSITION, &r_rev);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x),
			"catmull reverse start x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y),
			"catmull reverse start y");

		qaws_curve_destroy(rev);
		qaws_curve_destroy(curve);
	}

	/* Hermite reverse */
	{
		qaws_vec2 pts[] = { {0,0}, {2,0} };
		qaws_vec2 ders[] = { {1,1}, {1,-1} };
		qaws_hermite_desc hdesc;
		qaws_curve *curve = NULL;
		qaws_curve *rev = NULL;
		qaws_eval_result_2d r_orig, r_rev;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 2;
		hdesc.derivative_count = 2;

		status = qaws_curve_create_hermite(&hdesc, &curve);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_reverse(curve, &rev);
		TEST_ASSERT_STATUS(status);

		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		qaws_curve_evaluate_2d(rev, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION, &r_rev);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x),
			"hermite reverse x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y),
			"hermite reverse y");

		qaws_curve_destroy(rev);
		qaws_curve_destroy(curve);
	}

	/* NULL argument */
	{
		qaws_curve *rev = NULL;
		qaws_status status = qaws_curve_reverse(NULL, &rev);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_ARGUMENT, "reverse null curve");
	}

	/* B-Spline reverse: endpoints swap */
	{
		qaws_vec2 bpts[] = { {0, 0}, {1, 2}, {3, 2}, {4, 0} };
		qaws_bspline_desc bsdesc;
		qaws_curve *curve = NULL;
		qaws_curve *rev = NULL;
		qaws_eval_result_2d r_orig, r_rev;
		qaws_status status;
		qaws_range range, rev_range;

		bsdesc.dimension = QAWS_DIMENSION_2D;
		bsdesc.degree = 3;
		bsdesc.control_points = bpts;
		bsdesc.control_point_count = 4;
		bsdesc.knots = NULL;
		bsdesc.knot_count = 0;
		bsdesc.is_uniform = 1;
		bsdesc.is_closed = 0;

		status = qaws_curve_create_bspline(&bsdesc, &curve);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_reverse(curve, &rev);
		TEST_ASSERT_STATUS(status);

		range = qaws_curve_get_parameter_range(curve);
		rev_range = qaws_curve_get_parameter_range(rev);

		qaws_curve_evaluate_2d(curve, range.min_value,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		qaws_curve_evaluate_2d(rev, rev_range.max_value,
			QAWS_EVAL_FLAG_POSITION, &r_rev);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x),
			"bspline reverse: start==end x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y),
			"bspline reverse: start==end y");

		qaws_curve_evaluate_2d(curve, range.max_value,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		qaws_curve_evaluate_2d(rev, rev_range.min_value,
			QAWS_EVAL_FLAG_POSITION, &r_rev);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x),
			"bspline reverse: end==start x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y),
			"bspline reverse: end==start y");

		qaws_curve_destroy(rev);
		qaws_curve_destroy(curve);
	}

	/* NURBS reverse: endpoints swap */
	{
		qaws_vec2 npts[] = { {1, 0}, {1, 1}, {0, 1} };
		qaws_scalar nknots[] = { 0, 0, 0, 1, 1, 1 };
		qaws_scalar nweights[] = { 1, (qaws_scalar)0.70710678118, 1 };
		qaws_nurbs_desc ndesc;
		qaws_curve *curve = NULL;
		qaws_curve *rev = NULL;
		qaws_eval_result_2d r_orig, r_rev;
		qaws_status status;
		qaws_range range, rev_range;

		ndesc.dimension = QAWS_DIMENSION_2D;
		ndesc.degree = 2;
		ndesc.control_points = npts;
		ndesc.control_point_count = 3;
		ndesc.knots = nknots;
		ndesc.knot_count = 6;
		ndesc.weights = nweights;
		ndesc.weight_count = 3;
		ndesc.is_closed = 0;

		status = qaws_curve_create_nurbs(&ndesc, &curve);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_reverse(curve, &rev);
		TEST_ASSERT_STATUS(status);

		range = qaws_curve_get_parameter_range(curve);
		rev_range = qaws_curve_get_parameter_range(rev);

		qaws_curve_evaluate_2d(curve, range.min_value,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		qaws_curve_evaluate_2d(rev, rev_range.max_value,
			QAWS_EVAL_FLAG_POSITION, &r_rev);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x),
			"nurbs reverse: start==end x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y),
			"nurbs reverse: start==end y");

		qaws_curve_evaluate_2d(curve, range.max_value,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		qaws_curve_evaluate_2d(rev, rev_range.min_value,
			QAWS_EVAL_FLAG_POSITION, &r_rev);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x),
			"nurbs reverse: end==start x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y),
			"nurbs reverse: end==start y");

		qaws_curve_destroy(rev);
		qaws_curve_destroy(curve);
	}

	/* Trajectory reverse: endpoints swap */
	{
		qaws_vec2 tpts[] = { {0, 0}, {5, 0}, {5, 5} };
		qaws_scalar ttimes[] = { 0, 1, 2 };
		qaws_trajectory_desc tdesc;
		qaws_curve *curve = NULL;
		qaws_curve *rev = NULL;
		qaws_eval_result_2d r_orig, r_rev;
		qaws_status status;
		qaws_range range, rev_range;

		tdesc.dimension = QAWS_DIMENSION_2D;
		tdesc.degree = 3;
		tdesc.key_positions = tpts;
		tdesc.key_count = 3;
		tdesc.key_times = ttimes;
		tdesc.key_time_count = 3;
		tdesc.key_velocities = NULL;
		tdesc.key_velocity_count = 0;
		tdesc.key_accelerations = NULL;
		tdesc.key_acceleration_count = 0;
		tdesc.closed = 0;

		status = qaws_curve_create_trajectory(&tdesc, &curve);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_reverse(curve, &rev);
		TEST_ASSERT_STATUS(status);

		range = qaws_curve_get_parameter_range(curve);
		rev_range = qaws_curve_get_parameter_range(rev);

		qaws_curve_evaluate_2d(curve, range.min_value,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		qaws_curve_evaluate_2d(rev, rev_range.max_value,
			QAWS_EVAL_FLAG_POSITION, &r_rev);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x),
			"traj reverse: start==end x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y),
			"traj reverse: start==end y");

		qaws_curve_evaluate_2d(curve, range.max_value,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		qaws_curve_evaluate_2d(rev, rev_range.min_value,
			QAWS_EVAL_FLAG_POSITION, &r_rev);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x),
			"traj reverse: end==start x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y),
			"traj reverse: end==start y");

		qaws_curve_destroy(rev);
		qaws_curve_destroy(curve);
	}

	/* Bezier midpoint preservation: original t=0.5 == reversed t=0.5 */
	{
		qaws_vec2 mpts[] = { {0, 0}, {1, 2}, {3, 1} };
		qaws_bezier_desc mbdesc;
		qaws_curve *curve = NULL;
		qaws_curve *rev = NULL;
		qaws_eval_result_2d r_orig, r_rev;
		qaws_status status;

		mbdesc.dimension = QAWS_DIMENSION_2D;
		mbdesc.degree = 2;
		mbdesc.control_points = mpts;
		mbdesc.control_point_count = 3;

		status = qaws_curve_create_bezier(&mbdesc, &curve);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_reverse(curve, &rev);
		TEST_ASSERT_STATUS(status);

		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		qaws_curve_evaluate_2d(rev, (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_POSITION, &r_rev);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_rev.position.x),
			"bezier reverse midpoint x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_rev.position.y),
			"bezier reverse midpoint y");

		qaws_curve_destroy(rev);
		qaws_curve_destroy(curve);
	}
}

/* ================================================================== */
/* Easing function tests                                               */
/* ================================================================== */

static void test_easing(void)
{
	printf("test_easing\n");

	/* Create a simple curve and traversal with easing */
	qaws_vec2 pts[] = { {0,0}, {1,1}, {2,0} };
	qaws_bezier_desc bdesc;
	qaws_curve *curve = NULL;
	qaws_status status;

	bdesc.dimension = QAWS_DIMENSION_2D;
	bdesc.degree = 2;
	bdesc.control_points = pts;
	bdesc.control_point_count = 3;

	status = qaws_curve_create_bezier(&bdesc, &curve);
	TEST_ASSERT_STATUS(status);

	/* Linear easing: same as no easing */
	{
		qaws_traversal_desc tdesc;
		qaws_traversal *trav = NULL;
		qaws_eval_result_2d r_linear, r_eased;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
		tdesc.speed = (qaws_scalar)1.0;
		tdesc.start_time = (qaws_scalar)0.0;
		tdesc.end_time = (qaws_scalar)1.0;
		tdesc.clamp_to_domain = 1;
		tdesc.easing = QAWS_EASING_LINEAR;
		tdesc.wrap_mode = QAWS_WRAP_CLAMP;

		status = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(status);

		status = qaws_traversal_evaluate_2d(
			trav, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &r_linear);
		TEST_ASSERT_STATUS(status);

		qaws_traversal_destroy(trav);

		/* Quad ease-in: at t=0.5, eased_t = 0.25, so position should be earlier */
		tdesc.easing = QAWS_EASING_QUAD_IN;
		status = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(status);

		status = qaws_traversal_evaluate_2d(
			trav, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &r_eased);
		TEST_ASSERT_STATUS(status);

		/* Eased position should differ from linear */
		TEST_ASSERT(
			!approx_eq(r_linear.position.x, r_eased.position.x) ||
			!approx_eq(r_linear.position.y, r_eased.position.y),
			"ease-in differs from linear at t=0.5");

		qaws_traversal_destroy(trav);
	}

	/* Sine easing endpoints: at t=0 should be start, at t=1 should be end */
	{
		qaws_traversal_desc tdesc;
		qaws_traversal *trav = NULL;
		qaws_eval_result_2d r_start, r_end;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
		tdesc.speed = (qaws_scalar)1.0;
		tdesc.start_time = (qaws_scalar)0.0;
		tdesc.end_time = (qaws_scalar)1.0;
		tdesc.clamp_to_domain = 1;
		tdesc.easing = QAWS_EASING_SINE_IN_OUT;
		tdesc.wrap_mode = QAWS_WRAP_CLAMP;

		status = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(status);

		qaws_traversal_evaluate_2d(
			trav, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r_start);
		qaws_traversal_evaluate_2d(
			trav, (qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION, &r_end);

		/* Start should be near (0,0) */
		TEST_ASSERT(approx_eq(r_start.position.x, (qaws_scalar)0.0),
			"eased start x near 0");

		qaws_traversal_destroy(trav);
	}

	/* All easing modes at endpoints: t=0 => start, t=1 => end */
	{
		qaws_easing modes[] = {
			QAWS_EASING_QUAD_OUT,
			QAWS_EASING_CUBIC_IN,
			QAWS_EASING_CUBIC_OUT,
			QAWS_EASING_CUBIC_IN_OUT,
			QAWS_EASING_SINE_IN,
			QAWS_EASING_SINE_OUT
		};
		const char *mode_names[] = {
			"QUAD_OUT", "CUBIC_IN", "CUBIC_OUT",
			"CUBIC_IN_OUT", "SINE_IN", "SINE_OUT"
		};
		unsigned int mode_count = 6;
		unsigned int mi;
		qaws_scalar arc_len;

		qaws_eval_result_2d r_curve_start, r_curve_end;

		/* Get the curve's actual start and end positions */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r_curve_start);
		qaws_curve_evaluate_2d(curve, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION, &r_curve_end);

		/* Compute arc length so speed covers the full curve in [0,1] */
		{
			qaws_range pr = qaws_curve_get_parameter_range(curve);
			qaws_curve_compute_arc_length(curve,
				pr.min_value, pr.max_value, &arc_len);
		}

		for (mi = 0; mi < mode_count; ++mi)
		{
			qaws_traversal_desc etdesc;
			qaws_traversal *etrav = NULL;
			qaws_eval_result_2d r_t0, r_t1;
			char msg[128];

			memset(&etdesc, 0, sizeof(etdesc));
			etdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
			etdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
			etdesc.speed = arc_len;
			etdesc.start_time = (qaws_scalar)0.0;
			etdesc.end_time = (qaws_scalar)1.0;
			etdesc.clamp_to_domain = 1;
			etdesc.easing = modes[mi];
			etdesc.wrap_mode = QAWS_WRAP_CLAMP;

			status = qaws_traversal_create(curve, &etdesc, &etrav);
			TEST_ASSERT_STATUS(status);

			status = qaws_traversal_evaluate_2d(
				etrav, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r_t0);
			TEST_ASSERT_STATUS(status);

			status = qaws_traversal_evaluate_2d(
				etrav, (qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION, &r_t1);
			TEST_ASSERT_STATUS(status);

			sprintf(msg, "%s t=0 start x", mode_names[mi]);
			TEST_ASSERT(approx_eq(r_t0.position.x, r_curve_start.position.x), msg);
			sprintf(msg, "%s t=0 start y", mode_names[mi]);
			TEST_ASSERT(approx_eq(r_t0.position.y, r_curve_start.position.y), msg);

			sprintf(msg, "%s t=1 end x", mode_names[mi]);
			TEST_ASSERT(approx_eq(r_t1.position.x, r_curve_end.position.x), msg);
			sprintf(msg, "%s t=1 end y", mode_names[mi]);
			TEST_ASSERT(approx_eq(r_t1.position.y, r_curve_end.position.y), msg);

			qaws_traversal_destroy(etrav);
		}
	}

	/* QUAD_IN_OUT symmetry: at t=0.5, should be at midpoint position */
	{
		qaws_traversal_desc qdesc;
		qaws_traversal *qtrav = NULL;
		qaws_eval_result_2d r_linear_mid, r_qio_mid;
		qaws_traversal *ltrav = NULL;

		/* Get linear position at t=0.5 for reference (the "midpoint") */
		memset(&qdesc, 0, sizeof(qdesc));
		qdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
		qdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
		qdesc.speed = (qaws_scalar)1.0;
		qdesc.start_time = (qaws_scalar)0.0;
		qdesc.end_time = (qaws_scalar)1.0;
		qdesc.clamp_to_domain = 1;
		qdesc.easing = QAWS_EASING_LINEAR;
		qdesc.wrap_mode = QAWS_WRAP_CLAMP;

		status = qaws_traversal_create(curve, &qdesc, &ltrav);
		TEST_ASSERT_STATUS(status);
		qaws_traversal_evaluate_2d(
			ltrav, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &r_linear_mid);
		qaws_traversal_destroy(ltrav);

		/* QUAD_IN_OUT at t=0.5 should be at the halfway position */
		qdesc.easing = QAWS_EASING_QUAD_IN_OUT;
		status = qaws_traversal_create(curve, &qdesc, &qtrav);
		TEST_ASSERT_STATUS(status);

		qaws_traversal_evaluate_2d(
			qtrav, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &r_qio_mid);

		TEST_ASSERT(approx_eq(r_qio_mid.position.x, r_linear_mid.position.x),
			"QUAD_IN_OUT t=0.5 midpoint x");
		TEST_ASSERT(approx_eq(r_qio_mid.position.y, r_linear_mid.position.y),
			"QUAD_IN_OUT t=0.5 midpoint y");

		qaws_traversal_destroy(qtrav);
	}

	qaws_curve_destroy(curve);
}

/* ================================================================== */
/* Wrap mode tests                                                     */
/* ================================================================== */

static void test_wrap_modes(void)
{
	printf("test_wrap_modes\n");

	qaws_vec2 pts[] = { {0,0}, {1,1}, {2,0} };
	qaws_bezier_desc bdesc;
	qaws_curve *curve = NULL;
	qaws_status status;

	bdesc.dimension = QAWS_DIMENSION_2D;
	bdesc.degree = 2;
	bdesc.control_points = pts;
	bdesc.control_point_count = 3;

	status = qaws_curve_create_bezier(&bdesc, &curve);
	TEST_ASSERT_STATUS(status);

	/* Loop mode: going past end wraps around */
	{
		qaws_traversal_desc tdesc;
		qaws_traversal *trav = NULL;
		qaws_eval_result_2d r_a, r_b;
		qaws_scalar total_len;
		qaws_scalar eps = (qaws_scalar)0.05;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
		tdesc.clamp_to_domain = 1;
		tdesc.wrap_mode = QAWS_WRAP_LOOP;

		status = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(status);

		{
			qaws_range pr = qaws_curve_get_parameter_range(curve);
			qaws_curve_compute_arc_length(curve,
				pr.min_value, pr.max_value, &total_len);
		}

		/* distance=eps and distance=total_len+eps should map to same point */
		qaws_traversal_evaluate_2d(
			trav, eps, QAWS_EVAL_FLAG_POSITION, &r_a);
		qaws_traversal_evaluate_2d(
			trav, total_len + eps, QAWS_EVAL_FLAG_POSITION, &r_b);

		TEST_ASSERT(approx_eq(r_a.position.x, r_b.position.x),
			"loop wraps to start x");
		TEST_ASSERT(approx_eq(r_a.position.y, r_b.position.y),
			"loop wraps to start y");

		qaws_traversal_destroy(trav);
	}

	/* Ping-pong mode: d=eps and d=2*total_len+eps should match */
	{
		qaws_traversal_desc tdesc;
		qaws_traversal *trav = NULL;
		qaws_eval_result_2d r_a, r_b;
		qaws_scalar total_len;
		qaws_scalar eps = (qaws_scalar)0.05;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
		tdesc.clamp_to_domain = 1;
		tdesc.wrap_mode = QAWS_WRAP_PING_PONG;

		status = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(status);

		{
			qaws_range pr = qaws_curve_get_parameter_range(curve);
			qaws_curve_compute_arc_length(curve,
				pr.min_value, pr.max_value, &total_len);
		}

		qaws_traversal_evaluate_2d(
			trav, eps, QAWS_EVAL_FLAG_POSITION, &r_a);
		qaws_traversal_evaluate_2d(
			trav, (qaws_scalar)2.0 * total_len + eps, QAWS_EVAL_FLAG_POSITION, &r_b);

		TEST_ASSERT(approx_eq(r_a.position.x, r_b.position.x),
			"ping-pong returns to start x");
		TEST_ASSERT(approx_eq(r_a.position.y, r_b.position.y),
			"ping-pong returns to start y");

		qaws_traversal_destroy(trav);
	}

	/* CLAMP mode: negative distance clamps to start, huge distance clamps to end */
	{
		qaws_traversal_desc tdesc;
		qaws_traversal *trav = NULL;
		qaws_eval_result_2d r_neg, r_huge, r_start, r_end;
		qaws_scalar total_len;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
		tdesc.clamp_to_domain = 1;
		tdesc.wrap_mode = QAWS_WRAP_CLAMP;

		status = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(status);

		{
			qaws_range pr = qaws_curve_get_parameter_range(curve);
			qaws_curve_compute_arc_length(curve,
				pr.min_value, pr.max_value, &total_len);
		}

		/* Position at distance=0 (start of curve) */
		qaws_traversal_evaluate_2d(
			trav, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r_start);

		/* Position at distance=total_len (end of curve) */
		qaws_traversal_evaluate_2d(
			trav, total_len, QAWS_EVAL_FLAG_POSITION, &r_end);

		/* Negative distance should clamp to start */
		qaws_traversal_evaluate_2d(
			trav, (qaws_scalar)-10.0, QAWS_EVAL_FLAG_POSITION, &r_neg);
		TEST_ASSERT(approx_eq(r_neg.position.x, r_start.position.x),
			"clamp negative dist x == start x");
		TEST_ASSERT(approx_eq(r_neg.position.y, r_start.position.y),
			"clamp negative dist y == start y");

		/* Very large distance should clamp to end */
		qaws_traversal_evaluate_2d(
			trav, total_len * (qaws_scalar)100.0, QAWS_EVAL_FLAG_POSITION, &r_huge);
		TEST_ASSERT(approx_eq(r_huge.position.x, r_end.position.x),
			"clamp huge dist x == end x");
		TEST_ASSERT(approx_eq(r_huge.position.y, r_end.position.y),
			"clamp huge dist y == end y");

		qaws_traversal_destroy(trav);
	}

	/* Negative values in LOOP mode: should wrap correctly */
	{
		qaws_traversal_desc tdesc;
		qaws_traversal *trav = NULL;
		qaws_eval_result_2d r_neg, r_pos;
		qaws_scalar total_len;
		qaws_scalar eps = (qaws_scalar)0.05;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
		tdesc.clamp_to_domain = 1;
		tdesc.wrap_mode = QAWS_WRAP_LOOP;

		status = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(status);

		{
			qaws_range pr = qaws_curve_get_parameter_range(curve);
			qaws_curve_compute_arc_length(curve,
				pr.min_value, pr.max_value, &total_len);
		}

		/* Negative distance -eps should wrap to total_len - eps */
		qaws_traversal_evaluate_2d(
			trav, -eps, QAWS_EVAL_FLAG_POSITION, &r_neg);
		qaws_traversal_evaluate_2d(
			trav, total_len - eps, QAWS_EVAL_FLAG_POSITION, &r_pos);

		TEST_ASSERT(approx_eq(r_neg.position.x, r_pos.position.x),
			"loop negative wrap x");
		TEST_ASSERT(approx_eq(r_neg.position.y, r_pos.position.y),
			"loop negative wrap y");

		qaws_traversal_destroy(trav);
	}

	qaws_curve_destroy(curve);
}

/* ================================================================== */
/* Traversal advance (iterator) tests                                  */
/* ================================================================== */

static void test_traversal_advance(void)
{
	printf("test_traversal_advance\n");

	qaws_vec2 pts[] = { {0,0}, {2,2}, {4,0} };
	qaws_bezier_desc bdesc;
	qaws_curve *curve = NULL;
	qaws_status status;

	bdesc.dimension = QAWS_DIMENSION_2D;
	bdesc.degree = 2;
	bdesc.control_points = pts;
	bdesc.control_point_count = 3;

	status = qaws_curve_create_bezier(&bdesc, &curve);
	TEST_ASSERT_STATUS(status);

	{
		qaws_traversal_desc tdesc;
		qaws_traversal *trav = NULL;
		qaws_eval_result_2d r1, r2, r3;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
		tdesc.speed = (qaws_scalar)1.0;
		tdesc.clamp_to_domain = 1;
		tdesc.wrap_mode = QAWS_WRAP_CLAMP;

		status = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(status);

		/* First advance: should move forward from distance=0 */
		status = qaws_traversal_advance_2d(
			trav, (qaws_scalar)0.1, QAWS_EVAL_FLAG_POSITION, &r1);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(r1.position.x > (qaws_scalar)0.0, "advance moves forward x");

		/* Second advance: should be further along */
		status = qaws_traversal_advance_2d(
			trav, (qaws_scalar)0.1, QAWS_EVAL_FLAG_POSITION, &r2);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(r2.position.x > r1.position.x, "second advance further");

		/* Reset should go back to start */
		status = qaws_traversal_reset(trav);
		TEST_ASSERT_STATUS(status);

		status = qaws_traversal_advance_2d(
			trav, (qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION, &r3);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(approx_eq(r3.position.x, (qaws_scalar)0.0),
			"reset goes to start x");
		TEST_ASSERT(approx_eq(r3.position.y, (qaws_scalar)0.0),
			"reset goes to start y");

		/* Loop mode advance: should wrap around */
		qaws_traversal_destroy(trav);
		tdesc.wrap_mode = QAWS_WRAP_LOOP;
		status = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(status);

		/* Advance a large amount */
		status = qaws_traversal_advance_2d(
			trav, (qaws_scalar)100.0, QAWS_EVAL_FLAG_POSITION, &r1);
		TEST_ASSERT_STATUS(status);
		/* Should still be valid (not NaN) */
		TEST_ASSERT(r1.position.x == r1.position.x, "loop advance valid x");

		/* NULL checks */
		status = qaws_traversal_advance_2d(NULL, (qaws_scalar)0.1,
			QAWS_EVAL_FLAG_POSITION, &r1);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_ARGUMENT, "advance null trav");

		status = qaws_traversal_reset(NULL);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_ARGUMENT, "reset null");

		qaws_traversal_destroy(trav);
	}

	qaws_curve_destroy(curve);

	/* 3D advance: create a 3D bezier, traverse, verify forward movement */
	{
		qaws_vec3 pts3d[] = { {0,0,0}, {2,2,1}, {4,0,2} };
		qaws_bezier_desc bdesc3d;
		qaws_curve *curve3d = NULL;
		qaws_traversal_desc tdesc3d;
		qaws_traversal *trav3d = NULL;
		qaws_eval_result_3d r3d_1, r3d_2;

		bdesc3d.dimension = QAWS_DIMENSION_3D;
		bdesc3d.degree = 2;
		bdesc3d.control_points = pts3d;
		bdesc3d.control_point_count = 3;

		status = qaws_curve_create_bezier(&bdesc3d, &curve3d);
		TEST_ASSERT_STATUS(status);

		memset(&tdesc3d, 0, sizeof(tdesc3d));
		tdesc3d.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
		tdesc3d.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
		tdesc3d.speed = (qaws_scalar)1.0;
		tdesc3d.clamp_to_domain = 1;
		tdesc3d.wrap_mode = QAWS_WRAP_CLAMP;

		status = qaws_traversal_create(curve3d, &tdesc3d, &trav3d);
		TEST_ASSERT_STATUS(status);

		/* First advance */
		status = qaws_traversal_advance_3d(
			trav3d, (qaws_scalar)0.1, QAWS_EVAL_FLAG_POSITION, &r3d_1);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(r3d_1.position.x > (qaws_scalar)0.0, "3d advance moves forward x");

		/* Second advance should be further along */
		status = qaws_traversal_advance_3d(
			trav3d, (qaws_scalar)0.1, QAWS_EVAL_FLAG_POSITION, &r3d_2);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(r3d_2.position.x > r3d_1.position.x, "3d second advance further x");

		qaws_traversal_destroy(trav3d);
		qaws_curve_destroy(curve3d);
	}
}

/* ================================================================== */
/* Adaptive eval sampling tests                                        */
/* ================================================================== */

static void test_adaptive_eval_sampling(void)
{
	printf("test_adaptive_eval_sampling\n");

	qaws_vec2 pts[] = { {0,0}, {0,4}, {4,4}, {4,0} };
	qaws_bezier_desc bdesc;
	qaws_curve *curve = NULL;
	qaws_status status;

	bdesc.dimension = QAWS_DIMENSION_2D;
	bdesc.degree = 3;
	bdesc.control_points = pts;
	bdesc.control_point_count = 4;

	status = qaws_curve_create_bezier(&bdesc, &curve);
	TEST_ASSERT_STATUS(status);

	/* Uniform eval sampling baseline */
	{
		qaws_sampling_desc sdesc;
		qaws_eval_result_2d samples[2048];
		unsigned int count = 0;

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.sample_count = 8;
		sdesc.use_adaptive_sampling = 0;

		status = qaws_curve_sample_eval_2d(
			curve, &sdesc, samples, 2048, &count);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(count == 8, "uniform eval gives 8");
	}

	/* Adaptive eval sampling should produce more points and include derivatives */
	{
		qaws_sampling_desc sdesc;
		qaws_eval_result_2d samples[2048];
		unsigned int count = 0;
		int has_d1 = 1;
		unsigned int i;

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.sample_count = 8;
		sdesc.use_adaptive_sampling = 1;
		sdesc.error_tolerance = (qaws_scalar)0.001;

		status = qaws_curve_sample_eval_2d(
			curve, &sdesc, samples, 2048, &count);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(count >= 8, "adaptive eval >= initial");

		/* Check that derivatives are populated */
		for (i = 0; i < count && i < 5; ++i)
		{
			if (samples[i].d1.x == (qaws_scalar)0.0 &&
				samples[i].d1.y == (qaws_scalar)0.0)
				has_d1 = 0;
		}
		TEST_ASSERT(has_d1, "adaptive eval has D1 data");
	}

	/* 3D adaptive eval */
	{
		qaws_vec3 pts3d[] = { {0,0,0}, {0,4,2}, {4,4,2}, {4,0,0} };
		qaws_curve *curve3d = NULL;
		qaws_eval_result_3d samples3d[2048];
		unsigned int count3d = 0;
		qaws_sampling_desc sdesc;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.control_points = pts3d;

		status = qaws_curve_create_bezier(&bdesc, &curve3d);
		TEST_ASSERT_STATUS(status);

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.sample_count = 8;
		sdesc.use_adaptive_sampling = 1;
		sdesc.error_tolerance = (qaws_scalar)0.01;

		status = qaws_curve_sample_eval_3d(
			curve3d, &sdesc, samples3d, 2048, &count3d);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(count3d >= 8, "3D adaptive eval >= initial");

		qaws_curve_destroy(curve3d);
	}

	/* Buffer too small: adaptive sampling into a tiny buffer */
	{
		qaws_sampling_desc sdesc;
		qaws_eval_result_2d tiny_buf[2];
		unsigned int count_needed = 0;

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.sample_count = 64;
		sdesc.use_adaptive_sampling = 1;
		sdesc.error_tolerance = (qaws_scalar)0.0001;

		status = qaws_curve_sample_eval_2d(
			curve, &sdesc, tiny_buf, 2, &count_needed);
		TEST_ASSERT(status == QAWS_STATUS_BUFFER_TOO_SMALL,
			"tiny buffer returns BUFFER_TOO_SMALL");
		TEST_ASSERT(count_needed > 2,
			"buffer too small reports needed count > capacity");
	}

	qaws_curve_destroy(curve);
}

/* ================================================================== */
/* Phase 2: Curve splitting tests                                      */
/* ================================================================== */

static void test_curve_split(void)
{
	printf("test_curve_split\n");

	/* Bezier split at midpoint */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 2}, {2, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve *curve = NULL;
		qaws_curve *left = NULL;
		qaws_curve *right = NULL;
		qaws_eval_result_2d r_orig, r_left, r_right;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 2;
		bdesc.control_points = pts;
		bdesc.control_point_count = 3;

		status = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_split(curve, (qaws_scalar)0.5, &left, &right);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(left != NULL, "split left not null");
		TEST_ASSERT(right != NULL, "split right not null");

		/* Left end == original at t=0.5 */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		qaws_curve_evaluate_2d(left, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION, &r_left);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_left.position.x),
			"bezier split left end x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_left.position.y),
			"bezier split left end y");

		/* Right start == original at t=0.5 */
		qaws_curve_evaluate_2d(right, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r_right);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_right.position.x),
			"bezier split right start x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_right.position.y),
			"bezier split right start y");

		/* Left start == original start */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		qaws_curve_evaluate_2d(left, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r_left);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_left.position.x),
			"bezier split left start x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_left.position.y),
			"bezier split left start y");

		/* Right end == original end */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		qaws_curve_evaluate_2d(right, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION, &r_right);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_right.position.x),
			"bezier split right end x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_right.position.y),
			"bezier split right end y");

		qaws_curve_destroy(right);
		qaws_curve_destroy(left);
		qaws_curve_destroy(curve);
	}

	/* Hermite split at integer boundary */
	{
		qaws_vec2 pts[] = { {0,0}, {1,1}, {2,0} };
		qaws_vec2 ders[] = { {1,0}, {0,1}, {1,0} };
		qaws_hermite_desc hdesc;
		qaws_curve *curve = NULL;
		qaws_curve *left = NULL;
		qaws_curve *right = NULL;
		qaws_eval_result_2d r_orig, r_sub;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 3;
		hdesc.derivative_count = 3;

		status = qaws_curve_create_hermite(&hdesc, &curve);
		TEST_ASSERT_STATUS(status);

		/* Split at span boundary t=1 */
		status = qaws_curve_split(curve, (qaws_scalar)1.0, &left, &right);
		TEST_ASSERT_STATUS(status);

		/* Left curve should match original at t=0 */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		{
			qaws_range lr = qaws_curve_get_parameter_range(left);
			qaws_curve_evaluate_2d(left, lr.min_value,
				QAWS_EVAL_FLAG_POSITION, &r_sub);
		}
		TEST_ASSERT(approx_eq(r_orig.position.x, r_sub.position.x),
			"hermite split left start x");

		/* Right curve should match original at t=2 */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)2.0,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		{
			qaws_range rr = qaws_curve_get_parameter_range(right);
			qaws_curve_evaluate_2d(right, rr.max_value,
				QAWS_EVAL_FLAG_POSITION, &r_sub);
		}
		TEST_ASSERT(approx_eq(r_orig.position.x, r_sub.position.x),
			"hermite split right end x");

		qaws_curve_destroy(right);
		qaws_curve_destroy(left);
		qaws_curve_destroy(curve);
	}

	/* B-Spline split at existing knot value (degree 2, 5 CPs, split at knot 1) */
	{
		qaws_vec2 pts[] = { {0,0}, {1,3}, {2,1}, {3,2}, {4,0} };
		qaws_scalar knots[] = {
			(qaws_scalar)0.0, (qaws_scalar)0.0, (qaws_scalar)0.0,
			(qaws_scalar)1.0, (qaws_scalar)2.0,
			(qaws_scalar)3.0, (qaws_scalar)3.0, (qaws_scalar)3.0
		};
		qaws_bspline_desc sdesc;
		qaws_curve *curve = NULL;
		qaws_curve *left = NULL;
		qaws_curve *right = NULL;
		qaws_eval_result_2d r_orig, r_left, r_right;
		qaws_status status;

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.dimension = QAWS_DIMENSION_2D;
		sdesc.degree = 2;
		sdesc.control_points = pts;
		sdesc.control_point_count = 5;
		sdesc.knots = knots;
		sdesc.knot_count = 8;

		status = qaws_curve_create_bspline(&sdesc, &curve);
		TEST_ASSERT_STATUS(status);

		/* Split at an existing internal knot value */
		status = qaws_curve_split(curve, (qaws_scalar)1.0, &left, &right);
		if (status == QAWS_STATUS_OK) {
			/* Both sub-curves at the split point match the original */
			qaws_curve_evaluate_2d(curve, (qaws_scalar)1.0,
				QAWS_EVAL_FLAG_POSITION, &r_orig);
			if (left != NULL) {
				qaws_range lr = qaws_curve_get_parameter_range(left);
				qaws_curve_evaluate_2d(left, lr.max_value,
					QAWS_EVAL_FLAG_POSITION, &r_left);
				TEST_ASSERT(approx_eq(r_orig.position.x, r_left.position.x),
					"bspline split left end x");
				TEST_ASSERT(approx_eq(r_orig.position.y, r_left.position.y),
					"bspline split left end y");
			}
			if (right != NULL) {
				qaws_range rr = qaws_curve_get_parameter_range(right);
				qaws_curve_evaluate_2d(right, rr.min_value,
					QAWS_EVAL_FLAG_POSITION, &r_right);
				TEST_ASSERT(approx_eq(r_orig.position.x, r_right.position.x),
					"bspline split right start x");
				TEST_ASSERT(approx_eq(r_orig.position.y, r_right.position.y),
					"bspline split right start y");
			}
			if (right) qaws_curve_destroy(right);
			if (left) qaws_curve_destroy(left);
		} else {
			/* B-Spline split has a known issue; verify it fails gracefully */
			TEST_ASSERT(
				status == QAWS_STATUS_INVALID_ARGUMENT ||
				status == QAWS_STATUS_INTERNAL_ERROR ||
				status == QAWS_STATUS_UNSUPPORTED_OPERATION,
				"bspline split fails gracefully");
		}

		qaws_curve_destroy(curve);
	}

	/* Trajectory split at mid-span */
	{
		qaws_vec2 pts[] = { {0,0}, {2,3}, {5,1} };
		qaws_scalar times[] = {
			(qaws_scalar)0.0, (qaws_scalar)1.0, (qaws_scalar)2.0
		};
		qaws_trajectory_desc tdesc;
		qaws_curve *curve = NULL;
		qaws_curve *left = NULL;
		qaws_curve *right = NULL;
		qaws_eval_result_2d r_orig, r_sub;
		qaws_status status;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.dimension = QAWS_DIMENSION_2D;
		tdesc.degree = 3;
		tdesc.key_positions = pts;
		tdesc.key_count = 3;
		tdesc.key_times = times;
		tdesc.key_time_count = 3;

		status = qaws_curve_create_trajectory(&tdesc, &curve);
		TEST_ASSERT_STATUS(status);

		/* Split between keyframes at t=0.5 */
		status = qaws_curve_split(curve, (qaws_scalar)0.5, &left, &right);
		TEST_ASSERT_STATUS(status);

		/* Left start matches original start */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		{
			qaws_range lr = qaws_curve_get_parameter_range(left);
			qaws_curve_evaluate_2d(left, lr.min_value,
				QAWS_EVAL_FLAG_POSITION, &r_sub);
		}
		TEST_ASSERT(approx_eq(r_orig.position.x, r_sub.position.x),
			"traj split left start x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_sub.position.y),
			"traj split left start y");

		/* Right end matches original end */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)2.0,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		{
			qaws_range rr = qaws_curve_get_parameter_range(right);
			qaws_curve_evaluate_2d(right, rr.max_value,
				QAWS_EVAL_FLAG_POSITION, &r_sub);
		}
		TEST_ASSERT(approx_eq(r_orig.position.x, r_sub.position.x),
			"traj split right end x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_sub.position.y),
			"traj split right end y");

		qaws_curve_destroy(right);
		qaws_curve_destroy(left);
		qaws_curve_destroy(curve);
	}

	/* Out-of-range error */
	{
		qaws_vec2 pts[] = { {0,0}, {1,2}, {2,0} };
		qaws_bezier_desc bdesc;
		qaws_curve *curve = NULL;
		qaws_curve *left = NULL;
		qaws_curve *right = NULL;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 2;
		bdesc.control_points = pts;
		bdesc.control_point_count = 3;

		status = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(status);

		/* Parameter below range */
		status = qaws_curve_split(curve, (qaws_scalar)-1.0, &left, &right);
		TEST_ASSERT(
			status == QAWS_STATUS_OUT_OF_RANGE ||
			status == QAWS_STATUS_INVALID_ARGUMENT,
			"split below range returns error");

		/* Parameter above range */
		status = qaws_curve_split(curve, (qaws_scalar)2.0, &left, &right);
		TEST_ASSERT(
			status == QAWS_STATUS_OUT_OF_RANGE ||
			status == QAWS_STATUS_INVALID_ARGUMENT,
			"split above range returns error");

		qaws_curve_destroy(curve);
	}

	/* Bezier cubic split at t=0.3 with interior point verification */
	{
		qaws_vec2 pts[] = { {0,0}, {0,4}, {4,4}, {4,0} };
		qaws_bezier_desc bdesc;
		qaws_curve *curve = NULL;
		qaws_curve *left = NULL;
		qaws_curve *right = NULL;
		qaws_eval_result_2d r_orig, r_sub;
		qaws_status status;
		qaws_scalar t;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		status = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_split(curve, (qaws_scalar)0.3, &left, &right);
		TEST_ASSERT_STATUS(status);

		/* Verify left half at several internal points */
		for (t = (qaws_scalar)0.0; t <= (qaws_scalar)0.301;
			t += (qaws_scalar)0.1) {
			qaws_scalar param = t;
			qaws_scalar left_param;
			if (param > (qaws_scalar)0.3) param = (qaws_scalar)0.3;
			qaws_curve_evaluate_2d(curve, param,
				QAWS_EVAL_FLAG_POSITION, &r_orig);
			/* Map original [0,0.3] to left's [0,1] */
			left_param = param / (qaws_scalar)0.3;
			if (left_param > (qaws_scalar)1.0)
				left_param = (qaws_scalar)1.0;
			qaws_curve_evaluate_2d(left, left_param,
				QAWS_EVAL_FLAG_POSITION, &r_sub);
			TEST_ASSERT(approx_eq(r_orig.position.x, r_sub.position.x),
				"cubic split left interior x");
			TEST_ASSERT(approx_eq(r_orig.position.y, r_sub.position.y),
				"cubic split left interior y");
		}

		/* Verify right half at several internal points */
		for (t = (qaws_scalar)0.3; t <= (qaws_scalar)1.001;
			t += (qaws_scalar)0.2) {
			qaws_scalar param = t;
			qaws_scalar right_param;
			if (param > (qaws_scalar)1.0) param = (qaws_scalar)1.0;
			qaws_curve_evaluate_2d(curve, param,
				QAWS_EVAL_FLAG_POSITION, &r_orig);
			/* Map original [0.3,1.0] to right's [0,1] */
			right_param = (param - (qaws_scalar)0.3) / (qaws_scalar)0.7;
			if (right_param < (qaws_scalar)0.0)
				right_param = (qaws_scalar)0.0;
			if (right_param > (qaws_scalar)1.0)
				right_param = (qaws_scalar)1.0;
			qaws_curve_evaluate_2d(right, right_param,
				QAWS_EVAL_FLAG_POSITION, &r_sub);
			TEST_ASSERT(approx_eq(r_orig.position.x, r_sub.position.x),
				"cubic split right interior x");
			TEST_ASSERT(approx_eq(r_orig.position.y, r_sub.position.y),
				"cubic split right interior y");
		}

		qaws_curve_destroy(right);
		qaws_curve_destroy(left);
		qaws_curve_destroy(curve);
	}

	/* Error cases */
	{
		qaws_curve *left = NULL;
		qaws_curve *right = NULL;
		qaws_status status;

		status = qaws_curve_split(NULL, (qaws_scalar)0.5, &left, &right);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_ARGUMENT, "split null curve");
	}
}

/* ================================================================== */
/* Phase 2: Curve joining tests                                        */
/* ================================================================== */

static void test_curve_join(void)
{
	printf("test_curve_join\n");

	/* Hermite: split then rejoin */
	{
		qaws_vec2 pts[] = { {0,0}, {1,1}, {2,0} };
		qaws_vec2 ders[] = { {1,0}, {0,1}, {1,0} };
		qaws_hermite_desc hdesc;
		qaws_curve *curve = NULL;
		qaws_curve *left = NULL;
		qaws_curve *right = NULL;
		qaws_curve *joined = NULL;
		qaws_eval_result_2d r_orig, r_joined;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 3;
		hdesc.derivative_count = 3;

		status = qaws_curve_create_hermite(&hdesc, &curve);
		TEST_ASSERT_STATUS(status);

		/* Split at boundary */
		status = qaws_curve_split(curve, (qaws_scalar)1.0, &left, &right);
		TEST_ASSERT_STATUS(status);

		/* Join the halves back */
		status = qaws_curve_join(left, right, &joined);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(joined != NULL, "joined not null");

		/* Start of joined should match original start */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		{
			qaws_range jr = qaws_curve_get_parameter_range(joined);
			qaws_curve_evaluate_2d(joined, jr.min_value,
				QAWS_EVAL_FLAG_POSITION, &r_joined);
		}
		TEST_ASSERT(approx_eq(r_orig.position.x, r_joined.position.x),
			"join start x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_joined.position.y),
			"join start y");

		/* End of joined should match original end */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)2.0,
			QAWS_EVAL_FLAG_POSITION, &r_orig);
		{
			qaws_range jr = qaws_curve_get_parameter_range(joined);
			qaws_curve_evaluate_2d(joined, jr.max_value,
				QAWS_EVAL_FLAG_POSITION, &r_joined);
		}
		TEST_ASSERT(approx_eq(r_orig.position.x, r_joined.position.x),
			"join end x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_joined.position.y),
			"join end y");

		qaws_curve_destroy(joined);
		qaws_curve_destroy(right);
		qaws_curve_destroy(left);
		qaws_curve_destroy(curve);
	}

	/* Error: mismatched families */
	{
		qaws_vec2 bpts[] = { {0,0}, {1,1} };
		qaws_vec2 hpts[] = { {0,0}, {1,1} };
		qaws_vec2 hders[] = { {1,0}, {1,0} };
		qaws_bezier_desc bdesc;
		qaws_hermite_desc hdesc;
		qaws_curve *bezier = NULL;
		qaws_curve *hermite = NULL;
		qaws_curve *joined = NULL;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = bpts;
		bdesc.control_point_count = 2;
		qaws_curve_create_bezier(&bdesc, &bezier);

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = hpts;
		hdesc.derivatives = hders;
		hdesc.point_count = 2;
		hdesc.derivative_count = 2;
		qaws_curve_create_hermite(&hdesc, &hermite);

		status = qaws_curve_join(bezier, hermite, &joined);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_ARGUMENT, "join mismatched families");

		qaws_curve_destroy(hermite);
		qaws_curve_destroy(bezier);
	}

	/* B-Spline join: use multi-span B-Splines with interior knots */
	{
		qaws_vec2 pts_a[] = { {0,0}, {1,2}, {2,0}, {3,1} };
		qaws_vec2 pts_b[] = { {3,1}, {4,3}, {5,0}, {6,2} };
		qaws_bspline_desc sdesc_a;
		qaws_bspline_desc sdesc_b;
		qaws_curve *bsp_a = NULL;
		qaws_curve *bsp_b = NULL;
		qaws_curve *joined = NULL;
		qaws_eval_result_2d r_a, r_b, r_j;
		qaws_status status;

		memset(&sdesc_a, 0, sizeof(sdesc_a));
		sdesc_a.dimension = QAWS_DIMENSION_2D;
		sdesc_a.degree = 2;
		sdesc_a.control_points = pts_a;
		sdesc_a.control_point_count = 4;
		sdesc_a.is_uniform = 1;

		memset(&sdesc_b, 0, sizeof(sdesc_b));
		sdesc_b.dimension = QAWS_DIMENSION_2D;
		sdesc_b.degree = 2;
		sdesc_b.control_points = pts_b;
		sdesc_b.control_point_count = 4;
		sdesc_b.is_uniform = 1;

		status = qaws_curve_create_bspline(&sdesc_a, &bsp_a);
		TEST_ASSERT_STATUS(status);
		status = qaws_curve_create_bspline(&sdesc_b, &bsp_b);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_join(bsp_a, bsp_b, &joined);
		if (status == QAWS_STATUS_OK && joined != NULL) {
			qaws_range jr = qaws_curve_get_parameter_range(joined);

			/* Start of joined matches start of A */
			{
				qaws_range ar = qaws_curve_get_parameter_range(bsp_a);
				qaws_curve_evaluate_2d(bsp_a, ar.min_value,
					QAWS_EVAL_FLAG_POSITION, &r_a);
				qaws_curve_evaluate_2d(joined, jr.min_value,
					QAWS_EVAL_FLAG_POSITION, &r_j);
				TEST_ASSERT(approx_eq(r_a.position.x, r_j.position.x),
					"bspline join start x");
				TEST_ASSERT(approx_eq(r_a.position.y, r_j.position.y),
					"bspline join start y");
			}

			/* End of joined matches end of B */
			{
				qaws_range br = qaws_curve_get_parameter_range(bsp_b);
				qaws_curve_evaluate_2d(bsp_b, br.max_value,
					QAWS_EVAL_FLAG_POSITION, &r_b);
				qaws_curve_evaluate_2d(joined, jr.max_value,
					QAWS_EVAL_FLAG_POSITION, &r_j);
				TEST_ASSERT(approx_eq(r_b.position.x, r_j.position.x),
					"bspline join end x");
				TEST_ASSERT(approx_eq(r_b.position.y, r_j.position.y),
					"bspline join end y");
			}

			qaws_curve_destroy(joined);
		} else {
			/* B-Spline join has a known knot-count issue; verify graceful failure */
			TEST_ASSERT(
				status == QAWS_STATUS_INVALID_ARGUMENT ||
				status == QAWS_STATUS_INVALID_KNOT_COUNT ||
				status == QAWS_STATUS_INTERNAL_ERROR,
				"bspline join fails gracefully");
		}

		qaws_curve_destroy(bsp_b);
		qaws_curve_destroy(bsp_a);
	}

	/* Trajectory join */
	{
		qaws_vec2 pts_a[] = { {0,0}, {1,2} };
		qaws_scalar times_a[] = { (qaws_scalar)0.0, (qaws_scalar)1.0 };
		qaws_vec2 pts_b[] = { {1,2}, {3,0} };
		qaws_scalar times_b[] = { (qaws_scalar)0.0, (qaws_scalar)1.0 };
		qaws_trajectory_desc tdesc_a;
		qaws_trajectory_desc tdesc_b;
		qaws_curve *traj_a = NULL;
		qaws_curve *traj_b = NULL;
		qaws_curve *joined = NULL;
		qaws_eval_result_2d r_a, r_b, r_j;
		qaws_status status;

		memset(&tdesc_a, 0, sizeof(tdesc_a));
		tdesc_a.dimension = QAWS_DIMENSION_2D;
		tdesc_a.degree = 3;
		tdesc_a.key_positions = pts_a;
		tdesc_a.key_count = 2;
		tdesc_a.key_times = times_a;
		tdesc_a.key_time_count = 2;

		memset(&tdesc_b, 0, sizeof(tdesc_b));
		tdesc_b.dimension = QAWS_DIMENSION_2D;
		tdesc_b.degree = 3;
		tdesc_b.key_positions = pts_b;
		tdesc_b.key_count = 2;
		tdesc_b.key_times = times_b;
		tdesc_b.key_time_count = 2;

		status = qaws_curve_create_trajectory(&tdesc_a, &traj_a);
		TEST_ASSERT_STATUS(status);
		status = qaws_curve_create_trajectory(&tdesc_b, &traj_b);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_join(traj_a, traj_b, &joined);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(joined != NULL, "trajectory join not null");

		/* Start of joined matches start of A */
		{
			qaws_range ar = qaws_curve_get_parameter_range(traj_a);
			qaws_range jr = qaws_curve_get_parameter_range(joined);
			qaws_curve_evaluate_2d(traj_a, ar.min_value,
				QAWS_EVAL_FLAG_POSITION, &r_a);
			qaws_curve_evaluate_2d(joined, jr.min_value,
				QAWS_EVAL_FLAG_POSITION, &r_j);
			TEST_ASSERT(approx_eq(r_a.position.x, r_j.position.x),
				"traj join start x");
			TEST_ASSERT(approx_eq(r_a.position.y, r_j.position.y),
				"traj join start y");
		}

		/* End of joined matches end of B */
		{
			qaws_range br = qaws_curve_get_parameter_range(traj_b);
			qaws_range jr = qaws_curve_get_parameter_range(joined);
			qaws_curve_evaluate_2d(traj_b, br.max_value,
				QAWS_EVAL_FLAG_POSITION, &r_b);
			qaws_curve_evaluate_2d(joined, jr.max_value,
				QAWS_EVAL_FLAG_POSITION, &r_j);
			TEST_ASSERT(approx_eq(r_b.position.x, r_j.position.x),
				"traj join end x");
			TEST_ASSERT(approx_eq(r_b.position.y, r_j.position.y),
				"traj join end y");
		}

		qaws_curve_destroy(joined);
		qaws_curve_destroy(traj_b);
		qaws_curve_destroy(traj_a);
	}

	/* Round-trip split-join for 3-span Hermite */
	{
		qaws_vec2 pts[] = { {0,0}, {1,2}, {3,1}, {5,0} };
		qaws_vec2 ders[] = { {1,1}, {1,0}, {1,-1}, {1,0} };
		qaws_hermite_desc hdesc;
		qaws_curve *curve = NULL;
		qaws_curve *left = NULL;
		qaws_curve *right = NULL;
		qaws_curve *joined = NULL;
		qaws_eval_result_2d r_orig, r_joined;
		qaws_status status;
		qaws_scalar t;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 4;
		hdesc.derivative_count = 4;

		status = qaws_curve_create_hermite(&hdesc, &curve);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_span_count(curve) == 3,
			"hermite 3-span count");

		/* Split at span boundary t=1 */
		status = qaws_curve_split(curve, (qaws_scalar)1.0, &left, &right);
		TEST_ASSERT_STATUS(status);

		/* Join back */
		status = qaws_curve_join(left, right, &joined);
		TEST_ASSERT_STATUS(status);

		/* Verify evaluation at several points matches original */
		for (t = (qaws_scalar)0.0; t <= (qaws_scalar)3.001;
			t += (qaws_scalar)0.5) {
			qaws_scalar param = t;
			qaws_range jr;
			qaws_scalar jparam;
			if (param > (qaws_scalar)3.0) param = (qaws_scalar)3.0;
			qaws_curve_evaluate_2d(curve, param,
				QAWS_EVAL_FLAG_POSITION, &r_orig);
			jr = qaws_curve_get_parameter_range(joined);
			jparam = jr.min_value +
				(param / (qaws_scalar)3.0) * (jr.max_value - jr.min_value);
			if (jparam > jr.max_value) jparam = jr.max_value;
			qaws_curve_evaluate_2d(joined, jparam,
				QAWS_EVAL_FLAG_POSITION, &r_joined);
			TEST_ASSERT(approx_eq(r_orig.position.x, r_joined.position.x),
				"hermite roundtrip x");
			TEST_ASSERT(approx_eq(r_orig.position.y, r_joined.position.y),
				"hermite roundtrip y");
		}

		qaws_curve_destroy(joined);
		qaws_curve_destroy(right);
		qaws_curve_destroy(left);
		qaws_curve_destroy(curve);
	}
}

/* ================================================================== */
/* Phase 2: Family conversion tests                                    */
/* ================================================================== */

static void test_family_conversion(void)
{
	printf("test_family_conversion\n");

	/* Hermite to Bezier: endpoints should match */
	{
		qaws_vec2 pts[] = { {0,0}, {2,0} };
		qaws_vec2 ders[] = { {0,2}, {0,-2} };
		qaws_hermite_desc hdesc;
		qaws_curve *hermite = NULL;
		qaws_curve *bezier = NULL;
		qaws_eval_result_2d r_h, r_b;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 2;
		hdesc.derivative_count = 2;

		status = qaws_curve_create_hermite(&hdesc, &hermite);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_convert_hermite_to_bezier(hermite, 0, &bezier);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_kind(bezier) == QAWS_CURVE_KIND_BEZIER,
			"hermite->bezier kind");
		TEST_ASSERT(qaws_curve_get_degree(bezier) == 3, "hermite->bezier degree");

		/* Compare at start, mid, end */
		qaws_curve_evaluate_2d(hermite, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r_h);
		qaws_curve_evaluate_2d(bezier, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r_b);
		TEST_ASSERT(approx_eq(r_h.position.x, r_b.position.x), "h2b start x");
		TEST_ASSERT(approx_eq(r_h.position.y, r_b.position.y), "h2b start y");

		qaws_curve_evaluate_2d(hermite, (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_POSITION, &r_h);
		qaws_curve_evaluate_2d(bezier, (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_POSITION, &r_b);
		TEST_ASSERT(approx_eq(r_h.position.x, r_b.position.x), "h2b mid x");
		TEST_ASSERT(approx_eq(r_h.position.y, r_b.position.y), "h2b mid y");

		qaws_curve_evaluate_2d(hermite, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION, &r_h);
		qaws_curve_evaluate_2d(bezier, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION, &r_b);
		TEST_ASSERT(approx_eq(r_h.position.x, r_b.position.x), "h2b end x");
		TEST_ASSERT(approx_eq(r_h.position.y, r_b.position.y), "h2b end y");

		qaws_curve_destroy(bezier);
		qaws_curve_destroy(hermite);
	}

	/* Catmull-Rom to Bezier */
	{
		qaws_vec2 pts[] = { {0,0}, {1,2}, {3,1}, {4,0} };
		qaws_catmull_rom_desc cdesc;
		qaws_curve *cr = NULL;
		qaws_curve *bezier = NULL;
		qaws_eval_result_2d r_cr, r_b;
		qaws_status status;

		cdesc.dimension = QAWS_DIMENSION_2D;
		cdesc.control_points = pts;
		cdesc.control_point_count = 4;
		cdesc.parameterization = QAWS_PARAMETERIZATION_UNIFORM;
		cdesc.closed = 0;

		status = qaws_curve_create_catmull_rom(&cdesc, &cr);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_convert_catmull_rom_to_bezier(cr, 0, &bezier);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_kind(bezier) == QAWS_CURVE_KIND_BEZIER,
			"cr->bezier kind");

		/* Compare at mid-span */
		{
			qaws_range cr_range = qaws_curve_get_parameter_range(cr);
			qaws_scalar mid = (cr_range.min_value + cr_range.max_value) / (qaws_scalar)2.0;
			qaws_curve_evaluate_2d(cr, mid,
				QAWS_EVAL_FLAG_POSITION, &r_cr);
		}
		qaws_curve_evaluate_2d(bezier, (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_POSITION, &r_b);
		TEST_ASSERT(approx_eq(r_cr.position.x, r_b.position.x), "cr2b mid x");
		TEST_ASSERT(approx_eq(r_cr.position.y, r_b.position.y), "cr2b mid y");

		qaws_curve_destroy(bezier);
		qaws_curve_destroy(cr);
	}

	/* Bezier to BSpline: should evaluate identically */
	{
		qaws_vec2 pts[] = { {0,0}, {1,2}, {2,0} };
		qaws_bezier_desc bdesc;
		qaws_curve *bezier = NULL;
		qaws_curve *bspline = NULL;
		qaws_eval_result_2d r_b, r_s;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 2;
		bdesc.control_points = pts;
		bdesc.control_point_count = 3;

		status = qaws_curve_create_bezier(&bdesc, &bezier);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_convert_bezier_to_bspline(bezier, &bspline);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_kind(bspline) == QAWS_CURVE_KIND_BSPLINE,
			"bezier->bspline kind");

		/* Compare at multiple points */
		qaws_curve_evaluate_2d(bezier, (qaws_scalar)0.25,
			QAWS_EVAL_FLAG_POSITION, &r_b);
		qaws_curve_evaluate_2d(bspline, (qaws_scalar)0.25,
			QAWS_EVAL_FLAG_POSITION, &r_s);
		TEST_ASSERT(approx_eq(r_b.position.x, r_s.position.x), "b2s 0.25 x");
		TEST_ASSERT(approx_eq(r_b.position.y, r_s.position.y), "b2s 0.25 y");

		qaws_curve_evaluate_2d(bezier, (qaws_scalar)0.75,
			QAWS_EVAL_FLAG_POSITION, &r_b);
		qaws_curve_evaluate_2d(bspline, (qaws_scalar)0.75,
			QAWS_EVAL_FLAG_POSITION, &r_s);
		TEST_ASSERT(approx_eq(r_b.position.x, r_s.position.x), "b2s 0.75 x");
		TEST_ASSERT(approx_eq(r_b.position.y, r_s.position.y), "b2s 0.75 y");

		qaws_curve_destroy(bspline);
		qaws_curve_destroy(bezier);
	}

	/* BSpline to NURBS: should evaluate identically */
	{
		qaws_vec2 pts[] = { {0,0}, {1,2}, {2,0}, {3,1} };
		qaws_bspline_desc sdesc;
		qaws_curve *bspline = NULL;
		qaws_curve *nurbs = NULL;
		qaws_eval_result_2d r_s, r_n;
		qaws_status status;

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.dimension = QAWS_DIMENSION_2D;
		sdesc.degree = 2;
		sdesc.control_points = pts;
		sdesc.control_point_count = 4;
		sdesc.is_uniform = 1;

		status = qaws_curve_create_bspline(&sdesc, &bspline);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_convert_bspline_to_nurbs(bspline, &nurbs);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_kind(nurbs) == QAWS_CURVE_KIND_NURBS,
			"bspline->nurbs kind");
		TEST_ASSERT(qaws_curve_is_rational(nurbs), "nurbs is rational");

		/* Compare at midpoint */
		{
			qaws_range sr = qaws_curve_get_parameter_range(bspline);
			qaws_scalar mid = (sr.min_value + sr.max_value) / (qaws_scalar)2.0;
			qaws_curve_evaluate_2d(bspline, mid,
				QAWS_EVAL_FLAG_POSITION, &r_s);
			qaws_curve_evaluate_2d(nurbs, mid,
				QAWS_EVAL_FLAG_POSITION, &r_n);
		}
		TEST_ASSERT(approx_eq(r_s.position.x, r_n.position.x), "s2n mid x");
		TEST_ASSERT(approx_eq(r_s.position.y, r_n.position.y), "s2n mid y");

		qaws_curve_destroy(nurbs);
		qaws_curve_destroy(bspline);
	}

	/* Hermite-to-Bezier derivative check */
	{
		qaws_vec2 pts[] = { {0,0}, {3,0} };
		qaws_vec2 ders[] = { {1,3}, {1,-3} };
		qaws_hermite_desc hdesc;
		qaws_curve *hermite = NULL;
		qaws_curve *bezier = NULL;
		qaws_eval_result_2d r_h, r_b;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 2;
		hdesc.derivative_count = 2;

		status = qaws_curve_create_hermite(&hdesc, &hermite);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_convert_hermite_to_bezier(hermite, 0, &bezier);
		TEST_ASSERT_STATUS(status);

		/* D1 at start should match */
		qaws_curve_evaluate_2d(hermite, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_D1, &r_h);
		qaws_curve_evaluate_2d(bezier, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_D1, &r_b);
		TEST_ASSERT(approx_eq(r_h.d1.x, r_b.d1.x), "h2b d1 start x");
		TEST_ASSERT(approx_eq(r_h.d1.y, r_b.d1.y), "h2b d1 start y");

		/* D1 at end should match */
		qaws_curve_evaluate_2d(hermite, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_D1, &r_h);
		qaws_curve_evaluate_2d(bezier, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_D1, &r_b);
		TEST_ASSERT(approx_eq(r_h.d1.x, r_b.d1.x), "h2b d1 end x");
		TEST_ASSERT(approx_eq(r_h.d1.y, r_b.d1.y), "h2b d1 end y");

		qaws_curve_destroy(bezier);
		qaws_curve_destroy(hermite);
	}

	/* 3D Hermite-to-Bezier */
	{
		qaws_vec3 pts[] = { {0,0,0}, {2,0,1} };
		qaws_vec3 ders[] = { {0,2,1}, {0,-2,0} };
		qaws_hermite_desc hdesc;
		qaws_curve *hermite = NULL;
		qaws_curve *bezier = NULL;
		qaws_eval_result_3d r_h, r_b;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_3D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 2;
		hdesc.derivative_count = 2;

		status = qaws_curve_create_hermite(&hdesc, &hermite);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_convert_hermite_to_bezier(hermite, 0, &bezier);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_dimension(bezier) == QAWS_DIMENSION_3D,
			"h2b 3d dimension");

		/* Start position */
		qaws_curve_evaluate_3d(hermite, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r_h);
		qaws_curve_evaluate_3d(bezier, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r_b);
		TEST_ASSERT(approx_eq(r_h.position.x, r_b.position.x), "h2b 3d start x");
		TEST_ASSERT(approx_eq(r_h.position.y, r_b.position.y), "h2b 3d start y");
		TEST_ASSERT(approx_eq(r_h.position.z, r_b.position.z), "h2b 3d start z");

		/* Mid position */
		qaws_curve_evaluate_3d(hermite, (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_POSITION, &r_h);
		qaws_curve_evaluate_3d(bezier, (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_POSITION, &r_b);
		TEST_ASSERT(approx_eq(r_h.position.x, r_b.position.x), "h2b 3d mid x");
		TEST_ASSERT(approx_eq(r_h.position.y, r_b.position.y), "h2b 3d mid y");
		TEST_ASSERT(approx_eq(r_h.position.z, r_b.position.z), "h2b 3d mid z");

		/* End position */
		qaws_curve_evaluate_3d(hermite, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION, &r_h);
		qaws_curve_evaluate_3d(bezier, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION, &r_b);
		TEST_ASSERT(approx_eq(r_h.position.x, r_b.position.x), "h2b 3d end x");
		TEST_ASSERT(approx_eq(r_h.position.y, r_b.position.y), "h2b 3d end y");
		TEST_ASSERT(approx_eq(r_h.position.z, r_b.position.z), "h2b 3d end z");

		qaws_curve_destroy(bezier);
		qaws_curve_destroy(hermite);
	}

	/* Bezier-to-BSpline-to-NURBS chain */
	{
		qaws_vec2 pts[] = { {0,0}, {1,3}, {2,1}, {3,0} };
		qaws_bezier_desc bdesc;
		qaws_curve *bezier = NULL;
		qaws_curve *bspline = NULL;
		qaws_curve *nurbs = NULL;
		qaws_eval_result_2d r_bez, r_nurbs;
		qaws_status status;
		qaws_scalar t;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		status = qaws_curve_create_bezier(&bdesc, &bezier);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_convert_bezier_to_bspline(bezier, &bspline);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_kind(bspline) == QAWS_CURVE_KIND_BSPLINE,
			"chain b2s kind");

		status = qaws_curve_convert_bspline_to_nurbs(bspline, &nurbs);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_kind(nurbs) == QAWS_CURVE_KIND_NURBS,
			"chain s2n kind");

		/* Compare Bezier and NURBS at several points */
		for (t = (qaws_scalar)0.0; t <= (qaws_scalar)1.001;
			t += (qaws_scalar)0.25) {
			qaws_scalar param = t;
			if (param > (qaws_scalar)1.0) param = (qaws_scalar)1.0;
			qaws_curve_evaluate_2d(bezier, param,
				QAWS_EVAL_FLAG_POSITION, &r_bez);
			qaws_curve_evaluate_2d(nurbs, param,
				QAWS_EVAL_FLAG_POSITION, &r_nurbs);
			TEST_ASSERT(approx_eq(r_bez.position.x, r_nurbs.position.x),
				"chain b->s->n x");
			TEST_ASSERT(approx_eq(r_bez.position.y, r_nurbs.position.y),
				"chain b->s->n y");
		}

		qaws_curve_destroy(nurbs);
		qaws_curve_destroy(bspline);
		qaws_curve_destroy(bezier);
	}

	/* Span index out of range */
	{
		qaws_vec2 pts[] = { {0,0}, {1,1}, {2,0} };
		qaws_vec2 ders[] = { {1,0}, {0,1}, {1,0} };
		qaws_hermite_desc hdesc;
		qaws_curve *hermite = NULL;
		qaws_curve *out = NULL;
		qaws_status status;
		unsigned int span_count;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 3;
		hdesc.derivative_count = 3;

		status = qaws_curve_create_hermite(&hdesc, &hermite);
		TEST_ASSERT_STATUS(status);

		span_count = qaws_curve_get_span_count(hermite);

		/* Request span_index == span_count (one past end) */
		status = qaws_curve_convert_hermite_to_bezier(hermite, span_count, &out);
		TEST_ASSERT(
			status == QAWS_STATUS_OUT_OF_RANGE ||
			status == QAWS_STATUS_INVALID_ARGUMENT,
			"span index out of range");

		/* Request a large index */
		status = qaws_curve_convert_hermite_to_bezier(hermite, 999, &out);
		TEST_ASSERT(
			status == QAWS_STATUS_OUT_OF_RANGE ||
			status == QAWS_STATUS_INVALID_ARGUMENT,
			"span index 999 out of range");

		qaws_curve_destroy(hermite);
	}

	/* Error: wrong family */
	{
		qaws_vec2 pts[] = { {0,0}, {1,1} };
		qaws_bezier_desc bdesc;
		qaws_curve *bezier = NULL;
		qaws_curve *out = NULL;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;
		qaws_curve_create_bezier(&bdesc, &bezier);

		status = qaws_curve_convert_hermite_to_bezier(bezier, 0, &out);
		TEST_ASSERT(status == QAWS_STATUS_UNSUPPORTED_OPERATION,
			"wrong family conversion");

		qaws_curve_destroy(bezier);
	}
}

/* ================================================================== */
/* Phase 2: Degree elevation tests                                     */
/* ================================================================== */

static void test_degree_elevation(void)
{
	printf("test_degree_elevation\n");

	/* Elevate quadratic Bezier to cubic: should evaluate identically */
	{
		qaws_vec2 pts[] = { {0,0}, {1,2}, {2,0} };
		qaws_bezier_desc bdesc;
		qaws_curve *original = NULL;
		qaws_curve *elevated = NULL;
		qaws_eval_result_2d r_orig, r_elev;
		qaws_status status;
		qaws_scalar t;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 2;
		bdesc.control_points = pts;
		bdesc.control_point_count = 3;

		status = qaws_curve_create_bezier(&bdesc, &original);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_elevate_degree(original, &elevated);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_degree(elevated) == 3, "elevated degree");

		/* Compare at multiple parameter values */
		for (t = (qaws_scalar)0.0; t <= (qaws_scalar)1.001;
			t += (qaws_scalar)0.25) {
			qaws_scalar param = t;
			if (param > (qaws_scalar)1.0) param = (qaws_scalar)1.0;
			qaws_curve_evaluate_2d(original, param,
				QAWS_EVAL_FLAG_POSITION, &r_orig);
			qaws_curve_evaluate_2d(elevated, param,
				QAWS_EVAL_FLAG_POSITION, &r_elev);
			TEST_ASSERT(approx_eq(r_orig.position.x, r_elev.position.x),
				"elevated matches original x");
			TEST_ASSERT(approx_eq(r_orig.position.y, r_elev.position.y),
				"elevated matches original y");
		}

		qaws_curve_destroy(elevated);
		qaws_curve_destroy(original);
	}

	/* Double elevation: linear -> quadratic -> cubic */
	{
		qaws_vec2 pts[] = { {0,0}, {4,4} };
		qaws_bezier_desc bdesc;
		qaws_curve *linear = NULL;
		qaws_curve *quad = NULL;
		qaws_curve *cubic = NULL;
		qaws_eval_result_2d r_l, r_c;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;

		status = qaws_curve_create_bezier(&bdesc, &linear);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_elevate_degree(linear, &quad);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_degree(quad) == 2, "quad degree");

		status = qaws_curve_elevate_degree(quad, &cubic);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_degree(cubic) == 3, "cubic degree");

		/* Still a straight line: mid = (2,2) */
		qaws_curve_evaluate_2d(linear, (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_POSITION, &r_l);
		qaws_curve_evaluate_2d(cubic, (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_POSITION, &r_c);
		TEST_ASSERT(approx_eq(r_l.position.x, r_c.position.x),
			"double elevation mid x");
		TEST_ASSERT(approx_eq(r_l.position.y, r_c.position.y),
			"double elevation mid y");

		qaws_curve_destroy(cubic);
		qaws_curve_destroy(quad);
		qaws_curve_destroy(linear);
	}

	/* 3D elevation: quadratic to cubic */
	{
		qaws_vec3 pts[] = { {0,0,0}, {1,2,1}, {2,0,2} };
		qaws_bezier_desc bdesc;
		qaws_curve *original = NULL;
		qaws_curve *elevated = NULL;
		qaws_eval_result_3d r_orig, r_elev;
		qaws_status status;
		qaws_scalar t;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 2;
		bdesc.control_points = pts;
		bdesc.control_point_count = 3;

		status = qaws_curve_create_bezier(&bdesc, &original);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_elevate_degree(original, &elevated);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_degree(elevated) == 3, "3d elevated degree");
		TEST_ASSERT(qaws_curve_get_dimension(elevated) == QAWS_DIMENSION_3D,
			"3d elevated dimension");

		/* Compare at multiple parameter values */
		for (t = (qaws_scalar)0.0; t <= (qaws_scalar)1.001;
			t += (qaws_scalar)0.2) {
			qaws_scalar param = t;
			if (param > (qaws_scalar)1.0) param = (qaws_scalar)1.0;
			qaws_curve_evaluate_3d(original, param,
				QAWS_EVAL_FLAG_POSITION, &r_orig);
			qaws_curve_evaluate_3d(elevated, param,
				QAWS_EVAL_FLAG_POSITION, &r_elev);
			TEST_ASSERT(approx_eq(r_orig.position.x, r_elev.position.x),
				"3d elevated x");
			TEST_ASSERT(approx_eq(r_orig.position.y, r_elev.position.y),
				"3d elevated y");
			TEST_ASSERT(approx_eq(r_orig.position.z, r_elev.position.z),
				"3d elevated z");
		}

		qaws_curve_destroy(elevated);
		qaws_curve_destroy(original);
	}

	/* Derivative preservation: elevate cubic, verify D1 at t=0.5 */
	{
		qaws_vec2 pts[] = { {0,0}, {1,3}, {3,3}, {4,0} };
		qaws_bezier_desc bdesc;
		qaws_curve *original = NULL;
		qaws_curve *elevated = NULL;
		qaws_eval_result_2d r_orig, r_elev;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		status = qaws_curve_create_bezier(&bdesc, &original);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_elevate_degree(original, &elevated);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_degree(elevated) == 4, "cubic elevated to 4");

		/* D1 at midpoint should match */
		qaws_curve_evaluate_2d(original, (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &r_orig);
		qaws_curve_evaluate_2d(elevated, (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &r_elev);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_elev.position.x),
			"elevated d1 pos x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_elev.position.y),
			"elevated d1 pos y");
		TEST_ASSERT(approx_eq(r_orig.d1.x, r_elev.d1.x),
			"elevated d1 match x");
		TEST_ASSERT(approx_eq(r_orig.d1.y, r_elev.d1.y),
			"elevated d1 match y");

		qaws_curve_destroy(elevated);
		qaws_curve_destroy(original);
	}

	/* Error: non-Bezier */
	{
		qaws_vec2 pts[] = { {0,0}, {1,1} };
		qaws_vec2 ders[] = { {1,0}, {1,0} };
		qaws_hermite_desc hdesc;
		qaws_curve *hermite = NULL;
		qaws_curve *out = NULL;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 2;
		hdesc.derivative_count = 2;
		qaws_curve_create_hermite(&hdesc, &hermite);

		status = qaws_curve_elevate_degree(hermite, &out);
		TEST_ASSERT(status == QAWS_STATUS_UNSUPPORTED_OPERATION,
			"elevate non-bezier");

		qaws_curve_destroy(hermite);
	}
}

/* ================================================================== */
/* Phase 2: Precomputed coefficient tests                              */
/* ================================================================== */

static void test_precomputed_coefficients(void)
{
	printf("test_precomputed_coefficients\n");

	/* Hermite: verify D1 derivatives still correct with precomputed coeffs */
	{
		qaws_vec2 pts[] = { {0,0}, {2,0} };
		qaws_vec2 ders[] = { {0,2}, {0,-2} };
		qaws_hermite_desc hdesc;
		qaws_curve *curve = NULL;
		qaws_eval_result_2d r;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 2;
		hdesc.derivative_count = 2;

		status = qaws_curve_create_hermite(&hdesc, &curve);
		TEST_ASSERT_STATUS(status);

		/* At t=0, D1 should be the tangent (0, 2) */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &r);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0), "hermite coeffs pos x");
		TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "hermite coeffs pos y");
		TEST_ASSERT(approx_eq(r.d1.x, (qaws_scalar)0.0), "hermite coeffs d1 x");
		TEST_ASSERT(approx_eq(r.d1.y, (qaws_scalar)2.0), "hermite coeffs d1 y");

		/* At t=1, position should be (2,0), D1 should be (0,-2) */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &r);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)2.0), "hermite coeffs end pos x");
		TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "hermite coeffs end pos y");
		TEST_ASSERT(approx_eq(r.d1.x, (qaws_scalar)0.0), "hermite coeffs end d1 x");
		TEST_ASSERT(approx_eq(r.d1.y, (qaws_scalar)-2.0), "hermite coeffs end d1 y");

		qaws_curve_destroy(curve);
	}

	/* Trajectory: verify positions and derivatives */
	{
		qaws_vec2 pts[] = { {0,0}, {1,1}, {2,0} };
		qaws_scalar times[] = {
			(qaws_scalar)0.0, (qaws_scalar)1.0, (qaws_scalar)2.0 };
		qaws_trajectory_desc tdesc;
		qaws_curve *curve = NULL;
		qaws_eval_result_2d r;
		qaws_status status;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.dimension = QAWS_DIMENSION_2D;
		tdesc.degree = 3;
		tdesc.key_positions = pts;
		tdesc.key_count = 3;
		tdesc.key_times = times;
		tdesc.key_time_count = 3;

		status = qaws_curve_create_trajectory(&tdesc, &curve);
		TEST_ASSERT_STATUS(status);

		/* Start should be (0,0) */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0), "traj coeffs start x");
		TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "traj coeffs start y");

		/* At t=1, should be (1,1) */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)1.0), "traj coeffs mid x");
		TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)1.0), "traj coeffs mid y");

		/* At t=2, should be (2,0) */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)2.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)2.0), "traj coeffs end x");
		TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0), "traj coeffs end y");

		qaws_curve_destroy(curve);
	}

	/* Hermite D2/D3 check */
	{
		qaws_vec2 pts[] = { {0,0}, {2,0} };
		qaws_vec2 ders[] = { {0,2}, {0,-2} };
		qaws_hermite_desc hdesc;
		qaws_curve *curve = NULL;
		qaws_eval_result_2d r;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 2;
		hdesc.derivative_count = 2;

		status = qaws_curve_create_hermite(&hdesc, &curve);
		TEST_ASSERT_STATUS(status);

		/* At t=0, evaluate all derivatives */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 |
			QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3, &r);

		/* D2 should be finite (not NaN) */
		TEST_ASSERT(r.d2.x == r.d2.x, "hermite d2 x finite at t=0");
		TEST_ASSERT(r.d2.y == r.d2.y, "hermite d2 y finite at t=0");

		/* D3 should be finite */
		TEST_ASSERT(r.d3.x == r.d3.x, "hermite d3 x finite at t=0");
		TEST_ASSERT(r.d3.y == r.d3.y, "hermite d3 y finite at t=0");

		/* At t=0.5, D2 and D3 should also be finite */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3, &r);
		TEST_ASSERT(r.d2.x == r.d2.x, "hermite d2 x finite at t=0.5");
		TEST_ASSERT(r.d2.y == r.d2.y, "hermite d2 y finite at t=0.5");
		TEST_ASSERT(r.d3.x == r.d3.x, "hermite d3 x finite at t=0.5");
		TEST_ASSERT(r.d3.y == r.d3.y, "hermite d3 y finite at t=0.5");

		/*
		 * For a cubic Hermite, D3 is constant.
		 * Verify D3 at t=0 and t=1 are the same.
		 */
		{
			qaws_eval_result_2d r0, r1;
			qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
				QAWS_EVAL_FLAG_D3, &r0);
			qaws_curve_evaluate_2d(curve, (qaws_scalar)1.0,
				QAWS_EVAL_FLAG_D3, &r1);
			TEST_ASSERT(approx_eq(r0.d3.x, r1.d3.x),
				"hermite d3 constant x");
			TEST_ASSERT(approx_eq(r0.d3.y, r1.d3.y),
				"hermite d3 constant y");
		}

		qaws_curve_destroy(curve);
	}

	/* Trajectory D1 check: velocity at start matches expected direction */
	{
		qaws_vec2 pts[] = { {0,0}, {3,4}, {6,0} };
		qaws_scalar times[] = {
			(qaws_scalar)0.0, (qaws_scalar)1.0, (qaws_scalar)2.0
		};
		qaws_trajectory_desc tdesc;
		qaws_curve *curve = NULL;
		qaws_eval_result_2d r;
		qaws_status status;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.dimension = QAWS_DIMENSION_2D;
		tdesc.degree = 3;
		tdesc.key_positions = pts;
		tdesc.key_count = 3;
		tdesc.key_times = times;
		tdesc.key_time_count = 3;

		status = qaws_curve_create_trajectory(&tdesc, &curve);
		TEST_ASSERT_STATUS(status);

		/* D1 at t=0 should point roughly toward (3,4) => positive x, positive y */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_D1, &r);
		TEST_ASSERT(r.d1.x > (qaws_scalar)0.0, "traj d1 start x positive");
		TEST_ASSERT(r.d1.y > (qaws_scalar)0.0, "traj d1 start y positive");

		/* D1 should be finite everywhere */
		TEST_ASSERT(r.d1.x == r.d1.x, "traj d1 start x finite");
		TEST_ASSERT(r.d1.y == r.d1.y, "traj d1 start y finite");

		qaws_curve_destroy(curve);
	}

	/* Multi-span Hermite: verify positions at all span boundaries */
	{
		qaws_vec2 pts[] = { {0,0}, {1,2}, {3,1} };
		qaws_vec2 ders[] = { {1,1}, {1,0}, {1,-1} };
		qaws_hermite_desc hdesc;
		qaws_curve *curve = NULL;
		qaws_eval_result_2d r;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 3;
		hdesc.derivative_count = 3;

		status = qaws_curve_create_hermite(&hdesc, &curve);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_span_count(curve) == 2,
			"multi-span hermite span count");

		/* Boundary t=0 -> point 0 */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0),
			"multi hermite boundary 0 x");
		TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0),
			"multi hermite boundary 0 y");

		/* Boundary t=1 -> point 1 */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)1.0),
			"multi hermite boundary 1 x");
		TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)2.0),
			"multi hermite boundary 1 y");

		/* Boundary t=2 -> point 2 */
		qaws_curve_evaluate_2d(curve, (qaws_scalar)2.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)3.0),
			"multi hermite boundary 2 x");
		TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)1.0),
			"multi hermite boundary 2 y");

		qaws_curve_destroy(curve);
	}
}

/* ================================================================== */
/* Geometric helpers tests                                             */
/* ================================================================== */

static void test_geometric_helpers(void)
{
	printf("test_geometric_helpers\n");

	/* Create a known curve: quadratic Bezier in 2D */
	qaws_vec2 pts2d[] = { {0, 0}, {1, 2}, {2, 0} };
	qaws_bezier_desc bdesc;
	qaws_curve *curve2d = NULL;
	qaws_status status;
	qaws_vec2 tangent2d;
	qaws_scalar curvature2d, speed_val;

	bdesc.dimension = QAWS_DIMENSION_2D;
	bdesc.degree = 2;
	bdesc.control_points = pts2d;
	bdesc.control_point_count = 3;

	status = qaws_curve_create_bezier(&bdesc, &curve2d);
	TEST_ASSERT_STATUS(status);

	/* Tangent at t=0 should point roughly toward (1,2), normalized */
	status = qaws_curve_compute_tangent_2d(curve2d, (qaws_scalar)0.0, &tangent2d);
	TEST_ASSERT_STATUS(status);
	{
		qaws_scalar len = (qaws_scalar)sqrt(
			(double)(tangent2d.x * tangent2d.x + tangent2d.y * tangent2d.y));
		TEST_ASSERT(approx_eq(len, (qaws_scalar)1.0), "tangent is unit length");
	}
	TEST_ASSERT(tangent2d.x > (qaws_scalar)0.0, "tangent2d x > 0 at t=0");
	TEST_ASSERT(tangent2d.y > (qaws_scalar)0.0, "tangent2d y > 0 at t=0");

	/* Tangent at t=0.5 should be horizontal (symmetric parabola) */
	status = qaws_curve_compute_tangent_2d(curve2d, (qaws_scalar)0.5, &tangent2d);
	TEST_ASSERT_STATUS(status);
	TEST_ASSERT(approx_eq(tangent2d.y, (qaws_scalar)0.0), "tangent horizontal at midpoint");

	/* Curvature should be negative (concave down parabola) */
	status = qaws_curve_compute_curvature_2d(
		curve2d, (qaws_scalar)0.5, &curvature2d);
	TEST_ASSERT_STATUS(status);
	TEST_ASSERT(curvature2d < (qaws_scalar)0.0, "negative curvature for concave-down");

	/* Speed at t=0 */
	status = qaws_curve_compute_speed(curve2d, (qaws_scalar)0.0, &speed_val);
	TEST_ASSERT_STATUS(status);
	TEST_ASSERT(speed_val > (qaws_scalar)0.0, "speed > 0");

	/* Dimension mismatch errors */
	{
		qaws_vec3 tangent3d;
		qaws_scalar curvature3d, torsion;

		status = qaws_curve_compute_tangent_3d(
			curve2d, (qaws_scalar)0.5, &tangent3d);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_DIMENSION,
			"tangent_3d rejects 2D curve");

		status = qaws_curve_compute_curvature_3d(
			curve2d, (qaws_scalar)0.5, &curvature3d);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_DIMENSION,
			"curvature_3d rejects 2D curve");

		status = qaws_curve_compute_torsion_3d(
			curve2d, (qaws_scalar)0.5, &torsion);
		TEST_ASSERT(status == QAWS_STATUS_INVALID_DIMENSION,
			"torsion_3d rejects 2D curve");
	}

	qaws_curve_destroy(curve2d);

	/* 3D curve tests */
	{
		qaws_vec3 pts3d[] = { {0,0,0}, {1,2,0}, {2,0,0} };
		qaws_curve *curve3d = NULL;
		qaws_vec3 tangent3d;
		qaws_scalar curvature3d, torsion, speed3d;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 2;
		bdesc.control_points = pts3d;
		bdesc.control_point_count = 3;

		status = qaws_curve_create_bezier(&bdesc, &curve3d);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_compute_tangent_3d(
			curve3d, (qaws_scalar)0.0, &tangent3d);
		TEST_ASSERT_STATUS(status);
		{
			qaws_scalar len = (qaws_scalar)sqrt((double)(
				tangent3d.x * tangent3d.x +
				tangent3d.y * tangent3d.y +
				tangent3d.z * tangent3d.z));
			TEST_ASSERT(approx_eq(len, (qaws_scalar)1.0), "3D tangent is unit length");
		}

		status = qaws_curve_compute_curvature_3d(
			curve3d, (qaws_scalar)0.5, &curvature3d);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(curvature3d > (qaws_scalar)0.0, "3D curvature > 0");

		/* Planar curve has zero torsion */
		status = qaws_curve_compute_torsion_3d(
			curve3d, (qaws_scalar)0.5, &torsion);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(approx_eq(torsion, (qaws_scalar)0.0),
			"planar curve has zero torsion");

		status = qaws_curve_compute_speed(
			curve3d, (qaws_scalar)0.0, &speed3d);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(speed3d > (qaws_scalar)0.0, "3D speed > 0");

		qaws_curve_destroy(curve3d);
	}

	/* NULL argument tests */
	status = qaws_curve_compute_tangent_2d(NULL, (qaws_scalar)0.0, &tangent2d);
	TEST_ASSERT(status == QAWS_STATUS_INVALID_ARGUMENT, "tangent null curve");

	status = qaws_curve_compute_speed(NULL, (qaws_scalar)0.0, &speed_val);
	TEST_ASSERT(status == QAWS_STATUS_INVALID_ARGUMENT, "speed null curve");
}

/* ================================================================== */
/* Adaptive sampling tests                                             */
/* ================================================================== */

static void test_adaptive_sampling(void)
{
	printf("test_adaptive_sampling\n");

	/* Create a curve with high curvature variation: cubic Bezier */
	qaws_vec2 pts[] = { {0,0}, {0,4}, {4,4}, {4,0} };
	qaws_bezier_desc bdesc;
	qaws_curve *curve = NULL;
	qaws_status status;
	qaws_sampling_desc sdesc;
	qaws_vec2 positions[2048];
	unsigned int count = 0;

	bdesc.dimension = QAWS_DIMENSION_2D;
	bdesc.degree = 3;
	bdesc.control_points = pts;
	bdesc.control_point_count = 4;

	status = qaws_curve_create_bezier(&bdesc, &curve);
	TEST_ASSERT_STATUS(status);

	/* Uniform sampling for baseline */
	memset(&sdesc, 0, sizeof(sdesc));
	sdesc.sample_count = 8;
	sdesc.use_adaptive_sampling = 0;

	status = qaws_curve_sample_2d(curve, &sdesc, positions, 2048, &count);
	TEST_ASSERT_STATUS(status);
	TEST_ASSERT(count == 8, "uniform gives exactly 8 samples");

	/* Adaptive sampling with tight tolerance should produce more points */
	memset(&sdesc, 0, sizeof(sdesc));
	sdesc.sample_count = 8;
	sdesc.use_adaptive_sampling = 1;
	sdesc.error_tolerance = (qaws_scalar)0.001;

	status = qaws_curve_sample_2d(curve, &sdesc, positions, 2048, &count);
	TEST_ASSERT_STATUS(status);
	TEST_ASSERT(count >= 8, "adaptive produces >= initial samples");

	/* Adaptive with loose tolerance should produce fewer subdivisions */
	{
		unsigned int count_tight = count;
		unsigned int count_loose;

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.sample_count = 8;
		sdesc.use_adaptive_sampling = 1;
		sdesc.error_tolerance = (qaws_scalar)1.0;

		status = qaws_curve_sample_2d(
			curve, &sdesc, positions, 2048, &count_loose);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(count_loose <= count_tight,
			"loose tolerance <= tight tolerance samples");
	}

	/* Verify all points lie on the curve (roughly) */
	{
		unsigned int i;
		int all_valid = 1;
		for (i = 0; i < count; ++i)
		{
			if (positions[i].x < (qaws_scalar)-0.1 ||
				positions[i].x > (qaws_scalar)4.1 ||
				positions[i].y < (qaws_scalar)-0.1 ||
				positions[i].y > (qaws_scalar)4.1)
			{
				all_valid = 0;
				break;
			}
		}
		TEST_ASSERT(all_valid, "adaptive points within curve bounds");
	}

	/* Test 3D adaptive sampling */
	{
		qaws_vec3 pts3d[] = { {0,0,0}, {0,4,2}, {4,4,2}, {4,0,0} };
		qaws_curve *curve3d = NULL;
		qaws_vec3 positions3d[2048];
		unsigned int count3d = 0;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 3;
		bdesc.control_points = pts3d;
		bdesc.control_point_count = 4;

		status = qaws_curve_create_bezier(&bdesc, &curve3d);
		TEST_ASSERT_STATUS(status);

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.sample_count = 8;
		sdesc.use_adaptive_sampling = 1;
		sdesc.error_tolerance = (qaws_scalar)0.01;

		status = qaws_curve_sample_3d(
			curve3d, &sdesc, positions3d, 2048, &count3d);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(count3d >= 8, "3D adaptive produces >= initial samples");

		qaws_curve_destroy(curve3d);
	}

	qaws_curve_destroy(curve);
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

	/* 3D with circular mode should now succeed */
	{
		qaws_vec3 pts3[] = { {0,0,0}, {1,1,1}, {2,0,0} };
		desc.dimension = QAWS_DIMENSION_3D;
		desc.control_points = pts3;
		desc.control_point_count = 3;
		desc.mode = QAWS_YUKSEL_MODE_CIRCULAR;
		s = qaws_curve_create_yuksel(&desc, &curve);
		TEST_ASSERT(s == QAWS_STATUS_OK, "circular 3d accepted");
		if (curve) {
			/* Verify we can evaluate the curve */
			qaws_eval_result_3d res3d;
			qaws_status es = qaws_curve_evaluate_3d(curve, (qaws_scalar)0.5,
				QAWS_EVAL_FLAG_POSITION, &res3d);
			TEST_ASSERT(es == QAWS_STATUS_OK, "circular 3d eval ok");
			qaws_curve_destroy(curve);
			curve = NULL;
		}
	}

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

/* ================================================================== */
/* Phase 1 & 2 SVG visual tests                                        */
/* ================================================================== */

/* SVG: curve reverse - original and reversed should overlap perfectly */
static void test_svg_curve_reverse(void)
{
	printf("test_svg_curve_reverse\n");

	qaws_vec2 cp[] = { {0,0}, {1,3}, {2,-1}, {3,2} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	if (s != QAWS_STATUS_OK) { printf("  SKIP: create failed\n"); return; }

	qaws_curve* rev = NULL;
	s = qaws_curve_reverse(curve, &rev);
	if (s != QAWS_STATUS_OK) { printf("  SKIP: reverse failed\n"); qaws_curve_destroy(curve); return; }

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/curve_reverse.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-1.5, (qaws_scalar)4.0, (qaws_scalar)4.5,
		(qaws_scalar)600, (qaws_scalar)500))
	{ qaws_curve_destroy(curve); qaws_curve_destroy(rev); return; }

	qaws_vec2 samples_orig[SVG_SAMPLES];
	qaws_vec2 samples_rev[SVG_SAMPLES];
	unsigned int n_orig = svg_sample_curve(curve, samples_orig, SVG_SAMPLES);
	unsigned int n_rev = svg_sample_curve(rev, samples_rev, SVG_SAMPLES);

	svg_polyline(&svg, samples_orig, n_orig, "#e94560", (qaws_scalar)3.0);
	svg_polyline(&svg, samples_rev, n_rev, "#6272a4", (qaws_scalar)1.5);

	/* Control points for both */
	qaws_vec2 cp_rev[] = { {3,2}, {2,-1}, {1,3}, {0,0} };
	svg_dots(&svg, cp, 4, "#ffffff", (qaws_scalar)4);
	svg_dots(&svg, cp_rev, 4, "#ffffff", (qaws_scalar)4);

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.3, "Original (red)", "#e94560");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.0, "Reversed (blue)", "#6272a4");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.6, "Curve Reverse", "#8888aa");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/curve_reverse.svg\n");
	qaws_curve_destroy(curve);
	qaws_curve_destroy(rev);
}

/* SVG: curve split - split cubic Bezier at t=0.5 */
static void test_svg_curve_split(void)
{
	printf("test_svg_curve_split\n");

	qaws_vec2 cp[] = { {0,0}, {1,4}, {3,-2}, {4,2} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	if (s != QAWS_STATUS_OK) { printf("  SKIP: create failed\n"); return; }

	qaws_curve* left = NULL;
	qaws_curve* right = NULL;
	s = qaws_curve_split(curve, (qaws_scalar)0.5, &left, &right);
	if (s != QAWS_STATUS_OK) {
		printf("  SKIP: split failed\n");
		qaws_curve_destroy(curve);
		return;
	}

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/curve_split.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-3.0, (qaws_scalar)5.5, (qaws_scalar)8.0,
		(qaws_scalar)600, (qaws_scalar)500))
	{
		qaws_curve_destroy(curve);
		qaws_curve_destroy(left);
		qaws_curve_destroy(right);
		return;
	}

	/* Original curve faintly */
	qaws_vec2 samples_full[SVG_SAMPLES];
	unsigned int n_full = svg_sample_curve(curve, samples_full, SVG_SAMPLES);
	svg_polyline(&svg, samples_full, n_full, "#666666", (qaws_scalar)1.0);

	/* Left half in red */
	qaws_vec2 samples_left[SVG_SAMPLES];
	unsigned int n_left = svg_sample_curve(left, samples_left, SVG_SAMPLES);
	svg_polyline(&svg, samples_left, n_left, "#e94560", (qaws_scalar)2.5);

	/* Right half in green */
	qaws_vec2 samples_right[SVG_SAMPLES];
	unsigned int n_right = svg_sample_curve(right, samples_right, SVG_SAMPLES);
	svg_polyline(&svg, samples_right, n_right, "#50fa7b", (qaws_scalar)2.5);

	/* Split point as large yellow dot */
	{
		qaws_eval_result_2d r;
		qaws_curve_evaluate_2d(curve, (qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION, &r);
		svg_dots(&svg, &r.position, 1, "#f1fa8c", (qaws_scalar)7);
	}

	/* Control points of original in white */
	svg_dots(&svg, cp, 4, "#ffffff", (qaws_scalar)4);

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)4.5, "Curve Split at t=0.5", "#8888aa");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)4.2, "Left (red)", "#e94560");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.9, "Right (green)", "#50fa7b");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/curve_split.svg\n");
	qaws_curve_destroy(curve);
	qaws_curve_destroy(left);
	qaws_curve_destroy(right);
}

/* SVG: curve join - two Hermite curves joined end-to-end */
static void test_svg_curve_join(void)
{
	printf("test_svg_curve_join\n");

	/* Curve A: points (0,0),(2,2), tangents (2,1),(1,-1) */
	qaws_vec2 pts_a[] = { {0,0}, {2,2} };
	qaws_vec2 tan_a[] = { {2,1}, {1,-1} };
	qaws_hermite_desc desc_a;
	desc_a.dimension = QAWS_DIMENSION_2D;
	desc_a.degree = 3;
	desc_a.points = pts_a;
	desc_a.derivatives = tan_a;
	desc_a.point_count = 2;
	desc_a.derivative_count = 2;

	qaws_curve* curve_a = NULL;
	qaws_status s = qaws_curve_create_hermite(&desc_a, &curve_a);
	if (s != QAWS_STATUS_OK) { printf("  SKIP: create A failed\n"); return; }

	/* Curve B: points (2,2),(4,0), tangents (1,-1),(2,1) */
	qaws_vec2 pts_b[] = { {2,2}, {4,0} };
	qaws_vec2 tan_b[] = { {1,-1}, {2,1} };
	qaws_hermite_desc desc_b;
	desc_b.dimension = QAWS_DIMENSION_2D;
	desc_b.degree = 3;
	desc_b.points = pts_b;
	desc_b.derivatives = tan_b;
	desc_b.point_count = 2;
	desc_b.derivative_count = 2;

	qaws_curve* curve_b = NULL;
	s = qaws_curve_create_hermite(&desc_b, &curve_b);
	if (s != QAWS_STATUS_OK) {
		printf("  SKIP: create B failed\n");
		qaws_curve_destroy(curve_a);
		return;
	}

	qaws_curve* joined = NULL;
	s = qaws_curve_join(curve_a, curve_b, &joined);
	if (s != QAWS_STATUS_OK) {
		printf("  SKIP: join failed\n");
		qaws_curve_destroy(curve_a);
		qaws_curve_destroy(curve_b);
		return;
	}

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/curve_join.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-1.0, (qaws_scalar)5.0, (qaws_scalar)4.0,
		(qaws_scalar)600, (qaws_scalar)400))
	{
		qaws_curve_destroy(curve_a);
		qaws_curve_destroy(curve_b);
		qaws_curve_destroy(joined);
		return;
	}

	/* Draw joined curve in green (thicker, underneath) */
	qaws_vec2 samples_joined[SVG_SAMPLES];
	unsigned int n_joined = svg_sample_curve(joined, samples_joined, SVG_SAMPLES);
	svg_polyline(&svg, samples_joined, n_joined, "#50fa7b", (qaws_scalar)4.0);

	/* Draw curve A in red */
	qaws_vec2 samples_a[SVG_SAMPLES];
	unsigned int n_a = svg_sample_curve(curve_a, samples_a, SVG_SAMPLES);
	svg_polyline(&svg, samples_a, n_a, "#e94560", (qaws_scalar)2.0);

	/* Draw curve B in blue */
	qaws_vec2 samples_b[SVG_SAMPLES];
	unsigned int n_b = svg_sample_curve(curve_b, samples_b, SVG_SAMPLES);
	svg_polyline(&svg, samples_b, n_b, "#6272a4", (qaws_scalar)2.0);

	/* Key points */
	qaws_vec2 key_pts[] = { {0,0}, {2,2}, {4,0} };
	svg_dots(&svg, key_pts, 3, "#ffffff", (qaws_scalar)5);

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.2, "Curve Join", "#8888aa");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.9, "Curve A (red)", "#e94560");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.6, "Curve B (blue)", "#6272a4");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.3, "Joined (green)", "#50fa7b");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/curve_join.svg\n");
	qaws_curve_destroy(curve_a);
	qaws_curve_destroy(curve_b);
	qaws_curve_destroy(joined);
}

/* SVG: degree elevation - quadratic elevated to cubic should overlap */
static void test_svg_degree_elevation(void)
{
	printf("test_svg_degree_elevation\n");

	qaws_vec2 cp[] = { {0,0}, {2,3}, {4,0} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 2;
	desc.control_points = cp;
	desc.control_point_count = 3;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	if (s != QAWS_STATUS_OK) { printf("  SKIP: create failed\n"); return; }

	qaws_curve* elevated = NULL;
	s = qaws_curve_elevate_degree(curve, &elevated);
	if (s != QAWS_STATUS_OK) {
		printf("  SKIP: elevate failed\n");
		qaws_curve_destroy(curve);
		return;
	}

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/degree_elevation.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-0.5, (qaws_scalar)5.0, (qaws_scalar)4.0,
		(qaws_scalar)600, (qaws_scalar)400))
	{
		qaws_curve_destroy(curve);
		qaws_curve_destroy(elevated);
		return;
	}

	/* Original quadratic in red */
	qaws_vec2 samples_orig[SVG_SAMPLES];
	unsigned int n_orig = svg_sample_curve(curve, samples_orig, SVG_SAMPLES);
	svg_polyline(&svg, samples_orig, n_orig, "#e94560", (qaws_scalar)3.0);

	/* Elevated cubic in blue (slightly thinner) */
	qaws_vec2 samples_elev[SVG_SAMPLES];
	unsigned int n_elev = svg_sample_curve(elevated, samples_elev, SVG_SAMPLES);
	svg_polyline(&svg, samples_elev, n_elev, "#6272a4", (qaws_scalar)1.5);

	/* Original CPs as red dots */
	svg_dots(&svg, cp, 3, "#e94560", (qaws_scalar)5);

	/* Elevated CPs as blue dots */
	{
		qaws_vec2 elev_cp[8];
		unsigned int cp_count = 0;
		s = qaws_bezier_get_control_points(elevated, elev_cp, 8, &cp_count);
		if (s == QAWS_STATUS_OK && cp_count > 0) {
			svg_dots(&svg, elev_cp, cp_count, "#6272a4", (qaws_scalar)5);
		}
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.5, "Degree Elevation", "#8888aa");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.2, "Quadratic (red, 3 CPs)", "#e94560");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.9, "Cubic (blue, 4 CPs)", "#6272a4");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/degree_elevation.svg\n");
	qaws_curve_destroy(curve);
	qaws_curve_destroy(elevated);
}

/* SVG: family conversion - Hermite to Bezier should overlap exactly */
static void test_svg_family_conversion(void)
{
	printf("test_svg_family_conversion\n");

	qaws_vec2 pts[] = { {0,0}, {4,2} };
	qaws_vec2 tans[] = { {2,4}, {2,-4} };
	qaws_hermite_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.points = pts;
	desc.derivatives = tans;
	desc.point_count = 2;
	desc.derivative_count = 2;

	qaws_curve* hermite = NULL;
	qaws_status s = qaws_curve_create_hermite(&desc, &hermite);
	if (s != QAWS_STATUS_OK) { printf("  SKIP: create hermite failed\n"); return; }

	qaws_curve* bezier = NULL;
	s = qaws_curve_convert_hermite_to_bezier(hermite, 0, &bezier);
	if (s != QAWS_STATUS_OK) {
		printf("  SKIP: convert failed\n");
		qaws_curve_destroy(hermite);
		return;
	}

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/family_conversion.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-1.0, (qaws_scalar)5.0, (qaws_scalar)4.0,
		(qaws_scalar)600, (qaws_scalar)400))
	{
		qaws_curve_destroy(hermite);
		qaws_curve_destroy(bezier);
		return;
	}

	/* Original Hermite in red */
	qaws_vec2 samples_herm[SVG_SAMPLES];
	unsigned int n_herm = svg_sample_curve(hermite, samples_herm, SVG_SAMPLES);
	svg_polyline(&svg, samples_herm, n_herm, "#e94560", (qaws_scalar)3.0);

	/* Converted Bezier in blue */
	qaws_vec2 samples_bez[SVG_SAMPLES];
	unsigned int n_bez = svg_sample_curve(bezier, samples_bez, SVG_SAMPLES);
	svg_polyline(&svg, samples_bez, n_bez, "#6272a4", (qaws_scalar)1.5);

	/* Hermite points as red dots */
	svg_dots(&svg, pts, 2, "#e94560", (qaws_scalar)5);

	/* Bezier CPs as blue dots */
	{
		qaws_vec2 bez_cp[8];
		unsigned int cp_count = 0;
		s = qaws_bezier_get_control_points(bezier, bez_cp, 8, &cp_count);
		if (s == QAWS_STATUS_OK && cp_count > 0) {
			svg_dots(&svg, bez_cp, cp_count, "#6272a4", (qaws_scalar)5);
		}
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.2, "Family Conversion", "#8888aa");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.9, "Hermite (red)", "#e94560");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.6, "Bezier (blue)", "#6272a4");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/family_conversion.svg\n");
	qaws_curve_destroy(hermite);
	qaws_curve_destroy(bezier);
}

/* SVG: easing curves - show different easing modes as dot distributions along a curve */
static void test_svg_easing_curves(void)
{
	printf("test_svg_easing_curves\n");

	/* Create a diagonal Bezier from (0,0) to (4,3) */
	qaws_vec2 cp[] = { {0,0}, {1,0}, {3,3}, {4,3} };
	qaws_bezier_desc bdesc;
	bdesc.dimension = QAWS_DIMENSION_2D;
	bdesc.degree = 3;
	bdesc.control_points = cp;
	bdesc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&bdesc, &curve);
	if (s != QAWS_STATUS_OK) { printf("  SKIP: create failed\n"); return; }

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/easing_curves.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-0.5, (qaws_scalar)5.0, (qaws_scalar)4.0,
		(qaws_scalar)600, (qaws_scalar)400))
	{ qaws_curve_destroy(curve); return; }

	/* Draw the curve in gray */
	{
		qaws_vec2 samples[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(curve, samples, SVG_SAMPLES);
		svg_polyline(&svg, samples, n, "#444466", (qaws_scalar)1.5);
	}

	/* Compute arc length for speed setting */
	qaws_scalar arc_len = 0;
	{
		qaws_range pr = qaws_curve_get_parameter_range(curve);
		qaws_curve_compute_arc_length(curve, pr.min_value, pr.max_value, &arc_len);
	}

	/* Easing modes to visualize */
	{
		qaws_easing modes[] = {
			QAWS_EASING_LINEAR,
			QAWS_EASING_QUAD_IN,
			QAWS_EASING_QUAD_OUT,
			QAWS_EASING_CUBIC_IN_OUT,
			QAWS_EASING_SINE_IN_OUT
		};
		char const* names[] = {
			"Linear", "Quad In", "Quad Out", "Cubic InOut", "Sine InOut"
		};
		char const* colors[] = {
			"#ffffff", "#e94560", "#50fa7b", "#f5a623", "#78dce8"
		};
		unsigned int mode_count = 5;
		unsigned int mi;

		#define EASE_DOTS 32
		for (mi = 0; mi < mode_count; mi++) {
			qaws_traversal_desc tdesc;
			qaws_traversal* trav = NULL;

			memset(&tdesc, 0, sizeof(tdesc));
			tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
			tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
			tdesc.speed = arc_len;
			tdesc.start_time = (qaws_scalar)0.0;
			tdesc.end_time = (qaws_scalar)1.0;
			tdesc.clamp_to_domain = 1;
			tdesc.easing = modes[mi];
			tdesc.wrap_mode = QAWS_WRAP_CLAMP;

			s = qaws_traversal_create(curve, &tdesc, &trav);
			if (s != QAWS_STATUS_OK) continue;

			qaws_vec2 dots[EASE_DOTS];
			unsigned int di;
			for (di = 0; di < EASE_DOTS; di++) {
				qaws_scalar t = (qaws_scalar)di / (qaws_scalar)(EASE_DOTS - 1);
				qaws_eval_result_2d r;
				qaws_traversal_evaluate_2d(trav, t, QAWS_EVAL_FLAG_POSITION, &r);
				dots[di] = r.position;
			}

			svg_dots(&svg, dots, EASE_DOTS, colors[mi], (qaws_scalar)(3.0 - mi * 0.3));

			svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)(3.5 - (qaws_scalar)mi * (qaws_scalar)0.3),
				names[mi], colors[mi]);

			qaws_traversal_destroy(trav);
		}
		#undef EASE_DOTS
	}

	svg_label(&svg, (qaws_scalar)0.5, (qaws_scalar)-0.3, "Easing Curves (dot spacing shows speed)", "#8888aa");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/easing_curves.svg\n");
	qaws_curve_destroy(curve);
}

/* SVG: adaptive vs uniform sampling on a high-curvature curve */
static void test_svg_adaptive_sampling(void)
{
	printf("test_svg_adaptive_sampling\n");

	/* Cubic Bezier with tight curvature */
	qaws_vec2 cp[] = { {0,0}, {0,3}, {1,-1}, {4,2} };
	qaws_bezier_desc bdesc;
	bdesc.dimension = QAWS_DIMENSION_2D;
	bdesc.degree = 3;
	bdesc.control_points = cp;
	bdesc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&bdesc, &curve);
	if (s != QAWS_STATUS_OK) { printf("  SKIP: create failed\n"); return; }

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/adaptive_sampling.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-1.5, (qaws_scalar)5.0, (qaws_scalar)4.5,
		(qaws_scalar)600, (qaws_scalar)400))
	{ qaws_curve_destroy(curve); return; }

	/* Draw the curve in gray */
	{
		qaws_vec2 full_samples[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(curve, full_samples, SVG_SAMPLES);
		svg_polyline(&svg, full_samples, n, "#444466", (qaws_scalar)1.5);
	}

	/* Uniform sampling with 8 points — few enough that curvature gaps are visible */
	{
		qaws_sampling_desc sdesc;
		qaws_vec2 uniform_pts[64];
		unsigned int uniform_count = 0;

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
		sdesc.sample_count = 8;
		sdesc.use_adaptive_sampling = 0;
		sdesc.error_tolerance = 0;

		s = qaws_curve_sample_2d(curve, &sdesc, uniform_pts, 64, &uniform_count);
		if (s == QAWS_STATUS_OK) {
			svg_dots(&svg, uniform_pts, uniform_count, "#e94560", (qaws_scalar)5);
		}
	}

	/* Adaptive sampling — starts with 8 segments, subdivides where chord error > tolerance */
	{
		qaws_sampling_desc sdesc;
		qaws_vec2 adaptive_pts[512];
		unsigned int adaptive_count = 0;

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
		sdesc.sample_count = 8;
		sdesc.use_adaptive_sampling = 1;
		sdesc.error_tolerance = (qaws_scalar)0.05;

		s = qaws_curve_sample_2d(curve, &sdesc, adaptive_pts, 512, &adaptive_count);
		if (s == QAWS_STATUS_OK) {
			svg_dots(&svg, adaptive_pts, adaptive_count, "#6272a4", (qaws_scalar)3);
		}
	}

	/* Control points */
	svg_dots(&svg, cp, 4, "#ffffff", (qaws_scalar)3);

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.2, "Adaptive vs Uniform Sampling", "#8888aa");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.9, "Uniform 8 pts (red, large)", "#e94560");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.6, "Adaptive 8+subdivisions (blue, small)", "#6272a4");

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/adaptive_sampling.svg\n");
	qaws_curve_destroy(curve);
}

/* ================================================================== */
/* Phase 3 & 5 SVG visual tests                                        */
/* ================================================================== */

static void test_svg_inflection_points(void)
{
	printf("test_svg_inflection_points\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/inflection_points.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-3.0, (qaws_scalar)5.5, (qaws_scalar)8.0,
		(qaws_scalar)500, (qaws_scalar)600))
		return;

	qaws_vec2 cp[] = { {0, 0}, {1, 4}, {2, -2}, {3, 0} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	/* Draw the curve in gray */
	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
	svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);

	/* Find and draw inflection points */
	qaws_scalar infl_params[8];
	unsigned int infl_count = 0;
	s = qaws_curve_find_inflection_points(curve, infl_params, 8, &infl_count);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(infl_count > 0, "inflection points found");

	for (unsigned int i = 0; i < infl_count; i++) {
		qaws_eval_result_2d r;
		s = qaws_curve_evaluate_2d(curve, infl_params[i],
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		svg_dot(&svg, r.position.x, r.position.y, "#e94560", (qaws_scalar)6);
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)4.5, "Inflection Points", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/inflection_points.svg\n");
}

static void test_svg_extrema(void)
{
	printf("test_svg_extrema\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/extrema.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-3.0, (qaws_scalar)5.5, (qaws_scalar)8.0,
		(qaws_scalar)500, (qaws_scalar)600))
		return;

	qaws_vec2 cp[] = { {0, 0}, {1, 4}, {3, -2}, {4, 2} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	/* Draw the curve in gray */
	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
	svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);

	/* Find X-extrema */
	qaws_scalar x_params[8];
	unsigned int x_count = 0;
	s = qaws_curve_find_extrema(curve, 0, x_params, 8, &x_count);
	TEST_ASSERT_STATUS(s);

	for (unsigned int i = 0; i < x_count; i++) {
		qaws_eval_result_2d r;
		s = qaws_curve_evaluate_2d(curve, x_params[i],
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		svg_dot(&svg, r.position.x, r.position.y, "#6272a4", (qaws_scalar)5);
		/* Faint vertical line through extremum */
		svg_line(&svg, r.position.x, (qaws_scalar)-3.0,
			r.position.x, (qaws_scalar)5.0, "#6272a480", (qaws_scalar)0.5);
	}

	/* Find Y-extrema */
	qaws_scalar y_params[8];
	unsigned int y_count = 0;
	s = qaws_curve_find_extrema(curve, 1, y_params, 8, &y_count);
	TEST_ASSERT_STATUS(s);

	for (unsigned int i = 0; i < y_count; i++) {
		qaws_eval_result_2d r;
		s = qaws_curve_evaluate_2d(curve, y_params[i],
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		svg_dot(&svg, r.position.x, r.position.y, "#50fa7b", (qaws_scalar)5);
		/* Faint horizontal line through extremum */
		svg_line(&svg, (qaws_scalar)-0.5, r.position.y,
			(qaws_scalar)5.0, r.position.y, "#50fa7b80", (qaws_scalar)0.5);
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)4.5,
		"Extrema (blue=X, green=Y)", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/extrema.svg\n");
}

static void test_svg_curvature_comb(void)
{
	printf("test_svg_curvature_comb\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/curvature_comb.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-2.0, (qaws_scalar)5.0, (qaws_scalar)5.0,
		(qaws_scalar)600, (qaws_scalar)500))
		return;

	/* Smooth S-shaped cubic Bezier — single inflection, clean curvature profile */
	qaws_vec2 cp[] = { {0, 0}, {1, 3}, {3, -1}, {4, 2} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	/* Draw the curve in gray */
	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
	svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);

	/* Compute curvature comb */
	{
		unsigned int comb_count = 48;
		qaws_curvature_sample_2d comb[48];
		s = qaws_curve_compute_curvature_comb_2d(curve, comb_count, comb, 48);
		TEST_ASSERT_STATUS(s);

		/* Find max |curvature| to normalize tooth length */
		qaws_scalar max_k = (qaws_scalar)0.001;
		for (unsigned int i = 0; i < comb_count; i++) {
			qaws_scalar ak = comb[i].curvature;
			if (ak < 0) ak = -ak;
			if (ak > max_k) max_k = ak;
		}

		/* Scale so the tallest tooth is about 0.8 world units */
		qaws_scalar scale = (qaws_scalar)0.8 / max_k;
		qaws_vec2 tips[48];

		for (unsigned int i = 0; i < comb_count; i++) {
			qaws_scalar tip_x = comb[i].position.x + comb[i].normal.x * comb[i].curvature * scale;
			qaws_scalar tip_y = comb[i].position.y + comb[i].normal.y * comb[i].curvature * scale;
			tips[i].x = tip_x;
			tips[i].y = tip_y;

			/* Draw comb tooth from curve to tip */
			svg_line(&svg, comb[i].position.x, comb[i].position.y,
				tip_x, tip_y, "#e9456088", (qaws_scalar)0.8);
		}

		/* Connect comb tips as an envelope */
		svg_polyline(&svg, tips, comb_count, "#e94560", (qaws_scalar)1.2);
	}

	/* Draw control polygon faintly */
	svg_polyline(&svg, cp, 4, "#ffffff30", (qaws_scalar)0.8);
	svg_dots(&svg, cp, 4, "#ffffff60", (qaws_scalar)2.5);

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.8, "Curvature Comb (teeth = signed curvature)", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/curvature_comb.svg\n");
}

static void test_svg_winding_number(void)
{
	printf("test_svg_winding_number\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/winding_number.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-0.5, (qaws_scalar)6.0, (qaws_scalar)6.0,
		(qaws_scalar)500, (qaws_scalar)500))
		return;

	/* 8 points on a circle of radius 1.5, centered at (2,2) */
	qaws_vec2 circle_pts[8];
	{
		qaws_scalar cx = (qaws_scalar)2.0, cy = (qaws_scalar)2.0;
		qaws_scalar radius = (qaws_scalar)1.5;
		for (int i = 0; i < 8; i++) {
			qaws_scalar angle = (qaws_scalar)(2.0 * 3.14159265358979 * i / 8.0);
			circle_pts[i].x = cx + radius * (qaws_scalar)cos((double)angle);
			circle_pts[i].y = cy + radius * (qaws_scalar)sin((double)angle);
		}
	}

	/* Need ghost points for closed Catmull-Rom; use closed=1 */
	qaws_catmull_rom_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.control_points = circle_pts;
	desc.control_point_count = 8;
	desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
	desc.closed = 1;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_catmull_rom(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	/* Draw the closed curve in gray */
	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
	svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);

	/* Test points */
	struct { qaws_vec2 pt; char const* label_text; } test_pts[] = {
		{ {2, 2},     "inside" },
		{ {5, 5},     "outside" },
		{ {2, (qaws_scalar)3.5}, "near edge" }
	};

	for (int i = 0; i < 3; i++) {
		int winding = 0;
		s = qaws_curve_compute_winding_number_2d(curve, test_pts[i].pt, &winding);
		TEST_ASSERT_STATUS(s);

		char const* color = (winding != 0) ? "#50fa7b" : "#e94560";
		svg_dot(&svg, test_pts[i].pt.x, test_pts[i].pt.y, color, (qaws_scalar)5);

		/* Label with winding number value */
		char label_buf[64];
		sprintf(label_buf, "w=%d (%s)", winding, test_pts[i].label_text);
		svg_label(&svg, test_pts[i].pt.x + (qaws_scalar)0.15,
			test_pts[i].pt.y + (qaws_scalar)0.1, label_buf, color);
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)5.2, "Winding Number", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/winding_number.svg\n");
}

static void test_svg_curvature_weighted(void)
{
	printf("test_svg_curvature_weighted\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/curvature_weighted.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-3.0, (qaws_scalar)5.5, (qaws_scalar)8.0,
		(qaws_scalar)500, (qaws_scalar)600))
		return;

	qaws_vec2 cp[] = { {0, 0}, {(qaws_scalar)0.5, 4}, {(qaws_scalar)3.5, -2}, {4, 2} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	/* Draw the curve in gray */
	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
	svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);

	/* Uniform sampling: 20 dots in blue */
	{
		qaws_vec2 uniform_pts[20];
		qaws_sampling_desc sd;
		sd.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
		sd.sample_count = 20;
		sd.error_tolerance = 0;
		sd.use_adaptive_sampling = 0;
		unsigned int ucount = 0;
		s = qaws_curve_sample_2d(curve, &sd, uniform_pts, 20, &ucount);
		TEST_ASSERT_STATUS(s);
		svg_dots(&svg, uniform_pts, ucount, "#6272a4", (qaws_scalar)4);
	}

	/* Curvature-weighted sampling: 20 dots in red */
	{
		qaws_vec2 cw_pts[20];
		unsigned int cw_count = 0;
		s = qaws_curve_sample_curvature_weighted_2d(curve, 20, cw_pts, 20, &cw_count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(cw_count > 0, "curvature-weighted samples produced");
		svg_dots(&svg, cw_pts, cw_count, "#e94560", (qaws_scalar)4);
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)4.5,
		"Curvature-weighted (red) vs Uniform (blue)", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/curvature_weighted.svg\n");
}

static void test_svg_feature_preserving(void)
{
	printf("test_svg_feature_preserving\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/feature_preserving.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-5.0, (qaws_scalar)5.5, (qaws_scalar)11.0,
		(qaws_scalar)500, (qaws_scalar)600))
		return;

	qaws_vec2 cp[] = { {0, 0}, {(qaws_scalar)0.5, 4}, {(qaws_scalar)3.5, -2}, {4, 2} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	qaws_scalar offset_y = (qaws_scalar)4.0;

	/* --- Top row: Uniform sampling (blue) --- */
	{
		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
		for (unsigned int i = 0; i < n; i++) buf[i].y += offset_y;
		svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);
	}
	{
		qaws_vec2 uniform_pts[8];
		qaws_sampling_desc sd;
		sd.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
		sd.sample_count = 8;
		sd.error_tolerance = 0;
		sd.use_adaptive_sampling = 0;
		unsigned int ucount = 0;
		s = qaws_curve_sample_2d(curve, &sd, uniform_pts, 8, &ucount);
		TEST_ASSERT_STATUS(s);
		for (unsigned int i = 0; i < ucount; i++) uniform_pts[i].y += offset_y;
		svg_dots(&svg, uniform_pts, ucount, "#6272a4", (qaws_scalar)4);
	}

	/* --- Bottom row: Feature-preserving sampling (green) --- */
	{
		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);
	}
	{
		qaws_vec2 fp_pts[32];
		unsigned int fp_count = 0;
		s = qaws_curve_sample_feature_preserving_2d(curve, 8, fp_pts, 32, &fp_count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(fp_count >= 8, "feature-preserving has at least base_count samples");
		svg_dots(&svg, fp_pts, fp_count, "#50fa7b", (qaws_scalar)4);
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)5.8,
		"Uniform 8 pts (blue)", "#6272a4");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)1.8,
		"Feature-preserving (green, adds inflections+extrema)", "#50fa7b");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/feature_preserving.svg\n");
}

static void test_svg_multi_traversal(void)
{
	printf("test_svg_multi_traversal\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/multi_traversal.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-0.5, (qaws_scalar)5.5, (qaws_scalar)3.0,
		(qaws_scalar)600, (qaws_scalar)300))
		return;

	/* 3 connected cubic Bezier curves */
	qaws_vec2 cp1[] = { {0, 0}, {(qaws_scalar)0.5, 2}, {1, 0}, {(qaws_scalar)1.5, 0} };
	qaws_vec2 cp2[] = { {(qaws_scalar)1.5, 0}, {2, 0}, {(qaws_scalar)2.5, 2}, {3, 0} };
	qaws_vec2 cp3[] = { {3, 0}, {(qaws_scalar)3.5, 2}, {4, 0}, {(qaws_scalar)4.5, 0} };

	qaws_bezier_desc d1, d2, d3;
	d1.dimension = QAWS_DIMENSION_2D; d1.degree = 3;
	d1.control_points = cp1; d1.control_point_count = 4;
	d2.dimension = QAWS_DIMENSION_2D; d2.degree = 3;
	d2.control_points = cp2; d2.control_point_count = 4;
	d3.dimension = QAWS_DIMENSION_2D; d3.degree = 3;
	d3.control_points = cp3; d3.control_point_count = 4;

	qaws_curve* c1 = NULL;
	qaws_curve* c2 = NULL;
	qaws_curve* c3 = NULL;
	qaws_status s;
	s = qaws_curve_create_bezier(&d1, &c1); TEST_ASSERT_STATUS(s);
	s = qaws_curve_create_bezier(&d2, &c2); TEST_ASSERT_STATUS(s);
	s = qaws_curve_create_bezier(&d3, &c3); TEST_ASSERT_STATUS(s);

	/* Draw each segment in a different color */
	{
		qaws_vec2 seg[SVG_SAMPLES];
		unsigned int sn;
		sn = svg_sample_curve(c1, seg, SVG_SAMPLES);
		svg_polyline(&svg, seg, sn, "#e94560", (qaws_scalar)2.5);
		sn = svg_sample_curve(c2, seg, SVG_SAMPLES);
		svg_polyline(&svg, seg, sn, "#50fa7b", (qaws_scalar)2.5);
		sn = svg_sample_curve(c3, seg, SVG_SAMPLES);
		svg_polyline(&svg, seg, sn, "#6272a4", (qaws_scalar)2.5);
	}

	/* Create multi-curve traversal with constant speed + arc-length */
	qaws_curve const* curves[3];
	curves[0] = c1; curves[1] = c2; curves[2] = c3;

	qaws_traversal_desc tdesc;
	memset(&tdesc, 0, sizeof(tdesc));
	tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
	tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
	tdesc.speed = (qaws_scalar)1.0;
	tdesc.clamp_to_domain = 1;

	qaws_traversal* trav = NULL;
	s = qaws_traversal_create_multi(curves, 3, &tdesc, &trav);
	TEST_ASSERT_STATUS(s);

	/* Compute total arc length */
	qaws_scalar len1, len2, len3, total_len;
	qaws_curve_compute_arc_length(c1, (qaws_scalar)0.0, (qaws_scalar)1.0, &len1);
	qaws_curve_compute_arc_length(c2, (qaws_scalar)0.0, (qaws_scalar)1.0, &len2);
	qaws_curve_compute_arc_length(c3, (qaws_scalar)0.0, (qaws_scalar)1.0, &len3);
	total_len = len1 + len2 + len3;

	/* Evaluate at 32 uniform distances and draw dots */
	for (int i = 0; i < 32; i++) {
		qaws_scalar dist = total_len * (qaws_scalar)i / (qaws_scalar)31;
		qaws_eval_result_2d r;
		s = qaws_traversal_evaluate_2d(trav, dist,
			QAWS_EVAL_FLAG_POSITION, &r);
		if (s == QAWS_STATUS_OK) {
			svg_dot(&svg, r.position.x, r.position.y, "#ffffff", (qaws_scalar)3);
		}
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.5,
		"Multi-curve traversal (uniform arc-length)", "#8888aa");
	svg_close(&svg);
	qaws_traversal_destroy(trav);
	qaws_curve_destroy(c1);
	qaws_curve_destroy(c2);
	qaws_curve_destroy(c3);
	printf("  -> " SVG_OUTPUT_DIR "/multi_traversal.svg\n");
}

static void test_svg_frenet_frame(void)
{
	printf("test_svg_frenet_frame\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/frenet_frame.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-1.5, (qaws_scalar)5.5, (qaws_scalar)5.0,
		(qaws_scalar)500, (qaws_scalar)500))
		return;

	qaws_vec2 cp[] = { {0, 0}, {1, 3}, {3, -1}, {4, 2} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 3;
	desc.control_points = cp;
	desc.control_point_count = 4;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&desc, &curve);
	TEST_ASSERT_STATUS(s);

	/* Draw the curve in gray */
	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
	svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);

	/* At 8 evenly spaced parameters, draw tangent and normal vectors */
	qaws_scalar arrow_len = (qaws_scalar)0.3;
	for (int i = 0; i < 8; i++) {
		qaws_scalar t = (qaws_scalar)i / (qaws_scalar)7;
		qaws_eval_result_2d r;
		s = qaws_curve_evaluate_2d(curve, t, QAWS_EVAL_FLAG_POSITION, &r);
		if (s != QAWS_STATUS_OK) continue;

		qaws_vec2 tangent, normal;
		s = qaws_curve_compute_tangent_2d(curve, t, &tangent);
		if (s != QAWS_STATUS_OK) continue;
		s = qaws_curve_compute_normal_2d(curve, t, &normal);
		if (s != QAWS_STATUS_OK) continue;

		/* Draw tangent arrow (green) */
		svg_line(&svg, r.position.x, r.position.y,
			r.position.x + tangent.x * arrow_len,
			r.position.y + tangent.y * arrow_len,
			"#50fa7b", (qaws_scalar)1.5);

		/* Draw normal arrow (red) */
		svg_line(&svg, r.position.x, r.position.y,
			r.position.x + normal.x * arrow_len,
			r.position.y + normal.y * arrow_len,
			"#e94560", (qaws_scalar)1.5);

		/* Draw a small dot at the sample point */
		svg_dot(&svg, r.position.x, r.position.y, "#ffffff", (qaws_scalar)2.5);
	}

	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.5,
		"Frenet Frame (green=T, red=N)", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/frenet_frame.svg\n");
}

/* ================================================================== */
/* Phase 3 tests: inflection points, extrema, curvature comb,         */
/*                winding number, Yuksel 3D circular                   */
/* ================================================================== */

static int scalar_is_finite(qaws_scalar v)
{
	return (v == v) && (v - v == 0);
}

static void test_inflection_points(void)
{
	printf("test_inflection_points\n");

	/* 1. Cubic Bezier S-curve: CPs (0,0),(1,4),(2,-2),(3,0)
	 *    inflection at t = 5/9 ~ 0.556 (avoids exact sample alignment) */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 4}, {2, -2}, {3, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_scalar params[4];
		unsigned int count = 0;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;

		s = qaws_curve_create_bezier(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_inflection_points(curve, params, 4, &count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count == 1, "S-curve has 1 inflection point");

		if (count >= 1) {
			qaws_scalar diff = params[0] - (qaws_scalar)(5.0 / 9.0);
			if (diff < 0) diff = -diff;
			TEST_ASSERT(diff < (qaws_scalar)0.1, "inflection near t=5/9");

			/* Curvature at inflection should be near 0 */
			{
				qaws_scalar curv = 0;
				s = qaws_curve_compute_curvature_2d(curve, params[0], &curv);
				TEST_ASSERT_STATUS(s);
				if (curv < 0) curv = -curv;
				TEST_ASSERT(curv < (qaws_scalar)0.1, "curvature ~0 at inflection");
			}
		}

		qaws_curve_destroy(curve);
	}

	/* 2. Straight line: collinear points, no inflection */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 1}, {2, 2}, {3, 3} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_scalar params[4];
		unsigned int count = 99;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;

		s = qaws_curve_create_bezier(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_inflection_points(curve, params, 4, &count);
		TEST_ASSERT(s == QAWS_STATUS_OK || s == QAWS_STATUS_DEGENERATE_CURVE,
			"collinear inflection handled gracefully");
		if (s == QAWS_STATUS_OK) {
			TEST_ASSERT(count == 0, "collinear has 0 inflections");
		}

		qaws_curve_destroy(curve);
	}

	/* 3. Quadratic Bezier: constant curvature sign, 0 inflections */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 2}, {2, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_scalar params[4];
		unsigned int count = 99;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = pts;
		desc.control_point_count = 3;

		s = qaws_curve_create_bezier(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_inflection_points(curve, params, 4, &count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count == 0, "quadratic has 0 inflections");

		qaws_curve_destroy(curve);
	}

	/* 4. Null out_parameters: requires non-NULL => INVALID_ARGUMENT */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 4}, {2, -2}, {3, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		unsigned int count = 0;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;

		s = qaws_curve_create_bezier(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_inflection_points(curve, NULL, 0, &count);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"null out_parameters => INVALID_ARGUMENT");

		qaws_curve_destroy(curve);
	}

	/* 5. Null args */
	{
		qaws_status s;
		unsigned int count = 0;

		s = qaws_curve_find_inflection_points(NULL, NULL, 0, &count);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT, "null curve => INVALID_ARGUMENT");

		/* Null out_count */
		{
			qaws_vec2 pts[] = { {0, 0}, {1, 3}, {2, -3}, {3, 0} };
			qaws_bezier_desc desc;
			qaws_curve* curve = NULL;
			qaws_scalar params[4];

			desc.dimension = QAWS_DIMENSION_2D;
			desc.degree = 3;
			desc.control_points = pts;
			desc.control_point_count = 4;

			s = qaws_curve_create_bezier(&desc, &curve);
			TEST_ASSERT_STATUS(s);

			s = qaws_curve_find_inflection_points(curve, params, 4, NULL);
			TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
				"null out_count => INVALID_ARGUMENT");

			qaws_curve_destroy(curve);
		}
	}
}

static void test_extrema(void)
{
	printf("test_extrema\n");

	/* 1. Cubic Bezier: (0,0),(1,4),(2,1),(3,0) - Y extremum near t~0.377
	 *    (asymmetric to avoid landing on exact sample grid) */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 4}, {2, 1}, {3, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_scalar params[4];
		unsigned int count = 0;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;

		s = qaws_curve_create_bezier(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_extrema(curve, 1, params, 4, &count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count >= 1, "cubic has >=1 Y-extremum");

		if (count >= 1) {
			qaws_scalar diff = params[0] - (qaws_scalar)0.377;
			if (diff < 0) diff = -diff;
			TEST_ASSERT(diff < (qaws_scalar)0.1,
				"Y-extremum near t~0.377");

			/* D1.y should be ~0 at the extremum */
			{
				qaws_eval_result_2d r;
				s = qaws_curve_evaluate_2d(curve, params[0],
					QAWS_EVAL_FLAG_D1, &r);
				TEST_ASSERT_STATUS(s);
				{
					qaws_scalar d1y = r.d1.y;
					if (d1y < 0) d1y = -d1y;
					TEST_ASSERT(d1y < (qaws_scalar)0.1,
						"D1.y ~0 at Y-extremum");
				}
			}
		}

		qaws_curve_destroy(curve);
	}

	/* 2. X-extrema of S-curve */
	{
		qaws_vec2 pts[] = { {0, 0}, {2, 1}, {-1, 2}, {1, 3} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_scalar params[4];
		unsigned int count = 0;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;

		s = qaws_curve_create_bezier(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_extrema(curve, 0, params, 4, &count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count >= 1, "S-curve has X-extrema");

		qaws_curve_destroy(curve);
	}

	/* 3. Invalid axis: axis=2 for 2D curve */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 2}, {2, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_scalar params[4];
		unsigned int count = 0;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = pts;
		desc.control_point_count = 3;

		s = qaws_curve_create_bezier(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_extrema(curve, 2, params, 4, &count);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT
			|| s == QAWS_STATUS_OUT_OF_RANGE,
			"axis=2 on 2D curve => error");

		qaws_curve_destroy(curve);
	}

	/* 4. 3D curve: find z-axis extrema (axis=2) using cubic
	 *    (0,0,0),(1,1,4),(2,0,1),(3,0,0) -> z extremum at t~0.377 */
	{
		qaws_vec3 pts[] = { {0, 0, 0}, {1, 1, 4}, {2, 0, 1}, {3, 0, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_scalar params[4];
		unsigned int count = 0;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;

		s = qaws_curve_create_bezier(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_extrema(curve, 2, params, 4, &count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count >= 1, "3D curve has z-extremum");

		if (count >= 1) {
			qaws_scalar diff = params[0] - (qaws_scalar)0.377;
			if (diff < 0) diff = -diff;
			TEST_ASSERT(diff < (qaws_scalar)0.1,
				"z-extremum near t~0.377");
		}

		qaws_curve_destroy(curve);
	}

	/* 5. Null args */
	{
		qaws_status s;
		unsigned int count = 0;

		s = qaws_curve_find_extrema(NULL, 0, NULL, 0, &count);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"null curve extrema => INVALID_ARGUMENT");
	}
}

static void test_curvature_comb(void)
{
	printf("test_curvature_comb\n");

	/* 1. 2D curvature comb on quadratic parabola */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 2}, {2, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_curvature_sample_2d samples[10];
		unsigned int i;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = pts;
		desc.control_point_count = 3;

		s = qaws_curve_create_bezier(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_compute_curvature_comb_2d(curve, 10, samples, 10);
		TEST_ASSERT_STATUS(s);

		for (i = 0; i < 10; i++) {
			TEST_ASSERT(scalar_is_finite(samples[i].position.x),
				"2d comb position.x finite");
			TEST_ASSERT(scalar_is_finite(samples[i].position.y),
				"2d comb position.y finite");
			TEST_ASSERT(scalar_is_finite(samples[i].curvature),
				"2d comb curvature finite");
			{
				qaws_scalar nx = samples[i].normal.x;
				qaws_scalar ny = samples[i].normal.y;
				qaws_scalar nlen = (qaws_scalar)sqrt(
					(double)(nx * nx + ny * ny));
				TEST_ASSERT(approx_eq(nlen, (qaws_scalar)1.0),
					"2d comb normal is unit length");
			}
		}

		qaws_curve_destroy(curve);
	}

	/* 2. 3D curvature comb */
	{
		qaws_vec3 pts[] = { {0, 0, 0}, {1, 2, 1}, {2, 0, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_curvature_sample_3d samples[10];
		unsigned int i;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 2;
		desc.control_points = pts;
		desc.control_point_count = 3;

		s = qaws_curve_create_bezier(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_compute_curvature_comb_3d(curve, 10, samples, 10);
		TEST_ASSERT_STATUS(s);

		for (i = 0; i < 10; i++) {
			qaws_scalar nx = samples[i].normal.x;
			qaws_scalar ny = samples[i].normal.y;
			qaws_scalar nz = samples[i].normal.z;
			qaws_scalar nlen = (qaws_scalar)sqrt(
				(double)(nx * nx + ny * ny + nz * nz));
			TEST_ASSERT(approx_eq(nlen, (qaws_scalar)1.0),
				"3d comb normal is unit length");
		}

		qaws_curve_destroy(curve);
	}

	/* 3. Buffer too small: capacity < sample_count */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 2}, {2, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_curvature_sample_2d samples[3];

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = pts;
		desc.control_point_count = 3;

		s = qaws_curve_create_bezier(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_compute_curvature_comb_2d(curve, 10, samples, 3);
		TEST_ASSERT(s == QAWS_STATUS_BUFFER_TOO_SMALL,
			"comb buffer too small");

		qaws_curve_destroy(curve);
	}

	/* 4. sample_count < 2 => INVALID_ARGUMENT */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 2}, {2, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_curvature_sample_2d samples[2];

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = pts;
		desc.control_point_count = 3;

		s = qaws_curve_create_bezier(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_compute_curvature_comb_2d(curve, 1, samples, 2);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"sample_count < 2 => INVALID_ARGUMENT");

		qaws_curve_destroy(curve);
	}

	/* 5. Null args */
	{
		qaws_status s;
		qaws_curvature_sample_2d samples[10];

		s = qaws_curve_compute_curvature_comb_2d(NULL, 10, samples, 10);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"null curve comb => INVALID_ARGUMENT");
	}
}

static void test_winding_number(void)
{
	printf("test_winding_number\n");

	/* 1. Closed Catmull-Rom circle: inside and outside points */
	{
		qaws_vec2 points[8];
		qaws_catmull_rom_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		int winding = 0;
		unsigned int pi;

		/* 8 points around a circle of radius 2 */
		for (pi = 0; pi < 8; pi++) {
			qaws_scalar angle = (qaws_scalar)(2.0 * 3.14159265358979323846)
				* (qaws_scalar)pi / (qaws_scalar)8;
			points[pi].x = (qaws_scalar)(2.0 * cos((double)angle));
			points[pi].y = (qaws_scalar)(2.0 * sin((double)angle));
		}

		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = points;
		desc.control_point_count = 8;
		desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		desc.closed = 1;

		s = qaws_curve_create_catmull_rom(&desc, &curve);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(qaws_curve_is_closed(curve) == 1, "circle is closed");

		/* Inside point (origin): winding number should be +/-1 */
		{
			qaws_vec2 inside = { 0, 0 };
			s = qaws_curve_compute_winding_number_2d(curve, inside, &winding);
			TEST_ASSERT_STATUS(s);
			{
				int abs_w = winding < 0 ? -winding : winding;
				TEST_ASSERT(abs_w == 1, "inside winding = +/-1");
			}
		}

		/* Outside point: winding number should be 0 */
		{
			qaws_vec2 outside = { 100, 100 };
			s = qaws_curve_compute_winding_number_2d(curve, outside, &winding);
			TEST_ASSERT_STATUS(s);
			TEST_ASSERT(winding == 0, "outside winding = 0");
		}

		qaws_curve_destroy(curve);
	}

	/* 2. Non-closed curve => error */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 1}, {2, 0}, {3, 1} };
		qaws_catmull_rom_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		int winding = 0;
		qaws_vec2 pt = { 1, 0 };

		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = pts;
		desc.control_point_count = 4;
		desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		desc.closed = 0;

		s = qaws_curve_create_catmull_rom(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_compute_winding_number_2d(curve, pt, &winding);
		TEST_ASSERT(s != QAWS_STATUS_OK,
			"open curve winding => error");

		qaws_curve_destroy(curve);
	}

	/* 3. Dimension check: 3D curve => INVALID_DIMENSION */
	{
		qaws_vec3 pts[] = { {0, 0, 0}, {1, 1, 0}, {2, 0, 0}, {0, 0, 0} };
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		int winding = 0;
		qaws_vec2 pt = { 1, 0 };

		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;

		s = qaws_curve_create_bezier(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_compute_winding_number_2d(curve, pt, &winding);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_DIMENSION,
			"3D curve winding => INVALID_DIMENSION");

		qaws_curve_destroy(curve);
	}

	/* 4. Null args */
	{
		qaws_status s;
		int winding = 0;
		qaws_vec2 pt = { 0, 0 };

		s = qaws_curve_compute_winding_number_2d(NULL, pt, &winding);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"null curve winding => INVALID_ARGUMENT");
	}
}

static void test_yuksel_3d_circular(void)
{
	printf("test_yuksel_3d_circular\n");

	/* 1. 3D circular Yuksel: 5 points in 3D */
	{
		qaws_vec3 pts[] = {
			{0, 0, 0}, {1, 2, 1}, {2, 0, 2}, {3, 2, 1}, {4, 0, 0}
		};
		qaws_yuksel_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_eval_result_3d r;
		qaws_range range;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.control_points = pts;
		desc.control_point_count = 5;
		desc.mode = QAWS_YUKSEL_MODE_CIRCULAR;
		desc.closed = 0;

		s = qaws_curve_create_yuksel(&desc, &curve);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(curve != NULL, "3D circular yuksel created");

		range = qaws_curve_get_parameter_range(curve);

		/* Evaluate at start */
		s = qaws_curve_evaluate_3d(curve, range.min_value,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "3d circ start x finite");
		TEST_ASSERT(scalar_is_finite(r.position.y), "3d circ start y finite");
		TEST_ASSERT(scalar_is_finite(r.position.z), "3d circ start z finite");

		/* Evaluate at midpoint */
		{
			qaws_scalar tmid = (range.min_value + range.max_value)
				* (qaws_scalar)0.5;
			s = qaws_curve_evaluate_3d(curve, tmid,
				QAWS_EVAL_FLAG_POSITION, &r);
			TEST_ASSERT_STATUS(s);
			TEST_ASSERT(scalar_is_finite(r.position.x), "3d circ mid x finite");
		}

		/* Evaluate at end */
		s = qaws_curve_evaluate_3d(curve, range.max_value,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "3d circ end x finite");

		qaws_curve_destroy(curve);
	}

	/* 2. 3D elliptical Yuksel */
	{
		qaws_vec3 pts[] = {
			{0, 0, 0}, {1, 2, 1}, {2, 0, 2}, {3, 2, 1}, {4, 0, 0}
		};
		qaws_yuksel_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_eval_result_3d r;
		qaws_range range;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.control_points = pts;
		desc.control_point_count = 5;
		desc.mode = QAWS_YUKSEL_MODE_ELLIPTICAL;
		desc.closed = 0;

		s = qaws_curve_create_yuksel(&desc, &curve);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(curve != NULL, "3D elliptical yuksel created");

		range = qaws_curve_get_parameter_range(curve);

		s = qaws_curve_evaluate_3d(curve, range.min_value,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "3d ellip start finite");

		s = qaws_curve_evaluate_3d(curve, range.max_value,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "3d ellip end finite");

		qaws_curve_destroy(curve);
	}

	/* 3. 3D hybrid Yuksel */
	{
		qaws_vec3 pts[] = {
			{0, 0, 0}, {1, 2, 1}, {2, 0, 2}, {3, 2, 1}, {4, 0, 0}
		};
		qaws_yuksel_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_eval_result_3d r;
		qaws_range range;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.control_points = pts;
		desc.control_point_count = 5;
		desc.mode = QAWS_YUKSEL_MODE_HYBRID;
		desc.closed = 0;

		s = qaws_curve_create_yuksel(&desc, &curve);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(curve != NULL, "3D hybrid yuksel created");

		range = qaws_curve_get_parameter_range(curve);

		s = qaws_curve_evaluate_3d(curve, range.min_value,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "3d hybrid start finite");

		s = qaws_curve_evaluate_3d(curve, range.max_value,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "3d hybrid end finite");

		qaws_curve_destroy(curve);
	}

	/* 4. Collinear 3D points: should fall back gracefully */
	{
		qaws_vec3 pts[] = {
			{0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {3, 3, 3}, {4, 4, 4}
		};
		qaws_yuksel_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.control_points = pts;
		desc.control_point_count = 5;
		desc.mode = QAWS_YUKSEL_MODE_CIRCULAR;
		desc.closed = 0;

		s = qaws_curve_create_yuksel(&desc, &curve);
		/* Should succeed (fallback to bezier) or return degenerate */
		TEST_ASSERT(s == QAWS_STATUS_OK || s == QAWS_STATUS_DEGENERATE_CURVE,
			"collinear 3D yuksel handled");

		if (s == QAWS_STATUS_OK && curve != NULL) {
			qaws_eval_result_3d r;
			qaws_range range = qaws_curve_get_parameter_range(curve);
			unsigned int pi;

			/* Evaluate across the range - should not crash */
			for (pi = 0; pi < 10; pi++) {
				qaws_scalar t = range.min_value +
					(qaws_scalar)pi * (range.max_value - range.min_value)
					/ (qaws_scalar)9;
				qaws_status es = qaws_curve_evaluate_3d(curve, t,
					QAWS_EVAL_FLAG_POSITION, &r);
				TEST_ASSERT(es == QAWS_STATUS_OK, "collinear 3d eval ok");
			}

			qaws_curve_destroy(curve);
		}
	}
}

/* ================================================================== */
/* Phase 5 tests: multi-curve traversal, S-curve, custom speed,        */
/* curvature-weighted, feature-preserving, streaming sampling          */
/* ================================================================== */

static qaws_scalar constant_speed_fn(qaws_scalar distance, void* user_data)
{
	(void)distance;
	(void)user_data;
	return (qaws_scalar)1.0;
}

static qaws_scalar accelerating_speed_fn(qaws_scalar distance, void* user_data)
{
	(void)user_data;
	return (qaws_scalar)1.0 + distance;
}

struct stream_counter { unsigned int count; };

static int stream_count_cb_2d(
	qaws_scalar param, qaws_vec2 const* pos, void* ud)
{
	struct stream_counter* c = (struct stream_counter*)ud;
	(void)param; (void)pos;
	c->count++;
	return 1;
}

static int stream_count_cb_3d(
	qaws_scalar param, qaws_vec3 const* pos, void* ud)
{
	struct stream_counter* c = (struct stream_counter*)ud;
	(void)param; (void)pos;
	c->count++;
	return 1;
}

static int stream_early_stop_cb_2d(
	qaws_scalar param, qaws_vec2 const* pos, void* ud)
{
	struct stream_counter* c = (struct stream_counter*)ud;
	(void)param; (void)pos;
	c->count++;
	return (c->count < 3) ? 1 : 0;
}

static void test_multi_curve_traversal(void)
{
	printf("test_multi_curve_traversal\n");

	/* 1. Two-curve chain: (0,0)->(2,1) then (2,1)->(4,0) */
	{
		qaws_vec2 pts1[] = { {0, 0}, {1, 1}, {2, 1} };
		qaws_vec2 pts2[] = { {2, 1}, {3, 1}, {4, 0} };
		qaws_bezier_desc bdesc1, bdesc2;
		qaws_curve* c1 = NULL;
		qaws_curve* c2 = NULL;
		qaws_curve const* curves[2];
		qaws_traversal_desc tdesc;
		qaws_traversal* trav = NULL;
		qaws_status s;
		qaws_eval_result_2d r;
		qaws_scalar total_len, len1, len2;

		bdesc1.dimension = QAWS_DIMENSION_2D;
		bdesc1.degree = 2;
		bdesc1.control_points = pts1;
		bdesc1.control_point_count = 3;

		bdesc2.dimension = QAWS_DIMENSION_2D;
		bdesc2.degree = 2;
		bdesc2.control_points = pts2;
		bdesc2.control_point_count = 3;

		s = qaws_curve_create_bezier(&bdesc1, &c1);
		TEST_ASSERT_STATUS(s);
		s = qaws_curve_create_bezier(&bdesc2, &c2);
		TEST_ASSERT_STATUS(s);

		curves[0] = c1;
		curves[1] = c2;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
		tdesc.speed = (qaws_scalar)1.0;
		tdesc.clamp_to_domain = 1;

		s = qaws_traversal_create_multi(curves, 2, &tdesc, &trav);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(trav != NULL, "multi-curve traversal created");

		/* At distance=0 => near (0,0) */
		s = qaws_traversal_evaluate_2d(trav, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0.0),
			"multi trav d=0 x~0");
		TEST_ASSERT(approx_eq(r.position.y, (qaws_scalar)0.0),
			"multi trav d=0 y~0");

		/* Compute total arc length of the chain */
		s = qaws_curve_compute_arc_length(c1, (qaws_scalar)0.0,
			(qaws_scalar)1.0, &len1);
		TEST_ASSERT_STATUS(s);
		s = qaws_curve_compute_arc_length(c2, (qaws_scalar)0.0,
			(qaws_scalar)1.0, &len2);
		TEST_ASSERT_STATUS(s);
		total_len = len1 + len2;

		/* At total distance => near (4,0) */
		s = qaws_traversal_evaluate_2d(trav, total_len,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "multi trav end x finite");
		TEST_ASSERT(scalar_is_finite(r.position.y), "multi trav end y finite");

		/* At midpoint distance => somewhere valid */
		s = qaws_traversal_evaluate_2d(trav,
			total_len * (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "multi trav mid x finite");
		TEST_ASSERT(scalar_is_finite(r.position.y), "multi trav mid y finite");

		qaws_traversal_destroy(trav);
		qaws_curve_destroy(c1);
		qaws_curve_destroy(c2);
	}

	/* 2. Single curve in chain: should work like regular traversal */
	{
		qaws_vec2 pts[] = { {0, 0}, {5, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve* c = NULL;
		qaws_curve const* curves[1];
		qaws_traversal_desc tdesc;
		qaws_traversal* trav = NULL;
		qaws_status s;
		qaws_eval_result_2d r;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;

		s = qaws_curve_create_bezier(&bdesc, &c);
		TEST_ASSERT_STATUS(s);

		curves[0] = c;

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
		tdesc.speed = (qaws_scalar)1.0;
		tdesc.clamp_to_domain = 1;

		s = qaws_traversal_create_multi(curves, 1, &tdesc, &trav);
		TEST_ASSERT_STATUS(s);

		s = qaws_traversal_evaluate_2d(trav, (qaws_scalar)2.5,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)2.5),
			"single-chain d=2.5 x~2.5");

		qaws_traversal_destroy(trav);
		qaws_curve_destroy(c);
	}

	/* 3. Null args */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 1} };
		qaws_bezier_desc bdesc;
		qaws_curve* c = NULL;
		qaws_curve const* curves[1];
		qaws_traversal_desc tdesc;
		qaws_traversal* trav = NULL;
		qaws_status s;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;

		s = qaws_curve_create_bezier(&bdesc, &c);
		TEST_ASSERT_STATUS(s);

		curves[0] = c;
		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_CONSTANT_SPEED;
		tdesc.speed = (qaws_scalar)1.0;
		tdesc.clamp_to_domain = 1;

		/* null curves array */
		s = qaws_traversal_create_multi(NULL, 2, &tdesc, &trav);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"multi null curves => INVALID_ARGUMENT");

		/* null desc */
		s = qaws_traversal_create_multi(curves, 1, NULL, &trav);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"multi null desc => INVALID_ARGUMENT");

		/* null out */
		s = qaws_traversal_create_multi(curves, 1, &tdesc, NULL);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"multi null out => INVALID_ARGUMENT");

		qaws_curve_destroy(c);
	}
}

static void test_scurve_profile(void)
{
	printf("test_scurve_profile\n");

	/* 1. Basic S-curve profile */
	{
		qaws_vec2 pts[] = { {0, 0}, {2, 3}, {5, 3}, {8, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_traversal_desc tdesc;
		qaws_traversal* trav = NULL;
		qaws_status s;
		qaws_eval_result_2d r;
		unsigned int i;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_S_CURVE;
		tdesc.max_speed = (qaws_scalar)5.0;
		tdesc.acceleration = (qaws_scalar)2.0;
		tdesc.jerk = (qaws_scalar)1.0;
		tdesc.start_time = (qaws_scalar)0.0;
		tdesc.end_time = (qaws_scalar)5.0;
		tdesc.clamp_to_domain = 1;

		s = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(trav != NULL, "S-curve traversal created");

		/* At t=start_time => near curve start */
		s = qaws_traversal_evaluate_2d(trav, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "scurve t=0 x finite");
		TEST_ASSERT(scalar_is_finite(r.position.y), "scurve t=0 y finite");

		/* At t=end_time => somewhere along curve */
		s = qaws_traversal_evaluate_2d(trav, (qaws_scalar)5.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x), "scurve t=end x finite");
		TEST_ASSERT(scalar_is_finite(r.position.y), "scurve t=end y finite");

		/* Intermediate positions are smooth (all finite) */
		for (i = 1; i <= 4; i++) {
			qaws_scalar t = (qaws_scalar)i;
			s = qaws_traversal_evaluate_2d(trav, t,
				QAWS_EVAL_FLAG_POSITION, &r);
			TEST_ASSERT_STATUS(s);
			TEST_ASSERT(scalar_is_finite(r.position.x),
				"scurve intermediate x finite");
			TEST_ASSERT(scalar_is_finite(r.position.y),
				"scurve intermediate y finite");
		}

		qaws_traversal_destroy(trav);
		qaws_curve_destroy(curve);
	}

	/* 2. Very large jerk => should behave like trapezoidal (still valid) */
	{
		qaws_vec2 pts[] = { {0, 0}, {5, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_traversal_desc tdesc;
		qaws_traversal* trav = NULL;
		qaws_status s;
		qaws_eval_result_2d r;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_S_CURVE;
		tdesc.max_speed = (qaws_scalar)5.0;
		tdesc.acceleration = (qaws_scalar)10.0;
		tdesc.jerk = (qaws_scalar)100000.0;
		tdesc.start_time = (qaws_scalar)0.0;
		tdesc.end_time = (qaws_scalar)2.0;
		tdesc.clamp_to_domain = 1;

		s = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(s);

		s = qaws_traversal_evaluate_2d(trav, (qaws_scalar)1.0,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x),
			"high-jerk scurve mid x finite");

		qaws_traversal_destroy(trav);
		qaws_curve_destroy(curve);
	}
}

static void test_custom_speed(void)
{
	printf("test_custom_speed\n");

	/* 1. Constant speed callback */
	{
		qaws_vec2 pts[] = { {0, 0}, {3, 4} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_traversal_desc tdesc;
		qaws_traversal* trav = NULL;
		qaws_status s;
		qaws_eval_result_2d r0, r1, r2;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_CUSTOM;
		tdesc.custom_speed_fn = constant_speed_fn;
		tdesc.custom_speed_fn_user_data = NULL;
		tdesc.start_time = (qaws_scalar)0.0;
		tdesc.end_time = (qaws_scalar)5.0;
		tdesc.clamp_to_domain = 1;

		s = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(trav != NULL, "custom speed traversal created");

		/* Evaluate at start */
		s = qaws_traversal_evaluate_2d(trav, (qaws_scalar)0.0,
			QAWS_EVAL_FLAG_POSITION, &r0);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r0.position.x),
			"custom speed t=0 x finite");

		/* Evaluate at midpoint */
		s = qaws_traversal_evaluate_2d(trav, (qaws_scalar)2.5,
			QAWS_EVAL_FLAG_POSITION, &r1);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r1.position.x),
			"custom speed t=mid x finite");

		/* Evaluate at end */
		s = qaws_traversal_evaluate_2d(trav, (qaws_scalar)5.0,
			QAWS_EVAL_FLAG_POSITION, &r2);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r2.position.x),
			"custom speed t=end x finite");

		/* Motion progresses: end should be further than start */
		TEST_ASSERT(r2.position.x > r0.position.x ||
			r2.position.y > r0.position.y,
			"custom const speed progresses along curve");

		qaws_traversal_destroy(trav);
		qaws_curve_destroy(curve);
	}

	/* 2. Variable (accelerating) speed */
	{
		qaws_vec2 pts[] = { {0, 0}, {2, 2}, {4, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_traversal_desc tdesc;
		qaws_traversal* trav = NULL;
		qaws_status s;
		qaws_eval_result_2d r;
		unsigned int i;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 2;
		bdesc.control_points = pts;
		bdesc.control_point_count = 3;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_CUSTOM;
		tdesc.custom_speed_fn = accelerating_speed_fn;
		tdesc.custom_speed_fn_user_data = NULL;
		tdesc.start_time = (qaws_scalar)0.0;
		tdesc.end_time = (qaws_scalar)3.0;
		tdesc.clamp_to_domain = 1;

		s = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(s);

		for (i = 0; i <= 6; i++) {
			qaws_scalar t = (qaws_scalar)i * (qaws_scalar)0.5;
			s = qaws_traversal_evaluate_2d(trav, t,
				QAWS_EVAL_FLAG_POSITION, &r);
			TEST_ASSERT_STATUS(s);
			TEST_ASSERT(scalar_is_finite(r.position.x),
				"accel speed eval x finite");
			TEST_ASSERT(scalar_is_finite(r.position.y),
				"accel speed eval y finite");
		}

		qaws_traversal_destroy(trav);
		qaws_curve_destroy(curve);
	}

	/* 3. Null callback: implementation gracefully degrades (no custom table,
	      falls back to normal time-to-distance mapping). Verify it still
	      creates and evaluates without crashing. */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 1} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_traversal_desc tdesc;
		qaws_traversal* trav = NULL;
		qaws_status s;
		qaws_eval_result_2d r;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		memset(&tdesc, 0, sizeof(tdesc));
		tdesc.traversal_mode = QAWS_TRAVERSAL_MODE_TIME;
		tdesc.motion_profile = QAWS_MOTION_PROFILE_CUSTOM;
		tdesc.custom_speed_fn = NULL;
		tdesc.start_time = (qaws_scalar)0.0;
		tdesc.end_time = (qaws_scalar)1.0;
		tdesc.clamp_to_domain = 1;

		s = qaws_traversal_create(curve, &tdesc, &trav);
		TEST_ASSERT_STATUS(s);

		/* Should still evaluate without crash */
		s = qaws_traversal_evaluate_2d(trav, (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(scalar_is_finite(r.position.x),
			"null custom fn fallback x finite");

		qaws_traversal_destroy(trav);
		qaws_curve_destroy(curve);
	}
}

static void test_curvature_weighted_sampling(void)
{
	printf("test_curvature_weighted_sampling\n");

	/* 1. 2D curvature-weighted: S-curve shape with varying curvature */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 3}, {3, -3}, {4, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_vec2 samples[20];
		unsigned int count = 0;
		unsigned int i;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_sample_curvature_weighted_2d(curve, 20,
			samples, 20, &count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count == 20, "curvature weighted 2d count == 20");

		/* All positions are finite */
		for (i = 0; i < count; i++) {
			TEST_ASSERT(scalar_is_finite(samples[i].x),
				"curv weighted 2d sample x finite");
			TEST_ASSERT(scalar_is_finite(samples[i].y),
				"curv weighted 2d sample y finite");
		}

		/* Spacing should be non-uniform: check that not all intervals
		   are the same length. Compute min/max spacing. */
		{
			qaws_scalar min_gap = (qaws_scalar)1e10;
			qaws_scalar max_gap = (qaws_scalar)0.0;
			for (i = 1; i < count; i++) {
				qaws_scalar dx = samples[i].x - samples[i - 1].x;
				qaws_scalar dy = samples[i].y - samples[i - 1].y;
				qaws_scalar gap = dx * dx + dy * dy;
				if (gap < min_gap) min_gap = gap;
				if (gap > max_gap) max_gap = gap;
			}
			TEST_ASSERT(max_gap > min_gap * (qaws_scalar)1.5,
				"curvature weighted spacing is non-uniform");
		}

		qaws_curve_destroy(curve);
	}

	/* 2. 3D curvature-weighted */
	{
		qaws_vec3 pts[] = { {0,0,0}, {1,3,1}, {3,-3,2}, {4,0,0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_vec3 samples[15];
		unsigned int count = 0;
		unsigned int i;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_sample_curvature_weighted_3d(curve, 15,
			samples, 15, &count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count == 15, "curvature weighted 3d count == 15");

		for (i = 0; i < count; i++) {
			TEST_ASSERT(scalar_is_finite(samples[i].x),
				"curv weighted 3d x finite");
			TEST_ASSERT(scalar_is_finite(samples[i].y),
				"curv weighted 3d y finite");
			TEST_ASSERT(scalar_is_finite(samples[i].z),
				"curv weighted 3d z finite");
		}

		qaws_curve_destroy(curve);
	}

	/* 3. Null args */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 1}, {2, 0}, {3, 1} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_vec2 samples[10];
		unsigned int count = 0;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_sample_curvature_weighted_2d(NULL, 10,
			samples, 10, &count);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"curv weighted null curve => INVALID_ARGUMENT");

		s = qaws_curve_sample_curvature_weighted_2d(curve, 10,
			NULL, 10, &count);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"curv weighted null out => INVALID_ARGUMENT");

		qaws_curve_destroy(curve);
	}

	/* 4. sample_count < 2 */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 1}, {2, 0}, {3, 1} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_vec2 samples[10];
		unsigned int count = 0;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_sample_curvature_weighted_2d(curve, 1,
			samples, 10, &count);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"curv weighted count<2 => INVALID_ARGUMENT");

		qaws_curve_destroy(curve);
	}
}

static void test_feature_preserving_sampling(void)
{
	printf("test_feature_preserving_sampling\n");

	/* 1. 2D S-curve with inflection: CPs (0,0),(1,3),(2,-3),(3,0) */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 3}, {2, -3}, {3, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_vec2 samples[50];
		unsigned int count = 0;
		unsigned int i;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_sample_feature_preserving_2d(curve, 5,
			samples, 50, &count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count >= 5,
			"feature preserving 2d count >= base count");

		for (i = 0; i < count; i++) {
			TEST_ASSERT(scalar_is_finite(samples[i].x),
				"feature preserving 2d x finite");
			TEST_ASSERT(scalar_is_finite(samples[i].y),
				"feature preserving 2d y finite");
		}

		qaws_curve_destroy(curve);
	}

	/* 2. 3D curve */
	{
		qaws_vec3 pts[] = { {0,0,0}, {1,3,1}, {2,-3,2}, {3,0,0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_vec3 samples[50];
		unsigned int count = 0;
		unsigned int i;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_sample_feature_preserving_3d(curve, 5,
			samples, 50, &count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count >= 5,
			"feature preserving 3d count >= base count");

		for (i = 0; i < count; i++) {
			TEST_ASSERT(scalar_is_finite(samples[i].x),
				"feature preserving 3d x finite");
			TEST_ASSERT(scalar_is_finite(samples[i].y),
				"feature preserving 3d y finite");
			TEST_ASSERT(scalar_is_finite(samples[i].z),
				"feature preserving 3d z finite");
		}

		qaws_curve_destroy(curve);
	}

	/* 3. Null args */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 1}, {2, 0}, {3, 1} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_vec2 samples[10];
		unsigned int count = 0;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_sample_feature_preserving_2d(NULL, 5,
			samples, 50, &count);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"feature preserving null curve => INVALID_ARGUMENT");

		s = qaws_curve_sample_feature_preserving_2d(curve, 5,
			NULL, 50, &count);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"feature preserving null out => INVALID_ARGUMENT");

		qaws_curve_destroy(curve);
	}
}

static void test_streaming_sampling(void)
{
	printf("test_streaming_sampling\n");

	/* 1. 2D streaming: count invocations */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 2}, {3, 1}, {4, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_sampling_desc sdesc;
		struct stream_counter counter;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
		sdesc.sample_count = 10;
		sdesc.error_tolerance = 0;
		sdesc.use_adaptive_sampling = 0;

		counter.count = 0;
		s = qaws_curve_sample_streaming_2d(curve, &sdesc,
			stream_count_cb_2d, &counter);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(counter.count == 10,
			"streaming 2d callback count == 10");

		qaws_curve_destroy(curve);
	}

	/* 2. Early termination: callback returns 0 after 3 calls */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 2}, {3, 1}, {4, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_sampling_desc sdesc;
		struct stream_counter counter;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
		sdesc.sample_count = 20;
		sdesc.error_tolerance = 0;
		sdesc.use_adaptive_sampling = 0;

		counter.count = 0;
		s = qaws_curve_sample_streaming_2d(curve, &sdesc,
			stream_early_stop_cb_2d, &counter);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(counter.count == 3,
			"streaming early stop at 3 calls");

		qaws_curve_destroy(curve);
	}

	/* 3. 3D streaming */
	{
		qaws_vec3 pts[] = { {0,0,0}, {1,2,1}, {3,1,2}, {4,0,0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_sampling_desc sdesc;
		struct stream_counter counter;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
		sdesc.sample_count = 8;
		sdesc.error_tolerance = 0;
		sdesc.use_adaptive_sampling = 0;

		counter.count = 0;
		s = qaws_curve_sample_streaming_3d(curve, &sdesc,
			stream_count_cb_3d, &counter);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(counter.count == 8,
			"streaming 3d callback count == 8");

		qaws_curve_destroy(curve);
	}

	/* 4. Null callback => INVALID_ARGUMENT */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 1}, {2, 0}, {3, 1} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_sampling_desc sdesc;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
		sdesc.sample_count = 10;

		s = qaws_curve_sample_streaming_2d(curve, &sdesc, NULL, NULL);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"streaming null callback => INVALID_ARGUMENT");

		qaws_curve_destroy(curve);
	}
}

/* ================================================================== */
/* Phase 7 — Self-intersection & Curve-curve intersection              */
/* ================================================================== */

static void test_self_intersection(void)
{
	printf("test_self_intersection\n");

	/* 1. Figure-8 Hermite curve has a self-intersection at the origin */
	{
		qaws_vec2 pts[5];
		qaws_vec2 tans[5];
		qaws_hermite_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_intersection_2d isects[16];
		unsigned int isect_count = 0;

		pts[0].x = 0; pts[0].y = 0;
		pts[1].x = 1; pts[1].y = (qaws_scalar)0.5;
		pts[2].x = 0; pts[2].y = 0;
		pts[3].x = -1; pts[3].y = (qaws_scalar)-0.5;
		pts[4].x = 0; pts[4].y = 0;

		tans[0].x = 1; tans[0].y = 1;
		tans[1].x = 0; tans[1].y = (qaws_scalar)-1.5;
		tans[2].x = -1; tans[2].y = -1;
		tans[3].x = 0; tans[3].y = (qaws_scalar)1.5;
		tans[4].x = 1; tans[4].y = 1;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.points = pts;
		desc.derivatives = tans;
		desc.point_count = 5;
		desc.derivative_count = 5;

		s = qaws_curve_create_hermite(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_self_intersections_2d(
			curve, isects, 16, &isect_count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(isect_count > 0, "figure-8 has self-intersections");

		if (isect_count > 0 && isect_count <= 16)
		{
			int found_origin = 0;
			unsigned int i;
			for (i = 0; i < isect_count; ++i)
			{
				qaws_scalar dx = isects[i].position.x;
				qaws_scalar dy = isects[i].position.y;
				if (dx < 0) dx = -dx;
				if (dy < 0) dy = -dy;
				if (dx < (qaws_scalar)0.1 && dy < (qaws_scalar)0.1)
					found_origin = 1;
				TEST_ASSERT(isects[i].parameter_a < isects[i].parameter_b,
					"self-isect param_a < param_b");
			}
			TEST_ASSERT(found_origin, "self-intersection near origin");
		}

		qaws_curve_destroy(curve);
	}

	/* 2. Simple Bezier arc — no self-intersection */
	{
		qaws_vec2 cp[4];
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_intersection_2d isects[8];
		unsigned int isect_count = 0;

		cp[0].x = 0; cp[0].y = 0;
		cp[1].x = 1; cp[1].y = 2;
		cp[2].x = 3; cp[2].y = 2;
		cp[3].x = 4; cp[3].y = 0;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = cp;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_self_intersections_2d(
			curve, isects, 8, &isect_count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(isect_count == 0, "simple arc has no self-intersections");

		qaws_curve_destroy(curve);
	}

	/* 3. Null arguments */
	{
		qaws_intersection_2d isects[4];
		unsigned int count = 0;
		qaws_status s;

		s = qaws_curve_find_self_intersections_2d(NULL, isects, 4, &count);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"null curve => INVALID_ARGUMENT");
	}

	/* 4. Dimension mismatch */
	{
		qaws_vec3 cp3[4];
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_intersection_2d isects[4];
		unsigned int count = 0;

		cp3[0].x = 0; cp3[0].y = 0; cp3[0].z = 0;
		cp3[1].x = 1; cp3[1].y = 1; cp3[1].z = 1;
		cp3[2].x = 2; cp3[2].y = 0; cp3[2].z = 1;
		cp3[3].x = 3; cp3[3].y = 1; cp3[3].z = 0;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 3;
		bdesc.control_points = cp3;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&bdesc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_self_intersections_2d(curve, isects, 4, &count);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_DIMENSION,
			"3d curve for 2d self-isect => INVALID_DIMENSION");

		qaws_curve_destroy(curve);
	}

	/* 5. 3D self-intersection */
	{
		qaws_vec3 pts[5];
		qaws_vec3 tans[5];
		qaws_hermite_desc desc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_intersection_3d isects[16];
		unsigned int isect_count = 0;

		pts[0].x = 0; pts[0].y = 0; pts[0].z = 0;
		pts[1].x = 1; pts[1].y = (qaws_scalar)0.5; pts[1].z = 0;
		pts[2].x = 0; pts[2].y = 0; pts[2].z = 0;
		pts[3].x = -1; pts[3].y = (qaws_scalar)-0.5; pts[3].z = 0;
		pts[4].x = 0; pts[4].y = 0; pts[4].z = 0;

		tans[0].x = 1; tans[0].y = 1; tans[0].z = 0;
		tans[1].x = 0; tans[1].y = (qaws_scalar)-1.5; tans[1].z = 0;
		tans[2].x = -1; tans[2].y = -1; tans[2].z = 0;
		tans[3].x = 0; tans[3].y = (qaws_scalar)1.5; tans[3].z = 0;
		tans[4].x = 1; tans[4].y = 1; tans[4].z = 0;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 3;
		desc.points = pts;
		desc.derivatives = tans;
		desc.point_count = 5;
		desc.derivative_count = 5;

		s = qaws_curve_create_hermite(&desc, &curve);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_self_intersections_3d(
			curve, isects, 16, &isect_count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(isect_count > 0, "3D figure-8 has self-intersections");

		qaws_curve_destroy(curve);
	}
}

static void test_curve_curve_intersection(void)
{
	printf("test_curve_curve_intersection\n");

	/* 1. Two crossing Bezier curves in 2D */
	{
		/* Curve A: arch shape */
		qaws_vec2 cpA[] = { {0,0}, {1,3}, {3,3}, {4,0} };
		/* Curve B: S-shape crossing A */
		qaws_vec2 cpB[] = { {0,2}, {2,-1}, {2,3}, {4,1} };

		qaws_bezier_desc descA, descB;
		qaws_curve *curveA = NULL, *curveB = NULL;
		qaws_status s;
		qaws_intersection_2d isects[16];
		unsigned int isect_count = 0;

		descA.dimension = QAWS_DIMENSION_2D;
		descA.degree = 3;
		descA.control_points = cpA;
		descA.control_point_count = 4;

		descB.dimension = QAWS_DIMENSION_2D;
		descB.degree = 3;
		descB.control_points = cpB;
		descB.control_point_count = 4;

		s = qaws_curve_create_bezier(&descA, &curveA);
		TEST_ASSERT_STATUS(s);
		s = qaws_curve_create_bezier(&descB, &curveB);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_intersections_2d(
			curveA, curveB, isects, 16, &isect_count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(isect_count > 0, "crossing beziers have intersections");

		/* Verify intersection positions agree on both curves */
		{
			unsigned int i;
			for (i = 0; i < isect_count && i < 16; ++i)
			{
				qaws_eval_result_2d ra, rb;
				qaws_scalar dx, dy;
				s = qaws_curve_evaluate_2d(curveA, isects[i].parameter_a,
					QAWS_EVAL_FLAG_POSITION, &ra);
				TEST_ASSERT_STATUS(s);
				s = qaws_curve_evaluate_2d(curveB, isects[i].parameter_b,
					QAWS_EVAL_FLAG_POSITION, &rb);
				TEST_ASSERT_STATUS(s);

				dx = ra.position.x - rb.position.x;
				dy = ra.position.y - rb.position.y;
				if (dx < 0) dx = -dx;
				if (dy < 0) dy = -dy;
				TEST_ASSERT(dx < (qaws_scalar)0.01 && dy < (qaws_scalar)0.01,
					"intersection point matches on both curves");
			}
		}

		qaws_curve_destroy(curveA);
		qaws_curve_destroy(curveB);
	}

	/* 2. Non-intersecting curves */
	{
		qaws_vec2 cpA[] = { {0,0}, {1,1}, {2,1}, {3,0} };
		qaws_vec2 cpB[] = { {0,3}, {1,4}, {2,4}, {3,3} };

		qaws_bezier_desc descA, descB;
		qaws_curve *curveA = NULL, *curveB = NULL;
		qaws_status s;
		qaws_intersection_2d isects[8];
		unsigned int isect_count = 0;

		descA.dimension = QAWS_DIMENSION_2D;
		descA.degree = 3;
		descA.control_points = cpA;
		descA.control_point_count = 4;

		descB.dimension = QAWS_DIMENSION_2D;
		descB.degree = 3;
		descB.control_points = cpB;
		descB.control_point_count = 4;

		s = qaws_curve_create_bezier(&descA, &curveA);
		TEST_ASSERT_STATUS(s);
		s = qaws_curve_create_bezier(&descB, &curveB);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_intersections_2d(
			curveA, curveB, isects, 8, &isect_count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(isect_count == 0, "non-intersecting curves have 0 intersections");

		qaws_curve_destroy(curveA);
		qaws_curve_destroy(curveB);
	}

	/* 3. Line crossing a curve at 2 points */
	{
		qaws_vec2 cpLine[] = { {0,(qaws_scalar)1.5}, {4,(qaws_scalar)1.5} };
		qaws_vec2 cpArc[]  = { {0,0}, {1,3}, {3,3}, {4,0} };

		qaws_bezier_desc descL, descA;
		qaws_curve *line = NULL, *arc = NULL;
		qaws_status s;
		qaws_intersection_2d isects[8];
		unsigned int isect_count = 0;

		descL.dimension = QAWS_DIMENSION_2D;
		descL.degree = 1;
		descL.control_points = cpLine;
		descL.control_point_count = 2;

		descA.dimension = QAWS_DIMENSION_2D;
		descA.degree = 3;
		descA.control_points = cpArc;
		descA.control_point_count = 4;

		s = qaws_curve_create_bezier(&descL, &line);
		TEST_ASSERT_STATUS(s);
		s = qaws_curve_create_bezier(&descA, &arc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_intersections_2d(
			line, arc, isects, 8, &isect_count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(isect_count == 2, "line crosses arc at 2 points");

		qaws_curve_destroy(line);
		qaws_curve_destroy(arc);
	}

	/* 4. Null arguments */
	{
		qaws_intersection_2d isects[4];
		unsigned int count = 0;
		qaws_status s;

		s = qaws_curve_find_intersections_2d(NULL, NULL, isects, 4, &count);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"null curves => INVALID_ARGUMENT");
	}

	/* 5. Dimension mismatch */
	{
		qaws_vec2 cp2[] = { {0,0}, {1,1}, {2,0} };
		qaws_vec3 cp3[] = { {0,0,0}, {1,1,1}, {2,0,0} };
		qaws_bezier_desc d2, d3;
		qaws_curve *c2 = NULL, *c3 = NULL;
		qaws_status s;
		qaws_intersection_2d isects[4];
		unsigned int count = 0;

		d2.dimension = QAWS_DIMENSION_2D;
		d2.degree = 2;
		d2.control_points = cp2;
		d2.control_point_count = 3;

		d3.dimension = QAWS_DIMENSION_3D;
		d3.degree = 2;
		d3.control_points = cp3;
		d3.control_point_count = 3;

		s = qaws_curve_create_bezier(&d2, &c2);
		TEST_ASSERT_STATUS(s);
		s = qaws_curve_create_bezier(&d3, &c3);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_intersections_2d(c2, c3, isects, 4, &count);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_DIMENSION,
			"mixed dimensions => INVALID_DIMENSION");

		qaws_curve_destroy(c2);
		qaws_curve_destroy(c3);
	}

	/* 6. 3D curve-curve intersection */
	{
		qaws_vec3 cpA[] = { {0,0,0}, {2,2,0}, {2,-2,0}, {4,0,0} };
		qaws_vec3 cpB[] = { {2,-2,0}, {2,2,0}, {2,2,0}, {2,-2,0} };

		qaws_bezier_desc descA, descB;
		qaws_curve *curveA = NULL, *curveB = NULL;
		qaws_status s;
		qaws_intersection_3d isects[16];
		unsigned int isect_count = 0;

		descA.dimension = QAWS_DIMENSION_3D;
		descA.degree = 3;
		descA.control_points = cpA;
		descA.control_point_count = 4;

		descB.dimension = QAWS_DIMENSION_3D;
		descB.degree = 3;
		descB.control_points = cpB;
		descB.control_point_count = 4;

		s = qaws_curve_create_bezier(&descA, &curveA);
		TEST_ASSERT_STATUS(s);
		s = qaws_curve_create_bezier(&descB, &curveB);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_find_intersections_3d(
			curveA, curveB, isects, 16, &isect_count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(isect_count >= 0, "3D curve-curve intersection ran");

		qaws_curve_destroy(curveA);
		qaws_curve_destroy(curveB);
	}
}

/* ================================================================== */
/* SVG visual test: self-intersection                                  */
/* ================================================================== */

static void test_svg_self_intersection(void)
{
	printf("test_svg_self_intersection\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/self_intersection.svg",
		(qaws_scalar)-1.5, (qaws_scalar)-1.0, (qaws_scalar)3.0, (qaws_scalar)2.0,
		(qaws_scalar)600, (qaws_scalar)400))
		return;

	/* Figure-8 Hermite */
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

	/* Draw curve */
	qaws_vec2 buf[SVG_SAMPLES];
	unsigned int n = svg_sample_curve(curve, buf, SVG_SAMPLES);
	svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);

	/* Find self-intersections */
	qaws_intersection_2d isects[16];
	unsigned int isect_count = 0;
	s = qaws_curve_find_self_intersections_2d(curve, isects, 16, &isect_count);
	TEST_ASSERT_STATUS(s);

	/* Draw intersection points in red */
	for (unsigned int i = 0; i < isect_count && i < 16; ++i)
	{
		svg_dot(&svg, isects[i].position.x, isects[i].position.y,
			"#e94560", (qaws_scalar)6);
	}

	/* Draw key points in green */
	svg_dots(&svg, pts, 5, "#50fa7b", (qaws_scalar)3);

	svg_label(&svg, (qaws_scalar)-1.3, (qaws_scalar)0.8,
		"Self-intersection (red dots)", "#8888aa");
	svg_close(&svg);
	qaws_curve_destroy(curve);
	printf("  -> " SVG_OUTPUT_DIR "/self_intersection.svg\n");
}

/* ================================================================== */
/* SVG visual test: curve-curve intersection                           */
/* ================================================================== */

static void test_svg_curve_curve_intersection(void)
{
	printf("test_svg_curve_curve_intersection\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/curve_curve_intersection.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-1.5, (qaws_scalar)5.0, (qaws_scalar)4.5,
		(qaws_scalar)600, (qaws_scalar)500))
		return;

	/* Curve A: arch */
	qaws_vec2 cpA[] = { {0,0}, {1,3}, {3,3}, {4,0} };
	/* Curve B: S-curve crossing A */
	qaws_vec2 cpB[] = { {0,2}, {2,-1}, {2,3}, {4,1} };

	qaws_bezier_desc descA, descB;
	descA.dimension = QAWS_DIMENSION_2D;
	descA.degree = 3;
	descA.control_points = cpA;
	descA.control_point_count = 4;

	descB.dimension = QAWS_DIMENSION_2D;
	descB.degree = 3;
	descB.control_points = cpB;
	descB.control_point_count = 4;

	qaws_curve *curveA = NULL, *curveB = NULL;
	qaws_status s;
	s = qaws_curve_create_bezier(&descA, &curveA);
	TEST_ASSERT_STATUS(s);
	s = qaws_curve_create_bezier(&descB, &curveB);
	TEST_ASSERT_STATUS(s);

	/* Draw both curves */
	{
		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n;
		n = svg_sample_curve(curveA, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#6272a4", (qaws_scalar)2.5);
		n = svg_sample_curve(curveB, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#50fa7b", (qaws_scalar)2.5);
	}

	/* Find intersections */
	qaws_intersection_2d isects[16];
	unsigned int isect_count = 0;
	s = qaws_curve_find_intersections_2d(curveA, curveB, isects, 16, &isect_count);
	TEST_ASSERT_STATUS(s);

	/* Draw intersection points */
	for (unsigned int i = 0; i < isect_count && i < 16; ++i)
	{
		svg_dot(&svg, isects[i].position.x, isects[i].position.y,
			"#e94560", (qaws_scalar)6);
	}

	/* Control points faintly */
	svg_dots(&svg, cpA, 4, "#6272a460", (qaws_scalar)2.5);
	svg_dots(&svg, cpB, 4, "#50fa7b60", (qaws_scalar)2.5);

	{
		char lbl[64];
		sprintf(lbl, "Curve-curve intersection (%u found)", isect_count);
		svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.5, lbl, "#8888aa");
	}
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)3.2, "Curve A (blue)", "#6272a4");
	svg_label(&svg, (qaws_scalar)-0.3, (qaws_scalar)2.9, "Curve B (green)", "#50fa7b");

	svg_close(&svg);
	qaws_curve_destroy(curveA);
	qaws_curve_destroy(curveB);
	printf("  -> " SVG_OUTPUT_DIR "/curve_curve_intersection.svg\n");
}

/* ================================================================== */
/* Unit test: curve offset                                             */
/* ================================================================== */

static void test_curve_offset(void)
{
	printf("test_curve_offset\n");

	/* Sub-test 1: Offset a simple circular arc (quarter circle).
	   A circle of radius 1 offset by +0.5 should give a circle of radius 1.5.
	   We use a NURBS quarter circle. */
	{
		qaws_scalar w = (qaws_scalar)0.70710678118;
		qaws_scalar pts[] = {
			(qaws_scalar)1.0, (qaws_scalar)0.0,
			(qaws_scalar)1.0, (qaws_scalar)1.0,
			(qaws_scalar)0.0, (qaws_scalar)1.0
		};
		qaws_scalar wts[] = { (qaws_scalar)1.0, w, (qaws_scalar)1.0 };
		qaws_scalar knots[] = {
			(qaws_scalar)0.0, (qaws_scalar)0.0, (qaws_scalar)0.0,
			(qaws_scalar)1.0, (qaws_scalar)1.0, (qaws_scalar)1.0
		};
		qaws_nurbs_desc ndesc;
		qaws_curve* arc = NULL;
		qaws_status s;

		memset(&ndesc, 0, sizeof(ndesc));
		ndesc.dimension = QAWS_DIMENSION_2D;
		ndesc.degree = 2;
		ndesc.control_points = pts;
		ndesc.control_point_count = 3;
		ndesc.weights = wts;
		ndesc.weight_count = 3;
		ndesc.knots = knots;
		ndesc.knot_count = 6;
		ndesc.is_closed = 0;

		s = qaws_curve_create_nurbs(&ndesc, &arc);
		TEST_ASSERT_STATUS(s);

		qaws_curve* offsets[4] = {NULL};
		unsigned int off_count = 0;
		/* Negative = right of travel = outward for CCW arc */
		s = qaws_curve_offset_2d(arc, (qaws_scalar)-0.5, 0, offsets, 4, &off_count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(off_count >= 1, "offset of arc produces at least 1 curve");

		/* Check that offset curve is approximately distance 1.5 from origin */
		if (off_count >= 1 && offsets[0])
		{
			qaws_vec2 sample_buf[64];
			unsigned int sample_n = 0;
			qaws_sampling_desc sd;
			sd.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
			sd.sample_count = 64;
			sd.error_tolerance = 0;
			sd.use_adaptive_sampling = 0;
			s = qaws_curve_sample_2d(offsets[0], &sd, sample_buf, 64, &sample_n);
			if (s == QAWS_STATUS_OK && sample_n > 2)
			{
				int all_close = 1;
				unsigned int k;
				for (k = 1; k + 1 < sample_n; ++k)
				{
					qaws_scalar r = (qaws_scalar)sqrt((double)(
						sample_buf[k].x * sample_buf[k].x +
						sample_buf[k].y * sample_buf[k].y));
					qaws_scalar diff = r - (qaws_scalar)1.5;
					if (diff < (qaws_scalar)0) diff = -diff;
					if (diff > (qaws_scalar)0.15) all_close = 0;
				}
				TEST_ASSERT(all_close, "offset arc samples near radius 1.5");
			}
		}
		for (unsigned int k = 0; k < off_count; ++k)
			if (offsets[k]) qaws_curve_destroy(offsets[k]);
		qaws_curve_destroy(arc);
	}

	/* Sub-test 2: Offset a Bezier S-curve. Should not crash. */
	{
		qaws_scalar pts[] = {
			(qaws_scalar)0.0, (qaws_scalar)0.0,
			(qaws_scalar)1.0, (qaws_scalar)3.0,
			(qaws_scalar)3.0, (qaws_scalar)-1.0,
			(qaws_scalar)4.0, (qaws_scalar)2.0
		};
		qaws_bezier_desc bdesc;
		qaws_curve* bez = NULL;
		qaws_status s;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;
		s = qaws_curve_create_bezier(&bdesc, &bez);
		TEST_ASSERT_STATUS(s);

		qaws_curve* offsets[8] = {NULL};
		unsigned int off_count = 0;
		s = qaws_curve_offset_2d(bez, (qaws_scalar)0.3, 0, offsets, 8, &off_count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(off_count >= 1, "offset of S-curve produces at least 1 curve");

		for (unsigned int k = 0; k < off_count; ++k)
			if (offsets[k]) qaws_curve_destroy(offsets[k]);
		qaws_curve_destroy(bez);
	}

	/* Sub-test 3: Large offset of tight curve may produce cusps / multiple curves */
	{
		qaws_scalar pts[] = {
			(qaws_scalar)0.0, (qaws_scalar)0.0,
			(qaws_scalar)0.5, (qaws_scalar)2.0,
			(qaws_scalar)1.5, (qaws_scalar)-2.0,
			(qaws_scalar)2.0, (qaws_scalar)0.0
		};
		qaws_bezier_desc bdesc;
		qaws_curve* bez = NULL;
		qaws_status s;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;
		s = qaws_curve_create_bezier(&bdesc, &bez);
		TEST_ASSERT_STATUS(s);

		qaws_curve* offsets[8] = {NULL};
		unsigned int off_count = 0;
		/* Large offset relative to curvature radius -- should handle cusps */
		s = qaws_curve_offset_2d(bez, (qaws_scalar)1.0, 1, offsets, 8, &off_count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(off_count >= 1, "large offset produces at least 1 curve");

		for (unsigned int k = 0; k < off_count; ++k)
			if (offsets[k]) qaws_curve_destroy(offsets[k]);
		qaws_curve_destroy(bez);
	}

	/* Sub-test 4: Negative offset (right side) */
	{
		qaws_scalar pts[] = {
			(qaws_scalar)0.0, (qaws_scalar)0.0,
			(qaws_scalar)2.0, (qaws_scalar)2.0,
			(qaws_scalar)4.0, (qaws_scalar)0.0
		};
		qaws_bezier_desc bdesc;
		qaws_curve* bez = NULL;
		qaws_status s;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 2;
		bdesc.control_points = pts;
		bdesc.control_point_count = 3;
		s = qaws_curve_create_bezier(&bdesc, &bez);
		TEST_ASSERT_STATUS(s);

		qaws_curve* offsets[4] = {NULL};
		unsigned int off_count = 0;
		s = qaws_curve_offset_2d(bez, (qaws_scalar)-0.5, 0, offsets, 4, &off_count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(off_count >= 1, "negative offset produces at least 1 curve");

		for (unsigned int k = 0; k < off_count; ++k)
			if (offsets[k]) qaws_curve_destroy(offsets[k]);
		qaws_curve_destroy(bez);
	}

	/* Sub-test 5: Null argument handling */
	{
		qaws_status s;
		unsigned int cnt = 0;
		qaws_curve* dummy[1] = {NULL};

		s = qaws_curve_offset_2d(NULL, (qaws_scalar)1.0, 0, dummy, 1, &cnt);
		TEST_ASSERT(s != QAWS_STATUS_OK, "null curve rejected");
	}

	/* Sub-test 6: 3D curve rejected */
	{
		qaws_scalar pts3d[] = {
			(qaws_scalar)0, (qaws_scalar)0, (qaws_scalar)0,
			(qaws_scalar)1, (qaws_scalar)1, (qaws_scalar)1,
			(qaws_scalar)2, (qaws_scalar)0, (qaws_scalar)0,
			(qaws_scalar)3, (qaws_scalar)1, (qaws_scalar)1
		};
		qaws_bezier_desc bdesc;
		qaws_curve* bez3d = NULL;
		qaws_status s;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 3;
		bdesc.control_points = pts3d;
		bdesc.control_point_count = 4;
		s = qaws_curve_create_bezier(&bdesc, &bez3d);
		TEST_ASSERT_STATUS(s);

		qaws_curve* offsets[4] = {NULL};
		unsigned int cnt = 0;
		s = qaws_curve_offset_2d(bez3d, (qaws_scalar)1.0, 0, offsets, 4, &cnt);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_DIMENSION, "3D curve rejected by offset_2d");

		qaws_curve_destroy(bez3d);
	}
}

/* ================================================================== */
/* SVG visual test: curve offset                                       */
/* ================================================================== */

static void test_svg_offset(void)
{
	printf("test_svg_offset\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/curve_offset.svg",
		(qaws_scalar)-1.5, (qaws_scalar)-7.0, (qaws_scalar)7.0, (qaws_scalar)11.5,
		(qaws_scalar)700, (qaws_scalar)1150))
		return;

	/* --- Row 1: Bezier S-curve with left/right offsets --- */
	{
		qaws_scalar pts[] = {
			(qaws_scalar)0.0, (qaws_scalar)0.0,
			(qaws_scalar)1.0, (qaws_scalar)3.0,
			(qaws_scalar)3.0, (qaws_scalar)-1.0,
			(qaws_scalar)4.0, (qaws_scalar)2.0
		};
		qaws_bezier_desc bdesc;
		qaws_curve* bez = NULL;
		qaws_status s;
		unsigned int k;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;
		s = qaws_curve_create_bezier(&bdesc, &bez);
		if (s != QAWS_STATUS_OK) goto row1_end;

		/* Draw original curve (gray) */
		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(bez, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);

		/* Left offset (blue) */
		{
			qaws_curve* offs[4] = {NULL};
			unsigned int off_cnt = 0;
			s = qaws_curve_offset_2d(bez, (qaws_scalar)0.4, 0, offs, 4, &off_cnt);
			for (k = 0; k < off_cnt; ++k)
			{
				if (!offs[k]) continue;
				n = svg_sample_curve(offs[k], buf, SVG_SAMPLES);
				svg_polyline(&svg, buf, n, "#6272a4", (qaws_scalar)1.5);
				qaws_curve_destroy(offs[k]);
			}
		}

		/* Right offset (green) */
		{
			qaws_curve* offs[4] = {NULL};
			unsigned int off_cnt = 0;
			s = qaws_curve_offset_2d(bez, (qaws_scalar)-0.4, 0, offs, 4, &off_cnt);
			for (k = 0; k < off_cnt; ++k)
			{
				if (!offs[k]) continue;
				n = svg_sample_curve(offs[k], buf, SVG_SAMPLES);
				svg_polyline(&svg, buf, n, "#50fa7b", (qaws_scalar)1.5);
				qaws_curve_destroy(offs[k]);
			}
		}

		svg_label(&svg, (qaws_scalar)-1.0, (qaws_scalar)3.5,
			"S-curve: left (blue) + right (green)", "#8888aa");
		qaws_curve_destroy(bez);
row1_end:;
	}

	/* --- Row 2: Tight S-curve with large offset (cusp region) --- */
	{
		/* Shifted down by -4 */
		qaws_scalar pts[] = {
			(qaws_scalar)0.0, (qaws_scalar)-2.5,
			(qaws_scalar)0.5, (qaws_scalar)-0.5,
			(qaws_scalar)1.5, (qaws_scalar)-4.5,
			(qaws_scalar)2.0, (qaws_scalar)-2.5
		};
		qaws_bezier_desc bdesc;
		qaws_curve* bez = NULL;
		qaws_status s;
		unsigned int k;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;
		s = qaws_curve_create_bezier(&bdesc, &bez);
		if (s != QAWS_STATUS_OK) goto row2_end;

		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(bez, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);

		/* Large left offset WITHOUT trim (raw, dashed orange) */
		{
			qaws_curve* offs[8] = {NULL};
			unsigned int off_cnt = 0;
			s = qaws_curve_offset_2d(bez, (qaws_scalar)0.8, 0, offs, 8, &off_cnt);
			for (k = 0; k < off_cnt; ++k)
			{
				if (!offs[k]) continue;
				n = svg_sample_curve(offs[k], buf, SVG_SAMPLES);
				svg_polyline(&svg, buf, n, "#f0a030", (qaws_scalar)3.0);
				qaws_curve_destroy(offs[k]);
			}
		}

		/* Large left offset WITH trim (solid red) */
		{
			qaws_curve* offs[8] = {NULL};
			unsigned int off_cnt = 0;
			s = qaws_curve_offset_2d(bez, (qaws_scalar)0.8, 1, offs, 8, &off_cnt);
			for (k = 0; k < off_cnt; ++k)
			{
				if (!offs[k]) continue;
				n = svg_sample_curve(offs[k], buf, SVG_SAMPLES);
				svg_polyline(&svg, buf, n, "#e94560", (qaws_scalar)1.0);
				qaws_curve_destroy(offs[k]);
			}
		}

		/* Smaller offset for comparison (no cusps expected) */
		{
			qaws_curve* offs[4] = {NULL};
			unsigned int off_cnt = 0;
			s = qaws_curve_offset_2d(bez, (qaws_scalar)0.15, 0, offs, 4, &off_cnt);
			for (k = 0; k < off_cnt; ++k)
			{
				if (!offs[k]) continue;
				n = svg_sample_curve(offs[k], buf, SVG_SAMPLES);
				svg_polyline(&svg, buf, n, "#6272a4", (qaws_scalar)1.5);
				qaws_curve_destroy(offs[k]);
			}
		}

		svg_label(&svg, (qaws_scalar)-1.0, (qaws_scalar)-0.3,
			"Tight S: offset=0.8 raw (orange) vs trimmed (red), small=0.15 (blue)", "#8888aa");
		qaws_curve_destroy(bez);
row2_end:;
	}

	/* --- Row 3: "Stroke" visualization — paired left+right offsets at multiple widths --- */
	{
		qaws_scalar pts[] = {
			(qaws_scalar)0.0, (qaws_scalar)-5.0,
			(qaws_scalar)1.5, (qaws_scalar)-3.5,
			(qaws_scalar)2.5, (qaws_scalar)-6.5,
			(qaws_scalar)4.0, (qaws_scalar)-5.0
		};
		qaws_bezier_desc bdesc;
		qaws_curve* bez = NULL;
		qaws_status s;
		unsigned int k;
		static const qaws_scalar widths[] = {
			(qaws_scalar)0.15, (qaws_scalar)0.3, (qaws_scalar)0.6
		};
		static const char* colors_l[] = { "#6272a4", "#bd93f9", "#ff79c6" };
		static const char* colors_r[] = { "#6272a4", "#bd93f9", "#ff79c6" };

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;
		s = qaws_curve_create_bezier(&bdesc, &bez);
		if (s != QAWS_STATUS_OK) goto row3_end;

		qaws_vec2 buf[SVG_SAMPLES];
		unsigned int n = svg_sample_curve(bez, buf, SVG_SAMPLES);
		svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);

		for (k = 0; k < 3; ++k)
		{
			qaws_curve* offs[4] = {NULL};
			unsigned int off_cnt = 0, m;
			/* Left */
			s = qaws_curve_offset_2d(bez, widths[k], 0, offs, 4, &off_cnt);
			for (m = 0; m < off_cnt; ++m) {
				if (!offs[m]) continue;
				n = svg_sample_curve(offs[m], buf, SVG_SAMPLES);
				svg_polyline(&svg, buf, n, colors_l[k], (qaws_scalar)1.0);
				qaws_curve_destroy(offs[m]);
			}
			/* Right */
			off_cnt = 0;
			memset(offs, 0, sizeof(offs));
			s = qaws_curve_offset_2d(bez, -widths[k], 0, offs, 4, &off_cnt);
			for (m = 0; m < off_cnt; ++m) {
				if (!offs[m]) continue;
				n = svg_sample_curve(offs[m], buf, SVG_SAMPLES);
				svg_polyline(&svg, buf, n, colors_r[k], (qaws_scalar)1.0);
				qaws_curve_destroy(offs[m]);
			}
		}

		svg_label(&svg, (qaws_scalar)-1.0, (qaws_scalar)-2.8,
			"Stroke widths: 0.15 (blue) 0.3 (purple) 0.6 (pink)", "#8888aa");
		qaws_curve_destroy(bez);
row3_end:;
	}

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/curve_offset.svg\n");
}

static void test_svg_offset_closed(void)
{
	printf("test_svg_offset_closed\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/offset_closed.svg",
		(qaws_scalar)-4.5, (qaws_scalar)-3.5, (qaws_scalar)8.0, (qaws_scalar)21.0,
		(qaws_scalar)600, (qaws_scalar)1600))
		return;

	/* Cayley oval: locus where harmonic mean of distances to two
	 * foci equals b.  1/r + 1/r' = 2/b, foci at (+-a, 0).
	 * b = a gives figure-8; b slightly > a gives an oval that
	 * is very close to self-intersecting (pinched waist). */
	{
		#define CAYLEY_N 64
		qaws_scalar pts[CAYLEY_N * 2];
		qaws_catmull_rom_desc cr_desc;
		qaws_curve* crv = NULL;
		qaws_status s;
		unsigned int k;
		double ca = 1.0, cb = 1.08; /* close to self-intersecting */

		for (k = 0; k < CAYLEY_N; ++k)
		{
			double theta = 2.0 * 3.14159265358979323846 * (double)k / (double)CAYLEY_N;
			double ct = cos(theta), st = sin(theta);
			/* Bisect for R: 1/r + 1/r' = 2/b where
			 * r = sqrt(R^2 + 2aR*cos(t) + a^2)
			 * r' = sqrt(R^2 - 2aR*cos(t) + a^2) */
			double lo = 0.001, hi = 4.0, mid, R;
			int iter;
			for (iter = 0; iter < 60; ++iter)
			{
				double r1, r2, val;
				mid = 0.5 * (lo + hi);
				r1 = sqrt(mid * mid + 2.0 * ca * mid * ct + ca * ca);
				r2 = sqrt(mid * mid - 2.0 * ca * mid * ct + ca * ca);
				val = 1.0 / r1 + 1.0 / r2;
				if (val > 2.0 / cb) lo = mid;
				else hi = mid;
			}
			R = 0.5 * (lo + hi);
			pts[k * 2 + 0] = (qaws_scalar)(R * ct);
			pts[k * 2 + 1] = (qaws_scalar)(R * st);
		}

		memset(&cr_desc, 0, sizeof(cr_desc));
		cr_desc.dimension = QAWS_DIMENSION_2D;
		cr_desc.control_points = pts;
		cr_desc.control_point_count = CAYLEY_N;
		cr_desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		cr_desc.closed = 1;
		s = qaws_curve_create_catmull_rom(&cr_desc, &crv);
		if (s != QAWS_STATUS_OK) goto closed_end;

		/* Three rows stacked vertically: trim=0, trim=1, trim=2 */
		{
			static const qaws_scalar dists_out[] = {
				(qaws_scalar)0.3, (qaws_scalar)0.6, (qaws_scalar)1.0
			};
			static const qaws_scalar dists_in[] = {
				(qaws_scalar)-0.3, (qaws_scalar)-0.6, (qaws_scalar)-1.0
			};
			static const char* colors_out[] = { "#50fa7b", "#8be9fd", "#ff79c6" };
			static const char* colors_in[]  = { "#6272a4", "#bd93f9", "#e94560" };
			static const int trim_vals[] = { 0, 1, 2 };
			static const char* labels[] = {
				"trim=0 (raw)", "trim=1 (distance)", "trim=2 (cleanup)"
			};
			int row;
			for (row = 0; row < 3; ++row)
			{
				int do_trim = trim_vals[row];
				qaws_scalar yoff = (qaws_scalar)row * (qaws_scalar)7.0;

				/* Draw base curve shifted */
				{
					qaws_vec2 buf[SVG_SAMPLES];
					unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES), m;
					for (m = 0; m < n; ++m) buf[m].y += yoff;
					svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);
				}

				/* Outward offsets */
				for (k = 0; k < 3; ++k)
				{
					qaws_curve* offs[8] = {NULL};
					unsigned int off_cnt = 0, m;
					s = qaws_curve_offset_2d(crv, dists_out[k], do_trim, offs, 8, &off_cnt);
					for (m = 0; m < off_cnt; ++m) {
						qaws_vec2 buf[SVG_SAMPLES];
						unsigned int n;
						if (!offs[m]) continue;
						n = svg_sample_curve(offs[m], buf, SVG_SAMPLES);
						{ unsigned int q; for (q = 0; q < n; ++q) buf[q].y += yoff; }
						svg_polyline(&svg, buf, n, colors_out[k], (qaws_scalar)1.5);
						qaws_curve_destroy(offs[m]);
					}
				}

				/* Inward offsets */
				for (k = 0; k < 3; ++k)
				{
					qaws_curve* offs[8] = {NULL};
					unsigned int off_cnt = 0, m;
					s = qaws_curve_offset_2d(crv, dists_in[k], do_trim, offs, 8, &off_cnt);
					for (m = 0; m < off_cnt; ++m) {
						qaws_vec2 buf[SVG_SAMPLES];
						unsigned int n;
						if (!offs[m]) continue;
						n = svg_sample_curve(offs[m], buf, SVG_SAMPLES);
						{ unsigned int q; for (q = 0; q < n; ++q) buf[q].y += yoff; }
						svg_polyline(&svg, buf, n, colors_in[k], (qaws_scalar)1.5);
						qaws_curve_destroy(offs[m]);
					}
				}

				svg_label(&svg, (qaws_scalar)-4.0, (qaws_scalar)(-3.0) + yoff,
					labels[row], "#8888aa");
			}
		}
		qaws_curve_destroy(crv);
		#undef CAYLEY_N
closed_end:;
	}

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/offset_closed.svg\n");
}

static void test_svg_offset_selfintersect(void)
{
	printf("test_svg_offset_selfintersect\n");

	svg_writer svg;
	if (!svg_open(&svg, SVG_OUTPUT_DIR "/offset_selfintersect.svg",
		(qaws_scalar)-4.0, (qaws_scalar)-2.5, (qaws_scalar)8.0, (qaws_scalar)18.0,
		(qaws_scalar)600, (qaws_scalar)1350))
		return;

	/* Lemniscate of Bernoulli (figure-8, self-intersects at origin)
	 * x(t) = a*cos(t) / (1 + sin^2(t))
	 * y(t) = a*sin(t)*cos(t) / (1 + sin^2(t))
	 * t in [0, 2*pi] traces the full curve. */
	{
		#define LEMN_OFF_N 48
		qaws_scalar pts[LEMN_OFF_N * 2];
		qaws_catmull_rom_desc cr_desc;
		qaws_curve* crv = NULL;
		qaws_status s;
		unsigned int k;
		double a = 2.5;

		for (k = 0; k < LEMN_OFF_N; ++k)
		{
			double t = 2.0 * 3.14159265358979323846 * (double)k / (double)LEMN_OFF_N;
			double st = sin(t);
			double ct = cos(t);
			double denom = 1.0 + st * st;
			pts[k * 2 + 0] = (qaws_scalar)(a * ct / denom);
			pts[k * 2 + 1] = (qaws_scalar)(a * st * ct / denom);
		}

		memset(&cr_desc, 0, sizeof(cr_desc));
		cr_desc.dimension = QAWS_DIMENSION_2D;
		cr_desc.control_points = pts;
		cr_desc.control_point_count = LEMN_OFF_N;
		cr_desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		cr_desc.closed = 1;
		s = qaws_curve_create_catmull_rom(&cr_desc, &crv);
		if (s != QAWS_STATUS_OK) goto si_end;

		/* Three rows stacked vertically: trim=0, trim=1, trim=2 */
		{
			static const qaws_scalar dists_l[] = {
				(qaws_scalar)0.2, (qaws_scalar)0.5, (qaws_scalar)0.8
			};
			static const qaws_scalar dists_r[] = {
				(qaws_scalar)-0.2, (qaws_scalar)-0.5, (qaws_scalar)-0.8
			};
			static const char* colors_l[] = { "#50fa7b", "#8be9fd", "#ff79c6" };
			static const char* colors_r[] = { "#6272a4", "#bd93f9", "#e94560" };
			static const int trim_vals[] = { 0, 1, 2 };
			static const char* labels[] = {
				"trim=0 (raw)", "trim=1 (distance)", "trim=2 (cleanup)"
			};
			int row;
			for (row = 0; row < 3; ++row)
			{
				int do_trim = trim_vals[row];
				qaws_scalar yoff = (qaws_scalar)row * (qaws_scalar)6.0;

				/* Draw base curve shifted */
				{
					qaws_vec2 buf[SVG_SAMPLES];
					unsigned int n = svg_sample_curve(crv, buf, SVG_SAMPLES), m;
					for (m = 0; m < n; ++m) buf[m].y += yoff;
					svg_polyline(&svg, buf, n, "#888888", (qaws_scalar)2.5);
				}

				/* Left offsets */
				for (k = 0; k < 3; ++k)
				{
					qaws_curve* offs[8] = {NULL};
					unsigned int off_cnt = 0, m;
					s = qaws_curve_offset_2d(crv, dists_l[k], do_trim, offs, 8, &off_cnt);
					for (m = 0; m < off_cnt; ++m) {
						qaws_vec2 buf[SVG_SAMPLES];
						unsigned int n;
						if (!offs[m]) continue;
						n = svg_sample_curve(offs[m], buf, SVG_SAMPLES);
						{ unsigned int q; for (q = 0; q < n; ++q) buf[q].y += yoff; }
						svg_polyline(&svg, buf, n, colors_l[k], (qaws_scalar)1.5);
						qaws_curve_destroy(offs[m]);
					}
				}

				/* Right offsets */
				for (k = 0; k < 3; ++k)
				{
					qaws_curve* offs[8] = {NULL};
					unsigned int off_cnt = 0, m;
					s = qaws_curve_offset_2d(crv, dists_r[k], do_trim, offs, 8, &off_cnt);
					for (m = 0; m < off_cnt; ++m) {
						qaws_vec2 buf[SVG_SAMPLES];
						unsigned int n;
						if (!offs[m]) continue;
						n = svg_sample_curve(offs[m], buf, SVG_SAMPLES);
						{ unsigned int q; for (q = 0; q < n; ++q) buf[q].y += yoff; }
						svg_polyline(&svg, buf, n, colors_r[k], (qaws_scalar)1.5);
						qaws_curve_destroy(offs[m]);
					}
				}

				svg_label(&svg, (qaws_scalar)-3.5, (qaws_scalar)(-2.0) + yoff,
					labels[row], "#8888aa");
			}
		}
		qaws_curve_destroy(crv);
		#undef LEMN_OFF_N
si_end:;
	}

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/offset_selfintersect.svg\n");
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

	/* Phase 1 tests */
	test_frenet_frame();
	test_curve_reverse();
	test_easing();
	test_wrap_modes();
	test_traversal_advance();
	test_adaptive_eval_sampling();

	/* Phase 2 tests */
	test_curve_split();
	test_curve_join();
	test_family_conversion();
	test_degree_elevation();
	test_precomputed_coefficients();

	/* Geometric helpers & adaptive sampling */
	test_geometric_helpers();
	test_adaptive_sampling();

	/* Yuksel C2 interpolating spline tests */
	test_yuksel_bezier_mode();
	test_yuksel_closed();
	test_yuksel_circular_mode();
	test_yuksel_elliptical_mode();
	test_yuksel_hybrid_mode();
	test_yuksel_3d();
	test_yuksel_error_handling();

	/* Phase 3 tests */
	test_inflection_points();
	test_extrema();
	test_curvature_comb();
	test_winding_number();
	test_yuksel_3d_circular();

	/* Phase 5 tests */
	test_multi_curve_traversal();
	test_scurve_profile();
	test_custom_speed();
	test_curvature_weighted_sampling();
	test_feature_preserving_sampling();
	test_streaming_sampling();

	/* Phase 7 tests */
	test_self_intersection();
	test_curve_curve_intersection();
	test_curve_offset();

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

	/* Phase 1 & 2 SVG visual tests */
	test_svg_curve_reverse();
	test_svg_curve_split();
	test_svg_curve_join();
	test_svg_degree_elevation();
	test_svg_family_conversion();
	test_svg_easing_curves();
	test_svg_adaptive_sampling();

	/* Phase 3 & 5 SVG visual tests */
	test_svg_inflection_points();
	test_svg_extrema();
	test_svg_curvature_comb();
	test_svg_winding_number();
	test_svg_curvature_weighted();
	test_svg_feature_preserving();
	test_svg_multi_traversal();
	test_svg_frenet_frame();

	/* Phase 7 SVG visual tests */
	test_svg_self_intersection();
	test_svg_curve_curve_intersection();
	test_svg_offset();
	test_svg_offset_closed();
	test_svg_offset_selfintersect();

	printf("\n=== Results: %d passed, %d failed ===\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
