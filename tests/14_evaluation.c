#include "test_common.h"

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

int test_14_evaluation_main(void)
{
	g_pass = 0; g_fail = 0;

	test_arc_length();
	test_inspection_flags();
	test_bounds_3d();
	test_precomputed_coefficients();

	printf("14_evaluation: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
