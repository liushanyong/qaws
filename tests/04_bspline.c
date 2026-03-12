#include "test_common.h"

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

int test_04_bspline_main(void)
{
	g_pass = 0; g_fail = 0;
	test_bspline();
	test_bspline_custom_knots();
	return g_fail > 0 ? 1 : 0;
}
