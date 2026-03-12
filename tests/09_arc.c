#include "test_common.h"

static void test_arc(void)
{
	printf("test_arc\n");

	/* Single quarter-circle arc */
	{
		qaws_arc_segment seg;
		qaws_arc_desc desc;
		qaws_curve *crv = NULL;
		qaws_eval_result_2d r;
		qaws_scalar dist;
		qaws_range rng;

		memset(&seg, 0, sizeof(seg));
		seg.center[0] = 0; seg.center[1] = 0;
		seg.radius = 1;
		seg.angle_start = 0;
		seg.angle_end = (qaws_scalar)1.5707963268; /* pi/2 */

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.segments = &seg;
		desc.segment_count = 1;

		TEST_ASSERT_STATUS(qaws_curve_create_arc(&desc, &crv));
		TEST_ASSERT(crv != NULL, "arc created");
		TEST_ASSERT(qaws_curve_get_kind(crv) == QAWS_CURVE_KIND_ARC, "arc kind");

		/* Arc length should be pi/2 */
		rng = qaws_curve_get_parameter_range(crv);
		TEST_ASSERT(approx_eq(rng.max_value, (qaws_scalar)1.5707963268),
			"arc parameter range = arc length");

		/* t=0 (local_t=0): (1, 0) */
		qaws_curve_evaluate_2d(crv, 0, QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, 1) && approx_eq(r.position.y, 0),
			"arc start point");

		/* t=max (local_t=1): (0, 1) */
		qaws_curve_evaluate_2d(crv, rng.max_value, QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, 0) && approx_eq(r.position.y, 1),
			"arc end point");

		/* Midpoint should be on unit circle */
		qaws_curve_evaluate_2d(crv, rng.max_value * (qaws_scalar)0.5,
			QAWS_EVAL_FLAG_POSITION, &r);
		dist = (qaws_scalar)sqrt((double)(r.position.x * r.position.x + r.position.y * r.position.y));
		TEST_ASSERT(approx_eq(dist, 1), "arc midpoint on circle");

		/* All sampled points on unit circle */
		{
			unsigned int i;
			int all_on = 1;
			for (i = 0; i <= 20; ++i)
			{
				qaws_scalar t = rng.max_value * (qaws_scalar)i / (qaws_scalar)20;
				qaws_curve_evaluate_2d(crv, t, QAWS_EVAL_FLAG_POSITION, &r);
				dist = (qaws_scalar)sqrt((double)(r.position.x * r.position.x + r.position.y * r.position.y));
				if (!approx_eq(dist, 1)) all_on = 0;
			}
			TEST_ASSERT(all_on, "arc all samples on unit circle");
		}

		/* Derivatives at start: D1 should be tangent (0, sweep*r) direction */
		qaws_curve_evaluate_2d(crv, 0,
			QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3, &r);
		TEST_ASSERT(r.valid_flags & QAWS_EVAL_FLAG_D1, "arc d1 valid");
		TEST_ASSERT(r.valid_flags & QAWS_EVAL_FLAG_D2, "arc d2 valid");
		TEST_ASSERT(r.valid_flags & QAWS_EVAL_FLAG_D3, "arc d3 valid");
		/* At theta=0 on unit circle: D1 direction is (0, 1) * sweep * r */
		TEST_ASSERT(approx_eq(r.d1.x, 0), "arc d1.x at start");
		TEST_ASSERT(r.d1.y > 0, "arc d1.y positive at start");
		/* D2 points inward: (-1, 0) direction * sweep^2 * r */
		TEST_ASSERT(r.d2.x < 0, "arc d2.x negative at start");

		/* is_rational should return 1 */
		TEST_ASSERT(qaws_curve_is_rational(crv) == 1, "arc is_rational");

		qaws_curve_destroy(crv);
	}

	/* Multi-segment: two quarter arcs forming a half circle */
	{
		qaws_scalar pi_half = (qaws_scalar)1.5707963268;
		qaws_arc_segment segs[2];
		qaws_arc_desc desc;
		qaws_curve *crv = NULL;
		qaws_range rng;
		qaws_eval_result_2d r;

		memset(segs, 0, sizeof(segs));
		segs[0].center[0] = 0; segs[0].center[1] = 0;
		segs[0].radius = 1;
		segs[0].angle_start = 0;
		segs[0].angle_end = pi_half;

		segs[1].center[0] = 0; segs[1].center[1] = 0;
		segs[1].radius = 1;
		segs[1].angle_start = pi_half;
		segs[1].angle_end = pi_half * 2;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.segments = segs;
		desc.segment_count = 2;

		TEST_ASSERT_STATUS(qaws_curve_create_arc(&desc, &crv));
		TEST_ASSERT(qaws_curve_get_span_count(crv) == 2, "arc 2-span count");

		rng = qaws_curve_get_parameter_range(crv);
		/* Total arc length = pi */
		TEST_ASSERT(approx_eq(rng.max_value, pi_half * 2), "arc half-circle length");

		/* End point should be (-1, 0) */
		qaws_curve_evaluate_2d(crv, rng.max_value, QAWS_EVAL_FLAG_POSITION, &r);
		TEST_ASSERT(approx_eq(r.position.x, -1), "arc half end x");
		TEST_ASSERT(approx_eq(r.position.y, 0), "arc half end y");

		qaws_curve_destroy(crv);
	}

	/* Error: zero radius */
	{
		qaws_arc_segment seg;
		qaws_arc_desc desc;
		qaws_curve *crv = NULL;
		memset(&seg, 0, sizeof(seg));
		seg.radius = 0;
		seg.angle_start = 0;
		seg.angle_end = 1;
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.segments = &seg;
		desc.segment_count = 1;
		TEST_ASSERT(qaws_curve_create_arc(&desc, &crv) != QAWS_STATUS_OK,
			"arc zero radius rejected");
	}

	/* 3D arc */
	{
		qaws_arc_segment seg3d;
		qaws_arc_desc desc3d;
		qaws_curve *crv3d = NULL;
		qaws_eval_result_3d r3d;
		qaws_scalar dist;

		memset(&seg3d, 0, sizeof(seg3d));
		seg3d.center[0] = 0; seg3d.center[1] = 0; seg3d.center[2] = 0;
		seg3d.radius = 1;
		seg3d.angle_start = 0;
		seg3d.angle_end = (qaws_scalar)1.5707963268;
		seg3d.axis_u[0] = 1; seg3d.axis_u[1] = 0; seg3d.axis_u[2] = 0;
		seg3d.axis_v[0] = 0; seg3d.axis_v[1] = 0; seg3d.axis_v[2] = 1;

		memset(&desc3d, 0, sizeof(desc3d));
		desc3d.dimension = QAWS_DIMENSION_3D;
		desc3d.segments = &seg3d;
		desc3d.segment_count = 1;

		TEST_ASSERT_STATUS(qaws_curve_create_arc(&desc3d, &crv3d));
		qaws_curve_evaluate_3d(crv3d, 0, QAWS_EVAL_FLAG_POSITION, &r3d);
		TEST_ASSERT(approx_eq(r3d.position.x, 1) && approx_eq(r3d.position.z, 0),
			"3d arc start");

		/* End: (0, 0, 1) in XZ plane */
		{
			qaws_range rng = qaws_curve_get_parameter_range(crv3d);
			qaws_curve_evaluate_3d(crv3d, rng.max_value,
				QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &r3d);
			TEST_ASSERT(approx_eq(r3d.position.x, 0), "3d arc end x");
			TEST_ASSERT(approx_eq(r3d.position.z, 1), "3d arc end z");
			TEST_ASSERT(approx_eq(r3d.position.y, 0), "3d arc y=0 plane");
			/* 3D arc on circle: all points at unit distance from center */
			dist = (qaws_scalar)sqrt((double)(
				r3d.position.x * r3d.position.x +
				r3d.position.y * r3d.position.y +
				r3d.position.z * r3d.position.z));
			TEST_ASSERT(approx_eq(dist, 1), "3d arc end on unit sphere");
		}

		qaws_curve_destroy(crv3d);
	}
}

int test_09_arc_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_arc();

	printf("09_arc: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
