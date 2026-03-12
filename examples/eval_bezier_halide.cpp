#include "../src/qaws/qaws_halide.h"
#include "../src/qaws/core/qaws_decasteljau_core.h"

#include <Halide.h>
#include <cstdio>

int main()
{
    /* Control points as ImageParam (8 floats: 4 CPs x 2D) */
    Halide::ImageParam cp_buf(Halide::Float(32), 1, "cp");

    /* Evaluation parameters */
    Halide::Param<float> t_start("t_start");
    Halide::Param<float> t_step("t_step");

    Halide::Func bezier("bezier");
    Halide::Var i("i");

    /* Build control points from buffer */
    qaws_vec2 p0 = { cp_buf(0), cp_buf(1) };
    qaws_vec2 p1 = { cp_buf(2), cp_buf(3) };
    qaws_vec2 p2 = { cp_buf(4), cp_buf(5) };
    qaws_vec2 p3 = { cp_buf(6), cp_buf(7) };

    /* Parametric evaluation: t = t_start + i * t_step */
    qaws_scalar t = t_start + Halide::cast<float>(i) * t_step;

    /* Evaluate cubic Bezier */
    qaws_vec2 pos = qaws_bezier3_eval_2d(p0, p1, p2, p3, t);

    /* Output as tuple (x, y) */
    bezier(i) = Halide::Tuple(pos.x, pos.y);

    /* Schedule: vectorize across samples */
    bezier.vectorize(i, 8);

    printf("Halide Bezier generator built successfully.\n");
    return 0;
}
