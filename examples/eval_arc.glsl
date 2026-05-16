#version 450
#extension GL_GOOGLE_include_directive : enable

/*
 * eval_arc.glsl — Arc evaluation compute shader
 *
 * Evaluates a 3-D circular arc at a uniform set of parameter values
 * and writes the resulting XYZ points to an output buffer.
 *
 * Uses the qaws arc evaluation core (qaws_arc_core.h) which provides the
 * same formula as the C runtime but as an inline GLSL function.
 *
 * Buffer layout:
 *   Params (push constant):
 *     theta_start  — arc start angle in radians
 *     theta_range  — arc angular sweep in radians (signed)
 *     sample_count — number of output points
 *
 *   Binding 0 (input)  : center[3] followed by axis_u[3] and axis_v[3]
 *                         (radius-scaled, like qaws_arc_segment)
 *   Binding 1 (output) : results[], flat XYZ float array
 */

#include "../src/qaws/qaws_glsl.h"
#include "../src/qaws/core/qaws_arc_core.h"

layout(std430, binding = 0) readonly buffer ArcParams {
    /* [0..2]  center    xyz  */
    /* [3..5]  axis_u    xyz  (radius-scaled cos direction) */
    /* [6..8]  axis_v    xyz  (radius-scaled sin direction) */
    float arc_params[];
};

layout(std430, binding = 1) writeonly buffer OutBuffer {
    float results[];
};

layout(push_constant) uniform Params {
    float theta_start;
    float theta_range;
    uint  sample_count;
};

layout(local_size_x = 256) in;

void main()
{
    uint idx = gl_GlobalInvocationID.x;
    if (idx >= sample_count)
        return;

    qaws_vec3 center = qaws_vec3(arc_params[0], arc_params[1], arc_params[2]);
    qaws_vec3 axis_u = qaws_vec3(arc_params[3], arc_params[4], arc_params[5]);
    qaws_vec3 axis_v = qaws_vec3(arc_params[6], arc_params[7], arc_params[8]);

    float t = float(idx) / float(sample_count - 1u);

    /* qaws_arc_eval_3d is inlined from qaws_arc_core.h */
    qaws_vec3 pos = qaws_arc_eval_3d(center, axis_u, axis_v,
                                     theta_start, theta_range, t);

    results[idx * 3u + 0u] = pos.x;
    results[idx * 3u + 1u] = pos.y;
    results[idx * 3u + 2u] = pos.z;
}
