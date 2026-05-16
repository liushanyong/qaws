/*
 * eval_arc.hlsl — Arc evaluation compute shader (HLSL / Direct3D 12)
 *
 * Evaluates a 3-D circular arc at a uniform set of parameter values.
 * Uses the qaws arc evaluation core (qaws_arc_core.h).
 *
 * Buffer layout:
 *   ArcParams (t0): flat float array — center[3], axis_u[3], axis_v[3]
 *   Results   (u0): output XYZ float array
 *
 * Constant buffer b0:
 *   theta_start  — arc start angle in radians
 *   theta_range  — angular sweep in radians (signed)
 *   sample_count — number of output points
 */

#include "../src/qaws/qaws_hlsl.h"
#include "../src/qaws/core/qaws_arc_core.h"

StructuredBuffer<float>   ArcParams : register(t0);
RWStructuredBuffer<float> Results   : register(u0);

cbuffer Params : register(b0) {
    float theta_start;
    float theta_range;
    uint  sample_count;
};

[numthreads(256, 1, 1)]
void CSMain(uint3 id : SV_DispatchThreadID)
{
    if (id.x >= sample_count)
        return;

    qaws_vec3 center = { ArcParams[0], ArcParams[1], ArcParams[2] };
    qaws_vec3 axis_u = { ArcParams[3], ArcParams[4], ArcParams[5] };
    qaws_vec3 axis_v = { ArcParams[6], ArcParams[7], ArcParams[8] };

    float t = (float)id.x / (float)(sample_count - 1);

    qaws_vec3 pos = qaws_arc_eval_3d(center, axis_u, axis_v,
                                     theta_start, theta_range, t);

    Results[id.x * 3 + 0] = pos.x;
    Results[id.x * 3 + 1] = pos.y;
    Results[id.x * 3 + 2] = pos.z;
}
