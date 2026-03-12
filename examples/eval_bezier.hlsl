#include "../src/qaws/qaws_hlsl.h"
#include "../src/qaws/core/qaws_decasteljau_core.h"

StructuredBuffer<float> cp : register(t0);
RWStructuredBuffer<float> results : register(u0);

cbuffer Params : register(b0) {
    uint sample_count;
};

[numthreads(256, 1, 1)]
void CSMain(uint3 id : SV_DispatchThreadID)
{
    if (id.x >= sample_count)
        return;

    float t = (float)id.x / (float)(sample_count - 1);

    qaws_vec2 p0 = { cp[0], cp[1] };
    qaws_vec2 p1 = { cp[2], cp[3] };
    qaws_vec2 p2 = { cp[4], cp[5] };
    qaws_vec2 p3 = { cp[6], cp[7] };

    qaws_vec2 pos = qaws_bezier3_eval_2d(p0, p1, p2, p3, t);

    results[id.x * 2]     = pos.x;
    results[id.x * 2 + 1] = pos.y;
}
