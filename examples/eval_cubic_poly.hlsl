#include "../src/qaws/qaws_hlsl.h"
#include "../src/qaws/core/qaws_cubic_poly_core.h"

StructuredBuffer<float> span_coeffs : register(t0);
RWTexture1D<float2> results : register(u0);

cbuffer Params : register(b0) {
    uint span_count;
    uint sample_count;
};

[numthreads(256, 1, 1)]
void CSMain(uint3 id : SV_DispatchThreadID)
{
    if (id.x >= sample_count)
        return;

    float global_t = (float)id.x / (float)(sample_count - 1) * (float)span_count;
    uint span = min((uint)global_t, span_count - 1);
    float local_t = global_t - (float)span;

    /* Load pre-computed a,b,c,d coefficients for this span.
       Layout: span_count * 8 floats (2D), stored as [ax,bx,cx,dx, ay,by,cy,dy] */
    uint base = span * 8;
    qaws_vec2 a = { span_coeffs[base + 0], span_coeffs[base + 4] };
    qaws_vec2 b = { span_coeffs[base + 1], span_coeffs[base + 5] };
    qaws_vec2 c = { span_coeffs[base + 2], span_coeffs[base + 6] };
    qaws_vec2 d = { span_coeffs[base + 3], span_coeffs[base + 7] };

    qaws_eval_2d result = qaws_cubic_eval_2d(a, b, c, d, local_t, 1);
    results[id.x] = float2(result.position.x, result.position.y);
}
