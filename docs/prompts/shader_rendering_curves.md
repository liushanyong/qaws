# 提示词文档：使用 Qaws 实现曲线的着色器渲染

## 背景

**Qaws**（قوس，阿拉伯语"弧/弓/曲线"）是一个无依赖的 C11 参数曲线与曲面库，支持四种后端：
**C · HLSL · GLSL · Halide**。

核心求值数学代码位于 `src/qaws/core/` 下的纯头文件中，同一份源码可在所有后端编译，
因此同样的曲线求值逻辑可以无缝用于 GPU 计算着色器。

---

## 核心概念

### 后端选择

| 后端    | 宏值 | 适用场景                  | 标量类型          |
|---------|:----:|---------------------------|-------------------|
| C       | 0    | CPU 运行时、工具、SIMD批处理 | `float` / `double` |
| HLSL    | 1    | DirectX 计算/像素着色器    | `float`           |
| GLSL    | 2    | Vulkan / OpenGL 计算着色器 | `float`           |
| Halide  | 3    | AOT 图像处理管线           | `Halide::Expr`    |

后端由编译器宏自动检测（`__HLSL_VERSION`、`GL_core_profile`、`HALIDE_HALIDERUNTIME_H`），
也可通过引导头文件强制指定。

### 引导头文件

```c
// 在着色器文件顶部包含对应引导头，自动设置 QAWS_BACKEND
#include "qaws_hlsl.h"   // DirectX HLSL
#include "qaws_glsl.h"   // Vulkan / OpenGL GLSL
```

---

## 着色器渲染完整流程

### 第一步：CPU 端准备数据

在 C/C++ 宿主代码中调用 `prepare` 函数，将控制点打包成 GPU 友好的平坦数组，
再通过 D3D12 / Vulkan API 上传到 GPU Buffer：

```c
#include "qaws.h"
#include "qaws_prepare.h"

// Hermite 曲线：生成 a/b/c/d 系数（每段 4 组值，每组按维度展开）
float a[N*2], b[N*2], c[N*2], d[N*2];
qaws_hermite_prepare_2d(points, tangents, count, a, b, c, d);

// 通过图形 API 将 a,b,c,d 上传为 StructuredBuffer（HLSL）或 SSBO（GLSL）
```

> **提示**：对于 Bézier，控制点可直接上传，无需 prepare 步骤。

---

### 第二步：编写着色器

#### GLSL 示例（Vulkan Compute Shader）— Cubic Bézier

```glsl
#version 450
#extension GL_GOOGLE_include_directive : enable

#include "../src/qaws/qaws_glsl.h"
#include "../src/qaws/core/qaws_decasteljau_core.h"

layout(std430, binding = 0) readonly buffer CPBuffer { float cp[]; };
layout(std430, binding = 1) buffer OutBuffer        { float results[]; };

layout(push_constant) uniform Params {
    uint sample_count;
};

layout(local_size_x = 256) in;
void main()
{
    uint idx = gl_GlobalInvocationID.x;
    if (idx >= sample_count) return;

    float t = float(idx) / float(sample_count - 1u);

    qaws_vec2 p0 = qaws_vec2(cp[0], cp[1]);
    qaws_vec2 p1 = qaws_vec2(cp[2], cp[3]);
    qaws_vec2 p2 = qaws_vec2(cp[4], cp[5]);
    qaws_vec2 p3 = qaws_vec2(cp[6], cp[7]);

    qaws_vec2 pos = qaws_bezier3_eval_2d(p0, p1, p2, p3, t);

    results[idx * 2u]      = pos.x;
    results[idx * 2u + 1u] = pos.y;
}
```

#### HLSL 示例（DirectX Compute Shader）— Cubic Bézier

```hlsl
#include "../src/qaws/qaws_hlsl.h"
#include "../src/qaws/core/qaws_decasteljau_core.h"

StructuredBuffer<float>   cp      : register(t0);
RWStructuredBuffer<float> results : register(u0);

cbuffer Params : register(b0) {
    uint sample_count;
};

[numthreads(256, 1, 1)]
void CSMain(uint3 id : SV_DispatchThreadID)
{
    if (id.x >= sample_count) return;

    float t = (float)id.x / (float)(sample_count - 1);

    qaws_vec2 p0 = { cp[0], cp[1] };
    qaws_vec2 p1 = { cp[2], cp[3] };
    qaws_vec2 p2 = { cp[4], cp[5] };
    qaws_vec2 p3 = { cp[6], cp[7] };

    qaws_vec2 pos = qaws_bezier3_eval_2d(p0, p1, p2, p3, t);

    results[id.x * 2]     = pos.x;
    results[id.x * 2 + 1] = pos.y;
}
```

---

## 支持的曲线类型

| 曲线类型              | 核心头文件                          | 主要函数                         |
|-----------------------|-------------------------------------|----------------------------------|
| Bézier（de Casteljau）| `core/qaws_decasteljau_core.h`      | `qaws_bezier3_eval_2d/3d`        |
| Hermite / Catmull-Rom | `core/qaws_cubic_poly_core.h`       | `qaws_cubic_eval_2d/3d`          |
| B-Spline / NURBS      | `core/qaws_bspline_eval_core.h` + `core/qaws_bspline_basis_core.h` | `qaws_bspline_eval_2d/3d` |
| 圆弧（Arc）           | `core/qaws_arc_core.h`              | `qaws_arc_eval_3d`               |

---

## 示例文件速查

仓库 `examples/` 目录下提供了完整可运行的着色器示例：

| 文件                      | 后端  | 说明                          |
|---------------------------|-------|-------------------------------|
| `eval_bezier.hlsl`        | HLSL  | Cubic Bézier 计算着色器        |
| `eval_bezier.comp`        | GLSL  | Cubic Bézier 计算着色器（Vulkan） |
| `eval_cubic_poly.hlsl`    | HLSL  | Hermite / CatRom 三次多项式    |
| `eval_cubic_poly.comp`    | GLSL  | Hermite / CatRom 三次多项式    |
| `eval_bspline.hlsl`       | HLSL  | B-Spline（含节点区间搜索）     |
| `eval_bspline.comp`       | GLSL  | B-Spline（含节点区间搜索）     |
| `eval_arc.hlsl`           | HLSL  | 三维圆弧求值                   |
| `eval_bezier_halide.cpp`  | Halide| Cubic Bézier AOT 生成器        |

---

## 平台抽象宏（跨后端兼容）

| 宏                     | C (float)         | HLSL              | GLSL              |
|------------------------|-------------------|-------------------|-------------------|
| `qaws_scalar`          | `float`           | `float`           | `float`           |
| `QAWS_INLINE`          | `static inline`   | `inline`          | （空）            |
| `QAWS_OUT`             | （空）            | `out`             | `out`             |
| `QAWS_INOUT`           | （空）            | `inout`           | `inout`           |
| `qaws_lerp(a,b,t)`     | C 函数            | `lerp`            | `mix`             |
| `QAWS_ATAN2(y,x)`      | `atan2f`          | `atan2`           | `atan(y,x)`       |
| `QAWS_FMOD(a,b)`       | `fmodf`           | `fmod`            | `mod`             |

---

## 关键约束与注意事项

1. **无指针、无动态内存**：GPU 后端（HLSL/GLSL）只包含纯数学内联函数，完全 shader-safe，不依赖堆分配。
2. **每线程一个采样点**：Dispatch `ceil(sample_count / 256)` 个线程组，每线程计算参数 `t` 对应的曲线坐标。
3. **B-Spline 需要额外节点向量**：比 Bézier 多一步 `qaws_find_span()` 节点区间搜索，需将节点向量也上传到 GPU。
4. **GPU 精度固定为 float**：HLSL/GLSL 后端标量类型固定为 `float`；若需 `double` 精度，只能在 C 后端使用（通过 `QAWS_SCALAR_IS_FLOAT=OFF` 构建）。
5. **控制点布局**：平坦数组按 `[x0, y0, x1, y1, ...]`（2D）或 `[x0, y0, z0, x1, y1, z1, ...]`（3D）排列。
6. **结果布局**：输出 Buffer 同样按分量交错存储，读回后按步长 2（2D）或 3（3D）拆分。

---

## 参考文档

- [`docs/backends.md`](../backends.md) — 多后端架构详解
- [`docs/evaluation.md`](../evaluation.md) — 求值标志、批处理求值、core 头文件
- [`docs/feature_per_backends.md`](../feature_per_backends.md) — 各后端函数支持矩阵
- [`examples/`](../../examples/) — 完整可运行着色器示例
