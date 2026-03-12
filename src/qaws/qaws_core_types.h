#ifndef QAWS_CORE_TYPES_H
#define QAWS_CORE_TYPES_H

#include "qaws_platform.h"

/*
 * On the C backend, qaws_types.h already defines qaws_vec2/vec3 via typedef.
 * This header is for shader/Halide consumers who include the core headers
 * directly without the full C runtime. The definitions are binary-compatible.
 *
 * When both headers are present (C backend), skip the re-definitions.
 */
/* vec2/vec3: skip when C backend already has them from qaws_types.h */
#if QAWS_BACKEND != QAWS_BACKEND_C || !defined(QAWS_TYPES_H)

QAWS_TYPE_DEF struct {
    qaws_scalar x, y;
} qaws_vec2;

QAWS_TYPE_DEF struct {
    qaws_scalar x, y, z;
} qaws_vec3;

#endif /* QAWS_BACKEND != QAWS_BACKEND_C || !defined(QAWS_TYPES_H) */

/* eval_2d/eval_3d: lightweight result for core headers (always defined) */
#ifndef QAWS_EVAL_2D_DEFINED
#define QAWS_EVAL_2D_DEFINED

QAWS_TYPE_DEF struct {
    qaws_vec2 position;
    qaws_vec2 d1;
    qaws_vec2 d2;
} qaws_eval_2d;

QAWS_TYPE_DEF struct {
    qaws_vec3 position;
    qaws_vec3 d1;
    qaws_vec3 d2;
} qaws_eval_3d;

#endif /* QAWS_EVAL_2D_DEFINED */

#ifndef QAWS_CORE_MAX_DEGREE
# define QAWS_CORE_MAX_DEGREE 16
#endif

#define QAWS_CORE_MAX_POINTS (QAWS_CORE_MAX_DEGREE + 1)

#endif /* QAWS_CORE_TYPES_H */
