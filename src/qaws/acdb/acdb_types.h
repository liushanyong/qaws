#ifndef QAWS_ACDB_TYPES_H
#define QAWS_ACDB_TYPES_H

/*
 * acdb_types.h — AcDb-compatible type definitions for qaws
 *
 * Types are aligned with ObjectARX AcDb naming conventions.
 * Geometric primitives reuse the AcGe types from acge_types.h so that
 * there is a single set of structs across both modules.
 *
 * Revised from CADplatformer-bak sources/core/acdb/acdb_types.h:
 *   - Removed GLM dependency (replaced with acge simple structs).
 *   - Removed bgfx / Windows headers.
 *   - Uses qaws_scalar where interaction with the qaws C API is needed.
 */

#ifdef __cplusplus

#include <cstdint>
#include <string>
#include <vector>
#include <memory>
#include <functional>

#include "acge/acge_types.h"   /* AcGePoint3d, AcGeVector3d, etc. */

namespace acdb {

/* -------------------------------------------------------------------------
 * Reuse the AcGe geometric types
 * ------------------------------------------------------------------------- */
using AcGePoint2d  = acge::AcGePoint2d;
using AcGeVector2d = acge::AcGeVector2d;
using AcGePoint3d  = acge::AcGePoint3d;
using AcGeVector3d = acge::AcGeVector3d;

/* -------------------------------------------------------------------------
 * Constants
 * ------------------------------------------------------------------------- */
static constexpr double kPi    = 3.14159265358979323846;
static constexpr double kTwoPi = 2.0 * kPi;

/* -------------------------------------------------------------------------
 * Return status codes (Acad::ErrorStatus naming)
 * ------------------------------------------------------------------------- */
enum class ErrorStatus {
    eOk = 0,
    eNotImplementedYet,
    eInvalidInput,
    eNullObjectId,
    eNullObjectPointer,
    eInvalidObjectId,
    eObjectAlreadyInDb,
    eKeyNotFound,
    eDuplicateKey,
    eOutOfRange,
    eFileNotFound,
    eFileAccessError,
    eEndOfFile,
    eInvalidDwgVersion,
    eNotOpenForRead,
    eNotOpenForWrite,
    eWrongObjectType,
    eGeneralError,
    eQawsError,   /* Error returned by an underlying qaws function. */
};

/* -------------------------------------------------------------------------
 * Object open mode
 * ------------------------------------------------------------------------- */
enum class OpenMode {
    kForRead,
    kForWrite,
    kForNotify,
};

/* -------------------------------------------------------------------------
 * Color representation
 * ------------------------------------------------------------------------- */
struct AcCmColor {
    uint8_t r = 255, g = 255, b = 255, a = 255;
    AcCmColor() = default;
    AcCmColor(uint8_t r_, uint8_t g_, uint8_t b_, uint8_t a_ = 255)
        : r(r_), g(g_), b(b_), a(a_) {}
    bool operator==(const AcCmColor& o) const
    { return r==o.r && g==o.g && b==o.b && a==o.a; }
    bool operator!=(const AcCmColor& o) const { return !(*this == o); }
};

/* -------------------------------------------------------------------------
 * Line weight (matching ObjectARX AcDb::LineWeight)
 * ------------------------------------------------------------------------- */
enum class LineWeight : int16_t {
    kLnWt000  =   0,
    kLnWt005  =   5,
    kLnWt009  =   9,
    kLnWt013  =  13,
    kLnWt015  =  15,
    kLnWt018  =  18,
    kLnWt020  =  20,
    kLnWt025  =  25,
    kLnWt030  =  30,
    kLnWt035  =  35,
    kLnWt040  =  40,
    kLnWt050  =  50,
    kLnWt053  =  53,
    kLnWt060  =  60,
    kLnWt070  =  70,
    kLnWt080  =  80,
    kLnWt090  =  90,
    kLnWt100  = 100,
    kLnWt106  = 106,
    kLnWt120  = 120,
    kLnWt140  = 140,
    kLnWt158  = 158,
    kLnWt200  = 200,
    kLnWt211  = 211,
    kLnWtByLayer   = -1,
    kLnWtByBlock   = -2,
    kLnWtByLwDefault = -3,
};

/* -------------------------------------------------------------------------
 * Simple 4×4 affine transform (row-major)
 * ------------------------------------------------------------------------- */
struct AcGeMatrix3d {
    double m[4][4] = {
        {1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,1}
    };
};

} /* namespace acdb */

#endif /* __cplusplus */
#endif /* QAWS_ACDB_TYPES_H */
