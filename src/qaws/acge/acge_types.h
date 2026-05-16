#ifndef QAWS_ACGE_TYPES_H
#define QAWS_ACGE_TYPES_H

/*
 * acge_types.h — AcGe-compatible type definitions
 *
 * Naming conventions follow ObjectARX AcGe.
 * Uses only C++ standard library; zero external dependencies.
 * Scalar type is always `double` (AcGe geometry precision).
 *
 * qaws types are included via extern "C" so that the core arc/b-spline
 * evaluation headers can be used from C++ translation units.
 */

#ifdef __cplusplus

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#ifdef __cplusplus
extern "C" {
#endif
#include "qaws_types.h"
#include "qaws_platform.h"
#ifdef __cplusplus
}
#endif

namespace acge {

/* -------------------------------------------------------------------------
 * Fundamental scalar type
 * ------------------------------------------------------------------------- */
using real = double;

/* -------------------------------------------------------------------------
 * Mathematical constants
 * ------------------------------------------------------------------------- */
static constexpr real kPi       = 3.14159265358979323846264338327950288;
static constexpr real kTwoPi    = real(2.0) * kPi;
static constexpr real EPSILON_ACGE = real(1.49011611938e-08);
static constexpr real INFTY     = std::numeric_limits<real>::infinity();

/* -------------------------------------------------------------------------
 * 2-D / 3-D point and vector types
 *
 * Simple aggregates — no GLM dependency.
 * Operators are provided for convenience.
 * ------------------------------------------------------------------------- */
struct AcGePoint2d {
    real x = 0.0, y = 0.0;
    AcGePoint2d() = default;
    AcGePoint2d(real x_, real y_) : x(x_), y(y_) {}
    AcGePoint2d  operator+(const AcGePoint2d& v) const { return {x + v.x, y + v.y}; }
    AcGePoint2d  operator-(const AcGePoint2d& v) const { return {x - v.x, y - v.y}; }
    AcGePoint2d  operator*(real s)               const { return {x * s, y * s}; }
    AcGePoint2d& operator+=(const AcGePoint2d& v) { x += v.x; y += v.y; return *this; }
    AcGePoint2d& operator-=(const AcGePoint2d& v) { x -= v.x; y -= v.y; return *this; }
    real dot(const AcGePoint2d& v) const { return x*v.x + y*v.y; }
    real length() const { return std::sqrt(x*x + y*y); }
    real lengthSqrd() const { return x*x + y*y; }
};
using AcGeVector2d = AcGePoint2d;

struct AcGePoint3d {
    real x = 0.0, y = 0.0, z = 0.0;
    AcGePoint3d() = default;
    AcGePoint3d(real x_, real y_, real z_) : x(x_), y(y_), z(z_) {}
    AcGePoint3d  operator+(const AcGePoint3d& v) const { return {x+v.x, y+v.y, z+v.z}; }
    AcGePoint3d  operator-(const AcGePoint3d& v) const { return {x-v.x, y-v.y, z-v.z}; }
    AcGePoint3d  operator*(real s)               const { return {x*s, y*s, z*s}; }
    AcGePoint3d& operator+=(const AcGePoint3d& v) { x+=v.x; y+=v.y; z+=v.z; return *this; }
    AcGePoint3d& operator-=(const AcGePoint3d& v) { x-=v.x; y-=v.y; z-=v.z; return *this; }
    real dot(const AcGePoint3d& v)   const { return x*v.x + y*v.y + z*v.z; }
    AcGePoint3d cross(const AcGePoint3d& v) const {
        return {y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x};
    }
    real length()    const { return std::sqrt(x*x + y*y + z*z); }
    real lengthSqrd() const { return x*x + y*y + z*z; }
    AcGePoint3d normal() const {
        real len = length();
        return (len > EPSILON_ACGE) ? AcGePoint3d{x/len, y/len, z/len} : AcGePoint3d{};
    }
};
using AcGeVector3d = AcGePoint3d;

/* -------------------------------------------------------------------------
 * Conversion helpers: acge ↔ qaws vector types
 * ------------------------------------------------------------------------- */
inline qaws_vec2 to_qaws(const AcGePoint2d& p)
{
    return {static_cast<qaws_scalar>(p.x), static_cast<qaws_scalar>(p.y)};
}

inline qaws_vec3 to_qaws(const AcGePoint3d& p)
{
    return {static_cast<qaws_scalar>(p.x),
            static_cast<qaws_scalar>(p.y),
            static_cast<qaws_scalar>(p.z)};
}

inline AcGePoint2d from_qaws(const qaws_vec2& v)
{
    return {static_cast<real>(v.x), static_cast<real>(v.y)};
}

inline AcGePoint3d from_qaws(const qaws_vec3& v)
{
    return {static_cast<real>(v.x),
            static_cast<real>(v.y),
            static_cast<real>(v.z)};
}

} /* namespace acge */

#endif /* __cplusplus */

#endif /* QAWS_ACGE_TYPES_H */
