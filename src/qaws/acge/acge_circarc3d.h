#ifndef QAWS_ACGE_CIRCARC3D_H
#define QAWS_ACGE_CIRCARC3D_H

/*
 * acge_circarc3d.h — AcGeCircArc3d: 3D circular arc
 *
 * Curve evaluation is delegated to qaws_arc_eval_3d / qaws_arc_eval_full_3d
 * from the qaws core arc header.
 *
 * The arc lies in a plane defined by two orthogonal unit vectors
 * (m_axis_u, m_axis_v); m_axis_u is the "cos" direction and m_axis_v is the
 * "sin" direction (i.e. 90° CCW when viewed from m_normal).
 *
 * ObjectARX lineage: AcGeEntity3d -> AcGeCurve3d -> AcGeCircArc3d
 * Reference: ObjectARX 2026 AcGeCircArc3d (gecarc3.h)
 */

#ifdef __cplusplus

#include "acge_types.h"

extern "C" {
#include "core/qaws_arc_core.h"
}

namespace acge {

class AcGeCircArc3d
{
public:
    /* ------------------------------------------------------------------ */
    /* Data members                                                        */
    /* ------------------------------------------------------------------ */
    AcGePoint3d  m_center;
    real         m_radius   = 1.0;
    real         m_startAng = 0.0;     /* radians */
    real         m_endAng   = kTwoPi;  /* radians */
    AcGeVector3d m_refVec   = {1.0, 0.0, 0.0}; /* angle=0 direction (axis_u unit) */
    AcGeVector3d m_normal   = {0.0, 0.0, 1.0}; /* arc-plane normal */

    /* ------------------------------------------------------------------ */
    /* Constructors                                                        */
    /* ------------------------------------------------------------------ */
    AcGeCircArc3d() = default;

    AcGeCircArc3d(const AcGePoint3d& center, const AcGeVector3d& normal,
                  real radius)
        : m_center(center), m_radius(radius), m_normal(normal.normal()) {}

    AcGeCircArc3d(const AcGePoint3d& center, const AcGeVector3d& normal,
                  const AcGeVector3d& refVec,
                  real radius, real startAng, real endAng)
        : m_center(center), m_radius(radius)
        , m_startAng(startAng), m_endAng(endAng)
        , m_refVec(refVec.normal()), m_normal(normal.normal()) {}

    /* ------------------------------------------------------------------ */
    /* Internal helpers: build qaws axis vectors                          */
    /* ------------------------------------------------------------------ */
private:
    /** axis_u = radius * m_refVec */
    qaws_vec3 axisU() const
    {
        return {static_cast<qaws_scalar>(m_radius * m_refVec.x),
                static_cast<qaws_scalar>(m_radius * m_refVec.y),
                static_cast<qaws_scalar>(m_radius * m_refVec.z)};
    }

    /** axis_v = radius * (m_normal × m_refVec), the 90° CCW direction. */
    qaws_vec3 axisV() const
    {
        AcGeVector3d v = m_normal.cross(m_refVec);  /* CCW when viewed from normal */
        return {static_cast<qaws_scalar>(m_radius * v.x),
                static_cast<qaws_scalar>(m_radius * v.y),
                static_cast<qaws_scalar>(m_radius * v.z)};
    }

public:
    /* ------------------------------------------------------------------ */
    /* Core evaluation — delegates to qaws_arc_eval_3d                   */
    /* ------------------------------------------------------------------ */

    /**
     * Evaluate at angle @p param (radians).
     * Converts the angle to a normalized t ∈ [0,1] and calls qaws.
     */
    AcGePoint3d pointAt(real param) const
    {
        qaws_vec3 center_v = to_qaws(m_center);
        real range = m_endAng - m_startAng;
        real t = (std::abs(range) > EPSILON_ACGE)
                     ? (param - m_startAng) / range
                     : real(0);
        qaws_vec3 r = qaws_arc_eval_3d(
            center_v, axisU(), axisV(),
            static_cast<qaws_scalar>(m_startAng),
            static_cast<qaws_scalar>(range),
            static_cast<qaws_scalar>(t));
        return from_qaws(r);
    }

    /**
     * Evaluate at normalized parameter t ∈ [0, 1].
     */
    AcGePoint3d pointAtParam(real t) const
    {
        qaws_vec3 center_v = to_qaws(m_center);
        real range = m_endAng - m_startAng;
        qaws_vec3 r = qaws_arc_eval_3d(
            center_v, axisU(), axisV(),
            static_cast<qaws_scalar>(m_startAng),
            static_cast<qaws_scalar>(range),
            static_cast<qaws_scalar>(t));
        return from_qaws(r);
    }

    /**
     * Evaluate position and first derivative at normalized t ∈ [0, 1].
     * Uses qaws_arc_eval_full_3d to obtain position and tangent in one call.
     */
    void evalFull(real t, AcGePoint3d& pos, AcGeVector3d& d1) const
    {
        qaws_vec3 center_v = to_qaws(m_center);
        real range = m_endAng - m_startAng;
        qaws_eval_3d r = qaws_arc_eval_full_3d(
            center_v, axisU(), axisV(),
            static_cast<qaws_scalar>(m_startAng),
            static_cast<qaws_scalar>(range),
            static_cast<qaws_scalar>(t),
            0x1 | 0x2);   /* POSITION | D1 */
        pos = from_qaws(r.position);
        d1  = from_qaws(r.d1);
    }

    /* ------------------------------------------------------------------ */
    /* Convenience accessors                                               */
    /* ------------------------------------------------------------------ */
    AcGePoint3d startPoint() const { return pointAt(m_startAng); }
    AcGePoint3d endPoint()   const { return pointAt(m_endAng);   }

    real length()   const { return std::abs(m_endAng - m_startAng) * m_radius; }
    bool isClosed() const
    {
        return std::abs(m_endAng - m_startAng - kTwoPi) < EPSILON_ACGE ||
               std::abs(m_endAng - m_startAng + kTwoPi) < EPSILON_ACGE;
    }
    bool isDegenerate(real tol = EPSILON_ACGE) const
    {
        return m_radius <= tol || std::abs(m_endAng - m_startAng) <= tol;
    }

    /* ------------------------------------------------------------------ */
    /* Sampling                                                            */
    /* ------------------------------------------------------------------ */
    /**
     * Sample @p n evenly-spaced points via qaws_arc_eval_3d.
     */
    void getSamplePoints(int n, std::vector<AcGePoint3d>& pts) const
    {
        if (n < 2) n = 64;
        pts.clear();
        pts.reserve(static_cast<size_t>(n));

        qaws_vec3 center_v = to_qaws(m_center);
        qaws_vec3 au = axisU(), av = axisV();
        real range = m_endAng - m_startAng;

        for (int i = 0; i < n; ++i) {
            real t = static_cast<real>(i) / static_cast<real>(n - 1);
            qaws_vec3 p = qaws_arc_eval_3d(
                center_v, au, av,
                static_cast<qaws_scalar>(m_startAng),
                static_cast<qaws_scalar>(range),
                static_cast<qaws_scalar>(t));
            pts.push_back(from_qaws(p));
        }
    }

    /* ------------------------------------------------------------------ */
    /* Angle accessors                                                     */
    /* ------------------------------------------------------------------ */
    void setAngles(real start, real end) { m_startAng = start; m_endAng = end; }
    void setNormal(const AcGeVector3d& n) { m_normal = n.normal(); }
    void setRefVec(const AcGeVector3d& v) { m_refVec = v.normal(); }

    /* ------------------------------------------------------------------ */
    /* Transformation helpers                                             */
    /* ------------------------------------------------------------------ */
    void translateBy(const AcGeVector3d& offset) { m_center += offset; }

    void rotateBy(real angle, const AcGeVector3d& axis, const AcGePoint3d& pivot)
    {
        /* Rodrigues' rotation formula applied to m_center and m_refVec. */
        auto rodrigues = [&](const AcGeVector3d& v, real a,
                              const AcGeVector3d& k) -> AcGeVector3d
        {
            real c = std::cos(a), s = std::sin(a);
            return v * c + k.cross(v) * s + k * (k.dot(v) * (real(1) - c));
        };

        AcGeVector3d k = axis.normal();
        AcGeVector3d toCenter = m_center - pivot;
        m_center = pivot + rodrigues(toCenter, angle, k);
        m_refVec = rodrigues(m_refVec, angle, k);
        m_normal = rodrigues(m_normal, angle, k);
        m_startAng += angle;
        m_endAng   += angle;
    }
};

} /* namespace acge */

#endif /* __cplusplus */
#endif /* QAWS_ACGE_CIRCARC3D_H */
