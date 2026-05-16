#ifndef QAWS_ACGE_CIRCARC2D_H
#define QAWS_ACGE_CIRCARC2D_H

/*
 * acge_circarc2d.h — AcGeCircArc2d: 2D circular arc
 *
 * Curve evaluation is delegated to qaws_arc_eval_2d / qaws_arc_eval_full_2d
 * from the qaws core arc header, replacing the manual cos/sin arithmetic
 * that was used in the original CADplatformer-bak implementation.
 *
 * ObjectARX lineage: AcGeEntity2d -> AcGeCurve2d -> AcGeCircArc2d
 * Reference: ObjectARX 2026 AcGeCircArc2d (gecarc2.h)
 */

#ifdef __cplusplus

#include "acge_types.h"

extern "C" {
#include "core/qaws_arc_core.h"
}

namespace acge {

class AcGeCircArc2d
{
public:
    /* ------------------------------------------------------------------ */
    /* Data members (public, matching ObjectARX layout)                   */
    /* ------------------------------------------------------------------ */
    AcGePoint2d m_center;
    real        m_radius   = 1.0;
    real        m_startAng = 0.0;          /* radians */
    real        m_endAng   = kTwoPi;       /* radians */
    AcGeVector2d m_refVec  = {1.0, 0.0};  /* angle=0 direction */

    /* ------------------------------------------------------------------ */
    /* Constructors                                                        */
    /* ------------------------------------------------------------------ */
    AcGeCircArc2d() = default;

    AcGeCircArc2d(const AcGePoint2d& center, real radius)
        : m_center(center), m_radius(radius) {}

    AcGeCircArc2d(const AcGePoint2d& center, real radius,
                  real startAng, real endAng,
                  const AcGeVector2d& refVec = {1.0, 0.0})
        : m_center(center), m_radius(radius)
        , m_startAng(startAng), m_endAng(endAng)
        , m_refVec(refVec) {}

    /** Construct from three points: start, interior, end. */
    AcGeCircArc2d(const AcGePoint2d& startPt,
                  const AcGePoint2d& interior,
                  const AcGePoint2d& endPt)
    {
        fromThreePoints(startPt, interior, endPt);
    }

    /* ------------------------------------------------------------------ */
    /* Core evaluation — uses qaws_arc_eval_2d                            */
    /* ------------------------------------------------------------------ */

    /**
     * Evaluate the arc at angle @p param (radians).
     *
     * The angle is mapped to a normalized t ∈ [0,1] and forwarded to
     * qaws_arc_eval_2d so that all derivative/rounding logic lives in one
     * place (the qaws core header).
     */
    AcGePoint2d pointAt(real param) const
    {
        qaws_vec2 center_v = to_qaws(m_center);
        /* axis_u and axis_v are radius-scaled orthogonal directions. */
        qaws_vec2 axis_u = {static_cast<qaws_scalar>(m_radius * m_refVec.x),
                            static_cast<qaws_scalar>(m_radius * m_refVec.y)};
        /* axis_v is 90° CCW from axis_u (for a standard CCW circle). */
        qaws_vec2 axis_v = {static_cast<qaws_scalar>(-m_radius * m_refVec.y),
                            static_cast<qaws_scalar>( m_radius * m_refVec.x)};

        real range = m_endAng - m_startAng;
        real t = (std::abs(range) > EPSILON_ACGE)
                     ? (param - m_startAng) / range
                     : real(0);

        qaws_vec2 r = qaws_arc_eval_2d(
            center_v, axis_u, axis_v,
            static_cast<qaws_scalar>(m_startAng),
            static_cast<qaws_scalar>(range),
            static_cast<qaws_scalar>(t));
        return from_qaws(r);
    }

    /**
     * Evaluate arc at normalized parameter t ∈ [0, 1].
     * Wraps qaws_arc_eval_2d directly without angle-to-t conversion.
     */
    AcGePoint2d pointAtParam(real t) const
    {
        qaws_vec2 center_v = to_qaws(m_center);
        qaws_vec2 axis_u = {static_cast<qaws_scalar>(m_radius * m_refVec.x),
                            static_cast<qaws_scalar>(m_radius * m_refVec.y)};
        qaws_vec2 axis_v = {static_cast<qaws_scalar>(-m_radius * m_refVec.y),
                            static_cast<qaws_scalar>( m_radius * m_refVec.x)};
        real range = m_endAng - m_startAng;
        qaws_vec2 r = qaws_arc_eval_2d(
            center_v, axis_u, axis_v,
            static_cast<qaws_scalar>(m_startAng),
            static_cast<qaws_scalar>(range),
            static_cast<qaws_scalar>(t));
        return from_qaws(r);
    }

    /**
     * Evaluate position and first derivative at normalized t ∈ [0, 1].
     * Uses qaws_arc_eval_full_2d — avoids redundant trig.
     */
    void evalFull(real t, AcGePoint2d& pos, AcGeVector2d& d1) const
    {
        qaws_vec2 center_v = to_qaws(m_center);
        qaws_vec2 axis_u = {static_cast<qaws_scalar>(m_radius * m_refVec.x),
                            static_cast<qaws_scalar>(m_radius * m_refVec.y)};
        qaws_vec2 axis_v = {static_cast<qaws_scalar>(-m_radius * m_refVec.y),
                            static_cast<qaws_scalar>( m_radius * m_refVec.x)};
        real range = m_endAng - m_startAng;
        qaws_eval_2d r = qaws_arc_eval_full_2d(
            center_v, axis_u, axis_v,
            static_cast<qaws_scalar>(m_startAng),
            static_cast<qaws_scalar>(range),
            static_cast<qaws_scalar>(t),
            0x1 | 0x2);   /* QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 */
        pos = from_qaws(r.position);
        d1  = from_qaws(r.d1);
    }

    /* ------------------------------------------------------------------ */
    /* Convenience accessors                                               */
    /* ------------------------------------------------------------------ */
    AcGePoint2d startPoint() const { return pointAt(m_startAng); }
    AcGePoint2d endPoint()   const { return pointAt(m_endAng);   }

    /** Arc length = |Δθ| · r */
    real length() const { return std::abs(m_endAng - m_startAng) * m_radius; }

    /** Angle of angle parameter p mapped to the reference frame. */
    real paramOf(const AcGePoint2d& p) const
    {
        AcGePoint2d delta = p - m_center;
        return std::atan2(delta.y, delta.x);
    }

    bool isClosed() const
    {
        return std::abs(m_endAng - m_startAng - kTwoPi) < EPSILON_ACGE ||
               std::abs(m_endAng - m_startAng + kTwoPi) < EPSILON_ACGE;
    }

    bool isClockWise() const
    {
        real delta = m_endAng - m_startAng;
        while (delta >  kPi) delta -= kTwoPi;
        while (delta < -kPi) delta += kTwoPi;
        return delta < 0.0;
    }

    bool isDegenerate(real tol = EPSILON_ACGE) const
    {
        return m_radius <= tol ||
               std::abs(m_endAng - m_startAng) <= tol;
    }

    /* ------------------------------------------------------------------ */
    /* Sampling                                                            */
    /* ------------------------------------------------------------------ */
    /**
     * Sample @p n evenly-spaced points on the arc.
     * Each point is computed via qaws_arc_eval_2d.
     */
    void getSamplePoints(int n, std::vector<AcGePoint2d>& pts) const
    {
        if (n < 2) n = 64;
        pts.clear();
        pts.reserve(static_cast<size_t>(n));

        qaws_vec2 center_v = to_qaws(m_center);
        qaws_vec2 axis_u = {static_cast<qaws_scalar>(m_radius * m_refVec.x),
                            static_cast<qaws_scalar>(m_radius * m_refVec.y)};
        qaws_vec2 axis_v = {static_cast<qaws_scalar>(-m_radius * m_refVec.y),
                            static_cast<qaws_scalar>( m_radius * m_refVec.x)};
        real range = m_endAng - m_startAng;

        for (int i = 0; i < n; ++i) {
            real t = static_cast<real>(i) / static_cast<real>(n - 1);
            qaws_vec2 p = qaws_arc_eval_2d(
                center_v, axis_u, axis_v,
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

    /* ------------------------------------------------------------------ */
    /* Transformation helpers                                             */
    /* ------------------------------------------------------------------ */
    void translateBy(const AcGeVector2d& offset) { m_center += offset; }

    void rotateBy(real angle, const AcGePoint2d& pivot)
    {
        real c = std::cos(angle), s = std::sin(angle);
        AcGePoint2d o = m_center - pivot;
        m_center = pivot + AcGePoint2d{o.x*c - o.y*s, o.x*s + o.y*c};
        AcGeVector2d rv = m_refVec;
        m_refVec = {rv.x*c - rv.y*s, rv.x*s + rv.y*c};
        m_startAng += angle;
        m_endAng   += angle;
    }

    void scaleBy(real factor, const AcGePoint2d& pivot)
    {
        m_center = pivot + (m_center - pivot) * factor;
        m_radius *= std::abs(factor);
    }

private:
    /** Fit a circular arc through three non-collinear points. */
    void fromThreePoints(const AcGePoint2d& s,
                         const AcGePoint2d& m,
                         const AcGePoint2d& e)
    {
        /* Circumcenter formula. */
        real ax = m.x - s.x, ay = m.y - s.y;
        real bx = m.x - e.x, by = m.y - e.y;
        real cross = ax * by - ay * bx;
        if (std::abs(cross) < EPSILON_ACGE) {
            /* Collinear — degenerate arc. */
            m_center = s; m_radius = 0.0; return;
        }
        real aSq = ax*ax + ay*ay;
        real bSq = bx*bx + by*by;
        real cx  = (ay * bSq - by * aSq) / (real(2) * cross);
        real cy  = (bx * aSq - ax * bSq) / (real(2) * cross);
        m_center = AcGePoint2d{m.x + cy, m.y - cx};
        AcGePoint2d toStart = s - m_center;
        m_radius   = toStart.length();
        m_startAng = std::atan2(toStart.y, toStart.x);
        AcGePoint2d toEnd = e - m_center;
        m_endAng   = std::atan2(toEnd.y, toEnd.x);
        m_refVec   = {1.0, 0.0};
    }
};

} /* namespace acge */

#endif /* __cplusplus */
#endif /* QAWS_ACGE_CIRCARC2D_H */
