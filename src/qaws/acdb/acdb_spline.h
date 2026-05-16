#ifndef QAWS_ACDB_SPLINE_H
#define QAWS_ACDB_SPLINE_H

/*
 * acdb_spline.h — AcDbSpline: NURBS / B-spline curve entity
 *
 * Revised from CADplatformer-bak sources/core/acdb/acdb_spline.h:
 *   - evaluateBSpline delegates to the qaws NURBS/B-spline runtime instead
 *     of the embedded De Boor implementation.
 *   - A qaws_curve is built lazily and cached; any setter that changes the
 *     spline geometry invalidates the cache.
 *   - worldDraw uses acgi::AcGiGeometry (no bgfx).
 *
 * ObjectARX lineage: AcDbObject -> AcDbEntity -> AcDbCurve -> AcDbSpline
 */

#ifdef __cplusplus

#include "acdb_curve.h"
#include <vector>

extern "C" {
#include "qaws_types.h"   /* qaws_curve* forward type */
}

namespace acdb {

class AcDbSpline : public AcDbCurve
{
public:
    AcDbSpline() = default;

    /** Construct from fit points (degree-3 interpolating spline). */
    explicit AcDbSpline(const std::vector<AcGePoint3d>& fitPoints,
                        int degree = 3);

    /** Construct from explicit NURBS data. */
    AcDbSpline(int degree, bool rational, bool closed,
               const std::vector<AcGePoint3d>& controlPoints,
               const std::vector<double>& knots,
               const std::vector<double>& weights = {});

    ~AcDbSpline() override;

    /* Disallow accidental copies (they would share the cached qaws_curve). */
    AcDbSpline(const AcDbSpline&) = delete;
    AcDbSpline& operator=(const AcDbSpline&) = delete;

    /* Identity */
    const char* className() const override { return "AcDbSpline"; }
    bool isKindOf(const char* name) const override;

    /* Accessors */
    int  degree()     const { return m_degree;   }
    bool isRational() const { return m_rational; }
    bool isClosed()   const override { return m_closed; }

    const std::vector<AcGePoint3d>& fitPoints()     const { return m_fitPoints; }
    const std::vector<AcGePoint3d>& controlPoints() const { return m_controlPoints; }
    const std::vector<double>&      knots()         const { return m_knots; }
    const std::vector<double>&      weights()       const { return m_weights; }

    /* Mutators — invalidate the qaws curve cache. */
    void setFitPoints(const std::vector<AcGePoint3d>& pts)
    { m_fitPoints = pts; invalidateCache(); }
    void addFitPoint(const AcGePoint3d& pt)
    { m_fitPoints.push_back(pt); invalidateCache(); }
    void setControlPoints(const std::vector<AcGePoint3d>& pts)
    { m_controlPoints = pts; invalidateCache(); }
    void setKnots(const std::vector<double>& k)
    { m_knots = k; invalidateCache(); }
    void setWeights(const std::vector<double>& w)
    { m_weights = w; invalidateCache(); }
    void setDegree(int d)
    { m_degree = d; invalidateCache(); }
    void setClosed(bool c)
    { m_closed = c; invalidateCache(); }
    void setRational(bool r)
    { m_rational = r; invalidateCache(); }

    /* Curve protocol */
    ErrorStatus getStartParam(double& p) const override;
    ErrorStatus getEndParam(double& p)   const override;
    ErrorStatus getStartPoint(AcGePoint3d& pt) const override;
    ErrorStatus getEndPoint(AcGePoint3d& pt)   const override;
    ErrorStatus getPointAtParam(double param, AcGePoint3d& pt) const override;

    /**
     * Sample @p n points using qaws_curve_sample_3d.
     * Prefer control-point spline when available; fall back to fit points.
     */
    void getSamplePoints(int n, std::vector<AcGePoint3d>& pts) const override;

    ErrorStatus getExtents(AcGePoint3d& minPt, AcGePoint3d& maxPt) const override;
    ErrorStatus copyFrom(const AcDbCurve* src) override;
    AcDbCurve*  clone() const override;

    /** Draw the spline via acgi::AcGiGeometry::polyline(). */
    bool worldDraw(acgi::AcGiGeometry& geom) const override;

    /**
     * Evaluate the B-spline / NURBS at parameter @p t.
     *
     * For NURBS (rational=true) uses qaws_curve_create_nurbs;
     * otherwise uses qaws_curve_create_bspline.
     * The cached qaws_curve is reused across subsequent calls.
     *
     * @param t  Parameter in the curve's natural domain
     *           [knots[degree], knots[n_cp]] (or [0,1] for fit-point splines).
     */
    AcGePoint3d evaluateSpline(double t) const;

private:
    int  m_degree   = 3;
    bool m_rational = false;
    bool m_closed   = false;

    std::vector<AcGePoint3d> m_fitPoints;
    std::vector<AcGePoint3d> m_controlPoints;
    std::vector<double>      m_knots;
    std::vector<double>      m_weights;

    /* Lazy-initialised qaws_curve cache. */
    mutable qaws_curve* m_qawsCurve = nullptr;
    mutable bool        m_dirty     = true;

    /** Destroy the cached curve and mark dirty. */
    void invalidateCache();

    /** Build (or rebuild) the qaws_curve from current data. */
    bool buildCache() const;

    /** Return the active point source: control points if available,
     *  otherwise fit points. */
    const std::vector<AcGePoint3d>& activePoints() const
    {
        return m_controlPoints.empty() ? m_fitPoints : m_controlPoints;
    }
};

} /* namespace acdb */

#endif /* __cplusplus */
#endif /* QAWS_ACDB_SPLINE_H */
