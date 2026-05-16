#ifndef QAWS_ACDB_CURVE_H
#define QAWS_ACDB_CURVE_H

/*
 * acdb_curve.h — AcDbCurve: abstract base class for all curve entities
 *
 * Revised from CADplatformer-bak sources/core/acdb/acdb_curve.h:
 *   - Removed AcGiWorldDraw dependency (replaced with acgi forward reference).
 *   - Removed GLM dependency.
 *   - Uses AcGePoint3d from acge_types.h.
 *
 * ObjectARX lineage: AcDbObject -> AcDbEntity -> AcDbCurve
 */

#ifdef __cplusplus

#include <cstring>
#include <vector>
#include "acdb_types.h"
#include "acgi/acgi_geometry.h"

namespace acdb {

/* -------------------------------------------------------------------------
 * AcDbCurve — abstract base for all parameterised curve entities
 * ------------------------------------------------------------------------- */
class AcDbCurve
{
public:
    AcDbCurve() = default;
    virtual ~AcDbCurve() = default;

    /* Identity */
    virtual const char* className() const { return "AcDbCurve"; }
    virtual bool isKindOf(const char* name) const
    {
        return std::strcmp(name, "AcDbCurve") == 0;
    }

    /* ---- Curve protocol ---- */

    virtual bool isClosed()   const { return false; }
    virtual bool isPeriodic() const { return false; }

    virtual ErrorStatus getStartParam(double& param) const
    { (void)param; return ErrorStatus::eNotImplementedYet; }
    virtual ErrorStatus getEndParam(double& param) const
    { (void)param; return ErrorStatus::eNotImplementedYet; }

    virtual ErrorStatus getStartPoint(AcGePoint3d& pt) const
    { (void)pt; return ErrorStatus::eNotImplementedYet; }
    virtual ErrorStatus getEndPoint(AcGePoint3d& pt) const
    { (void)pt; return ErrorStatus::eNotImplementedYet; }

    virtual ErrorStatus getPointAtParam(double param, AcGePoint3d& pt) const
    { (void)param; (void)pt; return ErrorStatus::eNotImplementedYet; }

    virtual ErrorStatus getParamAtPoint(const AcGePoint3d& pt, double& param) const
    { (void)pt; (void)param; return ErrorStatus::eNotImplementedYet; }

    virtual ErrorStatus getDistAtParam(double param, double& dist) const
    { (void)param; (void)dist; return ErrorStatus::eNotImplementedYet; }

    virtual ErrorStatus getArea(double& area) const
    { (void)area; return ErrorStatus::eNotImplementedYet; }

    virtual ErrorStatus getExtents(AcGePoint3d& /*minPt*/, AcGePoint3d& /*maxPt*/) const
    { return ErrorStatus::eNotImplementedYet; }

    virtual ErrorStatus copyFrom(const AcDbCurve* /*src*/)
    { return ErrorStatus::eNotImplementedYet; }

    virtual AcDbCurve* clone() const { return nullptr; }

    /**
     * Sample @p numSamples points on the curve.
     * The default implementation calls getPointAtParam over a uniform
     * parameter grid between getStartParam and getEndParam.
     */
    virtual void getSamplePoints(int numSamples,
                                 std::vector<AcGePoint3d>& pts) const
    {
        pts.clear();
        if (numSamples < 2) return;
        double pStart = 0.0, pEnd = 1.0;
        getStartParam(pStart);
        getEndParam(pEnd);
        pts.reserve(static_cast<size_t>(numSamples));
        for (int i = 0; i < numSamples; ++i) {
            double t = pStart + (pEnd - pStart) *
                       static_cast<double>(i) / static_cast<double>(numSamples - 1);
            AcGePoint3d pt;
            getPointAtParam(t, pt);
            pts.push_back(pt);
        }
    }

    /**
     * Render the curve via the AcGiGeometry drawing interface.
     * The default implementation samples the curve and draws a polyline.
     */
    virtual bool worldDraw(acgi::AcGiGeometry& geom) const
    {
        std::vector<AcGePoint3d> pts;
        getSamplePoints(0, pts);
        if (pts.size() < 2) return true;

        std::vector<float> flat;
        flat.reserve(pts.size() * 3);
        for (const auto& p : pts) {
            flat.push_back(static_cast<float>(p.x));
            flat.push_back(static_cast<float>(p.y));
            flat.push_back(static_cast<float>(p.z));
        }
        geom.polyline(static_cast<int>(pts.size()), flat.data());
        return true;
    }

    /* ---- Entity traits (color, linetype, etc.) ---- */
    AcCmColor color()     const { return m_color;     }
    LineWeight lineWeight() const { return m_lineWeight; }
    void setColor(const AcCmColor& c) { m_color = c; }
    void setLineWeight(LineWeight w)  { m_lineWeight = w; }

protected:
    AcCmColor  m_color      = {255, 255, 255, 255};
    LineWeight m_lineWeight = LineWeight::kLnWtByLayer;
};

} /* namespace acdb */

#endif /* __cplusplus */
#endif /* QAWS_ACDB_CURVE_H */
