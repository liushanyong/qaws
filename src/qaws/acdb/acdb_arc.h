#ifndef QAWS_ACDB_ARC_H
#define QAWS_ACDB_ARC_H

/*
 * acdb_arc.h — AcDbArc: arc entity backed by the qaws arc curve
 *
 * Revised from CADplatformer-bak sources/core/acdb/acdb_arc.h:
 *   - arc sampling / evaluation delegate to qaws_curve_create_arc and
 *     qaws_curve_sample_3d, replacing the manual cos/sin loop.
 *   - worldDraw uses the acgi abstract interface (no bgfx).
 *
 * ObjectARX lineage: AcDbObject -> AcDbEntity -> AcDbCurve -> AcDbArc
 */

#ifdef __cplusplus

#include "acdb_curve.h"
#include <vector>

namespace acdb {

class AcDbArc : public AcDbCurve
{
public:
    AcDbArc();
    AcDbArc(const AcGePoint3d& center,
            double radius,
            double startAngle,
            double endAngle);
    AcDbArc(const AcGePoint3d& center,
            const AcGeVector3d& normal,
            double radius,
            double startAngle,
            double endAngle);
    ~AcDbArc() override = default;

    /* Identity */
    const char* className() const override { return "AcDbArc"; }
    bool isKindOf(const char* name) const override;

    /* Accessors */
    AcGePoint3d  center()     const { return m_center; }
    void         setCenter(const AcGePoint3d& c) { m_center = c; }
    double       radius()     const { return m_radius; }
    void         setRadius(double r) { m_radius = r; }
    double       startAngle() const { return m_startAngle; }
    void         setStartAngle(double a) { m_startAngle = a; }
    double       endAngle()   const { return m_endAngle; }
    void         setEndAngle(double a) { m_endAngle = a; }
    AcGeVector3d normal()     const { return m_normal; }
    void         setNormal(const AcGeVector3d& n) { m_normal = n; }

    double totalAngle() const;   /* |endAngle - startAngle|, wrapping CCW */
    double arcLength()  const;

    /* Curve protocol */
    bool isClosed() const override { return false; }
    ErrorStatus getStartParam(double& p) const override;
    ErrorStatus getEndParam(double& p)   const override;
    ErrorStatus getStartPoint(AcGePoint3d& pt) const override;
    ErrorStatus getEndPoint(AcGePoint3d& pt)   const override;
    ErrorStatus getPointAtParam(double param, AcGePoint3d& pt) const override;

    /**
     * Sample @p n points on the arc using qaws_curve_sample_3d.
     * Falls back to direct arc math if qaws curve creation fails.
     */
    void getSamplePoints(int n, std::vector<AcGePoint3d>& pts) const override;

    ErrorStatus getExtents(AcGePoint3d& minPt, AcGePoint3d& maxPt) const override;
    ErrorStatus copyFrom(const AcDbCurve* src) override;
    AcDbCurve*  clone() const override;

    /** Draw arc via acgi::AcGiGeometry::circularArc(). */
    bool worldDraw(acgi::AcGiGeometry& geom) const override;

private:
    AcGePoint3d  m_center     = {0, 0, 0};
    double       m_radius     = 0.0;
    double       m_startAngle = 0.0;
    double       m_endAngle   = kTwoPi;
    AcGeVector3d m_normal     = {0, 0, 1};
};

} /* namespace acdb */

#endif /* __cplusplus */
#endif /* QAWS_ACDB_ARC_H */
