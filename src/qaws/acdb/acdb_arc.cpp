/*
 * acdb_arc.cpp — AcDbArc implementation
 *
 * Revised from CADplatformer-bak sources/core/acdb/private/acdb_arc.cpp:
 *   - getSamplePoints delegates to qaws_curve_create_arc +
 *     qaws_curve_sample_3d instead of a manual for-loop with cos/sin.
 *   - worldDraw uses acgi::AcGiGeometry (no bgfx).
 *
 * The qaws arc API expects axis_u / axis_v (radius-scaled orthogonal plane
 * vectors).  For the standard XY-plane arc (normal = {0,0,1}) these are
 * {radius,0,0} and {0,radius,0}.  For a tilted arc we compute them from the
 * normal vector.
 */

#include "acdb_arc.h"

#include <cstring>
#include <cmath>
#include <algorithm>

extern "C" {
#include "qaws_arc.h"
#include "qaws_curve.h"
#include "qaws_sampling.h"
#include "qaws_inspect.h"
#include "core/qaws_arc_core.h"
}

namespace acdb {

/* -------------------------------------------------------------------------
 * Helper: build qaws_arc_desc from the entity fields
 * ------------------------------------------------------------------------- */
static void buildArcDesc(const AcGePoint3d& center,
                         const AcGeVector3d& normal,
                         double radius,
                         double startAngle,
                         double endAngle,
                         qaws_arc_segment& seg,
                         qaws_arc_desc&    desc)
{
    /* Plane vectors: axis_u is the "cos" direction, axis_v the "sin" direction. */
    /* For a given normal n, pick axis_u perpendicular to n.                     */
    AcGeVector3d nUnit = normal.normal();

    /* Gram–Schmidt: find a vector not parallel to nUnit. */
    AcGeVector3d up = (std::abs(nUnit.z) < 0.9) ?
                          AcGeVector3d{0,0,1} : AcGeVector3d{1,0,0};
    AcGeVector3d axU = nUnit.cross(up);
    double lenU = axU.length();
    if (lenU > 1e-10) {
        axU.x /= lenU; axU.y /= lenU; axU.z /= lenU;
    } else {
        axU = {1,0,0};
    }
    /* axis_v = normal × axis_u  (CCW orientation) */
    AcGeVector3d axV = nUnit.cross(axU);

    seg.center[0] = static_cast<qaws_scalar>(center.x);
    seg.center[1] = static_cast<qaws_scalar>(center.y);
    seg.center[2] = static_cast<qaws_scalar>(center.z);
    seg.radius     = static_cast<qaws_scalar>(radius);
    seg.angle_start = static_cast<qaws_scalar>(startAngle);
    seg.angle_end   = static_cast<qaws_scalar>(endAngle);
    /* axis_u and axis_v stored as unit vectors; qaws scales by radius via
     * its own "radius-scaled axis" convention — so we pass unit vectors and
     * put the radius into seg.radius. */
    seg.axis_u[0] = static_cast<qaws_scalar>(axU.x);
    seg.axis_u[1] = static_cast<qaws_scalar>(axU.y);
    seg.axis_u[2] = static_cast<qaws_scalar>(axU.z);
    seg.axis_v[0] = static_cast<qaws_scalar>(axV.x);
    seg.axis_v[1] = static_cast<qaws_scalar>(axV.y);
    seg.axis_v[2] = static_cast<qaws_scalar>(axV.z);

    desc.dimension     = QAWS_DIMENSION_3D;
    desc.segments      = &seg;
    desc.segment_count = 1;
}

/* -------------------------------------------------------------------------
 * AcDbArc implementation
 * ------------------------------------------------------------------------- */
AcDbArc::AcDbArc() = default;

AcDbArc::AcDbArc(const AcGePoint3d& c, double r, double sa, double ea)
    : m_center(c), m_radius(r), m_startAngle(sa), m_endAngle(ea) {}

AcDbArc::AcDbArc(const AcGePoint3d& c, const AcGeVector3d& n,
                 double r, double sa, double ea)
    : m_center(c), m_radius(r), m_startAngle(sa), m_endAngle(ea), m_normal(n) {}

bool AcDbArc::isKindOf(const char* name) const
{
    return std::strcmp(name, "AcDbArc") == 0 ||
           AcDbCurve::isKindOf(name);
}

double AcDbArc::totalAngle() const
{
    double a = m_endAngle - m_startAngle;
    if (a < 0.0) a += kTwoPi;
    return a;
}

double AcDbArc::arcLength() const { return totalAngle() * m_radius; }

ErrorStatus AcDbArc::getStartParam(double& p) const
{ p = m_startAngle; return ErrorStatus::eOk; }

ErrorStatus AcDbArc::getEndParam(double& p) const
{ p = m_endAngle; return ErrorStatus::eOk; }

ErrorStatus AcDbArc::getStartPoint(AcGePoint3d& pt) const
{
    AcGePoint3d dummy;
    return getPointAtParam(m_startAngle, pt);
}

ErrorStatus AcDbArc::getEndPoint(AcGePoint3d& pt) const
{
    return getPointAtParam(m_endAngle, pt);
}

ErrorStatus AcDbArc::getPointAtParam(double param, AcGePoint3d& pt) const
{
    /* Direct evaluation via the qaws arc core (single point). */
    AcGeVector3d nUnit = m_normal.normal();
    AcGeVector3d up    = (std::abs(nUnit.z) < 0.9) ?
                             AcGeVector3d{0,0,1} : AcGeVector3d{1,0,0};
    AcGeVector3d axU   = nUnit.cross(up);
    double lenU = axU.length();
    if (lenU > 1e-10) { axU.x/=lenU; axU.y/=lenU; axU.z/=lenU; }
    else               { axU = {1,0,0}; }
    AcGeVector3d axV = nUnit.cross(axU);

    double range = m_endAngle - m_startAngle;
    double t     = (std::abs(range) > 1e-15) ?
                       (param - m_startAngle) / range : 0.0;

    /* qaws_arc_eval_3d expects radius-scaled axis vectors. */
    qaws_vec3 c  = { (qaws_scalar)m_center.x, (qaws_scalar)m_center.y, (qaws_scalar)m_center.z };
    qaws_vec3 au = { (qaws_scalar)(m_radius*axU.x), (qaws_scalar)(m_radius*axU.y), (qaws_scalar)(m_radius*axU.z) };
    qaws_vec3 av = { (qaws_scalar)(m_radius*axV.x), (qaws_scalar)(m_radius*axV.y), (qaws_scalar)(m_radius*axV.z) };

    {
        qaws_scalar theta = (qaws_scalar)m_startAngle + (qaws_scalar)t * (qaws_scalar)range;
        qaws_scalar ct = (qaws_scalar)std::cos((double)theta);
        qaws_scalar st = (qaws_scalar)std::sin((double)theta);
        pt.x = (double)(c.x + ct*au.x + st*av.x);
        pt.y = (double)(c.y + ct*au.y + st*av.y);
        pt.z = (double)(c.z + ct*au.z + st*av.z);
    }
    return ErrorStatus::eOk;
}

/* getSamplePoints — uses qaws_curve_sample_3d */
void AcDbArc::getSamplePoints(int n, std::vector<AcGePoint3d>& pts) const
{
    pts.clear();
    if (n < 2) n = 64;

    qaws_arc_segment seg{};
    qaws_arc_desc    desc{};
    buildArcDesc(m_center, m_normal, m_radius,
                 m_startAngle, m_endAngle, seg, desc);

    qaws_curve* curve = nullptr;
    qaws_status st = qaws_curve_create_arc(&desc, &curve);
    if (st != QAWS_STATUS_OK || !curve) {
        /* Fallback: direct cos/sin sampling. */
        double ta = totalAngle();
        pts.reserve(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) {
            double a = m_startAngle + ta * static_cast<double>(i) / (n - 1);
            AcGePoint3d pt;
            getPointAtParam(a, pt);
            pts.push_back(pt);
        }
        return;
    }

    std::vector<qaws_vec3> buf(static_cast<size_t>(n + 1));
    unsigned int got = 0;
    qaws_sampling_desc sdesc{};
    sdesc.traversal_mode        = QAWS_TRAVERSAL_MODE_PARAMETER;
    sdesc.sample_count          = static_cast<unsigned int>(n);
    sdesc.use_adaptive_sampling = 0;

    qaws_curve_sample_3d(curve, &sdesc, buf.data(),
                         static_cast<unsigned int>(buf.size()), &got);
    qaws_curve_destroy(curve);

    pts.reserve(got);
    for (unsigned int i = 0; i < got; ++i) {
        pts.push_back({(double)buf[i].x, (double)buf[i].y, (double)buf[i].z});
    }
}

ErrorStatus AcDbArc::getExtents(AcGePoint3d& minPt, AcGePoint3d& maxPt) const
{
    std::vector<AcGePoint3d> pts;
    getSamplePoints(64, pts);
    if (pts.empty()) return ErrorStatus::eOutOfRange;
    minPt = maxPt = pts[0];
    for (const auto& p : pts) {
        if (p.x < minPt.x) minPt.x = p.x;
        if (p.y < minPt.y) minPt.y = p.y;
        if (p.z < minPt.z) minPt.z = p.z;
        if (p.x > maxPt.x) maxPt.x = p.x;
        if (p.y > maxPt.y) maxPt.y = p.y;
        if (p.z > maxPt.z) maxPt.z = p.z;
    }
    return ErrorStatus::eOk;
}

ErrorStatus AcDbArc::copyFrom(const AcDbCurve* src)
{
    auto* a = dynamic_cast<const AcDbArc*>(src);
    if (!a) return ErrorStatus::eWrongObjectType;
    m_center = a->m_center; m_radius = a->m_radius;
    m_startAngle = a->m_startAngle; m_endAngle = a->m_endAngle;
    m_normal = a->m_normal;
    m_color = a->m_color; m_lineWeight = a->m_lineWeight;
    return ErrorStatus::eOk;
}

AcDbCurve* AcDbArc::clone() const
{
    auto* c = new AcDbArc();
    c->copyFrom(this);
    return c;
}

bool AcDbArc::worldDraw(acgi::AcGiGeometry& geom) const
{
    float center[3] = {
        static_cast<float>(m_center.x),
        static_cast<float>(m_center.y),
        static_cast<float>(m_center.z)
    };
    float normal[3] = {
        static_cast<float>(m_normal.x),
        static_cast<float>(m_normal.y),
        static_cast<float>(m_normal.z)
    };
    float sweep = static_cast<float>(m_endAngle - m_startAngle);
    if (sweep < 0.f) sweep += static_cast<float>(kTwoPi);
    geom.circularArc(center,
                     static_cast<float>(m_radius),
                     static_cast<float>(m_startAngle),
                     sweep,
                     normal);
    return true;
}

} /* namespace acdb */
