/*
 * acdb_spline.cpp — AcDbSpline implementation
 *
 * Revised from CADplatformer-bak sources/core/acdb/private/acdb_spline.cpp:
 *   - evaluateSpline (was evaluateBSpline) now builds a qaws_curve via
 *     qaws_curve_create_nurbs (rational) or qaws_curve_create_bspline
 *     (non-rational) and calls qaws_curve_evaluate_3d, replacing the
 *     embedded De Boor algorithm.
 *   - getSamplePoints uses qaws_curve_sample_3d.
 *   - worldDraw uses acgi::AcGiGeometry (no bgfx dependency).
 */

#include "acdb_spline.h"

#include <cstring>
#include <cmath>
#include <algorithm>

extern "C" {
#include "qaws_bspline.h"
#include "qaws_nurbs.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include "qaws_sampling.h"
#include "qaws_inspect.h"
}

namespace acdb {

/* -------------------------------------------------------------------------
 * Constructors / destructor
 * ------------------------------------------------------------------------- */
AcDbSpline::AcDbSpline(const std::vector<AcGePoint3d>& fitPoints, int degree)
    : m_degree(degree), m_fitPoints(fitPoints) {}

AcDbSpline::AcDbSpline(int degree, bool rational, bool closed,
                       const std::vector<AcGePoint3d>& controlPoints,
                       const std::vector<double>& knots,
                       const std::vector<double>& weights)
    : m_degree(degree), m_rational(rational), m_closed(closed)
    , m_controlPoints(controlPoints), m_knots(knots), m_weights(weights) {}

AcDbSpline::~AcDbSpline()
{
    invalidateCache();
}

bool AcDbSpline::isKindOf(const char* name) const
{
    return std::strcmp(name, "AcDbSpline") == 0 ||
           AcDbCurve::isKindOf(name);
}

/* -------------------------------------------------------------------------
 * Cache management
 * ------------------------------------------------------------------------- */
void AcDbSpline::invalidateCache()
{
    if (m_qawsCurve) {
        qaws_curve_destroy(m_qawsCurve);
        m_qawsCurve = nullptr;
    }
    m_dirty = true;
}

/**
 * Build a qaws_curve from the current spline data.
 *
 * For rational splines (NURBS) uses qaws_curve_create_nurbs.
 * For non-rational splines uses qaws_curve_create_bspline.
 * For fit-point-only splines (no control points / knots), constructs a
 * uniform B-spline through the fit points as a best approximation.
 */
bool AcDbSpline::buildCache() const
{
    if (!m_dirty && m_qawsCurve) return true;
    if (m_qawsCurve) {
        qaws_curve_destroy(m_qawsCurve);
        m_qawsCurve = nullptr;
    }
    m_dirty = false;

    const auto& pts = activePoints();
    int n = static_cast<int>(pts.size());
    if (n < 2) return false;

    int p = std::min(m_degree, n - 1);

    /* Flatten control points to qaws_scalar array [x0,y0,z0, x1,y1,z1,...] */
    std::vector<qaws_scalar> cpFlat;
    cpFlat.reserve(static_cast<size_t>(n) * 3);
    for (const auto& pt : pts) {
        cpFlat.push_back(static_cast<qaws_scalar>(pt.x));
        cpFlat.push_back(static_cast<qaws_scalar>(pt.y));
        cpFlat.push_back(static_cast<qaws_scalar>(pt.z));
    }

    /* Build or convert knot vector. */
    std::vector<qaws_scalar> knotsF;
    if (!m_knots.empty()) {
        knotsF.reserve(m_knots.size());
        for (double k : m_knots)
            knotsF.push_back(static_cast<qaws_scalar>(k));
    }
    /* If no knots supplied, qaws will use its default clamped uniform knots. */

    qaws_status st = QAWS_STATUS_OK;

    if (m_rational && !m_weights.empty()) {
        /* NURBS path. */
        std::vector<qaws_scalar> wF;
        wF.reserve(m_weights.size());
        for (double w : m_weights)
            wF.push_back(static_cast<qaws_scalar>(w));

        qaws_nurbs_desc desc{};
        desc.dimension           = QAWS_DIMENSION_3D;
        desc.degree              = static_cast<unsigned int>(p);
        desc.control_points      = cpFlat.data();
        desc.control_point_count = static_cast<unsigned int>(n);
        desc.knots               = knotsF.empty() ? nullptr : knotsF.data();
        desc.knot_count          = static_cast<unsigned int>(knotsF.size());
        desc.weights             = wF.data();
        desc.weight_count        = static_cast<unsigned int>(wF.size());
        desc.is_closed           = m_closed ? 1 : 0;

        st = qaws_curve_create_nurbs(&desc, &m_qawsCurve);
    } else {
        /* Non-rational B-spline path. */
        qaws_bspline_desc desc{};
        desc.dimension           = QAWS_DIMENSION_3D;
        desc.degree              = static_cast<unsigned int>(p);
        desc.control_points      = cpFlat.data();
        desc.control_point_count = static_cast<unsigned int>(n);
        desc.knots               = knotsF.empty() ? nullptr : knotsF.data();
        desc.knot_count          = static_cast<unsigned int>(knotsF.size());
        desc.is_uniform          = knotsF.empty() ? 1 : 0;
        desc.is_closed           = m_closed ? 1 : 0;

        st = qaws_curve_create_bspline(&desc, &m_qawsCurve);
    }

    if (st != QAWS_STATUS_OK || !m_qawsCurve) {
        m_qawsCurve = nullptr;
        return false;
    }
    return true;
}

/* -------------------------------------------------------------------------
 * evaluateSpline — replaces the custom De Boor loop
 * ------------------------------------------------------------------------- */
AcGePoint3d AcDbSpline::evaluateSpline(double t) const
{
    if (!buildCache() || !m_qawsCurve) return {0, 0, 0};

    qaws_eval_result_3d r{};
    qaws_status st = qaws_curve_evaluate_3d(
        m_qawsCurve,
        static_cast<qaws_scalar>(t),
        QAWS_EVAL_FLAG_POSITION,
        &r);

    if (st != QAWS_STATUS_OK) return {0, 0, 0};
    return {static_cast<double>(r.position.x),
            static_cast<double>(r.position.y),
            static_cast<double>(r.position.z)};
}

/* -------------------------------------------------------------------------
 * Curve protocol
 * ------------------------------------------------------------------------- */
ErrorStatus AcDbSpline::getStartParam(double& p) const
{
    if (!buildCache() || !m_qawsCurve) { p = 0.0; return ErrorStatus::eQawsError; }
    qaws_range range = qaws_curve_get_parameter_range(m_qawsCurve);
    p = static_cast<double>(range.min_value);
    return ErrorStatus::eOk;
}

ErrorStatus AcDbSpline::getEndParam(double& p) const
{
    if (!buildCache() || !m_qawsCurve) { p = 1.0; return ErrorStatus::eQawsError; }
    qaws_range range = qaws_curve_get_parameter_range(m_qawsCurve);
    p = static_cast<double>(range.max_value);
    return ErrorStatus::eOk;
}

ErrorStatus AcDbSpline::getStartPoint(AcGePoint3d& pt) const
{
    double pStart = 0.0;
    getStartParam(pStart);
    pt = evaluateSpline(pStart);
    return ErrorStatus::eOk;
}

ErrorStatus AcDbSpline::getEndPoint(AcGePoint3d& pt) const
{
    double pEnd = 1.0;
    getEndParam(pEnd);
    pt = evaluateSpline(pEnd);
    return ErrorStatus::eOk;
}

ErrorStatus AcDbSpline::getPointAtParam(double param, AcGePoint3d& pt) const
{
    pt = evaluateSpline(param);
    return ErrorStatus::eOk;
}

/* getSamplePoints — uses qaws_curve_sample_3d */
void AcDbSpline::getSamplePoints(int n, std::vector<AcGePoint3d>& pts) const
{
    pts.clear();
    const auto& src = activePoints();
    if (src.size() < 2) return;
    if (n < 2) n = std::max(32, static_cast<int>(src.size()) * 4);

    if (!buildCache() || !m_qawsCurve) return;

    std::vector<qaws_vec3> buf(static_cast<size_t>(n + 1));
    unsigned int got = 0;
    qaws_sampling_desc sdesc{};
    sdesc.traversal_mode        = QAWS_TRAVERSAL_MODE_PARAMETER;
    sdesc.sample_count          = static_cast<unsigned int>(n);
    sdesc.use_adaptive_sampling = 0;

    qaws_curve_sample_3d(m_qawsCurve, &sdesc,
                         buf.data(),
                         static_cast<unsigned int>(buf.size()),
                         &got);
    pts.reserve(got);
    for (unsigned int i = 0; i < got; ++i) {
        pts.push_back({static_cast<double>(buf[i].x),
                       static_cast<double>(buf[i].y),
                       static_cast<double>(buf[i].z)});
    }
}

ErrorStatus AcDbSpline::getExtents(AcGePoint3d& minPt, AcGePoint3d& maxPt) const
{
    if (buildCache() && m_qawsCurve) {
        qaws_vec3 qMin{}, qMax{};
        qaws_status st = qaws_curve_compute_bounds_3d(m_qawsCurve, &qMin, &qMax);
        if (st == QAWS_STATUS_OK) {
            minPt = {(double)qMin.x, (double)qMin.y, (double)qMin.z};
            maxPt = {(double)qMax.x, (double)qMax.y, (double)qMax.z};
            return ErrorStatus::eOk;
        }
    }
    /* Fallback: sample-based bounding box. */
    std::vector<AcGePoint3d> sample;
    getSamplePoints(64, sample);
    if (sample.empty()) return ErrorStatus::eOutOfRange;
    minPt = maxPt = sample[0];
    for (const auto& p : sample) {
        if (p.x < minPt.x) minPt.x = p.x;
        if (p.y < minPt.y) minPt.y = p.y;
        if (p.z < minPt.z) minPt.z = p.z;
        if (p.x > maxPt.x) maxPt.x = p.x;
        if (p.y > maxPt.y) maxPt.y = p.y;
        if (p.z > maxPt.z) maxPt.z = p.z;
    }
    return ErrorStatus::eOk;
}

ErrorStatus AcDbSpline::copyFrom(const AcDbCurve* src)
{
    auto* s = dynamic_cast<const AcDbSpline*>(src);
    if (!s) return ErrorStatus::eWrongObjectType;
    m_degree       = s->m_degree;
    m_rational     = s->m_rational;
    m_closed       = s->m_closed;
    m_fitPoints    = s->m_fitPoints;
    m_controlPoints = s->m_controlPoints;
    m_knots        = s->m_knots;
    m_weights      = s->m_weights;
    m_color        = s->m_color;
    m_lineWeight   = s->m_lineWeight;
    invalidateCache();
    return ErrorStatus::eOk;
}

AcDbCurve* AcDbSpline::clone() const
{
    auto* c = new AcDbSpline();
    c->copyFrom(this);
    return c;
}

bool AcDbSpline::worldDraw(acgi::AcGiGeometry& geom) const
{
    std::vector<AcGePoint3d> sample;
    getSamplePoints(0, sample);
    if (sample.size() < 2) return true;

    std::vector<float> flat;
    flat.reserve(sample.size() * 3);
    for (const auto& p : sample) {
        flat.push_back(static_cast<float>(p.x));
        flat.push_back(static_cast<float>(p.y));
        flat.push_back(static_cast<float>(p.z));
    }
    geom.polyline(static_cast<int>(sample.size()), flat.data());
    return true;
}

} /* namespace acdb */
