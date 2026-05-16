#ifndef QAWS_ACGI_GEOMETRY_H
#define QAWS_ACGI_GEOMETRY_H

/*
 * acgi_geometry.h — AcGiGeometry-compatible drawing interface
 *
 * Revised from CADplatformer-bak sources/world/render/acgi_geometry.h:
 *   - Removed bgfx dependency.
 *   - Removed Windows.h and platform-specific headers.
 *   - Provides a pure abstract C++ interface that any rendering backend
 *     can implement (e.g. bgfx, OpenGL, software rasteriser, test stubs).
 *
 * Naming conventions follow ObjectARX AcGiGeometry:
 *   polyline()       — draw a polyline from a flat xyz float array
 *   circularArc()    — draw an arc by center, radius, angles
 *   circle()         — draw a full circle
 *   worldLine()      — draw a single line segment
 *   mesh()           — draw a wireframe mesh (rows × cols)
 *   shell()          — draw a face-indexed shell
 *   ellipticalArc()  — draw an elliptical arc
 *
 * Reference: ObjectARX 2026 AcGiGeometry.html
 */

#ifdef __cplusplus

#include <cstddef>
#include <cstdint>
#include <vector>

namespace acgi {

/* -------------------------------------------------------------------------
 * AcGiSubEntityTraits — render property bag (ObjectARX naming)
 * ------------------------------------------------------------------------- */

/** Line end-cap style. */
enum class AcGiLineEndStyle : uint8_t {
    kEndStyleFlat     = 0,
    kEndStyleSquare   = 1,
    kEndStyleRound    = 2,
};

/** Line join style. */
enum class AcGiLineJoinStyle : uint8_t {
    kJoinStyleNone  = 0,
    kJoinStyleBevel = 1,
    kJoinStyleMiter = 2,
};

/** Fill type. */
enum class AcGiFillType : uint8_t {
    kAcGiFillNever  = 0,
    kAcGiFillAlways = 1,
};

struct AcGiSubEntityTraits {
    /* Standard DWG line types. */
    enum LineType : uint8_t {
        kContinuous = 0, kDashed, kHidden, kCenter,
        kPhantom, kDot, kDashDot, kDashDotDot
    };

    /* Setters */
    void setTrueColor(uint32_t rgba)       { m_color = rgba; }
    void setLineWeight(uint8_t px)         { m_lineWeight = px; }
    void setLineType(LineType lt)          { m_lineType = lt; }
    void setLineTypeScale(float s)         { m_lineTypeScale = s; }
    void setFillType(AcGiFillType ft)      { m_fillType = ft; }
    void setEndStyle(AcGiLineEndStyle s)   { m_cap = s; }
    void setJoinStyle(AcGiLineJoinStyle j) { m_join = j; }
    void setVisibility(bool b)             { m_visible = b; }

    /* Getters */
    uint32_t          trueColor()     const { return m_color; }
    uint8_t           lineWeight()    const { return m_lineWeight; }
    LineType          lineType()      const { return m_lineType; }
    float             lineTypeScale() const { return m_lineTypeScale; }
    AcGiFillType      fillType()      const { return m_fillType; }
    AcGiLineEndStyle  endStyle()      const { return m_cap; }
    AcGiLineJoinStyle joinStyle()     const { return m_join; }
    bool              visible()       const { return m_visible; }

private:
    uint32_t          m_color         = 0xFFFFFFFF; /* RGBA white */
    uint8_t           m_lineWeight    = 1;
    LineType          m_lineType      = kContinuous;
    float             m_lineTypeScale = 1.0f;
    AcGiFillType      m_fillType      = AcGiFillType::kAcGiFillNever;
    AcGiLineEndStyle  m_cap           = AcGiLineEndStyle::kEndStyleRound;
    AcGiLineJoinStyle m_join          = AcGiLineJoinStyle::kJoinStyleMiter;
    bool              m_visible       = true;
};

/* -------------------------------------------------------------------------
 * AcGiGeometry — abstract drawing interface
 *
 * Rendering backends (e.g. bgfx, OpenGL, SVG, test stubs) derive from this
 * class and implement the virtual primitives.
 * ------------------------------------------------------------------------- */
class AcGiGeometry
{
public:
    virtual ~AcGiGeometry() = default;

    /**
     * Draw a polyline from a flat array of 3-D points.
     * @param numPoints  Number of vertices.
     * @param xyz        Flat float array [x0,y0,z0, x1,y1,z1, ...].
     */
    virtual void polyline(int numPoints, const float* xyz) = 0;

    /**
     * Draw a 3-D circular arc.
     * @param center      World-space center [x,y,z].
     * @param radius      Arc radius.
     * @param startAngle  Start angle in radians.
     * @param sweepAngle  Angular sweep in radians (may be negative for CW).
     * @param normal      Arc-plane normal [nx,ny,nz] (unit vector).
     */
    virtual void circularArc(const float center[3],
                             float radius,
                             float startAngle,
                             float sweepAngle,
                             const float normal[3]) = 0;

    /**
     * Draw a full 3-D circle.
     * @param center  World-space center [x,y,z].
     * @param radius  Circle radius.
     * @param normal  Circle-plane normal [nx,ny,nz] (unit vector).
     */
    virtual void circle(const float center[3],
                        float radius,
                        const float normal[3]) = 0;

    /**
     * Draw a single line segment.
     * @param start  Start point [x,y,z].
     * @param end    End point [x,y,z].
     */
    virtual void worldLine(const float start[3], const float end[3]) = 0;

    /**
     * Draw a parametric elliptical arc.
     * @param center     Arc center [x,y,z].
     * @param majorAxis  Major semi-axis vector [x,y,z].
     * @param minorAxis  Minor semi-axis vector [x,y,z].
     * @param startAngle Start eccentric anomaly in radians.
     * @param endAngle   End eccentric anomaly in radians.
     */
    virtual void ellipticalArc(const float center[3],
                               const float majorAxis[3],
                               const float minorAxis[3],
                               float startAngle,
                               float endAngle) = 0;

    /**
     * Draw a wireframe mesh.
     * @param rows  Number of rows.
     * @param cols  Number of columns.
     * @param xyz   Flat float array of (rows * cols) points [x,y,z,...].
     */
    virtual void mesh(int rows, int cols, const float* xyz) = 0;

    /**
     * Draw a face-indexed shell (wireframe).
     * @param numPoints  Number of vertices.
     * @param xyz        Flat float array of vertices.
     * @param numFaces   Number of face entries.
     * @param faces      Face list: each face is [count, i0, i1, ..., iN].
     */
    virtual void shell(int numPoints, const float* xyz,
                       int numFaces,  const int* faces) = 0;
};

/* -------------------------------------------------------------------------
 * AcGiWorldDraw — associates traits with a geometry drawer
 *
 * Matches the ObjectARX AcGiWorldDraw interface (simplified).
 * ------------------------------------------------------------------------- */
class AcGiWorldDraw
{
public:
    virtual ~AcGiWorldDraw() = default;

    virtual AcGiGeometry&       geometry()       = 0;
    virtual AcGiSubEntityTraits& subEntityTraits() = 0;
};

} /* namespace acgi */

#endif /* __cplusplus */
#endif /* QAWS_ACGI_GEOMETRY_H */
