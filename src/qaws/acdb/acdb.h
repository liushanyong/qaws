#ifndef QAWS_ACDB_H
#define QAWS_ACDB_H

/*
 * acdb.h — AcDb CAD entity library (curve subset) for qaws
 *
 * Revised from CADplatformer-bak sources/core/acdb/ to use qaws as the
 * underlying curve evaluation backend.
 *
 * Include this umbrella header from C++ translation units only.
 * The qaws library (libqaws) must be linked separately.
 * The acdb .cpp files must be compiled and linked into the application.
 *
 * Available entities (curve subset):
 *   AcDbCurve   — abstract base
 *   AcDbArc     — circular arc  (backed by qaws_curve_create_arc)
 *   AcDbSpline  — NURBS / B-spline (backed by qaws_curve_create_nurbs/bspline)
 *
 * Graphics interface:
 *   acgi::AcGiGeometry  — abstract drawing interface (see acgi/acgi.h)
 */

#ifdef __cplusplus

#include "acdb_types.h"
#include "acdb_curve.h"
#include "acdb_arc.h"
#include "acdb_spline.h"

#endif /* __cplusplus */
#endif /* QAWS_ACDB_H */
