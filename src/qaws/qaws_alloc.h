#ifndef QAWS_ALLOC_H
#define QAWS_ALLOC_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_bezier.h"
#include "qaws_hermite.h"
#include "qaws_catmull_rom.h"
#include "qaws_bspline.h"
#include "qaws_nurbs.h"
#include "qaws_trajectory.h"
#include "qaws_surface_bezier.h"
#include "qaws_surface_bspline.h"
#include "qaws_surface_nurbs.h"

/*
 * Allocator-aware creation functions (_ex variants).
 *
 * These behave identically to their non-_ex counterparts but route
 * all internal allocations through the provided qaws_allocator.
 * Pass NULL for allocator to use the default malloc/free.
 *
 * Curves created with a custom allocator MUST be destroyed with
 * qaws_curve_destroy (which reads the stored allocator). The
 * allocator must remain valid until the curve is destroyed.
 */

/* --- Curve creation _ex variants --- */

qaws_status qaws_curve_create_bezier_ex(
	qaws_bezier_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

qaws_status qaws_curve_create_hermite_ex(
	qaws_hermite_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

qaws_status qaws_curve_create_catmull_rom_ex(
	qaws_catmull_rom_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

qaws_status qaws_curve_create_bspline_ex(
	qaws_bspline_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

qaws_status qaws_curve_create_nurbs_ex(
	qaws_nurbs_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

qaws_status qaws_curve_create_trajectory_ex(
	qaws_trajectory_desc const* desc,
	qaws_allocator const* allocator,
	qaws_curve** out_curve);

/* --- Surface creation _ex variants --- */

qaws_status qaws_surface_create_bezier_ex(
	qaws_surface_bezier_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

qaws_status qaws_surface_create_bspline_ex(
	qaws_surface_bspline_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

qaws_status qaws_surface_create_nurbs_ex(
	qaws_surface_nurbs_desc const* desc,
	qaws_allocator const* allocator,
	qaws_surface** out_surface);

#endif /* QAWS_ALLOC_H */
