#include "qaws_export.h"
#include "qaws_bspline.h"
#include "internal/qaws_internal_basis.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* -------------------------------------------------------------------------- */
/*  Helpers                                                                    */
/* -------------------------------------------------------------------------- */

static qaws_scalar scalar_sqrt(qaws_scalar x)
{
#if QAWS_SCALAR_IS_FLOAT
	return sqrtf(x);
#else
	return sqrt(x);
#endif
}

static qaws_scalar scalar_floor(qaws_scalar x)
{
#if QAWS_SCALAR_IS_FLOAT
	return floorf(x);
#else
	return floor(x);
#endif
}

/* Solve A x = b via Gaussian elimination with partial pivoting.
   A is size x size (row-major, modified in place).
   b is size x dim_count (column-major per dimension: b[d * size + i]).
   Solution is written back into b.
   Returns QAWS_STATUS_OK on success, QAWS_STATUS_NUMERICAL_FAILURE if singular. */
static qaws_status solve_linear_system(
	qaws_scalar* A,
	qaws_scalar* b,
	unsigned int size,
	unsigned int dim_count)
{
	unsigned int i, j, k, d;

	for (i = 0; i < size; i++) {
		/* Partial pivoting: find row with largest absolute value in column i */
		unsigned int max_row = i;
		qaws_scalar max_val = A[i * size + i];
		qaws_scalar tmp;

		if (max_val < (qaws_scalar)0.0)
			max_val = -max_val;

		for (k = i + 1; k < size; k++) {
			qaws_scalar v = A[k * size + i];
			if (v < (qaws_scalar)0.0)
				v = -v;
			if (v > max_val) {
				max_val = v;
				max_row = k;
			}
		}

		if (max_val < (qaws_scalar)1e-12)
			return QAWS_STATUS_NUMERICAL_FAILURE;

		/* Swap rows i and max_row in A and b */
		if (max_row != i) {
			for (j = 0; j < size; j++) {
				tmp = A[i * size + j];
				A[i * size + j] = A[max_row * size + j];
				A[max_row * size + j] = tmp;
			}
			for (d = 0; d < dim_count; d++) {
				tmp = b[d * size + i];
				b[d * size + i] = b[d * size + max_row];
				b[d * size + max_row] = tmp;
			}
		}

		/* Eliminate below */
		for (k = i + 1; k < size; k++) {
			qaws_scalar factor = A[k * size + i] / A[i * size + i];
			A[k * size + i] = (qaws_scalar)0.0;
			for (j = i + 1; j < size; j++)
				A[k * size + j] -= factor * A[i * size + j];
			for (d = 0; d < dim_count; d++)
				b[d * size + k] -= factor * b[d * size + i];
		}
	}

	/* Back substitution */
	for (d = 0; d < dim_count; d++) {
		for (i = size; i-- > 0; ) {
			qaws_scalar sum = b[d * size + i];
			for (j = i + 1; j < size; j++)
				sum -= A[i * size + j] * b[d * size + j];
			b[d * size + i] = sum / A[i * size + i];
		}
	}

	return QAWS_STATUS_OK;
}

/* -------------------------------------------------------------------------- */
/*  B-spline least-squares fitting                                             */
/* -------------------------------------------------------------------------- */

qaws_status qaws_curve_fit_bspline(
	qaws_bspline_fit_desc const* desc,
	qaws_curve** out_curve)
{
	unsigned int m;            /* data_point_count */
	unsigned int n;            /* control_point_count */
	unsigned int p;            /* degree */
	unsigned int dim_count;
	unsigned int knot_count;
	unsigned int sys_size;     /* n - 2: number of unknowns */
	qaws_scalar const* data;

	qaws_scalar* params  = NULL;  /* m parameter values */
	qaws_scalar* knots   = NULL;  /* n + p + 1 knots */
	qaws_scalar* NtN     = NULL;  /* sys_size x sys_size normal matrix */
	qaws_scalar* rhs     = NULL;  /* dim_count x sys_size right-hand side */
	qaws_scalar* basis   = NULL;  /* p + 1 basis function values */
	qaws_scalar* ctrl    = NULL;  /* n * dim_count control point scalars */
	qaws_status status = QAWS_STATUS_OK;

	unsigned int i, j, k, d;

	/* --- Validation ---------------------------------------------------- */

	if (!desc || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_curve = NULL;

	if (desc->dimension != QAWS_DIMENSION_2D &&
	    desc->dimension != QAWS_DIMENSION_3D)
		return QAWS_STATUS_INVALID_DIMENSION;

	if (!desc->data_points)
		return QAWS_STATUS_INVALID_ARGUMENT;

	m = desc->data_point_count;
	n = desc->control_point_count;
	p = desc->degree;
	dim_count = (desc->dimension == QAWS_DIMENSION_2D) ? 2u : 3u;
	data = (qaws_scalar const*)desc->data_points;

	if (p < 1)
		return QAWS_STATUS_INVALID_DEGREE;

	if (m < p + 1)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (n <= p)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;

	if (n > m)
		return QAWS_STATUS_INVALID_CONTROL_POINT_COUNT;

	sys_size = n - 2;

	/* --- Step 1: Parameter assignment ---------------------------------- */

	params = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)m);
	if (!params) { status = QAWS_STATUS_ALLOCATION_FAILURE; goto cleanup; }

	if (desc->parameters) {
		memcpy(params, desc->parameters, sizeof(qaws_scalar) * (size_t)m);
	} else {
		/* Chord-length parameterization */
		qaws_scalar total_length = (qaws_scalar)0.0;
		params[0] = (qaws_scalar)0.0;
		for (i = 1; i < m; i++) {
			qaws_scalar dist_sq = (qaws_scalar)0.0;
			for (d = 0; d < dim_count; d++) {
				qaws_scalar diff = data[i * dim_count + d]
				                 - data[(i - 1) * dim_count + d];
				dist_sq += diff * diff;
			}
			total_length += scalar_sqrt(dist_sq);
			params[i] = total_length;
		}
		if (total_length > (qaws_scalar)0.0) {
			for (i = 1; i < m; i++)
				params[i] /= total_length;
		}
		params[m - 1] = (qaws_scalar)1.0; /* ensure exact endpoint */
	}

	/* --- Step 2: Knot placement (averaging) ---------------------------- */

	knot_count = n + p + 1;
	knots = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)knot_count);
	if (!knots) { status = QAWS_STATUS_ALLOCATION_FAILURE; goto cleanup; }

	/* First p+1 knots = 0 */
	for (i = 0; i <= p; i++)
		knots[i] = (qaws_scalar)0.0;

	/* Last p+1 knots = 1 */
	for (i = 0; i <= p; i++)
		knots[knot_count - 1 - i] = (qaws_scalar)1.0;

	/* Internal knots via averaging */
	{
		qaws_scalar d_float = (qaws_scalar)(m - 1) / (qaws_scalar)(n - p);
		for (j = 1; j <= n - p - 1; j++) {
			qaws_scalar jd = (qaws_scalar)j * d_float;
			unsigned int ii = (unsigned int)scalar_floor(jd);
			qaws_scalar alpha = jd - (qaws_scalar)ii;
			/* Piegl & Tiller: knot[p+j] = (1-alpha)*t[i-1] + alpha*t[i] */
			knots[p + j] = ((qaws_scalar)1.0 - alpha) * params[ii - 1]
			             + alpha * params[ii];
		}
	}

	/* --- Step 3: Build normal equations N^T N Q = N^T R ---------------- */

	basis = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)(p + 1));
	if (!basis) { status = QAWS_STATUS_ALLOCATION_FAILURE; goto cleanup; }

	if (sys_size > 0) {
		NtN = (qaws_scalar*)calloc(
			(size_t)sys_size * (size_t)sys_size, sizeof(qaws_scalar));
		if (!NtN) { status = QAWS_STATUS_ALLOCATION_FAILURE; goto cleanup; }

		rhs = (qaws_scalar*)calloc(
			(size_t)dim_count * (size_t)sys_size, sizeof(qaws_scalar));
		if (!rhs) { status = QAWS_STATUS_ALLOCATION_FAILURE; goto cleanup; }

		/* For each interior data point k = 1..m-2, compute basis and
		   accumulate into the normal equations. */
		for (k = 1; k < m - 1; k++) {
			unsigned int span = qaws_internal_find_knot_span(
				knots, knot_count, p, n, params[k]);

			qaws_internal_bspline_basis(
				knots, knot_count, p, span, params[k], basis);

			/* Compute R_k = D_k - N_{0,p}(t_k)*D_0 - N_{n-1,p}(t_k)*D_{m-1}
			   for each dimension. */
			for (d = 0; d < dim_count; d++) {
				qaws_scalar r_kd = data[k * dim_count + d];

				/* Subtract endpoint contributions */
				for (j = 0; j <= p; j++) {
					int cp_idx = (int)span - (int)p + (int)j;
					if (cp_idx == 0) {
						r_kd -= basis[j] * data[0 * dim_count + d];
					} else if (cp_idx == (int)(n - 1)) {
						r_kd -= basis[j] * data[(m - 1) * dim_count + d];
					}
				}

				/* Accumulate into rhs: for each interior CP index i in [1..n-2] */
				for (j = 0; j <= p; j++) {
					int cp_idx = (int)span - (int)p + (int)j;
					if (cp_idx >= 1 && cp_idx <= (int)(n - 2)) {
						unsigned int ii = (unsigned int)(cp_idx - 1);
						rhs[d * sys_size + ii] += basis[j] * r_kd;
					}
				}
			}

			/* Accumulate into NtN: for each pair of interior CP indices */
			for (i = 0; i <= p; i++) {
				int ci = (int)span - (int)p + (int)i;
				if (ci < 1 || ci > (int)(n - 2))
					continue;
				for (j = 0; j <= p; j++) {
					int cj = (int)span - (int)p + (int)j;
					if (cj < 1 || cj > (int)(n - 2))
						continue;
					NtN[(unsigned int)(ci - 1) * sys_size + (unsigned int)(cj - 1)]
						+= basis[i] * basis[j];
				}
			}
		}

		/* --- Step 4: Solve ------------------------------------------------- */

		status = solve_linear_system(NtN, rhs, sys_size, dim_count);
		if (status != QAWS_STATUS_OK)
			goto cleanup;
	}

	/* --- Step 5: Assemble control points ------------------------------- */

	ctrl = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (size_t)(n * dim_count));
	if (!ctrl) { status = QAWS_STATUS_ALLOCATION_FAILURE; goto cleanup; }

	/* Q_0 = D_0 (first data point) */
	for (d = 0; d < dim_count; d++)
		ctrl[0 * dim_count + d] = data[0 * dim_count + d];

	/* Q_{n-1} = D_{m-1} (last data point) */
	for (d = 0; d < dim_count; d++)
		ctrl[(n - 1) * dim_count + d] = data[(m - 1) * dim_count + d];

	/* Interior control points from solved system */
	for (i = 0; i < sys_size; i++) {
		for (d = 0; d < dim_count; d++)
			ctrl[(i + 1) * dim_count + d] = rhs[d * sys_size + i];
	}

	/* --- Create B-spline curve ----------------------------------------- */
	{
		qaws_bspline_desc bspline_desc;
		memset(&bspline_desc, 0, sizeof(bspline_desc));
		bspline_desc.dimension = desc->dimension;
		bspline_desc.degree = p;
		bspline_desc.control_points = ctrl;
		bspline_desc.control_point_count = n;
		bspline_desc.knots = knots;
		bspline_desc.knot_count = knot_count;
		bspline_desc.is_uniform = 0;
		bspline_desc.is_closed = 0;

		status = qaws_curve_create_bspline(&bspline_desc, out_curve);
	}

cleanup:
	free(params);
	free(knots);
	free(NtN);
	free(rhs);
	free(basis);
	free(ctrl);
	return status;
}
