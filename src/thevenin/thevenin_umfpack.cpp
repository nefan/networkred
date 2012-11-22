/*
 * This file is part of networkred.
 *
 * Copyright (c) 2006, Timothy A. Davis.
 * File derived from:
 * UMFPACK
 * http://www.suitesparse.com
 *
 * Modified by Stefan Sommer (shso@elektro.dtu.dk), November 2012
 *
 * networkred is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * networkred is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with networkred.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "thevenin_umfpack.hpp"

#include "umfpack.h"
#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <math.h>
#include <float.h>

#include <omp.h>

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MATCH(s1,s2) (strcmp ((s1), (s2)) == 0)
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif
#define Long SuiteSparse_long


/* ========================================================================== */
/* === UMFPACK ============================================================== */
/* ========================================================================== */

void mexFunction
(
    int nargout,		/* number of outputs */
    mxArray *pargout [ ],	/* output arguments */
    int nargin,			/* number of inputs */
    const mxArray *pargin [ ]	/* input arguments */
)
{

    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL], dx, dz, dexp ;
    double *Lx, *Lz, *Ux, *Uz, *Ax, *Az, *Bx, *Bz, *Xx, *Xz, *User_Control,
	*p, *q, *Out_Info, *p1, *p2, *p3, *p4, *Ltx, *Ltz, *Rs, *Px, *Qx ;
    void *Symbolic, *Numeric ;
    Long *Lp, *Li, *Up, *Ui, *Ap, *Ai, *P, *Q, do_solve, lnz, unz, nn, i,
	transpose, size, do_info, do_numeric, *Front_npivcol, op, k, *Rp, *Ri,
	*Front_parent, *Chain_start, *Chain_maxrows, *Chain_maxcols, nz, status,
	nfronts, nchains, *Ltp, *Ltj, *Qinit, print_level, status2, no_scale,
	*Front_1strow, *Front_leftmostdesc, n_row, n_col, n_inner, sys,
	ignore1, ignore2, ignore3, A_is_complex, B_is_complex, X_is_complex,
	*Pp, *Pi, *Qp, *Qi, do_recip, do_det ;
    mxArray *Amatrix, *Bmatrix, *User_Control_struct, *User_Qinit ;
    char *xoperator, *operation ;
    mxComplexity Atype, Xtype ;
    char warning [200] ;
    int info_details ;

    /* ---------------------------------------------------------------------- */
    /* define the memory manager and printf functions for UMFPACK and AMD */ 
    /* ---------------------------------------------------------------------- */

    /* with these settings, the UMFPACK mexFunction can use ../Lib/libumfpack.a
     * and ../Lib/libamd.a, instead compiling UMFPACK and AMD specifically for
     * the MATLAB mexFunction. */

    amd_malloc = mxMalloc ;
    amd_free = mxFree ;
    amd_calloc = mxCalloc ;
    amd_realloc = mxRealloc ;

    amd_printf = mexPrintf ;

    User_Control_struct = (mxArray *) NULL ;
    User_Qinit = (mxArray *) NULL ;

    do_info = 0 ;
    do_solve = FALSE ;
    do_numeric = TRUE ;
    transpose = FALSE ;
    no_scale = FALSE ;
    do_det = FALSE ;

	/* ------------------------------------------------------------------ */
	/* LU factorization */
	/* ------------------------------------------------------------------ */

	/*

	    with no scaling:
	*/

	Amatrix = (mxArray *) pargin [0] ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    if (nargin != 2)
    {
        mexErrMsgTxt (
            "Usage: [L,U,P,Q,timing] = thevenin_umfpack(A,ordering) ") ;
    }

    n_row = mxGetM (Amatrix) ;
    n_col = mxGetN (Amatrix) ;
    nn = MAX (n_row, n_col) ;
    n_inner = MIN (n_row, n_col) ;

    if (!mxIsSparse (Amatrix))
    {
	mexErrMsgTxt ("input matrix A must be sparse") ;
    }
    Atype = mxCOMPLEX;

    int uordering = (int)(mxGetScalar(pargin[1]));
    mxAssert(uordering == 0 || uordering == 1, "invalid ordering option\n");
    int ordering;
    if (uordering == 0) {
        mexPrintf("umfpack with no ordering\n");
        ordering = UMFPACK_ORDERING_NONE;
    } else {
        mexPrintf("umfpack with CHOLMOD ordering\n");
        ordering = UMFPACK_ORDERING_CHOLMOD;
    }

    /* The real/complex status of A determines which version to use, */
    /* (umfpack_dl_* or umfpack_zl_*). */
    A_is_complex = mxIsComplex (Amatrix) ;
    mxAssert(A_is_complex,"go complex, dude :-)");
    Ap = (Long *) mxGetJc (Amatrix) ;
    Ai = (Long *) mxGetIr (Amatrix) ;
    Ax = mxGetPr (Amatrix) ;
    Az = mxGetPi (Amatrix) ;

    /* ---------------------------------------------------------------------- */
    /* set the Control parameters */
    /* ---------------------------------------------------------------------- */

	umfpack_zl_defaults (Control) ;

	/* turn off scaling for [L, U, P, Q] = umfpack (A) ;
	 * ignoring the input value of Control (24) for the usage
	 * [L, U, P, Q] = umfpack (A, Control) ; */
	Control [UMFPACK_SCALE] = UMFPACK_SCALE_NONE ;
	Control [UMFPACK_ORDERING] = UMFPACK_ORDERING_NONE ;

    print_level = (Long) Control [UMFPACK_PRL] ;

	Qinit = (Long *) NULL ;

    /* ---------------------------------------------------------------------- */
    /* perform the symbolic factorization */
    /* ---------------------------------------------------------------------- */

    double t0symbolic = omp_get_wtime();
	status = umfpack_zl_qsymbolic (n_row, n_col, Ap, Ai, Ax, Az,
	    Qinit, &Symbolic, Control, Info) ;
    double t1symbolic = omp_get_wtime();
    mxAssert(status == UMFPACK_OK,"symbolic factorization failed");
    // mexPrintf("UMFPACK: symbolic factorization took %f ms.\n",(t1symbolic-t0symbolic)/1e-3);

    /* ---------------------------------------------------------------------- */
    /* report the Symbolic object */
    /* ---------------------------------------------------------------------- */

	(void) umfpack_zl_report_symbolic (Symbolic, Control) ;

    /* ---------------------------------------------------------------------- */
    /* perform numeric factorization, or just return symbolic factorization */
    /* ---------------------------------------------------------------------- */

    double t0factor = omp_get_wtime();
    for (int i=0; i<10; i++)
        status = umfpack_zl_numeric (Ap, Ai, Ax, Az, Symbolic, &Numeric,
                Control, Info) ;
    double t1factor = omp_get_wtime();
    double time = (t1factor-t0factor)/10;
    // mexPrintf("UMFPACK: factorization took %f ms.\n",time/1e-3);

	/* ------------------------------------------------------------------ */
	/* free the symbolic factorization */
	/* ------------------------------------------------------------------ */

    umfpack_zl_free_symbolic (&Symbolic) ;

	/* ------------------------------------------------------------------ */
	/* report the Numeric object */
	/* ------------------------------------------------------------------ */

	mxAssert(status == UMFPACK_OK,"numeric factorization failed");

    (void) umfpack_zl_report_numeric (Numeric, Control) ;

	/* ------------------------------------------------------------------ */
	/* return the solution, determinant, or the factorization */
	/* ------------------------------------------------------------------ */

    /* -------------------------------------------------------------- */
    /* get L, U, P, Q, and r */
    /* -------------------------------------------------------------- */

    status = umfpack_zl_get_lunz (&lnz, &unz, &ignore1, &ignore2,
            &ignore3, Numeric) ;

	    /* avoid malloc of zero-sized arrays */
	    lnz = MAX (lnz, 1) ;
	    unz = MAX (unz, 1) ;

	    /* get temporary space, for the *** ROW *** form of L */
	    Ltp = (Long *) mxMalloc ((n_row+1) * sizeof (Long)) ;
	    Ltj = (Long *) mxMalloc (lnz * sizeof (Long)) ;
	    Ltx = (double *) mxMalloc (lnz * sizeof (double)) ;
        Ltz = (double *) mxMalloc (lnz * sizeof (double)) ;

	    /* create permanent copy of the output matrix U */
	    pargout [1] = mxCreateSparse (n_inner, n_col, unz, Atype) ;
	    Up = (Long *) mxGetJc (pargout [1]) ;
	    Ui = (Long *) mxGetIr (pargout [1]) ;
	    Ux = mxGetPr (pargout [1]) ;
	    Uz = mxGetPi (pargout [1]) ;

	    /* temporary space for the integer permutation vectors */
	    P = (Long *) mxMalloc (n_row * sizeof (Long)) ;
	    Q = (Long *) mxMalloc (n_col * sizeof (Long)) ;

	    /* get Lt, U, P, Q, and Rs from the numeric object */
		status = umfpack_zl_get_numeric (Ltp, Ltj, Ltx, Ltz, Up, Ui, Ux,
		    Uz, P, Q, (double *) NULL, (double *) NULL,
		    &do_recip, Rs, Numeric) ;
		umfpack_zl_free_numeric (&Numeric) ;

	    /* create sparse permutation matrix for P */
	    pargout [2] = mxCreateSparse (n_row, n_row, n_row, mxREAL) ;
	    Pp = (Long *) mxGetJc (pargout [2]) ;
	    Pi = (Long *) mxGetIr (pargout [2]) ;
	    Px = mxGetPr (pargout [2]) ;
	    for (k = 0 ; k < n_row ; k++)
	    {
		Pp [k] = k ;
		Px [k] = 1 ;
		Pi [P [k]] = k ;
	    }
	    Pp [n_row] = n_row ;

	    /* create sparse permutation matrix for Q */
	    pargout [3] = mxCreateSparse (n_col, n_col, n_col, mxREAL) ;
	    Qp = (Long *) mxGetJc (pargout [3]) ;
	    Qi = (Long *) mxGetIr (pargout [3]) ;
	    Qx = mxGetPr (pargout [3]) ;
	    for (k = 0 ; k < n_col ; k++)
	    {
		Qp [k] = k ;
		Qx [k] = 1 ;
		Qi [k] = Q [k] ;
	    }
	    Qp [n_col] = n_col ;

	    /* permanent copy of L */
	    pargout [0] = mxCreateSparse (n_row, n_inner, lnz, Atype) ;
	    Lp = (Long *) mxGetJc (pargout [0]) ;
	    Li = (Long *) mxGetIr (pargout [0]) ;
	    Lx = mxGetPr (pargout [0]) ;
	    Lz = mxGetPi (pargout [0]) ;

	    /* convert L from row form to column form */
		/* non-conjugate array transpose */
	        status = umfpack_zl_transpose (n_inner, n_row, Ltp, Ltj, Ltx,
		    Ltz, (Long *) NULL, (Long *) NULL, Lp, Li, Lx, Lz,
		    FALSE) ;

        /* timing */
        pargout[4] = mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(pargout[4]) = time;

	    mxFree (Ltp) ;
	    mxFree (Ltj) ;
	    mxFree (Ltx) ;
	    if (Ltz) mxFree (Ltz) ;

}
