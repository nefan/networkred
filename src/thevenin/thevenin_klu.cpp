/*
 * This file is part of networkred.
 *
 * Copyright (c) 2006, Timothy A. Davis.
 * File derived from:
 * KLU
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

#include "thevenin_klu.hpp"
#include "thevenin_klu.cljmex.hpp"

#include "klu.h"
#include "klu_cholmod.h"
#define Long SuiteSparse_long

#include <omp.h>

cljmex_start()

    double *X, *B, *Xz, *Xx, *Bx, *Bz, *Lx, *Ux, *Rs, *Offx, *Wx,
        *Uz, *Lz, *Offz, *Wz, *W, *Xi, *Bi ;
    Long *Lp, *Li, *Up, *Ui, *P, *Q, *R, *Rp, *Ri, *Offp, *Offi ;
    mxArray *L_matlab, *U_matlab, *p_matlab, *q_matlab, *R_matlab, *F_matlab,
        *r_matlab, *field ;
    const mxArray *B_matlab = NULL, *opts_matlab ;
    klu_l_symbolic *Symbolic ;
    klu_l_numeric *Numeric ;
    klu_l_common Common ;
    Long k, nrhs = 0, symmetric,
        do_transpose = 0, p, pend, nblocks,
        R1 [2], chunk, nr, i, j, block, k1, k2, nk, bn = 0, ordering ;
    static const char *LUnames [ ] = { "L", "U", "p", "q", "R", "F", "r" } ;

    Long n = A.rows;
    Long nz = A.nz;

    // Set KLU options

    klu_l_defaults (&Common) ;

    /* memory management routines */
    Common.malloc_memory  = mxMalloc ;
    Common.calloc_memory  = mxCalloc ;
    Common.free_memory    = mxFree ;
    Common.realloc_memory = mxRealloc ;

    if (Common.ordering < 0 || Common.ordering > 4)
    {
        mexErrMsgTxt ("invalid ordering option") ;
    }
    ordering = Common.ordering ;

    /* ordering option 3,4 becomes KLU option 3, with symmetric 0 or 1 */
    symmetric = (Common.ordering == 4) ;
    if (symmetric) Common.ordering = 3 ;
    Common.user_order = klu_l_cholmod ;
    Common.user_data = &symmetric ;
    Common.tol = 0; // no pivoting

    // Factorize

    double t0symbolic = omp_get_wtime();
    Symbolic = klu_l_analyze (n, A.p, A.i, &Common) ;
    double t1symbolic = omp_get_wtime();
    if (Symbolic == (klu_l_symbolic *) NULL)
    {
        mexErrMsgTxt ("klu symbolic analysis failed") ;
    }
    // mexPrintf("KLU: symbolic factorization took %f ms.\n",(t1symbolic-t0symbolic)/1e-3);

    /* ------------------------------------------------------------------ */
    /* factorize */
    /* ------------------------------------------------------------------ */

    double t0factor = omp_get_wtime();
    for (int i=0; i<NITER; i++)
        Numeric = klu_zl_factor (A.p, A.i, (double*)A.v, Symbolic, &Common) ;
    double t1factor = omp_get_wtime();
    double time = (t1factor-t0factor)/NITER;
    // mexPrintf("KLU: factorization took %f ms.\n",time/1e-3);
    if (Common.status != KLU_OK)
    {
        mexErrMsgTxt ("klu numeric factorization failed") ;
    }

    // LU = klu (A) usage; extract factorization

    /* sort the row indices in each column of L and U */
    klu_zl_sort (Symbolic, Numeric, &Common) ;

    /* L */
    L_matlab = mxCreateSparse (n, n, Numeric->lnz,
            mxCOMPLEX) ;
    Lp = (Long *) mxGetJc (L_matlab) ;
    Li = (Long *) mxGetIr (L_matlab) ;
    Lx = mxGetPr (L_matlab) ;
    Lz = mxGetPi (L_matlab) ;

    /* U */
    U_matlab = mxCreateSparse (n, n, Numeric->unz,
            mxCOMPLEX) ;
    Up = (Long *) mxGetJc (U_matlab) ;
    Ui = (Long *) mxGetIr (U_matlab) ;
    Ux = mxGetPr (U_matlab) ;
    Uz = mxGetPi (U_matlab) ;

    /* p */
    p_matlab = mxCreateNumericMatrix (1, n, mx_int, mxREAL) ;
    P = (Long *) mxGetData (p_matlab) ;

    /* q */
    q_matlab = mxCreateNumericMatrix (1, n, mx_int, mxREAL) ;
    Q = (Long *) mxGetData (q_matlab) ;

    /* R, as a sparse diagonal matrix */
    R_matlab = mxCreateSparse (n, n, n+1, mxREAL) ;
    Rp = (Long *) mxGetJc (R_matlab) ;
    Ri = (Long *) mxGetIr (R_matlab) ;
    Rs = mxGetPr (R_matlab) ;
    for (k = 0 ; k <= n ; k++)
    {
        Rp [k] = k ;
        Ri [k] = k ;
    }

    /* F, off diagonal blocks */
    F_matlab = mxCreateSparse (n, n, Numeric->nzoff,
            mxCOMPLEX) ;
    Offp = (Long *) mxGetJc (F_matlab) ;
    Offi = (Long *) mxGetIr (F_matlab) ;
    Offx = mxGetPr (F_matlab) ;
    Offz = mxGetPi (F_matlab) ;

    /* r, block boundaries */
    nblocks = Symbolic->nblocks ;
    r_matlab = mxCreateNumericMatrix (1, nblocks+1, mx_int, mxREAL) ;
    R = (Long *) mxGetData (r_matlab) ;

    /* extract the LU factorization from KLU Numeric and Symbolic objects */
    klu_zl_extract (Numeric, Symbolic, Lp, Li, Lx, Lz, Up, Ui, Ux, Uz,
            Offp, Offi, Offx, Offz, P, Q, Rs, R, &Common) ;

    /* fix p and q for 1-based indexing */
    for (k = 0 ; k < n ; k++)
    {
        P [k]++ ;
        Q [k]++ ;
    }

    /* fix r for 1-based indexing */
    for (k = 0 ; k <= nblocks ; k++)
    {
        R [k]++ ;
    }

    // Output

    /* create output LU struct */
    pargout [0] = mxCreateStructMatrix (1, 1, 7, LUnames) ;
    mxSetFieldByNumber (pargout [0], 0, 0, L_matlab) ;
    mxSetFieldByNumber (pargout [0], 0, 1, U_matlab) ;
    mxSetFieldByNumber (pargout [0], 0, 2, p_matlab) ;
    mxSetFieldByNumber (pargout [0], 0, 3, q_matlab) ;
    mxSetFieldByNumber (pargout [0], 0, 4, R_matlab) ;
    mxSetFieldByNumber (pargout [0], 0, 5, F_matlab) ;
    mxSetFieldByNumber (pargout [0], 0, 6, r_matlab) ;

    // Clean up

    klu_l_free_symbolic (&Symbolic, &Common) ;
    klu_l_free_numeric (&Numeric, &Common) ;

cljmex_end()
