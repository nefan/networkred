/*
 * This file is part of networkred.
 *
 * Copyright (C) 2012, Technical University of Denmark
 * https://github.com/nefan/networkred
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
#include "thevenin_getz.cljmex.hpp"

#define CS_LONG
#define CS_COMPLEX
#include "cs.h"
#define Long SuiteSparse_long
#include "ssparse.h"

#include <omp.h>
#include <fstream>
#include <string>

cljmex_start()

    // Setup data structures

    Long *P = (Long*)mxMalloc(n * sizeof(Long));
    for (int i=0; i<n; i++)
        P[i] = Pone.x[i]-1;
    Long *Pinv = cs_pinv(P,n);
    mxAssert(Pinv != NULL, "pinv allocation failed");

    Long *Q = (Long*)mxMalloc(n * sizeof(Long));
    for (int i=0; i<n; i++)
        Q[i] = Qone.x[i]-1;
    Long *Qinv = cs_pinv(Q,n);
    mxAssert(Qinv != NULL, "qinv allocation failed");

    CLJ_CS_TO_CS_CL(Lcs,L);
    CLJ_CS_TO_CS_CL(Ucs,U);

    cs_complex_t *B = (cs_complex_t*)mxMalloc (Link.nz * sizeof (cs_complex_t)) ;
    Long *Bi = (Long*)mxMalloc (Link.nz * sizeof (Long)) ;
    for (int c=0; c<M; c++) {
        for (int k=Link.p[c]; k<Link.p[c+1]; k++) {
            B[k] = Link.v[k]/R.x[Pinv[Link.i[k]]];
            Bi[k] = Pinv[Link.i[k]];
        }
    }
    CLJ_CS_TO_CS_CL(Bcs,Link);
    Bcs.i = Bi;
    Bcs.x = B;

    Long *LinkTQi = (Long*)mxMalloc (LinkT.nz * sizeof (Long)) ;
    for (int c=0; c<M; c++) {
        for (int k=LinkT.p[c]; k<LinkT.p[c+1]; k++) {
            LinkTQi[k] = Qinv[LinkT.i[k]];
        }
    }
    CLJ_CS_TO_CS_CL(LinkTcs,LinkT);
    LinkTcs.i = LinkTQi;

    CS_INT *reachL = (CS_INT*)mxMalloc ((1 + n) * M * sizeof (CS_INT)) ;
    for (int i=0; i<(1+n)*M; i++)
        reachL[i] = (CS_INT)reachL_double.x[i];
    CS_INT *reachU = (CS_INT*)mxMalloc ((1 + n) * M * sizeof (CS_INT)) ;
    for (int i=0; i<(1+n)*M; i++)
        reachU[i] = (CS_INT)reachU_double.x[i];

    
    // Global results, workspace, etc.

    cs_complex_t *vs = (cs_complex_t*)mxMalloc (M * sizeof (cs_complex_t)) ;
    cs_complex_t *S = (cs_complex_t*)mxMalloc (M * M * sizeof (cs_complex_t)) ;
    // zero
    for (int i=0; i<M; i++)
        vs[i] = 0;
    for (int i=0; i<M*M; i++)
        S[i] = 0;

    cs_complex_t *global_x = (cs_complex_t*)mxMalloc(n * sizeof(cs_complex_t) * M);
    cs_complex_t *global_y = (cs_complex_t*)mxMalloc(n * sizeof(cs_complex_t) * M);
    Long *global_topx = (Long*)mxMalloc(M * sizeof(Long));
    Long *global_topy = (Long*)mxMalloc(M * sizeof(Long));

    // zero
    for (int i=0; i<n*M; i++)
        global_x[i] = global_y[i] = 0;

    // Compute

    double timebacksolve = 0;
    double timeips = 0;

    for (int iter=0; iter<NITER; iter++) {

        double t0backsolve = omp_get_wtime();

#pragma omp parallel for schedule(static)
        for (int K=0; K<M; K++) {

            cs_complex_t *x = &global_x[n*K];
            cs_complex_t *y = &global_y[n*K];
            Long *topx = &global_topx[K];
            Long *topy = &global_topy[K];

            // do backwards solve
            Long *topyi = reachL+(1+n)*K;
            Long *yi = topyi+1;
            Long *topxi = reachU+(1+n)*K;
            Long *xi = topxi+1;

            *topy = ssp_cl_splsolve(&Lcs,&Bcs,K,topyi,y,NULL);
            *topx = ssp_cl_sputsolve(&Ucs,&LinkTcs,K,topxi,x,NULL);
        }

        double t1backsolve = omp_get_wtime();
        timebacksolve += (t1backsolve-t0backsolve)/NITER;

        double t0ips = omp_get_wtime();

#pragma omp parallel for schedule(static)
        for (int c=0; c<M; c++) {

            int start_r, end_r;
            if (diagOnly) {
                start_r = c;
                end_r = c+1;
            } else {
                start_r = 0;
                end_r = M;
            }

            for (int r=start_r; r<end_r; r++) {
                // get data
                cs_complex_t *x = &global_x[n*r];
                cs_complex_t *y = &global_y[n*c];
                Long *topx = &global_topx[r];

                Long *topxi = reachU+(1+n)*r;
                Long *xi = topxi+1;
                Long *yi = &global_topy[c]+1;

                // inner product
                CS_ENTRY ip = ssp_spdot(xi+*topx,x,y,n-*topx);

                // Schur complement
                cs_complex_t Urc = 0;
                for (int i=D.p[c]; i<D.p[c+1]; i++) {
                    if (D.i[i] == r) {
                        Urc = D.v[i];
                        break;
                    }
                }

                cs_complex_t Src = Urc-ip;
                // store value
                S[M*c+r] = Src;

                // Thevenin imp
                if (c == r) {
                    cs_complex_t v = cs_complex_t(1,0)/Src;
                    vs[c] = v;
                }
            }

        }

        double t1ips = omp_get_wtime();
        timeips += (t1ips-t0ips)/NITER;

    }

    // Clean up
    cs_free(Pinv);
    cs_free(Qinv);

cljmex_end()
