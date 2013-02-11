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

#include "../networkred_common.h"
#include "../kronred/kronred.hpp"
#include "../kronred/kronred_symbolic.hpp"
#include "cljmex.hpp"
#undef printf // mex.h defines printf to be mexPrintf

#include <cstring>
#include <assert.h>


/*
 * Small Kron reduction example
 */

// print 4x4 sparse matrix
void print44Sparse(cljmexInt *Mp, cljmexInt *Mi, double *Mx) {
    const int Mn = 4;
    double Mdense[Mn][Mn];

    for (int c=0; c<Mn; c++)
        for (int r=0; r<Mn; r++)
            Mdense[c][r] = 0;
    for (int c=0; c<Mn; c++)
        for (int ri=Mp[c]; ri<Mp[c+1]; ri++)
            Mdense[c][Mi[ri]] = Mx[ri];

    printf("[ %f %f %f %f\n",Mdense[0][0],Mdense[1][0],Mdense[2][0],Mdense[3][0]);
    printf("  %f %f %f %f\n",Mdense[0][1],Mdense[1][1],Mdense[2][1],Mdense[3][1]);
    printf("  %f %f %f %f\n",Mdense[0][2],Mdense[1][2],Mdense[2][2],Mdense[3][2]);
    printf("  %f %f %f %f ]\n",Mdense[0][3],Mdense[1][3],Mdense[2][3],Mdense[3][3]);
}

// print 3x3 sparse matrix
void print33Sparse(cljmexInt *Mp, cljmexInt *Mi, double *Mx) {
    const int Mn = 3;
    double Mdense[Mn][Mn];

    for (int c=0; c<Mn; c++)
        for (int r=0; r<Mn; r++)
            Mdense[c][r] = 0;
    for (int c=0; c<Mn; c++)
        for (int ri=Mp[c]; ri<Mp[c+1]; ri++)
            Mdense[c][Mi[ri]] = Mx[ri];

    printf("[ %f %f %f\n",Mdense[0][0],Mdense[1][0],Mdense[2][0]);
    printf("  %f %f %f\n",Mdense[0][1],Mdense[1][1],Mdense[2][1]);
    printf("  %f %f %f ]\n",Mdense[0][2],Mdense[1][2],Mdense[2][2]);
}

int main(int argc, char* argv[]) {

    // define sparse 4x4 matrix
    // M = [    1   0.1         0.4 
    //        0.1     2   0.2      
    //              0.2     3   0.3
    //        0.4         0.3     4 ]
    cljmexInt Mp[] = {0,3,6,9,12}; // column indicies
    cljmexInt Mi[] =  {0,1,3,0,1,2,1,2,3,0,2,3};
    double Mx[] =  {1.0,0.1,0.4,0.1,2.0,0.2,0.2,3.0,0.3,0.4,0.3,4};
    double Mz[] =  {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0};
    const int Mn = 4;
    const int Mnz = Mp[Mn];

    // remove node 0 from the network
    //
    // symbolic analysis
    cljmexComplex_sparse_matrix M, Msymbolic;
    int Nkr; 
    std::vector<int> index;
    EncodedOps eOps1, eOps2;

    const int cutDegree = Mn; // reduce all nvc nodes
    const int Nsweeps = 1; // only one sweep necessary with this cutDegree
    const int Nnvc = 1;
    M.p = Mp; M.i = Mi; // setup cljmex sparse matrix
    M.x = Mx; M.z = Mz;
    M.rows = Mn; M.cols = Mn;
    M.nz = Mnz;
    // analyze
    kronred_symbolic(M, Nnvc, cutDegree, Nsweeps, 
        false, true,
        Msymbolic, Nkr, index, eOps1, eOps2);

    // actual reduction
    double *Mxtemp = new double[Msymbolic.nz+1]; // temp arrays
    double *Mztemp = new double[Msymbolic.nz+1]; // +1 for 'trash' location
    // data from M is in Msymbolic already. On subsequent calls to
    // kronred(...), data can be copied into Msymbolic without
    // actually doing the symbolic analysis
    std::memcpy(Mxtemp,Msymbolic.x,Msymbolic.nz*sizeof(double));
    std::memcpy(Mztemp,Msymbolic.z,Msymbolic.nz*sizeof(double));
    kronred(Mxtemp,Mztemp,eOps1);

    // extract result
    cljmexComplex_sparse_matrix Mkr;
    Mkr.p = new cljmexInt[3+1];
    Mkr.i = new cljmexInt[3*3]; // allocate sufficent space
    Mkr.x = new double[3*3]; // allocate sufficent space
    Mkr.z = new double[3*3]; // allocate sufficent space
    // get result from Mxtemp,Mztemp into Mkr
    kronred_extract(Msymbolic,Mxtemp,Mztemp,index,Mkr);
    assert(Mkr.rows == 3);
    assert(Mkr.cols == 3);

    // output
    printf("M before reduction:\n");
    print44Sparse(Mp,Mi,Mx);
    printf("\n");
    printf("M after reduction:\n");
    print33Sparse(Mkr.p,Mkr.i,Mkr.x);
    printf("\n");

    // clean up
    delete Mxtemp;
    delete Mztemp;
    delete Msymbolic.p;
    delete Msymbolic.i;
    delete Msymbolic.x;
    delete Msymbolic.z;
    delete Mkr.p;
    delete Mkr.i;
    delete Mkr.x;
    delete Mkr.z;

    return 0;
}
