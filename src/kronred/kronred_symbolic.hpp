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

#ifndef _KRONRED_SYMBOLIC_H
#define _KRONRED_SYMBOLIC_H

#include "kronred_common.hpp"
#include "kronred.hpp"

// kronred symbolic analysis
void kronred_symbolic(const cljmexComplex_sparse_matrix& M, 
        const int Nnvc, const int cutDegree, const int Nsweeps, 
        const bool split, const bool updateEntireMatrix,
        cljmexComplex_sparse_matrix& Msymbolic,
        int& Nkr, std::vector<int>& index,
        EncodedOps& eOps1, EncodedOps& eOps2);

#endif // _KRONRED_SYMBOLIC_H
