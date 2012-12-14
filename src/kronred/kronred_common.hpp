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

#include "networkred_common.h"
#include "cljmex.hpp"
#include <vector>

extern int linear_index( int row, int col, long *Mp, long *Mi);

int kronred_extract(const cljmexComplex_sparse_matrix& M, 
        double *Mx, double *Mz,
        std::vector<int>& index,
        const cljmexComplex_sparse_matrix& Mkr);
