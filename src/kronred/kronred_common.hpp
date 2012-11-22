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

#include "kronred.hpp"


extern int linear_index( int row, int col, long *Mp, long *Mi);

extern int kronred_analyze (
    long *Mp,
    long *Mi,
    double *Mx,
    double *Mz,
    int nz,
    double *Mmx,
    double *Mmz,
    double *ops1,
    int Nops,
    int *ops2
);

// extract from M to Mkr
extern int kronred_extract
(
    double *Mmx,
    double *Mmz,
    long *Mp,
    long *Mi,
    int Nkr,
    double *index,
    double *filter,
    long *Mkrp,
    long *Mkri,
    double *Mkrx,
    double *Mkrz
);
