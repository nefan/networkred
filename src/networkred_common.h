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

// MEX debugging
#undef NDEBUG

// convenient macros
#define CLJ_CS_TO_CS_CL(cs_s,clj_s) \
    cs_cl cs_s; \
    cs_s.m = clj_s.rows; cs_s.n = clj_s.cols; \
    cs_s.p = clj_s.p; cs_s.i = clj_s.i; \
    cs_s.x = clj_s.v; \
    cs_s.nzmax = L.nz; cs_s.nz = -1; \

