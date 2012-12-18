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

#ifndef _KRONRED_H
#define _KRONRED_H

#include "networkred_common.h"

// operations
#include <vector>
#include <array>
typedef std::pair<int,int> Op;
typedef std::vector<Op> Ops;
#define NoOp Op(0,-1)
#define NewOp(node1,node2) Op(node1,node2)
#define StartNodeOp(node) Op(0,node)

typedef std::array<int,4> EncodedOp;
typedef std::vector<EncodedOp> EncodedOps;

// kron reduction
void kronred(double *Mx, double *Mz, EncodedOps& ops);

#endif // _KRONRED_H
