#
# This file is part of networkred.
#
# Copyright (C) 2012, Technical University of Denmark
# https://github.com/nefan/networkred
#
# networkred is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# networkred is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with networkred.  If not, see <http://www.gnu.org/licenses/>.
# 

# USE_CUDA = true
# USE_AVX_SIMD = true

# directory for ssparse, klu, umfpack, etc.
SPARSELIBDIR=$(HOME)/alien
KLUINTERNAL=$(HOME)/alien/KLU/Include
KLUUSER=$(HOME)/alien/KLU/User
ATLAS=$(HOME)/alien/atlas/lib
OPENBLAS=/usr/lib/openblas-base
CUDALIB=/usr/local/cuda/lib64

# cljmex
CLJMEXHOME=$(HOME)/projects/cljmex
CLJMEXSRC=$(CLJMEXHOME)/src
CLJMEXINCLUDE=$(CLJMEXHOME)/include

# includes
INCLUDEPATH=$(SPARSELIBDIR)/include $(CLJMEXINCLUDE) 
INC=$(foreach d, $(INCLUDEPATH), -I'$d')

# debugging
MEXFLAGS+=-g
CFLAGS+=-UNDEBUG

# compilers and flags
MATLAB=matlab
MEX=mex -O # -lmwlapack -lmwblas
MEXFLAGS+=-largeArrayDims # matlab sparse stuff

CC=g++
CFLAGS+=-std=c++11
CFLAGS+=-O3
ifdef USE_AVX_SIMD
	CFLAGS+=-DUSE_AVX_SIMD
	CFLAGS+=-mavx # if avx supported
endif
# CFLAGS+=-msse4 # if sse4
CFLAGS+=-fPIC # make it work with matlab
CFLAGS+=-fopenmp # openmp
LDFLAGS=-fPIC # make it work with matlab
LDFLAGS+=-lgomp -fopenmp  # openmp

