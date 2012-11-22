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

include ../../networkred_config.mk

SOURCES=kronred.cpp
CLJMEX=kronred.mexa64
HEADERS=kronred.hpp
OBJECTS=$(SOURCES:.cpp=.mexa64) kronred_common.o kronred_factor.o

ifdef USE_AVX_SIMD
	CFLAGS+=-DUSE_AVX_SIMD
endif

# libraries
CLIBS+=$(CXXLIBS)

INCLUDEPATH+=$(KLUINTERNAL) /usr/local/MATLAB/R2012a/extern/include/ ..
INC=$(foreach d, $(INCLUDEPATH), -I'$d')

#-----------------------------------------------------------------------------
# Make rules:
#-----------------------------------------------------------------------------
all: $(OBJECTS)

%.o: %.c
	$(CC) -c $(CFLAGS) $(INC) $*.c -o $@
%.o: %.cpp
	$(CC) -c $(CFLAGS) $(INC) $*.cpp -o $@
%.s: %.c
	$(CC) -S $(CFLAGS) -fverbose-asm $(INC) $*.c -o $@


%.cljmex.hpp: %.cljmex
	clojure -cp $(CLJMEXSRC) $*.cljmex > $@
$(CLJMEX): %.mexa64: %.cpp %.cljmex.hpp

$(OBJECTS): $(HEADERS) Makefile

%.mexa64: %.cpp kronred_factor.o kronred_common.o
	$(MEX) $(MEXFLAGS) CFLAGS='$$CFLAGS $(CFLAGS)' LDFLAGS='$$LDFLAGS $(LDFLAGS)' CXXLIBS='$$CXXLIBS $(CXXLIBS)' $(INC) $(@:.mexa64=.cpp) kronred_factor.o kronred_common.o

clean:
	rm -f ${OBJECTS} *.cljmex.hpp *.s
