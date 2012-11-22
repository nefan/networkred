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

SOURCES=
MEXSOURCES=thevenin_getz.cpp thevenin_reach.cpp thevenin_umfpack.cpp thevenin_klu.cpp
CLJMEX=thevenin_getz.mexa64 thevenin_reach.mexa64 thevenin_klu.mexa64
HEADERS=thevenin_umfpack.hpp thevenin_klu.hpp 
MEXOBJECTS=$(MEXSOURCES:.cpp=.mexa64)
OBJECTS=$(SOURCES:.cpp=.o)

EXTERNAL=$(KLUUSER)/klu_l_cholmod.c
EXTERNALOBJECTS=$(EXTERNAL:.c=.o)

LDFLAGS+=-L$(HOME)/alien/lib -L$(HOME)/alien/KLU/User -L$(ATLAS) -L$(CUDALIB) 
CXXLIBS+=-lklu
CXXLIBS+=-lumfpack -lcholmod -lsuitesparseconfig -lamd -lcolamd -lbtf -lrt
CXXLIBS+=-lsatlas -lf77blas -lgfortran
CXXLIBS+=-lssparse -lcxsparse 
CLIBS+=$(CXXLIBS)

INCLUDEPATH+=$(KLUUSER) ..
INC=$(foreach d, $(INCLUDEPATH), -I'$d')

CXXFLAGS+=$(CFLAGS)

#-----------------------------------------------------------------------------
# Make rules:
#-----------------------------------------------------------------------------
all: $(OBJECTS) $(MEXOBJECTS)

$(EXTERNALOBJECTS): $(EXTERNAL)
%.o: %.c
	$(CC) -c -fpermissive $(CFLAGS) $(INC) $(@:.o=.c) -o $@
%.o: %.cpp
	$(CC) -c -fpermissive $(CFLAGS) $(CXXFLAGS) $(INC) $(@:.o=.cpp) -o $@

$(OBJECTS): $(HEADERS) $(EXTERNALOBJECTS)
$(MEXOBJECTS): $(HEADERS) $(EXTERNALOBJECTS)

%.cljmex.hpp: %.cljmex
	clojure -cp $(CLJMEXSRC) $*.cljmex > $@
$(CLJMEX): %.mexa64: %.cpp %.cljmex.hpp

%.mexa64: %.c
	$(MEX) $(MEXFLAGS) CFLAGS='$$CFLAGS $(CFLAGS)' LDFLAGS='$$LDFLAGS $(LDFLAGS)' CLIBS='$$CLIBS $(CLIBS)' $(INC) $(@:.mexa64=.c)
%.mexa64: %.cpp $(EXTERNALOBJECTS)
	$(MEX) $(MEXFLAGS) CFLAGS='$$CFLAGS $(CFLAGS)' CXXFLAGS='$$CXXFLAGS $(CXXFLAGS)' LDFLAGS='$$LDFLAGS $(LDFLAGS)' CXXLIBS='$$CXXLIBS $(CXXLIBS)' $(INC) $(@:.mexa64=.cpp) $(EXTERNALOBJECTS)

clean:
	rm -f ${OBJECTS} ${MEXOBJECTS} ${EXTERNALOBJECTS} *.cljmex.hpp
