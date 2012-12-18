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

SOURCES=hello_networkred.cpp
HEADERS=../networkred_common.h
OBJECTS=$(SOURCES:.cpp=.o)

# libraries
libdir := ../../lib
CLIBS+=$(CXXLIBS)

CFLAGS+=-DNOMEX
CXXFLAGS+=$(CFLAGS)

# INCLUDEPATH+=
INCLUDEPATH+=..
INC=$(foreach d, $(INCLUDEPATH), -I'$d')

#-----------------------------------------------------------------------------
# Make rules:
#-----------------------------------------------------------------------------
all: $(OBJECTS) hello_networkred

%.o: %.c
	$(CC) -c $(CFLAGS) $(INC) $*.c -o $@
%.o: %.cpp
	$(CC) -c $(CXXFLAGS) $(INC) $*.cpp -o $@


$(OBJECTS): $(HEADERS) Makefile

hello_networkred: hello_networkred.o $(libdir)/libnetworkred.a
	$(CC) hello_networkred.o -L$(libdir) -lnetworkred -o hello_networkred

clean:
	rm -f ${OBJECTS} hello_networkred
