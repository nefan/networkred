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

# directories
kronred := kronred
thevenin := thevenin
schur := schur
examples := examples
libdir := ../lib
bindir := ../bin

# library
NETWORKRED_LIB = kronred/kronred.o kronred/kronred_common.o kronred/kronred_factor.o kronred/kronred_symbolic.o

#-----------------------------------------------------------------------------
# Make rules:
#-----------------------------------------------------------------------------
all: $(kronred) $(thevenin) $(schur) lib bin

$(kronred) $(thevenin) $(schur) $(examples): force_look
	$(MAKE) --directory=$@

lib: $(kronred) $(thevenin) $(schur)
	find . -type f -iname \*.mexa64 | xargs -I{} cp {} $(libdir)
	ar cr $(libdir)/libnetworkred.a $(NETWORKRED_LIB)
	ranlib $(libdir)/libnetworkred.a


bin: $(kronred) $(thevenin) $(schur) $(examples)
	cp examples/hello_networkred $(bindir)

clean:
	$(MAKE) --directory=$(kronred) clean
	$(MAKE) --directory=$(thevenin) clean
	$(MAKE) --directory=$(schur) clean
	$(MAKE) --directory=$(examples) clean
	rm $(libdir)/*.mexa64 $(libdir)/*.a $(bindir)/*

force_look:
	@true

.PHONY: all
