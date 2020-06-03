CC := gcc
FC := gfortran
F2PY := f2py3.6


CF := "-Ofast -funroll-all-loops -fchecking"

FF := --opt="-Ofast -funroll-all-loops -fopt-info -fchecking \
-save-temps -dA -std=f2008"

GFF := -Ofast -funroll-all-loops -ffast-math -fchecking


SRCDIR := src/
LIBDIR = libs/


MODULE := FortranFunctions


export CFLAGS="$(CF)"


.PHONY : build clear
all: build | clear


build:
	@echo "Building $(MODULE) module"
	$(F2PY) -c $(FF) -m $(MODULE) $(SRCDIR)FortranFunctions.f90

clear:
	@mv $(MODULE).* $(LIBDIR)$
