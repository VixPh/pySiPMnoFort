CC := gcc
FC := gfortran
F2PY := python3 -m numpy.f2py


CF := "-Ofast -fomit-frame-pointer -march=native"

FF := --opt="-Ofast -std=f2008 -fomit-frame-pointer -march=native -mavx -mavx2"


SRCDIR := src/
LIBDIR = libs/


MODULE := FortranFunctions


export CFLAGS="$(CF)"


.PHONY : build clear
all: build | clear


build:
	@echo "Building $(MODULE) module"
	$(F2PY) -c $(FF)  --no-wrap-functions -m $(MODULE) $(SRCDIR)FortranFunctions.f90

clear:
	@mv $(MODULE).* $(LIBDIR)
