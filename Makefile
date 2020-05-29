CC := gcc
CF := "-Ofast -funroll-all-loops -ffast-math -fchecking -pipe"
FF := --opt="-Ofast -funroll-all-loops -ffast-math -fopt-info-optall -fchecking -save-temps -dA"
F2PY := f2py3.6

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
