# -- compilation flags

FC = gfortran

FFLAGS = $(STANDARD) $(WARNINGS) $(DEBUG) $(OPTIMISE) $(LIBRARIES)

STANDARD  = -std=f2003
WARNINGS  = $(F1) $(F2) $(F3) $(F4)
DEBUG     = $(F5)
OPTIMISE  = $(F6)
LIBRARIES = $(F7)
LINKING	  =

F1	= -Wall -pedantic -Wextra -Waliasing -Wsurprising -Wfunction-elimination
F2	= -fmodule-private
F3	= -Wimplicit-procedure -Wimplicit-interface
F4	= -Warray-bounds -Warray-temporaries -Wrealloc-lhs-all -fbounds-check
F5	= -g -fcheck=all -fdump-tree-original
F6	= -O2
F7	= -fall-intrinsics -ffree-line-length-none -fopenmp

# -- source files
object = rsg.o \
		 input_data.o target_data.o \
		 gridr.o sturmian.o \
		 V12me.o \
		 wigner.o vnc_module.o vnc.o \
		 state_class.o states_module.o \
		 form.o hartree_fock.o \
		 intp.o core_wavefunctions.o \
		 ovlpste1me_module.o one_electron_func_module.o one_electron_func.o \
		 main.o

# -- commands
make:
	make main

main: $(object)
	$(FC) $(FFLAGS) -o $@ $^

all: $(object)
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(F4) $(F5) $(F6) $(F7) -c -o $@ $<

%.o: %.f
	$(FC) $(F4) $(F5) $(F7) -c -o $@ $<

clean:
	rm -f *.o *.mod *.d *.original
