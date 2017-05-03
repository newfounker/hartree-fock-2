# -- compilation flags

FC := gfortran

FFLAGS := \
	# -Wall -pedantic -Wextra -Waliasing -Wsurprising -Wfunction-elimination \
	-pedantic -Wfunction-elimination -Waliasing \
	-Warray-bounds -Warray-temporaries -Wrealloc-lhs-all -fbounds-check \
	-Wimplicit-procedure -Wimplicit-interface \
	-g -fcheck=all -fdump-tree-original \
	-fmodule-private -fall-intrinsics -ffree-line-length-none -fno-align-commons \
	-fopenmp -O2

# -- object files
OBJ := rsg.o \
		 input_data.o target_data.o \
		 gridr.o sturmian.o \
		 V12me.o \
		 wigner.o vnc_module.o vnc.o \
		 state_class.o states_module.o \
		 form.o hartree_fock.o \
		 intp.o core_wavefunctions.o \
		 ovlpste1me_module.o one_electron_func_module.o one_electron_func.o \
		 main.o

# -- make commands
make:
	make main

main: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

all: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.mod *.d *.original
