# -- unoptimised compilation
	# -Warray-bounds -Warray-temporaries -Wrealloc-lhs-all -fbounds-check \
	# -Wimplicit-procedure -Wimplicit-interface \
	# -fopenmp -O3 -ffast-math \
	# -fall-intrinsics -ffree-line-length-none
	# -Wall -pedantic -Wextra -Waliasing -Wsurprising -Wfunction-elimination \
	# -fmodule-private -fall-intrinsics -ffree-line-length-none
	# -g -fcheck=all -fdump-tree-original \

# -- optimised ifortran flags
IFLAGS = -Warray-bounds -Warray-temporaries -Wrealloc-lhs-all -fbounds-check \
	-Wimplicit-procedure -Wimplicit-interface \
	-fall-intrinsics -ffree-line-length-none \
	-fopenmp -O3 -ffast-math -xHost -ipo -mcmodel=medium -shared-intel

FC = ifort
FFLAGS = $(IFLAGS)

# -- optimised gfortran compilations
GFLAGS = -Warray-bounds -Warray-temporaries -Wrealloc-lhs-all -fbounds-check \
	-Wimplicit-procedure -Wimplicit-interface \
	-fall-intrinsics -ffree-line-length-none
	-fopenmp -O3 -ffast-math -mtune=native -fwhole-program -combine \
	-mcmodel=medium

FC = gfortran
FFLAGS = $(GFLAGS)

# -- object files
OBJ = rsg.o \
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
