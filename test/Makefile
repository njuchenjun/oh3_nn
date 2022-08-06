
 FC = ifort
 FFLAGS = -O3 -132 -qopenmp
 LDFLAGS = -mkl -static-intel

#MKLROOT = /work1/soft/intel/compilers_and_libraries_2017.1.132/linux/mkl
#LDFLAGS = -L$(MKLROOT)/lib/intel64 -lmkl_rt -Wl,-rpath $(MKLROOT)/lib/intel64
#LDFLAGS =  $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
#           $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a \
#           -Wl,--start-group \
#           $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
#           $(MKLROOT)/lib/intel64/libmkl_core.a \
#           $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
#           -Wl,--end-group -lpthread -lm

#### end of macro definition ####

OBJS = pes.o test2.o

run:	a.exe
	./a.exe

a.exe:	$(OBJS)
	$(FC) $(FFLAGS) $+ -o $@ $(LDFLAGS)

%.o:	%.F90
	$(FC) $(FFLAGS) -c $< -o $@
%.o:	%.f90
	$(FC) $(FFLAGS) -c $< -o $@
%.o:	%.F
	$(FC) $(FFLAGS) -c $< -o $@
%.o:	%.f
	$(FC) $(FFLAGS) -c $< -o $@
clean:
	rm -f *.o *.exe *.mod fort.*

pes.o:	../pes/energy_oh3.F
	$(FC) $(FFLAGS) -save -c $< -o $@

