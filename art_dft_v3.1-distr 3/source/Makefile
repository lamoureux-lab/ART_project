OBJ = art.o defs.o  random.o initialize.o min_converge.o \
	utils.o lanczos.o find_saddle.o diis.o initialize_potential.o \
	saddle_converge.o read_parameters.o restart.o \
	calcforce.o neighbour.o siesta_force.o siesta_min.o calcfo_sw.o \
	gaussian_force.o gaussian_min.o

#SYS = MacIntel
#include ${SYS}.arch

FORTRAN = ifort -g -O2
FLAGS =  -O3 -cpp
LLIBS = -mkl=sequential

all : artdft

artdft : ${OBJ} ${SIESTA_pot}
	$(FORTRAN)  $(OBJ) $(SIESTA_pot) $(FFLAGS)  -o  artdft $(LLIBS)

silicon :$(OBJ) $(SW_pot)
	$(FORTRAN)   $(OBJ) $(SW_pot) $(FFLAGS)  -o  silicon $(LLIBS)

defs.o: defs.f90
	$(FORTRAN)  -c $(FLAGS) defs.f90

initialize_potential.o: defs.o initialize_potential.f90
	$(FORTRAN) -c $(FLAGS) initialize_potential.f90

read_parameters.o: defs.o read_parameters.f90
	$(FORTRAN)  -c $(FLAGS) read_parameters.f90

neib.o : neib.f90
	$(FORTRAN)  -c $(FLAGS) neib.f90

random.o: random.f90
	$(FORTRAN)  -c $(FLAGS) random.f90

art.o: art.f90 lanczos.o defs.o utils.o
	$(FORTRAN) -c $(FLAGS) art.f90

calcfo_sw.o: defs.o calcfo_sw.f90
	$(FORTRAN) -c $(FLAGS) calcfo_sw.f90

siesta_force.o: defs.o siesta_force.f90
	$(FORTRAN) -c $(FLAGS) siesta_force.f90

gaussian_force.o: defs.o gaussian_force.f90
	$(FORTRAN) -c $(FLAGS) gaussian_force.f90

find_saddle.o: defs.o find_saddle.f90
	$(FORTRAN) -c $(FLAGS) find_saddle.f90

initialize.o: defs.o find_saddle.o utils.o initialize.f90
	$(FORTRAN) -c $(FLAGS) initialize.f90

lanczos.o: defs.o lanczos.f90 random.o
	$(FORTRAN) -c $(FLAGS) lanczos.f90

restart.o: defs.o restart.f90
	$(FORTRAN) -c $(FLAGS) restart.f90

siesta_min.o: defs.o siesta_min.f90
	$(FORTRAN) -c $(FLAGS) siesta_min.f90

gaussian_min.o: defs.o gaussian_min.f90
	$(FORTRAN) -c $(FLAGS) gaussian_min.f90

min_converge.o: defs.o min_converge.f90
	$(FORTRAN) -c $(FLAGS) min_converge.f90

neighbour.o: defs.o neighbour.f90 calcfo_sw.o
	$(FORTRAN) -c $(FLAGS) neighbour.f90

saddle_converge.o: defs.o saddle_converge.f90
	$(FORTRAN) -c $(FLAGS) saddle_converge.f90

calcforce.o : defs.o calcforce.f90
	$(FORTRAN) -c $(FLAGS) calcforce.f90

utils.o: defs.o utils.f90
	$(FORTRAN) -c $(FLAGS) utils.f90

diis.o:defs.o diis.f90
	$(FORTRAN) -c $(FLAGS) diis.f90

clean:
	rm -f *.o *.mod siestart silicon artdft
