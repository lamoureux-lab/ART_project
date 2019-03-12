OBJ = art.o defs.o  random.o initialize.o min_converge.o \
	utils.o lanczos.o find_saddle.o diis.o initialize_potential.o \
	saddle_converge.o read_parameters.o restart.o \
	calcforce.o \

#SYS = MacIntel
#include ${SYS}.arch

FORTRAN = ifort -g 
FLAGS =  -O3 -cpp
LLIBS = -mkl=sequential
TRACE = -traceback

all : art

art : ${OBJ} 
	$(FORTRAN) $(TRACE)  $(OBJ) $(FFLAGS)  -o  art $(LLIBS)

defs.o: defs.f90
	$(FORTRAN) $(TRACE)  -c $(FLAGS) defs.f90

initialize_potential.o: defs.o initialize_potential.f90
	$(FORTRAN) $(TRACE) -c $(FLAGS) initialize_potential.f90

read_parameters.o: defs.o read_parameters.f90
	$(FORTRAN) $(TRACE)  -c $(FLAGS) read_parameters.f90

random.o: random.f90
	$(FORTRAN) $(TRACE)  -c $(FLAGS) random.f90

art.o: art.f90 lanczos.o defs.o utils.o
	$(FORTRAN) $(TRACE) -c $(FLAGS) art.f90

find_saddle.o: defs.o find_saddle.f90
	$(FORTRAN) $(TRACE) -c $(FLAGS) find_saddle.f90

initialize.o: defs.o find_saddle.o utils.o initialize.f90
	$(FORTRAN) $(TRACE) -c $(FLAGS) initialize.f90

lanczos.o: defs.o lanczos.f90 random.o
	$(FORTRAN) $(TRACE) -c $(FLAGS) lanczos.f90

restart.o: defs.o restart.f90
	$(FORTRAN) $(TRACE) -c $(FLAGS) restart.f90

min_converge.o: defs.o min_converge.f90
	$(FORTRAN) $(TRACE) -c $(FLAGS) min_converge.f90

saddle_converge.o: defs.o saddle_converge.f90
	$(FORTRAN) $(TRACE) -c $(FLAGS) saddle_converge.f90

calcforce.o : defs.o calcforce.f90
	$(FORTRAN) $(TRACE) -c $(FLAGS) calcforce.f90

utils.o: defs.o utils.f90
	$(FORTRAN) $(TRACE) -c $(FLAGS) utils.f90

diis.o:defs.o diis.f90
	$(FORTRAN) $(TRACE) -c $(FLAGS) diis.f90

clean:
	rm -f *.o *.mod art
