#!/bin/csh

prctl --unaligned=silent mpirun -np 16 ~/Siesta2_Altix/Exec/siesta_16cpu_201 < art2siesta >& siesta2art

