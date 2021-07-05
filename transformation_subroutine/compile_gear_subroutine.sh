#!/bin/bash
export OLD_LD=${LD_LIBRARY_PATH}
g++ -c -O3 -std=c++17 -fPIC -I "${EIGENPATH}" -o initial_state_gear.o initial_state_gear.cpp
g++ -c -O3 -std=c++17 -fPIC -I "${EIGENPATH}" -o transformation_umat.o transformation_umat.cpp
g++ -c -O3 -std=c++17 -fPIC -I "${EIGENPATH}" -o utilities.o utilities.cpp
gfortran -c -fPIC getpartinfo.f
ld -r initial_state_gear.o transformation_umat.o getpartinfo.o utilities.o -o gear_subroutine.o
export LD_LIBRARY_PATH=${OLD_LD}
