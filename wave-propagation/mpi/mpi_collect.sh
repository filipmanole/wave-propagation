#!/bin/bash
collect mpirun -np 1 ~/wave-propagation/mpi/acoustics ~/wave-propagation/mpi/input
collect mpirun -np 2 ~/wave-propagation/mpi/acoustics ~/wave-propagation/mpi/input
collect mpirun -np 4 ~/wave-propagation/mpi/acoustics ~/wave-propagation/mpi/input
collect mpirun -np 6 ~/wave-propagation/mpi/acoustics ~/wave-propagation/mpi/input
collect mpirun -np 8 ~/wave-propagation/mpi/acoustics ~/wave-propagation/mpi/input
