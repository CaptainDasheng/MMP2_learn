#!/bin/bash 

mpi=$1    # parrallel compiler % 
path=$2   # root of path of parrallel compiler % 

# if intel mpi % 
if [ "$mpi" == 'intel' ]; then
   echo $mpi, $path
fi

# if opempi % 
if [ "$mpi" == 'openp' ]; then
   echo $mpi, $path
fi

# if mpich % 
if [ "$mpi" == 'mpich' ]; then
   echo $mpi, $path
fi

# if mavpich % 
if [ "$mpi" == 'mpich' ]; then
   echo $mpi, $path
fi
