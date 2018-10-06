#! /bin/bash 
#$ -S /bin/bash 
#$ -V 
#$ -N Job1 
#$ -pe impi3 36 
#$ -j y
#$ -cwd 

export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=/opt/openmpi-1.8.4/lib:$LD_LIBRARY_PATH
export PATH=/opt/openmpi-1.8.4/bin:$PATH
source  /opt/intel/composer_xe_2013_sp1.0.080/bin/compilervars.sh intel64
# /opt/intel/openmpi-1.8.4/bin/mpirun -np ${NSLOTS}  /home/xgzhao/upload/vasp5.4/vasp.5.4.4/vasp5.4_openmpi/vasp_std   > vasp.log 2 >& 1
