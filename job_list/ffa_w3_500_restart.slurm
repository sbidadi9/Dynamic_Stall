#!/bin/bash
#SBATCH --nodes=24
#SBATCH --time=24:00:00
#SBATCH --account=storm
#SBATCH --job-name=ffa_w3_500_restart
#SBATCH --qos=high
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shreyas.bidadi@nrel.gov

module purge
module restore
source /nopt/nrel/apps/cpu_stack/env_cpe23.sh
module load PrgEnv-intel
module load cray-python 
module list

export EXAWIND_MANAGER=${HOME}/exawind-manager
cd ${EXAWIND_MANAGER} && source shortcut.sh && cd -
quick-activate /home/sbidadi/exawind-manager/environments/uaero
spack load nalu-wind
which naluX

ranks_per_node=104
mpi_ranks=$(expr $SLURM_JOB_NUM_NODES \* $ranks_per_node)

export LD_PRELOAD=/nopt/nrel/apps/cray-mpich-stall/libs_mpich_nrel_intel/libmpi_intel.so.12
ldd `which naluX` | grep mpi

export MPICH_VERSION_DISPLAY=1
export MPICH_OFI_CQ_STALL=1
export MPICH_OFI_CQ_STALL_USECS=3
export MPICH_OFI_CXI_COUNTER_REPORT=2
export MPICH_ENV_DISPLAY=1

export FI_MR_CACHE_MONITOR=memhooks
export MPICH_SMP_SINGLE_COPY_MODE=NONE
export FI_CXI_RX_MATCH_MODE=software

nalu_exec=naluX

echo "Job name       = $SLURM_JOB_NAME"
echo "Num. nodes     = $SLURM_JOB_NUM_NODES"
echo "Num. MPI Ranks = $mpi_ranks"
echo "Num. threads   = $OMP_NUM_THREADS"
echo "Working dir    = $PWD"

grids=/scratch/sbidadi/Dynamic_Stall/final/nalu_runs/ffa_w3_500/
echo ${grids}
list_of_cases=($(ls ${grids}))

for idx in {0..5};

   do

      echo " "
      echo "Job: "$idx &
      echo "Case: "${list_of_cases[$idx]} &
      echo "Directory contents of " ${list_of_cases[$idx]} &
      ls -alh ${grids}/${list_of_cases[$idx]} &

      date &
      
      echo "Yaml file: " $grids${list_of_cases[$idx]}/*_restart.yaml

      srun -u -N4 -n384 --ntasks-per-node=96 --distribution=block:cyclic --cpu_bind=map_cpu:0,52,13,65,26,78,39,91,1,53,14,66,27,79,40,92,2,54,15,67,28,80,41,93,3,55,16,68,29,81,42,94,4,56,17,69,30,82,43,95,5,57,18,70,31,83,44,96,6,58,19,71,32,84,45,97,7,59,20,72,33,85,46,98,8,60,21,73,34,86,47,99,9,61,22,74,35,87,48,100,10,62,23,75,36,88,49,101,11,63,24,76,37,89,50,102,12,64,25,77,38,90,51,103 ${nalu_exec} -i $grids${list_of_cases[$idx]}/*_restart.yaml -o $grids${list_of_cases[$idx]}/log$idx.out &

      sleep 1

   done

wait
