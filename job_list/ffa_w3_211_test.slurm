#!/bin/bash
#SBATCH --nodes=6
#SBATCH --time=00:20:00
#SBATCH --account=storm
#SBATCH --job-name=ffa_w3_211_test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shreyas.bidadi@nrel.gov

module purge
export SPACK_MANAGER=/scratch/sbidadi/spack-manager
source ${SPACK_MANAGER}/start.sh
spack-start
module restore && module swap PrgEnv-cray PrgEnv-intel
quick-activate /scratch/sbidadi/spack-manager/environments/example-env
spack load exawind

ranks_per_node=104
mpi_ranks=$(expr $SLURM_JOB_NUM_NODES \* $ranks_per_node)
export OMP_NUM_THREADS=1  # Max hardware threads = 4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

nalu_exec=naluX

echo "Job name       = $SLURM_JOB_NAME"
echo "Num. nodes     = $SLURM_JOB_NUM_NODES"
echo "Num. MPI Ranks = $mpi_ranks"
echo "Num. threads   = $OMP_NUM_THREADS"
echo "Working dir    = $PWD"

grids=/scratch/sbidadi/Dynamic_Stall/final/nalu_runs/ffa_w3_211_test/
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
      
      echo "Yaml file: " $grids${list_of_cases[$idx]}/*.yaml
      srun -u -N1 -n52 --ntasks-per-node=52 --distribution=cyclic:cyclic --cpu_bind=cores ${nalu_exec} -i $grids${list_of_cases[$idx]}/*.yaml -o $grids${list_of_cases[$idx]}/log$idx.out &

      sleep 1

   done

wait
