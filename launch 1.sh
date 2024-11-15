#!/bin/bash
#SBATCH -D /users/addh496/sharedscratch/CaNS/run # Working Directory

#SBATCH -J cans_cha					  # Name of Job
#SBATCH --nodes=2			 	 	  # Number of Nodes
#SBATCH --ntasks=32				 	  # Total Number of Processors
##SBATCH --ntasks-per-node=16			 	  # Number tasks per node
#! The hyperion nodes have 24*2 CPUs (cores) each.
#SBATCH --time=01:30:00				 	  # Total cpu time
#SBATCH --mem=8GB                                  	  # Expected memory usage (0 means use all available memory)
#SBATCH --exclusive
##SBATCH --mail-user=giorgio.cavallazzi@city.ac.uk	  # Email to send the notifications
#SBATCH --mail-type=ALL
#SBATCH --output=R-%x.%j.out
#!SBATCH --qos=standard
#!SBATCH --partition=test

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
#mpi_tasks_per_node=$SLURM_TASKS_PER_NODE

#! Modify the settings below to specify the application's environment, location 
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
flight env activate gridware

#! Insert additional module load commands after this line if needed:
module add mpi/openmpi/4.1.1
module add fftw/3.3.10

#! Full path to application executable: 
application="/users/addh496/sharedscratch/CaNS/run/cans"

#! Run options for the application:
options=""

#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"  # The value of SLURM_SUBMIT_DIR sets workdir to the directory
                             # in which sbatch is run.

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 32:
export OMP_NUM_THREADS=1

#! Number of MPI tasks to be started by the application per node and in total (do not change):
np=$SLURM_NTASKS

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
#export I_MPI_PIN_ORDER=compact # Adjacent domains have minimal sharing of caches/sockets

#! Uncomment one choice for CMD below (add mpirun/mpiexec options if necessary):

#! Choose this for a MPI code (possibly using OpenMP) using OpenMPI:
CMD="mpirun -np $np $application $options"

###############################################################
### You should not have to change anything below this line ####
###############################################################

cd $workdir
echo -e "Changed directory to `pwd`.\n"

JOBID=$SLURM_JOB_ID

eval $CMD 
