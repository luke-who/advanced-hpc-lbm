# Add any `module load` or `export` commands that your code needs to
# compile and run to this file.

# /* set OMP env variables */
export OMP_NUM_THREADS=14
# set OMP_NUM_THREADS #This command displays the current setting of the OMP_NUM_THREADS environment variable
# export OMP_SCHEDULE=static #pragma omp parallel for uses static schedule by default, which evenly divides the total loop iterations between all threads
# export OMP_PROC_BIND=true
export OMP_PROC_BIND=close
# export OMP_PROC_BIND=spread
export OMP_PLACES=cores

# /* set MPI env variables */
# export I_MPI_DEBUG=4

# /* load necessary modules */
# module load languages/gcc/9.3.0
module load languages/intel/2020-u4 #The Intel compiler is available in the `languages/intel/2020-u4` module. This module puts advixe-cl in your path, so there will be no need to source anything
# module load intel/2017.01
# module load vtune/2017.1.132-GCC-5.4.0-2.26
# module load tools/valgrind/3.15.0

