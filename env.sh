# Add any `module load` or `export` commands that your code needs to
# compile and run to this file.

################################
### Set OpenMP env variables ###
################################
export OMP_NUM_THREADS=14
# set OMP_NUM_THREADS # This command displays the current setting of the OMP_NUM_THREADS environment variable
# export OMP_SCHEDULE=static # Static schedule divides the iterations into chunks and assigns chunks to threads in a round-robin manner. 'pragma omp parallel for' uses static schedule by default, which evenly divides the total loop iterations between all threads.

#----- Thread Pinning -----#
# export OMP_PROC_BIND=true # Threads won’t move, and follow proc_bind clauses or else the implementation default pinning
export OMP_PROC_BIND=close # Threads are assigned to places close to the master thread. If OMP_NUM_THREADS==ncores: thread 0 will pin to core 0; thread 1 will pin to core 1; etc
# export OMP_PROC_BIND=spread # Threads are assigned to places “sparsely”. If OMP_NUM_THREADS==ncores: thread 0 will pin to socket 0 core0; thread 1 will pin to socket 1 core 0; thread 2 will pin to socket 0 core 1; etc.

#----- Places -----#
export OMP_PLACES=cores # Each place is a single core(containing one or more hardware threads).


#############################
### Set MPI env variables ###
#############################
# export I_MPI_DEBUG=4


################################
#### Load necessary modules ####
################################
# module load languages/gcc/9.3.0
module load languages/intel/2020-u4 #The Intel compiler is available in the `languages/intel/2020-u4` module. This module puts advixe-cl in your path, so there will be no need to source anything
# module load intel/2017.01
# module load vtune/2017.1.132-GCC-5.4.0-2.26
# module load tools/valgrind/3.15.0

