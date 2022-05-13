# HPC Coursework
[![Intel C/C++ compiler (icc)](https://img.shields.io/badge/icc-v19.1.3.304-blue?&logo=Intel)](https://www.openmp.org/specifications/)
[![OpenMP](https://img.shields.io/badge/OpemMP-v4.0-blue?&logo=Intel)](https://www.openmp.org/specifications/)
[![Intel MPI](https://img.shields.io/badge/MPI-2019%20Update%209-informational?&logo=intel)](https://www.intel.com/content/www/us/en/developer/articles/release-notes/mpi-library-release-notes-linux.html)
-------------------------------------------------------------------------------------------------------------------------------------
Coursework for the Advanced High Performance Computing class.

* Source code is in the [d2q9-bgk.c](d2q9-bgk.c) file
* Results checking scripts are in the [check/](check/) directory

## Calculating % of achieved memory bandwidth
The formula for calculating Bandwidth of the 256x256 input is as follows:

```math
Bandwidth = \frac{size\_of\_one\_grid * (timestep(params.maxIters))}{Elapsed\ Compute\ Time}
```

```math
    = \frac{256*256 * 9 (speeds) * 4 (bytes\_per\_float) * 2(cells+tmp\_cells) * 2(read+write) * 80000(timestep)}{35.7s}
```

```math
    = \frac{2.36(MiBytes) * 4 * 80000}{35.7s} = 21.25GB/s
```
 The maximum achievable L2 memory bandwidth is 84.88 GB/s, from here we can calculate the fraction of STREAM bandwidth = $`\frac{21.25}{84.88} = 25\%`$.

## Compiling and running

### Set env variables
To set environment variables for OpenMP and MPI,
    
    source env.sh

see examples in [env.sh](env.sh)

To compile type `make`. Editing the values for `CC` and `CFLAGS` in the [Makefile](Makefile) can be used to enable different compiler options or use a different compiler. These can also be passed on the command line:

    $ make CFLAGS="-O3 -fopenmp -DDEBUG"

Run all three make arguments(all check clean)

    make .PHONY

Input parameter and obstacle files are all specified on the command line of the `d2q9-bgk` executable.

### Serial & OpenMP Usage:

For OpenMP to set the number of cores to use, choose one of the following methods(note SLURM calls a core a "cpu").:

1. Set it in the SLURM job script

        #SBATCH --ntasks-per-node 1
        #SBATCH --cpus-per-task 28
    or 

        #SBATCH --ntasks-per-node 28
2. Set Environment variables `OMP_NUM_THREADS=28` in environment
3. API calls `#include <omp.h>` and use `omp_set_num_threads(num_threads);` function in your code right before the upcoming parallel regions `#pragma omp parallel`

4. Use Clauses `#pragma omp parallel num_threads(28)` to use all 28 cores

***Note any one of these methods above will override another when used together, depending on the order of execution***

Then run the executable `d2q9-bgk` with input parameter and obstacle files:

    $ ./d2q9-bgk <paramfile> <obstaclefile>
eg:

    $ ./d2q9-bgk input_256x256.params obstacles_256x256.dat

### MPI Usage:

For MPI to set the number of cores to use:

Set it in the SLURM job script, change the number of cores with `--nodes` & `--ntasks-per-node` accordingly. On BC4 each node has 28 cores

    #SBATCH --nodes 4
    #SBATCH --ntasks-per-node 28
    #SBATCH --cpus-per-task=1
    
Use `mpirun` to run the executable `d2q9-bgk`:

    $ mpirun -np <num_cores> ./d2q9-bgk <paramfile> <obstaclefile>
Alternatively use `srun` in SLURM:
    
    $ srun --mpi=pmi2 -n <num_cores> ./d2q9-bgk input_128x128.params obstacles_128x128.dat
eg:

`mpirun`:
    
    $ mpirun -np 112 ./d2q9-bgk input_128x128.params obstacles_128x128.dat
`srun`:
    
    $ srun --mpi=pmi2 -n 112 ./d2q9-bgk input_128x128.params obstacles_128x128.dat

see more examples in [job_submit_d2q9-bgk](job_submit_d2q9-bgk)

## Using OpenMP


## Checking results

An automated result checking function is provided that requires you to load a particular Python module (`module load languages/anaconda2/5.0.1`). Running `make check` will check the output file (average velocities and final state) against some reference results. By default, it should look something like this:
```
[ab12345@bc4login3 advanced-hpc-lbm]$ make check
python check/check.py --ref-av-vels-file=check/128x128.av_vels.dat --ref-final-state-file=check/128x128.final_state.dat --av-vels-file=./av_vels.dat --final-state-file=./final_state.dat
Total difference in av_vels : 3.370462974994E-02
Biggest difference (at step 39619) : -1.788416140000E-06
  1.317970920354E-02 vs. 1.317792078740E-02 = -0.014%

Total difference in final_state : 1.069222150899E-02
Biggest difference (at coord (3,5)) : 7.545896299962E-07
  3.334492444992E-02 vs. 3.334567903955E-02 = 0.0023%

Both tests passed!
```

This script takes both the reference results and the results to check (both average velocities and final state). This is also specified in the makefile and can be changed like the other options:

    $ make check REF_AV_VELS_FILE=check/128x256.av_vels.dat REF_FINAL_STATE_FILE=check/128x256.final_state.dat
    python check/check.py --ref-av-vels-file=check/128x256.av_vels.dat --ref-final-state-file=check/128x256.final_state.dat --av-vels-file=./av_vels.dat --final-state-file=./final_state.dat
    ...

All the options for this script can be examined by passing the --help flag to it.

    $ python check/check.py --help
    usage: check.py [-h] [--tolerance TOLERANCE] --ref-av-vels-file
                    REF_AV_VELS_FILE --ref-final-state-file REF_FINAL_STATE_FILE
    ...


## Running on BlueCrystal Phase 4

When you wish to submit a job to the queuing system on BlueCrystal, you should use the job submission script provided.

    $ sbatch job_submit_d2q9-bgk
    
This will dispatch a job to the queue, which you can monitor using the
`squeue` command:

    $ squeue -u $USER

Or use the one below to submit & monitor in one step, watch it run in real time:

    $ make run
 
Cancel Jobs: 
    
    $ make cancel
    
When finished, the output from your job will be in a file called
`d2q9-bgk.out`:

    $ less d2q9-bgk.out
Alternatively: 

    $ make cat
    
If you wish to run a different set of input parameters, you should
modify `job_submit_d2q9-bgk` to update the value assigned to `options`.

## Checking submission content

Before handing in the coursework, you can use the `check_submission.sh` script to make sure that your code builds in a clean environment. This will reduce the chances of the automarker failing to build or run your code.

To use the script, simply run it from the directory containing the files you intend to submit:

    $ /path/to/check_submission.sh

The script will:

1. Unload all the modules currently loaded.
2. Load your modules and environment variables specified in `env.sh`.
3. Use `make` to build your code and verify that an executable with the expected name is produced.

If the submission checking script prints any errors, you should try to address those before you hand in. 

Note that `check_submission.sh` does _not_ run your code, and so you _cannot_ verify that the results produced by your application validate just by running this script. You should check the correctness of your results separately, e.g. using `make check`.


# Serial, OpenMP & MPI output for sample inputs

## Serial
Single core, run times were taken on a Phase 4 node using the base (gcc) compiler and base compiler flags as found in the Makefile:

- 128x128
```
$ ./d2q9-bgk  input_128x128.params obstacles_128x128.dat
==done==
Reynolds number:		9.751927375793E+00
Elapsed time:			38.387577 (s)
Elapsed user CPU time:		38.388736 (s)
Elapsed system CPU time:	0.003000 (s)
```

- 128x256
```
$ ./d2q9-bgk  input_128x256.params obstacles_128x256.dat
==done==
Reynolds number:		3.715003967285E+01
Elapsed time:			77.446019 (s)
Elapsed user CPU time:		77.450619 (s)
Elapsed system CPU time:	0.003000 (s)
```

- 256x256
```
$ ./d2q9-bgk  input_256x256.params obstacles_256x256.dat
==done==
Reynolds number:		1.005141162872E+01
Elapsed time:			309.040200 (s)
Elapsed user CPU time:		309.061111 (s)
Elapsed system CPU time:	0.004000 (s)
```

- 1024x1024
```
$ ./d2q9-bgk  input_1024x1024.params obstacles_1024x1024.dat
==done==
Reynolds number:		3.375851392746E+00
Elapsed time:			1287.501875 (s)
Elapsed user CPU time:		1287.568113 (s)
Elapsed system CPU time:	0.029001 (s)
```
## OpemMP

OpenMP Running with 28 cores, run times were taken on a Phase 4 node using the Intel (icc) compiler and certain compiler flags as found in the Makefile:

- 128x128
```

```

- 128x256
```

```

- 256x256
```

```

- 1024x1024
```

```

## MPI

MPI Running with 112 cores, run times were taken on a Phase 4 node using the Intel (icc) compiler and certain compiler flags as found in the Makefile:

- 128x128
```
[ab12345@bc4login3 advanced-hpc-lbm]$ mpirun -print-rank-map -np 112 ./d2q9-bgk input_128x128.params obstacles_128x128.dat
Running on host compute106.bc4.acrc.priv
Time is Fri May 13 00:37:29 BST 2022
Directory is /user/home/az16408/advanced-hpc-lbm
Slurm job ID is 10338478
This job runs on the following machines:
compute[106-109]
(compute106:0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27)
(compute107:28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55)
(compute108:56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83)
(compute109:84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111)
==done==
Reynolds number:                9.764906883240E+00
Elapsed Init time:                      0.000307 (s)
Elapsed Compute time:                   0.837770 (s)
Elapsed Collate time:                   0.001810 (s)
Elapsed Total time:                     0.839887 (s)
$ make check
python check/check.py --ref-av-vels-file=check/128x128.av_vels.dat --ref-final-state-file=check/128x128.final_state.dat --av-vels-file=./av_vels.dat --final-state-file=./final_state.dat
Total difference in av_vels : 3.370462974994E-02
Biggest difference (at step 39619) : -1.788416140000E-06
  1.317970920354E-02 vs. 1.317792078740E-02 = -0.014%

Total difference in final_state : 1.069222150899E-02
Biggest difference (at coord (3,5)) : 7.545896299962E-07
  3.334492444992E-02 vs. 3.334567903955E-02 = 0.0023%

Both tests passed!
```

- 128x256
```
$ mpirun -np 112 ./d2q9-bgk input_128x256.params obstacles_128x256.dat
Running on host compute106.bc4.acrc.priv
Time is Fri May 13 14:47:30 BST 2022
Directory is /user/home/az16408/advanced-hpc-lbm
Slurm job ID is 10339316
This job runs on the following machines:
compute[106-109]
==done==
Reynolds number:                3.718533325195E+01
Elapsed Init time:                      0.014253 (s)
Elapsed Compute time:                   1.394733 (s)
Elapsed Collate time:                   0.002658 (s)
Elapsed Total time:                     1.411644 (s)
$ make check
python check/check.py --ref-av-vels-file=check/128x256.av_vels.dat --ref-final-state-file=check/128x256.final_state.dat --av-vels-file=./av_vels.dat --final-state-file=./final_state.dat
Total difference in av_vels : 8.179887892935E-03
Biggest difference (at step 68) : -3.011414789994E-08
  5.813845782541E-04 vs. 5.813544641062E-04 = -0.0052%

Total difference in final_state : 4.311647071313E-02
Biggest difference (at coord (76,165)) : 1.389606820001E-06
  3.316890075803E-02 vs. 3.317029036485E-02 = 0.0042%

Both tests passed!
```

- 256x256
```
$ mpirun -np 112 ./d2q9-bgk input_256x256.params obstacles_256x256.dat
Running on host compute106.bc4.acrc.priv
Time is Fri May 13 14:39:29 BST 2022
Directory is /user/home/az16408/advanced-hpc-lbm
Slurm job ID is 10339310
This job runs on the following machines:
compute[106-109]
==done==
Reynolds number:                1.008056449890E+01
Elapsed Init time:                      0.000864 (s)
Elapsed Compute time:                   3.305648 (s)
Elapsed Collate time:                   0.005390 (s)
Elapsed Total time:                     3.311902 (s)
$ make check
python check/check.py --ref-av-vels-file=check/256x256.av_vels.dat --ref-final-state-file=check/256x256.final_state.dat --av-vels-file=./av_vels.dat --final-state-file=./final_state.dat
Total difference in av_vels : 1.740173914993E-01
Biggest difference (at step 79520) : -4.954269670000E-06
  1.360202953219E-02 vs. 1.359707526252E-02 = -0.036%

Total difference in final_state : 8.806575322403E-02
Biggest difference (at coord (138,159)) : 1.448792590002E-06
  3.330772370100E-02 vs. 3.330917249359E-02 = 0.0043%

Both tests passed!
```

- 1024x1024
```
$ mpirun -np 112 ./d2q9-bgk input_1024x1024.params obstacles_1024x1024.dat
Running on host compute106.bc4.acrc.priv
Time is Fri May 13 14:42:47 BST 2022
Directory is /user/home/az16408/advanced-hpc-lbm
Slurm job ID is 10339311
This job runs on the following machines:
compute[106-109]
==done==
Reynolds number:                3.377239465714E+00
Elapsed Init time:                      0.009041 (s)
Elapsed Compute time:                   6.620043 (s)
Elapsed Collate time:                   0.056015 (s)
Elapsed Total time:                     6.685099 (s)
$ make check
python check/check.py --ref-av-vels-file=check/1024x1024.av_vels.dat --ref-final-state-file=check/1024x1024.final_state.dat --av-vels-file=./av_vels.dat --final-state-file=./final_state.dat
Total difference in av_vels : 3.487499986143E-03
Biggest difference (at step 326) : -2.434364265000E-07
  3.626788966358E-04 vs. 3.624354602093E-04 = -0.067%

Total difference in final_state : 9.898915688519E-01
Biggest difference (at coord (1,1021)) : 1.978754159998E-06
  3.326061367989E-02 vs. 3.326259243405E-02 = 0.0059%

Both tests passed!
```

# Visualisation

You can view the final state of the simulation by creating a .png image file using a provided Gnuplot script:

    $ gnuplot final_state.plt
