# Add any `module load` or `export` commands that your code needs to
# compile and run to this file.


module load languages/intel/2020-u4

module load languages/anaconda2/5.0.1

export OMP_PLACES=cores
export OMP_PROC_BIND=true
export OMP_NUM_THREADS=28


alias job='sbatch job_submit_d2q9-bgk'

alias o='less d2q9-bgk.out'

alias q='squeue -u $USER'

alias qcancel='scancel -u $USER'

alias n='nano job_submit_d2q9-bgk'

alias m='nano Makefile'

alias c='clear'

alias plot='gnuplot final_state.plt'