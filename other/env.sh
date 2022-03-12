module load languages/intel/2018-u3

module load languages/anaconda2/5.0.1


alias job='sbatch job_submit_d2q9-bgk'

alias o='less d2q9-bgk.out'

alias q='squeue -u $USER'

alias qcancel='scancel -u $USER'

alias n='nano job_submit_d2q9-bgk'

alias m='nano Makefile'

alias c='clear'

alias plot='gnuplot final_state.plt'