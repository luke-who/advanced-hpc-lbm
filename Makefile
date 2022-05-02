
   
# Makefile

EXE=d2q9-bgk

#*******Set compiler*******#
# CC=gcc
# CC=icc
CC=mpiicc

REPORT = -qopt-report=5 -qopt-report-phase=vec
# OPTMLEVEL = -O0
OPTMLEVEL = -Ofast -xAVX2
# OPTMLEVEL = -Ofast -xAVX #-fast #note -fast would make Intel Advisor stop working

# TARGET_PLATFORM = -mtune=native
TARGET_PLATFORM = -xHOST #-xbroadwell
Profile_Generate = -pg
Unroll_loops = -funroll-all-loops
#*******Set compiler flags*******#
# CFLAGS= -std=c11 -Wall $(OPTMLEVEL) $(TARGET_PLATFORM) -g #-qopenmp #NOALIAS
# CFLAGS= -std=c11 -Wall $(OPTMLEVEL) $(TARGET_PLATFORM) -pg -qopenmp $(REPORT) #$(Profile_Generate) 
CFLAGS= -std=c11 -Wall $(OPTMLEVEL) $(TARGET_PLATFORM) -pg $(REPORT) #$(Profile_Generate) 

Intel_advisor = -Wl,-u__poll -Wl,-udlclose -Wl,-udlopen

LIBS = -lm

#*******Compute final state & av_velocity and output to files*******#
FINAL_STATE_FILE=./final_state.dat
AV_VELS_FILE=./av_vels.dat

DIM=128
# DIM=256
# DIM=1024
#*******Set the referenced final state & av_velocity files*******#
REF_FINAL_STATE_FILE=check/$(DIM)x$(DIM).final_state.dat
REF_AV_VELS_FILE=check/$(DIM)x$(DIM).av_vels.dat



all: $(EXE)

$(EXE): $(EXE).c
	$(CC) $(CFLAGS) $^ $(LIBS) -o $@

check:
	python check/check.py --ref-av-vels-file=$(REF_AV_VELS_FILE) --ref-final-state-file=$(REF_FINAL_STATE_FILE) --av-vels-file=$(AV_VELS_FILE) --final-state-file=$(FINAL_STATE_FILE)

.PHONY: all check clean

run:
	sbatch job_submit_d2q9-bgk
	watch -n 1 squeue -u az16408
	# while true; do squeue -u az16408;date ; sleep 1; done

cancel:
	scancel -u az16408
cat:
	cat d2q9-bgk.out
clean:
	rm -f $(EXE) $(EXE).out $(EXE).optrpt