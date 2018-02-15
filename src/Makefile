CC = g++
CFLAGS = -Wall -march=native -Os -g -Wno-write-strings

# This part is if MKL is available. It is possible to switch using the
# -DMKL for compiling and adding the necessary linking arguments. 
# Below is a suggestion.
#MKLROOT = /opt/intel/composer_xe_2011_sp1.11.339

#Dynamic linking
#MKL = -I$(MKLROOT)/mkl/include/ -lm -L$(MKLROOT)/mkl/lib/intel64 -L$(MKLROOT)/compiler/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

MKL = -lm -L/mnt/data/Programs/anaconda/pkgs/mkl-rt-11.1-p0/lib/ -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

#Static linking
#MKL =  -Wl,--start-group  $(MKLROOT)/mkl/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/mkl/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -L$(MKLROOT)/compiler/lib/intel64 -liomp5 -lpthread -lm

LINK = -DMKL $(MKL)

#LINK = -lblas -llapack -lfftw3 -lfftw3f

all: timing benchmark

timing_dr_double: timing_dr_slp.cpp read_problem.cpp tools.cpp dr_admm_slp.cpp
	$(CC) $(CFLAGS) -DDOUBLE -o timing_dr_double timing_dr_slp.cpp read_problem.cpp tools.cpp dr_admm_slp.cpp $(LINK)

timing_dr_single: timing_dr_slp.cpp read_problem.cpp tools.cpp dr_admm_slp.cpp
	$(CC) $(CFLAGS) -DSINGLE -o timing_dr_single timing_dr_slp.cpp read_problem.cpp tools.cpp dr_admm_slp.cpp $(LINK)

timing_admm_double: timing_admm_slp.cpp read_problem.cpp tools.cpp dr_admm_slp.cpp
	$(CC) $(CFLAGS) -DDOUBLE -o timing_admm_double timing_admm_slp.cpp read_problem.cpp tools.cpp dr_admm_slp.cpp $(LINK)

timing_admm_kmax_double: timing_admm_slp_kmax.cpp read_problem.cpp tools.cpp
	$(CC) $(CFLAGS) -DDOUBLE -o timing_admm_kmax_double timing_admm_slp_kmax.cpp read_problem.cpp tools.cpp dr_admm_slp.cpp $(LINK)

timing_admm_kmax_single: timing_admm_slp_kmax.cpp read_problem.cpp tools.cpp
	$(CC) $(CFLAGS) -DSINGLE -o timing_admm_kmax_single timing_admm_slp_kmax.cpp read_problem.cpp tools.cpp dr_admm_slp.cpp $(LINK)

timing_admm_single: timing_admm_slp.cpp read_problem.cpp tools.cpp dr_admm_slp.cpp
	$(CC) $(CFLAGS) -DSINGLE -o timing_admm_single timing_admm_slp.cpp read_problem.cpp tools.cpp dr_admm_slp.cpp $(LINK)


benchmark_signal_dr_double: benchmark_signal_dr_slp.cpp read_problem.cpp write_solution.cpp tools.cpp  dr_admm_slp.cpp
	$(CC) $(CFLAGS)  -DDOUBLE -o benchmark_signal_dr_double benchmark_signal_dr_slp.cpp read_problem.cpp write_solution.cpp tools.cpp  dr_admm_slp.cpp $(LINK)

benchmark_signal_admm_double: benchmark_signal_admm_slp.cpp read_problem.cpp write_solution.cpp tools.cpp dr_admm_slp.cpp
	$(CC) $(CFLAGS)  -DDOUBLE -o benchmark_signal_admm_double benchmark_signal_admm_slp.cpp read_problem.cpp write_solution.cpp tools.cpp dr_admm_slp.cpp $(LINK)

benchmark: benchmark_signal_dr_double benchmark_signal_admm_double

timing: timing_dr_double timing_dr_single timing_admm_double timing_admm_single timing_admm_kmax_double timing_admm_kmax_single

clean:
	rm -f *.o
	find . -perm +100 -type f -delete
	rm -f *.*~
