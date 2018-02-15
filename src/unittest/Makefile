
GTEST = /home/tlj/Programs/googletest/googletest

# This part is if MKL is available. It is possible to switch using the
# -DMKL for compiling and adding the necessary linking arguments. 
# Below is a suggestion.
#MKLROOT = /opt/intel/composer_xe_2011_sp1.11.339

#Dynamic linking
#MKL = -I$(MKLROOT)/mkl/include/ -lm -L$(MKLROOT)/mkl/lib/intel64 -L$(MKLROOT)/compiler/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

#Static linking
#MKL =  -Wl,--start-group  $(MKLROOT)/mkl/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/mkl/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -L$(MKLROOT)/compiler/lib/intel64 -liomp5 -lpthread -lm

MKL = -lm -L/mnt/data/Programs/anaconda/pkgs/mkl-rt-11.1-p0/lib/ -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

LINK = -DMKL $(MKL)


#Include both single and double precision fftw libraries for convenience
#LINK = -lblas -llapack -lfftw3 -lfftw3f

all:
	g++ -Wall -g -Wno-write-strings -DDOUBLE -isystem $(GTEST)/include -pthread  ../dr_admm_slp.cpp ../tools.cpp test_dr_admm_slp.cpp libgtest.a -o test_dr_admm_slp_double $(LINK) 

	g++ -Wall -g -Wno-write-strings -DSINGLE -isystem $(GTEST)/include -pthread   ../dr_admm_slp.cpp ../tools.cpp test_dr_admm_slp.cpp libgtest.a -o test_dr_admm_slp_single $(LINK) 
