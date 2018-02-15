#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#include "ftype.hpp"
#include "dr_admm_slp.hpp"
#include "read_problem.hpp"


int main(int argc, char **argv)
{ 
	isignal sig;

	char *s = "../data/problem570250_60";
	int i, j;
	double l1, l2;
	double t1, t2, t=0.0;
	struct timeval tim;
	int rep = 100;
    int n;

    FTYPE * z0;

	/* get signal */
	if(read_problem(s, &sig) ){
	      
      n = sig.m + sig.n;

      slp_dr *s = new slp_dr(sig.m, sig.n, 1, sig.Y, sig.gamma);
      s->set_epsilon(1e-6);
      s->set_verbose(0);

      z0 = vector(n, "z0");

      /* warm up */
      for( i = 0 ; i < 100; ++i){
        for( j = 0 ; j < n ; j++ )
          z0[j] = 0;
        
        s->set_signal(sig.Y);
        dr(s, z0);
      }

      /* run timing */
      gettimeofday(&tim, NULL);
      l1 = tim.tv_sec + (tim.tv_usec/1000000.0);
      for( i = 0 ; i < rep; ++i){
        gettimeofday(&tim, NULL);
        t1 = tim.tv_sec + (tim.tv_usec/1000000.0);

        for( j = 0 ; j < n ; j++ )
          z0[j] = 0;
        
        s->set_signal(sig.Y);

        dr(s, z0);

        gettimeofday(&tim, NULL);
        t2 = tim.tv_sec + (tim.tv_usec/1000000.0);
        t += t2 - t1;
      }
      gettimeofday(&tim, NULL);
      l2 = tim.tv_sec + (tim.tv_usec/1000000.0);
      printf("Number of iterations : %2d \n", s->k);
      printf("Average time : %2.5f [s]\n", (l2-l1)/rep);
      printf("Average time via individual times: %2.5f [s]\n", t/rep);
      
      free_problem(&sig);
      free(z0);
      delete(s);

      return 0;
	}
	else
      return 1;
}
