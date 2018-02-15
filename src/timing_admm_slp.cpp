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

    FTYPE * y;
    FTYPE * u;

	/* get signal */
	if(read_problem(s, &sig) ){
	      
      n = sig.m + sig.n;

      y = vector(n, "y");
      u = vector(n, "u");

      slp_admm *s = new slp_admm(sig.m, sig.n, 1, sig.Y, sig.gamma);
      s->set_epsilon(1e-6);

      /* warm up */
      for( i = 0 ; i < 100; ++i ){
        for( j = 0 ; j < n ; j++ ){
          y[j] = 0;
          u[j] = 0;
        }

        s->set_signal(sig.Y);
        admm(s, y, u);
      }
      
      /* run timing */
      gettimeofday(&tim, NULL);
      l1 = tim.tv_sec + (tim.tv_usec/1000000.0);
      for( i = 0 ; i < rep; ++i ){
        gettimeofday(&tim, NULL);
        t1 = tim.tv_sec + (tim.tv_usec/1000000.0);
      
        for( j = 0 ; j < n ; j++ ){
          y[j] = 0;
          u[j] = 0;
        }
        
        s->set_signal(sig.Y);

        admm(s, y, u);
      
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
      free(y);
      free(u);
      
      delete(s);

      return 0;
	}
	else
      return 1;
}
