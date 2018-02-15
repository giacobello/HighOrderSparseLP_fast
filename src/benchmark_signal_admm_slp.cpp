#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#include "ftype.hpp"
#include "dr_admm_slp.hpp"
#include "read_problem.hpp"
#include "write_solution.hpp"



int main(int argc, char **argv)
{ 
    isignal sig;
    solution sol;

	if( argc != 5 ){
		printf("Number of input should be 3: problem_name solution_name no_of_frames gs \n");
		return 0;
	}

	char s [100];
	sprintf(s, "../data/%s", argv[1]);
	char v [100];
	sprintf(v, "../data/%s", argv[2]);
	
	
	char f [100];
	int no_of_frames;
	sscanf(argv[3], "%d", &no_of_frames);

    int gs;
    sscanf(argv[4], "%d", &gs);

	int rpt = 100;

    FTYPE * y;
    FTYPE * u;
    FTYPE * y0;
    FTYPE * u0;

	int i=1;
    int n;
	int r;
    int j;
	double l1, l2;
	struct timeval tim;


	/* init variables based on the first problem */
	sprintf(f,"%s_%d", s, i);
	if( read_problem(f, &sig) ) {
      
      n = sig.m + sig.n;
      y = vector(n, "y");
      u = vector(n, "u");
      y0 = vector(n, "u");
      u0 = vector(n, "y");

      init_solution(sig.n, &sol);
      /* run over all problems in the signal */
      for( i = 1 ; i <= no_of_frames ; i++ ){

        if( i == 1){  //start at zero initialization
          for( j = 0 ; j < n ; j++ ){
            y0[j] = 0;
            u0[j] = 0;
          }
        }
        else{ //warm-start
          for( j = 0 ; j < n ; j++ ){
            y0[j] = y[j];
            u0[j] = u[j];
          }
        }
        //printf("Run problem instance %d\n",i);
        
        /* get signal */
        sprintf(f,"%s_%d",s,i);
        if(read_problem(f, &sig)){
          
          slp_admm *s = new slp_admm(sig.m, sig.n, gs, sig.Y, sig.gamma);
          s->set_epsilon(1e-6);
          s->set_kmax(100);

          gettimeofday(&tim, NULL);
          l1 = tim.tv_sec + (tim.tv_usec/1000000.0);			
          for( r = 0 ; r < rpt; ++r){

            for( j = 0 ; j < n ; j++ ){
              y[j] = y0[j];
              u[j] = u0[j];
            }

            s->set_signal(sig.Y);
            admm(s, y, u);
          }
          
          gettimeofday(&tim, NULL);
          l2 = tim.tv_sec + (tim.tv_usec/1000000.0);
          
          sprintf(f, "%s_%d",v,i);
          *sol.kp = s->k;
          *sol.status = 1;
          for( j = 0 ; j < sig.n ; ++j)
            sol.a[j] = (1/sig.gamma)*y[j];

          write_solution(f, sig.n, sol, (l2-l1)/rpt);
          //printf("Average time : %2.5f [s]\n", (l2-l1)/rpt);
          free_problem(&sig);
          delete s;
        }
        else
          free_solution(&sol);
          return 0;
      }
      free(y);
      free(u);
	}
	else 
      return 0;
}
