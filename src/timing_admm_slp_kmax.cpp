#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <algorithm>

#include "ftype.hpp"
#include "dr_admm_slp.hpp"
#include "read_problem.hpp"


double min(double *x, int n){

  double m = x[0];
  for(int i = 1 ; i < n ; ++i){
    if(x[i] < m)
      m = x[i];
  }
  return m;
}

  /* run timing */
double timing(slp_admm * s, FTYPE *x, FTYPE *y, FTYPE *u, int n){

  double l1, l2;
  struct timeval tim;
  int i, j;
  int rpt = 1;
  double tau = -1;

  while( tau < 0.2){
    rpt *= 10; 

    gettimeofday(&tim, NULL);
    l1 = tim.tv_sec + (tim.tv_usec/1000000.0);

    for( i = 0 ; i < rpt; ++i ){
      
      for( j = 0 ; j < n ; j++ ){
	y[j] = 0;
	u[j] = 0;
      }

      s->set_signal(x);      
      admm(s, y, u);
    }
    
    gettimeofday(&tim, NULL);
    l2 = tim.tv_sec + (tim.tv_usec/1000000.0);

    tau = l2-l1;
  }

  return tau/rpt;
}

int main(int argc, char **argv)
{ 
  isignal sig;
  
  char *s = "../data/problem570250_60";
  int i, ii;
  int n;
  int T = 4;
 
  FTYPE * y;
  FTYPE * u;
  double *t;

  int K = 6;
  double k[6] = {5, 10, 15, 20, 25, 30};


  t = (double*)malloc(T*sizeof(double));
  
  /* get signal */
  if(read_problem(s, &sig) ){
    
    n = sig.m + sig.n;
    
    y = vector(n, "y");
    u = vector(n, "u");
    
    slp_admm *s = new slp_admm(sig.m, sig.n, 1, sig.Y, sig.gamma);
    s->set_epsilon(-1);
    
    for(ii = 0 ; ii < K ; ++ii){
      s->set_kmax(k[ii]);
      
      for(i = 0 ; i < T ; ++i){
	t[i] = timing(s, sig.Y, y, u, n);
      }
      
      printf("Number of iterations : %2d \n", s->k);
      printf("Average time : %2.5f [ms]\n", min(t, T)*1e3);
    }

    free_problem(&sig);
    free(y);
    free(u);
    free(t);

    delete s;
    
    return 0;
  }
  else
    return 1;
}
