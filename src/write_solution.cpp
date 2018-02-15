#include "write_solution.hpp"


/*writes the solution and timing to the file f*/
int write_solution(char *fn, int n, solution sol, FTYPE time){
	FILE *f;
	int i;

	f = fopen(fn,"w");
	if( f != NULL){
		fprintf( f, "%d\n", *sol.kp);
		fprintf( f, "%d\n", *sol.status);
		fprintf( f, "%lf\n", (double) time);
		for( i = 0 ; i < n ; i++)
			fprintf( f, "%lf\n", (double) sol.a[i]);
		
		fclose(f);
		return 1;
	}
	else{
		printf("Could not open file %s in write_solution\n",fn);
		return 0;
	}
}


/* Allocates and initializes the variables */
void init_solution(int n, solution *sol){

	/* allocate memory for solution */
	sol->a = (FTYPE*) malloc(n*sizeof(FTYPE));
	sol->kp = (INT*) malloc(sizeof(FTYPE));
	sol->status = (INT*) malloc(sizeof(FTYPE));

}

void free_solution(solution *sol){
  free(sol->a);
  free(sol->kp);
  free(sol->status);
}
