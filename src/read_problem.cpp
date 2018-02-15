#include "read_problem.hpp"


void init_problem(int m, int n, isignal *sig){
		sig->Y = (FTYPE*) malloc(m*n*sizeof(FTYPE));
		sig->y = (FTYPE*) malloc(m*sizeof(FTYPE));
}


void free_problem(isignal *sig){
	free(sig->Y);
	free(sig->y);
}
	

int read_problem(char *s, isignal *sig){
	FILE *f;
	int i;
	double r;

	// open the file
	f = fopen(s, "r");
	if( f != NULL ){
		// first read m, n and gamma
		if(fscanf(f, "%d", &(sig->m)) != 1)
			return 0;
		if(fscanf(f, "%d", &(sig->n)) != 1)
			return 0;
		if(fscanf(f, "%lf", &r) != 1)
			 return 0;

		sig->gamma = (FTYPE)r;


		// allocate space for X and x
		init_problem(sig->m, sig->n, sig);

		// read in X and x
		for( i = 0 ; i < sig->m*sig->n ; ++i ){
			if(fscanf(f,"%lf", &r) != 1) return 0; else sig->Y[i]=(FTYPE)r;
			 }
		for( i = 0 ; i < sig->m ; ++i ){
			if(fscanf(f, "%lf", &r) != 1) return 0; else sig->y[i]=(FTYPE)r;

		}
		// close the file
		fclose(f);
		
		return 1;
	}
	else 
		printf("Could not load file %s\n", s);
		return 0;

}
