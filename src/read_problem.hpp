#ifndef __READ_PROBLEM_H__
#define __READ_PROBLEM_H__

#include <stdio.h>
#include <stdlib.h>

#include "ftype.hpp"

typedef struct{
	FTYPE *y;
	FTYPE *Y;
	INT m;
	INT n;
	FTYPE gamma;
} isignal;


void init_problem(int m, int n, isignal *sig);


void free_problem(isignal *sig);
	

int read_problem(char *s, isignal *sig);

#endif
