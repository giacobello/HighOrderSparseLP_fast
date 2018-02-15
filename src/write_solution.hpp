#ifndef __WRITE_SOLUTION_H__
#define __WRITE_SOLUTION_H__

#include <stdio.h>
#include <stdlib.h>

#include "ftype.hpp"

typedef struct{
	FTYPE *a;
	INT *kp;
	INT *status;
} solution;


int write_solution(char *fn, int n, solution sol, FTYPE time);

void init_solution(int n, solution *sol);

void free_solution(solution *sol);

#endif
