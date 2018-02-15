/*  
  Copyright 2014-2015 Tobias L. Jensen <tlj@es.aau.dk>
  Department of Electronic Systems, Aalborg University, Denmark

  This file is part of slp_sm.
    
  slp_sm is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  slp_sm is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with slp_sm.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "tools.hpp"

#define MIN(A,B) A < B ? A : B

/* Print the vector a of size, id first with char array s */
void print_vector(int size, FTYPE *a, char *s){
	int i;

	printf("%s\n",s);
	
	for( i = 0 ; i < size ; ++i )
		printf("%d %.4f\n",i,a[i]);
}


void print_matrix(int size_m, int size_n, FTYPE *X, char *s){
	int i,j;

	printf("%s\n",s);
	
	for( i = 0 ; i < size_m ; ++i ){
		for( j = 0 ; j < size_n ; ++j )
			printf("%.4f  ",X[i+j*size_m]);
		printf("\n");
	}
}

/* calculate the stepsize */
FTYPE tools_stepsize(INT size, FTYPE *dx, FTYPE *x){
	INT i;
	FTYPE alpha = 1.0;

	for( i = 0 ; i < size ; ++i ){
		if( dx[i] < 0 )
			alpha = MIN(alpha,-x[i]/dx[i]);
	}
	return alpha;
}

/* calculate the stepsize with a scaling factor*/
FTYPE tools_stepsize_scale(INT size, FTYPE scale, FTYPE *dx, FTYPE *x){
	INT i;
	FTYPE alpha = 1.0;

	for( i = 0 ; i < size ; ++i ){
		if( dx[i] < 0 )
			alpha = MIN(alpha,-scale*x[i]/dx[i]);
	}
	return alpha;
}


void tools_axpy(INT size, FTYPE alpha, FTYPE *x, FTYPE *y){
	INT one = 1;

	AXPY (&size, &alpha, x, &one, y, &one);
}


void tools_axpy(INT size, FTYPE alpha, FTYPE *x, int incx, FTYPE *y, int incy){

	AXPY (&size, &alpha, x, &incx, y, &incy);
}


void tools_copy(INT size, FTYPE *x, FTYPE *y){
	INT one = 1;

	COPY (&size, x, &one, y, &one);
}


void tools_scal(INT size, FTYPE alpha, FTYPE *x){
	INT one = 1;
 
	SCAL (&size, &alpha, x, &one);
}

/* add alpha to x, x_i = x_i + alpha */
void tools_add(INT size, FTYPE alpha, FTYPE *x){
	INT i;
 
	for( i = 0 ; i < size ; ++i){
		x[i] += alpha;
	}
}

/* element-wise complex multiplication c[i] = a[i]*b[i] */
void tools_mulz(INT size, FTYPE *a, FTYPE *b, FTYPE *c){

  #ifdef MKL
  VMUL (size, a, b, c);
  #else
  for( int i = 0 ; i < size*2 ; i+=2 ){
    c[i] = a[i]*b[i] - a[i+1]*b[i+1];
    c[i+1] = a[i]*b[i+1] + a[i+1]*b[i];
  }
  #endif
}

/* element-wise complex addition c[step*i] = a[step*i] + b[step*i] */
void tools_addz(INT size, FTYPE *a, FTYPE *b, FTYPE *c, INT step){

  for( int i = 0 ; i < size*2 ; i+=2*step ){
    c[i] = a[i] + b[i];
    c[i+1] = a[i+1] + b[i+1];
  }

}

/* element-wise complex substraction c[step*i] = a[step*i] - b[step*i] */
void tools_subz(INT size, FTYPE *a, FTYPE *b, FTYPE *c, INT step){

  for( int i = 0 ; i < size*2 ; i+=2*step ){
    c[i] = a[i] - b[i];
    c[i+1] = a[i+1] - b[i+1];
  }

}

/* set x_i to alpha, x_i = alpha */
void tools_set( INT size, FTYPE alpha, FTYPE *x){
	INT i;
 
	for( i = 0 ; i < size ; ++i){
		x[i] = alpha;
	}
}

FTYPE tools_dot(INT size, FTYPE *x, FTYPE *y){
	INT one = 1;

	return DOT (&size, x, &one, y, &one);

}

FTYPE tools_dot(INT size, FTYPE *x, int incx, FTYPE *y, int incy){

	return DOT (&size, x, &incx, y, &incy);
}

FTYPE tools_nrm2(INT size, FTYPE *x){
	INT one = 1;

	return NRM2 (&size, x, &one);
}

FTYPE tools_asum(INT size, FTYPE *x){
	INT one = 1;

	return ASUM (&size, x, &one);
}

FTYPE tools_sum(INT size, FTYPE *x){
	INT i = 1;
	FTYPE s = 0.0;

	for( i = 0 ; i < size ; ++i)
		s += x[i];
	
	return s;
}

/* calculates y_i = 1/x_i */
void tools_m1(INT size, FTYPE *x, FTYPE *y){
	INT i;

for( i = 0 ; i < size ; ++i)
	y[i] = 1.0/x[i];
}

/* calculates the Hadamard product z_i = x_i*y_i */
void tools_hp(INT size, FTYPE *x, FTYPE *y, FTYPE *z){
	INT i;

	for( i = 0 ; i < size ; ++i)
		z[i] = x[i]*y[i];
}

/* calculates the elementwise abs y_i = |x_i|*/
void tools_ea(INT size, FTYPE *x, FTYPE *y){
	INT i;

	for( i = 0 ; i < size ; ++i)
		y[i] = fabs(x[i]);
}

/* calculates the max abs of the elements in x, ||x||_\infty */
FTYPE tools_inf(INT size, FTYPE *x){
	INT one = 1;
	/*INT i;
	FTYPE max;

	max = fabs(x[0]);
	for( i = 1 ; i < size ; ++i){
		if( fabs(x[i]) > max){
			max = fabs(x[i]);
			}*/

	return fabs(x[ AMAX ( &size, x, &one) - 1 ]);
}

/* calculates the multiplication with a diagonal Z = X*diag(d) */
void tools_md(INT size_m, INT size_n, FTYPE *X, FTYPE *d, FTYPE *Z){
	INT i; /*,j,k;*/
	INT one=1;

	/*k = 0;
	for( i = 0 ; i < size_n ; ++i ){
		for( j = 0 ; j < size_m ; ++j){
			Z[k] = X[k]*d[i];
			++k;
		}
		}*/


	for( i = 0 ; i < size_n ; ++i ){
		COPY (&size_m, &X[i*size_m], &one, &Z[i*size_m], &one);
		SCAL (&size_m, &d[i], &Z[i*size_m], &one);
	}

}

/* calculates the multiplication with a diagonal Z = X'*diag(d), transpose */
void tools_mdt(INT size_m, INT size_n, FTYPE *X, FTYPE *d, FTYPE *Z){
	/*INT ij,k;*/
	INT j;

	/*k = 0;
	for( i = 0 ; i < size_n ; ++i ){
		for( j = 0 ; j < size_m ; ++j){
			Z[k] = X[k]*d[j];
			++k;
		}
	}*/		

	for( j = 0 ; j < size_m ; ++j ){
		COPY (&size_n, &X[j], &size_m, &Z[j], &size_m);
		SCAL (&size_n, &d[j], &Z[j], &size_m);
	}

}


/* calculates the addition of a diagonal X = X + diag(d) */
void tools_ad(INT size, FTYPE *d, FTYPE *X){
	INT one = 1;
	FTYPE done = 1.0;
	INT sizep1 = size + 1;

	/*INT i,k;
	k = 0;
	for( i = 0 ; i < size ; ++i){
		X[k] = X[k]+d[i];
		k += size + 1; 
	}*/

	AXPY ( &size, &done, d, &one, X, &sizep1);
}

/* calculates the addition of a scalar times identity X = X + alpha I */
void tools_ads(INT size, FTYPE alpha, FTYPE *X){
	INT i,k;
	
	k = 0;
	for( i = 0 ; i < size ; ++i){
		X[k] = X[k] + alpha;
		k += size + 1; 
	}
}

/* custom gemvp - own reference implementation */
void tools_gemvp(char NT, INT m, INT n, FTYPE *A, FTYPE *x, FTYPE *b){
	INT i,j;
	FTYPE s;

	if( NT == 'T'){ /*Transpose*/
		for( i = 0 ; i < n ; ++i){
			s = 0;
			for( j = 0 ; j < m ; ++j){
				s += A[j+i*m] * x[j]; 
			}
			b[i] = s;
		}
	}
	else{
		for( j = 0 ; j < m ; ++j )
			b[j] = 0;

		for( i = 0 ; i < n ; ++i){
			s = x[i];
			for( j = 0 ; j < m ; ++j){
				b[j] += A[j+i*m] * s; 
			}
		}
	}
}

void tools_gemv(CONST char transa, CONST INT m, CONST INT n, 
                 CONST FTYPE alpha, CONST FTYPE *a, CONST INT lda, 
                 CONST FTYPE *x, CONST INT incx, CONST FTYPE beta, 
                 FTYPE *y, CONST INT incy){

  GEMV (&transa, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);

}

void tools_gemm(CONST char transa, CONST char transb, CONST INT m, 
                CONST INT n, CONST INT k,	CONST FTYPE alpha, CONST FTYPE *a, 
                CONST INT lda, CONST FTYPE *b, CONST INT ldb,
                CONST FTYPE beta, FTYPE *c, CONST INT ldc){
  
   GEMM (&transa, &transb, &m, &n, &k, &alpha, a, 
         &lda, b, &ldb, &beta, c, &ldc);
}

void tools_gesdd(CONST char jobz, CONST INT m, CONST INT n, FTYPE* a, 
             CONST INT lda, FTYPE* s, FTYPE* u, CONST INT ldu,
             FTYPE* vt, CONST INT ldvt, FTYPE* work,
								 CONST INT lwork, INTA* iwork, INTA *info){
	
	GESDD (&jobz, &m, &n, a, &lda, s, u, &ldu, 
           vt, &ldvt, work, &lwork, iwork, info);
}

void tools_potrf(CONST char uplo, CONST INT n, FTYPE* a, 
									CONST INT lda, INTA* info){

	POTRF (&uplo, &n, a, &lda, info);
}

void tools_potrs(CONST char uplo, CONST INT n, CONST INT nrhs, 
             CONST FTYPE* a, CONST INT lda, FTYPE* b,
             CONST INT ldb, INTA* info ){

  POTRS (&uplo, &n, &nrhs, a, &lda, b, &ldb, info);
}
