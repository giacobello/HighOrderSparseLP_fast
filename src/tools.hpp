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


#ifndef __TOOLS_H__
#define __TOOLS_H__

#include <stdio.h>
#include <math.h>

#include "ftype.hpp"

void print_vector(int size, FTYPE *a, char *s);

void print_matrix(int size_m, int size_n, FTYPE *X, char *s);

FTYPE tools_stepsize(INT size, FTYPE *dx, FTYPE *x);

FTYPE tools_stepsize_scale(INT size, FTYPE scale, FTYPE *dx, FTYPE *x);

void tools_axpy(INT size, FTYPE alpha, FTYPE *x, FTYPE *y);

void tools_axpy(INT size, FTYPE alpha, FTYPE *x, int incx, FTYPE *y, int incy);

void tools_copy(INT size, FTYPE *x, FTYPE *y);

void tools_scal(INT size, FTYPE alpha, FTYPE *x);

FTYPE tools_asum(INT size, FTYPE *x);

FTYPE tools_sum(INT size, FTYPE *x);

void tools_add(INT size, FTYPE alpha, FTYPE *x);

void tools_mulz(INT size, FTYPE *a, FTYPE *b, FTYPE *c);

void tools_addz(INT size, FTYPE *a, FTYPE *b, FTYPE *c, INT step);

void tools_subz(INT size, FTYPE *a, FTYPE *b, FTYPE *c, INT step);

void tools_set(INT size, FTYPE alpha, FTYPE *x);

FTYPE tools_dot(INT size, FTYPE *x, FTYPE *y);

FTYPE tools_dot(INT size, FTYPE *x, int incx, FTYPE *y, int incy);

FTYPE tools_nrm2(INT size, FTYPE *x);

void tools_m1(INT size, FTYPE *x, FTYPE *y);

void tools_ea(INT size, FTYPE *x, FTYPE *y);

FTYPE tools_inf(INT size, FTYPE *x);

void tools_hp(INT size, FTYPE *x, FTYPE *y, FTYPE *z);

void tools_md(INT size_m, INT size_n, FTYPE *X, FTYPE *d, FTYPE *Z);

void tools_mdt(INT size_m, INT size_n, FTYPE *X, FTYPE *d, FTYPE *Z);

void tools_ad(INT size, FTYPE *d, FTYPE *X);

void tools_ads(INT size, FTYPE alpha, FTYPE *X);

void tools_gemvp(char NT, INT m, INT n, FTYPE *A, FTYPE *x, FTYPE *b);

void tools_gemv(CONST char transa, CONST INT m, CONST INT n, 
                CONST FTYPE alpha, CONST FTYPE *a, CONST INT lda, 
                CONST FTYPE *x, CONST INT incx, CONST FTYPE beta, 
                FTYPE *y, CONST INT incy);

void tools_gemm(CONST char transa, CONST char transb, CONST INT m, 
                CONST INT n, CONST INT k,	CONST FTYPE alpha, CONST FTYPE *a, 
                CONST INT lda, CONST FTYPE *b, CONST INT ldb,
                CONST FTYPE beta, FTYPE *c, CONST INT ldc);

void tools_gesdd(CONST char jobz, CONST INT m, CONST INT n, FTYPE* a, 
                 CONST INT lda, FTYPE* s, FTYPE* u, CONST INT ldu,
                 FTYPE* vt, CONST INT ldvt, FTYPE* work,
                 CONST INT lwork, INTA* iwork, INTA *info);

void tools_potrf(CONST char uplo, CONST INT n, 
                 FTYPE* a, CONST INT lda, INTA* info);

void tools_potrs(CONST char uplo, CONST INT n, CONST INT nrhs, 
                 CONST FTYPE* a, CONST INT lda, FTYPE* b,
                 CONST INT ldb, INTA* info);


#endif
