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


#ifndef __DR_SLP_H__
#define __DR_SLP_H__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <complex.h>
#include <fftw3.h>

#include "ftype.hpp"
#include "tools.hpp"

#define MIN(A,B) (A < B ? A : B)
#define MAX(A,B) (A < B ? B : A)

void adjointfilter(FTYPE * x, FTYPE * y, int k, FTYPE *a, int l);
void filter(FTYPE * a, int n, FTYPE * x, int k, FTYPE *y);
void autocorr(FTYPE * x, int k, FTYPE * r, int n);
void S(FTYPE * x, int n, FTYPE t, FTYPE * y);
void levinson(int n, FTYPE * a, FTYPE * b, FTYPE * y);
FTYPE norm1(int n, FTYPE * x);
FTYPE * vector(int n, char * s);


class fastautocorr{
  FFTW_PLAN p1, p2;
  FFTW_COMPLEX *X;
  FTYPE *X2, *xp, *a;
  int N, Np;
  public:
   fastautocorr(int k);
   ~fastautocorr();
   void autocorr(FTYPE * x, int k, FTYPE * r, int n);  
};

class fftfilter{
    FFTW_PLAN p1, p2, p3;
    FFTW_COMPLEX *H, *X, *Y;
    FTYPE *yp, *ap, *xp;
    int i, N, Np;
    int n, k;
    FTYPE s;

 public:
    fftfilter(int n, FTYPE * x, int k);
    ~fftfilter();
    void set_signal(FTYPE *x);
    void filter(FTYPE *a, FTYPE *y);
};

class fftadjointfilter{
    FFTW_PLAN p1, p2, p3;
    FFTW_COMPLEX *H, *X, *Y;
    FTYPE *yp, *ap, *xp;
    int i, N, Np;
    int l, k;
    FTYPE s;

 public:
    fftadjointfilter(int l, FTYPE * x, int k);
    ~fftadjointfilter();
    void set_signal(FTYPE *x);
    void filter(FTYPE *y, FTYPE *a);
};


class flevinson{
  FFTW_PLAN pF0, pF1, pFb, pT1, pT2, pH1, pH2, pP;
  FFTW_COMPLEX *F0, *F1, *F, *Fb, *T1, *T2, *H1, *H2, *P;
  FTYPE *t1, *t2, *h1, *h2, *z, *f, *p;
  FTYPE s, inv_delta_k;
  int i, N, Np, Npp, n;

 public:
  FTYPE *r;

  flevinson(int nr){
    n = nr - 1;
    N = nr;

    // Find a suffiently large Np that makes the Toeplitz matrices T0,T1 cyclic
    //Np = 2;
    //while(1){
    //  if ( Np >= n+1 )
    //    break;
    //  Np *= 2;
    //}
    //Np *= 2;
    Np = 2*(n+1);

    Npp = Np/2 + 1;

    z = vector(n, "z");
    r = vector(n, "r");

    F0 = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * Npp);
    F1 = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * Npp);
    Fb = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * Npp);

    T1 = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * Npp);
    T2 = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * Npp);

    H1 = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * Npp);
    H2 = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * Npp);

    P = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * Npp);


    f = (FTYPE*) FFTW_MALLOC(sizeof(FTYPE) * Np);
    t1 = (FTYPE*) FFTW_MALLOC(sizeof(FTYPE) * Np);
    t2 = (FTYPE*) FFTW_MALLOC(sizeof(FTYPE) * Np);

    h1 = (FTYPE*) FFTW_MALLOC(sizeof(FTYPE) * Np);
    h2 = (FTYPE*) FFTW_MALLOC(sizeof(FTYPE) * Np);
    
    p = (FTYPE*) FFTW_MALLOC(sizeof(FTYPE) * Np);

    pF0 = FFTW_PLAN_DFT_R2C_1D(Np, f, F0, FFTW_MEASURE);
    pF1 = FFTW_PLAN_DFT_R2C_1D(Np, f, F1, FFTW_MEASURE);
    pFb = FFTW_PLAN_DFT_R2C_1D(Np, f, Fb, FFTW_MEASURE);

    pT1 = FFTW_PLAN_DFT_C2R_1D(Np, T1, t1, FFTW_MEASURE);
    pT2 = FFTW_PLAN_DFT_R2C_1D(Np, t2, T2, FFTW_MEASURE);

    pH1 = FFTW_PLAN_DFT_C2R_1D(Np, H1, h1, FFTW_MEASURE);
    pH2 = FFTW_PLAN_DFT_R2C_1D(Np, h2, H2, FFTW_MEASURE);

    pP = FFTW_PLAN_DFT_C2R_1D(Np, P, p, FFTW_MEASURE);

    s = 1.0/Np;
  }

  ~flevinson(){
    free(r);
    free(z);

    FFTW_FREE(F0);
    FFTW_FREE(F1);
    FFTW_FREE(f);
    FFTW_FREE(Fb);

    FFTW_FREE(t1);
    FFTW_FREE(t2);
    FFTW_FREE(T1);
    FFTW_FREE(T2);
    FFTW_FREE(P);

    FFTW_FREE(h1);
    FFTW_FREE(h2);
    FFTW_FREE(H1);
    FFTW_FREE(H2);
    FFTW_FREE(p);

    FFTW_DESTROY_PLAN(pF0);
    FFTW_DESTROY_PLAN(pF1);
    FFTW_DESTROY_PLAN(pFb);

    FFTW_DESTROY_PLAN(pT1);
    FFTW_DESTROY_PLAN(pT2);
    FFTW_DESTROY_PLAN(pH1);
    FFTW_DESTROY_PLAN(pH2);

    FFTW_DESTROY_PLAN(pP);

  }

  void update_mu(FTYPE * mu);
  void solve(FTYPE * b, FTYPE * x);
};


/* Subclass this abstract class to define what the ADMM solves
 */
class spec_admm {
 public:
  int N;
  int kmax;
  int k;
  FTYPE rho;
  int verbose;
  FTYPE epsilon;
  FTYPE * t1;
  FTYPE * t2;
  FTYPE * z;
  FTYPE * yp;

  virtual void proxth(FTYPE * z, FTYPE * y) =0;
  virtual void PQ(FTYPE * x) =0;
  virtual FTYPE g(FTYPE * z) =0;

  void set_rho(FTYPE r){rho = r;};
  void set_verbose(int r){verbose = r;};
  void set_kmax(int r){kmax = r;};
  void set_epsilon(FTYPE r){epsilon = r;};
};


/* Subclass this abstract class to define what the DR solves
 */
class spec_dr {
 public:
  int N;
  int kmax;
  int k;
  FTYPE rho;
  int verbose;
  FTYPE epsilon;
  FTYPE t;
  FTYPE * zp;
  FTYPE * u;
  FTYPE * y;
  
  virtual void proxth(FTYPE * z, FTYPE * y) =0;
  virtual void PQ(FTYPE * x) =0;
  virtual FTYPE g(FTYPE * z) =0;

  void set_rho(FTYPE r){rho = r;};
  void set_verbose(int r){verbose = r;};
  void set_t(FTYPE r){t = r;};
  void set_kmax(int r){kmax = r;};
  void set_epsilon(FTYPE r){epsilon = r;};

};

/*
   This class defines functions for solving the least 1-norm problem

   minimize   ||x||_1
   subject to A*x == b

   using the Douglas-Rachford splitting method dr.
*/

class leastnorm : public spec_dr{

 public:
  FTYPE * A;
  FTYPE * b;
  int m, n;

  leastnorm(int, int, FTYPE *, FTYPE *);
  virtual ~leastnorm();

  void proxth(FTYPE * z, FTYPE * y){S(z, n, 1.0, y);};
  void PQ(FTYPE * x){
    FTYPE * temp;
    FTYPE * C;
    int info;
    int one = 1;

    temp = vector(m, "temp");

    C = vector(m*m, "C");
    
    /* Calc x <= x + A'*(A*A')\(b-A*x) */
    
    // b - Ax
    tools_copy(m, b, temp);
    tools_gemv('N', m, n, -1.0, A, m, x, one, 1.0, temp, one);
      

    //C = A*A'
    tools_gemm('N', 'T', m, m, n, 
               1.0, A, m, A, m, 0.0, C, m);
    
    // Factor
    tools_potrf('L', m, C, m, &info);

    //Solve
    tools_potrs('L', m, 1, C, m, temp, m, &info);

    // y <- x + A' temp
    tools_gemv('T', m, n, 1.0, A, m, temp, one, 1.0, x, one);

    free(temp);
    free(C);
  };

  FTYPE g(FTYPE * z){return tools_asum(n, z);};

};

/*
   This class defines functions for solving the sparse linear prediction problem

   minimize   ||x - X\alpha||_1 + \gamma||\alpha||_1

   using the Douglas-Rachford splitting method dr.
*/

class slp_dr : public spec_dr{

 public:
  FTYPE * x;
  FTYPE * X1;
  FTYPE gamma;
  FTYPE * r;
  FTYPE * a;
  FTYPE s;
  FTYPE * temp;
  int m, n;

  int solver;

  fftfilter * h;
  fftadjointfilter * ha;
  flevinson * l;
  fastautocorr * f;

  slp_dr(int, int, int, FTYPE *, FTYPE);
  virtual ~slp_dr();

  void set_signal(FTYPE * x);

  void proxth(FTYPE * z, FTYPE * y){
    // in Matlab notation 
    // y = proxth = @(xb) [Sm(z(1:n), t*gamma); x-Sm(x-z(n+1:end), t*1)];
    S(z, n, t*gamma, y);
    tools_copy(m, x, &y[n]);
    tools_axpy(m, -1.0, &z[n], &y[n]);
    S(&y[n], m, t, &y[n]);
    tools_axpy(m, -1.0, x, &y[n]);
    tools_scal(m, -1.0, &y[n]);
  };


  void PQ(FTYPE * x){
    /* Calc x <= [ I ] 
                 [ X ] (I + X'*X)^-1 (x1 + X'x2) 
    */
    
    // temp = X'*x2
    //adjointfilter(X1, &x[n], m, temp, n);
    ha->filter(&x[n], temp);

    // temp = x1 + X'x2
    tools_axpy(n, 1.0, x, temp);

    // x <= (I + X'*X)^-1 (x1 + X'x2) 
    if(solver == 0){
      levinson(n, a, temp, x);
      tools_scal(n, 1/s, x); //since the coeffients a are normalized
    }
    else{
      l->solve(temp, x);
    }
    
    //filter(x, n-1, X1, m, temp);
    h->filter(x, temp);
    
    tools_copy(m, temp, &x[n]);
  };

  FTYPE g(FTYPE * z){
    h->filter(z, temp);
    tools_axpy(m, -1.0, x, temp);
    return tools_asum(m, temp) + gamma*tools_asum(n, z);
  };
  
};


/*
   This class defines functions for solving the sparse linear prediction problem

   minimize   ||x - X\alpha||_1 + \gamma||\alpha||_1

   using the ADMM method.
*/

class slp_admm : public spec_admm{

 public:
  FTYPE * x;
  FTYPE * X1;
  FTYPE gamma;
  FTYPE * r;
  FTYPE * a;
  FTYPE s;
  int m, n;
  int solver;

  fftfilter * h;
  fftadjointfilter * ha;
  flevinson * l;
  fastautocorr * f;

  slp_admm(int, int, int, FTYPE *, FTYPE);
  virtual ~slp_admm();

  void set_signal(FTYPE * x);

  void proxth(FTYPE * z, FTYPE * y){
    // y = S(z, 1/rho)
    S(z, N, 1/rho, y);
  };


  void PQ(FTYPE * z){
    /* Calc a <= (X^TX+gamma^2)^{-1} ( X^T(x-z_2) + gamma z_1 )
       */
    FTYPE * z1 = z;
    FTYPE * z2 = z + n;
      
    // temp = x1 + X'x2
    tools_copy(m, x, t1);
    tools_axpy(m, -1.0, z2, t1);
    
    //adjointfilter
    ha->filter(t1, t2);
    
    tools_axpy(n, gamma, z1, t2);

    // x <= (gamma^2 I + X'*X)^-1 t2
    if(solver == 0){
      levinson(n, a, t2, z);
      tools_scal(n, 1/s, z); //since the coeffients a are normalized
    }
    else{
      l->solve(t2, z);
    }
    
    /* z  = [gamma alpha; e] */
    h->filter(z1, t1);
    tools_copy(m, x, z2);
    tools_axpy(m, -1.0, t1, z2);

    tools_scal(n, gamma, z);

  };

  FTYPE g(FTYPE * z){
    tools_copy(n, z, t1);
    tools_scal(n, 1/gamma, t1);
    h->filter(t1, t2);
    tools_axpy(m, -1.0, x, t2);
    return tools_asum(m, t2) + gamma*tools_asum(n, t1);
  };
  
};


void dr(spec_dr * s, FTYPE * z);
void admm(spec_admm * s, FTYPE * y, FTYPE * u);
#endif //__DR_SLP_H__
