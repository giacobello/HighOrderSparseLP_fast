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


#include "dr_admm_slp.hpp"


/* To use this, subclass class spec, implement the neccessary setup
   and call the function. Initialized using z */
void dr(spec_dr * s, FTYPE * z){
  int k;
  FTYPE nek;

  //unpacking - makes it more readable
  FTYPE * u = s->u;
  FTYPE * zp = s->zp;
  FTYPE * y = s->y;
  int kmax = s->kmax;
  int n = s->N;
  FTYPE epsilon = s->epsilon;
  FTYPE rho = s->rho;
  int verbose = s->verbose;

  for( k = 0 ; k < kmax ; ++k ){
    // u = proxth(z);
    s->proxth(z, u);  

    // y = PQ(2*u - z);
    tools_copy(n, u, y);
    tools_scal(n, 2.0, y);
    tools_axpy(n, -1.0, z, y);
    s->PQ(y);

    // zp = z + rho*(y - u);
    tools_axpy(n, -1.0, u, y);
    tools_copy(n, z, zp);
    tools_axpy(n, rho, y, zp);

    // evaluate objective
    if(verbose)
      printf("k=%d, h(u)=%.4e \n", k, s->g(u));

    tools_axpy(n, -1.0, zp, z);
    nek = tools_nrm2(n, z);
    
    // test stopping condition
    if(nek*nek <= epsilon*n){
      break;
    }

    // update via a copy
    tools_copy(n, zp, z);
  }

  
  s->k = k;

}

/* To use this, subclass class spec, implement the neccessary setup
   and call the function. Initialized using z 

   Special version for problems on the form 

   minimize f(x)
   subject x in Q

*/

void admm(spec_admm * s, FTYPE * y, FTYPE * u){
  int k;
  FTYPE nrk = 0;
  FTYPE nsk = 0;

  //unpacking - makes it more readable
  FTYPE * z = s->z;
  FTYPE * yp = s->yp;
  FTYPE * t1 = s->t1;
  FTYPE rho = s->rho;
  int kmax = s->kmax;
  int n = s->N;
  FTYPE epsilon = s->epsilon;
  int verbose = s->verbose;

  if(verbose){
    printf("---- ADMM -------------- \n");
    printf("n = %d\n", n);
  }

  for( k = 0 ; k < kmax ; ++k ){
    // z^(k+1) = PC(y^(k) - u^(k));
    tools_copy(n, y, z);
    tools_axpy(n, -1.0, u, z);
    s->PQ(z);

    // y^(k+1) = proxth(z^(k+1) + u^(k));
    tools_copy(n, u, t1);
    tools_axpy(n, 1.0, z, t1);
    s->proxth(t1, yp);
    
    // up = u + z^(k+1) - y^(k+1);
    tools_axpy(n, 1.0, z, u);
    tools_axpy(n, -1.0, yp, u);
    
    // Calculate stopping criterias
    tools_copy(n, z, t1);
    tools_axpy(n, -1.0, yp, t1);
    nrk = tools_nrm2(n, t1);
    
    tools_copy(n, y, t1);
    tools_axpy(n, -1.0, yp, t1);
    nsk = tools_nrm2(n, t1);

    //update
    tools_copy(n, yp, y);

    
    if(verbose){
      printf("k=%d, f(y)=%.4e primal %f <= %f, dual %f <= %f\n", 
             k, s->g(y), nrk*nrk, epsilon*n, rho*nsk*nsk, epsilon*n);
    }
    
    // test stopping condition
    if(nrk*nrk <= epsilon*n && rho*nsk*nsk <= epsilon*n){
      break;
    }

  }
  
  s->k = k;

}

leastnorm::leastnorm(int mr, int nr, FTYPE * Ar, FTYPE * br){

  m = mr;
  n = nr;

  N = n;

  zp = vector(N, "spec");
  u = vector(N, "spec");
  y = vector(N, "spec");

  A = Ar;
  b = br;
  t = 0.1;
  rho = 1.8;
  kmax = 2000;
  verbose = 0;
  epsilon = 1e-4;
}

leastnorm::~leastnorm(){
  free(zp);
  free(u);
  free(y);
}


slp_dr::slp_dr(int mr, int nr, int sol, FTYPE * xr, FTYPE gammar){
  
  m = mr;
  n = nr;
  X1 = xr;
  gamma = gammar;
  solver = sol;

  N = m + n;

  zp = vector(N, "slp");
  u = vector(N, "slp");
  y = vector(N, "slp");
  temp = vector(N, "slp");
  r = vector(n, "slp");
  a = vector(n, "slp");
  x = vector(m, "slp");

  tools_copy(m-1, X1+1, x);
  x[m-1] = 0.0;

  h = new fftfilter(n-1, X1, m);
  ha = new fftadjointfilter(n, X1, m);
  

  /* length of x is k, autocorr order is n, r has then the length n+1 */
  if(solver==1){
    f = new fastautocorr(m);
    f->autocorr(x, m, r, n-1);
  }
  else
    autocorr(x, m, r, n-1);

  tools_copy(n, r, a);
  a[0] += 1.0;
    
  if(solver==1){
    l = new flevinson(n);
    l->update_mu(a);
  }
  else{
    /* calculate the normalized coefficients a 
       of the symmetric Toeplitz matrix */
    s = a[0];
    tools_scal(n, 1/s, a);
  }
  
  t = 0.1;
  rho = 1.8;
  kmax = 2000;
  verbose = 0;
  epsilon = 1e-4;
}

void slp_dr::set_signal(FTYPE * xr){
    X1 = xr;
    tools_copy(m-1, X1+1, x);
    x[m-1] = 0.0;

    h->set_signal(X1);
    ha->set_signal(X1);
    
    /* length of x is k, autocorr order is n, r has then the length n+1 */
    if(solver==1)
      f->autocorr(x, m, r, n-1);
    else
      autocorr(x, m, r, n-1);

    tools_copy(n, r, a);
    a[0] += 1.0;
    
    if(solver==1)
      l->update_mu(a);
    else{
      s = a[0];
      tools_scal(n, 1/s, a);
    }

    /* calculate the normalized coefficients a of the symmetric Toeplitz matrix */
}

slp_dr::~slp_dr(){
  free(zp);
  free(u);
  free(y);
  free(temp);
  free(a);
  free(r);
  free(x);
  delete h;
  delete ha;
  if(solver==1){
    delete l;
    delete f;
  }


}

slp_admm::slp_admm(int mr, int nr, int sol, FTYPE * xr, FTYPE gammar){
  
  m = mr;
  n = nr;
  X1 = xr;
  gamma = gammar;
  solver = sol;

  N = m + n;

  yp = vector(N, "slp");
  z = vector(N, "slp");
  t1 = vector(N, "slp");
  t2 = vector(N, "slp");
  r = vector(n, "slp");
  a = vector(n, "slp");
  x = vector(m, "slp");

  tools_copy(m-1, X1+1, x);
  x[m-1] = 0.0;

  h = new fftfilter(n-1, X1, m);
  ha = new fftadjointfilter(n, X1, m);
  
  /* length of x is k, autocorr order is n, r has then the length n+1 */
  if(solver==1){
    f = new fastautocorr(m);
    f->autocorr(x, m, r, n-1);
  }
  else
    autocorr(x, m, r, n-1);

  tools_copy(n, r, a);
  a[0] += gamma*gamma;
    
  if(solver==1){
    l = new flevinson(n);
    l->update_mu(a);
  }
  else{
    /* calculate the normalized coefficients a of the symmetric Toeplitz matrix */
    s = a[0];
    tools_scal(n, 1/s, a);
  }
  
  rho = 100.0;
  kmax = 2000;
  verbose = 0;
  epsilon = 1e-4;
}

void slp_admm::set_signal(FTYPE * xr){
    X1 = xr;
    tools_copy(m-1, X1+1, x);
    x[m-1] = 0.0;

    h->set_signal(X1);
    ha->set_signal(X1);
    
    /* length of x is k, autocorr order is n, r has then the length n+1 */
    if(solver==1)
      f->autocorr(x, m, r, n-1);
    else
      autocorr(x, m, r, n-1);

    tools_copy(n, r, a);
    a[0] += gamma*gamma;
    
    if(solver==1)
      l->update_mu(a);
    else{
      s = a[0];
      tools_scal(n, 1/s, a);
    }

    /* calculate the normalized coefficients a of the symmetric Toeplitz matrix */
}

slp_admm::~slp_admm(){
  free(yp);
  free(z);
  free(t1);
  free(t2);
  free(a);
  free(r);
  free(x);
  delete h;
  delete ha;
  if(solver==1){
    delete l;
    delete f;
  }
}

// Update with a new set of coefficients
void flevinson::update_mu(FTYPE * mu){
  FTYPE delta_k, gamma_kp1;
  int k, l;

  delta_k = mu[0];
  gamma_kp1 = -mu[1]/mu[0];
  r[0] = gamma_kp1;
  delta_k = (1-gamma_kp1*gamma_kp1)*delta_k;

  for( k = 1 ; k < n ; k++ ){
    gamma_kp1 = -(tools_dot(k, r, mu+1) + mu[k+1])/delta_k;

    for( l = 0 ; l < k  ; l++ )
      z[l] = r[l] + gamma_kp1 * r[k-l-1];
      
    tools_copy(k, z, r+1);
    r[0] = gamma_kp1;

    delta_k = delta_k*(1-gamma_kp1*gamma_kp1);
  }


  // --------------------------
  f[0] = 1.0;
  for( l = 0 ; l < n ; l++)
    f[l+1] = r[n-l-1];
  for( l = n +1 ; l < Np ; l++)
    f[l] = 0;
  
  FFTW_EXECUTE(pF0);

  f[0] = 1.0;
  for( l = 1 ; l < Np ; l++)
    f[l] = 0;
  tools_copy(n, r, f + Np - n);
  FFTW_EXECUTE(pF1);

  //for( i = 0 ; i < Np ; i++)
  //  printf("%d: %5.2f +j%5.2f\n", i, creal(F1p[i]), cimag(F1p[i]));

  inv_delta_k = 1.0/delta_k;
}


void flevinson::solve(FTYPE * b, FTYPE * x){
  int i;

  tools_copy(N, b, f);
  FFTW_EXECUTE(pFb);

  // --------------------------------------
  tools_mulz(Npp, (FTYPE*)Fb, (FTYPE*)F0, (FTYPE*)T1);

  FFTW_EXECUTE(pT1);

  tools_copy(N, t1, t2);
  tools_scal(N, s, t2);
  for( i = N ; i < Np ; i++ )
      t2[i] = 0.0;

  FFTW_EXECUTE(pT2);
  tools_mulz(Npp, (FTYPE*)F1, (FTYPE*)T2, (FTYPE*)T1);
  
  // --------------------------------
  // --------------------------------------
  tools_mulz(Npp, (FTYPE*)Fb, (FTYPE*)F1, (FTYPE*)H1);
  FFTW_EXECUTE(pH1);

  tools_copy(N, h1+N, h2);
  tools_scal(N, s, h2);
  for( i = N ; i < Np ; i++ )
      h2[i] = 0.0;

  FFTW_EXECUTE(pH2);
  tools_mulz(Npp, (FTYPE*)F0, (FTYPE*)H2, (FTYPE*)H1);

  tools_subz(Npp, (FTYPE*)T1, (FTYPE*)H1, (FTYPE*)P, 2);
  //for( i = 0 ; i < Npp ; i+=2 ){
  //  P[i][0] = T1[i][0] - H1[i][0];
  //  P[i][1] = T1[i][1] - H1[i][1];
  //}

  tools_addz(Npp, (FTYPE*)(T1+1), (FTYPE*)(H1+1), (FTYPE*)(P+1), 2);
  //for( i = 1 ; i < Npp ; i+=2 ){
  //  P[i][0] = T1[i][0] + H1[i][0];
  //  P[i][1] = T1[i][1] + H1[i][1];
  //}

  FFTW_EXECUTE(pP);

  tools_copy(N, p, x);
  tools_scal(N, s*inv_delta_k, x);
}


/* Adjoint FIR filter with all-zero initialization
   If filter is y <- Xa
   then the adjoint filter implements a <- X^T y

  x has length k
  y has lenght k
  a has length l 
*/
void adjointfilter(FTYPE * x, FTYPE * y, int k, FTYPE *a, int l){
  
  int i, j;

  for( i = 0 ; i < l ; i++ ){
    a[i] = 0.0;
    
    for( j = i ; j < k ; j++ ){
      a[i] += y[j]*x[j-i];
    }
  }
}

/* FIR filter with all-zero initialization
  Identical to Matlab y = filter(a, 1, x) 
  a has length n+1 (order n)
  x has length k
  y has lenght k
*/
/*void filter(FTYPE * a, int n, FTYPE * x, int k, FTYPE *y){
  int i, j;
  int m;

  for( i = 0 ; i < k ; i++ ){
    m = MIN(i, n);

    y[i] = 0.0;
    for( j = 0 ; j <= m ; j++ ){
      y[i] += a[j]*x[i-j];
    }
    
  }
  }*/

void filter(FTYPE * a, int n, FTYPE * x, int k, FTYPE *y){
  
  FFTW_PLAN p1, p2, p3;
  FFTW_COMPLEX *H, *X, *Y;
  FTYPE *yp, *ap, *xp;
  int i, N, Np;
  FTYPE s;

  N = 2;
  while(1){
    if ( N >= k && N>=n+1)
      break;
    N = N*2;
  }

  Np = N/2 + 1;
  s = 1.0/N;

  H = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * Np);
  X = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * Np);
  Y = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * Np);
  yp = (FTYPE*) FFTW_MALLOC(sizeof(FTYPE) * N);
  ap = (FTYPE*) FFTW_MALLOC(sizeof(FTYPE) * N);
  xp = (FTYPE*) FFTW_MALLOC(sizeof(FTYPE) * N);

  tools_copy(n+1, a, ap);
  for( i = n+1 ; i < N ; i++ )
    ap[i] = 0.0;

  tools_copy(k, x, xp);
  for( i = k ; i < N ; i++ )
    xp[i] = 0.0;

  p1 = FFTW_PLAN_DFT_R2C_1D(N, ap, H, FFTW_MEASURE);
  p2 = FFTW_PLAN_DFT_R2C_1D(N, xp, X, FFTW_MEASURE);
  p3 = FFTW_PLAN_DFT_C2R_1D(N, Y, yp, FFTW_MEASURE);

  FFTW_EXECUTE(p1);
  FFTW_EXECUTE(p2);

  tools_mulz(Np, (FTYPE*)H, (FTYPE*)X, (FTYPE*)Y);
  //for( i = 0 ; i < Np ; i++ )
  //  Y[i] = H[i]*X[i];

  FFTW_EXECUTE(p3);

  //for( i = 0 ; i < Np ; i++)
  //  printf("%d: %5.2f +j%5.2f\n", i, creal(H[i]), cimag(H[i]));

  for( i = 0 ; i < k ; i++)
    y[i] = s*yp[i];

  FFTW_DESTROY_PLAN(p1);
  FFTW_DESTROY_PLAN(p2);
  FFTW_DESTROY_PLAN(p3);
  FFTW_FREE(H); 
  FFTW_FREE(X); 
  FFTW_FREE(Y); 
  FFTW_FREE(yp); 
  FFTW_FREE(ap); 
  FFTW_FREE(xp);  
}

fftfilter::fftfilter(int nr, FTYPE * x, int kr){

  n = nr;
  k = kr;

  N = 2;
  while(1){
    if ( N >= 2*k-1 && N>=2*(n+1)-1)
      break;
    N = N*2;
  }

  Np = N/2 + 1;
  s = 1.0/N;

  H = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * Np);
  X = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * Np);
  Y = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * Np);
  yp = (FTYPE*) FFTW_MALLOC(sizeof(FTYPE) * N);
  ap = (FTYPE*) FFTW_MALLOC(sizeof(FTYPE) * N);
  xp = (FTYPE*) FFTW_MALLOC(sizeof(FTYPE) * N);

  p1 = FFTW_PLAN_DFT_R2C_1D(N, ap, H, FFTW_MEASURE);
  p2 = FFTW_PLAN_DFT_R2C_1D(N, xp, X, FFTW_MEASURE);
  p3 = FFTW_PLAN_DFT_C2R_1D(N, Y, yp, FFTW_MEASURE);

  set_signal(x);
}

void fftfilter::set_signal(FTYPE * x){
  tools_copy(k, x, xp);
  for( i = k ; i < N ; i++ )
    xp[i] = 0.0;

  FFTW_EXECUTE(p2);
}

fftfilter::~fftfilter(){
  FFTW_DESTROY_PLAN(p1);
  FFTW_DESTROY_PLAN(p2);
  FFTW_DESTROY_PLAN(p3);
  FFTW_FREE(H); 
  FFTW_FREE(X); 
  FFTW_FREE(Y); 
  FFTW_FREE(yp); 
  FFTW_FREE(ap); 
  FFTW_FREE(xp);  
}

void fftfilter::filter(FTYPE * a, FTYPE *y){

  tools_copy(n+1, a, ap);
  for( i = n+1 ; i < N ; i++ )
    ap[i] = 0.0;

  FFTW_EXECUTE(p1);

  tools_mulz(Np, (FTYPE*)H, (FTYPE*)X, (FTYPE*)Y);

  FFTW_EXECUTE(p3);


  for( i = 0 ; i < k ; i++)
    y[i] = s*yp[i];

}

fftadjointfilter::fftadjointfilter(int lr, FTYPE * x, int kr){

  l = lr;
  k = kr;

  N = 2;
  while(1){
    if ( N >= 2*k-1 && N>=2*l-1 )
      break;
    N = N*2;
  }

  Np = N/2 + 1;
  s = 1.0/N;

  H = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * Np);
  X = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * Np);
  Y = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * Np);
  yp = (FTYPE*) FFTW_MALLOC(sizeof(FTYPE) * N);
  ap = (FTYPE*) FFTW_MALLOC(sizeof(FTYPE) * N);
  xp = (FTYPE*) FFTW_MALLOC(sizeof(FTYPE) * N);

  p1 = FFTW_PLAN_DFT_C2R_1D(N, H, ap, FFTW_MEASURE);
  p2 = FFTW_PLAN_DFT_R2C_1D(N, xp, X, FFTW_MEASURE);
  p3 = FFTW_PLAN_DFT_R2C_1D(N, yp, Y, FFTW_MEASURE);

  set_signal(x);
}

void fftadjointfilter::set_signal(FTYPE * x){
  int i;

  for( i = 0 ; i < k ; i++ )
    xp[i] = x[k-1-i]; //flipped coefficients
  for( i = k ; i < N ; i++ )
    xp[i] = 0.0;

  FFTW_EXECUTE(p2);
}

fftadjointfilter::~fftadjointfilter(){
  FFTW_DESTROY_PLAN(p1);
  FFTW_DESTROY_PLAN(p2);
  FFTW_DESTROY_PLAN(p3);
  FFTW_FREE(H); 
  FFTW_FREE(X); 
  FFTW_FREE(Y); 
  FFTW_FREE(yp); 
  FFTW_FREE(ap); 
  FFTW_FREE(xp);  
}

void fftadjointfilter::filter(FTYPE *y, FTYPE * a){

  tools_copy(k, y, yp);
  for( i = k ; i < N ; i++ )
    yp[i] = 0.0;

  FFTW_EXECUTE(p3);

  tools_mulz(Np, (FTYPE*)X, (FTYPE*)Y, (FTYPE*)H);


  //for( i = 0 ; i < N ; i++ )
  //  printf("%d: %5.2f +j%5.2f\n", i, creal(H[i]), cimag(H[i]));

  FFTW_EXECUTE(p1);

  for( i = k - 1; i < k + l - 1; i++ )
    a[i-k+1] = s*ap[i];
}

/* Autocorrelation - biased - calculated the direct way */
/* length of x is k, autocorr order is n, r has then the length n+1 */
void autocorr(FTYPE * x, int k, FTYPE * r, int n){
  
  int i;

  for( i = 0 ; i <= n ; i++ ){
    
    r[i] = tools_dot(k-i, x, x+i);
    //r[i] = 0.0;
    //for(int j = 0 ; j < k - i; j++ ){
    //  r[i] += x[j]*x[j+i];
    //}
  }

}

fastautocorr::fastautocorr(int k){


  N = 2;
  while(1){
    if ( N >= 2*k-1)
      break;
    N = N*2;
  }

  Np = N/2 + 1;

  X = (FFTW_COMPLEX*) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * Np);
  X2 = (FTYPE*) FFTW_MALLOC(sizeof(FTYPE) * Np);
  a = (FTYPE*) FFTW_MALLOC(sizeof(FTYPE) * Np);
  xp = (FTYPE*) FFTW_MALLOC(sizeof(FTYPE) * N);

  p1 = FFTW_PLAN_DFT_R2C_1D(N, xp, X, FFTW_MEASURE);
  p2 = FFTW_PLAN_DFT_R2R_1D(Np, X2, a, FFTW_REDFT00, FFTW_MEASURE);
}

fastautocorr::~fastautocorr(){
  FFTW_DESTROY_PLAN(p1);
  FFTW_DESTROY_PLAN(p2);
  FFTW_FREE(X); 
  FFTW_FREE(X2); 
  FFTW_FREE(a); 
  FFTW_FREE(xp); 
}

/* Autocorrelation - biased - calculated using FFTs */
/* length of x is k, autocorr order is n, r has then the length n+1 */
void fastautocorr::autocorr(FTYPE * x, int k, FTYPE * r, int n){

  tools_copy(k, x, xp);
  for( int i = k ; i < N ; i++ )
    xp[i] = 0.0;

  FFTW_EXECUTE(p1);

  for( int i = 0 ; i < Np ; i++){
    X2[i] = X[i][0]*X[i][0] + X[i][1]*X[i][1];
  }

  FFTW_EXECUTE(p2);

  tools_copy(n+1, a, r);
  tools_scal(n+1, (FTYPE) 1.0/N, r);
}



/* Softhreshold of the vector x, length n using parameter t
   output is in y */
void S(FTYPE * x, int n, FTYPE t, FTYPE * y){

  int k;
  
  
  for( k = 0 ; k < n ; ++k ){
    if(x[k] > t)
      y[k] = x[k] - t;
    else if(x[k] < -t)
      y[k] = x[k] + t;
    else
      y[k] = 0.0;
  }
}

/* Levinson algorithm for unit diagonal symmetric Toeplitz matrix */
/* See S4.7.2 in Matrix Computations, G.H Golub and C.F. Van Loan, 4ed */
void levinson(int n, FTYPE * a, FTYPE * b, FTYPE * x){

  register FTYPE beta, alpha=0.0, mu, ax, ay;
  register FTYPE *y, *z;

  int k, l;
  
  z = vector(n, "z - levinson");
  y = vector(n, "y - levinson");

  /* 1) Initializations */
  beta = 1.0;
    
  x[0] = b[0];

  if( 0 < n-1 ){
    y[0] = -a[1];
    alpha = -a[1];
  }

  /* 2) Loop */
  for( k = 1 ; k < n ; k++ ){
    beta = (1.0 - alpha*alpha)*beta;
    
    
    ax = tools_dot(k, a+1, 1, x, -1);
    /*
    ax = 0.0;
    for( l = 1 ; l <= k ; l++ )
       ax += a[l] * x[k-l];
    */

    mu = (b[k] - ax)/beta;

    //tools_axpy(k, mu, y, -1, x, 1);
    for( l = 0 ; l < k ; ++l )
      x[l] += mu * y[k-1-l];
    

    x[k] = mu;
    
    if( k < n-1 ){
      /*
      ay = 0.0;
      for( l = 1 ; l <= k ; l++ )
        ay += a[l] * y[k-l];
      */
      ay = tools_dot(k, a+1, 1, y, -1);

      alpha = -(a[k+1] + ay)/beta;

      for( l = 0 ; l < k ; l++ )
        z[l] = y[l] + alpha * y[k-1-l];
      
      /* copy */
      tools_copy(k, z, y);
      
      y[k] = alpha;

    }
  }

  /* Free the temporary variable, not the solution */
  free(z);
  free(y);
}

FTYPE norm1(int n, FTYPE * x){
  int i;
  FTYPE v;

  v = 0;
  for( i = 0 ; i < n ; i++ )
    v += fabs(x[i]);

  return v;
}

FTYPE * vector(int n, char * s){
  FTYPE * x;

  x = (FTYPE*)malloc(n*sizeof(FTYPE));
  if(x == NULL){
    printf("%s - malloc failed\n", s);
    exit(1);
  }

  return x;
}
