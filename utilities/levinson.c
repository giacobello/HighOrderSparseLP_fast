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

#include <mex.h>
#include <stdlib.h>
#include <stdio.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  register double *a, *b, *x, *y, *z, beta, alpha, mu, ax, ay;
  mxArray *zp;
  int n, k, l;
  
  if(nrhs != 3)
    printf("Should contain 3 input parameters but has %i\n", nrhs);

  /* Pass input */
  zp = (mxArray*) prhs[0]; 
  n = (int)mxGetScalar(zp);

  zp = (mxArray*) prhs[1];
  a = mxGetPr(zp);
  
  zp = (mxArray*) prhs[2];
  b = mxGetPr(zp);
  
  /* Allocate memory and assign output pointer */
  plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL);
  x = mxGetPr(plhs[0]);

  y = malloc(n*sizeof(double));
  if(y == NULL){
    printf("levinson - malloc failed\n");
    exit(1);
  }

  z = malloc(n*sizeof(double));
  if(z == NULL){
    printf("levinson - malloc failed\n");
    exit(1);
  }

  /* Levinson algorithm for unit diagonal symmetric Toeplitz matrix */
  /* See S4.7.2 in Matrix Computations, G.H Golub and C.F. Van Loan, 4ed */

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
    
    ax = 0.0;
    for( l = 1 ; l <= k ; l++ )
      ax += a[l] * x[k-l];

    mu = (b[k] - ax)/beta;

    for( l = 0 ; l < k ; l++ )
      x[l] += mu * y[k-1-l];

    x[k] = mu;
    
    if( k < n-1 ){
      ay = 0.0;
      for( l = 1 ; l <= k ; l++ )
        ay += a[l] * y[k-l];

      alpha = -(a[k+1] + ay)/beta;

      for( l = 0 ; l < k ; l++ )
        z[l] = y[l] + alpha * y[k-1-l];
      
      /* copy */
      for( l = 0 ; l < k ; l++ )
        y[l] = z[l];
      
      y[k] = alpha;

    }
  }

  /* Free the temporary variable, not the solution */
  free(y);
  free(z);
}
