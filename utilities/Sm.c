#include <mex.h>
#include <stdlib.h>
#include <stdio.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  register double *x, *y, t;
  mxArray *zp;
  int n, k;
  
  if(nrhs != 2){
    printf("Should contain 2 input parameters but has %i\n", nrhs);
    exit(1);
  }
    
  /* Pass input */
  zp = (mxArray*) prhs[0];
  x = mxGetPr(zp);

  n = mxGetM(zp);

  zp = (mxArray*) prhs[1];
  t = (double)mxGetScalar(zp);
   
  /* Allocate memory and assign output pointer */
  plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
  y = mxGetPr(plhs[0]);

  /* Soft thresholding */
  for( k = 0 ; k < n ; ++k ){
    if(x[k] > t)
      y[k] = x[k] - t;
    else if(x[k] < -t)
      y[k] = x[k] + t;
    else
      y[k] = 0.0;
  }
}
