/*
function [row,column,value] = mc_B(T,m,n);
        input: triangle inequalities T (m,4) as a long vector
               m ... number of inequalities
                  
        output: [row,column,value] triple specifying the sparse matrix B
                for the inequality constraints B(X) <= e
 
       last modified: july 2019
*/

// TODO: you need to add diagonal elem

#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{
  
//   /* call with 2 inputs and 3 outputs */
//   if (nlhs != 3 || nrhs != 2) {
//       mexErrMsgTxt("input-output arguments inconsistent.");
//   }
  
  // input
  double *T = mxGetPr(prhs[0]);
  long m = (long)mxGetScalar(prhs[1]);
  long n = (long)mxGetScalar(prhs[2]);
  
  // output -> matrix in [sparse format row,column,value ]
  plhs[0] = mxCreateDoubleMatrix(6*m, 1, mxREAL); // row
  plhs[1] = mxCreateDoubleMatrix(6*m, 1, mxREAL); // column
  plhs[2] = mxCreateDoubleMatrix(6*m, 1, mxREAL); // value
  
  double *row = mxGetPr(plhs[0]);
  double *column = mxGetPr(plhs[1]);
  double *value = mxGetPr(plhs[2]);
  
  int plane, l;
  int i, j, k, type;
  
  for (plane = 0; plane < m; ++plane){
      i = (int)T[plane];
      j = (int)T[m + plane];
      k = (int)T[2*m + plane];
      type = (int)T[3*m + plane];
      
      
      
      // rows
      for (l = 6*plane; l < 6*(plane+1); ++l){
        row[l] = plane+1;
      }    
        
      // column
      column[6*plane] = (j-1)*n+i;
      column[6*plane+1] = (i-1)*n+j;
      column[6*plane+2] = (j-1)*n+k;
      column[6*plane+3] = (k-1)*n+j;
      column[6*plane+4] = (k-1)*n+i;
      column[6*plane+5] = (i-1)*n+k;
      
      // check type of triagle inequality
      if (type == 1){
        for (l = 6*plane; l < 6*(plane+1); ++l){
            value[l] = -0.5;
        } 
      }
      else if (type == 2){
        value[6*plane] = -0.5;
        value[6*plane+1] = -0.5;
        value[6*plane+2] = 0.5;
        value[6*plane+3] = 0.5;
        value[6*plane+4] = 0.5;
        value[6*plane+5] = 0.5;
      }  
      else if (type == 3){
        value[6*plane] = 0.5;
        value[6*plane+1] = 0.5;
        value[6*plane+2] = 0.5;
        value[6*plane+3] = 0.5;
        value[6*plane+4] = -0.5;
        value[6*plane+5] = -0.5;
      }
      else if (type == 4){
        value[6*plane] = 0.5;
        value[6*plane+1] = 0.5;
        value[6*plane+2] = -0.5;
        value[6*plane+3] = -0.5;
        value[6*plane+4] = 0.5;
        value[6*plane+5] = 0.5;
      }
  } 
}