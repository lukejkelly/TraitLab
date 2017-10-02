/*
 * fastBi2De.c - convert binary arrays to decimal vectors.
 * Adapted from a similar function by written by Giovanni Motta.
 *
 * The calling syntax is:
 *
 *		P_d = fastBi2De(P_b, 'left-msb'),
 *
 * where
 *	- P_b is a double array whose rows correspond to binary vectors with
 *      the most significant bit on the left-hand side ('left-msb' in
 *      bi2de function). We note that P_b has M rows and L columns when
 *      passed to C. However, in C it is stacked by columns as a vector of
 *      length M * L, hence the awkward indexing used below.
 *  - P_d is a double vector whose entries are the decimal representations
 *      of the rows in P_b.
 *  - 'left-msb' is ignored.
 *
 * Uncomment the corresponding lines to include the parallelised version.
 *
 * This is a MEX-file for MATLAB.
 */

// Header files.
#include "mex.h"
// #include <omp.h>

/*
 * Computational routines.
 */

// // Parallelised binary-to-decimal function.
// void fastBi2De_par(double *P_b, int M, int L, double *P_d)
// {
//   // Dynamic threading.
//   omp_set_dynamic(1);
//
//   // Allocate loop variables.
//   int i, j;
//
//   #pragma omp parallel for shared(P_b, M, L, P_d) private(i, j)
//   for(i = 0; i < M; i++)  // Loop over rows of P_b/P_d.
//   {
//     for(j = 0; j < L; j++) // Loop over columns to populate P_d(i).
//     {
//       // Use bit-shifting to convert binary to decimal.
//       P_d[i] += ((double) ( ((int) P_b[i + j * M]) << (L - 1 - j) ) );
//     }
//   }
// }

// Serial binary-to-decimal function.
void fastBi2De_ser(double *P_b, int M, int L, double *P_d)
{
  for(int i = 0; i < M; i++) // Loop over rows of P_b/P_d.
  {
    for(int j = 0; j < L; j++) // Loop over columns to populate P_d(i).
    {
      // Use bit-shifting to convert binary to decimal.
      P_d[i] += ((double) ( ((int) P_b[i + j * M]) << (L - 1 - j) ) );
    }
  }
}

/*
 * The gateway function.
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Checking inputs and outputs.

  // Checking correct number of input arguments.
  if(nrhs != 2)
  {
    mexErrMsgIdAndTxt("SDLT:Inputs:nrhs", "Two inputs required.");
  }

  // Checking correct number of output arguments.
  if(nlhs != 1)
  {
    mexErrMsgIdAndTxt("SDLT:Outputs:nlhs", "One output required.");
  }

  // P_b is an M x L binary array, where M is often 2^L - 1.
	if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
	{
		mexErrMsgIdAndTxt("SDLT:Inputs:notDouble", "P_b must be type double.");
	}

  // Declaring input variable.
  double *P_b;  // Binary array.

 	// Declaring output variable.
  double *P_d;  // Integer vector stored as double.

  // Declaring associated variables.
  int M;        // Scalar number of rows in P_b.
  int L;        // Scalar number of columns in P_b.

  // Getting pointers to input data.
  P_b = mxGetPr(prhs[0]);

  // Dimensions of P_b.
  M = mxGetM(prhs[0]);
  L = mxGetN(prhs[0]);

  // Preparing output data.
  plhs[0] = mxCreateDoubleMatrix(M, 1, mxREAL);

  // Assigning a pointer to output vector.
  P_d = mxGetPr(plhs[0]);

  // Calling the computational routine.
  // if(L <= 30)
  // {
    fastBi2De_ser(P_b, M, L, P_d);
  // }
  // else
  // {
  //   fastBi2De_par(P_b, M, L, P_d);
  // }
}
