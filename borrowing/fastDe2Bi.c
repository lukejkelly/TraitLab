/*
 * fastDe2Bi_c_par.c - convert decimal vectors to binary arrays.
 * Adapted from a similar function written by Giovanni Motta.
 *
 * The calling syntax is:
 *
 *		P_b = fastDe2Bi(P_d, L, 'left-msb'),
 *
 * where
 *	- P_b is a double array whose rows correspond to binary vectors with
 *      the most significant bit on the left-hand side ('left-msb' in
 *      de2bi function). We note that P_b has M rows and L columns when
 *      passed to C. However, in C it is stacked by columns as a vector of
 *      length M * L, hence the awkward indexing used below.
 *  - P_d is a double vector whose entries are the decimal representations
 *      of the rows in P_b.
 *  - L is the length of each binary vector / row of P_b.
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
// void fastDe2Bi_par(double *P_b, int M, double *L, double *P_d)
// {
//   // Dynamic threading.
//   omp_set_dynamic(1);
//
//   // Allocate loop variables.
//   int i, j;
//
//   #pragma omp parallel for shared(P_b, M, L, P_d) private(i, j)
//   for(i = 0; i < M; i++) // Loop over rows of P_b/P_d.
//   {
//     // Storing P_d[i] as an int so that we can bitshift it.
//     int k = (int) P_d[i];
//
//     for(j = 0; j < ((int) *L); j++) // Looping over columns of P_b(i, :).
//     {
//       // Use bit-shifting and bitwise AND to convert decimal to binary.
//       P_b[i + (((int) *L) - j - 1) * M] = ((double) ((k >> j) & 0x1));
//     }
//   }
// }

// Serial binary-to-decimal function.
void fastDe2Bi_ser(double *P_b, int M, double *L, double *P_d)
{
  for(int i = 0; i < M; i++) // Looping over rows of P_b/P_d.
  {
    // Storing P_d[i] as an int so that we can bitshift it.
    int k = (int) P_d[i];

    for(int j = 0; j < ((int) *L); j++) // Looping over columns of P_b(i, :).
    {
      // Use bit-shifting and bitwise AND to convert decimal to binary.
      P_b[i + (((int) *L) - j - 1) * M] = ((double) ((k >> j) & 0x1));
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
  if(nrhs != 3)
  {
    mexErrMsgIdAndTxt("SDLT:Inputs:nrhs",
      "Three inputs required: P_d, L and 'left-msb'.");
  }

  // Checking correct number of output arguments.
  if(nlhs != 1)
  {
    mexErrMsgIdAndTxt("SDLT:Outputs:nlhs", "One output required: P_b.");
  }

  // P_d is an M x 1 decimal array, where M is often 2^L - 1.
	if(!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
	{
		mexErrMsgIdAndTxt("SDLT:Inputs:notDouble", "P_d must be type double.");
	}

	if(!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
    || mxGetNumberOfElements(prhs[1]) != 1)
	{
		mexErrMsgIdAndTxt("SDLT:Inputs:notScalar",
      "L must be a scalar of type double.");
	}

  // Declaring input variables.
  double *P_d;  // Vector of integers.
  double *L;    // Scalar.

 	// Declaring output variable.
  double *P_b;  // Binary array stored as double.

  // Declaring associated variable.
  int m, n, M;  // Scalar number of rows in P_b.

  // Getting pointers to input data.
  P_d = mxGetPr(prhs[0]);
  L   = mxGetPr(prhs[1]);

  // Number of rows in P_b.
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  M = (m > n) ? m : n;

  // Preparing output data.
  plhs[0] = mxCreateDoubleMatrix(M, *L, mxREAL);

  // Assigning a pointer to output array.
  P_b = mxGetPr(plhs[0]);

  // Calling the computational routine.
  // if(*L <= 12)
  // {
    fastDe2Bi_ser(P_b, M, L, P_d);
  // }
  // else
  // {
  //   fastDe2Bi_par(P_b, M, L, P_d);
  // }
}
