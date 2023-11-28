/*
 * ==========================================================================
 *
 *       Filename:  mex_almssn_compute_pgU.c
 *
 *    Description:  compute pgU = U'R + R'U
 *
 *        Version:  1.0
 *        Created:  04/14/2021 06:30:05 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Haoyang Liu (@liuhy), liuhaoyang@pku.edu.cn
 *   Organization:  BICMR, PKU
 *
 * ==========================================================================
 */

#include "mex.h"
#include "blas.h"

/** pgU = mex_almssn_compute_pgU(U, R)
 * 
 *  compute pgU = U'R + R'U with a single call of DSYR2K
 *
 *  Compiled with: mex -O mex_almssn_compute_pgU.c -lmwblas
 *
 *  Note: only lower triangular part of pgU is stored/referenced
 */


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    if (nrhs != 2){
        mexErrMsgTxt("Argument error: two inputs needed.");
        return;
    }

    if (nlhs != 1){
        mexErrMsgTxt("Argument error: one output needed.");
    }

    mwSignedIndex n, k, n1, k1;
    double *U, *R;

    /* dimension check */
    k  = mxGetM(prhs[0]);
    n  = mxGetN(prhs[0]);
    k1 = mxGetM(prhs[1]);
    n1 = mxGetN(prhs[1]);

    if (n != n1 || k != k1){
        mexErrMsgTxt("Argument error: matrix dimension must agree.");
        return;
    }

    /* get U = prhs[0], R = prhs[1] */
    U = mxGetPr(prhs[0]);
    R = mxGetPr(prhs[1]);

    /* prepare outputs */
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
    double *pgU = mxGetPr(plhs[0]);

    /* prepare inputs */
    char uplo = 'L', trans = 'T';
    double one = 1, zero  = 0;

    /* call dsyr2k */
    dsyr2k(&uplo, &trans, &n, &k, &one, U, &k, R, &k, &zero, pgU, &n);

}
