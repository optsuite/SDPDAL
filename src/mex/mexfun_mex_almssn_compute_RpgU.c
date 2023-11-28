/*
 * ==========================================================================
 *
 *       Filename:  mex_almssn_compute_RpgU.c
 *
 *    Description:  compute RpgU = R * pgU
 *
 *        Version:  1.0
 *        Created:  04/14/2021 07:02:05 PM
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

/** RpgU = mex_almssn_compute_RpgU(R, pgU)
 * 
 *  compute RpgU = R * pgU with a single call of DSYMM
 *
 *  Compiled with: mex -O mex_almssn_compute_RpgU.c -lmwblas
 *
 *  Note: on entry, only lower triangular part of pgU is stored/referenced
 */


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    if (nrhs != 2){
        mexErrMsgTxt("Argument error: two inputs needed.");
        return;
    }

    if (nlhs != 1){
        mexErrMsgTxt("Argument error: one output needed.");
    }

    mwSignedIndex n, k, n1, n2;
    double *pgU, *R;

    /* dimension check */
    k  = mxGetM(prhs[0]);
    n  = mxGetN(prhs[0]);
    n1 = mxGetM(prhs[1]);
    n2 = mxGetN(prhs[1]);

    if (n != n1){
        mexErrMsgTxt("Argument error: matrix dimension must agree.");
        return;
    }

    /* get R = prhs[0], pgU = prhs[1] */
    R   = mxGetPr(prhs[0]);
    pgU = mxGetPr(prhs[1]);

    /* prepare outputs */
    plhs[0] = mxCreateDoubleMatrix(k, n2, mxREAL);
    double *RpgU = mxGetPr(plhs[0]);

    /* prepare inputs */
    char uplo = 'U', side = 'R';
    double one = 1, zero  = 0;

    /* call dsymm */
    dsymm(&side, &uplo, &k, &n2, &one, pgU, &n, R, &k, &zero, RpgU, &k);

}
