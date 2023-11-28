/*
 * ===========================================================================
 *
 *       Filename:  mex_theta_applyA_full.c
 *
 *    Description:  MEX interface for computing A(R'R) or A(U'R+R'U) for theta
 *                  use full matrix representation
 *
 *        Version:  1.0
 *        Created:  06/26/2021 14:46:57 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Haoyang Liu (@liuhy), liuhaoyang@pku.edu.cn
 *   Organization:  BICMR, Peking University
 *      Copyright:  Copyright (c) 2021, Haoyang Liu
 *
 * ===========================================================================
 */
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "mex.h"

/* RHS = mex_theta_applyA_full(E, X)
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    mwSize n, p, num_E;
    n = mxGetN(prhs[1]);
    num_E = mxGetM(prhs[0]);

    double *pE = mxGetPr(prhs[0]);
    double *pX = mxGetPr(prhs[1]);

    /* create output vector */
    plhs[0] = mxCreateDoubleMatrix(num_E, 1, mxREAL);
    double *g2 = mxGetPr(plhs[0]);

    /* main loop */
    for (mwIndex edge_i = 0; edge_i < num_E; ++edge_i)
    {
        mwIndex i = (mwIndex)(pE[edge_i]) - 1;
        mwIndex j = (mwIndex)(pE[edge_i + num_E]) - 1;

        g2[edge_i] = 2 * pX[j * n + i];

    }
}

