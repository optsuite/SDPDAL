/*
 * ===========================================================================
 *
 *       Filename:  mex_theta_applyAT.c
 *
 *    Description:  MEX interface for computing A'(y) for theta
 *
 *        Version:  1.0
 *        Created:  05/19/2020 16:06:00 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Haoyang Liu (@liuhy), liuhaoyang@pku.edu.cn
 *   Organization:  BICMR, Peking University
 *      Copyright:  Copyright (c) 2020, Haoyang Liu
 *
 * ===========================================================================
 */
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "mex.h"

/* LHS = mex_theta_applyAT(Ref_E, y)
 * */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    mwSize n, num_E;
    n = mxGetM(prhs[0]);
    num_E = mxGetM(prhs[1]);

    double *y = mxGetPr(prhs[1]);

    /* export pattern */
    mwIndex *ir_ref = mxGetIr(prhs[0]);
    mwIndex *jc_ref = mxGetJc(prhs[0]);

    /* create symmetric sparse mat */
    plhs[0] = mxCreateSparse(n, n, num_E, mxREAL);
    /* copy pattern */
    memcpy(mxGetIr(plhs[0]), ir_ref, num_E * sizeof(mwIndex));
    memcpy(mxGetJc(plhs[0]), jc_ref, (n + 1) * sizeof(mwIndex));

    double *g2 = mxGetPr(plhs[0]);

    /* copy value */
    memcpy(g2, y, num_E * sizeof(double));
}

