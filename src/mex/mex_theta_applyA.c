/*
 * ===========================================================================
 *
 *       Filename:  mex_theta_applyA.c
 *
 *    Description:  MEX interface for computing A(R'R) or A(U'R+R'U) for theta
 *
 *        Version:  1.0
 *        Created:  05/19/2020 15:28:57 PM
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

double cblas_ddot(long, double *, long, double *, long);
double ddot_(long*, double *, long *, double *, long*);

/* RHS = mex_theta_applyA(E, R)
 * RHS = mex_theta_applyA(E, R, U)
 * */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    mwSize n, p, num_E;
    n = mxGetN(prhs[1]);
    p = mxGetM(prhs[1]);
    num_E = mxGetM(prhs[0]);

    double *pE = mxGetPr(prhs[0]);
    double *pR = mxGetPr(prhs[1]);
    double *pU = NULL;
    if (nrhs > 2)
        pU = mxGetPr(prhs[2]);

#ifdef _DEBUG
    printf("(n, p, num_E) = (%ld, %ld, %ld)\n", n, p, num_E, sigma, tol);
#endif

    /* create output vector */
    plhs[0] = mxCreateDoubleMatrix(num_E, 1, mxREAL);
    double *g2 = mxGetPr(plhs[0]);

    /* main loop */
    if (pU == NULL){
        for (mwIndex edge_i = 0; edge_i < num_E; ++edge_i)
        {
            mwIndex i = (mwIndex)(pE[edge_i]) - 1;
            mwIndex j = (mwIndex)(pE[edge_i + num_E]) - 1;

            double g2i = cblas_ddot(p, pR + i * p, 1, pR + j * p, 1);
            g2[edge_i] = 2 * g2i;

        }
    } else {
        for (mwIndex edge_i = 0; edge_i < num_E; ++edge_i){
            mwIndex i = (mwIndex)(pE[edge_i]) - 1;
            mwIndex j = (mwIndex)(pE[edge_i + num_E]) - 1;

            double g2i = cblas_ddot(p, pR + i * p, 1, pU + j * p, 1);
            g2i += cblas_ddot(p, pU + i * p, 1, pR + j * p, 1);
            g2[edge_i] = 2 * g2i;
        }
    }
}

double cblas_ddot(long n, double *x, long incx, double *y, long incy){
    return ddot_(&n, x, &incx, y, &incy);
}
