#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "mex.h"

double cblas_ddot(long, double *, long, double *, long);
double ddot_(long*, double *, long *, double *, long*);

/**
 * LHS = mex_maxcut_applyA(A, R)    -> b = A(R'R)
 * LHS = mex_maxcut_applyA(A, R, U) -> b = A(U'R + R'U)
 *
 * A is 1 by m cell array
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    mwSize n, p, m;

    if (nrhs != 2 && nrhs != 3){
        mexErrMsgTxt("Invalid input arguments.");
        return;
    }

    if (nlhs > 1){
        mexErrMsgTxt("Too many output arguments.");
        return;
    }

    if (!mxIsCell(prhs[0])){
        mexErrMsgTxt("A is not cell array.");
        return;
    }
    n = mxGetN(prhs[1]);
    p = mxGetM(prhs[1]);
    m = mxGetN(prhs[0]);

    double *pR = mxGetPr(prhs[1]);
    double *pU = NULL;
    if (nrhs > 2)
        pU = mxGetPr(prhs[2]);

#ifdef _DEBUG
    printf("(n, p, m) = (%ld, %ld, %ld)\n", n, p, m);
#endif

    /* create output vector */
    plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
    double *g2 = mxGetPr(plhs[0]);

    /* main loop */
    if (pU == NULL){
        for (mwIndex ic = 0; ic < m; ++ic)
        {
            /* A{ic} is 3 x 3 matrix */
            mxArray *Ai = mxGetCell(prhs[0], ic);
            mwSize mAi = mxGetM(Ai);
            mwSize nAi = mxGetN(Ai);
            double *pAi = mxGetPr(Ai);

            if (mAi != 3 || nAi != 3){
                mexErrMsgTxt("Each element of A should be 3 x 3.");
                return;
            }

            for (mwIndex i = 0; i < 3; ++i){
                mwIndex ii = (mwIndex)(pAi[3 * i + 0]) - 1;
                mwIndex jj = (mwIndex)(pAi[3 * i + 1]) - 1;
                double  vv = pAi[3 * i + 2];

                double g2i = cblas_ddot(p, pR + ii * p, 1, pR + jj * p, 1);
                g2[ic] += vv * g2i;
            }

        }
    } else {
        for (mwIndex ic = 0; ic < m; ++ic)
        {
            /* A{ic} is 3 x 3 matrix */
            mxArray *Ai = mxGetCell(prhs[0], ic);
            mwSize mAi = mxGetM(Ai);
            mwSize nAi = mxGetN(Ai);
            double *pAi = mxGetPr(Ai);

            if (mAi != 3 || nAi != 3){
                mexErrMsgTxt("Each element of A should be 3 x 3.");
                return;
            }

            for (mwIndex i = 0; i < 3; ++i){
                mwIndex ii = (mwIndex)(pAi[3 * i + 0]) - 1;
                mwIndex jj = (mwIndex)(pAi[3 * i + 1]) - 1;
                double  vv = pAi[3 * i + 2];

                double g2i = cblas_ddot(p, pR + ii * p, 1, pU + jj * p, 1);
                g2i += cblas_ddot(p, pU + ii * p, 1, pR + jj * p, 1);
                g2[ic] += vv * g2i;
            }

        }
    }
}

double cblas_ddot(long n, double *x, long incx, double *y, long incy){
    return ddot_(&n, x, &incx, y, &incy);
}
