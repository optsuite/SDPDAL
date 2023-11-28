/*
 * ==========================================================================
 *
 *       Filename:  mex_cutting_plane_mc.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/31/2021 09:41:48 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Haoyang Liu (@liuhy), liuhaoyang@pku.edu.cn
 *   Organization:  BICMR, PKU
 *
 * ==========================================================================
 */

#include <stdlib.h>
#include <string.h>
#include "mex.h"

void insert(mwIndex, mwIndex, mwIndex, double, double, double, double, double*, double*,
        mwIndex *, double*, mwIndex);
void get_max_elem(const double*, mwIndex, double*, mwIndex*);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    if (nlhs != 2 || (nrhs != 2 && nrhs != 3)){
        mexErrMsgTxt("Usage: [pos, sgn] = mex_cutting_plane_mc(X, K[ ,eps]).");
        return;
    }

    mwSize n = mxGetM(prhs[0]);
    double *pK = mxGetPr(prhs[1]);
    double *pX = mxGetPr(prhs[0]);
    mwIndex K = (mwIndex)(*pK);
    double bnd = 1;

    if (nrhs == 3){
        double *eps = mxGetPr(prhs[2]);
        bnd = (1 + *eps);
    }

    double *bpos = (double*)malloc(4 * K * sizeof(double));
    double *bsgn = (double*)malloc(3 * K * sizeof(double));

    mwIndex ncuts = 0;
    double v, vmax = 0;

    for (mwIndex i = 0; i < n; ++i){
        for (mwIndex j = i+1; j < n; ++j){
            for (mwIndex k = j+1; k < n; ++k){
                double Xij = pX[i * n + j], Xik = pX[i * n + k], Xjk = pX[j * n + k];
                v = Xij + Xik + Xjk;
                if (ncuts < K ? v < -bnd : v < vmax)
                    insert(i, j, k, v, 1, 1, 1, bpos, bsgn, &ncuts, &vmax, K);

                v = Xij - Xik - Xjk;
                if (ncuts < K ? v < -bnd : v < vmax)
                    insert(i, j, k, v, 1, -1, -1, bpos, bsgn, &ncuts, &vmax, K);

                v = -Xij + Xik - Xjk;
                if (ncuts < K ? v < -bnd : v < vmax)
                    insert(i, j, k, v, -1, 1, -1, bpos, bsgn, &ncuts, &vmax, K);

                v = -Xij - Xik + Xjk;
                if (ncuts < K ? v < -bnd : v < vmax)
                    insert(i, j, k, v, -1, -1, 1, bpos, bsgn, &ncuts, &vmax, K);

            }
        }
    }

    plhs[0] = mxCreateDoubleMatrix(4, ncuts, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(3, ncuts, mxREAL);

    double *pos = mxGetPr(plhs[0]);
    double *sgn = mxGetPr(plhs[1]);

    memcpy(pos, bpos, 4 * ncuts * sizeof(double));
    memcpy(sgn, bsgn, 3 * ncuts * sizeof(double));

    free(bpos);
    free(bsgn);
}

void insert(mwIndex i, mwIndex j, mwIndex k, double v,
        double s1, double s2, double s3, double*pos, double*sgn,
        mwIndex *ncuts, double *vmax, mwIndex K){
    if (*ncuts < K){
        pos[*ncuts * 4] = i + 1;
        pos[*ncuts * 4 + 1] = j + 1;
        pos[*ncuts * 4 + 2] = k + 1;
        pos[*ncuts * 4 + 3] = v;

        sgn[*ncuts * 3 + 0] = s1;
        sgn[*ncuts * 3 + 1] = s2;
        sgn[*ncuts * 3 + 2] = s3;
        *ncuts += 1;
    } else {
        mwIndex imax;
        get_max_elem(pos, K, vmax, &imax);

        if ( v < *vmax){
            pos[imax * 4] = i + 1;
            pos[imax * 4 + 1] = j + 1;
            pos[imax * 4 + 2] = k + 1;
            pos[imax * 4 + 3] = v;

            sgn[imax * 3 + 0] = s1;
            sgn[imax * 3 + 1] = s2;
            sgn[imax * 3 + 2] = s3;
        }

    }
}

void get_max_elem(const double *pos, mwIndex K, double *vmax, mwIndex *imax){
    *imax = 0;
    *vmax = pos[3];
    for (mwIndex i = 1; i < K; ++i){
        if (pos[i * 4 + 3] > *vmax){
            *vmax = pos[i * 4 + 3];
            *imax = i;
        }
    }
}
