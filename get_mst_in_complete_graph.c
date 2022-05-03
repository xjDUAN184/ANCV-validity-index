#include "mex.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define INF DBL_MAX

void prim(double *dist, double *ret, double *out_dist, int N) {

    double *lowcost = (double *)malloc(sizeof(double) * N);
    if (lowcost == NULL) {
        mexErrMsgTxt("Do not have enough memory.");
    }
    int *start = (int *)malloc(sizeof(int) * N);
    if (start == NULL) {
        mexErrMsgTxt("Do not have enough memory.");
    }

    int i, j, minid;
    double min;

    for (i = 0; i < N; i++) {
        lowcost[i] = dist[i];
        start[i] = 0;
    }

    lowcost[0] = -1;
    start[0] = -1;
    int p = 0;

    for (i = 1; i < N; i++) {
        min = INF;
        minid = -1;
        for (j = 1; j < N; j++) {
            if (lowcost[j] < min && lowcost[j] != -1) {
                min = lowcost[j];
                minid = j;
            }
        }
        ret[p] = start[minid] + 1;
        ret[N + p - 1] = minid + 1;
        out_dist[p] = min;
        p++;

        lowcost[minid] = -1;

        for (j = 1; j < N; j++) {
            if (minid != j && dist[minid * N + j] < lowcost[j]) {
                lowcost[j] = dist[minid * N + j];
                start[j] = minid;
            }
        }
    }

    if (lowcost != NULL)
        free(lowcost);
    if (start != NULL)
        free(start);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *ret, *output, *input;

    if (nrhs > 1) {
        mexErrMsgTxt("Too few input.");
    }

    input = mxGetPr(prhs[0]);
    int N = mxGetN(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(N - 1, 2, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(N - 1, 1, mxREAL);
    ret = mxGetPr(plhs[0]);
    output = mxGetPr(plhs[1]);

    prim(input, ret, output, N);
}
