#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double euclid_dist(double *x, double *y, int features) {
    double sum = 0.0;
    int i;
    for (i = 0; i < features; i++) {
        sum += (x[i] - y[i]) * (x[i] - y[i]);
    }
    return sqrt(sum);
}

double *get_euclid_dist(double *data, double *dist, int N, int features) {
    int len = N * (N - 1) / 2;
    int i, j, k = 0;
    for (i = 0; i < N - 1; i++) {
        for (j = i + 1; j < N; j++) {
            dist[k] = i + 1;
            dist[k + len] = j + 1;
            dist[k + 2 * len] =
                euclid_dist(&data[i * features], &data[j * features], features);
            k++;
        }
    }
    return dist;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *output, *input;

    input = mxGetPr(prhs[0]);
    int N = mxGetN(prhs[0]);
    int features = mxGetM(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(N * (N - 1) / 2, 3, mxREAL);
    output = mxGetPr(plhs[0]);

    get_euclid_dist(input, output, N, features);
}
