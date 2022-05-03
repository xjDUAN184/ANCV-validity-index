#include "mex.h"

void full(double *in, double *out, int len, int N) {
    int i;
    for (i = 0; i < len; i++) {
        int x = in[i] - 1;
        int y = in[i + len] - 1;
        double d = in[i + 2 * len];
        out[x * N + y] = d;
        out[y * N + x] = d;
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *output, *input;

    input = mxGetPr(prhs[0]);
    int len = mxGetM(prhs[0]);
    int N = input[2 * len - 1];
    plhs[0] = mxCreateDoubleMatrix(N, N, mxREAL);
    output = mxGetPr(plhs[0]);

    full(input, output, len, N);
}