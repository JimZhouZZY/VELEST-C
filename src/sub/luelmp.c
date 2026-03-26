/*
 * LU elimination solver (back-substitution phase)
 * Solves Ax=b where A is lower-triangular (from LUDECP)
 * Diagonal elements stored in reciprocal form
 * from Fortran LUELMP subroutine
 */

#include <stddef.h>

void luelmp(double a[], double b[], int n, double x[]) {
    const double zero = 0.0;
    double t;
    int ip = 0;
    int iw = 0;

    if (a == NULL || b == NULL || x == NULL || n <= 0) {
        return;
    }

    for (int i = 1; i <= n; ++i) {
        t = b[i - 1];
        int im1 = i - 1;
        if (iw != 0) {
            ip = ip + iw - 1;
            for (int k = iw; k <= im1; ++k) {
                t = t - a[ip] * x[k - 1];
                ip += 1;
            }
        } else {
            if (t != zero) {
                iw = i;
            }
            ip = ip + im1;
        }
        x[i - 1] = t * a[ip];
        ip += 1;
    }

    int n1 = n + 1;
    for (int i = 1; i <= n; ++i) {
        int ii = n1 - i;
        ip -= 1;
        int is = ip;
        int iq = ii + 1;
        t = x[ii - 1];
        if (n >= iq) {
            int kk = n;
            for (int k = iq; k <= n; ++k) {
                t = t - a[is] * x[kk - 1];
                kk -= 1;
                is = is - kk;
            }
        }
        x[ii - 1] = t * a[is];
    }
}
