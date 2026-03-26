#include <math.h>
#include <stddef.h>

void ludecp(const double *a, double *ul, int n, double *d1, double *d2,
            int *ier) {
    const double one = 1.0f;
    const double four = 4.0f;
    const double sixtn = 16.0f;
    const double sixth = 0.0625f;
    double rn, x;
    int ip, i, iq, ir, j, k, ip1;

    if (a == NULL || ul == NULL || d1 == NULL || d2 == NULL || ier == NULL ||
        n <= 0) {
        *ier = -1;
        return;
    }

    *d1 = one;
    *d2 = 0.0;
    rn = one / ((double)n * sixtn);
    ip = 0;
    *ier = 0;

    for (i = 0; i < n; ++i) {
        iq = ip;
        ir = 0;
        for (j = 0; j <= i; ++j) {
            x = a[ip];

            if (j != 0) {
                for (k = iq; k <= ip1; ++k) {
                    x -= ul[k] * ul[ir];
                    ir += 1;
                }
            }

            if (i != j) {
                ul[ip] = x * ul[ir];
            } else {
                *d1 *= x;
                if (a[ip] + x * rn <= a[ip]) {
                    *ier = 129;
                    return;
                }

                while (fabs(*d1) > one) {
                    *d1 *= sixth;
                    *d2 += four;
                }
                while (fabs(*d1) < sixth) {
                    *d1 *= sixtn;
                    *d2 -= four;
                }

                ul[ip] = one / sqrtf(x);
            }

            ip1 = ip;
            ip += 1;
            ir += 1;
        }
    }
}
