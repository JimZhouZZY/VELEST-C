#include <stddef.h>

void matrmult(const double *a, int m, int p, const double *b, int p1, int n,
              double *c, int m1, int n1) {
    int i;
    int j;
    int k;

    if (a == NULL || b == NULL || c == NULL || m <= 0 || p <= 0 || n <= 0 ||
        p1 <= 0 || m1 <= 0 || n1 <= 0) {
        return;
    }

    if (m != m1 || p != p1 || n != n1) {
        return;
    }

    for (i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j) {
            double s = 0.0;
            for (k = 0; k < p; ++k) {
                s += a[i * p + k] * b[k * n + j];
            }
            c[i * n + j] = s;
        }
    }
}
