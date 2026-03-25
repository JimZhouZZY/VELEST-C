#include <stddef.h>

void matrtran(const double *a, int n, int m, double *at) {
    int i;
    int j;

    if (a == NULL || at == NULL || n <= 0 || m <= 0) {
        return;
    }

    for (i = 0; i < n; ++i) {
        for (j = 0; j < m; ++j) {
            at[j * n + i] = a[i * m + j];
        }
    }
}
