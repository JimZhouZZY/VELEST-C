#include <stddef.h>

void spreadb(const double *r, int m, double *spread) {
    int i;
    int j;

    if (r == NULL || m <= 0 || spread == NULL) {
        return;
    }

    *spread = 0.0;
    for (i = 0; i < m; ++i) {
        for (j = 0; j < m; ++j) {
            double diff = (double)(i - j);
            double val = r[i * m + j];
            *spread += diff * diff * val * val;
        }
    }
}
