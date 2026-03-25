#include <stddef.h>

void spreadd(const double *r, int m, double *spread) {
    int i;
    int j;

    if (r == NULL || m <= 0 || spread == NULL) {
        return;
    }

    *spread = 0.0f;
    for (i = 0; i < m; ++i) {
        for (j = 0; j < m; ++j) {
            double val = r[i * m + j];
            if (i == j) {
                double diff = val - 1.0f;
                *spread += diff * diff;
            } else {
                *spread += val * val;
            }
        }
    }
}
