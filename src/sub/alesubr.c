#include <math.h>

int alesubr(const double *sv, int m, double *ale_out) {
    double alesum;
    int i;
    int izero;

    if (sv == 0 || ale_out == 0 || m <= 0) {
        return -1;
    }

    alesum = 0.0;
    izero = 0;

    for (i = 0; i < m; ++i) {
        double ratio = sv[i] / sv[0];
        if (ratio <= 0.0) {
            izero += 1;
        } else {
            alesum += log10f(ratio);
        }
    }

    if (m - izero <= 0) {
        *ale_out = 10.0 * (double)izero;
    } else {
        *ale_out = -(alesum / (double)(m - izero)) + 10.0 * (double)izero;
    }

    return 0;
}