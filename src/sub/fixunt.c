#include <stddef.h>

void fixunt(double *b, int neqs, int nshot, int nl, int ksta,
            const double *scale, const double *vdamp, int itotmodels,
            int inltot, int nlfirst) {
    int i, j, k, l, m;

    if (b == NULL || scale == NULL || vdamp == NULL) {
        return;
    }

    i = -1;

    if (neqs > 0) {
        for (j = 0; j < neqs; ++j) {
            for (k = 0; k < 4; ++k) {
                ++i;
                b[i] *= scale[k];
            }
        }
    }

    if (nshot > 0) {
        for (j = 0; j < nshot; ++j) {
            ++i;
            b[i] *= scale[6];
        }
    }

    if (scale[5] == 0.0f) {
        return;
    }

    l = 0;
    m = 1;
    for (j = 0; j < nl; ++j) {
        ++i;
        ++l;
        if (l > nlfirst) {
            ++m;
            l = 1;
        }
        if (m < 1 || m > itotmodels || l < 1 || l > inltot) {
            return;
        }
        double vd = vdamp[(m - 1) * inltot + (l - 1)];
        if (vd != 0.0f) {
            b[i] *= scale[5] / vd;
        }
    }

    if (ksta > 0) {
        for (j = 0; j < ksta; ++j) {
            ++i;
            b[i] *= scale[4];
        }
    }
}
