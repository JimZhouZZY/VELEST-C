#include <math.h>

void tiddid(int jl, int nl, const double v[], const double vsq[], const double thk[],
            double tid[], double did[]) {
    int j1 = jl + 1;
    for (int m = j1; m < nl; ++m) {
        tid[m] = 0.0f;
        did[m] = 0.0f;

        double tid1 = 0.0f;
        double tid2 = 0.0f;
        double did1 = 0.0f;
        double did2 = 0.0f;

        for (int l = 0; l < m; ++l) {
            if (vsq[m] <= vsq[l]) {
                tid[m] = 100000.0f;
                did[m] = 100000.0f;
                break;
            }

            double sqt = sqrtf(vsq[m] - vsq[l]);
            double tim = thk[l] * sqt / (v[l] * v[m]);
            double dimm = thk[l] * v[l] / sqt;

            if (l < jl) {
                tid1 += tim;
                did1 += dimm;
            } else {
                tid2 += tim;
                did2 += dimm;
            }
        }

        if (tid[m] != 100000.0f) {
            tid[m] = tid1 + 2.0f * tid2;
            did[m] = did1 + 2.0f * did2;
        }
    }
}
