#include <math.h>

extern void tiddid(int jl, int nl, const double v[], const double vsq[], const double thk[],
                   double tid[], double did[]);

void refract(int nl, const double v[], const double vsq[], const double thk[],
             int jl, double tkj, double delta,
             int *kk, double *tref, double *didjkk, double *xovmax) {
    double tid[100] = {0.0};
    double did[100] = {0.0};
    double tinj[100] = {0.0};
    double didj[100] = {0.0};
    double tr[100] = {0.0};

    tiddid(jl, nl, v, vsq, thk, tid, did);

    *tref = 100000.0;
    int j1 = jl + 1;

    for (int m = j1; m < nl; ++m) {
        if (tid[m] == 100000.0) {
            tr[m] = 100000.0;
        } else {
            double sqt = sqrtf(vsq[m] - vsq[jl]);
            tinj[m] = tid[m] - tkj * sqt / (v[m] * v[jl]);
            didj[m] = did[m] - tkj * v[jl] / sqt;
            tr[m] = tinj[m] + delta / v[m];
            if (didj[m] > delta) {
                tr[m] = 100000.0;
            }
        }

        if (tr[m] < *tref) {
            *tref = tr[m];
            *kk = m;
            *didjkk = didj[*kk];
        }
    }

    if (*tref == 100000.0) {
        *didjkk = 100000.0;
        *xovmax = 100000.0;
        *kk = -1;
        return;
    }

    int m = jl + 1;
    while (m < nl && tid[m] == 100000.0) {
        m++;
    }
    int lx = m;

    if (jl == 0) {
        *xovmax = tinj[lx] * v[lx] * v[0] / (v[lx] - v[0]);
        return;
    }

    m = jl;
    while (1) {
        tid[m] = 0.0;
        for (int l = 0; l < m; ++l) {
            if (vsq[m] <= vsq[l]) {
                tid[m] = 100000.0;
                break;
            }
            double sqt = sqrtf(vsq[m] - vsq[l]);
            double tim = thk[l] * sqt / (v[l] * v[m]);
            tid[m] += tim;
        }

        m--;
        if (tid[m + 1] < 100000.0 || m == 0) {
            break;
        }
    }

    if (tid[m + 1] < 100000.0) {
        int jx = m + 1;
        *xovmax = (tinj[lx] - tid[jx]) * v[lx] * v[jx] / (v[lx] - v[jx]);
    } else {
        *xovmax = tinj[lx] * v[lx] * v[0] / (v[lx] - v[0]);
    }
}
