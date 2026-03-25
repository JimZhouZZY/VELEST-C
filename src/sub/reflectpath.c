#include <math.h>
#include <stdlib.h>

void reflectpath(double xe, double ye, double ze, double xr, double yr, double zr,
                 double delta, int nl, const double v[], const double vsq[], const double thk[],
                 int jl, double tkj, double ain, int mll,
                 double rp[3][200], int *nrp, double *reflerr) {
    (void)nl;

    double d1 = (xr - xe) / delta;
    double d2 = (yr - ye) / delta;

    *nrp = jl + 1 + (mll - jl) * 2 + 2;

    double salpha = ain;

    rp[0][0] = xe;
    rp[1][0] = ye;
    rp[2][0] = ze;

    double tng = tanf(asinf(salpha));
    double deljl = tng * (thk[jl] - tkj);

    rp[0][1] = rp[0][0] + deljl * d1;
    rp[1][1] = rp[1][0] + deljl * d2;
    rp[2][1] = rp[2][0] + (thk[jl] - tkj);

    int m = jl;
    int ld = mll - jl;
    int j = 1;

    for (int i = 0; i < ld; ++i) {
        ++j;
        m = jl + i + 1;
        tng = v[m] * salpha / sqrtf(vsq[jl] - vsq[m] * salpha * salpha);
        rp[0][j] = rp[0][j - 1] + thk[m] * tng * d1;
        rp[1][j] = rp[1][j - 1] + thk[m] * tng * d2;
        rp[2][j] = rp[2][j - 1] + thk[m];
    }

    if (m != mll) {
        abort();
    }

    ++j;
    deljl = tng * thk[m];
    rp[0][j] = rp[0][j - 1] + deljl * d1;
    rp[1][j] = rp[1][j - 1] + deljl * d2;
    rp[2][j] = rp[2][j - 1] - thk[m];

    for (int i = j + 1; i <= (*nrp - 1); ++i) {
        --m;
        tng = v[m] * salpha / sqrtf(vsq[jl] - vsq[m] * salpha * salpha);
        rp[0][i] = rp[0][i - 1] + thk[m] * tng * d1;
        rp[1][i] = rp[1][i - 1] + thk[m] * tng * d2;
        rp[2][i] = rp[2][i - 1] - thk[m];
    }

    int last = *nrp - 1;
    double dx = rp[0][last] - xr;
    double dy = rp[1][last] - yr;
    double dz = rp[2][last] - zr;
    *reflerr = sqrtf(dx * dx + dy * dy + dz * dz) * 1000.0f;

    rp[0][last] = xr;
    rp[1][last] = yr;
    rp[2][last] = zr;
}
