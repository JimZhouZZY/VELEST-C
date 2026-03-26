#include <math.h>

void refpath(double xe, double ye, double ze,
             double xr, double yr, double zr,
             double delta, int nl,
             double v[], double vsq[], double thk[],
             double tkj, int jl, int kk, double didjkk,
             double rp[3][200], int *nrp, double *refrerr) {
    (void)nl;

    double thkjl = thk[jl];
    double d1 = (xr - xe) / delta;
    double d2 = (yr - ye) / delta;

    int nrpd = kk - jl + 1;
    *nrp = nrpd + kk;

    thk[jl] = thkjl - tkj;

    rp[0][0] = xe;
    rp[1][0] = ye;
    rp[2][0] = ze;

    for (int i = 1; i < nrpd; ++i) {
        int m = jl + i - 1;
        double tng = v[m] / sqrtf(vsq[kk] - vsq[m]);
        rp[0][i] = rp[0][i - 1] + thk[m] * tng * d1;
        rp[1][i] = rp[1][i - 1] + thk[m] * tng * d2;
        rp[2][i] = rp[2][i - 1] + thk[m];
    }

    rp[0][nrpd] = rp[0][nrpd - 1] + (delta - didjkk) * d1;
    rp[1][nrpd] = rp[1][nrpd - 1] + (delta - didjkk) * d2;
    rp[2][nrpd] = rp[2][nrpd - 1];

    thk[jl] = thkjl;

    int nrp1 = *nrp - 1;
    int nrpd2 = nrpd + 1;

    for (int i = nrpd2; i < *nrp; ++i) {
        int m = (kk - 1) - (i - nrpd2);
        double tng = v[m] / sqrtf(vsq[kk] - vsq[m]);
        rp[0][i] = rp[0][i - 1] + thk[m] * tng * d1;
        rp[1][i] = rp[1][i - 1] + thk[m] * tng * d2;
        rp[2][i] = rp[2][i - 1] - thk[m];
    }

    int last = *nrp - 1;
    double dx = rp[0][last] - xr;
    double dy = rp[1][last] - yr;
    double dz = rp[2][last] - zr;
    *refrerr = sqrtf(dx * dx + dy * dy + dz * dz) * 1000.0;

    rp[0][last] = xr;
    rp[1][last] = yr;
    rp[2][last] = zr;
}
