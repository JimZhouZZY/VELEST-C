#include <stddef.h>

void strpath(double xe, double ye, double ze, double xr, double yr, double zr,
             double *rp, int *nrp) {
    if (rp == NULL || nrp == NULL) {
        return;
    }

    *nrp = 2;
    rp[0] = xe;
    rp[1] = ye;
    rp[2] = ze;
    rp[3] = xr;
    rp[4] = yr;
    rp[5] = zr;
}
