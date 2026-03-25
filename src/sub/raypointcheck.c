#include <stdio.h>
#include <stddef.h>
#include <string.h>

typedef void (*chtop_fn_t)(double xx, double yy, double *zk,
                          const char *topo1file, const char *topo2file);

void raypointcheck(const double *rp, int nrp, int inrpmax, const char *staname,
                   int isingle, chtop_fn_t chtop_fn, const char *topo1file,
                   const char *topo2file) {
    int j;
    double zzz, dzzz;

    if (rp == NULL || nrp <= 0 || inrpmax <= 0 || staname == NULL ||
        chtop_fn == NULL) {
        return;
    }

    for (j = 0; j < nrp; ++j) {
        if (rp[j * 3 + 2] < 0.0f) {
            chtop_fn(-rp[j * 3 + 0], rp[j * 3 + 1], &zzz, topo1file,
                     topo2file);
            dzzz = rp[j * 3 + 2] - zzz;
            if (dzzz < 0.0f) {
                fprintf(stderr,
                        " ray in the air... ! rp3=%6.3f ZZ=%6.3f dz=%6.3f rp# "
                        "=%2d nrp=%2d STN=%.4s i %4d\n",
                        rp[j * 3 + 2], zzz, dzzz, j + 1, nrp, staname, isingle);
            }
        }
    }
}
