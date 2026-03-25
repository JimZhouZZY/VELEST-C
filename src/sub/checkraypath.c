#include <stddef.h>

#include <stddef.h>

void checkraypath(double *rp, int *nrp, int inrpmax) {
    if (rp == NULL || nrp == NULL || inrpmax <= 0) {
        return;
    }

    if (rp[0] == rp[3] && rp[1] == rp[4] && rp[2] == rp[5]) {
        rp[2] -= 0.0001f;
        return;
    }
}
