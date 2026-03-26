#include <math.h>
#include <stddef.h>

void setunt(int nitt, int invertratio, int nsinv, int *icount,
            double xythet, double stathet, double othet, double vthet,
            double zthet, double *scale) {
    if (scale == NULL || icount == NULL) {
        return;
    }

    *icount = nitt % invertratio;

    scale[0] = 1.0f;
    scale[6] = 1.0f;
    scale[1] = sqrtf(othet / xythet);
    scale[2] = scale[1];
    scale[3] = sqrtf(othet / zthet);
    scale[4] = 0.0;

    if (nsinv != 0 && *icount == 0) {
        scale[4] = sqrtf(othet / stathet);
    }

    scale[5] = 0.0;
    if (*icount == 0) {
        scale[5] = sqrtf(othet / vthet);
    }
}
