/*
 * Givens rotation - compute rotation parameters
 * IMSL routine LSVG1, nucleus for LSV decomposition
 * from Fortran LSVG1 subroutine
 */

#include <stddef.h>
#include <math.h>

void lsvg1(double a, double b, double *dcos, double *dsin, double *sig) {
    double aa, bb;
    double doubleaa, doubleb, doubler;
    
    if (fabs(a) > fabs(b)) {
        /* Use A as main component; double precision for numerical stability */
        aa = fabs(a + a);
        doubleb = (double)b;
        doubleaa = (double)aa;
        doubler = sqrt(0.25 + (doubleb / doubleaa) * (doubleb / doubleaa));
        *sig = aa * (double)doubler;
        *dcos = a / (*sig);
        *dsin = b / (*sig);
    } else {
        /* Use B as main component or handle zero case */
        if (b == 0.0f) {
            /* Both A and B are effectively zero */
            *sig = 0.0f;
            *dcos = 0.0f;
            *dsin = 1.0f;
        } else {
            /* B is main component */
            bb = fabs(b + b);
            *sig = bb * sqrtf(0.25f + (a / bb) * (a / bb));
            *dcos = a / (*sig);
            *dsin = b / (*sig);
        }
    }
}
