/*
 * Reject observation with large residual from further calculation
 * Re-normalize weights after rejection
 * from Fortran REJECTOBS subroutine
 */

#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>

/* External global variables from vel_com.inc */
extern int nitt;
extern int isingle;
extern int nobswithw0;
extern int nvar;
extern bool single_turbo;
extern double w[][100];
extern int knobs[];
extern int kpwt[][100];
extern double sphase[][100];
extern double swtfac;
extern double res[][100];

void rejectobs(int i, int nobs, int *iresflag) {
    double wsum, xkndw;
    int knobst, iii, iobswt;
    
    /* Check if weight already zero */
    if (w[nobs - 1][i] == 0.0) {
        return;
    }
    
    /* Check if sufficient observations remain */
    if ((knobs[i] - nobswithw0) == nvar) {
        return;
    }
    
    /* Re-normalize weights for this event */
    nobswithw0++;
    wsum = 0.0;
    knobst = knobs[i];
    
    for (iii = 0; iii < knobst; iii++) {
        iobswt = kpwt[iii][i];
        if (iobswt < 4 && w[iii][i] != 0.0) {
            /* Apply phase-specific weight scaling */
            if (sphase[iii][i] == 1.0f || sphase[iii][i] == 2.0f) {
                w[iii][i] = swtfac * 1.0f / (double)(1 << (iobswt * 2));  /* S-phase */
            } else {
                w[iii][i] = 1.0f / (double)(1 << (iobswt * 2));  /* P or M phase */
            }
        } else {
            w[iii][i] = 0.0;  /* Weight 4 or invalid */
        }
        
        if (iii == (nobs - 1)) {
            w[iii][i] = 0.0;  /* Current observation always rejected */
        }
        wsum += w[iii][i];
    }
    
    /* Normalize weights */
    xkndw = (double)(knobst - nobswithw0) / wsum;
    for (iii = 0; iii < knobst; iii++) {
        w[iii][i] = w[iii][i] * xkndw;
    }
    
    /* Output warning messages */
    if (!single_turbo) {
        fprintf(stderr, "WARNING:\n");
        fprintf(stderr, " Iteration# %3d Event#%4d Res#%3d =%.2f; "
                "ABS > 2.0 ---> weight set to zero !\n",
                nitt, isingle, nobs, res[nobs - 1][i]);
    }
    printf(" Iteration# %3d Event#%4d Res#%3d =%.2f; "
           "ABS > 2.0 ---> weight set to zero !\n",
           nitt, isingle, nobs, res[nobs - 1][i]);
    
    if (!single_turbo) {
        fprintf(stderr, "knobs(i)   = %d\n", knobst);
        fprintf(stderr, "nobswithw0 = %d\n", nobswithw0);
    }
    
    /* Signal to inhibit pre-stopping in OUTPUT */
    *iresflag = 1;
}
