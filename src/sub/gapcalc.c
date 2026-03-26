/*
 * Calculate azimuthal gap for earthquake location
 * Determines largest gap between stations around hypocenter
 * from Fortran GAPCALC subroutine
 */

#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "../include/globals.h"

/* External global variables from vel_com.inc */
extern int knobs[IEQ];
extern bool single_turbo;
extern double w[IST][IEQ];
extern int istm[IST][IEQ];
extern double e[5][IEQ];
extern double x[IST][3];
extern int igap[IEQ];

/* External function */
extern void sorti(int arr[], int n);

void gapcalc(int i) {
    int nofgaps, j, ig;
    double xstn, ystn, xhyp, yhyp, dxstnhyp, dystnhyp;
    int iga[MAXOBSPEREVENT];
    
    nofgaps = 0;
    
    if (i < 0 || i >= IEQ) {
        return;
    }

    /* Check array size */
    if (knobs[i] > MAXOBSPEREVENT) {
        if (!single_turbo) {
            fprintf(stderr, "WARNING: Event# %d Nobs = %d > %d\n",
                    i, knobs[i], MAXOBSPEREVENT);
        }
        printf("Event# %d Nobs = %d > %d\n", i, knobs[i], MAXOBSPEREVENT);
        return;
    }
    
    /* Loop through observations and compute bearings */
    for (j = 0; j < knobs[i]; j++) {
        if (w[j][i] > 0.0) {
            nofgaps++;
            xstn = x[istm[j][i]][0];
            ystn = x[istm[j][i]][1];
            xhyp = e[1][i];  /* X coordinate of event */
            yhyp = e[2][i];  /* Y coordinate of event */
            
            dxstnhyp = fabs(xstn - xhyp);
            dystnhyp = fabs(ystn - yhyp);
            
            if (dxstnhyp > 0.0001f || dystnhyp > 0.0001f) {
                /* Compute bearing (degrees from north, clockwise) */
                iga[nofgaps - 1] = (int)(57.296f * atan2f(xstn - xhyp, ystn - yhyp));
            } else {
                iga[nofgaps - 1] = 359;
            }
            
            if (iga[nofgaps - 1] < 0) {
                iga[nofgaps - 1] = iga[nofgaps - 1] + 360;
            }
        }
    }
    
    if (nofgaps > 0) {
        sorti(iga, nofgaps);
    } else {
        printf("WARNING: Event-# %d has zero observations!\n", i);
        if (!single_turbo) {
            fprintf(stderr, "WARNING: Event-# %d has zero observations!\n", i);
        }
        return;
    }
    
    /* Calculate maximum gap */
    igap[i] = 0;
    ig = iga[0] - iga[nofgaps - 1];
    if (ig < 0) ig = ig + 360;
    if (ig > igap[i]) igap[i] = ig;
    
    for (j = 1; j < nofgaps; j++) {
        ig = iga[j] - iga[j - 1];
        if (ig < 0) ig = ig + 360;
        if (ig > igap[i]) igap[i] = ig;
    }
}
