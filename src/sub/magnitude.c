/*
 * Calculate local magnitude for an earthquake event
 * Uses seismic amplitudes from multiple stations
 * from Fortran MAGNITUDE subroutine
 */

#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/* External global variables from vel_com.inc */
extern int nitt;
extern int knobs[];
extern double kpwt[][100];
extern char smn[][5];
extern double e[][658];
extern double d[][3][100];
extern int iyr[];
extern double amx[];
extern double prx[];
extern double xmagni[];
extern double xmagnitude;
extern double sdxmagnitude;
extern int nmag;
extern double delta;

/* External functions */
extern void stinp(int itime, char *staname, char *ifilt, int *iseis,
                  int *iscon, int *isdmp, int *isamp, double *cormag, int *ier);
extern void muk(double delta, double depth, char *ifilt, int iseis,
                double sconst, double sdampf, double voltgain, double cormag,
                double amx, double prx, double *xmag);

void magnitude(int ievent, int itime) {
    int iseis, iscon, isdmp, isamp, ier, i;
    double cormag, sconst, sdampf, voltgain, xmag;
    char ifilt[3];
    char staname[5];
    
    nmag = 0;
    xmagnitude = 0.0;
    sdxmagnitude = 0.0;
    
    /* Loop over all observations for this event */
    for (i = 0; i < knobs[ievent]; i++) {
        /* Skip observations with weight > 4 (not readings) */
        if (kpwt[i][ievent] > 4.0f) {
            continue;
        }
        
        /* Get station name from observation array */
        strncpy(staname, smn[i], 4);
        staname[4] = '\0';
        
        /* Get station instrument parameters */
        stinp(itime, staname, ifilt, &iseis, &iscon, &isdmp, &isamp,
              &cormag, &ier);
        
        if (ier < 0) {
            continue;  /* Seismic parameters not found for this station */
        }
        
        /* Update filter type if after 1984 and DE type */
        if (iyr[ievent] >= 1984 && ifilt[0] == 'D' && ifilt[1] == 'E') {
            ifilt[0] = 'A';
            ifilt[1] = 'D';
        }
        
        /* Extract station parameter scales */
        sconst = (double)iscon / 10.0;
        sdampf = (double)isdmp / 10.0;
        voltgain = (double)isamp;
        
        /* Calculate distance to station */
        double dx = e[ievent][2] - d[i][0][ievent];
        double dy = e[ievent][3] - d[i][1][ievent];
        delta = sqrtf(dx * dx + dy * dy);
        
        /* Calculate magnitude for this observation */
        xmag = 0.0;
        muk(delta, e[ievent][4], ifilt, iseis, sconst, sdampf,
            voltgain, cormag, amx[i], prx[i], &xmagni[i]);
        
        /* Skip invalid magnitude values */
        if (xmagni[i] == -13.0f) {
            continue;
        }
        
        /* Accumulate valid magnitudes */
        nmag++;
        xmag = xmagni[i];
        xmagnitude = xmagnitude + xmag;
        sdxmagnitude = sdxmagnitude + xmag * xmag;
    }
    
    /* Calculate mean and standard deviation */
    if (nmag != 0) {
        if (nmag >= 2) {
            sdxmagnitude = sqrtf((nmag * sdxmagnitude - xmagnitude * xmagnitude) /
                                 (nmag * (nmag - 1)));
        } else {
            sdxmagnitude = 0.0;
        }
        xmagnitude = xmagnitude / nmag;
    } else {
        xmagnitude = 0.0;
    }
}
