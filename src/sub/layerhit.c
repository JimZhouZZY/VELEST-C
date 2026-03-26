/*
 * Count ray hits in each layer and compute ray path statistics
 * from Fortran LAYERHIT subroutine
 */

#include <stddef.h>
#include <math.h>
#include "../include/globals.h"

/* External global variables from vel_com.inc */
extern int lmax;
extern int noheadwave;
extern double avhraylen;
extern double avvraylen;
extern int irefllayer[];
extern int irefrlayer[];
extern double refraylen[];
extern double hitlay[][3];
extern double h[];      /* Layer top depths */
extern double thk[];    /* Layer thicknesses */

void layerhit(double rp[3][200], int *nrpdeep, int nl, int nrp, int mll) {
    double rpmax, avraydepth;
    int j, jlay;
    
    /* Find deepest ray point (maximum depth) */
    rpmax = -999.0f;
    for (j = 0; j < nrp; j++) {
        if (rp[2][j] > rpmax) {
            rpmax = rp[2][j];
            *nrpdeep = j;  /* Index of deepest ray point (0-based) */
        }
    }
    
    rpmax = rpmax + 0.000001f;  /* Avoid exact layer boundary */
    
    /* Find deepest layer hit by this ray */
    lmax = 0;  /* Layer 0 by default */
    for (j = 1; j < nl; j++) {
        if (h[j] > rpmax) {
            lmax = j - 1;  /* Layer number (0-based) for depth rpmax */
            break;
        }
    }
    if (lmax == 0 && h[1] <= rpmax) {
        lmax = nl - 1;  /* Lowest layer if depth exceeds all layer tops */
    }
    
    /* Check if ray is horizontal (headwave) */
    if (*nrpdeep >= 0 && *nrpdeep < nrp - 1 &&
        rp[2][*nrpdeep] == rp[2][*nrpdeep + 1]) {
        irefrlayer[lmax]++;  /* Count refraction in layer LMAX */
        double hz_dist = sqrtf((rp[0][*nrpdeep] - rp[0][*nrpdeep + 1]) *
                              (rp[0][*nrpdeep] - rp[0][*nrpdeep + 1]) +
                              (rp[1][*nrpdeep] - rp[1][*nrpdeep + 1]) *
                              (rp[1][*nrpdeep] - rp[1][*nrpdeep + 1]));
        refraylen[lmax] = refraylen[lmax] + hz_dist;
    } else {
        /* Non-horizontal ray */
        if (mll == 0) {
            noheadwave++;
        } else if (lmax > 0) {
            irefllayer[lmax - 1]++;  /* Count reflection in layer above */
        }
    }
    
    /* Accumulate horizontal and vertical ray path lengths */
    double hz_tot = sqrtf((rp[0][0] - rp[0][nrp - 1]) * (rp[0][0] - rp[0][nrp - 1]) +
                         (rp[1][0] - rp[1][nrp - 1]) * (rp[1][0] - rp[1][nrp - 1]));
    avhraylen = avhraylen + hz_tot;
    
    double vz_tot = fabs(rp[2][*nrpdeep] - rp[2][nrp - 1]);
    avvraylen = avvraylen + vz_tot;
    
    /* Accumulate hit statistics for each ray segment */
    for (j = 1; j < nrp; j++) {
        avraydepth = (rp[2][j] + rp[2][j - 1]) / 2.0f;
        avraydepth = avraydepth + 0.000001f;  /* Avoid exact layer boundary */
        
        for (jlay = 0; jlay < nl; jlay++) {
            /* Check if ray segment is in layer jlay */
            if (avraydepth >= h[jlay] && avraydepth < (h[jlay] + thk[jlay])) {
                hitlay[jlay][0]++;  /* Count hits */
                
                double seg_hz = sqrtf((rp[0][j] - rp[0][j - 1]) *
                                     (rp[0][j] - rp[0][j - 1]) +
                                     (rp[1][j] - rp[1][j - 1]) *
                                     (rp[1][j] - rp[1][j - 1]));
                hitlay[jlay][1] = hitlay[jlay][1] + seg_hz;  /* Horizontal distance */
                
                double seg_vz = fabs(rp[2][j] - rp[2][j - 1]);
                hitlay[jlay][2] = hitlay[jlay][2] + seg_vz;  /* Vertical distance */
            }
        }
    }
}
