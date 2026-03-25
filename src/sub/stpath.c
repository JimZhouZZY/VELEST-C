/*
 * Store ray path points and compute travel time
 * Computes straight-line ray path from event to receiver
 * from Fortran STPATH subroutine
 */

#include <stddef.h>
#include <math.h>

void stpath(double xe, double ye, double ze, double xr, double yr, double zr,
            double rp[3][200], int *nrp_ptr, double *tt_ptr,
            int nl, int jl, double tkj, double v[], double d[], double thk[],
            double *sterr_ptr) {
    int nrp, jl1, j, j1, jb, jt;
    double tt;
    (void)nl;  /* Parameter not used; kept for API compatibility */
    
    nrp = jl + 1;  /* Number of ray points */
    
    /* Assign starting point (event location) */
    rp[0][0] = xe;      /* x of first point */
    rp[1][0] = ye;      /* y of first point */
    rp[2][0] = ze;      /* z of first point */
    
    /* If only 2 points (event and receiver), skip intermediate points */
    if (nrp == 2) {
        *nrp_ptr = nrp;
        *tt_ptr = tkj / v[jl];
        rp[0][nrp - 1] = xr;
        rp[1][nrp - 1] = yr;
        rp[2][nrp - 1] = zr;
        *sterr_ptr = 0.0;
        return;
    }
    
    /* Assign intermediate points at layer boundaries */
    jl1 = jl - 1;
    for (j = 0; j <= jl1; j++) {
        j1 = j + 1;
        rp[0][j1] = xe;
        rp[1][j1] = ye;
        jb = jl - j - 1;  /* 0-based layer index relative from bottom */
        rp[2][j1] = d[jb];
    }
    
    /* Check if receiver is at intermediate layer (compute error) */
    if (nrp > 2) {
        double dx = rp[0][nrp - 1] - xr;
        double dy = rp[1][nrp - 1] - yr;
        double dz = rp[2][nrp - 1] - zr;
        *sterr_ptr = sqrtf(dx * dx + dy * dy + dz * dz) * 1000.0;
    } else {
        *sterr_ptr = 0.0;
    }
    
    /* Assign ending point (receiver location) */
    rp[0][nrp - 1] = xr;
    rp[1][nrp - 1] = yr;
    rp[2][nrp - 1] = zr;
    
    /* Compute travel time */
    tt = tkj / v[jl];  /* Time in deepest layer */
    if (nrp == 2) {
        *nrp_ptr = nrp;
        *tt_ptr = tt;
        return;
    }
    
    /* Add travel time through intermediate layers */
    for (j = 0; j < jl1; j++) {
        jt = jl - j - 1;  /* 0-based layer index */
        tt = tt + thk[jt] / v[jt];
    }
    
    *nrp_ptr = nrp;
    *tt_ptr = tt;
}
