/*
 * Direct ray path - compute raypoint coordinates for direct seismic wave
 * Assigns raypoints in a flat layered earth model
 * from Fortran DIRPATH subroutine
 */

#include <stddef.h>
#include <math.h>

void dirpath(double xe, double ye, double ze, double xr, double yr, double zr,
             double delta, int nl, double v[], double vsq[], double thk[],
             int jl, double tkj, double salpha, double deljl,
             double rp[3][200], int *nrp, double *direrr) {
    double d1, d2, tng;
    int nrp1, i, m;
    
    (void)nl;  /* Parameter used for array dimensioning only in Fortran */
    
    /* Compute direction cosines for horizontal plane */
    d1 = (xr - xe) / delta;
    d2 = (yr - ye) / delta;
    
    /* Number of ray points = depth layers + 1 (surface) */
    *nrp = jl + 1;
    
    /* Assign first ray point (event location) */
    rp[0][0] = xe;
    rp[1][0] = ye;
    rp[2][0] = ze;
    
    /* Assign second ray point (bottom of event layer) */
    rp[0][1] = rp[0][0] + deljl * d1;
    rp[1][1] = rp[1][0] + deljl * d2;
    rp[2][1] = rp[2][0] - tkj;
    
    /* Calculate intermediate ray points at layer boundaries */
    nrp1 = *nrp;  /* Include up to surface layer for error calculation */
    for (i = 2; i < nrp1; i++) {
        m = jl - i;  /* Layer number (0-based) below ray point i */
        
        /* Tangent of incidence angle in layer m */
        double denom = sqrtf(vsq[jl] - vsq[m] * salpha * salpha);
        tng = v[m] * salpha / denom;
        
        /* Ray point coordinates */
        rp[0][i] = rp[0][i - 1] + thk[m] * tng * d1;
        rp[1][i] = rp[1][i - 1] + thk[m] * tng * d2;
        rp[2][i] = rp[2][i - 1] - thk[m];
    }
    
    /* Compute error (distance from final ray point to receiver) */
    double dx = rp[0][*nrp - 1] - xr;
    double dy = rp[1][*nrp - 1] - yr;
    double dz = rp[2][*nrp - 1] - zr;
    *direrr = sqrtf(dx * dx + dy * dy + dz * dz) * 1000.0f;
    
    /* Assign final ray point (receiver location) */
    rp[0][*nrp - 1] = xr;
    rp[1][*nrp - 1] = yr;
    rp[2][*nrp - 1] = zr;
}
