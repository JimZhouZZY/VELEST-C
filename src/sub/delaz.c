/*
 * Calculate geographic distance and azimuth between two points
 * Uses spherical trigonometry with Earth ellipsoid parameters
 * from Fortran DELAZ subroutine
 */

#include <stddef.h>
#include <math.h>

/* External global variables (common blocks) */
extern double rearth;   /* Earth radius */
extern double ellip;    /* Ellipsoid parameter */
extern double rlatc;    /* Flattening parameter (flat) */
extern double rad;      /* Conversion factor degrees to radians */

void delaz(double alat, double alon, double blat, double blon,
           double *del, double *dist, double *az) {
    double pi2 = 1.570796;  /* pi/2 in radians */
    double flat = rlatc;    /* Flattening parameter */
    double alatr, alonr, blatr, blonr;
    double tana, geoa, acol, tanb, geob, bcol;
    double diflon, cosdel, delr, top, den, azr, colat, radius;
    
    /* Convert latitude/longitude from degrees to radians */
    alatr = alat * rad;
    alonr = alon * rad;
    blatr = blat * rad;
    blonr = blon * rad;
    
    /* Convert geodetic latitudes to geocentric colatitudes */
    tana = flat * tan(alatr);
    geoa = atan(tana);
    acol = pi2 - geoa;
    
    tanb = flat * tan(blatr);
    geob = atan(tanb);
    bcol = pi2 - geob;
    
    /* Calculate angular distance (delta) */
    diflon = blonr - alonr;
    cosdel = sin(acol) * sin(bcol) * cos(diflon) + cos(acol) * cos(bcol);
    delr = acos(cosdel);
    
    /* Calculate azimuth from A to B */
    top = sin(diflon);
    den = (sin(acol) / tan(bcol)) - cos(diflon) * cos(acol);
    azr = atan2(top, den);
    
    /* Convert angular measurements back to degrees */
    *del = (double)(delr / rad);
    *az = (double)(azr / rad);
    if (*az < 0.0f) {
        *az = 360.0f + *az;
    }
    
    /* Compute distance in kilometers */
    colat = pi2 - (alatr + blatr) / 2.0;
    radius = 6378.163 * (1.0 + 3.35278e-3 * (1.0 / 3.0 - cos(colat) * cos(colat)));
    *dist = (double)(delr * radius);
}
