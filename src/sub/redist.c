/*
 * Convert rectangular coordinates (km) to lat/lon
 * Inverse of coordinate system transformation
 * from Fortran REDIST subroutine
 */

#include <math.h>
#include <stdio.h>

/* External global variables from GEO_COORSYSTEM common block */
extern double olat, olon;
extern double rearth, ellip, rlatc, rad;
extern double aa, bb, bc, sint, cost, rotate;
extern int icoordsystem;

void redist(double xkm, double ykm, double *xlat, double *xlon) {
    double xx, yy, y, x, q, lat, yp, bcl, p, lon;
    double lat1, lat2, lat3, clat1 = 0.0;
    
    xx = (double)xkm;
    yy = (double)ykm;
    
    /* Undo rotation to original coordinate orientation */
    x = xx * cost - yy * sint;
    y = yy * cost + xx * sint;
    
    /* Check for division by zero */
    if (fabs((double)aa) < 1e-7) {
        printf("REDIST: AA=%g BB=%g COS(LAT1)=%g >> DIVISION BY ZERO\n",
               aa, bb, clat1);
        return;
    }
    
    /* Convert Y coordinate to latitude */
    q = y / aa;
    lat = (q + olat) / 60.0;
    *xlat = (double)(q + olat - 60.0 * lat);
    yp = 60.0 * lat + (*xlat);
    
    lat1 = atan(rlatc * tan(yp * rad / 60.0));
    lat2 = atan(rlatc * tan(olat * rad / 60.0));
    lat3 = (lat1 + lat2) / 2.0;
    clat1 = cos(lat3);
    bcl = bb * clat1;
    
    if (fabs((double)bcl) < 1e-6) {
        printf("REDIST: AA=%g BB=%g COS(LAT1)=%g >> DIVISION BY ZERO\n",
               aa, bb, clat1);
        return;
    }
    
    /* Convert X coordinate to longitude */
    p = x / (bb * clat1);
    lon = (p + olon) / 60.0;
    *xlon = (double)(p + olon - 60.0 * lon);
    
    /* Combine degrees and minutes */
    *xlat = (double)(lat + (*xlat) / 60.0);
    *xlon = (double)(lon + (*xlon) / 60.0);
}
