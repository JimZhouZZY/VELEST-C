#include <math.h>
#include <stdio.h>

extern double olat, olon;
extern double rlatc, rad;
extern double aa, bb, sint, cost;

void dist(double xlat, double xlon, double *xkm, double *ykm) {
    fprintf(stderr ,"DEBUG dist inputs:\n");
    fprintf(stderr ,"  xlat=%.10f xlon=%.10f\n", xlat, xlon);

    fprintf(stderr ,"DEBUG globals:\n");
    fprintf(stderr ,"  olat=%.10f olon=%.10f\n", olat, olon);
    fprintf(stderr ,"  rlatc=%.15f rad=%.15f\n", rlatc, rad);
    fprintf(stderr ,"  aa=%.10f bb=%.10f\n", aa, bb);
    fprintf(stderr ,"  cost=%.10f sint=%.10f\n", cost, sint);
    double q, yp, xx;
    double lat1, lat2, lat3;

    if (!xkm || !ykm) {
        return;
    }

    q = 60.0 * xlat - olat;
    yp = q + olat;
    lat1 = atan(rlatc * tan(rad * yp / 60.0));
    lat2 = atan(rlatc * tan(rad * olat / 60.0));
    lat3 = (lat2 + lat1) / 2.0;

    xx = 60.0 * xlon - olon;
    q = q * aa;
    xx = xx * bb * cos(lat3);

    yp = cost * q - sint * xx;
    xx = cost * xx + sint * q;
    q = yp;

    *xkm = xx;
    *ykm = q;
    fprintf(stderr, "DEBUG intermediate variables:\n");
    fprintf(stderr, "  q=%.10f yp=%.10f xx=%.10f\n", q, yp, xx);
    fprintf(stderr, "  lat1=%.15f lat2=%.15f lat3=%.15f\n", lat1, lat2, lat3);
}
