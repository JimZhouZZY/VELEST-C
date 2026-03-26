#include <math.h>
#include <stdio.h>

extern double olat;
extern double olon;
extern double rearth;
extern double ellip;
extern double rlatc;
extern double rad;
extern double aa;
extern double bb;
extern double bc;
extern double sint;
extern double cost;
extern double rotate;
extern int icoordsystem;

void setorg(double orlat, double orlon, double rrotate, int ifil) {
    (void)icoordsystem;

    rotate = rrotate;

    if (orlat == 0.0 && orlon == 0.0) {
        olat = 46.95240f;
        olon = -7.439583f;
    } else {
        olat = orlat;
        olon = orlon;
    }

    olat *= 60.0;
    olon *= 60.0;

    rad = 0.017453292;

    rearth = 6378.135;
    ellip = 298.26;

    double phi = olat * rad / 60.0;
    double beta = phi - sin(phi * 2.0) / ellip;
    rlatc = tan(beta) / tan(phi);

    if (ifil > 0) {
        FILE *out = stdout;
        fprintf(out, "\n\n");
        fprintf(out, "SHORT DISTANCE CONVERSION on ELLIPSOIDE of WORLD GEODETIC SYSTEM 1972 (WGS72)\n");
        fprintf(out, "================================================================================\n\n");
        fprintf(out, " Radius at equator (REARTH)= %10.5f  km\n", rearth);
        fprintf(out, "   1. / (ellipticity)      = %10.3f\n\n", ellip);
        fprintf(out, "Origin of cartesian coordinates [degrees]:\n");
        if (orlat == 0.0 && orlon == 0.0) {
            fprintf(out, " (Origin = city of BERNE, Switzerland)\n");
        }
        fprintf(out, " %12.7f N     %12.7f W\n\n", olat / 60.0, olon / 60.0);
        fprintf(out, " Rotation angle (in degr.) clockwise from\n");
        fprintf(out, "North  rotate= %6.1f\n\n", rotate);
    }

    double lat1 = atan(rlatc * tan(olat * rad / 60.0));
    double lat2 = atan(rlatc * tan((olat + 1.0f) * rad / 60.0));
    double dela = lat2 - lat1;
    double r = rearth * (1.0 - (sin(lat1) * sin(lat1)) / ellip);

    aa = (double)(r * dela);

    double delb = acos(sin(lat1) * sin(lat1) + cos(rad / 60.0) * cos(lat1) * cos(lat1));
    bc = (double)(r * delb);
    bb = (double)(r * delb / cos(lat1));

    if (ifil > 0) {
        FILE *out = stdout;
        fprintf(out, " Radius of sphere at OLAT = %10.3f km\n\n", r);
        fprintf(out, "Conversion of GEOGRAPHICAL LATITUDE to GEOCENTRICAL LATITUDE:\n");
        fprintf(out, "RLATC = TAN(GEOCENTR.LAT) / TAN(GEOGRAPH.LAT)\n");
        fprintf(out, " RLATC = %12.8f\n\n", rlatc);
        fprintf(out, "          Short distance conversions:\n");
        fprintf(out, "          one min lat ~ %7.4f km \n", aa);
        fprintf(out, "          one min lon ~ %7.4f km \n\n\n", bc);
    }

    sint = sinf(rotate * (double)rad);
    cost = cosf(rotate * (double)rad);
}
