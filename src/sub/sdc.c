#include<stdio.h>

void dist(double xlat, double xlon, double *xkm, double *ykm);
void redist(double xkm, double ykm, double *xlat, double *xlon);

void sdc(double *x, double *y, double xlat, double xlon, int i) {
    if (!x || !y) {
        return;
    }
    if (i != -1 && i != 1) {
        return;
    }

    if (i == -1) {
        dist(xlat, xlon, x, y);
    } else {
        redist(xlat, xlon, x, y);
    }
    fprintf(stderr, "DEBUG: sdc x0 = %f \n", xlat);
    fprintf(stderr, "DEBUG: sdc y0 = %f \n", xlon);
    fprintf(stderr, "DEBUG: sdc i = %d \n", i);
    fprintf(stderr, "DEBUG: sdc x = %f \n", *x);
    fprintf(stderr, "DEBUG: sdc y = %f \n", *y);
}
