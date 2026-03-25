/*
 * LSV G matrix core computation
 * Rotates a 2D vector (X, Y) by angle defined by DCOS, DSIN
 * from Fortran LSVG2 subroutine
 */

#include <stddef.h>

void lsvg2(double *dcos, double *dsin, double *x, double *y) {
    double xr = (*dcos) * (*x) + (*dsin) * (*y);
    *y = -(*dsin) * (*x) + (*dcos) * (*y);
    *x = xr;
}
