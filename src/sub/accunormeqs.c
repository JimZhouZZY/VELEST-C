#include <stdlib.h>

void outer(double *g, const double *rowofa, int nvar, const double *w);
void rside(double *rhs, const double *rowofa, int nvar, const double *res, const double *w);

void accunormeqs(double *rowofa, int nvar, double res, double w,
                 double *g, double *rhs) {
    outer(g, rowofa, nvar, &w);
    rside(rhs, rowofa, nvar, &res, &w);
}