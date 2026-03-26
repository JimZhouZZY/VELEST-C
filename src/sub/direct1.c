#include <math.h>
#include <stdio.h>

void direct1(int nl, const double v[], const double vsq[], const double thk[],
             int jl, double tkj, double delta, double depth,
             double *tdir, double *u, double *x) {
    (void)nl;

    double del = 0.0;

    int lmax = jl;
    double tklmax = tkj;
    double vlmax = v[jl];

    for (int l = 0; l < jl; ++l) {
        if (v[l] > vlmax) {
            lmax = l;
            tklmax = thk[l];
            vlmax = v[l];
        }
    }

    if (tklmax <= 0.05f) {
        tklmax = 0.05f;
    }

    double ua = (v[jl] / vlmax) * delta / sqrtf(delta * delta + depth * depth);
    double ub = (v[jl] / vlmax) * delta / sqrtf(delta * delta + tklmax * tklmax);

    double uasq = ua * ua;
    double ubsq = ub * ub;
    if (ubsq >= 1.0f) ubsq = 0.99999f;
    if (uasq >= 1.0f) uasq = 0.99999f;

    double xa = tkj * ua / sqrtf(1.0f - uasq);
    double xb = (lmax == jl) ? delta : (tkj * ub / sqrtf(1.0f - ubsq));

    double dela = xa;
    double delb = xb;

    for (int l = 0; l < jl; ++l) {
        dela += thk[l] * ua / sqrtf(vsq[jl] / vsq[l] - uasq);
        double ubdiv = sqrtf(vsq[jl] / vsq[l] - ubsq);
        if (ubdiv <= 1.0e-20f) {
            fprintf(stdout, "WARNING: subr. direct1 near-singular denominator\n");
            ubdiv = 1.0e-20f;
        }
        delb += thk[l] * ub / ubdiv;
    }

    double usq = 0.0;
    for (int kount = 0; kount < 25; ++kount) {
        if ((delb - dela) < 0.02f) {
            *x = 0.5f * (xa + xb);
            *u = *x / sqrtf((*x) * (*x) + tkj * tkj);
            usq = (*u) * (*u);
            break;
        }

        *x = xa + (delta - dela) * (xb - xa) / (delb - dela);
        *u = *x / sqrtf((*x) * (*x) + tkj * tkj);
        usq = (*u) * (*u);

        del = *x;
        for (int l = 0; l < jl; ++l) {
            del += thk[l] * (*u) / sqrtf(vsq[jl] / vsq[l] - usq);
        }

        double xtest = del - delta;
        if (fabs(xtest) < 0.02f) {
            break;
        }
        if (xtest < 0.0) {
            xa = *x;
            dela = del;
        } else {
            xb = *x;
            delb = del;
        }
    }

    if (del == 0.0) {
        del = *x;
    }

    *tdir = sqrtf((*x) * (*x) + tkj * tkj) / v[jl];
    for (int l = 0; l < jl; ++l) {
        *tdir += thk[l] * v[jl] / (vsq[l] * sqrtf(vsq[jl] / vsq[l] - usq));
    }
    *tdir -= (*u / v[jl]) * (del - delta);
}
