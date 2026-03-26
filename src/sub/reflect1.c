#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void reflect1(int nl, const double v[], const double vsq[], const double thk[],
              int jl, double tkj, double delta, double z, int mll,
              double *trefl, double *ain, int *ierr,
              double *dtddrefl, double *dtdhrefl) {
    (void)z;

    *ierr = 0;
    *dtddrefl = 0.0;
    *dtdhrefl = 0.0;

    double div[100] = {0.0};

    if (jl > mll) {
        *ierr = -1;
        fprintf(stdout, "WARNING: subr. REFLECT1 >>> hypocenter is below reflector!\n");
    }
    if (mll >= nl) {
        abort();
    }

    for (int m = 0; m <= mll; ++m) {
        double fac = 1.0f;
        if (jl <= m) {
            fac = 2.0f;
            if (jl == m) {
                fac = (2.0f * thk[m] - tkj) / thk[m];
            }
        }
        div[m] = fac * thk[m];
    }

    double depth = 0.0;
    for (int i = 0; i <= mll; ++i) {
        depth += div[i];
    }

    int lmax = mll;
    double tklmax = div[mll];
    double vlmax = v[mll];

    for (int l = 0; l < mll; ++l) {
        if (v[l] > vlmax) {
            lmax = l;
            tklmax = div[l];
            vlmax = v[l];
        }
    }

    if (tklmax <= 0.05f) {
        tklmax = 0.05f;
    }

    double ua = (v[mll] / vlmax) * delta / sqrtf(delta * delta + depth * depth);
    double ub = (v[mll] / vlmax) * delta / sqrtf(delta * delta + tklmax * tklmax);

    double uasq = ua * ua;
    double ubsq = ub * ub;
    if (ubsq >= 1.0f) ubsq = 0.99999f;
    if (uasq >= 1.0f) uasq = 0.99999f;

    double xa = div[mll] * ua / sqrtf(1.0f - uasq);
    double xb = (lmax == mll) ? delta : div[mll] * ub / sqrtf(1.0f - ubsq);

    double dela = xa;
    double delb = xb;

    for (int l = 0; l < mll; ++l) {
        dela += div[l] * ua / sqrtf(vsq[mll] / vsq[l] - uasq);
        double ubdiv = sqrtf(vsq[mll] / vsq[l] - ubsq);
        if (ubdiv <= 1.e-20f) {
            ubdiv = 1.e-20f;
        }
        delb += div[l] * ub / ubdiv;
    }

    double del = 0.0;
    double x = 0.0;
    double u = 0.0;
    double usq = 0.0;

    for (int kount = 0; kount < 25; ++kount) {
        if ((delb - dela) < 0.02f) {
            x = 0.5f * (xa + xb);
            u = x / sqrtf(x * x + div[mll] * div[mll]);
            usq = u * u;
            break;
        }

        x = xa + (delta - dela) * (xb - xa) / (delb - dela);
        u = x / sqrtf(x * x + div[mll] * div[mll]);
        usq = u * u;

        del = x;
        for (int l = 0; l < mll; ++l) {
            del += div[l] * u / sqrtf(vsq[mll] / vsq[l] - usq);
        }

        double xtest = del - delta;
        if (fabs(xtest) < 0.02f) {
            break;
        }

        if (xtest < 0.0) {
            xa = x;
            dela = del;
        } else {
            xb = x;
            delb = del;
        }
    }

    if (del == 0.0) {
        del = x;
    }

    double tdir = sqrtf(x * x + div[mll] * div[mll]) / v[mll];
    for (int l = 0; l < mll; ++l) {
        tdir += div[l] * v[mll] / (vsq[l] * sqrtf(vsq[mll] / vsq[l] - usq));
    }
    tdir -= (u / v[mll]) * (del - delta);

    *trefl = tdir;
    *ain = u;

    if (mll > jl) {
        for (int i = mll - 1; i >= jl; --i) {
            *ain = (*ain) * v[i] / v[i + 1];
        }
    }
}
