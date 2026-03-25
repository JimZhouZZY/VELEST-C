#include <math.h>
#include <stdio.h>

void reflect(int nl, const double v[], const double vsq[], const double thk[],
             int *jl, double *tkj, double delta, double z, int mll,
             double *trefl, double *ain, int *ierr,
             double *dtddrefl, double *dtdhrefl) {
    double vi[100] = {0.0f};
    double vsqu[100] = {0.0f};
    double div[100] = {0.0f};

    *ierr = 0;
    *dtddrefl = 0.0f;
    *dtdhrefl = 0.0f;

    const double test4 = 0.01f;

    for (int i = 0; i < nl; ++i) {
        vi[i] = 1.0f / v[i];
    }

    double sum = 0.0f;
    *jl = 0;
    while (*jl < nl) {
        sum += thk[*jl];
        if (z <= sum) {
            break;
        }
        (*jl)++;
    }

    if (*jl >= nl) {
        *jl = nl - 1;
        *tkj = z - sum;
    } else {
        *tkj = z - sum + thk[*jl];
    }

    if (*jl > mll) {
        *ierr = -1;
        fprintf(stdout, "WARNING: subr. REFLECT >>> hypocenter is below reflector!\n");
    }
    if (mll >= nl) {
        *ierr = 50;
        return;
    }

    int nit = 0;
    double zbig = *tkj;
    sum = 0.0f;

    for (int m = 0; m <= mll; ++m) {
        double fac = 1.0f;
        if (*jl <= m) {
            fac = 2.0f;
            if (*jl == m) {
                fac = (2.0f * thk[m] - zbig) / thk[m];
            }
        }
        div[m] = fac * thk[m];
        vsqu[m] = vsq[mll] / vsq[m];
        sum += div[m];
    }

    double a = atan2f(delta, sum);
    double sina = 0.0f;
    double sina2 = 0.0f;
    double da = 0.0f, da2 = 0.0f;
    double dddas = 0.0f;

    while (1) {
        sina = sinf(a);
        if (sina == 1.0f) {
            sina = 0.9999999f;
        }
        sina2 = sina * sina;

        double dellit = 0.0f;
        dddas = 0.0f;

        for (int m = 0; m <= mll; ++m) {
            double sqr = vsqu[m] - sina2;
            if (sqr <= 0.0f) {
                *ierr = 50;
                return;
            }
            double ddda = div[m] / sqrtf(sqr);
            dellit += ddda;
            dddas += vsqu[m] * ddda / sqr;
        }

        dellit *= sina;
        dddas *= cosf(a);
        da = dellit - delta;
        da2 = da * da;

        if (da2 <= test4) {
            break;
        }

        nit++;
        if (nit >= 50) {
            *ierr = 50;
            return;
        }

        a = a - da / dddas;
        if (a >= 1.570796f) {
            a = 1.5f;
        }
    }

    double tt = 0.0f;
    double dtda = 0.0f;

    for (int m = 0; m <= mll; ++m) {
        double ti = div[m] * vi[m] / sqrtf(1.0f - sina2 / vsqu[m]);
        double dti = ti / (vsqu[m] - sina2);
        tt += ti;
        dtda += dti;
    }

    dtda = dtda * sina * cosf(a);
    double dtdd = dtda / dddas;
    double dtdh = -(1.0f - v[*jl] * sina * dtdd) / (v[*jl] * sqrtf(1.0f - sina2));

    double anin = sina;
    *ain = anin;

    if (mll > *jl) {
        for (int i = mll - 1; i >= *jl; --i) {
            *ain = (*ain) * v[i] / v[i + 1];
        }
    }

    *trefl = tt;
    *dtddrefl = dtdd;
    *dtdhrefl = dtdh;
}
