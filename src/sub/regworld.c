#include <math.h>
#include <string.h>

extern int indfe[];
extern int lt50[];
extern int lt25[];
extern char irname[];
extern double xnrlon[];

void regworld(int ityp, double cord1, double cord2, char *place, int *nreg) {
    const double galat = 46.7706f;
    const double galon = 8.09428f;

    double gelat = 0.0;
    double gelon = 0.0;

    if (ityp <= 0) {
        gelat = cord1;
        gelon = cord2;
    } else {
        double pi = 3.141593f;
        double pirad = pi / 180.0;
        double galor = galon * pirad;
        double galar = galat * pirad;
        double xazr = cord2 * pirad;
        double xdr = cord1 * pirad;

        double q = sinf(galar) * cosf(xdr) + cosf(galar) * sinf(xdr) * cosf(xazr);
        if (fabs(q) > 1.0f) {
            place[0] = '\0';
            return;
        }

        double gelac = asinf(q);
        double xlatr = gelac;
        double gelor = 0.0;

        if (cosf(gelac) != 0.0 && cosf(galar) != 0.0) {
            double ck1 = (cosf(xdr) - sinf(galar) * sinf(xlatr)) / (cosf(galar) * cosf(xlatr));
            double ck2 = (sinf(xazr) * sinf(xdr)) / cosf(xlatr);

            if (fabs(ck1) <= 1.0f && fabs(ck2) <= 1.0f) {
                if (ck1 < 0.0 && ck2 >= 0.0) gelor = galor + acosf(ck1);
                else if (ck1 >= 0.0 && ck2 < 0.0) gelor = galor + asinf(ck2);
                else if (ck1 >= 0.0 && ck2 >= 0.0) gelor = galor + acosf(ck1);
                else gelor = galor - pi - asinf(ck2);
            }
        }

        double alon = fabs(gelor);
        while (alon >= 2.0f * pi) {
            alon -= 2.0f * pi;
        }
        gelor = copysignf(alon, gelor);
        if (fabs(gelor) > pi) gelor -= copysignf(2.0f * pi, gelor);

        gelon = gelor / pirad;
        gelat = gelac / pirad;
    }

    memset(place, ' ', 32);
    place[32] = '\0';

    double xlat = gelat;
    double xlon = gelon;
    if (xlon < 0.0) xlon += 360.0;

    int mhemis = 0;
    if (xlat >= 0.0) {
        if (xlon > 180.0) {
            mhemis = 3600;
            xlon = 360.0 - xlon;
        }
    } else {
        if (xlon > 180.0) {
            mhemis = 10800;
            xlon = 360.0 - xlon;
        } else {
            mhemis = 7200;
        }
    }

    int latin = (fabs(xlat) < 1.0f) ? 0 : (int)fabs(xlat);
    if (latin > 89) latin = 89;

    int nlo = (int)fabs(xlon);
    int i0 = mhemis + 40 * latin;

    double buf[40];
    for (int i = 0; i < 40; ++i) buf[i] = xnrlon[i0 + i];

    int nval = (int)(buf[0] * 1000.0 + 1.1f);
    int idx = 1;
    for (; idx < nval; ++idx) {
        if (nlo < (int)buf[idx]) break;
    }

    double xreg = buf[idx - 1] - (int)buf[idx - 1];
    *nreg = (int)(xreg * 1000.0 + 0.1f);
    int i1 = indfe[*nreg];
    int i2 = indfe[*nreg + 1];
    int n = i2 - i1;
    if (n > 32) n = 32;
    memcpy(place, irname + i1, (size_t)n);
    if (n < 32) place[n] = '\0';
}

void matrinv(int n, const double *a, double *b) {
    b[0] = 1.0f / a[0];
    if (n <= 1) return;

    int nn = n * n;
    for (int i = 1; i < nn; ++i) {
        b[i] = 0.0;
    }

    int mm = 0;
    int kn = 0;
    for (int m = 2; m <= n; ++m) {
        int k = m - 1;
        mm = mm + n + 1;
        kn = kn + n;

        double ek = a[mm - 1];
        int mi = m - n;
        for (int i = 1; i <= k; ++i) {
            mi = mi + n;
            int ij = i - n;
            int jm = kn;
            for (int j = 1; j <= k; ++j) {
                ij = ij + n;
                jm = jm + 1;
                ek = ek - a[mi - 1] * b[ij - 1] * a[jm - 1];
            }
        }

        b[mm - 1] = 1.0f / ek;
        mi = m - n;
        int im = kn;
        for (int i = 1; i <= k; ++i) {
            im = im + 1;
            int ij = i - n;
            int mj = kn;
            for (int j = 1; j <= k; ++j) {
                ij = ij + n;
                mj = mj + 1;
                b[im - 1] = b[im - 1] - b[ij - 1] * a[mj - 1] * b[mm - 1];
            }
            mi = mi + n;
            b[mi - 1] = b[im - 1];
        }

        im = kn;
        for (int i = 1; i <= k; ++i) {
            im = im + 1;
            int mj = m - n;
            int ij = i - n;
            for (int j = 1; j <= k; ++j) {
                mj = mj + n;
                ij = ij + n;
                b[ij - 1] = b[ij - 1] + b[im - 1] * b[mj - 1] * ek;
            }
        }
    }
}
