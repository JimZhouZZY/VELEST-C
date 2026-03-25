#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern void layers(int nx, int ny, int nz, const double x[], const double y[], const double z[],
                   const double vel[][100][100], double xe, double ye, double xr, double yr,
                   int nl, double v[], double vsq[]);
extern void stpath(double xe, double ye, double ze, double xr, double yr, double zr,
                   double rp[3][200], int *nrp, double *tt, int nl, int jl, double tkj,
                   double v[], double d[], double thk[], double *sterr);
extern void direct1(int nl, const double v[], const double vsq[], const double thk[],
                    int jl, double tkj, double delta, double depth, double *tdir, double *u, double *x);
extern void dirpath(double xe, double ye, double ze, double xr, double yr, double zr,
                    double delta, int nl, double v[], double vsq[], double thk[],
                    int jl, double tkj, double salpha, double deljl,
                    double rp[3][200], int *nrp, double *direrr);
extern void reflect1(int nl, const double v[], const double vsq[], const double thk[],
                     int jl, double tkj, double delta, double z, int mll,
                     double *trefl, double *ain, int *ierr, double *dtddrefl, double *dtdhrefl);
extern void reflectpath(double xe, double ye, double ze, double xr, double yr, double zr,
                        double delta, int nl, const double v[], const double vsq[], const double thk[],
                        int jl, double tkj, double ain, int mll,
                        double rp[3][200], int *nrp, double *reflerr);
extern void refract(int nl, const double v[], const double vsq[], const double thk[],
                    int jl, double tkj, double delta, int *kk, double *tref, double *didjkk, double *xovmax);
extern void refpath(double xe, double ye, double ze, double xr, double yr, double zr,
                    double delta, int nl, double v[], double vsq[], double thk[],
                    double tkj, int jl, int kk, double didjkk,
                    double rp[3][200], int *nrp, double *refrerr);

void raypath(int nx, int ny, int nz,
             const double x[], const double y[], const double z[], const double vel[][100][100],
             int nl, double thk[], double d[], double v[], double vsq[],
             double xe, double ye, double ze, double xr, double yr, double zr,
             double rp[3][200], int *nrp, int *nrtn, int *jl,
             double *tkj, int itype, double *ttt, int *mll,
             double *sterr, double *direrr, double *refrerr, double *reflerr,
             double *dtddrefl, double *dtdhrefl) {
retry:
    *sterr = 0.0f;
    *direrr = 0.0f;
    *refrerr = 0.0f;
    *reflerr = 0.0f;

    if (ze < d[0] || zr < d[0]) {
        abort();
    }

    double depth = ze - zr;
    double delta = sqrtf((xr - xe) * (xr - xe) + (yr - ye) * (yr - ye));

    int l = nl - 1;
    while (l >= 0 && d[l] >= ze) {
        l--;
    }
    *jl = l;

    d[0] = zr;
    thk[0] = d[1] - zr;
    *tkj = ze - d[*jl];

    if (delta < 0.05f) {
        stpath(xe, ye, ze, xr, yr, zr, rp, nrp, ttt, nl, *jl, *tkj, v, d, thk, sterr);
        *nrtn = 1;
        return;
    }

    if (itype != 1) {
        layers(nx, ny, nz, x, y, z, vel, xe, ye, xr, yr, nl, v, vsq);
    }

    if (*jl == nl - 1) {
        double tdir = 0.0f;
        double salpha = 0.0f;
        double deljl = 0.0f;
        direct1(nl, v, vsq, thk, *jl, *tkj, delta, depth, &tdir, &salpha, &deljl);
        dirpath(xe, ye, ze, xr, yr, zr, delta, nl, v, vsq, thk, *jl, *tkj, salpha, deljl, rp, nrp, direrr);
        *nrtn = 2;
        *ttt = tdir;
        return;
    }

    if (*mll > 0) {
        double trefl = 0.0f;
        double ain = 0.0f;
        int ierr = 0;
        reflect1(nl, v, vsq, thk, *jl, *tkj, delta, depth, *mll, &trefl, &ain, &ierr, dtddrefl, dtdhrefl);
        if (ierr != 0) {
            *mll = 0;
            goto retry;
        }
        reflectpath(xe, ye, ze, xr, yr, zr, delta, nl, v, vsq, thk, *jl, *tkj, ain, *mll, rp, nrp, reflerr);
        *ttt = trefl;
        *nrtn = 8;
        return;
    }

    int kk = -1;
    double tref = 0.0f, didjkk = 0.0f, xovmax = 0.0f;
    refract(nl, v, vsq, thk, *jl, *tkj, delta, &kk, &tref, &didjkk, &xovmax);

    if (delta > xovmax) {
        refpath(xe, ye, ze, xr, yr, zr, delta, nl, v, vsq, thk, *tkj, *jl, kk, didjkk, rp, nrp, refrerr);
        *nrtn = 3;
        *ttt = tref;
        return;
    }

    if (*jl == 0) {
        double tdir = sqrtf((*tkj) * (*tkj) + delta * delta) / v[0];
        if (tref < tdir) {
            refpath(xe, ye, ze, xr, yr, zr, delta, nl, v, vsq, thk, *tkj, *jl, kk, didjkk, rp, nrp, refrerr);
            *nrtn = 4;
            *ttt = tref;
        } else {
            stpath(xe, ye, ze, xr, yr, zr, rp, nrp, ttt, nl, *jl, *tkj, v, d, thk, sterr);
            *nrtn = 5;
            *ttt = tdir;
        }
        return;
    }

    double tdir = 0.0f;
    double salpha = 0.0f;
    double deljl = 0.0f;
    direct1(nl, v, vsq, thk, *jl, *tkj, delta, depth, &tdir, &salpha, &deljl);

    if (tref < tdir) {
        refpath(xe, ye, ze, xr, yr, zr, delta, nl, v, vsq, thk, *tkj, *jl, kk, didjkk, rp, nrp, refrerr);
        *ttt = tref;
        *nrtn = 6;
    } else {
        dirpath(xe, ye, ze, xr, yr, zr, delta, nl, v, vsq, thk, *jl, *tkj, salpha, deljl, rp, nrp, direrr);
        *ttt = tdir;
        *nrtn = 7;
    }
    
    /*
    printf("DEBUG raypath: nrtn=%d, ttt=%.6f, nrp=%d, jl=%d, tkj=%.6f\n", *nrtn, *ttt, *nrp, *jl, *tkj);
    printf("DEBUG raypath: delta=%.6f, depth=%.6f, xe=%.6f, ye=%.6f, ze=%.6f\n", delta, depth, xe, ye, ze);
    printf("DEBUG raypath: xr=%.6f, yr=%.6f, zr=%.6f\n", xr, yr, zr);
    */
}
