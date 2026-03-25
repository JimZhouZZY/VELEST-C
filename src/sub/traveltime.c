#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define IEQ 658
#define INSHOT 50
#define ITOTMODELS 2
#define INLTOT 100
#define IST 650
#define INRPMAX (2 * INLTOT + 2)

extern int isingle, neqs, nitt, nsp, ireflector, itopo, ifixsolution;
extern int nplay[INLTOT], iphase[IST][IEQ], istm[IST][IEQ], map2[INSHOT], map1[IST], kpwt[IST][IEQ], iain[IST];
extern int jl;
extern double tkj, delta;
extern double e[5][IEQ], d[IST][3][IEQ], x[IST][3], pt[IST][IEQ], ptcor[IST], stcor[IST], sphase[IST][IEQ];
extern double v[INLTOT], vsq[INLTOT], h[INLTOT], thk[INLTOT], vp[ITOTMODELS][INLTOT], hp[ITOTMODELS][INLTOT], thkp[ITOTMODELS][INLTOT];
extern double dtdr[3], dtdv[INLTOT], tctime[IST][IEQ], tcalc[IST], res[IST][IEQ], w[IST][IEQ];
extern double sterr, direrr, refrerr, reflerr, zmin, zmininput, vpvs;
extern char stn[IST][5];
extern char topo1file[81], topo2file[81];

extern void resetstatis(void);
extern void chtop(double xx, double yy, double *zk, const char *topo1file, const char *topo2file);
extern void raypath(int nx, int ny, int nz,
                    const double xarr[], const double yarr[], const double zarr[], const double vel[][100][100],
                    int nl, double thk[], double d[], double v[], double vsq[],
                    double xe, double ye, double ze, double xr, double yr, double zr,
                    double rp[3][200], int *nrp, int *nrtn, int *jl,
                    double *tkj, int itype, double *ttt, int *mll,
                    double *sterr, double *direrr, double *refrerr, double *reflerr,
                    double *dtddrefl, double *dtdhrefl);
extern void checkraypath(double *rp, int *nrp, int inrpmax);
extern void raypointcheck(const double *rp, int nrp, int inrpmax, const char *staname,
                          int isingle, void (*chtop_fn)(double, double, double *, const char *, const char *),
                          const char *topo1file, const char *topo2file);
extern void bendray(double rp[3][200], int nrp, const char *staname, double vtop, double *ttt);
extern void layerhit(double rp[3][200], int *nrpdeep, int nl, int nrp, int mll);
extern void travderiv(const char *raytype, int nl, int mll,
                      const double *v1, const double *vsq1,
                      const double rp[3][200], int nrp,
                      double *x2, double *y2, double *z2, double *ss,
                      double r1, double r2, int i, int nobs);
extern void resisave(int nrp, int nrpdeep, double rp[3][200], int nobs, int i, int k1, int mll);

void traveltime(int i, int nobs, int iresflag) {
    (void)iresflag;
    static int nittdone = -1;
    static int ifirst = 0;

    int ie = i - 1;
    int io = nobs - 1;
    //fprintf(stderr, "DEBUG traveltime: i=%d, nobs=%d, ie=%d, io=%d\n", i, nobs, ie, io);
    //fprintf(stderr, "DEBUG e: e[0]=%d, e[1]=%d, e[2]=%d, e[3]=%d, e[4]=%d\n", 
    //        e[0][ie], e[1][ie], e[2][ie], e[3][ie], e[4][ie]);

    if (isingle > 0 && kpwt[io][ie] == 5) {
        if (w[io][ie] != 0.0f) {
            fprintf(stderr, "TRAVELTIME>>> error in obs-wt!\n");
            exit(1);
        }
        res[io][ie] = 0.0f;
        for (int ii = 0; ii < 3; ++ii) dtdr[ii] = 0.0f;
        return;
    }

    if (isingle == 0) {
        if (i == 1 && nobs == 1 && nitt > nittdone) ifirst = 1;
        if (i == 1 && nobs == 1 && nitt == nittdone) ifirst = 2;
        if (nitt > nittdone || ifirst == 2) {
            nittdone = nitt;
            ifirst = 1;
            resetstatis();
        }
    }

    int do_s = 1;
    int k2 = iphase[io][ie];
    
    //fprintf(stderr, "DEBUG traveltime: i=%d, nobs=%d, ie=%d, io=%d, k2=%d\n", i, nobs, ie, io, k2);
    //fprintf(stderr, "DEBUG e: e[0]=%d, e[1]=%d, e[2]=%d, e[3]=%d, e[4]=%d\n", 
    //        e[0][ie], e[1][ie], e[2][ie], e[3][ie], e[4][ie]);
again:
    {
        int k2idx = k2 - 1;
        int nl = nplay[k2idx];
        int k1 = istm[io][ie];
        double v1[INLTOT], vsq1[INLTOT];
        for (int ii = 0; ii < nl; ++ii) {
            v[ii] = vp[k2idx][ii];
            v1[ii] = vp[k2idx][ii];
            if (nsp == 3 && sphase[io][ie] == 1.0f) v[ii] = vp[k2idx][ii] / vpvs;
            if (nsp == 3 && sphase[io][ie] == 2.0f && do_s == 1) v[ii] = vp[k2idx][ii] / vpvs;
            vsq[ii] = v[ii] * v[ii];
            vsq1[ii] = v1[ii] * v1[ii];
            h[ii] = hp[k2idx][ii];
            thk[ii] = thkp[k2idx][ii];
        }

        double ttt = 0.0f;
        double tkh = x[k1][2] - h[0];
        double r1 = e[1][ie] - x[k1][0];
        double r2 = e[2][ie] - x[k1][1];

        if (itopo > 0) {
            if (e[3][ie] < 0.0f) chtop(-e[1][ie], e[2][ie], &zmin, topo1file, topo2file);
            else zmin = zmininput;
        }

        if (i <= neqs) {
            if (e[3][ie] <= zmin) e[3][ie] = (isingle == 0) ? (zmin + 0.1f) : (zmin + 0.011f);
        } else {
            if (e[3][ie] <= zmin) e[3][ie] = zmin + 0.011f;
        }
        if (ifixsolution > 0 && e[3][ie] <= 0.0f) e[3][ie] = zmin + 0.001f;

        delta = sqrtf(r1 * r1 + r2 * r2);

        int mll = 0;
        if (ireflector > 0 && sphase[io][ie] == -1.0f) mll = ireflector;

        double rp[3][200];
        int nrp = 0, nrtn = 0;
        int nrpdeep = 0;
        double ster = 0.0f, direr = 0.0f, refrer = 0.0f, refler = 0.0f;
        double dtddrefl = 0.0f, dtdhrefl = 0.0f;

        static const double d1[1] = {1.0f};
        static double vel_dummy[1][100][100];
        raypath(1, 1, 1, d1, d1, d1, vel_dummy, nl, thk, h, v, vsq,
                e[1][ie], e[2][ie], e[3][ie], x[k1][0], x[k1][1], x[k1][2],
                rp, &nrp, &nrtn, &jl, &tkj, 1, &ttt, &mll,
                &ster, &direr, &refrer, &refler, &dtddrefl, &dtdhrefl);

        checkraypath(&rp[0][0], &nrp, INRPMAX);
        if (itopo == 2) raypointcheck(&rp[0][0], nrp, INRPMAX, stn[k1], isingle, chtop, topo1file, topo2file);
        if (itopo == 3) bendray(rp, nrp, stn[k1], v[0], &ttt);

        if (sphase[io][ie] != 2.0f || do_s == 0) {
            if (isingle == 0) {
                sterr += ster;
                direrr += direr;
                refrerr += refrer;
                reflerr += refler;
            }
            if (isingle == 0) layerhit(rp, &nrpdeep, nl, nrp, mll);
        }

        double x2[INRPMAX], y2[INRPMAX], z2[INRPMAX], ss[INRPMAX];
        if (nrtn == 1 || nrtn == 2 || nrtn == 5 || nrtn == 7) {
            travderiv("direct", nl, mll, v1, vsq1, rp, nrp, x2, y2, z2, ss, r1, r2, i, nobs);
        } else if (nrtn == 3 || nrtn == 4 || nrtn == 6) {
            travderiv("refracted", nl, mll, v1, vsq1, rp, nrp, x2, y2, z2, ss, r1, r2, i, nobs);
        } else if (nrtn == 8) {
            travderiv("reflected", nl, mll, v1, vsq1, rp, nrp, x2, y2, z2, ss, r1, r2, i, nobs);
        } else {
            fprintf(stderr, "TRAVELTIME>>> illegal nrtn from raytracer!!!\n");
            exit(1);
        }

        static double dtdr_s[3], dtdv_s[INLTOT], ttt_s;
        if (sphase[io][ie] == 2.0f && do_s == 1) {
            do_s = 0;
            ttt_s = ttt;
            for (int ii = 0; ii < 3; ++ii) dtdr_s[ii] = dtdr[ii];
            for (int ii = 0; ii < nl; ++ii) dtdv_s[ii] = dtdv[ii];
            if (nsp == 2) k2 = k2 - 1;
            goto again;
        }
        if (sphase[io][ie] == 2.0f) {
            for (int ii = 0; ii < 3; ++ii) dtdr[ii] = dtdr_s[ii] - dtdr[ii];
            ttt = ttt_s - ttt;
        }

        double pobs = (sphase[io][ie] != 2.0f) ? (pt[io][ie] - e[0][ie]) : pt[io][ie];
        double extrat1 = ptcor[k1];
        double extrat2 = 0.0f;

        extern int nshcor;
        if (nsp == 2 && sphase[io][ie] == 1.0f) extrat1 = stcor[k1];
        if (nsp == 2 && sphase[io][ie] == 2.0f) extrat1 = stcor[k1] - ptcor[k1];
        if (nsp == 3 && sphase[io][ie] == 1.0f) extrat1 = ptcor[k1] * vpvs;
        if (nsp == 3 && sphase[io][ie] == 2.0f) extrat1 = ptcor[k1] * vpvs - ptcor[k1];

        if (nshcor != 0 && i > neqs) {
            int kk1 = map2[i - neqs - 1];
            if (kk1 != 0) {
                for (int j = 0; j < IST; ++j) {
                    if (kk1 == map1[j]) {
                        extrat2 = ptcor[j];
                        if (nsp == 2 && sphase[io][ie] == 1.0f) extrat2 = stcor[j];
                        if (nsp == 3 && sphase[io][ie] == 1.0f) extrat2 = ptcor[j] * vpvs;
                        break;
                    }
                }
            }
        }

        res[io][ie] = pobs - (ttt + extrat1 + extrat2);

        if (isingle != 0) {
            tcalc[io] = ttt;
            double num = sqrtf((rp[0][1] - rp[0][0]) * (rp[0][1] - rp[0][0]) +
                              (rp[1][1] - rp[1][0]) * (rp[1][1] - rp[1][0]));
            double den = sqrtf(num * num + (rp[2][1] - rp[2][0]) * (rp[2][1] - rp[2][0]));
            double takeoff = (den > 0.0f) ? (57.296f * asinf(num / den)) : 0.0f;
            if ((rp[2][1] - rp[2][0]) < 0.0f) takeoff = 180.0f - takeoff;
            iain[io] = (int)lroundf(takeoff);
        }

        if (isingle == 0) {
            resisave(nrp, nrpdeep, rp, nobs, i, k1, mll);
        }

        tctime[io][ie] = ttt;
        h[0] = x[k1][2] - tkh;
        thk[0] = thk[0] + tkh;
    }
}
