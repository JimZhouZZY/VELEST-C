#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define IEQ 658
#define INSHOT 50
#define ITOTMODELS 2
#define INLTOT 100
#define IST 650
#define MAXOBSPEREVENT 180
#define ILIP (4 * IEQ + INSHOT + INLTOT + IST)
#define INVA (ILIP - 1)

extern int nitt, nvar, istopflag, nittc, icount, ibackups, iabort;
extern int legs, neqs, nshot, nsinv, invertratio, ksta, nltot, iturbo, ifixsolution;
extern int irayout, idrvout, icnvout, ismpout, isingle, ialeout, idspout, irflout, irfrout, iresout;
extern int iresolcalc, ittmax, nsp, icoordsystem, nmag;
extern int knobs[IEQ], nobswithw0, istm[IST][IEQ], iyr[IEQ], imo[IEQ], iday[IEQ], ihr[IEQ], imin[IEQ];
extern int igap[IEQ], isconstrain[3], iconstrain[IEQ];
extern bool single_turbo;
extern double xythet, stathet, othet, vthet, zthet, delmin, davar1, xmsqrs1, steplen, spread;
extern double zmin, zmininput, vpvs;
extern double e[5][IEQ], pt[IST][IEQ], rms[IEQ], avres[IEQ], b[], rht[], g[];
extern double w[IST][IEQ];
extern double scale[7], ale[IEQ], sv[4], xmagni[IST], amx[IST], prx[IST], emag[IEQ];
extern double dtdr[3], gg2[MAXOBSPEREVENT][4], ggt[4][MAXOBSPEREVENT], ggg[4][MAXOBSPEREVENT];
extern double gtg[4][4], ggti[4][4], drm[MAXOBSPEREVENT][MAXOBSPEREVENT];
extern char headerline[3][81], smpline[81], regionname[33];
extern char rayfile[81], drvfile[81], cnvfile[81], smpfile[81], velfile[81];
extern char alefile[81], dsprfile[81], rflfile[81], rfrfile[81], resfile[81];
extern char regnamfile[81], regkoordfile[81];

extern void cputimer(double *cpusec);
extern void datetime(char *ctime, int ctime_len);
extern void inputparam(void);
extern void setunt(int nitt, int invertratio, int nsinv, int *icount,
                   double xythet, double stathet, double othet, double vthet, double zthet, double *scale);
extern void detnofunknowns(void);
extern void inputdata(int i);
extern void traveltime(int i, int nobs, int iresflag);
extern void setupmatrixg(int i, int l);
extern void singularvalues(int i);
extern void actualstations(void);
extern void storeg(int src_col, int dst_col);
extern void rmsdatvar(void);
extern void checksolution(int *istopflag, int *better);
extern void dampg(void);
extern void ludecp(const double *a, double *ul, int n, double *d1, double *d2, int *ier);
extern void luelmp(const double *a, const double *b, int n, double *x);
extern void fixunt(double *b, int neqs, int nshot, int nltot, int ksta, double *scale,
                   double vdamp[ITOTMODELS][INLTOT], int itotmodels, int inltot, int nplay1);
extern void steplengthdamp(double *damp);
extern void adjustmodel(double damp);
extern void steplengthcalc(void);
extern void nittoutput(double damp);
extern void backup(void);
extern void geoko(double *x, double *y, double lat, double lon, int dir);
extern void sdc(double *x, double *y, double lat, double lon, int dir);
extern void region(int ityp, double c1, double c2, char *place, int *nreg, const char *regnamfile, const char *regkoordfile);
extern void finalhypocout(void);
extern void finalstaresi(void);
extern void resolcovar(double datvar);
extern void statisticsout(void);
extern void timeclear(int *iyr, int *imo, int *iday, int *ihr, int *imin, double *sec, int *itime);
extern void magnitude(int i, int itime);
extern void statislout(void);
extern void matrtran(const double *a, int n, int m, double *at);
extern void matrmult(const double *a, int m, int p, const double *b, int p1, int n, double *c, int m1, int n1);
extern void matrinv(int n, const double *a, double *b);

extern double vdamp[ITOTMODELS][INLTOT];
extern int nplay[INLTOT];
extern int nreg;
extern double xmagnitude;

static int juliam_calc(int iyrv, int imo_v, int idy_v, int ihr_v, int imn_v) {
    static const int kmo[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
    int ky = iyrv;
    int km = (imo_v <= 0) ? 1 : imo_v;
    int kd = kmo[km - 1] + idy_v;
    int ky4 = ky / 4;
    int ky1 = ky / 100;
    int ky0 = ky / 1000;
    int kl = (ky4 - ky1 + ky0);
    int leap = (ky4 * 4 == ky && (ky1 * 100 != ky || ky0 * 1000 == ky)) ? 1 : 0;
    if (leap && km < 3) kl -= 1;
    int out = 365 * ky + kd + kl;
    out = out * 24 + ihr_v;
    out = out * 60 + imn_v;
    return out;
}

int main(void) {
    int i4enabletest = 123456789;
    (void)i4enabletest;

    int k = 0;
    int iresflag = 0;
    int better = 1;
    int ier = 0;
    int itime = 0;
    int iminold = 0;
    double damp = 0.0f;
    double cpusec = 0.0f;
    double cpmintot = 0.0f;
    double d1 = 0.0f, d2 = 0.0f;
    double xlat = 0.0f, xlon = 0.0f;
    char ctime[21] = {0};

    cputimer(&cpusec);
    datetime(ctime, 20);
    snprintf(headerline[0], 81, " >>> Start of program VELEST at %s <<<", ctime);
    fprintf(stdout, "%s\n\n", headerline[0]);
    snprintf(headerline[1], 81, "::::::::::::::::::::::::::::::::::::::");
    snprintf(headerline[2], 81, " V E L E S T  - Version : January 3, 1995");

    nitt = 0;
    nittc = 0;
    nvar = 0;
    istopflag = 0;

    inputparam();

single_event_loop:
    if (k == -1) {
        datetime(ctime, 20);
        cputimer(&cpusec);
        fprintf(stdout, " >>>   End of program VELEST at %s <<< [ CPU-sec : %7.1f]\n", ctime, cpusec);
        return 0;
    }

    nitt = 0;
    setunt(1, invertratio, nsinv, &icount, xythet, stathet, othet, vthet, zthet, scale);
    detnofunknowns();

    for (int j = 0; j < INVA * ILIP / 2; ++j) g[j] = 0.0f;
    for (int kk = 0; kk < nvar; ++kk) rht[kk] = 0.0f;

    for (int i = 1; i <= legs; ++i) {
        avres[i - 1] = 0.0f;
        rms[i - 1] = 0.0f;
        inputdata(i);
        k = knobs[i - 1];
        if (k == -1) goto single_event_loop;

        for (int l = 1; l <= k; ++l) {
            traveltime(i, l, iresflag);
            setupmatrixg(i, l);
        }

        if (isingle != 0) {
            if ((knobs[i - 1] - nobswithw0) < nvar) {
                iabort = 1;
                goto after_iterations;
            }
            singularvalues(i);
        }
    }

    actualstations();
    for (int col = 1; col <= 4; ++col) storeg(col, col);
    if (scale[5] != 0.0f) {
        int i1 = 4 * neqs + nshot + 1;
        int j1 = i1 + nltot - 1;
        int ii = 4;
        for (int col = i1; col <= j1; ++col) {
            ii++;
            storeg(col, ii);
        }
    }

    rmsdatvar();
    if (ittmax == 0) goto final_solution;
    checksolution(&istopflag, &better);

iteration_loop:
    datetime(ctime, 20);
    cputimer(&cpusec);
    nittc++;
    nitt++;

    setunt(nitt, invertratio, nsinv, &icount, xythet, stathet, othet, vthet, zthet, scale);
    detnofunknowns();
    dampg();

    for (int kk = 0; kk < nvar; ++kk) b[kk] = rht[kk];

    {
        int diag_last = (nvar > 0) ? ((nvar * (nvar + 1)) / 2 - 1) : 0;
        fprintf(stderr,
                "DEBUG LUDECP(pre): nitt=%d nvar=%d g0=%g gdiag_last=%g rht0=%g\n",
                nitt, nvar,
                (double)g[0], (double)g[diag_last], (double)rht[0]);
    }

    ludecp(g, g, nvar, &d1, &d2, &ier);
    if (ier != 0) {
        fprintf(stderr, "WARNING: error in ludecp ier=%d\n", ier);
        return 1;
    }

    luelmp(g, b, nvar, b);
    fixunt(b, neqs, nshot, nltot, ksta, scale, vdamp, ITOTMODELS, INLTOT, nplay[0]);
    steplengthdamp(&damp);
    adjustmodel(damp);
    steplengthcalc();

forward_after_adjust:
    if (ibackups < 4) setunt(nitt + 1, invertratio, nsinv, &icount, xythet, stathet, othet, vthet, zthet, scale);
    else setunt(nitt, invertratio, nsinv, &icount, xythet, stathet, othet, vthet, zthet, scale);

    detnofunknowns();
    for (int j = 0; j < INVA * ILIP / 2; ++j) g[j] = 0.0f;
    for (int kk = 0; kk < nvar; ++kk) rht[kk] = 0.0f;

    for (int i = 1; i <= legs; ++i) {
        avres[i - 1] = 0.0f;
        rms[i - 1] = 0.0f;
        int kkobs = knobs[i - 1];
        for (int l = 1; l <= kkobs; ++l) {
            traveltime(i, l, iresflag);
            setupmatrixg(i, l);
            if (isingle != 0) {
                gg2[l - 1][0] = 1.0f * w[l - 1][0];
                gg2[l - 1][1] = dtdr[0] * w[l - 1][0];
                gg2[l - 1][2] = dtdr[1] * w[l - 1][0];
                gg2[l - 1][3] = dtdr[2] * w[l - 1][0];
            }
        }

        if (isingle != 0) {
            if ((knobs[i - 1] - nobswithw0) < nvar) {
                iabort = 1;
                goto after_iterations;
            }
            singularvalues(i);
        }
    }

    for (int col = 1; col <= 4; ++col) storeg(col, col);
    if (scale[5] != 0.0f) {
        int i1 = 4 * neqs + nshot + 1;
        int j1 = i1 + nltot - 1;
        int ii = 4;
        for (int col = i1; col <= j1; ++col) {
            ii++;
            storeg(col, ii);
        }
    }

    setunt(nitt, invertratio, nsinv, &icount, xythet, stathet, othet, vthet, zthet, scale);
    detnofunknowns();
    rmsdatvar();
    checksolution(&istopflag, &better);

    if (ibackups < 4) {
        if (better) {
            if (ibackups < 4) ibackups = 0;
        } else {
            if (ibackups == 0 && !single_turbo) nittoutput(damp);
            if (ibackups < 4) {
                backup();
                ibackups++;
                steplengthcalc();
                goto forward_after_adjust;
            }
        }
    }

    if (!single_turbo) nittoutput(damp);

    if (ibackups == 4) goto final_solution;
    if (steplen > 0.0f && steplen < delmin) goto final_solution;
    if (nitt == ittmax) goto final_solution;
    if (isingle != 0 && istopflag == 1) goto final_solution;

    ibackups = 0;
    goto iteration_loop;

final_solution:
    cputimer(&cpusec);
    nittc++;

    if (iresolcalc > 0) {
        setunt(nitt + 1, invertratio, nsinv, &icount, xythet, stathet, othet, vthet, zthet, scale);
        detnofunknowns();
        dampg();
        {
            int diag_last = (nvar > 0) ? ((nvar * (nvar + 1)) / 2 - 1) : 0;
            fprintf(stderr,
                    "DEBUG LUDECP(pre,resol): nitt=%d nvar=%d g0=%g gdiag_last=%g\n",
                    nitt, nvar,
                    (double)g[0], (double)g[diag_last]);
        }
        ludecp(g, g, nvar, &d1, &d2, &ier);
        if (ier != 0) {
            fprintf(stderr, "WARNING: error in ludecp ier=%d\n", ier);
            return 1;
        }
    }

    if (isingle != 0 && icoordsystem == 2) {
        for (int i = 0; i < legs; ++i) {
            if (iyr[i] < 100) iminold = juliam_calc(iyr[i] + 1900, imo[i], iday[i], ihr[i], imin[i]);
            else iminold = juliam_calc(iyr[i], imo[i], iday[i], ihr[i], imin[i]);
            timeclear(&iyr[i], &imo[i], &iday[i], &ihr[i], &imin[i], &e[0][i], &itime);
            for (int j = 0; j < knobs[i]; ++j) pt[j][i] = pt[j][i] + (double)((iminold - itime) * 60);
            magnitude(i + 1, itime);
            emag[i] = xmagnitude;
        }
    }

    if (isingle != 0) {
        if (icoordsystem == 2) {
            region(1, -e[1][0], e[2][0], regionname, &nreg, regnamfile, regkoordfile);
        } else {
            sdc(&xlat, &xlon, e[1][0], e[2][0], 1);
            region(2, xlat, -xlon, regionname, &nreg, regnamfile, regkoordfile);
        }
    }

    finalhypocout();
    if (!single_turbo) finalstaresi();
    if (iresolcalc > 0) resolcovar(davar1);

    if (ialeout == 1) {
        FILE *fp = fopen(alefile, "a");
        if (fp) {
            fprintf(fp, "  ALE%10.4f%10.4f%10.4f        .1\n", -e[1][0], e[2][0], ale[0]);
            fclose(fp);
        }
    }
    if (idspout == 1) {
        FILE *fp = fopen(dsprfile, "a");
        if (fp) {
            fprintf(fp, "  DSP%10.4f%10.4f%10.4f        .1\n", -e[1][0], e[2][0], spread);
            fclose(fp);
        }
    }

    if (isingle == 0) statisticsout();

    cputimer(&cpusec);
    cpmintot = cpusec / 60.0f;
    (void)cpmintot;

    if (isingle == 0) {
        datetime(ctime, 20);
        fprintf(stdout, " >>>   End of program VELEST at %s <<< [ CPU-sec : %7.1f]\n", ctime, cpusec);
    }

after_iterations:
    if (isingle > 0) {
        datetime(ctime, 20);
        if (iabort != 1) {
            fprintf(stdout, " GAP = %3d   NITT = %2d   D-Spread =%5.2f\n", igap[0], nitt, spread);
        }

        matrtran(&gg2[0][0], MAXOBSPEREVENT, 4, &ggt[0][0]);
        matrmult(&ggt[0][0], 4, MAXOBSPEREVENT, &gg2[0][0], MAXOBSPEREVENT, 4, &gtg[0][0], 4, 4);
        gtg[0][0] += othet;
        gtg[1][1] += xythet;
        gtg[2][2] += xythet;
        gtg[3][3] += zthet;

        double gtg_v[16];
        double ggti_v[16];
        for (int r = 0; r < 4; ++r) for (int c = 0; c < 4; ++c) gtg_v[r * 4 + c] = gtg[r][c];
        matrinv(4, gtg_v, ggti_v);
        for (int r = 0; r < 4; ++r) for (int c = 0; c < 4; ++c) ggti[r][c] = ggti_v[r * 4 + c];

        matrmult(&ggti[0][0], 4, 4, &ggt[0][0], 4, MAXOBSPEREVENT, &ggg[0][0], 4, MAXOBSPEREVENT);
        matrmult(&gg2[0][0], MAXOBSPEREVENT, 4, &ggg[0][0], 4, MAXOBSPEREVENT, &drm[0][0], MAXOBSPEREVENT, MAXOBSPEREVENT);

        statislout();

        iabort = 0;
        isingle = isingle + 1;
        nitt = 0;
        nvar = 0;
        icount = 0;
        istopflag = 0;
        iresflag = 0;
        nittc = 0;
        ibackups = 0;
        ale[0] = 0.0f;
        for (int ii = 0; ii < 3; ++ii) isconstrain[ii] = 0;
        iconstrain[0] = 0;
        emag[0] = 0.0f;
        nmag = 0;
        for (int ii = 0; ii < knobs[0]; ++ii) {
            xmagni[ii] = 0.0f;
            amx[ii] = 0.0f;
            prx[ii] = 0.0f;
            istm[ii][0] = 0;
        }
        ifixsolution = 0;
        zmin = zmininput;
        goto single_event_loop;
    }

    return 0;
}
