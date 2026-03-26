#include <math.h>
#include <stdio.h>
#include <string.h>
#include "../include/globals.h"

extern int iabort, ifixsolution, icoordsystem, nvar, nobswithw0, nreg, nmag, nitt, nsp;
extern int iyr[IEQ], imo[IEQ], iday[IEQ], ihr[IEQ], imin[IEQ], igap[IEQ], iain[IST], idelta[IST], kpwt[IST][IEQ];
extern int knobs[IEQ], istm[IST][IEQ];
extern double e[5][IEQ], emag[IEQ], rms[IEQ], ale[IEQ], spread, steplen, avres[IEQ];
extern double s[2 * (4 * IEQ + 50 + 100 + 650 - 1)], pt[IST][IEQ], tcalc[IST], res[IST][IEQ], w[IST][IEQ];
extern double d[IST][3][IEQ], ptcor[IST], stcor[IST], sphase[IST][IEQ], drm[MAXOBSPEREVENT][MAXOBSPEREVENT];
extern double amx[IST], prx[IST], xmagni[IST], xmagnitude, sdxmagnitude, vpvs;
extern char smn[IST][IEQ][5], prmk[IST][2], regionname[33];

extern void geoko(double *x, double *y, double lat, double lon, int dir);
extern void sdc(double *x, double *y, double lat, double lon, int dir);
extern void datetime(char *dattim, int dattim_len);

void statislout(void) {
    if ((knobs[0] - nobswithw0) < nvar && iabort == 0) {
        iabort = 1;
        fprintf(stderr, "knobs(i)-nobswithw0 < nvar !!!\nEvent cannot be located!!!\n");
    }
    if (iabort == 1) {
        fprintf(stdout, " ERROR: insufficient data to locate the quake!\n");
        return;
    }

    if (ifixsolution == 0) fprintf(stdout, "0 DATE  ORIGIN   TIME   LAT      LON     DEPTH  MAG  NO  DM GAP  RMS   ALE D-SPR\n");
    if (ifixsolution == 1) fprintf(stdout, "0 DATE  ORIGIN TIME   LAT       LON     *DEPTH* MAG  NO  DM GAP  RMS   ALE D-SPR\n");
    if (ifixsolution == 9) fprintf(stdout, "0 DATE  ORIGIN TIME  *LAT*     *LON*    *DEPTH* MAG  NO  DM GAP  RMS   ALE D-SPR\n");

    double sec = e[0][0];
    int nin = imin[0];
    if (sec < 0.0) { sec += 60.0; nin -= 1; }
    if (sec > 60.0) { sec -= 60.0; nin += 1; }
    if (nin < 0) { nin += 60; ihr[0] -= 1; }

    double xlat = 0.0, xlon = 0.0;
    if (icoordsystem == 2) {
        geoko(&xlat, &xlon, -e[1][0], e[2][0], 1);
        xlon = -xlon;
    } else {
        sdc(&xlat, &xlon, e[1][0], e[2][0], 1);
    }
    char cns = 'N', cew = 'W';
    if (xlat < 0.0) { cns = 'S'; xlat = -xlat; }
    if (xlon < 0.0) { cew = 'E'; xlon = -xlon; }

    int idmin = 999;
    double aar = 0.0;
    for (int k = 0; k < knobs[0]; ++k) {
        double delta = sqrtf((e[1][0] - d[k][0][0]) * (e[1][0] - d[k][0][0]) +
                            (e[2][0] - d[k][1][0]) * (e[2][0] - d[k][1][0]));
        idelta[k] = (int)lroundf(delta);
        if (w[k][0] > 0.0) {
            if (idelta[k] < idmin) idmin = idelta[k];
            aar += fabs(res[k][0]) * w[k][0];
        }
    }
    aar /= (double)(knobs[0] - nobswithw0);

    fprintf(stdout,
            " %02d%02d%02d %2d:%02d:%6.3f %7.4f%c%8.4f%c %7.3f  %3.1f  %2d %3d %3d%5.2f%6.2f %5.2f\n",
            iyr[0], imo[0], iday[0], ihr[0], nin, sec,
            xlat, cns, xlon, cew, e[3][0], emag[0], knobs[0] - nobswithw0, idmin, igap[0], rms[0], ale[0], spread);

    double erh = sqrtf(s[1] * s[1] + s[2] * s[2]);
    double erx = s[1], ery = s[2], erz = s[3];
    int js = 4;
    if ((rms[0] < 1.0f) && (erh <= 15.0f)) js = 3;
    if ((rms[0] < 0.60f) && (erh <= 8.0f) && (erz <= 15.0f)) js = 2;
    if ((rms[0] < 0.30f) && (erh <= 2.0f) && (erz <= 6.0f)) js = 1;
    int jd = 4;
    int no = knobs[0] - nobswithw0;
    double ofd = e[3][0];
    double tfd = 2.0f * e[3][0];
    if (ofd < 10.0) ofd = 10.0;
    if (tfd < 30.0) tfd = 30.0;
    if ((igap[0] <= 180) || ((no >= 4) && (idmin <= 100))) jd = 3;
    if ((igap[0] <= 135) || ((no >= 5) && (idmin <= tfd))) jd = 2;
    if ((igap[0] <= 90) || ((no >= 6) && (idmin <= ofd))) jd = 1;
    static const char classv[4] = {'A', 'B', 'C', 'D'};
    int jav = (js + jd + 1) / 2;

    int knobs1 = 0;
    for (int i = 0; i < knobs[0]; ++i) if (kpwt[i][0] < 5) knobs1++;

    fprintf(stdout, " %5.1f%5.1f%5.1f %c %c/%c%6.2f %2d%3d%6.2f%6.2f%3d  %3.1f  %3.1f %3d\n",
            erx, ery, erz, classv[jav - 1], classv[js - 1], classv[jd - 1], steplen, 0,
            knobs1, avres[0], aar, nmag, xmagnitude, sdxmagnitude, nitt);

    if (icoordsystem == 2) {
        fprintf(stdout, "0 %s NR:%4d %-32s CH-COORD.:%9.3f /%9.3f KM\n",
                (nreg >= 1000) ? "L+T" : "F-E", nreg, regionname, -e[1][0], e[2][0]);
    } else {
        fprintf(stdout, "0 F-E NR:%4d %-32s\n", nreg, regionname);
    }

    fprintf(stdout, "0 STN  DIST AZM AIN PRMK HRMN  P-SEC  TPOBS  TPCAL  -TSCOR  P-RES   P-WT IMP STURES\n");
    fprintf(stdout, "        AMX PRX     SRMK XMAG  S-SEC  TSOBS  TSCAL  -TSCOR  S-RES   S-WT IMP STURES\n");

    for (int k = 0; k < knobs[0]; ++k) {
        char char1 = 'P';
        char clay = ' ';
        if (sphase[k][0] == 1.0f) char1 = 'S';
        if (sphase[k][0] == 2.0f) char1 = '-';
        if (sphase[k][0] == -1.0f) { char1 = 'P'; clay = 'M'; }

        double xstn = d[k][0][0], ystn = d[k][1][0], xhyp = e[1][0], yhyp = e[2][0];
        double azi = 57.296f * atan2f(xhyp - xstn, yhyp - ystn);
        if (azi < 0.0) azi += 360.0;
        azi = fmodf(azi + 180.0, 360.0);
        int iazi = (int)lroundf(azi);

        double tobs = pt[k][0] - e[0][0];
        if (tobs < 0.0) pt[k][0] += 60.0;
        if (tobs > 60.0) pt[k][0] -= 60.0;
        tobs = pt[k][0] - e[0][0];

        if (tcalc[k] < 0.0) tcalc[k] += 60.0;
        if (tcalc[k] > 60.0) tcalc[k] -= 60.0;

        double tcorr = ptcor[istm[k][0]];
        if (nsp == 2 && sphase[k][0] == 1.0f) tcorr = stcor[istm[k][0]];
        if (nsp == 3 && sphase[k][0] == 1.0f) tcorr *= vpvs;
        if (sphase[k][0] == 2.0f) tcorr = 0.0;

        double studres = res[k][0] / sqrtf(fmaxf(1.0e-6f, 1.0f - drm[k][k]));
        if (studres > 999.0f) studres = 999.999f;

        fprintf(stdout,
                "  %4s %4d%4d%4d %c%c%c%d%c%2d%2d%7.3f%7.3f%7.3f%7.3f %7.3f %6.2f %6.4f%7.3f\n",
                smn[k][0], idelta[k], iazi, iain[k], prmk[k][0], char1, prmk[k][1],
                kpwt[k][0], clay, ihr[0], nin, pt[k][0], tobs, tcalc[k], tcorr,
                res[k][0], w[k][0], drm[k][k], studres);
    }

    char ctime[21] = {0};
    datetime(ctime, 20);
    fprintf(stdout, "  $$$   VELEST-Version ETH-11FEB92 located at: %s\n\n", ctime);
}
