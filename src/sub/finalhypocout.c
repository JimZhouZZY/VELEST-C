#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include "../include/globals.h"

extern bool single_turbo;
extern int isingle;
extern int legs;
extern int icoordsystem;
extern int ismpout;
extern int icnvout;
extern int nreg;
extern int nobswithw0;
extern int istm[IST][IEQ];
extern int kpwt[IST][IEQ];
extern int iyr[IEQ], imo[IEQ], iday[IEQ], ihr[IEQ], imin[IEQ];
extern int knobs[IEQ], igap[IEQ];
extern double e[5][IEQ];
extern double emag[IEQ], rms[IEQ];
extern double x[IST][3];
extern double pt[IST][IEQ];
extern double sphase[IST][IEQ];
extern char smn[IST][IEQ][5];
extern char smpline[81];
extern char regionname[33];

extern void geoko(double *x, double *y, double lat, double lon, int dir);
extern void sdc(double *x, double *y, double lat, double lon, int dir);
extern void gapcalc(int i);

extern FILE *fm_cnvfile;

void finalhypocout(void) {
    double tt[IST] = {0.0};
    char phzz[IST] = {0};

    if (!single_turbo) {
        fprintf(stdout, "\n");
    }

    for (int i = 0; i < legs; ++i) {
        double xlat = 0.0;
        double xlon = 0.0;
        double xxx = 0.0;
        char cns = 'N';
        char cew = 'W';

        if (icoordsystem == 2) {
            geoko(&xlat, &xlon, -e[1][i], e[2][i], 1);
            xlon = -xlon;
            xxx = -e[1][i];
        } else {
            sdc(&xlat, &xlon, e[1][i], e[2][i], 1);
            xxx = e[1][i];
        }

        if (xlat < 0.0) {
            cns = 'S';
            xlat = -xlat;
        }
        if (xlon < 0.0) {
            cew = 'E';
            xlon = -xlon;
        }

        double sec = e[0][i];
        int nin = imin[i];
        while (sec < 0.0) {
            sec += 60.0;
            nin--;
        }
        while (sec >= 60.0) {
            sec -= 60.0;
            nin++;
        }
        if (nin < 0) {
            nin += 60;
            ihr[i]--;
        }

        gapcalc(i);

        if (!single_turbo) {
            fprintf(stdout,
                    " %3d %02d%02d%02d %02d%02d%6.2f %7.4f%c %8.4f%c %6.2f%5.2f%4d%6.3f %7.2f%7.2f%7.2f\n",
                    i + 1, iyr[i], imo[i], iday[i], ihr[i], nin, sec,
                    xlat, cns, xlon, cew, e[3][i], emag[i], knobs[i] - nobswithw0,
                    rms[i], xxx, e[2][i], e[3][i]);
        }

        if (isingle != 0) {
            fprintf(stdout,
                    " %02d%02d%02d %02d%02d%6.2f %7.4f%c %8.4f%c %6.2f%5.2f%4d%6.3f %7.2f%7.2f%6.2f\n",
                    iyr[i], imo[i], iday[i], ihr[i], nin, sec,
                    xlat, cns, xlon, cew, e[3][i], emag[i], knobs[i] - nobswithw0,
                    rms[i], xxx, e[2][i], e[3][i]);
        }

        if (icnvout != 0) {
            fprintf(fm_cnvfile,
                    "%02d%02d%02d %02d%02d %5.2f %7.4f%c %8.4f%c %7.2f%7.2f    %3d     %5.2f\n",
                    iyr[i], imo[i], iday[i], ihr[i], nin, sec,
                    xlat, cns, xlon, cew, e[3][i], emag[i], igap[i], rms[i]);

            imin[i] = nin;
            for (int j = 0; j < knobs[i] && j < IST; ++j) {
                tt[j] = pt[j][i] - e[0][i];
            }
            e[0][i] = sec;

            for (int j = 0; j < knobs[i] && j < IST; ++j) {
                phzz[j] = 'P';
                if (sphase[j][i] == 1.0f) phzz[j] = 'S';
                if (sphase[j][i] == -1.0f) phzz[j] = 'M';
                if (sphase[j][i] == 2.0f) phzz[j] = '-';
                fprintf(fm_cnvfile, "%4.4s%c%1d%6.2f", smn[j][i], phzz[j], kpwt[j][i], tt[j]);
                // 如果需要每行输出 6 个记录（对应 Fortran 的 6(a4,a1,i1,f6.2)）
                if ((j + 1) % 6 == 0 || j == knobs[i] - 1) {
                    fprintf(fm_cnvfile, "\n");
                } else {
                    fprintf(fm_cnvfile, " "); // 空格分隔
                }
            }
            fprintf(fm_cnvfile, "\n"); // Fortran 的 write(7,*)
        }

        if (isingle != 0) {
            if (!single_turbo) {
                fprintf(stdout, " Event# %3d GAP = %3d\n", i + 1, igap[i]);
            }
            fprintf(stdout, " %s   Nr.: %4d\n", regionname, nreg);
        }

        if (!single_turbo) {
            fprintf(stdout, "\n");
        }
    }
}
