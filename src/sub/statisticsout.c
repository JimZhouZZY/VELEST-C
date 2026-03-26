#include <math.h>
#include <stdio.h>
#include <string.h>
#include "../include/globals.h"

extern int neqs, nltot, nsta, noheadwave;
extern int igap[IEQ], ihypoclayer[INLTOT], irefrlayer[INLTOT], irefllayer[INLTOT];
extern double e[5][IEQ], h[INLTOT], thk[INLTOT], v[INLTOT], depthsofinput[IEQ], emag[IEQ];
extern double refraylen[INLTOT], hitlay[INLTOT][3], avhraylen, avvraylen, sterr, direrr, refrerr, reflerr;
extern double stnazires[IST][8];
extern int knobs[IEQ];
extern char stn[IST][5];

void statisticsout(void) {
    int depthnri[105];
    int depthnrl[105];
    int magnr[49];
    int nobsnr[49];
    memset(depthnri, 0, sizeof(depthnri));
    memset(depthnrl, 0, sizeof(depthnrl));
    memset(magnr, 0, sizeof(magnr));
    memset(nobsnr, 0, sizeof(nobsnr));

    for (int i = 0; i < neqs; ++i) {
        for (int ii = 0; ii < nltot - 1; ++ii) {
            if (e[3][i] >= h[ii] && e[3][i] < h[ii + 1]) ihypoclayer[ii] += 1;
        }
        if (e[3][i] >= h[nltot - 1]) ihypoclayer[nltot - 1] += 1;
    }

    double rlen = 0.0;
    for (int i = 0; i < nltot; ++i) rlen += refraylen[i];

    for (int i = 0; i < nltot; ++i) {
        refraylen[i] = (rlen > 1.0e-5f) ? (100.0 * refraylen[i] / rlen) : 0.0;
        if (hitlay[i][0] >= 1.0f) {
            hitlay[i][1] /= hitlay[i][0];
            hitlay[i][2] /= hitlay[i][0];
        }
    }

    fprintf(stdout, "\n\n\n\nRAY-STATISTICS FOR  L A S T  ITERATION\n");
    fprintf(stdout, "--------------------------------------\n\n");
    fprintf(stdout, " nlay   top ..... bottom     velocity   NHYP NREF %%len  NHIT xy-km  z-km  RFLX\n\n");

    int irefr = 0, isour = 0, irefl = 0;
    for (int i = 0; i < nltot - 1; ++i) {
        int nhit = (int)lroundf(hitlay[i][0]);
        fprintf(stdout,
                "  %2d   %6.2f...%6.2f km   %5.2f km/s %4d %4d %5.1f %5d %5.1f %5.1f %4d\n",
                i + 1, h[i], h[i] + thk[i], v[i], ihypoclayer[i], irefrlayer[i], refraylen[i],
                nhit, hitlay[i][1], hitlay[i][2], irefllayer[i]);
        irefr += irefrlayer[i];
        isour += ihypoclayer[i];
        irefl += irefllayer[i];
    }
    int nhit = (int)lroundf(hitlay[nltot - 1][0]);
    fprintf(stdout,
            "  %2d   %6.2f...      km   %5.2f km/s %4d %4d %5.1f %5d %5.1f %5.1f\n",
            nltot, h[nltot - 1], v[nltot - 1], ihypoclayer[nltot - 1], irefrlayer[nltot - 1],
            refraylen[nltot - 1], nhit, hitlay[nltot - 1][1], hitlay[nltot - 1][2]);
    irefr += irefrlayer[nltot - 1];
    isour += ihypoclayer[nltot - 1];
    irefl += irefllayer[nltot - 1];

    int itotal = noheadwave + irefr + irefl;
    fprintf(stdout, "\n Total nr of events was %4d\n\n", isour);
    fprintf(stdout, " Total nr of refracted rays = %5d\n", irefr);
    fprintf(stdout, " Total nr of reflected rays = %5d\n", irefl);
    fprintf(stdout, " Total nr of   other   rays = %5d\n", noheadwave);
    fprintf(stdout, "                               ------\n");
    fprintf(stdout, " Total nr of    all    rays = %5d\n\n", itotal);

    double err = (noheadwave > 0) ? ((sterr + direrr) / noheadwave) : 0.0;
    fprintf(stdout, " Straight and direct rays : %7.2f meters\n", err);
    err = (irefr > 0) ? (refrerr / irefr) : 0.0;
    fprintf(stdout, " Refracted           rays : %7.2f meters\n", err);
    err = (irefl > 0) ? (reflerr / irefl) : 0.0;
    fprintf(stdout, " Reflected           rays : %7.2f meters\n\n", err);

    if (itotal > 0) {
        avhraylen /= itotal;
        avvraylen /= itotal;
    }
    fprintf(stdout, "ALL RAYS TOGETHER:\n");
    fprintf(stdout, "Average horizontal ray length = %6.1f km   (Hypocenter --> Station)\n", avhraylen);
    fprintf(stdout, "Average  vertical  ray length = %6.1f km   (Deepest ray-point --> Station)\n\n", avvraylen);

    fprintf(stdout, "... and some more STATISTICS \n---------------------------- \n GAP of final epicenters:\n\n");
    int mini = 361, maxi = -1, iavgap = 0;
    for (int i = 0; i < neqs; ++i) {
        if (igap[i] > maxi) maxi = igap[i];
        if (igap[i] < mini) mini = igap[i];
        iavgap += igap[i];
    }
    if (neqs > 0) iavgap = (int)lroundf((double)iavgap / (double)neqs);
    fprintf(stdout, "GAPs were between %3d and %3d\n", mini, maxi);
    fprintf(stdout, "      (average GAP was %3d)\n\n", iavgap);

    int lesseq1 = 0, mge6 = 0, nobslesseq1 = 0, nobsmge6 = 0;
    for (int i = 0; i < neqs; ++i) {
        int mag = (int)lroundf(emag[i] * 10.0);
        if (mag <= 10) { lesseq1++; nobslesseq1 += knobs[i]; }
        else if (mag < 60) { magnr[mag - 10]++; nobsnr[mag - 10] += knobs[i]; }
        else { mge6++; nobsmge6 += knobs[i]; }
    }

    const char *cstari = "*********1*********2*********3*********4*********5>";
    fprintf(stdout, "\n MAGNITUDES of INPUT-DATA:\n\nMagnitude (# of events)   ***average number of obs***\n\n");
    if (lesseq1 > 0) nobslesseq1 = (int)lroundf((double)nobslesseq1 / (double)lesseq1);
    if (nobslesseq1 > 50) nobslesseq1 = 51;
    fprintf(stdout, "MAG<= 1.0 (%3d) %.*s\n", lesseq1, (nobslesseq1 > 0) ? nobslesseq1 : 1, (nobslesseq1 > 0) ? cstari : " ");

    for (int i = 0; i < 49; ++i) {
        double xmag = (10.0 + (i + 1)) / 10.0;
        if (magnr[i] > 0) nobsnr[i] = (int)lroundf((double)nobsnr[i] / (double)magnr[i]);
        if (nobsnr[i] > 50) nobsnr[i] = 51;
        fprintf(stdout, "MAG = %3.1f (%3d) %.*s\n", xmag, magnr[i], (nobsnr[i] > 0) ? nobsnr[i] : 1, (nobsnr[i] > 0) ? cstari : " ");
    }

    if (mge6 > 0) nobsmge6 = (int)lroundf((double)nobsmge6 / (double)mge6);
    if (nobsmge6 > 50) nobsmge6 = 51;
    fprintf(stdout, "MAG>= 6.0 (%3d) %.*s\n\n\n", mge6, (nobsmge6 > 0) ? nobsmge6 : 1, (nobsmge6 > 0) ? cstari : " ");

    int mge5i = 0, mge5l = 0;
    for (int i = 0; i < neqs; ++i) {
        int idepi = (int)lroundf(depthsofinput[i]);
        if (idepi < 100) depthnri[idepi + 5]++; else mge5i++;
        int idepl = (int)lroundf(e[3][i]);
        if (idepl < 100) depthnrl[idepl + 5]++; else mge5l++;
    }

    const char *cstarl = "*********1*********2*********3*********4*********5>";
    const char *cstari2 = ".........1.........2.........3.........4.........5>";
    fprintf(stdout, " DEPTHs of INPUT-DATA and of LAST ITERATION:\n ===========================================\n\n      Depth       # of events\n\n");
    for (int i = 0; i < 105; ++i) {
        int depth = i - 5;
        int ni = depthnri[i] > 50 ? 51 : depthnri[i];
        int nl = depthnrl[i] > 50 ? 51 : depthnrl[i];
        fprintf(stdout, "DEPTH ( input ) = %4d km :  %.*s\n", depth, (ni > 0) ? ni : 1, (ni > 0) ? cstari2 : " ");
        fprintf(stdout, "DEPTH (last_IT) = %4d km :  %.*s\n", depth, (nl > 0) ? nl : 1, (nl > 0) ? cstarl : " ");
    }
    if (mge5i > 50) mge5i = 51;
    if (mge5l > 50) mge5l = 51;
    fprintf(stdout, "DEPTH ( input ) > 100. km :  %.*s\n", (mge5i > 0) ? mge5i : 1, (mge5i > 0) ? cstari2 : " ");
    fprintf(stdout, "DEPTH (last_IT) > 100. km :  %.*s\n\n\n", (mge5l > 0) ? mge5l : 1, (mge5l > 0) ? cstarl : " ");

    fprintf(stdout, "Residuals of the stations according to the azimuth:\n");
    fprintf(stdout, "(RES  = total average residual at station)\n");
    fprintf(stdout, "(RES1 = average residual of rays from 1st quadrant)\n");
    fprintf(stdout, "(RES2 = average residual of rays from 2nd quadrant)\n");
    fprintf(stdout, "(RES3 = average residual of rays from 3rd quadrant)\n");
    fprintf(stdout, "(RES4 = average residual of rays from 4th quadrant)\n\n");
    fprintf(stdout, "Stn#  Stn     RES          RES1          RES2         RES3         RES4\n");

    for (int i = 0; i < nsta; ++i) {
        double res1 = (stnazires[i][1] > 0.0) ? stnazires[i][0] / stnazires[i][1] : 0.0;
        double res2 = (stnazires[i][3] > 0.0) ? stnazires[i][2] / stnazires[i][3] : 0.0;
        double res3 = (stnazires[i][5] > 0.0) ? stnazires[i][4] / stnazires[i][5] : 0.0;
        double res4 = (stnazires[i][7] > 0.0) ? stnazires[i][6] / stnazires[i][7] : 0.0;
        double tot = stnazires[i][1] + stnazires[i][3] + stnazires[i][5] + stnazires[i][7];
        if (tot > 0.0) {
            double res0 = (stnazires[i][1] * res1 + stnazires[i][3] * res2 + stnazires[i][5] * res3 + stnazires[i][7] * res4) / tot;
            fprintf(stdout, "%3d   %4s %7.2f(%4d)%7.2f(%4d)%7.2f(%4d)%7.2f(%4d)%7.2f(%4d)\n",
                    i + 1, stn[i], res0, (int)lroundf(tot),
                    res1, (int)lroundf(stnazires[i][1]),
                    res2, (int)lroundf(stnazires[i][3]),
                    res3, (int)lroundf(stnazires[i][5]),
                    res4, (int)lroundf(stnazires[i][7]));
        } else {
            fprintf(stdout, "%3d   %4s -.-\n", i + 1, stn[i]);
        }
    }

    fprintf(stdout, "\n\n");
}
