#include <math.h>
#include <stdio.h>

extern int nitt;
extern int nrotheres, nrrefrres, nrreflres;
extern double avotheres, abotheres, avrefrres, abrefrres, avreflres, abreflres;
extern FILE *fm_ptr;

void avresistatist(void) {
    static int ifirstrun = 0;
    static double oldres = 0.0;

    int nrtotres;
    double abtotres = 0.0;
    double avtotres = 0.0;
    double proz = 0.0;

    FILE *logfp = (fm_ptr != NULL) ? fm_ptr : stdout;

    fprintf(logfp, "\n");

    if (nrotheres > 0) {
        avotheres = avotheres / (double)nrotheres;
        abotheres = abotheres / (double)nrotheres;
    }

    if (nrrefrres > 0) {
        avrefrres = avrefrres / (double)nrrefrres;
        abrefrres = abrefrres / (double)nrrefrres;
    }

    if (nrreflres > 0) {
        avreflres = avreflres / (double)nrreflres;
        abreflres = abreflres / (double)nrreflres;
    }

    nrtotres = nrotheres + nrrefrres + nrreflres;

    if (nrtotres > 0) {
        abtotres = (abotheres * (double)nrotheres + abrefrres * (double)nrrefrres +
                abreflres * (double)nrreflres) /
                   (double)nrtotres;
        avtotres = (avotheres * (double)nrotheres + avrefrres * (double)nrrefrres +
                avreflres * (double)nrreflres) /
                   (double)nrtotres;
    }

    fprintf(logfp, " After%3d iterations we got:\n", nitt);
    fprintf(logfp, "Average absolute & unweighted [and mean] residual of\n");
    fprintf(logfp, " %5d straight and direct rays =%9.5f [%9.5f]\n", nrotheres,
            abotheres, avotheres);
    fprintf(logfp, " %5d refracted           rays =%9.5f [%9.5f]\n", nrrefrres,
            abrefrres, avrefrres);
    fprintf(logfp, " %5d reflected           rays =%9.5f [%9.5f]\n", nrreflres,
            abreflres, avreflres);
    fprintf(logfp, "\n");

    if (ifirstrun != 10000001) {
        ifirstrun = 10000001;
        proz = 0.0;
    } else if (fabs(oldres) > 1.0e-10f) {
        proz = 100.0 * (abtotres - oldres) / oldres;
    }

    fprintf(logfp, " %5d ALL                 RAYS =%9.5f [%9.5f]     %7.2f %%\n",
            nrtotres, abtotres, avtotres, proz);

    oldres = abtotres;
    fprintf(logfp, "\n");
}
