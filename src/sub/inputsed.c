#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void timeclear(int *iyr, int *imo, int *iday, int *ihr, int *imin,
                      double *sec, int *itime);

void inputsed(FILE *fp, int *nobs,
              char sta[][5], int *iyr, int *imo, int *iday, int *ihr, int *imin, double sec[],
              char rmk1[], char rmk2[], char cphase[], int iwt[], double amx[], double prx[],
              double *xlat, double *xlon, double *xmagni, double *depth, double *origtime,
              int *iswt, int *ifixsolution, int *ievnr, char *eventtype) {
    (void)xmagni;
    const int maxobsperevent = 180;
    char cline[128];

    int iyr1[180] = {0}, imo1[180] = {0}, iday1[180] = {0}, ihr1[180] = {0}, kmin1[180] = {0};
    int itime1[180] = {0};

    *ifixsolution = 0;
    *nobs = 0;

    int j = 0;

    while (1) {
        if (fgets(cline, sizeof(cline), fp) == NULL) {
            *nobs = -1;
            return;
        }

        if (strstr(cline, "INST") || strstr(cline, "SED") || strstr(cline, "BOL") || strstr(cline, "EVENT")) {
            continue;
        }

        j++;

        if (strncmp(cline, "    ", 4) == 0) {
            if (sscanf(cline + 16, "%1d%1d%lf%lf%lf", iswt, ifixsolution, depth, xlat, xlon) < 5) {
                *nobs = -1;
                return;
            }
            *nobs = j - 1;
            break;
        }

        if (j > maxobsperevent - 1) {
            *nobs = maxobsperevent - 1;
            break;
        }

        char clayP = ' ';
        char clayS = ' ';

        if (sscanf(cline,
                   "%4s%1c%1c%1c%1d%1c%2d%2d%2d%2d%2d%lf%*7c%lf%*2c%1c%1d%*4c%lf%lf",
                   sta[j - 1], &rmk1[j - 1], &cphase[j - 1], &rmk2[j - 1], &iwt[j - 1], &clayP,
                   &iyr1[j - 1], &imo1[j - 1], &iday1[j - 1], &ihr1[j - 1], &kmin1[j - 1], &sec[j - 1],
                   &sec[j], &clayS, &iwt[j], &amx[j - 1], &prx[j - 1]) < 17) {
            fprintf(stderr, "INPUTSED read-error! input-line is:\n%s\n", cline);
            exit(1);
        }

        if (j == 1) {
            *eventtype = cline[63];
            char tmp[5] = {0};
            strncpy(tmp, cline + 75, 4);
            *ievnr = atoi(tmp);
        }

        if (iwt[j - 1] == 9) {
            j--;
            continue;
        }

        if (clayP == 'm' || clayP == 'M') {
            cphase[j - 1] = 'M';
        }

        if (strncmp(cline + 19, "     ", 5) == 0) {
            iwt[j - 1] = 5;
        }

        if (clayS != ' ') {
            cphase[j] = 'S';
            strcpy(sta[j], sta[j - 1]);
            rmk1[j] = ' ';
            rmk2[j] = ' ';
            iyr1[j] = iyr1[j - 1];
            imo1[j] = imo1[j - 1];
            iday1[j] = iday1[j - 1];
            ihr1[j] = ihr1[j - 1];
            kmin1[j] = kmin1[j - 1];
            if (iwt[j - 1] == 5) {
                amx[j] = amx[j - 1];
                prx[j] = prx[j - 1];
            }
            j++;
        }
    }

    int jjmin1 = -1;
    for (int idx = 0; idx < *nobs; ++idx) {
        if (iwt[idx] < 5) {
            jjmin1 = idx;
            break;
        }
    }
    if (jjmin1 < 0) {
        return;
    }

    int itime = 0;
    timeclear(&iyr1[jjmin1], &imo1[jjmin1], &iday1[jjmin1], &ihr1[jjmin1], &kmin1[jjmin1], &sec[jjmin1], &itime1[jjmin1]);
    itime = itime1[jjmin1];
    int jjmin = jjmin1;

    for (int idx = 0; idx < *nobs; ++idx) {
        if (iwt[idx] < 5) {
            timeclear(&iyr1[idx], &imo1[idx], &iday1[idx], &ihr1[idx], &kmin1[idx], &sec[idx], &itime1[idx]);
            if (itime1[idx] < itime) {
                itime = itime1[idx];
                jjmin = idx;
            }
        }
    }

    for (int idx = 0; idx < *nobs; ++idx) {
        if (iwt[idx] < 5) {
            sec[idx] = (double)(itime1[idx] - itime) * 60.0f + sec[idx];
        }
    }

    *origtime = sec[jjmin];
    *iyr = iyr1[jjmin];
    *imo = imo1[jjmin];
    *iday = iday1[jjmin];
    *ihr = ihr1[jjmin];
    *imin = kmin1[jjmin];
}
