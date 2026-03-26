#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void timeclear(int *iyr, int *imo, int *iday, int *ihr, int *imin, double *sec, int *itime);
extern void casefold(char *s);

void inputsed_new(FILE *fp, int *nobs,
                  char sta[][5], int *iyr, int *imo, int *iday, int *ihr, int *imin, double sec[],
                  char rmk1[], char rmk2[], char cphase[], int iwt[], double amx[], double prx[],
                  double *xlat, double *xlon, double *xmagni, double *depth, double *origtime,
                  int *iswt, int *ifixsolution, int *ievnr, char *eventtype, int itrial) {
    (void)xmagni;

    const int maxobsperevent = 180;
    int iyr1[180] = {0}, imo1[180] = {0}, iday1[180] = {0}, ihr1[180] = {0}, kmin1[180] = {0};
    int itime1[180] = {0};
    int itime_o = 0;

    char cline[128];
    *ifixsolution = 0;
    *nobs = 0;
    int j = 0;

    int ev_year = 0, ev_month = 0, ev_day = 0, ev_hour = 0, ev_min = 0;

    while (fgets(cline, sizeof(cline), fp) != NULL) {
        if (strstr(cline, "SED") || strstr(cline, "BOL") || strstr(cline, "EVENT")) {
            continue;
        }

        if (strstr(cline, "INST")) {
            sscanf(cline + 28, "%d", iswt);
            continue;
        }

        if (strstr(cline, "TRIAL")) {
            sscanf(cline, "%lf%lf%lf%d%d%d%d%d%lf%d",
                   xlat, xlon, depth, iyr, imo, iday, ihr, imin, origtime, ifixsolution);
            *xlon = -*xlon;
            *iyr = *iyr - 1900;
            if (*ifixsolution == 4) *ifixsolution = 9;
            timeclear(iyr, imo, iday, ihr, imin, origtime, &itime_o);

            if (!fgets(cline, sizeof(cline), fp)) break;
            sscanf(cline + 2, "%d %d %d %d %d", &ev_year, &ev_month, &ev_day, &ev_hour, &ev_min);
            *eventtype = cline[54];
            sscanf(cline + 70, "%d", ievnr);
            continue;
        }

        if (strstr(cline, "SKIP")) {
            *nobs = j;
            goto process_event;
        }

        if (j >= maxobsperevent) {
            *nobs = maxobsperevent;
            goto process_event;
        }

        char phase_id[9] = {0}, amx_type[9] = {0};
        int use_flag = 0;
        int amx_flag = 0;

        sscanf(cline, "%4s%*4c%8s%1c%1c%lf%d%lf%lf%8s%d",
               sta[j], phase_id, &rmk1[j], &rmk2[j], &sec[j], &use_flag,
               &prx[j], &amx[j], amx_type, &amx_flag);

        casefold(phase_id);
        casefold(amx_type);
        casefold(&rmk1[j]);
        casefold(&rmk2[j]);

        int found = 0;
        if (strcmp(phase_id, "P") == 0) {
            found = 1; cphase[j] = 'P';
        }
        if (strcmp(phase_id, "S") == 0) {
            found = 1; cphase[j] = 'S';
        }
        if (strcmp(phase_id, "PMP") == 0) {
            found = 1; cphase[j] = 'M';
        }
        if (strcmp(phase_id, "S-P") == 0) {
            found = 1; cphase[j] = '-';
        }

        if (!found) {
            continue;
        }

        if (rmk1[j] == 'I') iwt[j] = 0;
        if (rmk1[j] == 'E') iwt[j] = 1;
        if (rmk1[j] == 'Q') iwt[j] = 2;
        if (rmk1[j] == ' ') iwt[j] = 5;
        if (use_flag == 0) iwt[j] = 4;

        iyr1[j] = ev_year;
        imo1[j] = ev_month;
        iday1[j] = ev_day;
        ihr1[j] = ev_hour;
        kmin1[j] = ev_min;

        j++;
    }

    *nobs = -1;
    return;

process_event:
    {
        int jjmin1 = -1;
        for (int i = 0; i < *nobs; ++i) {
            if (iwt[i] < 5) {
                jjmin1 = i;
                break;
            }
        }
        if (jjmin1 < 0) return;

        int itime = 0;
        timeclear(&iyr1[jjmin1], &imo1[jjmin1], &iday1[jjmin1], &ihr1[jjmin1], &kmin1[jjmin1], &sec[jjmin1], &itime1[jjmin1]);
        itime = itime1[jjmin1];
        int jjmin = jjmin1;

        for (int i = 0; i < *nobs; ++i) {
            if (iwt[i] < 5) {
                timeclear(&iyr1[i], &imo1[i], &iday1[i], &ihr1[i], &kmin1[i], &sec[i], &itime1[i]);
                if (itime1[i] < itime) {
                    itime = itime1[i];
                    jjmin = i;
                }
            }
        }

        if (itrial > 0) {
            for (int i = 0; i < *nobs; ++i) {
                if (iwt[i] < 5) sec[i] = (itime1[i] - itime) * 60.0 + sec[i];
            }
            *origtime = sec[jjmin];
            *iyr = iyr1[jjmin];
            *imo = imo1[jjmin];
            *iday = iday1[jjmin];
            *ihr = ihr1[jjmin];
            *imin = kmin1[jjmin];
        } else {
            for (int i = 0; i < *nobs; ++i) {
                if (iwt[i] < 5) {
                    sec[i] = (itime1[i] - itime_o) * 60.0 + sec[i];
                    if (sec[i] < 0.0) sec[i] = 0.0;
                }
            }
        }
    }
}
