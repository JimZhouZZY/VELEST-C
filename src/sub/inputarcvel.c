#include <stdio.h>
#include <stdlib.h>
#include <string.h>


extern void timeclear(int *iyr, int *imo, int *iday, int *ihr, int *imin,
                      double *sec, int *itime);

static int is_blank_line(const char *line) {
    while (*line) {
        if (*line != ' ' && *line != '\t' && *line != '\n' && *line != '\r') {
            return 0;  // 不是空行
        }
        line++;
    }
    return 1;  // 全是空白
}

void inputarcvel(FILE *fp, int *nobs,
                 char sta[][5], int *iyr, int *imo, int *iday, int *ihr, int *imin, double sec[],
                 char rmk1[], char rmk2[], char cphase[], int iwt[], double amx[], double prx[],
                 double *xlat, double *xlon, double *xmagni, double *depth, double *origtime,
                 int *i1, int *i2, int *i3, char *eventtype) {
    (void)rmk1;
    (void)rmk2;
    (void)amx;
    (void)prx;
    (void)i2;
    (void)i3;

    const int maxobsperevent = 180;
    char cline[128];
    double ttime[180];
    char cns;
    char cew;
    int itime;

    *eventtype = 'L';
    *nobs = -1;

    while (fgets(cline, sizeof(cline), fp) != NULL) {
        if (cline[0] == '\n' || cline[0] == '\0') {
            continue;
        }
        break;
    }

    if (feof(fp)) {
        *nobs = -1;
        return;
    }

    /*
    if (sscanf(cline,
               "%2d%2d%2d %2d%2d %lf %lf%c %lf%c %lf %lf",
               iyr, imo, iday, ihr, imin, origtime,
               xlat, &cns, xlon, &cew, depth, xmagni) < 12) {
        fprintf(stderr, "INPUTARCVEL read error, line: %s\n", cline);
        exit(1);
    }
    */
   char buf[16];

    /* YY */
    fprintf(stderr, "DEBUG: cline = '%s'\n", cline);
    memcpy(buf, cline + 0, 2);
    buf[2] = '\0';
    *iyr = atoi(buf);
    fprintf(stderr, "DEBUG: iyr = '%d'\n", *iyr);
    if (*iyr < 100) {
        if (*iyr < 50) *iyr += 2000;  // TILL 2050
        else *iyr += 1900;            // 20..99 → 1920..1999
    }

    /* MM */
    memcpy(buf, cline+2, 2);
    buf[2] = '\0';
    fprintf(stderr, "DEBUG: buf = '%s'\n", buf);
    *imo = atoi(buf);

    /* DD */
    memcpy(buf, cline + 4, 2);
    buf[2] = '\0';
    *iday = atoi(buf);

    /* HH */
    memcpy(buf, cline + 7, 2);
    buf[2] = '\0';
    *ihr = atoi(buf);

    /* MIN */
    memcpy(buf, cline + 9, 2);
    buf[2] = '\0';
    *imin = atoi(buf);

    /* SS.SS */
    memcpy(buf, cline + 12, 6);
    buf[6] = '\0';
    *origtime = atof(buf);

    /* LAT */
    sscanf(cline + 18, "%lf%c", xlat, &cns);

    /* LON */
    sscanf(cline + 27, "%lf%c", xlon, &cew);

    /* DEPTH & MAGNITUDE */
    sscanf(cline + 38, "%lf %lf", depth, xmagni);

    fprintf(stderr, "DEBUG: Read event header:\n");
    fprintf(stderr, "  Date: %04d-%02d-%02d %02d:%02d\n", *iyr, *imo, *iday, *ihr, *imin);
    fprintf(stderr, "  Origin time: %.3f\n", *origtime);
    fprintf(stderr, "  Latitude: %.4f %c\n", *xlat, cns);
    fprintf(stderr, "  Longitude: %.4f %c\n", *xlon, cew);
    fprintf(stderr, "  Depth: %.3f\n", *depth);
    fprintf(stderr,"  Magnitude: %.3f\n", *xmagni);

    *i1 = 0;

if (cns == 'S') {
    *xlat = -*xlat;
}
if (cew == 'E') {
    *xlon = -*xlon;
}

    timeclear(iyr, imo, iday, ihr, imin, origtime, &itime);

    int j1 = 0;
    for (int j = 0; j < maxobsperevent; ++j) {
        char sta4[5] = {0};
        char phase = '\0';
        int wt = 0;
        double tt = 0.0;

        if (fgets(cline, sizeof(cline), fp) == NULL) {
            break;
        }

        if ( is_blank_line(cline)) {
            break;
        }

        if (sscanf(cline, "  %4s  %c   %d   %lf", sta4, &phase, &wt, &tt) < 4) {
            continue;
        }
        fprintf(stderr, "DEBUG: sta4='%s', phase='%c', wt=%d, tt=%.3f\n", sta4, phase, wt, tt);

        if (sta4[0] == '\0' || sta4[0] == '\n') {
            break;
        }

        strncpy(sta[j], sta4, 4);
        sta[j][4] = '\0';
        cphase[j] = phase;
        iwt[j] = wt;
        ttime[j] = tt;
        j1 = j + 1;
    }

    *nobs = j1;
    for (int j = 0; j < *nobs; ++j) {
        sec[j] = *origtime + ttime[j];
    }
}
