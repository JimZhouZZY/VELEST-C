#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

extern void timeclear(int *iyr, int *imo, int *iday, int *ihr, int *imin,
                      double *sec, int *itime);

static int is_blank_line(const char *s) {
    if (s == NULL) return 1;
    for (const unsigned char *p = (const unsigned char *)s; *p; ++p) {
        if (!isspace(*p)) return 0;
    }
    return 1;
}

static void slice_copy(char *dst, size_t dst_size, const char *src, int start, int len) {
    if (dst == NULL || dst_size == 0) return;
    int i = 0;
    for (; i < len && (size_t)i + 1 < dst_size && src[start + i] != '\0'; ++i) {
        dst[i] = src[start + i];
    }
    dst[i] = '\0';
}

static void rstrip_inplace(char *s) {
    if (s == NULL) return;
    size_t n = strlen(s);
    while (n > 0 && isspace((unsigned char)s[n - 1])) {
        s[--n] = '\0';
    }
}

void inputcnv(FILE *fp, int *nobs,
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
    char cline[256];
    double ttime[180] = {0.0};
    char cns = 'N';
    char cew = 'W';
    int itime;

    *eventtype = 'L';
    int j2 = -1;
    *nobs = 0;

read_event_header:
    while (fgets(cline, sizeof(cline), fp) != NULL) {
        if (is_blank_line(cline)) continue;
        if (strncmp(cline, "9999", 4) == 0) {
            *nobs = -1;
            return;
        }
        break;
    }

    if (feof(fp)) {
        *nobs = -1;
        return;
    }

    if (sscanf(cline,
               "%2d%2d%2d %2d%2d %5lf %7lf%c %8lf%c %7lf %5lf",
               iyr, imo, iday, ihr, imin, origtime,
               xlat, &cns, xlon, &cew, depth, xmagni) < 12) {
        fprintf(stderr, "INPUTCNV>>> read-error! input-line is:\n");
        fprintf(stderr, "%s", cline);
        exit(1);
    }

    *i1 = 0;
    if (cns == 'S') *xlat = -*xlat;
    if (cew == 'E') *xlon = -*xlon;

    timeclear(iyr, imo, iday, ihr, imin, origtime, &itime);

    while (1) {
        int j1 = j2 + 1;
        j2 = j1 + 5;

        if (fgets(cline, sizeof(cline), fp) == NULL) {
            break;
        }
        for (int k = 0; k < 6; ++k) {
            int j = j1 + k;
            if (j >= maxobsperevent) break;

            int off = 12 * k;
            char s4[8] = {0}, ph[4] = {0}, wt[8] = {0}, tt[16] = {0};
            slice_copy(s4, sizeof(s4), cline, off + 0, 4);
            slice_copy(ph, sizeof(ph), cline, off + 4, 1);
            slice_copy(wt, sizeof(wt), cline, off + 5, 1);
            slice_copy(tt, sizeof(tt), cline, off + 6, 6);

            rstrip_inplace(s4);

            if (s4[0] == '\0') {
                if (j == 0 && j1 == 0) {
                    goto read_event_header;
                }
                j2 = j - 1;
                goto done_obs;
            }

            strncpy(sta[j], s4, 4);
            sta[j][4] = '\0';
            cphase[j] = (ph[0] == '\0') ? ' ' : ph[0];
            iwt[j] = atoi(wt);
            ttime[j] = strtof(tt, NULL);
        }
    }

done_obs:
    *nobs = j2 + 1;
    if (*nobs < 0) {
        *nobs = -1;
        return;
    }

    for (int j = 0; j < *nobs; ++j) {
        sec[j] = *origtime + ttime[j];
    }
}
