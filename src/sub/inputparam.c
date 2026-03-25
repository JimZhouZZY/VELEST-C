#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define IEQ 658
#define INSHOT 50
#define ITOTMODELS 2
#define INLTOT 100
#define IST 650

extern int ieq, inshot, itotmodels, inltot;
extern int maxobsperevent;

extern int neqs, nshot, icoordsystem, itrial, ised;
extern int isingle, iresolcalc, itopo, lowveloclay;
extern int nsp, nmod;
extern int nsinv, nshcor, nshfix, iuseelev, iusestacorr;
extern int iturbo, icnvout, istaout, ismpout;
extern int irayout, idrvout, ialeout, idspout, irflout, irfrout, iresout;
extern int ittmax, invertratio;
extern int legs, nsta, ksta, nltot;
extern int ifixsolution;

extern int nplay[INLTOT], laysum[INLTOT];
extern int model[ITOTMODELS * IST], map1[IST], map2[INSHOT];
extern int ielev[IST];
extern int ireflector, lmax;

extern double zshift, ztrial, dmax, zmininput, zmin;
extern double delmin, veladj, zadj, swtfac, vpvs;
extern double othet, xythet, zthet, vthet, stathet, rmsmin;
extern double vp[ITOTMODELS][INLTOT], hp[ITOTMODELS][INLTOT], vdamp[ITOTMODELS][INLTOT];
extern double thkp[ITOTMODELS][INLTOT];
extern double x[IST][3], xla[IST], xlo[IST], ptcor[IST], stcor[IST];
extern double rotate;

extern char reflchar;
extern char stn[IST][5];
extern char fm[81];
extern char headerline[3][81];

extern char modelfilename[81], stationfilename[81], seismofilename[81], scratchfilename[81];
extern char phasefile[81], shotfile[81], topo1file[81], topo2file[81], regnamfile[81], regkoordfile[81];
extern char velfile[81], cnvfile[81], rayfile[81], outfile[81], smpfile[81], stafile[81];
extern char drvfile[81], alefile[81], dsprfile[81], rflfile[81], rfrfile[81], resfile[81];

extern bool single_turbo;
extern int iabort;
extern double olat, olon;

extern FILE *fm_ptr;
extern FILE *fm_cnvfile;

void setorg(double orlat, double orlon, double rrotate, int ifil);
void geoko(double *x, double *y, double lat, double lon, int dir);
void sdc(double *x, double *y, double lat, double lon, int dir);
void casefold(char *cn);
int maxii(int n, const int *nx, int *imax, int *jndex);
int openerror(const char *subr, const char *file_name, bool single_turbo, FILE *logfp);

static void rstrip(char *s) {
    if (!s) return;
    size_t n = strlen(s);
    while (n > 0 && (s[n - 1] == '\n' || s[n - 1] == '\r' || isspace((unsigned char)s[n - 1]))) {
        s[--n] = '\0';
    }
}

static void casefold_char(char *c) {
    if (*c >= 'A' && *c <= 'Z') *c += 'a' - 'A';
}

static void copy_trim(char *dst, size_t cap, const char *src) {
    if (!dst || cap == 0) return;
    if (!src) {
        dst[0] = '\0';
        return;
    }
    strncpy(dst, src, cap - 1);
    dst[cap - 1] = '\0';
    rstrip(dst);
}

static bool is_blank_str(const char *s) {
    if (!s) return true;
    while (*s) {
        if (!isspace((unsigned char)*s)) return false;
        ++s;
    }
    return true;
}

static void slice_copy(char *dst, size_t cap, const char *src, int start, int len) {
    if (!dst || cap == 0) return;
    if (!src || start < 0 || len <= 0) {
        dst[0] = '\0';
        return;
    }
    size_t src_len = strlen(src);
    if ((size_t)start >= src_len) {
        dst[0] = '\0';
        return;
    }
    int i = 0;
    for (; i < len && (size_t)i + 1 < cap && (size_t)(start + i) < src_len; ++i) {
        char ch = src[start + i];
        if (ch == '\0' || ch == '\n' || ch == '\r') break;
        dst[i] = ch;
    }
    dst[i] = '\0';
}

static int parse_first_int(const char *s, int *out) {
    if (!s || !out) return -1;
    while (*s && !isdigit((unsigned char)*s) && *s != '-' && *s != '+') s++;
    if (!*s) return -1;
    *out = (int)strtol(s, NULL, 10);
    return 0;
}

void inputparam(void) {
    FILE *fp = fopen("velest.cmn", "r");
    if (!fp) {
        openerror("inputparam", "control-input-file FOR010", single_turbo, fm_ptr);
        iabort = 1;
        return;
    }

    char card[256];
    char lines[32][256];
    int nline = 0;

    while (fgets(card, sizeof(card), fp)) {
        if (card[0] == '*') continue;
        if (nline >= 32) {
            fprintf(stderr, "INPUTPARAM>>> control-file not correct!\n");
            fclose(fp);
            iabort = 1;
            return;
        }
        copy_trim(lines[nline], sizeof(lines[nline]), card);
        nline++;
    }
    fclose(fp);

    if (!(nline == 24 || nline == 32)) {
        fprintf(stderr, "INPUTPARAM>>> control-file not correct!\n");
        iabort = 1;
        return;
    }

    char titleline[256] = {0};
    copy_trim(titleline, sizeof(titleline), lines[0]);

    double olat_local = 0.0f, olon_local = 0.0f;
    if (sscanf(lines[1], "%lf %lf %d %lf %d %lf %d", &olat_local, &olon_local, &icoordsystem, &zshift, &itrial, &ztrial, &ised) < 7) {
        fprintf(stderr, "INPUTPARAM>>> failed parsing line 2\n");
        iabort = 1;
        return;
    }
    if (sscanf(lines[2], "%d %d %lf", &neqs, &nshot, &rotate) < 3) {
        fprintf(stderr, "INPUTPARAM>>> failed parsing line 3\n");
        iabort = 1;
        return;
    }
    sscanf(lines[3], "%d %d", &isingle, &iresolcalc);
    sscanf(lines[4], "%lf %d %lf %lf %lf %d", &dmax, &itopo, &zmininput, &veladj, &zadj, &lowveloclay);
    sscanf(lines[5], "%d %lf %lf %d", &nsp, &swtfac, &vpvs, &nmod);
    sscanf(lines[6], "%lf %lf %lf %lf %lf", &othet, &xythet, &zthet, &vthet, &stathet);
    sscanf(lines[7], "%d %d %d %d %d", &nsinv, &nshcor, &nshfix, &iuseelev, &iusestacorr);
    sscanf(lines[8], "%d %d %d %d", &iturbo, &icnvout, &istaout, &ismpout);
    sscanf(lines[9], "%d %d %d %d %d %d %d", &irayout, &idrvout, &ialeout, &idspout, &irflout, &irfrout, &iresout);
    sscanf(lines[10], "%lf %d %d", &delmin, &ittmax, &invertratio);

    copy_trim(modelfilename, sizeof(modelfilename), lines[11]);
    copy_trim(stationfilename, sizeof(stationfilename), lines[12]);
    copy_trim(seismofilename, sizeof(seismofilename), lines[13]);
    copy_trim(regnamfile, sizeof(regnamfile), lines[14]);
    copy_trim(regkoordfile, sizeof(regkoordfile), lines[15]);
    copy_trim(topo1file, sizeof(topo1file), lines[16]);
    copy_trim(topo2file, sizeof(topo2file), lines[17]);
    copy_trim(phasefile, sizeof(phasefile), lines[18]);
    copy_trim(shotfile, sizeof(shotfile), lines[19]);
    // Debug print for what is read
    fprintf(stderr, "DEBUG: modelfilename = '%s'\n", modelfilename);
    fprintf(stderr, "DEBUG: stationfilename = '%s'\n", stationfilename);
    fprintf(stderr, "DEBUG: seismofilename = '%s'\n", seismofilename);
    fprintf(stderr, "DEBUG: regnamfile = '%s'\n", regnamfile);
    fprintf(stderr, "DEBUG: regkoordfile = '%s'\n", regkoordfile);
    fprintf(stderr, "DEBUG: topo1file = '%s'\n", topo1file);
    fprintf(stderr, "DEBUG: topo2file = '%s'\n", topo2file);
    fprintf(stderr, "DEBUG: phasefile = '%s'\n", phasefile);
    fprintf(stderr, "DEBUG: shotfile = '%s'\n", shotfile);
    copy_trim(outfile, sizeof(outfile), lines[20]);
    if (is_blank_str(outfile)) strcpy(outfile, "vel.out");
    copy_trim(velfile, sizeof(velfile), lines[21]);
    if (is_blank_str(velfile)) strcpy(velfile, "velout.vel");
    copy_trim(cnvfile, sizeof(cnvfile), lines[22]);
    if (is_blank_str(cnvfile)) strcpy(cnvfile, "velout.cnv");
    copy_trim(stafile, sizeof(stafile), lines[23]);
    if (is_blank_str(stafile)) strcpy(stafile, "velout.sta");
    fm_cnvfile = fopen(cnvfile, "w");
    // Debug print for what is read
    fprintf(stderr, "DEBUG: outfile = '%s'\n", outfile);
    fprintf(stderr, "DEBUG: velfile = '%s'\n", velfile);
    fprintf(stderr, "DEBUG: cnvfile = '%s'\n", cnvfile);
    fprintf(stderr, "DEBUG: stafile = '%s'\n", stafile);
    if (nline == 32) {
        copy_trim(smpfile, sizeof(smpfile), lines[24]); if (is_blank_str(smpfile)) strcpy(smpfile, "velout.smp");
        copy_trim(rayfile, sizeof(rayfile), lines[25]); if (is_blank_str(rayfile)) strcpy(rayfile, "velout.ray");
        copy_trim(drvfile, sizeof(drvfile), lines[26]); if (is_blank_str(drvfile)) strcpy(drvfile, "velout.drv");
        copy_trim(alefile, sizeof(alefile), lines[27]); if (is_blank_str(alefile)) strcpy(alefile, "velout.ale");
        copy_trim(dsprfile, sizeof(dsprfile), lines[28]); if (is_blank_str(dsprfile)) strcpy(dsprfile, "velout.dspr");
        copy_trim(rflfile, sizeof(rflfile), lines[29]); if (is_blank_str(rflfile)) strcpy(rflfile, "velout.rfl");
        copy_trim(rfrfile, sizeof(rfrfile), lines[30]); if (is_blank_str(rfrfile)) strcpy(rfrfile, "velout.rfr");
        copy_trim(resfile, sizeof(resfile), lines[31]); if (is_blank_str(resfile)) strcpy(resfile, "velout.res");

        // Debug print for what is read
        fprintf(stderr, "DEBUG: smpfile = '%s'\n", smpfile);
        fprintf(stderr, "DEBUG: rayfile = '%s'\n", rayfile);
        fprintf(stderr, "DEBUG: drvfile = '%s'\n", drvfile);
        fprintf(stderr, "DEBUG: alefile = '%s'\n", alefile);
        fprintf(stderr, "DEBUG: dsprfile = '%s'\n", dsprfile);
        fprintf(stderr, "DEBUG: rflfile = '%s'\n", rflfile);
        fprintf(stderr, "DEBUG: rfrfile = '%s'\n", rfrfile);
        fprintf(stderr, "DEBUG: resfile = '%s'\n", resfile);
    }

    single_turbo = (isingle == 1 && iturbo == 1);

    if (!single_turbo) {
        fm_ptr = fopen(outfile, "w");
        if (!fm_ptr) {
            openerror("inputparam", "main-output-file", single_turbo, fm_ptr);
            iabort = 1;
            return;
        }
        fprintf(fm_ptr, "%s\n\n%s\n%s\n\n", headerline[0], headerline[1], headerline[2]);
        if (!is_blank_str(titleline)) {
            fprintf(fm_ptr, "Title of this VELEST run:\n%s\n\n", titleline);
        }
    } else {
        fm_ptr = stdout;
    }

    if (nsp < 1) nsp = 1;
    if (nsp > 3) nsp = 3;
    nmod = (nsp == 2) ? 2 : 1;

    zmin = zmininput;
    rmsmin = 0.0f;
    if (invertratio <= 0) invertratio = 999;

    if (isingle != 1) {
        legs = neqs + nshot;
    } else {
        legs = 1;
        neqs = 1;
        nshot = 0;
    }

    if (icoordsystem != 2) {
    fprintf(stderr, "DEBUG olat_local = %lf \n", olat_local);
        setorg(olat_local, olon_local, rotate, (!single_turbo) ? 16 : 0);
    }

    if (dmax == 0.0f) dmax = 150.0f;

    fp = fopen(modelfilename, "r");
    if (!fp) {
        openerror("inputparam", "model-input-file (FOR010)", single_turbo, fm_ptr);
        iabort = 1;
        return;
    }

    ireflector = 0;
    reflchar = ' ';

    char line[256];
    if (!fgets(line, sizeof(line), fp)) {
        fclose(fp);
        fprintf(stderr, "INPUTPARAM>>> empty model file\n");
        iabort = 1;
        return;
    }

    for (int im = 0; im < nmod; ++im) {
        int nl = 0;
        do {
            if (!fgets(line, sizeof(line), fp)) {
                fclose(fp);
                fprintf(stderr, "INPUTPARAM>>> unexpected EOF in model file\n");
                iabort = 1;
                return;
            }
            rstrip(line);
        } while (is_blank_str(line));

        if (parse_first_int(line, &nl) != 0 || nl <= 0 || nl > INLTOT) {
            fclose(fp);
            fprintf(stderr, "INPUTPARAM>>> invalid nplay in model file: %s\n", line);
            iabort = 1;
            return;
        }
        nplay[im] = nl;

        for (int j = 0; j < nl; ++j) {
            if (!fgets(line, sizeof(line), fp)) {
                fclose(fp);
                fprintf(stderr, "INPUTPARAM>>> unexpected EOF in velocity layers\n");
                iabort = 1;
                return;
            }
            double vpj = 0.0f, hpj = 0.0f, vdj = 0.0f;
            if (sscanf(line, "%lf %lf %lf", &vpj, &hpj, &vdj) < 3) {
                fclose(fp);
                fprintf(stderr, "INPUTPARAM>>> bad velocity layer line: %s\n", line);
                iabort = 1;
                return;
            }
            vp[im][j] = vpj;
            hp[im][j] = hpj;
            vdamp[im][j] = vdj;

            char reflch = ' ';
            size_t ln = strlen(line);
            if (ln > 29) reflch = line[29];
            if (reflch == 'm' || reflch == 'M') {
                reflchar = reflch;
                ireflector = j + 1;
            }
        }
    }
    fclose(fp);

    fp = fopen(stationfilename, "r");
    if (!fp) {
        openerror("inputparam", "station-input-file (FOR010)", single_turbo, fm_ptr);
        iabort = 1;
        return;
    }

    if (!fgets(line, sizeof(line), fp)) {
        fclose(fp);
        fprintf(stderr, "INPUTPARAM>>> empty station file\n");
        iabort = 1;
        return;
    }
    copy_trim(fm, sizeof(fm), line);

    nsta = 0;
    while (fgets(line, sizeof(line), fp)) {
        if (is_blank_str(line)) break;
        if (nsta >= IST) break;

        char st[8] = {0}, lats[16] = {0}, lons[16] = {0}, elevs[16] = {0}, modes[16] = {0}, iccs[16] = {0}, ptcs[16] = {0}, stcs[16] = {0};
        char cns = 'N', cew = 'W';

        slice_copy(st, sizeof(st), line, 0, 4);
        slice_copy(lats, sizeof(lats), line, 6, 8);
        if ((int)strlen(line) > 11) cns = line[13];
        slice_copy(lons, sizeof(lons), line, 15, 9);
        if ((int)strlen(line) > 21) cew = line[23];
        slice_copy(elevs, sizeof(elevs), line, 25, 4);
        slice_copy(modes, sizeof(modes), line, 30, 1);
        slice_copy(iccs, sizeof(iccs), line, 32, 3);
        slice_copy(ptcs, sizeof(ptcs), line, 36, 5);
        slice_copy(stcs, sizeof(stcs), line, 43, 5);

        if (is_blank_str(st)) break;

        strncpy(stn[nsta], st, 4);
        stn[nsta][4] = '\0';

        xla[nsta] = strtof(lats, NULL);
        fprintf(stderr, "DEBUG: RAW %d lat = %s, lon = %s\n", nsta, lats, lons);
        xlo[nsta] = strtof(lons, NULL);
        ielev[nsta] = (int)strtol(elevs, NULL, 10);
        fprintf(stderr, "DEBUG: RAW1 %d lat = %f, lon = %f\n", nsta, xla[nsta], xlo[nsta]);

        int mode = 1;
        int icc = (int)strtol(iccs, NULL, 10);

        ptcor[nsta] = strtof(ptcs, NULL);
        stcor[nsta] = strtof(stcs, NULL);

        casefold_char(&cns);
        casefold_char(&cew);
        fprintf(stderr, "DEBUG: cns = %c, cew = %c\n", cns, cew);
        if (cns == 's') xla[nsta] = -xla[nsta];
        if (cew == 'e') xlo[nsta] = -xlo[nsta];

        if (nsp == 2) {
            model[2 * nsta] = mode;
            model[2 * nsta + 1] = mode + 1;
        } else {
            model[nsta] = mode;
        }

        double dx = 0.0f, dy = 0.0f;
        double dz = -(double)ielev[nsta] / 1000.0f;
        // Debug print for latitude and longitude
        fprintf(stderr, "DEBUG: Station %d lat = %f, lon = %f\n", nsta, xla[nsta], xlo[nsta]);
        if (icoordsystem == 2) {
            geoko(&dx, &dy, xla[nsta], -xlo[nsta], -1);
            dx = -dx;
        } else {
            sdc(&dx, &dy, xla[nsta], xlo[nsta], -1);
        }

        x[nsta][0] = dx;
        x[nsta][1] = dy;
        x[nsta][2] = (iuseelev == 1) ? dz : 0.0f;

        if (iusestacorr == 0) {
            ptcor[nsta] = 0.0f;
            stcor[nsta] = 0.0f;
        }

        map1[nsta] = icc;
        nsta++;
    }
    fclose(fp);

    if (!is_blank_str(seismofilename)) {
        FILE *fsp = fopen(seismofilename, "r");
        if (!fsp) {
            openerror("inputparam", "seismo-input-file (FOR010)", single_turbo, fm_ptr);
            iabort = 1;
            return;
        }
        fclose(fsp);
    }

    int ml = 0, jndex = 0;
    maxii(nsta, map1, &ml, &jndex);

    ksta = 0;
    if (nsinv != 0) {
        ksta = ml - 1;
        if (ksta < 0) ksta = 0;
        if (nsp == 2) ksta = 2 * ksta + 1;
    }

    laysum[0] = 1;
    if (nmod > 1) laysum[1] = nplay[0] + 1;
    for (int i = 2; i < nmod; ++i) {
        laysum[i] = laysum[i - 1] + nplay[i - 1];
    }

    nltot = 0;
    for (int i = 0; i < nmod; ++i) nltot += nplay[i];

    for (int i = 0; i < ITOTMODELS; ++i) {
        for (int j = 0; j < INLTOT; ++j) thkp[i][j] = 0.0f;
    }
    for (int i = 0; i < nmod; ++i) {
        if (nplay[i] <= 1) continue;
        for (int j = 0; j < nplay[i] - 1; ++j) {
            thkp[i][j] = hp[i][j + 1] - hp[i][j];
        }
    }
}
