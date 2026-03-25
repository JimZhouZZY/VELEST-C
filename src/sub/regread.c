#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern int trimlen(const char *s);

extern int indfe[];
extern int lt50[];
extern int lt25[];
extern char irname[];
extern double xnrlon[];

static void read_region_line(FILE *fp, int *nr, char *region, size_t region_size) {
    char line[256];
    if (fgets(line, sizeof(line), fp) == NULL) {
        fprintf(stderr, "REGREAD: unexpected EOF while reading region names\n");
        exit(1);
    }

    char temp[200] = {0};
    if (sscanf(line, "%d  %199[^\n]", nr, temp) < 1) {
        fprintf(stderr, "REGREAD: invalid region line: %s\n", line);
        exit(1);
    }

    strncpy(region, temp, region_size - 1);
    region[region_size - 1] = '\0';
}

void regread(const char *regnamfile, const char *regkoordfile) {
    int nstart = 0;
    int nr;
    int k;
    char cregion[64];

    FILE *fp_names = fopen(regnamfile, "r");
    if (fp_names == NULL) {
        fprintf(stderr, "REGREAD: cannot open %s\n", regnamfile);
        exit(1);
    }

    for (int i = 0; i < 729; ++i) {
        indfe[i] = nstart;
        read_region_line(fp_names, &nr, cregion, sizeof(cregion));
        if (nr != (i + 1)) {
            fprintf(stderr, "REGREAD sequence error at %d (%s), expected %d\n", nr, cregion, i + 1);
            exit(1);
        }
        int nc = trimlen(cregion);
        memcpy(irname + nstart, cregion, (size_t)nc);
        nstart += nc;
    }
    indfe[729] = nstart;

    for (int i = 0; i < 420; ++i) {
        lt25[i] = nstart;
        read_region_line(fp_names, &nr, cregion, sizeof(cregion));
        if (i <= 399) {
            k = (i + 1) + 999;
        } else {
            k = 1019 + (i - 400) * 20;
        }
        if (k != nr) {
            fprintf(stderr, "REGREAD sequence error at %d (%s), expected %d\n", nr, cregion, k);
            exit(1);
        }
        int nc = trimlen(cregion);
        memcpy(irname + nstart, cregion, (size_t)nc);
        nstart += nc;
    }
    lt25[420] = nstart;

    for (int i = 0; i < 110; ++i) {
        lt50[i] = nstart;
        read_region_line(fp_names, &nr, cregion, sizeof(cregion));
        if (i <= 99) {
            k = (i + 1) + 199;
        } else {
            k = 209 + (i - 100) * 10;
        }
        if (k != nr) {
            fprintf(stderr, "REGREAD sequence error at %d (%s), expected %d\n", nr, cregion, k);
            exit(1);
        }
        int nc = trimlen(cregion);
        memcpy(irname + nstart, cregion, (size_t)nc);
        nstart += nc;
    }
    lt50[110] = nstart;

    fclose(fp_names);

    FILE *fp_coords = fopen(regkoordfile, "r");
    if (fp_coords == NULL) {
        fprintf(stderr, "REGREAD: cannot open %s\n", regkoordfile);
        exit(1);
    }

    for (int j = 0; j < 14400; j += 8) {
        if (fscanf(fp_coords, "%lf%lf%lf%lf%lf%lf%lf%lf",
                   &xnrlon[j], &xnrlon[j + 1], &xnrlon[j + 2], &xnrlon[j + 3],
                   &xnrlon[j + 4], &xnrlon[j + 5], &xnrlon[j + 6], &xnrlon[j + 7]) != 8) {
            fprintf(stderr, "REGREAD: error reading region coordinate table\n");
            exit(1);
        }
    }

    fclose(fp_coords);
}
