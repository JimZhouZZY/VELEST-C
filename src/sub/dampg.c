#include "../include/globals.h"

extern int nvar;
extern double othet;
extern double g[IKVA * 2];

void dampg(void) {
    int j = 0;
    for (int k = 1; k <= nvar; ++k) {
        j += k;
        g[j - 1] += othet;
    }
}
