#include <stdio.h>
#include <stddef.h>

void freeunit(int *iunit) {
    int i;
    FILE *f;

    if (iunit == NULL) {
        return;
    }

    for (i = 10; i < 999; ++i) {
        f = fdopen(i, "r");
        if (f == NULL) {
            *iunit = i;
            return;
        }
    }

    *iunit = 10;
}
