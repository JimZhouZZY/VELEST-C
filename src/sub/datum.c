#include <stddef.h>

int trimlen(const char *t) {
    int len;

    if (t == NULL) {
        return 1;
    }

    for (len = 0; t[len]; ++len) {
    }

    while (len > 0 && t[len - 1] == ' ') {
        --len;
    }

    return len > 0 ? len : 1;
}

void datum(int itf, int *iyr, int *imo, int *idy, int *ihr, int *imn) {
    int kmo[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    int k, kh, id, l, iyr4, iyrh, iyrt, ld, i;

    if (iyr == NULL || imo == NULL || idy == NULL || ihr == NULL ||
        imn == NULL) {
        return;
    }

    k = itf / 60;
    *imn = itf - k * 60;
    kh = k / 24;
    *ihr = k - kh * 24;
    *iyr = kh / 365;

    while (1) {
        id = kh - *iyr * 365;
        l = 0;
        iyr4 = *iyr / 4;
        iyrh = *iyr / 100;
        iyrt = *iyr / 1000;
        ld = iyr4 - iyrh + iyrt;

        if (iyr4 * 4 == *iyr && (iyrh * 100 != *iyr || iyrt * 1000 == *iyr)) {
            l = 1;
        }

        id = id - ld + l;

        if (id > 0) {
            break;
        }

        if (id == 0 && *ihr == 0 && *imn == 0) {
            *idy = 0;
            *imo = 0;
            return;
        }

        *iyr = *iyr - 1;
    }

    kmo[1] = 28 + l;

    for (i = 0; i < 12; ++i) {
        id = id - kmo[i];
        if (id <= 0) {
            break;
        }
    }

    if (i == 12) {
        i = 11;
    }

    *idy = id + kmo[i];
    *imo = i + 1;
}
