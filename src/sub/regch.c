#include <stddef.h>
#include <string.h>

void regch(double x, double y, const int *lt25, int lt25_len,
           const char *irname, int irname_len, char *place, int place_len,
           int *nreg) {
    int nx, ny, nrpu, i1, i2;
    double dx, dy, xx, yy;

    if (lt25 == NULL || irname == NULL || place == NULL || nreg == NULL) {
        return;
    }

    *nreg = 0;

    if (x < 480.0f || x >= 865.0f || y <= 62.0f || y > 302.0f) {
        return;
    }

    dx = 17.5f;
    dy = 12.0f;
    xx = x - 480.0f;
    yy = 302.0f - y;
    nx = (int)(xx / dx);
    ny = (int)(yy / dy);
    nrpu = ny * 20 + nx;
    *nreg = nrpu;

    if (nx >= 20) {
        nrpu = 400 + ny;
        *nreg = ny * 20 + 19;
    }

    if (nrpu < 0 || nrpu + 1 >= lt25_len) {
        *nreg = 0;
        return;
    }

    i1 = lt25[nrpu];
    i2 = lt25[nrpu + 1];

    if (i1 < 0 || i2 > irname_len || i1 > i2) {
        *nreg = 0;
        return;
    }

    *nreg += 1000;

    int len = i2 - i1;
    if (len > place_len - 1) {
        len = place_len - 1;
    }
    strncpy(place, irname + i1, (size_t)len);
    place[len] = '\0';
}
