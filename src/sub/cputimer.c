#include <time.h>

void cputimer(double *cpusec) {
    clock_t clk;

    if (cpusec == NULL) {
        return;
    }

    clk = clock();
    *cpusec = (double)clk / (double)CLOCKS_PER_SEC;
}
