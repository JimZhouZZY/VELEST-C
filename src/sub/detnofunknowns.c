#include <stdio.h>
#include <stdbool.h>

extern int neqs;
extern int nshot;
extern int nltot;
extern int ksta;
extern int nstaeff;
extern int icount;
extern int nsinv;
extern int nvar;
extern int nvareff;
extern int kvar;
extern bool single_turbo;
extern FILE *fm_ptr;

void detnofunknowns(void) {
    int lip;
    int lipeff;

    if (icount != 0) {
        lip = 4 * neqs + nshot + 1;
        lipeff = lip;
    } else {
        lip = 4 * neqs + nshot + nltot + ksta + 1;
        lipeff = 4 * neqs + nshot + nltot + nstaeff + 1;

        if (nsinv == 0) {
            lip = 4 * neqs + nshot + nltot + 1;
        }
    }

    nvar = lip - 1;
    nvareff = lipeff - 1;
    kvar = nvar * lip / 2;

    FILE *out = fm_ptr ? fm_ptr : stdout;
    //FILE *out = stdout;
    if (!single_turbo) {
        fprintf(out, "\n");
        fprintf(out, " Number of unknowns (for array-indexing):  nvar = %4d\n", nvar);
        fprintf(out, " Number of effective unknowns         :  nvareff = %4d\n", nvareff);
        fprintf(out,
                " Number of elements on/below main diagonal of matrix G = At*A : kvar = %7d\n",
                kvar);
        fprintf(out, "\n");
    }
}
