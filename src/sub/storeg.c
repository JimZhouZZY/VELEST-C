#define IEQ 658
#define INSHOT 50
#define INLTOT 100
#define IST 650
#define ILIP (4 * IEQ + INSHOT + INLTOT + IST)
#define INVA (ILIP - 1)

extern int nvar;
extern double g[];
extern double gcopy[60][INVA * 2];

void storeg(int k, int l) {
    int ii = 0;
    int jj = 0;
    int k0 = k - 1;
    int l0 = l - 1;

    if (k <= 0 || l <= 0 || l > 60 || nvar <= 0) {
        return;
    }

    for (int i = 0; i < nvar; ++i) {
        for (int j = 0; j <= i; ++j) {
            jj++;
            if (i != k0 && j != k0) {
                continue;
            }
            ii++;
            gcopy[l0][ii - 1] = g[jj - 1];
        }
    }
}
