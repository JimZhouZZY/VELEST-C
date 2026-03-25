#define IEQ 658
#define INSHOT 50
#define INLTOT 100
#define IST 650
#define ILIP (4 * IEQ + INSHOT + INLTOT + IST)
#define INVA (ILIP - 1)
#define IKVA ((INVA * ILIP) / 2)

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
