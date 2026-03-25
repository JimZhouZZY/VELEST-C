#include <stddef.h>

void chtop(double xx, double yy, double *zk, const char *topo1file,
           const char *topo2file) {
    (void)xx;
    (void)yy;
    (void)topo1file;
    (void)topo2file;
    
    if (zk == NULL) {
        return;
    }

    *zk = 0.0f;
}
