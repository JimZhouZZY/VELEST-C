#include <stdbool.h>
#include <stdio.h>

int openerror(const char *subr, const char *file_name, bool single_turbo,
              FILE *logfp) {
    if (subr == NULL) {
        subr = "";
    }
    if (file_name == NULL) {
        file_name = "";
    }

    fprintf(stderr, "WARNING:\n");
    fprintf(stderr, "SUBROUTINE :%s    ERROR OPENING FILE: %s\n", subr,
            file_name);

    if (!single_turbo) {
        FILE *out = (logfp != NULL) ? logfp : stdout;
        fprintf(out, "WARNING:\n");
        fprintf(out, "SUBROUTINE :%s    ERROR OPENING FILE: %s\n", subr,
                file_name);
    }

    return -1;
}
