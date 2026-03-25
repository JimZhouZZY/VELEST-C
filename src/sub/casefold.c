#include <ctype.h>
#include <string.h>

void casefold(char *cn) {
    size_t i;
    size_t len;

    if (cn == NULL) {
        return;
    }

    len = strlen(cn);
    for (i = 0; i < len; ++i) {
        unsigned char c = (unsigned char)cn[i];
        if (c >= 'a' && c <= 'z') {
            cn[i] = (char)(c - 32);
        }
    }
}
