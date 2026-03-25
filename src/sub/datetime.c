#include <stdio.h>
#include <string.h>
#include <time.h>

void datetime(char *dattim, int dattim_len) {
    time_t current_time;
    struct tm *time_info;
    char result[30];

    if (dattim == NULL || dattim_len <= 0) {
        return;
    }

    current_time = time(NULL);
    time_info = localtime(&current_time);

    strftime(result, sizeof(result), "%a %b %d %H:%M:%S", time_info);
    strncpy(dattim, result, (size_t)(dattim_len - 1));
    dattim[dattim_len - 1] = '\0';
}
