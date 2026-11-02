#ifndef PERFORMANCE_H
#define PERFORMANCE_H 1

#include <stddef.h>

/* Optional process-level timing artifact. Calls are cheap no-ops when the
   command line did not request a timing document. */
int performance_open(const char *destination, const char *deck_path);
int performance_active(void);
double performance_now(void);
void performance_add(const char *phase, double elapsed_seconds);
void performance_dense_backend(size_t order);
void performance_finish(int success);
void performance_abandon(void);

#endif
