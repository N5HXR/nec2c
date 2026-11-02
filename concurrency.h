#ifndef NEC2C_CONCURRENCY_H
#define NEC2C_CONCURRENCY_H 1

#include <stddef.h>

typedef struct
{
  size_t jobs;
  size_t threads;
  size_t allowed_cpus;
  size_t available_memory_bytes;
  size_t estimated_worker_bytes;
  int jobs_automatic;
  int threads_automatic;
  int memory_limited;
} concurrency_plan;

/* Re-exec once with numerical-library thread variables established before
   shared-library constructors run. Call at the beginning of main(). */
int concurrency_bootstrap(int argc, char **argv);

/* "auto" is represented by zero. */
int concurrency_set_jobs(const char *value);
int concurrency_set_threads(const char *value);
size_t concurrency_requested_jobs(void);
size_t concurrency_requested_threads(void);
const char *concurrency_requested_jobs_name(void);
const char *concurrency_requested_threads_name(void);

/* Build a deterministic plan from the resources visible to this process.
   estimated_worker_bytes may be zero when no trustworthy estimate exists. */
concurrency_plan concurrency_plan_for(size_t frequency_count,
    size_t estimated_worker_bytes);

/* Apply the selected numerical-library thread count before the first backend
   call. This deliberately controls the common BLAS/OpenMP implementations. */
int concurrency_apply_threads(size_t threads);

#endif
