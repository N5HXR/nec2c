#include "concurrency.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void require(int condition, const char *message)
{
  if(!condition)
  {
    fprintf(stderr, "concurrency_test: %s\n", message);
    exit(EXIT_FAILURE);
  }
}

int main(void)
{
  concurrency_plan plan;
  const char *value;

  require(concurrency_set_jobs("auto"), "auto jobs rejected");
  require(concurrency_set_threads("auto"), "auto threads rejected");
  require(!concurrency_set_jobs("0"), "zero jobs accepted");
  require(!concurrency_set_threads("-1"), "negative threads accepted");
  require(!concurrency_set_jobs("1x"), "trailing jobs text accepted");

  plan = concurrency_plan_for(8, 0);
  require(plan.jobs >= 1 && plan.jobs <= 8, "automatic job count out of range");
  require(plan.threads >= 1 && plan.threads <= 4,
      "automatic thread count out of range");
  require(plan.allowed_cpus >= 1, "allowed CPU count is zero");
  require(plan.jobs * plan.threads <= plan.allowed_cpus,
      "automatic plan oversubscribes affinity");
  require(plan.jobs_automatic && plan.threads_automatic,
      "automatic request not reported");

  require(concurrency_set_jobs("1"), "one job rejected");
  require(concurrency_set_threads("1"), "one thread rejected");
  plan = concurrency_plan_for(32, 1);
  require(plan.jobs == 1 && plan.threads == 1,
      "explicit serial plan was not preserved");
  require(!plan.jobs_automatic && !plan.threads_automatic,
      "explicit request reported as automatic");

  require(concurrency_set_jobs("32"), "large explicit job count rejected");
  plan = concurrency_plan_for(2, 0);
  require(plan.jobs == 2,
      "explicit jobs were not limited to the available frequencies");

  require(concurrency_apply_threads(3), "could not apply thread environment");
  value = getenv("OPENBLAS_NUM_THREADS");
  require(value && !strcmp(value, "3"), "OpenBLAS thread environment mismatch");
  value = getenv("OMP_NUM_THREADS");
  require(value && !strcmp(value, "3"), "OpenMP thread environment mismatch");
  return EXIT_SUCCESS;
}
