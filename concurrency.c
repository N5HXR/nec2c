#define _POSIX_C_SOURCE 200809L

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "concurrency.h"

#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _WIN32
#include <sched.h>
#include <unistd.h>
#endif

static size_t requested_jobs;
static size_t requested_threads;

static int parse_count(const char *value, size_t *result)
{
  char *end;
  unsigned long parsed;
  if(!value || !strcmp(value, "auto"))
  {
    *result = 0;
    return 1;
  }
  errno = 0;
  parsed = strtoul(value, &end, 10);
  if(errno || !*value || *end || parsed == 0 || parsed > 65536UL)
    return 0;
  *result = (size_t)parsed;
  return 1;
}

int concurrency_set_jobs(const char *value)
{ return parse_count(value, &requested_jobs); }

int concurrency_set_threads(const char *value)
{ return parse_count(value, &requested_threads); }

size_t concurrency_requested_jobs(void) { return requested_jobs; }
size_t concurrency_requested_threads(void) { return requested_threads; }
const char *concurrency_requested_jobs_name(void)
{ return requested_jobs ? "explicit" : "auto"; }
const char *concurrency_requested_threads_name(void)
{ return requested_threads ? "explicit" : "auto"; }

static size_t allowed_cpu_count(void)
{
#ifndef _WIN32
#ifdef CPU_COUNT
  cpu_set_t set;
  if(sched_getaffinity(0, sizeof(set), &set) == 0)
  {
    int count = CPU_COUNT(&set);
    if(count > 0) return (size_t)count;
  }
#endif
  {
    long count = sysconf(_SC_NPROCESSORS_ONLN);
    if(count > 0) return (size_t)count;
  }
#endif
  return 1;
}

static size_t available_memory(void)
{
#ifndef _WIN32
#if defined(_SC_AVPHYS_PAGES) && defined(_SC_PAGESIZE)
  long pages = sysconf(_SC_AVPHYS_PAGES);
  long page_size = sysconf(_SC_PAGESIZE);
  if(pages > 0 && page_size > 0 &&
     (unsigned long)pages <= SIZE_MAX / (unsigned long)page_size)
    return (size_t)pages * (size_t)page_size;
#endif
#endif
  return 0;
}

static size_t minimum_size(size_t a, size_t b) { return a < b ? a : b; }

concurrency_plan concurrency_plan_for(size_t frequency_count,
    size_t estimated_worker_bytes)
{
  concurrency_plan plan;
  size_t cpu_jobs, memory_jobs = SIZE_MAX, reserve, usable;
  memset(&plan, 0, sizeof(plan));
  plan.allowed_cpus = allowed_cpu_count();
  plan.available_memory_bytes = available_memory();
  plan.estimated_worker_bytes = estimated_worker_bytes;
  plan.jobs_automatic = requested_jobs == 0;
  plan.threads_automatic = requested_threads == 0;
  if(frequency_count == 0) frequency_count = 1;

  if(estimated_worker_bytes && plan.available_memory_bytes)
  {
    reserve = plan.available_memory_bytes / 5;
    if(reserve < (size_t)1024 * 1024 * 1024)
      reserve = (size_t)1024 * 1024 * 1024;
    usable = reserve < plan.available_memory_bytes ?
        plan.available_memory_bytes - reserve : plan.available_memory_bytes / 2;
    memory_jobs = usable / estimated_worker_bytes;
    if(memory_jobs == 0) memory_jobs = 1;
  }

  if(requested_threads)
    plan.threads = requested_threads;
  else
    plan.threads = 1;

  cpu_jobs = plan.allowed_cpus / plan.threads;
  if(cpu_jobs == 0) cpu_jobs = 1;
  if(requested_jobs)
    plan.jobs = requested_jobs;
  else
  {
    plan.jobs = minimum_size(frequency_count, cpu_jobs);
    plan.jobs = minimum_size(plan.jobs, memory_jobs);
    if(plan.jobs == 0) plan.jobs = 1;
  }

  if(!requested_threads)
  {
    size_t available_per_job = plan.allowed_cpus / plan.jobs;
    if(available_per_job == 0) available_per_job = 1;
    /* M2 measurements found diminishing or negative returns above four BLAS
       threads. Prefer more independent frequency workers when available. */
    plan.threads = minimum_size((size_t)4, available_per_job);
  }
  if(requested_jobs && requested_jobs > memory_jobs)
    plan.memory_limited = 1;
  else if(!requested_jobs && plan.jobs == memory_jobs && memory_jobs < frequency_count)
    plan.memory_limited = 1;
  return plan;
}

int concurrency_apply_threads(size_t threads)
{
  static const char *names[] = {
    "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"
  };
  char value[32];
  size_t index;
  if(threads == 0 || snprintf(value, sizeof(value), "%lu",
      (unsigned long)threads) >= (int)sizeof(value)) return 0;
#ifndef _WIN32
  for(index = 0; index < sizeof(names) / sizeof(names[0]); ++index)
    if(setenv(names[index], value, 1) != 0) return 0;
  return 1;
#else
  (void)names; (void)index;
  return 0;
#endif
}
