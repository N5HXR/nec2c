#define _POSIX_C_SOURCE 200809L

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "performance.h"
#include "nec2c.h"
#include "dense_backend.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef _WIN32
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>
#else
#include <process.h>
#include <windows.h>
#endif

#define PERFORMANCE_PHASES 32
#define PERFORMANCE_NAME 64

typedef struct
{
  char name[PERFORMANCE_NAME];
  double seconds;
  size_t calls;
} performance_phase;

typedef struct
{
  int active;
  int written;
  int success;
  char destination[NEC2C_PATH_MAX * 4 + 1];
  char temporary[NEC2C_PATH_MAX * 4 + 64];
  char deck[NEC2C_PATH_MAX * 4 + 1];
  FILE *document;
  double wall_started;
  double cpu_started;
  performance_phase phases[PERFORMANCE_PHASES];
  size_t phase_count;
  int dense_backend_set;
  int dense_backend_mixed;
  char dense_backend_requested[16];
  char dense_backend_used[16];
  int fft_backend_set;
  int fft_backend_mixed;
  char fft_backend_requested[16];
  char fft_backend_used[16];
  int concurrency_set;
  size_t jobs, threads, allowed_cpus;
  int jobs_automatic, threads_automatic, memory_limited;
} performance_session;

static performance_session session;

static void json_string(FILE *stream, const char *value)
{
  const unsigned char *p = (const unsigned char *)(value ? value : "");
  fputc('"', stream);
  while(*p)
  {
    unsigned char c = *p++;
    switch(c)
    {
      case '"': fputs("\\\"", stream); break;
      case '\\': fputs("\\\\", stream); break;
      case '\b': fputs("\\b", stream); break;
      case '\f': fputs("\\f", stream); break;
      case '\n': fputs("\\n", stream); break;
      case '\r': fputs("\\r", stream); break;
      case '\t': fputs("\\t", stream); break;
      default:
        if(c < 0x20) fprintf(stream, "\\u%04x", (unsigned)c);
        else fputc((int)c, stream);
    }
  }
  fputc('"', stream);
}

static void json_number(FILE *stream, double value)
{
  char buffer[96], *comma;
  if(!isfinite(value)) { fputs("null", stream); return; }
  snprintf(buffer, sizeof(buffer), "%.17g", value);
  comma = strchr(buffer, ',');
  if(comma) *comma = '.';
  fputs(buffer, stream);
}

double performance_now(void)
{
  if(!session.active) return 0.0;
#ifndef _WIN32
  struct timespec value;
  if(clock_gettime(CLOCK_MONOTONIC, &value) != 0) return 0.0;
  return (double)value.tv_sec + 1.0e-9 * (double)value.tv_nsec;
#else
  LARGE_INTEGER value, frequency;
  if(!QueryPerformanceCounter(&value) || !QueryPerformanceFrequency(&frequency))
    return 0.0;
  return (double)value.QuadPart / (double)frequency.QuadPart;
#endif
}

static long peak_rss_kib(void)
{
#ifndef _WIN32
  struct rusage usage, children;
  long value;
  if(getrusage(RUSAGE_SELF, &usage) != 0) return -1;
  value = (long)usage.ru_maxrss;
  if(getrusage(RUSAGE_CHILDREN, &children) == 0 &&
     (long)children.ru_maxrss > value) value = (long)children.ru_maxrss;
#ifdef __APPLE__
  return value / 1024;
#else
  return value;
#endif
#else
  return -1;
#endif
}

static double process_cpu_seconds(void)
{
#ifndef _WIN32
  struct rusage self, children;
  double value;
  if(getrusage(RUSAGE_SELF, &self) != 0) return 0.0;
  value = self.ru_utime.tv_sec + 1.0e-6 * self.ru_utime.tv_usec +
      self.ru_stime.tv_sec + 1.0e-6 * self.ru_stime.tv_usec;
  if(getrusage(RUSAGE_CHILDREN, &children) == 0)
    value += children.ru_utime.tv_sec + 1.0e-6 * children.ru_utime.tv_usec +
        children.ru_stime.tv_sec + 1.0e-6 * children.ru_stime.tv_usec;
  return value;
#else
  return (double)clock() / (double)CLOCKS_PER_SEC;
#endif
}

static void write_document(void)
{
  size_t index;
  double wall, cpu;
  long rss;
  if(!session.active || session.written) return;
  session.written = 1;
  wall = performance_now() - session.wall_started;
  cpu = process_cpu_seconds() - session.cpu_started;
  rss = peak_rss_kib();
  fprintf(session.document,
      "{\"schema\":\"nec2c-timing\",\"schema_version\":2,"
      "\"producer\":{\"program\":\"nec2c-n5hxr\",\"program_version\":");
  json_string(session.document, PACKAGE_VERSION);
  fputs("},\"input\":{\"deck\":", session.document);
  json_string(session.document, session.deck);
  fputs("},\"dense_backend\":{\"requested\":", session.document);
  if(session.dense_backend_set)
    json_string(session.document, session.dense_backend_requested);
  else fputs("null", session.document);
  fputs(",\"used\":", session.document);
  if(session.dense_backend_set)
    json_string(session.document,
        session.dense_backend_mixed ? "mixed" : session.dense_backend_used);
  else fputs("null", session.document);
  fputs("},\"fft_backend\":{\"requested\":", session.document);
  if(session.fft_backend_set)
    json_string(session.document, session.fft_backend_requested);
  else fputs("null", session.document);
  fputs(",\"used\":", session.document);
  if(session.fft_backend_set)
    json_string(session.document,
        session.fft_backend_mixed ? "mixed" : session.fft_backend_used);
  else fputs("null", session.document);
  fputs("},\"concurrency\":", session.document);
  if(session.concurrency_set)
    fprintf(session.document,
        "{\"jobs_requested\":\"%s\",\"threads_requested\":\"%s\","
        "\"jobs_used\":%lu,\"threads_per_job\":%lu,\"allowed_cpus\":%lu,"
        "\"memory_limited\":%s}",
        session.jobs_automatic ? "auto" : "explicit",
        session.threads_automatic ? "auto" : "explicit",
        (unsigned long)session.jobs, (unsigned long)session.threads,
        (unsigned long)session.allowed_cpus,
        session.memory_limited ? "true" : "false");
  else fputs("null", session.document);
  fputs(",\"clock\":{\"wall\":\"monotonic\",\"cpu\":\"process\"},"
        "\"completion\":{\"status\":", session.document);
  json_string(session.document, session.success ? "ok" : "error");
  fputs("},\"process\":{\"wall_seconds\":", session.document);
  json_number(session.document, wall);
  fputs(",\"cpu_seconds\":", session.document);
  json_number(session.document, cpu);
  fputs(",\"peak_rss_kib\":", session.document);
  if(rss < 0) fputs("null", session.document); else fprintf(session.document, "%ld", rss);
  fputs("},\"phases\":[", session.document);
  for(index = 0; index < session.phase_count; ++index)
  {
    performance_phase *phase = &session.phases[index];
    if(index) fputc(',', session.document);
    fputs("{\"name\":", session.document);
    json_string(session.document, phase->name);
    fputs(",\"seconds\":", session.document);
    json_number(session.document, phase->seconds);
    fprintf(session.document, ",\"calls\":%lu}", (unsigned long)phase->calls);
  }
  fputs("]}\n", session.document);
  if(fflush(session.document) != 0 || fclose(session.document) != 0 ||
     rename(session.temporary, session.destination) != 0)
  {
    fprintf(stderr, "nec2c: unable to publish timing JSON %s: %s\n",
            session.destination, strerror(errno));
    remove(session.temporary);
  }
  session.document = NULL;
}

static void write_at_exit(void) { write_document(); }

int performance_open(const char *destination, const char *deck_path)
{
  long process_id;
  if(!destination || !*destination || strcmp(destination, "-") == 0 || session.active)
    return 0;
#ifndef _WIN32
  process_id = (long)getpid();
#else
  process_id = (long)_getpid();
#endif
  if(strlen(destination) >= sizeof(session.destination) ||
     strlen(deck_path ? deck_path : "") >= sizeof(session.deck)) return 0;
  memset(&session, 0, sizeof(session));
  snprintf(session.destination, sizeof(session.destination), "%s", destination);
  snprintf(session.temporary, sizeof(session.temporary), "%s.tmp.%ld",
           destination, process_id);
  snprintf(session.deck, sizeof(session.deck), "%s", deck_path ? deck_path : "");
  session.document = fopen(session.temporary, "w");
  if(!session.document) { memset(&session, 0, sizeof(session)); return 0; }
  session.active = 1;
  session.wall_started = performance_now();
  session.cpu_started = process_cpu_seconds();
  if(atexit(write_at_exit) != 0)
  {
    fclose(session.document);
    remove(session.temporary);
    memset(&session, 0, sizeof(session));
    return 0;
  }
  return 1;
}

int performance_active(void) { return session.active; }

void performance_add(const char *name, double seconds)
{
  size_t index;
  performance_phase *phase;
  if(!session.active || session.written || !name || !*name ||
     !isfinite(seconds) || seconds < 0.0) return;
  for(index = 0; index < session.phase_count; ++index)
    if(strcmp(session.phases[index].name, name) == 0) break;
  if(index == session.phase_count)
  {
    if(session.phase_count == PERFORMANCE_PHASES) return;
    phase = &session.phases[session.phase_count++];
    snprintf(phase->name, sizeof(phase->name), "%s", name);
  }
  else phase = &session.phases[index];
  phase->seconds += seconds;
  ++phase->calls;
}

void performance_dense_backend(size_t order)
{
  const char *requested, *used;
  if(!session.active || session.written) return;
  requested = dense_backend_requested_name();
  used = dense_backend_name_for_order(order);
  if(!session.dense_backend_set)
  {
    snprintf(session.dense_backend_requested,
             sizeof(session.dense_backend_requested), "%s", requested);
    snprintf(session.dense_backend_used,
             sizeof(session.dense_backend_used), "%s", used);
    session.dense_backend_set = 1;
  }
  else if(strcmp(session.dense_backend_requested, requested) != 0 ||
          strcmp(session.dense_backend_used, used) != 0)
    session.dense_backend_mixed = 1;
}

void performance_fft_backend(const char *requested, const char *used)
{
  if(!session.active || session.written || !requested || !used) return;
  if(!session.fft_backend_set)
  {
    snprintf(session.fft_backend_requested,
             sizeof(session.fft_backend_requested), "%s", requested);
    snprintf(session.fft_backend_used,
             sizeof(session.fft_backend_used), "%s", used);
    session.fft_backend_set = 1;
  }
  else if(strcmp(session.fft_backend_requested, requested) != 0 ||
          strcmp(session.fft_backend_used, used) != 0)
    session.fft_backend_mixed = 1;
}

void performance_concurrency(size_t jobs, size_t threads, size_t allowed_cpus,
    int jobs_automatic, int threads_automatic, int memory_limited)
{
  if(!session.active || session.written) return;
  session.concurrency_set = 1;
  session.jobs = jobs; session.threads = threads;
  session.allowed_cpus = allowed_cpus;
  session.jobs_automatic = jobs_automatic;
  session.threads_automatic = threads_automatic;
  session.memory_limited = memory_limited;
}

void performance_finish(int success)
{
  if(!session.active || session.written) return;
  session.success = success != 0;
  write_document();
}

void performance_abandon(void)
{
  if(!session.active) return;
  if(session.document) fclose(session.document);
  remove(session.temporary);
  memset(&session, 0, sizeof(session));
}
