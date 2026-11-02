#define _POSIX_C_SOURCE 200809L

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "performance.h"
#include "nec2c.h"

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
  clock_t cpu_started;
  performance_phase phases[PERFORMANCE_PHASES];
  size_t phase_count;
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
  struct rusage usage;
  if(getrusage(RUSAGE_SELF, &usage) != 0) return -1;
#ifdef __APPLE__
  return (long)(usage.ru_maxrss / 1024);
#else
  return (long)usage.ru_maxrss;
#endif
#else
  return -1;
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
  cpu = (double)(clock() - session.cpu_started) / (double)CLOCKS_PER_SEC;
  rss = peak_rss_kib();
  fprintf(session.document,
      "{\"schema\":\"nec2c-timing\",\"schema_version\":1,"
      "\"producer\":{\"program\":\"nec2c-n5hxr\",\"program_version\":");
  json_string(session.document, PACKAGE_VERSION);
  fputs("},\"input\":{\"deck\":", session.document);
  json_string(session.document, session.deck);
  fputs("},\"clock\":{\"wall\":\"monotonic\",\"cpu\":\"process\"},"
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
  session.cpu_started = clock();
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
