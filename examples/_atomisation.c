#line 0 "atomisation-cpp.c"
#line 0 "<built-in>"
#line 0 "<command-line>"
#line 1 "/usr/include/stdc-predef.h"
#line 0 "<command-line>"
#line 1 "atomisation-cpp.c"
#if _XOPEN_SOURCE < 700
  #undef _XOPEN_SOURCE
  #define _XOPEN_SOURCE 700
#endif
#if _GNU_SOURCE
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#endif



#line 1 "/home/gongminjiang/basilisk/src/common.h"
#line 1 "/home/gongminjiang/basilisk/src/ast/std/stdlib.h"
#include <stdlib.h>
#line 2 "/home/gongminjiang/basilisk/src/common.h"
#line 1 "/home/gongminjiang/basilisk/src/ast/std/stdio.h"
#include <stdio.h>
#line 3 "/home/gongminjiang/basilisk/src/common.h"
#line 1 "/home/gongminjiang/basilisk/src/ast/std/stddef.h"
#include <stddef.h>
#line 4 "/home/gongminjiang/basilisk/src/common.h"
#line 1 "/home/gongminjiang/basilisk/src/ast/std/stdbool.h"
#include <stdbool.h>
#line 5 "/home/gongminjiang/basilisk/src/common.h"
#line 1 "/home/gongminjiang/basilisk/src/ast/std/stdarg.h"
#include <stdarg.h>
#line 6 "/home/gongminjiang/basilisk/src/common.h"
#line 1 "/home/gongminjiang/basilisk/src/ast/std/string.h"
#include <string.h>
#line 7 "/home/gongminjiang/basilisk/src/common.h"
#line 1 "/home/gongminjiang/basilisk/src/ast/std/float.h"
#include <float.h>
#line 8 "/home/gongminjiang/basilisk/src/common.h"
#line 1 "/home/gongminjiang/basilisk/src/ast/std/limits.h"
#include <limits.h>
#line 9 "/home/gongminjiang/basilisk/src/common.h"
#line 1 "/home/gongminjiang/basilisk/src/ast/std/math.h"
#include <math.h>
#line 10 "/home/gongminjiang/basilisk/src/common.h"
#line 1 "/home/gongminjiang/basilisk/src/ast/std/time.h"
#include <time.h>
#line 11 "/home/gongminjiang/basilisk/src/common.h"
#line 1 "/home/gongminjiang/basilisk/src/ast/std/sys/time.h"
#include <sys/time.h>
#line 12 "/home/gongminjiang/basilisk/src/common.h"
#line 1 "/home/gongminjiang/basilisk/src/ast/std/sys/resource.h"
#include <sys/resource.h>
#line 13 "/home/gongminjiang/basilisk/src/common.h"

#if _OPENMP
# include <omp.h>
# define OMP(x) _Pragma(#x)
#elif 1

# define OMP(x)

# include <mpi.h>
static int mpi_rank, mpi_npe;
# define tid() mpi_rank
# define pid() mpi_rank
# define npe() mpi_npe

#else

# define OMP(x)

#endif
#line 46 "/home/gongminjiang/basilisk/src/common.h"
#undef HUGE
#define HUGE ((double)1e30)

#define _NVARMAX 65536
#define is_constant(v) ((v).i >= _NVARMAX)
#define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : HUGE)

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define sq(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define sign(x) ((x) > 0 ? 1 : -1)
#define noise() (1. - 2.*rand()/(double)RAND_MAX)
#define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))

#define unmap(x,y)

#define trash(x)


#define systderr stderr
#define systdout stdout

#if 1
FILE * qstderr (void);
FILE * qstdout (void);
FILE * ferr = NULL, * fout = NULL;
#define not_mpi_compatible()\
do {\
  if (npe() > 1) {\
    fprintf (ferr, "%s() is not compatible with MPI (yet)\n", __func__);\
    exit (1);\
  }\
} while(0)\

#line 80

# define system(command) (pid() == 0 ? system(command) : 0)
#else
# define qstderr() stderr
# define qstdout() stdout
# define ferr stderr
# define fout stdout
# define not_mpi_compatible()
#endif



static inline void qassert (const char * file, int line, const char * cond) {
  fprintf (ferr, "%s:%d: Assertion `%s' failed.\n", file, line, cond);
  abort();
}
#line 104 "/home/gongminjiang/basilisk/src/common.h"
#define sysmalloc malloc
#define syscalloc calloc
#define sysrealloc realloc
#define sysfree free
#define systrdup strdup

#if MTRACE

struct {
  FILE * fp;
  size_t total, max;
  size_t overhead, maxoverhead;
  size_t nr;
  size_t startrss, maxrss;
  char * fname;
} pmtrace;

typedef struct {
  char * func, * file;
  size_t max, total;
  int line, id;
} pmfunc;

typedef struct {
  size_t id, size;
} pmdata;

static pmfunc * pmfuncs = NULL;
static int pmfuncn = 0;

static int pmfunc_index (const char * func, const char * file, int line)
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++)
    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))
      return p->id;
  pmfuncn++;
  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));
  p = &pmfuncs[pmfuncn - 1];
  memset (p, 0, sizeof(pmfunc));
  p->func = systrdup(func);
  p->file = systrdup(file);
  p->line = line;
  p->id = pmfuncn;
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "@ %d %s %s %d\n", pmfuncn, func, file, line);
  return pmfuncn;
}

static void pmfunc_trace (pmfunc * f, char c)
{
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "%c %d %ld %ld %ld",
      c, f->id, pmtrace.nr, pmtrace.total, f->total);
#if _GNU_SOURCE
  if (pmtrace.nr % 1 == 0) {
    struct rusage usage;
    getrusage (RUSAGE_SELF, &usage);
    if (pmtrace.fp)
      fprintf (pmtrace.fp, " %ld", usage.ru_maxrss*1024);
    if (!pmtrace.nr)
      pmtrace.startrss = usage.ru_maxrss;
    if (usage.ru_maxrss > pmtrace.maxrss)
      pmtrace.maxrss = usage.ru_maxrss;
  }
#endif
  if (pmtrace.fp)
    fputc ('\n', pmtrace.fp);
  pmtrace.nr++;
}

static void * pmfunc_alloc (pmdata * d, size_t size,
       const char * func, const char * file, int line,
       char c)
{
  if (!(d != NULL)) qassert ("/home/gongminjiang/basilisk/src/common.h", 179, "d != NULL");
  OMP (omp critical)
  {
    d->id = pmfunc_index(func, file, line);
    d->size = size;
    pmfunc * f = &pmfuncs[d->id - 1];
    f->total += size;
    if (f->total > f->max)
      f->max = f->total;
    pmtrace.total += size;
    pmtrace.overhead += sizeof(pmdata);
    if (pmtrace.total > pmtrace.max) {
      pmtrace.max = pmtrace.total;
      pmtrace.maxoverhead = pmtrace.overhead;
    }
    pmfunc_trace (f, c);
  }
  return ((char *)d) + sizeof(pmdata);
}

static void * pmfunc_free (void * ptr, char c)
{
  if (!ptr)
    return ptr;
  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));
  if (d->id < 1 || d->id > pmfuncn) {
    fputs ("*** MTRACE: ERROR!: corrupted free()", ferr);
    if (d->size == 0)
      fputs (", possible double free()", ferr);
    else
      fputs (", not traced?", ferr);
    fputs (", aborting...\n", ferr);
    abort();
    return ptr;
  }
  else
  OMP (omp critical)
  {
    pmfunc * f = &pmfuncs[d->id - 1];
    if (f->total < d->size) {
      fprintf (ferr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        f->total, d->size);
      abort();
    }
    else
      f->total -= d->size;
    if (pmtrace.total < d->size) {
      fprintf (ferr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        pmtrace.total, d->size);
      abort();
    }
    else {
      pmtrace.total -= d->size;
      pmtrace.overhead -= sizeof(pmdata);
    }
    d->id = 0;
    d->size = 0;
    pmfunc_trace (f, c);
  }
  return d;
}

static void * pmalloc (size_t size,
         const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysmalloc (sizeof(pmdata) + size),
         size, func, file, line, '+');
}

static void * pcalloc (size_t nmemb, size_t size,
         const char * func, const char * file, int line)
{
  void * p = pmalloc (nmemb*size, func, file, line);
  return memset (p, 0, nmemb*size);
}

static void * prealloc (void * ptr, size_t size,
   const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysrealloc (pmfunc_free(ptr, '<'),
           sizeof(pmdata) + size),
         size, func, file, line, '>');
}

static void pfree (void * ptr,
     const char * func, const char * file, int line)
{
  sysfree (pmfunc_free (ptr, '-'));
}

static char * pstrdup (const char * s,
         const char * func, const char * file, int line)
{
  char * d = (char *) pmalloc (strlen(s) + 1, func, file, line);
  return strcpy (d, s);
}

#if MTRACE < 3
static int pmaxsort (const void * a, const void * b) {
  const pmfunc * p1 = a, * p2 = b;
  return p1->max < p2->max;
}
#endif

static int ptotalsort (const void * a, const void * b) {
  const pmfunc * p1 = (const pmfunc *) a, * p2 = (const pmfunc *) b;
  return p1->total < p2->total;
}

static void pmfuncs_free()
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++) {
    sysfree (p->func);
    sysfree (p->file);
  }
  sysfree (pmfuncs);
}

void pmuntrace (void)
{
#if MTRACE < 3
  fprintf (ferr,
    "*** MTRACE: max resident  set size: %10ld bytes\n"
    "*** MTRACE: max traced memory size: %10ld bytes"
    " (tracing overhead %.1g%%)\n"
    "%10s    %20s   %s\n",
    pmtrace.maxrss*1024,
    pmtrace.max, pmtrace.maxoverhead*100./pmtrace.max,
    "max bytes", "function", "file");
  qsort (pmfuncs, pmfuncn, sizeof(pmfunc), pmaxsort);
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn && p->max > 0; i++, p++)
    fprintf (ferr, "%10ld    %20s   %s:%d\n",
      p->max, p->func, p->file, p->line);

  if (pmtrace.fp) {
    char * fname = pmtrace.fname, * s;
    while ((s = strchr(fname,'/')))
      fname = s + 1;

    fputs ("load(\"`echo $BASILISK`/mtrace.plot\")\n", pmtrace.fp);
    fprintf (pmtrace.fp,
      "plot '%s' u 3:($6-%g) w l t 'ru_maxrss - %.3g',"
      "total(\"%s\") w l t 'total'",
      fname,
      pmtrace.startrss*1024.,
      pmtrace.startrss*1024.,
      fname);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->max > 0.01*pmtrace.max; i++, p++)
      fprintf (pmtrace.fp,
        ",func(\"%s\",%d) w l t '%s'",
        fname, p->id, p->func);
    fputc ('\n', pmtrace.fp);
    fprintf (ferr,
      "*** MTRACE: To get a graph use: tail -n 2 %s | gnuplot -persist\n",
      fname);
    fclose (pmtrace.fp);
    pmtrace.fp = NULL;
    sysfree (pmtrace.fname);
  }
#endif

  if (pmtrace.total > 0) {
    qsort (pmfuncs, pmfuncn, sizeof(pmfunc), ptotalsort);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->total > 0; i++, p++)
      fprintf (ferr, "%s:%d: error: %ld bytes leaked here\n",
        p->file, p->line, p->total);
    pmfuncs_free();
    exit(1);
  }
  else {
#if MTRACE < 3
    fputs ("*** MTRACE: No memory leaks\n", ferr);
#endif
    pmfuncs_free();
  }
}

#else
# define pmalloc(s,func,file,line) malloc(s)
# define pcalloc(n,s,func,file,line) calloc(n,s)
# define prealloc(p,s,func,file,line) realloc(p,s)
# define pfree(p,func,file,line) free(p)
# define pstrdup(s,func,file,line) strdup(s)
#endif







typedef struct {
  void * p;
  long max, len;
} Array;

Array * array_new()
{
  Array * a = ((Array *) pmalloc ((1)*sizeof(Array),__func__,__FILE__,__LINE__));
  a->p = NULL;
  a->max = a->len = 0;
  return a;
}

void array_free (Array * a)
{
  pfree (a->p,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
}

void array_append (Array * a, void * elem, size_t size)
{
  if (a->len + size >= a->max) {
    a->max += max (size, 4096);
    a->p = prealloc (a->p, a->max,__func__,__FILE__,__LINE__);
  }
  memcpy (((char *)a->p) + a->len, elem, size);
  a->len += size;
}

void * array_shrink (Array * a)
{
  void * p = prealloc (a->p, a->len,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
  return p;
}



#if TRACE == 1
#include <extrae_user_events.h>

typedef struct {
  Array index, stack;
  extrae_type_t type;
} Trace;

Trace trace_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000010,
};

Trace trace_mpi_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000011,
};

static int lookup_func (Array * a, const char * func)
{
  for (int i = 0; i < a->len/sizeof(char *); i++) {
    char * s = ((char **)a->p)[i];
    if (!strcmp (func, s))
      return i + 1;
  }
  char * s = pstrdup (func,__func__,__FILE__,__LINE__);
  array_append (a, &s, sizeof(char *));
  return a->len;
}

static void trace_push (Trace * t, const char * func)
{
  int value = lookup_func (&t->index, func);
  Extrae_eventandcounters (t->type, value);
  array_append (&t->stack, &value, sizeof(int));
}

static void trace_pop (Trace * t, const char * func)
{
  if (!(t->stack.len > 0)) qassert ("/home/gongminjiang/basilisk/src/common.h", 451, "t->stack.len > 0");
  t->stack.len -= sizeof(int);
  int value = t->stack.len > 0 ?
    ((int *)t->stack.p)[t->stack.len/sizeof(int) - 1] : 0;
  Extrae_eventandcounters (t->type, value);
}

static void trace_define (Trace * t, char * description)
{
  if (t->index.len > 0) {
    extrae_value_t values[t->index.len/sizeof(char *) + 1];
    char * names[t->index.len/sizeof(char *) + 1],
      ** func = (char **) t->index.p;
    names[0] = "OTHER";
    values[0] = 0;
    unsigned len = 1;
    for (int i = 0; i < t->index.len/sizeof(char *); i++, func++) {
      names[len] = *func;
      values[len++] = i + 1;
    }
    Extrae_define_event_type (&t->type, description, &len, values, names);
  }
}

static void trace_free (Trace * t)
{
  char ** func = (char **) t->index.p;
  for (int i = 0; i < t->index.len/sizeof(char *); i++, func++)
    pfree (*func,__func__,__FILE__,__LINE__);
  pfree (t->index.p,__func__,__FILE__,__LINE__);
  pfree (t->stack.p,__func__,__FILE__,__LINE__);
}

static void trace_off()
{
  trace_define (&trace_func, "Basilisk functions");
  trace_define (&trace_mpi_func, "Basilisk functions (MPI-related)");
  trace_free (&trace_func);
  trace_free (&trace_mpi_func);
}






# define tracing(func, file, line) trace_push (&trace_func, func)
# define end_tracing(func, file, line) trace_pop (&trace_func, func)

#elif TRACE

typedef struct {
  char * func, * file;
  int line, calls;
  double total, self;
#if 1
  double min, max;
#endif
} TraceIndex;

struct {
  Array stack, index;
  double t0;
} Trace = {
  {NULL, 0, 0}, {NULL, 0, 0},
  -1
};

static void trace_add (const char * func, const char * file, int line,
         double total, double self)
{
  TraceIndex * t = (TraceIndex *) Trace.index.p;
  int i, len = Trace.index.len/sizeof(TraceIndex);
  for (i = 0; i < len; i++, t++)
    if (t->line == line && !strcmp (func, t->func) && !strcmp (file, t->file))
      break;
  if (i == len) {
    TraceIndex t = {pstrdup(func,__func__,__FILE__,__LINE__), pstrdup(file,__func__,__FILE__,__LINE__), line, 1, total, self};
    array_append (&Trace.index, &t, sizeof(TraceIndex));
  }
  else
    t->calls++, t->total += total, t->self += self;
}

static void tracing (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  if (Trace.t0 < 0)
    Trace.t0 = tv.tv_sec + tv.tv_usec/1e6;
  double t[2] = {(tv.tv_sec - Trace.t0) + tv.tv_usec/1e6, 0.};
  array_append (&Trace.stack, t, 2*sizeof(double));




}

static void end_tracing (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  double te = (tv.tv_sec - Trace.t0) + tv.tv_usec/1e6;
  double * t = (double *) Trace.stack.p;
  if (!(Trace.stack.len >= 2*sizeof(double))) qassert ("/home/gongminjiang/basilisk/src/common.h", 555, "Trace.stack.len >= 2*sizeof(double)");
  t += Trace.stack.len/sizeof(double) - 2;
  Trace.stack.len -= 2*sizeof(double);
  double dt = te - t[0];




  trace_add (func, file, line, dt, dt - t[1]);
  if (Trace.stack.len >= 2*sizeof(double)) {
    t -= 2;
    t[1] += dt;
  }
}

static int compar_self (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  return t1->self < t2->self;
}

#if 1
static int compar_func (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  if (t1->line != t2->line)
    return t1->line < t2->line;
  return strcmp (t1->file, t2->file);
}
#endif

void trace_print (FILE * fp, double threshold)
{
  int i, len = Trace.index.len/sizeof(TraceIndex);
  double total = 0.;
  TraceIndex * t;
  Array * index = array_new();
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    array_append (index, t, sizeof(TraceIndex)), total += t->self;
#if 1
  qsort (index->p, len, sizeof(TraceIndex), compar_func);
  double tot[len], self[len], min[len], max[len];
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    tot[i] = t->total, self[i] = t->self;
  MPI_Reduce (self, min, len, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce (self, max, len, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? self : MPI_IN_PLACE,
       self, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? tot : MPI_IN_PLACE,
       tot, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  total = 0.;
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    t->total = tot[i]/npe(), t->self = self[i]/npe(),
      t->max = max[i], t->min = min[i], total += t->self;
#endif
  qsort (index->p, len, sizeof(TraceIndex), compar_self);
  fprintf (fp, "   calls    total     self   %% total   function\n");
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    if (t->self*100./total > threshold) {
      fprintf (fp, "%8d   %6.2f   %6.2f     %4.1f%%",
        t->calls, t->total, t->self, t->self*100./total);
#if 1
      fprintf (fp, " (%4.1f%% - %4.1f%%)", t->min*100./total, t->max*100./total);
#endif
      fprintf (fp, "   %s():%s:%d\n", t->func, t->file, t->line);
    }
  fflush (fp);
  array_free (index);
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    t->calls = t->total = t->self = 0.;
}

static void trace_off()
{
  trace_print (fout, 0.);

  int i, len = Trace.index.len/sizeof(TraceIndex);
  TraceIndex * t;
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    pfree (t->func,__func__,__FILE__,__LINE__), pfree (t->file,__func__,__FILE__,__LINE__);

  pfree (Trace.index.p,__func__,__FILE__,__LINE__);
  Trace.index.p = NULL;
  Trace.index.len = Trace.index.max = 0;

  pfree (Trace.stack.p,__func__,__FILE__,__LINE__);
  Trace.stack.p = NULL;
  Trace.stack.len = Trace.stack.max = 0;
}

#else
# define tracing(...)
# define end_tracing(...)
#endif



#if _OPENMP

#define tid() omp_get_thread_num()
#define pid() 0
#define npe() omp_get_num_threads()
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_array(v,type,op,elem)

#elif 1

static bool in_prof = false;
static double prof_start, _prof;
#define prof_start(name)\
  if (!(!in_prof)) qassert ("/home/gongminjiang/basilisk/src/common.h", 665, "!in_prof"); in_prof = true;\
  prof_start = MPI_Wtime();\

#line 667

#define prof_stop()\
  if (!(in_prof)) qassert ("/home/gongminjiang/basilisk/src/common.h", 669, "in_prof"); in_prof = false;\
  _prof = MPI_Wtime();\
  mpi_time += _prof - prof_start;\

#line 672


#if FAKE_MPI
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_array(v,type,op,elem)
#else
     
int mpi_all_reduce0 (void *sendbuf, void *recvbuf, int count,
       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{tracing("mpi_all_reduce0","/home/gongminjiang/basilisk/src/common.h",679);
  { int _ret= MPI_Allreduce (sendbuf, recvbuf, count, datatype, op, comm);end_tracing("mpi_all_reduce0","/home/gongminjiang/basilisk/src/common.h",682);return _ret;}
end_tracing("mpi_all_reduce0","/home/gongminjiang/basilisk/src/common.h",683);}
#define mpi_all_reduce(v,type,op) {\
  prof_start ("mpi_all_reduce");\
  union { int a; float b; double c;} global;\
  mpi_all_reduce0 (&(v), &global, 1, type, op, MPI_COMM_WORLD);\
  memcpy (&(v), &global, sizeof (v));\
  prof_stop();\
}\

#line 691

#define mpi_all_reduce_array(v,type,op,elem) {\
  prof_start ("mpi_all_reduce");\
  type global[elem], tmp[elem];\
  for (int i = 0; i < elem; i++)\
    tmp[i] = (v)[i];\
  MPI_Datatype datatype;\
  if (!strcmp(#type, "double")) datatype = MPI_DOUBLE;\
  else if (!strcmp(#type, "int")) datatype = MPI_INT;\
  else if (!strcmp(#type, "long")) datatype = MPI_LONG;\
  else if (!strcmp(#type, "bool")) datatype = MPI_C_BOOL;\
  else {\
    fprintf (stderr, "unknown reduction type '%s'\n", #type);\
    fflush (stderr);\
    abort();\
  }\
  mpi_all_reduce0 (tmp, global, elem, datatype, op, MPI_COMM_WORLD);\
  for (int i = 0; i < elem; i++)\
    (v)[i] = global[i];\
  prof_stop();\
}\

#line 712


#endif

#define QFILE FILE

FILE * qstderr (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "log-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systderr;
  }
  return fp;
}

FILE * qstdout (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "out-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systdout;
  }
  return fp;
}

static void finalize (void)
{
  MPI_Finalize();
}

void mpi_init()
{
  int initialized;
  MPI_Initialized (&initialized);
  if (!initialized) {
    MPI_Init (NULL, NULL);
    MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    atexit (finalize);
  }
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_npe);
  srand (mpi_rank + 1);
  if (ferr == NULL) {
    if (mpi_rank > 0) {
      ferr = fopen ("/dev/null", "w");
      fout = fopen ("/dev/null", "w");
    }
    else {
      ferr = systderr;
      fout = systdout;
    }
    char * etrace = getenv ("MALLOC_TRACE"), name[80];
    if (etrace && mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      setenv ("MALLOC_TRACE", name, 1);
    }
#if MTRACE == 1
    etrace = getenv ("MTRACE");
    if (!etrace)
      etrace = "mtrace";
    if (mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      pmtrace.fp = fopen (name, "w");
      pmtrace.fname = systrdup(name);
    }
    else {
      pmtrace.fp = fopen (etrace, "w");
      pmtrace.fname = systrdup(etrace);
    }
#endif
  }
}

#else

#define tid() 0
#define pid() 0
#define npe() 1
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_array(v,type,op,elem)

#endif

#define OMP_PARALLEL() OMP(omp parallel)

#define NOT_UNUSED(x) (void)(x)

#define VARIABLES ;
#define _index(a,m) (a.i)
#define val(a,k,l,m) data(k,l,m)[_index(a,m)]

double _val_higher_dimension = 0.;
#line 823 "/home/gongminjiang/basilisk/src/common.h"
#if (_GNU_SOURCE || __APPLE__) && !_OPENMP && !_CADNA
double undefined;
# if __APPLE__
# include <stdint.h>
# include "fp_osx.h"
# endif
# define enable_fpe(flags) feenableexcept (flags)
# define disable_fpe(flags) fedisableexcept (flags)
static void set_fpe (void) {
  int64_t lnan = 0x7ff0000000000001;
  if (!(sizeof (int64_t) == sizeof (double))) qassert ("/home/gongminjiang/basilisk/src/common.h", 833, "sizeof (int64_t) == sizeof (double)");
  memcpy (&undefined, &lnan, sizeof (double));
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#else
# define undefined ((double) DBL_MAX)
# define enable_fpe(flags)
# define disable_fpe(flags)
static void set_fpe (void) {}
#endif


typedef struct {
  long n;
  long tn;
  int depth;
  int maxdepth;
} Grid;
Grid * grid = NULL;

double X0 = 0., Y0 = 0., Z0 = 0.;

double L0 = 1.;


int N = 64;




typedef struct { int i; } scalar;

typedef struct {
  scalar x;

  scalar y;




} vector;

typedef struct {
  scalar * x;

  scalar * y;




} vectorl;

typedef struct {
  vector x;

  vector y;




} tensor;

struct { int x, y, z; } Period = {false, false, false};

typedef struct {
  double x, y, z;
} coord;

OMP(omp declare reduction (+ : coord :
      omp_out.x += omp_in.x,
      omp_out.y += omp_in.y,
      omp_out.z += omp_in.z))
#line 917 "/home/gongminjiang/basilisk/src/common.h"
void normalize (coord * n)
{
  double norm = 0.;
  
    norm += sq(n->x);
    
#line 921
norm += sq(n->y);
  norm = sqrt(norm);
  
    n->x /= norm;
    
#line 924
n->y /= norm;
}

struct _origin { double x, y, z; };

void origin (struct _origin p) {
  X0 = p.x; Y0 = p.y; Z0 = p.z;
}

void size (double L) {
  L0 = L;
}

double zero (double s0, double s1, double s2) { return 0.; }






  enum { right, left, top, bottom };



int nboundary = 2*2;



#define dirichlet(expr) (2.*(expr) - val(_s,0,0,0))
#define dirichlet_homogeneous() (- val(_s,0,0,0))
#define dirichlet_face(expr) (expr)
#define dirichlet_face_homogeneous() (0.)
#define neumann(expr) (Delta*(expr) + val(_s,0,0,0))
#define neumann_homogeneous() (val(_s,0,0,0))

double * _constant = NULL;
size_t datasize = 0;
typedef struct _Point Point;

#line 1 "/home/gongminjiang/basilisk/src/grid/boundaries.h"


typedef struct _Boundary Boundary;

struct _Boundary {
  void (* destroy) (Boundary * b);
  void (* level) (const Boundary * b, scalar * list, int l);

  void (* restriction) (const Boundary * b, scalar * list, int l);
};

static Boundary ** boundaries = NULL;

void add_boundary (Boundary * b) {
  int len = 0;
  if (boundaries) {
    Boundary ** i = boundaries;
    while (*i++) len++;
  }
  boundaries = (Boundary * *) prealloc (boundaries, (len + 2)*sizeof(Boundary *),__func__,__FILE__,__LINE__);
  boundaries[len] = b;
  boundaries[len+1] = NULL;
}

void free_boundaries() {
  if (!boundaries)
    return;
  Boundary ** i = boundaries, * b;
  while ((b = *i++))
    if (b->destroy)
      b->destroy (b);
    else
      pfree (b,__func__,__FILE__,__LINE__);
  pfree (boundaries,__func__,__FILE__,__LINE__);
  boundaries = NULL;
}
#line 47 "/home/gongminjiang/basilisk/src/grid/boundaries.h"
typedef struct {
  Boundary parent;
  int d;
} BoxBoundary;
#line 964 "/home/gongminjiang/basilisk/src/common.h"



typedef struct {
  double (** boundary) (Point, Point, scalar, void *);
  double (** boundary_homogeneous) (Point, Point, scalar, void *);
  double (* gradient) (double, double, double);
  void (* delete) (scalar);
  char * name;
  struct {
    int x;

    int y;




  } d;
  vector v;
  int face;
  bool nodump, freed;
  int block;
  scalar * depends;

  
#line 19 "/home/gongminjiang/basilisk/src/grid/stencils.h"
bool input, output;
  int width;
  int dirty;
  
#line 18 "/home/gongminjiang/basilisk/src/grid/multigrid-common.h"
void (* prolongation) (Point, scalar);
  void (* restriction) (Point, scalar);
  
#line 9 "/home/gongminjiang/basilisk/src/grid/tree-common.h"
void (* refine) (Point, scalar);
  
#line 97
void (* coarsen) (Point, scalar);
  
#line 28 "/home/gongminjiang/basilisk/src/vof.h"
scalar * tracers, c;
  bool inverse;
  
#line 82 "/home/gongminjiang/basilisk/src/fractions.h"
vector n;
  
#line 21 "/home/gongminjiang/basilisk/src/iforce.h"
scalar phi;
  
#line 456 "/home/gongminjiang/basilisk/src/heights.h"
vector height;
  
#line 22 "/home/gongminjiang/basilisk/src/tension.h"
double sigma;

#line 987 "/home/gongminjiang/basilisk/src/common.h"
} _Attributes;

static _Attributes * _attribute = NULL;






int list_len (scalar * list)
{
  if (!list) return 0;
  int ns = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ ns++;}}
  return ns;
}

scalar * list_append (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  list[len] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_prepend (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  for (int i = len; i >= 1; i--)
    list[i] = list[i-1];
  list[0] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_add (scalar * list, scalar s)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){
    if (t.i == s.i)
      return list;}}
  return list_append (list, s);
}

int list_lookup (scalar * l, scalar s)
{
  if (l != NULL)
    {scalar*_i=(scalar*)( l);if(_i)for(scalar s1=*_i;(&s1)->i>=0;s1=*++_i){
      if (s1.i == s.i)
 return true;}}
  return false;
}

scalar * list_copy (scalar * l)
{
  scalar * list = NULL;
  if (l != NULL)
    {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      list = list_append (list, s);}}
  return list;
}

scalar * list_concat (scalar * l1, scalar * l2)
{
  scalar * l3 = list_copy (l1);
  {scalar*_i=(scalar*)( l2);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    l3 = list_append (l3, s);}}
  return l3;
}

void list_print (scalar * l, FILE * fp)
{
  int i = 0;
  {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    fprintf (fp, "%s%s", i++ == 0 ? "{" : ",", _attribute[s.i].name);}}
  fputs (i > 0 ? "}\n" : "{}\n", fp);
}

int vectors_len (vector * list)
{
  if (!list) return 0;
  int nv = 0;
  {vector*_i=(vector*)( list);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ nv++;}}
  return nv;
}

vector * vectors_append (vector * list, vector v)
{
  int len = vectors_len (list);
  list = (vector *) prealloc (list, (len + 2)*sizeof(vector),__func__,__FILE__,__LINE__);
  list[len] = v;
  list[len + 1] = (vector){{-1}};
  return list;
}

vector * vectors_add (vector * list, vector v)
{
  {vector*_i=(vector*)( list);if(_i)for(vector w=*_i;(&w)->x.i>=0;w=*++_i){ {
    bool id = true;
    
      if (w.x.i != v.x.i)
 id = false;
      
#line 1088
if (w.y.i != v.y.i)
 id = false;
    if (id)
      return list;
  }}}
  return vectors_append (list, v);
}

vector * vectors_copy (vector * l)
{
  vector * list = NULL;
  if (l != NULL)
    {vector*_i=(vector*)( l);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
      list = vectors_append (list, v);}}
  return list;
}

vector * vectors_from_scalars (scalar * s)
{
  vector * list = NULL;
  while (s->i >= 0) {
    vector v;
     {
      if (!(s->i >= 0)) qassert ("/home/gongminjiang/basilisk/src/common.h", 1111, "s->i >= 0");
      v.x = *s++;
    } 
#line 1110
{
      if (!(s->i >= 0)) qassert ("/home/gongminjiang/basilisk/src/common.h", 1111, "s->i >= 0");
      v.y = *s++;
    }
    list = vectors_append (list, v);
  }
  return list;
}

int tensors_len (tensor * list)
{
  if (!list) return 0;
  int nt = 0;
  {tensor*_i=(tensor*)( list);if(_i)for(tensor t=*_i;(&t)->x.x.i>=0;t=*++_i){ nt++;}}
  return nt;
}

tensor * tensors_append (tensor * list, tensor t)
{
  int len = tensors_len (list);
  list = (tensor *) prealloc (list, (len + 2)*sizeof(tensor),__func__,__FILE__,__LINE__);
  list[len] = t;
  list[len + 1] = (tensor){{{-1}}};
  return list;
}

tensor * tensors_from_vectors (vector * v)
{
  tensor * list = NULL;
  while (v->x.i >= 0) {
    tensor t;
     {
      if (!(v->x.i >= 0)) qassert ("/home/gongminjiang/basilisk/src/common.h", 1142, "v->x.i >= 0");
      t.x = *v++;
    } 
#line 1141
{
      if (!(v->y.i >= 0)) qassert ("/home/gongminjiang/basilisk/src/common.h", 1142, "v->x.i >= 0");
      t.y = *v++;
    }
    list = tensors_append (list, t);
  }
  return list;
}

static inline bool is_vertex_scalar (scalar s)
{
  
    if (_attribute[s.i].d.x != -1)
      return false;
    
#line 1153
if (_attribute[s.i].d.y != -1)
      return false;
  return true;
}

scalar * all = NULL;
scalar * baseblock = NULL;



scalar (* init_scalar) (scalar, const char *);
scalar (* init_vertex_scalar) (scalar, const char *);
vector (* init_vector) (vector, const char *);
tensor (* init_tensor) (tensor, const char *);
vector (* init_face_vector) (vector, const char *);





typedef struct _Event Event;
typedef int (* Expr) (int *, double *, Event *);

struct _Event {
  int last, nexpr;
  int (* action) (const int, const double, Event *);
  Expr expr[3];
  int * arrayi;
  double * arrayt;
  char * file;
  int line;
  char * name;
  double t;
  int i, a;
  void * data;
  Event * next;
};

static Event * Events = NULL;

int iter = 0, inext = 0;
double t = 0, tnext = 0;
void init_events (void);
void event_register (Event event);
static void _init_solver (void);

void init_solver()
{
  Events = pmalloc (sizeof (Event),__func__,__FILE__,__LINE__);
  Events[0].last = 1;
  _attribute = pcalloc (datasize/sizeof(double), sizeof (_Attributes),__func__,__FILE__,__LINE__);
  int n = datasize/sizeof(double);
  all = (scalar *) pmalloc (sizeof (scalar)*(n + 1),__func__,__FILE__,__LINE__);
  baseblock = (scalar *) pmalloc (sizeof (scalar)*(n + 1),__func__,__FILE__,__LINE__);
  for (int i = 0; i < n; i++)
    baseblock[i].i = all[i].i = i;
  baseblock[n].i = all[n].i = -1;
#if _CADNA
  cadna_init (-1);
#endif
#if 1
  mpi_init();
#elif MTRACE == 1
  char * etrace = getenv ("MTRACE");
  pmtrace.fp = fopen (etrace ? etrace : "mtrace", "w");
  pmtrace.fname = systrdup (etrace ? etrace : "mtrace");
#endif
}



#if 1
static double mpi_time = 0.;
#endif

typedef struct {
  clock_t c;
  struct timeval tv;
  double tm;
} timer;

timer timer_start (void)
{
  timer t;
  t.c = clock();
  gettimeofday (&t.tv, NULL);
#if 1
  t.tm = mpi_time;
#endif
  return t;
}

double timer_elapsed (timer t)
{
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  return ((tvend.tv_sec - t.tv.tv_sec) +
   (tvend.tv_usec - t.tv.tv_usec)/1e6);
}



const vector zerof = {{_NVARMAX+0},{_NVARMAX+1}};
const vector unityf = {{_NVARMAX+2},{_NVARMAX+3}};
const scalar unity = {_NVARMAX+4};
const scalar zeroc = {_NVARMAX+5};



        vector fm = {{_NVARMAX+2},{_NVARMAX+3}};
        scalar cm = {_NVARMAX+4};
#line 1277 "/home/gongminjiang/basilisk/src/common.h"
static FILE ** qpopen_pipes = NULL;

FILE * qpopen (const char * command, const char * type)
{
  if (pid() > 0)
    return fopen ("/dev/null", type);
  FILE * fp = popen (command, type);
  if (fp) {
    FILE ** i = qpopen_pipes;
    int n = 0;
    while (i && *i) { n++; i++; }
    qpopen_pipes = (FILE * *) prealloc (qpopen_pipes, (n + 2)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    qpopen_pipes[n] = fp;
    qpopen_pipes[n+1] = NULL;
  }
  return fp;
}

int qpclose (FILE * fp)
{
  if (pid() > 0)
    return fclose (fp);
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i == fp)
      *i = (FILE *) 1;
    i++;
  }
  return pclose (fp);
}

static void qpclose_all()
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i != (FILE *) 1)
      pclose (*i);
    i++;
  }
  pfree (qpopen_pipes,__func__,__FILE__,__LINE__);
  qpopen_pipes = NULL;
}






FILE * lfopen (const char * name, const char * mode)
{
  char fname[80];
  sprintf (fname, "%s-%d", name, pid());
  return fopen (fname, mode);
}



void * matrix_new (int n, int p, size_t size)
{
  void ** m = ((void * *) pmalloc ((n)*sizeof(void *),__func__,__FILE__,__LINE__));
  char * a = ((char *) pmalloc ((n*p*size)*sizeof(char),__func__,__FILE__,__LINE__));
  for (int i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return m;
}

double matrix_inverse (double ** m, int n, double pivmin)
{
  int indxc[n], indxr[n], ipiv[n];
  int i, icol = 0, irow = 0, j, k, l, ll;
  double big, dum, pivinv, minpiv = HUGE;

  for (j = 0; j < n; j++)
    ipiv[j] = -1;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
 for (k = 0; k < n; k++) {
   if (ipiv[k] == -1) {
     if (fabs (m[j][k]) >= big) {
       big = fabs (m[j][k]);
       irow = j;
       icol = k;
     }
   }
 }
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++)
 do { double __tmp = m[irow][l]; m[irow][l] = m[icol][l]; m[icol][l] = __tmp; } while(0);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin)
      return 0.;
    if (fabs (m[icol][icol]) < minpiv)
      minpiv = fabs (m[icol][icol]);
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
 dum = m[ll][icol];
 m[ll][icol] = 0.0;
 for (l = 0; l < n; l++)
   m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
 do { double __tmp = m[k][indxr[l]]; m[k][indxr[l]] = m[k][indxc[l]]; m[k][indxc[l]] = __tmp; } while(0);
  }
  return minpiv;
}

void matrix_free (void * m)
{
  pfree (((void **) m)[0],__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
}



typedef void (* free_solver_func) (void);

static Array * free_solver_funcs = NULL;

void free_solver_func_add (free_solver_func func)
{
  if (!free_solver_funcs)
    free_solver_funcs = array_new();
  array_append (free_solver_funcs, &func, sizeof(free_solver_func));
}



static char * display_defaults = NULL;

struct _display {
  const char * commands;
  bool overwrite;
};

static void free_display_defaults() {
  pfree (display_defaults,__func__,__FILE__,__LINE__);
}

void display (struct _display p)
{
  if (display_defaults == NULL)
    free_solver_func_add (free_display_defaults);
  if (p.overwrite) {
    pfree (display_defaults,__func__,__FILE__,__LINE__);
    display_defaults = pmalloc (strlen(p.commands) + 2,__func__,__FILE__,__LINE__);
    strcpy (display_defaults, "@");
    strcat (display_defaults, p.commands);
  }
  else {
    if (!display_defaults)
      display_defaults = pstrdup ("@",__func__,__FILE__,__LINE__);
    display_defaults =
      prealloc (display_defaults,
        strlen(display_defaults) + strlen(p.commands) + 1,__func__,__FILE__,__LINE__);
    strcat (display_defaults, p.commands);
  }
}







#line 1 "/home/gongminjiang/basilisk/src/grid/stencils.h"
#line 17 "/home/gongminjiang/basilisk/src/grid/stencils.h"










typedef struct {
  const char * fname;
  int line;
  int first;
  int face;
  bool vertex;
} ForeachData;


#define foreach_stencil() {\
  static ForeachData _loop = {\
    __FILE__, __LINE__,\
    1, 0, 0\
  };\
  if (baseblock) for (scalar s = baseblock[0], * i = baseblock;\
  s.i >= 0; i++, s = *i) {\
    _attribute[s.i].input = _attribute[s.i].output = false;\
    _attribute[s.i].width = 0;\
  }\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0}; NOT_UNUSED (point);\

#line 48


#define end_foreach_stencil()\
  end_stencil (&_loop);\
  _loop.first = 0;\
}\

#line 54


#define foreach_vertex_stencil() foreach_stencil() _loop.vertex = true;
#define end_foreach_vertex_stencil() end_foreach_stencil()

#define foreach_face_stencil() foreach_stencil()
#define end_foreach_face_stencil() end_foreach_stencil()

#define foreach_visible_stencil(...) foreach_stencil()
#define end_foreach_visible_stencil(...) end_foreach_stencil()

#define _stencil_is_face_x() { _loop.face |= (1 << 0);
#define end__stencil_is_face_x() }
#define _stencil_is_face_y() { _loop.face |= (1 << 1);
#define end__stencil_is_face_y() }
#define _stencil_is_face_z() { _loop.face |= (1 << 2);
#define end__stencil_is_face_z() }

void stencil_val (Point p, scalar s, int i, int j, int k,
    const char * file, int line, bool overflow);
void stencil_val_a (Point p, scalar s, int i, int j, int k, bool input,
      const char * file, int line);

#define _stencil_val(a,_i,_j,_k)\
  stencil_val (point, a, _i, _j, _k, __FILE__, __LINE__, false)\

#line 79

#define _stencil_val_o(a,_i,_j,_k)\
  stencil_val (point, a, _i, _j, _k, __FILE__, __LINE__, true)\

#line 82

#define _stencil_val_a(a,_i,_j,_k)\
  stencil_val_a (point, a, _i, _j, _k, false, __FILE__, __LINE__)\

#line 85

#define _stencil_val_r(a,_i,_j,_k)\
  stencil_val_a (point, a, _i, _j, _k, true, __FILE__, __LINE__)\

#line 88


#define _stencil_fine(a,_i,_j,_k) _stencil_val(a,_i,_j,_k)
#define _stencil_fine_a(a,_i,_j,_k) _stencil_val_a(a,_i,_j,_k)
#define _stencil_fine_r(a,_i,_j,_k) _stencil_val_r(a,_i,_j,_k)

#define _stencil_coarse(a,_i,_j,_k) _stencil_val(a,_i,_j,_k)
#define _stencil_coarse_a(a,_i,_j,_k) _stencil_val_a(a,_i,_j,_k)
#define _stencil_coarse_r(a,_i,_j,_k) _stencil_val_r(a,_i,_j,_k)

#define r_assign(x)
#define _assign(x)

#define _stencil_neighbor(i,j,k)
#define _stencil_child(i,j,k)
#define _stencil_aparent(i,j,k)
#define _stencil_aparent_a(i,j,k)
#define _stencil_aparent_r(i,j,k)

#define _stencil_neighborp(i,j,k) neighborp(i,j,k)

int _stencil_nop;
#define _stencil_val_higher_dimension (_stencil_nop = 1)
#define _stencil__val_constant(a,_i,_j,_k) (_stencil_nop = 1)

typedef void _stencil_undefined;

#define o_stencil -2







static inline bool scalar_is_dirty (scalar s)
{
  if (_attribute[s.i].dirty)
    return true;
  scalar * depends = _attribute[s.i].depends;
  {scalar*_i=(scalar*)( depends);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
    if (_attribute[d.i].dirty)
      return true;}}
  return false;
}




static inline bool scalar_depends_from (scalar a, scalar b)
{
  scalar * depends = _attribute[a.i].depends;
  {scalar*_i=(scalar*)( depends);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (s.i == b.i)
      return true;}}
  return false;
}







void boundary_internal (scalar * list, const char * fname, int line);
void (* boundary_face) (vectorl);







void end_stencil (ForeachData * loop)
{
  scalar * listc = NULL, * dirty = NULL;
  vectorl listf = {NULL};
  bool flux = false;




  {scalar*_i=(scalar*)( baseblock);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    bool write = _attribute[s.i].output, read = _attribute[s.i].input;




    {





      if (read && scalar_is_dirty (s)) {





 if (_attribute[s.i].face) {
   if (_attribute[s.i].width > 0)
     listc = list_append (listc, s);
   else if (!write) {
     scalar sn = _attribute[s.i].v.x.i >= 0 ? _attribute[s.i].v.x : s;
     
       if (_attribute[s.i].v.x.i == s.i) {




  if (_attribute[sn.i].boundary[left] || _attribute[sn.i].boundary[right])
    listc = list_append (listc, s);
  else if (_attribute[s.i].dirty != 2) {
    listf.x = list_append (listf.x, s);
    flux = true;
  }
       }
       
#line 194
if (_attribute[s.i].v.y.i == s.i) {




  if (_attribute[sn.i].boundary[bottom] || _attribute[sn.i].boundary[top])
    listc = list_append (listc, s);
  else if (_attribute[s.i].dirty != 2) {
    listf.y = list_append (listf.y, s);
    flux = true;
  }
       }
   }
 }





 else if (_attribute[s.i].width > 0)
   listc = list_append (listc, s);
      }





      if (write) {
 if (2 > 1 && !loop->vertex && loop->first) {
   bool vertex = true;
   
     if (_attribute[s.i].d.x != -1)
       vertex = false;
     
#line 225
if (_attribute[s.i].d.y != -1)
       vertex = false;
   if (vertex)
     fprintf (ferr,
       "%s:%d: warning: vertex scalar '%s' should be assigned with"
       " a foreach_vertex() loop\n",
       loop->fname, loop->line, _attribute[s.i].name);
 }
 if (_attribute[s.i].face) {
   if (loop->face == 0 && loop->first)
     fprintf (ferr,
       "%s:%d: warning: face vector '%s' should be assigned with"
       " a foreach_face() loop\n",
       loop->fname, loop->line, _attribute[s.i].name);
 }
 else if (loop->face) {
   if (_attribute[s.i].v.x.i < 0) {
     int d = 1, i = 0;
      {
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.x.i = s.i;
  _attribute[s.i].boundary[left] = _attribute[s.i].boundary[right] = NULL;





       }
       d *= 2, i++;
     } 
#line 243
{
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.y.i = s.i;
  _attribute[s.i].boundary[bottom] = _attribute[s.i].boundary[top] = NULL;





       }
       d *= 2, i++;
     }
     if (!_attribute[s.i].face && loop->first)
       fprintf (ferr,
         "%s:%d: warning: scalar '%s' should be assigned with "
         "a foreach_face(x|y|z) loop\n",
         loop->fname, loop->line, _attribute[s.i].name);
   }
   else {
     char * name = NULL;
     if (_attribute[s.i].name) {
       name = pstrdup (_attribute[s.i].name,__func__,__FILE__,__LINE__);
       char * s = name + strlen(name) - 1;
       while (s != name && *s != '.') s--;
       if (s != name) *s = '\0';
     }
     struct { int x, y, z; } input, output;
     vector v = _attribute[s.i].v;

     
       input.x = _attribute[v.x.i].input, output.x = _attribute[v.x.i].output;
       
#line 273
input.y = _attribute[v.y.i].input, output.y = _attribute[v.y.i].output;

     init_face_vector (v, name);


     
       _attribute[v.x.i].input = input.x, _attribute[v.x.i].output = output.x;
       
#line 279
_attribute[v.y.i].input = input.y, _attribute[v.y.i].output = output.y;





     pfree (name,__func__,__FILE__,__LINE__);
   }
 }
 else if (loop->vertex) {
   bool vertex = true;
   
     if (_attribute[s.i].d.x != -1)
       vertex = false;
     
#line 291
if (_attribute[s.i].d.y != -1)
       vertex = false;
   if (!vertex) {
     char * name = NULL;
     if (_attribute[s.i].name) name = pstrdup (_attribute[s.i].name,__func__,__FILE__,__LINE__);
     init_vertex_scalar (s, name);
     
       _attribute[s.i].v.x.i = -1;
       
#line 298
_attribute[s.i].v.y.i = -1;




     pfree (name,__func__,__FILE__,__LINE__);
   }
 }





 dirty = list_append (dirty, s);
 {scalar*_i=(scalar*)( baseblock);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
   if (scalar_depends_from (d, s))
     dirty = list_append (dirty, d);}}
      }
    }
  }}}




  if (flux) {
#line 335 "/home/gongminjiang/basilisk/src/grid/stencils.h"
    boundary_face (listf);
    
      pfree (listf.x,__func__,__FILE__,__LINE__);
      
#line 337
pfree (listf.y,__func__,__FILE__,__LINE__);
  }




  if (listc) {






    boundary_internal (listc, loop->fname, loop->line);
    pfree (listc,__func__,__FILE__,__LINE__);
  }





  if (dirty) {






    {scalar*_i=(scalar*)( dirty);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = true;}}
    pfree (dirty,__func__,__FILE__,__LINE__);
  }
}
#line 1453 "/home/gongminjiang/basilisk/src/common.h"
#line 14 "atomisation-cpp.c"
#line 1 "grid/quadtree.h"
#line 1 "/home/gongminjiang/basilisk/src/grid/quadtree.h"


#line 1 "grid/tree.h"
#line 1 "/home/gongminjiang/basilisk/src/grid/tree.h"
#line 1 "grid/mempool.h"
#line 1 "/home/gongminjiang/basilisk/src/grid/mempool.h"





typedef struct _Pool Pool;

struct _Pool {
  Pool * next;
};

typedef struct {
  char * first, * lastb;
  size_t size;
  size_t poolsize;
  Pool * pool, * last;
} Mempool;

typedef struct {
  char * next;
} FreeBlock;

Mempool * mempool_new (size_t poolsize, size_t size)
{

  if (!(poolsize % 8 == 0)) qassert ("/home/gongminjiang/basilisk/src/grid/mempool.h", 26, "poolsize % 8 == 0");
  if (!(size >= sizeof(FreeBlock))) qassert ("/home/gongminjiang/basilisk/src/grid/mempool.h", 27, "size >= sizeof(FreeBlock)");


  poolsize = min(1 << 20, poolsize + sizeof(Pool));
  Mempool * m = ((Mempool *) pcalloc (1, sizeof(Mempool),__func__,__FILE__,__LINE__));
  m->poolsize = poolsize;
  m->size = size;
  return m;
}

void mempool_destroy (Mempool * m)
{
  Pool * p = m->pool;
  while (p) {
    Pool * next = p->next;
    pfree (p,__func__,__FILE__,__LINE__);
    p = next;
  }
  pfree (m,__func__,__FILE__,__LINE__);
}

void * mempool_alloc (Mempool * m)
{
  if (!m->first) {

    Pool * p = (Pool *) pmalloc (m->poolsize,__func__,__FILE__,__LINE__);
    p->next = NULL;
    if (m->last)
      m->last->next = p;
    else
      m->pool = p;
    m->last = p;
    m->first = m->lastb = ((char *)m->last) + sizeof(Pool);
    FreeBlock * b = (FreeBlock *) m->first;
    b->next = NULL;
  }
  void * ret = m->first;
  FreeBlock * b = (FreeBlock *) ret;
  char * next = b->next;
  if (!next) {
    m->lastb += m->size;
    next = m->lastb;
    if (next + m->size > ((char *) m->last) + m->poolsize)
      next = NULL;
    else {
      FreeBlock * b = (FreeBlock *) next;
      b->next = NULL;
    }
  }
  m->first = next;
#if TRASH
  double * v = (double *) ret;
  for (int i = 0; i < m->size/sizeof(double); i++)
    v[i] = undefined;
#endif
  return ret;
}

void * mempool_alloc0 (Mempool * m)
{
  void * ret = mempool_alloc (m);
  memset (ret, 0, m->size);
  return ret;
}

void mempool_free (Mempool * m, void * p)
{
#if TRASH
  double * v = (double *) p;
  for (int i = 0; i < m->size/sizeof(double); i++)
    v[i] = undefined;
#endif
  FreeBlock * b = (FreeBlock *) p;
  b->next = m->first;
  m->first = (char *) p;
}
#line 2 "/home/gongminjiang/basilisk/src/grid/tree.h"




#line 1 "grid/memindex/range.h"
#line 1 "/home/gongminjiang/basilisk/src/grid/memindex/range.h"
#line 15 "/home/gongminjiang/basilisk/src/grid/memindex/range.h"
typedef struct {
  void ** p;
  int size;
} Memalloc;

typedef struct {
  int start, end;
} Memrange;
#line 34 "/home/gongminjiang/basilisk/src/grid/memindex/range.h"
void memrange_alloc (Memrange * r, Memalloc * mem, int i)
{
  if (r->start == r->end) {
    r->start = i;
    r->end = i + 1;
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = pcalloc (1, m->size,__func__,__FILE__,__LINE__);
      *m->p = (char *)(*m->p) - i*m->size;
    }
  }
  else if (i >= r->end) {
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = prealloc ((char *)(*m->p) + r->start*m->size,
         m->size*(i + 1 - r->start),__func__,__FILE__,__LINE__);
      *m->p = (char *)(*m->p) - r->start*m->size;
      memset ((char *)(*m->p) + r->end*m->size, 0, (i - r->end + 1)*m->size);
    }
    r->end = i + 1;
  }
  else if (i < r->start) {
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = prealloc ((char *)(*m->p) + r->start*m->size, m->size*(r->end - i),__func__,__FILE__,__LINE__);
      memmove ((char *)(*m->p) + (r->start - i)*m->size, *m->p,
        m->size*(r->end - r->start));
      memset ((char *)(*m->p), 0, (r->start - i)*m->size);
      *m->p = (char *)(*m->p) - i*m->size;
    }
    r->start = i;
  }
}
#line 73 "/home/gongminjiang/basilisk/src/grid/memindex/range.h"
bool memrange_free (Memrange * r, Memalloc * mem, int i)
{
  if (i == r->start) {
    if (i == r->end - 1) {
      for (Memalloc * m = mem; m->p; m++) {
 pfree ((char *)(*m->p) + r->start*m->size,__func__,__FILE__,__LINE__);
 *m->p = NULL;
      }
      r->start = r->end = 0;
      return true;
    }
    else {
      for (i = i + 1; i < r->end &&
      !*(void **)((char *)(*mem->p) + i*mem->size); i++);
      for (Memalloc * m = mem; m->p; m++) {
 memmove ((char *)(*m->p) + r->start*m->size,
   (char *)(*m->p) + i*m->size, m->size*(r->end - i));
 *m->p = prealloc ((char *)(*m->p) + r->start*m->size,
    m->size*(r->end - i),__func__,__FILE__,__LINE__);
 *m->p = (char *)(*m->p) - i*m->size;
      }
      r->start = i;
    }
  }
  else if (i == r->end - 1) {
    for (i = i - 1; i >= r->start &&
    !*(void **)((char *)(*mem->p) + i*mem->size); i--);
    r->end = i + 1;
    for (Memalloc * m = mem; m->p; m++) {
      *m->p = prealloc ((char *)(*m->p) + r->start*m->size,
         m->size*(r->end - r->start),__func__,__FILE__,__LINE__);
      *m->p = (char *)(*m->p) - r->start*m->size;
    }
  }
  else {
    if (!(i > r->start && i < r->end)) qassert ("/home/gongminjiang/basilisk/src/grid/memindex/range.h", 108, "i > r->start && i < r->end");
    for (Memalloc * m = mem; m->p; m++)
      memset ((char *)(*m->p) + i*m->size, 0, m->size);
  }
  return false;
}







struct _Memindex {
  Memrange r1;

  Memrange * r2;







  char *** b;



};
#line 171 "/home/gongminjiang/basilisk/src/grid/memindex/range.h"
struct _Memindex * mem_new (int len)
{
  struct _Memindex * m = pcalloc (1, sizeof (struct _Memindex),__func__,__FILE__,__LINE__);
  return m;
}





void mem_destroy (struct _Memindex * m, int len)
{

  for (int i = m->r1.start; i < m->r1.end; i++)
    if (m->b[i]) {






      pfree (m->b[i] + m->r2[i].start,__func__,__FILE__,__LINE__);
    }
  if (m->b) {
    pfree (m->r2 + m->r1.start,__func__,__FILE__,__LINE__);



  }

  if (m->b)
    pfree (m->b + m->r1.start,__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
}
#line 218 "/home/gongminjiang/basilisk/src/grid/memindex/range.h"
void mem_assign (struct _Memindex * m, int i, int j, int len, void * b)
{
  Memalloc mem[] = {{(void **)&m->b, sizeof(char **)},
      {(void **)&m->r2, sizeof(Memrange)},
      {NULL}};
  memrange_alloc (&m->r1, mem, i);
  Memalloc mem1[] = {{(void **)&m->b[i], sizeof(char *)},
       {NULL}};
  memrange_alloc (&m->r2[i], mem1, j);
  ((m)->b[i][j]) = b;
}
#line 259 "/home/gongminjiang/basilisk/src/grid/memindex/range.h"
void mem_free (struct _Memindex * m, int i, int j, int len)
{
  Memalloc mem[] = {{(void **)&m->b[i], sizeof(char *)},
      {NULL}};
  if (memrange_free (&m->r2[i], mem, j)) {
    Memalloc mem[] = {{(void **)&m->b, sizeof(char **)},
        {(void **)&m->r2, sizeof(Memrange)},
        {NULL}};
    memrange_free (&m->r1, mem, i);
  }
}
#line 305 "/home/gongminjiang/basilisk/src/grid/memindex/range.h"
#define foreach_mem(_m, _len, _i) {\
  Point point = {0};\
  for (point.i = max(Period.x*2, (_m)->r1.start);\
       point.i < min(_len - Period.x*2, (_m)->r1.end);\
       point.i += _i)\
    if ((_m)->b[point.i])\
      for (point.j = max(Period.y*2, (_m)->r2[point.i].start);\
    point.j < min(_len - Period.y*2, (_m)->r2[point.i].end);\
    point.j += _i)\
 if ((_m)->b[point.i][point.j]) {\

#line 315

#define end_foreach_mem() }}
#line 7 "/home/gongminjiang/basilisk/src/grid/tree.h"
#line 24 "/home/gongminjiang/basilisk/src/grid/tree.h"
typedef struct {
  unsigned short flags;

  unsigned short neighbors;
  int pid;
} Cell;

enum {
  active = 1 << 0,
  leaf = 1 << 1,
  border = 1 << 2,
  vertex = 1 << 3,
  user = 4,

  face_x = 1 << 0

  , face_y = 1 << 1




};

#define is_active(cell) ((cell).flags & active)
#define is_leaf(cell) ((cell).flags & leaf)
#define is_coarse() ((cell).neighbors > 0)
#define is_border(cell) ((cell).flags & border)
#define is_local(cell) ((cell).pid == pid())
#define is_vertex(cell) ((cell).flags & vertex)



typedef struct {
  int i;

  int j;




} IndexLevel;

typedef struct {
  IndexLevel * p;
  int n, nm;
} CacheLevel;

typedef struct {
  int i;

  int j;




  int level, flags;
} Index;

typedef struct {
  Index * p;
  int n, nm;
} Cache;



typedef struct {
  struct _Memindex * m;
  Mempool * pool;
  long nc;
  int len;
} Layer;

static size_t _size (size_t depth)
{
  return (1 << depth) + 2*2;
}

static size_t poolsize (size_t depth, size_t size)
{




  return sq(_size(depth))*size;



}

static Layer * new_layer (int depth)
{
  Layer * l = ((Layer *) pmalloc ((1)*sizeof(Layer),__func__,__FILE__,__LINE__));
  l->len = _size (depth);
  if (depth == 0)
    l->pool = NULL;
  else {
    size_t size = sizeof(Cell) + datasize;


    l->pool = mempool_new (poolsize (depth, size), (1 << 2)*size);
  }
  l->m = mem_new (l->len);
  l->nc = 0;
  return l;
}

static void destroy_layer (Layer * l)
{
  if (l->pool)
    mempool_destroy (l->pool);
  mem_destroy (l->m, l->len);
  pfree (l,__func__,__FILE__,__LINE__);
}



typedef struct {
  Grid g;
  Layer ** L;

  Cache leaves;
  Cache faces;
  Cache vertices;
  Cache refined;
  CacheLevel * active;
  CacheLevel * prolongation;
  CacheLevel * boundary;

  CacheLevel * restriction;

  bool dirty;
} Tree;



struct _Point {

  int i;

  int j;




  int level;
#ifdef foreach_block
  int l;
  #define _BLOCK_INDEX , point.l
#else
  #define _BLOCK_INDEX
#endif
};
static Point last_point;



static void cache_level_append (CacheLevel * c, Point p)
{
  if (c->n >= c->nm) {
    c->nm += 128;
    c->p = (IndexLevel *) prealloc (c->p, (c->nm)*sizeof(IndexLevel),__func__,__FILE__,__LINE__);
  }
  c->p[c->n].i = p.i;

  c->p[c->n].j = p.j;




  c->n++;
}

static void cache_level_shrink (CacheLevel * c)
{
  if (c->nm > (c->n/128 + 1)*128) {
    c->nm = (c->n/128 + 1)*128;
    if (!(c->nm > c->n)) qassert ("/home/gongminjiang/basilisk/src/grid/tree.h", 200, "c->nm > c->n");
    c->p = (IndexLevel *) prealloc (c->p, sizeof (Index)*c->nm,__func__,__FILE__,__LINE__);
  }
}

static void cache_append (Cache * c, Point p, unsigned short flags)
{
  if (c->n >= c->nm) {
    c->nm += 128;
    c->p = (Index *) prealloc (c->p, (c->nm)*sizeof(Index),__func__,__FILE__,__LINE__);
  }
  c->p[c->n].i = p.i;

  c->p[c->n].j = p.j;




  c->p[c->n].level = p.level;
  c->p[c->n].flags = flags;
  c->n++;
}

void cache_shrink (Cache * c)
{
  cache_level_shrink ((CacheLevel *)c);
}
#line 243 "/home/gongminjiang/basilisk/src/grid/tree.h"
#define allocated(k,l,n) (((point.i+k) >= (((Tree *)grid)->L[point.level]->m)->r1.start && (point.i+k) < (((Tree *)grid)->L[point.level]->m->r1.end) && (((Tree *)grid)->L[point.level]->m)->b[point.i+k] && (point.j+l) >= (((Tree *)grid)->L[point.level]->m)->r2[point.i+k].start && (point.j+l) < (((Tree *)grid)->L[point.level]->m)->r2[point.i+k].end && (((Tree *)grid)->L[point.level]->m)->b[point.i+k][point.j+l])\
                               )\

#line 245

#define NEIGHBOR(k,l,n) (((((Tree *)grid)->L[point.level]->m)->b[point.i+k][point.j+l])\
                            )\

#line 248

#define PARENT(k,l,n) (((((Tree *)grid)->L[point.level-1]->m)->b[(point.i+2)/2+k][(point.j+2)/2+l])\
                                                    )\

#line 251

#define allocated_child(k,l,n) (level < depth() &&\
         ((2*point.i-2 +k) >= (((Tree *)grid)->L[point.level+1]->m)->r1.start && (2*point.i-2 +k) < (((Tree *)grid)->L[point.level+1]->m->r1.end) && (((Tree *)grid)->L[point.level+1]->m)->b[2*point.i-2 +k] && (2*point.j-2 +l) >= (((Tree *)grid)->L[point.level+1]->m)->r2[2*point.i-2 +k].start && (2*point.j-2 +l) < (((Tree *)grid)->L[point.level+1]->m)->r2[2*point.i-2 +k].end && (((Tree *)grid)->L[point.level+1]->m)->b[2*point.i-2 +k][2*point.j-2 +l])\
\
                             )\

#line 256

#define CHILD(k,l,n) (((((Tree *)grid)->L[point.level+1]->m)->b[2*point.i-2 +k][2*point.j-2 +l])\
                                                )\

#line 259

#line 284 "/home/gongminjiang/basilisk/src/grid/tree.h"
#define CELL(m) (*((Cell *)(m)))


#define depth() (grid->depth)
#define aparent(k,l,n) CELL(PARENT(k,l,n))
#define child(k,l,n) CELL(CHILD(k,l,n))


#define cell CELL(NEIGHBOR(0,0,0))
#define neighbor(k,l,n) CELL(NEIGHBOR(k,l,n))
#define neighborp(l,m,n) (Point) {\
    point.i + l,\
\
    point.j + m,\
\
\
\
\
    point.level\
    _BLOCK_INDEX\
}\

#line 305



#define data(k,l,n) ((double *) (NEIGHBOR(k,l,n) + sizeof(Cell)))
#define fine(a,k,p,n) ((double *) (CHILD(k,p,n) + sizeof(Cell)))[_index(a,n)]
#define coarse(a,k,p,n) ((double *) (PARENT(k,p,n) + sizeof(Cell)))[_index(a,n)]

#define POINT_VARIABLES\
  VARIABLES\
  int level = point.level; NOT_UNUSED(level);\
\
\
\
  struct { int x, y; } child = {\
    2*((point.i+2)%2)-1, 2*((point.j+2)%2)-1\
  };\
\
\
\
\
\
  NOT_UNUSED(child);\
  Point parent = point; NOT_UNUSED(parent);\
  parent.level--;\
  parent.i = (point.i + 2)/2;\
\
  parent.j = (point.j + 2)/2;\
\

#line 341


#line 1 "grid/foreach_cell.h"
#line 1 "/home/gongminjiang/basilisk/src/grid/foreach_cell.h"
#line 66 "/home/gongminjiang/basilisk/src/grid/foreach_cell.h"
#define foreach_cell_root(root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point = {0};\
\
\
\
    struct { int l, i, j, stage; } stack[20];\
\
\
\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
 POINT_VARIABLES;\
\

#line 89

#define end_foreach_cell_root()\
        if (point.level < grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 1; };\
          { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
        }\
        break;\
      }\
\
\
\
      case 1: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 2; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; }; break;\
      case 2: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 3; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; }; break;\
      case 3: { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; }; break;\
\
      }\
    }\
  }\

#line 123


#define foreach_cell() {\
\
\
\
  Point root = {2,2,0};\
\
\
\
  foreach_cell_root (root)\

#line 134

#define end_foreach_cell() end_foreach_cell_root() }

#define foreach_cell_all() {\
  Point root = {0};\
  for (root.i = 2*Period.x; root.i <= 2*(2 - Period.x); root.i++)\
\
    for (root.j = 2*Period.y; root.j <= 2*(2 - Period.y); root.j++)\
\
\
\
\
 foreach_cell_root (root)\

#line 147

#define end_foreach_cell_all() end_foreach_cell_root() }

#define foreach_cell_post_root(condition, root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point = {0};\
\
\
\
    struct { int l, i, j, stage; } stack[20];\
\
\
\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
        POINT_VARIABLES;\
 if (point.level == grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 8; };\
 }\
 else {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 1; };\
   if (condition)\
     { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
 }\
 break;\
      }\
\
\
\
\
\
\
\
      case 1:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 2; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 2:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 3; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
 break;\
      case 3:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 4; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; };\
 break;\
\
      default: {\
        POINT_VARIABLES;\
\

#line 244

#define end_foreach_cell_post_root()\
      }\
      }\
    }\
  }\

#line 250


#define foreach_cell_post(condition)\
  {\
\
\
\
    Point root = {2,2,0};\
\
\
\
    foreach_cell_post_root(condition, root)\

#line 262

#define end_foreach_cell_post() end_foreach_cell_post_root() }

#define foreach_cell_post_all(condition) {\
  Point root = {0};\
  for (root.i = 0; root.i <= 2*2; root.i++)\
\
    for (root.j = 0; root.j <= 2*2; root.j++)\
\
\
\
\
 foreach_cell_post_root (condition, root)\

#line 275

#define end_foreach_cell_post_all() end_foreach_cell_post_root() }

#define foreach_leaf() foreach_cell()\
  if (is_leaf (cell)) {\
    if (is_active(cell) && is_local(cell)) {\

#line 281

#define end_foreach_leaf() } continue; } end_foreach_cell()
#line 344 "/home/gongminjiang/basilisk/src/grid/tree.h"
#line 361 "/home/gongminjiang/basilisk/src/grid/tree.h"
#define foreach_child() {\
  int _i = 2*point.i - 2, _j = 2*point.j - 2;\
  point.level++;\
  for (int _k = 0; _k < 2; _k++) {\
    point.i = _i + _k;\
    for (int _l = 0; _l < 2; _l++) {\
      point.j = _j + _l;\
      POINT_VARIABLES;\

#line 369

#define end_foreach_child()\
    }\
  }\
  point.i = (_i + 2)/2; point.j = (_j + 2)/2;\
  point.level--;\
}\

#line 376

#define foreach_child_break() _k = _l = 2
#line 407 "/home/gongminjiang/basilisk/src/grid/tree.h"
#define is_refined_check() ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) &&\
    point.i > 0 && point.i < (1 << level) + 2*2 - 1\
\
    && point.j > 0 && point.j < (1 << level) + 2*2 - 1\
\
\
\
\
    )\

#line 416


#define foreach_cache(_cache) {\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.i = 2;\
\
  point.j = 2;\
\
\
\
\
  int _k; unsigned short _flags; NOT_UNUSED(_flags);\
  OMP(omp for schedule(static))\
  for (_k = 0; _k < _cache.n; _k++) {\
    point.i = _cache.p[_k].i;\
\
    point.j = _cache.p[_k].j;\
\
\
\
\
    point.level = _cache.p[_k].level;\
    _flags = _cache.p[_k].flags;\
    POINT_VARIABLES;\

#line 442

#define end_foreach_cache() } } }

#define foreach_cache_level(_cache,_l) {\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.i = 2;\
\
  point.j = 2;\
\
\
\
\
  point.level = _l;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 0; _k < _cache.n; _k++) {\
    point.i = _cache.p[_k].i;\
\
    point.j = _cache.p[_k].j;\
\
\
\
\
    POINT_VARIABLES;\

#line 468

#define end_foreach_cache_level() } } }

#define foreach_boundary_level(_l) {\
  if (_l <= depth()) {\
    { if (((Tree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _boundary = ((Tree *)grid)->boundary[_l];\
    foreach_cache_level (_boundary,_l)\

#line 476

#define end_foreach_boundary_level() end_foreach_cache_level(); }}



#define foreach_boundary(_b) {\
  for (int _l = depth(); _l >= 0; _l--)\
    foreach_boundary_level(_l) {\
      if ((- cell.pid - 1) == _b)\
 for (int _d = 0; _d < 2; _d++) {\
   for (int _i = -1; _i <= 1; _i += 2) {\
     if (_d == 0) ig = _i; else if (_d == 1) jg = _i; else kg = _i;\
     if (allocated(-ig,-jg,-kg) &&\
  is_leaf (neighbor(-ig,-jg,-kg)) &&\
  !(neighbor(-ig,-jg,-kg).pid < 0) &&\
  is_local(neighbor(-ig,-jg,-kg))) {\
       point.i -= ig; x -= ig*Delta/2.;\
\
       point.j -= jg; y -= jg*Delta/2.;\
\
\
\
\

#line 499

#define end_foreach_boundary()\
       point.i += ig; x += ig*Delta/2.;\
\
       point.j += jg; y += jg*Delta/2.;\
\
\
\
\
            }\
   }\
   ig = jg = kg = 0;\
 }\
    } end_foreach_boundary_level(); }\

#line 513


#define foreach_halo(_name,_l) {\
  if (_l <= depth()) {\
    { if (((Tree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _cache = ((Tree *)grid)->_name[_l];\
    foreach_cache_level (_cache, _l)\

#line 520

#define end_foreach_halo() end_foreach_cache_level(); }}

#line 1 "grid/neighbors.h"
#line 1 "/home/gongminjiang/basilisk/src/grid/neighbors.h"
#line 17 "/home/gongminjiang/basilisk/src/grid/neighbors.h"
#define foreach_neighbor(_s) {\
  int _nn = _s + 0 ? _s + 0 : 2;\
  int _i = point.i, _j = point.j;\
  for (int _k = - _nn; _k <= _nn; _k++) {\
    point.i = _i + _k;\
    for (int _l = - _nn; _l <= _nn; _l++) {\
      point.j = _j + _l;\
      POINT_VARIABLES;\

#line 25

#define end_foreach_neighbor()\
    }\
  }\
  point.i = _i; point.j = _j;\
}\

#line 31

#define foreach_neighbor_break() _k = _l = _nn + 1
#line 524 "/home/gongminjiang/basilisk/src/grid/tree.h"

static inline bool has_local_children (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    if (is_local(cell))
      return true;end_foreach_child()}
  return false;
}

static inline void cache_append_face (Point point, unsigned short flags)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  Tree * q = ((Tree *)grid);
  cache_append (&q->faces, point, flags);

  if (!is_vertex(cell)) {
    cache_append (&q->vertices, point, 0);
    cell.flags |= vertex;
  }
  
    if ((flags & face_y) && !is_vertex(neighbor(1,0,0))) {
      cache_append (&q->vertices, neighborp(1,0,0), 0);
      neighbor(1,0,0).flags |= vertex;
    }
    
#line 543
if ((flags & face_x) && !is_vertex(neighbor(0,1,0))) {
      cache_append (&q->vertices, neighborp(0,1,0), 0);
      neighbor(0,1,0).flags |= vertex;
    }
#line 557 "/home/gongminjiang/basilisk/src/grid/tree.h"
}



static void update_cache_f (void)
{
  Tree * q = ((Tree *)grid);

  {foreach_cache (q->vertices)
    if (level <= depth() && allocated(0,0,0))
      cell.flags &= ~vertex;end_foreach_cache();}


  q->leaves.n = q->faces.n = q->vertices.n = 0;
  for (int l = 0; l <= depth(); l++)
    q->active[l].n = q->prolongation[l].n =
      q->boundary[l].n = q->restriction[l].n = 0;

  const unsigned short fboundary = 1 << user;
  {foreach_cell() {



    if (is_local(cell) && is_active(cell)) {


      cache_level_append (&q->active[level], point);
    }
#line 601 "/home/gongminjiang/basilisk/src/grid/tree.h"
    if (!(cell.pid < 0)) {

      {foreach_neighbor (2)
 if (allocated(0,0,0) && (cell.pid < 0) && !(cell.flags & fboundary)) {
   cache_level_append (&q->boundary[level], point);
   cell.flags |= fboundary;
 }end_foreach_neighbor()}
    }

    else if (level > 0 && is_local(aparent(0,0,0)))
      cache_level_append (&q->restriction[level], point);

    if (is_leaf (cell)) {
      if (is_local(cell)) {
 cache_append (&q->leaves, point, 0);

 unsigned short flags = 0;
 
   if ((neighbor(-1,0,0).pid < 0) || (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) ||
       is_leaf(neighbor(-1,0,0)))
     flags |= face_x;
   
#line 619
if ((neighbor(0,-1,0).pid < 0) || (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) ||
       is_leaf(neighbor(0,-1,0)))
     flags |= face_y;
 if (flags)
   cache_append (&q->faces, point, flags);
 
   if ((neighbor(1,0,0).pid < 0) || (!is_leaf(neighbor(1,0,0)) && !neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0) ||
       (!is_local(neighbor(1,0,0)) && is_leaf(neighbor(1,0,0))))
     cache_append (&q->faces, neighborp(1,0,0), face_x);
   
#line 625
if ((neighbor(0,1,0).pid < 0) || (!is_leaf(neighbor(0,1,0)) && !neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0) ||
       (!is_local(neighbor(0,1,0)) && is_leaf(neighbor(0,1,0))))
     cache_append (&q->faces, neighborp(0,1,0), face_y);

 for (int i = 0; i <= 1; i++)

   for (int j = 0; j <= 1; j++)




       if (!is_vertex(neighbor(i,j,k))) {
  cache_append (&q->vertices, neighborp(i,j,k), 0);
  neighbor(i,j,k).flags |= vertex;
       }

        if (cell.neighbors > 0)
   cache_level_append (&q->prolongation[level], point);
      }
      else if (!(cell.pid < 0) || is_local(aparent(0,0,0))) {

 unsigned short flags = 0;
 
   if (allocated(-1,0,0) &&
       is_local(neighbor(-1,0,0)) && (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))
     flags |= face_x;
   
#line 648
if (allocated(0,-1,0) &&
       is_local(neighbor(0,-1,0)) && (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))
     flags |= face_y;
 if (flags)
   cache_append_face (point, flags);
 
   if (allocated(1,0,0) && is_local(neighbor(1,0,0)) &&
       (!is_leaf(neighbor(1,0,0)) && !neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0))
     cache_append_face (neighborp(1,0,0), face_x);
   
#line 654
if (allocated(0,1,0) && is_local(neighbor(0,1,0)) &&
       (!is_leaf(neighbor(0,1,0)) && !neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0))
     cache_append_face (neighborp(0,1,0), face_y);
      }

      continue;

    }
  }end_foreach_cell();}


  cache_shrink (&q->leaves);
  cache_shrink (&q->faces);
  cache_shrink (&q->vertices);
  for (int l = 0; l <= depth(); l++) {
    cache_level_shrink (&q->active[l]);
    cache_level_shrink (&q->prolongation[l]);
    cache_level_shrink (&q->boundary[l]);
    cache_level_shrink (&q->restriction[l]);
}

  q->dirty = false;


  for (int l = depth(); l >= 0; l--)
    {foreach_boundary_level (l)
      cell.flags &= ~fboundary;end_foreach_boundary_level();}



  grid->n = q->leaves.n;

#if !1
  grid->tn = grid->n;
  grid->maxdepth = grid->depth;
#endif
}

#define foreach() { if (((Tree *)grid)->dirty) update_cache_f(); }; foreach_cache(((Tree *)grid)->leaves)
#define end_foreach() end_foreach_cache()

#define foreach_face_generic()\
  { if (((Tree *)grid)->dirty) update_cache_f(); };\
  foreach_cache(((Tree *)grid)->faces) 
#line 716

#define end_foreach_face_generic() end_foreach_cache()

#define is_face_x() { int ig = -1; VARIABLES; if (_flags & face_x) {
#define end_is_face_x() }}


#define is_face_y() { int jg = -1; VARIABLES; if (_flags & face_y) {
#define end_is_face_y() }}






#define foreach_vertex()\
  { if (((Tree *)grid)->dirty) update_cache_f(); };\
  foreach_cache(((Tree *)grid)->vertices) {\
    x -= Delta/2.;\
\
    y -= Delta/2.;\
\
\
\
\

#line 742

#define end_foreach_vertex() } end_foreach_cache()
#line 734 "/home/gongminjiang/basilisk/src/grid/tree.h"
#define foreach_level(l) {\
  if (l <= depth()) {\
    { if (((Tree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _active = ((Tree *)grid)->active[l];\
    foreach_cache_level (_active,l)\

#line 739

#define end_foreach_level() end_foreach_cache_level(); }}

#define foreach_coarse_level(l) foreach_level(l) if (!is_leaf(cell)) {
#define end_foreach_coarse_level() } end_foreach_level()

#define foreach_level_or_leaf(l) {\
  for (int _l1 = l; _l1 >= 0; _l1--)\
    foreach_level(_l1)\
      if (_l1 == l || is_leaf (cell)) {\

#line 749

#define end_foreach_level_or_leaf() } end_foreach_level(); }

#if TRASH
# undef trash
# define trash(list) reset(list, undefined)
#endif

void reset (void * alist, double val)
{
  scalar * list = (scalar *) alist;
  Tree * q = ((Tree *)grid);

  for (int l = 0; l <= depth(); l++) {
    Layer * L = q->L[l];
    {foreach_mem (L->m, L->len, 1) {
      point.level = l;
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
 if (!is_constant(s))
   for (int b = 0; b < _attribute[s.i].block; b++)
     data(0,0,0)[s.i + b] = val;
      }}}
    }end_foreach_mem();}
  }
}

static CacheLevel * cache_level_resize (CacheLevel * name, int a)
{
  for (int i = 0; i <= depth() - a; i++)
    pfree (name[i].p,__func__,__FILE__,__LINE__);
  pfree (name,__func__,__FILE__,__LINE__);
  return ((CacheLevel *) pcalloc (depth() + 1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
}

static void update_depth (int inc)
{
  Tree * q = ((Tree *)grid);
  grid->depth += inc;
  q->L = &(q->L[-1]);
  q->L = (Layer * *) prealloc (q->L, (grid->depth + 2)*sizeof(Layer *),__func__,__FILE__,__LINE__);
  q->L = &(q->L[1]);
  if (inc > 0)
    q->L[grid->depth] = new_layer (grid->depth);
  q->active = cache_level_resize (q->active, inc);
  q->prolongation = cache_level_resize (q->prolongation, inc);
  q->boundary = cache_level_resize (q->boundary, inc);
  q->restriction = cache_level_resize (q->restriction, inc);
}
#line 823 "/home/gongminjiang/basilisk/src/grid/tree.h"
typedef void (* PeriodicFunction) (struct _Memindex *, int, int, int, void *);

static void periodic_function (struct _Memindex * m, int i, int j, int len, void * b,
          PeriodicFunction f)
{
  f(m, i, j, len, b);
  if (Period.x) {
    int nl = len - 2*2;
    for (int l = - 1; l <= 1; l += 2)
      for (int n = i + l*nl; n >= 0 && n < len; n += l*nl)
 f(m, n, j, len, b);
    if (Period.y)
      for (int l = - 1; l <= 1; l += 2)
 for (int n = j + l*nl; n >= 0 && n < len; n += l*nl) {
   f(m, i, n, len, b);
   for (int o = - 1; o <= 1; o += 2)
     for (int p = i + o*nl; p >= 0 && p < len; p += o*nl)
       f(m, p, n, len, b);
 }
  }
  else if (Period.y) {
    int nl = len - 2*2;
    for (int l = - 1; l <= 1; l += 2)
      for (int n = j + l*nl; n >= 0 && n < len; n += l*nl)
 f(m, i, n, len, b);
  }
}

static void assign_periodic (struct _Memindex * m, int i, int j, int len, void * b)
{
  periodic_function (m, i, j, len, b, mem_assign);
}

static void free_periodic (struct _Memindex * m, int i, int j, int len)
{
  periodic_function (m, i, j, len, NULL, (PeriodicFunction) mem_free);
}
#line 938 "/home/gongminjiang/basilisk/src/grid/tree.h"
static void alloc_children (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (point.level == grid->depth)
    update_depth (+1);
  else if (allocated_child(0,0,0))
    return;


  Layer * L = ((Tree *)grid)->L[point.level + 1];
  L->nc++;
  size_t len = sizeof(Cell) + datasize;
  char * b = (char *) mempool_alloc0 (L->pool);
  int i = 2*point.i - 2;
  for (int k = 0; k < 2; k++, i++) {




    int j = 2*point.j - 2;
    for (int l = 0; l < 2; l++, j++) {
      assign_periodic (L->m, i, j, L->len, b);
      b += len;
    }
#line 971 "/home/gongminjiang/basilisk/src/grid/tree.h"
  }

  int pid = cell.pid;
  {foreach_child() {
    cell.pid = pid;
#if TRASH
    {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      val(s,0,0,0) = undefined;}}
#endif
  }end_foreach_child()}
}
#line 1000 "/home/gongminjiang/basilisk/src/grid/tree.h"
static void free_children (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  Layer * L = ((Tree *)grid)->L[point.level + 1];
  int i = 2*point.i - 2, j = 2*point.j - 2;
  if (!(((L->m)->b[i][j]))) qassert ("/home/gongminjiang/basilisk/src/grid/tree.h", 1005, "mem_data (L->m,i,j)");
  mempool_free (L->pool, ((L->m)->b[i][j]));
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++)
      free_periodic (L->m, i + k, j + l, L->len);
  if (--L->nc == 0) {
    destroy_layer (L);
    if (!(point.level + 1 == grid->depth)) qassert ("/home/gongminjiang/basilisk/src/grid/tree.h", 1012, "point.level + 1 == grid->depth");
    update_depth (-1);
  }
}
#line 1041 "/home/gongminjiang/basilisk/src/grid/tree.h"
void increment_neighbors (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  ((Tree *)grid)->dirty = true;
  if (cell.neighbors++ == 0)
    alloc_children (point);
  {foreach_neighbor (2/2)
    if (cell.neighbors++ == 0)
      alloc_children (point);end_foreach_neighbor()}
  cell.neighbors--;
}

void decrement_neighbors (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  ((Tree *)grid)->dirty = true;
  {foreach_neighbor (2/2)
    if (allocated(0,0,0)) {
      cell.neighbors--;
      if (cell.neighbors == 0)
 free_children (point);
    }end_foreach_neighbor()}
  if (cell.neighbors) {
    int pid = cell.pid;
    {foreach_child() {
      cell.flags = 0;
      cell.pid = pid;
    }end_foreach_child()}
  }
}

void realloc_scalar (int size)
{

  Tree * q = ((Tree *)grid);
  size_t oldlen = sizeof(Cell) + datasize;
  size_t newlen = oldlen + size;
  datasize += size;

  Layer * L = q->L[0];
  {foreach_mem (L->m, L->len, 1) {




    char * p = (char *) prealloc (((L->m)->b[point.i][point.j]),
     newlen*sizeof(char),__func__,__FILE__,__LINE__);
    assign_periodic (L->m, point.i, point.j, L->len, p);





  }end_foreach_mem();}

  for (int l = 1; l <= depth(); l++) {
    Layer * L = q->L[l];
    Mempool * oldpool = L->pool;
    L->pool = mempool_new (poolsize (l, newlen), (1 << 2)*newlen);
    {foreach_mem (L->m, L->len, 2) {
      char * new = (char *) mempool_alloc (L->pool);







      for (int k = 0; k < 2; k++)
 for (int o = 0; o < 2; o++) {
   memcpy (new, ((L->m)->b[point.i + k][point.j + o]), oldlen);
   assign_periodic (L->m, point.i + k, point.j + o, L->len, new);
   new += newlen;
 }
#line 1124 "/home/gongminjiang/basilisk/src/grid/tree.h"
    }end_foreach_mem();}
    mempool_destroy (oldpool);
  }
}



#define VN v.x
#define VT v.y
#define VR v.z




#if 1
# define disable_fpe_for_mpi() disable_fpe (FE_DIVBYZERO|FE_INVALID)
# define enable_fpe_for_mpi() enable_fpe (FE_DIVBYZERO|FE_INVALID)
#else
# define disable_fpe_for_mpi()
# define enable_fpe_for_mpi()
#endif

static inline void no_restriction (Point point, scalar s);

static bool normal_neighbor (Point point, scalar * scalars, vector * vectors)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  for (int k = 1; k <= 2; k++)
    {
      for (int i = -k; i <= k; i += 2*k)
 if ((allocated(i,0,0) && !(neighbor(i,0,0).pid < 0))) {
   Point neighbor = neighborp(i,0,0);
   int id = (- cell.pid - 1);
   {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    
       val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s, NULL);}}
   {vector*_i=(vector*)( vectors);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
     {
       scalar vn = VN;
       val(v.x,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.x, NULL);

       scalar vt = VT;
       val(v.y,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.y, NULL);





     }}}
   return true;
 }
      
#line 1152
for (int i = -k; i <= k; i += 2*k)
 if ((allocated(0,i,0) && !(neighbor(0,i,0).pid < 0))) {
   Point neighbor = neighborp(0,i,0);
   int id = (- cell.pid - 1);
   {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    
       val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s, NULL);}}
   {vector*_i=(vector*)( vectors);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
     {
       scalar vn = VN;
       val(v.y,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.y, NULL);

       scalar vt = VT;
       val(v.x,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.x, NULL);





     }}}
   return true;
 }}
  return false;
}

static bool diagonal_neighbor_2D (Point point,
      scalar * scalars, vector * vectors)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  for (int k = 1; k <= 2; k++)



      for (int i = -k; i <= k; i += 2*k)
 for (int j = -k; j <= k; j += 2*k)
   if (allocated(i,j,0) && (allocated(i,j,0) && !(neighbor(i,j,0).pid < 0)) &&
       allocated(i,0,0) && (neighbor(i,0,0).pid < 0) &&
       allocated(0,j,0) && (neighbor(0,j,0).pid < 0)) {
     Point n = neighborp(i,j,0),
       n1 = neighborp(i,0,0), n2 = neighborp(0,j,0);
     int id1 = (- neighbor(i,0,0).pid - 1), id2 = (- neighbor(0,j,0).pid - 1);
     {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      
  val(s,0,0,0) = (_attribute[s.i].boundary[id1](n,n1,s,NULL) +
         _attribute[s.i].boundary[id2](n,n2,s,NULL) -
         val(s,i,j,0));}}
     {vector*_i=(vector*)( vectors);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
       {
  scalar vt = VT, vn = VN;
  val(v.x,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.x,NULL) +
    _attribute[vn.i].boundary[id2](n,n2,v.x,NULL) -
    val(v.x,i,j,0));
  val(v.y,0,0,0) = (_attribute[vn.i].boundary[id1](n,n1,v.y,NULL) +
    _attribute[vt.i].boundary[id2](n,n2,v.y,NULL) -
    val(v.y,i,j,0));






       }}}
     return true;
   }

  return false;
}

static bool diagonal_neighbor_3D (Point point,
      scalar * scalars, vector * vectors)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
#line 1266 "/home/gongminjiang/basilisk/src/grid/tree.h"
  return false;
}



static Point tangential_neighbor_x (Point point, bool * zn)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  for (int k = 1; k <= 2; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(0,j,0) && !(neighbor(0,j,0).pid < 0)) || (allocated(-1,j,0) && !(neighbor(-1,j,0).pid < 0))) {
 *zn = false;
 return neighborp(0,j,0);
      }







    }
  return (Point){.level = -1};
}

#line 1271
static Point tangential_neighbor_y (Point point, bool * zn)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  for (int k = 1; k <= 2; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(j,0,0) && !(neighbor(j,0,0).pid < 0)) || (allocated(j,-1,0) && !(neighbor(j,-1,0).pid < 0))) {
 *zn = false;
 return neighborp(j,0,0);
      }







    }
  return (Point){.level = -1};
}


static inline bool is_boundary_point (Point point) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return (cell.pid < 0);
}

static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  disable_fpe_for_mpi();
  scalar * scalars = NULL;
  vector * vectors = NULL, * faces = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].v.x.i == s.i) {
 if (_attribute[s.i].face)
   faces = vectors_add (faces, _attribute[s.i].v);
 else
   vectors = vectors_add (vectors, _attribute[s.i].v);
      }
      else if (_attribute[s.i].v.x.i < 0 && _attribute[s.i].boundary[0])
 scalars = list_add (scalars, s);
    }}}

  {foreach_boundary_level (l) {
    if (!normal_neighbor (point, scalars, vectors) &&
 !diagonal_neighbor_2D (point, scalars, vectors) &&
 !diagonal_neighbor_3D (point, scalars, vectors)) {

      {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){

   val(s,0,0,0) = undefined;}}
      {vector*_i=(vector*)( vectors);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){

   {
     val(v.x,0,0,0) = undefined;
     
#line 1323
val(v.y,0,0,0) = undefined;}}}
    }
    if (faces) {
      int id = (- cell.pid - 1);
      
 for (int i = -1; i <= 1; i += 2) {

   if ((allocated(i,0,0) && !(neighbor(i,0,0).pid < 0))) {
     Point neighbor = neighborp(i,0,0);
     {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ {
       scalar vn = VN;
       if (_attribute[vn.i].boundary[id])
 
    val(v.x,(i + 1)/2,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.x, NULL);
     }}}
   }

   else if (i == -1) {

     bool zn;
     Point neighbor = tangential_neighbor_x (point, &zn);
     if (neighbor.level >= 0) {
       int id = is_boundary_point (neighbor) ?
  (- neighbor(-1,0,0).pid - 1) : (- cell.pid - 1);
       {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ {

  scalar vt = VT;



 
    val(v.x,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.x, NULL);
       }}}
     }
     else

       {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
 
    val(v.x,0,0,0) = 0.;}}
   }

 }
 
#line 1328
for (int i = -1; i <= 1; i += 2) {

   if ((allocated(0,i,0) && !(neighbor(0,i,0).pid < 0))) {
     Point neighbor = neighborp(0,i,0);
     {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ {
       scalar vn = VN;
       if (_attribute[vn.i].boundary[id])
 
    val(v.y,0,(i + 1)/2,0) = _attribute[vn.i].boundary[id](neighbor, point, v.y, NULL);
     }}}
   }

   else if (i == -1) {

     bool zn;
     Point neighbor = tangential_neighbor_y (point, &zn);
     if (neighbor.level >= 0) {
       int id = is_boundary_point (neighbor) ?
  (- neighbor(0,-1,0).pid - 1) : (- cell.pid - 1);
       {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ {

  scalar vt = VT;



 
    val(v.y,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.y, NULL);
       }}}
     }
     else

       {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
 
    val(v.y,0,0,0) = 0.;}}
   }

 }
    }
  }end_foreach_boundary_level();}

  pfree (scalars,__func__,__FILE__,__LINE__);
  pfree (vectors,__func__,__FILE__,__LINE__);
  pfree (faces,__func__,__FILE__,__LINE__);
  enable_fpe_for_mpi();
}



#undef VN
#undef VT
#define VN _attribute[s.i].v.x
#define VT _attribute[s.i].v.y

static double masked_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double sum = 0., n = 0.;
  {foreach_child()
    if (!(cell.pid < 0) && val(s,0,0,0) != HUGE)
      sum += val(s,0,0,0), n++;end_foreach_child()}
  return n ? sum/n : HUGE;
}


static double masked_average_x (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double sum = 0., n = 0.;
  {foreach_child()
    if (child.x < 0 && (!(cell.pid < 0) || !(neighbor(1,0,0).pid < 0)) &&
 val(s,1,0,0) != HUGE)
      sum += val(s,1,0,0), n++;end_foreach_child()}
  return n ? sum/n : HUGE;
}

#line 1391
static double masked_average_y (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double sum = 0., n = 0.;
  {foreach_child()
    if (child.y < 0 && (!(cell.pid < 0) || !(neighbor(0,1,0).pid < 0)) &&
 val(s,0,1,0) != HUGE)
      sum += val(s,0,1,0), n++;end_foreach_child()}
  return n ? sum/n : HUGE;
}

static void masked_boundary_restriction (const Boundary * b,
      scalar * list, int l)
{
  scalar * scalars = NULL;
  vector * faces = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].v.x.i == s.i && _attribute[s.i].face)
 faces = vectors_add (faces, _attribute[s.i].v);
      else
 scalars = list_add (scalars, s);
    }}}

  {foreach_halo (restriction, l) {
    {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      val(s,0,0,0) = masked_average (parent, s);}}
    {vector*_i=(vector*)( faces);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
      { {
 double average = masked_average_x (parent, v.x);
 if ((neighbor(-1,0,0).pid < 0))
   val(v.x,0,0,0) = average;
 if ((neighbor(1,0,0).pid < 0))
   val(v.x,1,0,0) = average;
      } 
#line 1418
{
 double average = masked_average_y (parent, v.y);
 if ((neighbor(0,-1,0).pid < 0))
   val(v.y,0,0,0) = average;
 if ((neighbor(0,1,0).pid < 0))
   val(v.y,0,1,0) = average;
      }}}}
  }end_foreach_halo();}

  pfree (scalars,__func__,__FILE__,__LINE__);
  pfree (faces,__func__,__FILE__,__LINE__);
}
#line 1454 "/home/gongminjiang/basilisk/src/grid/tree.h"
static void free_cache (CacheLevel * c)
{
  for (int l = 0; l <= depth(); l++)
    pfree (c[l].p,__func__,__FILE__,__LINE__);
  pfree (c,__func__,__FILE__,__LINE__);
}

void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  Tree * q = ((Tree *)grid);
  pfree (q->leaves.p,__func__,__FILE__,__LINE__);
  pfree (q->faces.p,__func__,__FILE__,__LINE__);
  pfree (q->vertices.p,__func__,__FILE__,__LINE__);
  pfree (q->refined.p,__func__,__FILE__,__LINE__);


  Layer * L = q->L[0];
  {foreach_mem (L->m, L->len, 1) {



    pfree (((L->m)->b[point.i][point.j]),__func__,__FILE__,__LINE__);



  }end_foreach_mem();}
  for (int l = 0; l <= depth(); l++)
    destroy_layer (q->L[l]);
  q->L = &(q->L[-1]);
  pfree (q->L,__func__,__FILE__,__LINE__);
  free_cache (q->active);
  free_cache (q->prolongation);
  free_cache (q->boundary);
  free_cache (q->restriction);
  pfree (q,__func__,__FILE__,__LINE__);
  grid = NULL;
}

static void refine_level (int depth);

     
void init_grid (int n)
{tracing("init_grid","/home/gongminjiang/basilisk/src/grid/tree.h",1498);

  if (!(sizeof(Cell) % 8 == 0)) qassert ("/home/gongminjiang/basilisk/src/grid/tree.h", 1501, "sizeof(Cell) % 8 == 0");

  free_grid();
  int depth = 0;
  while (n > 1) {
    if (n % 2) {
      fprintf (ferr, "tree: N must be a power-of-two\n");
      exit (1);
    }
    n /= 2;
    depth++;
  }
  Tree * q = ((Tree *) pcalloc (1, sizeof(Tree),__func__,__FILE__,__LINE__));
  grid = (Grid *) q;
  grid->depth = 0;


  q->L = ((Layer * *) pmalloc ((2)*sizeof(Layer *),__func__,__FILE__,__LINE__));

  q->L[0] = NULL; q->L = &(q->L[1]);

  Layer * L = new_layer (0);
  q->L[0] = L;
#line 1537 "/home/gongminjiang/basilisk/src/grid/tree.h"
  for (int i = Period.x*2; i < L->len - Period.x*2; i++)
    for (int j = Period.y*2; j < L->len - Period.y*2; j++)
      assign_periodic (L->m, i, j, L->len,
         (char *) pcalloc (1, sizeof(Cell) + datasize,__func__,__FILE__,__LINE__));
  CELL(((L->m)->b[2][2])).flags |= leaf;
  if (pid() == 0)
    CELL(((L->m)->b[2][2])).flags |= active;
  for (int k = - 2*(1 - Period.x); k <= 2*(1 - Period.x); k++)
    for (int l = -2*(1 - Period.y); l <= 2*(1 - Period.y); l++)
      CELL(((L->m)->b[2 +k][2 +l])).pid =
 (k < 0 ? -1 - left :
  k > 0 ? -1 - right :
  l > 0 ? -1 - top :
  l < 0 ? -1 - bottom :
  0);
  CELL(((L->m)->b[2][2])).pid = 0;
#line 1575 "/home/gongminjiang/basilisk/src/grid/tree.h"
  q->active = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->prolongation = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->boundary = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->restriction = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->dirty = true;
  N = 1 << depth;
#if 1
  void mpi_boundary_new();
  mpi_boundary_new();
#endif

  Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,__LINE__));
  b->level = box_boundary_level;
  b->restriction = masked_boundary_restriction;
  add_boundary (b);
  refine_level (depth);
  reset (all, 0.);
  { if (((Tree *)grid)->dirty) update_cache_f(); };
end_tracing("init_grid","/home/gongminjiang/basilisk/src/grid/tree.h",1593);}


void check_two_one (void)
{
  {foreach_leaf()
    if (level > 0)
      for (int k = -1; k <= 1; k++)
 for (int l = -1; l <= 1; l++) {

   int i = (point.i + 2)/2 + k;
   int j = (point.j + 2)/2 + l;
   double x = ((i - 2 + 0.5)*(1./(1 << point.level))*2. - 0.5);
   double y = ((j - 2 + 0.5)*(1./(1 << point.level))*2. - 0.5);
   if (x > -0.5 && x < 0.5 && y > -0.5 && y < 0.5 &&
       !(aparent(k,l,0).flags & active)) {
     FILE * fp = fopen("check_two_one_loc", "w");
     fprintf (fp,
       "# %d %d\n"
       "%g %g\n%g %g\n",
       k, l,
       (((point.i - 2) + 0.5)*(1./(1 << point.level)) - 0.5),
       (((point.j - 2) + 0.5)*(1./(1 << point.level)) - 0.5),
       x, y);
     fclose (fp);





     if (!(false)) qassert ("/home/gongminjiang/basilisk/src/grid/tree.h", 1623, "false");
   }
 }end_foreach_leaf();}
}


struct _locate { double x, y, z; };

Point locate (struct _locate p)
{
  for (int l = depth(); l >= 0; l--) {
    Point point = {0};int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
    point.level = l;
    int n = 1 << point.level;
    point.i = (p.x - X0)/L0*n + 2;

    point.j = (p.y - Y0)/L0*n + 2;




    if (point.i >= 0 && point.i < n + 2*2

 && point.j >= 0 && point.j < n + 2*2




 ) {
      if (allocated(0,0,0) && is_local(cell) && is_leaf(cell))
 return point;
    }
    else
      break;
  }
  Point point = {0};int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  point.level = -1;
  return point;
}



bool tree_is_full()
{
  { if (((Tree *)grid)->dirty) update_cache_f(); };
  return (grid->tn == 1L << grid->maxdepth*2);
}

#line 1 "grid/tree-common.h"
#line 1 "/home/gongminjiang/basilisk/src/grid/tree-common.h"



#line 1 "grid/multigrid-common.h"
#line 1 "/home/gongminjiang/basilisk/src/grid/multigrid-common.h"


#line 1 "grid/cartesian-common.h"
#line 1 "/home/gongminjiang/basilisk/src/grid/cartesian-common.h"
#line 1 "grid/events.h"
#line 1 "/home/gongminjiang/basilisk/src/grid/events.h"







static void event_error (Event * ev, const char * s)
{
  fprintf (ferr, "%s:%d: error: %s\n", ev->file, ev->line, s);
  exit (1);
}

static void init_event (Event * ev)
{
  if (ev->arrayi || ev->arrayt) {
    ev->i = ev->t = -1;
    if (ev->arrayi)
      ev->i = ev->arrayi[0];
    else
      ev->t = ev->arrayt[0];
    ev->a = 1;
    ev->expr[1] = NULL;
  }
  else {
    if (ev->nexpr > 0) {
      Expr init = NULL, cond = NULL, inc = NULL;
      for (int j = 0; j < ev->nexpr; j++) {
 int i = -123456; double t = i;
 (* ev->expr[j]) (&i, &t, ev);
 if (i == -123456 && t == -123456) {

   if (cond)
     event_error (ev, "events can only use a single condition");
   cond = ev->expr[j];
 }
 else {

   int i1 = i; double t1 = t;
   (* ev->expr[j]) (&i1, &t1, ev);
   if (i1 == i && t1 == t) {


     if (init)
       event_error (ev, "events can only use a single initialisation");
     init = ev->expr[j];
   }
   else {

     if (inc)
       event_error (ev, "events can only use a single increment");
     inc = ev->expr[j];
   }
 }
      }
      ev->expr[0] = init;
      ev->expr[1] = cond;
      ev->expr[2] = inc;
      ev->nexpr = 0;
    }
    ev->i = ev->t = -1;
    if (ev->expr[0]) {
      (* ev->expr[0]) (&ev->i, &ev->t, ev);
      if (ev->i == 1234567890 || ev->t == 1234567890) {
 ev->i = 1234567890; ev->t = -1;
      }
    }
    else if (ev->expr[2]) {
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (ev->i != -1)
 ev->i = 0;
      if (ev->t != -1)
 ev->t = 0;
    }
  }
}

enum { event_done, event_alive, event_stop };

static int event_finished (Event * ev)
{
  ev->t = ev->i = -1;
  return event_done;
}

void event_register (Event event) {
  if (!(Events)) qassert ("/home/gongminjiang/basilisk/src/grid/events.h", 87, "Events");
  if (!(!event.last)) qassert ("/home/gongminjiang/basilisk/src/grid/events.h", 88, "!event.last");
  int n = 0, parent = -1;
  for (Event * ev = Events; !ev->last; ev++) {
    if (!strcmp (event.name, ev->name)) {
      if (!(parent < 0)) qassert ("/home/gongminjiang/basilisk/src/grid/events.h", 92, "parent < 0");
      parent = n;
    }
    n++;
  }
  if (parent < 0) {
    Events = (Event *) prealloc (Events, (n + 2)*sizeof(Event),__func__,__FILE__,__LINE__);
    Events[n] = event;
    Events[n].next = NULL;
    Events[n + 1].last = true;
    init_event (&Events[n]);
  }
  else {
    Event * ev = ((Event *) pcalloc (1, sizeof(Event),__func__,__FILE__,__LINE__));
    *ev = Events[parent];
    Events[parent] = event;
    Events[parent].next = ev;
    init_event (&Events[parent]);
  }
}

static int event_cond (Event * ev, int i, double t)
{
  if (!ev->expr[1])
    return true;
  return (* ev->expr[1]) (&i, &t, ev);
}
#line 131 "/home/gongminjiang/basilisk/src/grid/events.h"
static int event_do (Event * ev, bool action)
{
  if ((iter > ev->i && t > ev->t) || !event_cond (ev, iter, t))
    return event_finished (ev);
  if (iter == ev->i || fabs (t - ev->t) <= 1e-9) {
    if (action) {
      bool finished = false;
      for (Event * e = ev; e; e = e->next) {



 if ((* e->action) (iter, t, e))
   finished = true;
      }
      if (finished) {
 event_finished (ev);
 return event_stop;
      }
    }
    if (ev->arrayi) {
      ev->i = ev->arrayi[ev->a++];
      if (ev->i < 0)
 return event_finished (ev);
    }
    if (ev->arrayt) {
      ev->t = ev->arrayt[ev->a++];
      if (ev->t < 0)
 return event_finished (ev);
    }
    else if (ev->expr[2]) {
      int i0 = ev->i;
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (i0 == -1 && ev->i != i0)
 ev->i += iter + 1;
      if (!event_cond (ev, iter + 1, ev->t))
 return event_finished (ev);
    }
    else if (ev->expr[0] && !ev->expr[1])
      return event_finished (ev);
  }
  return event_alive;
}

static void end_event_do (bool action)
{




  for (Event * ev = Events; !ev->last; ev++)
    if (ev->i == 1234567890 && action)
      for (Event * e = ev; e; e = e->next) {



 e->action (iter, t, e);
      }
}

int events (bool action)
{





  if (iter == 0)
    for (Event * ev = Events; !ev->last; ev++)
      init_event (ev);

  int cond = 0, cond1 = 0;
  inext = 1234567890; tnext = HUGE;
  for (Event * ev = Events; !ev->last && !cond; ev++)
    if (ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond = 1;
  for (Event * ev = Events; !ev->last; ev++) {
    int status = event_do (ev, action);
    if (status == event_stop) {
      end_event_do (action);
      return 0;
    }
    if (status == event_alive && ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond1 = 1;
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
    if (ev->i > iter && ev->i < inext)
      inext = ev->i;
  }
  if ((!cond || cond1) && (tnext != HUGE || inext != 1234567890)) {
    inext = iter + 1;
    return 1;
  }
  end_event_do (action);
  return 0;
}

void event (const char * name)
{
  for (Event * ev = Events; !ev->last; ev++)
    if (!strcmp (ev->name, name))
      for (Event * e = ev; e; e = e->next) {



 (* e->action) (0, 0, e);
      }
}

double dtnext (double dt)
{
  if (tnext != HUGE && tnext > t) {
    unsigned int n = (tnext - t)/dt;
    if (!(n < INT_MAX)) qassert ("/home/gongminjiang/basilisk/src/grid/events.h", 245, "n < INT_MAX");
    if (n == 0)
      dt = tnext - t;
    else {
      double dt1 = (tnext - t)/n;
      if (dt1 > dt + 1e-9)
 dt = (tnext - t)/(n + 1);
      else if (dt1 < dt)
 dt = dt1;
      tnext = t + dt;
    }
  }
  else
    tnext = t + dt;
  return dt;
}
#line 2 "/home/gongminjiang/basilisk/src/grid/cartesian-common.h"

void (* debug) (Point);

#define _val_constant(a,k,l,m) ((const double) _constant[a.i -_NVARMAX])
#define diagonalize(a)
#define val_diagonal(a,k,l,m) ((k) == 0 && (l) == 0 && (m) == 0)

#undef VARIABLES
#define VARIABLES\
  double Delta = L0*(1./(1 << point.level));\
  double Delta_x = Delta;\
\
  double Delta_y = Delta;\
\
\
\
\
\
  double x = (ig/2. + (point.i - 2) + 0.5)*Delta + X0; NOT_UNUSED(x);\
\
  double y = (jg/2. + (point.j - 2) + 0.5)*Delta + Y0;\
\
\
\
 NOT_UNUSED(y);\
\
\
\
  double z = 0.;\
\
  NOT_UNUSED(z);\
\
  NOT_UNUSED(Delta);\
  NOT_UNUSED(Delta_x);\
\
  NOT_UNUSED(Delta_y);\
\
\
\
\
\
  ;\

#line 44


#line 1 "grid/fpe.h"
#line 1 "/home/gongminjiang/basilisk/src/grid/fpe.h"


#include <signal.h>
#include <unistd.h>

static int gdb()
{
  if (last_point.level >= 0) {
    debug (last_point);
    fputc ('\n', ferr);
    fflush (ferr);
  }
  char command[80];
  sprintf (command, "exec xterm -e 'gdb -p %d' & xterm -e 'gnuplot plot -'",
    getpid());
  return system (command);
}

static void caught_abort (int sig)
{
  fprintf (ferr, "Caught signal %d (Aborted)\n", sig);
  gdb();
}

static void caught_fpe (int sig)
{
  fprintf (ferr, "Caught signal %d (Floating Point Exception)\n", sig);
  gdb();
  exit (1);
}

static void caught_segfault (int sig)
{
  fprintf (ferr, "Caught signal %d (Segmentation Fault)\n", sig);
  gdb();
  exit (2);
}

void catch_fpe (void)
{
  struct sigaction act;
  act.sa_handler = caught_fpe;
  sigemptyset (&act.sa_mask);
  act.sa_flags = 0;
  last_point.level = -1;
  sigaction (8, &act, NULL);
  act.sa_handler = caught_segfault;
  sigaction (11, &act, NULL);
  act.sa_handler = caught_abort;
  act.sa_flags = SA_RESETHAND;
  sigaction (6, &act, NULL);
}
#line 47 "/home/gongminjiang/basilisk/src/grid/cartesian-common.h"

#define end_foreach_face()



static void init_block_scalar (scalar sb, const char * name, const char * ext,
          int n, int block)
{
  char bname[strlen(name) + strlen(ext) + 10];
  if (n == 0) {
    sprintf (bname, "%s%s", name, ext);
    _attribute[sb.i].block = block;
    init_scalar (sb, bname);
    baseblock = list_append (baseblock, sb);
  }
  else {
    sprintf (bname, "%s%d%s", name, n, ext);
    _attribute[sb.i].block = - n;
    init_scalar (sb, bname);
  }
  all = list_append (all, sb);
}

scalar new_block_scalar (const char * name, const char * ext, int block)
{
  int nvar = datasize/sizeof(double);

  scalar s = {0};
  while (s.i < nvar) {
    int n = 0;
    scalar sb = s;
    while (sb.i < nvar && n < block && _attribute[sb.i].freed)
      n++, sb.i++;
    if (n >= block) {
      for (sb.i = s.i, n = 0; n < block; n++, sb.i++)
 init_block_scalar (sb, name, ext, n, block);
      trash (((scalar []){s, {-1}}));
      return s;
    }
    s.i = sb.i + 1;
  }


  s = (scalar){nvar};
  if (!(nvar + block <= _NVARMAX)) qassert ("/home/gongminjiang/basilisk/src/grid/cartesian-common.h", 91, "nvar + block <= _NVARMAX");
  _attribute = (_Attributes *) prealloc (_attribute, (nvar + block)*sizeof(_Attributes),__func__,__FILE__,__LINE__);
  memset (&_attribute[nvar], 0, block*sizeof (_Attributes));
  for (int n = 0; n < block; n++, nvar++) {
    scalar sb = (scalar){nvar};
    init_block_scalar (sb, name, ext, n, block);
  }

  realloc_scalar (block*sizeof(double));
  trash (((scalar []){s, {-1}}));
  return s;
}

scalar new_scalar (const char * name)
{
  return new_block_scalar (name, "", 1);
}

scalar new_vertex_scalar (const char * name)
{
  scalar s = new_block_scalar (name, "", 1);
  init_vertex_scalar (s, NULL);
  return s;
}

scalar new_block_vertex_scalar (const char * name, int block)
{
  scalar s = new_block_scalar (name, "", block);
  for (int i = 0; i < block; i++) {
    scalar sb = {s.i + i};
    init_vertex_scalar (sb, NULL);
  }
  return s;
}

static vector alloc_block_vector (const char * name, int block)
{
  vector v;
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  
    v.x = new_block_scalar (name, ext.x, block);
    
#line 131
v.y = new_block_scalar (name, ext.y, block);
  return v;
}

vector new_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_vector (v, NULL);
  return v;
}

vector new_face_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_face_vector (v, NULL);
  return v;
}

vector new_block_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    
      vb.x.i = v.x.i + i;
      
#line 155
vb.y.i = v.y.i + i;
    init_vector (vb, NULL);
    
      _attribute[vb.x.i].block = - i;
      
#line 158
_attribute[vb.y.i].block = - i;
  }
  
    _attribute[v.x.i].block = block;
    
#line 161
_attribute[v.y.i].block = block;
  return v;
}

vector new_block_face_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    
      vb.x.i = v.x.i + i;
      
#line 171
vb.y.i = v.y.i + i;
    init_face_vector (vb, NULL);
    
      _attribute[vb.x.i].block = - i;
      
#line 174
_attribute[vb.y.i].block = - i;
  }
  
    _attribute[v.x.i].block = block;
    
#line 177
_attribute[v.y.i].block = block;
  return v;
}

tensor new_tensor (const char * name)
{
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {"%s.x", "%s.y", "%s.z"};
  tensor t;
   {
    sprintf (cname, ext.x, name);
    t.x = new_vector (cname);
  } 
#line 186
{
    sprintf (cname, ext.y, name);
    t.y = new_vector (cname);
  }
  init_tensor (t, NULL);
  return t;
}

tensor new_symmetric_tensor (const char * name)
{
  char cname[strlen(name) + 5];
  struct { char * x, * y, * z; } ext = {"%s.x.x", "%s.y.y", "%s.z.z"};
  tensor t;
   {
    sprintf (cname, ext.x, name);
    t.x.x = new_scalar(cname);
  } 
#line 199
{
    sprintf (cname, ext.y, name);
    t.y.y = new_scalar(cname);
  }

    sprintf (cname, "%s.x.y", name);
    t.x.y = new_scalar(cname);
    t.y.x = t.x.y;
#line 219 "/home/gongminjiang/basilisk/src/grid/cartesian-common.h"
  init_tensor (t, NULL);
  return t;
}

static int nconst = 0;

void init_const_scalar (scalar s, const char * name, double val)
{
  if (s.i - _NVARMAX >= nconst) {
    nconst = s.i - _NVARMAX + 1;
    _constant = (double *) prealloc (_constant, (nconst)*sizeof(double),__func__,__FILE__,__LINE__);
  }
  _constant[s.i - _NVARMAX] = val;
}

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = (scalar){i + _NVARMAX};
  init_const_scalar (s, name, val);
  return s;
}

void init_const_vector (vector v, const char * name, double * val)
{
  
    init_const_scalar (v.x, name, *val++);
    
#line 244
init_const_scalar (v.y, name, *val++);
}

vector new_const_vector (const char * name, int i, double * val)
{
  vector v;
  
    v.x.i = _NVARMAX + i++;
    
#line 251
v.y.i = _NVARMAX + i++;
  init_const_vector (v, name, val);
  return v;
}

void scalar_clone (scalar a, scalar b)
{
  char * name = _attribute[a.i].name;
  double (** boundary) (Point, Point, scalar, void *) = _attribute[a.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar, void *) =
    _attribute[a.i].boundary_homogeneous;
  if (!(_attribute[b.i].block > 0 && _attribute[a.i].block == _attribute[b.i].block)) qassert ("/home/gongminjiang/basilisk/src/grid/cartesian-common.h", 262, "b.block > 0 && a.block == b.block");
  pfree (_attribute[a.i].depends,__func__,__FILE__,__LINE__);
  _attribute[a.i] = _attribute[b.i];
  _attribute[a.i].name = name;
  _attribute[a.i].boundary = boundary;
  _attribute[a.i].boundary_homogeneous = boundary_homogeneous;
  for (int i = 0; i < nboundary; i++) {
    _attribute[a.i].boundary[i] = _attribute[b.i].boundary[i];
    _attribute[a.i].boundary_homogeneous[i] = _attribute[b.i].boundary_homogeneous[i];
  }
  _attribute[a.i].depends = list_copy (_attribute[b.i].depends);
}

scalar * list_clone (scalar * l)
{
  scalar * list = NULL;
  int nvar = datasize/sizeof(double), map[nvar];
  for (int i = 0; i < nvar; i++)
    map[i] = -1;
  {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    scalar c = _attribute[s.i].block > 1 ? new_block_scalar("c", "", _attribute[s.i].block) :
      new_scalar("c");
    scalar_clone (c, s);
    map[s.i] = c.i;
    list = list_append (list, c);
  }}}
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    {
      if (_attribute[s.i].v.x.i >= 0 && map[_attribute[s.i].v.x.i] >= 0)
 _attribute[s.i].v.x.i = map[_attribute[s.i].v.x.i];
      
#line 290
if (_attribute[s.i].v.y.i >= 0 && map[_attribute[s.i].v.y.i] >= 0)
 _attribute[s.i].v.y.i = map[_attribute[s.i].v.y.i];}}}
  return list;
}

void delete (scalar * list)
{
  if (all == NULL)
    return;

  {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){ {
    for (int i = 0; i < _attribute[f.i].block; i++) {
      scalar fb = {f.i + i};
      if (_attribute[f.i].delete)
 _attribute[f.i].delete (fb);
      pfree (_attribute[fb.i].name,__func__,__FILE__,__LINE__); _attribute[fb.i].name = NULL;
      pfree (_attribute[fb.i].boundary,__func__,__FILE__,__LINE__); _attribute[fb.i].boundary = NULL;
      pfree (_attribute[fb.i].boundary_homogeneous,__func__,__FILE__,__LINE__); _attribute[fb.i].boundary_homogeneous = NULL;
      pfree (_attribute[fb.i].depends,__func__,__FILE__,__LINE__); _attribute[fb.i].depends = NULL;
      _attribute[fb.i].freed = true;
    }
  }}}

  if (list == all) {
    all[0].i = -1;
    baseblock[0].i = -1;
    return;
  }

  trash (list);
  {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){ {
    if (_attribute[f.i].block > 0) {
      scalar * s = all;
      for (; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
 for (; s[_attribute[f.i].block].i >= 0; s++)
   s[0] = s[_attribute[f.i].block];
 s->i = -1;
      }
      for (s = baseblock; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
 for (; s[1].i >= 0; s++)
   s[0] = s[1];
 s->i = -1;
      }
    }
  }}}
}

void free_solver()
{
  if (!(_val_higher_dimension == 0.)) qassert ("/home/gongminjiang/basilisk/src/grid/cartesian-common.h", 341, "_val_higher_dimension == 0.");

  if (free_solver_funcs) {
    free_solver_func * a = (free_solver_func *) free_solver_funcs->p;
    for (int i = 0; i < free_solver_funcs->len/sizeof(free_solver_func); i++)
      a[i] ();
    array_free (free_solver_funcs);
  }

  delete (all);
  pfree (all,__func__,__FILE__,__LINE__); all = NULL;
  pfree (baseblock,__func__,__FILE__,__LINE__); baseblock = NULL;
  for (Event * ev = Events; !ev->last; ev++) {
    Event * e = ev->next;
    while (e) {
      Event * next = e->next;
      pfree (e,__func__,__FILE__,__LINE__);
      e = next;
    }
  }

  pfree (Events,__func__,__FILE__,__LINE__); Events = NULL;
  pfree (_attribute,__func__,__FILE__,__LINE__); _attribute = NULL;
  pfree (_constant,__func__,__FILE__,__LINE__); _constant = NULL;
  free_grid();
  qpclose_all();
#if TRACE
  trace_off();
#endif
#if MTRACE
  pmuntrace();
#endif
#if _CADNA
  cadna_end();
#endif
}



void (* boundary_level) (scalar *, int l);
void (* boundary_face) (vectorl);




void boundary_flux (vector * list) __attribute__ ((deprecated));

void boundary_flux (vector * list)
{
  vectorl list1 = {NULL};
  {vector*_i=(vector*)( list);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
    {
      list1.x = list_append (list1.x, v.x);
      
#line 393
list1.y = list_append (list1.y, v.y);}}}
  boundary_face (list1);
  
    pfree (list1.x,__func__,__FILE__,__LINE__);
    
#line 396
pfree (list1.y,__func__,__FILE__,__LINE__);
}

static scalar * list_add_depends (scalar * list, scalar s)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){
    if (t.i == s.i)
      return list;}}
  scalar * list1 = list;
  {scalar*_i=(scalar*)( _attribute[s.i].depends);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
    if (_attribute[d.i].dirty)
      list1 = list_add_depends (list1, d);}}
  return list_append (list1, s);
}

     
void boundary_internal (scalar * list, const char * fname, int line)
{tracing("boundary_internal","/home/gongminjiang/basilisk/src/grid/cartesian-common.h",412);
  if (list == NULL)
    {end_tracing("boundary_internal","/home/gongminjiang/basilisk/src/grid/cartesian-common.h",415);return;}
  scalar * listc = NULL;
  vectorl listf = {NULL};
  bool flux = false;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].block > 0) {
      if (scalar_is_dirty (s)) {
 if (_attribute[s.i].face && _attribute[s.i].dirty != 2)
   {
     if (_attribute[s.i].v.x.i == s.i)
       listf.x = list_add (listf.x, s), flux = true;
     
#line 424
if (_attribute[s.i].v.y.i == s.i)
       listf.y = list_add (listf.y, s), flux = true;}
 if (!is_constant(cm) && _attribute[cm.i].dirty)
   listc = list_add_depends (listc, cm);
 if (_attribute[s.i].face != 2)
   listc = list_add_depends (listc, s);
      }




    }}}
  if (flux) {
    boundary_face (listf);
    
      pfree (listf.x,__func__,__FILE__,__LINE__);
      
#line 439
pfree (listf.y,__func__,__FILE__,__LINE__);
  }
  if (listc) {
    boundary_level (listc, -1);
    {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = false;}}
    pfree (listc,__func__,__FILE__,__LINE__);
  }
end_tracing("boundary_internal","/home/gongminjiang/basilisk/src/grid/cartesian-common.h",447);}

void cartesian_boundary_level (scalar * list, int l)
{
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, l); };
}

void cartesian_boundary_face (vectorl list)
{
  
    {scalar*_i=(scalar*)( list.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
    
#line 457
{scalar*_i=(scalar*)( list.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
}

static double symmetry (Point point, Point neighbor, scalar s, void * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return val(s,0,0,0);
}

static double antisymmetry (Point point, Point neighbor, scalar s, void * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return -val(s,0,0,0);
}

double (* default_scalar_bc[]) (Point, Point, scalar, void *) = {
  symmetry, symmetry, symmetry, symmetry, symmetry, symmetry
};

scalar cartesian_init_scalar (scalar s, const char * name)
{

  char * pname;
  if (name) {
    pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
    pname = pstrdup (name,__func__,__FILE__,__LINE__);
  }
  else
    pname = _attribute[s.i].name;
  int block = _attribute[s.i].block;
  double (** boundary) (Point, Point, scalar, void *) = _attribute[s.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar, void *) =
    _attribute[s.i].boundary_homogeneous;

  _attribute[s.i] = (const _Attributes){0};
  _attribute[s.i].name = pname;
  _attribute[s.i].block = block == 0 ? 1 : block;

  _attribute[s.i].boundary = boundary ? boundary :
    (double (**)(Point, Point, scalar, void *))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  _attribute[s.i].boundary_homogeneous = boundary_homogeneous ? boundary_homogeneous :
    (double (**)(Point, Point, scalar, void *))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  for (int b = 0; b < nboundary; b++)
    _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] =
      b < 2*2 ? default_scalar_bc[b] : symmetry;
  _attribute[s.i].gradient = NULL;
   {
    _attribute[s.i].d.x = 0;
    _attribute[s.i].v.x.i = -1;
  } 
#line 504
{
    _attribute[s.i].d.y = 0;
    _attribute[s.i].v.y.i = -1;
  }
  _attribute[s.i].face = false;
  return s;
}

scalar cartesian_init_vertex_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  
    _attribute[s.i].d.x = -1;
    
#line 516
_attribute[s.i].d.y = -1;
  for (int d = 0; d < nboundary; d++)
    _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = NULL;
  return s;
}

double (* default_vector_bc[]) (Point, Point, scalar, void *) = {
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry
};

vector cartesian_init_vector (vector v, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
   {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_scalar (v.x, cname);
    }
    else
      init_scalar (v.x, NULL);
    _attribute[v.x.i].v = v;
  } 
#line 531
{
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_scalar (v.y, cname);
    }
    else
      init_scalar (v.y, NULL);
    _attribute[v.y.i].v = v;
  }

  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] =
      d < 2*2 ? default_vector_bc[d] : antisymmetry;
  return v;
}

vector cartesian_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
   {
    _attribute[v.x.i].d.x = -1;
    _attribute[v.x.i].face = true;
  } 
#line 551
{
    _attribute[v.y.i].d.y = -1;
    _attribute[v.y.i].face = true;
  }
  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
  return v;
}

tensor cartesian_init_tensor (tensor t, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
   {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_vector (t.x, cname);
    }
    else
      init_vector (t.x, NULL);
  } 
#line 563
{
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_vector (t.y, cname);
    }
    else
      init_vector (t.y, NULL);
  }






    for (int b = 0; b < nboundary; b++) {
      _attribute[t.x.x.i].boundary[b] = _attribute[t.y.x.i].boundary[b] =
 _attribute[t.x.x.i].boundary_homogeneous[b] = _attribute[t.y.y.i].boundary_homogeneous[b] =
 b < 2*2 ? default_scalar_bc[b] : symmetry;
      _attribute[t.x.y.i].boundary[b] = _attribute[t.y.y.i].boundary[b] =
 _attribute[t.x.y.i].boundary_homogeneous[b] = _attribute[t.y.x.i].boundary_homogeneous[b] =
 b < 2*2 ? default_vector_bc[b] : antisymmetry;
    }



  return t;
}

struct OutputCells {
  FILE * fp;
  coord c;
  double size;
};

void output_cells (struct OutputCells p)
{
  if (!p.fp) p.fp = fout;
  {foreach() {
    bool inside = true;
    coord o = {x,y,z};
    
      if (inside && p.size > 0. &&
   (o.x > p.c.x + p.size || o.x < p.c.x - p.size))
 inside = false;
      
#line 605
if (inside && p.size > 0. &&
   (o.y > p.c.y + p.size || o.y < p.c.y - p.size))
 inside = false;
    if (inside) {
      Delta /= 2.;



      fprintf (p.fp, "%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n\n",
        x - Delta, y - Delta,
        x - Delta, y + Delta,
        x + Delta, y + Delta,
        x + Delta, y - Delta,
        x - Delta, y - Delta);
#line 633 "/home/gongminjiang/basilisk/src/grid/cartesian-common.h"
    }
  }end_foreach();}
  fflush (p.fp);
}


static void output_cells_internal (FILE * fp)
{
  output_cells ((struct OutputCells){fp});
}


static char * replace_ (const char * vname)
{
  char * name = pstrdup (vname,__func__,__FILE__,__LINE__), * c = name;
  while (*c != '\0') {
    if (*c == '.')
      *c = '_';
    c++;
  }
  return name;
}

static void debug_plot (FILE * fp, const char * name, const char * cells,
   const char * stencil)
{
  char * vname = replace_ (name);
  fprintf (fp,
    "  load 'debug.plot'\n"
    "  v=%s\n"




    "  plot '%s' w l lc 0, "
    "'%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 1 title columnhead(3+3*v)",





    vname, cells, stencil);
  pfree (vname,__func__,__FILE__,__LINE__);
}

void cartesian_debug (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  char name[80] = "cells";
  if (pid() > 0)
    sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells ((struct OutputCells){fp, (coord){x,y,z}, 4.*Delta});
  fclose (fp);

  char stencil[80] = "stencil";
  if (pid() > 0)
    sprintf (stencil, "stencil-%d", pid());
  fp = fopen (stencil, "w");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){



    fprintf (fp, "x y %s ", _attribute[v.i].name);}}



  fputc ('\n', fp);
#line 712 "/home/gongminjiang/basilisk/src/grid/cartesian-common.h"
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++) {
 {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){ {
   fprintf (fp, "%g %g ",
     x + k*Delta + _attribute[v.i].d.x*Delta/2.,
     y + l*Delta + _attribute[v.i].d.y*Delta/2.);
   if (allocated(k,l,0))
     fprintf (fp, "%g ", val(v,k,l,0));
   else
     fputs ("n/a ", fp);
 }}}
 fputc ('\n', fp);
      }
#line 742 "/home/gongminjiang/basilisk/src/grid/cartesian-common.h"
  fclose (fp);

  fp = fopen ("debug.plot", "w");
  fprintf (fp,
    "set term x11\n"
    "set size ratio -1\n"
    "set key outside\n");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    char * name = replace_ (_attribute[s.i].name);
    fprintf (fp, "%s = %d\n", name, s.i);
    pfree (name,__func__,__FILE__,__LINE__);
  }}}
  fclose (fp);

  fprintf (ferr, "Last point stencils can be displayed using (in gnuplot)\n");
  debug_plot (ferr, _attribute[0].name, name, stencil);
  fflush (ferr);

  fp = fopen ("plot", "w");
  debug_plot (fp, _attribute[0].name, name, stencil);
  fclose (fp);
}

void cartesian_methods()
{
  init_scalar = cartesian_init_scalar;
  init_vertex_scalar = cartesian_init_vertex_scalar;
  init_vector = cartesian_init_vector;
  init_tensor = cartesian_init_tensor;
  init_face_vector = cartesian_init_face_vector;
  boundary_level = cartesian_boundary_level;
  boundary_face = cartesian_boundary_face;
  debug = cartesian_debug;
}

tensor init_symmetric_tensor (tensor t, const char * name)
{
  return init_tensor (t, name);
}

struct _interpolate {
  scalar v;
  double x, y, z;
};

static double interpolate_linear (Point point, struct _interpolate p)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  scalar v = p.v;







  x = (p.x - x)/Delta - _attribute[v.i].d.x/2.;
  y = (p.y - y)/Delta - _attribute[v.i].d.y/2.;
  int i = sign(x), j = sign(y);
  x = fabs(x); y = fabs(y);

  return ((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +
   (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y);
#line 816 "/home/gongminjiang/basilisk/src/grid/cartesian-common.h"
}

     
double interpolate (struct _interpolate p)
{tracing("interpolate","/home/gongminjiang/basilisk/src/grid/cartesian-common.h",819);
  scalar v = p.v;
  boundary_internal ((scalar *)((scalar[]){v,{-1}}), "/home/gongminjiang/basilisk/src/grid/cartesian-common.h", 822);
  Point point = locate ((struct _locate){p.x, p.y, p.z});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (point.level < 0)
    {end_tracing("interpolate","/home/gongminjiang/basilisk/src/grid/cartesian-common.h",825);return HUGE;}
  { double _ret= interpolate_linear (point, p);end_tracing("interpolate","/home/gongminjiang/basilisk/src/grid/cartesian-common.h",826);return _ret;}
end_tracing("interpolate","/home/gongminjiang/basilisk/src/grid/cartesian-common.h",827);}

     
void interpolate_array (scalar * list, coord * a, int n, double * v, bool linear)
{tracing("interpolate_array","/home/gongminjiang/basilisk/src/grid/cartesian-common.h",830);
  boundary_internal ((scalar *)list, "/home/gongminjiang/basilisk/src/grid/cartesian-common.h", 832);
  int j = 0;
  for (int i = 0; i < n; i++) {
    Point point = locate ((struct _locate){a[i].x, a[i].y, a[i].z});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
    if (point.level >= 0) {
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 v[j++] = !linear ? val(s,0,0,0) :
   interpolate_linear (point,
         (struct _interpolate){s, a[i].x, a[i].y, a[i].z});}}
    }
    else
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 v[j++] = HUGE;}}
  }
#if 1
  if (pid() == 0)
    MPI_Reduce (MPI_IN_PLACE, v, n*list_len(list), MPI_DOUBLE,
  MPI_MIN, 0, MPI_COMM_WORLD);
  else
    MPI_Reduce (v, v, n*list_len(list), MPI_DOUBLE,
  MPI_MIN, 0, MPI_COMM_WORLD);
#endif
end_tracing("interpolate_array","/home/gongminjiang/basilisk/src/grid/cartesian-common.h",854);}



typedef int bid;

bid new_bid()
{
  int b = nboundary++;
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    _attribute[s.i].boundary = (double (**)(Point, Point, scalar, void *))
      prealloc (_attribute[s.i].boundary, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
    _attribute[s.i].boundary_homogeneous = (double (**)(Point, Point, scalar, void *))
      prealloc (_attribute[s.i].boundary_homogeneous, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  }}}
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    if (_attribute[s.i].v.x.i < 0)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] = symmetry;
    else if (_attribute[s.i].v.x.i == s.i) {
      vector v = _attribute[s.i].v;
      
 _attribute[v.y.i].boundary[b] = _attribute[v.y.i].boundary_homogeneous[b] = symmetry;
 
#line 875
_attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] = symmetry;
      _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] =
 _attribute[v.x.i].face ? NULL : antisymmetry;
    }
  }}}
  return b;
}



static double periodic_bc (Point point, Point neighbor, scalar s, void * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return HUGE;
}

static void periodic_boundary (int d)
{

  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (is_vertex_scalar (s))
      _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = NULL;
    else
      _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = periodic_bc;}}

  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (_attribute[s.i].face) {
      vector v = _attribute[s.i].v;
      _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
    }}}

  default_scalar_bc[d] = periodic_bc;
  default_vector_bc[d] = periodic_bc;
}

void periodic (int dir)
{



    if (!(dir <= bottom)) qassert ("/home/gongminjiang/basilisk/src/grid/cartesian-common.h", 914, "dir <= bottom");




  int c = dir/2;
  periodic_boundary (2*c);
  periodic_boundary (2*c + 1);
  (&Period.x)[c] = true;
}


double getvalue (Point point, scalar s, int i, int j, int k)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  return val(s,i,j,k);
}

void default_stencil (Point p, scalar * list)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].input = true, _attribute[s.i].width = 2;}}
}




static void write_stencil_index (int * index)
{
  fprintf (qstderr(), "[%d", index[0]);
  for (int d = 1; d < 2; d++)
    fprintf (qstderr(), ",%d", index[d]);
  fputs ("]", qstderr());
}

void stencil_val (Point p, scalar s, int i, int j, int k,
    const char * file, int line, bool overflow)
{
  if (is_constant(s) || s.i < 0)
    return;
  int index[] = {i, j, k};
  for (int d = 0; d < 2; d++)
    index[d] += (&p.i)[d];
  bool central = true;
  for (int d = 0; d < 2; d++) {
    if (!overflow && (index[d] > 2 || index[d] < - 2)) {
      fprintf (qstderr(), "%s:%d: error: stencil overflow: %s",
        file, line, _attribute[s.i].name);
      write_stencil_index (index);
      fprintf (qstderr(), "\n");
      fflush (qstderr());
      abort();
    }
    if (index[d] != 0)
      central = false;
  }
  if (central) {
    if (!_attribute[s.i].output)
      _attribute[s.i].input = true;
  }
  else {
    _attribute[s.i].input = true;
    int d = 0;
     {
      if ((!_attribute[s.i].face || _attribute[s.i].v.x.i != s.i) && abs(index[d]) > _attribute[s.i].width)
 _attribute[s.i].width = abs(index[d]);
      d++;
    } 
#line 976
{
      if ((!_attribute[s.i].face || _attribute[s.i].v.y.i != s.i) && abs(index[d]) > _attribute[s.i].width)
 _attribute[s.i].width = abs(index[d]);
      d++;
    }
  }
}

void stencil_val_a (Point p, scalar s, int i, int j, int k, bool input,
      const char * file, int line)
{
  if (is_constant(s) || s.i < 0)
    abort();
  int index[] = {i, j, k};
  for (int d = 0; d < 2; d++)
    index[d] += (&p.i)[d];
  for (int d = 0; d < 2; d++)
    if (index[d] != 0) {
      fprintf (qstderr(), "%s:%d: error: illegal write: %s",
        file, line, _attribute[s.i].name);
      write_stencil_index (index);
      fprintf (qstderr(), "\n");
      fflush (qstderr());
      abort();
    }
  if (input && !_attribute[s.i].output)
    _attribute[s.i].input = true;
  _attribute[s.i].output = true;
}
#line 4 "/home/gongminjiang/basilisk/src/grid/multigrid-common.h"

#ifndef foreach_level_or_leaf
# define foreach_level_or_leaf foreach_level
# define end_foreach_level_or_leaf end_foreach_level
#endif

#ifndef foreach_coarse_level
# define foreach_coarse_level foreach_level
# define end_foreach_coarse_level end_foreach_level
#endif










void (* restriction) (scalar *);

static inline void restriction_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double sum = 0.;
  {foreach_child()
    sum += val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 2);
}

static inline void restriction_volume_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

#line 35
if(!is_constant(cm)){{
  double sum = 0.;
  {foreach_child()
    sum += val(cm,0,0,0)*val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 2)/(val(cm,0,0,0) + 1e-30);
}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

#line 35
{
  double sum = 0.;
  {foreach_child()
    sum += _const_cm*val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 2)/(_const_cm + 1e-30);
}}

#line 40
}

static inline void face_average (Point point, vector v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
   {




      val(v.x,0,0,0) = (fine(v.x,0,0,0) + fine(v.x,0,1,0))/2.;
      val(v.x,1,0,0) = (fine(v.x,2,0,0) + fine(v.x,2,1,0))/2.;






  } 
#line 44
{




      val(v.y,0,0,0) = (fine(v.y,0,0,0) + fine(v.y,1,0,0))/2.;
      val(v.y,0,1,0) = (fine(v.y,0,2,0) + fine(v.y,1,2,0))/2.;






  }
}

static inline void restriction_face (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  face_average (point, _attribute[s.i].v);
}

static inline void restriction_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  for (int i = 0; i <= 1; i++) {
    val(s,i,0,0) = fine(s,2*i,0,0);

    val(s,i,1,0) = fine(s,2*i,2,0);





  }
}

static inline void no_restriction (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;}

static inline void no_data (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = HUGE;end_foreach_child()}
}

void wavelet (scalar s, scalar w)
{
  restriction (((scalar[]){s,{-1}}));
  for (int l = grid->maxdepth - 1; l >= 0; l--) {
    {foreach_coarse_level (l) {
      {foreach_child()
        val(w,0,0,0) = val(s,0,0,0);end_foreach_child()}
      _attribute[s.i].prolongation (point, s);
      {foreach_child() {
        double sp = val(s,0,0,0);
        val(s,0,0,0) = val(w,0,0,0);

        val(w,0,0,0) -= sp;
      }end_foreach_child()}
    }end_foreach_coarse_level();}
    boundary_level (((scalar[]){w,{-1}}), l + 1);
  }

  {foreach_level(0)
    val(w,0,0,0) = val(s,0,0,0);end_foreach_level();}
  boundary_level (((scalar[]){w,{-1}}), 0);
}

void inverse_wavelet (scalar s, scalar w)
{
  {foreach_level(0)
    val(s,0,0,0) = val(w,0,0,0);end_foreach_level();}
  boundary_level (((scalar[]){s,{-1}}), 0);
  for (int l = 0; l <= grid->maxdepth - 1; l++) {
    {foreach_coarse_level (l) {
      _attribute[s.i].prolongation (point, s);
      {foreach_child()
        val(s,0,0,0) += val(w,0,0,0);end_foreach_child()}
    }end_foreach_coarse_level();}
    boundary_level (((scalar[]){s,{-1}}), l + 1);
  }
}

static inline double bilinear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;



    return (9.*coarse(s,0,0,0) +
     3.*(coarse(s,child.x,0,0) + coarse(s,0,child.y,0)) +
     coarse(s,child.x,child.y,0))/16.;
#line 140 "/home/gongminjiang/basilisk/src/grid/multigrid-common.h"
}

static inline void refine_bilinear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = bilinear (point, s);end_foreach_child()}
}

static inline double quadratic (double a, double b, double c)
{
  return (30.*a + 5.*b - 3.*c)/32.;
}

static inline double biquadratic (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;



  return
    quadratic (quadratic (coarse(s,0,0,0),
     coarse(s,child.x,0,0),
     coarse(s,-child.x,0,0)),
        quadratic (coarse(s,0,child.y,0),
     coarse(s,child.x,child.y,0),
     coarse(s,-child.x,child.y,0)),
        quadratic (coarse(s,0,-child.y,0),
     coarse(s,child.x,-child.y,0),
     coarse(s,-child.x,-child.y,0)));




}

static inline double biquadratic_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;



  return (36.*val(s,0,0,0) + 18.*(val(s,-1,0,0) + val(s,0,-1,0)) - 6.*(val(s,1,0,0) + val(s,0,1,0)) +
   9.*val(s,-1,-1,0) - 3.*(val(s,1,-1,0) + val(s,-1,1,0)) + val(s,1,1,0))/64.;




}

static inline void refine_biquadratic (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = biquadratic (point, s);end_foreach_child()}
}

static inline void refine_linear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

#line 194
if(!is_constant(cm)){{
  coord g;
  if (_attribute[s.i].gradient)
    {
      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
      
#line 198
g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));}
  else
    {
      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
      
#line 201
g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val(cm,0,0,0), sum = val(cm,0,0,0)*(1 << 2);
  {foreach_child() {
    val(s,0,0,0) = sc;
    
      val(s,0,0,0) += child.x*g.x*val(cm,-child.x,0,0)/cmc;
      
#line 207
val(s,0,0,0) += child.y*g.y*val(cm,0,-child.y,0)/cmc;
    sum -= val(cm,0,0,0);
  }end_foreach_child()}
  if (!(fabs(sum) < 1e-10)) qassert ("/home/gongminjiang/basilisk/src/grid/multigrid-common.h", 210, "fabs(sum) < 1e-10");
}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

#line 194
{
  coord g;
  if (_attribute[s.i].gradient)
    {
      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
      
#line 198
g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));}
  else
    {
      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
      
#line 201
g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*_const_cm, sum = _const_cm*(1 << 2);
  {foreach_child() {
    val(s,0,0,0) = sc;
    
      val(s,0,0,0) += child.x*g.x*_const_cm/cmc;
      
#line 207
val(s,0,0,0) += child.y*g.y*_const_cm/cmc;
    sum -= _const_cm;
  }end_foreach_child()}
  if (!(fabs(sum) < 1e-10)) qassert ("/home/gongminjiang/basilisk/src/grid/multigrid-common.h", 210, "fabs(sum) < 1e-10");
}}

#line 211
}

static inline void refine_reset (Point point, scalar v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    val(v,0,0,0) = 0.;end_foreach_child()}
}

static inline void refine_injection (Point point, scalar v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double val = val(v,0,0,0);
  {foreach_child()
    val(v,0,0,0) = val;end_foreach_child()}
}

static scalar multigrid_init_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  _attribute[s.i].prolongation = refine_bilinear;
  _attribute[s.i].restriction = restriction_average;
  return s;
}

static scalar multigrid_init_vertex_scalar (scalar s, const char * name)
{
  s = cartesian_init_vertex_scalar (s, name);
  _attribute[s.i].restriction = restriction_vertex;
  return s;
}

static vector multigrid_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  
    _attribute[v.y.i].restriction = no_restriction;
    
#line 245
_attribute[v.x.i].restriction = no_restriction;
  _attribute[v.x.i].restriction = restriction_face;
  return v;
}

void multigrid_debug (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  cartesian_debug (point);

  FILE * plot = fopen ("plot", "a");
  if (point.level > 0) {
    char name[80] = "coarse";
    if (pid() > 0)
      sprintf (name, "coarse-%d", pid());
    FILE * fp = fopen (name, "w");
#line 271 "/home/gongminjiang/basilisk/src/grid/multigrid-common.h"
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      for (int k = 0; k <= 1; k++)
 for (int l = 0; l <= 1; l++) {
   {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){
     fprintf (fp, "%g %g %g ",
       xc + k*child.x*Delta*2. + _attribute[v.i].d.x*Delta,
       yc + l*child.y*Delta*2. + _attribute[v.i].d.y*Delta,
       coarse(v,k*child.x,l*child.y,0));}}
   fputc ('\n', fp);
 }
      fprintf (ferr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
      fprintf (plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
#line 302 "/home/gongminjiang/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }

  if (is_coarse()) {
    char name[80] = "fine";
    if (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
#line 324 "/home/gongminjiang/basilisk/src/grid/multigrid-common.h"
      double xf = x - Delta/4., yf = y - Delta/4.;
      for (int k = -2; k <= 3; k++)
 for (int l = -2; l <= 3; l++) {
   {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){ {
     fprintf (fp, "%g %g ",
       xf + k*Delta/2. + _attribute[v.i].d.x*Delta/4.,
       yf + l*Delta/2. + _attribute[v.i].d.y*Delta/4.);
     if (allocated_child(k,l,0))
       fprintf (fp, "%g ", fine(v,k,l,0));
     else
       fputs ("n/a ", fp);
   }}}
   fputc ('\n', fp);
 }
      fprintf (ferr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
      fprintf (plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
#line 362 "/home/gongminjiang/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }
  fflush (ferr);
  fclose (plot);
}

static void multigrid_restriction (scalar * list)
{
  scalar * listdef = NULL, * listc = NULL, * list2 = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant (s) && _attribute[s.i].block > 0) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
     list2 = list_add (list2, _attribute[s.i].v.x);
     
#line 381
list2 = list_add (list2, _attribute[s.i].v.y);}
 else
   list2 = list_add (list2, s);
      }
    }}}

  if (listdef || listc) {
    for (int l = depth() - 1; l >= 0; l--) {
      {foreach_coarse_level(l) {
 {scalar*_i=(scalar*)( listdef);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     restriction_average (point, s);}}
 {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
  
     _attribute[s.i].restriction (point, s);
 }}}
      }end_foreach_coarse_level();}
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,__LINE__);
    pfree (listc,__func__,__FILE__,__LINE__);
    pfree (list2,__func__,__FILE__,__LINE__);
  }
}

void multigrid_methods()
{
  cartesian_methods();
  debug = multigrid_debug;
  init_scalar = multigrid_init_scalar;
  init_vertex_scalar = multigrid_init_vertex_scalar;
  init_face_vector = multigrid_init_face_vector;
  restriction = multigrid_restriction;
}







void subtree_size (scalar size, bool leaves)
{




  foreach_stencil()
    {_stencil_val_a(size,0,0,0);  }end_foreach_stencil();




  {
#line 428
foreach()
    val(size,0,0,0) = 1;end_foreach();}





  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar[]){size,{-1}}), depth()); };
  for (int l = depth() - 1; l >= 0; l--) {
    {foreach_coarse_level(l) {
      double sum = !leaves;
      {foreach_child()
 sum += val(size,0,0,0);end_foreach_child()}
      val(size,0,0,0) = sum;
    }end_foreach_coarse_level();}
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar[]){size,{-1}}), l); };
  }
}
#line 5 "/home/gongminjiang/basilisk/src/grid/tree-common.h"




#line 21 "/home/gongminjiang/basilisk/src/grid/tree-common.h"
int refine_cell (Point point, scalar * list, int flag, Cache * refined)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  int nr = 0;


  if (level > 0)
    for (int k = 0; k != 2*child.x; k += child.x)

      for (int l = 0; l != 2*child.y; l += child.y)




   if (aparent(k,l,m).pid >= 0 && is_leaf(aparent(k,l,m))) {
     Point p = point;


     p.level = point.level - 1;
     p.i = (point.i + 2)/2 + k;
     do { if (p.i < 2) p.i += 1 << p.level; else if (p.i >= 2 + (1 << p.level)) p.i -= 1 << p.level; } while(0);

       p.j = (point.j + 2)/2 + l;
       do { if (p.j < 2) p.j += 1 << p.level; else if (p.j >= 2 + (1 << p.level)) p.j -= 1 << p.level; } while(0);





     nr += refine_cell (p, list, flag, refined);
     aparent(k,l,m).flags |= flag;
   }



  increment_neighbors (point);

  int cflag = is_active(cell) ? (active|leaf) : leaf;
  {foreach_child()
    cell.flags |= cflag;end_foreach_child()}


  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (is_local(cell) || _attribute[s.i].face)
      _attribute[s.i].refine (point, s);}}


  cell.flags &= ~leaf;

#if 1
  if (is_border(cell)) {
    {foreach_child() {
      bool bord = false;
      {foreach_neighbor() {
 if (!is_local(cell) || (level > 0 && !is_local(aparent(0,0,0)))) {
   bord = true; foreach_neighbor_break();
 }
 if (is_refined_check())
   {foreach_child()
     if (!is_local(cell)) {
       bord = true; foreach_child_break();
     }end_foreach_child()}
 if (bord)
   foreach_neighbor_break();
      }end_foreach_neighbor()}
      if (bord)
 cell.flags |= border;
    }end_foreach_child()}
    if (refined)
      cache_append (refined, point, cell.flags);
    nr++;
  }
#endif
  return nr;
}





bool coarsen_cell (Point point, scalar * list)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;



  int pid = cell.pid;
  {foreach_child()
    if (cell.neighbors || (cell.pid < 0 && cell.pid != pid))
      return false;end_foreach_child()}



  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    _attribute[s.i].restriction (point, s);
    if (_attribute[s.i].coarsen)
      _attribute[s.i].coarsen (point, s);
  }}}


  cell.flags |= leaf;


  decrement_neighbors (point);

#if 1
  if (!is_local(cell)) {
    cell.flags &= ~(active|border);
    {foreach_neighbor(1)
      if (cell.neighbors)
 {foreach_child()
   if (allocated(0,0,0) && is_local(cell) && !is_border(cell))
     cell.flags |= border;end_foreach_child()}end_foreach_neighbor()}
  }
#endif

  return true;
}

void coarsen_cell_recursive (Point point, scalar * list)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;


  {foreach_child()
    if (cell.neighbors)
      {foreach_neighbor(1)
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
   coarsen_cell_recursive (point, list);end_foreach_neighbor()}end_foreach_child()}

  if (!(coarsen_cell (point, list))) qassert ("/home/gongminjiang/basilisk/src/grid/tree-common.h", 148, "coarsen_cell (point, list)");
}

void mpi_boundary_refine (scalar *);
void mpi_boundary_coarsen (int, int);
void mpi_boundary_update (scalar *);

typedef struct {
  int nc, nf;
} astats;

struct Adapt {
  scalar * slist;
  double * max;
  int maxlevel;
  int minlevel;
  scalar * list;
};

     
astats adapt_wavelet (struct Adapt p)
{tracing("adapt_wavelet","/home/gongminjiang/basilisk/src/grid/tree-common.h",168);
  scalar * list = p.list;

  if (is_constant(cm)) {
    if (list == NULL || list == all)
      list = list_copy (all);
    boundary_internal ((scalar *)list, "/home/gongminjiang/basilisk/src/grid/tree-common.h", 175);
    restriction (p.slist);
  }
  else {
    if (list == NULL || list == all) {
      list = list_copy (((scalar[]){cm, fm.x,fm.y,{-1}}));
      {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 list = list_add (list, s);}}
    }
    boundary_internal ((scalar *)list, "/home/gongminjiang/basilisk/src/grid/tree-common.h", 184);
    scalar * listr = list_concat (p.slist, ((scalar[]){cm,{-1}}));
    restriction (listr);
    pfree (listr,__func__,__FILE__,__LINE__);
  }

  astats st = {0, 0};
  scalar * listc = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].restriction != no_restriction)
      listc = list_add (listc, s);}}


  if (p.minlevel < 1)
    p.minlevel = 1;
  ((Tree *)grid)->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);
  {foreach_cell() {
    if (is_active(cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
 if (cell.flags & too_coarse) {
   cell.flags &= ~too_coarse;
   refine_cell (point, listc, refined, &((Tree *)grid)->refined);
   st.nf++;
 }
 continue;
      }
      else {
 if (cell.flags & refined) {

   cell.flags &= ~too_coarse;
   continue;
 }

 bool local = is_local(cell);
 if (!local)
   {foreach_child()
     if (is_local(cell)) {
       local = true; foreach_child_break();
     }end_foreach_child()}
 if (local) {
   int i = 0;
   static const int just_fine = 1 << (user + 3);
   {scalar*_i=(scalar*)( p.slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
     double max = p.max[i++], sc[1 << 2];
     int c = 0;
     {foreach_child()
       sc[c++] = val(s,0,0,0);end_foreach_child()}
     _attribute[s.i].prolongation (point, s);
     c = 0;
     {foreach_child() {
       double e = fabs(sc[c] - val(s,0,0,0));
       if (e > max && level < p.maxlevel) {
  cell.flags &= ~too_fine;
  cell.flags |= too_coarse;
       }
       else if ((e <= max/1.5 || level > p.maxlevel) &&
         !(cell.flags & (too_coarse|just_fine))) {
  if (level >= p.minlevel)
    cell.flags |= too_fine;
       }
       else if (!(cell.flags & too_coarse)) {
  cell.flags &= ~too_fine;
  cell.flags |= just_fine;
       }
       val(s,0,0,0) = sc[c++];
     }end_foreach_child()}
   }}}
   {foreach_child() {
     cell.flags &= ~just_fine;
     if (!is_leaf(cell)) {
       cell.flags &= ~too_coarse;
       if (level >= p.maxlevel)
  cell.flags |= too_fine;
     }
     else if (!is_active(cell))
       cell.flags &= ~too_coarse;
   }end_foreach_child()}
 }
      }
    }
    else
      continue;
  }end_foreach_cell();}
  mpi_boundary_refine (listc);



  for (int l = depth(); l >= 0; l--) {
    {foreach_cell()
      if (!(cell.pid < 0)) {
 if (level == l) {
   if (!is_leaf(cell)) {
     if (cell.flags & refined)

       cell.flags &= ~(refined|too_fine);
     else if (cell.flags & too_fine) {
       if (is_local(cell) && coarsen_cell (point, listc))
  st.nc++;
       cell.flags &= ~too_fine;
     }
   }
   if (cell.flags & too_fine)
     cell.flags &= ~too_fine;
   else if (level > 0 && (aparent(0,0,0).flags & too_fine))
     aparent(0,0,0).flags &= ~too_fine;
   continue;
 }
 else if (is_leaf(cell))
   continue;
      }end_foreach_cell();}
    mpi_boundary_coarsen (l, too_fine);
  }
  pfree (listc,__func__,__FILE__,__LINE__);

  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (list);

  if (list != p.list)
    pfree (list,__func__,__FILE__,__LINE__);

  {end_tracing("adapt_wavelet","/home/gongminjiang/basilisk/src/grid/tree-common.h",308);return st;}
end_tracing("adapt_wavelet","/home/gongminjiang/basilisk/src/grid/tree-common.h",309);}
#line 331 "/home/gongminjiang/basilisk/src/grid/tree-common.h"
static void refine_level (int depth)
{
  int refined;
  do {
    refined = 0;
    ((Tree *)grid)->refined.n = 0;
    {foreach_leaf()
      if (level < depth) {
 refine_cell (point, NULL, 0, &((Tree *)grid)->refined);
 refined++;
 continue;
      }end_foreach_leaf();}
    mpi_all_reduce (refined, MPI_INT, MPI_SUM);
    if (refined) {
      mpi_boundary_refine (NULL);
      mpi_boundary_update (NULL);
    }
  } while (refined);
}
#line 376 "/home/gongminjiang/basilisk/src/grid/tree-common.h"
     
static void halo_face (vectorl vl)
{tracing("halo_face","/home/gongminjiang/basilisk/src/grid/tree-common.h",377);
  
    {scalar*_i=(scalar*)( vl.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
    
#line 380
{scalar*_i=(scalar*)( vl.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}

  for (int l = depth() - 1; l >= 0; l--)
    {foreach_halo (prolongation, l)
      {
        if (vl.x) {
#line 395 "/home/gongminjiang/basilisk/src/grid/tree-common.h"
   if ((!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))
     {scalar*_i=(scalar*)( vl.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = (fine(s,0,0,0) + fine(s,0,1,0))/2.;}}
   if ((!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0))
     {scalar*_i=(scalar*)( vl.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,1,0,0) = (fine(s,2,0,0) + fine(s,2,1,0))/2.;}}
#line 411 "/home/gongminjiang/basilisk/src/grid/tree-common.h"
 }
        
#line 386
if (vl.y) {
#line 395 "/home/gongminjiang/basilisk/src/grid/tree-common.h"
   if ((!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))
     {scalar*_i=(scalar*)( vl.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = (fine(s,0,0,0) + fine(s,1,0,0))/2.;}}
   if ((!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0))
     {scalar*_i=(scalar*)( vl.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,1,0) = (fine(s,0,2,0) + fine(s,1,2,0))/2.;}}
#line 411 "/home/gongminjiang/basilisk/src/grid/tree-common.h"
 }}end_foreach_halo();}
end_tracing("halo_face","/home/gongminjiang/basilisk/src/grid/tree-common.h",412);}



static scalar tree_init_scalar (scalar s, const char * name)
{
  s = multigrid_init_scalar (s, name);
  _attribute[s.i].refine = _attribute[s.i].prolongation;
  return s;
}

static void prolongation_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  fine(s,1,1,0) = (val(s,0,0,0) + val(s,1,0,0) + val(s,0,1,0) + val(s,1,1,0))/4.;





  for (int i = 0; i <= 1; i++) {
    for (int j = 0; j <= 1; j++)





      if (allocated_child(2*i,2*j,0))
 fine(s,2*i,2*j,0) = val(s,i,j,0);


    
      if (neighbor(i,0,0).neighbors) {

 fine(s,2*i,1,0) = (val(s,i,0,0) + val(s,i,1,0))/2.;
#line 456 "/home/gongminjiang/basilisk/src/grid/tree-common.h"
      }
      
#line 444
if (neighbor(0,i,0).neighbors) {

 fine(s,1,2*i,0) = (val(s,0,i,0) + val(s,1,i,0))/2.;
#line 456 "/home/gongminjiang/basilisk/src/grid/tree-common.h"
      }
  }
}

static scalar tree_init_vertex_scalar (scalar s, const char * name)
{
  s = multigrid_init_vertex_scalar (s, name);
  _attribute[s.i].refine = _attribute[s.i].prolongation = prolongation_vertex;
  return s;
}


static void refine_face_x (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  vector v = _attribute[s.i].v;

  if (!(!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(-1,0,0)))) {
    double g1 = (val(v.x,0,+1,0) - val(v.x,0,-1,0))/8.;
    for (int j = 0; j <= 1; j++)
      fine(v.x,0,j,0) = val(v.x,0,0,0) + (2*j - 1)*g1;
  }
  if (!(!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0) && neighbor(1,0,0).neighbors &&
      (is_local(cell) || is_local(neighbor(1,0,0)))) {
    double g1 = (val(v.x,1,+1,0) - val(v.x,1,-1,0))/8.;
    for (int j = 0; j <= 1; j++)
      fine(v.x,2,j,0) = val(v.x,1,0,0) + (2*j - 1)*g1;
  }
  if (is_local(cell)) {
    double g1 = (val(v.x,0,+1,0) - val(v.x,0,-1,0) + val(v.x,1,+1,0) - val(v.x,1,-1,0))/16.;
    for (int j = 0; j <= 1; j++)
      fine(v.x,1,j,0) = (val(v.x,0,0,0) + val(v.x,1,0,0))/2. + (2*j - 1)*g1;
  }
#line 514 "/home/gongminjiang/basilisk/src/grid/tree-common.h"
}

#line 468
static void refine_face_y (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  vector v = _attribute[s.i].v;

  if (!(!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(0,-1,0)))) {
    double g1 = (val(v.y,+1,0,0) - val(v.y,-1,0,0))/8.;
    for (int j = 0; j <= 1; j++)
      fine(v.y,j,0,0) = val(v.y,0,0,0) + (2*j - 1)*g1;
  }
  if (!(!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0) && neighbor(0,1,0).neighbors &&
      (is_local(cell) || is_local(neighbor(0,1,0)))) {
    double g1 = (val(v.y,+1,1,0) - val(v.y,-1,1,0))/8.;
    for (int j = 0; j <= 1; j++)
      fine(v.y,j,2,0) = val(v.y,0,1,0) + (2*j - 1)*g1;
  }
  if (is_local(cell)) {
    double g1 = (val(v.y,+1,0,0) - val(v.y,-1,0,0) + val(v.y,+1,1,0) - val(v.y,-1,1,0))/16.;
    for (int j = 0; j <= 1; j++)
      fine(v.y,j,1,0) = (val(v.y,0,0,0) + val(v.y,0,1,0))/2. + (2*j - 1)*g1;
  }
#line 514 "/home/gongminjiang/basilisk/src/grid/tree-common.h"
}

void refine_face (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  vector v = _attribute[s.i].v;
  
    _attribute[v.x.i].prolongation (point, v.x);
    
#line 520
_attribute[v.y.i].prolongation (point, v.y);
}

void refine_face_solenoidal (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  refine_face (point, s);

  if (is_local(cell)) {

    vector v = _attribute[s.i].v;
    double d[1 << 2], p[1 << 2];
    int i = 0;
    {foreach_child() {
      d[i] = 0.;
      
 d[i] += val(v.x,1,0,0) - val(v.x,0,0,0);
 
#line 535
d[i] += val(v.y,0,1,0) - val(v.y,0,0,0);
      i++;
    }end_foreach_child()}

    p[0] = 0.;
    p[1] = (3.*d[3] + d[0])/4. + d[2]/2.;
    p[2] = (d[3] + 3.*d[0])/4. + d[2]/2.;
    p[3] = (d[3] + d[0])/2. + d[2];
    fine(v.x,1,1,0) += p[1] - p[0];
    fine(v.x,1,0,0) += p[3] - p[2];
    fine(v.y,0,1,0) += p[0] - p[2];
    fine(v.y,1,1,0) += p[1] - p[3];
#line 574 "/home/gongminjiang/basilisk/src/grid/tree-common.h"
  }

}

vector tree_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  
    _attribute[v.x.i].restriction = _attribute[v.x.i].refine = no_restriction;
    
#line 582
_attribute[v.y.i].restriction = _attribute[v.y.i].refine = no_restriction;
  _attribute[v.x.i].restriction = restriction_face;
  _attribute[v.x.i].refine = refine_face;
  
    _attribute[v.x.i].prolongation = refine_face_x;
    
#line 586
_attribute[v.y.i].prolongation = refine_face_y;
  return v;
}

     
static void tree_boundary_level (scalar * list, int l)
{tracing("tree_boundary_level","/home/gongminjiang/basilisk/src/grid/tree-common.h",591);
  int depth = l < 0 ? depth() : l;

  if (tree_is_full()) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, depth); };
    {end_tracing("tree_boundary_level","/home/gongminjiang/basilisk/src/grid/tree-common.h",597);return;}
  }

  scalar * listdef = NULL, * listc = NULL, * list2 = NULL, * vlist = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant (s)) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
     list2 = list_add (list2, _attribute[s.i].v.x);
     
#line 611
list2 = list_add (list2, _attribute[s.i].v.y);}
 else {
   list2 = list_add (list2, s);
   if (_attribute[s.i].restriction == restriction_vertex)
     vlist = list_add (vlist, s);
 }
      }
    }}}

  if (vlist)






    {foreach_vertex () {
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || (!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) ||
   (!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) || (!is_leaf (neighbor(-1,-1,0)) && neighbor(-1,-1,0).neighbors && neighbor(-1,-1,0).pid >= 0)) {

 {scalar*_i=(scalar*)( vlist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   val(s,0,0,0) = is_vertex (child(0,0,0)) ? fine(s,0,0,0) : HUGE;}}
      }
      else
 {
   if (child.y == 1 &&
       ((!is_leaf(cell) && !cell.neighbors && cell.pid >= 0) || (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))) {

     {scalar*_i=(scalar*)( vlist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = is_vertex(neighbor(0,-1,0)) && is_vertex(neighbor(0,1,0)) ?
  (val(s,0,-1,0) + val(s,0,1,0))/2. : HUGE;}}
   }
   
#line 636
if (child.x == 1 &&
       ((!is_leaf(cell) && !cell.neighbors && cell.pid >= 0) || (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))) {

     {scalar*_i=(scalar*)( vlist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
       val(s,0,0,0) = is_vertex(neighbor(-1,0,0)) && is_vertex(neighbor(1,0,0)) ?
  (val(s,-1,0,0) + val(s,1,0,0))/2. : HUGE;}}
   }}
    }end_foreach_vertex();}
#line 676 "/home/gongminjiang/basilisk/src/grid/tree-common.h"
  pfree (vlist,__func__,__FILE__,__LINE__);

  if (listdef || listc) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, list2, depth); };
    for (int l = depth - 1; l >= 0; l--) {
      {foreach_coarse_level(l) {
 {scalar*_i=(scalar*)( listdef);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   restriction_average (point, s);}}
 {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   _attribute[s.i].restriction (point, s);}}
      }end_foreach_coarse_level();}
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,__LINE__);
    pfree (listc,__func__,__FILE__,__LINE__);
    pfree (list2,__func__,__FILE__,__LINE__);
  }

  scalar * listr = NULL;
  vector * listf = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant (s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].face)
 listf = vectors_add (listf, _attribute[s.i].v);
      else
 listr = list_add (listr, s);
    }}}

  if (listr || listf) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, 0); };
    for (int i = 0; i < depth; i++) {
      {foreach_halo (prolongation, i) {
 {scalar*_i=(scalar*)( listr);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
          _attribute[s.i].prolongation (point, s);}}
 {vector*_i=(vector*)( listf);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
   {
     _attribute[v.x.i].prolongation (point, v.x);
     
#line 712
_attribute[v.y.i].prolongation (point, v.y);}}}
      }end_foreach_halo();}
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, i + 1); };
    }
    pfree (listr,__func__,__FILE__,__LINE__);
    pfree (listf,__func__,__FILE__,__LINE__);
  }
end_tracing("tree_boundary_level","/home/gongminjiang/basilisk/src/grid/tree-common.h",719);}

double treex (Point point) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (level == 0)
    return 0;

  double i = 2*child.x - child.y;
  if (i <= 1 && i >= -1) i = -i;




  return treex(parent) + i/(1 << 2*(level - 1));
}

double treey (Point point) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (level == 0)
    return 0;
  return treey(parent) + 4./(1 << 2*(level - 1));
}

void output_tree (FILE * fp)
{
  {foreach_cell()
    if (cell.neighbors)
      {foreach_child()
 if (is_local(cell))
   fprintf (fp, "%g %g\n%g %g\n\n",
     treex(parent), treey(parent), treex(point), treey(point));end_foreach_child()}end_foreach_cell();}
}

     
void tree_check()
{tracing("tree_check","/home/gongminjiang/basilisk/src/grid/tree-common.h",751);


  long nleaves = 0, nactive = 0;
  {foreach_cell_all() {
    if (is_leaf(cell)) {
      if (!(cell.pid >= 0)) qassert ("/home/gongminjiang/basilisk/src/grid/tree-common.h", 758, "cell.pid >= 0");
      nleaves++;
    }
    if (is_local(cell))
      if (!(is_active(cell) || (!is_leaf(cell) && !cell.neighbors && cell.pid >= 0))) qassert ("/home/gongminjiang/basilisk/src/grid/tree-common.h", 762, "is_active(cell) || is_prolongation(cell)");
    if (is_active(cell))
      nactive++;

    int neighbors = 0;
    {foreach_neighbor(1)
      if (allocated(0,0,0) && (!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 neighbors++;end_foreach_neighbor()}
    if (!(cell.neighbors == neighbors)) qassert ("/home/gongminjiang/basilisk/src/grid/tree-common.h", 770, "cell.neighbors == neighbors");


    if (!cell.neighbors)
      if (!(!allocated_child(0,0,0))) qassert ("/home/gongminjiang/basilisk/src/grid/tree-common.h", 774, "!allocated_child(0)");
  }end_foreach_cell_all();}


  long reachable = 0;
  {foreach_cell() {
    if (is_active(cell))
      reachable++;
    else
      continue;
  }end_foreach_cell();}
  if (!(nactive == reachable)) qassert ("/home/gongminjiang/basilisk/src/grid/tree-common.h", 785, "nactive == reachable");


  reachable = 0;
  {foreach_cell()
    if (is_leaf(cell)) {
      reachable++;
      continue;
    }end_foreach_cell();}
  if (!(nleaves == reachable)) qassert ("/home/gongminjiang/basilisk/src/grid/tree-common.h", 794, "nleaves == reachable");
end_tracing("tree_check","/home/gongminjiang/basilisk/src/grid/tree-common.h",795);}

     
static void tree_restriction (scalar * list) {tracing("tree_restriction","/home/gongminjiang/basilisk/src/grid/tree-common.h",798);
  boundary_internal ((scalar *)list, "/home/gongminjiang/basilisk/src/grid/tree-common.h", 799);
  if (tree_is_full())
    multigrid_restriction (list);
end_tracing("tree_restriction","/home/gongminjiang/basilisk/src/grid/tree-common.h",802);}

void tree_methods()
{
  multigrid_methods();
  init_scalar = tree_init_scalar;
  init_vertex_scalar = tree_init_vertex_scalar;
  init_face_vector = tree_init_face_vector;
  boundary_level = tree_boundary_level;
  boundary_face = halo_face;
  restriction = tree_restriction;
}
#line 1672 "/home/gongminjiang/basilisk/src/grid/tree.h"


void tree_periodic (int dir)
{
  int depth = grid ? depth() : -1;
  if (grid)
    free_grid();
  periodic (dir);
  if (depth >= 0)
    init_grid (1 << depth);
}


#if 1
#line 1 "grid/tree-mpi.h"
#line 1 "/home/gongminjiang/basilisk/src/grid/tree-mpi.h"

int debug_iteration = -1;

void debug_mpi (FILE * fp1);

typedef struct {
  CacheLevel * halo;
  void * buf;
  MPI_Request r;
  int depth;
  int pid;
  int maxdepth;
} Rcv;

typedef struct {
  Rcv * rcv;
  char * name;
  int npid;
} RcvPid;

typedef struct {
  RcvPid * rcv, * snd;
} SndRcv;

typedef struct {
  Boundary parent;

  SndRcv mpi_level, mpi_level_root, restriction;
  Array * send, * receive;
} MpiBoundary;

static void cache_level_init (CacheLevel * c)
{
  c->p = NULL;
  c->n = c->nm = 0;
}

static void rcv_append (Point point, Rcv * rcv)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (level > rcv->depth) {
    rcv->halo = (CacheLevel *) prealloc (rcv->halo, (level + 1)*sizeof(CacheLevel),__func__,__FILE__,__LINE__);
    for (int j = rcv->depth + 1; j <= level; j++)
      cache_level_init (&rcv->halo[j]);
    rcv->depth = level;
  }
  cache_level_append (&rcv->halo[level], point);
  if (level > rcv->maxdepth)
    rcv->maxdepth = level;
}

void rcv_print (Rcv * rcv, FILE * fp, const char * prefix)
{
  for (int l = 0; l <= rcv->depth; l++)
    if (rcv->halo[l].n > 0)
      {foreach_cache_level(rcv->halo[l], l)
 fprintf (fp, "%s%g %g %g %d %d\n", prefix, x, y, z, rcv->pid, level);end_foreach_cache_level();}
}

static void rcv_free_buf (Rcv * rcv)
{
  if (rcv->buf) {
    prof_start ("rcv_pid_receive");
    MPI_Wait (&rcv->r, MPI_STATUS_IGNORE);
    pfree (rcv->buf,__func__,__FILE__,__LINE__);
    rcv->buf = NULL;
    prof_stop();
  }
}

static void rcv_destroy (Rcv * rcv)
{
  rcv_free_buf (rcv);
  for (int i = 0; i <= rcv->depth; i++)
    if (rcv->halo[i].n > 0)
      pfree (rcv->halo[i].p,__func__,__FILE__,__LINE__);
  pfree (rcv->halo,__func__,__FILE__,__LINE__);
}

static RcvPid * rcv_pid_new (const char * name)
{
  RcvPid * r = ((RcvPid *) pcalloc (1, sizeof(RcvPid),__func__,__FILE__,__LINE__));
  r->name = pstrdup (name,__func__,__FILE__,__LINE__);
  return r;
}

static Rcv * rcv_pid_pointer (RcvPid * p, int pid)
{
  if (!(pid >= 0 && pid < npe())) qassert ("/home/gongminjiang/basilisk/src/grid/tree-mpi.h", 88, "pid >= 0 && pid < npe()");

  int i;
  for (i = 0; i < p->npid; i++)
    if (pid == p->rcv[i].pid)
      break;

  if (i == p->npid) {
    p->rcv = (Rcv *) prealloc (p->rcv, (++p->npid)*sizeof(Rcv),__func__,__FILE__,__LINE__);
    Rcv * rcv = &p->rcv[p->npid-1];
    rcv->pid = pid;
    rcv->depth = rcv->maxdepth = 0;
    rcv->halo = ((CacheLevel *) pmalloc ((1)*sizeof(CacheLevel),__func__,__FILE__,__LINE__));
    rcv->buf = NULL;
    cache_level_init (&rcv->halo[0]);
  }
  return &p->rcv[i];
}

static void rcv_pid_append (RcvPid * p, int pid, Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  rcv_append (point, rcv_pid_pointer (p, pid));
}

static void rcv_pid_append_pids (RcvPid * p, Array * pids)
{

  for (int i = 0; i < p->npid; i++) {
    int pid = p->rcv[i].pid, j, * a;
    for (j = 0, a = pids->p; j < pids->len/sizeof(int); j++,a++)
      if (*a == pid)
 break;
    if (j == pids->len/sizeof(int))
      array_append (pids, &pid, sizeof(int));
  }
}

void rcv_pid_write (RcvPid * p, const char * name)
{
  for (int i = 0; i < p->npid; i++) {
    Rcv * rcv = &p->rcv[i];
    char fname[80];
    sprintf (fname, "%s-%d-%d", name, pid(), rcv->pid);
    FILE * fp = fopen (fname, "w");
    rcv_print (rcv, fp, "");
    fclose (fp);
  }
}

static void rcv_pid_print (RcvPid * p, FILE * fp, const char * prefix)
{
  for (int i = 0; i < p->npid; i++)
    rcv_print (&p->rcv[i], fp, prefix);
}

static void rcv_pid_destroy (RcvPid * p)
{
  for (int i = 0; i < p->npid; i++)
    rcv_destroy (&p->rcv[i]);
  pfree (p->rcv,__func__,__FILE__,__LINE__);
  pfree (p->name,__func__,__FILE__,__LINE__);
  pfree (p,__func__,__FILE__,__LINE__);
}

static Boundary * mpi_boundary = NULL;






void debug_mpi (FILE * fp1);

static void apply_bc (Rcv * rcv, scalar * list, scalar * listv,
        vector * listf, int l, MPI_Status s)
{
  double * b = rcv->buf;
  {foreach_cache_level(rcv->halo[l], l) {
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      memcpy (&val(s,0,0,0), b, sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    }}}
    {vector*_i=(vector*)( listf);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
      { {
 memcpy (&val(v.x,0,0,0), b, sizeof(double)*_attribute[v.x.i].block);
 b += _attribute[v.x.i].block;
 if (*b != HUGE && allocated(1,0,0))
   memcpy (&val(v.x,1,0,0), b, sizeof(double)*_attribute[v.x.i].block);
 b += _attribute[v.x.i].block;
      } 
#line 171
{
 memcpy (&val(v.y,0,0,0), b, sizeof(double)*_attribute[v.y.i].block);
 b += _attribute[v.y.i].block;
 if (*b != HUGE && allocated(0,1,0))
   memcpy (&val(v.y,0,1,0), b, sizeof(double)*_attribute[v.y.i].block);
 b += _attribute[v.y.i].block;
      }}}}
    {scalar*_i=(scalar*)( listv);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)







          {
     if (*b != HUGE && allocated(i,j,0))
       memcpy (&val(s,i,j,0), b, sizeof(double)*_attribute[s.i].block);
     b += _attribute[s.i].block;
          }

    }}}
  }end_foreach_cache_level();}
  size_t size = b - (double *) rcv->buf;
  pfree (rcv->buf,__func__,__FILE__,__LINE__);
  rcv->buf = NULL;

  int rlen;
  MPI_Get_count (&s, MPI_DOUBLE, &rlen);
  if (rlen != size) {
    fprintf (ferr,
      "rlen (%d) != size (%ld), %d receiving from %d at level %d\n"
      "Calling debug_mpi(NULL)...\n"
      "Aborting...\n",
      rlen, size, pid(), rcv->pid, l);
    fflush (ferr);
    debug_mpi (NULL);
    MPI_Abort (MPI_COMM_WORLD, -2);
  }
}
#line 234 "/home/gongminjiang/basilisk/src/grid/tree-mpi.h"
static void mpi_recv_check (void * buf, int count, MPI_Datatype datatype,
       int source, int tag,
       MPI_Comm comm, MPI_Status * status,
       const char * name)
{
#line 269 "/home/gongminjiang/basilisk/src/grid/tree-mpi.h"
  int errorcode = MPI_Recv (buf, count, datatype, source, tag, comm, status);
  if (errorcode != MPI_SUCCESS) {
    char string[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_string (errorcode, string, &resultlen);
    fprintf (ferr,
      "ERROR MPI_Recv \"%s\" (count = %d, source = %d, tag = %d):\n%s\n"
      "Calling debug_mpi(NULL)...\n"
      "Aborting...\n",
      name, count, source, tag, string);
    fflush (ferr);
    debug_mpi (NULL);
    MPI_Abort (MPI_COMM_WORLD, -1);
  }





}

     
static int mpi_waitany (int count, MPI_Request array_of_requests[], int *indx,
   MPI_Status *status)
{tracing("mpi_waitany","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",291);
  { int _ret= MPI_Waitany (count, array_of_requests, indx, status);end_tracing("mpi_waitany","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",294);return _ret;}
end_tracing("mpi_waitany","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",295);}

static int list_lenb (scalar * list) {
  int len = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    len += _attribute[s.i].block;}}
  return len;
}

static int vectors_lenb (vector * list) {
  int len = 0;
  {vector*_i=(vector*)( list);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
    len += _attribute[v.x.i].block;}}
  return len;
}

static void rcv_pid_receive (RcvPid * m, scalar * list, scalar * listv,
        vector * listf, int l)
{
  if (m->npid == 0)
    return;

  prof_start ("rcv_pid_receive");

  int len = list_lenb (list) + 2*2*vectors_lenb (listf) +
    (1 << 2)*list_lenb (listv);

  MPI_Request r[m->npid];
  Rcv * rrcv[m->npid];
  int nr = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      if (!(!rcv->buf)) qassert ("/home/gongminjiang/basilisk/src/grid/tree-mpi.h", 328, "!rcv->buf");
      rcv->buf = pmalloc (sizeof (double)*rcv->halo[l].n*len,__func__,__FILE__,__LINE__);






      MPI_Irecv (rcv->buf, rcv->halo[l].n*len, MPI_DOUBLE, rcv->pid,
   (l), MPI_COMM_WORLD, &r[nr]);
      rrcv[nr++] = rcv;






    }
  }


  if (nr > 0) {
    int i;
    MPI_Status s;
    mpi_waitany (nr, r, &i, &s);
    while (i != MPI_UNDEFINED) {
      Rcv * rcv = rrcv[i];
      if (!(l <= rcv->depth && rcv->halo[l].n > 0)) qassert ("/home/gongminjiang/basilisk/src/grid/tree-mpi.h", 355, "l <= rcv->depth && rcv->halo[l].n > 0");
      if (!(rcv->buf)) qassert ("/home/gongminjiang/basilisk/src/grid/tree-mpi.h", 356, "rcv->buf");
      apply_bc (rcv, list, listv, listf, l, s);
      mpi_waitany (nr, r, &i, &s);
    }
  }

  prof_stop();
}

     
static void rcv_pid_wait (RcvPid * m)
{tracing("rcv_pid_wait","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",366);

  for (int i = 0; i < m->npid; i++)
    rcv_free_buf (&m->rcv[i]);
end_tracing("rcv_pid_wait","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",371);}

static void rcv_pid_send (RcvPid * m, scalar * list, scalar * listv,
     vector * listf, int l)
{
  if (m->npid == 0)
    return;

  prof_start ("rcv_pid_send");

  int len = list_lenb (list) + 2*2*vectors_lenb (listf) +
    (1 << 2)*list_lenb (listv);


  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      if (!(!rcv->buf)) qassert ("/home/gongminjiang/basilisk/src/grid/tree-mpi.h", 388, "!rcv->buf");
      rcv->buf = pmalloc (sizeof (double)*rcv->halo[l].n*len,__func__,__FILE__,__LINE__);
      double * b = rcv->buf;
      {foreach_cache_level(rcv->halo[l], l) {
 {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
   memcpy (b, &val(s,0,0,0), sizeof(double)*_attribute[s.i].block);
   b += _attribute[s.i].block;
 }}}
 {vector*_i=(vector*)( listf);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
   { {
     memcpy (b, &val(v.x,0,0,0), sizeof(double)*_attribute[v.x.i].block);
     b += _attribute[v.x.i].block;
     if (allocated(1,0,0))
       memcpy (b, &val(v.x,1,0,0), sizeof(double)*_attribute[v.x.i].block);
     else
       *b = HUGE;
     b += _attribute[v.x.i].block;
   } 
#line 397
{
     memcpy (b, &val(v.y,0,0,0), sizeof(double)*_attribute[v.y.i].block);
     b += _attribute[v.y.i].block;
     if (allocated(0,1,0))
       memcpy (b, &val(v.y,0,1,0), sizeof(double)*_attribute[v.y.i].block);
     else
       *b = HUGE;
     b += _attribute[v.y.i].block;
   }}}}
 {scalar*_i=(scalar*)( listv);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
   for (int i = 0; i <= 1; i++)
     for (int j = 0; j <= 1; j++)
#line 418 "/home/gongminjiang/basilisk/src/grid/tree-mpi.h"
       {
  if (allocated(i,j,0))
    memcpy (b, &val(s,i,j,0), sizeof(double)*_attribute[s.i].block);
  else
    *b = HUGE;
  b += _attribute[s.i].block;
       }

 }}}
      }end_foreach_cache_level();}





      MPI_Isend (rcv->buf, (b - (double *) rcv->buf),
   MPI_DOUBLE, rcv->pid, (l), MPI_COMM_WORLD,
   &rcv->r);
    }
  }

  prof_stop();
}

static void rcv_pid_sync (SndRcv * m, scalar * list, int l)
{
  scalar * listr = NULL, * listv = NULL;
  vector * listf = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].block > 0) {
      if (_attribute[s.i].face)
 listf = vectors_add (listf, _attribute[s.i].v);
      else if (_attribute[s.i].restriction == restriction_vertex)
 listv = list_add (listv, s);
      else
 listr = list_add (listr, s);
    }}}
  rcv_pid_send (m->snd, listr, listv, listf, l);
  rcv_pid_receive (m->rcv, listr, listv, listf, l);
  rcv_pid_wait (m->snd);
  pfree (listr,__func__,__FILE__,__LINE__);
  pfree (listf,__func__,__FILE__,__LINE__);
  pfree (listv,__func__,__FILE__,__LINE__);
}

static void snd_rcv_destroy (SndRcv * m)
{
  rcv_pid_destroy (m->rcv);
  rcv_pid_destroy (m->snd);
}

static void snd_rcv_init (SndRcv * m, const char * name)
{
  char s[strlen(name) + 5];
  strcpy (s, name);
  strcat (s, ".rcv");
  m->rcv = rcv_pid_new (s);
  strcpy (s, name);
  strcat (s, ".snd");
  m->snd = rcv_pid_new (s);
}

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
  snd_rcv_destroy (&m->mpi_level);
  snd_rcv_destroy (&m->mpi_level_root);
  snd_rcv_destroy (&m->restriction);
  array_free (m->send);
  array_free (m->receive);
  pfree (m,__func__,__FILE__,__LINE__);
}

     
static void mpi_boundary_level (const Boundary * b, scalar * list, int l)
{tracing("mpi_boundary_level","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",492);
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->mpi_level, list, l);
  rcv_pid_sync (&m->mpi_level_root, list, l);
end_tracing("mpi_boundary_level","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",497);}

     
static void mpi_boundary_restriction (const Boundary * b, scalar * list, int l)
{tracing("mpi_boundary_restriction","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",500);
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->restriction, list, l);
end_tracing("mpi_boundary_restriction","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",504);}

void mpi_boundary_new()
{
  mpi_boundary = (Boundary *) ((MpiBoundary *) pcalloc (1, sizeof(MpiBoundary),__func__,__FILE__,__LINE__));
  mpi_boundary->destroy = mpi_boundary_destroy;
  mpi_boundary->level = mpi_boundary_level;
  mpi_boundary->restriction = mpi_boundary_restriction;
  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;
  snd_rcv_init (&mpi->mpi_level, "mpi_level");
  snd_rcv_init (&mpi->mpi_level_root, "mpi_level_root");
  snd_rcv_init (&mpi->restriction, "restriction");
  mpi->send = array_new();
  mpi->receive = array_new();
  add_boundary (mpi_boundary);
}

static FILE * fopen_prefix (FILE * fp, const char * name, char * prefix)
{
  if (fp) {
    sprintf (prefix, "%s-%d ", name, pid());
    return fp;
  }
  else {
    strcpy (prefix, "");
    char fname[80];
    if (debug_iteration >= 0)
      sprintf (fname, "%s-%d-%d", name, debug_iteration, pid());
    else
      sprintf (fname, "%s-%d", name, pid());
    return fopen (fname, "w");
  }
}

void debug_mpi (FILE * fp1)
{
  void output_cells_internal (FILE * fp);

  char prefix[80];
  FILE * fp;


  if (fp1 == NULL) {
    char name[80];
    sprintf (name, "halo-%d", pid()); remove (name);
    sprintf (name, "cells-%d", pid()); remove (name);
    sprintf (name, "faces-%d", pid()); remove (name);
    sprintf (name, "vertices-%d", pid()); remove (name);
    sprintf (name, "neighbors-%d", pid()); remove (name);
    sprintf (name, "mpi-level-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-level-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-level-root-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-level-root-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-restriction-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-restriction-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-border-%d", pid()); remove (name);
    sprintf (name, "exterior-%d", pid()); remove (name);
    sprintf (name, "depth-%d", pid()); remove (name);
    sprintf (name, "refined-%d", pid()); remove (name);
  }


  fp = fopen_prefix (fp1, "halo", prefix);
  for (int l = 0; l < depth(); l++)
    {foreach_halo (prolongation, l)
      {foreach_child()
        fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);end_foreach_child()}end_foreach_halo();}
  if (!fp1)
    fclose (fp);

  if (!fp1) {
    fp = fopen_prefix (fp1, "cells", prefix);
    output_cells_internal (fp);
    fclose (fp);
  }

  fp = fopen_prefix (fp1, "faces", prefix);
  {foreach_face_generic(){is_face_x(){
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);}end_is_face_x()
#line 581
is_face_y(){
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);}end_is_face_y()}end_foreach_face_generic();}
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "vertices", prefix);
  {foreach_vertex()
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);end_foreach_vertex();}
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "neighbors", prefix);
  {foreach() {
    int n = 0;
    {foreach_neighbor(1)
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 n++;end_foreach_neighbor()}
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, cell.neighbors);
    if (!(cell.neighbors == n)) qassert ("/home/gongminjiang/basilisk/src/grid/tree-mpi.h", 599, "cell.neighbors == n");
  }end_foreach();}
  if (!fp1)
    fclose (fp);

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;

  fp = fopen_prefix (fp1, "mpi-level-rcv", prefix);
  rcv_pid_print (mpi->mpi_level.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-root-rcv", prefix);
  rcv_pid_print (mpi->mpi_level_root.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-restriction-rcv", prefix);
  rcv_pid_print (mpi->restriction.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-snd", prefix);
  rcv_pid_print (mpi->mpi_level.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-root-snd", prefix);
  rcv_pid_print (mpi->mpi_level_root.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-restriction-snd", prefix);
  rcv_pid_print (mpi->restriction.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-border", prefix);
  {foreach_cell() {
    if (is_border(cell))
      fprintf (fp, "%s%g %g %g %d %d %d\n",
        prefix, x, y, z, level, cell.neighbors, cell.pid);
    else
      continue;
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "exterior", prefix);
  {foreach_cell() {
    if (!is_local(cell))
      fprintf (fp, "%s%g %g %g %d %d %d %d\n",
        prefix, x, y, z, level, cell.neighbors,
        cell.pid, cell.flags & leaf);






  }end_foreach_cell();}
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "depth", prefix);
  fprintf (fp, "depth: %d %d\n", pid(), depth());
  fprintf (fp, "======= mpi_level.snd ======\n");
  RcvPid * snd = mpi->mpi_level.snd;
  for (int i = 0; i < snd->npid; i++)
    fprintf (fp, "%d %d %d\n", pid(), snd->rcv[i].pid, snd->rcv[i].maxdepth);
  fprintf (fp, "======= mpi_level.rcv ======\n");
  snd = mpi->mpi_level.rcv;
  for (int i = 0; i < snd->npid; i++)
    fprintf (fp, "%d %d %d\n", pid(), snd->rcv[i].pid, snd->rcv[i].maxdepth);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "refined", prefix);
  {foreach_cache (((Tree *)grid)->refined)
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level);end_foreach_cache();}
  if (!fp1)
    fclose (fp);
}

static void snd_rcv_free (SndRcv * p)
{
  char name[strlen(p->rcv->name) + 1];
  strcpy (name, p->rcv->name);
  rcv_pid_destroy (p->rcv);
  p->rcv = rcv_pid_new (name);
  strcpy (name, p->snd->name);
  rcv_pid_destroy (p->snd);
  p->snd = rcv_pid_new (name);
}

static bool is_root (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
    {foreach_child()
      if (is_local(cell))
 return true;end_foreach_child()}
  return false;
}


static bool is_local_prolongation (Point point, Point p)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

  struct { int x, y; } dp = {p.i - point.i, p.j - point.j};



   {
    if (dp.x == 0 && ((!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) || (!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(dp.x,0,0)) && neighbor(dp.x,0,0).neighbors && neighbor(dp.x,0,0).pid >= 0))
      return true;
  } 
#line 713
{
    if (dp.y == 0 && ((!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) || (!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(0,dp.y,0)) && neighbor(0,dp.y,0).neighbors && neighbor(0,dp.y,0).pid >= 0))
      return true;
  }
  return false;
}



static void append_pid (Array * pids, int pid)
{
  for (int i = 0, * p = (int *) pids->p; i < pids->len/sizeof(int); i++, p++)
    if (*p == pid)
      return;
  array_append (pids, &pid, sizeof(int));
}

static int locals_pids (Point point, Array * pids)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (is_leaf(cell)) {
    if (is_local(cell)) {
      Point p = point;
      {foreach_neighbor(1) {
 if ((cell.pid >= 0 && cell.pid != pid()) &&
     ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || is_local_prolongation (point, p)))
   append_pid (pids, cell.pid);
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
   {foreach_child()
     if ((cell.pid >= 0 && cell.pid != pid()))
       append_pid (pids, cell.pid);end_foreach_child()}
      }end_foreach_neighbor()}
    }
  }
  else
    {foreach_neighbor(1) {
      if ((cell.pid >= 0 && cell.pid != pid()))
 append_pid (pids, cell.pid);
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 {foreach_child()
   if ((cell.pid >= 0 && cell.pid != pid()))
     append_pid (pids, cell.pid);end_foreach_child()}
    }end_foreach_neighbor()}
  return pids->len/sizeof(int);
}

static int root_pids (Point point, Array * pids)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    if ((cell.pid >= 0 && cell.pid != pid()))
      append_pid (pids, cell.pid);end_foreach_child()}
  return pids->len/sizeof(int);
}







static void rcv_pid_row (RcvPid * m, int l, int * row)
{
  for (int i = 0; i < npe(); i++)
    row[i] = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0)
      row[rcv->pid] = rcv->halo[l].n;
  }
}

void check_snd_rcv_matrix (SndRcv * sndrcv, const char * name)
{
  int maxlevel = depth();
  mpi_all_reduce (maxlevel, MPI_INT, MPI_MAX);
  int * row = ((int *) pmalloc ((npe())*sizeof(int),__func__,__FILE__,__LINE__));
  for (int l = 0; l <= maxlevel; l++) {
    int status = 0;
    if (pid() == 0) {


      int ** send = matrix_new (npe(), npe(), sizeof(int));
      int ** receive = matrix_new (npe(), npe(), sizeof(int));
      rcv_pid_row (sndrcv->snd, l, row);
      MPI_Gather (row, npe(), MPI_INT, &send[0][0], npe(), MPI_INT, 0,
    MPI_COMM_WORLD);
      rcv_pid_row (sndrcv->rcv, l, row);
      MPI_Gather (row, npe(), MPI_INT, &receive[0][0], npe(), MPI_INT, 0,
    MPI_COMM_WORLD);

      int * astatus = ((int *) pmalloc ((npe())*sizeof(int),__func__,__FILE__,__LINE__));
      for (int i = 0; i < npe(); i++)
 astatus[i] = 0;
      for (int i = 0; i < npe(); i++)
 for (int j = 0; j < npe(); j++)
   if (send[i][j] != receive[j][i]) {
     fprintf (ferr, "%s: %d sends    %d to   %d at level %d\n",
       name, i, send[i][j], j, l);
     fprintf (ferr, "%s: %d receives %d from %d at level %d\n",
       name, j, receive[j][i], i, l);
     fflush (ferr);
     for (int k = i - 2; k <= i + 2; k++)
       if (k >= 0 && k < npe())
  astatus[k] = 1;
     for (int k = j - 2; k <= j + 2; k++)
       if (k >= 0 && k < npe())
  astatus[k] = 1;
   }
      MPI_Scatter (astatus, 1, MPI_INT, &status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      pfree (astatus,__func__,__FILE__,__LINE__);

      matrix_free (send);
      matrix_free (receive);
    }
    else {
      rcv_pid_row (sndrcv->snd, l, row);
      MPI_Gather (row, npe(), MPI_INT, NULL, npe(), MPI_INT, 0, MPI_COMM_WORLD);
      rcv_pid_row (sndrcv->rcv, l, row);
      MPI_Gather (row, npe(), MPI_INT, NULL, npe(), MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Scatter (NULL, 1, MPI_INT, &status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    if (status) {
      fprintf (ferr,
        "check_snd_rcv_matrix \"%s\" failed\n"
        "Calling debug_mpi(NULL)...\n"
        "Aborting...\n",
        name);
      fflush (ferr);
      debug_mpi (NULL);
      MPI_Abort (MPI_COMM_WORLD, -3);
    }
  }
  pfree (row,__func__,__FILE__,__LINE__);
}

static bool has_local_child (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child()
    if (is_local(cell))
      return true;end_foreach_child()}
  return false;
}

     
void mpi_boundary_update_buffers()
{tracing("mpi_boundary_update_buffers","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",858);
  if (npe() == 1)
    {end_tracing("mpi_boundary_update_buffers","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",861);return;}

  prof_start ("mpi_boundary_update_buffers");

  MpiBoundary * m = (MpiBoundary *) mpi_boundary;
  SndRcv * mpi_level = &m->mpi_level;
  SndRcv * mpi_level_root = &m->mpi_level_root;
  SndRcv * restriction = &m->restriction;

  snd_rcv_free (mpi_level);
  snd_rcv_free (mpi_level_root);
  snd_rcv_free (restriction);

  static const unsigned short used = 1 << user;
  {foreach_cell() {
    if (is_active(cell) && !is_border(cell))



      continue;

    if (cell.neighbors) {

      Array pids = {NULL, 0, 0};
      int n = locals_pids (point, &pids);
      if (n) {
 {foreach_child()
   if (is_local(cell))
     for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
       rcv_pid_append (mpi_level->snd, *p, point);end_foreach_child()}
 pfree (pids.p,__func__,__FILE__,__LINE__);
      }

      bool locals = false;
      if (is_leaf(cell)) {
 if ((cell.pid >= 0 && cell.pid != pid())) {
   Point p = point;
   {foreach_neighbor(1)
     if ((is_local(cell) &&
   ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || is_local_prolongation (point, p))) ||
  is_root(point)) {
       locals = true; foreach_neighbor_break();
     }end_foreach_neighbor()}
 }
      }
      else
 {foreach_neighbor(1)
   if (is_local(cell) || is_root(point)) {
     locals = true; foreach_neighbor_break();
   }end_foreach_neighbor()}
      if (locals)
 {foreach_child()
   if ((cell.pid >= 0 && cell.pid != pid()))
            rcv_pid_append (mpi_level->rcv, cell.pid, point),
       cell.flags |= used;end_foreach_child()}


      if (!is_leaf(cell)) {

 if (is_local(cell)) {
   Array pids = {NULL, 0, 0};

   int n = root_pids (point, &pids);
   if (n) {
     {foreach_neighbor()
       for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
  if (cell.pid >= 0 && cell.pid != *p)
    rcv_pid_append (mpi_level_root->snd, *p, point);end_foreach_neighbor()}

     for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
       rcv_pid_append (restriction->snd, *p, point);
     pfree (pids.p,__func__,__FILE__,__LINE__);
   }
 }

 else if ((cell.pid >= 0 && cell.pid != pid())) {
   bool root = false;
   {foreach_child()
     if (is_local(cell)) {
       root = true; foreach_child_break();
     }end_foreach_child()}
   if (root) {
     int pid = cell.pid;
     {foreach_neighbor()
       if ((cell.pid >= 0 && cell.pid != pid()))
  rcv_pid_append (mpi_level_root->rcv, pid, point),
    cell.flags |= used;end_foreach_neighbor()}

     rcv_pid_append (restriction->rcv, pid, point);
   }
 }
      }
    }


    if (level > 0) {
      if (is_local(cell)) {

 Array pids = {NULL, 0, 0};
 if ((aparent(0,0,0).pid >= 0 && aparent(0,0,0).pid != pid()))
   append_pid (&pids, aparent(0,0,0).pid);
 int n = root_pids (parent, &pids);
 if (n) {
   for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
     rcv_pid_append (restriction->snd, *p, point);
   pfree (pids.p,__func__,__FILE__,__LINE__);
 }
      }
      else if ((cell.pid >= 0 && cell.pid != pid())) {

 if (is_local(aparent(0,0,0)) || has_local_child (parent))
   rcv_pid_append (restriction->rcv, cell.pid, point);
      }
    }
  }end_foreach_cell();}





  static const unsigned short keep = 1 << (user + 1);
  for (int l = depth(); l >= 0; l--)
    {foreach_cell()
      if (level == l) {
 if (level > 0 && (cell.pid < 0 || is_local(cell) || (cell.flags & used)))
   aparent(0,0,0).flags |= keep;
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) && !(cell.flags & keep))
   coarsen_cell (point, NULL);
 cell.flags &= ~(used|keep);
 continue;
      }end_foreach_cell();}


  m->send->len = m->receive->len = 0;
  rcv_pid_append_pids (mpi_level->snd, m->send);
  rcv_pid_append_pids (mpi_level_root->snd, m->send);
  rcv_pid_append_pids (mpi_level->rcv, m->receive);
  rcv_pid_append_pids (mpi_level_root->rcv, m->receive);

  prof_stop();
#line 1015 "/home/gongminjiang/basilisk/src/grid/tree-mpi.h"
end_tracing("mpi_boundary_update_buffers","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",1015);}

     
void mpi_boundary_refine (scalar * list)
{tracing("mpi_boundary_refine","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",1018);
  prof_start ("mpi_boundary_refine");

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;


  Array * snd = mpi->send;
  MPI_Request r[2*snd->len/sizeof(int)];
  int nr = 0;
  for (int i = 0, * dest = snd->p; i < snd->len/sizeof(int); i++,dest++) {
    int len = ((Tree *)grid)->refined.n;
    MPI_Isend (&((Tree *)grid)->refined.n, 1, MPI_INT, *dest,
        (128), MPI_COMM_WORLD, &r[nr++]);
    if (len > 0)
      MPI_Isend (((Tree *)grid)->refined.p, sizeof(Index)/sizeof(int)*len,
   MPI_INT, *dest, (128), MPI_COMM_WORLD, &r[nr++]);
  }



  Array * rcv = mpi->receive;
  Cache rerefined = {NULL, 0, 0};
  for (int i = 0, * source = rcv->p; i < rcv->len/sizeof(int); i++,source++) {
    int len;
    mpi_recv_check (&len, 1, MPI_INT, *source, (128),
      MPI_COMM_WORLD, MPI_STATUS_IGNORE,
      "mpi_boundary_refine (len)");
    if (len > 0) {
      Index p[len];
      mpi_recv_check (p, sizeof(Index)/sizeof(int)*len,
        MPI_INT, *source, (128),
        MPI_COMM_WORLD, MPI_STATUS_IGNORE,
        "mpi_boundary_refine (p)");
      Cache refined = {p, len, len};
      {foreach_cache (refined)
 if (level <= depth() && allocated(0,0,0)) {
   if (is_leaf(cell)) {
     bool neighbors = false;
     {foreach_neighbor()
       if (allocated(0,0,0) && (is_active(cell) || is_local(aparent(0,0,0)))) {
  neighbors = true; foreach_neighbor_break();
       }end_foreach_neighbor()}

     if (neighbors)
       refine_cell (point, list, 0, &rerefined);
   }
 }end_foreach_cache();}
    }
  }


  if (nr)
    MPI_Waitall (nr, r, MPI_STATUSES_IGNORE);


  pfree (((Tree *)grid)->refined.p,__func__,__FILE__,__LINE__);
  ((Tree *)grid)->refined = rerefined;

  prof_stop();



  mpi_all_reduce (rerefined.n, MPI_INT, MPI_SUM);
  if (rerefined.n)
    mpi_boundary_refine (list);
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].dirty = true;}}
end_tracing("mpi_boundary_refine","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",1086);}

static void check_depth()
{
#line 1121 "/home/gongminjiang/basilisk/src/grid/tree-mpi.h"
}

typedef struct {
  int refined, leaf;
} Remote;



     
void mpi_boundary_coarsen (int l, int too_fine)
{tracing("mpi_boundary_coarsen","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",1130);
  if (npe() == 1)
    {end_tracing("mpi_boundary_coarsen","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",1133);return;}

  check_depth();

  if (!(sizeof(Remote) == sizeof(double))) qassert ("/home/gongminjiang/basilisk/src/grid/tree-mpi.h", 1137, "sizeof(Remote) == sizeof(double)");

  scalar  remote=new_scalar("remote");
  {foreach_cell() {
    if (level == l) {
      if (is_local(cell)) {
 ((Remote *)&val(remote,0,0,0))->refined = (!is_leaf (cell) && cell.neighbors && cell.pid >= 0);
 ((Remote *)&val(remote,0,0,0))->leaf = is_leaf(cell);
      }
      else {
 ((Remote *)&val(remote,0,0,0))->refined = true;
 ((Remote *)&val(remote,0,0,0))->leaf = false;
      }
      continue;
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}
  mpi_boundary_level (mpi_boundary, ((scalar[]){remote,{-1}}), l);

  {foreach_cell() {
    if (level == l) {
      if (!is_local(cell)) {
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) && !((Remote *)&val(remote,0,0,0))->refined)
   coarsen_cell_recursive (point, NULL);
 else if (is_leaf(cell) && cell.neighbors && ((Remote *)&val(remote,0,0,0))->leaf) {
   int pid = cell.pid;
   {foreach_child()
     cell.pid = pid;end_foreach_child()}
 }
      }
      continue;
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  check_depth();

  if (l > 0) {
    {foreach_cell() {
      if (level == l) {
 val(remote,0,0,0) = is_local(cell) ? cell.neighbors : 0;
 continue;
      }
      if (is_leaf(cell))
 continue;
    }end_foreach_cell();}
    mpi_boundary_level (mpi_boundary, ((scalar[]){remote,{-1}}), l);
    {foreach_cell() {
      if (level == l)
 if (!is_local(cell) && is_local(aparent(0,0,0)) && val(remote,0,0,0)) {
   aparent(0,0,0).flags &= ~too_fine;
   continue;
 }
      if (is_leaf(cell))
 continue;
    }end_foreach_cell();}
  }delete((scalar*)((scalar[]){remote,{-1}}));
end_tracing("mpi_boundary_coarsen","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",1196);}

static void flag_border_cells()
{
  {foreach_cell() {
    if (is_active(cell)) {
      short flags = cell.flags & ~border;
      {foreach_neighbor() {
 if (!is_local(cell) || (level > 0 && !is_local(aparent(0,0,0)))) {
   flags |= border; foreach_neighbor_break();
 }

 if (is_refined_check())
   {foreach_child()
     if (!is_local(cell)) {
       flags |= border; foreach_child_break();
     }end_foreach_child()}
 if (flags & border)
   foreach_neighbor_break();
      }end_foreach_neighbor()}
      cell.flags = flags;
    }
    else {
      cell.flags &= ~border;

    }
    if (is_leaf(cell)) {
      if (cell.neighbors) {
 {foreach_child()
   cell.flags &= ~border;end_foreach_child()}
 if (is_border(cell)) {
   bool remote = false;
   {foreach_neighbor (2/2)
     if (!is_local(cell)) {
       remote = true; foreach_neighbor_break();
     }end_foreach_neighbor()}
   if (remote)
     {foreach_child()
       cell.flags |= border;end_foreach_child()}
 }
      }
      continue;
    }
  }end_foreach_cell();}
}

static int balanced_pid (long index, long nt, int nproc)
{
  long ne = max(1, nt/nproc), nr = nt % nproc;
  int pid = index < nr*(ne + 1) ?
    index/(ne + 1) :
    nr + (index - nr*(ne + 1))/ne;
  return min(nproc - 1, pid);
}


     
void mpi_partitioning()
{tracing("mpi_partitioning","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",1253);
  prof_start ("mpi_partitioning");

  long nt = 0;
  
#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
{
#line 1258
foreach ()
    nt++;end_foreach();}
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif



  
#line 1262
long i = 0;
  ((Tree *)grid)->dirty = true;
  {foreach_cell_post (is_active (cell))
    if (is_active (cell)) {
      if (is_leaf (cell)) {
 cell.pid = balanced_pid (i++, nt, npe());
 if (cell.neighbors > 0) {
   int pid = cell.pid;
   {foreach_child()
     cell.pid = pid;end_foreach_child()}
 }
 if (!is_local(cell))
   cell.flags &= ~active;
      }
      else {
 cell.pid = child(0,0,0).pid;
 bool inactive = true;
 {foreach_child()
   if (is_active(cell)) {
     inactive = false; foreach_child_break();
   }end_foreach_child()}
 if (inactive)
   cell.flags &= ~active;
      }
    }end_foreach_cell_post();}

  flag_border_cells();

  prof_stop();

  mpi_boundary_update_buffers();
end_tracing("mpi_partitioning","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",1293);}

void restore_mpi (FILE * fp, scalar * list1)
{
  long index = 0, nt = 0, start = ftell (fp);
  scalar  size=new_scalar("size"), * list = list_concat (((scalar[]){size,{-1}}), list1);;
  long offset = sizeof(double)*list_len(list);


  static const unsigned short set = 1 << user;
  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector[]){fm,{{-1},{-1}}});
  {foreach_cell()
    if (balanced_pid (index, nt, npe()) <= pid()) {
      unsigned flags;
      if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
 fprintf (ferr, "restore(): error: expecting 'flags'\n");
 exit (1);
      }
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
 double val;
 if (fread (&val, sizeof(double), 1, fp) != 1) {
   fprintf (ferr, "restore(): error: expecting scalar\n");
   exit (1);
 }
 if (s.i != INT_MAX)
   val(s,0,0,0) = val;
      }}}
      if (level == 0)
 nt = val(size,0,0,0);
      cell.pid = balanced_pid (index, nt, npe());
      cell.flags |= set;
      if (!(flags & leaf) && is_leaf(cell)) {
 if (balanced_pid (index + val(size,0,0,0) - 1, nt, npe()) < pid()) {
   fseek (fp, (sizeof(unsigned) + offset)*(val(size,0,0,0) - 1), SEEK_CUR);
   index += val(size,0,0,0);
   continue;
 }
 refine_cell (point, listm, 0, NULL);
      }
      index++;
      if (is_leaf(cell))
 continue;
    }end_foreach_cell();}


  fseek (fp, start, SEEK_SET);
  index = 0;
  {foreach_cell() {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (ferr, "restore(): error: expecting 'flags'\n");
      exit (1);
    }
    if (cell.flags & set)
      fseek (fp, offset, SEEK_CUR);
    else {
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
 double val;
 if (fread (&val, sizeof(double), 1, fp) != 1) {
   fprintf (ferr, "restore(): error: expecting a scalar\n");
   exit (1);
 }
 if (s.i != INT_MAX)
   val(s,0,0,0) = val;
      }}}
      cell.pid = balanced_pid (index, nt, npe());
      if (is_leaf(cell) && cell.neighbors) {
 int pid = cell.pid;
 {foreach_child()
   cell.pid = pid;end_foreach_child()}
      }
    }
    if (!(flags & leaf) && is_leaf(cell)) {
      bool locals = false;
      {foreach_neighbor(1)
 if ((cell.flags & set) && (is_local(cell) || is_root(point))) {
   locals = true; foreach_neighbor_break();
 }end_foreach_neighbor()}
      if (locals)
 refine_cell (point, listm, 0, NULL);
      else {
 fseek (fp, (sizeof(unsigned) + offset)*(val(size,0,0,0) - 1), SEEK_CUR);
 index += val(size,0,0,0);
 continue;
      }
    }
    index++;
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}


  {foreach_cell_post (is_active (cell)) {
    cell.flags &= ~set;
    if (is_active (cell)) {
      if (is_leaf (cell)) {
 if (cell.neighbors > 0) {
   int pid = cell.pid;
   {foreach_child()
     cell.pid = pid;end_foreach_child()}
 }
 if (!is_local(cell))
   cell.flags &= ~active;
      }
      else if (!is_local(cell)) {
 bool inactive = true;
 {foreach_child()
   if (is_active(cell)) {
     inactive = false; foreach_child_break();
   }end_foreach_child()}
 if (inactive)
   cell.flags &= ~active;
      }
    }
  }end_foreach_cell_post();}

  flag_border_cells();

  mpi_boundary_update (list);
  pfree (list,__func__,__FILE__,__LINE__);delete((scalar*)((scalar[]){size,{-1}}));
}
#line 1435 "/home/gongminjiang/basilisk/src/grid/tree-mpi.h"
     
double z_indexing (scalar index, bool leaves)
{tracing("z_indexing","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",1436);



  scalar  size=new_scalar("size");
  subtree_size (size, leaves);






  double maxi = -1.;
  if (pid() == 0)
    {foreach_level(0)
      maxi = val(size,0,0,0) - 1.;end_foreach_level();}




  {foreach_level(0)
    val(index,0,0,0) = 0;end_foreach_level();}
  for (int l = 0; l < depth(); l++) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar[]){index,{-1}}), l); };
    {foreach_cell() {
      if (level == l) {
 if (is_leaf(cell)) {
   if (is_local(cell) && cell.neighbors) {
     int i = val(index,0,0,0);
     {foreach_child()
       val(index,0,0,0) = i;end_foreach_child()}
   }
 }
 else {
   bool loc = is_local(cell);
   if (!loc)
     {foreach_child()
       if (is_local(cell)) {
  loc = true; foreach_child_break();
       }end_foreach_child()}
   if (loc) {
     int i = val(index,0,0,0) + !leaves;
     {foreach_child() {
       val(index,0,0,0) = i;
       i += val(size,0,0,0);
     }end_foreach_child()}
   }
 }
 continue;
      }
      if (is_leaf(cell))
 continue;
    }end_foreach_cell();}
  }
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar[]){index,{-1}}), depth()); };

  {delete((scalar*)((scalar[]){size,{-1}}));{end_tracing("z_indexing","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",1493);return maxi;}}delete((scalar*)((scalar[]){size,{-1}}));
end_tracing("z_indexing","/home/gongminjiang/basilisk/src/grid/tree-mpi.h",1494);}
#line 1687 "/home/gongminjiang/basilisk/src/grid/tree.h"
#line 1 "grid/balance.h"
#line 1 "/home/gongminjiang/basilisk/src/grid/balance.h"


typedef struct {
  short leaf, prolongation;
  int pid;
} NewPid;



#if TRASH
# define is_newpid() (!isnan(val(newpid,0,0,0)) && ((NewPid *)&val(newpid,0,0,0))->pid > 0)
#else
# define is_newpid() (((NewPid *)&val(newpid,0,0,0))->pid > 0)
#endif

Array * linear_tree (size_t size, scalar newpid)
{
  const unsigned short sent = 1 << user, next = 1 << (user + 1);
  Array * a = array_new();

  {foreach_cell_post_all (true)
    if (level > 0 && (cell.flags & (sent|next)))
      aparent(0,0,0).flags |= next;end_foreach_cell_post_all();}

  bool empty = true;
  {foreach_cell_all() {
    if (cell.flags & sent) {
      array_append (a, &cell, size);
      cell.flags &= ~sent;
      empty = false;
    }
    else {
      if (cell.pid >= 0 && ((NewPid *)&val(newpid,0,0,0))->leaf)
 if (!(is_leaf(cell))) qassert ("/home/gongminjiang/basilisk/src/grid/balance.h", 34, "is_leaf(cell)");
      if (is_refined_check()) {


 bool prolo = false;
 {foreach_child()
   if (((NewPid *)&val(newpid,0,0,0))->prolongation)
     prolo = true;end_foreach_child()}
 if (prolo) {

   cell.flags |= leaf;
   array_append (a, &cell, sizeof(Cell));
   cell.flags &= ~leaf;
 }
 else
   array_append (a, &cell, sizeof(Cell));
      }
      else
 array_append (a, &cell, sizeof(Cell));
    }
    if (cell.flags & next)
      cell.flags &= ~next;
    else
      continue;
  }end_foreach_cell_all();}

  if (empty)
    a->len = 0;
  return a;
}

#define foreach_tree(t, size, list)\
{\
  const unsigned short _sent = 1 << user, _next = 1 << (user + 1);\
  scalar * _list = list;\
  char * _i = (char *) (t)->p;\
  foreach_cell_all() {\
    Cell * c = (Cell *) _i;\
    if (c->flags & _sent) {\
      _i += size;\

#line 74


#define end_foreach_tree()\
    }\
    else\
      _i += sizeof(Cell);\
    if (c->flags & _next) {\
      if (!(c->neighbors)) qassert ("/home/gongminjiang/basilisk/src/grid/balance.h", 81, "c->neighbors");\
      if (!(c->flags & leaf) && is_leaf(cell) &&\
   (!is_newpid() || !((NewPid *)&val(newpid,0,0,0))->leaf))\
\
 refine_cell (point, _list, 0, NULL);\
      else if (!cell.neighbors)\
\
 alloc_children (point);\
    }\
    else\
      continue;\
  } end_foreach_cell_all();\
}\

#line 94


Array * neighborhood (scalar newpid, int nextpid, FILE * fp)
{
  const unsigned short sent = 1 << user;
  {foreach_cell() {

    bool root = false;
    if ((!is_local(cell) || ((NewPid *)&val(newpid,0,0,0))->pid - 1 != nextpid) && (!is_leaf (cell) && cell.neighbors && cell.pid >= 0)) {
      {foreach_child()
 if (is_local(cell) && ((NewPid *)&val(newpid,0,0,0))->pid - 1 == nextpid) {
   root = true; foreach_child_break();
 }end_foreach_child()}
      if (root && cell.pid != nextpid) {
 {foreach_neighbor()
   if (cell.pid != nextpid && is_newpid()) {
     if (fp)
       fprintf (fp, "%g %g %g %d %d root\n",
         x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
     cell.flags |= sent;
   }end_foreach_neighbor()}
      }
    }

    if ((is_local(cell) && ((NewPid *)&val(newpid,0,0,0))->pid - 1 == nextpid) || root) {
      {foreach_neighbor(1)
 if (cell.neighbors && cell.pid != nextpid)
   {foreach_child()
     if (cell.pid != nextpid && is_newpid()) {
       if (fp)
  fprintf (fp, "%g %g %g %d %d nextpid\n",
    x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
       cell.flags |= sent;
     }end_foreach_child()}end_foreach_neighbor()}
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  return linear_tree (sizeof(Cell) + datasize, newpid);
}

static void send_tree (Array * a, int to, MPI_Request * r)
{
  MPI_Isend (&a->len, 1, MPI_LONG, to, (256), MPI_COMM_WORLD, &r[0]);
  if (a->len > 0) {
    MPI_Isend (a->p, a->len, MPI_BYTE, to, (256), MPI_COMM_WORLD, &r[1]);
    ((Tree *)grid)->dirty = true;
  }
}

static void receive_tree (int from, scalar newpid, FILE * fp)
{
  Array a;
  mpi_recv_check (&a.len, 1, MPI_LONG, from, (256),
    MPI_COMM_WORLD, MPI_STATUS_IGNORE, "receive_tree (len)");
  if (a.len > 0) {
    a.p = pmalloc (a.len,__func__,__FILE__,__LINE__);
    if (fp)
      fprintf (fp, "receiving %ld from %d\n", a.len, from);
    mpi_recv_check (a.p, a.len, MPI_BYTE, from, (256),
      MPI_COMM_WORLD, MPI_STATUS_IGNORE, "receive_tree (p)");

    {foreach_tree (&a, sizeof(Cell) + datasize, NULL) {
      memcpy (((char *)&cell) + sizeof(Cell), ((char *)c) + sizeof(Cell),
       datasize);
      if (!(((NewPid *)&val(newpid,0,0,0))->pid > 0)) qassert ("/home/gongminjiang/basilisk/src/grid/balance.h", 160, "NEWPID()->pid > 0");
      if (fp)
 fprintf (fp, "%g %g %g %d %d %d %d %d %d recv\n",
   x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid,
   c->flags & leaf,
   cell.flags & leaf, from, ((NewPid *)&val(newpid,0,0,0))->leaf);
    }end_foreach_tree();}
    pfree (a.p,__func__,__FILE__,__LINE__);
    ((Tree *)grid)->dirty = true;
  }
}

static void wait_tree (Array * a, MPI_Request * r)
{
  MPI_Wait (&r[0], MPI_STATUS_IGNORE);
  if (a->len > 0)
    MPI_Wait (&r[1], MPI_STATUS_IGNORE);
}

static void check_flags()
{







}

struct {
  int min;
  bool leaves;

  int npe;
} mpi = {
  1,
  true
};

     
bool balance()
{tracing("balance","/home/gongminjiang/basilisk/src/grid/balance.h",201);
  if (npe() == 1)
    {end_tracing("balance","/home/gongminjiang/basilisk/src/grid/balance.h",204);return false;}

  if (!(sizeof(NewPid) == sizeof(double))) qassert ("/home/gongminjiang/basilisk/src/grid/balance.h", 206, "sizeof(NewPid) == sizeof(double)");

  check_flags();

  long nl = 0, nt = 0;
  {foreach_cell() {
    if (is_local(cell)) {
      nt++;
      if (is_leaf(cell))
 nl++;
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  grid->n = grid->tn = nl;
  grid->maxdepth = depth();
  long nmin = nl, nmax = nl;

  mpi_all_reduce (nmax, MPI_LONG, MPI_MAX);
  mpi_all_reduce (nmin, MPI_LONG, MPI_MIN);
  mpi_all_reduce (grid->tn, MPI_LONG, MPI_SUM);
  mpi_all_reduce (grid->maxdepth, MPI_INT, MPI_MAX);
  if (mpi.leaves)
    nt = grid->tn;
  else
    mpi_all_reduce (nt, MPI_LONG, MPI_SUM);

  long ne = max(1, nt/npe());

  if (ne < mpi.min) {
    mpi.npe = max(1, nt/mpi.min);
    ne = max(1, nt/mpi.npe);
  }
  else
    mpi.npe = npe();

  if (nmax - nmin <= 1)
    {end_tracing("balance","/home/gongminjiang/basilisk/src/grid/balance.h",244);return false;}

  scalar  newpid=new_scalar("newpid");
  double zn = z_indexing (newpid, mpi.leaves);
  if (pid() == 0)
    if (!(zn + 1 == nt)) qassert ("/home/gongminjiang/basilisk/src/grid/balance.h", 249, "zn + 1 == nt");

  FILE * fp = NULL;
#line 261 "/home/gongminjiang/basilisk/src/grid/balance.h"
  bool next = false, prev = false;
  {foreach_cell_all() {
    if (is_local(cell)) {
      int pid = balanced_pid (val(newpid,0,0,0), nt, mpi.npe);
      pid = clamp (pid, cell.pid - 1, cell.pid + 1);
      if (pid == pid() + 1)
 next = true;
      else if (pid == pid() - 1)
 prev = true;
      ((NewPid *)&val(newpid,0,0,0))->pid = pid + 1;
      ((NewPid *)&val(newpid,0,0,0))->leaf = is_leaf(cell);
      ((NewPid *)&val(newpid,0,0,0))->prolongation = (!is_leaf(cell) && !cell.neighbors && cell.pid >= 0);
      if (fp)
 fprintf (fp, "%g %g %d %d newpid\n", x, y, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
    }
    else
      val(newpid,0,0,0) = 0;
  }end_foreach_cell_all();}
  for (int l = 0; l <= depth(); l++)
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, ((scalar[]){newpid,{-1}}), l); };
#line 305 "/home/gongminjiang/basilisk/src/grid/balance.h"
  Array * anext = next ? neighborhood (newpid, pid() + 1, fp) : array_new();
  Array * aprev = prev ? neighborhood (newpid, pid() - 1, fp) : array_new();

  if (fp)
    fflush (fp);

  check_flags();


  MPI_Request rprev[2], rnext[2];
  if (pid() > 0)
    send_tree (aprev, pid() - 1, rprev);
  if (pid() < npe() - 1)
    send_tree (anext, pid() + 1, rnext);


  if (pid() < npe() - 1)
    receive_tree (pid() + 1, newpid, fp);
  if (pid() > 0)
    receive_tree (pid() - 1, newpid, fp);


  if (pid() > 0)
    wait_tree (aprev, rprev);
  array_free (aprev);
  if (pid() < npe() - 1)
    wait_tree (anext, rnext);
  array_free (anext);

  if (fp)
    fflush (fp);


  int pid_changed = false;
  {foreach_cell_all() {
    if (cell.pid >= 0) {
      if (is_newpid()) {
 if (fp)
   fprintf (fp, "%g %g %g %d %d %d %d %d new\n",
     x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid,
     is_leaf(cell), cell.neighbors, ((NewPid *)&val(newpid,0,0,0))->leaf);
 if (cell.pid != ((NewPid *)&val(newpid,0,0,0))->pid - 1) {
   cell.pid = ((NewPid *)&val(newpid,0,0,0))->pid - 1;
   cell.flags &= ~(active|border);
   if (is_local(cell))
     cell.flags |= active;
   pid_changed = true;
 }
 if (((NewPid *)&val(newpid,0,0,0))->leaf && !is_leaf(cell) && cell.neighbors)
   coarsen_cell_recursive (point, NULL);
      }
      else if (level > 0 && ((NewPid *)&coarse(newpid,0,0,0))->leaf)
 cell.pid = aparent(0,0,0).pid;
    }

    if (!cell.neighbors && allocated_child(0,0,0)) {
      if (fp)
 fprintf (fp, "%g %g %g %d %d freechildren\n",
   x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
      free_children (point);
    }
  }end_foreach_cell_all();}

  if (((Tree *)grid)->dirty || pid_changed) {


    {foreach_cell_post (!is_leaf (cell))
      if (!is_leaf(cell) && !is_local(cell)) {
 unsigned short flags = cell.flags & ~active;
 {foreach_child()
   if (is_active(cell)) {
     flags |= active; foreach_child_break();
   }end_foreach_child()}
 cell.flags = flags;
      }end_foreach_cell_post();}

    flag_border_cells();
    pid_changed = true;
  }

  if (fp)
    fclose (fp);

  mpi_all_reduce (pid_changed, MPI_INT, MPI_MAX);
  if (pid_changed)
    mpi_boundary_update_buffers();

  {delete((scalar*)((scalar[]){newpid,{-1}}));{end_tracing("balance","/home/gongminjiang/basilisk/src/grid/balance.h",392);return pid_changed;}}delete((scalar*)((scalar[]){newpid,{-1}}));
end_tracing("balance","/home/gongminjiang/basilisk/src/grid/balance.h",393);}

void mpi_boundary_update (scalar * list)
{
  mpi_boundary_update_buffers();
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].dirty = true;}}
  grid->tn = 0;
  boundary_internal ((scalar *)list, "/home/gongminjiang/basilisk/src/grid/balance.h", 401);
  while (balance());
}
#line 1688 "/home/gongminjiang/basilisk/src/grid/tree.h"
#else
void mpi_boundary_refine (scalar * list){}
void mpi_boundary_coarsen (int a, int b){}
void mpi_boundary_update (scalar * list) {
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].dirty = true;}}
  boundary_internal ((scalar *)list, "/home/gongminjiang/basilisk/src/grid/tree.h", 1694);
}
#endif
#line 4 "/home/gongminjiang/basilisk/src/grid/quadtree.h"

void quadtree_methods() {
  tree_methods();
}
#line 15 "atomisation-cpp.c"
#line 1 "atomisation.c"
#line 14 "atomisation.c"
#line 1 "navier-stokes/centered.h"
#line 1 "/home/gongminjiang/basilisk/src/navier-stokes/centered.h"
#line 27 "/home/gongminjiang/basilisk/src/navier-stokes/centered.h"
#line 1 "./run.h"
#line 1 "/home/gongminjiang/basilisk/src/run.h"
#line 9 "/home/gongminjiang/basilisk/src/run.h"
double dt = 1.;

#line 1 "./utils.h"
#line 1 "/home/gongminjiang/basilisk/src/utils.h"







double DT = 1e10, CFL = 0.5;




struct {

  long nc;

  long tnc;

  double t;

  double speed;

  timer gt;
} perf;





void update_perf() {
  perf.nc += grid->n;
  perf.tnc += grid->tn;
  perf.t = timer_elapsed (perf.gt);
  perf.speed = perf.tnc/perf.t;
}






typedef struct {
  double cpu;
  double real;
  double speed;
  double min;
  double avg;
  double max;
  size_t tnc;
  long mem;
} timing;






timing timer_timing (timer t, int i, size_t tnc, double * mpi)
{
  timing s;
#if 1
  s.avg = mpi_time - t.tm;
#endif
  clock_t end = clock();
  s.cpu = ((double) (end - t.c))/CLOCKS_PER_SEC;
  s.real = timer_elapsed (t);
  if (tnc == 0) {
    double n = 0;
    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:n)){
#line 69
foreach() n++;end_foreach();mpi_all_reduce_array(&n,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
    
#line 70
s.tnc = n;
    tnc = n*i;
  }
  else
    s.tnc = tnc;
#if _GNU_SOURCE
  struct rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  s.mem = usage.ru_maxrss;
#else
  s.mem = 0;
#endif
#if 1
  if (mpi)
    MPI_Allgather (&s.avg, 1, MPI_DOUBLE, mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  s.max = s.min = s.avg;
  mpi_all_reduce (s.max, MPI_DOUBLE, MPI_MAX);
  mpi_all_reduce (s.min, MPI_DOUBLE, MPI_MIN);
  mpi_all_reduce (s.avg, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.real, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.mem, MPI_LONG, MPI_SUM);
  s.real /= npe();
  s.avg /= npe();
  s.mem /= npe();
#else
  s.min = s.max = s.avg = 0.;
#endif
  s.speed = s.real > 0. ? tnc/s.real : -1.;
  return s;
}




void timer_print (timer t, int i, size_t tnc)
{
  timing s = timer_timing (t, i, tnc, NULL);
  fprintf (fout,
    "\n# " "Quadtree"
    ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
    i, s.cpu, s.real, s.speed, (int) (datasize/sizeof(double)));
#if 1
  fprintf (fout,
    "# %d procs, MPI: min %.2g (%.2g%%) "
    "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
    npe(),
    s.min, 100.*s.min/s.real,
    s.avg, 100.*s.avg/s.real,
    s.max, 100.*s.max/s.real);
#endif
}







typedef struct {
  double avg, rms, max, volume;
} norm;

norm normf (scalar f)
{
  double avg = 0., rms = 0., max = 0., volume = 0.;
  foreach_stencil(
)
    {_stencil_val(f,0,0,0);_stencil_val(cm,0,0,0); {   
      _stencil_val(f,0,0,0);   
         
_stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0); 
       
    
#line 143
}       }end_foreach_stencil();
  
#line 135
if(!is_constant(cm)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(+:volume)
   reduction(+:rms) reduction(+:avg)reduction(max:max)){
#line 135
foreach(
)
    if (val(f,0,0,0) != HUGE && (sq(Delta)*val(cm,0,0,0)) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > max) max = v;
      volume += (sq(Delta)*val(cm,0,0,0));
      avg += (sq(Delta)*val(cm,0,0,0))*v;
      rms += (sq(Delta)*val(cm,0,0,0))*sq(v);
    }end_foreach();mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&rms,double,MPI_SUM,1);mpi_all_reduce_array(&avg,double,MPI_SUM,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 143
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(+:volume)
   reduction(+:rms) reduction(+:avg)reduction(max:max)){
#line 135
foreach(
)
    if (val(f,0,0,0) != HUGE && (sq(Delta)*_const_cm) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > max) max = v;
      volume += (sq(Delta)*_const_cm);
      avg += (sq(Delta)*_const_cm)*v;
      rms += (sq(Delta)*_const_cm)*sq(v);
    }end_foreach();mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&rms,double,MPI_SUM,1);mpi_all_reduce_array(&avg,double,MPI_SUM,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 143
}
  norm n;
  n.avg = volume ? avg/volume : 0.;
  n.rms = volume ? sqrt(rms/volume) : 0.;
  n.max = max;
  n.volume = volume;
  return n;
}





typedef struct {
  double min, max, sum, stddev, volume;
} stats;

stats statsf (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
  foreach_stencil(
)
    {_stencil_val(cm,0,0,0); _stencil_val(f,0,0,0); {
_stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0);_stencil_val(f,0,0,0); 
       _stencil_val(cm,0,0,0);_stencil_val(f,0,0,0); 
       _stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
_stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
         
         
    
#line 171
}      }end_foreach_stencil();
  
#line 163
if(!is_constant(cm)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:min)
   reduction(max:max) reduction(+:volume) reduction(+:sum2)reduction(+:sum)){
#line 163
foreach(
)
    if ((sq(Delta)*val(cm,0,0,0)) > 0. && val(f,0,0,0) != HUGE) {
      volume += (sq(Delta)*val(cm,0,0,0));
      sum += (sq(Delta)*val(cm,0,0,0))*val(f,0,0,0);
      sum2 += (sq(Delta)*val(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > max) max = val(f,0,0,0);
      if (val(f,0,0,0) < min) min = val(f,0,0,0);
    }end_foreach();mpi_all_reduce_array(&min,double,MPI_MIN,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&sum2,double,MPI_SUM,1);mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 171
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:min)
   reduction(max:max) reduction(+:volume) reduction(+:sum2)reduction(+:sum)){
#line 163
foreach(
)
    if ((sq(Delta)*_const_cm) > 0. && val(f,0,0,0) != HUGE) {
      volume += (sq(Delta)*_const_cm);
      sum += (sq(Delta)*_const_cm)*val(f,0,0,0);
      sum2 += (sq(Delta)*_const_cm)*sq(val(f,0,0,0));
      if (val(f,0,0,0) > max) max = val(f,0,0,0);
      if (val(f,0,0,0) < min) min = val(f,0,0,0);
    }end_foreach();mpi_all_reduce_array(&min,double,MPI_MIN,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&sum2,double,MPI_SUM,1);mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 171
}
  stats s;
  s.min = min, s.max = max, s.sum = sum, s.volume = volume;
  if (volume > 0.)
    sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
}
#line 187 "/home/gongminjiang/basilisk/src/utils.h"
static double generic_limiter (double r, double beta)
{
  double v1 = min (r, beta), v2 = min (beta*r, 1.);
  v1 = max (0., v1);
  return max (v1, v2);
}

double minmod (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.)*(s1 - s0);
}

double superbee (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 2.)*(s1 - s0);
}

double sweby (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.5)*(s1 - s0);
}
#line 213 "/home/gongminjiang/basilisk/src/utils.h"
double theta = 1.3;

double minmod2 (double s0, double s1, double s2)
{
  if (s0 < s1 && s1 < s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 < d1) d1 = d2;
    return min(d1, d3);
  }
  if (s0 > s1 && s1 > s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 > d1) d1 = d2;
    return max(d1, d3);
  }
  return 0.;
}
#line 237 "/home/gongminjiang/basilisk/src/utils.h"
void gradients (scalar * f, vector * g)
{
  if (!(list_len(f) == vectors_len(g))) qassert ("/home/gongminjiang/basilisk/src/utils.h", 239, "list_len(f) == vectors_len(g)");
  foreach_stencil() {
    scalar s; vector v;
    {vector*_i0=g;scalar*_i1= f;if(_i0)for(v=*_i0,s=*_i1;_i0->x.i>= 0;v=*++_i0,s=*++_i1){ {
      if (_attribute[s.i].gradient)
 { {





     _stencil_val_a(v.x,0,0,0);_stencil_val(s,-1,0,0); _stencil_val(s,0,0,0); _stencil_val(s,1,0,0);   
 } 
#line 244
{





     _stencil_val_a(v.y,0,0,0);_stencil_val(s,0,-1,0); _stencil_val(s,0,0,0); _stencil_val(s,0,1,0);   
 }}
      else
 { {





     _stencil_val_a(v.x,0,0,0);_stencil_val(s,1,0,0); _stencil_val(s,-1,0,0);   
 } 
#line 253
{





     _stencil_val_a(v.y,0,0,0);_stencil_val(s,0,1,0); _stencil_val(s,0,-1,0);   
 }}
    }}}
  }end_foreach_stencil();
  {
#line 240
foreach() {
    scalar s; vector v;
    {vector*_i0=g;scalar*_i1= f;if(_i0)for(v=*_i0,s=*_i1;_i0->x.i>= 0;v=*++_i0,s=*++_i1){ {
      if (_attribute[s.i].gradient)
 { {





     val(v.x,0,0,0) = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0))/Delta;
 } 
#line 244
{





     val(v.y,0,0,0) = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0))/Delta;
 }}
      else
 { {





     val(v.x,0,0,0) = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);
 } 
#line 253
{





     val(v.y,0,0,0) = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);
 }}
    }}}
  }end_foreach();}
}
#line 280 "/home/gongminjiang/basilisk/src/utils.h"
void vorticity (const vector u, scalar omega)
{
  foreach_stencil()
    {_stencil_val_a(omega,0,0,0);_stencil_val(fm.x,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,0,0,0);
        _stencil_val(fm.x,1,0,0);_stencil_val(u.y,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,-1,0,0);
_stencil_val(fm.y,0,1,0); _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,0,0);
        _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,-1,0); _stencil_val(fm.y,0,1,0);_stencil_val(u.x,0,1,0);_stencil_val(cm,0,0,0);      
             
#line 286
}end_foreach_stencil();
  
#line 282
if(!is_constant(fm.x) && !is_constant(cm)){{foreach()
    val(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +
        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +
        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*val(cm,0,0,0)*Delta + 0.);end_foreach();}}else if(is_constant(fm.x) && !is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  {
#line 282
foreach()
    val(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +
        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -
        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +
        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*val(cm,0,0,0)*Delta + 0.);end_foreach();}}else if(!is_constant(fm.x) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  {
#line 282
foreach()
    val(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +
        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +
        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*_const_cm*Delta + 0.);end_foreach();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  {
#line 282
foreach()
    val(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +
        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -
        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +
        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*_const_cm*Delta + 0.);end_foreach();}}
}





double change (scalar s, scalar sn)
{
  double max = 0.;
  foreach_stencil() {
_stencil_val(cm,0,0,0); {     
       _stencil_val(sn,0,0,0);_stencil_val(s,0,0,0);   
       
  
    }
       
    
#line 302
_stencil_val_a(sn,0,0,0); _stencil_val(s,0,0,0); 
  }end_foreach_stencil();
  
#line 296
if(!is_constant(cm)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:max)){
#line 296
foreach() {
    if ((sq(Delta)*val(cm,0,0,0)) > 0.) {
      double ds = fabs (val(s,0,0,0) - val(sn,0,0,0));
      if (ds > max)
 max = ds;
    }
    val(sn,0,0,0) = val(s,0,0,0);
  }end_foreach();mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 303
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:max)){
#line 296
foreach() {
    if ((sq(Delta)*_const_cm) > 0.) {
      double ds = fabs (val(s,0,0,0) - val(sn,0,0,0));
      if (ds > max)
 max = ds;
    }
    val(sn,0,0,0) = val(s,0,0,0);
  }end_foreach();mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 303
}
  return max;
}





scalar lookup_field (const char * name)
{
  if (name)
    {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (!strcmp (_attribute[s.i].name, name))
 return s;}}
  return (scalar){-1};
}

vector lookup_vector (const char * name)
{
  if (name) {
    char component[strlen(name) + 3];
    strcpy (component, name);
    strcat (component, ".x");
    {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (!strcmp (_attribute[s.i].name, component))
 return _attribute[s.i].v;}}
  }
  return (vector){{-1}};
}
#line 340 "/home/gongminjiang/basilisk/src/utils.h"
#define foreach_segment(_S,_p) {\
  coord t = {(_S)[1].x - (_S)[0].x, (_S)[1].y - (_S)[0].y};\
  double norm = sqrt(sq(t.x) + sq(t.y));\
  if (!(norm > 0.)) qassert ("/home/gongminjiang/basilisk/src/utils.h", 343, "norm > 0.");\
  t.x = t.x/norm + 1e-6, t.y = t.y/norm - 1.5e-6;\
  double alpha = ((_S)[0].x*((_S)[1].y - (_S)[0].y) -\
    (_S)[0].y*((_S)[1].x - (_S)[0].x))/norm;\
  foreach()\
    if (fabs(t.y*x - t.x*y - alpha) < 0.708*Delta) {\
      coord _o = {x,y}, _p[2];\
      int _n = 0;\
 if (t.x)\
   for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {\
     _p[_n].x = _o.x + _i*Delta/2.;\
     double a = (_p[_n].x - (_S)[0].x)/t.x;\
     _p[_n].y = (_S)[0].y + a*t.y;\
     if (fabs(_p[_n].y - _o.y) <= Delta/2.) {\
       a = clamp (a, 0., norm);\
       _p[_n].x = (_S)[0].x + a*t.x, _p[_n].y = (_S)[0].y + a*t.y;\
       if (fabs(_p[_n].x - _o.x) <= Delta/2. &&\
    fabs(_p[_n].y - _o.y) <= Delta/2.)\
  _n++;\
     }\
   }\
\
 if (t.y)\
   for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {\
     _p[_n].y = _o.y + _i*Delta/2.;\
     double a = (_p[_n].y - (_S)[0].y)/t.y;\
     _p[_n].x = (_S)[0].x + a*t.x;\
     if (fabs(_p[_n].x - _o.x) <= Delta/2.) {\
       a = clamp (a, 0., norm);\
       _p[_n].y = (_S)[0].y + a*t.y, _p[_n].x = (_S)[0].x + a*t.x;\
       if (fabs(_p[_n].y - _o.y) <= Delta/2. &&\
    fabs(_p[_n].x - _o.x) <= Delta/2.)\
  _n++;\
     }\
   }\
\
      if (_n == 2) {\

#line 380

#line 410 "/home/gongminjiang/basilisk/src/utils.h"
#define end_foreach_segment() } } end_foreach(); }




void fields_stats()
{
  fprintf (ferr, "# t = %g, fields = {", t);
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    fprintf (ferr, " %s", _attribute[s.i].name);}}
  fputs (" }\n", ferr);
  fprintf (ferr, "# %12s: %12s %12s %12s %12s\n",
    "name", "min", "avg", "stddev", "max");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    stats ss = statsf (s);
    fprintf (ferr, "# %12s: %12g %12g %12g %12g\n",
      _attribute[s.i].name, ss.min, ss.sum/ss.volume, ss.stddev, ss.max);
  }}}
}

#line 1 "./output.h"
#line 1 "/home/gongminjiang/basilisk/src/output.h"
#line 37 "/home/gongminjiang/basilisk/src/output.h"
struct OutputField {
  scalar * list;
  FILE * fp;
  int n;
  bool linear;
  double box[2][2];
};

     
void output_field (struct OutputField p)
{tracing("output_field","/home/gongminjiang/basilisk/src/output.h",46);
  if (!p.list) p.list = all;
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = fout;
  p.n++;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }

  boundary_internal ((scalar *)p.list, "/home/gongminjiang/basilisk/src/output.h", 58);
  int len = list_len(p.list);
  double Delta = 0.999999*(p.box[1][0] - p.box[0][0])/(p.n - 1);
  int ny = (p.box[1][1] - p.box[0][1])/Delta + 1;
  double ** field = (double **) matrix_new (p.n, ny, len*sizeof(double));
  for (int i = 0; i < p.n; i++) {
    double x = Delta*i + p.box[0][0];
    for (int j = 0; j < ny; j++) {
      double y = Delta*j + p.box[0][1];
      if (p.linear) {
 int k = 0;
 {scalar*_i=(scalar*)( p.list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   field[i][len*j + k++] = interpolate ((struct _interpolate){s, x, y});}}
      }
      else {
 Point point = locate ((struct _locate){x, y});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
 int k = 0;
 {scalar*_i=(scalar*)( p.list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   field[i][len*j + k++] = point.level >= 0 ? val(s,0,0,0) : HUGE;}}
      }
    }
  }

  if (pid() == 0) {
#if 1
    MPI_Reduce (MPI_IN_PLACE, field[0], len*p.n*ny, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif
    fprintf (p.fp, "# 1:x 2:y");
    int i = 3;
    {scalar*_i=(scalar*)( p.list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      fprintf (p.fp, " %d:%s", i++, _attribute[s.i].name);}}
    fputc('\n', p.fp);
    for (int i = 0; i < p.n; i++) {
      double x = Delta*i + p.box[0][0];
      for (int j = 0; j < ny; j++) {
 double y = Delta*j + p.box[0][1];

 fprintf (p.fp, "%g %g", x, y);
 int k = 0;
 {scalar*_i=(scalar*)( p.list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   fprintf (p.fp, " %g", field[i][len*j + k++]);}}
 fputc ('\n', p.fp);
      }
      fputc ('\n', p.fp);
    }
    fflush (p.fp);
  }
#if 1
  else
    MPI_Reduce (field[0], NULL, len*p.n*ny, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (field);
end_tracing("output_field","/home/gongminjiang/basilisk/src/output.h",113);}
#line 141 "/home/gongminjiang/basilisk/src/output.h"
struct OutputMatrix {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
};

     
void output_matrix (struct OutputMatrix p)
{tracing("output_matrix","/home/gongminjiang/basilisk/src/output.h",149);
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = fout;
  if (p.linear) {
    scalar f = p.f;
    boundary_internal ((scalar *)((scalar[]){f,{-1}}), "/home/gongminjiang/basilisk/src/output.h", 155);
  }
  float fn = p.n;
  float Delta = (float) L0/fn;
  fwrite (&fn, sizeof(float), 1, p.fp);
  for (int j = 0; j < p.n; j++) {
    float yp = (float) (Delta*j + X0 + Delta/2.);
    fwrite (&yp, sizeof(float), 1, p.fp);
  }
  for (int i = 0; i < p.n; i++) {
    float xp = (float) (Delta*i + X0 + Delta/2.);
    fwrite (&xp, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = (float)(Delta*j + Y0 + Delta/2.), v;
      if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
 if (!(point.level >= 0)) qassert ("/home/gongminjiang/basilisk/src/output.h", 173, "point.level >= 0");
 v = val(p.f,0,0,0);
      }
      fwrite (&v, sizeof(float), 1, p.fp);
    }
  }
  fflush (p.fp);
end_tracing("output_matrix","/home/gongminjiang/basilisk/src/output.h",180);}
#line 189 "/home/gongminjiang/basilisk/src/output.h"
typedef void (* colormap) (double cmap[127][3]);

void jet (double cmap[127][3])
{
  for (int i = 0; i < 127; i++) {
    cmap[i][0] =
      i <= 46 ? 0. :
      i >= 111 ? -0.03125*(i - 111) + 1. :
      i >= 78 ? 1. :
      0.03125*(i - 46);
    cmap[i][1] =
      i <= 14 || i >= 111 ? 0. :
      i >= 79 ? -0.03125*(i - 111) :
      i <= 46 ? 0.03125*(i - 14) :
      1.;
    cmap[i][2] =
      i >= 79 ? 0. :
      i >= 47 ? -0.03125*(i - 79) :
      i <= 14 ? 0.03125*(i - 14) + 1.:
      1.;
  }
}

void cool_warm (double cmap[127][3])
{






  static double basemap[33][3] = {
    {0.2298057, 0.298717966, 0.753683153},
    {0.26623388, 0.353094838, 0.801466763},
    {0.30386891, 0.406535296, 0.84495867},
    {0.342804478, 0.458757618, 0.883725899},
    {0.38301334, 0.50941904, 0.917387822},
    {0.424369608, 0.558148092, 0.945619588},
    {0.46666708, 0.604562568, 0.968154911},
    {0.509635204, 0.648280772, 0.98478814},
    {0.552953156, 0.688929332, 0.995375608},
    {0.596262162, 0.726149107, 0.999836203},
    {0.639176211, 0.759599947, 0.998151185},
    {0.681291281, 0.788964712, 0.990363227},
    {0.722193294, 0.813952739, 0.976574709},
    {0.761464949, 0.834302879, 0.956945269},
    {0.798691636, 0.849786142, 0.931688648},
    {0.833466556, 0.860207984, 0.901068838},
    {0.865395197, 0.86541021, 0.865395561},
    {0.897787179, 0.848937047, 0.820880546},
    {0.924127593, 0.827384882, 0.774508472},
    {0.944468518, 0.800927443, 0.726736146},
    {0.958852946, 0.769767752, 0.678007945},
    {0.96732803, 0.734132809, 0.628751763},
    {0.969954137, 0.694266682, 0.579375448},
    {0.966811177, 0.650421156, 0.530263762},
    {0.958003065, 0.602842431, 0.481775914},
    {0.943660866, 0.551750968, 0.434243684},
    {0.923944917, 0.49730856, 0.387970225},
    {0.89904617, 0.439559467, 0.343229596},
    {0.869186849, 0.378313092, 0.300267182},
    {0.834620542, 0.312874446, 0.259301199},
    {0.795631745, 0.24128379, 0.220525627},
    {0.752534934, 0.157246067, 0.184115123},
    {0.705673158, 0.01555616, 0.150232812}
  };

  for (int i = 0; i < 127; i++) {
    double x = i*(32 - 1e-10)/(127 - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

void gray (double cmap[127][3])
{
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = i/(127 - 1.);
}

void randomap (double cmap[127][3])
{
  srand(0);
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (noise() + 1.)/2.;
}

void blue_white_red (double cmap[127][3])
{
  for (int i = 0; i < (127 + 1)/2; i++) {
    cmap[i][0] = i/((127 - 1)/2.);
    cmap[i][1] = i/((127 - 1)/2.);
    cmap[i][2] = 1.;
  }
  for (int i = 0; i < (127 - 1)/2; i++) {
    cmap[i + (127 + 1)/2][0] = 1.;
    cmap[i + (127 + 1)/2][1] = cmap[(127 - 3)/2 - i][1];
    cmap[i + (127 + 1)/2][2] = cmap[(127 - 3)/2 - i][1];
  }
}





typedef struct {
  unsigned char r, g, b;
} color;

color colormap_color (double cmap[127][3],
        double val, double min, double max)
{
  color c;
  if (val == HUGE) {
    c.r = c.g = c.b = 0;
    return c;
  }
  int i;
  double coef;
  if (max != min)
    val = (val - min)/(max - min);
  else
    val = 0.;
  if (val <= 0.) i = 0, coef = 0.;
  else if (val >= 1.) i = 127 - 2, coef = 1.;
  else {
    i = val*(127 - 1);
    coef = val*(127 - 1) - i;
  }
  if (!(i < 127 - 1)) qassert ("/home/gongminjiang/basilisk/src/output.h", 321, "i < NCMAP - 1");
  unsigned char * c1 = (unsigned char *) &c;
  for (int j = 0; j < 3; j++)
    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
  return c;
}
#line 340 "/home/gongminjiang/basilisk/src/output.h"
static const char * extension (const char * file, const char * ext) {
  int len = strlen(file);
  return len > 4 && !strcmp (file + len - 4, ext) ? file + len - 4 : NULL;
}

static const char * is_animation (const char * file) {
  const char * ext;
  if ((ext = extension (file, ".mp4")) ||
      (ext = extension (file, ".ogv")) ||
      (ext = extension (file, ".gif")))
    return ext;
  return NULL;
}

static struct {
  FILE ** fp;
  char ** names;
  int n;
} open_image_data = {NULL, NULL, 0};

static void open_image_cleanup()
{
  for (int i = 0; i < open_image_data.n; i++) {
    qpclose (open_image_data.fp[i]);
    pfree (open_image_data.names[i],__func__,__FILE__,__LINE__);
  }
  pfree (open_image_data.fp,__func__,__FILE__,__LINE__);
  pfree (open_image_data.names,__func__,__FILE__,__LINE__);
  open_image_data.fp = NULL;
  open_image_data.names = NULL;
  open_image_data.n = 0;
}

static FILE * open_image_lookup (const char * file)
{
  for (int i = 0; i < open_image_data.n; i++)
    if (!strcmp (file, open_image_data.names[i]))
      return open_image_data.fp[i];
  return NULL;
}

static bool which (const char * command)
{
  char * s = getenv ("PATH");
  if (!s)
    return false;
  char path[strlen(s) + 1];
  strcpy (path, s);
  s = strtok (path, ":");
  while (s) {
    char f[strlen(s) + strlen(command) + 2];
    strcpy (f, s);
    strcat (f, "/");
    strcat (f, command);
    FILE * fp = fopen (f, "r");
    if (fp) {
      fclose (fp);
      return true;
    }
    s = strtok (NULL, ":");
  }
  return false;
}

static FILE * ppm_fallback (const char * file, const char * mode)
{
  char filename[strlen(file) + 5];
  strcpy (filename, file);
  strcat (filename, ".ppm");
  FILE * fp = fopen (filename, mode);
  if (!fp) {
    perror (file);

    MPI_Abort (MPI_COMM_WORLD, 1);

    exit (1);
  }
  return fp;
}

FILE * open_image (const char * file, const char * options)
{
  if (!(pid() == 0)) qassert ("/home/gongminjiang/basilisk/src/output.h", 422, "pid() == 0");
  const char * ext;
  if ((ext = is_animation (file))) {
    FILE * fp = open_image_lookup (file);
    if (fp)
      return fp;

    int len = strlen ("ppm2???    ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "ppm2"); strcat (command, ext + 1);

    static int has_ffmpeg = -1;
    if (has_ffmpeg < 0) {
      if (which (command) && (which ("ffmpeg") || which ("avconv")))
 has_ffmpeg = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find '%s' or 'ffmpeg'/'avconv'\n"
   "  falling back to raw PPM outputs\n", command);
 has_ffmpeg = false;
      }
    }
    if (!has_ffmpeg)
      return ppm_fallback (file, "a");

    static bool added = false;
    if (!added) {
      free_solver_func_add (open_image_cleanup);
      added = true;
    }
    open_image_data.n++;
    open_image_data.names = (char * *) prealloc (open_image_data.names, (open_image_data.n)*sizeof(char *),__func__,__FILE__,__LINE__);
    open_image_data.names[open_image_data.n - 1] = pstrdup (file,__func__,__FILE__,__LINE__);

    if (options) {
      strcat (command, " ");
      strcat (command, options);
    }
    strcat (command, !strcmp (ext, ".mp4") ? " " : " > ");
    strcat (command, file);
    open_image_data.fp = (FILE * *) prealloc (open_image_data.fp, (open_image_data.n)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    return open_image_data.fp[open_image_data.n - 1] = qpopen (command, "w");
  }
  else {
    static int has_convert = -1;
    if (has_convert < 0) {
      if (which ("convert"))
 has_convert = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find 'convert'\n"
   "  falling back to raw PPM outputs\n");
 has_convert = false;
      }
    }
    if (!has_convert)
      return ppm_fallback (file, "w");

    int len = strlen ("convert ppm:-   ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "convert ppm:- ");
    if (options) {
      strcat (command, options);
      strcat (command, " ");
    }
    strcat (command, file);
    return qpopen (command, "w");
  }
}

void close_image (const char * file, FILE * fp)
{
  if (!(pid() == 0)) qassert ("/home/gongminjiang/basilisk/src/output.h", 496, "pid() == 0");
  if (is_animation (file)) {
    if (!open_image_lookup (file))
      fclose (fp);
  }
  else if (which ("convert"))
    qpclose (fp);
  else
    fclose (fp);
}
#line 571 "/home/gongminjiang/basilisk/src/output.h"
struct OutputPPM {
  scalar f;
  FILE * fp;
  int n;
  char * file;
  double min, max, spread, z;
  bool linear;
  double box[2][2];
  scalar mask;
  colormap map;
  char * opt;
};

     
void output_ppm (struct OutputPPM p)
{tracing("output_ppm","/home/gongminjiang/basilisk/src/output.h",585);

  if (p.n == 0) p.n = N;
  if (p.min == 0 && p.max == 0) {
    stats s = statsf (p.f);
    if (p.spread < 0.)
      p.min = s.min, p.max = s.max;
    else {
      double avg = s.sum/s.volume, spread = (p.spread ? p.spread : 5.)*s.stddev;
      p.min = avg - spread; p.max = avg + spread;
    }
  }
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }
  if (!p.map)
    p.map = jet;
  if (p.linear) {
    scalar f = p.f, mask = p.mask;
    if (mask.i)
      boundary_internal ((scalar *)((scalar[]){f, mask,{-1}}), "/home/gongminjiang/basilisk/src/output.h", 608);
    else
      boundary_internal ((scalar *)((scalar[]){f,{-1}}), "/home/gongminjiang/basilisk/src/output.h", 610);
  }

  double fn = p.n;
  double Delta = (p.box[1][0] - p.box[0][0])/fn;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;
  if (ny % 2) ny++;

  color ** ppm = (color **) matrix_new (ny, p.n, sizeof(color));
  double cmap[127][3];
  p.map (cmap);
  OMP_PARALLEL() {
    OMP(omp for schedule(static))
      for (int j = 0; j < ny; j++) {
 double yp = Delta*j + p.box[0][1] + Delta/2.;
 for (int i = 0; i < p.n; i++) {
   double xp = Delta*i + p.box[0][0] + Delta/2., v;
   if (p.mask.i) {
     if (p.linear) {
       double m = interpolate ((struct _interpolate){p.mask, xp, yp, p.z});
       if (m < 0.)
  v = HUGE;
       else
  v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
     }
     else {
       Point point = locate ((struct _locate){xp, yp, p.z});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
       if (point.level < 0 || val(p.mask,0,0,0) < 0.)
  v = HUGE;
       else
  v = val(p.f,0,0,0);
     }
   }
   else if (p.linear)
     v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
   else {
     Point point = locate ((struct _locate){xp, yp, p.z});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
     v = point.level >= 0 ? val(p.f,0,0,0) : HUGE;
   }
   ppm[ny - 1 - j][i] = colormap_color (cmap, v, p.min, p.max);
 }
      }
  }

  if (pid() == 0) {
#if 1
    MPI_Reduce (MPI_IN_PLACE, ppm[0], 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif
    if (!p.fp) p.fp = fout;
    if (p.file)
      p.fp = open_image (p.file, p.opt);

    fprintf (p.fp, "P6\n%u %u 255\n", p.n, ny);
    fwrite (((void **) ppm)[0], sizeof(color), ny*p.n, p.fp);

    if (p.file)
      close_image (p.file, p.fp);
    else
      fflush (p.fp);
  }
#if 1
  else
    MPI_Reduce (ppm[0], NULL, 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (ppm);
end_tracing("output_ppm","/home/gongminjiang/basilisk/src/output.h",678);}
#line 710 "/home/gongminjiang/basilisk/src/output.h"
struct OutputGRD {
  scalar f;
  FILE * fp;
  double Delta;
  bool linear;
  double box[2][2];
  scalar mask;
};

     
void output_grd (struct OutputGRD p)
{tracing("output_grd","/home/gongminjiang/basilisk/src/output.h",720);

  if (!p.fp) p.fp = fout;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
    if (p.Delta == 0) p.Delta = L0/N;
  }
  if (p.linear) {
    scalar f = p.f, mask = p.mask;
    if (mask.i)
      boundary_internal ((scalar *)((scalar[]){f, mask,{-1}}), "/home/gongminjiang/basilisk/src/output.h", 733);
    else
      boundary_internal ((scalar *)((scalar[]){f,{-1}}), "/home/gongminjiang/basilisk/src/output.h", 735);
  }

  double Delta = p.Delta;
  int nx = (p.box[1][0] - p.box[0][0])/Delta;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;


  fprintf (p.fp, "ncols          %d\n", nx);
  fprintf (p.fp, "nrows          %d\n", ny);
  fprintf (p.fp, "xllcorner      %g\n", p.box[0][0]);
  fprintf (p.fp, "yllcorner      %g\n", p.box[0][1]);
  fprintf (p.fp, "cellsize       %g\n", Delta);
  fprintf (p.fp, "nodata_value   -9999\n");


  for (int j = ny-1; j >= 0; j--) {
    double yp = Delta*j + p.box[0][1] + Delta/2.;
    for (int i = 0; i < nx; i++) {
      double xp = Delta*i + p.box[0][0] + Delta/2., v;
      if (p.mask.i) {
 if (p.linear) {
   double m = interpolate ((struct _interpolate){p.mask, xp, yp});
   if (m < 0.)
     v = HUGE;
   else
     v = interpolate ((struct _interpolate){p.f, xp, yp});
 }
 else {
   Point point = locate ((struct _locate){xp, yp});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
   if (point.level < 0 || val(p.mask,0,0,0) < 0.)
     v = HUGE;
   else
     v = val(p.f,0,0,0);
 }
      }
      else if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
 v = point.level >= 0 ? val(p.f,0,0,0) : HUGE;
      }
      if (v == HUGE)
 fprintf (p.fp, "-9999 ");
      else
 fprintf (p.fp, "%f ", v);
    }
    fprintf (p.fp, "\n");
  }

  fflush (p.fp);
end_tracing("output_grd","/home/gongminjiang/basilisk/src/output.h",786);}
#line 813 "/home/gongminjiang/basilisk/src/output.h"
struct OutputGfs {
  FILE * fp;
  scalar * list;
  double t;
  char * file;
  bool translate;
};

static char * replace (const char * input, int target, int with,
         bool translate)
{
  if (translate) {
    if (!strcmp (input, "u.x"))
      return pstrdup ("U",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.y"))
      return pstrdup ("V",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.z"))
      return pstrdup ("W",__func__,__FILE__,__LINE__);
  }
  char * name = pstrdup (input,__func__,__FILE__,__LINE__), * i = name;
  while (*i != '\0') {
    if (*i == target)
      *i = with;
    i++;
  }
  return name;
}

     
void output_gfs (struct OutputGfs p)
{tracing("output_gfs","/home/gongminjiang/basilisk/src/output.h",842);
  char * fname = p.file;

#if 1



  FILE * fp = p.fp;
  if (p.file == NULL) {
    long pid = getpid();
    MPI_Bcast (&pid, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    fname = ((char *) pmalloc ((80)*sizeof(char),__func__,__FILE__,__LINE__));
    snprintf (fname, 80, ".output-%ld", pid);
    p.fp = NULL;
  }
#endif

  bool opened = false;
  if (p.fp == NULL) {
    if (fname == NULL)
      p.fp = fout;
    else if (!(p.fp = fopen (fname, "w"))) {
      perror (fname);
      exit (1);
    }
    else
      opened = true;
  }

  scalar * list = p.list ? p.list : list_copy (all);

  restriction (list);
  fprintf (p.fp,
    "1 0 GfsSimulation GfsBox GfsGEdge { binary = 1"
    " x = %g y = %g ",
    0.5 + X0/L0, 0.5 + Y0/L0);




  if (list != NULL && list[0].i != -1) {
    scalar s = list[0];
    char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
    fprintf (p.fp, "variables = %s", name);
    pfree (name,__func__,__FILE__,__LINE__);
    for (int i = 1; i < list_len(list); i++) {
      scalar s = list[i];
      if (_attribute[s.i].name) {
 char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
 fprintf (p.fp, ",%s", name);
 pfree (name,__func__,__FILE__,__LINE__);
      }
    }
    fprintf (p.fp, " ");
  }
  fprintf (p.fp, "} {\n");
  fprintf (p.fp, "  Time { t = %g }\n", t);
  if (L0 != 1.)
    fprintf (p.fp, "  PhysicalParams { L = %g }\n", L0);
  fprintf (p.fp, "  VariableTracerVOF f\n");
  fprintf (p.fp, "}\nGfsBox { x = 0 y = 0 z = 0 } {\n");

#if 1
  long header;
  if ((header = ftell (p.fp)) < 0) {
    perror ("output_gfs(): error in header");
    exit (1);
  }
  int cell_size = sizeof(unsigned) + sizeof(double);
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (_attribute[s.i].name)
      cell_size += sizeof(double);}}
  scalar index = new_scalar("index");
  size_t total_size = header + (z_indexing (index, false) + 1)*cell_size;
#endif



  {foreach_cell() {
#if 1
    if (is_local(cell))
#endif
    {
#if 1
      if (fseek (p.fp, header + val(index,0,0,0)*cell_size, SEEK_SET) < 0) {
 perror ("output_gfs(): error while seeking");
 exit (1);
      }
#endif
      unsigned flags =
 level == 0 ? 0 :



      child.x == -1 && child.y == -1 ? 0 :
 child.x == -1 && child.y == 1 ? 1 :
 child.x == 1 && child.y == -1 ? 2 :
 3;
#line 951 "/home/gongminjiang/basilisk/src/output.h"
      if (is_leaf(cell))
 flags |= (1 << 4);
      fwrite (&flags, sizeof (unsigned), 1, p.fp);
      double a = -1;
      fwrite (&a, sizeof (double), 1, p.fp);
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 if (_attribute[s.i].name) {
   if (_attribute[s.i].v.x.i >= 0) {




     if (_attribute[s.i].v.x.i == s.i) {
       s = _attribute[s.i].v.y;
       a = is_local(cell) && val(s,0,0,0) != HUGE ? val(s,0,0,0) : (double) DBL_MAX;
     }
     else if (_attribute[s.i].v.y.i == s.i) {
       s = _attribute[s.i].v.x;
       a = is_local(cell) && val(s,0,0,0) != HUGE ? - val(s,0,0,0) : (double) DBL_MAX;
     }





   }
   else
     a = is_local(cell) && val(s,0,0,0) != HUGE ? val(s,0,0,0) : (double) DBL_MAX;
   fwrite (&a, sizeof (double), 1, p.fp);
 }}}
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

#if 1
  delete (((scalar[]){index,{-1}}));
  if (!pid() && fseek (p.fp, total_size, SEEK_SET) < 0) {
    perror ("output_gfs(): error while finishing");
    exit (1);
  }
  if (!pid())
#endif
    fputs ("}\n", p.fp);
  fflush (p.fp);

  if (!p.list)
    pfree (list,__func__,__FILE__,__LINE__);
  if (opened)
    fclose (p.fp);

#if 1
  if (p.file == NULL) {
    MPI_Barrier (MPI_COMM_WORLD);
    if (pid() == 0) {
      if (fp == NULL)
 fp = fout;
      p.fp = fopen (fname, "r");
      size_t l;
      unsigned char buffer[8192];
      while ((l = fread (buffer, 1, 8192, p.fp)) > 0)
 fwrite (buffer, 1, l, fp);
      fflush (fp);
      remove (fname);
    }
    pfree (fname,__func__,__FILE__,__LINE__);
  }
#endif
end_tracing("output_gfs","/home/gongminjiang/basilisk/src/output.h",1019);}
#line 1043 "/home/gongminjiang/basilisk/src/output.h"
struct Dump {
  char * file;
  scalar * list;
  FILE * fp;
  bool unbuffered;
};

struct DumpHeader {
  double t;
  long len;
  int i, depth, npe, version;
  coord n;
};

static const int dump_version =

  170901;

static scalar * dump_list (scalar * lista)
{
  scalar * list = is_constant(cm) ? NULL : list_concat (((scalar[]){cm,{-1}}), NULL);
  {scalar*_i=(scalar*)( lista);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!_attribute[s.i].face && !_attribute[s.i].nodump && s.i != cm.i)
      list = list_add (list, s);}}
  return list;
}

static void dump_header (FILE * fp, struct DumpHeader * header, scalar * list)
{
  if (fwrite (header, sizeof(struct DumpHeader), 1, fp) < 1) {
    perror ("dump(): error while writing header");
    exit (1);
  }
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    unsigned len = strlen(_attribute[s.i].name);
    if (fwrite (&len, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing len");
      exit (1);
    }
    if (fwrite (_attribute[s.i].name, sizeof(char), len, fp) < len) {
      perror ("dump(): error while writing s.name");
      exit (1);
    }
  }}}
  double o[4] = {X0,Y0,Z0,L0};
  if (fwrite (o, sizeof(double), 4, fp) < 4) {
    perror ("dump(): error while writing coordinates");
    exit (1);
  }
}

#if !1
     
void dump (struct Dump p)
{tracing("dump","/home/gongminjiang/basilisk/src/output.h",1096);
  FILE * fp = p.fp;
  char def[] = "dump", * file = p.file ? p.file : p.fp ? NULL : def;

  char * name = NULL;
  if (file) {
    name = (char *) pmalloc (strlen(file) + 2,__func__,__FILE__,__LINE__);
    strcpy (name, file);
    if (!p.unbuffered)
      strcat (name, "~");
    if ((fp = fopen (name, "w")) == NULL) {
      perror (name);
      exit (1);
    }
  }
  if (!(fp)) qassert ("/home/gongminjiang/basilisk/src/output.h", 1112, "fp");

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar  size=new_scalar("size");
  scalar * list = list_concat (((scalar[]){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
          dump_version };
  dump_header (fp, &header, list);

  subtree_size (size, false);

  {foreach_cell() {
    unsigned flags = is_leaf(cell) ? leaf : 0;
    if (fwrite (&flags, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing flags");
      exit (1);
    }
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (fwrite (&val(s,0,0,0), sizeof(double), 1, fp) < 1) {
 perror ("dump(): error while writing scalars");
 exit (1);
      }}}
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  pfree (list,__func__,__FILE__,__LINE__);
  if (file) {
    fclose (fp);
    if (!p.unbuffered)
      rename (name, file);
    pfree (name,__func__,__FILE__,__LINE__);
  }delete((scalar*)((scalar[]){size,{-1}}));
end_tracing("dump","/home/gongminjiang/basilisk/src/output.h",1145);}
#else
     
void dump (struct Dump p)
{tracing("dump","/home/gongminjiang/basilisk/src/output.h",1148);
  FILE * fp = p.fp;
  char def[] = "dump", * file = p.file ? p.file : p.fp ? NULL : def;

  if (fp != NULL || file == NULL) {
    fprintf (ferr, "dump(): must specify a file name when using MPI\n");
    exit(1);
  }

  char name[strlen(file) + 2];
  strcpy (name, file);
  if (!p.unbuffered)
    strcat (name, "~");
  FILE * fh = fopen (name, "w");
  if (fh == NULL) {
    perror (name);
    exit (1);
  }

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar  size=new_scalar("size");
  scalar * list = list_concat (((scalar[]){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
          dump_version };







  if (pid() == 0)
    dump_header (fh, &header, list);

  scalar index = {-1};

  index = new_scalar("index");
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
  int sizeofheader = sizeof(header) + 4*sizeof(double);
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    sizeofheader += sizeof(unsigned) + sizeof(char)*strlen(_attribute[s.i].name);}}
  long pos = pid() ? 0 : sizeofheader;

  subtree_size (size, false);

  {foreach_cell() {

    if (is_local(cell)) {
      long offset = sizeofheader + val(index,0,0,0)*cell_size;
      if (pos != offset) {
 fseek (fh, offset, SEEK_SET);
 pos = offset;
      }
      unsigned flags = is_leaf(cell) ? leaf : 0;
      fwrite (&flags, 1, sizeof(unsigned), fh);
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 fwrite (&val(s,0,0,0), 1, sizeof(double), fh);}}
      pos += cell_size;
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}

  delete (((scalar[]){index,{-1}}));

  pfree (list,__func__,__FILE__,__LINE__);
  fclose (fh);
  if (!p.unbuffered && pid() == 0)
    rename (name, file);delete((scalar*)((scalar[]){size,{-1}}));
end_tracing("dump","/home/gongminjiang/basilisk/src/output.h",1219);}
#endif

     
bool restore (struct Dump p)
{tracing("restore","/home/gongminjiang/basilisk/src/output.h",1223);
  FILE * fp = p.fp;
  char * file = p.file;
  if (file && (fp = fopen (file, "r")) == NULL)
    {end_tracing("restore","/home/gongminjiang/basilisk/src/output.h",1228);return false;}
  if (!(fp)) qassert ("/home/gongminjiang/basilisk/src/output.h", 1229, "fp");

  struct DumpHeader header;
  if (fread (&header, sizeof(header), 1, fp) < 1) {
    fprintf (ferr, "restore(): error: expecting header\n");
    exit (1);
  }


  init_grid (1);
  {foreach_cell() {
    cell.pid = pid();
    cell.flags |= active;
  }end_foreach_cell();}
  ((Tree *)grid)->dirty = true;
#line 1264 "/home/gongminjiang/basilisk/src/output.h"
  bool restore_all = (p.list == all);
  scalar * list = dump_list (p.list ? p.list : all);
  if (header.version == 161020) {
    if (header.len - 1 != list_len (list)) {
      fprintf (ferr,
        "restore(): error: the list lengths don't match: "
        "%ld (file) != %d (code)\n",
        header.len - 1, list_len (list));
      exit (1);
    }
  }
  else {
    if (header.version != dump_version) {
      fprintf (ferr,
        "restore(): error: file version mismatch: "
        "%d (file) != %d (code)\n",
        header.version, dump_version);
      exit (1);
    }

    scalar * input = NULL;
    for (int i = 0; i < header.len; i++) {
      unsigned len;
      if (fread (&len, sizeof(unsigned), 1, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting len\n");
 exit (1);
      }
      char name[len + 1];
      if (fread (name, sizeof(char), len, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting s.name\n");
 exit (1);
      }
      name[len] = '\0';

      if (i > 0) {
 bool found = false;
 {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   if (!strcmp (_attribute[s.i].name, name)) {
     input = list_append (input, s);
     found = true; break;
   }}}
 if (!found) {
   if (restore_all) {
     scalar s = new_scalar("s");
     pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
     _attribute[s.i].name = pstrdup (name,__func__,__FILE__,__LINE__);
     input = list_append (input, s);
   }
   else
     input = list_append (input, (scalar){INT_MAX});
 }
      }
    }
    pfree (list,__func__,__FILE__,__LINE__);
    list = input;

    double o[4];
    if (fread (o, sizeof(double), 4, fp) < 4) {
      fprintf (ferr, "restore(): error: expecting coordinates\n");
      exit (1);
    }
    origin ((struct _origin){o[0], o[1], o[2]});
    size (o[3]);
  }
#line 1339 "/home/gongminjiang/basilisk/src/output.h"
  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector[]){fm,{{-1},{-1}}});

  restore_mpi (fp, list);
#line 1369 "/home/gongminjiang/basilisk/src/output.h"
  scalar * other = NULL;
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!list_lookup (list, s) && !list_lookup (listm, s))
      other = list_append (other, s);}}
  reset (other, 0.);
  pfree (other,__func__,__FILE__,__LINE__);

  pfree (list,__func__,__FILE__,__LINE__);
  if (file)
    fclose (fp);


  while (iter < header.i && events (false))
    iter = inext;
  events (false);
  while (t < header.t && events (false))
    t = tnext;
  t = header.t;
  events (false);

  {end_tracing("restore","/home/gongminjiang/basilisk/src/output.h",1389);return true;}
end_tracing("restore","/home/gongminjiang/basilisk/src/output.h",1390);}
#line 431 "/home/gongminjiang/basilisk/src/utils.h"
#line 12 "/home/gongminjiang/basilisk/src/run.h"

     
void run (void)
{tracing("run","/home/gongminjiang/basilisk/src/run.h",14);
  iter = 0, t = 0., dt = 1.;
  init_grid (N);

  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  while (events (true)) {





    update_perf();
    iter = inext, t = tnext;
  }




  timer_print (perf.gt, iter, perf.tnc);

  free_grid();
end_tracing("run","/home/gongminjiang/basilisk/src/run.h",37);}




static int defaults_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0);*ip=i;*tp=t;return ret;}      static int defaults(const int i,const double t,Event *_ev){tracing("defaults","/home/gongminjiang/basilisk/src/run.h",42); {
  display ((struct _display){"box();"});
}{end_tracing("defaults","/home/gongminjiang/basilisk/src/run.h",44);return 0;}end_tracing("defaults","/home/gongminjiang/basilisk/src/run.h",44);}





static int cleanup_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t = 1234567890);*ip=i;*tp=t;return ret;}      static int cleanup(const int i,const double t,Event *_ev){tracing("cleanup","/home/gongminjiang/basilisk/src/run.h",50); {
  display ((struct _display){"", true});
}{end_tracing("cleanup","/home/gongminjiang/basilisk/src/run.h",52);return 0;}end_tracing("cleanup","/home/gongminjiang/basilisk/src/run.h",52);}
#line 28 "/home/gongminjiang/basilisk/src/navier-stokes/centered.h"
#line 1 "./timestep.h"
#line 1 "/home/gongminjiang/basilisk/src/timestep.h"

double timestep (const vector u, double dtmax)
{
  static double previous = 0.;
  dtmax /= CFL;
  foreach_face_stencil(){_stencil_is_face_x(){
    {_stencil_val(u.x,0,0,0); {   
      _stencil_val(u.x,0,0,0); 




_stencil_val(cm,0,0,0);   




       

         
    
#line 16
}   }}end__stencil_is_face_x()
#line 6
_stencil_is_face_y(){
    {_stencil_val(u.y,0,0,0); {   
      _stencil_val(u.y,0,0,0); 




_stencil_val(cm,0,0,0);   




       

         
    
#line 16
}   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#line 6
if(!is_constant(cm)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(min:dtmax)){
#line 6
foreach_face_generic(){is_face_x(){
    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.x,0,0,0));




      dt *= val(cm,0,0,0);

      if (dt < dtmax) dtmax = dt;
    }}end_is_face_x()
#line 6
is_face_y(){
    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.y,0,0,0));




      dt *= val(cm,0,0,0);

      if (dt < dtmax) dtmax = dt;
    }}end_is_face_y()}end_foreach_face_generic();mpi_all_reduce_array(&dtmax,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 16
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(min:dtmax)){
#line 6
foreach_face_generic(){is_face_x(){
    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.x,0,0,0));




      dt *= _const_cm;

      if (dt < dtmax) dtmax = dt;
    }}end_is_face_x()
#line 6
is_face_y(){
    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.y,0,0,0));




      dt *= _const_cm;

      if (dt < dtmax) dtmax = dt;
    }}end_is_face_y()}end_foreach_face_generic();mpi_all_reduce_array(&dtmax,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 16
}
  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}
#line 29 "/home/gongminjiang/basilisk/src/navier-stokes/centered.h"
#line 1 "./bcg.h"
#line 1 "/home/gongminjiang/basilisk/src/bcg.h"
#line 11 "/home/gongminjiang/basilisk/src/bcg.h"
void tracer_fluxes (scalar f,
      vector uf,
      vector flux,
      double dt,
              scalar src)
{





  vector  g=new_vector("g");
  gradients (((scalar[]){f,{-1}}), ((vector[]){g,{{-1},{-1}}}));




  foreach_face_stencil(){_stencil_is_face_x(){ {        







    _stencil_val(fm.x,0,0,0);_stencil_val(uf.x,0,0,0);              
    
    _stencil_val(g.x,o_stencil,0,0); _stencil_val(src,-1,0,0);_stencil_val(src,0,0,0);_stencil_val(f, o_stencil,0,0);





_stencil_val(fm.y,o_stencil,0,0); _stencil_val(fm.y,o_stencil,1,0); {     
       _stencil_val(fm.y,o_stencil,1,0);_stencil_val(fm.y,o_stencil,0,0); _stencil_val(uf.y,o_stencil,1,0);_stencil_val(uf.y,o_stencil,0,0);         
       _stencil_val(f,o_stencil,-1,0);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,1,0);
        
    }





      
#line 58 "/home/gongminjiang/basilisk/src/bcg.h"
    _stencil_val_a(flux.x,0,0,0);_stencil_val(uf.x,0,0,0);  
  }}end__stencil_is_face_x()
#line 28
_stencil_is_face_y(){ {        







    _stencil_val(fm.y,0,0,0);_stencil_val(uf.y,0,0,0);              
    
    _stencil_val(g.y,0,o_stencil,0); _stencil_val(src,0,-1,0);_stencil_val(src,0,0,0);_stencil_val(f,0, o_stencil,0);





_stencil_val(fm.x,0,o_stencil,0); _stencil_val(fm.x,1,o_stencil,0); {     
       _stencil_val(fm.x,1,o_stencil,0);_stencil_val(fm.x,0,o_stencil,0); _stencil_val(uf.x,1,o_stencil,0);_stencil_val(uf.x,0,o_stencil,0);         
       _stencil_val(f,-1,o_stencil,0);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,1,o_stencil,0);
        
    }





      
#line 58 "/home/gongminjiang/basilisk/src/bcg.h"
    _stencil_val_a(flux.y,0,0,0);_stencil_val(uf.y,0,0,0);  
  }}end__stencil_is_face_y()}end_foreach_face_stencil();




  
#line 28
if(!is_constant(fm.x) && !is_constant(src)){{foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(val(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val(src,0,0,0) + val(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val(fm.y,i,0,0) + val(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/gongminjiang/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(val(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val(src,0,0,0) + val(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val(fm.x,0,i,0) + val(fm.x,1,i,0));
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/gongminjiang/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(src)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);




  {
#line 28
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(_const_fm.x*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val(src,0,0,0) + val(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (_const_fm.y && _const_fm.y) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(_const_fm.y + _const_fm.y);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/gongminjiang/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(_const_fm.y*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val(src,0,0,0) + val(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (_const_fm.x && _const_fm.x) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(_const_fm.x + _const_fm.x);
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/gongminjiang/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(src)){double _const_src=_constant[src.i-_NVARMAX];NOT_UNUSED(_const_src);




  {
#line 28
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(val(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val(fm.y,i,0,0) + val(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/gongminjiang/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(val(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val(fm.x,0,i,0) + val(fm.x,1,i,0));
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/gongminjiang/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_src=_constant[src.i-_NVARMAX];NOT_UNUSED(_const_src);




  {
#line 28
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(_const_fm.x*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (_const_fm.y && _const_fm.y) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(_const_fm.y + _const_fm.y);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/gongminjiang/basilisk/src/bcg.h"
    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 28
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(_const_fm.y*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (_const_fm.x && _const_fm.x) {
      double vn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(_const_fm.x + _const_fm.x);
      double fyy = vn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }
#line 58 "/home/gongminjiang/basilisk/src/bcg.h"
    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}delete((scalar*)((vector[]){g,{{-1},{-1}}}));
}






struct Advection {
  scalar * tracers;
  vector u;
  double dt;
  scalar * src;
};

void advection (struct Advection p)
{




  scalar * lsrc = p.src;
  if (!lsrc)
    {scalar*_i=(scalar*)( p.tracers);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      lsrc = list_append (lsrc, zeroc);}}
  if (!(list_len(p.tracers) == list_len(lsrc))) qassert ("/home/gongminjiang/basilisk/src/bcg.h", 84, "list_len(p.tracers) == list_len(lsrc)");

  scalar f, src;
  {scalar*_i0=lsrc;scalar*_i1= p.tracers;if(_i0)for(src=*_i0,f=*_i1;_i0->i>= 0;src=*++_i0,f=*++_i1){ {
    vector  flux=new_face_vector("flux");
    tracer_fluxes (f, p.u, flux, p.dt, src);

    foreach_stencil()
      {
        {_stencil_val_r(f,0,0,0);_stencil_val(flux.x,0,0,0); _stencil_val(flux.x,1,0,0);_stencil_val(cm,0,0,0);   }
        
#line 93
{_stencil_val_r(f,0,0,0);_stencil_val(flux.y,0,0,0); _stencil_val(flux.y,0,1,0);_stencil_val(cm,0,0,0);   }}end_foreach_stencil();

    
#line 91
if(!is_constant(cm)){{foreach()
      {
        val(f,0,0,0) += p.dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val(cm,0,0,0));
        
#line 93
val(f,0,0,0) += p.dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val(cm,0,0,0));}end_foreach();}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

    {
#line 91
foreach()
      {
        val(f,0,0,0) += p.dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*_const_cm);
        
#line 93
val(f,0,0,0) += p.dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*_const_cm);}end_foreach();}}delete((scalar*)((vector[]){flux,{{-1},{-1}}}));



  }}}

  if (!p.src)
    pfree (lsrc,__func__,__FILE__,__LINE__);
}
#line 30 "/home/gongminjiang/basilisk/src/navier-stokes/centered.h"



#line 1 "./viscosity.h"
#line 1 "/home/gongminjiang/basilisk/src/viscosity.h"
#line 1 "./poisson.h"
#line 1 "/home/gongminjiang/basilisk/src/poisson.h"
#line 32 "/home/gongminjiang/basilisk/src/poisson.h"
void mg_cycle (scalar * a, scalar * res, scalar * da,
        void (* relax) (scalar * da, scalar * res,
          int depth, void * data),
        void * data,
        int nrelax, int minlevel, int maxlevel)
{




  restriction (res);





  minlevel = min (minlevel, maxlevel);
  for (int l = minlevel; l <= maxlevel; l++) {




    if (l == minlevel)
      {foreach_level_or_leaf (l)
 {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     val(s,0,0,0) = 0.;}}end_foreach_level_or_leaf();}





    else
      {foreach_level (l)
 {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     val(s,0,0,0) = bilinear (point, s);}}end_foreach_level();}





    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
    }
  }




  foreach_stencil() {
    scalar s, ds;
    {scalar*_i0= da;scalar*_i1= a;if(_i0)for(ds=*_i0,s=*_i1;_i0->i>= 0;ds=*++_i0,s=*++_i1){
     
 {_stencil_val_r(s,0,0,0); _stencil_val(ds,0,0,0); }}}
  }end_foreach_stencil();




  {
#line 84
foreach() {
    scalar s, ds;
    {scalar*_i0= da;scalar*_i1= a;if(_i0)for(ds=*_i0,s=*_i1;_i0->i>= 0;ds=*++_i0,s=*++_i1){
     
 val(s,0,0,0) += val(ds,0,0,0);}}
  }end_foreach();}
}
#line 102 "/home/gongminjiang/basilisk/src/poisson.h"
int NITERMAX = 100, NITERMIN = 1;
double TOLERANCE = 1e-3;




typedef struct {
  int i;
  double resb, resa;
  double sum;
  int nrelax;
  int minlevel;
} mgstats;
#line 125 "/home/gongminjiang/basilisk/src/poisson.h"
struct MGSolve {
  scalar * a, * b;
  double (* residual) (scalar * a, scalar * b, scalar * res,
         void * data);
  void (* relax) (scalar * da, scalar * res, int depth,
    void * data);
  void * data;

  int nrelax;
  scalar * res;
  int minlevel;
  double tolerance;
};

mgstats mg_solve (struct MGSolve p)
{





  scalar * da = list_clone (p.a), * res = p.res;
  if (!res)
    res = list_clone (p.b);






  for (int b = 0; b < nboundary; b++)
    {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b];}}




  mgstats s = {0};
  double sum = 0.;
  foreach_stencil ()
    {scalar*_i=(scalar*)( p.b);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      { _stencil_val(s,0,0,0); }}}end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:sum)){
#line 164
foreach ()
    {scalar*_i=(scalar*)( p.b);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      sum += val(s,0,0,0);}}end_foreach();mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
  
#line 167
s.sum = sum;
  s.nrelax = p.nrelax > 0 ? p.nrelax : 4;




  double resb;
  resb = s.resb = s.resa = p.residual (p.a, p.b, res, p.data);






  if (p.tolerance == 0.)
    p.tolerance = TOLERANCE;
  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > p.tolerance);
       s.i++) {
    mg_cycle (p.a, res, da, p.relax, p.data,
       s.nrelax,
       p.minlevel,
       grid->maxdepth);
    s.resa = p.residual (p.a, p.b, res, p.data);
#line 199 "/home/gongminjiang/basilisk/src/poisson.h"
    if (s.resa > p.tolerance) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
 s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
 s.nrelax--;
    }







    resb = s.resa;
  }
  s.minlevel = p.minlevel;




  if (s.resa > p.tolerance) {
    scalar v = p.a[0];
    fprintf (ferr,
      "WARNING: convergence for %s not reached after %d iterations\n"
      "  res: %g sum: %g nrelax: %d\n", _attribute[v.i].name,
      s.i, s.resa, s.sum, s.nrelax), fflush (ferr);
  }




  if (!p.res)
    delete (res), pfree (res,__func__,__FILE__,__LINE__);
  delete (da), pfree (da,__func__,__FILE__,__LINE__);

  return s;
}
#line 258 "/home/gongminjiang/basilisk/src/poisson.h"
struct Poisson {
  scalar a, b;
          vector alpha;
          scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;



};





static void relax (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = (struct Poisson *) data;
          vector alpha = p->alpha;
          scalar lambda = p->lambda;
#line 296 "/home/gongminjiang/basilisk/src/poisson.h"
  scalar c = a;






  if(!is_constant(lambda) && !is_constant(alpha.x)){{foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - val(lambda,0,0,0)*sq(Delta);
     {
      n += val(alpha.x,1,0,0)*val(a,1,0,0) + val(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val(alpha.x,1,0,0) + val(alpha.x,0,0,0);
    } 
#line 305
{
      n += val(alpha.y,0,1,0)*val(a,0,1,0) + val(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val(alpha.y,0,1,0) + val(alpha.y,0,0,0);
    }
#line 319 "/home/gongminjiang/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}else if(is_constant(lambda) && !is_constant(alpha.x)){double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);






  {
#line 303
foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - _const_lambda*sq(Delta);
     {
      n += val(alpha.x,1,0,0)*val(a,1,0,0) + val(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val(alpha.x,1,0,0) + val(alpha.x,0,0,0);
    } 
#line 305
{
      n += val(alpha.y,0,1,0)*val(a,0,1,0) + val(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val(alpha.y,0,1,0) + val(alpha.y,0,0,0);
    }
#line 319 "/home/gongminjiang/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}else if(!is_constant(lambda) && is_constant(alpha.x)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);






  {
#line 303
foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - val(lambda,0,0,0)*sq(Delta);
     {
      n += _const_alpha.x*val(a,1,0,0) + _const_alpha.x*val(a,-1,0,0);
      d += _const_alpha.x + _const_alpha.x;
    } 
#line 305
{
      n += _const_alpha.y*val(a,0,1,0) + _const_alpha.y*val(a,0,-1,0);
      d += _const_alpha.y + _const_alpha.y;
    }
#line 319 "/home/gongminjiang/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}else {double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);






  {
#line 303
foreach_level_or_leaf (l) {
    double n = - sq(Delta)*val(b,0,0,0), d = - _const_lambda*sq(Delta);
     {
      n += _const_alpha.x*val(a,1,0,0) + _const_alpha.x*val(a,-1,0,0);
      d += _const_alpha.x + _const_alpha.x;
    } 
#line 305
{
      n += _const_alpha.y*val(a,0,1,0) + _const_alpha.y*val(a,0,-1,0);
      d += _const_alpha.y + _const_alpha.y;
    }
#line 319 "/home/gongminjiang/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}}
#line 338 "/home/gongminjiang/basilisk/src/poisson.h"
}






static double residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = (struct Poisson *) data;
          vector alpha = p->alpha;
          scalar lambda = p->lambda;
  double maxres = 0.;


  vector  g=new_face_vector("g");
  foreach_face_stencil(){_stencil_is_face_x(){
    {_stencil_val_a(g.x,0,0,0); _stencil_val(alpha.x,0,0,0);_stencil_val(a,0,0,0); _stencil_val(a,0 -1,0,0);  }}end__stencil_is_face_x()
#line 355
_stencil_is_face_y(){
    {_stencil_val_a(g.y,0,0,0); _stencil_val(alpha.y,0,0,0);_stencil_val(a,0,0,0); _stencil_val(a,0,0 -1,0);  }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#line 355
if(!is_constant(alpha.x)){{foreach_face_generic(){is_face_x(){
    val(g.x,0,0,0) = val(alpha.x,0,0,0)*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta);}end_is_face_x()
#line 355
is_face_y(){
    val(g.y,0,0,0) = val(alpha.y,0,0,0)*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta);}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  {
#line 355
foreach_face_generic(){is_face_x(){
    val(g.x,0,0,0) = _const_alpha.x*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta);}end_is_face_x()
#line 355
is_face_y(){
    val(g.y,0,0,0) = _const_alpha.y*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta);}end_is_face_y()}end_foreach_face_generic();}}
  foreach_stencil () {
    _stencil_val_a(res,0,0,0); _stencil_val(b,0,0,0); _stencil_val(lambda,0,0,0);_stencil_val(a,0,0,0);  
    
      {_stencil_val_r(res,0,0,0);_stencil_val(g.x,1,0,0); _stencil_val(g.x,0,0,0);   }
      
#line 360
{_stencil_val_r(res,0,0,0);_stencil_val(g.y,0,1,0); _stencil_val(g.y,0,0,0);   }






_stencil_val(res,0,0,0);
      {_stencil_val(res,0,0,0);   }






        
  
#line 369
}end_foreach_stencil();
  
#line 357
if(!is_constant(lambda)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 357
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - val(lambda,0,0,0)*val(a,0,0,0);
    
      val(res,0,0,0) -= (val(g.x,1,0,0) - val(g.x,0,0,0))/Delta;
      
#line 360
val(res,0,0,0) -= (val(g.y,0,1,0) - val(g.y,0,0,0))/Delta;






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 369
}else {double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 357
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - _const_lambda*val(a,0,0,0);
    
      val(res,0,0,0) -= (val(g.x,1,0,0) - val(g.x,0,0,0))/Delta;
      
#line 360
val(res,0,0,0) -= (val(g.y,0,1,0) - val(g.y,0,0,0))/Delta;






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 369
}
#line 387 "/home/gongminjiang/basilisk/src/poisson.h"
  {delete((scalar*)((vector[]){g,{{-1},{-1}}}));return maxres;}delete((scalar*)((vector[]){g,{{-1},{-1}}}));
}
#line 399 "/home/gongminjiang/basilisk/src/poisson.h"
mgstats poisson (struct Poisson p)
{






  if (!p.alpha.x.i)
    p.alpha = unityf;
  if (!p.lambda.i)
    p.lambda = zeroc;




  vector alpha = p.alpha;
  scalar lambda = p.lambda;
  restriction (((scalar[]){alpha.x,alpha.y,lambda,{-1}}));





  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  scalar a = p.a, b = p.b;




  mgstats s = mg_solve ((struct MGSolve){((scalar[]){a,{-1}}), ((scalar[]){b,{-1}}), residual, relax,
   &p, p.nrelax, p.res, .minlevel = max(1, p.minlevel)});




  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}
#line 461 "/home/gongminjiang/basilisk/src/poisson.h"
struct Project {
  vector uf;
  scalar p;
  vector alpha;
  double dt;
  int nrelax;
};

     
mgstats project (struct Project q)
{tracing("project","/home/gongminjiang/basilisk/src/poisson.h",470);
  vector uf = q.uf;
  scalar p = q.p;
          vector alpha = q.alpha.x.i ? q.alpha : unityf;
  double dt = q.dt ? q.dt : 1.;
  int nrelax = q.nrelax ? q.nrelax : 4;






  scalar  div=new_scalar("div");
  foreach_stencil() {
    _stencil_val_a(div,0,0,0);  
    
      {_stencil_val_r(div,0,0,0); _stencil_val(uf.x,1,0,0); _stencil_val(uf.x,0,0,0);  }
      
#line 487
{_stencil_val_r(div,0,0,0); _stencil_val(uf.y,0,1,0); _stencil_val(uf.y,0,0,0);  }
    _stencil_val_r(div,0,0,0);  
  }end_foreach_stencil();
  {
#line 484
foreach() {
    val(div,0,0,0) = 0.;
    
      val(div,0,0,0) += val(uf.x,1,0,0) - val(uf.x,0,0,0);
      
#line 487
val(div,0,0,0) += val(uf.y,0,1,0) - val(uf.y,0,0,0);
    val(div,0,0,0) /= dt*Delta;
  }end_foreach();}
#line 500 "/home/gongminjiang/basilisk/src/poisson.h"
  mgstats mgp = poisson ((struct Poisson){p, div, alpha,
    .tolerance = TOLERANCE/sq(dt), .nrelax = nrelax});




  foreach_face_stencil(){_stencil_is_face_x(){
    {_stencil_val_r(uf.x,0,0,0);_stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0 -1,0,0);   }}end__stencil_is_face_x()
#line 506
_stencil_is_face_y(){
    {_stencil_val_r(uf.y,0,0,0);_stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,0 -1,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();




  
#line 506
if(!is_constant(alpha.x)){{foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) -= dt*val(alpha.x,0,0,0)*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()
#line 506
is_face_y(){
    val(uf.y,0,0,0) -= dt*val(alpha.y,0,0,0)*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);




  {
#line 506
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) -= dt*_const_alpha.x*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()
#line 506
is_face_y(){
    val(uf.y,0,0,0) -= dt*_const_alpha.y*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()}end_foreach_face_generic();}}

  {delete((scalar*)((scalar[]){div,{-1}}));{end_tracing("project","/home/gongminjiang/basilisk/src/poisson.h",509);return mgp;}}delete((scalar*)((scalar[]){div,{-1}}));
end_tracing("project","/home/gongminjiang/basilisk/src/poisson.h",510);}
#line 2 "/home/gongminjiang/basilisk/src/viscosity.h"

struct Viscosity {
  vector u;
  vector mu;
  scalar rho;
  double dt;
  int nrelax;
  scalar * res;
};
#line 25 "/home/gongminjiang/basilisk/src/viscosity.h"
static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
          vector mu = p->mu;
          scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0])));




  vector w = u;


  if(!is_constant(rho) && !is_constant(mu.x)){{foreach_level_or_leaf (l) {
    
      val(w.x,0,0,0) = (dt/val(rho,0,0,0)*(2.*val(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)
#line 58 "/home/gongminjiang/basilisk/src/viscosity.h"
      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).x + dt/val(rho,0,0,0)*(2.*val(mu.x,1,0,0) + 2.*val(mu.x,0,0,0)

          + val(mu.y,0,1,0) + val(mu.y,0,0,0)




        ));
      
#line 41
val(w.y,0,0,0) = (dt/val(rho,0,0,0)*(2.*val(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val(mu.x,1,0,0)*(val(u.y,1,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)
#line 58 "/home/gongminjiang/basilisk/src/viscosity.h"
      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).y + dt/val(rho,0,0,0)*(2.*val(mu.y,0,1,0) + 2.*val(mu.y,0,0,0)

          + val(mu.x,1,0,0) + val(mu.x,0,0,0)




        ));
  }end_foreach_level_or_leaf();}}else if(is_constant(rho) && !is_constant(mu.x)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);


  {
#line 39
foreach_level_or_leaf (l) {
    
      val(w.x,0,0,0) = (dt/_const_rho*(2.*val(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)
#line 58 "/home/gongminjiang/basilisk/src/viscosity.h"
      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).x + dt/_const_rho*(2.*val(mu.x,1,0,0) + 2.*val(mu.x,0,0,0)

          + val(mu.y,0,1,0) + val(mu.y,0,0,0)




        ));
      
#line 41
val(w.y,0,0,0) = (dt/_const_rho*(2.*val(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val(mu.x,1,0,0)*(val(u.y,1,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)
#line 58 "/home/gongminjiang/basilisk/src/viscosity.h"
      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).y + dt/_const_rho*(2.*val(mu.y,0,1,0) + 2.*val(mu.y,0,0,0)

          + val(mu.x,1,0,0) + val(mu.x,0,0,0)




        ));
  }end_foreach_level_or_leaf();}}else if(!is_constant(rho) && is_constant(mu.x)){struct{double x,y;}_const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX]};NOT_UNUSED(_const_mu);


  {
#line 39
foreach_level_or_leaf (l) {
    
      val(w.x,0,0,0) = (dt/val(rho,0,0,0)*(2.*_const_mu.x*val(u.x,1,0,0) + 2.*_const_mu.x*val(u.x,-1,0,0)

      + _const_mu.y*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - _const_mu.y*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)
#line 58 "/home/gongminjiang/basilisk/src/viscosity.h"
      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).x + dt/val(rho,0,0,0)*(2.*_const_mu.x + 2.*_const_mu.x

          + _const_mu.y + _const_mu.y




        ));
      
#line 41
val(w.y,0,0,0) = (dt/val(rho,0,0,0)*(2.*_const_mu.y*val(u.y,0,1,0) + 2.*_const_mu.y*val(u.y,0,-1,0)

      + _const_mu.x*(val(u.y,1,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - _const_mu.x*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)
#line 58 "/home/gongminjiang/basilisk/src/viscosity.h"
      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).y + dt/val(rho,0,0,0)*(2.*_const_mu.y + 2.*_const_mu.y

          + _const_mu.x + _const_mu.x




        ));
  }end_foreach_level_or_leaf();}}else {double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);struct{double x,y;}_const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX]};NOT_UNUSED(_const_mu);


  {
#line 39
foreach_level_or_leaf (l) {
    
      val(w.x,0,0,0) = (dt/_const_rho*(2.*_const_mu.x*val(u.x,1,0,0) + 2.*_const_mu.x*val(u.x,-1,0,0)

      + _const_mu.y*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - _const_mu.y*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)
#line 58 "/home/gongminjiang/basilisk/src/viscosity.h"
      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).x + dt/_const_rho*(2.*_const_mu.x + 2.*_const_mu.x

          + _const_mu.y + _const_mu.y




        ));
      
#line 41
val(w.y,0,0,0) = (dt/_const_rho*(2.*_const_mu.y*val(u.y,0,1,0) + 2.*_const_mu.y*val(u.y,0,-1,0)

      + _const_mu.x*(val(u.y,1,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - _const_mu.x*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)
#line 58 "/home/gongminjiang/basilisk/src/viscosity.h"
      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.}).y + dt/_const_rho*(2.*_const_mu.y + 2.*_const_mu.y

          + _const_mu.x + _const_mu.x




        ));
  }end_foreach_level_or_leaf();}}
#line 85 "/home/gongminjiang/basilisk/src/viscosity.h"
}

static double residual_viscosity (scalar * a, scalar * b, scalar * resl,
      void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
          vector mu = p->mu;
          scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0]))), res = (*((vector *)&(resl[0])));
  double maxres = 0.;
#line 104 "/home/gongminjiang/basilisk/src/viscosity.h"
  boundary_internal ((scalar *)((vector[]){u,{{-1},{-1}}}), "/home/gongminjiang/basilisk/src/viscosity.h", 104);

   {
    vector  taux=new_face_vector("taux");
    foreach_face_stencil()_stencil_is_face_x(){
      {_stencil_val_a(taux.x,0,0,0);_stencil_val(mu.x,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,-1,0,0);   }}end__stencil_is_face_x()end_foreach_face_stencil();
    
#line 108
if(!is_constant(mu.x)){{foreach_face_generic()is_face_x(){
      val(taux.x,0,0,0) = 2.*val(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta;}end_is_face_x()end_foreach_face_generic();}}else {struct{double x,y;}_const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX]};NOT_UNUSED(_const_mu);
    {
#line 108
foreach_face_generic()is_face_x(){
      val(taux.x,0,0,0) = 2.*_const_mu.x*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta;}end_is_face_x()end_foreach_face_generic();}}

      foreach_face_stencil()_stencil_is_face_y(){
 {_stencil_val_a(taux.y,0,0,0); _stencil_val(mu.y,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,0,-1,0);
_stencil_val(u.y,1,-1,0); _stencil_val(u.y,1,0,0);
_stencil_val(u.y,-1,-1,0); _stencil_val(u.y,-1,0,0);    
        
      
#line 114
}}end__stencil_is_face_y()end_foreach_face_stencil();

      
#line 111
if(!is_constant(mu.x)){{foreach_face_generic()is_face_y(){
 val(taux.y,0,0,0) = val(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta;}end_is_face_y()end_foreach_face_generic();}}else {struct{double x,y;}_const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX]};NOT_UNUSED(_const_mu);

      {
#line 111
foreach_face_generic()is_face_y(){
 val(taux.y,0,0,0) = _const_mu.y*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta;}end_is_face_y()end_foreach_face_generic();}}







    foreach_stencil () {   
      
      
 { _stencil_val(taux.x,1,0,0); _stencil_val(taux.x,0,0,0);  }
 
#line 125
{ _stencil_val(taux.y,0,1,0); _stencil_val(taux.y,0,0,0);  }
      _stencil_val_a(res.x,0,0,0); _stencil_val(r.x,0,0,0);_stencil_val(u.x,0,0,0);_stencil_val(rho,0,0,0);
_stencil_val(res.x,0,0,0);
 {_stencil_val(res.x,0,0,0);   }     
          
    
#line 129
}end_foreach_stencil();







    
#line 122
if(!is_constant(rho)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 122
foreach () {
      double d = 0.;
      
 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
 
#line 125
d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.}).x*val(u.x,0,0,0) + dt/val(rho,0,0,0)*d/Delta;
      if (fabs (val(res.x,0,0,0)) > maxres)
 maxres = fabs (val(res.x,0,0,0));
    }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 129
}else {double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);







    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 122
foreach () {
      double d = 0.;
      
 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
 
#line 125
d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.}).x*val(u.x,0,0,0) + dt/_const_rho*d/Delta;
      if (fabs (val(res.x,0,0,0)) > maxres)
 maxres = fabs (val(res.x,0,0,0));
    }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 129
}delete((scalar*)((vector[]){taux,{{-1},{-1}}}));
  } 
#line 106
{
    vector  taux=new_face_vector("taux");
    foreach_face_stencil()_stencil_is_face_y(){
      {_stencil_val_a(taux.y,0,0,0);_stencil_val(mu.y,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,-1,0);   }}end__stencil_is_face_y()end_foreach_face_stencil();
    
#line 108
if(!is_constant(mu.y)){{foreach_face_generic()is_face_y(){
      val(taux.y,0,0,0) = 2.*val(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta;}end_is_face_y()end_foreach_face_generic();}}else {struct{double x,y;}_const_mu={_constant[mu.y.i-_NVARMAX],_constant[mu.x.i-_NVARMAX]};NOT_UNUSED(_const_mu);
    {
#line 108
foreach_face_generic()is_face_y(){
      val(taux.y,0,0,0) = 2.*_const_mu.y*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta;}end_is_face_y()end_foreach_face_generic();}}

      foreach_face_stencil()_stencil_is_face_x(){
 {_stencil_val_a(taux.x,0,0,0); _stencil_val(mu.x,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,-1,0,0);
_stencil_val(u.x,-1,1,0); _stencil_val(u.x,0,1,0);
_stencil_val(u.x,-1,-1,0); _stencil_val(u.x,0,-1,0);    
        
      
#line 114
}}end__stencil_is_face_x()end_foreach_face_stencil();

      
#line 111
if(!is_constant(mu.y)){{foreach_face_generic()is_face_x(){
 val(taux.x,0,0,0) = val(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta;}end_is_face_x()end_foreach_face_generic();}}else {struct{double x,y;}_const_mu={_constant[mu.y.i-_NVARMAX],_constant[mu.x.i-_NVARMAX]};NOT_UNUSED(_const_mu);

      {
#line 111
foreach_face_generic()is_face_x(){
 val(taux.x,0,0,0) = _const_mu.x*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta;}end_is_face_x()end_foreach_face_generic();}}







    foreach_stencil () {   
      
      
 { _stencil_val(taux.y,0,1,0); _stencil_val(taux.y,0,0,0);  }
 
#line 125
{ _stencil_val(taux.x,1,0,0); _stencil_val(taux.x,0,0,0);  }
      _stencil_val_a(res.y,0,0,0); _stencil_val(r.y,0,0,0);_stencil_val(u.y,0,0,0);_stencil_val(rho,0,0,0);
_stencil_val(res.y,0,0,0);
 {_stencil_val(res.y,0,0,0);   }     
          
    
#line 129
}end_foreach_stencil();







    
#line 122
if(!is_constant(rho)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 122
foreach () {
      double d = 0.;
      
 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
 
#line 125
d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.}).y*val(u.y,0,0,0) + dt/val(rho,0,0,0)*d/Delta;
      if (fabs (val(res.y,0,0,0)) > maxres)
 maxres = fabs (val(res.y,0,0,0));
    }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 129
}else {double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);







    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)){
#line 122
foreach () {
      double d = 0.;
      
 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
 
#line 125
d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.}).y*val(u.y,0,0,0) + dt/_const_rho*d/Delta;
      if (fabs (val(res.y,0,0,0)) > maxres)
 maxres = fabs (val(res.y,0,0,0));
    }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 129
}delete((scalar*)((vector[]){taux,{{-1},{-1}}}));
  }
#line 159 "/home/gongminjiang/basilisk/src/viscosity.h"
  return maxres;
}



     
mgstats viscosity (struct Viscosity p)
{tracing("viscosity","/home/gongminjiang/basilisk/src/viscosity.h",165);
  vector u = p.u,  r=new_vector("r");
  foreach_stencil()
    {
      {_stencil_val_a(r.x,0,0,0); _stencil_val(u.x,0,0,0); }
      
#line 170
{_stencil_val_a(r.y,0,0,0); _stencil_val(u.y,0,0,0); }}end_foreach_stencil();
  {
#line 168
foreach()
    {
      val(r.x,0,0,0) = val(u.x,0,0,0);
      
#line 170
val(r.y,0,0,0) = val(u.y,0,0,0);}end_foreach();}

  vector mu = p.mu;
  scalar rho = p.rho;
  restriction (((scalar[]){mu.x,mu.y,rho,{-1}}));

  { mgstats _ret= mg_solve ((struct MGSolve){(scalar *)((vector[]){u,{{-1},{-1}}}), (scalar *)((vector[]){r,{{-1},{-1}}}),
     residual_viscosity, relax_viscosity, &p, p.nrelax, p.res});delete((scalar*)((vector[]){r,{{-1},{-1}}}));{end_tracing("viscosity","/home/gongminjiang/basilisk/src/viscosity.h",177);return _ret;}}delete((scalar*)((vector[]){r,{{-1},{-1}}}));
end_tracing("viscosity","/home/gongminjiang/basilisk/src/viscosity.h",178);}

     
mgstats viscosity_explicit (struct Viscosity p)
{tracing("viscosity_explicit","/home/gongminjiang/basilisk/src/viscosity.h",181);
  vector u = p.u,  r=new_vector("r");
  mgstats mg = {0};
  mg.resb = residual_viscosity ((scalar *)((vector[]){u,{{-1},{-1}}}), (scalar *)((vector[]){u,{{-1},{-1}}}), (scalar *)((vector[]){r,{{-1},{-1}}}), &p);
  foreach_stencil()
    {
      {_stencil_val_r(u.x,0,0,0); _stencil_val(r.x,0,0,0); }
      
#line 188
{_stencil_val_r(u.y,0,0,0); _stencil_val(r.y,0,0,0); }}end_foreach_stencil();
  {
#line 186
foreach()
    {
      val(u.x,0,0,0) += val(r.x,0,0,0);
      
#line 188
val(u.y,0,0,0) += val(r.y,0,0,0);}end_foreach();}
  {delete((scalar*)((vector[]){r,{{-1},{-1}}}));{end_tracing("viscosity_explicit","/home/gongminjiang/basilisk/src/viscosity.h",189);return mg;}}delete((scalar*)((vector[]){r,{{-1},{-1}}}));
end_tracing("viscosity_explicit","/home/gongminjiang/basilisk/src/viscosity.h",190);}
#line 34 "/home/gongminjiang/basilisk/src/navier-stokes/centered.h"
#line 44 "/home/gongminjiang/basilisk/src/navier-stokes/centered.h"
scalar  p={0};
vector  u={{1},{2}},  g={{3},{4}};
scalar  pf={5};
vector  uf={{6},{7}};
#line 70 "/home/gongminjiang/basilisk/src/navier-stokes/centered.h"
        vector mu = {{_NVARMAX+0},{_NVARMAX+1}}, a = {{_NVARMAX+0},{_NVARMAX+1}}, alpha = {{_NVARMAX+2},{_NVARMAX+3}};
        scalar rho = {_NVARMAX+4};
mgstats mgp, mgpf, mgu;
bool stokes = false;
#line 91
static double _boundary0(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann ((val(a.x,1,0,0)*val(fm.x,1,0,0)/val(alpha.x,1,0,0))));}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann ((_const_a.x*val(fm.x,1,0,0)/val(alpha.x,1,0,0))));}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann ((val(a.x,1,0,0)*_const_fm.x/val(alpha.x,1,0,0))));}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann ((_const_a.x*_const_fm.x/val(alpha.x,1,0,0))));}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann ((val(a.x,1,0,0)*val(fm.x,1,0,0)/_const_alpha.x)));}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann ((_const_a.x*val(fm.x,1,0,0)/_const_alpha.x)));}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann ((val(a.x,1,0,0)*_const_fm.x/_const_alpha.x)));}}}else {struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann ((_const_a.x*_const_fm.x/_const_alpha.x)));}}}}static double _boundary0_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann_homogeneous ());}}
static double _boundary1(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann (- (val(a.x,0,0,0)*val(fm.x,0,0,0)/val(alpha.x,0,0,0))));}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann (- (_const_a.x*val(fm.x,0,0,0)/val(alpha.x,0,0,0))));}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann (- (val(a.x,0,0,0)*_const_fm.x/val(alpha.x,0,0,0))));}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann (- (_const_a.x*_const_fm.x/val(alpha.x,0,0,0))));}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann (- (val(a.x,0,0,0)*val(fm.x,0,0,0)/_const_alpha.x)));}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann (- (_const_a.x*val(fm.x,0,0,0)/_const_alpha.x)));}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann (- (val(a.x,0,0,0)*_const_fm.x/_const_alpha.x)));}}}else {struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann (- (_const_a.x*_const_fm.x/_const_alpha.x)));}}}}static double _boundary1_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann_homogeneous ());}}








static double _boundary2(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann ((val(a.y,0,1,0)*val(fm.y,0,1,0)/val(alpha.y,0,1,0))));}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann ((_const_a.y*val(fm.y,0,1,0)/val(alpha.y,0,1,0))));}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann ((val(a.y,0,1,0)*_const_fm.y/val(alpha.y,0,1,0))));}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann ((_const_a.y*_const_fm.y/val(alpha.y,0,1,0))));}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann ((val(a.y,0,1,0)*val(fm.y,0,1,0)/_const_alpha.y)));}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann ((_const_a.y*val(fm.y,0,1,0)/_const_alpha.y)));}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann ((val(a.y,0,1,0)*_const_fm.y/_const_alpha.y)));}}}else {struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann ((_const_a.y*_const_fm.y/_const_alpha.y)));}}}}static double _boundary2_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann_homogeneous ());}}
static double _boundary3(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann (- (val(a.y,0,0,0)*val(fm.y,0,0,0)/val(alpha.y,0,0,0))));}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann (- (_const_a.y*val(fm.y,0,0,0)/val(alpha.y,0,0,0))));}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann (- (val(a.y,0,0,0)*_const_fm.y/val(alpha.y,0,0,0))));}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann (- (_const_a.y*_const_fm.y/val(alpha.y,0,0,0))));}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann (- (val(a.y,0,0,0)*val(fm.y,0,0,0)/_const_alpha.y)));}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann (- (_const_a.y*val(fm.y,0,0,0)/_const_alpha.y)));}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann (- (val(a.y,0,0,0)*_const_fm.y/_const_alpha.y)));}}}else {struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann (- (_const_a.y*_const_fm.y/_const_alpha.y)));}}}}static double _boundary3_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann_homogeneous ());}}
#line 126 "/home/gongminjiang/basilisk/src/navier-stokes/centered.h"
static int defaults_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0);*ip=i;*tp=t;return ret;}      static int defaults_0(const int i,const double t,Event *_ev){tracing("defaults_0","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",126);
{

  CFL = 0.8;




  _attribute[p.i].nodump = _attribute[pf.i].nodump = true;




  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {
    vector alphav = alpha;
    foreach_face_stencil(){_stencil_is_face_x(){
      {_stencil_val_a(alphav.x,0,0,0); _stencil_val(fm.x,0,0,0); }}end__stencil_is_face_x()
#line 145
_stencil_is_face_y(){
      {_stencil_val_a(alphav.y,0,0,0); _stencil_val(fm.y,0,0,0); }}end__stencil_is_face_y()}end_foreach_face_stencil();
    
#line 145
if(!is_constant(fm.x)){{foreach_face_generic(){is_face_x(){
      val(alphav.x,0,0,0) = val(fm.x,0,0,0);}end_is_face_x()
#line 145
is_face_y(){
      val(alphav.y,0,0,0) = val(fm.y,0,0,0);}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
    {
#line 145
foreach_face_generic(){is_face_x(){
      val(alphav.x,0,0,0) = _const_fm.x;}end_is_face_x()
#line 145
is_face_y(){
      val(alphav.y,0,0,0) = _const_fm.y;}end_is_face_y()}end_foreach_face_generic();}}
  }






  _attribute[uf.x.i].refine = refine_face_solenoidal;
#line 173 "/home/gongminjiang/basilisk/src/navier-stokes/centered.h"
}{end_tracing("defaults_0","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",173);return 0;}end_tracing("defaults_0","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",173);}





static int default_display_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0);*ip=i;*tp=t;return ret;}      static int default_display(const int i,const double t,Event *_ev){tracing("default_display","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",179);
  display ((struct _display){"squares (color = 'u.x', spread = -1);"});{end_tracing("default_display","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",180);return 0;}end_tracing("default_display","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",180);}





double dtmax;

static int init_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0);*ip=i;*tp=t;return ret;}      static int init(const int i,const double t,Event *_ev){tracing("init","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",188);
{
  trash (((vector[]){uf,{{-1},{-1}}}));
  foreach_face_stencil(){_stencil_is_face_x(){
    {_stencil_val_a(uf.x,0,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);  }}end__stencil_is_face_x()
#line 191
_stencil_is_face_y(){
    {_stencil_val_a(uf.y,0,0,0); _stencil_val(fm.y,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);  }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#line 191
if(!is_constant(fm.x)){{foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = val(fm.x,0,0,0)*((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.);}end_is_face_x()
#line 191
is_face_y(){
    val(uf.y,0,0,0) = val(fm.y,0,0,0)*((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.);}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  {
#line 191
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = _const_fm.x*((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.);}end_is_face_x()
#line 191
is_face_y(){
    val(uf.y,0,0,0) = _const_fm.y*((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.);}end_is_face_y()}end_foreach_face_generic();}}




  event ("properties");





  dtmax = DT;
  event ("stability");
}{end_tracing("init","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",205);return 0;}end_tracing("init","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",205);}
#line 214 "/home/gongminjiang/basilisk/src/navier-stokes/centered.h"
static int set_dtmax_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int set_dtmax(const int i,const double t,Event *_ev){tracing("set_dtmax","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",214); dtmax = DT;{end_tracing("set_dtmax","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",214);return 0;}end_tracing("set_dtmax","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",214);}

static int stability_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int stability(const int i,const double t,Event *_ev){tracing("stability","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",216); {
  dt = dtnext (stokes ? dtmax : timestep (uf, dtmax));
}{end_tracing("stability","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",218);return 0;}end_tracing("stability","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",218);}







static int vof_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}static int vof(const int i,const double t,Event *_ev){;return 0;}
static int tracer_advection_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}static int tracer_advection(const int i,const double t,Event *_ev){;return 0;}
static int tracer_diffusion_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}static int tracer_diffusion(const int i,const double t,Event *_ev){;return 0;}






static int properties_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}static int properties(const int i,const double t,Event *_ev){;return 0;}
#line 247 "/home/gongminjiang/basilisk/src/navier-stokes/centered.h"
void prediction()
{
  vector du;
   {
    scalar s = new_scalar("s");
    du.x = s;
  } 
#line 250
{
    scalar s = new_scalar("s");
    du.y = s;
  }

  if (_attribute[u.x.i].gradient)
    {
    
#line 256
foreach_stencil()
      { {





   _stencil_val_a(du.x,0,0,0);_stencil_val(u.x,-1,0,0); _stencil_val(u.x,0,0,0); _stencil_val(u.x,1,0,0);   
      } 
#line 257
{





   _stencil_val_a(du.y,0,0,0);_stencil_val(u.y,0,-1,0); _stencil_val(u.y,0,0,0); _stencil_val(u.y,0,1,0);   
      }}end_foreach_stencil();{
#line 256
foreach()
      { {





   val(du.x,0,0,0) = _attribute[u.x.i].gradient (val(u.x,-1,0,0), val(u.x,0,0,0), val(u.x,1,0,0))/Delta;
      } 
#line 257
{





   val(du.y,0,0,0) = _attribute[u.y.i].gradient (val(u.y,0,-1,0), val(u.y,0,0,0), val(u.y,0,1,0))/Delta;
      }}end_foreach();}}
  else
    {
    
#line 266
foreach_stencil()
      { {





   _stencil_val_a(du.x,0,0,0);_stencil_val(u.x,1,0,0); _stencil_val(u.x,-1,0,0);   
    } 
#line 267
{





   _stencil_val_a(du.y,0,0,0);_stencil_val(u.y,0,1,0); _stencil_val(u.y,0,-1,0);   
    }}end_foreach_stencil();{
#line 266
foreach()
      { {





   val(du.x,0,0,0) = (val(u.x,1,0,0) - val(u.x,-1,0,0))/(2.*Delta);
    } 
#line 267
{





   val(du.y,0,0,0) = (val(u.y,0,1,0) - val(u.y,0,-1,0))/(2.*Delta);
    }}end_foreach();}}

  trash (((vector[]){uf,{{-1},{-1}}}));
  foreach_face_stencil(){_stencil_is_face_x(){ {       
     _stencil_val(u.x,-1,0,0);_stencil_val(u.x,0,0,0);     
    
    _stencil_val_a(uf.x,0,0,0);_stencil_val(u.x, o_stencil,0,0);_stencil_val(g.x,0,0,0); _stencil_val(g.x,-1,0,0);_stencil_val(du.x,o_stencil,0,0);

_stencil_val(fm.y,o_stencil,0,0); _stencil_val(fm.y,o_stencil,1,0); {        
       _stencil_val(u.x,o_stencil,-1,0);_stencil_val(u.x, o_stencil,0,0);_stencil_val(u.x, o_stencil,0,0); _stencil_val(u.x,o_stencil,1,0);_stencil_val(u.y, o_stencil,0,0);
      _stencil_val_r(uf.x,0,0,0);_stencil_val(u.y,o_stencil,0,0);  
    }        

      







    
#line 293
_stencil_val_r(uf.x,0,0,0); _stencil_val(fm.x,0,0,0); 
  }}end__stencil_is_face_x()
#line 277
_stencil_is_face_y(){ {       
     _stencil_val(u.y,0,-1,0);_stencil_val(u.y,0,0,0);     
    
    _stencil_val_a(uf.y,0,0,0);_stencil_val(u.y,0, o_stencil,0);_stencil_val(g.y,0,0,0); _stencil_val(g.y,0,-1,0);_stencil_val(du.y,0,o_stencil,0);

_stencil_val(fm.x,0,o_stencil,0); _stencil_val(fm.x,1,o_stencil,0); {        
       _stencil_val(u.y,-1,o_stencil,0);_stencil_val(u.y,0, o_stencil,0);_stencil_val(u.y,0, o_stencil,0); _stencil_val(u.y,1,o_stencil,0);_stencil_val(u.x,0, o_stencil,0);
      _stencil_val_r(uf.y,0,0,0);_stencil_val(u.x,0,o_stencil,0);  
    }        

      







    
#line 293
_stencil_val_r(uf.y,0,0,0); _stencil_val(fm.y,0,0,0); 
  }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#line 277
if(!is_constant(fm.x)){{foreach_face_generic(){is_face_x(){ {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {
      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);
    }







    val(uf.x,0,0,0) *= val(fm.x,0,0,0);
  }}end_is_face_x()
#line 277
is_face_y(){ {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {
      double fyy = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fyy/(2.*Delta);
    }







    val(uf.y,0,0,0) *= val(fm.y,0,0,0);
  }}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  {
#line 277
foreach_face_generic(){is_face_x(){ {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

    if (_const_fm.y && _const_fm.y) {
      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);
    }







    val(uf.x,0,0,0) *= _const_fm.x;
  }}end_is_face_x()
#line 277
is_face_y(){ {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

    if (_const_fm.x && _const_fm.x) {
      double fyy = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fyy/(2.*Delta);
    }







    val(uf.y,0,0,0) *= _const_fm.y;
  }}end_is_face_y()}end_foreach_face_generic();}}

  delete ((scalar *)((vector[]){du,{{-1},{-1}}}));
}
#line 308 "/home/gongminjiang/basilisk/src/navier-stokes/centered.h"
static int advection_term_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int advection_term(const int i,const double t,Event *_ev){tracing("advection_term","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",308);
{
  if (!stokes) {
    prediction();
    mgpf = project ((struct Project){uf, pf, alpha, dt/2., mgpf.nrelax});
    advection ((struct Advection){(scalar *)((vector[]){u,{{-1},{-1}}}), uf, dt, (scalar *)((vector[]){g,{{-1},{-1}}})});
  }
}{end_tracing("advection_term","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",315);return 0;}end_tracing("advection_term","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",315);}







static void correction (double dt)
{
  foreach_stencil()
    {
      {_stencil_val_r(u.x,0,0,0);_stencil_val(g.x,0,0,0);  }
      
#line 327
{_stencil_val_r(u.y,0,0,0);_stencil_val(g.y,0,0,0);  }}end_foreach_stencil();
  {
#line 325
foreach()
    {
      val(u.x,0,0,0) += dt*val(g.x,0,0,0);
      
#line 327
val(u.y,0,0,0) += dt*val(g.y,0,0,0);}end_foreach();}
}
#line 337 "/home/gongminjiang/basilisk/src/navier-stokes/centered.h"
static int viscous_term_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int viscous_term(const int i,const double t,Event *_ev){tracing("viscous_term","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",337);
{
  if (constant(mu.x) != 0.) {
    correction (dt);
    mgu = viscosity ((struct Viscosity){u, mu, rho, dt, mgu.nrelax});
    correction (-dt);
  }




  if (!is_constant(a.x)) {
    vector af = a;
    trash (((vector[]){af,{{-1},{-1}}}));
    foreach_face_stencil(){_stencil_is_face_x(){
      {_stencil_val_a(af.x,0,0,0);  }}end__stencil_is_face_x()
#line 351
_stencil_is_face_y(){
      {_stencil_val_a(af.y,0,0,0);  }}end__stencil_is_face_y()}end_foreach_face_stencil();
    {
#line 351
foreach_face_generic(){is_face_x(){
      val(af.x,0,0,0) = 0.;}end_is_face_x()
#line 351
is_face_y(){
      val(af.y,0,0,0) = 0.;}end_is_face_y()}end_foreach_face_generic();}
  }
}{end_tracing("viscous_term","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",354);return 0;}end_tracing("viscous_term","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",354);}
#line 373 "/home/gongminjiang/basilisk/src/navier-stokes/centered.h"
static int acceleration_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int acceleration(const int i,const double t,Event *_ev){tracing("acceleration","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",373);
{
  trash (((vector[]){uf,{{-1},{-1}}}));
  foreach_face_stencil(){_stencil_is_face_x(){
    {_stencil_val_a(uf.x,0,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);_stencil_val(a.x,0,0,0);    }}end__stencil_is_face_x()
#line 376
_stencil_is_face_y(){
    {_stencil_val_a(uf.y,0,0,0); _stencil_val(fm.y,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);_stencil_val(a.y,0,0,0);    }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#line 376
if(!is_constant(fm.x) && !is_constant(a.x)){{foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = val(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val(a.x,0,0,0));}end_is_face_x()
#line 376
is_face_y(){
    val(uf.y,0,0,0) = val(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val(a.y,0,0,0));}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(a.x)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  {
#line 376
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = _const_fm.x*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val(a.x,0,0,0));}end_is_face_x()
#line 376
is_face_y(){
    val(uf.y,0,0,0) = _const_fm.y*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val(a.y,0,0,0));}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(a.x)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);
  {
#line 376
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = val(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*_const_a.x);}end_is_face_x()
#line 376
is_face_y(){
    val(uf.y,0,0,0) = val(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*_const_a.y);}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);
  {
#line 376
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = _const_fm.x*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*_const_a.x);}end_is_face_x()
#line 376
is_face_y(){
    val(uf.y,0,0,0) = _const_fm.y*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*_const_a.y);}end_is_face_y()}end_foreach_face_generic();}}
}{end_tracing("acceleration","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",378);return 0;}end_tracing("acceleration","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",378);}
#line 387 "/home/gongminjiang/basilisk/src/navier-stokes/centered.h"
void centered_gradient (scalar p, vector g)
{





  vector  gf=new_face_vector("gf");
  foreach_face_stencil(){_stencil_is_face_x(){
    {_stencil_val_a(gf.x,0,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(a.x,0,0,0); _stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);   }}end__stencil_is_face_x()
#line 395
_stencil_is_face_y(){
    {_stencil_val_a(gf.y,0,0,0); _stencil_val(fm.y,0,0,0);_stencil_val(a.y,0,0,0); _stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#line 395
if(!is_constant(fm.x) && !is_constant(a.x) && !is_constant(alpha.x)){{foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*val(a.x,0,0,0) - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 395
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*val(a.y,0,0,0) - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(a.x) && !is_constant(alpha.x)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  {
#line 395
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*val(a.x,0,0,0) - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 395
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*val(a.y,0,0,0) - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(a.x) && !is_constant(alpha.x)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);
  {
#line 395
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*_const_a.x - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 395
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*_const_a.y - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && is_constant(a.x) && !is_constant(alpha.x)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);
  {
#line 395
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*_const_a.x - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 395
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*_const_a.y - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && !is_constant(a.x) && is_constant(alpha.x)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  {
#line 395
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*val(a.x,0,0,0) - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 395
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*val(a.y,0,0,0) - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(a.x) && is_constant(alpha.x)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  {
#line 395
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*val(a.x,0,0,0) - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 395
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*val(a.y,0,0,0) - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(a.x) && is_constant(alpha.x)){struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  {
#line 395
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*_const_a.x - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 395
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*_const_a.y - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  {
#line 395
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*_const_a.x - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 395
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*_const_a.y - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()}end_foreach_face_generic();}}





  trash (((vector[]){g,{{-1},{-1}}}));
  foreach_stencil()
    {
      {_stencil_val_a(g.x,0,0,0);_stencil_val(gf.x,0,0,0); _stencil_val(gf.x,1,0,0);_stencil_val(fm.x,0,0,0); _stencil_val(fm.x,1,0,0);      }
      
#line 405
{_stencil_val_a(g.y,0,0,0);_stencil_val(gf.y,0,0,0); _stencil_val(gf.y,0,1,0);_stencil_val(fm.y,0,0,0); _stencil_val(fm.y,0,1,0);      }}end_foreach_stencil();
  
#line 403
if(!is_constant(fm.x)){{foreach()
    {
      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(val(fm.x,0,0,0) + val(fm.x,1,0,0) + 0.);
      
#line 405
val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(val(fm.y,0,0,0) + val(fm.y,0,1,0) + 0.);}end_foreach();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  {
#line 403
foreach()
    {
      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(_const_fm.x + _const_fm.x + 0.);
      
#line 405
val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(_const_fm.y + _const_fm.y + 0.);}end_foreach();}}delete((scalar*)((vector[]){gf,{{-1},{-1}}}));
}






static int projection_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int projection(const int i,const double t,Event *_ev){tracing("projection","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",413);
{
  mgp = project ((struct Project){uf, p, alpha, dt, mgp.nrelax});
  centered_gradient (p, g);




  correction (dt);
}{end_tracing("projection","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",422);return 0;}end_tracing("projection","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",422);}





static int end_timestep_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}static int end_timestep(const int i,const double t,Event *_ev){;return 0;}
#line 438 "/home/gongminjiang/basilisk/src/navier-stokes/centered.h"
static int adapt_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int adapt(const int i,const double t,Event *_ev){tracing("adapt","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",438); {






  event ("properties");
}{end_tracing("adapt","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",446);return 0;}end_tracing("adapt","/home/gongminjiang/basilisk/src/navier-stokes/centered.h",446);}
#line 15 "atomisation.c"
#line 1 "two-phase.h"
#line 1 "/home/gongminjiang/basilisk/src/two-phase.h"
#line 13 "/home/gongminjiang/basilisk/src/two-phase.h"
#line 1 "vof.h"
#line 1 "/home/gongminjiang/basilisk/src/vof.h"
#line 27 "/home/gongminjiang/basilisk/src/vof.h"





#line 1 "fractions.h"
#line 1 "/home/gongminjiang/basilisk/src/fractions.h"
#line 12 "/home/gongminjiang/basilisk/src/fractions.h"
#line 1 "geometry.h"
#line 1 "/home/gongminjiang/basilisk/src/geometry.h"
#line 28 "/home/gongminjiang/basilisk/src/geometry.h"
double line_alpha (double c, coord n)
{
  double alpha, n1, n2;

  n1 = fabs (n.x); n2 = fabs (n.y);
  if (n1 > n2)
    do { double __tmp = n1; n1 = n2; n2 = __tmp; } while(0);

  c = clamp (c, 0., 1.);
  double v1 = n1/2.;
  if (c <= v1/n2)
    alpha = sqrt (2.*c*n1*n2);
  else if (c <= 1. - v1/n2)
    alpha = c*n2 + v1;
  else
    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;

  return alpha - (n.x + n.y)/2.;
}
#line 133 "/home/gongminjiang/basilisk/src/geometry.h"
double line_area (double nx, double ny, double alpha)
{
  double a, v, area;

  alpha += (nx + ny)/2.;
  if (nx < 0.) {
    alpha -= nx;
    nx = - nx;
  }
  if (ny < 0.) {
    alpha -= ny;
    ny = - ny;
  }

  if (alpha <= 0.)
    return 0.;

  if (alpha >= nx + ny)
    return 1.;

  if (nx < 1e-10)
    area = alpha/ny;
  else if (ny < 1e-10)
    area = alpha/nx;
  else {
    v = sq(alpha);

    a = alpha - nx;
    if (a > 0.)
      v -= a*a;

    a = alpha - ny;
    if (a > 0.)
      v -= a*a;

    area = v/(2.*nx*ny);
  }

  return clamp (area, 0., 1.);
}
#line 237 "/home/gongminjiang/basilisk/src/geometry.h"
double rectangle_fraction (coord n, double alpha, coord a, coord b)
{
  coord n1;
   {
    alpha -= n.x*(b.x + a.x)/2.;
    n1.x = n.x*(b.x - a.x);
  } 
#line 240
{
    alpha -= n.y*(b.y + a.y)/2.;
    n1.y = n.y*(b.y - a.y);
  }
  return line_area(n1.x, n1.y, alpha);
}
#line 262 "/home/gongminjiang/basilisk/src/geometry.h"
int facets (coord n, double alpha, coord p[2])
{
  int i = 0;
  for (double s = -0.5; s <= 0.5; s += 1.)
    {
      if (fabs (n.y) > 1e-4 && i < 2) {
 double a = (alpha - s*n.x)/n.y;
 if (a >= -0.5 && a <= 0.5) {
   p[i].x = s;
   p[i++].y = a;
 }
      }
      
#line 267
if (fabs (n.x) > 1e-4 && i < 2) {
 double a = (alpha - s*n.y)/n.x;
 if (a >= -0.5 && a <= 0.5) {
   p[i].y = s;
   p[i++].x = a;
 }
      }}
  return i;
}
#line 352 "/home/gongminjiang/basilisk/src/geometry.h"
double line_length_center (coord m, double alpha, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
    
#line 358
if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }

  p->x = p->y = p->z = 0.;

  if (alpha <= 0. || alpha >= n.x + n.y)
    return 0.;

  
    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }
    
#line 369
if (n.y < 1e-4) {
      p->y = 0.;
      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }

  if (alpha >= n.x) {
    p->x += 1.;
    p->y += (alpha - n.x)/n.y;
  }
  else
    p->x += alpha/n.x;

  double ax = p->x, ay = p->y;
  if (alpha >= n.y) {
    p->y += 1.;
    ay -= 1.;
    p->x += (alpha - n.y)/n.x;
    ax -= (alpha - n.y)/n.x;
  }
  else {
    p->y += alpha/n.y;
    ay -= alpha/n.y;
  }

   {
    p->x /= 2.;
    p->x = clamp (p->x, 0., 1.);
    if (m.x < 0.)
      p->x = 1. - p->x;
    p->x -= 0.5;
  } 
#line 394
{
    p->y /= 2.;
    p->y = clamp (p->y, 0., 1.);
    if (m.y < 0.)
      p->y = 1. - p->y;
    p->y -= 0.5;
  }

  return sqrt (ax*ax + ay*ay);
}
#line 482 "/home/gongminjiang/basilisk/src/geometry.h"
void line_center (coord m, double alpha, double a, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
    
#line 488
if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }

  p->z = 0.;
  if (alpha <= 0.) {
    p->x = p->y = -0.5;
    return;
  }

  if (alpha >= n.x + n.y) {
    p->x = p->y = 0.;
    return;
  }

  
    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = sign(m.y)*(a/2. - 0.5);
      return;
    }
    
#line 505
if (n.y < 1e-4) {
      p->y = 0.;
      p->x = sign(m.x)*(a/2. - 0.5);
      return;
    }

  p->x = p->y = cube(alpha);

   {
    double b = alpha - n.x;
    if (b > 0.) {
      p->x -= sq(b)*(alpha + 2.*n.x);
      p->y -= cube(b);
    }
  } 
#line 513
{
    double b = alpha - n.y;
    if (b > 0.) {
      p->y -= sq(b)*(alpha + 2.*n.y);
      p->x -= cube(b);
    }
  }

   {
    p->x /= 6.*sq(n.x)*n.y*a;
    p->x = sign(m.x)*(p->x - 0.5);
  } 
#line 521
{
    p->y /= 6.*sq(n.y)*n.x*a;
    p->y = sign(m.y)*(p->y - 0.5);
  }
}
#line 13 "/home/gongminjiang/basilisk/src/fractions.h"






#line 1 "myc2d.h"
#line 1 "/home/gongminjiang/basilisk/src/myc2d.h"





coord mycs (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  int ix;
  double c_t,c_b,c_r,c_l;
  double mx0,my0,mx1,my1,mm1,mm2;


  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);
  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);
  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);
  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);



  mx0 = 0.5*(c_l-c_r);
  my0 = 0.5*(c_b-c_t);


  if (fabs(mx0) <= fabs(my0)) {
    my0 = my0 > 0. ? 1. : -1.;
    ix = 1;
  }
  else {
    mx0 = mx0 > 0. ? 1. : -1.;
    ix = 0;
  }


  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);
  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);
  mx1 = mm1 - mm2 + 1.e-30;
  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);
  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);
  my1 = mm1 - mm2 + 1.e-30;


  if (ix) {
    mm1 = fabs(my1);
    mm1 = fabs(mx1)/mm1;
    if (mm1 > fabs(mx0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }
  else {
    mm1 = fabs(mx1);
    mm1 = fabs(my1)/mm1;
    if (mm1 > fabs(my0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }



  mm1 = fabs(mx0) + fabs(my0);
  coord n = {mx0/mm1, my0/mm1};

  return n;
}
#line 13 "/home/gongminjiang/basilisk/src/fractions.h"






#line 1 "myc2d.h"
#line 1 "/home/gongminjiang/basilisk/src/myc2d.h"





static void _stencil_mycs (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;   
  
  
   


_stencil_val(c,-1,1,0); _stencil_val(c,0,1,0); _stencil_val(c,1,1,0); 


     
#line 14
_stencil_val(c,-1,-1,0); _stencil_val(c,0,-1,0); _stencil_val(c,1,-1,0); 
     _stencil_val(c,1,-1,0); _stencil_val(c,1,0,0); _stencil_val(c,1,1,0); 
     _stencil_val(c,-1,-1,0); _stencil_val(c,-1,0,0); _stencil_val(c,-1,1,0);
            
      
   
            
     
       
    



   
    





_stencil_val(c,-1,-1,0);_stencil_val(c,-1,0,0); _stencil_val(c,-1,1,0); 


     
  


      
#line 35
_stencil_val(c,1,-1,0);_stencil_val(c,1,0,0); _stencil_val(c,1,1,0);  
     
        _stencil_val(c,-1,-1,0);_stencil_val(c,0,-1,0); _stencil_val(c,1,-1,0);
       _stencil_val(c,-1,1,0);_stencil_val(c,0,1,0); _stencil_val(c,1,1,0);    
        
        
    
      
       
       
   
    
      
       


  
        
        
    
      
       
       
   



      
  

  
#line 64
return ;
}
#line 20 "/home/gongminjiang/basilisk/src/fractions.h"
#line 41 "/home/gongminjiang/basilisk/src/fractions.h"
void fraction_refine (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;





  double cc = val(c,0,0,0);
  if (cc <= 0. || cc >= 1.)
    {foreach_child()
      val(c,0,0,0) = cc;end_foreach_child()}
  else {




    coord n = mycs (point, c);
    double alpha = line_alpha (cc, n);






    {foreach_child() {
      static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
      coord nc;
      
 nc.x = child.x*n.x;
 
#line 69
nc.y = child.y*n.y;
      val(c,0,0,0) = rectangle_fraction (nc, alpha, a, b);
    }end_foreach_child()}
  }
}











static void alpha_refine (Point point, scalar alpha)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  vector n = _attribute[alpha.i].n;
  double alphac = 2.*val(alpha,0,0,0);
  coord m;
  
    m.x = val(n.x,0,0,0);
    
#line 91
m.y = val(n.y,0,0,0);
  {foreach_child() {
    val(alpha,0,0,0) = alphac;
    
      val(alpha,0,0,0) -= child.x*m.x/2.;
      
#line 95
val(alpha,0,0,0) -= child.y*m.y/2.;
  }end_foreach_child()}
}
#line 121 "/home/gongminjiang/basilisk/src/fractions.h"
struct Fractions {
  scalar Phi;
  scalar c;
  vector s;
  double val;
};

     
void fractions (struct Fractions a)
{tracing("fractions","/home/gongminjiang/basilisk/src/fractions.h",129);
  scalar Phi = a.Phi;
  scalar c = a.c;
  vector   s=(a.s).x.i?(a.s):new_face_vector("s");
  double val = a.val;
#line 145 "/home/gongminjiang/basilisk/src/fractions.h"
  vector p;
  p.x = s.y; p.y = s.x;
#line 155 "/home/gongminjiang/basilisk/src/fractions.h"
  foreach_face_stencil(){_stencil_is_face_y(){ {





_stencil_val(Phi,0,0,0);_stencil_val(Phi,1,0,0);{ {






      _stencil_val_a(p.x,0,0,0);_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,0,0); _stencil_val(Phi,1,0,0);
_stencil_val(Phi,0,0,0);
 {_stencil_val_a(p.x,0,0,0); _stencil_val(p.x,0,0,0);   }     
         
    
#line 171
}
      








{_stencil_val_a(p.x,0,0,0);_stencil_val(Phi,0,0,0); _stencil_val(Phi,1,0,0);       }}





           
#line 180 "/home/gongminjiang/basilisk/src/fractions.h"
    
  
}}end__stencil_is_face_y()
#line 155
_stencil_is_face_x(){ {





_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,1,0);{ {






      _stencil_val_a(p.y,0,0,0);_stencil_val(Phi,0,0,0);_stencil_val(Phi,0,0,0); _stencil_val(Phi,0,1,0);
_stencil_val(Phi,0,0,0);
 {_stencil_val_a(p.y,0,0,0); _stencil_val(p.y,0,0,0);   }     
         
    
#line 171
}
      








{_stencil_val_a(p.y,0,0,0);_stencil_val(Phi,0,0,0); _stencil_val(Phi,0,1,0);       }}





           
#line 180 "/home/gongminjiang/basilisk/src/fractions.h"
    
  
}}end__stencil_is_face_x()}end_foreach_face_stencil();
#line 155 "/home/gongminjiang/basilisk/src/fractions.h"
  {foreach_face_generic(){is_face_y(){ {





    if ((val(Phi,0,0,0) - val)*(val(Phi,1,0,0) - val) < 0.) {






      val(p.x,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,1,0,0));
      if (val(Phi,0,0,0) < val)
 val(p.x,0,0,0) = 1. - val(p.x,0,0,0);
    }
#line 180 "/home/gongminjiang/basilisk/src/fractions.h"
    else
      val(p.x,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,1,0,0) > val);
  }}end_is_face_y()
#line 155
is_face_x(){ {





    if ((val(Phi,0,0,0) - val)*(val(Phi,0,1,0) - val) < 0.) {






      val(p.y,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,0,1,0));
      if (val(Phi,0,0,0) < val)
 val(p.y,0,0,0) = 1. - val(p.y,0,0,0);
    }
#line 180 "/home/gongminjiang/basilisk/src/fractions.h"
    else
      val(p.y,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,0,1,0) > val);
  }}end_is_face_x()}end_foreach_face_generic();}
#line 205 "/home/gongminjiang/basilisk/src/fractions.h"
  scalar s_z = c;
  foreach_stencil()

  {    
#line 240 "/home/gongminjiang/basilisk/src/fractions.h"
    
    
     { 
_stencil_val(p.y,0,0,0); _stencil_val(p.y,1,0,0);  
       
       
    
#line 245
} 
#line 242
{ 
_stencil_val(p.x,0,0,0); _stencil_val(p.x,0,1,0);  
       
       
    
#line 245
}





{
      {_stencil_val_a(s_z,0,0,0); _stencil_val(p.x,0,0,0); } 
{      





      
   






      
      for (int i = 0; i <= 1; i++)
 {
   {_stencil_val(p.x,0,i,0); _stencil_val(p.x,0,i,0); {       
     _stencil_val(p.x,0,i,0);_stencil_val(Phi,0,i,0); 
          
     
   }      }
   
#line 270
{_stencil_val(p.y,i,0,0); _stencil_val(p.y,i,0,0); {       
     _stencil_val(p.y,i,0,0);_stencil_val(Phi,i,0,0); 
          
     
   }      }}








{
 {_stencil_val_a(s_z,0,0,0);_stencil_val(p.x,0,0,0); _stencil_val(p.y,0,0,0);   }
{
 {_stencil_val_a(s_z,0,0,0);     } 
{



 _stencil_val_a(s_z,0,0,0);  

      }}}
#line 283 "/home/gongminjiang/basilisk/src/fractions.h"
         
          
      
    







}}





       
    
  
#line 295
}end_foreach_stencil();
  {
#line 206
foreach()

  {
#line 240 "/home/gongminjiang/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
     {
      n.x = val(p.y,0,0,0) - val(p.y,1,0,0);
      nn += fabs(n.x);
    } 
#line 242
{
      n.y = val(p.x,0,0,0) - val(p.x,0,1,0);
      nn += fabs(n.y);
    }





    if (nn == 0.)
      val(s_z,0,0,0) = val(p.x,0,0,0);
    else {





      
 n.x /= nn;
 
#line 260
n.y /= nn;






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
   if (val(p.x,0,i,0) > 0. && val(p.x,0,i,0) < 1.) {
     double a = sign(val(Phi,0,i,0) - val)*(val(p.x,0,i,0) - 0.5);
     alpha += n.x*a + n.y*(i - 0.5);
     ni++;
   }
   
#line 270
if (val(p.y,i,0,0) > 0. && val(p.y,i,0,0) < 1.) {
     double a = sign(val(Phi,i,0,0) - val)*(val(p.y,i,0,0) - 0.5);
     alpha += n.y*a + n.x*(i - 0.5);
     ni++;
   }}
#line 283 "/home/gongminjiang/basilisk/src/fractions.h"
      if (ni == 0)
 val(s_z,0,0,0) = max (val(p.x,0,0,0), val(p.y,0,0,0));
      else if (ni != 4)
 val(s_z,0,0,0) = line_area (n.x, n.y, alpha/ni);
      else {



 val(s_z,0,0,0) = 0.;

      }
    }
  }end_foreach();}{if(!(a.s).x.i)delete((scalar*)((vector[]){s,{{-1},{-1}}}));}
#line 347 "/home/gongminjiang/basilisk/src/fractions.h"
end_tracing("fractions","/home/gongminjiang/basilisk/src/fractions.h",347);}
#line 391 "/home/gongminjiang/basilisk/src/fractions.h"
coord youngs_normal (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  coord n;
  double nn = 0.;
  if (!(2 == 2)) qassert ("/home/gongminjiang/basilisk/src/fractions.h", 395, "dimension == 2");
   {
    n.x = (val(c,-1,1,0) + 2.*val(c,-1,0,0) + val(c,-1,-1,0) -
    val(c,+1,1,0) - 2.*val(c,+1,0,0) - val(c,+1,-1,0));
    nn += fabs(n.x);
  } 
#line 396
{
    n.y = (val(c,1,-1,0) + 2.*val(c,0,-1,0) + val(c,-1,-1,0) -
    val(c,1,+1,0) - 2.*val(c,0,+1,0) - val(c,-1,+1,0));
    nn += fabs(n.y);
  }

  if (nn > 0.)
    {
      n.x /= nn;
      
#line 404
n.y /= nn;}
  else
    n.x = 1.;
  return n;
}





coord facet_normal (Point point, scalar c, vector s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (s.x.i >= 0) {
    coord n;
    double nn = 0.;
     {
      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
      nn += fabs(n.x);
    } 
#line 419
{
      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
      nn += fabs(n.y);
    }
    if (nn > 0.)
      {
 n.x /= nn;
 
#line 425
n.y /= nn;}
    else
      {
 n.x = 1./2;
 
#line 428
n.y = 1./2;}
    return n;
  }
  return mycs (point, c);
}






#line 414
static void _stencil_facet_normal (Point point, scalar c, vector s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (s.x.i >= 0) {    
    
    
     { 
_stencil_val(s.x,0,0,0); _stencil_val(s.x,1,0,0);  
       
       
    
#line 422
} 
#line 419
{ 
_stencil_val(s.y,0,0,0); _stencil_val(s.y,0,1,0);  
       
       
    
#line 422
}
      
   
       
    
       
  
    return ;
  } 
_stencil_mycs (point, c);
  
#line 431
return;
}
#line 441 "/home/gongminjiang/basilisk/src/fractions.h"
     
void reconstruction (const scalar c, vector n, scalar alpha)
{tracing("reconstruction","/home/gongminjiang/basilisk/src/fractions.h",442);
  foreach_stencil() {





_stencil_val(c,0,0,0); _stencil_val(c,0,0,0);{ {
      _stencil_val_a(alpha,0,0,0);  
      
 {_stencil_val_a(n.x,0,0,0);  }
 
#line 453
{_stencil_val_a(n.y,0,0,0);  }
    } 
{  






       _stencil_mycs (point, c);
      
 {_stencil_val_a(n.x,0,0,0);  }
 
#line 464
{_stencil_val_a(n.y,0,0,0);  }
      _stencil_val_a(alpha,0,0,0);_stencil_val(c,0,0,0);    
    }}





          
    
  
#line 467
}end_foreach_stencil();
  {
#line 444
foreach() {





    if (val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) {
      val(alpha,0,0,0) = 0.;
      
 val(n.x,0,0,0) = 0.;
 
#line 453
val(n.y,0,0,0) = 0.;
    }
    else {






      coord m = mycs (point, c);
      
 val(n.x,0,0,0) = m.x;
 
#line 464
val(n.y,0,0,0) = m.y;
      val(alpha,0,0,0) = line_alpha (val(c,0,0,0), m);
    }
  }end_foreach();}
#line 476 "/home/gongminjiang/basilisk/src/fractions.h"
  
    _attribute[n.x.i].refine = _attribute[n.x.i].prolongation = refine_injection;
    
#line 477
_attribute[n.y.i].refine = _attribute[n.y.i].prolongation = refine_injection;




  _attribute[alpha.i].n = n;
  _attribute[alpha.i].refine = _attribute[alpha.i].prolongation = alpha_refine;

end_tracing("reconstruction","/home/gongminjiang/basilisk/src/fractions.h",485);}
#line 505 "/home/gongminjiang/basilisk/src/fractions.h"
struct OutputFacets {
  scalar c;
  FILE * fp;
  vector s;
};

     
void output_facets (struct OutputFacets p)
{tracing("output_facets","/home/gongminjiang/basilisk/src/fractions.h",512);
  scalar c = p.c;
  vector s = p.s;
  if (!p.fp) p.fp = fout;
  if (!s.x.i) s.x.i = -1;

  foreach_stencil()
    {_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {  
       _stencil_facet_normal (point, c, s);     
      _stencil_val(c,0,0,0); 

            
           
        
     
 
#line 538 "/home/gongminjiang/basilisk/src/fractions.h"
    }        }end_foreach_stencil();

  {
#line 519
foreach()
    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha = line_alpha (val(c,0,0,0), n);

      coord segment[2];
      if (facets (n, alpha, segment) == 2)
 fprintf (p.fp, "%g %g\n%g %g\n\n",
   x + segment[0].x*Delta, y + segment[0].y*Delta,
   x + segment[1].x*Delta, y + segment[1].y*Delta);
#line 538 "/home/gongminjiang/basilisk/src/fractions.h"
    }end_foreach();}

  fflush (p.fp);
end_tracing("output_facets","/home/gongminjiang/basilisk/src/fractions.h",541);}







     
double interface_area (scalar c)
{tracing("interface_area","/home/gongminjiang/basilisk/src/fractions.h",550);
  double area = 0.;
  foreach_stencil ()
    {_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {   
       _stencil_mycs (point, c);     
      _stencil_val(c,0,0,0); 
          
    }        }end_foreach_stencil();
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:area)){
#line 553
foreach ()
    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = mycs (point, c), p;
      double alpha = line_alpha (val(c,0,0,0), n);
      area += pow(Delta, 2 - 1)*line_length_center(n,alpha,&p);
    }end_foreach();mpi_all_reduce_array(&area,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
  
#line 559
{end_tracing("interface_area","/home/gongminjiang/basilisk/src/fractions.h",559);return area;}
end_tracing("interface_area","/home/gongminjiang/basilisk/src/fractions.h",560);}
#line 36 "/home/gongminjiang/basilisk/src/vof.h"
#line 44 "/home/gongminjiang/basilisk/src/vof.h"
extern scalar * interfaces;
extern vector uf;
extern double dt;








static double vof_concentration_gradient_x (Point point, scalar c, scalar t)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  static const double cmin = 0.5;
  double cl = val(c,-1,0,0), cc = val(c,0,0,0), cr = val(c,1,0,0);
  if (_attribute[t.i].inverse)
    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  if (cc >= cmin && _attribute[t.i].gradient != zero) {
    if (cr >= cmin) {
      if (cl >= cmin) {
 if (_attribute[t.i].gradient)
   return _attribute[t.i].gradient (val(t,-1,0,0)/cl, val(t,0,0,0)/cc, val(t,1,0,0)/cr)/Delta;
 else
   return (val(t,1,0,0)/cr - val(t,-1,0,0)/cl)/(2.*Delta);
      }
      else
 return (val(t,1,0,0)/cr - val(t,0,0,0)/cc)/Delta;
    }
    else if (cl >= cmin)
      return (val(t,0,0,0)/cc - val(t,-1,0,0)/cl)/Delta;
  }
  return 0.;
}

#line 55
static double vof_concentration_gradient_y (Point point, scalar c, scalar t)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  static const double cmin = 0.5;
  double cl = val(c,0,-1,0), cc = val(c,0,0,0), cr = val(c,0,1,0);
  if (_attribute[t.i].inverse)
    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  if (cc >= cmin && _attribute[t.i].gradient != zero) {
    if (cr >= cmin) {
      if (cl >= cmin) {
 if (_attribute[t.i].gradient)
   return _attribute[t.i].gradient (val(t,0,-1,0)/cl, val(t,0,0,0)/cc, val(t,0,1,0)/cr)/Delta;
 else
   return (val(t,0,1,0)/cr - val(t,0,-1,0)/cl)/(2.*Delta);
      }
      else
 return (val(t,0,1,0)/cr - val(t,0,0,0)/cc)/Delta;
    }
    else if (cl >= cmin)
      return (val(t,0,0,0)/cc - val(t,0,-1,0)/cl)/Delta;
  }
  return 0.;
}









#line 55
static void _stencil_vof_concentration_gradient_x (Point point, scalar c, scalar t)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;           
  
   _stencil_val(c,1,0,0); _stencil_val(c,0,0,0); _stencil_val(c,-1,0,0); 


{
{ {
{ {
 if (_attribute[t.i].gradient)
   {_stencil_val(t,-1,0,0); _stencil_val(t,0,0,0); _stencil_val(t,1,0,0);  }
 else
   {_stencil_val(t,1,0,0); _stencil_val(t,-1,0,0);  }
      }
 
{_stencil_val(t,1,0,0); _stencil_val(t,0,0,0);  }}
         
      
    
#line 71
}
      
{_stencil_val(t,0,0,0); _stencil_val(t,-1,0,0);  }}
       
        
  
#line 74
} 
  
                  
         
  
#line 75
return ;
}

#line 55
static void _stencil_vof_concentration_gradient_y (Point point, scalar c, scalar t)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;           
  
   _stencil_val(c,0,1,0); _stencil_val(c,0,0,0); _stencil_val(c,0,-1,0); 


{
{ {
{ {
 if (_attribute[t.i].gradient)
   {_stencil_val(t,0,-1,0); _stencil_val(t,0,0,0); _stencil_val(t,0,1,0);  }
 else
   {_stencil_val(t,0,1,0); _stencil_val(t,0,-1,0);  }
      }
 
{_stencil_val(t,0,1,0); _stencil_val(t,0,0,0);  }}
         
      
    
#line 71
}
      
{_stencil_val(t,0,0,0); _stencil_val(t,0,-1,0);  }}
       
        
  
#line 74
} 
  
                  
         
  
#line 75
return ;
}






static void vof_concentration_refine (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;

#line 84
if(!is_constant(cm)){{
  scalar f = _attribute[s.i].c;
  if (val(cm,0,0,0) == 0. || (!_attribute[s.i].inverse && val(f,0,0,0) <= 0.) || (_attribute[s.i].inverse && val(f,0,0,0) >= 1.))
    {foreach_child()
      val(s,0,0,0) = 0.;end_foreach_child()}
  else {
    coord g;
    
      g.x = Delta*vof_concentration_gradient_x (point, f, s);
      
#line 92
g.y = Delta*vof_concentration_gradient_y (point, f, s);
    double sc = _attribute[s.i].inverse ? val(s,0,0,0)/(1. - val(f,0,0,0)) : val(s,0,0,0)/val(f,0,0,0), cmc = 4.*val(cm,0,0,0);
    {foreach_child() {
      val(s,0,0,0) = sc;
      
 val(s,0,0,0) += child.x*g.x*val(cm,-child.x,0,0)/cmc;
 
#line 97
val(s,0,0,0) += child.y*g.y*val(cm,0,-child.y,0)/cmc;
      val(s,0,0,0) *= _attribute[s.i].inverse ? 1. - val(f,0,0,0) : val(f,0,0,0);
    }end_foreach_child()}
  }
}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

#line 84
{
  scalar f = _attribute[s.i].c;
  if (_const_cm == 0. || (!_attribute[s.i].inverse && val(f,0,0,0) <= 0.) || (_attribute[s.i].inverse && val(f,0,0,0) >= 1.))
    {foreach_child()
      val(s,0,0,0) = 0.;end_foreach_child()}
  else {
    coord g;
    
      g.x = Delta*vof_concentration_gradient_x (point, f, s);
      
#line 92
g.y = Delta*vof_concentration_gradient_y (point, f, s);
    double sc = _attribute[s.i].inverse ? val(s,0,0,0)/(1. - val(f,0,0,0)) : val(s,0,0,0)/val(f,0,0,0), cmc = 4.*_const_cm;
    {foreach_child() {
      val(s,0,0,0) = sc;
      
 val(s,0,0,0) += child.x*g.x*_const_cm/cmc;
 
#line 97
val(s,0,0,0) += child.y*g.y*_const_cm/cmc;
      val(s,0,0,0) *= _attribute[s.i].inverse ? 1. - val(f,0,0,0) : val(f,0,0,0);
    }end_foreach_child()}
  }
}}

#line 101
}





static int defaults_1_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0);*ip=i;*tp=t;return ret;}      static int defaults_1(const int i,const double t,Event *_ev){tracing("defaults_1","/home/gongminjiang/basilisk/src/vof.h",107);
{
  {scalar*_i=(scalar*)( interfaces);if(_i)for(scalar c=*_i;(&c)->i>=0;c=*++_i){ {
    _attribute[c.i].refine = _attribute[c.i].prolongation = fraction_refine;
    _attribute[c.i].dirty = true;
    scalar * tracers = _attribute[c.i].tracers;
    {scalar*_i=(scalar*)( tracers);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){ {
      _attribute[t.i].restriction = restriction_volume_average;
      _attribute[t.i].refine = _attribute[t.i].prolongation = vof_concentration_refine;
      _attribute[t.i].dirty = true;
      _attribute[t.i].c = c;
    }}}
  }}}
}{end_tracing("defaults_1","/home/gongminjiang/basilisk/src/vof.h",120);return 0;}end_tracing("defaults_1","/home/gongminjiang/basilisk/src/vof.h",120);}






static int defaults_2_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0);*ip=i;*tp=t;return ret;}      static int defaults_2(const int i,const double t,Event *_ev){tracing("defaults_2","/home/gongminjiang/basilisk/src/vof.h",127);
{
  {scalar*_i=(scalar*)( interfaces);if(_i)for(scalar c=*_i;(&c)->i>=0;c=*++_i){ {
    scalar * tracers = _attribute[c.i].tracers;
    {scalar*_i=(scalar*)( tracers);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){
      _attribute[t.i].depends = list_add (_attribute[t.i].depends, c);}}
  }}}
}{end_tracing("defaults_2","/home/gongminjiang/basilisk/src/vof.h",134);return 0;}end_tracing("defaults_2","/home/gongminjiang/basilisk/src/vof.h",134);}





static int stability_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int stability_0(const int i,const double t,Event *_ev){tracing("stability_0","/home/gongminjiang/basilisk/src/vof.h",140); {
  if (CFL > 0.5)
    CFL = 0.5;
}{end_tracing("stability_0","/home/gongminjiang/basilisk/src/vof.h",143);return 0;}end_tracing("stability_0","/home/gongminjiang/basilisk/src/vof.h",143);}
#line 157 "/home/gongminjiang/basilisk/src/vof.h"

static void sweep_x (scalar c, scalar cc, scalar * tcl)
{
  vector  n=new_vector("n");
  scalar  alpha=new_scalar("alpha"),  flux=new_scalar("flux");
  double cfl = 0.;
#line 171 "/home/gongminjiang/basilisk/src/vof.h"
  scalar * tracers = _attribute[c.i].tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    {scalar*_i=(scalar*)( tracers);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){ {
      scalar gf = new_scalar("gf"), flux = new_scalar("flux");
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }}}




    foreach_stencil() {
      scalar t, gf;
      {scalar*_i0=gfl;scalar*_i1= tracers;if(_i0)for(gf=*_i0,t=*_i1;_i0->i>= 0;gf=*++_i0,t=*++_i1){
 {_stencil_val_a(gf,0,0,0); _stencil_vof_concentration_gradient_x (point, c, t); }}}
    }end_foreach_stencil();




    {
#line 182
foreach() {
      scalar t, gf;
      {scalar*_i0=gfl;scalar*_i1= tracers;if(_i0)for(gf=*_i0,t=*_i1;_i0->i>= 0;gf=*++_i0,t=*++_i1){
 val(gf,0,0,0) = vof_concentration_gradient_x (point, c, t);}}
    }end_foreach();}
  }






  reconstruction (c, n, alpha);
  foreach_face_stencil()_stencil_is_face_x(){ {       






    _stencil_val(fm.x,0,0,0); _stencil_val(uf.x,0,0,0);     
    







_stencil_val(fm.x,0,0,0);_stencil_val(cm,0,0,0);
      {_stencil_val(fm.x,0,0,0);_stencil_val(cm,0,0,0);    }              
       
      
      







         
#line 225 "/home/gongminjiang/basilisk/src/vof.h"
    
_stencil_val(alpha, o_stencil,0,0);_stencil_val_higher_dimension;_stencil_val(n.y, o_stencil,0,0);_stencil_val(n.x,o_stencil,0,0);
#line 225
_stencil_val(c, o_stencil,0,0);_stencil_val(c, o_stencil,0,0);_stencil_val(c,o_stencil,0,0);





    


_stencil_val_a(flux,0,0,0);_stencil_val(uf.x,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c, o_stencil,0,0);


{ {       
 _stencil_val(gf,o_stencil,0,0);_stencil_val(t, o_stencil,0,0);
 _stencil_val_a(tflux,0,0,0);_stencil_val(uf.x,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_x()end_foreach_face_stencil();
  
#line 195
if(!is_constant(fm.x) && !is_constant(cm)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_x(){ {






    double un = val(uf.x,0,0,0)*dt/(Delta*val(fm.x,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*val(fm.x,0,0,0)*s/(val(cm,0,0,0) + 0.) > cfl)
      cfl = un*val(fm.x,0,0,0)*s/(val(cm,0,0,0) + 0.);
#line 225 "/home/gongminjiang/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_x()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 253
}else if(is_constant(fm.x) && !is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_x(){ {






    double un = val(uf.x,0,0,0)*dt/(Delta*_const_fm.x + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*_const_fm.x*s/(val(cm,0,0,0) + 0.) > cfl)
      cfl = un*_const_fm.x*s/(val(cm,0,0,0) + 0.);
#line 225 "/home/gongminjiang/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_x()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 253
}else if(!is_constant(fm.x) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_x(){ {






    double un = val(uf.x,0,0,0)*dt/(Delta*val(fm.x,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*val(fm.x,0,0,0)*s/(_const_cm + 0.) > cfl)
      cfl = un*val(fm.x,0,0,0)*s/(_const_cm + 0.);
#line 225 "/home/gongminjiang/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_x()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 253
}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_x(){ {






    double un = val(uf.x,0,0,0)*dt/(Delta*_const_fm.x + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*_const_fm.x*s/(_const_cm + 0.) > cfl)
      cfl = un*_const_fm.x*s/(_const_cm + 0.);
#line 225 "/home/gongminjiang/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), _val_higher_dimension}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_x()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 253
}
  delete (gfl); pfree (gfl,__func__,__FILE__,__LINE__);




  if (cfl > 0.5 + 1e-6)
    fprintf (ferr,
      "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n",
      cfl - 0.5), fflush (ferr);
#line 282 "/home/gongminjiang/basilisk/src/vof.h"
  foreach_stencil() {
    _stencil_val_r(c,0,0,0);_stencil_val(flux,0,0,0); _stencil_val(flux,1,0,0); _stencil_val(cc,0,0,0);_stencil_val(uf.x,1,0,0); _stencil_val(uf.x,0,0,0);_stencil_val(cm,0,0,0);     
    scalar t, tc, tflux;
    {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
      {_stencil_val_r(t,0,0,0);_stencil_val(tflux,0,0,0); _stencil_val(tflux,1,0,0); _stencil_val(tc,0,0,0);_stencil_val(uf.x,1,0,0); _stencil_val(uf.x,0,0,0);_stencil_val(cm,0,0,0);     }}}
  }end_foreach_stencil();
#line 282 "/home/gongminjiang/basilisk/src/vof.h"
  if(!is_constant(cm)){{foreach() {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,1,0,0) + val(cc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/(val(cm,0,0,0)*Delta);
    scalar t, tc, tflux;
    {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,1,0,0) + val(tc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/(val(cm,0,0,0)*Delta);}}
  }end_foreach();}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
#line 282 "/home/gongminjiang/basilisk/src/vof.h"
  {foreach() {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,1,0,0) + val(cc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/(_const_cm*Delta);
    scalar t, tc, tflux;
    {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,1,0,0) + val(tc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/(_const_cm*Delta);}}
  }end_foreach();}}
#line 305 "/home/gongminjiang/basilisk/src/vof.h"
  delete (tfluxl); pfree (tfluxl,__func__,__FILE__,__LINE__);delete((scalar*)((scalar[]){flux,alpha,n.x,n.y,{-1}}));
}

#line 158
static void sweep_y (scalar c, scalar cc, scalar * tcl)
{
  vector  n=new_vector("n");
  scalar  alpha=new_scalar("alpha"),  flux=new_scalar("flux");
  double cfl = 0.;
#line 171 "/home/gongminjiang/basilisk/src/vof.h"
  scalar * tracers = _attribute[c.i].tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    {scalar*_i=(scalar*)( tracers);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){ {
      scalar gf = new_scalar("gf"), flux = new_scalar("flux");
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }}}




    foreach_stencil() {
      scalar t, gf;
      {scalar*_i0=gfl;scalar*_i1= tracers;if(_i0)for(gf=*_i0,t=*_i1;_i0->i>= 0;gf=*++_i0,t=*++_i1){
 {_stencil_val_a(gf,0,0,0); _stencil_vof_concentration_gradient_y (point, c, t); }}}
    }end_foreach_stencil();




    {
#line 182
foreach() {
      scalar t, gf;
      {scalar*_i0=gfl;scalar*_i1= tracers;if(_i0)for(gf=*_i0,t=*_i1;_i0->i>= 0;gf=*++_i0,t=*++_i1){
 val(gf,0,0,0) = vof_concentration_gradient_y (point, c, t);}}
    }end_foreach();}
  }






  reconstruction (c, n, alpha);
  foreach_face_stencil()_stencil_is_face_y(){ {       






    _stencil_val(fm.y,0,0,0); _stencil_val(uf.y,0,0,0);     
    







_stencil_val(fm.y,0,0,0);_stencil_val(cm,0,0,0);
      {_stencil_val(fm.y,0,0,0);_stencil_val(cm,0,0,0);    }              
       
      
      







         
#line 225 "/home/gongminjiang/basilisk/src/vof.h"
    
_stencil_val(alpha,0, o_stencil,0);_stencil_val_higher_dimension;_stencil_val(n.x,0, o_stencil,0);_stencil_val(n.y,0,o_stencil,0);
#line 225
_stencil_val(c,0, o_stencil,0);_stencil_val(c,0, o_stencil,0);_stencil_val(c,0,o_stencil,0);





    


_stencil_val_a(flux,0,0,0);_stencil_val(uf.y,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c,0, o_stencil,0);


{ {       
 _stencil_val(gf,0,o_stencil,0);_stencil_val(t,0, o_stencil,0);
 _stencil_val_a(tflux,0,0,0);_stencil_val(uf.y,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_y()end_foreach_face_stencil();
  
#line 195
if(!is_constant(fm.y) && !is_constant(cm)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_y(){ {






    double un = val(uf.y,0,0,0)*dt/(Delta*val(fm.y,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*val(fm.y,0,0,0)*s/(val(cm,0,0,0) + 0.) > cfl)
      cfl = un*val(fm.y,0,0,0)*s/(val(cm,0,0,0) + 0.);
#line 225 "/home/gongminjiang/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_y()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 253
}else if(is_constant(fm.y) && !is_constant(cm)){struct{double x,y;}_const_fm={_constant[fm.y.i-_NVARMAX],_constant[fm.x.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_y(){ {






    double un = val(uf.y,0,0,0)*dt/(Delta*_const_fm.y + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*_const_fm.y*s/(val(cm,0,0,0) + 0.) > cfl)
      cfl = un*_const_fm.y*s/(val(cm,0,0,0) + 0.);
#line 225 "/home/gongminjiang/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_y()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 253
}else if(!is_constant(fm.y) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_y(){ {






    double un = val(uf.y,0,0,0)*dt/(Delta*val(fm.y,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*val(fm.y,0,0,0)*s/(_const_cm + 0.) > cfl)
      cfl = un*val(fm.y,0,0,0)*s/(_const_cm + 0.);
#line 225 "/home/gongminjiang/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_y()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 253
}else {struct{double x,y;}_const_fm={_constant[fm.y.i-_NVARMAX],_constant[fm.x.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_y(){ {






    double un = val(uf.y,0,0,0)*dt/(Delta*_const_fm.y + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*_const_fm.y*s/(_const_cm + 0.) > cfl)
      cfl = un*_const_fm.y*s/(_const_cm + 0.);
#line 225 "/home/gongminjiang/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.x,0,i,0), _val_higher_dimension}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_y()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 253
}
  delete (gfl); pfree (gfl,__func__,__FILE__,__LINE__);




  if (cfl > 0.5 + 1e-6)
    fprintf (ferr,
      "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n",
      cfl - 0.5), fflush (ferr);
#line 282 "/home/gongminjiang/basilisk/src/vof.h"
  foreach_stencil() {
    _stencil_val_r(c,0,0,0);_stencil_val(flux,0,0,0); _stencil_val(flux,0,1,0); _stencil_val(cc,0,0,0);_stencil_val(uf.y,0,1,0); _stencil_val(uf.y,0,0,0);_stencil_val(cm,0,0,0);     
    scalar t, tc, tflux;
    {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
      {_stencil_val_r(t,0,0,0);_stencil_val(tflux,0,0,0); _stencil_val(tflux,0,1,0); _stencil_val(tc,0,0,0);_stencil_val(uf.y,0,1,0); _stencil_val(uf.y,0,0,0);_stencil_val(cm,0,0,0);     }}}
  }end_foreach_stencil();
#line 282 "/home/gongminjiang/basilisk/src/vof.h"
  if(!is_constant(cm)){{foreach() {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,1,0) + val(cc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/(val(cm,0,0,0)*Delta);
    scalar t, tc, tflux;
    {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,0,1,0) + val(tc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/(val(cm,0,0,0)*Delta);}}
  }end_foreach();}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
#line 282 "/home/gongminjiang/basilisk/src/vof.h"
  {foreach() {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,1,0) + val(cc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/(_const_cm*Delta);
    scalar t, tc, tflux;
    {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,0,1,0) + val(tc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/(_const_cm*Delta);}}
  }end_foreach();}}
#line 305 "/home/gongminjiang/basilisk/src/vof.h"
  delete (tfluxl); pfree (tfluxl,__func__,__FILE__,__LINE__);delete((scalar*)((scalar[]){flux,alpha,n.x,n.y,{-1}}));
}






void vof_advection (scalar * interfaces, int i)
{
  {scalar*_i=(scalar*)( interfaces);if(_i)for(scalar c=*_i;(&c)->i>=0;c=*++_i){ {
#line 325 "/home/gongminjiang/basilisk/src/vof.h"
    scalar  cc=new_scalar("cc"), * tcl = NULL, * tracers = _attribute[c.i].tracers;
    {scalar*_i=(scalar*)( tracers);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){ {
      scalar tc = new_scalar("tc");
      tcl = list_append (tcl, tc);

      if (_attribute[t.i].refine != vof_concentration_refine) {
 _attribute[t.i].refine = _attribute[t.i].prolongation = vof_concentration_refine;
 _attribute[t.i].restriction = restriction_volume_average;
 _attribute[t.i].dirty = true;
 _attribute[t.i].c = c;
      }

    }}}
    foreach_stencil() {
      _stencil_val_a(cc,0,0,0);_stencil_val(c,0,0,0);    
      scalar t, tc;
      {scalar*_i0= tcl;scalar*_i1= tracers;if(_i0)for(tc=*_i0,t=*_i1;_i0->i>= 0;tc=*++_i0,t=*++_i1){ {
 if (_attribute[t.i].inverse)
   {_stencil_val_a(tc,0,0,0); _stencil_val(c,0,0,0); _stencil_val(t,0,0,0); _stencil_val(c,0,0,0);       }
 else
   {_stencil_val_a(tc,0,0,0); _stencil_val(c,0,0,0); _stencil_val(t,0,0,0);_stencil_val(c,0,0,0);      }
      }}}
    }end_foreach_stencil();
    {
#line 338
foreach() {
      val(cc,0,0,0) = (val(c,0,0,0) > 0.5);
      scalar t, tc;
      {scalar*_i0= tcl;scalar*_i1= tracers;if(_i0)for(tc=*_i0,t=*_i1;_i0->i>= 0;tc=*++_i0,t=*++_i1){ {
 if (_attribute[t.i].inverse)
   val(tc,0,0,0) = val(c,0,0,0) < 0.5 ? val(t,0,0,0)/(1. - val(c,0,0,0)) : 0.;
 else
   val(tc,0,0,0) = val(c,0,0,0) > 0.5 ? val(t,0,0,0)/val(c,0,0,0) : 0.;
      }}}
    }end_foreach();}






    void (* sweep[2]) (scalar, scalar, scalar *);
    int d = 0;
    
      sweep[d++] = sweep_x;
      
#line 357
sweep[d++] = sweep_y;
    for (d = 0; d < 2; d++)
      sweep[(i + d) % 2] (c, cc, tcl);
    delete (tcl), pfree (tcl,__func__,__FILE__,__LINE__);delete((scalar*)((scalar[]){cc,{-1}}));
  }}}
}

static int vof_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int vof_0(const int i,const double t,Event *_ev){tracing("vof_0","/home/gongminjiang/basilisk/src/vof.h",364);
  vof_advection (interfaces, i);{end_tracing("vof_0","/home/gongminjiang/basilisk/src/vof.h",365);return 0;}end_tracing("vof_0","/home/gongminjiang/basilisk/src/vof.h",365);}
#line 14 "/home/gongminjiang/basilisk/src/two-phase.h"

scalar  f={8}, * interfaces = ((scalar[]){{8},{-1}});
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;





vector  alphav={{9},{10}};
scalar  rhov={11};

static int defaults_3_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0);*ip=i;*tp=t;return ret;}      static int defaults_3(const int i,const double t,Event *_ev){tracing("defaults_3","/home/gongminjiang/basilisk/src/two-phase.h",25); {
  alpha = alphav;
  rho = rhov;





  if (mu1 || mu2)
    mu = new_face_vector("mu");




  display ((struct _display){"draw_vof (c = 'f');"});
}{end_tracing("defaults_3","/home/gongminjiang/basilisk/src/two-phase.h",40);return 0;}end_tracing("defaults_3","/home/gongminjiang/basilisk/src/two-phase.h",40);}
#line 64 "/home/gongminjiang/basilisk/src/two-phase.h"
static int tracer_advection_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int tracer_advection_0(const int i,const double t,Event *_ev){tracing("tracer_advection_0","/home/gongminjiang/basilisk/src/two-phase.h",64);
{
#line 90 "/home/gongminjiang/basilisk/src/two-phase.h"
  _attribute[f.i].prolongation = refine_bilinear;
  _attribute[f.i].dirty = true;

}{end_tracing("tracer_advection_0","/home/gongminjiang/basilisk/src/two-phase.h",93);return 0;}end_tracing("tracer_advection_0","/home/gongminjiang/basilisk/src/two-phase.h",93);}

static int properties_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int properties_0(const int i,const double t,Event *_ev){tracing("properties_0","/home/gongminjiang/basilisk/src/two-phase.h",95);
{
  foreach_face_stencil(){_stencil_is_face_x(){ {    
     _stencil_val(f,-1,0,0);_stencil_val(f,0,0,0);
    _stencil_val_a(alphav.x,0,0,0); _stencil_val(fm.x,0,0,0);     
    if (mu1 || mu2) {
      vector muv = mu;
      _stencil_val_a(muv.x,0,0,0); _stencil_val(fm.x,0,0,0);     
    }
  }}end__stencil_is_face_x()
#line 97
_stencil_is_face_y(){ {    
     _stencil_val(f,0,-1,0);_stencil_val(f,0,0,0);
    _stencil_val_a(alphav.y,0,0,0); _stencil_val(fm.y,0,0,0);     
    if (mu1 || mu2) {
      vector muv = mu;
      _stencil_val_a(muv.y,0,0,0); _stencil_val(fm.y,0,0,0);     
    }
  }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#line 97
if(!is_constant(fm.x)){{foreach_face_generic(){is_face_x(){ {
    double ff = (val(f,0,0,0) + val(f,-1,0,0))/2.;
    val(alphav.x,0,0,0) = val(fm.x,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.x,0,0,0) = val(fm.x,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  }}end_is_face_x()
#line 97
is_face_y(){ {
    double ff = (val(f,0,0,0) + val(f,0,-1,0))/2.;
    val(alphav.y,0,0,0) = val(fm.y,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.y,0,0,0) = val(fm.y,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  }}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  {
#line 97
foreach_face_generic(){is_face_x(){ {
    double ff = (val(f,0,0,0) + val(f,-1,0,0))/2.;
    val(alphav.x,0,0,0) = _const_fm.x/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.x,0,0,0) = _const_fm.x*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  }}end_is_face_x()
#line 97
is_face_y(){ {
    double ff = (val(f,0,0,0) + val(f,0,-1,0))/2.;
    val(alphav.y,0,0,0) = _const_fm.y/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.y,0,0,0) = _const_fm.y*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  }}end_is_face_y()}end_foreach_face_generic();}}

  foreach_stencil()
    {_stencil_val_a(rhov,0,0,0); _stencil_val(cm,0,0,0);_stencil_val(f,0,0,0);     }end_foreach_stencil();

  
#line 106
if(!is_constant(cm)){{foreach()
    val(rhov,0,0,0) = val(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2);end_foreach();}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

  {
#line 106
foreach()
    val(rhov,0,0,0) = _const_cm*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2);end_foreach();}}


  _attribute[f.i].prolongation = fraction_refine;
  _attribute[f.i].dirty = true;

}{end_tracing("properties_0","/home/gongminjiang/basilisk/src/two-phase.h",113);return 0;}end_tracing("properties_0","/home/gongminjiang/basilisk/src/two-phase.h",113);}
#line 16 "atomisation.c"
#line 1 "tension.h"
#line 1 "/home/gongminjiang/basilisk/src/tension.h"
#line 15 "/home/gongminjiang/basilisk/src/tension.h"
#line 1 "iforce.h"
#line 1 "/home/gongminjiang/basilisk/src/iforce.h"
#line 20 "/home/gongminjiang/basilisk/src/iforce.h"








static int defaults_4_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0);*ip=i;*tp=t;return ret;}

      static int defaults_4(const int i,const double t,Event *_ev){tracing("defaults_4","/home/gongminjiang/basilisk/src/iforce.h",30); {
  if (is_constant(a.x)) {
    a = new_face_vector("a");
    foreach_face_stencil(){_stencil_is_face_x(){
      {_stencil_val_a(a.x,0,0,0);  }}end__stencil_is_face_x()
#line 33
_stencil_is_face_y(){
      {_stencil_val_a(a.y,0,0,0);  }}end__stencil_is_face_y()}end_foreach_face_stencil();
    
#line 33
if(!is_constant(a.x)){{foreach_face_generic(){is_face_x(){
      val(a.x,0,0,0) = 0.;}end_is_face_x()
#line 33
is_face_y(){
      val(a.y,0,0,0) = 0.;}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX]};NOT_UNUSED(_const_a);
    {
#line 33
foreach_face_generic(){is_face_x(){
      _const_a.x = 0.;}end_is_face_x()
#line 33
is_face_y(){
      _const_a.y = 0.;}end_is_face_y()}end_foreach_face_generic();}}
  }
}{end_tracing("defaults_4","/home/gongminjiang/basilisk/src/iforce.h",36);return 0;}end_tracing("defaults_4","/home/gongminjiang/basilisk/src/iforce.h",36);}






static int acceleration_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int acceleration_0(const int i,const double t,Event *_ev){tracing("acceleration_0","/home/gongminjiang/basilisk/src/iforce.h",43);
{





  scalar * list = NULL;
  {scalar*_i=(scalar*)( interfaces);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
    if (_attribute[f.i].phi.i) {
      list = list_add (list, f);






      foreach_stencil()
 {_stencil_val_a(f,0,0,0);_stencil_val(f,0,0,0);     }end_foreach_stencil();






      {
#line 60
foreach()
 val(f,0,0,0) = clamp (val(f,0,0,0), 0., 1.);end_foreach();}
    }}}
#line 72 "/home/gongminjiang/basilisk/src/iforce.h"
  {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){ {
    _attribute[f.i].prolongation = _attribute[p.i].prolongation;
    _attribute[f.i].dirty = true;
  }}}
#line 86 "/home/gongminjiang/basilisk/src/iforce.h"
  vector ia = a;
  foreach_face_stencil(){_stencil_is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,-1,0,0); _stencil_val(fm.x,0,0,0); {
#line 99 "/home/gongminjiang/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,-1,0,0);
   
#line 104
_stencil_val(phi,-1,0,0); 
#line 103
_stencil_val(phi,0,0,0);
   
#line 103
_stencil_val(phi,0,0,0); 
#line 102
_stencil_val(phi,-1,0,0);_stencil_val(phi,0,0,0); 
#line 101
_stencil_val(phi,-1,0,0);_stencil_val(phi,0,0,0);

 



_stencil_val_r(ia.x,0,0,0); _stencil_val(alpha.x,0,0,0);_stencil_val(fm.x,0,0,0);_stencil_val(f,0,0,0); _stencil_val(f,-1,0,0);  
      }     }}}}end__stencil_is_face_x()
#line 87
_stencil_is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,0,-1,0); _stencil_val(fm.y,0,0,0); {
#line 99 "/home/gongminjiang/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,0,-1,0);
   
#line 104
_stencil_val(phi,0,-1,0); 
#line 103
_stencil_val(phi,0,0,0);
   
#line 103
_stencil_val(phi,0,0,0); 
#line 102
_stencil_val(phi,0,-1,0);_stencil_val(phi,0,0,0); 
#line 101
_stencil_val(phi,0,-1,0);_stencil_val(phi,0,0,0);

 



_stencil_val_r(ia.y,0,0,0); _stencil_val(alpha.y,0,0,0);_stencil_val(fm.y,0,0,0);_stencil_val(f,0,0,0); _stencil_val(f,0,-1,0);  
      }     }}}}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#line 87
if(!is_constant(fm.x) && !is_constant(alpha.x)){{foreach_face_generic(){is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,-1,0,0) && val(fm.x,0,0,0) > 0.) {
#line 99 "/home/gongminjiang/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < HUGE && val(phi,-1,0,0) < HUGE) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < HUGE ? val(phi,0,0,0) :
   val(phi,-1,0,0) < HUGE ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += val(alpha.x,0,0,0)/val(fm.x,0,0,0)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      }}}}end_is_face_x()
#line 87
is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,0,-1,0) && val(fm.y,0,0,0) > 0.) {
#line 99 "/home/gongminjiang/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < HUGE && val(phi,0,-1,0) < HUGE) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < HUGE ? val(phi,0,0,0) :
   val(phi,0,-1,0) < HUGE ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += val(alpha.y,0,0,0)/val(fm.y,0,0,0)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      }}}}end_is_face_y()}end_foreach_face_generic();}}else if(is_constant(fm.x) && !is_constant(alpha.x)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  {
#line 87
foreach_face_generic(){is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,-1,0,0) && _const_fm.x > 0.) {
#line 99 "/home/gongminjiang/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < HUGE && val(phi,-1,0,0) < HUGE) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < HUGE ? val(phi,0,0,0) :
   val(phi,-1,0,0) < HUGE ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += val(alpha.x,0,0,0)/_const_fm.x*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      }}}}end_is_face_x()
#line 87
is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,0,-1,0) && _const_fm.y > 0.) {
#line 99 "/home/gongminjiang/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < HUGE && val(phi,0,-1,0) < HUGE) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < HUGE ? val(phi,0,0,0) :
   val(phi,0,-1,0) < HUGE ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += val(alpha.y,0,0,0)/_const_fm.y*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      }}}}end_is_face_y()}end_foreach_face_generic();}}else if(!is_constant(fm.x) && is_constant(alpha.x)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  {
#line 87
foreach_face_generic(){is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,-1,0,0) && val(fm.x,0,0,0) > 0.) {
#line 99 "/home/gongminjiang/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < HUGE && val(phi,-1,0,0) < HUGE) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < HUGE ? val(phi,0,0,0) :
   val(phi,-1,0,0) < HUGE ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += _const_alpha.x/val(fm.x,0,0,0)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      }}}}end_is_face_x()
#line 87
is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,0,-1,0) && val(fm.y,0,0,0) > 0.) {
#line 99 "/home/gongminjiang/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < HUGE && val(phi,0,-1,0) < HUGE) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < HUGE ? val(phi,0,0,0) :
   val(phi,0,-1,0) < HUGE ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += _const_alpha.y/val(fm.y,0,0,0)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      }}}}end_is_face_y()}end_foreach_face_generic();}}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  {
#line 87
foreach_face_generic(){is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,-1,0,0) && _const_fm.x > 0.) {
#line 99 "/home/gongminjiang/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < HUGE && val(phi,-1,0,0) < HUGE) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < HUGE ? val(phi,0,0,0) :
   val(phi,-1,0,0) < HUGE ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += _const_alpha.x/_const_fm.x*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      }}}}end_is_face_x()
#line 87
is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,0,-1,0) && _const_fm.y > 0.) {
#line 99 "/home/gongminjiang/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < HUGE && val(phi,0,-1,0) < HUGE) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < HUGE ? val(phi,0,0,0) :
   val(phi,0,-1,0) < HUGE ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += _const_alpha.y/_const_fm.y*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      }}}}end_is_face_y()}end_foreach_face_generic();}}






  {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){ {
    _attribute[f.i].prolongation = fraction_refine;
    _attribute[f.i].dirty = true;
  }}}






  {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){ {
    scalar phi = _attribute[f.i].phi;
    delete (((scalar[]){phi,{-1}}));
    _attribute[f.i].phi.i = 0;
  }}}
  pfree (list,__func__,__FILE__,__LINE__);
}{end_tracing("acceleration_0","/home/gongminjiang/basilisk/src/iforce.h",131);return 0;}end_tracing("acceleration_0","/home/gongminjiang/basilisk/src/iforce.h",131);}
#line 16 "/home/gongminjiang/basilisk/src/tension.h"
#line 1 "curvature.h"
#line 1 "/home/gongminjiang/basilisk/src/curvature.h"
#line 12 "/home/gongminjiang/basilisk/src/curvature.h"
static void curvature_restriction (Point point, scalar kappa)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double k = 0., s = 0.;
  {foreach_child()
    if (val(kappa,0,0,0) != HUGE)
      k += val(kappa,0,0,0), s++;end_foreach_child()}
  val(kappa,0,0,0) = s ? k/s : HUGE;
}







static void curvature_prolongation (Point point, scalar kappa)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  {foreach_child() {
    double sk = 0., s = 0.;
    for (int i = 0; i <= 1; i++)

      for (int j = 0; j <= 1; j++)




   if (coarse(kappa,child.x*i,child.y*j,child.z*k) != HUGE)
     sk += coarse(kappa,child.x*i,child.y*j,child.z*k), s++;
    val(kappa,0,0,0) = s ? sk/s : HUGE;
  }end_foreach_child()}
}
#line 66 "/home/gongminjiang/basilisk/src/curvature.h"
#line 1 "heights.h"
#line 1 "/home/gongminjiang/basilisk/src/heights.h"
#line 29 "/home/gongminjiang/basilisk/src/heights.h"
static inline double height (double H) {
  return H > 20./2. ? H - 20. : H < -20./2. ? H + 20. : H;
}

static inline int orientation (double H) {
  return fabs(H) > 20./2.;
}
#line 49 "/home/gongminjiang/basilisk/src/heights.h"
static void half_column (Point point, scalar c, vector h, vector cs, int j)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;






  const int complete = -1;

   {







    double S = val(c,0,0,0), H = S, ci, a;







    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {




      if (val(h.x,0,0,0) == 300.)
 state.s = complete, state.h = HUGE;




      else {
 int s = (val(h.x,0,0,0) + 20./2.)/100.;
 state.h = val(h.x,0,0,0) - 100.*s;
 state.s = s - 1;
      }





      if (state.s != complete)
 S = state.s, H = state.h;
    }
#line 109 "/home/gongminjiang/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? val(c,i*j,0,0) : val(cs.x,(i - 2)*j,0,0);
      H += ci;




      if (S > 0. && S < 1.) {
 S = ci;
 if (ci <= 0. || ci >= 1.) {







   H -= i*ci;
   break;
 }
      }
#line 138 "/home/gongminjiang/basilisk/src/heights.h"
      else if (S >= 1. && ci <= 0.) {
 H = (H - 0.5)*j + (j == -1)*20.;
 S = complete;
 break;
      }
      else if (S <= 0. && ci >= 1.) {
 H = (i + 0.5 - H)*j + (j == 1)*20.;
 S = complete;
 break;
      }
#line 156 "/home/gongminjiang/basilisk/src/heights.h"
      else if (S == ci && modf(H, &a))
 break;
    }





    if (j == -1) {







      if (S != complete && ((val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) ||
       (S > 0. && S < 1.)))
 val(h.x,0,0,0) = 300.;
      else if (S == complete)
 val(h.x,0,0,0) = H;
      else





 val(h.x,0,0,0) = H + 100.*(1. + (S >= 1.));
    }
    else {
#line 195 "/home/gongminjiang/basilisk/src/heights.h"
      if (state.s != complete ||
   (S == complete && fabs(height(H)) < fabs(height(state.h))))
 state.s = S, state.h = H;





      if (state.s != complete)
 val(h.x,0,0,0) = HUGE;
      else
 val(h.x,0,0,0) = (state.h > 1e10 ? HUGE : state.h);
    }
  } 
#line 59
{







    double S = val(c,0,0,0), H = S, ci, a;







    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {




      if (val(h.y,0,0,0) == 300.)
 state.s = complete, state.h = HUGE;




      else {
 int s = (val(h.y,0,0,0) + 20./2.)/100.;
 state.h = val(h.y,0,0,0) - 100.*s;
 state.s = s - 1;
      }





      if (state.s != complete)
 S = state.s, H = state.h;
    }
#line 109 "/home/gongminjiang/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? val(c,0,i*j,0) : val(cs.y,0,(i - 2)*j,0);
      H += ci;




      if (S > 0. && S < 1.) {
 S = ci;
 if (ci <= 0. || ci >= 1.) {







   H -= i*ci;
   break;
 }
      }
#line 138 "/home/gongminjiang/basilisk/src/heights.h"
      else if (S >= 1. && ci <= 0.) {
 H = (H - 0.5)*j + (j == -1)*20.;
 S = complete;
 break;
      }
      else if (S <= 0. && ci >= 1.) {
 H = (i + 0.5 - H)*j + (j == 1)*20.;
 S = complete;
 break;
      }
#line 156 "/home/gongminjiang/basilisk/src/heights.h"
      else if (S == ci && modf(H, &a))
 break;
    }





    if (j == -1) {







      if (S != complete && ((val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) ||
       (S > 0. && S < 1.)))
 val(h.y,0,0,0) = 300.;
      else if (S == complete)
 val(h.y,0,0,0) = H;
      else





 val(h.y,0,0,0) = H + 100.*(1. + (S >= 1.));
    }
    else {
#line 195 "/home/gongminjiang/basilisk/src/heights.h"
      if (state.s != complete ||
   (S == complete && fabs(height(H)) < fabs(height(state.h))))
 state.s = S, state.h = H;





      if (state.s != complete)
 val(h.y,0,0,0) = HUGE;
      else
 val(h.y,0,0,0) = (state.h > 1e10 ? HUGE : state.h);
    }
  }
}
#line 222 "/home/gongminjiang/basilisk/src/heights.h"
static void column_propagation (vector h)
{
  foreach_stencil ()
    for (int i = -2; i <= 2; i++)
      {
 {_stencil_val(h.x,i,0,0);
_stencil_val(h.x,i,0,0);_stencil_val(h.x,0,0,0);
   {_stencil_val_a(h.x,0,0,0); _stencil_val(h.x,i,0,0);   }      
       
#line 229
}
 
#line 227
{_stencil_val(h.y,0,i,0);
_stencil_val(h.y,0,i,0);_stencil_val(h.y,0,0,0);
   {_stencil_val_a(h.y,0,0,0); _stencil_val(h.y,0,i,0);   }      
       
#line 229
}}end_foreach_stencil();
  
#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
{
#line 224
foreach ()
    for (int i = -2; i <= 2; i++)
      {
 if (fabs(height(val(h.x,i,0,0))) <= 3.5 &&
     fabs(height(val(h.x,i,0,0)) + i) < fabs(height(val(h.x,0,0,0))))
   val(h.x,0,0,0) = val(h.x,i,0,0) + i;
 
#line 227
if (fabs(height(val(h.y,0,i,0))) <= 3.5 &&
     fabs(height(val(h.y,0,i,0)) + i) < fabs(height(val(h.y,0,0,0))))
   val(h.y,0,0,0) = val(h.y,0,i,0) + i;}end_foreach();}
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif


#line 230
}
#line 289 "/home/gongminjiang/basilisk/src/heights.h"

static void refine_h_x (Point point, scalar h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;




  bool complete = true;
  {foreach_child() {
    for (int i = -2; i <= 2; i++)
      if (allocated(i,0,0) &&
   !(!is_leaf(neighbor(i,0,0)) && !neighbor(i,0,0).neighbors && neighbor(i,0,0).pid >= 0) && !(neighbor(i,0,0).pid < 0) &&
   fabs(height(val(h,i,0,0))) <= 3.5 &&
   fabs(height(val(h,i,0,0)) + i) < fabs(height(val(h,0,0,0))))
 val(h,0,0,0) = val(h,i,0,0) + i;
    if (val(h,0,0,0) == HUGE)
      complete = false;
  }end_foreach_child()}
  if (complete)
    return;
#line 317 "/home/gongminjiang/basilisk/src/heights.h"
  int ori = orientation(val(h,0,0,0));

  for (int i = -1; i <= 1; i++)
    if (val(h,0,i,0) == HUGE || orientation(val(h,0,i,0)) != ori)
      return;

  double h0 = (30.*height(val(h,0,0,0)) + height(val(h,0,1,0)) + height(val(h,0,-1,0)))/16.
    + 20.*ori;
  double dh = (height(val(h,0,1,0)) - height(val(h,0,-1,0)))/4.;
  {foreach_child()
    if (val(h,0,0,0) == HUGE)
      val(h,0,0,0) = h0 + dh*child.y - child.x/2.;end_foreach_child()}
#line 352 "/home/gongminjiang/basilisk/src/heights.h"
}

#line 290
static void refine_h_y (Point point, scalar h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;




  bool complete = true;
  {foreach_child() {
    for (int i = -2; i <= 2; i++)
      if (allocated(0,i,0) &&
   !(!is_leaf(neighbor(0,i,0)) && !neighbor(0,i,0).neighbors && neighbor(0,i,0).pid >= 0) && !(neighbor(0,i,0).pid < 0) &&
   fabs(height(val(h,0,i,0))) <= 3.5 &&
   fabs(height(val(h,0,i,0)) + i) < fabs(height(val(h,0,0,0))))
 val(h,0,0,0) = val(h,0,i,0) + i;
    if (val(h,0,0,0) == HUGE)
      complete = false;
  }end_foreach_child()}
  if (complete)
    return;
#line 317 "/home/gongminjiang/basilisk/src/heights.h"
  int ori = orientation(val(h,0,0,0));

  for (int i = -1; i <= 1; i++)
    if (val(h,i,0,0) == HUGE || orientation(val(h,i,0,0)) != ori)
      return;

  double h0 = (30.*height(val(h,0,0,0)) + height(val(h,1,0,0)) + height(val(h,-1,0,0)))/16.
    + 20.*ori;
  double dh = (height(val(h,1,0,0)) - height(val(h,-1,0,0)))/4.;
  {foreach_child()
    if (val(h,0,0,0) == HUGE)
      val(h,0,0,0) = h0 + dh*child.x - child.y/2.;end_foreach_child()}
#line 352 "/home/gongminjiang/basilisk/src/heights.h"
}






     
void heights (scalar c, vector h)
{tracing("heights","/home/gongminjiang/basilisk/src/heights.h",360);
  vector  s=new_vector("s");
  
    for (int i = 0; i < nboundary; i++)
      _attribute[s.x.i].boundary[i] = _attribute[c.i].boundary[i];
    
#line 364
for (int i = 0; i < nboundary; i++)
      _attribute[s.y.i].boundary[i] = _attribute[c.i].boundary[i];





  restriction (((scalar[]){c,{-1}}));
  for (int j = -1; j <= 1; j += 2) {





    {foreach_level(0)
      {
        val(h.x,0,0,0) = HUGE;
        
#line 380
val(h.y,0,0,0) = HUGE;}end_foreach_level();}

    for (int l = 1; l <= depth(); l++) {




      {foreach_level (l)
 {
   val(s.x,0,0,0) = val(c,2*j,0,0);
   
#line 389
val(s.y,0,0,0) = val(c,0,2*j,0);}end_foreach_level();}
#line 399 "/home/gongminjiang/basilisk/src/heights.h"
      {foreach_level (l - 1)
 { {
   val(s.x,0,0,0) = val(c,j,0,0);
   val(s.x,j,0,0) = val(c,2*j,0,0);
        } 
#line 400
{
   val(s.y,0,0,0) = val(c,0,j,0);
   val(s.y,0,j,0) = val(c,0,2*j,0);
        }}end_foreach_level();}






      {foreach_halo (prolongation, l - 1)
 {
   _attribute[c.i].prolongation (point, s.x);
   
#line 412
_attribute[c.i].prolongation (point, s.y);}end_foreach_halo();}
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, (scalar *)((vector[]){s,{{-1},{-1}}}), l); };





      {foreach_level (l)
        half_column (point, c, h, s, j);end_foreach_level();}
    }
  }






   {
    _attribute[h.x.i].prolongation = no_data;
    _attribute[h.x.i].restriction = no_restriction;
    _attribute[h.x.i].dirty = true;
  } 
#line 429
{
    _attribute[h.y.i].prolongation = no_data;
    _attribute[h.y.i].restriction = no_restriction;
    _attribute[h.y.i].dirty = true;
  }




  column_propagation (h);






  
    _attribute[h.x.i].prolongation = refine_h_x;
    
#line 446
_attribute[h.y.i].prolongation = refine_h_y;delete((scalar*)((vector[]){s,{{-1},{-1}}}));
end_tracing("heights","/home/gongminjiang/basilisk/src/heights.h",447);}








#line 67 "/home/gongminjiang/basilisk/src/curvature.h"



static double kappa_y (Point point, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  int ori = orientation(val(h.y,0,0,0));
  for (int i = -1; i <= 1; i++)
    if (val(h.y,i,0,0) == HUGE || orientation(val(h.y,i,0,0)) != ori)
      return HUGE;
  double hx = (val(h.y,1,0,0) - val(h.y,-1,0,0))/2.;
  double hxx = (val(h.y,1,0,0) + val(h.y,-1,0,0) - 2.*val(h.y,0,0,0))/Delta;
  return hxx/pow(1. + sq(hx), 3/2.);
}

#line 70
static double kappa_x (Point point, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  int ori = orientation(val(h.x,0,0,0));
  for (int i = -1; i <= 1; i++)
    if (val(h.x,0,i,0) == HUGE || orientation(val(h.x,0,i,0)) != ori)
      return HUGE;
  double hx = (val(h.x,0,1,0) - val(h.x,0,-1,0))/2.;
  double hxx = (val(h.x,0,1,0) + val(h.x,0,-1,0) - 2.*val(h.x,0,0,0))/Delta;
  return hxx/pow(1. + sq(hx), 3/2.);
}


static coord normal_y (Point point, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  coord n = {HUGE, HUGE, HUGE};
  if (val(h.y,0,0,0) == HUGE)
    return n;
  int ori = orientation(val(h.y,0,0,0));
  if (val(h.y,-1,0,0) != HUGE && orientation(val(h.y,-1,0,0)) == ori) {
    if (val(h.y,1,0,0) != HUGE && orientation(val(h.y,1,0,0)) == ori)
      n.x = (val(h.y,-1,0,0) - val(h.y,1,0,0))/2.;
    else
      n.x = val(h.y,-1,0,0) - val(h.y,0,0,0);
  }
  else if (val(h.y,1,0,0) != HUGE && orientation(val(h.y,1,0,0)) == ori)
    n.x = val(h.y,0,0,0) - val(h.y,1,0,0);
  else
    return n;
  double nn = (ori ? -1. : 1.)*sqrt(1. + sq(n.x));
  n.x /= nn;
  n.y = 1./nn;
  return n;
}

#line 82
static coord normal_x (Point point, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  coord n = {HUGE, HUGE, HUGE};
  if (val(h.x,0,0,0) == HUGE)
    return n;
  int ori = orientation(val(h.x,0,0,0));
  if (val(h.x,0,-1,0) != HUGE && orientation(val(h.x,0,-1,0)) == ori) {
    if (val(h.x,0,1,0) != HUGE && orientation(val(h.x,0,1,0)) == ori)
      n.y = (val(h.x,0,-1,0) - val(h.x,0,1,0))/2.;
    else
      n.y = val(h.x,0,-1,0) - val(h.x,0,0,0);
  }
  else if (val(h.x,0,1,0) != HUGE && orientation(val(h.x,0,1,0)) == ori)
    n.y = val(h.x,0,0,0) - val(h.x,0,1,0);
  else
    return n;
  double nn = (ori ? -1. : 1.)*sqrt(1. + sq(n.y));
  n.y /= nn;
  n.x = 1./nn;
  return n;
}
#line 179 "/home/gongminjiang/basilisk/src/curvature.h"
static double height_curvature (Point point, scalar c, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;






  typedef struct {
    double n;
    double (* kappa) (Point, vector);
  } NormKappa;
  struct { NormKappa x, y, z; } n;
  
    n.x.n = val(c,1,0,0) - val(c,-1,0,0), n.x.kappa = kappa_x;
    
#line 193
n.y.n = val(c,0,1,0) - val(c,0,-1,0), n.y.kappa = kappa_y;
  double (* kappaf) (Point, vector) = NULL; NOT_UNUSED (kappaf);




  if (fabs(n.x.n) < fabs(n.y.n))
    do { NormKappa __tmp = n.x; n.x = n.y; n.y = __tmp; } while(0);
#line 211 "/home/gongminjiang/basilisk/src/curvature.h"
  double kappa = HUGE;
  
    if (kappa == HUGE) {
      kappa = n.x.kappa (point, h);
      if (kappa != HUGE) {
 kappaf = n.x.kappa;
 if (n.x.n < 0.)
   kappa = - kappa;
      }
    }
    
#line 213
if (kappa == HUGE) {
      kappa = n.y.kappa (point, h);
      if (kappa != HUGE) {
 kappaf = n.y.kappa;
 if (n.y.n < 0.)
   kappa = - kappa;
      }
    }

  if (kappa != HUGE) {




    if (fabs(kappa) > 1./Delta)
      kappa = sign(kappa)/Delta;
#line 247 "/home/gongminjiang/basilisk/src/curvature.h"
  }

  return kappa;
}
#line 179 "/home/gongminjiang/basilisk/src/curvature.h"
static void _stencil_height_curvature (Point point, scalar c, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;        
      
     
      






  
  
  
    { _stencil_val(c,1,0,0); _stencil_val(c,-1,0,0);     }
    
#line 193
{ _stencil_val(c,0,1,0); _stencil_val(c,0,-1,0);     }                      
  
      




     
#line 211 "/home/gongminjiang/basilisk/src/curvature.h"
     
     
           
      
         
   
 
      
       
     

  
        




       
#line 247 "/home/gongminjiang/basilisk/src/curvature.h"
   

  return ;
}






coord height_normal (Point point, scalar c, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;






  typedef struct {
    double n;
    coord (* normal) (Point, vector);
  } NormNormal;
  struct { NormNormal x, y, z; } n;
  
    n.x.n = val(c,1,0,0) - val(c,-1,0,0), n.x.normal = normal_x;
    
#line 271
n.y.n = val(c,0,1,0) - val(c,0,-1,0), n.y.normal = normal_y;




  if (fabs(n.x.n) < fabs(n.y.n))
    do { NormNormal __tmp = n.x; n.x = n.y; n.y = __tmp; } while(0);
#line 288 "/home/gongminjiang/basilisk/src/curvature.h"
  coord normal = {HUGE, HUGE, HUGE};
  
    if (normal.x == HUGE)
      normal = n.x.normal (point, h);
    
#line 290
if (normal.y == HUGE)
      normal = n.y.normal (point, h);

  return normal;
}
#line 330 "/home/gongminjiang/basilisk/src/curvature.h"
#line 1 "parabola.h"
#line 1 "/home/gongminjiang/basilisk/src/parabola.h"
#line 1 "utils.h"
#line 2 "/home/gongminjiang/basilisk/src/parabola.h"






typedef struct {
  coord o;

  coord m;
  double ** M, rhs[3], a[3];
#line 21 "/home/gongminjiang/basilisk/src/parabola.h"
} ParabolaFit;

static void parabola_fit_init (ParabolaFit * p, coord o, coord m)
{
  
    p->o.x = o.x;
    
#line 26
p->o.y = o.y;

  
    p->m.x = m.x;
    
#line 29
p->m.y = m.y;
  normalize (&p->m);
  int n = 3;
#line 65 "/home/gongminjiang/basilisk/src/parabola.h"
  p->M = (double **) matrix_new (n, n, sizeof(double));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      p->M[i][j] = 0.;
    p->rhs[i] = 0.;
  }
}

static void parabola_fit_add (ParabolaFit * p, coord m, double w)
{

  double x1 = m.x - p->o.x, y1 = m.y - p->o.y;
  double x = p->m.y*x1 - p->m.x*y1;
  double y = p->m.x*x1 + p->m.y*y1;
  double x2 = w*x*x, x3 = x2*x, x4 = x3*x;
  p->M[0][0] += x4;
  p->M[1][0] += x3; p->M[1][1] += x2;
  p->M[2][1] += w*x; p->M[2][2] += w;
  p->rhs[0] += x2*y; p->rhs[1] += w*x*y; p->rhs[2] += w*y;
#line 111 "/home/gongminjiang/basilisk/src/parabola.h"
}

static double parabola_fit_solve (ParabolaFit * p)
{

  p->M[0][1] = p->M[1][0];
  p->M[0][2] = p->M[2][0] = p->M[1][1];
  p->M[1][2] = p->M[2][1];
  double pivmin = matrix_inverse (p->M, 3, 1e-10);
  if (pivmin) {
    p->a[0] = p->M[0][0]*p->rhs[0] + p->M[0][1]*p->rhs[1] + p->M[0][2]*p->rhs[2];
    p->a[1] = p->M[1][0]*p->rhs[0] + p->M[1][1]*p->rhs[1] + p->M[1][2]*p->rhs[2];
  }
  else
    p->a[0] = p->a[1] = 0.;
#line 158 "/home/gongminjiang/basilisk/src/parabola.h"
  matrix_free (p->M);
  return pivmin;
}

static double parabola_fit_curvature (ParabolaFit * p,
          double kappamax, double * kmax)
{
  double kappa;

  double dnm = 1. + sq(p->a[1]);
  kappa = - 2.*p->a[0]/pow(dnm, 3/2.);
  if (kmax)
    *kmax = fabs (kappa);
#line 190 "/home/gongminjiang/basilisk/src/parabola.h"
  if (fabs (kappa) > kappamax) {
    if (kmax)
      *kmax = kappamax;
    return kappa > 0. ? kappamax : - kappamax;
  }
  return kappa;
}
#line 331 "/home/gongminjiang/basilisk/src/curvature.h"






static int independents (coord * p, int n)
{
  if (n < 2)
    return n;
  int ni = 1;
  for (int j = 1; j < n; j++) {
    bool depends = false;
    for (int i = 0; i < j && !depends; i++) {
      double d2 = 0.;
      
 d2 += sq(p[i].x - p[j].x);
 
#line 347
d2 += sq(p[i].y - p[j].y);
      depends = (d2 < sq(0.5));
    }
    ni += !depends;
  }
  return ni;
}






static double height_curvature_fit (Point point, scalar c, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;





  coord ip[2 == 2 ? 6 : 27];
  int n = 0;




   {





    int n1 = 0, n2 = 0;

    for (int i = -1; i <= 1; i++)
      if (val(h.y,i,0,0) != HUGE) {
 if (orientation(val(h.y,i,0,0))) n1++; else n2++;
      }







    int ori = (n1 > n2);







    for (int i = -1; i <= 1; i++)
      if (val(h.y,i,0,0) != HUGE && orientation(val(h.y,i,0,0)) == ori)
 ip[n].x = i, ip[n++].y = height(val(h.y,i,0,0));






  } 
#line 373
{





    int n1 = 0, n2 = 0;

    for (int i = -1; i <= 1; i++)
      if (val(h.x,0,i,0) != HUGE) {
 if (orientation(val(h.x,0,i,0))) n1++; else n2++;
      }







    int ori = (n1 > n2);







    for (int i = -1; i <= 1; i++)
      if (val(h.x,0,i,0) != HUGE && orientation(val(h.x,0,i,0)) == ori)
 ip[n].y = i, ip[n++].x = height(val(h.x,0,i,0));






  }





  if (independents (ip, n) < (2 == 2 ? 3 : 9))
    return HUGE;





  coord m = mycs (point, c), fc;
  double alpha = line_alpha (val(c,0,0,0), m);
  double area = line_length_center(m,alpha,&fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);

  NOT_UNUSED(area);
  parabola_fit_add (&fit, fc, .1);
#line 438 "/home/gongminjiang/basilisk/src/curvature.h"
  for (int i = 0; i < n; i++)
    parabola_fit_add (&fit, ip[i], 1.);
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;



  return kappa;
}







#line 360
static void _stencil_height_curvature_fit (Point point, scalar c, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;          





  
  




   {      





    

    for (int i = -1; i <= 1; i++)
      {_stencil_val(h.y,i,0,0); {
_stencil_val(h.y,i,0,0);  
   
      
#line 384
}   }     







    







    for (int i = -1; i <= 1; i++)
      {_stencil_val(h.y,i,0,0);_stencil_val(h.y,i,0,0);
 {_stencil_val(h.y,i,0,0);     }       }






  } 
#line 373
{      





    

    for (int i = -1; i <= 1; i++)
      {_stencil_val(h.x,0,i,0); {
_stencil_val(h.x,0,i,0);  
   
      
#line 384
}   }     







    







    for (int i = -1; i <= 1; i++)
      {_stencil_val(h.x,0,i,0);_stencil_val(h.x,0,i,0);
 {_stencil_val(h.x,0,i,0);     }       }






  }    
    





             





   _stencil_mycs (point, c);     
  _stencil_val(c,0,0,0);          
  
                 
  

  
  
#line 438 "/home/gongminjiang/basilisk/src/curvature.h"
     
    
  
  



  return ;
}






static double centroids_curvature_fit (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;





  coord m = mycs (point, c), fc;
  double alpha = line_alpha (val(c,0,0,0), m);
  line_length_center(m,alpha,&fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);





  coord r = {x,y,z};
  {foreach_neighbor(1)
    if (val(c,0,0,0) > 0. && val(c,0,0,0) < 1.) {
      coord m = mycs (point, c), fc;
      double alpha = line_alpha (val(c,0,0,0), m);
      double area = line_length_center(m,alpha,&fc);
      coord rn = {x,y,z};
      
 fc.x += (rn.x - r.x)/Delta;
 
#line 478
fc.y += (rn.y - r.y)/Delta;
      parabola_fit_add (&fit, fc, area);
    }end_foreach_neighbor()}
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;



  return kappa;
}







#line 453
static void _stencil_centroids_curvature_fit (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;   





   _stencil_mycs (point, c);     
  _stencil_val(c,0,0,0);    
  
     
  





  
  {foreach_neighbor(1)
    {_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {   
       _stencil_mycs (point, c);     
      _stencil_val(c,0,0,0);      
      
      
       
       
      
    }      }end_foreach_neighbor()}       
  
  



  return ;
}
#line 500 "/home/gongminjiang/basilisk/src/curvature.h"
static inline bool interfacial (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (val(c,0,0,0) >= 1.) {
    for (int i = -1; i <= 1; i += 2)
      {
 if (val(c,i,0,0) <= 0.)
   return true;
 
#line 505
if (val(c,0,i,0) <= 0.)
   return true;}
  }
  else if (val(c,0,0,0) <= 0.) {
    for (int i = -1; i <= 1; i += 2)
      {
 if (val(c,i,0,0) >= 1.)
   return true;
 
#line 511
if (val(c,0,i,0) >= 1.)
   return true;}
  }
  else
    return true;
  return false;
}
#line 500 "/home/gongminjiang/basilisk/src/curvature.h"
static void _stencil_interfacial (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
_stencil_val(c,0,0,0);{ {
    for (int i = -1; i <= 1; i += 2)
      {
 {_stencil_val(c,i,0,0); 
      }
 
#line 505
{_stencil_val(c,0,i,0); 
      }}
  } 
{_stencil_val(c,0,0,0);{ {
    for (int i = -1; i <= 1; i += 2)
      {
 {_stencil_val(c,i,0,0); 
      }
 
#line 511
{_stencil_val(c,0,i,0); 
      }}
  } 
    
}   
  
#line 515
}}
     
  
  
#line 516
return ;
}
#line 530 "/home/gongminjiang/basilisk/src/curvature.h"
typedef struct {
  int h;
  int f;
  int a;
  int c;
} cstats;

struct Curvature {
  scalar c, kappa;
  double sigma;
  bool add;
};

     
cstats curvature (struct Curvature p)
{tracing("curvature","/home/gongminjiang/basilisk/src/curvature.h",544);
  scalar c = p.c, kappa = p.kappa;
  double sigma = p.sigma ? p.sigma : 1.;
  int sh = 0, f = 0, sa = 0, sc = 0;
  vector ch = _attribute[c.i].height,   h=(ch).x.i?(ch):new_vector("h");
  if (!ch.x.i)
    heights (c, h);






  _attribute[kappa.i].refine = _attribute[kappa.i].prolongation = curvature_prolongation;
  _attribute[kappa.i].restriction = curvature_restriction;






  scalar  k=new_scalar("k");
  scalar_clone (k, kappa);

  foreach_stencil() {




_stencil_interfacial (point, c);{
      {_stencil_val_a(k,0,0,0);  } 





{_stencil_val_a(k,0,0,0); _stencil_height_curvature (point, c, h);{
       
{_stencil_val_a(k,0,0,0); _stencil_height_curvature_fit (point, c, h);
          }}    
    
#line 584
}}




     





    
  
#line 585
}end_foreach_stencil();

  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(+:f)reduction(+:sh)){
#line 569
foreach() {




    if (!interfacial (point, c))
      val(k,0,0,0) = HUGE;





    else if ((val(k,0,0,0) = height_curvature (point, c, h)) != HUGE)
      sh++;
    else if ((val(k,0,0,0) = height_curvature_fit (point, c, h)) != HUGE)
      f++;
  }end_foreach();mpi_all_reduce_array(&f,int,MPI_SUM,1);mpi_all_reduce_array(&sh,int,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}

  
#line 587
foreach_stencil () { 





    
_stencil_val(k,0,0,0);{
      { _stencil_val(k,0,0,0); } 
{_stencil_interfacial (point, c);{ {      





      
      {foreach_neighbor(1)
 {_stencil_val(k,0,0,0);
   { _stencil_val(k,0,0,0);  }   }end_foreach_neighbor()}




 


{ _stencil_centroids_curvature_fit (point, c);  }
         
    
      
    
#line 614
}
        
} 
    
#line 616
}}




{
      {_stencil_val_a(kappa,0,0,0);  } 
if (p.add)
      {_stencil_val_r(kappa,0,0,0);  }
    else
      {_stencil_val_a(kappa,0,0,0);  }}
       
    




       
    
  
#line 627
}end_foreach_stencil();

  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(+:sc)reduction(+:sa)){
#line 587
foreach () {





    double kf;
    if (val(k,0,0,0) < HUGE)
      kf = val(k,0,0,0);
    else if (interfacial (point, c)) {





      double sk = 0., a = 0.;
      {foreach_neighbor(1)
 if (val(k,0,0,0) < HUGE)
   sk += val(k,0,0,0), a++;end_foreach_neighbor()}
      if (a > 0.)
 kf = sk/a, sa++;
      else




 kf = centroids_curvature_fit (point, c), sc++;
    }
    else
      kf = HUGE;




    if (kf == HUGE)
      val(kappa,0,0,0) = HUGE;
    else if (p.add)
      val(kappa,0,0,0) += sigma*kf;
    else
      val(kappa,0,0,0) = sigma*kf;
  }end_foreach();mpi_all_reduce_array(&sc,int,MPI_SUM,1);mpi_all_reduce_array(&sa,int,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}

  
#line 629
{ cstats _ret= (cstats){sh, f, sa, sc};{delete((scalar*)((scalar[]){k,{-1}}));if(!(ch).x.i)delete((scalar*)((vector[]){h,{{-1},{-1}}}));}{end_tracing("curvature","/home/gongminjiang/basilisk/src/curvature.h",629);return _ret;}}{delete((scalar*)((scalar[]){k,{-1}}));if(!(ch).x.i)delete((scalar*)((vector[]){h,{{-1},{-1}}}));}
end_tracing("curvature","/home/gongminjiang/basilisk/src/curvature.h",630);}
#line 649 "/home/gongminjiang/basilisk/src/curvature.h"

static double pos_x (Point point, vector h, coord * G, coord * Z)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (fabs(height(val(h.x,0,0,0))) > 1.)
    return HUGE;
  coord o = {x, y, z};
  o.x += height(val(h.x,0,0,0))*Delta;
  double pos = 0.;
  
    pos += (o.x - Z->x)*G->x;
    
#line 658
pos += (o.y - Z->y)*G->y;
  return pos;
}

#line 650
static double pos_y (Point point, vector h, coord * G, coord * Z)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (fabs(height(val(h.y,0,0,0))) > 1.)
    return HUGE;
  coord o = {x, y, z};
  o.y += height(val(h.y,0,0,0))*Delta;
  double pos = 0.;
  
    pos += (o.y - Z->y)*G->y;
    
#line 658
pos += (o.x - Z->x)*G->x;
  return pos;
}







static double height_position (Point point, scalar f, vector h,
          coord * G, coord * Z)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;






  typedef struct {
    double n;
    double (* pos) (Point, vector, coord *, coord *);
  } NormPos;
  struct { NormPos x, y, z; } n;
  
    n.x.n = val(f,1,0,0) - val(f,-1,0,0), n.x.pos = pos_x;
    
#line 683
n.y.n = val(f,0,1,0) - val(f,0,-1,0), n.y.pos = pos_y;




  if (fabs(n.x.n) < fabs(n.y.n))
    do { NormPos __tmp = n.x; n.x = n.y; n.y = __tmp; } while(0);
#line 700 "/home/gongminjiang/basilisk/src/curvature.h"
  double pos = HUGE;
  
    if (pos == HUGE)
      pos = n.x.pos (point, h, G, Z);
    
#line 702
if (pos == HUGE)
      pos = n.y.pos (point, h, G, Z);

  return pos;
}








#line 668
static void _stencil_height_position (Point point, scalar f, vector h,
          coord * G, coord * Z)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;        
          
     
      






  
  
  
    { _stencil_val(f,1,0,0); _stencil_val(f,-1,0,0);     }
    
#line 683
{ _stencil_val(f,0,1,0); _stencil_val(f,0,-1,0);     }                
    




     
#line 700 "/home/gongminjiang/basilisk/src/curvature.h"
  
     
          
      

  return ;
}
#line 717 "/home/gongminjiang/basilisk/src/curvature.h"
struct Position {
  scalar f, pos;
  coord G, Z;
  bool add;
};

void position (struct Position p)
{
  scalar f = p.f, pos = p.pos;
  coord * G = &p.G, * Z = &p.Z;






  _attribute[pos.i].refine = _attribute[pos.i].prolongation = curvature_prolongation;
  _attribute[pos.i].restriction = curvature_restriction;


  vector fh = _attribute[f.i].height,   h=(fh).x.i?(fh):new_vector("h");
  if (!fh.x.i)
    heights (f, h);
  foreach_stencil() {
_stencil_interfacial (point, f);{ {  
       _stencil_height_position (point, f, h, G, Z); 
{      





  _stencil_mycs (point, f);     
 _stencil_val(f,0,0,0); 
 
  
 
         
      }
         
      
#line 756
if (p.add)
 {_stencil_val_r(pos,0,0,0);  }
      else
 {_stencil_val_a(pos,0,0,0);  }
    }
      
{_stencil_val_a(pos,0,0,0);  }}
     
    
  
#line 763
}end_foreach_stencil();
  {
#line 740
foreach() {
    if (interfacial (point, f)) {
      double hp = height_position (point, f, h, G, Z);
      if (hp == HUGE) {





 coord n = mycs (point, f), o = {x,y,z}, c;
 double alpha = line_alpha (val(f,0,0,0), n);
 line_length_center(n,alpha,&c);
 hp = 0.;
 
   hp += (o.x + Delta*c.x - Z->x)*G->x;
   
#line 754
hp += (o.y + Delta*c.y - Z->y)*G->y;
      }
      if (p.add)
 val(pos,0,0,0) += hp;
      else
 val(pos,0,0,0) = hp;
    }
    else
      val(pos,0,0,0) = HUGE;
  }end_foreach();}{if(!(fh).x.i)delete((scalar*)((vector[]){h,{{-1},{-1}}}));}
}
#line 17 "/home/gongminjiang/basilisk/src/tension.h"





#line 36 "/home/gongminjiang/basilisk/src/tension.h"
static int stability_1_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int stability_1(const int i,const double t,Event *_ev){tracing("stability_1","/home/gongminjiang/basilisk/src/tension.h",36);
{





  double amin = HUGE, amax = -HUGE, dmin = HUGE;
  foreach_face_stencil (){_stencil_is_face_x(){
    {_stencil_val(fm.x,0,0,0); {
_stencil_val(alpha.x,0,0,0);_stencil_val(fm.x,0,0,0); { _stencil_val(alpha.x,0,0,0);_stencil_val(fm.x,0,0,0); }
_stencil_val(alpha.x,0,0,0);_stencil_val(fm.x,0,0,0); { _stencil_val(alpha.x,0,0,0);_stencil_val(fm.x,0,0,0); }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_x()
#line 44
_stencil_is_face_y(){
    {_stencil_val(fm.y,0,0,0); {
_stencil_val(alpha.y,0,0,0);_stencil_val(fm.y,0,0,0); { _stencil_val(alpha.y,0,0,0);_stencil_val(fm.y,0,0,0); }
_stencil_val(alpha.y,0,0,0);_stencil_val(fm.y,0,0,0); { _stencil_val(alpha.y,0,0,0);_stencil_val(fm.y,0,0,0); }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_y()}end_foreach_face_stencil();
  
#line 44
if(!is_constant(fm.x) && !is_constant(alpha.x)){
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:dmin) reduction(max:amax)reduction(min:amin)){
#line 44
foreach_face_generic (){is_face_x(){
    if (val(fm.x,0,0,0) > 0.) {
      if (val(alpha.x,0,0,0)/val(fm.x,0,0,0) > amax) amax = val(alpha.x,0,0,0)/val(fm.x,0,0,0);
      if (val(alpha.x,0,0,0)/val(fm.x,0,0,0) < amin) amin = val(alpha.x,0,0,0)/val(fm.x,0,0,0);
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_x()
#line 44
is_face_y(){
    if (val(fm.y,0,0,0) > 0.) {
      if (val(alpha.y,0,0,0)/val(fm.y,0,0,0) > amax) amax = val(alpha.y,0,0,0)/val(fm.y,0,0,0);
      if (val(alpha.y,0,0,0)/val(fm.y,0,0,0) < amin) amin = val(alpha.y,0,0,0)/val(fm.y,0,0,0);
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_y()}end_foreach_face_generic();mpi_all_reduce_array(&dmin,double,MPI_MIN,1);mpi_all_reduce_array(&amax,double,MPI_MAX,1);mpi_all_reduce_array(&amin,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 49
}else if(is_constant(fm.x) && !is_constant(alpha.x)){struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:dmin) reduction(max:amax)reduction(min:amin)){
#line 44
foreach_face_generic (){is_face_x(){
    if (_const_fm.x > 0.) {
      if (val(alpha.x,0,0,0)/_const_fm.x > amax) amax = val(alpha.x,0,0,0)/_const_fm.x;
      if (val(alpha.x,0,0,0)/_const_fm.x < amin) amin = val(alpha.x,0,0,0)/_const_fm.x;
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_x()
#line 44
is_face_y(){
    if (_const_fm.y > 0.) {
      if (val(alpha.y,0,0,0)/_const_fm.y > amax) amax = val(alpha.y,0,0,0)/_const_fm.y;
      if (val(alpha.y,0,0,0)/_const_fm.y < amin) amin = val(alpha.y,0,0,0)/_const_fm.y;
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_y()}end_foreach_face_generic();mpi_all_reduce_array(&dmin,double,MPI_MIN,1);mpi_all_reduce_array(&amax,double,MPI_MAX,1);mpi_all_reduce_array(&amin,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 49
}else if(!is_constant(fm.x) && is_constant(alpha.x)){struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:dmin) reduction(max:amax)reduction(min:amin)){
#line 44
foreach_face_generic (){is_face_x(){
    if (val(fm.x,0,0,0) > 0.) {
      if (_const_alpha.x/val(fm.x,0,0,0) > amax) amax = _const_alpha.x/val(fm.x,0,0,0);
      if (_const_alpha.x/val(fm.x,0,0,0) < amin) amin = _const_alpha.x/val(fm.x,0,0,0);
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_x()
#line 44
is_face_y(){
    if (val(fm.y,0,0,0) > 0.) {
      if (_const_alpha.y/val(fm.y,0,0,0) > amax) amax = _const_alpha.y/val(fm.y,0,0,0);
      if (_const_alpha.y/val(fm.y,0,0,0) < amin) amin = _const_alpha.y/val(fm.y,0,0,0);
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_y()}end_foreach_face_generic();mpi_all_reduce_array(&dmin,double,MPI_MIN,1);mpi_all_reduce_array(&amax,double,MPI_MAX,1);mpi_all_reduce_array(&amin,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 49
}else {struct{double x,y;}_const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);struct{double x,y;}_const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(min:dmin) reduction(max:amax)reduction(min:amin)){
#line 44
foreach_face_generic (){is_face_x(){
    if (_const_fm.x > 0.) {
      if (_const_alpha.x/_const_fm.x > amax) amax = _const_alpha.x/_const_fm.x;
      if (_const_alpha.x/_const_fm.x < amin) amin = _const_alpha.x/_const_fm.x;
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_x()
#line 44
is_face_y(){
    if (_const_fm.y > 0.) {
      if (_const_alpha.y/_const_fm.y > amax) amax = _const_alpha.y/_const_fm.y;
      if (_const_alpha.y/_const_fm.y < amin) amin = _const_alpha.y/_const_fm.y;
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_y()}end_foreach_face_generic();mpi_all_reduce_array(&dmin,double,MPI_MIN,1);mpi_all_reduce_array(&amax,double,MPI_MAX,1);mpi_all_reduce_array(&amin,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 49
}
  double rhom = (1./amin + 1./amax)/2.;





  double sigma = 0.;
  {scalar*_i=(scalar*)( interfaces);if(_i)for(scalar c=*_i;(&c)->i>=0;c=*++_i){
    sigma += _attribute[c.i].sigma;}}
  if (sigma) {
    double dt = sqrt (rhom*cube(dmin)/(3.14159265358979*sigma));
    if (dt < dtmax)
      dtmax = dt;
  }
}{end_tracing("stability_1","/home/gongminjiang/basilisk/src/tension.h",64);return 0;}end_tracing("stability_1","/home/gongminjiang/basilisk/src/tension.h",64);}







static int acceleration_1_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int acceleration_1(const int i,const double t,Event *_ev){tracing("acceleration_1","/home/gongminjiang/basilisk/src/tension.h",72);
{




  {scalar*_i=(scalar*)( interfaces);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
    if (_attribute[f.i].sigma) {





      scalar phi = _attribute[f.i].phi;
      if (phi.i)
 curvature ((struct Curvature){f, phi, _attribute[f.i].sigma, .add = true});
      else {
 phi = new_scalar("phi");
 curvature ((struct Curvature){f, phi, _attribute[f.i].sigma, .add = false});
 _attribute[f.i].phi = phi;
      }
    }}}
}{end_tracing("acceleration_1","/home/gongminjiang/basilisk/src/tension.h",94);return 0;}end_tracing("acceleration_1","/home/gongminjiang/basilisk/src/tension.h",94);}
#line 17 "atomisation.c"
#line 1 "tag.h"
#line 1 "/home/gongminjiang/basilisk/src/tag.h"
#line 13 "/home/gongminjiang/basilisk/src/tag.h"
static void restriction_tag (Point point, scalar t)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double min = HUGE;
  {foreach_child()
    if (val(t,0,0,0) < min)
      min = val(t,0,0,0);end_foreach_child()}
  val(t,0,0,0) = min;
}







static long lookup_tag (Array * a, double tag)
{
  long len = a->len/sizeof(double);
  double * p = (double *) a->p;
  if (tag == p[0])
    return 0;
  if (tag == p[len - 1])
    return len - 1;

  long s = 0, e = len - 1;
  while (s < e - 1) {
    long m = (s + e)/2;
    if (p[m] <= tag)
      s = m;
    else
      e = m;
  }
  return s;
}


static int compar_double (const void * p1, const void * p2)
{
  const double * a = p1, * b = p2;
  return *a > *b;
}







     
int tag (scalar t)
{tracing("tag","/home/gongminjiang/basilisk/src/tag.h",62);




  _attribute[t.i].restriction = restriction_tag;

  _attribute[t.i].refine = _attribute[t.i].prolongation = refine_injection;
  _attribute[t.i].dirty = true;
#line 82 "/home/gongminjiang/basilisk/src/tag.h"
  scalar  index=new_scalar("index");
  z_indexing (index, true);
  foreach_stencil()
    {_stencil_val_a(t,0,0,0);_stencil_val(t,0,0,0);_stencil_val(index,0,0,0);      }end_foreach_stencil();
  {
#line 84
foreach()
    val(t,0,0,0) = (val(t,0,0,0) != 0)*(val(index,0,0,0) + 1);end_foreach();}
#line 99 "/home/gongminjiang/basilisk/src/tag.h"
  bool changed;
  do {






    restriction (((scalar[]){t,{-1}}));





    changed = false;
    for (int l = 1; l <= grid->maxdepth; l++) {






      {foreach_level(l)
 if (coarse(t,0,0,0))
   val(t,0,0,0) = coarse(t,0,0,0);end_foreach_level();}
      boundary_level (((scalar[]){t,{-1}}), l);







      
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel  reduction(||:changed)){
#line 132
foreach_level (l)
        if (val(t,0,0,0)) {
   double min = val(t,0,0,0);
   {foreach_neighbor(1)
     if (val(t,0,0,0) && val(t,0,0,0) < min)
       min = val(t,0,0,0);end_foreach_neighbor()}






   
     for (int i = -1; i <= 2; i += 3)
       if ((!is_leaf (neighbor((2*i - 1)/3,0,0)) && neighbor((2*i - 1)/3,0,0).neighbors && neighbor((2*i - 1)/3,0,0).pid >= 0))
  for (int j = 0; j <= 1; j++)
    for (int k = 0; k <= 1; k++)
      if (fine(t,i,j,k) && fine(t,i,j,k) < min)
        min = fine(t,i,j,k);
     
#line 145
for (int i = -1; i <= 2; i += 3)
       if ((!is_leaf (neighbor(0,(2*i - 1)/3,0)) && neighbor(0,(2*i - 1)/3,0).neighbors && neighbor(0,(2*i - 1)/3,0).pid >= 0))
  for (int j = 0; j <= 1; j++)
    for (int k = 0; k <= 1; k++)
      if (fine(t,j,i,k) && fine(t,j,i,k) < min)
        min = fine(t,j,i,k);


   if (val(t,0,0,0) != min) {
     changed = true;
     val(t,0,0,0) = min;
   }
 }end_foreach_level();mpi_all_reduce_array(&changed,bool,MPI_LOR,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
      
#line 158
boundary_level (((scalar[]){t,{-1}}), l);
    }
  } while (changed);
#line 171 "/home/gongminjiang/basilisk/src/tag.h"
  Array * a = array_new();
  foreach_stencil ()
    {_stencil_val(t,0,0,0); {        







       _stencil_val(t,0,0,0);   
           
      
     
   
       
  
    
               
          
                






      
  
    
       
         
    





}   }end_foreach_stencil();
  
#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
{
#line 172
foreach ()
    if (val(t,0,0,0) > 0) {







      double tag = val(t,0,0,0), * ap = (double *) a->p;
      long s = -1;
      if (a->len == 0 || tag > ap[a->len/sizeof(double) - 1])
 s = a->len/sizeof(double);
      else if (tag < ap[0])
 s = 0;
      else {






 s = lookup_tag (a, tag) + 1;
 if (tag == ap[s - 1] || tag == ap[s])
   s = -1;
      }
      if (s >= 0) {





 array_append (a, &tag, sizeof(double)), ap = (double *) a->p;
 for (int i = a->len/sizeof(double) - 1; i > s; i--)
   ap[i] = ap[i-1];
 ap[s] = tag;
      }
    }end_foreach();}
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif

#line 222 "/home/gongminjiang/basilisk/src/tag.h"
  long lmax = a->len;
  mpi_all_reduce (lmax, MPI_LONG, MPI_MAX);
  a->p = prealloc (a->p, lmax,__func__,__FILE__,__LINE__);
  lmax /= sizeof(double);






  double * q = a->p;
  for (int i = a->len/sizeof(double); i < lmax; i++)
    q[i] = -1;
  double p[lmax*npe()];
  MPI_Allgather (a->p, lmax, MPI_DOUBLE, p, lmax, MPI_DOUBLE, MPI_COMM_WORLD);
  qsort (p, lmax*npe(), sizeof(double), compar_double);







  array_free (a);
  a = array_new();
  double last = -1;
  for (int i = 0; i < lmax*npe(); i++)
    if (p[i] != last) {
      array_append (a, &p[i], sizeof(double));
      last = p[i];
    }






  foreach_stencil()
    {_stencil_val(t,0,0,0);
      {_stencil_val_a(t,0,0,0); _stencil_val(t,0,0,0);     }   }end_foreach_stencil();






  {
#line 259
foreach()
    if (val(t,0,0,0) > 0)
      val(t,0,0,0) = lookup_tag (a, val(t,0,0,0)) + 1;end_foreach();}




  int n = a->len/sizeof(double);
  array_free (a);
  {delete((scalar*)((scalar[]){index,{-1}}));{end_tracing("tag","/home/gongminjiang/basilisk/src/tag.h",268);return n;}}delete((scalar*)((scalar[]){index,{-1}}));
end_tracing("tag","/home/gongminjiang/basilisk/src/tag.h",269);}
#line 278 "/home/gongminjiang/basilisk/src/tag.h"
struct RemoveDroplets {
  scalar f;
  int minsize;
  double threshold;
  bool bubbles;
};

void remove_droplets (struct RemoveDroplets p)
{
  scalar  d=new_scalar("d"), f = p.f;
  double threshold = p.threshold ? p.threshold : 1e-4;
  foreach_stencil()
    {_stencil_val_a(d,0,0,0); _stencil_val(f,0,0,0); _stencil_val(f,0,0,0);        }end_foreach_stencil();
  {
#line 289
foreach()
    val(d,0,0,0) = (p.bubbles ? 1. - val(f,0,0,0) : val(f,0,0,0)) > threshold;end_foreach();}
  int n = tag (d), size[n];
  for (int i = 0; i < n; i++)
    size[i] = 0;
  foreach_stencil ()
    {_stencil_val(d,0,0,0);
      { _stencil_val(d,0,0,0);  }   }end_foreach_stencil();
  
#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
{
#line 294
foreach ()
    if (val(d,0,0,0) > 0)
      size[((int) val(d,0,0,0)) - 1]++;end_foreach();}
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif


  
#line 298
MPI_Allreduce (MPI_IN_PLACE, size, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  int minsize = pow (p.minsize ? p.minsize : 3, 2);
  foreach_stencil()
    {_stencil_val(d,0,0,0); _stencil_val(d,0,0,0);
      {_stencil_val_a(f,0,0,0);  }         }end_foreach_stencil();
  {
#line 301
foreach()
    if (val(d,0,0,0) > 0 && size[((int) val(d,0,0,0)) - 1] < minsize)
      val(f,0,0,0) = p.bubbles;end_foreach();}delete((scalar*)((scalar[]){d,{-1}}));
}
#line 18 "atomisation.c"
#line 1 "view.h"
#line 1 "/home/gongminjiang/basilisk/src/view.h"
#line 67 "/home/gongminjiang/basilisk/src/view.h"
#line 1 "gl/framebuffer.h"
#line 1 "/home/gongminjiang/basilisk/src/gl/framebuffer.h"
typedef struct _framebuffer framebuffer;
typedef unsigned int fbdepth_t;

framebuffer * framebuffer_new (unsigned width, unsigned height);
void framebuffer_destroy (framebuffer * p);
unsigned char * framebuffer_image (framebuffer * p);
fbdepth_t * framebuffer_depth (framebuffer * p);
#line 68 "/home/gongminjiang/basilisk/src/view.h"
#line 1 "gl/trackball.h"
#line 1 "/home/gongminjiang/basilisk/src/gl/trackball.h"
#line 50 "/home/gongminjiang/basilisk/src/gl/trackball.h"
void
gl_trackball(float q[4], float p1x, float p1y, float p2x, float p2y);
#line 61 "/home/gongminjiang/basilisk/src/gl/trackball.h"
void
gl_add_quats(float *q1, float *q2, float *dest);





void
gl_build_rotmatrix(float m[4][4], float q[4]);






void
gl_axis_to_quat(float a[3], float phi, float q[4]);
#line 69 "/home/gongminjiang/basilisk/src/view.h"
#line 1 "gl/utils.h"
#line 1 "/home/gongminjiang/basilisk/src/gl/utils.h"
#line 1 "./gl/gl.h"
#line 1 "/home/gongminjiang/basilisk/src/gl/gl.h"
#line 120 "/home/gongminjiang/basilisk/src/gl/gl.h"
typedef unsigned int GLenum;
typedef int GLint;
typedef unsigned int GLuint;
typedef float GLfloat;
typedef double GLdouble;
typedef int GLsizei;
typedef unsigned int GLbitfield;
typedef unsigned char GLubyte;

void glBegin (GLenum mode);
void glEnd (void);
void glClear (GLbitfield mask);
void glClearColor (float red, float green, float blue, float alpha);
void glBindTexture (GLenum target, GLuint texture);
void glColor3f (GLfloat red, GLfloat green, GLfloat blue);
void glColorMaterial (GLenum face, GLenum mode);
void glDisable (GLenum cap);
void glEnable (GLenum cap);
void glFinish (void);
void glGetDoublev (GLenum pname, GLdouble * params);
void glGetIntegerv (GLenum pname, GLint * params);
void glHint (GLenum target, GLenum mode);
void glLightfv (GLenum light, GLenum pname, const GLfloat *params);
void glGetLightfv (GLenum light, GLenum pname, GLfloat * params);
void glLightModeli (GLenum pname, GLint param);
void glLineWidth (GLfloat width);
void glNormal3d (GLdouble nx, GLdouble ny, GLdouble nz);
void glOrtho (GLdouble left, GLdouble right, GLdouble bottom, GLdouble top,
       GLdouble nearVal, GLdouble farVal);
void glShadeModel (GLenum mode);
void glTexCoord1d (GLdouble s);
void glTexCoord2f (GLfloat s, GLfloat t);
void glTexImage1D (GLenum target, GLint level, GLint internalFormat,
     GLsizei width, GLint border, GLenum format, GLenum type,
     const void * data);
void glTexParameteri (GLenum target, GLenum pname, GLint param);
void glVertex3d (GLdouble x, GLdouble y, GLdouble z);
void glVertex3f (GLfloat x, GLfloat y, GLfloat z);

GLenum glGetError (void);
void glGetFloatv (GLenum pname, GLfloat * params);
void glMultMatrixf (const GLfloat * m);
void glLoadIdentity (void);
void glScalef (GLfloat x, GLfloat y, GLfloat z);
void glTranslatef (GLfloat x, GLfloat y, GLfloat z);
void glRotatef (GLfloat angle, GLfloat x, GLfloat y, GLfloat z);
void glMatrixMode (GLenum mode);
void glPopMatrix (void);
void glPushMatrix (void);
void glLoadMatrixd (const GLdouble * m);
#line 2 "/home/gongminjiang/basilisk/src/gl/utils.h"
#line 1 "/home/gongminjiang/basilisk/src/ast/std/stdio.h"
#include <stdio.h>
#line 3 "/home/gongminjiang/basilisk/src/gl/utils.h"
#line 1 "/home/gongminjiang/basilisk/src/ast/std/stdbool.h"
#include <stdbool.h>
#line 4 "/home/gongminjiang/basilisk/src/gl/utils.h"

void gl_write_image (FILE * fp, const GLubyte * buffer,
       unsigned width, unsigned height, unsigned samples);
void init_gl();

void matrix_multiply (float * m, const float * n);
void vector_multiply (float * v, const float * m);

typedef struct {
  float m[16], p[16];
  float n[6][3];
  float d[6];
  unsigned width;
} Frustum;

void gl_get_frustum (Frustum * f);
int sphere_in_frustum (double x, double y, double z, double r, Frustum * f);
float sphere_diameter (double x, double y, double z, double r, Frustum * f);
void gl_check_error();

int polygonize (const double val[8], double isolevel, double triangles[5][3][3]);
void gl_perspective (double fovy, double aspect, double zNear, double zFar);
int gl_project (float objx, float objy, float objz,
  const float modelMatrix[16],
  const float projMatrix[16],
  const GLint viewport[4],
  float *winx, float *winy, float *winz);

#line 1 "gl/parser.h"
#line 1 "/home/gongminjiang/basilisk/src/gl/parser.h"
typedef struct _Node Node;

struct _Node {
  char type;
  union {
    char * id;
    double (* func) (double);
    double value;
  } d;
  int s;
  Node * e[3];
};

Node * parse_node (char * code);
void free_node (Node * n);
void print_node (Node * n, FILE * fp);
void reset_node_type (Node * n, char type);
#line 33 "/home/gongminjiang/basilisk/src/gl/utils.h"
#line 70 "/home/gongminjiang/basilisk/src/view.h"



#line 1 "input.h"
#line 1 "/home/gongminjiang/basilisk/src/input.h"
#line 16 "/home/gongminjiang/basilisk/src/input.h"
struct InputPGM {

  scalar s;
  FILE * fp;

  double ox, oy, width;
};

void input_pgm (struct InputPGM p)
{
  scalar s = p.s;
  if (p.width == 0.) p.width = L0;

  char line[81];
  if (!fgets (line, 81, p.fp)) {
    fprintf (ferr, "input_pgm: could not read magic number\n");
    exit (1);
  }
  if (strcmp (line, "P2\n") && strcmp (line, "P5\n")) {
    fprintf (ferr, "input_pgm: magic number '%s' does not match PGM\n",
      line);
    exit (1);
  }
  int binary = !strcmp (line, "P5\n");
  if (!fgets (line, 81, p.fp)) {
    fprintf (ferr, "input_pgm: could not read width and height\n");
    exit (1);
  }
  int width, height;
  while (line[0] == '#' && fgets (line, 81, p.fp));
  if (line[0] == '#' || sscanf (line, "%d %d", &width, &height) != 2) {
    fprintf (ferr, "input_pgm: could not read width and height\n");
    exit (1);
  }
  if (!fgets (line, 81, p.fp)) {
    fprintf (ferr, "input_pgm: could not read maxval\n");
    exit (1);
  }
  int maxval;
  if (sscanf (line, "%d", &maxval) != 1) {
    fprintf (ferr, "input_pgm: could not read maxval\n");
    exit (1);
  }
  if (maxval < 256) {
    unsigned char * a = ((unsigned char *) pmalloc ((width*height)*sizeof(unsigned char),__func__,__FILE__,__LINE__));
    size_t n = 0;
    if (binary)
      n = fread (a, 1, width*height, p.fp);
    else {
      int v;
      while (n < width*height && fscanf (p.fp, "%d ", &v) == 1)
 a[n++] = v;
    }
    if (n != width*height) {
      fprintf (ferr, "input_pgm: read only %ld values\n", n);
      exit (1);
    }
    foreach_stencil() {          
      
{
 {_stencil_val_a(s,0,0,0);          }
 
{_stencil_val_a(s,0,0,0);  }}
                     
      
    
#line 79
}end_foreach_stencil();
    {
#line 73
foreach() {
      int i = (x - p.ox)*width/p.width, j = (y - p.oy)*width/p.width;
      if (i >= 0 && i < width && j >= 0 && j < height)
 val(s,0,0,0) = 1. - a[(height - 1 - j)*width + i]/(double)maxval;
      else
 val(s,0,0,0) = 0.;
    }end_foreach();}
    pfree (a,__func__,__FILE__,__LINE__);
  }
  else {
    unsigned short * a = ((unsigned short *) pmalloc ((width*height)*sizeof(unsigned short),__func__,__FILE__,__LINE__));
    size_t n = 0;
    if (binary)
      n = fread (a, 2, width*height, p.fp);
    else {
      int v;
      while (n < width*height && fscanf (p.fp, "%d ", &v) == 1)
 a[n++] = v;
    }
    if (n != width*height) {
      fprintf (ferr, "input_pgm: read only %ld values\n", n);
      exit (1);
    }
    foreach_stencil() {          
      
{
 {_stencil_val_a(s,0,0,0);          }
 
{_stencil_val_a(s,0,0,0);  }}
                     
      
    
#line 102
}end_foreach_stencil();
    {
#line 96
foreach() {
      int i = (x - p.ox)*width/p.width, j = (y - p.oy)*width/p.width;
      if (i >= 0 && i < width && j >= 0 && j < height)
 val(s,0,0,0) = 1. - a[(height - 1 - j)*width + i]/(double)maxval;
      else
 val(s,0,0,0) = 0.;
    }end_foreach();}
    pfree (a,__func__,__FILE__,__LINE__);
  }
}

static void next_char (FILE * fp, int target)
{
  int c = fgetc(fp), para = 0;
  while (c != EOF && (c != target || para > 0)) {
    if (c == '{') para++;
    if (c == '}') para--;
    c = fgetc(fp);
  }
  if (c != target) {
    fprintf (ferr, "input_gfs(): error: expecting '%c'\n", target);
    exit (1);
  }
}

static int next_string (FILE * fp, const char * target)
{
  int slen = strlen (target), para = 0;
  char s[slen + 1];
  s[slen] = '\0';
  int len = 0, c = fgetc (fp);
  while (c != EOF && len < slen) {
    if (c == '{') para++;
    if (c == '}') para--;
    s[len++] = c;
    c = fgetc (fp);
  }
  while (c != EOF && para >= 0) {
    if (!strcmp (s, target) && para == 0)
      break;
    if (c == '{') para++;
    if (c == '}') para--;
    for (int i = 0; i < slen - 1; i++)
      s[i] = s[i+1];
    s[slen - 1] = c;
    c = fgetc (fp);
  }
  if (strcmp (s, target))
    c = -1;
  return c;
}
#line 166 "/home/gongminjiang/basilisk/src/input.h"
     
void input_gfs (struct OutputGfs p)
{tracing("input_gfs","/home/gongminjiang/basilisk/src/input.h",167);
  not_mpi_compatible();

  bool opened = false;
  if (p.fp == NULL) {
    if (p.file == NULL)
      p.fp = stdin;
    else if (!(p.fp = fopen (p.file, "r"))) {
      perror (p.file);
      exit (1);
    }
    else
      opened = true;
  }
  bool input_all = (p.list == all);
  if (p.list == NULL) p.list = all;


  init_grid (1);


  next_char (p.fp, '{');

  char * s = ((char *) pmalloc ((1)*sizeof(char),__func__,__FILE__,__LINE__));
  int len = 0;
  int c = fgetc(p.fp);
  while (c != EOF && c != '}') {
    s[len++] = c;
    s = (char *) prealloc (s, (len + 1)*sizeof(char),__func__,__FILE__,__LINE__);
    s[len] = '\0';
    c = fgetc(p.fp);
  }
  if (c != '}') {
    fprintf (ferr, "input_gfs(): error: expecting '}'\n");
    exit (1);
  }

  char * s1 = strstr (s, "variables");
  if (!s1) {
    fprintf (ferr, "input_gfs(): error: expecting 'variables'\n");
    exit (1);
  }

  s1 = strstr (s1, "=");
  if (!s1) {
    fprintf (ferr, "input_gfs(): error: expecting '='\n");
    exit (1);
  }
  s1++;

  while (strchr (" \t", *s1))
    s1++;

  scalar * input = NULL;
  s1 = strtok (s1, ", \t");
  while (s1) {
    char * name = replace (s1, '_', '.', false);
    bool found = false;
    {scalar*_i=(scalar*)( p.list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (!is_constant(s) && _attribute[s.i].name && !strcmp (_attribute[s.i].name, name)) {
 input = list_append (input, s);
 found = true; break;
      }}}
    if (!found) {
      if (input_all) {
 scalar s = new_scalar("s");
 pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
 _attribute[s.i].name = pstrdup (name,__func__,__FILE__,__LINE__);
 input = list_append (input, s);
      }
      else
 input = list_append (input, (scalar){INT_MAX});
    }
    pfree (name,__func__,__FILE__,__LINE__);
    s1 = strtok (NULL, ", \t");
  }
  pfree (s,__func__,__FILE__,__LINE__);

  next_char (p.fp, '{');
  double t1 = 0.;
  if (next_string (p.fp, "Time") >= 0) {
    next_char (p.fp, '{');
    next_char (p.fp, 't');
    next_char (p.fp, '=');
    if (fscanf (p.fp, "%lf", &t1) != 1) {
      fprintf (ferr, "input_gfs(): error: expecting 't'\n");
      exit (1);
    }
    next_char (p.fp, '}');
    next_char (p.fp, '}');
  }

  if (next_string (p.fp, "Box") < 0) {
    fprintf (ferr, "input_gfs(): error: expecting 'GfsBox'\n");
    exit (1);
  }

  next_char (p.fp, '{');
  next_char (p.fp, '{');
  next_char (p.fp, '\n');

  scalar * listm = ((scalar[]){cm,fm.x,fm.y,{-1}});
  scalar * listr = !is_constant(cm) ? listm : NULL;
  NOT_UNUSED (listr);

  {foreach_cell() {
    unsigned flags;
    if (fread (&flags, sizeof (unsigned), 1, p.fp) != 1) {
      fprintf (ferr, "input_gfs(): error: expecting 'flags'\n");
      exit (1);
    }
    if (!(flags & (1 << 4)) && is_leaf(cell))
      refine_cell (point, listr, 0, NULL);
    double a;
    if (fread (&a, sizeof (double), 1, p.fp) != 1 || a != -1) {
      fprintf (ferr, "input_gfs(): error: expecting '-1'\n");
      exit (1);
    }
    {scalar*_i=(scalar*)( input);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      if (fread (&a, sizeof (double), 1, p.fp) != 1) {
 fprintf (ferr, "input_gfs(): error: expecting a scalar\n");
 exit (1);
      }
      if (s.i != INT_MAX) {
 if (_attribute[s.i].v.x.i >= 0) {



   if (_attribute[s.i].v.x.i == s.i) {
     s = _attribute[s.i].v.y;
     val(s,0,0,0) = a;
   }
   else if (_attribute[s.i].v.y.i == s.i) {
     s = _attribute[s.i].v.x;
     val(s,0,0,0) = - a;
   }





 }
 else
   val(s,0,0,0) = a;
      }
    }}}
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}
  {scalar*_i=(scalar*)( listm);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s))
      _attribute[s.i].dirty = true;}}
  {scalar*_i=(scalar*)( input);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s))
      _attribute[s.i].dirty = true;}}

  pfree (input,__func__,__FILE__,__LINE__);
  if (opened)
    fclose (p.fp);


  while (t < t1 && events (false))
    t = tnext;
  events (false);
end_tracing("input_gfs","/home/gongminjiang/basilisk/src/input.h",332);}
#line 371 "/home/gongminjiang/basilisk/src/input.h"
struct InputGRD {
  scalar s;
  FILE * fp;
  char * file;
  double nodatavalue;
  bool linear, periodic, zero;
  int smooth;
};

void input_grd (struct InputGRD p)
{
  scalar input = p.s;

  bool opened = false;
  if (p.fp == NULL) {
    if (p.file == NULL)
      p.fp = stdin;
    else if (!(p.fp = fopen (p.file, "r"))) {
      perror (p.file);
      exit (1);
    }
    else
      opened = true;
  }


  double DeltaGRD;
  int nx, ny;
  double XG0, YG0, ndv;


  char waste[100];
  if (fscanf (p.fp, "%s %d", waste, &nx) != 2) {
    fprintf (ferr, "input_grd(): error reading 'nx'\n");
    if (opened) fclose (p.fp);
    return;
  }
  if (fscanf (p.fp, "%s %d", waste, &ny) != 2) {
    fprintf (ferr, "input_grd(): error reading 'ny'\n");
    if (opened) fclose (p.fp);
    return;
  }
  if (fscanf (p.fp, "%s %lf", waste, &XG0) != 2) {
    fprintf (ferr, "input_grd(): error reading 'XG0'\n");
    if (opened) fclose (p.fp);
    return;
  }
  if (fscanf (p.fp, "%s %lf", waste, &YG0) != 2) {
    fprintf (ferr, "input_grd(): error reading 'YG0'\n");
    if (opened) fclose (p.fp);
    return;
  }
  if (fscanf (p.fp, "%s %lf", waste, &DeltaGRD) != 2) {
    fprintf (ferr, "input_grd(): error reading 'DeltaGRD'\n");
    if (opened) fclose (p.fp);
    return;
  }
  if (fscanf (p.fp, "%s %lf", waste, &ndv) != 2) {
    fprintf (ferr, "input_grd(): error reading 'ndv'\n");
    if (opened) fclose (p.fp);
    return;
  }


  if (!p.nodatavalue)
    p.nodatavalue = ndv;


  double * value = ((double *) pmalloc ((nx*ny)*sizeof(double),__func__,__FILE__,__LINE__));
  for (int i = ny - 1; i >= 0; i--)
    for (int j = 0 ; j < nx; j++) {
      if (fscanf (p.fp, "%lf ", &value[j + i*nx]) != 1) {
 fprintf (ferr, "input_grd(): error reading value %d,%d\n", i, j);
 if (opened) fclose (p.fp);
 pfree (value,__func__,__FILE__,__LINE__);
 return;
      }
      if (p.zero && value[j + i*nx] == ndv)
 value[j + i*nx] = 0.;
    }


  if (p.smooth > 0) {
    double * smoothed = ((double *) pmalloc ((nx*ny)*sizeof(double),__func__,__FILE__,__LINE__));
    for (int s = 0; s < p.smooth; s++) {
      for (int i = 0; i < ny; i++)
 for (int j = 0 ; j < nx; j++) {
   int n = 0;
   smoothed[j + i*nx] = 0.;
   for (int k = -1; k <= 1; k++)
     for (int l = -1; l <= 1; l++)
       if ((l != 0 || k != 0) &&
    i + k >= 0 && i + k < ny &&
    j + l >= 0 && j + l < nx &&
    value[j + l + (i + k)*nx] != ndv)
  smoothed[j + i*nx] += value[j + l + (i + k)*nx], n++;
   if (n == 0)
     smoothed[j + i*nx] = p.zero ? 0. : ndv;
   else
     smoothed[j + i*nx] /= n;
 }
      do { double * __tmp = value; value = smoothed; smoothed = __tmp; } while(0);
    }
    pfree (smoothed,__func__,__FILE__,__LINE__);
  }

  bool warning = false;
  foreach_stencil () {                   
     
  
           
           
  
     

    
    
{ {              

 
  
          
                 
           
                  
                   
         
                 
      

      
       





{
 {_stencil_val_a(input,0,0,0);  }
 
{_stencil_val_a(input,0,0,0);  }}
                               
                  
              
      
     
         
      
    
#line 512
} 
{
      _stencil_val_a(input,0,0,0);    
      
    }}
                   
    
  
#line 517
}end_foreach_stencil();
  
#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
{
#line 478
foreach () {
    if (p.periodic || _attribute[input.i].boundary[right] == periodic_bc) {
      if (x > XG0 + nx*DeltaGRD)
 x -= nx*DeltaGRD;
      else if (x < XG0)
 x += nx*DeltaGRD;
    }

    int j = (x - XG0 + DeltaGRD/2.)/DeltaGRD;
    int i = (y - YG0 + DeltaGRD/2.)/DeltaGRD;
    if (i >= 0 && i < ny && j >= 0 && j < nx) {
      double val;

      int j1 = (x - XG0)/DeltaGRD;
      int i1 = (y - YG0)/DeltaGRD;
      if (p.linear && i1 >= 0 && j1 >= 0 && i1 < ny - 1 && j1 < nx - 1 &&
   value[j1 + i1*nx] != ndv && value[j1 + 1 + i1*nx] != ndv &&
   value[j1 + (i1 + 1)*nx] != ndv && value[j1 + 1 + (i1 + 1)*nx] != ndv) {

 double dx = x - (j1*DeltaGRD + XG0);
 double dy = y - (i1*DeltaGRD + YG0);
 val = (value[j1 + i1*nx] +
        dx*(value[j1 + 1 + i1*nx] - value[j1 + i1*nx])/DeltaGRD +
        dy*(value[j1 + (i1 + 1)*nx] - value[j1 + i1*nx])/DeltaGRD +
        dx*dy*(value[j1 + i1*nx] + value[j1 + 1 + (i1 + 1)*nx] -
        value[j1 + (i1 + 1)*nx] - value[j1 + 1 + i1*nx])
        /sq(DeltaGRD));
      }
      else
 val = value[j + i*nx];
      if (val == ndv)
 val(input,0,0,0) = HUGE;
      else
 val(input,0,0,0) = val;
    }
    else {
      val(input,0,0,0) = HUGE;
      warning = true;
    }
  }end_foreach();}
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif

  
#line 518
pfree (value,__func__,__FILE__,__LINE__);

  if (warning)
    fprintf (ferr,
      "input_grd(): Warning: Raster data is not covering all"
      " the simulation area\n");

  if (opened)
    fclose (p.fp);
}
#line 74 "/home/gongminjiang/basilisk/src/view.h"







typedef struct {
  char * expr;
  scalar s;
} cexpr;

static scalar get_cexpr (cexpr * cache, const char * expr)
{
  cexpr * c = cache;
  while (c->expr) {
    if (!strcmp (c->expr, expr)) {


      cexpr tmp = *c;
      while ((c + 1)->expr)
 *c = *(c + 1), c++;
      *c = tmp;
      return c->s;
    }
    c++;
  }
  return (scalar){-1};
}

static cexpr * add_cexpr (cexpr * cache, int maxlen,
     const char * expr, scalar s)
{
  cexpr * c = cache;
  while (c->expr) c++;
  int len = c - cache;
  if (len < maxlen) {
    cache = prealloc (cache, sizeof(cexpr)*(len + 2),__func__,__FILE__,__LINE__);
    c = &cache[len];
  }
  else {

    c = cache;
    pfree (c->expr,__func__,__FILE__,__LINE__);
    scalar s = c->s;
    delete (((scalar[]){s,{-1}}));

    while ((c + 1)->expr)
      *c = *(c + 1), c++;
  }
  c->expr = pstrdup (expr,__func__,__FILE__,__LINE__);
  c->s = s;
  (c + 1)->expr = NULL;
  return cache;
}

static void free_cexpr (cexpr * cache)
{
  cexpr * c = cache;
  while (c->expr) {
    pfree (c->expr,__func__,__FILE__,__LINE__);
    scalar s = c->s;
    delete (((scalar[]){s,{-1}}));
    c++;
  }
  pfree (cache,__func__,__FILE__,__LINE__);
}






struct _bview {
  float tx, ty, sx, sy, sz;
  float quat[4];
  float fov;
  float tz, near, far;

  bool gfsview;
  bool reversed;

  float bg[3];
  float lc;
  float res;

  unsigned width, height, samples;

  framebuffer * fb;
  Frustum frustum;

  void (* map) (coord *);

  int ni;

  bool active;

  cexpr * cache;
  int maxlen;
};

typedef struct _bview bview;




bview * bview_new()
{
  bview * p = ((bview *) pcalloc (1, sizeof(bview),__func__,__FILE__,__LINE__));

  p->tx = p->ty = 0;
  p->sx = p->sy = p->sz = 1.;
  p->quat[0] = p->quat[1] = p->quat[2] = 0; p->quat[3] = 1;
  p->fov = 24.;
  gl_trackball (p->quat, 0.0, 0.0, 0.0, 0.0);


  p->bg[0] = 1; p->bg[1] = 1; p->bg[2] = 1;



  p->res = 1.;
  p->lc = 0.001;

  p->samples = 4;
  p->width = 600*p->samples, p->height = 600*p->samples;


  disable_fpe (FE_DIVBYZERO|FE_INVALID);

  p->fb = framebuffer_new (p->width, p->height);

  init_gl();
  p->active = false;

  enable_fpe (FE_DIVBYZERO|FE_INVALID);

  return p;
}




void bview_destroy (bview * p)
{
  framebuffer_destroy (p->fb);
  if (p->cache)
    free_cexpr (p->cache);
  pfree (p,__func__,__FILE__,__LINE__);
}




static bview * _view = NULL;






static void destroy_view()
{
  if (!(_view)) qassert ("/home/gongminjiang/basilisk/src/view.h", 237, "_view");
  bview_destroy (_view);
}

bview * get_view() {
  if (!_view) {
    _view = bview_new();
    free_solver_func_add (destroy_view);
  }
  return _view;
}




static void redraw() {
  bview * view = get_view();


  disable_fpe (FE_DIVBYZERO|FE_INVALID);

  glMatrixMode (0x1701);
  glLoadIdentity ();

  if (view->far <= view->near) {
    double max = 2.;
    gl_perspective (view->fov, view->width/(float)view->height, 1., 1. + 2.*max);

    glMatrixMode (0x1700);
    glLoadIdentity ();
    glTranslatef (view->tx, view->ty, - (1. + max));
  }
  else {
    gl_perspective (view->fov, view->width/(float)view->height,
      view->near, view->far);

    glMatrixMode (0x1700);
    glLoadIdentity ();
    glTranslatef (view->tx, view->ty, view->tz);
  }

  GLfloat m[4][4];
  gl_build_rotmatrix (m, view->quat);
  glMultMatrixf (&m[0][0]);

  if (view->gfsview) {
    m[0][0] = 0., m[0][1] = 0., m[0][2] = -1.;
    m[1][0] = 0., m[1][1] = -1., m[1][2] = 0.;
    m[2][0] = 1., m[2][1] = 0., m[2][2] = 0.;
    glMultMatrixf (&m[0][0]);
  }

  glScalef (view->sx/L0, view->sy/L0, view->sz/L0);

  glClearColor (view->bg[0], view->bg[1], view->bg[2], 0.);
  glClear (0x00004000|0x00000100);

  gl_get_frustum (&view->frustum);

  view->active = true;
  view->ni = 0;
}




bview * draw() {
  bview * view = get_view();
  if (!view->active)
    redraw();
  else


    disable_fpe (FE_DIVBYZERO|FE_INVALID);
  glMatrixMode (0x1701);
  glTranslatef (0, 0, - 1e-4);
  return view;
}







typedef void * pointer;
#line 331 "/home/gongminjiang/basilisk/src/view.h"
typedef struct {
  GLubyte a[4];
} RGBA;

static void compose_image_op (void * pin, void * pout, int * len,
         MPI_Datatype * dptr)
{
  RGBA * rin = pin, * out = pout;
  for (int i = 0; i < *len; i++,rin++,out++)
    if (out->a[3] == 0)
      *out = *rin;
}

     
static pointer compose_image (bview * view)
{tracing("compose_image","/home/gongminjiang/basilisk/src/view.h",345);
  unsigned char * image = framebuffer_image (view->fb);
  if (!(image)) qassert ("/home/gongminjiang/basilisk/src/view.h", 348, "image");
  if (npe() > 1) {
    MPI_Op op;
    MPI_Op_create (compose_image_op, true, &op);
    MPI_Datatype rgba;
    MPI_Type_contiguous (4, MPI_BYTE, &rgba);
    MPI_Type_commit (&rgba);
    int size = view->width*view->height;
    if (pid() == 0)
      MPI_Reduce (MPI_IN_PLACE, image, size, rgba, op, 0, MPI_COMM_WORLD);
    else
      MPI_Reduce (image, image, size, rgba, op, 0, MPI_COMM_WORLD);
    MPI_Op_free (&op);
    MPI_Type_free (&rgba);
  }
  {end_tracing("compose_image","/home/gongminjiang/basilisk/src/view.h",363);return image;}
end_tracing("compose_image","/home/gongminjiang/basilisk/src/view.h",364);}
#line 423 "/home/gongminjiang/basilisk/src/view.h"
#line 1 "vertexbuffer.h"
#line 1 "/home/gongminjiang/basilisk/src/vertexbuffer.h"
#line 14 "/home/gongminjiang/basilisk/src/vertexbuffer.h"
struct {

  Array * position, * normal, * color, * index;
  float modelview[16];
  int type;
  int dim;
  int vertex, nvertex;
  bool visible;


  int line_loop, lines, line_strip ;
  int quads, polygon, fan;
  int state;
} VertexBuffer = {
  .visible = false,
  .modelview = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 }
};

static void vertex_buffer_push_index (unsigned int i)
{
  i -= VertexBuffer.vertex;
  array_append (VertexBuffer.index, &i, sizeof(unsigned int));
}

void vertex_buffer_setup()
{
  VertexBuffer.nvertex = 0;
  VertexBuffer.type = -1;
  VertexBuffer.dim = -1;
  VertexBuffer.position = array_new();
  VertexBuffer.normal = array_new();
  VertexBuffer.color = array_new();
  VertexBuffer.index = array_new();
}

void vertex_buffer_free()
{
  array_free (VertexBuffer.position);
  VertexBuffer.position = NULL;
  array_free (VertexBuffer.normal);
  VertexBuffer.normal = NULL;
  array_free (VertexBuffer.color);
  VertexBuffer.color = NULL;
  array_free (VertexBuffer.index);
  VertexBuffer.index = NULL;
}

static void vertex_buffer_glBegin (int state)
{
  if (VertexBuffer.index) {

    glGetFloatv (0x0BA6, VertexBuffer.modelview);

    bview * view = get_view();

    float q[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,
      - view->tx, - view->ty, 3, 1 };
    matrix_multiply (q, VertexBuffer.modelview);
    for (int i = 0; i < 16; i++)
      VertexBuffer.modelview[i] = q[i];

    gl_build_rotmatrix ((float (*)[4])q, view->quat);
    do { float __tmp = q[1]; q[1] = q[4]; q[4] = __tmp; } while(0);
    do { float __tmp = q[2]; q[2] = q[8]; q[8] = __tmp; } while(0);
    do { float __tmp = q[6]; q[6] = q[9]; q[9] = __tmp; } while(0);
    matrix_multiply (q, VertexBuffer.modelview);
    for (int i = 0; i < 16; i++)
      VertexBuffer.modelview[i] = q[i];

    VertexBuffer.state = state;
    switch (state) {
    case 0x0002:
      VertexBuffer.line_loop = VertexBuffer.nvertex;
      break;
    case 0x0001:
      VertexBuffer.lines = VertexBuffer.nvertex;
      break;
    case 0x0003:
      VertexBuffer.line_strip = VertexBuffer.nvertex;
      break;
    case 0x0007:
      VertexBuffer.quads = VertexBuffer.nvertex;
      break;
    case 0x0009:
      VertexBuffer.polygon = VertexBuffer.nvertex;
      break;
    case 0x0006:
      VertexBuffer.fan = VertexBuffer.nvertex;
      break;
    default:
      fprintf (ferr, "glBegin (%d) not implemented yet\n", state);
      break;
    }
  }
  else
    glBegin (state);
}

static void vertex_buffer_glEnd()
{
  if (VertexBuffer.index) {
    int type = -1;
    switch (VertexBuffer.state) {

    case 0x0002:
      for (int i = VertexBuffer.line_loop; i < VertexBuffer.nvertex - 1; i++) {
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      vertex_buffer_push_index (VertexBuffer.nvertex - 1);
      vertex_buffer_push_index (VertexBuffer.line_loop);
      type = 0;
      break;

    case 0x0001:
      for (int i = VertexBuffer.lines; i < VertexBuffer.nvertex; i += 2) {
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      type = 0;
      break;

    case 0x0003:
      for (int i = VertexBuffer.line_strip; i < VertexBuffer.nvertex - 1; i++) {
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      type = 0;
      break;

    case 0x0007:
      for (int i = VertexBuffer.quads; i < VertexBuffer.nvertex; i += 4)
 for (int j = 1; j <= 2; j++) {
   vertex_buffer_push_index (i);
   vertex_buffer_push_index (i + j);
   vertex_buffer_push_index (i + j + 1);
 }
      type = 1;
      break;

    case 0x0009:
      for (int j = 1; j <= VertexBuffer.nvertex - VertexBuffer.polygon - 2;
    j++) {
 vertex_buffer_push_index (VertexBuffer.polygon);
 vertex_buffer_push_index (VertexBuffer.polygon + j);
 vertex_buffer_push_index (VertexBuffer.polygon + j + 1);
      }
      type = 1;
      break;

    case 0x0006:
      for (int i = VertexBuffer.fan + 1; i < VertexBuffer.nvertex - 1; i++) {
 vertex_buffer_push_index (VertexBuffer.fan);
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      type = 1;
      break;

    default:
      break;
    }
    VertexBuffer.state = 0;
    if (VertexBuffer.type >= 0 && type >= 0) {

      if (!(VertexBuffer.type == type)) qassert ("/home/gongminjiang/basilisk/src/vertexbuffer.h", 179, "VertexBuffer.type == type");
    }
    else
      VertexBuffer.type = type;
  }
  else
    glEnd();
}

static void vertex_buffer_glColor3f (float r, float g, float b)
{
  if (VertexBuffer.color) {
    struct { float x, y, z; } color = {r, g, b};
    array_append (VertexBuffer.color, &color, 3*sizeof(float));
  }
  else
    glColor3f (r, g, b);
}

static void vertex_buffer_glNormal3d (double nx, double ny, double nz)
{
  if (VertexBuffer.normal) {
    struct { float x, y, z; } normal = {nx, ny, nz};
    array_append (VertexBuffer.normal, &normal, 3*sizeof(float));
  }
  else
    glNormal3d (nx, ny, nz);
}

static void vertex_buffer_glVertex3d (double x, double y, double z)
{
  if (VertexBuffer.position) {
    if (VertexBuffer.dim < 3)
      VertexBuffer.dim = 3;
    float v[4] = {x, y, z, 1.};
    vector_multiply (v, VertexBuffer.modelview);
    array_append (VertexBuffer.position, v, 3*sizeof(float));
    VertexBuffer.nvertex++;
  }
  else
    glVertex3d (x, y, z);
}

static void vertex_buffer_glVertex2d (double x, double y)
{
  if (VertexBuffer.position) {
    if (VertexBuffer.dim < 2)
      VertexBuffer.dim = 2;
    float v[4] = {x, y, 0, 1.};
    vector_multiply (v, VertexBuffer.modelview);
    array_append (VertexBuffer.position, v, 3*sizeof(float));
    VertexBuffer.nvertex++;
  }
  else
    glVertex3d (x, y, 0.);
}
#line 424 "/home/gongminjiang/basilisk/src/view.h"






#line 1 "draw.h"
#line 1 "/home/gongminjiang/basilisk/src/draw.h"





#line 1 "gl/font.h"
#line 1 "/home/gongminjiang/basilisk/src/gl/font.h"
#line 27 "/home/gongminjiang/basilisk/src/gl/font.h"
#line 1 "/home/gongminjiang/basilisk/src/ast/std/stdio.h"
#include <stdio.h>
#line 28 "/home/gongminjiang/basilisk/src/gl/font.h"
#line 1 "gl/og_font.h"
#line 1 "/home/gongminjiang/basilisk/src/gl/og_font.h"




#line 1 "gl/gl.h"
#line 6 "/home/gongminjiang/basilisk/src/gl/og_font.h"

typedef struct tagSOG_StrokeVertex SOG_StrokeVertex;
struct tagSOG_StrokeVertex
{
    GLfloat X, Y;
};

typedef struct tagSOG_StrokeStrip SOG_StrokeStrip;
struct tagSOG_StrokeStrip
{
    int Number;
    const SOG_StrokeVertex *Vertices;
};

typedef struct tagSOG_StrokeChar SOG_StrokeChar;
struct tagSOG_StrokeChar
{
    GLfloat Right;
    int Number;
    const SOG_StrokeStrip* Strips;
};

typedef struct tagSOG_StrokeFont SOG_StrokeFont;
struct tagSOG_StrokeFont
{
    char *Name;
    int Quantity;
    GLfloat Height;
    const SOG_StrokeChar **Characters;
};
#line 29 "/home/gongminjiang/basilisk/src/gl/font.h"
#line 39 "/home/gongminjiang/basilisk/src/gl/font.h"
extern SOG_StrokeFont ogStrokeMonoRoman;
#line 48 "/home/gongminjiang/basilisk/src/gl/font.h"
static SOG_StrokeFont *oghStrokeByID( void *font )
{


    if( font == ((void *)0x0001) )
        return &ogStrokeMonoRoman;

    fprintf (ferr, "stroke font %p not found", font );
    return 0;
}
#line 83 "/home/gongminjiang/basilisk/src/gl/font.h"
void gl_StrokeCharacter( int character )
{
    void *fontID = ((void *)0x0001);
    const SOG_StrokeChar *schar;
    const SOG_StrokeStrip *strip;
    int i, j;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( !font ||
        ( 1 > character ) ||
        ( font->Quantity < character ) )
        return;

    schar = font->Characters[ character ];
    if( schar )
    {
        strip = schar->Strips;

        for( i = 0; i < schar->Number; i++, strip++ )
        {
            vertex_buffer_glBegin( 0x0003 );
            for( j = 0; j < strip->Number; j++ )
                vertex_buffer_glVertex2d( strip->Vertices[ j ].X, strip->Vertices[ j ].Y );
            vertex_buffer_glEnd( );
        }
        glTranslatef( schar->Right, 0.0, 0.0 );
    }
}
#line 147 "/home/gongminjiang/basilisk/src/gl/font.h"
void gl_StrokeString( const char *string )
{
    void *fontID = ((void *)0x0001);
    int i, j;
    float length = 0.0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );
    unsigned char c;

    if( font && string )





        while(( c = *string++ ))
       if( c < font->Quantity ) {
                if( c == '\n' )
                {
                    glTranslatef ( -length, -( float )( font->Height ), 0.0 );
                    length = 0.0;
                }
                else
                {
                    const SOG_StrokeChar *schar =
                        font->Characters[ c ];
                    if( schar )
                    {
                        const SOG_StrokeStrip *strip = schar->Strips;

                        for( i = 0; i < schar->Number; i++, strip++ )
                        {
                            vertex_buffer_glBegin( 0x0003 );

                            for( j = 0; j < strip->Number; j++ )
                                vertex_buffer_glVertex2d( strip->Vertices[ j ].X,
                                            strip->Vertices[ j ].Y);

                            vertex_buffer_glEnd( );
                        }

                        length += schar->Right;
                        glTranslatef( schar->Right, 0.0, 0.0 );
                    }
                }
     }
}
#line 226 "/home/gongminjiang/basilisk/src/gl/font.h"
float gl_StrokeWidth( int character )
{
    void *fontID = ((void *)0x0001);
    float ret = 0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( font &&
        ( 0 < character ) &&
        ( font->Quantity > character ) )
    {
        const SOG_StrokeChar *schar = font->Characters[ character ];
        if( schar )
            ret = schar->Right;
    }

    return ret;
}
#line 269 "/home/gongminjiang/basilisk/src/gl/font.h"
float gl_StrokeLength( const char *string )
{
    void *fontID = ((void *)0x0001);
    unsigned char c;
    float length = 0.0;
    float this_line_length = 0.0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( font && string )
        while(( c = *string++ ))
            if( c < font->Quantity )
            {
                if( c == '\n' )
                {
                    if( length < this_line_length )
                        length = this_line_length;
                    this_line_length = 0.0;
                }
                else
                {
                    const SOG_StrokeChar *schar =
                        font->Characters[ c ];
                    if( schar )
                        this_line_length += schar->Right;
                }
            }

    if( length < this_line_length )
        length = this_line_length;
    return length;
}
#line 321 "/home/gongminjiang/basilisk/src/gl/font.h"
GLfloat gl_StrokeHeight()
{
    void *fontID = ((void *)0x0001);
    GLfloat ret = 0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( font )
        ret = font->Height;

    return ret;
}
#line 7 "/home/gongminjiang/basilisk/src/draw.h"




void clear()
{
  bview * view = get_view();
  if (view->active)
    view->active = false;
  draw();
}
#line 49 "/home/gongminjiang/basilisk/src/draw.h"
struct _view_set {
  float tx, ty;
  float fov;
  float quat[4];
  float sx, sy, sz;
  unsigned width, height, samples;
  float bg[3];
  float theta, phi, psi;
  bool relative;
  float tz, near, far;
  float res;
  char * camera;
  void (* map) (coord *);
  int cache;
  float p1x, p1y, p2x, p2y;
  bview * view;
};

void view (struct _view_set p)
{
  bview * v = p.view ? p.view : get_view();
  if (p.fov) {
    if (p.relative)
      v->fov += (0.1 + 3.*v->fov)*p.fov;
    else
      v->fov = p.fov;
    v->fov = clamp(v->fov,0.01,100.);
  }
  for (int i = 0; i < 4; i++)
    if (p.quat[i]) {
      for (int j = 0; j < 4; j++)
 v->quat[j] = p.quat[j];
      break;
    }
  if (p.tx) v->tx = p.relative ? v->tx + p.tx*0.02*(0.01 + 3.*v->fov) : p.tx;
  if (p.ty) v->ty = p.relative ? v->ty + p.ty*0.02*(0.01 + 3.*v->fov) : p.ty;
  if (p.sx) v->sx = p.sx;
  if (p.sy) v->sy = p.sy;
  if (p.sz) v->sz = p.sz;
  if (p.bg[0] || p.bg[1] || p.bg[2])
    for (int i = 0; i < 3; i++)
      v->bg[i] = p.bg[i];

  if (p.camera) {
    v->gfsview = false;
    if (strlen(p.camera) >= 4 &&
 !strcmp (&p.camera[strlen(p.camera) - 4], ".gfv")) {
      FILE * fp = fopen (p.camera, "r");
      if (!fp) {
 perror (p.camera);
 exit (1);
      }
      char s[81];
      float q[4], fov;
      int nq = 0, nf = 0;
      while (fgets (s, 81, fp) && (!nq || !nf)) {
 if (!nq)
   nq = sscanf (s, "  q0 = %f q1 = %f q2 = %f q3 = %f",
         &q[0], &q[1], &q[2], &q[3]);
 if (!nf)
   nf = sscanf (s, "  fov = %f", &fov);
      }
      if (nq != 4 || nf != 1) {
 fprintf (ferr, "%s: not a valid gfv file\n", p.camera);
 exit (1);
      }
      for (int j = 0; j < 4; j++)
 v->quat[j] = q[j];
      v->fov = fov;
      v->gfsview = true;
    }
    else if (!strcmp (p.camera, "left"))
      gl_axis_to_quat ((float[]){0,1,0}, - 3.14159265358979/2., v->quat);
    else if (!strcmp (p.camera, "right"))
      gl_axis_to_quat ((float[]){0,1,0}, 3.14159265358979/2., v->quat);
    else if (!strcmp (p.camera, "top"))
      gl_axis_to_quat ((float[]){1,0,0}, - 3.14159265358979/2., v->quat);
    else if (!strcmp (p.camera, "bottom"))
      gl_axis_to_quat ((float[]){1,0,0}, 3.14159265358979/2., v->quat);
    else if (!strcmp (p.camera, "front"))
      gl_axis_to_quat ((float[]){0,0,1}, 0., v->quat);
    else if (!strcmp (p.camera, "back"))
      gl_axis_to_quat ((float[]){0,1,0}, 3.14159265358979, v->quat);
    else if (!strcmp (p.camera, "iso")) {
      gl_axis_to_quat ((float[]){0,1,0}, 3.14159265358979/4., v->quat);
      float q[4];
      gl_axis_to_quat ((float[]){1,0,0}, - 3.14159265358979/4., q);
      gl_add_quats(q, v->quat, v->quat);
    }
    else {
      fprintf (ferr, "view(): unknown camera '%s'\n", p.camera);
      exit (1);
    }
  }
  else if (p.theta || p.phi || p.psi) {
    v->gfsview = false;
    float q[4];
    gl_axis_to_quat ((float[]){1,0,0}, - p.phi, q);
    if (p.relative) {
      float q1[4];
      gl_axis_to_quat ((float[]){0,1,0}, p.theta, q1);
      gl_add_quats(q, q1, q1);
      float q2[4];
      gl_axis_to_quat ((float[]){0,0,1}, p.psi, q2);
      gl_add_quats(q1, q2, q2);
      gl_add_quats(q2, v->quat, v->quat);
    }
    else {
      gl_axis_to_quat ((float[]){0,1,0}, p.theta, v->quat);
      gl_add_quats(q, v->quat, v->quat);
      gl_axis_to_quat ((float[]){0,0,1}, p.psi, q);
      gl_add_quats(q, v->quat, v->quat);
    }
  }

  if (p.map)
    v->map = p.map;

  if (p.p1x || p.p1y || p.p2x || p.p2y) {
    float q[4];
    gl_trackball(q, p.p1x, p.p1y, p.p2x, p.p2y);
    gl_add_quats (q, v->quat, v->quat);
  }

  if (p.far > p.near) {
    v->tz = p.tz;
    v->far = p.far;
    v->near = p.near;
  }

  if (p.res)
    v->res = p.res;

  if ((p.width && p.width != v->width) ||
      (p.height && p.height != v->height) ||
      (p.samples && p.samples != v->samples)) {
    v->width = v->width/v->samples;
    v->height = v->height/v->samples;
    if (p.width) v->width = p.width;
    if (p.height) v->height = p.height;
    if (p.samples) v->samples = p.samples;
    v->width *= v->samples;
    v->height *= v->samples;
    framebuffer_destroy (v->fb);


    disable_fpe (FE_DIVBYZERO|FE_INVALID);

    v->fb = framebuffer_new (v->width, v->height);
    init_gl();

    enable_fpe (FE_DIVBYZERO|FE_INVALID);
  }

  if (p.cache > 0) {
    v->cache = pcalloc (1, sizeof (cexpr),__func__,__FILE__,__LINE__);
    v->maxlen = p.cache;
  }

  clear();
}







struct _translate {
  float x, y, z;
};

void begin_translate (struct _translate p)
{
  bview * view = draw();
  glMatrixMode (0x1700);
  glPushMatrix();
  glTranslatef (p.x, p.y, p.z);
  gl_get_frustum (&view->frustum);
}

void end_translate()
{
  bview * view = draw();
  glMatrixMode (0x1700);
  glPopMatrix();
  gl_get_frustum (&view->frustum);
}
#line 246 "/home/gongminjiang/basilisk/src/draw.h"
struct _mirror {
  coord n;
  double alpha;
};

void begin_mirror (struct _mirror p)
{
  bview * view = draw();
  glMatrixMode (0x1700);
  glPushMatrix();
  normalize (&p.n);
  GLfloat s[16], t[16];
  s[0] = 1. - 2.*p.n.x*p.n.x;
  s[1] = - 2.*p.n.x*p.n.y; s[2] = - 2.*p.n.x*p.n.z;
  s[3] = 0.;
  s[4] = s[1];
  s[5] = 1. - 2.*p.n.y*p.n.y; s[6] = - 2.*p.n.y*p.n.z;
  s[7] = 0.;
  s[8] = s[2]; s[9] = s[6]; s[10] = 1. - 2.*p.n.z*p.n.z;
  s[11] = 0.;
  s[12] = 0.; s[13] = 0.; s[14] = 0.;
  s[15] = 1.;

  t[0] = 1.; t[1] = 0.; t[2] = 0.; t[3] = 0.;
  t[4] = 0.; t[5] = 1.; t[6] = 0.; t[7] = 0.;
  t[8] = 0.; t[9] = 0.; t[10] = 1.; t[11] = 0.;
  t[12] = - 2.*p.n.x*p.alpha;
  t[13] = - 2.*p.n.y*p.alpha;
  t[14] = - 2.*p.n.z*p.alpha;
  t[15] = 1.;
  matrix_multiply (s, t);
  glMultMatrixf (s);
  gl_get_frustum (&view->frustum);
  view->reversed = !view->reversed;
}

void end_mirror() {
  end_translate();
  bview * view = draw();
  view->reversed = !view->reversed;
}







static void mapped_position (bview * view, coord * p, double * r)
{
  double x = p->x, y = p->y, z = p->z, rm = 0.;
  view->map (p);
  for (int i = -1; i <= 1; i += 2)
    for (int j = -1; j <= 1; j += 2)
      for (int k = -1; k <= 1; k += 2) {
 coord q = {x + i**r, y + j**r, z + k**r};
 view->map (&q);
 double pq = sq(p->x - q.x) + sq(p->y - q.y) + sq(p->z - q.z);
 if (pq > rm)
   rm = pq;
      }
  *r = sqrt (rm);
}

#define foreach_visible(view)\
foreach_cell() {\
\
  double _r = Delta*0.71;\
\
\
\
  coord _p = {x, y, z};\
  if ((view)->map)\
    mapped_position (view, &_p, &_r);\
  if (VertexBuffer.visible &&\
      !sphere_in_frustum (_p.x, _p.y, _p.z, _r, &(view)->frustum))\
    continue;\
  if (is_leaf(cell) ||\
      (VertexBuffer.visible &&\
       sphere_diameter (_p.x, _p.y, _p.z, _r/L0, &(view)->frustum)\
       < (view)->res)) {\
    if (is_active(cell) && is_local(cell)) {\

#line 328

#define end_foreach_visible()\
    }\
    continue;\
  }\
}\
end_foreach_cell();\

#line 335

#line 386 "/home/gongminjiang/basilisk/src/draw.h"
static bool _reversed = false;

static void begin_draw_lines (bview * view, float color[3], float lw)
{
  glMatrixMode (0x1701);
  glPushMatrix();
  glTranslatef (0., 0., view->lc*view->fov/24.);
  vertex_buffer_glColor3f (color[0], color[1], color[2]);
  glLineWidth (view->samples*(lw > 0. ? lw : 1.));
  _reversed = view->reversed;
  view->reversed = false;
}

static void end_draw_lines()
{
  glMatrixMode (0x1701);
  glPopMatrix();
  bview * view = draw();
  view->reversed = _reversed;
}

static inline double interp (Point point, coord p, scalar col) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  struct _interpolate _r = { col, x + p.x*Delta, y + p.y*Delta, z + p.z*Delta };
  return interpolate_linear (point, _r);
}

static double evaluate_expression (Point point, Node * n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  if (!(n)) qassert ("/home/gongminjiang/basilisk/src/draw.h", 414, "n");
  switch (n->type) {
  case '1': return n->d.value;
  case '+': return (evaluate_expression (point, n->e[0]) +
      evaluate_expression(point, n->e[1]));
  case '-': return (evaluate_expression (point, n->e[0]) -
      evaluate_expression(point, n->e[1]));
  case '*': return (evaluate_expression (point, n->e[0]) *
      evaluate_expression(point, n->e[1]));
  case '/': return (evaluate_expression (point, n->e[0]) /
      evaluate_expression(point, n->e[1]));
  case '^': return pow (evaluate_expression (point, n->e[0]),
   evaluate_expression(point, n->e[1]));
  case '>': return (evaluate_expression (point, n->e[0]) >
      evaluate_expression(point, n->e[1]));
  case '<': return (evaluate_expression (point, n->e[0]) <
      evaluate_expression(point, n->e[1]));
  case 'L': return (evaluate_expression (point, n->e[0]) <=
      evaluate_expression(point, n->e[1]));
  case 'G': return (evaluate_expression (point, n->e[0]) >=
      evaluate_expression(point, n->e[1]));
  case '=': return (evaluate_expression (point, n->e[0]) ==
      evaluate_expression(point, n->e[1]));
  case 'i': return (evaluate_expression (point, n->e[0]) !=
      evaluate_expression(point, n->e[1]));
  case 'O': return (evaluate_expression (point, n->e[0]) ||
      evaluate_expression(point, n->e[1]));
  case 'A': return (evaluate_expression (point, n->e[0]) &&
      evaluate_expression(point, n->e[1]));
  case '?': return (evaluate_expression (point, n->e[0]) ?
      evaluate_expression(point, n->e[1]) :
      evaluate_expression(point, n->e[2]));
  case 'm': return - evaluate_expression (point, n->e[0]);
  case 'f': return n->d.func (evaluate_expression (point, n->e[0]));
  case 'v': {
    scalar s = {n->s};
    int k[3] = {0,0,0};
    for (int i = 0; i < 3; i++)
      if (n->e[i])
 k[i] = evaluate_expression (point, n->e[i]);
    return val(s,k[0],k[1],k[2]);
  }
  case 'D': return Delta;
  case 'x': return x;
  case 'y': return y;
  case 'z': return z;
  default:
    fprintf (ferr, "unknown operation type '%c'\n", n->type);
    if (!(false)) qassert ("/home/gongminjiang/basilisk/src/draw.h", 462, "false");
  }
  return undefined;
}


#line 412
static void _stencil_evaluate_expression (Point point, Node * n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES; 
      
  switch (n->type) {
  case '1': return ;
  case '+': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 417
return  
;}
  case '-': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 419
return  
;}
  case '*': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 421
return  
;}
  case '/': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 423
return  
;}
  case '^': {_stencil_evaluate_expression (point, n->e[0]);
   _stencil_evaluate_expression(point, n->e[1]);
#line 425
return  
;}
  case '>': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 427
return  
;}
  case '<': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 429
return  
;}
  case 'L': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 431
return  
;}
  case 'G': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 433
return  
;}
  case '=': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 435
return  
;}
  case 'i': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 437
return  
;}
  case 'O': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 439
return  
;}
  case 'A': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 441
return  
;}
  case '?': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
      _stencil_evaluate_expression(point, n->e[2]);
#line 443
return   

;}
  case 'm': { _stencil_evaluate_expression (point, n->e[0]);return ;}
  case 'f': {_stencil_evaluate_expression (point, n->e[0]);return  ;}
  case 'v': {
    scalar s = {n->s};   
    
    for (int i = 0; i < 3; i++)
      if (n->e[i])
 { _stencil_evaluate_expression (point, n->e[i]); } 
_stencil_val(s,o_stencil,o_stencil,o_stencil);
    
#line 454
return;
  }
  case 'D': return ;
  case 'x': return ;
  case 'y': return ;
  case 'z': return ; 
     
    
        
  }
  return ;
}

static bool assemble_node (Node * n)
{
  if (n->type == 'v') {
    char * id = n->d.id;
    scalar s = lookup_field (id);
    if (s.i >= 0)
      n->s = s.i;
    else {
      n->s = -1;
      if (!strcmp (id, "Delta"))
 reset_node_type (n, 'D');
      else if (!strcmp (id, "x"))
 reset_node_type (n, 'x');
      else if (!strcmp (id, "y"))
 reset_node_type (n, 'y');
      else if (!strcmp (id, "z"))
 reset_node_type (n, 'z');
      else {
 typedef struct { char * name; double val; } Constant;
 static Constant constants[] = {
   {"pi", 3.14159265358979 },
   {"nodata", HUGE },
   {"HUGE", HUGE },
   { NULL },
 };
 Constant * p = constants;
 while (p->name) {
   if (!strcmp (p->name, id)) {
     reset_node_type (n, '1');
     n->d.value = p->val;
     break;
   }
   p++;
 }
 if (n->type == 'v') {
   fprintf (ferr, "unknown identifier '%s'\n", id);
   return false;
 }
      }
    }
  }
  for (int i = 0; i < 3; i++)
    if (n->e[i] && !assemble_node (n->e[i]))
      return false;
  return true;
}

static scalar compile_expression (char * expr, bool * isexpr)
{
  *isexpr = false;
  if (!expr)
    return (scalar){-1};

  bview * view = get_view();
  scalar s;
  if (view->cache && (s = get_cexpr (view->cache, expr)).i >= 0)
    return s;

  Node * node = parse_node (expr);
  if (node == NULL) {
    fprintf (ferr, "'%s': syntax error\n", expr);
    return (scalar){-1};
  }
  if (!assemble_node (node)) {
    free_node (node);
    return (scalar){-1};
  }
  if (node->type == 'v' && node->e[0] == NULL) {
    scalar s = {node->s};
    free_node (node);
    return s;
  }
  s = new_scalar("s");
  pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
  _attribute[s.i].name = pstrdup (expr,__func__,__FILE__,__LINE__);
  foreach_stencil()
    {_stencil_val_a(s,0,0,0); _stencil_evaluate_expression (point, node); }end_foreach_stencil();
  {
#line 542
foreach()
    val(s,0,0,0) = evaluate_expression (point, node);end_foreach();}
  restriction (((scalar[]){s,{-1}}));
  free_node (node);

  if (view->cache)
    view->cache = add_cexpr (view->cache, view->maxlen, expr, s);
  else
    *isexpr = true;
  return s;
}
#line 612 "/home/gongminjiang/basilisk/src/draw.h"
static void begin_colorized (float fc[3], bool constant_color,
        double cmap[127][3], bool use_texture)
{

  if (use_texture) {
    GLfloat texture[3*256];
    for (int i = 0; i < 256; i++) {
      color j = colormap_color (cmap, i/255., 0, 1);
      texture[3*i] = j.r/255.;
      texture[3*i + 1] = j.g/255.;
      texture[3*i + 2] = j.b/255.;
    }
    glTexImage1D (0x0DE0, 0, 0x1907, 256,0, 0x1907, 0x1406, texture);
    glTexParameteri (0x0DE0, 0x2801, 0x2601);
    glTexParameteri (0x0DE0, 0x2800, 0x2601);
    glTexParameteri (0x0DE0, 0x2802, 0x812F);
    glTexParameteri (0x0DE0, 0x2803, 0x812F);
    glEnable (0x0DE0);
  }
  if (constant_color)
    vertex_buffer_glColor3f (fc[0], fc[1], fc[2]);
}

static void end_colorized() {
  glDisable (0x0DE0);
}
#line 669 "/home/gongminjiang/basilisk/src/draw.h"
struct _draw_vof {
  char * c;
  char * s;
  bool edges;
  double larger;
  int filled;

  char * color;
  double min, max, spread;
  bool linear;
  colormap map;
  float fc[3], lc[3], lw;
  bool expr;
};







static bool cfilter (Point point, scalar c, double cmin)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  double cmin1 = 4.*cmin;
  if (val(c,0,0,0) <= cmin) {
    
      if (val(c,1,0,0) >= 1. - cmin1 || val(c,-1,0,0) >= 1. - cmin1)
 return true;
      
#line 695
if (val(c,0,1,0) >= 1. - cmin1 || val(c,0,-1,0) >= 1. - cmin1)
 return true;
    return false;
  }
  if (val(c,0,0,0) >= 1. - cmin) {
    
      if (val(c,1,0,0) <= cmin1 || val(c,-1,0,0) <= cmin1)
 return true;
      
#line 701
if (val(c,0,1,0) <= cmin1 || val(c,0,-1,0) <= cmin1)
 return true;
    return false;
  }
  int n = 0;
  double min = HUGE, max = - HUGE;
  {foreach_neighbor(1) {
    if (val(c,0,0,0) > cmin && val(c,0,0,0) < 1. - cmin && ++n >= (1 << 2))
      return true;
    if (val(c,0,0,0) > max) max = val(c,0,0,0);
    if (val(c,0,0,0) < min) min = val(c,0,0,0);
  }end_foreach_neighbor()}
  return max - min > 0.5;
}








#line 690
static void _stencil_cfilter (Point point, scalar c, double cmin)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;   
  
_stencil_val(c,0,0,0); {
    
      {_stencil_val(c,1,0,0); _stencil_val(c,-1,0,0); 
           }
      
#line 695
{_stencil_val(c,0,1,0); _stencil_val(c,0,-1,0); 
           } 
    
  }
_stencil_val(c,0,0,0); {
    
      {_stencil_val(c,1,0,0); _stencil_val(c,-1,0,0); 
       }
      
#line 701
{_stencil_val(c,0,1,0); _stencil_val(c,0,-1,0); 
       } 
    
  }          
     
       
  
  
  {
#line 707
foreach_neighbor(1) {
_stencil_val(c,0,0,0); _stencil_val(c,0,0,0);

_stencil_val(c,0,0,0); { _stencil_val(c,0,0,0); }
_stencil_val(c,0,0,0); { _stencil_val(c,0,0,0); } 
      
                  
       
       
  
#line 712
}end_foreach_neighbor()}
  return     ;
}

static void glvertex3d (bview * view, double x, double y, double z) {
  if (view->map) {
    coord p = {x, y, z};
    view->map (&p);
    vertex_buffer_glVertex3d (p.x, p.y, p.z);
  }
  else
    vertex_buffer_glVertex3d (x, y, z);
}


static void glvertex2d (bview * view, double x, double y) {
  if (view->map) {
    coord p = {x, y, 0.};
    view->map (&p);
    vertex_buffer_glVertex2d (p.x, p.y);
  }
  else
    vertex_buffer_glVertex2d (x, y);
}

static void glvertex_normal3d (bview * view, Point point, vector n,
          double xp, double yp, double zp)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;
  coord v = {(xp - x)/Delta, (yp - y)/Delta}, np;
  
    np.x = - interp (point, v, n.x);
    
#line 742
np.y = - interp (point, v, n.y);
  vertex_buffer_glNormal3d (np.x, np.y, 1.);
  glvertex3d (view, xp, yp, zp);
}


     
bool draw_vof (struct _draw_vof p)
{tracing("draw_vof","/home/gongminjiang/basilisk/src/draw.h",749);
  scalar c = lookup_field (p.c);
  if (c.i < 0) {
    fprintf (ferr, "draw_vof(): no field named '%s'\n", p.c);
    {end_tracing("draw_vof","/home/gongminjiang/basilisk/src/draw.h",754);return false;}
  }
  vector s = lookup_vector (p.s);

  scalar col = {-1}; if (p.color && strcmp (p.color, "level")) { col = compile_expression (p.color, &p.expr); if (col.i < 0) {end_tracing("draw_vof","/home/gongminjiang/basilisk/src/draw.h",758);return false;} boundary_internal ((scalar *)((scalar[]){col,{-1}}), "/home/gongminjiang/basilisk/src/draw.h", 758); } double cmap[127][3]; if (p.color) { if (p.min == 0 && p.max == 0) { if (col.i < 0) p.min = 0, p.max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (p.spread < 0.) p.min = s.min, p.max = s.max; else { if (!p.spread) p.spread = 5.; double spread = p.spread*s.stddev; p.min = avg - spread; p.max = avg + spread; } } } if (!p.map) p.map = jet; p.map (cmap); } if ((2 > 2 || p.linear) && !p.fc[0] && !p.fc[1] && !p.fc[2]) p.fc[0] = p.fc[1] = p.fc[2] = 1.;;

  double cmin = 1e-3;



  void (* prolongation) (Point, scalar) = _attribute[c.i].prolongation;
  if (prolongation != fraction_refine) {
    _attribute[c.i].prolongation = fraction_refine;
    _attribute[c.i].dirty = true;
  }


  bview * view = draw();

  if (p.filled) {
    vertex_buffer_glColor3f (p.fc[0], p.fc[1], p.fc[2]);
    vertex_buffer_glNormal3d (0, 0, view->reversed ? -1 : 1);
    foreach_visible_stencil (view) { 
_stencil_val(c,0,0,0); _stencil_val(c,0,0,0);{                             
 
 
 
 
 
 
 
        
{_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {     
  _stencil_facet_normal (point, c, s);              
 
   
        
  _stencil_val(c,0,0,0); _stencil_val(c,0,0,0);    
     
            
      
 
   
         
       
       
              
            
         
    
       
  
        
              
            
         
    
             
         
        
         
       
       
             
       
         
     
     
 
         
       
         
     
 
 
      }      }}
                   
      
    
#line 830
}end_foreach_visible_stencil();
    {
#line 776
foreach_visible (view) {
      if ((p.filled > 0 && val(c,0,0,0) >= 1.) || (p.filled < 0 && val(c,0,0,0) <= 0.)) {
 vertex_buffer_glBegin (0x0007);
 glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
 glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
 glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
 glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
 vertex_buffer_glEnd();
 view->ni++;
      }
      else if (val(c,0,0,0) > 0. && val(c,0,0,0) < 1.) {
 coord n = facet_normal (point, c, s), r = {1.,1.};
 if (p.filled < 0)
   {
     n.x = - n.x;
     
#line 790
n.y = - n.y;}
 double alpha = line_alpha (p.filled < 0. ? 1. - val(c,0,0,0) : val(c,0,0,0), n);
 alpha += (n.x + n.y)/2.;
 
   if (n.x < 0.) alpha -= n.x, n.x = - n.x, r.x = - 1.;
   
#line 794
if (n.y < 0.) alpha -= n.y, n.y = - n.y, r.y = - 1.;
 coord v[5];
 int nv = 0;
 if (alpha >= 0. && alpha <= n.x) {
   v[nv].x = alpha/n.x, v[nv++].y = 0.;
   if (alpha <= n.y)
     v[nv].x = 0., v[nv++].y = alpha/n.y;
   else if (alpha >= n.y && alpha - n.y <= n.x) {
     v[nv].x = (alpha - n.y)/n.x, v[nv++].y = 1.;
     v[nv].x = 0., v[nv++].y = 1.;
   }
   v[nv].x = 0., v[nv++].y = 0.;
 }
 else if (alpha >= n.x && alpha - n.x <= n.y) {
   v[nv].x = 1., v[nv++].y = (alpha - n.x)/n.y;
   if (alpha >= n.y && alpha - n.y <= n.x) {
     v[nv].x = (alpha - n.y)/n.x, v[nv++].y = 1.;
     v[nv].x = 0., v[nv++].y = 1.;
   }
   else if (alpha <= n.y)
     v[nv].x = 0., v[nv++].y = alpha/n.y;
   v[nv].x = 0., v[nv++].y = 0.;
   v[nv].x = 1., v[nv++].y = 0.;
 }
 vertex_buffer_glBegin (0x0009);
 if (r.x*r.y < 0.)
   for (int i = nv - 1; i >= 0; i--)
     glvertex2d (view, x + r.x*(v[i].x - 0.5)*Delta,
   y + r.y*(v[i].y - 0.5)*Delta);
 else
   for (int i = 0; i < nv; i++)
     glvertex2d (view, x + r.x*(v[i].x - 0.5)*Delta,
   y + r.y*(v[i].y - 0.5)*Delta);
 vertex_buffer_glEnd ();
 view->ni++;
      }
    }end_foreach_visible();}
  }
  else
    {begin_draw_lines (view, p.lc, p.lw); {
      vertex_buffer_glBegin (0x0001);
      foreach_visible_stencil (view)
 {_stencil_cfilter (point, c, cmin); {  
    _stencil_facet_normal (point, c, s);     
   _stencil_val(c,0,0,0); 
                 
     
     
     
    
         
 } }end_foreach_visible_stencil();
      {
#line 835
foreach_visible (view)
 if (cfilter (point, c, cmin)) {
   coord n = facet_normal (point, c, s);
   double alpha = line_alpha (val(c,0,0,0), n);
   coord segment[2];
   if (facets (n, alpha, segment) == 2) {
     glvertex2d (view, x + segment[0].x*Delta, y + segment[0].y*Delta);
     glvertex2d (view, x + segment[1].x*Delta, y + segment[1].y*Delta);
     view->ni++;
   }
 }end_foreach_visible();}
      vertex_buffer_glEnd ();
    }end_draw_lines();}
#line 899 "/home/gongminjiang/basilisk/src/draw.h"
  if (prolongation != fraction_refine) {
    _attribute[c.i].prolongation = prolongation;
    _attribute[c.i].dirty = true;
  }


  if (p.expr) delete(((scalar[]){col,{-1}}));
  {end_tracing("draw_vof","/home/gongminjiang/basilisk/src/draw.h",906);return true;}
end_tracing("draw_vof","/home/gongminjiang/basilisk/src/draw.h",907);}
#line 918 "/home/gongminjiang/basilisk/src/draw.h"
struct _isoline {
  char * phi;
  double val;
  int n;


  char * c;
  char * s;
  bool edges;
  double larger;
  int filled;

  char * color;
  double min, max, spread;
  bool linear;
  colormap map;
  float fc[3], lc[3], lw;
  bool expr;
};

     
bool isoline (struct _isoline p)
{tracing("isoline","/home/gongminjiang/basilisk/src/draw.h",939);

  if (!p.color) p.color = p.phi;
  scalar col = {-1}; if (p.color && strcmp (p.color, "level")) { col = compile_expression (p.color, &p.expr); if (col.i < 0) {end_tracing("isoline","/home/gongminjiang/basilisk/src/draw.h",943);return false;} boundary_internal ((scalar *)((scalar[]){col,{-1}}), "/home/gongminjiang/basilisk/src/draw.h", 943); } double cmap[127][3]; if (p.color) { if (p.min == 0 && p.max == 0) { if (col.i < 0) p.min = 0, p.max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (p.spread < 0.) p.min = s.min, p.max = s.max; else { if (!p.spread) p.spread = 5.; double spread = p.spread*s.stddev; p.min = avg - spread; p.max = avg + spread; } } } if (!p.map) p.map = jet; p.map (cmap); } if ((2 > 2 || p.linear) && !p.fc[0] && !p.fc[1] && !p.fc[2]) p.fc[0] = p.fc[1] = p.fc[2] = 1.;;
  scalar phi = col,  fiso=new_scalar("fiso");
  vector  siso=new_face_vector("siso");
  p.c = "fiso", p.s = "siso";
  struct _draw_vof a = *((struct _draw_vof *)&p.c);
  if (p.n < 2) {
    fractions ((struct Fractions){phi, fiso, siso, p.val});
    draw_vof (a);
  }
  else if (p.max > p.min) {
    double dv = (p.max - p.min)/(p.n - 1);
    for (p.val = p.min; p.val <= p.max; p.val += dv) {
      fractions ((struct Fractions){phi, fiso, siso, p.val});
      draw_vof (a);
    }
  }
  if (p.expr) delete(((scalar[]){col,{-1}}));



  {delete((scalar*)((scalar[]){siso.x,siso.y,fiso,{-1}}));{end_tracing("isoline","/home/gongminjiang/basilisk/src/draw.h",963);return true;}}delete((scalar*)((scalar[]){siso.x,siso.y,fiso,{-1}}));
end_tracing("isoline","/home/gongminjiang/basilisk/src/draw.h",964);}
#line 977 "/home/gongminjiang/basilisk/src/draw.h"
struct _cells {
  coord n;
  double alpha;
  float lc[3], lw;
};

     
bool cells (struct _cells p)
{tracing("cells","/home/gongminjiang/basilisk/src/draw.h",984);
  bview * view = draw();
  {begin_draw_lines (view, p.lc, p.lw); {

    {foreach_visible (view) {
      vertex_buffer_glBegin (0x0002);
      glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
      glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
      glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
      glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
      vertex_buffer_glEnd();
      view->ni++;
    }end_foreach_visible();}
#line 1011 "/home/gongminjiang/basilisk/src/draw.h"
  }end_draw_lines();}
  {end_tracing("cells","/home/gongminjiang/basilisk/src/draw.h",1012);return true;}
end_tracing("cells","/home/gongminjiang/basilisk/src/draw.h",1013);}






struct _vectors {
  char * u;
  double scale;
  float lc[3], lw;
};

     
bool vectors (struct _vectors p)
{tracing("vectors","/home/gongminjiang/basilisk/src/draw.h",1027);

  vector u;
  struct { char x, y, z; } index = {'x', 'y', 'z'};
   {
    char name[80];
    sprintf (name, "%s.%c", p.u, index.x);
    u.x = lookup_field (name);
  } 
#line 1032
{
    char name[80];
    sprintf (name, "%s.%c", p.u, index.y);
    u.y = lookup_field (name);
  }
  bview * view = draw();
  float res = view->res;
  if (view->res < 15*view->samples)
    view->res = 15*view->samples;
  {begin_draw_lines (view, p.lc, p.lw); {
    double scale = (p.scale ? p.scale : 1.)*view->res/view->samples;
    vertex_buffer_glBegin (0x0001);
    foreach_visible_stencil (view)
      {_stencil_val(u.x,0,0,0); {      
 _stencil_val(u.y,0,0,0);_stencil_val(u.x,0,0,0);                
                                  
              
 
 
 
 
 
 
 
      }   }end_foreach_visible_stencil();
    {
#line 1044
foreach_visible (view)
      if (val(u.x,0,0,0) != HUGE) {
 coord f = { scale*val(u.x,0,0,0), scale*val(u.y,0,0,0) };
 glvertex2d (view, x + f.x - (f.x - f.y/2.)/5.,
      y + f.y - (f.x/2. + f.y)/5.);
 glvertex2d (view, x + f.x, y + f.y);
 glvertex2d (view, x + f.x, y + f.y);
 glvertex2d (view, x + f.x - (f.x + f.y/2.)/5.,
      y + f.y + (f.x/2. - f.y)/5.);
 glvertex2d (view, x, y);
 glvertex2d (view, x + f.x, y + f.y);
 view->ni++;
      }end_foreach_visible();}
    vertex_buffer_glEnd();
  }end_draw_lines();}
  view->res = res;



  {end_tracing("vectors","/home/gongminjiang/basilisk/src/draw.h",1063);return true;}
end_tracing("vectors","/home/gongminjiang/basilisk/src/draw.h",1064);}
#line 1084 "/home/gongminjiang/basilisk/src/draw.h"
struct _squares {
  char * color;
  char * z;
  double min, max, spread;
  bool linear;
  colormap map;
  float fc[3], lc[3];
  bool expr;

  coord n;
  double alpha;
};

     
bool squares (struct _squares p)
{tracing("squares","/home/gongminjiang/basilisk/src/draw.h",1098);

  scalar Z = {-1};
  vector n;
  bool zexpr = false;
  if (p.z) {
    Z = compile_expression (p.z, &zexpr);
    if (Z.i < 0)
      {end_tracing("squares","/home/gongminjiang/basilisk/src/draw.h",1107);return false;}
    n = new_vector("n");
    foreach_stencil()
      {
        {_stencil_val_a(n.x,0,0,0);_stencil_val(Z,1,0,0); _stencil_val(Z,-1,0,0);   }
        
#line 1111
{_stencil_val_a(n.y,0,0,0);_stencil_val(Z,0,1,0); _stencil_val(Z,0,-1,0);   }}end_foreach_stencil();
    {
#line 1109
foreach()
      {
        val(n.x,0,0,0) = (val(Z,1,0,0) - val(Z,-1,0,0))/(2.*Delta_x);
        
#line 1111
val(n.y,0,0,0) = (val(Z,0,1,0) - val(Z,0,-1,0))/(2.*Delta_y);}end_foreach();}
  }

  scalar col = {-1}; if (p.color && strcmp (p.color, "level")) { col = compile_expression (p.color, &p.expr); if (col.i < 0) {end_tracing("squares","/home/gongminjiang/basilisk/src/draw.h",1114);return false;} boundary_internal ((scalar *)((scalar[]){col,{-1}}), "/home/gongminjiang/basilisk/src/draw.h", 1114); } double cmap[127][3]; if (p.color) { if (p.min == 0 && p.max == 0) { if (col.i < 0) p.min = 0, p.max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (p.spread < 0.) p.min = s.min, p.max = s.max; else { if (!p.spread) p.spread = 5.; double spread = p.spread*s.stddev; p.min = avg - spread; p.max = avg + spread; } } } if (!p.map) p.map = jet; p.map (cmap); } if ((2 > 2 || p.linear) && !p.fc[0] && !p.fc[1] && !p.fc[2]) p.fc[0] = p.fc[1] = p.fc[2] = 1.;;
  scalar f = col;

  bview * view = draw();
  glShadeModel (0x1D01);
  if (p.linear) {
    {begin_colorized (p.fc, !VertexBuffer.color || !p.color, cmap, !VertexBuffer.color && p.color && p.linear && col.i >= 0); {

      if (Z.i < 0) {
 vertex_buffer_glNormal3d (0, 0, view->reversed ? -1 : 1);
 foreach_visible_stencil (view)
   {_stencil_val(f,0,0,0); { 
     
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) {                  _stencil_val(f,1,-1,0); _stencil_val(f,-1,1,0); _stencil_val(f,1,1,0); _stencil_val(f,-1,-1,0); _stencil_val(f,0,-1,0); _stencil_val(f,0,1,0); _stencil_val(f,-1,0,0);_stencil_val(f,1,0,0);_stencil_val(f,0,0,0);     } else {              _stencil_val(f,1,-1,0); _stencil_val(f,-1,1,0); _stencil_val(f,1,1,0); _stencil_val(f,-1,-1,0); _stencil_val(f,0,-1,0); _stencil_val(f,0,1,0); _stencil_val(f,-1,0,0);_stencil_val(f,1,0,0);_stencil_val(f,0,0,0);                } }   



      
     
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) {            _stencil_val(f,0,-1,0); _stencil_val(f,-1,-1,0); _stencil_val(f,-1,0,0);_stencil_val(f,0,0,0);     } else {        _stencil_val(f,0,-1,0); _stencil_val(f,-1,-1,0); _stencil_val(f,-1,0,0);_stencil_val(f,0,0,0);                } }       
     
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) {            _stencil_val(f,0,-1,0); _stencil_val(f,1,-1,0); _stencil_val(f,1,0,0);_stencil_val(f,0,0,0);     } else {        _stencil_val(f,0,-1,0); _stencil_val(f,1,-1,0); _stencil_val(f,1,0,0);_stencil_val(f,0,0,0);                } }       
     
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) {            _stencil_val(f,0,1,0); _stencil_val(f,1,1,0); _stencil_val(f,1,0,0);_stencil_val(f,0,0,0);     } else {        _stencil_val(f,0,1,0); _stencil_val(f,1,1,0); _stencil_val(f,1,0,0);_stencil_val(f,0,0,0);                } }       
     
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) {            _stencil_val(f,0,1,0); _stencil_val(f,-1,1,0); _stencil_val(f,-1,0,0);_stencil_val(f,0,0,0);     } else {        _stencil_val(f,0,1,0); _stencil_val(f,-1,1,0); _stencil_val(f,-1,0,0);_stencil_val(f,0,0,0);                } }       
     
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) {            _stencil_val(f,0,-1,0); _stencil_val(f,-1,-1,0); _stencil_val(f,-1,0,0);_stencil_val(f,0,0,0);     } else {        _stencil_val(f,0,-1,0); _stencil_val(f,-1,-1,0); _stencil_val(f,-1,0,0);_stencil_val(f,0,0,0);                } }       
     
     
     
   }   }end_foreach_visible_stencil();
 {
#line 1124
foreach_visible (view)
   if (val(f,0,0,0) != HUGE) {
     vertex_buffer_glBegin (0x0006);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (4.*val(f,0,0,0) + 2.*(val(f,1,0,0) + val(f,-1,0,0) + val(f,0,1,0) + val(f,0,-1,0)) + val(f,-1,-1,0) + val(f,1,1,0) + val(f,-1,1,0) + val(f,1,-1,0))/16., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (4.*val(f,0,0,0) + 2.*(val(f,1,0,0) + val(f,-1,0,0) + val(f,0,1,0) + val(f,0,-1,0)) + val(f,-1,-1,0) + val(f,1,1,0) + val(f,-1,1,0) + val(f,1,-1,0))/16.; if (p.max > p.min) glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); else glTexCoord1d (0.); } }



      ;
     glvertex2d (view, x, y);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4.; if (p.max > p.min) glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,1,0,0) + val(f,1,-1,0) + val(f,0,-1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,1,0,0) + val(f,1,-1,0) + val(f,0,-1,0))/4.; if (p.max > p.min) glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,1,0,0) + val(f,1,1,0) + val(f,0,1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,1,0,0) + val(f,1,1,0) + val(f,0,1,0))/4.; if (p.max > p.min) glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,1,0) + val(f,0,1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,1,0) + val(f,0,1,0))/4.; if (p.max > p.min) glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4.; if (p.max > p.min) glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
     vertex_buffer_glEnd();
     view->ni++;
   }end_foreach_visible();}
      }
      else
 {foreach_leaf()
   if (val(f,0,0,0) != HUGE) {
     vertex_buffer_glBegin (0x0006);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (4.*val(f,0,0,0) + 2.*(val(f,1,0,0) + val(f,-1,0,0) + val(f,0,1,0) + val(f,0,-1,0)) + val(f,-1,-1,0) + val(f,1,1,0) + val(f,-1,1,0) + val(f,1,-1,0))/16., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (4.*val(f,0,0,0) + 2.*(val(f,1,0,0) + val(f,-1,0,0) + val(f,0,1,0) + val(f,0,-1,0)) + val(f,-1,-1,0) + val(f,1,1,0) + val(f,-1,1,0) + val(f,1,-1,0))/16.; if (p.max > p.min) glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); else glTexCoord1d (0.); } }

                                                    ;
     glvertex_normal3d (view, point, n, x, y, val(Z,0,0,0));
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4.; if (p.max > p.min) glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex_normal3d (view, point, n, x - Delta_x/2., y - Delta_y/2.,
          (val(Z,0,0,0) + val(Z,-1,0,0) + val(Z,-1,-1,0) + val(Z,0,-1,0))/4.);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,1,0,0) + val(f,1,-1,0) + val(f,0,-1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,1,0,0) + val(f,1,-1,0) + val(f,0,-1,0))/4.; if (p.max > p.min) glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex_normal3d (view, point, n, x + Delta_x/2., y - Delta_y/2.,
          (val(Z,0,0,0) + val(Z,1,0,0) + val(Z,1,-1,0) + val(Z,0,-1,0))/4.);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,1,0,0) + val(f,1,1,0) + val(f,0,1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,1,0,0) + val(f,1,1,0) + val(f,0,1,0))/4.; if (p.max > p.min) glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex_normal3d (view, point, n, x + Delta_x/2., y + Delta_y/2.,
          (val(Z,0,0,0) + val(Z,1,0,0) + val(Z,1,1,0) + val(Z,0,1,0))/4.);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,1,0) + val(f,0,1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,1,0) + val(f,0,1,0))/4.; if (p.max > p.min) glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex_normal3d (view, point, n, x - Delta_x/2., y + Delta_y/2.,
          (val(Z,0,0,0) + val(Z,-1,0,0) + val(Z,-1,1,0) + val(Z,0,1,0))/4.);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4.; if (p.max > p.min) glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex_normal3d (view, point, n, x - Delta_x/2., y - Delta_y/2.,
          (val(Z,0,0,0) + val(Z,-1,0,0) + val(Z,-1,-1,0) + val(Z,0,-1,0))/4.);
     vertex_buffer_glEnd();
     view->ni++;
   }end_foreach_leaf();}
#line 1199 "/home/gongminjiang/basilisk/src/draw.h"
    }end_colorized();}
  }
  else {

    vertex_buffer_glNormal3d (0, 0, view->reversed ? -1 : 1);
    vertex_buffer_glBegin (0x0007);
    foreach_visible_stencil (view)
      {_stencil_val(f,0,0,0); {
 if (p.color && (!p.linear || col.i < 0)) {               _stencil_val(col,0,0,0);     }       
 
 if (p.color && (!p.linear || col.i < 0)) {               _stencil_val(col,0,0,0);     }       
 
 if (p.color && (!p.linear || col.i < 0)) {               _stencil_val(col,0,0,0);     }       
 
 if (p.color && (!p.linear || col.i < 0)) {               _stencil_val(col,0,0,0);     }       
 
 
      }   }end_foreach_visible_stencil();
    {
#line 1205
foreach_visible (view)
      if (val(f,0,0,0) != HUGE) {
 if (p.color && (!p.linear || col.i < 0)) { color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
 glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
 if (p.color && (!p.linear || col.i < 0)) { color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
 glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
 if (p.color && (!p.linear || col.i < 0)) { color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
 glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
 if (p.color && (!p.linear || col.i < 0)) { color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
 glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
 view->ni++;
      }end_foreach_visible();}
    vertex_buffer_glEnd();
#line 1235 "/home/gongminjiang/basilisk/src/draw.h"
  }
  if (p.expr) delete (((scalar[]){col,{-1}}));

  if (zexpr) delete (((scalar[]){Z,{-1}}));
  if (p.z) delete ((scalar *)((vector[]){n,{{-1},{-1}}}));

  {end_tracing("squares","/home/gongminjiang/basilisk/src/draw.h",1241);return true;}
end_tracing("squares","/home/gongminjiang/basilisk/src/draw.h",1242);}
#line 1253 "/home/gongminjiang/basilisk/src/draw.h"
struct _box {
  bool notics;
  float lc[3], lw;
};

     
bool box (struct _box p)
{tracing("box","/home/gongminjiang/basilisk/src/draw.h",1259);
  bview * view = draw();
  {begin_draw_lines (view, p.lc, p.lw); {

    float height = 0.5*gl_StrokeHeight();
    float width = gl_StrokeWidth ('1'), scale = L0/(60.*width), length;
    float Z1 = 2 == 2 ? 0. : Z0;
    char label[80];

    glMatrixMode (0x1700);

    if (!p.notics) {
      int nt = 8;
      for (int i = 0; i <= nt; i++) {
 glPushMatrix();
 glTranslatef (X0 + i*L0/nt - height/2.*scale, Y0 - width/3.*scale, Z1);
 glRotatef (-90, 0, 0, 1);
 glScalef (scale, scale, 1.);
 sprintf (label, "%g", X0 + i*L0/nt);
 gl_StrokeString (label);
 glPopMatrix();

 glPushMatrix();
 sprintf (label, "%g", Y0 + i*L0/nt);
 length = gl_StrokeLength (label);
 glTranslatef (X0 - (length + width/3.)*scale,
        Y0 + i*L0/nt - height/2.*scale, Z1);
 glScalef (scale, scale, 1.);
 gl_StrokeString (label);
 glPopMatrix();
#line 1302 "/home/gongminjiang/basilisk/src/draw.h"
      }

      glPushMatrix();
      sprintf (label, "%g", X0 + L0/2.);
      length = gl_StrokeLength (label);
      glTranslatef (X0 + L0/2 - height*scale, Y0 - (length + 4.*width)*scale, Z1);
      glScalef (2.*scale, 2.*scale, 1.);
      gl_StrokeString ("X");
      glPopMatrix();


      glPushMatrix();
      sprintf (label, "%g", Y0 + L0/2.);
      length = gl_StrokeLength (label);
      glTranslatef (X0 - (length + 4.*width)*scale,
      Y0 + L0/2. - height*scale, Z1);
      glScalef (2.*scale, 2.*scale, 1.);
      gl_StrokeString ("Y");
      glPopMatrix();
#line 1333 "/home/gongminjiang/basilisk/src/draw.h"
    }


    {foreach_level (0) {
      vertex_buffer_glBegin (0x0002);
      glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
      glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
      glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
      glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
      vertex_buffer_glEnd ();
      view->ni++;
    }end_foreach_level();}
#line 1365 "/home/gongminjiang/basilisk/src/draw.h"
  }end_draw_lines();}
  {end_tracing("box","/home/gongminjiang/basilisk/src/draw.h",1366);return true;}
end_tracing("box","/home/gongminjiang/basilisk/src/draw.h",1367);}
#line 1380 "/home/gongminjiang/basilisk/src/draw.h"
struct _isosurface {
  char * f;
  double v;

  char * color;
  double min, max, spread;
  bool linear;
  colormap map;
  float fc[3], lc[3], lw;
  bool expr;
};

     
bool isosurface (struct _isosurface p)
{tracing("isosurface","/home/gongminjiang/basilisk/src/draw.h",1393);
#line 1454 "/home/gongminjiang/basilisk/src/draw.h"
  {end_tracing("isosurface","/home/gongminjiang/basilisk/src/draw.h",1454);return true;}
end_tracing("isosurface","/home/gongminjiang/basilisk/src/draw.h",1455);}
#line 1465 "/home/gongminjiang/basilisk/src/draw.h"
struct _travelling {
  double start, end;
  float tx, ty, quat[4], fov;
};




void travelling (struct _travelling p)
{
  static float tx, ty, quat[4], fov;
  static double told = -1.;
  if (told < p.start && t >= p.start) {
    bview * view = get_view();
    tx = view->tx, ty = view->ty, fov = view->fov;
    for (int i = 0; i < 4; i++)
      quat[i] = view->quat[i];
  }
  if (t >= p.start && t <= p.end)
    view ((struct _view_set){.tx = (!p.tx ? tx : ((t - p.start)*(p.tx) + (p.end - t)*(tx))/(p.end - p.start)), .ty = (!p.ty ? ty : ((t - p.start)*(p.ty) + (p.end - t)*(ty))/(p.end - p.start)),
   .fov = (!p.fov ? fov : ((t - p.start)*(p.fov) + (p.end - t)*(fov))/(p.end - p.start)),
   .quat = {(!p.quat[0] ? quat[0] : ((t - p.start)*(p.quat[0]) + (p.end - t)*(quat[0]))/(p.end - p.start)), (!p.quat[1] ? quat[1] : ((t - p.start)*(p.quat[1]) + (p.end - t)*(quat[1]))/(p.end - p.start)),
           (!p.quat[2] ? quat[2] : ((t - p.start)*(p.quat[2]) + (p.end - t)*(quat[2]))/(p.end - p.start)), (!p.quat[3] ? quat[3] : ((t - p.start)*(p.quat[3]) + (p.end - t)*(quat[3]))/(p.end - p.start))}});
  if (told < p.end && t >= p.end) {
    bview * view = get_view();
    tx = view->tx, ty = view->ty, fov = view->fov;
    for (int i = 0; i < 4; i++)
      quat[i] = view->quat[i];
  }
  told = t;
}
#line 1512 "/home/gongminjiang/basilisk/src/draw.h"
struct _draw_string {
  char * str;
  int pos;
  float size;
  float lc[3], lw;
};

     
bool draw_string (struct _draw_string p)
{tracing("draw_string","/home/gongminjiang/basilisk/src/draw.h",1520);
  bview * view = draw();

  glMatrixMode (0x1701);
  glPushMatrix();
  glLoadIdentity();

  glMatrixMode (0x1700);
  glPushMatrix();
  glLoadIdentity();

  vertex_buffer_glColor3f (p.lc[0], p.lc[1], p.lc[2]);
  glLineWidth (view->samples*(p.lw > 0. ? p.lw : 1.));

  float width = gl_StrokeWidth ('1'), height = gl_StrokeHeight();
  if (!p.size)
    p.size = 40;
  float hscale = 2./(p.size*width), vscale = hscale*view->width/view->height;
  float vmargin = width/2.*vscale;
  if (p.pos == 0)
    glTranslatef (-1., -1. + vmargin, 0.);
  else if (p.pos == 1)
    glTranslatef (-1., 1. - height*vscale, 0.);
  else if (p.pos == 2)
    glTranslatef (1. - strlen(p.str)*width*hscale, 1. - height*vscale, 0.);
  else
    glTranslatef (1. - strlen(p.str)*width*hscale, -1. + vmargin, 0.);
  glScalef (hscale, vscale, 1.);
  gl_StrokeString (p.str);

  glMatrixMode (0x1700);
  glPopMatrix();
  glMatrixMode (0x1701);
  glPopMatrix();

  {end_tracing("draw_string","/home/gongminjiang/basilisk/src/draw.h",1556);return true;}
end_tracing("draw_string","/home/gongminjiang/basilisk/src/draw.h",1557);}




struct _labels {
  char * f;
  float lc[3], lw;
};

     
bool labels (struct _labels p)
{tracing("labels","/home/gongminjiang/basilisk/src/draw.h",1568);

  bool expr = false;
  scalar f = compile_expression (p.f, &expr);
  if (f.i < 0)
    {end_tracing("labels","/home/gongminjiang/basilisk/src/draw.h",1574);return false;}
  bview * view = draw();
  float width = gl_StrokeWidth ('1'), height = gl_StrokeHeight();
  float res = view->res;
  if (view->res < 150*view->samples)
    view->res = 150*view->samples;
  {begin_draw_lines (view, p.lc, p.lw); {
    glMatrixMode (0x1700);
    foreach_visible_stencil (view)
      {_stencil_val(f,0,0,0); { 
 
  
_stencil_val(f,0,0,0);     
 
            
 
 
 
 
      
#line 1592
}   }end_foreach_visible_stencil();
    {
#line 1582
foreach_visible (view)
      if (val(f,0,0,0) != HUGE) {
 glPushMatrix();
 char s[80];
 sprintf (s, "%g", val(f,0,0,0));
 float scale = 0.8*Delta_x/(strlen(s)*width);
 glTranslatef (x - 0.4*Delta_x, y - scale*height/3., 0.);
 glScalef (scale, scale, 1.);
 gl_StrokeString (s);
 glPopMatrix();
      }end_foreach_visible();}
  }end_draw_lines();}
  view->res = res;
  if (expr) delete (((scalar[]){f,{-1}}));
  {end_tracing("labels","/home/gongminjiang/basilisk/src/draw.h",1596);return true;}




end_tracing("labels","/home/gongminjiang/basilisk/src/draw.h",1601);}







#line 1 "draw_json.h"
#line 1 "/home/gongminjiang/basilisk/src/draw_json.h"



int _view_set_json (void * q, char * s, int len) {
  struct _view_set * p = (struct _view_set *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"view_set\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"tx\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->tx);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"ty\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->ty);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fov\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->fov);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"quat\": { \"type\": \"pfloat\", \"cardinality\": 4, \"value\": [%f,%f,%f,%f] }", p->quat[0], p->quat[1], p->quat[2], p->quat[3]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"sx\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->sx);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"sy\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->sy);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"sz\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->sz);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"width\": { \"type\": \"punsigned\", \"cardinality\": 1, \"value\": \"%d\" }", p->width);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"height\": { \"type\": \"punsigned\", \"cardinality\": 1, \"value\": \"%d\" }", p->height);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"samples\": { \"type\": \"punsigned\", \"cardinality\": 1, \"value\": \"%d\" }", p->samples);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"bg\": { \"type\": \"pfloat\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->bg[0], p->bg[1], p->bg[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"theta\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->theta);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"phi\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->phi);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"psi\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->psi);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"relative\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->relative);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"tz\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->tz);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"near\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->near);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"far\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->far);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"res\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->res);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"camera\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->camera);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"cache\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", p->cache);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p1x\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->p1x);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p1y\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->p1y);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p2x\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->p2x);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p2y\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->p2y);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _translate_json (void * q, char * s, int len) {
  struct _translate * p = (struct _translate *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"translate\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"x\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->x);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"y\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->y);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"z\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->z);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _mirror_json (void * q, char * s, int len) {
  struct _mirror * p = (struct _mirror *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"mirror\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"n\": { \"type\": \"coord\", \"cardinality\": 3, \"value\": [%lf,%lf,%lf] }", p->n.x, p->n.y, p->n.z);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"alpha\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->alpha);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _draw_vof_json (void * q, char * s, int len) {
  struct _draw_vof * p = (struct _draw_vof *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"draw_vof\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"c\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->c);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"s\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->s);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"edges\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->edges);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"larger\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->larger);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"filled\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", p->filled);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->color);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->min);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->max);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->spread);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->linear);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->fc[0], p->fc[1], p->fc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _isoline_json (void * q, char * s, int len) {
  struct _isoline * p = (struct _isoline *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"isoline\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"phi\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->phi);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"val\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->val);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"n\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", p->n);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"c\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->c);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"s\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->s);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"edges\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->edges);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"larger\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->larger);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"filled\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", p->filled);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->color);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->min);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->max);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->spread);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->linear);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->fc[0], p->fc[1], p->fc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _cells_json (void * q, char * s, int len) {
  struct _cells * p = (struct _cells *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"cells\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"n\": { \"type\": \"coord\", \"cardinality\": 3, \"value\": [%lf,%lf,%lf] }", p->n.x, p->n.y, p->n.z);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"alpha\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->alpha);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _vectors_json (void * q, char * s, int len) {
  struct _vectors * p = (struct _vectors *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"vectors\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"u\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->u);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"scale\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->scale);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _squares_json (void * q, char * s, int len) {
  struct _squares * p = (struct _squares *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"squares\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->color);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"z\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->z);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->min);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->max);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->spread);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->linear);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->fc[0], p->fc[1], p->fc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"n\": { \"type\": \"coord\", \"cardinality\": 3, \"value\": [%lf,%lf,%lf] }", p->n.x, p->n.y, p->n.z);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"alpha\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->alpha);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _box_json (void * q, char * s, int len) {
  struct _box * p = (struct _box *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"box\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"notics\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->notics);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _isosurface_json (void * q, char * s, int len) {
  struct _isosurface * p = (struct _isosurface *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"isosurface\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"f\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->f);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"v\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->v);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->color);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->min);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->max);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->spread);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->linear);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->fc[0], p->fc[1], p->fc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _travelling_json (void * q, char * s, int len) {
  struct _travelling * p = (struct _travelling *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"travelling\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"start\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->start);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"end\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->end);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"tx\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->tx);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"ty\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->ty);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"quat\": { \"type\": \"pfloat\", \"cardinality\": 4, \"value\": [%f,%f,%f,%f] }", p->quat[0], p->quat[1], p->quat[2], p->quat[3]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fov\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->fov);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _draw_string_json (void * q, char * s, int len) {
  struct _draw_string * p = (struct _draw_string *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"draw_string\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"str\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->str);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"pos\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", p->pos);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"size\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->size);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _labels_json (void * q, char * s, int len) {
  struct _labels * p = (struct _labels *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"labels\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"f\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->f);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
#line 1610 "/home/gongminjiang/basilisk/src/draw.h"

struct {
  int (* json) (void * q, char * s, int len);
} bview_interface[] = {
  { _draw_vof_json },
  { _squares_json },
  { _cells_json },
  { _box_json },

  { _isoline_json },
  { _labels_json },
  { _vectors_json },



  { NULL }
};
#line 431 "/home/gongminjiang/basilisk/src/view.h"
#line 446 "/home/gongminjiang/basilisk/src/view.h"
struct _load {
  FILE * fp;
  char * file;
  Array * buf;
};

bool load (struct _load p);
#line 500 "/home/gongminjiang/basilisk/src/view.h"
struct _save {
  char * file, * format, * opt;
  FILE * fp;
  float lw;
  int sort, options;

  bview * view;
};

static void bview_draw (bview * view)
{
  if (!view->active)
    return;
  view->active = false;
  glFinish ();
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}

     
bool save (struct _save p)
{tracing("save","/home/gongminjiang/basilisk/src/view.h",519);
  char ppm[] = "ppm";
  if (!p.format) {
    p.format = ppm;
    if (p.file) {
      char * s = strchr (p.file, '.'), * dot = s;
      while (s) {
 dot = s;
 s = strchr (s + 1, '.');
      }
      if (dot)
 p.format = dot + 1;
    }
  }

  bview * view = p.view ? p.view : get_view();

  if (!strcmp (p.format, "png") ||
      !strcmp (p.format, "jpg") ||
      (p.file && is_animation (p.file))) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0) {
      FILE * fp = open_image (p.file, p.opt);
      if (!fp) {
 perror (p.file);
 {end_tracing("save","/home/gongminjiang/basilisk/src/view.h",546);return false;}
      }
      gl_write_image (fp, image, view->width, view->height, view->samples);
      close_image (p.file, fp);
    }
    {end_tracing("save","/home/gongminjiang/basilisk/src/view.h",551);return true;}
  }

  if (p.file && (p.fp = fopen (p.file, "w")) == NULL) {
    perror (p.file);
    {end_tracing("save","/home/gongminjiang/basilisk/src/view.h",556);return false;}
  }
  if (!p.fp)
    p.fp = fout;

  if (!strcmp (p.format, "ppm")) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0)
      gl_write_image (p.fp, image, view->width, view->height, view->samples);
  }

  else if (!strcmp (p.format, "bv")) {

    fprintf (ferr, "save(): error: the '%s' format is no longer supported\n",
      p.format);
#line 588 "/home/gongminjiang/basilisk/src/view.h"
  }

  else if (!strcmp (p.format, "gnu") ||
    !strcmp (p.format, "obj") ||
    !strcmp (p.format, "kml") ||
    !strcmp (p.format, "ps") ||
    !strcmp (p.format, "eps") ||
    !strcmp (p.format, "tex") ||
    !strcmp (p.format, "pdf") ||
    !strcmp (p.format, "svg") ||
    !strcmp (p.format, "pgf"))
    fprintf (ferr, "save(): error: the '%s' format is no longer supported\n",
      p.format);

  else {
    fprintf (ferr, "save(): unknown format '%s'\n", p.format);
    if (p.file) {
      fclose (p.fp);
      remove (p.file);
    }
    {end_tracing("save","/home/gongminjiang/basilisk/src/view.h",608);return false;}
  }

  fflush (p.fp);
  if (p.file)
    fclose (p.fp);

  {end_tracing("save","/home/gongminjiang/basilisk/src/view.h",615);return true;}
end_tracing("save","/home/gongminjiang/basilisk/src/view.h",616);}







static char * remove_blanks (char * line)
{
  while (strchr (" \t", *line)) line++;
  char * s = line, * cur = line;
  bool instring = false;
  while (*s != '\0' && *s != '#') {
    if (*s == '"')
      instring = !instring;
    if (instring || !strchr (" \t", *s))
      *cur++ = *s;
    s++;
  }
  *cur = '\0';
  return line;
}






#line 1 "draw_get.h"
#line 1 "/home/gongminjiang/basilisk/src/draw_get.h"

#line 1 "parse.h"
#line 1 "/home/gongminjiang/basilisk/src/parse.h"




enum ParamsType { pstring, pint, punsigned, pbool, pfloat, pdouble, pcolormap };

typedef struct {
  char * key;
  enum ParamsType type;
  void * val;
  int n;
} Params;

static bool atobool (char * s)
{
  if (!strcmp (s, "true"))
    return true;
  if (!strcmp (s, "false"))
    return false;
  return atoi (s) != 0;
}

static bool args (Params * p, char * val)
{
  static char * name[] = { "string", "int", "unsigned",
      "bool", "float", "double", "colormap" };
  switch (p->type) {

  case pstring:
    if (val[0] != '"') {
      fprintf (ferr, "expecting a string for '%s' got '%s'\n", p->key, val);
      return false;
    }
    if (val[strlen(val) - 1] != '"') {
      fprintf (ferr, "unterminated quoted string '%s'\n", val);
      return false;
    }
    val[strlen(val) - 1] = '\0';
    char * s = &val[1];
    int nc = 0;
    while (*s != '\0') {
      if (!strchr (" \t\n\r", *s))
 nc++;
      s++;
    }
    *((char **)p->val) = nc > 0 ? &val[1] : NULL;
    break;

  case pcolormap:
    if (!strcmp (val, "jet"))
      *((colormap *)p->val) = jet;
    else if (!strcmp (val, "cool_warm"))
      *((colormap *)p->val) = cool_warm;
    else if (!strcmp (val, "gray"))
      *((colormap *)p->val) = gray;
    else if (!strcmp (val, "randomap"))
      *((colormap *)p->val) = randomap;
    else {
      fprintf (ferr, "unknown colormap '%s'\n", val);
      return false;
    }
    break;

  case pint: case punsigned: case pbool: case pdouble: case pfloat:
    if (val[0] == '"') {
      fprintf (ferr, "expecting a %s for '%s' got %s\n",
        name[p->type], p->key, val);
      return false;
    }
    if (!p->n) {
      switch (p->type) {
      case pint: *((int *)p->val) = atoi(val); break;
      case punsigned: *((unsigned *)p->val) = atoi(val); break;
      case pbool: *((bool *)p->val) = atobool(val); break;
      case pfloat: *((float *)p->val) = atof(val); break;
      case pdouble: *((double *)p->val) = atof(val); break;
      default: if (!(false)) qassert ("/home/gongminjiang/basilisk/src/parse.h", 77, "false");
      }
    }
    else {
      if (val[0] != '{') {
 fprintf (ferr, "expecting an array for '%s' got %s\n", p->key, val);
 return false;
      }
      val++;
      int i = 0;
      char c = ',';
      while (i < p->n && c != '}') {
 char * s = strchr (val, ',');
 if (!s)
   s = strchr (val, '}');
 if (!s) {
   fprintf (ferr, "expecting an array for '%s' got %s\n", p->key, val);
   return false;
 }
 c = *s;
 *s++ = '\0';
 switch (p->type) {
 case pint: ((int *)p->val)[i++] = atoi (val); break;
 case punsigned: ((unsigned *)p->val)[i++] = atoi (val); break;
 case pbool: ((bool *)p->val)[i++] = atobool (val); break;
 case pfloat: ((float *)p->val)[i++] = atof (val); break;
 case pdouble: ((double *)p->val)[i++] = atof (val); break;
 default: if (!(false)) qassert ("/home/gongminjiang/basilisk/src/parse.h", 104, "false");
 }
 val = s;
      }
      if (c != '}') {
 fprintf (ferr, "expecting '}' for '%s' got %s\n", p->key, val);
 return false;
      }
    }
    break;

  default:
    if (!(false)) qassert ("/home/gongminjiang/basilisk/src/parse.h", 116, "false");
  }
  return true;
}

static char * find_comma (char * s)
{
  int par = 0;
  while (*s != '\0') {
    if (*s == ',' && par == 0) {
      *s = '\0';
      return s + 1;
    }
    if (*s == '{')
      par++;
    else if (*s == '}')
      par--;
    s++;
  }
  return NULL;
}

static char * mystrtok (char * str, const char * delim)
{
  static char * s = NULL;
  char * start = str ? str : s;
  bool string = false;
  s = start;
  while (*s != '\0') {
    if (*s == '"')
      string = !string;
    if (!string && strchr(delim, *s))
      break;
    s++;
  }
  if (*s != '\0')
    *s++ = '\0';
  return start;
}

int parse_params (Params * params)
{
  char * s;
  int i = 0, n = 0;
  Params * p = params;
  while (p->key) p++, n++;
  if (!(s = mystrtok (NULL, ");")) || s[0] == '\n')
    return false;
  while (s && *s != '\0') {
    char * next = find_comma (s), * key = s;
    if ((s = strchr (key, '='))) {
      s[0] = '\0', s++;
      i = -1;
      Params * p = params;
      while (p->key && strcmp(p->key, key)) p++;
      if (!p->key) {
 fprintf (ferr, "unknown key '%s'\n", key);
 return false;
      }
      if (!args (p, s))
 return false;
    }
    else {
      if (i < 0) {
 fprintf (ferr, "anonymous value '%s' after keys\n", key);
 return false;
      }
      if (i >= n) {
 fprintf (ferr, "too many parameters: '%s' %d %d\n", key, i, n);
 return false;
      }
      if (!args (&params[i], key))
 return false;
      i++;
    }
    s = next;
  }
  return true;
}
#line 3 "/home/gongminjiang/basilisk/src/draw_get.h"

bool _view_set_get (struct _view_set * p) {
  Params params[] = {
    {"tx", pfloat, &p->tx},
    {"ty", pfloat, &p->ty},
    {"fov", pfloat, &p->fov},
    {"quat", pfloat, p->quat, 4},
    {"sx", pfloat, &p->sx},
    {"sy", pfloat, &p->sy},
    {"sz", pfloat, &p->sz},
    {"width", punsigned, &p->width},
    {"height", punsigned, &p->height},
    {"samples", punsigned, &p->samples},
    {"bg", pfloat, p->bg, 3},
    {"theta", pfloat, &p->theta},
    {"phi", pfloat, &p->phi},
    {"psi", pfloat, &p->psi},
    {"relative", pbool, &p->relative},
    {"tz", pfloat, &p->tz},
    {"near", pfloat, &p->near},
    {"far", pfloat, &p->far},
    {"res", pfloat, &p->res},
    {"camera", pstring, &p->camera},
    {"cache", pint, &p->cache},
    {"p1x", pfloat, &p->p1x},
    {"p1y", pfloat, &p->p1y},
    {"p2x", pfloat, &p->p2x},
    {"p2y", pfloat, &p->p2y},
    {NULL}
  };
  return parse_params (params);
}

bool _translate_get (struct _translate * p) {
  Params params[] = {
    {"x", pfloat, &p->x},
    {"y", pfloat, &p->y},
    {"z", pfloat, &p->z},
    {NULL}
  };
  return parse_params (params);
}

bool _mirror_get (struct _mirror * p) {
  Params params[] = {
    {"n", pdouble, &p->n, 3},
    {"alpha", pdouble, &p->alpha},
    {NULL}
  };
  return parse_params (params);
}

bool _draw_vof_get (struct _draw_vof * p) {
  Params params[] = {
    {"c", pstring, &p->c},
    {"s", pstring, &p->s},
    {"edges", pbool, &p->edges},
    {"larger", pdouble, &p->larger},
    {"filled", pint, &p->filled},
    {"color", pstring, &p->color},
    {"min", pdouble, &p->min},
    {"max", pdouble, &p->max},
    {"spread", pdouble, &p->spread},
    {"linear", pbool, &p->linear},
    {"map", pcolormap, &p->map},
    {"fc", pfloat, p->fc, 3},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {"expr", pbool, &p->expr},
    {NULL}
  };
  return parse_params (params);
}

bool _isoline_get (struct _isoline * p) {
  Params params[] = {
    {"phi", pstring, &p->phi},
    {"val", pdouble, &p->val},
    {"n", pint, &p->n},
    {"c", pstring, &p->c},
    {"s", pstring, &p->s},
    {"edges", pbool, &p->edges},
    {"larger", pdouble, &p->larger},
    {"filled", pint, &p->filled},
    {"color", pstring, &p->color},
    {"min", pdouble, &p->min},
    {"max", pdouble, &p->max},
    {"spread", pdouble, &p->spread},
    {"linear", pbool, &p->linear},
    {"map", pcolormap, &p->map},
    {"fc", pfloat, p->fc, 3},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {"expr", pbool, &p->expr},
    {NULL}
  };
  return parse_params (params);
}

bool _cells_get (struct _cells * p) {
  Params params[] = {
    {"n", pdouble, &p->n, 3},
    {"alpha", pdouble, &p->alpha},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}

bool _vectors_get (struct _vectors * p) {
  Params params[] = {
    {"u", pstring, &p->u},
    {"scale", pdouble, &p->scale},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}

bool _squares_get (struct _squares * p) {
  Params params[] = {
    {"color", pstring, &p->color},
    {"z", pstring, &p->z},
    {"min", pdouble, &p->min},
    {"max", pdouble, &p->max},
    {"spread", pdouble, &p->spread},
    {"linear", pbool, &p->linear},
    {"map", pcolormap, &p->map},
    {"fc", pfloat, p->fc, 3},
    {"lc", pfloat, p->lc, 3},
    {"expr", pbool, &p->expr},
    {"n", pdouble, &p->n, 3},
    {"alpha", pdouble, &p->alpha},
    {NULL}
  };
  return parse_params (params);
}

bool _box_get (struct _box * p) {
  Params params[] = {
    {"notics", pbool, &p->notics},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}

bool _isosurface_get (struct _isosurface * p) {
  Params params[] = {
    {"f", pstring, &p->f},
    {"v", pdouble, &p->v},
    {"color", pstring, &p->color},
    {"min", pdouble, &p->min},
    {"max", pdouble, &p->max},
    {"spread", pdouble, &p->spread},
    {"linear", pbool, &p->linear},
    {"map", pcolormap, &p->map},
    {"fc", pfloat, p->fc, 3},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {"expr", pbool, &p->expr},
    {NULL}
  };
  return parse_params (params);
}

bool _travelling_get (struct _travelling * p) {
  Params params[] = {
    {"start", pdouble, &p->start},
    {"end", pdouble, &p->end},
    {"tx", pfloat, &p->tx},
    {"ty", pfloat, &p->ty},
    {"quat", pfloat, p->quat, 4},
    {"fov", pfloat, &p->fov},
    {NULL}
  };
  return parse_params (params);
}

bool _draw_string_get (struct _draw_string * p) {
  Params params[] = {
    {"str", pstring, &p->str},
    {"pos", pint, &p->pos},
    {"size", pfloat, &p->size},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}

bool _labels_get (struct _labels * p) {
  Params params[] = {
    {"f", pstring, &p->f},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}
#line 646 "/home/gongminjiang/basilisk/src/view.h"

bool process_line (char * line)
{
  if (line[0] == '\0')
    return true;
  char * buf = pstrdup (line,__func__,__FILE__,__LINE__);
  char * s = mystrtok (remove_blanks (line), "(");
  if (!s) {
    pfree (buf,__func__,__FILE__,__LINE__);
    return true;
  }

  if (!strcmp (s, "restore")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file) {
      bview * view = get_view();
      if (view->cache) {
 free_cexpr (view->cache);
 view->cache = pcalloc (1, sizeof (cexpr),__func__,__FILE__,__LINE__);
      }
      if (!restore ((struct Dump){.file = file, .list = all}))
 fprintf (ferr, "could not restore from '%s'\n", file);
      else {
 restriction (all);
 fields_stats();
 clear();
      }
    }
  }

  else if (!strcmp (s, "dump")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    dump ((struct Dump){.file = file});
  }

  else if (!strcmp (s, "input_gfs")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file) {
      input_gfs ((struct OutputGfs){.file = file, .list = all});
      restriction (all);
      fields_stats();
      clear();
    }
  }

  else if (!strcmp (s, "save")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file)
      save ((struct _save){.file = file});
  }

  else if (!strcmp (s, "load")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file)
      load ((struct _load){.file = file});
  }

  else if (!strcmp (s, "cells")) {
    struct _cells p = {{0}};
    if (!_cells_get (&p) || !cells (p))
      return false;
  }

  else if (!strcmp (s, "vectors")) {
    struct _vectors p = {0};
    if (!_vectors_get (&p) || !vectors (p))
      return false;
  }

  else if (!strcmp (s, "draw_vof")) {
    struct _draw_vof p = {0};
    if (!_draw_vof_get (&p) || !draw_vof (p))
      return false;
  }

  else if (!strcmp (s, "isoline")) {
    struct _isoline p = {0};
    if (!_isoline_get (&p) || !isoline (p))
      return false;
  }

  else if (!strcmp (s, "squares")) {
    struct _squares p = {0};
    if (!_squares_get (&p) || !squares (p))
      return false;
  }

  else if (!strcmp (s, "begin_translate")) {
    struct _translate p = {0};
    _translate_get (&p);
    begin_translate (p);
  }

  else if (!strcmp (s, "end_translate"))
    end_translate();

  else if (!strcmp (s, "begin_mirror")) {
    struct _mirror p = {{0}};
    _mirror_get (&p);
    begin_mirror (p);
  }

  else if (!strcmp (s, "end_mirror"))
    end_mirror();

  else if (!strcmp (s, "squares")) {
    struct _squares p = {0};
    if (!_squares_get (&p) || !squares (p))
      return false;
  }

  else if (!strcmp (s, "isosurface")) {
    struct _isosurface p = {0};
    if (!_isosurface_get (&p) || !isosurface (p))
      return false;
  }

  else if (!strcmp (s, "draw_string")) {
    struct _draw_string p = {0};
    if (!_draw_string_get (&p) || !draw_string (p))
      return false;
  }

  else if (!strcmp (s, "labels")) {
    struct _labels p = {0};
    if (!_labels_get (&p) || !labels (p))
      return false;
  }

  else if (!strcmp (s, "clear"))
    clear();

  else if (!strcmp (s, "box")) {
    struct _box p = {0};
    if (!_box_get (&p) || !box (p))
      return false;
  }

  else if (!strcmp (s, "view")) {
    struct _view_set p = {0};
    _view_set_get (&p);
    view (p);
  }

  else if (s[0] != '\n' && s[0] != '\0')
    fprintf (ferr, "load(): syntax error: '%s'\n", s);

  pfree (buf,__func__,__FILE__,__LINE__);
  return true;
}

bool load (struct _load p) {
  if (p.file) {
    p.fp = fopen (p.file, "r");
    if (!p.fp) {
      perror (p.file);
      return false;
    }
  }

  if (p.fp) {
    char line[256];
    while (fgets (line, 256, p.fp) && process_line (line));
  }
  else if (p.buf) {
    int i = 0;
    char * s = (char *) p.buf->p;
    while (i < p.buf->len) {
      char * start = s;
      while (i < p.buf->len && *s != '\n')
 s++, i++;
      if (*s == '\n' && ++s > start) {
 char line[s - start + 1];
 strncpy (line, start, s - start);
 line[s - start] = '\0';
 process_line (line);
      }
    }
  }
  return true;
}
#line 19 "atomisation.c"
#line 33 "atomisation.c"
int maxlevel = 10;
double uemax = 0.1;







scalar  f0={12};
static double _boundary4(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( dirichlet(val(f0,0,0,0)*(1. + 0.05*sin (10.*2.*3.14159265358979*t))));}}static double _boundary4_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( dirichlet_homogeneous());}}
static double _boundary5(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( dirichlet(0));}}static double _boundary5_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( dirichlet_homogeneous());}}



static double _boundary6(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann(0));}}static double _boundary6_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann_homogeneous());}}
static double _boundary7(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( val(f0,0,0,0));}}

static double _boundary8(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann(0));}}static double _boundary8_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( neumann_homogeneous());}}
static double _boundary9(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( dirichlet(0));}}static double _boundary9_homogeneous(Point point,Point neighbor,scalar _s,void *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);POINT_VARIABLES;return( dirichlet_homogeneous());}}





int main (int argc, char * argv[])
{_init_solver();
  if (argc > 1)
    maxlevel = atoi (argv[1]);
  if (argc > 2)
    uemax = atof (argv[2]);





  init_grid (64);
  origin ((struct _origin){0, -1.5, -1.5});
  size (3.);





  rho1 = 1., rho2 = 1./27.84;
  mu1 = 2.*1./12./5800*rho1, mu2 = 2.*1./12./5800*rho2;
  _attribute[f.i].sigma = 3e-5;

  run();
free_solver();}




static int init_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t = 0);*ip=i;*tp=t;return ret;}      static int init_0(const int i,const double t,Event *_ev){tracing("init_0","atomisation.c",87); {
  if (!restore ((struct Dump){.file = "restart"})) {





    do { int refined; do { boundary_internal ((scalar *)all, "atomisation.c", 94); refined = 0; ((Tree *)grid)->refined.n = 0; {foreach_leaf() if (x < 1.2*0.025 && sq(y) + sq(z) < 2.*sq(1./12.) && level < maxlevel) { refine_cell (point, all, 0, &((Tree *)grid)->refined); refined++; continue; }end_foreach_leaf();} mpi_all_reduce (refined, MPI_INT, MPI_SUM); if (refined) { mpi_boundary_refine (all); mpi_boundary_update (all); } } while (refined); } while(0);





    do { scalar  phi=new_vertex_scalar("phi"); foreach_vertex_stencil() {_stencil_val_a(phi,0,0,0);      }end_foreach_vertex_stencil(); {foreach_vertex() val(phi,0,0,0) = sq(1./12.) - sq(y) - sq(z);end_foreach_vertex();} fractions ((struct Fractions){phi, f0});delete((scalar*)((scalar[]){phi,{-1}})); } while(0);
    _attribute[f0.i].refine = _attribute[f0.i].prolongation = fraction_refine;
    restriction (((scalar[]){f0,{-1}}));




    foreach_stencil() {
      _stencil_val_a(f,0,0,0); _stencil_val(f0,0,0,0);   
      _stencil_val_a(u.x,0,0,0); _stencil_val(f,0,0,0); 
    }end_foreach_stencil();




    {
#line 107
foreach() {
      val(f,0,0,0) = val(f0,0,0,0)*(x < 0.025);
      val(u.x,0,0,0) = val(f,0,0,0);
    }end_foreach();}
  }
}{end_tracing("init_0","atomisation.c",112);return 0;}end_tracing("init_0","atomisation.c",112);}






static int logfile_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int logfile(const int i,const double t,Event *_ev){tracing("logfile","atomisation.c",119); {
  if (i == 0)
    fprintf (ferr,
      "t dt mgp.i mgpf.i mgu.i grid->tn perf.t perf.speed\n");
  fprintf (ferr, "%g %g %d %d %d %ld %g %g\n",
    t, dt, mgp.i, mgpf.i, mgu.i,
    grid->tn, perf.t, perf.speed);
}{end_tracing("logfile","atomisation.c",126);return 0;}end_tracing("logfile","atomisation.c",126);}




static int movie_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t += 1e-2);*ip=i;*tp=t;return ret;}      static int movie(const int i,const double t,Event *_ev){tracing("movie","atomisation.c",131);
{

  scalar  omega=new_scalar("omega");
  vorticity (u, omega);
  view ((struct _view_set){.tx = -0.5});
  clear();
  draw_vof ((struct _draw_vof){"f"});
  squares ((struct _squares){"omega", .linear = true, .spread = 10});
  box ((struct _box){0});
#line 151 "atomisation.c"
  save ((struct _save){"movie.mp4"});delete((scalar*)((scalar[]){omega,{-1}}));
}{end_tracing("movie","atomisation.c",152);return 0;}end_tracing("movie","atomisation.c",152);}





static int snapshot_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=( t <= 3.8);*ip=i;*tp=t;return ret;}static int snapshot_expr1(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=( t += 0.1);*ip=i;*tp=t;return ret;}static int snapshot_expr2(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t = 0.1);*ip=i;*tp=t;return ret;}      static int snapshot(const int i,const double t,Event *_ev){tracing("snapshot","atomisation.c",158); {
  char name[80];
  sprintf (name, "snapshot-%g", t);
  scalar  pid=new_scalar("pid");
  foreach_stencil()
    {_stencil_val_a(pid,0,0,0);     }end_foreach_stencil();
  {
#line 162
foreach()
    val(pid,0,0,0) = fmod(pid()*(npe() + 37), npe());end_foreach();}
  dump ((struct Dump){name});delete((scalar*)((scalar[]){pid,{-1}}));
}{end_tracing("snapshot","atomisation.c",165);return 0;}end_tracing("snapshot","atomisation.c",165);}
#line 177 "atomisation.c"
static int droplets_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t += 0.1);*ip=i;*tp=t;return ret;}      static int droplets(const int i,const double t,Event *_ev){tracing("droplets","atomisation.c",177);
{
  scalar  m=new_scalar("m");
  foreach_stencil()
    {_stencil_val_a(m,0,0,0); _stencil_val(f,0,0,0);   }end_foreach_stencil();
  {
#line 180
foreach()
    val(m,0,0,0) = val(f,0,0,0) > 1e-3;end_foreach();}
  int n = tag (m);
#line 191 "atomisation.c"
  double v[n];
  coord b[n];
  for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = b[j].z = 0.;
  foreach_stencil ()
    {_stencil_val(m,0,0,0); {    
       _stencil_val(m,0,0,0);
_stencil_val(cm,0,0,0);_stencil_val(f,0,0,0);   
        
      
      
 
#line 201
{_stencil_val(cm,0,0,0);_stencil_val(f,0,0,0);  }
 
#line 201
{_stencil_val(cm,0,0,0);_stencil_val(f,0,0,0);  }
    }   }end_foreach_stencil();
  
#line 195
if(!is_constant(cm)){
#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
{
#line 195
foreach ()
    if (val(m,0,0,0) > 0) {
      int j = val(m,0,0,0) - 1;
      v[j] += (sq(Delta)*val(cm,0,0,0))*val(f,0,0,0);
      coord p = {x,y,z};
      
 b[j].x += (sq(Delta)*val(cm,0,0,0))*val(f,0,0,0)*p.x;
 
#line 201
b[j].y += (sq(Delta)*val(cm,0,0,0))*val(f,0,0,0)*p.y;
    }end_foreach();}
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif

#line 202
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
{
#line 195
foreach ()
    if (val(m,0,0,0) > 0) {
      int j = val(m,0,0,0) - 1;
      v[j] += (sq(Delta)*_const_cm)*val(f,0,0,0);
      coord p = {x,y,z};
      
 b[j].x += (sq(Delta)*_const_cm)*val(f,0,0,0)*p.x;
 
#line 201
b[j].y += (sq(Delta)*_const_cm)*val(f,0,0,0)*p.y;
    }end_foreach();}
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif

#line 202
}






  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);






  for (int j = 0; j < n; j++)
    fprintf (fout, "%d %g %d %g %g %g\n", i, t,
      j, v[j], b[j].x/v[j], b[j].y/v[j]);
  fflush (fout);delete((scalar*)((scalar[]){m,{-1}}));
}{end_tracing("droplets","atomisation.c",221);return 0;}end_tracing("droplets","atomisation.c",221);}







static int adapt_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++);*ip=i;*tp=t;return ret;}      static int adapt_0(const int i,const double t,Event *_ev){tracing("adapt_0","atomisation.c",229); {
  adapt_wavelet ((struct Adapt){((scalar[]){f,u.x,u.y,{-1}}), (double[]){0.01,uemax,uemax,uemax}, maxlevel});
}{end_tracing("adapt_0","atomisation.c",231);return 0;}end_tracing("adapt_0","atomisation.c",231);}
#line 2 "ast/init_solver.h"

static void _init_solver (void)
{
  void init_solver();
  datasize=13*sizeof(double);init_solver();
  quadtree_methods();{  event_register((Event){0,1,defaults,{defaults_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/run.h",42,"defaults"});
  event_register((Event){0,1,defaults_0,{defaults_0_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/navier-stokes/centered.h",126,"defaults"});
  event_register((Event){0,1,default_display,{default_display_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/navier-stokes/centered.h",179,"default_display"});
  event_register((Event){0,1,init,{init_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/navier-stokes/centered.h",188,"init"});
  event_register((Event){0,1,defaults_1,{defaults_1_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/vof.h",107,"defaults"});
  event_register((Event){0,1,defaults_2,{defaults_2_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/vof.h",127,"defaults"});
  event_register((Event){0,1,defaults_3,{defaults_3_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/two-phase.h",25,"defaults"});
  event_register((Event){0,1,defaults_4,{defaults_4_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/iforce.h",30,"defaults"});
  event_register((Event){0,1,init_0,{init_0_expr0},((int *)0),((double *)0),"atomisation.c",87,"init"});
  event_register((Event){0,1,logfile,{logfile_expr0},((int *)0),((double *)0),"atomisation.c",119,"logfile"});
  event_register((Event){0,1,movie,{movie_expr0},((int *)0),((double *)0),"atomisation.c",131,"movie"});
  event_register((Event){0,3,snapshot,{snapshot_expr0,snapshot_expr1,snapshot_expr2},((int *)0),((double *)0),"atomisation.c",158,"snapshot"});
  event_register((Event){0,1,droplets,{droplets_expr0},((int *)0),((double *)0),"atomisation.c",177,"droplets"});


    

    init_const_vector((vector){{_NVARMAX+0},{_NVARMAX+1}},"zerof",(double[]) {0.,0.,0.});
  init_const_vector((vector){{_NVARMAX+2},{_NVARMAX+3}},"unityf",(double[]) {1.,1.,1.});
  init_const_scalar((scalar){_NVARMAX+4},"unity", 1.);
  init_const_scalar((scalar){_NVARMAX+5},"zeroc", 0.);
  event_register((Event){0,1,cleanup,{cleanup_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/run.h",50,"cleanup"});
  init_scalar((scalar){0},"p");
  init_vector((vector){{1},{2}},"u");
  init_vector((vector){{3},{4}},"g");
  init_scalar((scalar){5},"pf");
  init_face_vector((vector){{6},{7}},"uf");
  event_register((Event){0,1,set_dtmax,{set_dtmax_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/navier-stokes/centered.h",214,"set_dtmax"});
  event_register((Event){0,1,stability,{stability_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/navier-stokes/centered.h",216,"stability"});
  event_register((Event){0,1,vof,{vof_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/navier-stokes/centered.h",226,"vof"});
  event_register((Event){0,1,tracer_advection,{tracer_advection_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/navier-stokes/centered.h",227,"tracer_advection"});
  event_register((Event){0,1,tracer_diffusion,{tracer_diffusion_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/navier-stokes/centered.h",228,"tracer_diffusion"});
  event_register((Event){0,1,properties,{properties_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/navier-stokes/centered.h",235,"properties"});
  event_register((Event){0,1,advection_term,{advection_term_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/navier-stokes/centered.h",308,"advection_term"});
  event_register((Event){0,1,viscous_term,{viscous_term_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/navier-stokes/centered.h",337,"viscous_term"});
  event_register((Event){0,1,acceleration,{acceleration_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/navier-stokes/centered.h",373,"acceleration"});
  event_register((Event){0,1,projection,{projection_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/navier-stokes/centered.h",413,"projection"});
  event_register((Event){0,1,end_timestep,{end_timestep_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/navier-stokes/centered.h",428,"end_timestep"});
  event_register((Event){0,1,adapt,{adapt_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/navier-stokes/centered.h",438,"adapt"});
  event_register((Event){0,1,stability_0,{stability_0_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/vof.h",140,"stability"});
  event_register((Event){0,1,vof_0,{vof_0_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/vof.h",364,"vof"});
  init_scalar((scalar){8},"f");
  init_face_vector((vector){{9},{10}},"alphav");
  init_scalar((scalar){11},"rhov");
  event_register((Event){0,1,tracer_advection_0,{tracer_advection_0_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/two-phase.h",64,"tracer_advection"});
  event_register((Event){0,1,properties_0,{properties_0_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/two-phase.h",95,"properties"});
  event_register((Event){0,1,acceleration_0,{acceleration_0_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/iforce.h",43,"acceleration"});
  event_register((Event){0,1,stability_1,{stability_1_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/tension.h",36,"stability"});
  event_register((Event){0,1,acceleration_1,{acceleration_1_expr0},((int *)0),((double *)0),"/home/gongminjiang/basilisk/src/tension.h",72,"acceleration"});
  init_scalar((scalar){12},"f0");
  event_register((Event){0,1,adapt_0,{adapt_0_expr0},((int *)0),((double *)0),"atomisation.c",229,"adapt"});

#line 13
}  _attribute[p.i].dirty=1,_attribute[p.i].boundary[right]=_boundary0,_attribute[p.i].boundary_homogeneous[right]=_boundary0_homogeneous;
  _attribute[p.i].dirty=1,_attribute[p.i].boundary[left]=_boundary1,_attribute[p.i].boundary_homogeneous[left]=_boundary1_homogeneous;
  _attribute[p.i].dirty=1,_attribute[p.i].boundary[top]=_boundary2,_attribute[p.i].boundary_homogeneous[top]=_boundary2_homogeneous;
  _attribute[p.i].dirty=1,_attribute[p.i].boundary[bottom]=_boundary3,_attribute[p.i].boundary_homogeneous[bottom]=_boundary3_homogeneous;
  _attribute[u.x.i].dirty=1,_attribute[u.x.i].boundary[left]=_boundary4,_attribute[u.x.i].boundary_homogeneous[left]=_boundary4_homogeneous;
  _attribute[u.y.i].dirty=1,_attribute[u.y.i].boundary[left]=_boundary5,_attribute[u.y.i].boundary_homogeneous[left]=_boundary5_homogeneous;
  _attribute[p.i].dirty=1,_attribute[p.i].boundary[left]=_boundary6,_attribute[p.i].boundary_homogeneous[left]=_boundary6_homogeneous;
  _attribute[f.i].dirty=1,_attribute[f.i].boundary[left]=_boundary7,_attribute[f.i].boundary_homogeneous[left]=_boundary7;
  _attribute[u.x.i].dirty=1,_attribute[u.x.i].boundary[right]=_boundary8,_attribute[u.x.i].boundary_homogeneous[right]=_boundary8_homogeneous;
  _attribute[p.i].dirty=1,_attribute[p.i].boundary[right]=_boundary9,_attribute[p.i].boundary_homogeneous[right]=_boundary9_homogeneous;

  
#line 14
set_fpe();
}
