#ifndef _TIMER_H_
#define _TIMER_H_

#ifdef _WIN32

typedef unsigned long Timeval;

#else

#include "sys/time.h"

typedef struct timeval Timeval;

#endif

double timer_diff(Timeval* time1, Timeval* time2);
void timer_gettime(Timeval* out);

#endif

