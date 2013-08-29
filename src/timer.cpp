#include "timer.h"

#ifdef _WIN32

#define WIN32_LEAN_AND_MEAN
#include "windows.h"

double timer_diff(DWORD* time1, DWORD* time2)
{
	return (double)(*time2 - *time1) / 1000.0f;
}

void timer_gettime(DWORD* out)
{
	*out = GetTickCount();
}

#else

double timer_diff(struct timeval *time1, struct timeval *time2)
{
	double tmp_time, tmp_usec;

	tmp_time = time2->tv_sec - time1->tv_sec;
	tmp_usec = time2->tv_usec - time1->tv_usec;

	if (tmp_usec < 0.0)
	{
		tmp_usec += 1000000.0;
		if (tmp_time > 0.0)
			tmp_time -= 1.0;
	}
	return tmp_time + (tmp_usec / 1000000.0);
}

void timer_gettime(struct timeval* out)
{
	gettimeofday(out, 0);
}

#endif
