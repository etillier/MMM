#include <string>
#include "sighand.h"
#include "mmm_algorithm.h"

using namespace std;

static CMMMAlgorithm* s_algo_ptr = 0;
static string s_state_file;

#ifdef _WIN32

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

static BOOL handler(DWORD signum)
{
	if (signum == CTRL_BREAK_EVENT)
	{
		s_algo_ptr->stateSave(s_state_file);
		return TRUE;
	}
	else if (signum == CTRL_C_EVENT)
	{
		exit(0);
	}

	return FALSE;
}


void sighand_init(CMMMAlgorithm* algo, string& filename)
{
	s_state_file = filename;
	s_algo_ptr = algo;
	SetConsoleCtrlHandler( (PHANDLER_ROUTINE) handler, TRUE );
}


void sighand_restore()
{
	SetConsoleCtrlHandler( NULL, TRUE );
}


#else

#include <signal.h>

static void handler(int signum)
{
	if (signum == SIGUSR2)
	{
		s_algo_ptr->stateSave(s_state_file);
	}
}


void sighand_init(CMMMAlgorithm* algo, string& filename)
{
	s_state_file = filename;
	s_algo_ptr = algo;
	signal(SIGUSR2, handler);
}


void sighand_restore()
{
	signal(SIGUSR2, SIG_DFL);	
}

#endif
