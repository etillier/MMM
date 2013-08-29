#include <vector>
#include "tmports.h"
using namespace std;

static char g_portsfile[] =
"# Name	direction	bits	Handshake_from_circuit	Handshake_from_SUN \n"
"o_outdata	o	32	o_outdata_ready	i_outdata_want \n"
"i_indata	i	32	o_indata_want	i_indata_ready \n"
"o_dbg_outdata	o	64 \n"
"i_dbg_indata	i	32	o_dbg_indata_want	i_dbg_indata_have";

static char g_server[] = "skynet.eecg.toronto.edu";

static int g_p_indata;
static int g_p_outdata;
static int g_p_dbg_indata;
static int g_p_dbg_outdata;

bool hw_init()
{
	if (tm_init(g_server, g_portsfile))
		return false;

	g_p_outdata = tm_open("o_outdata", "r");
	g_p_indata = tm_open("i_indata", "w");

	g_p_dbg_outdata = tm_open("o_dbg_outdata", "r");
	g_p_dbg_indata = tm_open("i_dbg_indata", "w");

	if (g_p_outdata < 0 || g_p_indata < 0)
		return false;

	return true;
}


bool hw_write(vector<unsigned int>& buf)
{
	int bytes_to_write = (int)buf.size() * 4;

	int bytes_written = tm_write(g_p_indata, (char*)&buf[0], bytes_to_write);

	if (bytes_written != bytes_to_write)
		return false;

	return true;
}


bool hw_read(vector<unsigned short>& buf)
{
	buf.resize(2);

	if (tm_read(g_p_outdata, (char*)&buf[0], 4) != 4)
	{
		return false;
	}

	unsigned int n_cliques = buf[0];
	unsigned int global_maxsize = buf[1];
	int words_to_read = (int)((n_cliques * (global_maxsize+1) + 1) & 0xFFFFFFFE);

	if (words_to_read == 0)
		return true;
	
	buf.resize(2 + words_to_read);

	int bytes_to_read =  words_to_read * 2;
	int bytes_read = tm_read(g_p_outdata, (char*)&buf[2], bytes_to_read);
	if (bytes_read != bytes_to_read)
		return false;

	return true;
}


void hw_dbg_clear_counters()
{
	unsigned int cmd = 0x00000001;
	tm_write(g_p_dbg_indata, (char*)&cmd, sizeof(cmd));
}


void hw_dbg_get_counters(unsigned long long int* total_ticks, 
						 unsigned long long int* wu_ticks,
						 unsigned long long int* usage_ticks)
{
	unsigned int cmd = 0x00000000;
	tm_write(g_p_dbg_indata, (char*)&cmd, sizeof(cmd));
	tm_read(g_p_dbg_outdata, (char*)total_ticks, sizeof(unsigned long long int));

	cmd = 0x00010000;
	tm_write(g_p_dbg_indata, (char*)&cmd, sizeof(cmd));
	tm_read(g_p_dbg_outdata, (char*)wu_ticks, sizeof(unsigned long long int));

	cmd = 0x00020000;
	tm_write(g_p_dbg_indata, (char*)&cmd, sizeof(cmd));
	tm_read(g_p_dbg_outdata, (char*)usage_ticks, sizeof(unsigned long long int));
}


void hw_close()
{
	tm_close(g_p_outdata);
	tm_close(g_p_indata);
	tm_exit();
}


static const unsigned int CACHE_LINES = 32;
static const unsigned int CACHE_LINE_SIZE = 1024;

static struct
{
	bool valid;
	unsigned int tag;
} s_cache_lines[CACHE_LINES];

static unsigned int s_cache_hits = 0;
static unsigned int s_cache_misses = 0;

void hw_sim_cache_init()
{
	s_cache_hits = 0;
	s_cache_misses = 0;
}

void hw_sim_cache_clear()
{
	for (unsigned int i = 0; i < CACHE_LINES; i++)
	{
		s_cache_lines[i].valid = false;
	}
}

void hw_sim_cache_access(unsigned int bit_addr)
{
	unsigned int tag = bit_addr / (CACHE_LINES * CACHE_LINE_SIZE);
	unsigned int line = (bit_addr / CACHE_LINE_SIZE) & (CACHE_LINES-1);

	if (s_cache_lines[line].valid && s_cache_lines[line].tag == tag)
	{
		s_cache_hits++;
	}
	else
	{
		s_cache_misses++;
		s_cache_lines[line].valid = true;
		s_cache_lines[line].tag = tag;
	}
}

double hw_sim_cache_get_hit_ratio()
{
	if (s_cache_hits == 0 && s_cache_misses == 0)
		return 1.0;
	else
		return (double)s_cache_hits / (double)(s_cache_hits + s_cache_misses);
}