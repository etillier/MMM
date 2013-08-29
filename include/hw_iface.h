#ifndef _HW_IFACE_H_
#define _HW_IFACE_H_

#include <vector>
using namespace std;

bool hw_init();
bool hw_write(vector<unsigned int>& buf);
bool hw_read(vector<unsigned short>& buf);
void hw_close();

void hw_dbg_clear_counters();
void hw_dbg_get_counters(unsigned long long int* total_ticks, 
						 unsigned long long int* wu_ticks,
						 unsigned long long int* usage_ticks);

void hw_sim_cache_init();
void hw_sim_cache_clear();
void hw_sim_cache_access(unsigned int bit_addr);
double hw_sim_cache_get_hit_ratio();

#define HW_MAX_VERTS (1 << 12)
#define HW_MAX_PROBS (1 << 16)
#define HW_MAX_STACKSIZE (1 << 14)
#define HW_TICKS_PER_S 150000000

#endif
