#ifndef _SIGHAND_H_
#define _SIGHAND_H_

#include <string>
using namespace std;

void sighand_init(class CMMMAlgorithm* algo, string& filename);
void sighand_restore();

#endif
