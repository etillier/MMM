#ifndef _TMPORTS_H_
#define _TMPORTS_H_

int tm_init(const char *tmhost, char* portdefs);
void tm_exit();
int tm_open(const char* portname, const char* mode);
int tm_close(int port);
int tm_read(int port, char* buf, int nbytes);
int tm_write(int port, char* buf, int nbytes);
int tm_get_port_width(char *portname);

#endif
