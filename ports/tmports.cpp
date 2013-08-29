/* A library of routines that will talk to a design on the TM.
 * The TM design must be wrapped with a communications layer like the
 * one that portmux_gen creates.
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "../include/profiler.h"

#ifdef _WIN32
	#include <winsock.h>
	#define close closesocket
	#define read(a,b,c) recv(a,b,c,0)
	#define write(a,b,c) send(a,b,c,0)
	#define fdopen _fdopen
	#define strdup _strdup
#else
	#include <unistd.h>
	#include <ctype.h>
	#include <sys/socket.h>
	#include <netinet/in.h>
	#include <netdb.h>
#endif



static char Rcsid[] = "$Header: /jayar/i/tm4/src/ports/RCS/tmports.c,v 1.12 2010/07/15 17:52:18 drg Exp $";

#define DEFAULT_TMHOST	"andy.eecg"
#define NUM_TM4_CHIPS	4
#define TM_MAX_PACKET_DATA	16383
#define TM_MAX_PORTS		256
#define TM_PACKET_HLEN		4
#define TM_PORTA	3462
#define TM_PORTB	9100

struct tm_packet
{
	unsigned int h;	
	char data[TM_MAX_PACKET_DATA];
};

#define TM_PACKET_WRITE(p)	(((p).h & 0x80000000) != 0)
#define TM_PACKET_NEEDS_ACK(p)	(((p).h & 0x40000000) != 0)
#define TM_PACKET_CHIP(p)	(((p).h & 0x30000000) >> 28)
#define TM_PACKET_PORT(p)	(((p).h & 0x03FC0000) >> 18)
#define TM_PACKET_PORT_BYTES(p)	((((p).h & 0x0000F000) >> 10) | (((p).h & 0x00030000) >> 16))
#define TM_PACKET_LENGTH(p)	(((p).h & 0x00000FFF) | (((p).h & 0x0C000000) >> 14))

#define TM_SET_PACKET_WRITE(p, val)	(p).h = ((p).h & ~0x80000000) | ((val)<<31);
#define TM_SET_PACKET_NEEDS_ACK(p, val)	(p).h = ((p).h & ~0x40000000) | (((val)<<30) & 0x40000000);
#define TM_SET_PACKET_CHIP(p, val)	(p).h = ((p).h & ~0x30000000) | (((val)<<28) & 0x30000000);
#define TM_SET_PACKET_PORT(p, val)	(p).h = ((p).h & ~0x03FC0000) | (((val)<<18) & 0x03FC0000);
#define TM_SET_PACKET_PORT_BYTES(p,val) (p).h = ((p).h & ~0x0003F000) | (((val)<<10) & 0x0000F000) | (((val)<<16) & 0x00030000);;
#define TM_SET_PACKET_LENGTH(p, val)	(p).h = ((p).h & ~0x0C000FFF) | ((val) & 0x00000FFF) | (((val)<<14) & 0x0C000000);

static struct portinfo
{
	char *name;
	char iox;
	int nbits;
	int needs_ack;
	int chip;
	int port;
} portinfo[TM_MAX_PORTS];

static struct tm_packet	packet;
static int tm_sock;
static int nports;

static int l_trychip(char* portdefs);
static int l_socket_read(int fd, char *buf, int nbytes);
static int l_socket_write(int fd, char *buf, int nbytes);

/* Contact the tmumon process, and get the list of ports */

int tm_init(const char *tmhost, char* portdefs)
{
	int i;
	struct sockaddr_in serveraddr;
	struct hostent *hp;

	#ifdef _WIN32
	WSADATA wsadata;
	if (WSAStartup(MAKEWORD(2,2), &wsadata) != 0)
	{
		perror("tm_init: can't initialize WinSock");
		return(-1);
	}
	#endif

	if(tm_sock > 0)
	{
		close(tm_sock);
		for(i = 0; i<nports; i++)
		{
			free(portinfo[nports].name);
		}
		nports = 0;
	}

	tm_sock = socket(AF_INET, SOCK_STREAM, 0);
	if(tm_sock < 0)
	{
		perror("tm_init: opening stream socket");
		return(-1);
	}

	if((tmhost == (char *) NULL) || (tmhost[0] == '\0'))
	{
		if((tmhost = getenv("TM_SERVER")) == (char *) NULL)
			tmhost = DEFAULT_TMHOST;
	}

	serveraddr.sin_family = AF_INET;
	hp = gethostbyname(tmhost);
	if(hp == (struct hostent *) NULL)
	{
		fprintf(stderr, "tm_init: can't find address for '%s'\n", tmhost);
		return(-1);
	}

	memcpy(&serveraddr.sin_addr.s_addr, hp->h_addr, sizeof(serveraddr.sin_addr.s_addr));
	serveraddr.sin_port = htons(TM_PORTB);

	if(connect(tm_sock, (struct sockaddr *) &serveraddr, sizeof(serveraddr)))
	{
		perror("tm_init: can't connect to tmumon");
		return(-1);
	}

	return l_trychip(portdefs);
}


void tm_exit()
{
	int i;
	for (i = 0; i < nports; i++)
	{
		free(portinfo[i].name);
	}

	close(tm_sock);
}


static int l_trychip(char portdefs[])
{
	char* buf;
	char name[1024], mode[1024], handshake_from_circuit[1024], handshake_from_sun[1024]; 
	int nbits, nfields, needs_ack;
	int portno;


	portno = 0;
	buf = strtok(portdefs, "\n");
	while(buf != NULL)
	{
		if((buf[0] == '#') || (buf[0] == '\0') || (buf[0] == '\n'))
		{
			buf = strtok(NULL, "\n");
			continue;
		}

		nfields = sscanf(buf, "%s %s %d %s %s", name, mode, &nbits,
						handshake_from_circuit, handshake_from_sun);
		if(nfields == 3)
		{
			needs_ack = 0;
		}
		else if(nfields == 5)
		{
			needs_ack = 1;
		}
		else
		{
			fprintf(stderr, "tm_init: don't understand line '%s'", buf);
			return -1;
		}

		portinfo[nports].name = strdup(name);
		if(isupper(mode[0]))
			mode[0] = tolower(mode[0]);
		portinfo[nports].iox = mode[0];
		portinfo[nports].nbits = nbits;
		portinfo[nports].needs_ack = needs_ack;
		portinfo[nports].chip = 0;
		portinfo[nports].port = portno++;
		nports++;

		buf = strtok(NULL, "\n");
	}

	return 0;
}



int tm_open(const char *portname, const char *mode)
{
	int portno;

	if(tm_sock <= 0)
		return(-1);
	
	for(portno = 0; portno < nports; portno++)
	{
		if(!strcmp(portname, portinfo[portno].name))
		{
			switch(portinfo[portno].iox)
			{
			case 'i':
				if(mode[0] != 'w')
					return(-1);
				break;

			case 'o':
				if(mode[0] != 'r')
					return(-1);
				break;
			}
			
			return(portno);
		}
	}

	return(-1);
}


int tm_get_port_width(char *portname)
{
	int portno;

	if(tm_sock <= 0)
		return(-1);
	
	for(portno = 0; portno < nports; portno++)
	{
		if(!strcmp(portname, portinfo[portno].name))
		{
			return(portinfo[portno].nbits);
		}
	}

	return(-1);
}


int tm_write(int port, char *buf, int nbytes)
{
	int n, retval;
	int portbytes, i;
	char *packet_data_ptr, *itemptr;
	int save_nbytes = nbytes;
	int result_code;

	if(tm_sock <= 0)
		return(-1);

	if((port < 0) || (port >= nports))
		return(-1);

	if(portinfo[port].iox != 'i')
		return(-1);

	portbytes = (portinfo[port].nbits + 7) / 8;
	if((nbytes % portbytes) != 0)
		return(-1);
	
	while(nbytes > 0)
	{
		if(nbytes > TM_MAX_PACKET_DATA)
			n = (TM_MAX_PACKET_DATA / portbytes) * portbytes;
		else
			n = nbytes;

		packet.h = 0;
		TM_SET_PACKET_WRITE(packet, 1);
		TM_SET_PACKET_NEEDS_ACK(packet, portinfo[port].needs_ack);
		TM_SET_PACKET_CHIP(packet, portinfo[port].chip);
		TM_SET_PACKET_PORT(packet, portinfo[port].port);
		TM_SET_PACKET_PORT_BYTES(packet, (portinfo[port].nbits + 7) / 8);
		TM_SET_PACKET_LENGTH(packet, n);

		/* Convert the data to least significant byte first.  */

		if(ntohl(1) == 1)
		{
			packet_data_ptr = packet.data;
			for(itemptr = buf; itemptr < &buf[n]; itemptr += portbytes)
			{
				for(i=0; i<portbytes; i++)
				{
					packet_data_ptr[i] = itemptr[portbytes - 1 - i];
				}

				packet_data_ptr += portbytes;
			}
		}
		else
		{
			memcpy(packet.data, buf, n);
		}
	
		packet.h = htonl(packet.h);
		retval = l_socket_write(tm_sock, (char *) &packet, n + TM_PACKET_HLEN);
		if(retval != (n + TM_PACKET_HLEN))
			return(-1);

		retval = l_socket_read(tm_sock, (char *) &result_code, sizeof(result_code));
		if(retval != sizeof(result_code))
			return(-1);

		result_code = ntohl(result_code);
		if(result_code == -1)
			return(-1);

		nbytes -= result_code;
		buf += result_code;
	}

	return(save_nbytes - nbytes);
}


int tm_read(int port, char *buf, int nbytes)
{
	int n, retval;
	int portbytes, i;
	char *itemptr;
	int save_nbytes = nbytes;

	if(tm_sock <= 0)
		return(-1);

	if((port < 0) || (port >= nports))
		return(-1);

	if(portinfo[port].iox != 'o')
		return(-1);

	portbytes = (portinfo[port].nbits + 7) / 8;
	if((nbytes % portbytes) != 0)
		return(-1);
	
	while(nbytes > 0)
	{
		if(nbytes > TM_MAX_PACKET_DATA)
			n = (TM_MAX_PACKET_DATA / portbytes) * portbytes;
		else
			n = nbytes;

		packet.h = 0;
		TM_SET_PACKET_WRITE(packet, 0);
		TM_SET_PACKET_NEEDS_ACK(packet, portinfo[port].needs_ack);
		TM_SET_PACKET_CHIP(packet, portinfo[port].chip);
		TM_SET_PACKET_PORT(packet, portinfo[port].port);
		TM_SET_PACKET_PORT_BYTES(packet, (portinfo[port].nbits + 7) / 8);
		TM_SET_PACKET_LENGTH(packet, n);

		packet.h = htonl(packet.h);
		retval = l_socket_write(tm_sock, (char *) &packet, TM_PACKET_HLEN);
		if(retval != TM_PACKET_HLEN)
			return(-1);

		retval = l_socket_read(tm_sock, (char *) &packet, TM_PACKET_HLEN);
		if(retval != TM_PACKET_HLEN)
			return(-1);

		packet.h = ntohl(packet.h);
		if(packet.h == -1)
			return(-1);

		retval = l_socket_read(tm_sock, (char *) packet.data, n);
		if(retval != n)
			return(-1);

		/* Convert the data we just read to host byte order.  */

		if(ntohl(1) == 1)
		{
			for(itemptr = buf; itemptr < &buf[n]; itemptr += portbytes)
			{
				for(i=0; i<portbytes; i++)
				{
					itemptr[i] = packet.data[(itemptr - buf) + portbytes - 1 - i];
				}
			}
		}
		else
		{
			memcpy(buf, packet.data, n);
		}

		nbytes -= n;
		buf += n;
	}

	return(save_nbytes - nbytes);
}


int tm_close(int port)
{
	if(tm_sock <= 0)
		return(-1);

	if((port < 0) || (port >= nports))
		return(-1);

	return(0);
}


static int l_socket_read(int fd, char *buf, int nbytes)
{
   /* Read nbytes from a socket.  Handle partial reads. */

   int n, retval;
   char *cp;

   cp = buf;
   n = nbytes;

   while(n > 0)
   {
      retval = read(fd, cp, n);

      if(retval <= 0)
	  {
         return(-1);
      }

      n -= retval;
      cp += retval;
   }

   return(nbytes);
}


static int l_socket_write(int fd, char *buf, int nbytes)
{
   /* Write nbytes to a socket.  Handle partial writes. */

   int n, retval;
   char *cp;

   cp = buf;
   n = nbytes;

   while(n > 0)
   {
      retval = write(fd, cp, n);

      if(retval <= 0)
	  {
         return(-1);
      }

      n -= retval;
      cp += retval;
   }

   return(nbytes);
}

