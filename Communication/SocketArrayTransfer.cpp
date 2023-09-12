#include "SocketArrayTransfer.h"
#ifndef MSG_NOSIGNAL
# define MSG_NOSIGNAL 0
#endif

void recvArray(int socket, void *array,int SIZE,int emit)
{
	int bytes = 0;
	for(int i=0;i<SIZE;)
	{
		if(emit > SIZE-i) emit = SIZE - i;
		bytes = recv(socket, &((char*)array)[i], emit, MSG_NOSIGNAL);
		if(bytes<=0) throw CommException(); 
        i  += bytes;	
	}
}

void sendArray(int socket, void *array,int SIZE,int emit)
{
	int bytes = 0;
	for(int i=0;i<SIZE;)
	{
		if(emit > SIZE-i) emit = SIZE - i;
		bytes = send(socket, &((char*)array)[i], emit, MSG_NOSIGNAL);
		if(bytes < 0) throw CommException(); 
        i += bytes;
	}
}
