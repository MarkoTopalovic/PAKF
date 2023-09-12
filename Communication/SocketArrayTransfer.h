#if defined (_WIN32) || defined (_WIN64)

#include<winsock2.h>
#pragma comment(lib,"ws2_32.lib") 

#else
#include <sys/types.h>                  		
#include <sys/fcntl.h> 
#include <sys/socket.h> 
#include <netinet/in.h> 
#include <netdb.h> 
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CommException.h"


void sendArray(int socket, void *array,int SIZE,int emit);
void recvArray(int socket, void *array,int SIZE,int emit);


