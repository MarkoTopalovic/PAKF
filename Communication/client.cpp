// Bogdan
// novembar 2015
#include <string.h>
#include "SocketArrayTransfer.h"
        
#define BUF_SIZE    4096                    

#if defined (_WIN32) || defined (_WIN64)
	WSADATA wsa;
#endif
int server_socket;

extern "C"
{
	
#if defined (_WIN32) || defined (_WIN64)
void __cdecl CONNECT_TO_SERVIS(char*servername,int *SERVER_PORT)
#else 
void connect_to_servis_(char*servername,int *SERVER_PORT)
#endif
{
	int c;
	struct hostent *h;                      
	struct sockaddr_in channel;                   
	#if defined (_WIN32) || defined (_WIN64)
	  if (WSAStartup(MAKEWORD(2,2),&wsa) != 0)
		{
			printf("Failed. Error Code : %d",WSAGetLastError());
		}	
	#endif	
	 h = gethostbyname(servername);                   
	 if (!h) { printf("gethostbyname failed"); exit(0);}
	 server_socket = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP); 
	 if (server_socket <0) {printf("socket"); exit(0);} 
	 memset(&channel, 0, sizeof(channel)); 
	 channel.sin_family= AF_INET; 
	 memcpy(&channel.sin_addr.s_addr, h->h_addr, h->h_length); 
	 channel.sin_port= htons(*SERVER_PORT); 
	 c = connect(server_socket, (struct sockaddr *) &channel, sizeof(channel)); 
	 if (c < 0) {printf("connect failed"); exit(0);}

}


#if defined (_WIN32) || defined (_WIN64)
void __cdecl SEND_REQUEST(int *m,int *koliko,int *irows,int *jcols,double *vals,double *vektor) 
#else 
void send_request_(int *m,int *koliko,int *irows,int *jcols,double *vals,double *vektor) 
#endif
{  	  
	 sendArray(server_socket,m,1*sizeof(int),BUF_SIZE);
	 sendArray(server_socket,koliko,1*sizeof(int),BUF_SIZE);
	 sendArray(server_socket,irows,(*koliko)*sizeof(int),BUF_SIZE); 
	 sendArray(server_socket,jcols,(*koliko)*sizeof(int),BUF_SIZE);
	 sendArray(server_socket,vals,(*koliko)*sizeof(double),BUF_SIZE); 
	 sendArray(server_socket,vektor,(*m)*sizeof(double),BUF_SIZE);   

} 

#if defined (_WIN32) || defined (_WIN64)
void __cdecl RECV_SOLUTION(double *vektor,int*m) 
#else 
void recv_solution_(double *vektor,int*m) 
#endif
{
	   recvArray(server_socket,vektor,(*m)*sizeof(double),BUF_SIZE);    	
}

#if defined (_WIN32) || defined (_WIN64)
void __cdecl SEND_INT_ARRAY(int*m,int*ind) 
#else 
void send_int_array_(int*m,int*ind) 
#endif
{ 
	sendArray(server_socket,ind,(*m)*sizeof(int),BUF_SIZE);
}

#if defined (_WIN32) || defined (_WIN64)
void __cdecl SEND_DOUBLE_ARRAY(int*m,double*rhs) 
#else 
void send_double_array_(int*m,double*rhs) 
#endif
{ 
   sendArray(server_socket,rhs,(*m)*sizeof(double),BUF_SIZE);
}	


#if defined (_WIN32) || defined (_WIN64)
void __cdecl GETCHAR()
#else
void getchar_() 
#endif
{ getchar(); } 

}

