#include <iostream>
#include <exception>
#include <string>
using namespace std;

class CommException: public exception
{
	public:
	string msg;
	CommException(){};
	CommException(string msg){ this->msg = msg; }
	~CommException() throw(){};
	virtual const char* what() const throw(){return msg.c_str();}
};