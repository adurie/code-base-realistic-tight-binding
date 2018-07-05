#include <iostream>
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
using namespace std;

int main(){
	double a, b, c, d;
	a = 5;
	b = 3;
	c = -2;
	d = 47;
	SHFT(a,b,c,a + 4*(b-a));
	cout<<a<<endl;
	cout<<b<<endl;
	cout<<c<<endl;
	cout<<d<<endl;
	return 0;
}
