#include <iostream>
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#include <eigen3/Eigen/Dense>
#include <vector>
using namespace std;
using namespace Eigen;

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
	VectorXcd T(5);
	T << 3, 5, -3.2, -2313, 0.17181;
	cout<<T.transpose()<<endl;
	T.reverseInPlace();
	cout<<T.transpose()<<endl;
	vector<int> vec;
	int ii;
	for (int i = 1; i < 10; i++){
		ii = pow(i,2);
		vec.emplace_back(ii);
	}
	for (int i = 0; i < vec.size(); i++)
		cout<<vec[i]<<endl;

	return 0;
}
