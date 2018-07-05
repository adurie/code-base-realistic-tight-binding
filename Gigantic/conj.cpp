#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "cgls.h"

using namespace std;
using namespace Eigen;

double xargs(VectorXd& p){
	double f;
	f = pow((p(0)-1),4) + 2*pow((p(1)-1),4) + 3*pow((p(2)-3),4);
	return f;
}

void dxargs(VectorXd& p, VectorXd& xi){
	xi(0) = 4*pow((p(0)-1),3);
	xi(1) = 8*pow((p(1)-1),3);
	xi(2) = 12*pow((p(2)-3),3);
}

int main(){
	VectorXd p(3);
	p << 1.1, 2.2, 3.3;
	double ftol = 1e-5;
	double fret;
	int iter;
	frprmn(p, ftol, iter, fret, xargs, dxargs);
	cout<<p.transpose()<<endl;
	cout<<"number of iterations: "<<iter<<endl;
	cout<<"minimum of function: "<<fret<<endl;
	return 0;
}
