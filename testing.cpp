#include <cmath>
#include <iostream>
#include "cunningham_diamond.h"
#include <fstream>
#include <eigen3/Eigen/Dense>

double test(double x, double y, double a){
	return 1;
}

int main(){
	double a = 1;
	cout<<kspace(&test, 0, 5e-2, 0, a)/(4*M_PI*M_PI);
	return 0;
}
