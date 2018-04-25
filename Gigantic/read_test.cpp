#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <eigen3/Eigen/Dense>
using namespace std;
using namespace Eigen;

int main(){
	ifstream infile("iron_up_hr.dat");
	string line;
	double a, b, c, d, e, f, g;
	double eV_Ry = 0.073498618;
	complex<double> im;
	im = -1.;
	im = sqrt(im);
	MatrixXcd U(9,9), H(9,9);
	U.fill(0.);
	while (!infile.eof()) 
	{
		getline(infile, line);
		istringstream iss(line);
		iss >> a >> b >> c >> d >> e >> f >> g;
		if ((a == 0) && (b == 0) && (c == 0))
		{
			if (d == e)
				U(d-1,e-1) = (f + g*im)*eV_Ry;
		}
		if ((a == 1) && (b == 1) && (c == 1))
			H(d-1,e-1) = (f + g*im)*eV_Ry;
	}
	cout<<U.real()<<endl<<endl;
	cout<<H.real()<<endl;
	return 0;
}
