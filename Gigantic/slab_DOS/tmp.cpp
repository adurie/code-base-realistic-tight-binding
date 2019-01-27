#include <eigen3/Eigen/Dense>
#include <iostream>
#include <cmath>
using namespace std;
using namespace Eigen;

int main(){
	complex<double> im = -1;
	im = sqrt(im);
	Matrix3cd tmp;
	Matrix3cd I;
	I.setIdentity();
	complex<double> el1, el2, el3, el4;
	el1 = 3.+4.*im;
	el2 = 3.-4.*im;
	el3 = 8.-14.*im;
	el4 = 8.+14.*im;
	tmp << 1, el1, 7,
	       el2, 7, el3,
	       7, el4, 21; 
	I = I*el2;
	SelfAdjointEigenSolver<Matrix3cd> ES;
	ES.compute(tmp);
	Matrix3cd vecs;
	Matrix3d lam;
	lam = ES.eigenvalues().asDiagonal();
	vecs = ES.eigenvectors();
	cout<< (tmp+I).inverse() <<endl<<endl;
	Matrix3cd QVQ;
	QVQ = vecs*(lam.cast<complex<double>>()+I).inverse()*vecs.inverse();
	cout<< QVQ <<endl;
	cout<<QVQ.array() + 3<<endl;
	return 0;
}

