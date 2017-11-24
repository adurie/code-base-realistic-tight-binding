#include <iostream>
#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "cunningham_spawn.h"
#include "TBdynamic.h"

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
typedef Matrix<dcomp, 4, 4> dmat;
typedef Matrix<dcomp, 8, 8> ddmat;
typedef Matrix<dcomp, 16, 16> dddmat;
typedef Vector3d vec;

ddmat gs(ddmat &OM, ddmat &T)
{
	ddmat zero = ddmat::Zero();
	dddmat X,O;
	X << 	zero,	T.inverse(),
		-T.adjoint(),	OM*T.inverse();
	ComplexEigenSolver<dddmat> ces;
	ces.compute(X);
	O = ces.eigenvectors();
	ddmat b = O.topRightCorner(8, 8);
	ddmat d = O.bottomRightCorner(8, 8);
	ddmat GR;
	GR = b*d.inverse();
	return GR;
}

double greens(double k_x, double k_z, double a, dcomp omega, dmat &u, dmat &t_1,
		dmat &t_2, dmat &t_3, dmat &t_4, dmat &t_5, dmat &t_6, dmat &t_7, 
		dmat &t_8, dmat &t_9, dmat &t_10, dmat &t_11, dmat &t_12, dmat &t_13,
	  	dmat &t_14, dmat &t_15, dmat &t_16, dmat &t_17, dmat &t_18, vec &d_1,
	       	vec &d_2, vec &d_3, vec &d_4, vec &d_5, vec &d_6, vec &d_7, vec &d_8,
	       	vec &d_9, vec &d_10, vec &d_11, vec &d_12, vec &d_13, vec &d_14, vec &d_15,
	       	vec &d_16, vec &d_17, vec &d_18){

	dcomp i;
	i = -1.;
	i = sqrt(i);
	double k_y = 0;

	Vector3d K;
	K(0) = k_x;
	K(1) = k_y;
	K(2) = k_z;

	//construct diagonalised in-plane matrices
	dmat u_11, u_12, u_21, T_21;
	u_11 = u + t_3*exp(-i*d_3.dot(K))+ t_4*exp(-i*d_4.dot(K))+ t_9*exp(-i*d_9.dot(K)) + t_10*exp(-i*d_10.dot(K)) + 
		t_13*exp(-i*d_13.dot(K))+ t_14*exp(-i*d_14.dot(K))+ t_17*exp(-i*d_17.dot(K)) + t_18*exp(-i*d_18.dot(K));
	u_12 = t_1 + t_5*exp(-i*d_9.dot(K)) + t_7*exp(-i*d_14.dot(K)) + t_12*exp(-i*d_4.dot(K));
	u_21 = t_2 + t_8*exp(-i*d_13.dot(K)) + t_11*exp(-i*d_3.dot(K)) + t_6*exp(-i*d_10.dot(K));
	ddmat U, T, Tdagg, OM, GL, GR, GN, GRinv, GNinv;
	U << u_11, u_12, u_21, u_11;
	dmat zero = dmat::Zero();
	T_21 = t_7 + t_1*exp(-i*d_13.dot(K)) + t_5*exp(-i*d_3.dot(K)) + t_12*exp(-i*d_10.dot(K));
	T << t_15, zero, T_21, t_15;

      	ddmat I = ddmat::Identity();

	OM = omega*I-U;
	Tdagg = T.adjoint();

	GL = gs(OM, T);
	GR = gs(OM, Tdagg);
	GRinv = GR.inverse();
	GNinv = GRinv - T.adjoint()*GL*T;
	GN = GNinv.inverse();

	return imag(GN.trace());

}

int main(){

	cout<<"Name the data file\n";
	string Mydata;
	getline(cin, Mydata);
	ofstream Myfile;	
	Mydata += ".txt";
	Myfile.open( Mydata.c_str(),ios::trunc );

	Vector3d d_1, d_2, d_3, d_4, d_5, d_6, d_7, d_8, d_9;
	Vector3d d_10, d_11, d_12, d_13, d_14, d_15, d_16;
	Vector3d d_17, d_18;
	
	double a = 1.;

	//position vectors of nearest neighbours in fcc
	d_1 << a, a, 0;
	d_2 << -a, -a, 0;
	d_3 << a, 0, a;
	d_4 << -a, 0, -a;
	d_5 << 0, a, a;
	d_6 << 0, -a, -a;
	d_7 << -a, a, 0;
	d_8 << a, -a, 0;
	d_9 << -a, 0, a;
	d_10 << a, 0, -a;
	d_11 << 0, -a, a;
	d_12 << 0, a, -a;

	//position vectors of next nearest neighbours
	d_13 << 2*a, 0, 0;
	d_14 << -2*a, 0, 0;
	d_15 << 0, 2*a, 0;
	d_16 << 0, -2*a, 0;
	d_17 << 0, 0, 2*a;
	d_18 << 0, 0, -2*a;

	//initialise onsite and hopping matrices for each nn
	dmat u, E;
	dmat t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9;
	dmat t_10, t_11, t_12, t_13, t_14, t_15, t_16;
	dmat t_17, t_18;
	int dim = 4;
	u = TB(2, 0, 0, dim, d_1);
	t_1 = TB(2, 1, 0, dim, d_1);
	t_2 = TB(2, 1, 0, dim, d_2);
	t_3 = TB(2, 1, 0, dim, d_3);
	t_4 = TB(2, 1, 0, dim, d_4);
	t_5 = TB(2, 1, 0, dim, d_5);
	t_6 = TB(2, 1, 0, dim, d_6);
	t_7 = TB(2, 1, 0, dim, d_7);
	t_8 = TB(2, 1, 0, dim, d_8);
	t_9 = TB(2, 1, 0, dim, d_9);
	t_10 = TB(2, 1, 0, dim, d_10);
	t_11 = TB(2, 1, 0, dim, d_11);
	t_12 = TB(2, 1, 0, dim, d_12);

	t_13 = TB(2, 1, 1, dim, d_13);
	t_14 = TB(2, 1, 1, dim, d_14);
	t_15 = TB(2, 1, 1, dim, d_15);
	t_16 = TB(2, 1, 1, dim, d_16);
	t_17 = TB(2, 1, 1, dim, d_17);
	t_18 = TB(2, 1, 1, dim, d_18);

	
	dcomp i;
	i = -1.;
	i = sqrt(i);

	double start = -0.2;
	double end = 1.1;
	/* double end = start; */
	double step = 0.001;
	double result;

	for (double j = start; j<end + step; j=j+step){

		/* result = greens(0, 0, 2*a, j + 1e-4*i, u, t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9, */
		/* 		t_10, t_11, t_12, t_13, t_14, t_15, t_16, t_17, t_18, d_1, d_2, d_3, d_4, */
		/* 		d_5, d_6, d_7, d_8, d_9, d_10, d_11, d_12, d_13, d_14, d_15, d_16, d_17, d_18); */

		result = kspace(&greens, 3, 0.01, 2*a, j + 1e-4*i, u, t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9,
				t_10, t_11, t_12, t_13, t_14, t_15, t_16, t_17, t_18, d_1, d_2, d_3, d_4,
				d_5, d_6, d_7, d_8, d_9, d_10, d_11, d_12, d_13, d_14, d_15, d_16, d_17, d_18);

		cout<<100*(j-start+step)/(end-start+step)<<"% completed"<<endl;

		Myfile<<j<<" "<<-result*a*a/(M_PI*M_PI*M_PI)<<endl;
		/* Myfile<<j<<" "<<-0.5*result/M_PI<<endl; */
	}

	Myfile.close();
	return 0;
}
