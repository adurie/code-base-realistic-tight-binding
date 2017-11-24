#include <iostream>
#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include "TB.h"
/* #include "cunningham_spawn.h" */
#include "cunningham.h"


using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
typedef Matrix<dcomp, 9, 9> dmat;
typedef Matrix<dcomp, 18, 18> ddmat;
typedef Vector3d vec;

ddmat gs(ddmat &OM, ddmat &T)
{
	ddmat zero = ddmat::Zero();
	Matrix<dcomp, 36, 36> X,O;
	X << 	zero,	T.inverse(),
		-T.adjoint(),	OM*T.inverse();
	ComplexEigenSolver<Matrix<dcomp, 36, 36>> ces;
	ces.compute(X);
	O = ces.eigenvectors();
	ddmat b = O.topRightCorner(18, 18);
	ddmat d = O.bottomRightCorner(18, 18);
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
	Matrix<dcomp, 9, 9> u_11, u_12, u_21, T_21;
	u_11 = u + t_3*exp(i*d_3.dot(K))+ t_4*exp(i*d_4.dot(K))+ t_9*exp(i*d_9.dot(K)) + t_10*exp(i*d_10.dot(K)) + 
		t_13*exp(i*d_13.dot(K))+ t_14*exp(i*d_14.dot(K))+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	u_12 = t_1 + t_5*exp(i*d_9.dot(K)) + t_7*exp(i*d_14.dot(K)) + t_12*exp(i*d_4.dot(K));
	u_21 = t_2 + t_8*exp(i*d_13.dot(K)) + t_11*exp(i*d_3.dot(K)) + t_6*exp(i*d_10.dot(K));
	/* u_21 = u_12.adjoint(); */
	Matrix<complex<double>, 9, 9> zero = Matrix<complex<double>, 9, 9>::Zero();
	Matrix<dcomp, 18, 18> U, T, Tdagg, OM, GL, GR, GN, GRinv, GNinv;
	U << u_11, u_12, u_21, u_11;
	T_21 = t_7 + t_1*exp(i*d_13.dot(K)) + t_5*exp(i*d_3.dot(K)) + t_12*exp(i*d_10.dot(K));
	T << t_15, zero, T_21, t_15;

      	Matrix<complex<double>, 18, 18> I = Matrix<complex<double>, 18, 18>::Identity();
      	/* Matrix<complex<double>, 18, 18> I = Matrix<complex<double>, 18, 18>::Ones(); */
	/* cout<<U.trace()<<endl; */
	/* cout<<T.trace()<<endl; */
	/* cout<<U.determinant()<<endl; */
	/* cout<<T.determinant()<<endl; */
	/* SelfAdjointEigenSolver<Matrix<dcomp, 18, 18>> ces; */
	/* ComplexEigenSolver<Matrix<dcomp, 18, 18>> ces; */
	/* ces.compute(T); */
	/* cout<<ces.eigenvalues()<<endl; */

	OM = omega*I-U;
	Tdagg = T.adjoint();

	GL = gs(OM, T);
	GR = gs(OM, Tdagg);
	GRinv = GR.inverse();
	GNinv = GRinv - Tdagg*GL*T;
	GN = GNinv.inverse();
	/* cout<<GL.trace()<<endl; */
	/* cout<<GR.trace()<<endl; */

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
	
	double a = 1;
	
	//position vectors of nearest neighbours in fcc
	d_1 << a/2., a/2., 0;
	d_2 << -a/2., -a/2., 0;
	d_3 << a/2., 0, a/2.;
	d_4 << -a/2., 0, -a/2.;
	d_5 << 0, a/2., a/2.;
	d_6 << 0, -a/2., -a/2.;
	d_7 << -a/2., a/2., 0;
	d_8 << a/2., -a/2., 0;
	d_9 << -a/2., 0, a/2.;
	d_10 << a/2., 0, -a/2.;
	d_11 << 0, -a/2., a/2.;
	d_12 << 0, a/2., -a/2.;

	//position vectors of next nearest neighbours
	d_13 << a, 0, 0;
	d_14 << -a, 0, 0;
	d_15 << 0, a, 0;
	d_16 << 0, -a, 0;
	d_17 << 0, 0, a;
	d_18 << 0, 0, -a;

	//initialise onsite and hopping matrices for each nn
	Matrix<dcomp, 9, 9> u, E;
	Matrix<dcomp, 9, 9> t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9;
	Matrix<dcomp, 9, 9> t_10, t_11, t_12, t_13, t_14, t_15, t_16;
	Matrix<dcomp, 9, 9> t_17, t_18;
	u = TB(2, 0, 0, 8, d_1);
	t_1 = TB(2, 1, 0, 8, d_1);
	t_2 = TB(2, 1, 0, 8, d_2);
	t_3 = TB(2, 1, 0, 8, d_3);
	t_4 = TB(2, 1, 0, 8, d_4);
	t_5 = TB(2, 1, 0, 8, d_5);
	t_6 = TB(2, 1, 0, 8, d_6);
	t_7 = TB(2, 1, 0, 8, d_7);
	t_8 = TB(2, 1, 0, 8, d_8);
	t_9 = TB(2, 1, 0, 8, d_9);
	t_10 = TB(2, 1, 0, 8, d_10);
	t_11 = TB(2, 1, 0, 8, d_11);
	t_12 = TB(2, 1, 0, 8, d_12);

	t_13 = TB(2, 1, 1, 8, d_13);
	t_14 = TB(2, 1, 1, 8, d_14);
	t_15 = TB(2, 1, 1, 8, d_15);
	t_16 = TB(2, 1, 1, 8, d_16);
	t_17 = TB(2, 1, 1, 8, d_17);
	t_18 = TB(2, 1, 1, 8, d_18);

	dcomp i;
	i = -1.;
	i = sqrt(i);

	double result;

	double start = 0.35;
	double end = 0.5;
	double step = 0.0026;

	/* double start = 0; */
	/* double end = 1; */
	/* double step = 1/500.; */

	for (double j = start; j<end + step; j=j+step){

		/* result = greens(M_PI/2., M_PI/2.0, a, 0.57553 + 1e-4*i, u, t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9, */
		/* 		t_10, t_11, t_12, t_13, t_14, t_15, t_16, t_17, t_18, d_1, d_2, d_3, d_4, */
		/* 		d_5, d_6, d_7, d_8, d_9, d_10, d_11, d_12, d_13, d_14, d_15, d_16, d_17, d_18); */
		/* cout<<-0.5*result/M_PI<<endl; */

		/* result = greens(0, 0, a, j + 1e-4*i, u, t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9, */
		/* 		t_10, t_11, t_12, t_13, t_14, t_15, t_16, t_17, t_18, d_1, d_2, d_3, d_4, */
		/* 		d_5, d_6, d_7, d_8, d_9, d_10, d_11, d_12, d_13, d_14, d_15, d_16, d_17, d_18); */

		result = kspace(&greens, 3000, 0.001, 10, a, j + 1e-4*i, u, t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9,
				t_10, t_11, t_12, t_13, t_14, t_15, t_16, t_17, t_18, d_1, d_2, d_3, d_4,
				d_5, d_6, d_7, d_8, d_9, d_10, d_11, d_12, d_13, d_14, d_15, d_16, d_17, d_18);

		cout<<100*(j-start+step)/(end-start+step)<<"% completed"<<endl;

		Myfile<<j<<" "<<-result*a*a/(4.*M_PI*M_PI*M_PI)<<endl;
		/* Myfile<<j<<" "<<-0.5*result/M_PI<<endl; */
	}

	Myfile.close();
	return 0;
}
