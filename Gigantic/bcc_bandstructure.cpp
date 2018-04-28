#include <iostream>
#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "TB.h"
/* #include "TBdynamic.h" */

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
//calculates the bandstructure of fcc Cu

int main(){

	cout<<"Name the data file\n";
	string Mydata;
	getline(cin, Mydata);
	ofstream Myfile;	
	Mydata += ".txt";
	Myfile.open( Mydata.c_str(),ios::trunc );
	int ind = 1;

	Vector3d d_1, d_2, d_3, d_4, d_5, d_6, d_7, d_8;
	Vector3d d_13, d_14, d_15, d_16;
	Vector3d d_17, d_18;
	/* double a = 6.692; */
	double a = 1.;
	
	//position vectors of nearest neighbours in fcc
	d_1 << a, a, a;
	d_2 << -a, -a, -a;
	d_3 << -a, a, a;
	d_4 << a, -a, -a;
	d_5 << a, -a, a;
	d_6 << -a, a, -a;
	d_7 << a, a, -a;
	d_8 << -a, -a, a;

	//position vectors of next nearest neighbours
	d_13 << 2*a, 0, 0;
	d_14 << -2*a, 0, 0;
	d_15 << 0, 2*a, 0;
	d_16 << 0, -2*a, 0;
	d_17 << 0, 0, 2*a;
	d_18 << 0, 0, -2*a;

	//initialise onsite and hopping matrices for each nn
	Matrix<dcomp, 9, 9> u;
	Matrix<dcomp, 9, 9> t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8;
	Matrix<dcomp, 9, 9> t_13, t_14, t_15, t_16;
	Matrix<dcomp, 9, 9> t_17, t_18;
	u = TB(ind, 0, 0, 8, d_1);
	t_1 = TB(ind, 1, 0, 8, d_1);
	t_2 = TB(ind, 1, 0, 8, d_2);
	t_3 = TB(ind, 1, 0, 8, d_3);
	t_4 = TB(ind, 1, 0, 8, d_4);
	t_5 = TB(ind, 1, 0, 8, d_5);
	t_6 = TB(ind, 1, 0, 8, d_6);
	t_7 = TB(ind, 1, 0, 8, d_7);
	t_8 = TB(ind, 1, 0, 8, d_8);

	t_13 = TB(ind, 1, 1, 8, d_13);
	t_14 = TB(ind, 1, 1, 8, d_14);
	t_15 = TB(ind, 1, 1, 8, d_15);
	t_16 = TB(ind, 1, 1, 8, d_16);
	t_17 = TB(ind, 1, 1, 8, d_17);
	t_18 = TB(ind, 1, 1, 8, d_18);
	Myfile<<"P X Y"<<endl;

	dcomp i;
	i = -1.;
	i = sqrt(i);

	Matrix<dcomp, 9, 9> E;

	/* Matrix<dcomp, 18, 18> E; */
	/* Matrix<dcomp, 9, 9> u_11, u_12; */
	/* Matrix<dcomp, 18, 18> U, T; */
	/* Matrix<complex<double>, 9, 9> zero = Matrix<complex<double>, 9, 9>::Zero(); */

	double k_x, k_y, k_z, pi;
	Vector3d K;
	for (int k = 0; k < 251; k++)
	{
		if (k < 101){
			pi = M_PI*k/100.;
			k_x = pi;
			k_y = 0;
			k_z = 0;
		}
		if ((k > 100) && (k < 151)){
			pi = M_PI*(k-100)/100.;
			k_x = M_PI;
			k_y = pi;
			k_z = 0;
		}	
		if ((k > 150) && (k < 201)){
			pi = M_PI*(k-150)/100.;
			k_x = M_PI - pi;
			k_y = M_PI/2.;
			k_z = pi;
		}
		if ((k > 200) && (k < 251)){
			pi = M_PI*(k-200)/100.;
			k_x = M_PI/2.-pi;
			k_y = M_PI/2.-pi;
			k_z = M_PI/2.-pi;
		}
		K(0) = k_x;
		K(1) = k_y;
		K(2) = k_z;

		//fully diagonalised Hamiltonian
		E = u + t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
			+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
				+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K))
				+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
				+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
				+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
		SelfAdjointEigenSolver<Matrix<dcomp, 9, 9>> es;
		es.compute(E);
		Matrix<double, 9, 1> O;
		O = es.eigenvalues();

		/* //Hamiltonian firstly diagonalised in-plane */
		/* u_11 = u + t_3*exp(i*d_3.dot(K))+ t_4*exp(i*d_4.dot(K))+ t_9*exp(i*d_9.dot(K)) + t_10*exp(i*d_10.dot(K)) + */ 
		/* 	t_13*exp(i*d_13.dot(K))+ t_14*exp(i*d_14.dot(K))+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K)); */
		/* u_12 = t_1 + t_5*exp(i*d_9.dot(K)) + t_7*exp(i*d_14.dot(K)) + t_12*exp(i*d_4.dot(K)); */
		/* U << u_11, u_12, u_12.adjoint(), u_11; */
		/* T << t_15, zero, u_12.adjoint(), t_16; */
		/* E = U + T*exp(i*d_15.dot(K)) + T.adjoint()*exp(i*d_16.dot(K)); */
		/* SelfAdjointEigenSolver<Matrix<dcomp, 18, 18>> es; */
		/* es.compute(E); */
		/* Matrix<double, 18, 1> O; */
		/* O = es.eigenvalues(); */

		Myfile<<"A"<<" "<<k<<" "<<O(0)<<endl;
		Myfile<<"B"<<" "<<k<<" "<<O(1)<<endl;
		Myfile<<"C"<<" "<<k<<" "<<O(2)<<endl;
		Myfile<<"D"<<" "<<k<<" "<<O(3)<<endl;
		Myfile<<"E"<<" "<<k<<" "<<O(4)<<endl;
		Myfile<<"F"<<" "<<k<<" "<<O(5)<<endl;
		Myfile<<"G"<<" "<<k<<" "<<O(6)<<endl;
		Myfile<<"H"<<" "<<k<<" "<<O(7)<<endl;
		Myfile<<"I"<<" "<<k<<" "<<O(8)<<endl;

		/* Myfile<<"J"<<" "<<k<<" "<<O(9)<<endl; */
		/* Myfile<<"K"<<" "<<k<<" "<<O(10)<<endl; */
		/* Myfile<<"L"<<" "<<k<<" "<<O(11)<<endl; */
		/* Myfile<<"M"<<" "<<k<<" "<<O(12)<<endl; */
		/* Myfile<<"N"<<" "<<k<<" "<<O(13)<<endl; */
		/* Myfile<<"O"<<" "<<k<<" "<<O(14)<<endl; */
		/* Myfile<<"P"<<" "<<k<<" "<<O(15)<<endl; */
		/* Myfile<<"Q"<<" "<<k<<" "<<O(16)<<endl; */
		/* Myfile<<"R"<<" "<<k<<" "<<O(17)<<endl; */

	}

	Myfile.close();
	return 0;
}
