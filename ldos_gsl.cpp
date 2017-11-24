#include <iostream>
#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "TB.h"
#include "cunningham_spawn.h"
#include <gsl/gsl_integration.h>

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
typedef Matrix<dcomp, 9, 9> dmat;
typedef Vector3d vec;

struct a_struct{
	Vector3d d_1, d_2, d_3, d_4, d_5, d_6, d_7, d_8, d_9;
	Vector3d d_10, d_11, d_12, d_13, d_14, d_15, d_16;
	Vector3d d_17, d_18;
	
	double a = 1.;
	
	Matrix<dcomp, 9, 9> u, E;
	Matrix<dcomp, 9, 9> t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9;
	Matrix<dcomp, 9, 9> t_10, t_11, t_12, t_13, t_14, t_15, t_16;
	Matrix<dcomp, 9, 9> t_17, t_18;
	double j;
	};


double greens(double k_x, double k_z, double a, double k_y, dcomp omega, dmat &u, dmat &t_1,
		dmat &t_2, dmat &t_3, dmat &t_4, dmat &t_5, dmat &t_6, dmat &t_7, 
		dmat &t_8, dmat &t_9, dmat &t_10, dmat &t_11, dmat &t_12, dmat &t_13,
	  	dmat &t_14, dmat &t_15, dmat &t_16, dmat &t_17, dmat &t_18, vec &d_1,
	       	vec &d_2, vec &d_3, vec &d_4, vec &d_5, vec &d_6, vec &d_7, vec &d_8,
	       	vec &d_9, vec &d_10, vec &d_11, vec &d_12, vec &d_13, vec &d_14, vec &d_15,
	       	vec &d_16, vec &d_17, vec &d_18){

	dcomp i;
	i = -1.;
	i = sqrt(i);

	Vector3d K;
	K(0) = k_x;
	K(1) = k_y;
	K(2) = k_z;
	dmat E;
	E = u + t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K))
			+ t_9*exp(i*d_9.dot(K)) + t_10*exp(i*d_10.dot(K))
			+ t_11*exp(i*d_11.dot(K)) + t_12*exp(i*d_12.dot(K))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	dmat G;
      	Matrix<complex<double>, 9, 9> I = Matrix<complex<double>, 9, 9>::Identity();
	G = (omega*I - E).inverse();
	return (G.imag()).trace();
}

double temp(double k_y, void * alpha) {
	a_struct params = *(a_struct *) alpha;
	dcomp i, im, Ec;
	i = -1;
	i = sqrt(i);
	im = 1e-4*i;	
	double E = params.j;
	Ec = E + im;
	double result;
	result = kspace(&greens, 2, 1e-2, 2.*params.a, k_y, Ec, params.u, params.t_1, params.t_2, params.t_3,
		       	params.t_4, params.t_5, params.t_6, params.t_7, params.t_8, params.t_9,
			params.t_10, params.t_11, params.t_12, params.t_13, params.t_14, params.t_15,
		       	params.t_16, params.t_17, params.t_18, params.d_1, params.d_2, params.d_3, params.d_4,
			params.d_5, params.d_6, params.d_7, params.d_8, params.d_9, params.d_10, params.d_11,
		       	params.d_12, params.d_13, params.d_14, params.d_15, params.d_16, params.d_17, params.d_18);
	return result;
}

int main(){

	a_struct params;
	cout<<"Name the data file\n";
	string Mydata;
	getline(cin, Mydata);
	ofstream Myfile;	
	Mydata += ".txt";
	Myfile.open( Mydata.c_str(),ios::trunc );
	
	//position vectors of nearest neighbours in fcc
	params.d_1 << params.a, params.a, 0;
	params.d_2 << -params.a, -params.a, 0;
	params.d_3 << params.a, 0, params.a;
	params.d_4 << -params.a, 0, -params.a;
	params.d_5 << 0, params.a, params.a;
	params.d_6 << 0, -params.a, -params.a;
	params.d_7 << -params.a, params.a, 0;
	params.d_8 << params.a, -params.a, 0;
	params.d_9 << -params.a, 0, params.a;
	params.d_10 << params.a, 0, -params.a;
	params.d_11 << 0, -params.a, params.a;
	params.d_12 << 0, params.a, -params.a;

	//position vectors of next nearest neighbours
	params.d_13 << 2*params.a, 0, 0;
	params.d_14 << -2*params.a, 0, 0;
	params.d_15 << 0, 2*params.a, 0;
	params.d_16 << 0, -2*params.a, 0;
	params.d_17 << 0, 0, 2*params.a;
	params.d_18 << 0, 0, -2*params.a;

	//initialise onsite anparams.d hopping matrices for each nn
	params.u = TB(2, 0, 0, 8, params.d_1);
	params.t_1 = TB(2, 1, 0, 8, params.d_1);
	params.t_2 = TB(2, 1, 0, 8, params.d_2);
	params.t_3 = TB(2, 1, 0, 8, params.d_3);
	params.t_4 = TB(2, 1, 0, 8, params.d_4);
	params.t_5 = TB(2, 1, 0, 8, params.d_5);
	params.t_6 = TB(2, 1, 0, 8, params.d_6);
	params.t_7 = TB(2, 1, 0, 8, params.d_7);
	params.t_8 = TB(2, 1, 0, 8, params.d_8);
	params.t_9 = TB(2, 1, 0, 8, params.d_9);
	params.t_10 = TB(2, 1, 0, 8, params.d_10);
	params.t_11 = TB(2, 1, 0, 8, params.d_11);
	params.t_12 = TB(2, 1, 0, 8, params.d_12);

	params.t_13 = TB(2, 1, 1, 8, params.d_13);
	params.t_14 = TB(2, 1, 1, 8, params.d_14);
	params.t_15 = TB(2, 1, 1, 8, params.d_15);
	params.t_16 = TB(2, 1, 1, 8, params.d_16);
	params.t_17 = TB(2, 1, 1, 8, params.d_17);
	params.t_18 = TB(2, 1, 1, 8, params.d_18);

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	double result, error;
	size_t evals;

	gsl_function F;
	F.function = &temp;
	F.params = &params;
	const double start = 0;
	const double end = M_PI/(2.*params.a);
	double k_y = 0;
	for (params.j = 0.2; params.j<0.701; params.j=params.j+0.001){
		gsl_integration_qags(&F, start, end, 0, 1e-2, 1000, w, &result, &error);
		cout<<error<<endl;
		result = -result*2.*params.a*params.a*params.a/(M_PI*M_PI*M_PI*M_PI);
		Myfile<<params.j<<" "<<result<<endl;
	}


	Myfile.close();
	return 0;
}
