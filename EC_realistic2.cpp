#include <iostream>
#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "TB.h"
#include "cunningham_spawn.h"

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

dcomp greens(double k_x, double k_z, double a, dcomp omega, int N, dmat &u,
		dmat &t_1, dmat &t_3, dmat &t_4, dmat &t_5, dmat &t_7, 
		dmat &t_9, dmat &t_10, dmat &t_12, dmat &t_13,
	  	dmat &t_14, dmat &t_15, dmat &t_16, dmat &t_17, dmat &t_18, dmat &u_u, 
		dmat &tu_1, dmat &tu_4, dmat &tu_5,
	       	dmat &tu_7, dmat &tu_13, dmat &tu_14, dmat &tu_15, dmat &tu_16,
	       	dmat &tu_17, dmat &tu_18, dmat &u_d, dmat &td_1,
	       	dmat &td_4, dmat &td_5, dmat &td_7, dmat &td_13, 
		dmat &td_14, dmat &td_15, dmat &td_16, dmat &td_17, dmat &td_18, 
		vec &d_3, vec &d_4,
	       	vec &d_9, vec &d_10, vec &d_13, vec &d_14,
	       	vec &d_17, vec &d_18){

	dcomp i;
	i = -1.;
	i = sqrt(i);
	double k_y = 0;

	Vector3d K;
	K(0) = k_x;
	K(1) = k_y;
	K(2) = k_z;

	//construct diagonalised in-plane matrices
	//Cu fcc
	Matrix<dcomp, 9, 9> u_11, u_12;
	u_11 = u + t_3*exp(i*d_3.dot(K))+ t_4*exp(i*d_4.dot(K))+ t_9*exp(i*d_9.dot(K)) + t_10*exp(i*d_10.dot(K)) + 
		t_13*exp(i*d_13.dot(K))+ t_14*exp(i*d_14.dot(K))+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	u_12 = t_1 + t_5*exp(i*d_9.dot(K)) + t_7*exp(i*d_14.dot(K)) + t_12*exp(i*d_4.dot(K));
	Matrix<dcomp, 18, 18> U, T, OM;
	U << u_11, u_12, u_12.adjoint(), u_11;
	Matrix<complex<double>, 9, 9> zero = Matrix<complex<double>, 9, 9>::Zero();
	T << t_15, zero, u_12.adjoint(), t_16;

	//Co spin up bcc
	Matrix<dcomp, 9, 9> uu_11, uu_12;
	uu_11 = u_u + tu_13*exp(i*d_13.dot(K))+ tu_14*exp(i*d_14.dot(K))+ tu_17*exp(i*d_17.dot(K)) + tu_18*exp(i*d_18.dot(K));
	uu_12 = tu_1 + tu_5*exp(i*d_14.dot(K)) + tu_7*exp(i*d_18.dot(K)) + tu_4*exp(i*(d_18 + d_14).dot(K));
	Matrix<dcomp, 18, 18> Uu, Tu, OMu, GLu, GRu;
	Uu << uu_11, uu_12, uu_12.adjoint(), uu_11;
	Tu << tu_15, zero, uu_12.adjoint(), tu_16;

	//Co spin up bcc
	Matrix<dcomp, 9, 9> ud_11, ud_12;
	ud_11 = u_d + td_13*exp(i*d_13.dot(K))+ td_14*exp(i*d_14.dot(K))+ td_17*exp(i*d_17.dot(K)) + td_18*exp(i*d_18.dot(K));
	ud_12 = td_1 + td_5*exp(i*d_14.dot(K)) + td_7*exp(i*d_18.dot(K)) + td_4*exp(i*(d_18 + d_14).dot(K));
	Matrix<dcomp, 18, 18> Ud, Td, OMd, GLd, GRd;
	Ud << ud_11, ud_12, ud_12.adjoint(), ud_11;
	Td << td_15, zero, ud_12.adjoint(), td_16;

      	Matrix<complex<double>, 18, 18> I = Matrix<complex<double>, 18, 18>::Identity();
	ddmat Tudagg, Tddagg, Tdagg;
	Tudagg = Tu.adjoint();
	Tddagg = Td.adjoint();
	Tdagg = T.adjoint();

	OM = omega*I-U;
	OMu = omega*I-Uu;
	OMd = omega*I-Ud;

	GLu = gs(OMu, Tu);
	GLd = gs(OMd, Td);
	GRu = gs(OMu, Tudagg);
	GRd = gs(OMd, Tddagg);

	ddmat Rsigma_0_u, Rsigma_0_d, Rsigma_PI_u, Rsigma_PI_d;
	dcomp Fsigma;

//mobius transformation layer 2 from layer 1 to spacer thickness, N
	ddmat Tinv;
	Tinv = T.inverse();
	ddmat zero2 = ddmat::Zero();
	Matrix<dcomp, 36, 36> X,O,A,Oinv,OAOinv;
	X << 	zero2,	Tinv,
		-Tdagg,	OM*Tinv;
	ComplexEigenSolver<Matrix<dcomp, 36, 36>> ces;
	ces.compute(X);
	O = ces.eigenvectors();
	A = ces.eigenvalues().asDiagonal();
	Oinv = O.inverse();
	ddmat GNu, GNd, aa, b, c, d;
	OAOinv = O*A.array().pow(N+1).matrix()*Oinv;
	aa = OAOinv.topLeftCorner(18, 18);
	b = OAOinv.topRightCorner(18, 18);
	c = OAOinv.bottomLeftCorner(18, 18);
	d = OAOinv.bottomRightCorner(18, 18);
	GNu = (aa*GLu + b)*(c*GLu + d).inverse();
	GNd = (aa*GLd + b)*(c*GLd + d).inverse();
	Rsigma_0_u = (I-GRu*Tdagg*GNu*T);
	Rsigma_0_d = (I-GRd*Tdagg*GNd*T);
	Rsigma_PI_u = (I-GRd*Tdagg*GNu*T);
	Rsigma_PI_d = (I-GRu*Tdagg*GNd*T);
	Fsigma = (1./M_PI)*log((Rsigma_0_d*Rsigma_0_u*Rsigma_PI_u.inverse()*Rsigma_PI_d.inverse()).determinant());

	return Fsigma;
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
	Vector3d c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8;
	
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

	//position vectors of nearest neighbours in bcc
	c_1 << a, a, a;
	c_2 << -a, -a, -a;
	c_3 << a, -a, a;
	c_4 << -a, a, -a;
	c_5 << -a, a, a;
	c_6 << a, -a, -a;
	c_7 << a, a, -a;
	c_8 << -a, -a, a;

	//position vectors of next nearest neighbours
	d_13 << 2*a, 0, 0;
	d_14 << -2*a, 0, 0;
	d_15 << 0, 2*a, 0;
	d_16 << 0, -2*a, 0;
	d_17 << 0, 0, 2*a;
	d_18 << 0, 0, -2*a;

	//initialise onsite for fcc Cu
	Matrix<dcomp, 9, 9> u;
	u = TB(2, 0, 0, 8, d_1);

	//initialise nn hopping for fcc Cu
	Matrix<dcomp, 9, 9> t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9;
	Matrix<dcomp, 9, 9> t_10, t_11, t_12;
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

	//initialise next nn hopping for fcc Cu
	Matrix<dcomp, 9, 9> t_13, t_14, t_15, t_16, t_17, t_18;
	t_13 = TB(2, 1, 1, 8, d_13);
	t_14 = TB(2, 1, 1, 8, d_14);
	t_15 = TB(2, 1, 1, 8, d_15);
	t_16 = TB(2, 1, 1, 8, d_16);
	t_17 = TB(2, 1, 1, 8, d_17);
	t_18 = TB(2, 1, 1, 8, d_18);

	//initialise onsite for bcc Co spin up
	Matrix<dcomp, 9, 9> u_u;
	u_u = TB(0, 0, 0, 8, d_1);

	//initialise nn hopping for bcc Co spin up
	Matrix<dcomp, 9, 9> tu_1, tu_2, tu_3, tu_4, tu_5, tu_6, tu_7, tu_8;
	tu_1 = TB(0, 1, 0, 8, c_1);
	tu_2 = TB(0, 1, 0, 8, c_2);
	tu_3 = TB(0, 1, 0, 8, c_3);
	tu_4 = TB(0, 1, 0, 8, c_4);
	tu_5 = TB(0, 1, 0, 8, c_5);
	tu_6 = TB(0, 1, 0, 8, c_6);
	tu_7 = TB(0, 1, 0, 8, c_7);
	tu_8 = TB(0, 1, 0, 8, c_8);

	//initialise next nn hopping for bcc Co spin up
	Matrix<dcomp, 9, 9> tu_13, tu_14, tu_15, tu_16, tu_17, tu_18;
	tu_13 = TB(0, 1, 1, 8, d_13);
	tu_14 = TB(0, 1, 1, 8, d_14);
	tu_15 = TB(0, 1, 1, 8, d_15);
	tu_16 = TB(0, 1, 1, 8, d_16);
	tu_17 = TB(0, 1, 1, 8, d_17);
	tu_18 = TB(0, 1, 1, 8, d_18);

	//initialise onsite for bcc Co spin down 
	Matrix<dcomp, 9, 9> u_d;
	u_d = TB(1, 0, 0, 8, d_1);

	//initialise nn hopping for bcc Co spin down
	Matrix<dcomp, 9, 9> td_1, td_2, td_3, td_4, td_5, td_6, td_7, td_8;
	td_1 = TB(1, 1, 0, 8, c_1);
	td_2 = TB(1, 1, 0, 8, c_2);
	td_3 = TB(1, 1, 0, 8, c_3);
	td_4 = TB(1, 1, 0, 8, c_4);
	td_5 = TB(1, 1, 0, 8, c_5);
	td_6 = TB(1, 1, 0, 8, c_6);
	td_7 = TB(1, 1, 0, 8, c_7);
	td_8 = TB(1, 1, 0, 8, c_8);

	//initialise next nn hopping for bcc Co spin down
	Matrix<dcomp, 9, 9> td_13, td_14, td_15, td_16, td_17, td_18;
	td_13 = TB(1, 1, 1, 8, d_13);
	td_14 = TB(1, 1, 1, 8, d_14);
	td_15 = TB(1, 1, 1, 8, d_15);
	td_16 = TB(1, 1, 1, 8, d_16);
	td_17 = TB(1, 1, 1, 8, d_17);
	td_18 = TB(1, 1, 1, 8, d_18);

	dcomp i;
	i = -1.;
	i = sqrt(i);

	//number of principle layers of spacer
	const int N = 50;

	dcomp E = 0.;
	const double Ef = 0.0;
	const double kT = 8.617342857e-5*300/13.6058;
	dcomp result_complex;
	Myfile<<"N , Gamma"<<endl;
	for (int ii=0; ii < N ; ++ii)
	{
		result_complex = 0.;
		for (int j=0; j!=10; j++){
			E = Ef + (2.*j + 1.)*kT*M_PI*i;
			result_complex = result_complex + kspace(&greens, 2, 5e-2, 2*a, E, ii,
					u, t_1, t_3, t_4, t_5, t_7, t_9,
					t_10, t_12, t_13, t_14, t_15, t_16, t_17, t_18,
					u_u, tu_1, tu_4, tu_5, tu_7,
					tu_13, tu_14, tu_15, tu_16, tu_17, tu_18,
					u_d, td_1, td_4, td_5, td_7,
					td_13, td_14, td_15, td_16, td_17, td_18,
					d_3, d_4, d_9, d_10,
				       	d_13, d_14, d_17, d_18);

		}
		double result = real(result_complex);

		result *= kT/(4.*M_PI*M_PI);

		Myfile << ii+1 <<" ,  "<< -2.*M_PI*result << endl;
	}

	cout<<"finished!"<<endl;

			Myfile.close();
	return 0;
}
