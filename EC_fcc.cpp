#include <iostream>
#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "TBdynamic.h"
#include "cunningham_spawn_diamond.h"
//This program calculates the realistic exchange coupling in a Co/Cu/Co(001)
//trilayer. It does so for bcc Co and fcc Cu. Interatomic spacing is considered
//identical and the interfaces are abrupt.
//use adlayer for integer thickness n, mobius method for fractional n.
//each method is arranged in a block for simplicity.

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
typedef Matrix<dcomp, 9, 9> dmat;
typedef Matrix<dcomp, 18, 18> ddmat;
typedef Matrix<dcomp, 36, 36> dddmat;
typedef Vector3d vec;

//calculates g at the interfaces
ddmat gs(dmat &OM, dmat &t_1, dmat &t_2)
{
	dmat I = dmat::Identity();
	dmat t_2inv;
	t_2inv = t_2.inverse();
	dmat zero = dmat::Zero();
	ddmat a0, b0, c0, d0;
	a0 << zero, I, zero, zero;
	b0 << zero, zero, t_2inv, zero;
	c0 << zero, zero, -t_2.adjoint(), -t_1.adjoint();
	d0 << -t_1*t_2inv, I, OM*t_2inv, zero; 
	dddmat stack, O;
	stack << a0, b0, c0, d0;
	ComplexEigenSolver<Matrix<dcomp, 36, 36>> ces;
	ces.compute(stack);
	O = ces.eigenvectors();
	ddmat b = O.topRightCorner(18, 18);
	ddmat d = O.bottomRightCorner(18, 18);
	ddmat GR;
	GR = b*d.inverse();
	return GR;
}

//prepares SGF for RH surfaces
ddmat gr(dmat &OM, dmat &t_1, dmat &t_2)
{
	dmat I = dmat::Identity();
	dmat t_2adj, t_2adjinv;
	t_2adj = t_2.adjoint();
	t_2adjinv = t_2adj.inverse();
	dmat zero = dmat::Zero();
	ddmat a0, b0, c0, d0;
	a0 << zero, zero, I, zero;
	b0 << zero, t_2adjinv, zero, zero;
	c0 << -t_1, -t_2, zero, zero; 
	d0 << zero, OM*t_2adjinv, I, -t_1.adjoint()*t_2adjinv; 
	dddmat stack, O;
	stack << a0, b0, c0, d0;
	ComplexEigenSolver<Matrix<dcomp, 36, 36>> ces;
	ces.compute(stack);
	O = ces.eigenvectors();
	ddmat b = O.topRightCorner(18, 18);
	ddmat d = O.bottomRightCorner(18, 18);
	ddmat GR;
	GR = b*d.inverse();
	return GR;
}

/* ddmat gs(ddmat &OM, ddmat &T) */
/* { */
/* 	ddmat zero = ddmat::Zero(); */
/* 	Matrix<dcomp, 36, 36> X,O; */
/* 	X << 	zero,	T.inverse(), */
/* 		-T.adjoint(),	OM*T.inverse(); */
/* 	ComplexEigenSolver<Matrix<dcomp, 36, 36>> ces; */
/* 	ces.compute(X); */
/* 	O = ces.eigenvectors(); */
/* 	ddmat b = O.topRightCorner(18, 18); */
/* 	ddmat d = O.bottomRightCorner(18, 18); */
/* 	ddmat GR; */
/* 	GR = b*d.inverse(); */
/* 	return GR; */
/* } */

VectorXcd greens(double k_x, double k_z, double a, dcomp omega, int N, dmat &u,
		dmat &t_1, dmat &t_2, dmat &t_3, dmat &t_4, dmat &t_5, dmat &t_6, dmat &t_7, 
		dmat &t_8, dmat &t_9, dmat &t_10, dmat &t_11, dmat &t_12, dmat &t_13,
	  	dmat &t_14, dmat &t_15, dmat &t_16, dmat &t_17, dmat &t_18, dmat &u_u, 
		dmat &tu_1, dmat &tu_2, dmat &tu_3, dmat &tu_4, dmat &tu_5, dmat &tu_6,
	       	dmat &tu_7, dmat &tu_8, dmat &tu_9, dmat &tu_10, dmat &tu_11, dmat &tu_12, dmat &tu_13,
		dmat &tu_14, dmat &tu_15, dmat &tu_16, dmat &tu_17, dmat &tu_18, dmat &u_d,
		dmat &td_1, dmat &td_2, dmat &td_3, dmat &td_4, dmat &td_5, dmat &td_6,
	       	dmat &td_7, dmat &td_8, dmat &td_9, dmat &td_10, dmat &td_11, dmat &td_12, dmat &td_13,
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
	Matrix<dcomp, 9, 9> u_11, T_21;
	u_11 = u + t_3*exp(i*d_3.dot(K))+ t_4*exp(i*d_4.dot(K))+ t_9*exp(i*d_9.dot(K)) + t_10*exp(i*d_10.dot(K)) + 
		t_13*exp(i*d_13.dot(K))+ t_14*exp(i*d_14.dot(K))+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* u_12 = t_1 + t_5*exp(i*d_9.dot(K)) + t_7*exp(i*d_14.dot(K)) + t_12*exp(i*d_4.dot(K)); */
	/* u_21 = t_2 + t_8*exp(i*d_13.dot(K)) + t_11*exp(i*d_3.dot(K)) + t_6*exp(i*d_10.dot(K)); */
	Matrix<dcomp, 18, 18> T, GL, GR, GN, GRinv, GNinv;
	/* U << u_11, u_12, u_21, u_11; */
	Matrix<complex<double>, 9, 9> zero = Matrix<complex<double>, 9, 9>::Zero();
	T_21 = t_7 + t_1*exp(i*d_13.dot(K)) + t_5*exp(i*d_3.dot(K)) + t_12*exp(i*d_10.dot(K));
	T << t_15, zero, T_21, t_15;

	//Co spin up fcc
	Matrix<dcomp, 9, 9> uu_11, Tu_21;
	uu_11 = u_u + tu_3*exp(i*d_3.dot(K))+ tu_4*exp(i*d_4.dot(K))+ tu_9*exp(i*d_9.dot(K)) + tu_10*exp(i*d_10.dot(K)) + 
		tu_13*exp(i*d_13.dot(K))+ tu_14*exp(i*d_14.dot(K))+ tu_17*exp(i*d_17.dot(K)) + tu_18*exp(i*d_18.dot(K));
	/* uu_12 = tu_1 + tu_5*exp(i*d_9.dot(K)) + tu_7*exp(i*d_14.dot(K)) + tu_12*exp(i*d_4.dot(K)); */
	/* uu_21 = tu_2 + tu_8*exp(i*d_13.dot(K)) + tu_11*exp(i*d_3.dot(K)) + tu_6*exp(i*d_10.dot(K)); */
	Matrix<dcomp, 18, 18> Tu, GLu, GRu;
	/* Uu << uu_11, uu_12, uu_21, uu_11; */
	Tu_21 = tu_7 + tu_1*exp(i*d_13.dot(K)) + tu_5*exp(i*d_3.dot(K)) + tu_12*exp(i*d_10.dot(K));
	Tu << tu_15, zero, Tu_21, tu_15;

	//Co spin up fcc
	Matrix<dcomp, 9, 9> ud_11, Td_21;
	ud_11 = u_d + td_3*exp(i*d_3.dot(K))+ td_4*exp(i*d_4.dot(K))+ td_9*exp(i*d_9.dot(K)) + td_10*exp(i*d_10.dot(K)) + 
		td_13*exp(i*d_13.dot(K))+ td_14*exp(i*d_14.dot(K))+ td_17*exp(i*d_17.dot(K)) + td_18*exp(i*d_18.dot(K));
	/* ud_12 = td_1 + td_5*exp(i*d_9.dot(K)) + td_7*exp(i*d_14.dot(K)) + td_12*exp(i*d_4.dot(K)); */
	/* ud_21 = td_2 + td_8*exp(i*d_13.dot(K)) + td_11*exp(i*d_3.dot(K)) + td_6*exp(i*d_10.dot(K)); */
	Matrix<dcomp, 18, 18>  Td, GLd, GRd;
	/* Ud << ud_11, ud_12, ud_21, ud_11; */
	Td_21 = td_7 + td_1*exp(i*d_13.dot(K)) + td_5*exp(i*d_3.dot(K)) + td_12*exp(i*d_10.dot(K));
	Td << td_15, zero, Td_21, td_15;

      	Matrix<complex<double>, 18, 18> I = Matrix<complex<double>, 18, 18>::Identity();
	ddmat Tudagg, Tddagg, Tdagg;
	Tudagg = Tu.adjoint();
	Tddagg = Td.adjoint();
	Tdagg = T.adjoint();

	dmat small_I = dmat::Identity();
	dmat OM, OMu, OMd;
	OM = omega*small_I - u_11; 
	OMu = omega*small_I-uu_11;
	OMd = omega*small_I-ud_11;

	GLu = gs(OMu, Tu_21, tu_15);
	GLd = gs(OMd, Td_21, td_15);
	GRu = gr(OMu, Tu_21, tu_15);
	GRd = gr(OMd, Td_21, td_15);

	ddmat Rsigma_0_u, Rsigma_0_d, Rsigma_PI_u, Rsigma_PI_d;
	dcomp Fsigma;
	VectorXcd result(N);
	result.fill(0.);

//mobius transformation layer 2 from layer 1 to spacer thickness, N
	/* ddmat Tinv; */
	/* Tinv = T.inverse(); */
	ddmat zero2 = ddmat::Zero();

	/* Matrix<dcomp, 36, 36> X,O,Oinv,OAOinv; */
	/* Matrix<dcomp, 36, 1> A; */
	/* X << 	zero2,	Tinv, */
	/* 	-Tdagg,	OM*Tinv; */
	/* ComplexEigenSolver<Matrix<dcomp, 36, 36>> ces; */
	/* ces.compute(X); */
	/* O = ces.eigenvectors(); */
	/* A = ces.eigenvalues(); */
	/* Oinv = O.inverse(); */
	/* ddmat GNu, GNd, a1, b1, c1, d1, fu, fd, a2, b2, c2, d2, tmpu, tmpd; */
	/* Matrix<dcomp, 18, 1> A1, A2; */
	/* A1 = A.topLeftCorner(18, 1); */
	/* A2 = A.bottomLeftCorner(18, 1); */
	/* a1 = Oinv.topLeftCorner(18, 18); */
	/* b1 = Oinv.topRightCorner(18, 18); */
	/* c1 = Oinv.bottomLeftCorner(18, 18); */
	/* d1 = Oinv.bottomRightCorner(18, 18); */
	/* fu = (a1*GLu + b1)*(c1*GLu + d1).inverse(); */
	/* fd = (a1*GLd + b1)*(c1*GLd + d1).inverse(); */
	/* a2 = O.topLeftCorner(18, 18); */
	/* b2 = O.topRightCorner(18, 18); */
	/* c2 = O.bottomLeftCorner(18, 18); */
	/* d2 = O.bottomRightCorner(18, 18); */

	dmat t_15inv;
	t_15inv = t_15.inverse();
	ddmat GNu, GNd, a0, b0, c0, d0;
	a0 << zero, small_I, zero, zero;
	b0 << zero, zero, t_15inv, zero;
	c0 << zero, zero, -t_15.adjoint(), -T_21.adjoint();
	d0 << -T_21*t_15inv, small_I, OM*t_15inv, zero; 
	dddmat stack, powstack;
	stack << a0, b0, c0, d0;
	ComplexEigenSolver<Matrix<dcomp, 36, 36>> ces;
	ces.compute(stack);
	Matrix<dcomp, 36, 36>O,Oinv,OAOinv;
	Matrix<dcomp, 36, 1> A;
	O = ces.eigenvectors();
	A = ces.eigenvalues();
	Oinv = O.inverse();
	ddmat a1, b1, c1, d1, fu, fd, a2, b2, c2, d2, tmpu, tmpd;
	Matrix<dcomp, 18, 1> A1, A2;
	A1 = A.topLeftCorner(18, 1);
	A2 = A.bottomLeftCorner(18, 1);
	a1 = Oinv.topLeftCorner(18, 18);
	b1 = Oinv.topRightCorner(18, 18);
	c1 = Oinv.bottomLeftCorner(18, 18);
	d1 = Oinv.bottomRightCorner(18, 18);
	fu = (a1*GLu + b1)*(c1*GLu + d1).inverse();
	fd = (a1*GLd + b1)*(c1*GLd + d1).inverse();
	a2 = O.topLeftCorner(18, 18);
	b2 = O.topRightCorner(18, 18);
	c2 = O.bottomLeftCorner(18, 18);
	d2 = O.bottomRightCorner(18, 18);

	for (int it=0; it < N; ++it){
		tmpu = (A1.array().pow(it/10.).matrix()).asDiagonal()*fu*(A2.array().pow(-it/10.).matrix()).asDiagonal();
		tmpd = (A1.array().pow(it/10.).matrix()).asDiagonal()*fd*(A2.array().pow(-it/10.).matrix()).asDiagonal();
		GNu = (a2*tmpu + b2)*(c2*tmpu + d2).inverse();
		GNd = (a2*tmpd + b2)*(c2*tmpd + d2).inverse();
		Rsigma_0_u = (I-GRu*Tdagg*GNu*T);
		Rsigma_0_d = (I-GRd*Tdagg*GNd*T);
		Rsigma_PI_u = (I-GRd*Tdagg*GNu*T);
		Rsigma_PI_d = (I-GRu*Tdagg*GNd*T);
		Fsigma = (1./M_PI)*log((Rsigma_0_d*Rsigma_0_u*Rsigma_PI_u.inverse()*Rsigma_PI_d.inverse()).determinant());
		result[it] = Fsigma;
	}

/* //adlayer layer 2 from layer 1 to spacer thickness, N */
/* 	for (int it=0; it < N; ++it){ */
/* 		GLu = (OM -Tdagg*GLu*T).inverse(); */
/* 		GLd = (OM -Tdagg*GLd*T).inverse(); */
/* 		Rsigma_0_u = (I-GRu*Tdagg*GLu*T); */
/* 		Rsigma_0_d = (I-GRd*Tdagg*GLd*T); */
/* 		Rsigma_PI_u = (I-GRd*Tdagg*GLu*T); */
/* 		Rsigma_PI_d = (I-GRu*Tdagg*GLd*T); */
/* 		Fsigma = (1./M_PI)*log((Rsigma_0_d*Rsigma_0_u*Rsigma_PI_u.inverse()*Rsigma_PI_d.inverse()).determinant()); */
/* 		result[it] = Fsigma; */
/* 	} */

	return result;
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

	//initialise onsite for fcc Cu
	Matrix<dcomp, 9, 9> u;
	u = TB(2, 0, 0, 9, d_1);

	//initialise nn hopping for fcc Cu
	Matrix<dcomp, 9, 9> t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9;
	Matrix<dcomp, 9, 9> t_10, t_11, t_12;
	t_1 = TB(2, 1, 0, 9, d_1);
	t_2 = TB(2, 1, 0, 9, d_2);
	t_3 = TB(2, 1, 0, 9, d_3);
	t_4 = TB(2, 1, 0, 9, d_4);
	t_5 = TB(2, 1, 0, 9, d_5);
	t_6 = TB(2, 1, 0, 9, d_6);
	t_7 = TB(2, 1, 0, 9, d_7);
	t_8 = TB(2, 1, 0, 9, d_8);
	t_9 = TB(2, 1, 0, 9, d_9);
	t_10 = TB(2, 1, 0, 9, d_10);
	t_11 = TB(2, 1, 0, 9, d_11);
	t_12 = TB(2, 1, 0, 9, d_12);

	//initialise next nn hopping for fcc Cu
	Matrix<dcomp, 9, 9> t_13, t_14, t_15, t_16, t_17, t_18;
	t_13 = TB(2, 1, 1, 9, d_13);
	t_14 = TB(2, 1, 1, 9, d_14);
	t_15 = TB(2, 1, 1, 9, d_15);
	t_16 = TB(2, 1, 1, 9, d_16);
	t_17 = TB(2, 1, 1, 9, d_17);
	t_18 = TB(2, 1, 1, 9, d_18);

	//initialise onsite for fcc Co spin up
	Matrix<dcomp, 9, 9> u_u;
	u_u = TB(0, 0, 0, 9, d_1);

	//initialise nn hopping for fcc Co spin up
	Matrix<dcomp, 9, 9> tu_1, tu_2, tu_3, tu_4, tu_5, tu_6, tu_7, tu_8;
	Matrix<dcomp, 9, 9> tu_9, tu_10, tu_11, tu_12;
	tu_1 = TB(0, 1, 0, 9, d_1);
	tu_2 = TB(0, 1, 0, 9, d_2);
	tu_3 = TB(0, 1, 0, 9, d_3);
	tu_4 = TB(0, 1, 0, 9, d_4);
	tu_5 = TB(0, 1, 0, 9, d_5);
	tu_6 = TB(0, 1, 0, 9, d_6);
	tu_7 = TB(0, 1, 0, 9, d_7);
	tu_8 = TB(0, 1, 0, 9, d_8);
	tu_9 = TB(0, 1, 0, 9, d_9);
	tu_10 = TB(0, 1, 0, 9, d_10);
	tu_11 = TB(0, 1, 0, 9, d_11);
	tu_12 = TB(0, 1, 0, 9, d_12);

	//initialise next nn hopping for fcc Co spin up
	Matrix<dcomp, 9, 9> tu_13, tu_14, tu_15, tu_16, tu_17, tu_18;
	tu_13 = TB(0, 1, 1, 9, d_13);
	tu_14 = TB(0, 1, 1, 9, d_14);
	tu_15 = TB(0, 1, 1, 9, d_15);
	tu_16 = TB(0, 1, 1, 9, d_16);
	tu_17 = TB(0, 1, 1, 9, d_17);
	tu_18 = TB(0, 1, 1, 9, d_18);

	//initialise onsite for fcc Co spin down 
	Matrix<dcomp, 9, 9> u_d;
	u_d = TB(1, 0, 0, 9, d_1);

	//initialise nn hopping for fcc Co spin down
	Matrix<dcomp, 9, 9> td_1, td_2, td_3, td_4, td_5, td_6, td_7, td_8;
	Matrix<dcomp, 9, 9> td_9, td_10, td_11, td_12;
	td_1 = TB(1, 1, 0, 9, d_1);
	td_2 = TB(1, 1, 0, 9, d_2);
	td_3 = TB(1, 1, 0, 9, d_3);
	td_4 = TB(1, 1, 0, 9, d_4);
	td_5 = TB(1, 1, 0, 9, d_5);
	td_6 = TB(1, 1, 0, 9, d_6);
	td_7 = TB(1, 1, 0, 9, d_7);
	td_8 = TB(1, 1, 0, 9, d_8);
	td_9 = TB(1, 1, 0, 9, d_9);
	td_10 = TB(1, 1, 0, 9, d_10);
	td_11 = TB(1, 1, 0, 9, d_11);
	td_12 = TB(1, 1, 0, 9, d_12);

	//initialise next nn hopping for fcc Co spin down
	Matrix<dcomp, 9, 9> td_13, td_14, td_15, td_16, td_17, td_18;
	td_13 = TB(1, 1, 1, 9, d_13);
	td_14 = TB(1, 1, 1, 9, d_14);
	td_15 = TB(1, 1, 1, 9, d_15);
	td_16 = TB(1, 1, 1, 9, d_16);
	td_17 = TB(1, 1, 1, 9, d_17);
	td_18 = TB(1, 1, 1, 9, d_18);

	dcomp i;
	i = -1.;
	i = sqrt(i);

	//number of principle layers of spacer
	/* const int N = 50; */
	const int N = 200;

	dcomp E = 0.;
	/* const double Ef = 0.5805; */
	const double Ef = 0.57553;
	/* const double Ef = -0.038; */
	const double kT = 8.617342857e-5*316/13.6058;
	VectorXcd result_complex(N);
	result_complex.fill(0.);
	for (int j=0; j!=10; j++){
		E = Ef + (2.*j + 1.)*kT*M_PI*i;
		result_complex = result_complex + kspace(&greens, 2, 5e-2, 2*a, E, N,
				u, t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9,
				t_10, t_11, t_12, t_13, t_14, t_15, t_16, t_17, t_18,
				u_u, tu_1, tu_2, tu_3, tu_4, tu_5, tu_6, tu_7, tu_8, tu_9,
			       	tu_10, tu_11, tu_12, tu_13, tu_14, tu_15, tu_16, tu_17, tu_18,
				u_d, td_1, td_2, td_3, td_4, td_5, td_6, td_7, td_8, td_9,
			       	td_10, td_11, td_12, td_13, td_14, td_15, td_16, td_17, td_18,
				d_3, d_4, d_9, d_10,
			       	d_13, d_14, d_17, d_18);

	}
	VectorXd result = result_complex.real();

	result *= kT/(4.*M_PI*M_PI);
	Myfile<<"N , Gamma"<<endl;

	for (int ii=0; ii < N ; ++ii){
		/* Myfile << (ii+1)/10. <<" ,  "<< -2.*M_PI*result[ii] << endl; */
		Myfile << (ii)/10. <<" "<< 4.*M_PI*result[ii] << endl;
		/* Myfile << ii+1 <<" ,  "<< -2.*M_PI*result[ii] << endl; */
	}

	cout<<"finished!"<<endl;

			Myfile.close();
	return 0;
}
