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
	Matrix<dcomp, 9, 9> u_11, u_12, u_21, T_21;
	u_11 = u + t_3*exp(i*d_3.dot(K))+ t_4*exp(i*d_4.dot(K))+ t_9*exp(i*d_9.dot(K)) + t_10*exp(i*d_10.dot(K)) + 
		t_13*exp(i*d_13.dot(K))+ t_14*exp(i*d_14.dot(K))+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	u_12 = t_1 + t_5*exp(i*d_9.dot(K)) + t_7*exp(i*d_14.dot(K)) + t_12*exp(i*d_4.dot(K));
	u_21 = t_2 + t_8*exp(i*d_13.dot(K)) + t_11*exp(i*d_3.dot(K)) + t_6*exp(i*d_10.dot(K));
	Matrix<dcomp, 18, 18> U, OM, GL, GR, GN, GRinv, GNinv;
	U << u_11, u_12, u_21, u_11;
	Matrix<complex<double>, 9, 9> zero = Matrix<complex<double>, 9, 9>::Zero();
	T_21 = t_7 + t_1*exp(i*d_13.dot(K)) + t_5*exp(i*d_3.dot(K)) + t_12*exp(i*d_10.dot(K));

	//Co spin up fcc
	Matrix<dcomp, 9, 9> uu_11, uu_12, uu_21, Tu_21;
	uu_11 = u_u + tu_3*exp(i*d_3.dot(K))+ tu_4*exp(i*d_4.dot(K))+ tu_9*exp(i*d_9.dot(K)) + tu_10*exp(i*d_10.dot(K)) + 
		tu_13*exp(i*d_13.dot(K))+ tu_14*exp(i*d_14.dot(K))+ tu_17*exp(i*d_17.dot(K)) + tu_18*exp(i*d_18.dot(K));
	uu_12 = tu_1 + tu_5*exp(i*d_9.dot(K)) + tu_7*exp(i*d_14.dot(K)) + tu_12*exp(i*d_4.dot(K));
	uu_21 = tu_2 + tu_8*exp(i*d_13.dot(K)) + tu_11*exp(i*d_3.dot(K)) + tu_6*exp(i*d_10.dot(K));
	Matrix<dcomp, 18, 18> Uu, Tu_B, Tu_C, OMu, GLu, GRu_B, GRu_C;
	Uu << uu_11, uu_12, uu_21, uu_11;
	Tu_21 = tu_7 + tu_1*exp(i*d_13.dot(K)) + tu_5*exp(i*d_3.dot(K)) + tu_12*exp(i*d_10.dot(K));
	Tu_B << tu_15, zero, Tu_21, tu_15;
	Tu_C << tu_15, zero, uu_12, tu_15;

	//Co spin up fcc
	Matrix<dcomp, 9, 9> ud_11, ud_12, ud_21, Td_21;
	ud_11 = u_d + td_3*exp(i*d_3.dot(K))+ td_4*exp(i*d_4.dot(K))+ td_9*exp(i*d_9.dot(K)) + td_10*exp(i*d_10.dot(K)) + 
		td_13*exp(i*d_13.dot(K))+ td_14*exp(i*d_14.dot(K))+ td_17*exp(i*d_17.dot(K)) + td_18*exp(i*d_18.dot(K));
	ud_12 = td_1 + td_5*exp(i*d_9.dot(K)) + td_7*exp(i*d_14.dot(K)) + td_12*exp(i*d_4.dot(K));
	ud_21 = td_2 + td_8*exp(i*d_13.dot(K)) + td_11*exp(i*d_3.dot(K)) + td_6*exp(i*d_10.dot(K));
	Matrix<dcomp, 18, 18> Ud, Td_B, Td_C, OMd, GLd, GRd_B, GRd_C;
	Ud << ud_11, ud_12, ud_21, ud_11;
	Td_21 = td_7 + td_1*exp(i*d_13.dot(K)) + td_5*exp(i*d_3.dot(K)) + td_12*exp(i*d_10.dot(K));
	Td_B << td_15, zero, Td_21, td_15;
	Td_C << td_15, zero, ud_12, td_15;

	/* cout<<u_12<<endl<<endl; */
	/* cout<<T_21<<endl<<endl; */
	/* cout<<" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "<<endl; */

      	Matrix<complex<double>, 18, 18> I = Matrix<complex<double>, 18, 18>::Identity();
	ddmat Tu_Bdagg, Td_Bdagg, Tu_Cdagg, Td_Cdagg;
	Tu_Bdagg = Tu_B.adjoint();
	Td_Bdagg = Td_B.adjoint();
	Tu_Cdagg = Tu_C.adjoint();
	Td_Cdagg = Td_C.adjoint();

	dmat miniOM;
	dmat small_I = dmat::Identity();
	miniOM = omega*small_I - u_11; 
	OM = omega*I-U;
	OMu = omega*I-Uu;
	OMd = omega*I-Ud;

	GLu = gs(OMu, Tu_B);
	GLd = gs(OMd, Td_B);
	GRu_B = gs(OMu, Tu_Bdagg);
	GRd_B = gs(OMd, Td_Bdagg);
	GRu_C = (OMu -Tu_C*GRu_B*Tu_Cdagg).inverse();
	GRd_C = (OMd -Td_C*GRd_B*Td_Cdagg).inverse();

	/* GRu_C = gs(OMu, Tu_Cdagg); */
	/* GRd_C = gs(OMd, Td_Cdagg); */

	ddmat Rsigma_0_u, Rsigma_0_d, Rsigma_PI_u, Rsigma_PI_d;
	dcomp Fsigma;
	VectorXcd result(N);
	result.fill(0.);

//mobius transformation layer 2 from layer 1 to spacer thickness, N
	ddmat T_B, T_C;
	T_B << t_15, zero, T_21, t_15;
	T_C << t_15, zero, u_12, t_15;
	ddmat T_Binv, T_Cinv;
	T_Binv = T_B.inverse();
	T_Cinv = T_C.inverse();
	ddmat T_Bdagg, T_Cdagg;
	T_Bdagg = T_B.adjoint();
	T_Cdagg = T_C.adjoint();

	/* ddmat zero2 = ddmat::Zero(); */
	/* Matrix<dcomp, 36, 36> X_B, X_C; */
	/* Matrix<dcomp, 36, 1> A; */

	/* X_B << 	zero2,		T_Binv, */
	/* 	-T_Bdagg,	OM*T_Binv; */

	/* X_C << 	zero2,		T_Cinv, */
	/* 	-T_Cdagg,	OM*T_Cinv; */

	/* dmat t_15inv; */
	/* t_15inv = t_15.inverse(); */
	/* ddmat a0, b0, c0, d0, e0, f0; */
	/* a0 << zero, small_I, zero, zero; */
	/* b0 << zero, zero, t_15inv, zero; */
	/* c0 << zero, zero, -t_15.adjoint(), -T_21.adjoint(); */
	/* d0 << -T_21*t_15inv, small_I, miniOM*t_15inv, zero; */ 
	/* dddmat Xi_B, Xi_C; */
	/* Xi_B << a0, b0, c0, d0; */
	/* e0 << zero, zero, -t_15.adjoint(), -u_12.adjoint(); */
	/* f0 << -u_12*t_15inv, small_I, miniOM*t_15inv, zero; */ 
	/* Xi_C << a0, b0, e0, f0; */

	/* ComplexEigenSolver<Matrix<dcomp, 36, 36>> ces; */
	/* ces.compute(X); */
	/* O = ces.eigenvectors(); */
	/* A = ces.eigenvalues(); */
	/* Oinv = O.inverse(); */

	/* ddmat GNu, GNd, a1, b1, c1, d1; */
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

	/* dmat t_15inv; */
	/* t_15inv = t_15.inverse(); */
	/* ddmat GNu, GNd, a0, b0, c0, d0; */
	/* a0 << zero, small_I, zero, zero; */
	/* b0 << zero, zero, t_15inv, zero; */
	/* c0 << zero, zero, -t_15.adjoint(), -T_21.adjoint(); */
	/* d0 << -T_21*t_15inv, small_I, miniOM*t_15inv, zero; */ 
	/* dddmat stack; */
	/* stack << a0, b0, c0, d0; */
	/* ComplexEigenSolver<Matrix<dcomp, 36, 36>> ces; */
	/* ces.compute(stack); */
	/* Matrix<dcomp, 36, 36>O,Oinv; */
	/* Matrix<dcomp, 36, 1> A; */

	/* ddmat zero2 = ddmat::Zero(); */
	/* Matrix<dcomp, 36, 36> X; */
	/* X << 	zero2,	Tinv, */
	/* 	-Tdagg,	OM*Tinv; */
	/* ComplexEigenSolver<Matrix<dcomp, 36, 36>> ces2; */
	/* ces2.compute(X); */
	/* O = ces2.eigenvectors(); */

	/* O = ces.eigenvectors(); */
	/* A = ces.eigenvalues(); */
	/* /1* for (int it = 0; it < 36; it++) *1/ */
	/* 	/1* cout<<sqrt(real(A(it))*real(A(it)) + imag(A(it))*imag(A(it)))<<endl; *1/ */
	/* Oinv = O.inverse(); */
	/* ddmat a1, b1, c1, d1, fu, fd, a2, b2, c2, d2, tmpu, tmpd; */
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

	/* dddmat tmp; */
	/* tmp = Xi_C; */

	/* for (int it=1; it < N; ++it){ */
	/* 	a1 = tmp.topLeftCorner(18, 18); */
	/* 	b1 = tmp.topRightCorner(18, 18); */
	/* 	c1 = tmp.bottomLeftCorner(18, 18); */
	/* 	d1 = tmp.bottomRightCorner(18, 18); */
	/* 	GNu = (a1*GLu + b1)*(c1*GLu + d1).inverse(); */
	/* 	GNd = (a1*GLd + b1)*(c1*GLd + d1).inverse(); */
	/* 	if (it%2 == 0){ */
	/* 		tmp = tmp*Xi_C; */
	/* 		Rsigma_0_u = (I-GRu*T_Bdagg*GNu*T_B); */
	/* 		Rsigma_0_d = (I-GRd*T_Bdagg*GNd*T_B); */
	/* 		Rsigma_PI_u = (I-GRd*T_Bdagg*GNu*T_B); */
	/* 		Rsigma_PI_d = (I-GRu*T_Bdagg*GNd*T_B); */
	/* 	} */	
	/* 	if (it%2 == 1){ */
	/* 		tmp = tmp*Xi_B; */
	/* 		Rsigma_0_u = (I-GRu*T_Cdagg*GNu*T_C); */
	/* 		Rsigma_0_d = (I-GRd*T_Cdagg*GNd*T_C); */
	/* 		Rsigma_PI_u = (I-GRd*T_Cdagg*GNu*T_C); */
	/* 		Rsigma_PI_d = (I-GRu*T_Cdagg*GNd*T_C); */
	/* 	} */	
	/* 	Fsigma = (1./M_PI)*log((Rsigma_0_d*Rsigma_0_u*Rsigma_PI_u.inverse()*Rsigma_PI_d.inverse()).determinant()); */
	/* 	result[it] = Fsigma; */
	/* } */

	ddmat GNu, GNd;

	GNu = GLu;
	GNd = GLd;
	Rsigma_0_u = (I-GRu_B*T_Bdagg*GNu*T_B);
	Rsigma_0_d = (I-GRd_B*T_Bdagg*GNd*T_B);
	Rsigma_PI_u = (I-GRd_B*T_Bdagg*GNu*T_B);
	Rsigma_PI_d = (I-GRu_B*T_Bdagg*GNd*T_B);

	Fsigma = (1./M_PI)*log((Rsigma_0_d*Rsigma_0_u*Rsigma_PI_u.inverse()*Rsigma_PI_d.inverse()).determinant());
	result[0] = Fsigma;

//adlayer layer 2 from layer 1 to spacer thickness, N
	for (int it=1; it < N; ++it){
		if (it%2 == 0){
			GLu = (OM -T_Bdagg*GLu*T_B).inverse();
			GLd = (OM -T_Bdagg*GLd*T_B).inverse();
			Rsigma_0_u = (I-GRu_B*T_Bdagg*GLu*T_B);
			Rsigma_0_d = (I-GRd_B*T_Bdagg*GLd*T_B);
			Rsigma_PI_u = (I-GRd_B*T_Bdagg*GLu*T_B);
			Rsigma_PI_d = (I-GRu_B*T_Bdagg*GLd*T_B);
		}
		if (it%2 == 1){
			GLu = (OM -T_Cdagg*GLu*T_C).inverse();
			GLd = (OM -T_Cdagg*GLd*T_C).inverse();
			Rsigma_0_u = (I-GRu_C*T_Cdagg*GLu*T_C);
			Rsigma_0_d = (I-GRd_C*T_Cdagg*GLd*T_C);
			Rsigma_PI_u = (I-GRd_C*T_Cdagg*GLu*T_C);
			Rsigma_PI_d = (I-GRu_C*T_Cdagg*GLd*T_C);
		}
		Fsigma = (1./M_PI)*log((Rsigma_0_d*Rsigma_0_u*Rsigma_PI_u.inverse()*Rsigma_PI_d.inverse()).determinant());
		result[it] = Fsigma;
	}

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
	const int N = 30;

	dcomp E = 0.;
	/* const double Ef = 0.5805; */
	const double Ef = 0.57553;
	/* const double Ef = -0.038; */
	const double kT = 8.617342857e-5*315.79/13.6058;
	VectorXcd result_complex(N);
	E = Ef + kT*M_PI*i;
	cout<<E<<endl;
	result_complex = greens(-2.19911485751286, -0.314159265358979, a, E, N,
			u, t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9,
			t_10, t_11, t_12, t_13, t_14, t_15, t_16, t_17, t_18,
			u_u, tu_1, tu_2, tu_3, tu_4, tu_5, tu_6, tu_7, tu_8, tu_9,
		       	tu_10, tu_11, tu_12, tu_13, tu_14, tu_15, tu_16, tu_17, tu_18,
			u_d, td_1, td_2, td_3, td_4, td_5, td_6, td_7, td_8, td_9,
		       	td_10, td_11, td_12, td_13, td_14, td_15, td_16, td_17, td_18,
			d_3, d_4, d_9, d_10,
		       	d_13, d_14, d_17, d_18);

	VectorXd result = result_complex.real();

	/* result *= 1/(4.*M_PI*M_PI); */
	Myfile<<"N , Gamma"<<endl;

	for (int ii=0; ii < N ; ++ii){
		/* Myfile << (ii+1)/10. <<" ,  "<< -2.*M_PI*result[ii] << endl; */
		Myfile << (ii) <<" "<< result[ii] << endl;
		/* Myfile << ii+1 <<" ,  "<< -2.*M_PI*result[ii] << endl; */
	}

	cout<<"finished!"<<endl;

			Myfile.close();
	return 0;
}
