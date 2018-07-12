#include <iostream>
#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <vector>
#include "TB_sk.h"
#include <iomanip>
#include "cunningham_spawn_old.h"

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
		dmat &t_8, dmat &t_13, dmat &t_14, dmat &t_15, dmat &t_16, dmat &t_17, dmat &t_18,
		dmat &t_19, dmat &t_20, dmat &t_21, dmat &t_22, dmat &t_23, dmat &t_24, dmat &t_25, dmat &t_26, dmat &t_27, 
		dmat &t_28, dmat &t_29, dmat &t_30,
		vec &d_1, vec &d_2, vec &d_3, vec &d_4, vec &d_5, vec &d_6, vec &d_7, vec &d_8,
	       	vec &d_13, vec &d_14, vec &d_15, vec &d_16, vec &d_17, vec &d_18, vec &d_19, vec &d_20,
		vec &d_21, vec &d_22, vec &d_23, vec &d_24, vec &d_25, vec &d_26, vec &d_27, vec &d_28,
		vec &d_29, vec &d_30){

	dcomp i;
	i = -1.;
	i = sqrt(i);
	double k_y = 0;

	Vector3d K;
	K(0) = k_x;
	K(1) = k_y;
	K(2) = k_z;

	//construct diagonalised in-plane matrices
	Matrix<dcomp, 9, 9> u_11, u_12, u_21, T_21, T_11;
	u_11 = u + t_13*exp(i*d_13.dot(K))+ t_14*exp(i*d_14.dot(K))+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
		+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K)) + t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K));
	u_12 = t_1 + t_3*exp(i*d_14.dot(K)) + t_7*exp(i*d_18.dot(K)) + t_6*exp(i*(d_22).dot(K));
	u_21 = t_2 + t_4*exp(i*d_13.dot(K)) + t_8*exp(i*d_17.dot(K)) + t_5*exp(i*(d_13 + d_17).dot(K));
	Matrix<complex<double>, 9, 9> zero = Matrix<complex<double>, 9, 9>::Zero();
	Matrix<dcomp, 18, 18> U, T, Tdagg, OM, GL, GR, GN, GRinv, GNinv;
	U << u_11, u_12, u_21, u_11;
	/* T_21 = t_7 + t_1*exp(i*d_13.dot(K)) + t_5*exp(i*d_3.dot(K)) + t_12*exp(i*d_10.dot(K)); */
	T_11 = t_15 + t_19*exp(i*d_13.dot(K)) + t_23*exp(i*d_17.dot(K)) + t_26*exp(i*d_14.dot(K)) + t_30*exp(i*d_18.dot(K));
	/* T_21 = u_12; */
	T_21 = t_5 + t_1*exp(i*d_21.dot(K)) + t_7*exp(i*d_13.dot(K)) + t_3*exp(i*d_17.dot(K));
	T << T_11, zero, T_21, T_11;

      	Matrix<complex<double>, 18, 18> I = Matrix<complex<double>, 18, 18>::Identity();

	OM = omega*I-U;
	Tdagg = T.adjoint();

	GL = gs(OM, T);
	GR = gs(OM, Tdagg);
	GRinv = GR.inverse();
	GNinv = GRinv - Tdagg*GL*T;
	GN = GNinv.inverse();

	return imag(GN.trace());

}

Matrix<dcomp, 9, 9> convert(Matrix<dcomp, 9, 9> &t){
	Matrix<dcomp, 9, 9> result;
	//row ordering comes from inspecting the diagonalised onsite matrix
	//U, and comparing to Papaconstantopoulos
	vector<int> m;
	m.emplace_back(5);
	m.emplace_back(6);
	m.emplace_back(7);
	m.emplace_back(8);
	m.emplace_back(2);
	m.emplace_back(3);
	m.emplace_back(4);
	m.emplace_back(1);
	m.emplace_back(0);
	for (int i = 0; i<9; i++){
		for (int j = 0; j<9; j++){
			result(i,j) = t(m[i],m[j]);
		}
	}
	return result;
}

Matrix<dcomp, 9, 9> read(Vector3d &dvec, int ispin){
      string input;
      if (ispin == +1)
        input = "iron_up_hr.dat";
      if (ispin == -1)
        input = "iron_dn_hr.dat";
      ifstream infile(input);
      Matrix<dcomp, 9, 9> rt;
      rt.fill(0.);
      Vector3d a_1, a_2, a_3;
      a_1 << 1, 1, 1;
      a_2 << -1, 1, 1;
      a_3 << -1, -1, 1;
      Vector3d A;
  
      dcomp i;
      i = -1;
      i = sqrt(i);
      string line;
      double a, b, c, d, e, f, g;
      double eV_Ry = 0.073498618;
      /* double eV_Ry = 1; */
      /* double iron_fermi = 12.63; // obtained from Wannier90 scf.out */
      /* iron_fermi *= eV_Ry; */
      /* double silver_fermi = 0.4635; //this is lazy - amend to carry this through */
      /* double cshift = silver_fermi - iron_fermi; */
      complex<double> im;
      while (!infile.eof()) 
      {
	getline(infile, line);
	istringstream iss(line);
	iss >> a >> b >> c >> d >> e >> f >> g;
	A = a*a_1 + b*a_2 + c*a_3;
	if ((A(0) == dvec(0)) && (A(1) == dvec(1)) && (A(2) == dvec(2)))
		rt(d-1,e-1) = (f + g*i)*eV_Ry;
      }
      return rt;
}

int main(){

	cout<<"Name the data file\n";
	string Mydata;
	getline(cin, Mydata);
	ofstream Myfile;
	Mydata += ".txt";
	Myfile.open( Mydata.c_str(),ios::trunc );

	cout<<endl;

	int ispin = -1;

	Vector3d d_0, d_1, d_2, d_3, d_4, d_5, d_6, d_7, d_8;
	Vector3d d_13, d_14, d_15, d_16;
	Vector3d d_17, d_18;
	Vector3d d_19, d_20, d_21, d_22, d_23, d_24, d_25, d_26,
		 d_27, d_28, d_29, d_30;
	double a = 1.;
	
	//position of onsite potentials
	d_0 << 0, 0, 0;

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

	//position vectors of third nearest neighbours
	d_19 << 2*a, 2*a, 0;
	d_20 << -2*a, -2*a, 0;
	d_21 << 2*a, 0, 2*a;
	d_22 << -2*a, 0, -2*a;
	d_23 << 0, 2*a, 2*a;
	d_24 << 0, -2*a, -2*a;
	d_25 << 2*a, -2*a, 0;
	d_26 << -2*a, 2*a, 0;
	d_27 << -2*a, 0, 2*a;
	d_28 << 2*a, 0, -2*a;
	d_29 << 0, -2*a, 2*a;
	d_30 << 0, 2*a, -2*a;

	Matrix<dcomp, 9, 9> u;
	Matrix<dcomp, 9, 9> t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8;
	Matrix<dcomp, 9, 9> t_13, t_14, t_15, t_16;
	Matrix<dcomp, 9, 9> t_17, t_18;
	Matrix<dcomp, 9, 9> t_19, t_20, t_21, t_22, t_23, t_24, t_25,
		t_26, t_27, t_28, t_29, t_30;
	u = read(d_0, ispin);

	SelfAdjointEigenSolver<Matrix<dcomp, 9, 9>> uu;
	uu.compute(u);
	Matrix<double, 9, 9> lambdareal;
	Matrix<dcomp, 9, 9> Oo, Odagg, lambda;
	lambdareal = uu.eigenvalues().asDiagonal();
	lambda = lambdareal.cast<dcomp>();
	Oo = uu.eigenvectors();
	Odagg = Oo.adjoint();
	lambda = convert(lambda);

	dcomp i;
	i = -1.;
	i = sqrt(i);
	double s, p, d1, d2, sss1, sss2, pps1, pps2, ppp1, ppp2, dds1, dds2, ddp1, ddp2, ddd1, ddd2, sps1, sps2, sds1, sds2, pds1, pds2, pdp1, pdp2;
	double sss3, sps3, pps3, ppp3, sds3, pds3, pdp3, dds3, ddp3, ddd3;

	sds1 = -0.07158; dds1 = -0.04897; ddp1 = 0.02434; ddd1 = -0.00178;
	pps1 = 0.26892; ppp1 = -0.01859; sps1 = 0.16918;
	pds1 = -0.11882; pdp1 = 0.03462; 
	pdp2 = -0.01088; pds2 = -0.05257;
	sds2 = -0.02805; dds2 = -0.02267; ddp2 = -0.00468; ddd2 = 0.00209;
	ppp2 = 0.03060;
	pps2 = 0.16341; 
	sps2 = 0.06189; 
	sss1 = -0.118047;
	sss2 = -0.0227164;
	s = 1.33239;
	p = 1.94576;
	d1 = 0.846975;
	d2 = 0.79515;

//this block from Papa
      s =  1.14481; // on-site
      p =  1.80769;
      d1 =  0.5*(0.78456 + 0.75661);
      d2 =  0.5*(0.78456 + 0.75661);
      sss1 = -0.13243;   //  same atom hopping
      sps1 =  0.17278;
      pps1 =  0.25911;
      ppp1 =  0.02653;
      sds1 = -0.07145;
      pds1 = -0.09702;
      pdp1 =  0.02129;
      dds1 = -0.05266;
      ddp1 =  0.03276;
      ddd1 = -0.00286;
      sss2 = -0.03003;
      sps2 =  0.07159;
      pps2 =  0.18256;
      ppp2 =  0.03703;
      sds2 = -0.04075;
      pds2 = -0.06522;
      pdp2 = -0.00467;
      dds2 = -0.03396;
      ddp2 =  0.00581;
      ddd2 =  0.00114;
      sss3 =  0.01589;
      sps3 = -0.02306;
      pps3 = -0.04253;
      ppp3 =  0.01538;
      sds3 =  0.00016;
      pds3 =  0.00222;
      pdp3 = -0.00351;
      dds3 =  0.00233;
      ddp3 =  0.00013;
      ddd3 = -0.00060;

	VectorXd nn(10),nnn(20);
	nn<<sss1, sps1, pps1, ppp1, sds1, pds1, pdp1, dds1, ddp1, ddd1;
	nnn<<sss2, sps2, pps2, ppp2, sds2, pds2, pdp2, dds2, ddp2, ddd2, sss3, sps3, pps3, ppp3, sds3, pds3, pdp3, dds3, ddp3, ddd3;

	/* t_1 = read(d_1, ispin); */
	/* t_1 = Odagg*t_1*Oo; */
	/* t_1 = convert(t_1); */
	/* t_2 = read(d_2, ispin); */
	/* t_2 = Odagg*t_2*Oo; */
	/* t_2 = convert(t_2); */
	/* t_3 = read(d_3, ispin); */
	/* t_3 = Odagg*t_3*Oo; */
	/* t_3 = convert(t_3); */
	/* t_4 = read(d_4, ispin); */
	/* t_4 = Odagg*t_4*Oo; */
	/* t_4 = convert(t_4); */
	/* t_5 = read(d_5, ispin); */
	/* t_5 = Odagg*t_5*Oo; */
	/* t_5 = convert(t_5); */
	/* t_6 = read(d_6, ispin); */
	/* t_6 = Odagg*t_6*Oo; */
	/* t_6 = convert(t_6); */
	/* t_7 = read(d_7, ispin); */
	/* t_7 = Odagg*t_7*Oo; */
	/* t_7 = convert(t_7); */
	/* t_8 = read(d_8, ispin); */
	/* t_8 = Odagg*t_8*Oo; */
	/* t_8 = convert(t_8); */
	/* t_13 = read(d_13, ispin); */
	/* t_13 = Odagg*t_13*Oo; */
	/* t_13 = convert(t_13); */
	/* t_14 = read(d_14, ispin); */
	/* t_14 = Odagg*t_14*Oo; */
	/* t_14 = convert(t_14); */
	/* t_15 = read(d_15, ispin); */
	/* t_15 = Odagg*t_15*Oo; */
	/* t_15 = convert(t_15); */
	/* t_16 = read(d_16, ispin); */
	/* t_16 = Odagg*t_16*Oo; */
	/* t_16 = convert(t_16); */
	/* t_17 = read(d_17, ispin); */
	/* t_17 = Odagg*t_17*Oo; */
	/* t_17 = convert(t_17); */
	/* t_18 = read(d_18, ispin); */
	/* t_18 = Odagg*t_18*Oo; */
	/* t_18 = convert(t_18); */

	lambda.fill(0.);
	lambda(0,0) = s;
	lambda(1,1) = lambda(2,2) = lambda(3,3) = p;
	lambda(4,4) = lambda(5,5) = lambda(6,6) = d1;
	lambda(7,7) = lambda(8,8) = d2;
	t_1 = TB(0,1,0,9,d_1,nn,nnn);
	t_2 = TB(0,1,0,9,d_2,nn,nnn);
	t_3 = TB(0,1,0,9,d_3,nn,nnn);
	t_4 = TB(0,1,0,9,d_4,nn,nnn);
	t_5 = TB(0,1,0,9,d_5,nn,nnn);
	t_6 = TB(0,1,0,9,d_6,nn,nnn);
	t_7 = TB(0,1,0,9,d_7,nn,nnn);
	t_8 = TB(0,1,0,9,d_8,nn,nnn);
	t_13 = TB(0,1,1,9,d_13,nn,nnn);
	t_14 = TB(0,1,1,9,d_14,nn,nnn);
	t_15 = TB(0,1,1,9,d_15,nn,nnn);
	t_16 = TB(0,1,1,9,d_16,nn,nnn);
	t_17 = TB(0,1,1,9,d_17,nn,nnn);
	t_18 = TB(0,1,1,9,d_18,nn,nnn);
	t_19 = TB(0,1,2,9,d_19,nn,nnn);
	t_20 = TB(0,1,2,9,d_20,nn,nnn);
	t_21 = TB(0,1,2,9,d_21,nn,nnn);
	t_22 = TB(0,1,2,9,d_22,nn,nnn);
	t_23 = TB(0,1,2,9,d_23,nn,nnn);
	t_24 = TB(0,1,2,9,d_24,nn,nnn);
	t_25 = TB(0,1,2,9,d_25,nn,nnn);
	t_26 = TB(0,1,2,9,d_26,nn,nnn);
	t_27 = TB(0,1,2,9,d_27,nn,nnn);
	t_28 = TB(0,1,2,9,d_28,nn,nnn);
	t_29 = TB(0,1,2,9,d_29,nn,nnn);
	t_30 = TB(0,1,2,9,d_30,nn,nnn);
	d_1 = 0.5*d_1;
	d_2 = 0.5*d_2;
	d_3 = 0.5*d_3;
	d_4 = 0.5*d_4;
	d_5 = 0.5*d_5;
	d_6 = 0.5*d_6;
	d_7 = 0.5*d_7;
	d_8 = 0.5*d_8;
	d_13 = 0.5*d_13;
	d_14 = 0.5*d_14;
	d_15 = 0.5*d_15;
	d_16 = 0.5*d_16;
	d_17 = 0.5*d_17;
	d_18 = 0.5*d_18;
	d_20 = 0.5*d_20;
	d_21 = 0.5*d_21;
	d_22 = 0.5*d_22;
	d_23 = 0.5*d_23;
	d_24 = 0.5*d_24;
	d_25 = 0.5*d_25;
	d_26 = 0.5*d_26;
	d_27 = 0.5*d_27;
	d_28 = 0.5*d_28;
	d_29 = 0.5*d_29;
	d_30 = 0.5*d_30;

	double result;

	double start = 0.4;
	double end = 1.5;
	double step = 0.0026;

	for (double j = start; j<end + step; j=j+step){
		result = kspace(&greens, 1, 0.05, a, j + 1e-4*i, lambda, t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8,
				t_13, t_14, t_15, t_16, t_17, t_18, t_19, t_20, t_21, t_22, t_23, t_24, t_25, 
			       t_26, t_27, t_28, t_29, t_30, d_1, d_2, d_3, d_4,
				d_5, d_6, d_7, d_8, d_13, d_14, d_15, d_16, d_17, d_18, d_19, d_20, d_21, d_22, d_23, d_24, d_25, 
			       d_26, d_27, d_28, d_29, d_30);

		cout<<100*(j-start+step)/(end-start+step)<<"% completed"<<endl;

		Myfile<<j<<" "<<-result*(a)*(a)/(4.*M_PI*M_PI*M_PI)<<endl;
		/* Myfile<<j<<" "<<-0.5*result/M_PI<<endl; */
	}

	Myfile.close();
	return 0;
}
