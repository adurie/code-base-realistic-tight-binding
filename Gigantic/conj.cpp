#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "cgls.h"
#include <vector>
#include "TB_sk.h"
#include <iomanip>
#include <fstream>

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
typedef Vector3d v3;
/* typedef Matrix<double, 9, 1> v4; */
typedef VectorXd v4;
typedef Matrix<dcomp, 9, 9> m3;

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

double xargs2(VectorXd& nnn, v4& W_eigs, v4& W_eigs2, v4& W_eigs3, v4& W_eigs4,
	       	v4& W_eigs5, v4& W_eigs6, v4& W_eigs7, v4& W_eigs8,
	       	const Matrix<dcomp, 9, 9> lambda, const v3& d_13,
	       	const v3& d_14, const v3& d_15, 
		const v3& d_16, const v3& d_17, const v3& d_18){

	int sz = nnn.size();
	VectorXd nn(sz);
	Matrix<dcomp, 9, 9> t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8;
	Matrix<dcomp, 9, 9> t_13, t_14, t_15, t_16;
	Matrix<dcomp, 9, 9> t_17, t_18;
	/* t_1 = TB(0,1,0,9,d_1,nn,nnn); */
	/* t_2 = TB(0,1,0,9,d_2,nn,nnn); */
	/* t_3 = TB(0,1,0,9,d_3,nn,nnn); */
	/* t_4 = TB(0,1,0,9,d_4,nn,nnn); */
	/* t_5 = TB(0,1,0,9,d_5,nn,nnn); */
	/* t_6 = TB(0,1,0,9,d_6,nn,nnn); */
	/* t_7 = TB(0,1,0,9,d_7,nn,nnn); */
	/* t_8 = TB(0,1,0,9,d_8,nn,nnn); */
	t_13 = TB(0,1,1,9,d_13,nn,nnn);
	t_14 = TB(0,1,1,9,d_14,nn,nnn);
	t_15 = TB(0,1,1,9,d_15,nn,nnn);
	t_16 = TB(0,1,1,9,d_16,nn,nnn);
	t_17 = TB(0,1,1,9,d_17,nn,nnn);
	t_18 = TB(0,1,1,9,d_18,nn,nnn);

	dcomp i;
	i = -1.;
	i = sqrt(i);

	Matrix<dcomp, 9, 9> Ew;
	double b = 2.;
	Vector3d K;
	K << 0, 0, 0;
	double fret;
	fret = 0;
	Ew = lambda
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* Matrix4cd AA; */
	/* AA = Ew.topLeftCorner(4,4); */
	ComplexEigenSolver<MatrixXcd> CA(sz);
	CA.compute(Ew);
	VectorXd eigs;
	eigs = CA.eigenvalues().real();
	sort(eigs.data(), eigs.data()+eigs.size());
	for (int k = 0; k < eigs.size(); k++)
		fret += (eigs(k) - W_eigs(k))*(eigs(k) - W_eigs(k));

	K << 3.*M_PI/(2.*b), M_PI/(2.*b), M_PI/(2.*b);
	Ew = lambda
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs2;
	eigs2 = CA.eigenvalues().real();
	sort(eigs2.data(), eigs2.data()+eigs2.size());
	for (int k = 0; k < eigs2.size(); k++)
		fret += (eigs2(k) - W_eigs2(k))*(eigs2(k) - W_eigs2(k));

	K << M_PI/b, 0, 0;
	Ew = lambda
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs3;
	eigs3 = CA.eigenvalues().real();
	sort(eigs3.data(), eigs3.data()+eigs3.size());
	for (int k = 0; k < eigs3.size(); k++)
		fret += (eigs3(k) - W_eigs3(k))*(eigs3(k) - W_eigs3(k));

	K << 2.*M_PI/b, 0, 0;
	Ew = lambda
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs4;
	eigs4 = CA.eigenvalues().real();
	sort(eigs4.data(), eigs4.data()+eigs4.size());
	for (int k = 0; k < eigs4.size(); k++)
		fret += (eigs4(k) - W_eigs4(k))*(eigs4(k) - W_eigs4(k));

	K << M_PI/b, M_PI/b, M_PI/b;
	Ew = lambda
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs5;
	eigs5 = CA.eigenvalues().real();
	sort(eigs5.data(), eigs5.data()+eigs5.size());
	for (int k = 0; k < eigs5.size(); k++)
		fret += (eigs5(k) - W_eigs5(k))*(eigs5(k) - W_eigs5(k));

	K << 0, M_PI/b, M_PI/b;
	Ew = lambda
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs6;
	eigs6 = CA.eigenvalues().real();
	sort(eigs6.data(), eigs6.data()+eigs6.size());
	for (int k = 0; k < eigs6.size(); k++)
		fret += (eigs6(k) - W_eigs6(k))*(eigs6(k) - W_eigs6(k));

	K << M_PI/(2.*b), M_PI/b, M_PI/b;
	Ew = lambda
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs7;
	eigs7 = CA.eigenvalues().real();
	sort(eigs7.data(), eigs7.data()+eigs7.size());
	for (int k = 0; k < eigs7.size(); k++)
		fret += (eigs7(k) - W_eigs7(k))*(eigs7(k) - W_eigs7(k));

	K << 0, M_PI/(2.*b), M_PI/(2.*b);
	Ew = lambda
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs8;
	eigs8 = CA.eigenvalues().real();
	sort(eigs8.data(), eigs8.data()+eigs8.size());
	for (int k = 0; k < eigs8.size(); k++)
		fret += (eigs8(k) - W_eigs8(k))*(eigs8(k) - W_eigs8(k));

	return fret;

}
VectorXd Eigens(double x, double y, double z, const Matrix<dcomp, 9, 9> lambda, const v3& d_1,
	       	const v3& d_2, const v3& d_3, const v3& d_4, const v3& d_5, 
		const v3& d_6, const v3& d_7, const v3& d_8,
		const v3& d_13, const v3& d_14, const v3& d_15, const v3& d_16,
		const v3& d_17, const v3& d_18, const v3& d_19, const v3& d_20,
	       	const v3& d_21, const v3& d_22, const v3& d_23, const v3& d_24, 
		const v3& d_25,	const v3& d_26, const v3& d_27, const v3& d_28,
		const v3& d_29, const v3& d_30, const m3& t_1,
	       	const m3& t_2, const m3& t_3, const m3& t_4, const m3& t_5, 
		const m3& t_6, const m3& t_7, const m3& t_8,
		const m3& t_13, const m3& t_14, const m3& t_15, const m3& t_16,
		const m3& t_17, const m3& t_18, const m3& t_19, const m3& t_20,
	       	const m3& t_21, const m3& t_22, const m3& t_23, const m3& t_24, 
		const m3& t_25, const m3& t_26, const m3& t_27, const m3& t_28,
		const m3& t_29, const m3& t_30){
	Matrix<dcomp, 9, 9> Ew;
	dcomp i;
	i = -1.;
	i = sqrt(i);
	int sz = sqrt(t_1.size());

	/* double b = 2.; */
	Vector3d K;
	/* K << x/b, y/b, z/b; */
	K << x, y, z;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	ComplexEigenSolver<MatrixXcd> CA(sz);
	CA.compute(Ew);
	VectorXd eigs;
	eigs = CA.eigenvalues().real();
	sort(eigs.data(), eigs.data()+eigs.size());
	return eigs;
}

double xargs6(VectorXd& nnnn, const Matrix<dcomp, 9, 9> lambda, const v3& d_1,
	       	const v3& d_2, const v3& d_3, const v3& d_4, const v3& d_5, 
		const v3& d_6, const v3& d_7, const v3& d_8,
		const v3& d_13, const v3& d_14, const v3& d_15, const v3& d_16,
		const v3& d_17, const v3& d_18, const v3& d_19, const v3& d_20,
	       	const v3& d_21, const v3& d_22, const v3& d_23, const v3& d_24, 
		const v3& d_25,	const v3& d_26, const v3& d_27, const v3& d_28,
		const v3& d_29, const v3& d_30, const m3& wt_1,
	       	const m3& wt_2, const m3& wt_3, const m3& wt_4, const m3& wt_5, 
		const m3& wt_6, const m3& wt_7, const m3& wt_8,
		const m3& wt_13, const m3& wt_14, const m3& wt_15, const m3& wt_16,
		const m3& wt_17, const m3& wt_18, const m3& wt_19, const m3& wt_20,
	       	const m3& wt_21, const m3& wt_22, const m3& wt_23, const m3& wt_24, 
		const m3& wt_25, const m3& wt_26, const m3& wt_27, const m3& wt_28,
		const m3& wt_29, const m3& wt_30){

	int sz = nnnn.size()/3;
	VectorXd nn(sz), nnn(2*sz);
	nn = nnnn.head(sz);
	nnn = nnnn.tail(2*sz);
	Matrix<dcomp, 9, 9> t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8;
	Matrix<dcomp, 9, 9> t_13, t_14, t_15, t_16;
	Matrix<dcomp, 9, 9> t_17, t_18;
	Matrix<dcomp, 9, 9> t_19, t_20, t_21, t_22, t_23, t_24, t_25,
		t_26, t_27, t_28, t_29, t_30;
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
	Matrix<dcomp, 9, 9> u;
	u.fill(0.);
	u(0,0) = real(lambda(0,0));
	u(1,1) = u(2,2) = u(3,3) = real((lambda(1,1) + lambda(2,2) + lambda(3,3))/3.);
	u(4,4) = u(5,5) = u(6,6) = real((lambda(4,4) + lambda(5,5) + lambda(6,6))/3.);
	u(7,7) = u(8,8) = real((lambda(7,7) + lambda(8,8))/2.);

	double fret;
	fret = 0;
	VectorXd eigs, W_eigs;
	double x, y, z;
	/* int in = 6; */ 
	/* int it = 0; */
	double pi;
	double b = 2.;
	for (int k = 0; k < 351; k++)
	{
		if (k < 101){
			pi = 2*M_PI*k/100.;
			z = pi/b;
			y = 0;
			x = 0;
		}
		if ((k > 100) && (k < 201)){
			pi = M_PI*(k-100)/100.;
			x = (2*M_PI-pi)/b;
			y = pi/b;
			z = pi/b;
		}	
		if ((k > 200) && (k < 251)){
			pi = M_PI*(k-200)/50.;
			x = (M_PI - pi)/b;
			y = M_PI/b;
			z = M_PI/b;
		}
		if ((k > 250) && (k < 351)){
			pi = M_PI*(k-250)/100.;
			x = 0;
			y = (M_PI-pi)/b;
			z = (M_PI-pi)/b;
		}

	/* for (int xx = 0; xx < (in+1); xx++){ */
	/* 	for (int yy = 0; yy < xx; yy++){ */
	/* 		for (int zz = 0; zz < yy; zz++){ */
	/* 			x = 2.*M_PI*xx/(in*1.); */
	/* 			y = 2.*M_PI*yy/(in*1.); */
	/* 			z = 2.*M_PI*zz/(in*1.); */
	/* 			it++; */

				eigs = Eigens(x, y, z, u,
					d_1, d_2, d_3, d_4, d_5, d_6, d_7, d_8, 
					d_13, d_14, d_15, d_16, d_17, d_18, d_19, d_20,
					d_21, d_22, d_23, d_24, d_25, d_26, d_27, d_28,
					d_29, d_30,
					t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, 
					t_13, t_14, t_15, t_16, t_17, t_18, t_19, t_20,
					t_21, t_22, t_23, t_24, t_25, t_26, t_27, t_28,
					t_29, t_30);
				W_eigs = Eigens(x, y, z, lambda,
					d_1, d_2, d_3, d_4, d_5, d_6, d_7, d_8, 
					d_13, d_14, d_15, d_16, d_17, d_18, d_19, d_20,
					d_21, d_22, d_23, d_24, d_25, d_26, d_27, d_28,
					d_29, d_30,
					wt_1, wt_2, wt_3, wt_4, wt_5, wt_6, wt_7, wt_8, 
					wt_13, wt_14, wt_15, wt_16, wt_17, wt_18, wt_19, wt_20,
					wt_21, wt_22, wt_23, wt_24, wt_25, wt_26, wt_27, wt_28,
					wt_29, wt_30);
				for (int k = 0; k < eigs.size(); k++)
					fret += (eigs(k) - W_eigs(k))*(eigs(k) - W_eigs(k));
			/* } */
		/* } */
	}
	/* cout<<it<<endl; */

	return fret;

}
double xargs5(VectorXd& nnnn, v4& W_eigs, v4& W_eigs2, v4& W_eigs3, v4& W_eigs4,
	       	v4& W_eigs5, v4& W_eigs6, v4& W_eigs7, v4& W_eigs8,
	       	const Matrix<dcomp, 9, 9> lambda, const v3& d_1,
	       	const v3& d_2, const v3& d_3, const v3& d_4, const v3& d_5, 
		const v3& d_6, const v3& d_7, const v3& d_8,
		const v3& d_13, const v3& d_14, const v3& d_15, const v3& d_16,
		const v3& d_17, const v3& d_18, const v3& d_19, const v3& d_20,
	       	const v3& d_21, const v3& d_22, const v3& d_23, const v3& d_24, 
		const v3& d_25,	const v3& d_26, const v3& d_27, const v3& d_28,
		const v3& d_29, const v3& d_30){

	int sz = nnnn.size()/3;
	VectorXd nn(sz), nnn(2*sz);
	nn = nnnn.head(sz);
	nnn = nnnn.tail(2*sz);
	Matrix<dcomp, 9, 9> t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8;
	Matrix<dcomp, 9, 9> t_13, t_14, t_15, t_16;
	Matrix<dcomp, 9, 9> t_17, t_18;
	Matrix<dcomp, 9, 9> t_19, t_20, t_21, t_22, t_23, t_24, t_25,
		t_26, t_27, t_28, t_29, t_30;
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

	dcomp i;
	i = -1.;
	i = sqrt(i);

	Matrix<dcomp, 9, 9> Ew;
	double b = 2.;
	Vector3d K;
	K << 0, 0, 0;
	double fret;
	fret = 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* Matrix4cd AA; */
	/* AA = Ew.topLeftCorner(4,4); */
	ComplexEigenSolver<MatrixXcd> CA(sz);
	CA.compute(Ew);
	VectorXd eigs;
	eigs = CA.eigenvalues().real();
	sort(eigs.data(), eigs.data()+eigs.size());
	for (int k = 0; k < eigs.size(); k++)
		fret += (eigs(k) - W_eigs(k))*(eigs(k) - W_eigs(k));

	K << 3.*M_PI/(2.*b), M_PI/(2.*b), M_PI/(2.*b);
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs2;
	eigs2 = CA.eigenvalues().real();
	sort(eigs2.data(), eigs2.data()+eigs2.size());
	for (int k = 0; k < eigs2.size(); k++)
		fret += (eigs2(k) - W_eigs2(k))*(eigs2(k) - W_eigs2(k));

	K << M_PI/b, 0, 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs3;
	eigs3 = CA.eigenvalues().real();
	sort(eigs3.data(), eigs3.data()+eigs3.size());
	for (int k = 0; k < eigs3.size(); k++)
		fret += (eigs3(k) - W_eigs3(k))*(eigs3(k) - W_eigs3(k));

	K << 2.*M_PI/b, 0, 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs4;
	eigs4 = CA.eigenvalues().real();
	sort(eigs4.data(), eigs4.data()+eigs4.size());
	for (int k = 0; k < eigs4.size(); k++)
		fret += (eigs4(k) - W_eigs4(k))*(eigs4(k) - W_eigs4(k));

	K << M_PI/b, M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs5;
	eigs5 = CA.eigenvalues().real();
	sort(eigs5.data(), eigs5.data()+eigs5.size());
	for (int k = 0; k < eigs5.size(); k++)
		fret += (eigs5(k) - W_eigs5(k))*(eigs5(k) - W_eigs5(k));

	K << 0, M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs6;
	eigs6 = CA.eigenvalues().real();
	sort(eigs6.data(), eigs6.data()+eigs6.size());
	for (int k = 0; k < eigs6.size(); k++)
		fret += (eigs6(k) - W_eigs6(k))*(eigs6(k) - W_eigs6(k));

	K << M_PI/(2.*b), M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs7;
	eigs7 = CA.eigenvalues().real();
	sort(eigs7.data(), eigs7.data()+eigs7.size());
	for (int k = 0; k < eigs7.size(); k++)
		fret += (eigs7(k) - W_eigs7(k))*(eigs7(k) - W_eigs7(k));

	K << 0, M_PI/(2.*b), M_PI/(2.*b);
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs8;
	eigs8 = CA.eigenvalues().real();
	sort(eigs8.data(), eigs8.data()+eigs8.size());
	for (int k = 0; k < eigs8.size(); k++)
		fret += (eigs8(k) - W_eigs8(k))*(eigs8(k) - W_eigs8(k));

	return fret;

}
double xargs4(VectorXd& NNN, VectorXd& nnnn, v4& W_eigs, v4& W_eigs2, v4& W_eigs3, v4& W_eigs4,
	       	v4& W_eigs5, v4& W_eigs6, v4& W_eigs7, v4& W_eigs8,
	       	const Matrix<dcomp, 9, 9> lambda, const v3& d_1,
	       	const v3& d_2, const v3& d_3, const v3& d_4, const v3& d_5, 
		const v3& d_6, const v3& d_7, const v3& d_8,
		const v3& d_13, const v3& d_14, const v3& d_15, const v3& d_16,
		const v3& d_17, const v3& d_18, const v3& d_19, const v3& d_20,
	       	const v3& d_21, const v3& d_22, const v3& d_23, const v3& d_24, 
		const v3& d_25,	const v3& d_26, const v3& d_27, const v3& d_28,
		const v3& d_29, const v3& d_30){

	int sz = nnnn.size()/2;
	VectorXd nn(sz), nnn(2*sz);
	nn = nnnn.head(sz);
	nnn.head(sz) = nnnn.tail(sz);
	nnn.tail(sz) = NNN;
	Matrix<dcomp, 9, 9> t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8;
	Matrix<dcomp, 9, 9> t_13, t_14, t_15, t_16;
	Matrix<dcomp, 9, 9> t_17, t_18;
	Matrix<dcomp, 9, 9> t_19, t_20, t_21, t_22, t_23, t_24, t_25,
		t_26, t_27, t_28, t_29, t_30;
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

	dcomp i;
	i = -1.;
	i = sqrt(i);

	Matrix<dcomp, 9, 9> Ew;
	double b = 2.;
	Vector3d K;
	K << 0, 0, 0;
	double fret;
	fret = 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* Matrix4cd AA; */
	/* AA = Ew.topLeftCorner(4,4); */
	ComplexEigenSolver<MatrixXcd> CA(sz);
	CA.compute(Ew);
	VectorXd eigs;
	eigs = CA.eigenvalues().real();
	sort(eigs.data(), eigs.data()+eigs.size());
	for (int k = 0; k < eigs.size(); k++)
		fret += (eigs(k) - W_eigs(k))*(eigs(k) - W_eigs(k));

	K << 3.*M_PI/(2.*b), M_PI/(2.*b), M_PI/(2.*b);
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs2;
	eigs2 = CA.eigenvalues().real();
	sort(eigs2.data(), eigs2.data()+eigs2.size());
	for (int k = 0; k < eigs2.size(); k++)
		fret += (eigs2(k) - W_eigs2(k))*(eigs2(k) - W_eigs2(k));

	K << M_PI/b, 0, 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs3;
	eigs3 = CA.eigenvalues().real();
	sort(eigs3.data(), eigs3.data()+eigs3.size());
	for (int k = 0; k < eigs3.size(); k++)
		fret += (eigs3(k) - W_eigs3(k))*(eigs3(k) - W_eigs3(k));

	K << 2.*M_PI/b, 0, 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs4;
	eigs4 = CA.eigenvalues().real();
	sort(eigs4.data(), eigs4.data()+eigs4.size());
	for (int k = 0; k < eigs4.size(); k++)
		fret += (eigs4(k) - W_eigs4(k))*(eigs4(k) - W_eigs4(k));

	K << M_PI/b, M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs5;
	eigs5 = CA.eigenvalues().real();
	sort(eigs5.data(), eigs5.data()+eigs5.size());
	for (int k = 0; k < eigs5.size(); k++)
		fret += (eigs5(k) - W_eigs5(k))*(eigs5(k) - W_eigs5(k));

	K << 0, M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs6;
	eigs6 = CA.eigenvalues().real();
	sort(eigs6.data(), eigs6.data()+eigs6.size());
	for (int k = 0; k < eigs6.size(); k++)
		fret += (eigs6(k) - W_eigs6(k))*(eigs6(k) - W_eigs6(k));

	K << M_PI/(2.*b), M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs7;
	eigs7 = CA.eigenvalues().real();
	sort(eigs7.data(), eigs7.data()+eigs7.size());
	for (int k = 0; k < eigs7.size(); k++)
		fret += (eigs7(k) - W_eigs7(k))*(eigs7(k) - W_eigs7(k));

	K << 0, M_PI/(2.*b), M_PI/(2.*b);
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs8;
	eigs8 = CA.eigenvalues().real();
	sort(eigs8.data(), eigs8.data()+eigs8.size());
	for (int k = 0; k < eigs8.size(); k++)
		fret += (eigs8(k) - W_eigs8(k))*(eigs8(k) - W_eigs8(k));

	return fret;

}

double xargs3(VectorXd& nnnn, v4& W_eigs, v4& W_eigs2, v4& W_eigs3, v4& W_eigs4,
	       	v4& W_eigs5, v4& W_eigs6, v4& W_eigs7, v4& W_eigs8,
	       	const Matrix<dcomp, 9, 9> lambda, const v3& d_1,
	       	const v3& d_2, const v3& d_3, const v3& d_4, const v3& d_5, 
		const v3& d_6, const v3& d_7, const v3& d_8,
		const v3& d_13, const v3& d_14, const v3& d_15, const v3& d_16,
		const v3& d_17, const v3& d_18){

	int sz = nnnn.size()/2;
	VectorXd nn(sz), nnn(sz);
	nn = nnnn.head(sz);
	nnn = nnnn.tail(sz);
	Matrix<dcomp, 9, 9> t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8;
	Matrix<dcomp, 9, 9> t_13, t_14, t_15, t_16;
	Matrix<dcomp, 9, 9> t_17, t_18;
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

	dcomp i;
	i = -1.;
	i = sqrt(i);

	Matrix<dcomp, 9, 9> Ew;
	double b = 2.;
	Vector3d K;
	K << 0, 0, 0;
	double fret;
	fret = 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* Matrix4cd AA; */
	/* AA = Ew.topLeftCorner(4,4); */
	ComplexEigenSolver<MatrixXcd> CA(sz);
	CA.compute(Ew);
	VectorXd eigs;
	eigs = CA.eigenvalues().real();
	sort(eigs.data(), eigs.data()+eigs.size());
	for (int k = 0; k < eigs.size(); k++)
		fret += (eigs(k) - W_eigs(k))*(eigs(k) - W_eigs(k));

	K << 3.*M_PI/(2.*b), M_PI/(2.*b), M_PI/(2.*b);
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs2;
	eigs2 = CA.eigenvalues().real();
	sort(eigs2.data(), eigs2.data()+eigs2.size());
	for (int k = 0; k < eigs2.size(); k++)
		fret += (eigs2(k) - W_eigs2(k))*(eigs2(k) - W_eigs2(k));

	K << M_PI/b, 0, 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs3;
	eigs3 = CA.eigenvalues().real();
	sort(eigs3.data(), eigs3.data()+eigs3.size());
	for (int k = 0; k < eigs3.size(); k++)
		fret += (eigs3(k) - W_eigs3(k))*(eigs3(k) - W_eigs3(k));

	K << 2.*M_PI/b, 0, 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs4;
	eigs4 = CA.eigenvalues().real();
	sort(eigs4.data(), eigs4.data()+eigs4.size());
	for (int k = 0; k < eigs4.size(); k++)
		fret += (eigs4(k) - W_eigs4(k))*(eigs4(k) - W_eigs4(k));

	K << M_PI/b, M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs5;
	eigs5 = CA.eigenvalues().real();
	sort(eigs5.data(), eigs5.data()+eigs5.size());
	for (int k = 0; k < eigs5.size(); k++)
		fret += (eigs5(k) - W_eigs5(k))*(eigs5(k) - W_eigs5(k));

	K << 0, M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs6;
	eigs6 = CA.eigenvalues().real();
	sort(eigs6.data(), eigs6.data()+eigs6.size());
	for (int k = 0; k < eigs6.size(); k++)
		fret += (eigs6(k) - W_eigs6(k))*(eigs6(k) - W_eigs6(k));

	K << M_PI/(2.*b), M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs7;
	eigs7 = CA.eigenvalues().real();
	sort(eigs7.data(), eigs7.data()+eigs7.size());
	for (int k = 0; k < eigs7.size(); k++)
		fret += (eigs7(k) - W_eigs7(k))*(eigs7(k) - W_eigs7(k));

	K << 0, M_PI/(2.*b), M_PI/(2.*b);
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs8;
	eigs8 = CA.eigenvalues().real();
	sort(eigs8.data(), eigs8.data()+eigs8.size());
	for (int k = 0; k < eigs8.size(); k++)
		fret += (eigs8(k) - W_eigs8(k))*(eigs8(k) - W_eigs8(k));

	return fret;

}

double xargs(VectorXd& nn, v4& W_eigs, v4& W_eigs2, v4& W_eigs3, v4& W_eigs4,
	       	v4& W_eigs5, v4& W_eigs6, v4& W_eigs7, v4& W_eigs8,
	       	const Matrix<dcomp, 9, 9> lambda, const v3& d_1,
	       	const v3& d_2, const v3& d_3, const v3& d_4, const v3& d_5, 
		const v3& d_6, const v3& d_7, const v3& d_8){

	int sz = nn.size();
	VectorXd nnn(sz);
	Matrix<dcomp, 9, 9> t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8;
	Matrix<dcomp, 9, 9> t_13, t_14, t_15, t_16;
	Matrix<dcomp, 9, 9> t_17, t_18;
	t_1 = TB(0,1,0,9,d_1,nn,nnn);
	t_2 = TB(0,1,0,9,d_2,nn,nnn);
	t_3 = TB(0,1,0,9,d_3,nn,nnn);
	t_4 = TB(0,1,0,9,d_4,nn,nnn);
	t_5 = TB(0,1,0,9,d_5,nn,nnn);
	t_6 = TB(0,1,0,9,d_6,nn,nnn);
	t_7 = TB(0,1,0,9,d_7,nn,nnn);
	t_8 = TB(0,1,0,9,d_8,nn,nnn);
	/* t_13 = TB(0,1,1,9,d_13,nn,nnn); */
	/* t_14 = TB(0,1,1,9,d_14,nn,nnn); */
	/* t_15 = TB(0,1,1,9,d_15,nn,nnn); */
	/* t_16 = TB(0,1,1,9,d_16,nn,nnn); */
	/* t_17 = TB(0,1,1,9,d_17,nn,nnn); */
	/* t_18 = TB(0,1,1,9,d_18,nn,nnn); */

	dcomp i;
	i = -1.;
	i = sqrt(i);

	Matrix<dcomp, 9, 9> Ew;
	double b = 2.;
	Vector3d K;
	K << 0, 0, 0;
	double fret;
	fret = 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			/* + t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K)) */
			/* + t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K)) */
			/* + t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))); */
		;
	/* Matrix4cd AA; */
	/* AA = Ew.topLeftCorner(4,4); */
	ComplexEigenSolver<MatrixXcd> CA(sz);
	CA.compute(Ew);
	VectorXd eigs;
	eigs = CA.eigenvalues().real();
	sort(eigs.data(), eigs.data()+eigs.size());
	for (int k = 0; k < eigs.size(); k++)
		fret += (eigs(k) - W_eigs(k))*(eigs(k) - W_eigs(k));

	K << 3.*M_PI/(2.*b), M_PI/(2.*b), M_PI/(2.*b);
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
		;
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs2;
	eigs2 = CA.eigenvalues().real();
	sort(eigs2.data(), eigs2.data()+eigs2.size());
	for (int k = 0; k < eigs2.size(); k++)
		fret += (eigs2(k) - W_eigs2(k))*(eigs2(k) - W_eigs2(k));

	K << M_PI/b, 0, 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
		;
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs3;
	eigs3 = CA.eigenvalues().real();
	sort(eigs3.data(), eigs3.data()+eigs3.size());
	for (int k = 0; k < eigs3.size(); k++)
		fret += (eigs3(k) - W_eigs3(k))*(eigs3(k) - W_eigs3(k));

	K << 2.*M_PI/b, 0, 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
		;
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs4;
	eigs4 = CA.eigenvalues().real();
	sort(eigs4.data(), eigs4.data()+eigs4.size());
	for (int k = 0; k < eigs4.size(); k++)
		fret += (eigs4(k) - W_eigs4(k))*(eigs4(k) - W_eigs4(k));

	K << M_PI/b, M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
		;
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs5;
	eigs5 = CA.eigenvalues().real();
	sort(eigs5.data(), eigs5.data()+eigs5.size());
	for (int k = 0; k < eigs5.size(); k++)
		fret += (eigs5(k) - W_eigs5(k))*(eigs5(k) - W_eigs5(k));

	K << 0, M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
		;
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs6;
	eigs6 = CA.eigenvalues().real();
	sort(eigs6.data(), eigs6.data()+eigs6.size());
	for (int k = 0; k < eigs6.size(); k++)
		fret += (eigs6(k) - W_eigs6(k))*(eigs6(k) - W_eigs6(k));

	K << M_PI/(2.*b), M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
		;
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs7;
	eigs7 = CA.eigenvalues().real();
	sort(eigs7.data(), eigs7.data()+eigs7.size());
	for (int k = 0; k < eigs7.size(); k++)
		fret += (eigs7(k) - W_eigs7(k))*(eigs7(k) - W_eigs7(k));

	K << 0, M_PI/(2.*b), M_PI/(2.*b);
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
		;
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs8;
	eigs8 = CA.eigenvalues().real();
	sort(eigs8.data(), eigs8.data()+eigs8.size());
	for (int k = 0; k < eigs8.size(); k++)
		fret += (eigs8(k) - W_eigs8(k))*(eigs8(k) - W_eigs8(k));

	return fret;

}

template <typename func1, typename... Args>
void dxargs(VectorXd& p, VectorXd& xi, func1&& func, Args&&... params){
	double delta = 1e-4;
	int n;
	n = p.size();
	VectorXd pd(n), pn(n);
	pd = p;
	pn = p;
	for (int k = 0; k<p.size(); k++){
		pd(k) = p(k) + delta;
		pn(k) = p(k) - delta;
		xi(k) = (forward<func1>(func)(pd, forward<Args>(params)...) 
				- forward<func1>(func)(pn, forward<Args>(params)...))/(2.*delta); 
	}
}

int main(){
	int ispin = +1;

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
	cout<<endl;

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
	t_1 = read(d_1, ispin);
	t_1 = Odagg*t_1*Oo;
	t_1 = convert(t_1);
	t_2 = read(d_2, ispin);
	t_2 = Odagg*t_2*Oo;
	t_2 = convert(t_2);
	t_3 = read(d_3, ispin);
	t_3 = Odagg*t_3*Oo;
	t_3 = convert(t_3);
	t_4 = read(d_4, ispin);
	t_4 = Odagg*t_4*Oo;
	t_4 = convert(t_4);
	t_5 = read(d_5, ispin);
	t_5 = Odagg*t_5*Oo;
	t_5 = convert(t_5);
	t_6 = read(d_6, ispin);
	t_6 = Odagg*t_6*Oo;
	t_6 = convert(t_6);
	t_7 = read(d_7, ispin);
	t_7 = Odagg*t_7*Oo;
	t_7 = convert(t_7);
	t_8 = read(d_8, ispin);
	t_8 = Odagg*t_8*Oo;
	t_8 = convert(t_8);
	t_13 = read(d_13, ispin);
	t_13 = Odagg*t_13*Oo;
	t_13 = convert(t_13);
	t_14 = read(d_14, ispin);
	t_14 = Odagg*t_14*Oo;
	t_14 = convert(t_14);
	t_15 = read(d_15, ispin);
	t_15 = Odagg*t_15*Oo;
	t_15 = convert(t_15);
	t_16 = read(d_16, ispin);
	t_16 = Odagg*t_16*Oo;
	t_16 = convert(t_16);
	t_17 = read(d_17, ispin);
	t_17 = Odagg*t_17*Oo;
	t_17 = convert(t_17);
	t_18 = read(d_18, ispin);
	t_18 = Odagg*t_18*Oo;
	t_18 = convert(t_18);
	t_19 = read(d_19, ispin);
	t_19 = Odagg*t_19*Oo;
	t_19 = convert(t_19);
	t_20 = read(d_20, ispin);
	t_20 = Odagg*t_20*Oo;
	t_20 = convert(t_20);
	t_21 = read(d_21, ispin);
	t_21 = Odagg*t_21*Oo;
	t_21 = convert(t_21);
	t_22 = read(d_22, ispin);
	t_22 = Odagg*t_22*Oo;
	t_22 = convert(t_22);
	t_23 = read(d_23, ispin);
	t_23 = Odagg*t_23*Oo;
	t_23 = convert(t_23);
	t_24 = read(d_24, ispin);
	t_24 = Odagg*t_24*Oo;
	t_24 = convert(t_24);
	t_25 = read(d_25, ispin);
	t_25 = Odagg*t_25*Oo;
	t_25 = convert(t_25);
	t_26 = read(d_26, ispin);
	t_26 = Odagg*t_26*Oo;
	t_26 = convert(t_26);
	t_27 = read(d_27, ispin);
	t_27 = Odagg*t_27*Oo;
	t_27 = convert(t_27);
	t_28 = read(d_28, ispin);
	t_28 = Odagg*t_28*Oo;
	t_28 = convert(t_28);
	t_29 = read(d_29, ispin);
	t_29 = Odagg*t_29*Oo;
	t_29 = convert(t_29);
	t_30 = read(d_30, ispin);
	t_30 = Odagg*t_30*Oo;
	t_30 = convert(t_30);

	dcomp i;
	i = -1.;
	i = sqrt(i);

	Matrix<dcomp, 9, 9> Ew;
	double b = 2.;
	Vector3d K;
	K << 0, 0, 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			/* + t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K)) */
			/* + t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K)) */
			/* + t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))); */
		;
	/* Matrix4cd AA; */
	/* AA = Ew.topLeftCorner(4,4); */
	ComplexEigenSolver<MatrixXcd> CA(9);
	CA.compute(Ew);
	VectorXd eigs;
	eigs = CA.eigenvalues().real();
	sort(eigs.data(), eigs.data()+eigs.size());

	K << 3.*M_PI/(2.*b), M_PI/(2.*b), M_PI/(2.*b);
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
		;
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs2;
	eigs2 = CA.eigenvalues().real();
	sort(eigs2.data(), eigs2.data()+eigs2.size());

	K << M_PI/b, 0, 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
		;
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs3;
	eigs3 = CA.eigenvalues().real();
	sort(eigs3.data(), eigs3.data()+eigs3.size());

	K << 2.*M_PI/b, 0, 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
		;
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs4;
	eigs4 = CA.eigenvalues().real();
	sort(eigs4.data(), eigs4.data()+eigs4.size());

	K << M_PI/b, M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
		;
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs5;
	eigs5 = CA.eigenvalues().real();
	sort(eigs5.data(), eigs5.data()+eigs5.size());

	K << 0, M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
		;
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs6;
	eigs6 = CA.eigenvalues().real();
	sort(eigs6.data(), eigs6.data()+eigs6.size());

	K << M_PI/(2.*b), M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
		;
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs7;
	eigs7 = CA.eigenvalues().real();
	sort(eigs7.data(), eigs7.data()+eigs7.size());

	K << 0, M_PI/(2.*b), M_PI/(2.*b);
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
		;
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	VectorXd eigs8;
	eigs8 = CA.eigenvalues().real();
	sort(eigs8.data(), eigs8.data()+eigs8.size());

	double s, p, d1, d2, sss1, sss2, pps1, pps2, ppp1, ppp2, dds1, dds2, ddp1, ddp2, ddd1, ddd2, sps1, sps2, sds1, sds2, pds1, pds2, pdp1, pdp2;
	//intial guess
	sss2 = sps2 = pps2 = ppp2 = sds2 = pds2 = pdp2 = dds2 = ddp2 = ddd2 = 0.;
	sss1 = -0.118047;
	sps1 = 0.16918;
	pps1 = 0.26892; 
	ppp1 = -0.01859; 
	sds1 = -0.07158; dds1 = -0.04897; ddp1 = 0.02434; ddd1 = -0.00178;
	pds1 = -0.11882; pdp1 = 0.03462; 
	pdp2 = -0.01088; pds2 = -0.05257;
	sds2 = -0.02805; dds2 = -0.02267; ddp2 = -0.00468; ddd2 = 0.00209;
	ppp2 = 0.03060;
	pps2 = 0.16341; 
	sps2 = 0.06189; 
	VectorXd nn(10),nnn(10);
	nn<<sss1, sps1, pps1, ppp1, sds1, pds1, pdp1, dds1, ddp1, ddd1;
	nnn<<sss2, sps2, pps2, ppp2, sds2, pds2, pdp2, dds2, ddp2, ddd2;

	double ftol = 1e-5;
	double fret;
	int iter;
	cout<<"Find the SK potentials of the 1st neighbour terms:"<<endl;
	frprmn(nn, ftol, iter, fret, xargs, eigs, eigs2, eigs3,
		       	eigs4, eigs5, eigs6, eigs7, eigs8, lambda,
		d_1, d_2, d_3, d_4, d_5, d_6, d_7, d_8);
	cout<<"sss1 = "<<nn(0)<<"; sps1 = "<<nn(1)<<"; pps1 = "<<nn(2)<<"; ppp1 = "
		<<nn(3)<<"; sds1 = "<<nn(4)<<";"<<endl<<"pds1 = "<<nn(5)<<"; pdp1 = "<<nn(6)<<
		"; dds1 = "<<nn(7)<<"; ddp1 = "<<nn(8)<<"; ddd1 = "<<nn(9)<<";"<<endl;
	cout<<endl;
	cout<<"number of iterations: "<<iter<<endl;
	cout<<"minimum of function: "<<fret<<endl;
	cout<<endl<<endl;
	K << 0, 0, 0;
	Ew = lambda
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* Matrix4cd AA; */
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs = CA.eigenvalues().real();
	sort(eigs.data(), eigs.data()+eigs.size());

	K << 3.*M_PI/(2.*b), M_PI/(2.*b), M_PI/(2.*b);
	Ew = lambda
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs2 = CA.eigenvalues().real();
	sort(eigs2.data(), eigs2.data()+eigs2.size());

	K << M_PI/b, 0, 0;
	Ew = lambda
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs3 = CA.eigenvalues().real();
	sort(eigs3.data(), eigs3.data()+eigs3.size());

	K << 2.*M_PI/b, 0, 0;
	Ew = lambda
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs4 = CA.eigenvalues().real();
	sort(eigs4.data(), eigs4.data()+eigs4.size());

	K << M_PI/b, M_PI/b, M_PI/b;
	Ew = lambda
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs5 = CA.eigenvalues().real();
	sort(eigs5.data(), eigs5.data()+eigs5.size());

	K << 0, M_PI/b, M_PI/b;
	Ew = lambda
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs6 = CA.eigenvalues().real();
	sort(eigs6.data(), eigs6.data()+eigs6.size());

	K << M_PI/(2.*b), M_PI/b, M_PI/b;
	Ew = lambda
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs7 = CA.eigenvalues().real();
	sort(eigs7.data(), eigs7.data()+eigs7.size());

	K << 0, M_PI/(2.*b), M_PI/(2.*b);
	Ew = lambda
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs8 = CA.eigenvalues().real();
	sort(eigs8.data(), eigs8.data()+eigs8.size());
	cout<<"Now find the SK potentials of the 2nd neighbours:"<<endl;
	frprmn(nnn, ftol, iter, fret, xargs2, eigs, eigs2, eigs3,
		       	eigs4, eigs5, eigs6, eigs7, eigs8, lambda,
		d_13, d_14, d_15, d_16, d_17, d_18);
	cout<<"sss2 = "<<nnn(0)<<"; sps2 = "<<nnn(1)<<"; pps2 = "<<nnn(2)<<"; ppp2 = "
		<<nnn(3)<<"; sds2 = "<<nnn(4)<<";"<<endl<<"pds2 = "<<nnn(5)<<"; pdp2 = "<<nnn(6)<<
		"; dds2 = "<<nnn(7)<<"; ddp2 = "<<nnn(8)<<"; ddd2 = "<<nnn(9)<<";"<<endl;
	cout<<endl;
	cout<<"number of iterations: "<<iter<<endl;
	cout<<"minimum of function: "<<fret<<endl;
	cout<<endl<<endl;
	K << 0, 0, 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* Matrix4cd AA; */
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs = CA.eigenvalues().real();
	sort(eigs.data(), eigs.data()+eigs.size());

	K << 3.*M_PI/(2.*b), M_PI/(2.*b), M_PI/(2.*b);
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs2 = CA.eigenvalues().real();
	sort(eigs2.data(), eigs2.data()+eigs2.size());

	K << M_PI/b, 0, 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs3 = CA.eigenvalues().real();
	sort(eigs3.data(), eigs3.data()+eigs3.size());

	K << 2.*M_PI/b, 0, 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs4 = CA.eigenvalues().real();
	sort(eigs4.data(), eigs4.data()+eigs4.size());

	K << M_PI/b, M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs5 = CA.eigenvalues().real();
	sort(eigs5.data(), eigs5.data()+eigs5.size());

	K << 0, M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs6 = CA.eigenvalues().real();
	sort(eigs6.data(), eigs6.data()+eigs6.size());

	K << M_PI/(2.*b), M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs7 = CA.eigenvalues().real();
	sort(eigs7.data(), eigs7.data()+eigs7.size());

	K << 0, M_PI/(2.*b), M_PI/(2.*b);
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs8 = CA.eigenvalues().real();
	sort(eigs8.data(), eigs8.data()+eigs8.size());
	VectorXd nnnn(20);
	nnnn.head(10) = nn;
	nnnn.tail(10) = nnn;
	cout<<"Now find 1st and 2nd neigbour terms in one run:"<<endl;
	frprmn(nnnn, ftol, iter, fret, xargs3, eigs, eigs2, eigs3,
		       	eigs4, eigs5, eigs6, eigs7, eigs8, lambda,
			d_1, d_2, d_3, d_4, d_5, d_6, d_7, d_8, 
		d_13, d_14, d_15, d_16, d_17, d_18);
	cout<<"1st nearest neighbour SK terms:"<<endl;
	cout<<"sss1 = "<<nnnn(0)<<"; sps1 = "<<nnnn(1)<<"; pps1 = "<<nnnn(2)<<"; ppp1 = "
		<<nnnn(3)<<"; sds1 = "<<nnnn(4)<<";"<<endl<<"pds1 = "<<nnnn(5)<<"; pdp1 = "<<nnnn(6)<<
		"; dds1 = "<<nnnn(7)<<"; ddp1 = "<<nnnn(8)<<"; ddd1 = "<<nnnn(9)<<";"<<endl;
	cout<<endl;
	cout<<"2nd nearest neighbour SK terms:"<<endl;
	cout<<"sss2 = "<<nnnn(10)<<"; sps2 = "<<nnnn(11)<<"; pps2 = "<<nnnn(12)<<"; ppp2 = "
		<<nnnn(13)<<"; sds2 = "<<nnnn(14)<<";"<<endl<<"pds2 = "<<nnnn(15)<<"; pdp2 = "<<nnnn(16)<<
		"; dds2 = "<<nnnn(17)<<"; ddp2 = "<<nnnn(18)<<"; ddd2 = "<<nnnn(19)<<";"<<endl;
	cout<<endl;
	cout<<"number of iterations: "<<iter<<endl;
	cout<<"minimum of function: "<<fret<<endl;
	cout<<endl<<endl;
	K << 0, 0, 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* Matrix4cd AA; */
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs = CA.eigenvalues().real();
	sort(eigs.data(), eigs.data()+eigs.size());

	K << 3.*M_PI/(2.*b), M_PI/(2.*b), M_PI/(2.*b);
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs2 = CA.eigenvalues().real();
	sort(eigs2.data(), eigs2.data()+eigs2.size());

	K << M_PI/b, 0, 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs3 = CA.eigenvalues().real();
	sort(eigs3.data(), eigs3.data()+eigs3.size());

	K << 2.*M_PI/b, 0, 0;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs4 = CA.eigenvalues().real();
	sort(eigs4.data(), eigs4.data()+eigs4.size());

	K << M_PI/b, M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs5 = CA.eigenvalues().real();
	sort(eigs5.data(), eigs5.data()+eigs5.size());

	K << 0, M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs6 = CA.eigenvalues().real();
	sort(eigs6.data(), eigs6.data()+eigs6.size());

	K << M_PI/(2.*b), M_PI/b, M_PI/b;
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs7 = CA.eigenvalues().real();
	sort(eigs7.data(), eigs7.data()+eigs7.size());

	K << 0, M_PI/(2.*b), M_PI/(2.*b);
	Ew = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
			+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
			+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
			+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
			+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
			+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K));
	/* AA = Ew.topLeftCorner(4,4); */
	CA.compute(Ew);
	eigs8 = CA.eigenvalues().real();
	sort(eigs8.data(), eigs8.data()+eigs8.size());
	VectorXd NNN(10);
	NNN.fill(0.);
	cout<<"Now find 3rd nearest neighbour SK terms whilst keeping 1st & 2nd neighbour terms fixed:"<<endl;
	frprmn(NNN, ftol, iter, fret, xargs4, nnnn, eigs, eigs2, eigs3,
		       	eigs4, eigs5, eigs6, eigs7, eigs8, lambda,
			d_1, d_2, d_3, d_4, d_5, d_6, d_7, d_8, 
		d_13, d_14, d_15, d_16, d_17, d_18, d_19, d_20,
		d_21, d_22, d_23, d_24, d_25, d_26, d_27, d_28,
		d_29, d_30);
	cout<<"sss3 = "<<NNN(0)<<"; sps3 = "<<NNN(1)<<"; pps3 = "<<NNN(2)<<"; ppp3 = "
		<<NNN(3)<<"; sds3 = "<<NNN(4)<<";"<<endl<<"pds3 = "<<NNN(5)<<"; pdp3 = "<<NNN(6)<<
		"; dds3 = "<<NNN(7)<<"; ddp3 = "<<NNN(8)<<"; ddd3 = "<<NNN(9)<<";"<<endl;
	cout<<endl;
	cout<<"number of iterations: "<<iter<<endl;
	cout<<"minimum of function: "<<fret<<endl;
	cout<<endl<<endl;
	VectorXd NNNN(30);
	NNNN.head(20) = nnnn;
	NNNN.tail(10) = NNN;
	cout<<"Now, find up to 3rd neighbour terms in one run:"<<endl;
	frprmn(NNNN, ftol, iter, fret, xargs5, eigs, eigs2, eigs3,
		       	eigs4, eigs5, eigs6, eigs7, eigs8, lambda,
			d_1, d_2, d_3, d_4, d_5, d_6, d_7, d_8, 
		d_13, d_14, d_15, d_16, d_17, d_18, d_19, d_20,
		d_21, d_22, d_23, d_24, d_25, d_26, d_27, d_28,
		d_29, d_30);
	cout<<"1st nearest neighbour SK terms:"<<endl;
	cout<<"sss1 = "<<NNNN(0)<<"; sps1 = "<<NNNN(1)<<"; pps1 = "<<NNNN(2)<<"; ppp1 = "
		<<NNNN(3)<<"; sds1 = "<<NNNN(4)<<";"<<endl<<"pds1 = "<<NNNN(5)<<"; pdp1 = "<<NNNN(6)<<
		"; dds1 = "<<NNNN(7)<<"; ddp1 = "<<NNNN(8)<<"; ddd1 = "<<NNNN(9)<<";"<<endl;
	cout<<endl;
	cout<<"2nd nearest neighbour SK terms:"<<endl;
	cout<<"sss2 = "<<NNNN(10)<<"; sps2 = "<<NNNN(11)<<"; pps2 = "<<NNNN(12)<<"; ppp2 = "
		<<NNNN(13)<<"; sds2 = "<<NNNN(14)<<";"<<endl<<"pds2 = "<<NNNN(15)<<"; pdp2 = "<<NNNN(16)<<
		"; dds2 = "<<NNNN(17)<<"; ddp2 = "<<NNNN(18)<<"; ddd2 = "<<NNNN(19)<<";"<<endl;
	cout<<endl;
	cout<<"3rd nearest neighbour SK terms:"<<endl;
	cout<<"sss3 = "<<NNNN(20)<<"; sps3 = "<<NNNN(21)<<"; pps3 = "<<NNNN(22)<<"; ppp3 = "
		<<NNNN(23)<<"; sds3 = "<<NNNN(24)<<"; "<<endl<<"pds3 = "<<NNNN(25)<<"; pdp3 = "<<NNNN(26)<<
		"; dds3 = "<<NNNN(27)<<"; ddp3 = "<<NNNN(28)<<"; ddd3 = "<<NNNN(29)<<";"<<endl;
	cout<<endl;
	cout<<"number of iterations: "<<iter<<endl;
	cout<<"minimum of function: "<<fret<<endl;
	cout<<endl<<endl;
	cout<<"Finally, find up to 3rd neighbour terms in one run, using loads of k-points:"<<endl;
	frprmn(NNNN, ftol, iter, fret, xargs6, lambda,
			d_1, d_2, d_3, d_4, d_5, d_6, d_7, d_8, 
		d_13, d_14, d_15, d_16, d_17, d_18, d_19, d_20,
		d_21, d_22, d_23, d_24, d_25, d_26, d_27, d_28,
		d_29, d_30,
			t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, 
		t_13, t_14, t_15, t_16, t_17, t_18, t_19, t_20,
		t_21, t_22, t_23, t_24, t_25, t_26, t_27, t_28,
		t_29, t_30);
	cout<<"1st nearest neighbour SK terms:"<<endl;
	cout<<"sss1 = "<<NNNN(0)<<"; sps1 = "<<NNNN(1)<<"; pps1 = "<<NNNN(2)<<"; ppp1 = "
		<<NNNN(3)<<"; sds1 = "<<NNNN(4)<<";"<<endl<<"pds1 = "<<NNNN(5)<<"; pdp1 = "<<NNNN(6)<<
		"; dds1 = "<<NNNN(7)<<"; ddp1 = "<<NNNN(8)<<"; ddd1 = "<<NNNN(9)<<";"<<endl;
	cout<<endl;
	cout<<"2nd nearest neighbour SK terms:"<<endl;
	cout<<"sss2 = "<<NNNN(10)<<"; sps2 = "<<NNNN(11)<<"; pps2 = "<<NNNN(12)<<"; ppp2 = "
		<<NNNN(13)<<"; sds2 = "<<NNNN(14)<<";"<<endl<<"pds2 = "<<NNNN(15)<<"; pdp2 = "<<NNNN(16)<<
		"; dds2 = "<<NNNN(17)<<"; ddp2 = "<<NNNN(18)<<"; ddd2 = "<<NNNN(19)<<";"<<endl;
	cout<<endl;
	cout<<"3rd nearest neighbour SK terms:"<<endl;
	cout<<"sss3 = "<<NNNN(20)<<"; sps3 = "<<NNNN(21)<<"; pps3 = "<<NNNN(22)<<"; ppp3 = "
		<<NNNN(23)<<"; sds3 = "<<NNNN(24)<<"; "<<endl<<"pds3 = "<<NNNN(25)<<"; pdp3 = "<<NNNN(26)<<
		"; dds3 = "<<NNNN(27)<<"; ddp3 = "<<NNNN(28)<<"; ddd3 = "<<NNNN(29)<<";"<<endl;
	cout<<endl;
	cout<<"number of iterations: "<<iter<<endl;
	cout<<"minimum of function: "<<fret<<endl;
	return 0;
}
