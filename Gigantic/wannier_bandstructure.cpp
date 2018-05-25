#include <iostream>
#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <vector>
#include "TB_sk.h"
/* #include "TBdynamic.h" */

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
//calculates the bandstructure of fcc Cu
Matrix<double,10,1> sk_extraction(const Matrix<double,9,9> &t, const Vector3d &pos){
	Matrix<double,81,10> b;
	b.fill(0);

	  double x, y, z;
	  Vector3d X, Y, Z;
	  X << 1, 0, 0;
	  Y << 0, 1, 0;
	  Z << 0, 0, 1;
	  x = pos.dot(X)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	  y = pos.dot(Y)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	  z = pos.dot(Z)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 

      double xx=x*x;
      double xy=x*y;
      double yy=y*y;
      double yz=y*z;
      double zz=z*z;
      double zx=z*x;
      double xxyy=xx*yy;
      double yyzz=yy*zz;
      double zzxx=zz*xx;
      double r3=sqrt(3.);
      double f8=3.*zz-1.;
      double f1=xx+yy;
      double f2=xx-yy;
      double f3=zz-.5*f1;
	/*the following rh indices are;
	 * 0 = sss
	 * 1 = sps
	 * 2 = pps
	 * 3 = ppp
	 * 4 = ss
	 * 5 = ps
	 * 6 = pp
	 * 7 = ds
	 * 8 = dp
	 * 9 = dd */
	b(0,0) = 1;
	b(1,1) = x;
	b(2,1) = y;
	b(3,1) = z;
	b(4,4) = xy*r3;
	b(5,4) = yz*r3;
	b(6,4) = zx*r3;
	b(7,4) = 0.5*f2*r3;
	b(8,4) = 0.5*f8;
	b(9,1) = -x;
	b(10,2) = xx;
	b(10,3) = 1.-xx;
	b(11,2) = xy;
	b(11,3) = -xy;
	b(12,2) = zx;
	b(12,3) = -zx;
	b(13,5) = y*r3*xx;
	b(13,6) = y*(1.-2.*xx); 
	b(14,5) = r3*xy*z;
	b(14,6) = -2.*xy*z;
	b(15,5) = z*r3*xx;
	b(15,6) = z*(1.-2.*xx);
	b(16,5) = x*.5*r3*f2;
	b(16,6) = x*(1.-f2);
	b(17,5) = x*.5*f8;
	b(17,6) = -r3*zz*x;
	b(18,1) = -y;
	b(19,2) = xy;
	b(19,3) = -xy;
	b(20,2) = yy;
	b(20,3) = (1.-yy);
	b(21,2) = yz;
	b(21,3) = -yz;
	b(22,5) = r3*yy*x;
	b(22,6) = (1.-2.*yy)*x;
	b(23,5) = r3*yy*z;
	b(23,6) = (1.-2.*yy)*z;
	b(24,5) = r3*xy*z;
	b(24,6) = -2.*xy*z;
	b(25,5) = .5*r3*f2*y;
	b(25,6) = -(1.+f2)*y;
	b(26,5) = .5*f8*y;
	b(26,6) = -r3*zz*y;
	b(27,1) = -z;
	b(28,2) = zx;
	b(28,3) = -zx;
	b(29,2) = yz;
	b(29,3) = -yz;
	b(30,2) = zz;
	b(30,3) = (1.-zz);
	b(31,5) = r3*xy*z;
	b(31,6) = -2.*xy*z;
	b(32,5) = r3*zz*y;
	b(32,6) = (1.-2.*zz)*y;
	b(33,5) = r3*zz*x;
	b(33,6) = (1.-2.*zz)*x;
	b(34,5) = .5*r3*f2*z;
	b(34,6) = -f2*z;
	b(35,5) = .5*f8*z;
	b(35,6) = r3*f1*z;
	b(36,4) = xy*r3;
	b(37,5) = -y*r3*xx;
	b(37,6) = -y*(1.-2.*xx); 
	b(38,5) = -r3*yy*x;
	b(38,6) = -(1.-2.*yy)*x;
	b(39,5) = -r3*xy*z;
	b(39,6) = 2.*xy*z;
	b(40,7) = 3.*xxyy;
	b(40,8) = (f1-4.*xxyy);
	b(40,9) = (zz+xxyy);
	b(41,7) = 3.*yy*zx;
	b(41,8) = (1.-4.*yy)*zx;
	b(41,9) = (yy-1.)*zx;
	b(42,7) = 3.*xx*yz;
	b(42,8) = (1.-4.*xx)*yz;
	b(42,9) = (xx-1.)*yz;
	b(43,7) = 1.5*f2*xy;
	b(43,8) = -2.*f2*xy;
	b(43,9) = .5*f2*xy;
	b(44,7) = r3*f3*xy;
	b(44,8) = -2.*r3*zz*xy;
	b(44,9) = .5*r3*(1.+zz)*xy;
	b(45,4) = yz*r3;
	b(46,5) = -r3*xy*z;
	b(46,6) = 2.*xy*z;
	b(47,5) = -r3*yy*z;
	b(47,6) = -(1.-2.*yy)*z;
	b(48,5) = -r3*zz*y;
	b(48,6) = -(1.-2.*zz)*y;
	b(49,7) = 3.*yy*zx;
	b(49,8) = (1.-4.*yy)*zx;
	b(49,9) = (yy-1.)*zx;
	b(50,7) = 3.*yyzz;
	b(50,8) = (yy+zz-4.*yyzz);
	b(50,9) = (xx+yyzz);
	b(51,7) = 3.*zz*xy;
	b(51,8) = (1.-4.*zz)*xy;
	b(51,9) = (zz-1.)*xy;
	b(52,7) = 1.5*f2*yz;
	b(52,8) = -(1.+2.*f2)*yz;
	b(52,9) = (1.+.5*f2)*yz;
	b(53,7) = r3*f3*yz;
	b(53,8) = r3*(f1-zz)*yz;
	b(53,9) = -.5*r3*f1*yz;
	b(54,4) = zx*r3;
	b(55,5) = -z*r3*xx;
	b(55,6) = -z*(1.-2.*xx);
	b(56,5) = -r3*xy*z;
	b(56,6) = 2.*xy*z;
	b(57,5) = -r3*zz*x;
	b(57,6) = -(1.-2.*zz)*x;
	b(58,7) = 3.*xx*yz;
	b(58,8) = (1.-4.*xx)*yz;
	b(58,9) = (xx-1.)*yz;
	b(59,7) = 3.*zz*xy;
	b(59,8) = (1.-4.*zz)*xy;
	b(59,9) = (zz-1.)*xy;
	b(60,7) = 3.*zzxx;
	b(60,8) = (zz+xx-4.*zzxx);
	b(60,9) = (yy+zzxx);
	b(61,7) = 1.5*f2*zx;
	b(61,8) = (1.-2.*f2)*zx;
	b(61,9) = -(1.-.5*f2)*zx;
	b(62,7) = r3*f3*zx;
	b(62,8) = r3*(f1-zz)*zx;
	b(62,9) = -.5*r3*f1*zx;
	b(63,4) = 0.5*f2*r3;
	b(64,5) = -x*.5*r3*f2;
	b(64,6) = -x*(1.-f2);
	b(65,5) = -.5*r3*f2*y;
	b(65,6) = (1.+f2)*y;
	b(66,5) = -.5*r3*f2*z;
	b(66,6) = f2*z;
	b(67,7) = 1.5*f2*xy;
	b(67,8) = -2.*f2*xy;
	b(67,9) = .5*f2*xy;
	b(68,7) = 1.5*f2*yz;
	b(68,8) = -(1.+2.*f2)*yz;
	b(68,9) = (1.+.5*f2)*yz;
	b(69,7) = 1.5*f2*zx;
	b(69,8) = (1.-2.*f2)*zx;
	b(69,9) = -(1.-.5*f2)*zx;
	b(70,7) = .75*f2*f2;
	b(70,8) = (f1-f2*f2);
	b(70,9) = (zz+.25*f2*f2);
	b(71,7) = .5*f2*r3*f3;
	b(71,8) = -r3*zz*f2;
	b(71,9) = .25*r3*(1.+zz)*f2;
	b(72,4) = 0.5*f8;
	b(73,5) = -x*.5*f8;
	b(73,6) = r3*zz*x;
	b(74,5) = -.5*f8*y;
	b(74,6) = r3*zz*y;
	b(75,5) = -.5*f8*z;
	b(75,6) = -r3*f1*z;
	b(76,7) = r3*f3*xy;
	b(76,8) = -2.*r3*zz*xy;
	b(76,9) = .5*r3*(1.+zz)*xy;
	b(77,7) = r3*f3*yz;
	b(77,8) = r3*(f1-zz)*yz;
	b(77,9) = -.5*r3*f1*yz;
	b(78,7) = r3*f3*zx;
	b(78,8) = r3*(f1-zz)*zx;
	b(78,9) = -.5*r3*f1*zx;
	b(79,7) = .5*f2*r3*f3;
	b(79,8) = -r3*zz*f2;
	b(79,9) = .25*r3*(1.+zz)*f2;
	b(80,7) = f3*f3;
	b(80,8) = 3.*zz*f1;
	b(80,9) = .75*f1*f1;
	/* Matrix<double,10,1> test; */
        /* test(0) = -0.05276; */
        /* test(1) = -0.06725; */
        /* test(2) =  0.17199; */
        /* test(3) =  0.03242; */
        /* test(4) = -0.03035; */
        /* test(5) =  0.05202; */
        /* test(6) = -0.00092; */
        /* test(7) = -0.03367; */
        /* test(8) =  0.01000; */
        /* test(9) = -0.00075; */
        /* test(0) = -0.10605; */
        /* test(1) = -0.19326; */
        /* test(2) =  0.26932; */
        /* test(3) = -0.00672; */
        /* test(4) = -0.08018; */
        /* test(5) =  0.09915; */
        /* test(6) = -0.02467; */
        /* test(7) = -0.05669; */
        /* test(8) =  0.03738; */
        /* test(9) = -0.00669; */
	/* Matrix<double,81,1> test2; */
	/* test2 = b*test; */
	/* Matrix<double,9,9> test3; */
	/* test3.row(0) = test2.head(9); */
	/* test3.row(1) = test2.segment(9,9); */
	/* test3.row(2) = test2.segment(18,9); */
	/* test3.row(3) = test2.segment(27,9); */
	/* test3.row(4) = test2.segment(36,9); */
	/* test3.row(5) = test2.segment(45,9); */
	/* test3.row(6) = test2.segment(54,9); */
	/* test3.row(7) = test2.segment(63,9); */
	/* test3.row(8) = test2.tail(9); */
	Matrix<double,10,81> bdagg;
	bdagg = b.adjoint();
	Matrix<double,10,10> bbdg, bbdg_inv;
	bbdg = bdagg*b;
	bbdg_inv = bbdg.inverse();
	Matrix<double,81,1> tvec;
	tvec.head(9) = t.row(0);
	tvec.segment(9,9) = t.row(1);
	tvec.segment(18,9) = t.row(2);
	tvec.segment(27,9) = t.row(3);
	tvec.segment(36,9) = t.row(4);
	tvec.segment(45,9) = t.row(5);
	tvec.segment(54,9) = t.row(6);
	tvec.segment(63,9) = t.row(7);
	tvec.tail(9) = t.row(8);

	Matrix<double,10,81> tmp;
	tmp = bbdg_inv*bdagg;
	Matrix<double,10,1> V;
	V = tmp*tvec;

	return V;
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
	m.emplace_back(0);
	m.emplace_back(1);
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
	int ispin = -1;

	Vector3d d_0, d_1, d_2, d_3, d_4, d_5, d_6, d_7, d_8;
	Vector3d d_13, d_14, d_15, d_16;
	Vector3d d_17, d_18;
	/* double a = 6.692; */
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

	Matrix<dcomp, 9, 9> u;
	Matrix<dcomp, 9, 9> t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8;
	Matrix<dcomp, 9, 9> t_13, t_14, t_15, t_16;
	Matrix<dcomp, 9, 9> t_17, t_18;
	u = read(d_0, ispin);

	t_1 = read(d_1, ispin);
	t_3 = read(d_3, ispin);
	t_5 = read(d_5, ispin);
	t_13 = read(d_13, ispin);

	SelfAdjointEigenSolver<Matrix<dcomp, 9, 9>> uu;
	uu.compute(u);
	Matrix<double, 9, 9> lambdareal;
	Matrix<dcomp, 9, 9> Oo, Odagg, lambda;
	lambdareal = uu.eigenvalues().asDiagonal();
	lambda = lambdareal.cast<dcomp>();
	/* cout<<lambda.real()<<endl<<endl; */
	Oo = uu.eigenvectors();
	Odagg = Oo.adjoint();

	t_1 = Odagg*t_1*Oo;
	t_1 = convert(t_1);
	t_3 = Odagg*t_3*Oo;
	t_3 = convert(t_3);
	t_5 = Odagg*t_5*Oo;
	t_5 = convert(t_5);
	t_13 = Odagg*t_13*Oo;
	t_13 = convert(t_13);
	/* cout<<(Odagg*u*O).real()<<endl<<endl; */
	lambda = convert(lambda);
	Matrix<double,10,1> nn,nnn;
	Matrix<double,9,9> test;

	/* test = t_1.real(); */
	/* nn = sk_extraction(test, d_1); */

	/* test = t_3.real(); */
	/* nn = sk_extraction(test, d_3); */

	test = t_5.real();
	nn = sk_extraction(test, d_5);

	test = t_13.real();
	nnn = sk_extraction(test, d_13);

	//initialise onsite and hopping matrices for each nn

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

	/* Myfile<<"P X Y"<<endl; */
	Myfile<<"X Y"<<endl;

	dcomp i;
	i = -1.;
	i = sqrt(i);

	Matrix<dcomp, 9, 9> E;
	double b = 2.;

	/* Matrix<dcomp, 18, 18> E; */
	/* Matrix<dcomp, 9, 9> u_11, u_12; */
	/* Matrix<dcomp, 18, 18> U, T; */
	/* Matrix<complex<double>, 9, 9> zero = Matrix<complex<double>, 9, 9>::Zero(); */


	double k_x, k_y, k_z, pi;
	Vector3d K;
	for (int k = 0; k < 351; k++)
	{
		if (k < 101){
			pi = 2*M_PI*k/100.;
			k_z = pi/b;
			k_y = 0;
			k_x = 0;
		}
		if ((k > 100) && (k < 201)){
			pi = M_PI*(k-100)/100.;
			k_x = (2*M_PI-pi)/b;
			k_y = pi/b;
			k_z = pi/b;
		}	
		if ((k > 200) && (k < 251)){
			pi = M_PI*(k-200)/50.;
			k_x = (M_PI - pi)/b;
			k_y = M_PI/b;
			k_z = M_PI/b;
		}
		if ((k > 250) && (k < 351)){
			pi = M_PI*(k-250)/100.;
			k_x = 0;
			k_y = (M_PI-pi)/b;
			k_z = (M_PI-pi)/b;
		}
		K(0) = k_x;
		K(1) = k_y;
		K(2) = k_z;

		//fully diagonalised Hamiltonian
		E = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
			+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
				+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K))
				+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
				+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
				+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K)));

		/* if (k == 5) */
		/* 	cout<<E<<endl<<endl; */
		/* if (k == 15) */
		/* 	cout<<E<<endl<<endl; */
		/* if (k == 50) */
		/* 	cout<<E<<endl<<endl; */

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

		Myfile<<k<<" "<<O(0)<<endl;
		Myfile<<k<<" "<<O(1)<<endl;
		Myfile<<k<<" "<<O(2)<<endl;
		Myfile<<k<<" "<<O(3)<<endl;
		Myfile<<k<<" "<<O(4)<<endl;
		Myfile<<k<<" "<<O(5)<<endl;
		Myfile<<k<<" "<<O(6)<<endl;
		Myfile<<k<<" "<<O(7)<<endl;
		Myfile<<k<<" "<<O(8)<<endl;

		/* Myfile<<"A"<<" "<<k<<" "<<O(0)<<endl; */
		/* Myfile<<"B"<<" "<<k<<" "<<O(1)<<endl; */
		/* Myfile<<"C"<<" "<<k<<" "<<O(2)<<endl; */
		/* Myfile<<"D"<<" "<<k<<" "<<O(3)<<endl; */
		/* Myfile<<"E"<<" "<<k<<" "<<O(4)<<endl; */
		/* Myfile<<"F"<<" "<<k<<" "<<O(5)<<endl; */
		/* Myfile<<"G"<<" "<<k<<" "<<O(6)<<endl; */
		/* Myfile<<"H"<<" "<<k<<" "<<O(7)<<endl; */
		/* Myfile<<"I"<<" "<<k<<" "<<O(8)<<endl; */

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
