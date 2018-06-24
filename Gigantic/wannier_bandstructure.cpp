#include <iostream>
#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <vector>
#include "TB_sk.h"
#include <iomanip>
/* #include "TBdynamic.h" */
//Very good matching of SK potentials obtained through comparison of eigenvalues through interpolation scheme, for the 4x4 matrix s,p but not o good for the rest yet
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
	string Mydata, Mydata2, Mydata3, Mydata4;
	getline(cin, Mydata);
	ofstream Myfile, Myfile2, Myfile3, Myfile4;	
	Mydata2 = Mydata;
	/* Mydata3 = Mydata; */
	Mydata4 = Mydata;
	Mydata += ".txt";
	Mydata2 += ".dat";
	/* Mydata3 += "_catch.dat"; */
	Mydata4 += "_eigs.dat";
	Myfile.open( Mydata.c_str(),ios::trunc );
	Myfile2.open( Mydata2.c_str(),ios::trunc );
	/* Myfile3.open( Mydata3.c_str(),ios::trunc ); */
	Myfile4.open( Mydata4.c_str(),ios::trunc );

	cout<<endl;

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

	SelfAdjointEigenSolver<Matrix<dcomp, 9, 9>> uu;
	uu.compute(u);
	Matrix<double, 9, 9> lambdareal;
	Matrix<dcomp, 9, 9> Oo, Odagg, lambda;
	lambdareal = uu.eigenvalues().asDiagonal();
	lambda = lambdareal.cast<dcomp>();
	/* cout<<lambda.real()<<endl<<endl; */
	Oo = uu.eigenvectors();
	Odagg = Oo.adjoint();

	/* t_1 = read(d_1, ispin); */
	/* t_1 = Odagg*t_1*Oo; */
	/* t_1 = convert(t_1); */
	/* t_2 = read(d_2, ispin); */
	/* t_2 = Odagg*t_2*Oo; */
	/* t_2 = convert(t_2); */
	/* t_14 = read(d_14, ispin); */
	/* t_14 = Odagg*t_14*Oo; */
	/* t_14 = convert(t_14); */
	/* t_13 = read(d_13, ispin); */
	/* t_13 = Odagg*t_13*Oo; */
	/* t_13 = convert(t_13); */

	/* cout<<"t_1"<<endl; */
	/* cout<<t_1.real()<<endl<<endl; */
	/* cout<<"t_13"<<endl; */
	/* cout<<t_13.real()<<endl<<endl; */

	/* cout<<(Odagg*u*O).real()<<endl<<endl; */
	lambda = convert(lambda);
	/* cout<<lambda<<endl<<endl; */

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

	/* cout<<"t_15"<<endl<<t_15.real()<<endl<<endl; */
	/* cout<<"t_16"<<endl<<t_16.real()<<endl<<endl; */
	/* cout<<"t_17"<<endl<<t_17.real()<<endl<<endl; */
	/* cout<<"t_18"<<endl<<t_18.real()<<endl<<endl; */

	/* cout<<lambda<<endl<<endl; */
	/* cout<<fixed<<(t_13+t_14+t_15+t_16+t_17+t_18).real()<<endl<<endl; */
	/* cout<<fixed<<(t_1+t_2+t_3+t_4+t_5+t_6+t_7+t_8).real()<<endl<<endl; */


	//initialise onsite and hopping matrices for each nn

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

	Matrix<double,10,1> nn1,nn2,nn3,nn4,nn5,nn6,nn7,nn8,nnav,nnn1,nnn2,nnn3,nnn4,nnn5,nnn6,nnnav;
	Matrix<double,9,9> test;

	test = t_1.real();
	nn1 = sk_extraction(test, d_1);
	test = t_2.real();
	nn2 = sk_extraction(test, d_2);
	test = t_3.real();
	nn3 = sk_extraction(test, d_3);
	test = t_4.real();
	nn4 = sk_extraction(test, d_4);
	test = t_5.real();
	nn5 = sk_extraction(test, d_5);
	test = t_6.real();
	nn6 = sk_extraction(test, d_6);
	test = t_7.real();
	nn7 = sk_extraction(test, d_7);
	test = t_8.real();
	nn8 = sk_extraction(test, d_8);
	nnav = (nn1+nn2+nn3+nn4+nn5+nn6+nn7+nn8)/8.;

	test = t_13.real();
	nnn1 = sk_extraction(test, d_13);
	test = t_14.real();
	nnn2 = sk_extraction(test, d_14);
	test = t_15.real();
	nnn3 = sk_extraction(test, d_15);
	test = t_16.real();
	nnn4 = sk_extraction(test, d_16);
	test = t_17.real();
	nnn5 = sk_extraction(test, d_17);
	test = t_18.real();
	nnn6 = sk_extraction(test, d_18);
	nnnav = (nnn1+nnn2+nnn3+nnn4+nnn5+nnn6)/6.;
	cout<<"nn1 = "<<nn1.transpose()<<endl;
	cout<<"nn2 = "<<nn2.transpose()<<endl;
	cout<<"nn3 = "<<nn3.transpose()<<endl;
	cout<<"nn4 = "<<nn4.transpose()<<endl;
	cout<<"nn5 = "<<nn5.transpose()<<endl;
	cout<<"nn6 = "<<nn6.transpose()<<endl;
	cout<<"nn7 = "<<nn7.transpose()<<endl;
	cout<<"nn8 = "<<nn8.transpose()<<endl;
	cout<<"nnav = "<<nnav.transpose()<<endl;
	cout<<"nnn1 = "<<nnn1.transpose()<<endl;
	cout<<"nnn2 = "<<nnn2.transpose()<<endl;
	cout<<"nnn4 = "<<nnn3.transpose()<<endl;
	cout<<"nnn4 = "<<nnn4.transpose()<<endl;
	cout<<"nnn5 = "<<nnn5.transpose()<<endl;
	cout<<"nnn6 = "<<nnn6.transpose()<<endl;
	cout<<"nnnav= "<<nnnav.transpose()<<endl;

	dcomp i;
	i = -1.;
	i = sqrt(i);

	Matrix<dcomp, 9, 9> E_gamma, E_H, E_P, E_N;
	double k_x, k_y, k_z, pi;
	double b = 2.;
	Vector3d K;
	pi = 2*M_PI;
	K << 0, 0, 0;
	E_gamma = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K)));
	K << 0, 0, pi/b;
	E_H = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K)));
	K << M_PI/b, M_PI/b, M_PI/b;
	E_P = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K)));
	K << M_PI/b, M_PI/b, 0;
	E_N = lambda + (t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K))
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K)));
	double s, p, d1, d2, sss1, sss2, pps1, pps2, ppp1, ppp2, dds1, dds2, ddp1, ddp2, ddd1, ddd2, sps1, sps2, sds1, sds2, pds1, pds2, pdp1, pdp2;
	sss1 = real((1./16.)*(E_gamma(0,0)-E_H(0,0)));
	sss2 = real((1./24.)*(E_gamma(0,0)+E_H(0,0)-2.*E_P(0,0)));// very important - typo in intro of Papa p12
	s = 0.25*real((E_gamma(0,0) + E_H(0,0) + 2.*E_P(0,0)));
	p = 0.125*real(E_gamma(1,1) + E_H(1,1) + 2.*(E_N(2,2) + E_N(3,3) + E_N(1,1)));
	cout<<"onsite terms"<<endl;
	cout<<"s = "<<s<<endl;
	cout<<"p = "<<p<<endl;
	d1 = real(lambda(5,5));
	cout<<"d1 = "<<d1<<endl;
	d2 = 0.25*real(E_gamma(7,7) + E_H(7,7) + 2.*E_P(7,7));
	cout<<"d2 = "<<d2<<endl;
	double xx1, xy1, r3;
	r3 = 1./sqrt(3.);
	xx1 = (1./16.)*real(E_gamma(1,1) - E_H(1,1));
	xy1 = -0.125*real(E_N(2,3));
  
	/* //this block obtained from matching eigenvalues, see block below */
	/* sps1 = 0.20140; */
	/* pps1 = 0.00725; */
	/* ppp1 = 0.11224; */

	/* //this block obtained from matching eigenvalues, a second attempt */
	/* sss1 = -0.11807; */
	/* sps1 = -0.16907; */
	/* pps1 = 0.07499; */
	/* ppp1 = 0.07835; */

	/* pps1 = 2.*xy1 + xx1;//obtained from Papa formula in intro *1/ */
	/* ppp1 = xx1 - xy1;//obtained from Papa formula in intro *1/ */
	cout<<endl;
	cout<<"1st neighbours"<<endl;

	cout<<"sss = "<<sss1<<endl;
	/* cout<<"pps = "<<pps1<<endl; */
	/* cout<<"ppp = "<<ppp1<<endl; */

	double xyxy1;
	xyxy1 = (1./16.)*real(E_gamma(4,4)-E_H(6,6));
	double d2d21;
	d2d21 = (1./16.)*real(E_gamma(7,7)-E_H(7,7));
	double xyyz1;
	xyyz1 = -(1./8.)*real(E_N(5,6));
	ddd1 = (6./5.)*(xyxy1-xyyz1) - (3./5.)*d2d21;
	ddp1 = (3./5.)*(2.*d2d21 + xyxy1 - xyyz1);
	dds1 = xyxy1 + 2.*xyyz1;
	cout<<endl;
	cout<<"These are Papa eq. first nn SK potentials {"<<endl;
	cout<<"ddd = "<<ddd1<<endl<<"ddp = "<<ddp1<<endl<<"dds = "<<dds1<<endl;
	cout<<"}"<<endl;

	/* cout<<"dds = "<<dds1<<endl; */
	/* cout<<"ddp = "<<ddp1<<endl; */
	/* cout<<"ddd = "<<ddd1<<endl; */

	//this block obtained from matching eigenvalues, see block below accurate!
	ppp1 = -0.018593;
	pps1 = 0.2689;
	sps1 = -0.16918;

	cout<<"ppp = "<<ppp1<<endl;
	cout<<"pps = "<<pps1<<endl;
	cout<<"sps = "<<sps1<<endl;

	/* ppp2 = (1./16.)*real(E_gamma(2,1) + E_H(1,1) - 2.*E_N(1,1)); //obtained from Papa formula in intro */
	/* pps2 = (1./16.)*real(E_gamma(1,1) + E_H(1,1) + 2.*(-E_N(2,2) - E_N(3,3) + E_N(1,1)));//obtained from Papa formula in intro */
	cout<<endl;
	//this block obtained from matching eigenvalues, see block below accurate!
	ppp2 = 0.03060;
	pps2 = 0.16341; 
	sps2 = -0.06189; 

	cout<<"2nd neighbours"<<endl;
	cout<<"sss = "<<sss2<<endl;
	cout<<"pps = "<<pps2<<endl;
	cout<<"ppp = "<<ppp2<<endl;
	double delta;
	//magic 1.693948 below extracted from 5th element of onsite matrix (x2)
	/* delta = 1.693948 + real(-0.5*(E_gamma(4,4) + E_H(6,6)) - E_N(5,5) - E_N(6,6)); */
	ddd2 = 0.25*real(E_N(7,7) - E_P(7,7));
	double ddd2_alt;
	/* ddd2_alt = (1./8.)*real(E_gamma(4,4) + E_H(6,6) - 2.*delta); */
	delta = -4.*ddd2 + real(0.5*E_gamma(4,4) + 0.5*E_H(6,6));
	ddp2 = (1./8.)*real(2.*delta - 1.5*(E_N(5,5) + E_N(6,6)) - E_N(5,6));
	dds2 = (1./12.)*real(E_gamma(7,7) + E_H(7,7) + E_P(7,7) - (3./2.)*(E_N(5,5) + 2.*E_N(5,6) + E_N(6,6)));
	cout<<endl;
	cout<<"These are Papa eq. second nn SK potentials {"<<endl;
	cout<<"dds = "<<dds2<<endl;
	cout<<"ddp = "<<ddp2<<endl;
	cout<<"ddd = "<<ddd2<<endl;
	cout<<"}"<<endl;
	/* cout<<"ddd = "<<ddd2_alt<<endl; */
	cout<<endl;
	cout<<"delta = "<<delta<<endl;
	ddd2_alt = 0.25*real(E_gamma(4,4) + E_H(6,6) + 2.*E_N(5,5) + 2.*E_N(6,6) + 2*delta);
	cout<<"unknown potential = "<<ddd2_alt<<endl<<endl;

	/* cout<<lambda.real()<<endl; */

	/* K << 0, M_PI/b, M_PI/b; */
	K << 3.*M_PI/(2.*b), M_PI/(2.*b), M_PI/(2.*b);
	Matrix<dcomp, 9, 9> E_AA;
	/* cout<<exp(i*d_17.dot(K))<<" "<<exp(i*d_18.dot(K))<<endl; */
	/* E_AA = t_13; */
	E_AA = lambda 
		/* + t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K)) */
		/* + t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K)) */
		/* 	+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)) */
			+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
			+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
			+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
			;
	Matrix4cd AA, BB;
	AA = E_AA.topLeftCorner(4,4);
	/* cout<<AA<<endl<<endl; */
	/* ComplexEigenSolver<Matrix<dcomp,9,9>> CA; */
	ComplexEigenSolver<Matrix<dcomp,4,4>> CA, CB;
	/* CA.compute(E_AA); */
	CA.compute(AA);
	Vector4d eigs, eigs2;
	eigs = CA.eigenvalues().real();
	/* cout<<CA.eigenvalues().transpose()<<endl<<endl; */
	/* cout<<CA.eigenvalues().sum()<<endl<<endl; */

	Matrix<dcomp, 5, 5> CC, DD;
	CC = E_AA.bottomRightCorner(5,5);
	ComplexEigenSolver<Matrix<dcomp, 5, 5>> DA, DB;
	DA.compute(CC);
	Matrix<double, 5, 1> deigs, deigs2;
	deigs = DA.eigenvalues().real();

	double a1, a2, a3, a4, a5;
	double b1, b2, b3;
	dcomp c1, c2, c3, c4, c5;

	/* cout<<t_13.real()<<endl<<endl; */
	/* dds1 = -0.04253; */
	/* ddp1 = 0.01857; */
	/* ddd1 = 0.000316; */
	
	/* a1 = (1./9.)*(3.*dds1 + 2.*ddp1 + 4.*ddd1); */
	/* a2 = (1./9.)*(3.*dds1 - ddp1 - 2.*ddd1); */
	/* a3 = (sqrt(3.)/9.)*(ddp1 - ddd1); */
	/* a4 = (1./3.)*(ddp1 - ddd1); */
	/* a5 = (1./3.)*(2.*ddp1 + ddd1); */
	/* DD<< a1,  a2, a2,   0, -2.*a3, */ 
	/*      a2,  a1, a2, -a4,     a3, */
	/*      a2,  a2, a1,  a4,     a3, */
	/*       0, -a4, a4,  a5,      0, */
	/*  -2.*a3,  a3, a3,   0,     a5; */

	/* b1 = r3*sps1; */
	/* b2 = (1./3.)*(pps1 + 2.*ppp1); */
	/* b3 = (1./3.)*(pps1 - ppp1); */
	/* BB<<sss1, b1, b1, b1, */
	/*      -b1, b2, b3, b3, */
	/*      -b1, b3, b2, b3, */
	/*      -b1, b3, b3, b2; */
	/* Matrix<dcomp, 4, 5> BBBB; */
	/* Matrix<dcomp, 5, 4> CCCC; */
	/* Matrix<dcomp, 9, 9> MASTER; */
	/* ComplexEigenSolver<Matrix<dcomp, 9, 9>> EA, EB; */
	/* EA.compute(E_AA); */
	/* Matrix<double, 9, 1> feigs, feigs2; */
	/* feigs = EA.eigenvalues().real(); */
	/* cout<<feigs.transpose()<<endl<<endl; */

	/* sds1 = -0.02; */
	/* pds1 = 0.07; */
	/* pdp1 = -0.0001; */

	/* ddd2 = -0.00678; */

	/* cout<<t_13.topLeftCorner(4,4)<<endl<<endl; */

	BB<< s - 2.*sss2,                     0,           0,           0,
	               0, p + 2.*pps2 - 4.*ppp2,           0,           0,
	               0,                     0, p - 2.*pps2,           0,
	               0,                     0,           0, p - 2.*pps2;
	Matrix<dcomp, 4, 5> BBBB;
	Matrix<dcomp, 5, 4> CCCC;
	Matrix<dcomp, 9, 9> MASTER;
	ComplexEigenSolver<Matrix<dcomp, 9, 9>> EA, EB;
	EA.compute(E_AA);
	Matrix<double, 9, 1> feigs, feigs2;
	feigs = EA.eigenvalues().real();
	sort(feigs.data(), feigs.data()+feigs.size());
	cout<<feigs.transpose()<<endl<<endl;
	cout<<feigs.sum()<<" "<<feigs.prod()<<endl;
	//alternatively sds2 could be positive below
	sds2 = -0.02805; dds2 = -0.02267; ddp2 = -0.00468; ddd2 = 0.00209;
	/* sds2 = -0.02775; dds2 = -0.0224; ddp2 = -0.00488; ddd2 = 0.00279; */

	/* this computed at K=(0,pi,pi)*/

	/* for (double pot = -0.02824; pot <= -0.02804; pot += 0.00001){ */
	/* 	cout<<pot<<endl<<endl; */
	/* 	/1* if ((pot > -0.02800) && (pot < -0.02729)) *1/ */
	/* 	/1* 	pot = 0.02800; *1/ */
	/* for (double pot1 = -0.02305; pot1 <= -0.02249; pot1 += 0.00001){ */
	/* for (double pot2 = -0.00492; pot2 <= -0.00447; pot2 += 0.00001){ */
	/* for (double pot3 = 0.00209; pot3 <= 0.00226; pot3 += 0.00001){ */
	/* 	sds2 = pot; */
	/* 	dds2 = pot1; */
	/* 	ddp2 = pot2; */
	/* 	ddd2 = pot3; */
	/* 	a1 = d1 - 2.*ddd2; */
	/* 	a2 = d1 + 2.*ddd2 - 4.*ddp2; */
	/* 	a3 = d2 - 2*ddd2; */
	/* 	a4 = sqrt(3.)*(ddd2 - dds2); */
	/* 	a5 = d2 - 2*dds2; */
	/* 	DD<< a1, 0, 0,  0,  0, */ 
	/* 	     0, a2, 0,  0,  0, */
	/* 	     0, 0, a1,  0,  0, */
	/* 	     0, 0,  0, a3, a4, */
	/* 	     0, 0,  0, a4, a5; */
	/* 	c1 = 2.*sqrt(3.)*sds2; */
	/* 	c2 = -2.*sds2; */
	/* 	BBBB<<  0, 0, 0, c1, c2, */
	/* 	        0, 0, 0,  0,  0, */
	/* 	        0, 0, 0,  0,  0, */
	/* 	        0, 0, 0,  0,  0; */
	/* 	CCCC<<  0, 0, 0, 0, */
	/* 	        0, 0, 0, 0, */
	/* 		0, 0, 0, 0, */
	/* 	       c1, 0, 0, 0, */
	/* 	       c2, 0, 0, 0; */
	/* 	MASTER.topLeftCorner(4,4) = BB; */
	/* 	MASTER.bottomRightCorner(5,5) = DD; */
	/* 	MASTER.topRightCorner(4,5) = BBBB; */
	/* 	MASTER.bottomLeftCorner(5,4) = CCCC; */
	/* 	EB.compute(MASTER); */
	/* 	feigs2 = EB.eigenvalues().real(); */
	/* 	sort(feigs2.data(), feigs2.data()+feigs2.size()); */
	/* 	/1* cout<<eigs2.transpose()<<endl<<endl; *1/ */
	/* 	if ((abs(feigs(0)-feigs2(0)) < 1.6e-3) && (abs(feigs(1)-feigs2(1)) < 2.5e-3) */
	/* 		       	&& (abs(feigs(2)-feigs2(2)) < 2.5e-3) && (abs(feigs(3)-feigs2(3)) < 1.6e-3) */
	/* 		       	&& (abs(feigs(4)-feigs2(4)) < 1.6e-3) && (abs(feigs(5)-feigs2(5)) < 2.5e-3) */
	/* 		       	&& (abs(feigs(6)-feigs2(6)) < 2.5e-3) && (abs(feigs(7)-feigs2(7)) < 2.5e-3) */
	/* 		       	&& (abs(feigs(8)-feigs2(8)) < 2.5e-3) && (abs(feigs.sum()-feigs2.sum())<1e-5)){ */
	/* 		Myfile<<setprecision(9)<<"sds2 = "<<sds2<<"; " */
	/* 			<<"dds2 = "<<dds2<<"; "<<"ddp2 = "<<ddp2<<"; "<<"ddd2 = "<<ddd2<<";"<<endl; */
	/* 		Myfile2<<setprecision(9)<<sds2<<" "<<dds2<<" "<<ddp2<<" "<<ddd2<<" "<<abs(feigs.sum()-feigs2.sum())<<endl; */
	/* 		Myfile4<<feigs2.transpose()<<endl; */
	/* 		cout<<feigs.transpose()<<endl<<feigs2.transpose()<<endl<<endl; */
	/* 		/1* break; *1/ */
	/* 	} */
	/* } */
	/* } */
	/* } */
	/* } */

	/* //K << 3.*M_PI/(2.*b), M_PI/(2.*b), M_PI/(2.*b); */
	/* BB<< s        , -2.*i*sps2, 2.*i*sps2, 2.*i*sps2, */
	/*      2.*i*sps2,          p,         0,         0, */
	/*     -2.*i*sps2,          0,         p,         0, */
	/*     -2.*i*sps2,          0,         0,         p; */
	/* DD<< d1, 0, 0,  0,  0, */ 
	/*      0, d1, 0,  0,  0, */
	/*      0, 0, d1,  0,  0, */
	/*      0, 0,  0, d2,  0, */
	/*      0, 0,  0,  0, d2; */
	/* //-ve chosen, but it appears sign is irrelevant for pdp2 & pds2 */
	/* for (double pot = -0.0120; pot <= -0.0100; pot += 0.00001){ */
	/* for (double pot1 = -0.0540; pot1 <= -0.0500; pot1 += 0.00001){ */
	/* 	pdp2 = pot; */
	/* 	pds2 = pot1; */
	/* 	c1 = i*2.*pdp2; */
	/* 	c2 = i*pds2; */
	/* 	c3 = -i*sqrt(3.)*pds2; */
	/* 	BBBB<<  0,  0,   0,  0,     0, */
	/* 	       c1,  0,  c1, c3,    c2, */
	/* 	      -c1, c1,   0, c3,   -c2, */
	/* 	        0, c1, -c1,  0, 2.*c2; */
	/* 	CCCC<< 0, -c1,  c1,      0, */
	/* 	       0,   0, -c1,    -c1, */
	/* 	       0, -c1,   0,     c1, */
	/* 	       0, -c3, -c3,      0, */
	/* 	       0, -c2,  c2, -2.*c2; */
	/* 	MASTER.topLeftCorner(4,4) = BB; */
	/* 	MASTER.bottomRightCorner(5,5) = DD; */
	/* 	MASTER.topRightCorner(4,5) = BBBB; */
	/* 	MASTER.bottomLeftCorner(5,4) = CCCC; */
	/* 	EB.compute(MASTER); */
	/* 	feigs2 = EB.eigenvalues().real(); */
	/* 	sort(feigs2.data(), feigs2.data()+feigs2.size()); */
	/* 	/1* cout<<eigs2.transpose()<<endl<<endl; *1/ */
	/* 	if ((abs(feigs(0)-feigs2(0)) < 1.5e-5) && (abs(feigs(1)-feigs2(1)) < 1.5e-5) */
	/* 		       	&& (abs(feigs(2)-feigs2(2)) < 1.5e-5) && (abs(feigs(3)-feigs2(3)) < 1.5e-5) */
	/* 		       	&& (abs(feigs(4)-feigs2(4)) < 1.5e-5) && (abs(feigs(5)-feigs2(5)) < 1.5e-5) */
	/* 		       	&& (abs(feigs(6)-feigs2(6)) < 1.5e-5) && (abs(feigs(7)-feigs2(7)) < 1.5e-5) */
	/* 		       	&& (abs(feigs(8)-feigs2(8)) < 1.5e-5) && (abs(feigs.sum()-feigs2.sum())<3e-5)){ */
	/* 		Myfile<<setprecision(9)<<"pdp2 = "<<pdp2<<"; "<<"pds2 = "<<pds2<<";"<<endl; */
	/* 		Myfile2<<setprecision(9)<<pdp2<<" "<<pds2<<" "<<abs(feigs.sum()-feigs2.sum())<<endl; */
	/* 		Myfile4<<feigs2.transpose()<<endl; */
	/* 		cout<<feigs.transpose()<<endl<<feigs2.transpose()<<endl<<endl; */
	/* 		/1* break; *1/ */
	/* 	} */
	/* } */
	/* } */

	pdp2 = -0.01088; pds2 = -0.05257;

	/* sds2 = -0.01490; pds2 = 0.04464; pdp2 = -0.01124; dds2 = -0.00668; ddp2 = -0.00138; ddd2 = -0.01802; */
	//Given the apparent symmetry, for pot -ve, pot1 +ve and vice versa over same energy range. Also pot2 seems sign independent.
	//only -ve pot, +ve pot1 and -ve pot2
	/* for (double pot = -0.01504; pot <= -0.01474; pot += 0.00002){ */
	/* 	cout<<pot<<endl<<endl; */
	/* 	/1* if ((pot > -0.01881) && (pot < -0.01879)) *1/ */
	/* 	/1* 	pot = -0.01510; *1/ */
	/* 	/1* if (pot < 0){ *1/ */
	/* 	/1*        start = 0.076; *1/ */
	/* 	/1*        end = 0.11; *1/ */
	/* 	/1* } *1/ */
	/* 	/1* else{ *1/ */
	/* 	/1* 	start = -0.1; *1/ */
	/* 	/1* 	end = -0.076; *1/ */
	/* 	/1* } *1/ */
	/* for (double pot1 = 0.04456; pot1 <= 0.04468; pot1 += 0.00002){ */
	/* for (double pot2 = -0.01124; pot2 <= -0.01124; pot2 += 0.00002){ */
	/* for (double pot3 = -0.00676; pot3 <= -0.00642; pot3 += 0.00002){ */
	/* for (double pot4 = -0.00138; pot4 <= -0.00138; pot4 += 0.00002){ */
	/* for (double pot5 = -0.01802; pot5 <= -0.01802; pot5 += 0.00002){ */
	/* 	sds2 = pot; */
	/* 	pds2 = pot1; */
	/* 	pdp2 = pot2; */
	/* 	dds2 = pot3; */
	/* 	ddp2 = pot4; */
	/* 	ddd2 = pot5; */
	/* 	a1 = ddp2; */
	/* 	a2 = ddd2; */
	/* 	a3 = 0.75*dds2 + 0.25*ddd2; */
	/* 	a4 = 0.25*sqrt(3.)*(ddd2 - dds2); */
	/* 	a5 = 0.25*dds2 + 0.75*ddd2; */
	/* 	DD<< a1, 0, 0,  0,  0, */ 
	/* 	     0, a2, 0,  0,  0, */
	/* 	     0, 0, a1,  0,  0, */
	/* 	     0, 0,  0, a3, a4, */
	/* 	     0, 0,  0, a4, a5; */
	/* 	c1 = 0.5*sqrt(3.)*sds2; */
	/* 	c2 = 0.5*sqrt(3.)*pds2; */
	/* 	c3 = -0.5*sds2; */
	/* 	c4 = -0.5*pds2; */
	/* 	c5 = pdp2; */
	/* 	BBBB<<  0, 0,  0, c1, c3, */
	/* 	        0, 0,  0, c2, c4, */
	/* 	       c5, 0,  0,  0,  0, */
	/* 	        0, 0, c5,  0,  0; */
	/* 	CCCC<<  0,   0, -c5,   0, */
	/* 	        0,   0,   0,   0, */
	/* 		0,   0,   0, -c5, */
	/* 	       c1, -c2,   0,   0, */
	/* 	       c3, -c4,   0,   0; */
	/* 	MASTER.topLeftCorner(4,4) = BB; */
	/* 	MASTER.bottomRightCorner(5,5) = DD; */
	/* 	MASTER.topRightCorner(4,5) = BBBB; */
	/* 	MASTER.bottomLeftCorner(5,4) = CCCC; */
	/* 	EB.compute(MASTER); */
	/* 	feigs2 = EB.eigenvalues().real(); */
	/* 	sort(feigs2.data(), feigs2.data()+feigs2.size()); */
	/* 	/1* cout<<eigs2.transpose()<<endl<<endl; *1/ */
	/* 	if ((abs(feigs(0)-feigs2(0))/abs(feigs(0)) < 5e-2) && */
	/* 			(abs(feigs(4)-feigs2(4))/abs(feigs(4)) < 2.9e-2) */
	/* 			 && (abs(feigs2(2)) < 1e-4) */ 
	/* 			 && (abs(feigs(5)-feigs2(5))/abs(feigs(5)) < 1e-2) */
	/* 			 && (abs(feigs(6)-feigs2(6))/abs(feigs(6)) < 2.8e-2) */
	/* 			 && (abs(feigs(7)-feigs2(7))/abs(feigs(7)) < 2.8e-2) && */
	/* 			(abs(feigs(8)-feigs2(8))/abs(feigs(8)) < 1e-4)){ */
	/* 		Myfile<<setprecision(9)<<"sds2 = "<<sds2<<"; "<<"pds2 = "<<pds2<<"; "<< */
	/* 			"pdp2 = "<<pdp2<<"; "<<"dds2 = "<<dds2<<"; "<<"ddp2 = "<<ddp2<<"; "<<"ddd2 = "<<ddd2<<";"<<endl; */
	/* 		Myfile2<<setprecision(9)<<sds2<<" "<<pds2<<" "<<pdp2<<" "<<dds2<<" "<<ddp2<<" "<<ddd2<<endl; */
	/* 		Myfile4<<feigs2.transpose()<<endl; */
	/* 		cout<<feigs.transpose()<<endl<<feigs2.transpose()<<endl<<endl; */
	/* 		break; */
	/* 	} */
	/* 	/1* else if ((abs(feigs(0)-feigs2(0))/abs(feigs(0)) < 1e-1) && *1/ */
	/* 	/1* 		(abs(feigs(4)-feigs2(4))/abs(feigs(4)) < 1e-1) *1/ */
	/* 	/1* 		 && (abs(feigs(5)-feigs2(5))/abs(feigs(5)) < 1e-1) *1/ */
	/* 	/1* 		 && (abs(feigs(6)-feigs2(6))/abs(feigs(6)) < 1e-1) *1/ */
	/* 	/1* 		 && (abs(feigs(7)-feigs2(7))/abs(feigs(7)) < 1e-1) && *1/ */
	/* 	/1* 		(abs(feigs(8)-feigs2(8))/abs(feigs(8)) < 1e-1)){ *1/ */
	/* 	/1* 	Myfile3<<setprecision(9)<<sds2<<" "<<pds2<<" "<<pdp2<<" "<<dds2<<" "<<ddp2<<" "<<ddd2<<endl; *1/ */
	/* 	/1* 	break; *1/ */
	/* 	/1* } *1/ */
	/* } */
	/* } */
	/* } */
	/* } */
	/* } */
	/* } */

	/* for (double xx = 1.; xx >= 0.; xx -= 0.01){ */
	/* 	for (double yy = 1.; yy >= 0.; yy -= 0.01){ */
	/* 		for (double zz = 1.; zz >= 0.; zz -= 0.01){ */
	/* 			double result; */
	/* 			result = sqrt(3.)*(xx*xx - yy*yy)/(3.*zz*zz - 1.); */
	/* 			if (result < 0){ */
	/* 				if (((-result > 0.474) && (-result < 0.479)) || ((-result > 2.097) && (-result < 2.102))) */
	/* 					cout<<result<<" x = "<<xx<<" y = "<<yy<<" z = "<<zz<<endl; */
	/* 			} */
	/* 		} */
	/* 	} */
	/* } */


	/* //This block uses the bottom right 5x5 submatrix of t_13 to find the SK potentials of d orbitals */
	/* for (double pot = -0.02667; pot <= -0.02086; pot += 0.00001){ */
	/* for (double pot1 = 0.00043; pot1 <= 0.00328; pot1 += 0.00001){ */
	/* 	dds2 = pot; */
	/* 	ddp2 = pot1; */
	/* 	a1 = (sqrt(3.)/4.)*(ddd2-dds2); */
	/* 	DD<< ddp2,    0,    0,                     0,                     0, */ 
	/* 	        0, ddd2,    0,                     0,                     0, */
	/* 	        0,    0, ddp2,                     0,                     0, */
	/* 	        0,    0,    0, 0.75*dds2 + 0.25*ddd2,                    a1, */
	/* 	        0,    0,    0,                    a1, 0.25*dds2 + 0.75*ddd2; */
	/* 	DB.compute(DD); */
	/* 	deigs2 = DB.eigenvalues().real(); */
	/* 	/1* cout<<eigs2.transpose()<<endl<<endl; *1/ */
	/* 	if ((abs(deigs(0)-deigs2(0)) < 2.9e-3) && (abs(deigs(1)-deigs2(1)) < 2.9e-3) && (abs(deigs(2)-deigs2(2)) < 2.9e-3) && (abs(deigs(3)-deigs2(3)) < 2.9e-3) && (abs(deigs(4)-deigs2(4)) < 2.9e-3)){ */
	/* 		Myfile<<setprecision(9)<<dds2<<" "<<ddp2<<" "<<ddd2<<endl; */
	/* 		cout<<deigs2.transpose()<<endl; */
	/* 		break; */
	/* 	} */
	/* } */
	/* } */

	/* double start, end; */
	/* sds1 = -0.076; pds1 = 0.108; pdp1 = 0.019; dds1 = -0.026; ddp1 = 0.02; ddd1 = 0.011; */
	sds1 = -0.072; pds1 = 0.108; pdp1 = -0.02; dds1 = -0.036; ddp1 = 0.021; ddd1 = -0.009;
	cout<<"These are my first nn SK potentials {"<<endl;
	cout<<"ddd = "<<ddd1<<endl<<"ddp = "<<ddp1<<endl<<"dds = "<<dds1<<endl;
	cout<<"sds = "<<sds1<<endl<<"pds = "<<pds1<<endl<<"pdp = "<<pdp1<<endl;
	cout<<"}"<<endl;
	/* for (double pot = -0.074; pot <= -0.070; pot += 0.0001){ */
	/* 	/1* if (pot == -0.036) *1/ */
	/* 	/1* 	pot = 0.037; *1/ */
	/* 	/1* if (pot < 0){ *1/ */
	/* 	/1*        start = 0.076; *1/ */
	/* 	/1*        end = 0.11; *1/ */
	/* 	/1* } *1/ */
	/* 	/1* else{ *1/ */
	/* 	/1* 	start = -0.1; *1/ */
	/* 	/1* 	end = -0.076; *1/ */
	/* 	/1* } *1/ */
	/* for (double pot1 = 0.105; pot1 <= 0.109; pot1 += 0.0001){ */
	/* for (double pot2 = -0.02; pot2 <= -0.019; pot2 += 0.0001){ */
	/* for (double pot3 = -0.037; pot3 <= -0.034; pot3 += 0.0001){ */
	/* for (double pot4 = 0.020; pot4 <= 0.021; pot4 += 0.0001){ */
	/* for (double pot5 = -0.008; pot5 <= 0.011; pot5 += 0.0001){ */
	/* 	sds1 = pot; */
	/* 	pds1 = pot1; */
	/* 	pdp1 = pot2; */
	/* 	dds1 = pot3; */
	/* 	ddp1 = pot4; */
	/* 	ddd1 = pot5; */
	/* 	a1 = (1./9.)*(3.*dds1 + 2.*ddp1 + 4.*ddd1); */
	/* 	a2 = (1./9.)*(3.*dds1 - ddp1 - 2.*ddd1); */
	/* 	a3 = (sqrt(3.)/9.)*(ddp1 - ddd1); */
	/* 	a4 = (1./3.)*(ddp1 - ddd1); */
	/* 	a5 = (1./3.)*(2.*ddp1 + ddd1); */
	/* 	DD<< a1,  a2, a2,   0, -2.*a3, */ 
	/* 	     a2,  a1, a2, -a4,     a3, */
	/* 	     a2,  a2, a1,  a4,     a3, */
	/* 	      0, -a4, a4,  a5,      0, */
	/* 	 -2.*a3,  a3, a3,   0,     a5; */
	/* 	c1 = r3*sds1; */
	/* 	c2 = r3*(r3*pds1 + (1./3.)*pdp1); */
	/* 	c3 = (r3/3.)*(sqrt(3.)*pds1 - 2.*pdp1); */
	/* 	c4 = r3*pdp1; */
	/* 	c5 = (1./3.)*pdp1; */
	/* 	BBBB<< c1, c1, c1,  0,     0, */
	/* 	       c2, c3, c2, c4,   -c5, */
	/* 	       c2, c2, c3, c4,   -c5, */
	/* 	       c3, c2, c2,  0, 2.*c5; */
	/* 	CCCC<< c1, -c2, -c2,    -c3, */
	/* 	       c1, -c3, -c2,    -c2, */
	/* 	       c1, -c2, -c3,    -c2, */
	/* 	        0, -c4, -c4,      0, */
	/* 		0,  c5,  c5, -2.*c5; */
	/* 	MASTER.topLeftCorner(4,4) = BB; */
	/* 	MASTER.bottomRightCorner(5,5) = DD; */
	/* 	MASTER.topRightCorner(4,5) = BBBB; */
	/* 	MASTER.bottomLeftCorner(5,4) = CCCC; */
	/* 	EB.compute(MASTER); */
	/* 	feigs2 = EB.eigenvalues().real(); */
	/* 	/1* cout<<eigs2.transpose()<<endl<<endl; *1/ */
	/* 	if ((abs(feigs(0)-feigs2(0))/abs(feigs(0)) < 1e-2) && (abs(feigs(1)-feigs2(1))/abs(feigs(1)) < 1e-2) && */ 
	/* 			(abs(feigs(2)-feigs2(2))/abs(feigs(2)) < 1e-2) */
	/* 			 && (abs(feigs(7)-feigs2(7))/abs(feigs(7)) < 1e-2) && */
	/* 			(abs(feigs(8)-feigs2(8))/abs(feigs(8)) < 1e-2)){ */
	/* 		Myfile<<setprecision(9)<<sds1<<" "<<pds1<<" "<<pdp1<<" "<<dds1<<" "<<ddp1<<" "<<ddd1<<endl; */
	/* 		cout<<feigs.transpose()<<endl<<feigs2.transpose()<<endl<<endl; */
	/* 		break; */
	/* 	} */
	/* } */
	/* } */
	/* } */
	/* } */
	/* } */
	/* } */

	/* //This block uses the bottom right 5x5 submatrix of t_1 to find the SK potentials of d orbitals accurate */
	/* for (double pot = -0.04255; pot <= -0.04250; pot += 0.000001){ */
	/* for (double pot1 = 0.00030; pot1 <= 0.00033; pot1 += 0.000001){ */
	/* for (double pot2 = 0.01855; pot2 <= 0.01857; pot2 += 0.000001){ */
	/* 	dds1 = pot; */
	/* 	ddp1 = pot1; */
	/* 	ddd1 = pot2; */
	/* 	a1 = (1./9.)*(3.*dds1 + 2.*ddp1 + 4.*ddd1); */
	/* 	a2 = (1./9.)*(3.*dds1 - ddp1 - 2.*ddd1); */
	/* 	a3 = (sqrt(3.)/9.)*(ddp1 - ddd1); */
	/* 	a4 = (1./3.)*(ddp1 - ddd1); */
	/* 	a5 = (1./3.)*(2.*ddp1 + ddd1); */
	/* 	DD<< a1,  a2, a2,   0, -2.*a3, */ 
	/* 	     a2,  a1, a2, -a4,     a3, */
	/* 	     a2,  a2, a1,  a4,     a3, */
	/* 	      0, -a4, a4,  a5,      0, */
	/* 	 -2.*a3,  a3, a3,   0,     a5; */
	/* 	DB.compute(DD); */
	/* 	deigs2 = DB.eigenvalues().real(); */
	/* 	/1* cout<<eigs2.transpose()<<endl<<endl; *1/ */
	/* 	if ((abs(deigs(0)-deigs2(0)) < 7.1e-6) && (abs(deigs(1)-deigs2(1)) < 7.1e-6) && (abs(deigs(2)-deigs2(2)) < 7.1e-6) && (abs(deigs(3)-deigs2(3)) < 7.1e-6) && (abs(deigs(4)-deigs2(4)) < 7.1e-6)){ */
	/* 		Myfile<<setprecision(9)<<dds1<<" "<<ddp1<<" "<<ddd1<<endl; */
	/* 		break; */
	/* 	} */
	/* } */
	/* } */
	/* } */
	/* //This block uses the bottom right 5x5 submatrix of t_1 to find the SK potentials of d orbitals */
	/* for (double pot = -0.04255; pot <= -0.04250; pot += 0.000001){ */
	/* for (double pot1 = 0.01855; pot1 <= 0.01859; pot1 += 0.000001){ */
	/* for (double pot2 = 0.00030; pot2 <= 0.00032; pot2 += 0.000001){ */
	/* 	dds1 = pot; */
	/* 	ddp1 = pot1; */
	/* 	ddd1 = pot2; */
	/* 	a1 = (1./9.)*(3.*dds1 + 2.*ddp1 + 4.*ddd1); */
	/* 	a2 = (1./9.)*(3.*dds1 - ddp1 - 2.*ddd1); */
	/* 	a3 = (sqrt(3.)/9.)*(ddp1 - ddd1); */
	/* 	a4 = (1./3.)*(ddp1 - ddd1); */
	/* 	a5 = (1./3.)*(2.*ddp1 + ddd1); */
	/* 	DD<< a1,  a2, a2,   0, -2.*a3, */ 
	/* 	     a2,  a1, a2, -a4,     a3, */
	/* 	     a2,  a2, a1,  a4,     a3, */
	/* 	      0, -a4, a4,  a5,      0, */
	/* 	 -2.*a3,  a3, a3,   0,     a5; */
	/* 	DB.compute(DD); */
	/* 	deigs2 = DB.eigenvalues().real(); */
	/* 	/1* cout<<eigs2.transpose()<<endl<<endl; *1/ */
	/* 	if ((abs(deigs(0)-deigs2(0)) < 7.1e-6) && (abs(deigs(1)-deigs2(1)) < 7.1e-6) && (abs(deigs(2)-deigs2(2)) < 7.1e-6) && (abs(deigs(3)-deigs2(3)) < 7.1e-6) && (abs(deigs(4)-deigs2(4)) < 7.1e-6)){ */
	/* 		Myfile<<setprecision(9)<<dds1<<" "<<ddp1<<" "<<ddd1<<endl; */
	/* 		break; */
	/* 	} */
	/* } */
	/* } */
	/* } */

	/* //This block uses the top left 4x4 submatrix of t_1 to find the SK potentials of p, s orbitals and does so correctly */
	/* for (double pot = 0.26889233; pot <= 0.26894123; pot += 0.00000001){ */
	/* /1* for (double pot1 = -0.0185935; pot1 <= -0.0185932; pot1 += 0.00000001){ *1/ */
	/* for (double pot2 = -0.16919453; pot2 <= -0.16916656; pot2 += 0.00000001){ */
	/* 	pps1 = pot; */
	/* 	/1* ppp1 = pot1; *1/ */
	/* 	sps1 = pot2; */
	/* 	BB<<sss1   , r3*sps1                 , r3*sps1                 , r3*sps1                 , */
	/* 	   -r3*sps1, (1./3.)*(pps1 + 2.*ppp1), (1./3.)*(pps1 - ppp1)   , (1./3.)*(pps1 - ppp1)   , */
	/* 	   -r3*sps1, (1./3.)*(pps1 - ppp1)   , (1./3.)*(pps1 + 2.*ppp1), (1./3.)*(pps1 - ppp1)   , */
	/* 	   -r3*sps1, (1./3.)*(pps1 - ppp1)   , (1./3.)*(pps1 - ppp1)   , (1./3.)*(pps1 + 2.*ppp1); */
	/* 	CB.compute(BB); */
	/* 	eigs2 = CB.eigenvalues().real(); */
	/* 	/1* cout<<eigs2.transpose()<<endl<<endl; *1/ */
	/* 	if ((abs(eigs(0)-eigs2(0)) < 1.2349e-5) && (abs(eigs(1)-eigs2(1)) < 1.2349e-5) && (abs(eigs(2)-eigs2(2)) < 1.2349e-5) && (abs(eigs(3)-eigs2(3)) < 1.2349e-5)){ */
	/* 		Myfile<<setprecision(9)<<pps1<<" "<<ppp1<<" "<<sps1<<endl; */
	/* 		break; */
	/* 	} */
	/* } */
	/* } */

	//produces incorrect potentials
	/* for (double pot0 = -0.11810; pot0 <= -0.11805; pot0 += 0.00001){ */
	/* for (double pot = -0.16908; pot <= -0.16900; pot += 0.00001){ */
	/* for (double pot1 = 0.07498; pot1 <= 0.07500; pot1 += 0.00001){ */
	/* for (double pot2 = 0.07834; pot2 <= 0.07836; pot2 += 0.00001){ */
	/* 	sss1 = pot0; */
	/* 	sps1 = pot; */
	/* 	pps1 = pot1; */
	/* 	ppp1 = pot2; */
	/* 	BB<<s+5.65685*sss1, 0, 0, (5.65865*sqrt(3.)/3.)*i*sps1, */
	/* 		0, p+(5.65865/3.)*(pps1+2.*ppp1), 0, 0, */
	/* 		0, 0, p+(5.65865/3.)*(pps1+2.*ppp1), (4./3.)*(ppp1-pps1), */
	/* 		-(5.65865*sqrt(3.)/3.)*i*sps1, 0, (4./3.)*(ppp1-pps1), p+(5.65865/3.)*(pps1+2.*ppp1); */
	/* 	CB.compute(BB); */
	/* 	eigs2 = CB.eigenvalues().real(); */
	/* 	/1* cout<<eigs2.transpose()<<endl<<endl; *1/ */
	/* 	if ((abs(eigs(0)-eigs2(0)) < 5.77e-5) && (abs(eigs(1)-eigs2(1)) < 5.77e-5) && (abs(eigs(2)-eigs2(2)) < 5.77e-5) && (abs(eigs(3)-eigs2(3)) < 5.77e-5)){ */
	/* 		Myfile<<sss1<<" "<<sps1<<" "<<pps1<<" "<<ppp1<<endl; */
	/* 	} */
	/* } */
	/* }}} */

	//produces incorrect potentials
	/* for (double pot = -0.2019; pot < -0.2009; pot += 0.00001){ */
	/* for (double pot1 = 0.0067; pot1 < 0.0078; pot1 += 0.00001){ */
	/* for (double pot2 = 0.1117; pot2 < 0.1127; pot2 += 0.00001){ */
	/* 	sps1 = pot; */
	/* 	pps1 = pot1; */
	/* 	ppp1 = pot2; */
	/* 	BB<<s+4.*sss1, 0, (4.*sqrt(3.)/3.)*i*sps1, (4.*sqrt(3.)/3.)*i*sps1, */
	/* 		0, p+(4./3.)*(pps1+2.*ppp1), 0, 0, */
	/* 		-(4.*sqrt(3.)/3.)*i*sps1, 0, p+(4./3.)*(pps1+2.*ppp1), (4./3.)*(ppp1-pps1), */
	/* 		-(4.*sqrt(3.)/3.)*i*sps1, 0, (4./3.)*(ppp1-pps1), p+(4./3.)*(pps1+2.*ppp1); */
	/* 	CB.compute(BB); */
	/* 	eigs2 = CB.eigenvalues().real(); */
	/* 	/1* cout<<eigs2.transpose()<<endl<<endl; *1/ */
	/* 	if ((abs(eigs(0)-eigs2(0)) < 1.5e-5) && (abs(eigs(1)-eigs2(1)) < 1.5e-5) && (abs(eigs(2)-eigs2(2)) < 1.5e-5) && (abs(eigs(3)-eigs2(3)) < 1.5e-5)){ */
	/* 		Myfile<<sps1<<" "<<pps1<<" "<<ppp1<<endl<<endl; */
	/* 		break; */
	/* 	} */
	/* } */
	/* }} */

	/* //this block is used to obtain potentials by matching eigenvalues */
	/* for (double pot = -0.062; pot < -0.0618; pot += 0.00001){ */
	/* for (double pot1 = 0.1633; pot1 < 0.1635; pot1 += 0.00001){ */
	/* for (double pot2 = 0.0305; pot2 < 0.0307; pot2 += 0.00001){ */
	/* 	sps2 = pot; */
	/* 	pps2 = pot1; */
	/* 	ppp2 = pot2; */
	/* 	BB<<s+2.*sss2, 0, 2.*i*sps2, 2.*i*sps2, */
	/* 		0, p+2*pps2, 0, 0, */
	/* 		-2.*i*sps2, 0, p+2*ppp2, 0, */
	/* 		-2.*i*sps2, 0, 0, p+2*ppp2; */
	/* 	CB.compute(BB); */
	/* 	eigs2 = CB.eigenvalues().real(); */
	/* 	/1* cout<<eigs2.transpose()<<endl<<endl; *1/ */
	/* 	if ((abs(eigs(0)-eigs2(0)) < 1.1e-5) && (abs(eigs(1)-eigs2(1)) < 1.1e-5) && (abs(eigs(2)-eigs2(2)) < 1.1e-5) && (abs(eigs(3)-eigs2(3)) < 1.1e-5)){ */
	/* 		Myfile<<sps2<<" "<<pps2<<" "<<ppp2<<endl<<endl; */
	/* 		break; */
	/* 	} */
	/* } */
	/* }} */
	/* cout<<sps2<<endl<<endl; */
	
	Matrix<double,10,1>nn,nnn;
	nn<<sss1, sps1, pps1, ppp1, sds1, pds1, pdp1, dds1, ddp1, ddd1;
	nnn<<sss2, sps2, pps2, ppp2, sds2, pds2, pdp2, dds2, ddp2, ddd2;

	/* t_1 = TB(0,1,0,9,d_1,nn,nnn); */
	/* t_2 = TB(0,1,0,9,d_2,nn,nnn); */
	/* t_3 = TB(0,1,0,9,d_3,nn,nnn); */
	/* t_4 = TB(0,1,0,9,d_4,nn,nnn); */
	/* t_5 = TB(0,1,0,9,d_5,nn,nnn); */
	/* t_6 = TB(0,1,0,9,d_6,nn,nnn); */
	/* t_7 = TB(0,1,0,9,d_7,nn,nnn); */
	/* t_8 = TB(0,1,0,9,d_8,nn,nnn); */
	/* t_13 = TB(0,1,1,9,d_13,nn,nnn); */
	/* t_14 = TB(0,1,1,9,d_14,nn,nnn); */
	/* t_15 = TB(0,1,1,9,d_15,nn,nnn); */
	/* t_16 = TB(0,1,1,9,d_16,nn,nnn); */
	/* t_17 = TB(0,1,1,9,d_17,nn,nnn); */
	/* t_18 = TB(0,1,1,9,d_18,nn,nnn); */

	/* cout<<t_13.topLeftCorner(4,4)<<endl<<endl; */

	/* Myfile<<"P X Y"<<endl; */
	Myfile<<"X Y"<<endl;

	Matrix<dcomp, 9, 9> E;
	Matrix<dcomp, 4, 4> E_small;
	Matrix<dcomp, 5, 5> E_mid;

	/* Matrix<dcomp, 18, 18> E; */
	/* Matrix<dcomp, 9, 9> u_11, u_12; */
	/* Matrix<dcomp, 18, 18> U, T; */
	/* Matrix<complex<double>, 9, 9> zero = Matrix<complex<double>, 9, 9>::Zero(); */


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
		E = lambda
		       	/* + t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K)) */
			/* + t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K)) */
				/* + t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K)) */
				+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
				+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
				+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
				;

		/* E_mid = E.bottomRightCorner(5,5); */
		/* E_small = E.topLeftCorner(4,4); */

		/* if (k == 5) */
		/* 	cout<<E<<endl<<endl; */
		/* if (k == 15) */
		/* 	cout<<E<<endl<<endl; */
		/* if (k == 50) */
		/* 	cout<<E<<endl<<endl; */

		SelfAdjointEigenSolver<Matrix<dcomp, 9, 9>> es;
		/* SelfAdjointEigenSolver<Matrix<dcomp, 4, 4>> es; */
		/* SelfAdjointEigenSolver<Matrix<dcomp, 5, 5>> es; */
		/* es.compute(E_small); */
		/* es.compute(E_mid); */
		es.compute(E);
		Matrix<double, 9, 1> O;
		/* Matrix<double, 4, 1> O; */
		/* Matrix<double, 5, 1> O; */
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
	Myfile2.close();
	/* Myfile3.close(); */
	Myfile4.close();
	return 0;
}
