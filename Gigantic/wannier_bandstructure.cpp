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
	string Mydata;//, Mydata2, Mydata3, Mydata4;
	getline(cin, Mydata);
	ofstream Myfile;//, Myfile2, Myfile3, Myfile4;	
	/* Mydata2 = Mydata; */
	/* Mydata3 = Mydata; */
	/* Mydata4 = Mydata; */
	Mydata += ".txt";
	/* Mydata2 += ".dat"; */
	/* Mydata3 += "_catch.dat"; */
	/* Mydata4 += "_eigs.dat"; */
	Myfile.open( Mydata.c_str(),ios::trunc );
	/* Myfile2.open( Mydata2.c_str(),ios::trunc ); */
	/* Myfile3.open( Mydata3.c_str(),ios::trunc ); */
	/* Myfile4.open( Mydata4.c_str(),ios::trunc ); */

	cout<<endl;

	int ispin = -1;

	Vector3d d_0, d_1, d_2, d_3, d_4, d_5, d_6, d_7, d_8;
	Vector3d d_13, d_14, d_15, d_16;
	Vector3d d_17, d_18;
	Vector3d d_19, d_20, d_21, d_22, d_23, d_24, d_25, d_26,
		 d_27, d_28, d_29, d_30;
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
	K << 0, M_PI/b, M_PI/b;
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
  
	pps1 = 2.*xy1 + xx1;//obtained from Papa formula in intro */
	ppp1 = xx1 - xy1;//obtained from Papa formula in intro */
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
	cout<<"sss = "<<sss1<<endl<<"pps = "<<pps1<<endl<<"ppp = "<<ppp1<<endl;
	cout<<"ddd = "<<ddd1<<endl<<"ddp = "<<ddp1<<endl<<"dds = "<<dds1<<endl;
	cout<<"}"<<endl;

	/* cout<<"dds = "<<dds1<<endl; */
	/* cout<<"ddp = "<<ddp1<<endl; */
	/* cout<<"ddd = "<<ddd1<<endl; */

	//this block obtained from matching eigenvalues, see block below accurate!
	/* ppp1 = -0.018593; */
	/* pps1 = 0.2689; */
	/* sps1 = -0.16918; */

	/* cout<<"ppp = "<<ppp1<<endl; */
	/* cout<<"pps = "<<pps1<<endl; */
	/* cout<<"sps = "<<sps1<<endl; */

	ppp2 = (1./16.)*real(E_gamma(1,1) + E_H(1,1) - 2.*E_N(1,1)); //obtained from Papa formula in intro
	/* pps2 = (1./16.)*real(E_gamma(1,1) + E_H(1,1) + 2.*(-E_N(2,2) - E_N(3,3) + E_N(1,1)));//obtained from Papa formula in intro */
	cout<<endl;
	//this block obtained from matching eigenvalues, see block below accurate!
	/* ppp2 = 0.03060; */
	/* pps2 = 0.16341; */ 
	/* sps2 = 0.06189; */ 

	cout<<"Papa 2nd neighbours"<<endl;
	/* cout<<"sss = "<<sss2<<endl; */
	/* cout<<"pps = "<<pps2<<endl; */
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

	/* K << M_PI/b, 0, 0; */
	/* K << 0, M_PI/b, M_PI/b; */
	K << 3.*M_PI/(2.*b), M_PI/(2.*b), M_PI/(2.*b);
	Matrix<dcomp, 9, 9> E_AA;
	/* cout<<exp(i*d_17.dot(K))<<" "<<exp(i*d_18.dot(K))<<endl; */
	/* E_AA = t_1; */
	E_AA = lambda 
		+ t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
		+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
			+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K))
			/* + t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K)) */
			/* + t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K)) */
			/* + t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K)) */
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
	sort(eigs.data(), eigs.data()+eigs.size());
	/* cout<<CA.eigenvalues().transpose()<<endl<<endl; */
	/* cout<<CA.eigenvalues().sum()<<endl<<endl; */

	Matrix<dcomp, 5, 5> CC, DD;
	CC = E_AA.bottomRightCorner(5,5);
	ComplexEigenSolver<Matrix<dcomp, 5, 5>> DA, DB;
	DA.compute(CC);
	Matrix<double, 5, 1> deigs, deigs2;
	deigs = DA.eigenvalues().real();

	double a1, a2, a3, a4, a5;
	dcomp b1, b2, b3, b4, b5;
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

	/* //this computed at K=(0,pi,pi)*/
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

	/* sds1 = -0.07770; */
	/* //this computed at K=(0,pi,pi) */
	/* //taken sps to be -ve */
	/* BB<< s, 0,                   0,                   0, */
	/*      0, p,                   0,                   0, */
	/*      0, 0,                   p, (8./3.)*(ppp1-pps1), */
	/*      0, 0, (8./3.)*(ppp1-pps1),                   p; */
	/* for (double pot = -0.07771; pot <= -0.07768; pot += 0.000001){ */
	/* 	cout<<pot<<endl<<endl; */
	/* for (double pot1 = -0.09161; pot1 <= -0.09159; pot1 += 0.000001){ */
	/* for (double pot2 = -0.03316; pot2 <= -0.03307; pot2 += 0.000001){ */
	/* for (double pot3 = -0.04985; pot3 <= -0.04978; pot3 += 0.000001){ */
	/* 	sds1 = pot; */
	/* 	dds1 = pot1; */
	/* 	ddp1 = pot2; */
	/* 	ddd1 = pot3; */
	/* 	a1 = d1; */
	/* 	a2 = -(8./3.)*dds1 + (8./9.)*ddp1 + (16./9.)*ddd1; */
	/* 	a3 = d2; */
	/* 	a4 = (8./3.)*(ddp1 - ddd1); */
	/* 	a5 = (8./9.)*sqrt(3.)*(ddd1 - ddp1); */
	/* 	DD<< a1,  0, a2,  0,  0, */ 
	/* 	      0, a1,  0, a4, a5, */
	/* 	     a2,  0, a1,  0,  0, */
	/* 	      0, a4,  0, a3,  0, */
	/* 	      0, a5,  0,  0, a3; */
	/* 	c1 = -(8./3.)*sqrt(3.)*sds1; */
	/* 	BBBB<<  0, c1, 0, 0, 0, */
	/* 	        0,  0, 0, 0, 0, */
	/* 	        0,  0, 0, 0, 0, */
	/* 	        0,  0, 0, 0, 0; */
	/* 	CCCC<<  0, 0, 0, 0, */
	/* 	       c1, 0, 0, 0, */
	/* 		0, 0, 0, 0, */
	/* 	        0, 0, 0, 0, */
	/* 	        0, 0, 0, 0; */
	/* 	MASTER.topLeftCorner(4,4) = BB; */
	/* 	MASTER.bottomRightCorner(5,5) = DD; */
	/* 	MASTER.topRightCorner(4,5) = BBBB; */
	/* 	MASTER.bottomLeftCorner(5,4) = CCCC; */
	/* 	EB.compute(MASTER); */
	/* 	feigs2 = EB.eigenvalues().real(); */
	/* 	sort(feigs2.data(), feigs2.data()+feigs2.size()); */
	/* 	/1* cout<<eigs2.transpose()<<endl<<endl; *1/ */
	/* 	if ((abs(feigs(0)-feigs2(0)) < 4.9813e-5) && (abs(feigs(1)-feigs2(1)) < 4.9813e-5) */
	/* 		       	&& (abs(feigs(2)-feigs2(2)) < 4.9813e-5) && (abs(feigs(3)-feigs2(3)) < 4.9813e-5) */
	/* 		       	&& (abs(feigs(4)-feigs2(4)) < 4.9813e-5) && (abs(feigs(5)-feigs2(5)) < 4.9813e-5) */
	/* 		       	&& (abs(feigs(6)-feigs2(6)) < 4.9813e-5) && (abs(feigs(7)-feigs2(7)) < 4.9813e-5) */
	/* 		       	&& (abs(feigs(8)-feigs2(8)) < 4.9813e-5) && (abs(feigs.sum()-feigs2.sum())<3e-5)){ */
	/* 		Myfile<<setprecision(9)<<"sds1 = "<<sds1<<"; " */
	/* 			<<"dds1 = "<<dds1<<"; "<<"ddp1 = "<<ddp1<<"; "<<"ddd1 = "<<ddd1<<";"<<endl; */
	/* 		Myfile2<<setprecision(9)<<sds1<<" "<<dds1<<" "<<ddp1<<" "<<ddd1<<" "<<abs(feigs.sum()-feigs2.sum())<<endl; */
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

	/* double start, end; */
	/* sds1 = -0.076; pds1 = 0.108; pdp1 = 0.019; dds1 = -0.026; ddp1 = 0.02; ddd1 = 0.011; */
	/* sds1 = -0.072; pds1 = 0.108; pdp1 = -0.02; dds1 = -0.036; ddp1 = 0.021; ddd1 = -0.009; */
	/* cout<<"These are my first nn SK potentials {"<<endl; */
	/* cout<<"ddd = "<<ddd1<<endl<<"ddp = "<<ddp1<<endl<<"dds = "<<dds1<<endl; */
	/* cout<<"sds = "<<sds1<<endl<<"pds = "<<pds1<<endl<<"pdp = "<<pdp1<<endl; */
	/* cout<<"}"<<endl; */

	/* pds1 = -0.07884; pdp1 = -0.03462; */
	/* pds1 = 0.07884; pdp1 = 0.03462; */
	pds1 = -0.11882; pdp1 = 0.03462; 
	/* pds1 = 0.11882; pdp1 = -0.03462; */ 
	/* //K << M_PI/b, 0, 0; */
	/* b1 = (8./3.)*sqrt(3.)*i*sps1; */
	/* BB<< s, b1, 0, 0, */
	/*    -b1,  p, 0, 0, */
	/*      0,  0, p, 0, */
	/*      0,  0, 0, p; */
	/* DD<< d1,  0,  0,  0,  0, */
	/*       0, d1,  0,  0,  0, */
	/*       0,  0, d1,  0,  0, */
	/*       0,  0,  0, d2,  0, */
	/*       0,  0,  0,  0, d2; */
	/* for (double pot1 = -0.118822; pot1 <= -0.118812; pot1 += 0.0000001){ */
	/* 	cout<<pot1<<endl; */
	/* for (double pot2 = 0.034620; pot2 <= 0.034624; pot2 += 0.0000001){ */
	/* 	pds1 = pot1; */
	/* 	pdp1 = pot2; */
	/* 	c1 = i*(8./9.)*sqrt(3.)*(sqrt(3.)*pds1 + pdp1); */
	/* 	c2 = i*(8./3.)*sqrt(3.)*pdp1; */
	/* 	c3 = i*(8./3.)*pdp1; */
	/* 	BBBB<<  0,  0,  0,  0,  0, */
	/* 	        0,  0,  0, c2,-c3, */
	/* 	       c1,  0,  0,  0,  0, */
	/* 	        0,  0, c1,  0,  0; */
	/* 	CCCC<<  0,   0, -c1,   0, */
	/* 	        0,   0,   0,   0, */
	/* 	        0,   0,   0, -c1, */
	/* 	        0, -c2,   0,   0, */
	/* 		0,  c3,   0,   0; */
	/* 	MASTER.topLeftCorner(4,4) = BB; */
	/* 	MASTER.bottomRightCorner(5,5) = DD; */
	/* 	MASTER.topRightCorner(4,5) = BBBB; */
	/* 	MASTER.bottomLeftCorner(5,4) = CCCC; */
	/* 	EB.compute(MASTER); */
	/* 	feigs2 = EB.eigenvalues().real(); */
	/* 	sort(feigs2.data(), feigs2.data()+feigs2.size()); */
	/* 	/1* cout<<eigs2.transpose()<<endl<<endl; *1/ */
	/* 	if ((abs(feigs(0)-feigs2(0)) < 8e-6) && (abs(feigs(1)-feigs2(1)) < 8e-6) */
	/* 		       	&& (abs(feigs(2)-feigs2(2)) < 8e-6) && (abs(feigs(3)-feigs2(3)) < 8e-6) */
	/* 		       	&& (abs(feigs(4)-feigs2(4)) < 8e-6) && (abs(feigs(5)-feigs2(5)) < 8e-6) */
	/* 		       	&& (abs(feigs(6)-feigs2(6)) < 8e-6) && (abs(feigs(7)-feigs2(7)) < 8e-6) */
	/* 		       	&& (abs(feigs(8)-feigs2(8)) < 8e-6) && (abs(feigs.sum()-feigs2.sum())<3e-4)){ */
	/* 		Myfile<<setprecision(9)<<"pds1 = "<<pds1<<"; pdp1 = "<<pdp1<<"; "<<endl; */
	/* 		Myfile2<<setprecision(9)<<pds1<<" "<<pdp1<<" "<<abs(feigs.sum()-feigs2.sum())<<endl; */
	/* 		Myfile4<<feigs2.transpose()<<endl; */
	/* 		cout<<feigs.transpose()<<endl<<feigs2.transpose()<<endl<<endl; */
	/* 	} */
	/* } */
	/* } */

	pps1 = 0.26892; ppp1 = -0.01859; sps1 = 0.16918;
	/* //K << 3.*M_PI/(2.*b), M_PI/(2.*b), M_PI/(2.*b); */
	/* for (double pot = 0.2682; pot <= 0.277; pot += 0.0001){ */
	/* 	cout<<pot<<endl; */
	/* for (double pot1 = -0.0192; pot1 <= -0.012; pot1 += 0.0001){ */
	/* for (double pot2 = -0.1696; pot2 <= 0.175; pot2 += 0.0001){ */
	/* 	pps1 = pot; */
	/* 	ppp1 = pot1; */
	/* 	sps1 = pot2; */
	/* 	b1 = s - 2.*M_SQRT2*sss1; */
	/* 	b2 = i*2.*M_SQRT2*r3*sps1; */
	/* 	b3 = p - (2.*M_SQRT2/3.)*(pps1 + 2.*ppp1); */
	/* 	b4 = -(2.*M_SQRT2/3.)*(pps1 - ppp1); */
	/* 	BB<< b1, b2, -b2, -b2, */
	/* 	    -b2, b3,  b4,  b4, */
	/* 	     b2, b4,  b3, -b4, */
	/* 	     b2, b4, -b4,  b3; */
	/* 	CB.compute(BB); */
	/* 	eigs2 = CB.eigenvalues().real(); */
	/* 	sort(eigs2.data(), eigs2.data()+eigs2.size()); */
	/* 	/1* cout<<eigs2.transpose()<<endl<<endl; *1/ */
	/* 	if ((abs(eigs(0)-eigs2(0)) < 1e-3) && (abs(eigs(1)-eigs2(1)) < 1e-3) && (abs(eigs(2)-eigs2(2)) < 1e-3) && (abs(eigs(3)-eigs2(3)) < 1e-3)){ */
	/* 		Myfile<<setprecision(9)<<"pps1 = "<<pps1<<"; "<<"ppp1 = "<<ppp1<<"; "<<"sps1 = "<<sps1<<";"<<endl; */
	/* 		Myfile2<<setprecision(9)<<pps1<<" "<<ppp1<<" "<<sps1<<" "<<abs(eigs.sum()-eigs2.sum())<<endl; */
	/* 		Myfile4<<eigs2.transpose()<<endl; */
	/* 		cout<<eigs.transpose()<<endl<<eigs2.transpose()<<endl<<endl; */
	/* 	} */
	/* } */
	/* } */
	/* } */

	/* pds1 = 0.11758; pdp1 = -0.01924; dds1 = -0.03998; ddp1 = 0.02639; ddd1 = -0.00877; */

	/* sds1 = -0.08552; dds1 = -0.04549; ddp1 = 0.02456; ddd1 = -0.00497; */
	sds1 = -0.07158; dds1 = -0.04897; ddp1 = 0.02434; ddd1 = -0.00178;

	/* //K << 3.*M_PI/(2.*b), M_PI/(2.*b), M_PI/(2.*b); */
	/* b1 = s - 2.*M_SQRT2*sss1; */
	/* b2 = i*2.*M_SQRT2*r3*sps1; */
	/* b3 = p - (2.*M_SQRT2/3.)*(pps1 + 2.*ppp1); */
	/* b4 = -(2.*M_SQRT2/3.)*(pps1 - ppp1); */
	/* BB<< b1, b2, -b2, -b2, */
	/*     -b2, b3,  b4,  b4, */
	/*      b2, b4,  b3, -b4, */
	/*      b2, b4, -b4,  b3; */
	/* double f1,f2,f3,f4,f5,f6; */
	/* f1 = 2.*M_SQRT2; */
	/* f2 = f1/3.; */
	/* f3 = f2*2.; */
	/* f4 = f3/3.; */
	/* f5 = f4*2.; */
	/* f6 = f4/2.; */
	/* /1* //this in range but not pursued{ *1/ */
	/* 0.079	-0.056	-0.032	-0.013	*/
	/* /1* 0.105	0.002	0.028	0.035 *1/ */
	/* /1* } *1/ */
	/* for (double pot = 0.0854; pot <= 0.0863; pot += 0.0001){ */
	/* 	cout<<pot<<endl; */
	/* for (double pot1 = -0.0468; pot1 <= -0.0449; pot1 += 0.0001){ */
	/* for (double pot2 = 0.0236; pot2 <= 0.0259; pot2 += 0.0001){ */
	/* for (double pot3 = -0.0054; pot3 <= -0.0037; pot3 += 0.0001){ */
	/* /1* for (double pot = -0.07169; pot <= -0.07155; pot += 0.00001){ *1/ */
	/* /1* 	cout<<pot<<endl; *1/ */
	/* /1* for (double pot1 = -0.04899; pot1 <= -0.04892; pot1 += 0.00001){ *1/ */
	/* /1* for (double pot2 = 0.02433; pot2 <= 0.02440; pot2 += 0.00001){ *1/ */
	/* /1* for (double pot3 = -0.00182; pot3 <= -0.00176; pot3 += 0.00001){ *1/ */
	/* /1* { This block yields potentials in the range *1/ */
	/* /1* for (double pot = 0.06; pot <= 0.095; pot += 0.005){ *1/ */
	/* /1* 	cout<<pot<<endl; *1/ */
	/* /1* for (double pot1 = -0.065; pot1 <= 0.015; pot1 += 0.005){ *1/ */
	/* /1* for (double pot2 = -0.035; pot2 <= 0.03; pot2 += 0.005){ *1/ */
	/* /1* for (double pot3 = -0.015; pot3 <= 0.04; pot3 += 0.005){ *1/ */
	/* /1* but decided not to pursue this sign} *1/ */
	/* /1* for (double pot = -0.08558; pot <= -0.08546; pot += 0.00001){ *1/ */
	/* /1* 	cout<<pot<<endl; *1/ */
	/* /1* for (double pot1 = -0.04553; pot1 <= -0.04546; pot1 += 0.00001){ *1/ */
	/* /1* for (double pot2 = 0.02453; pot2 <= 0.02462; pot2 += 0.00001){ *1/ */
	/* /1* for (double pot3 = -0.00500; pot3 <= -0.00495; pot3 += 0.00001){ *1/ */
	/* /1* for (double pot1 = -0.0604; pot1 <= -0.0592; pot1 += 0.00001){ *1/ */
	/* /1* 	cout<<pot1<<endl; *1/ */
	/* /1* for (double pot2 = 0.0162; pot2 <= 0.0194; pot2 += 0.00001){ *1/ */
	/* /1* for (double pot3 = 0.0082; pot3 <= 0.0116; pot3 += 0.00001){ *1/ */
	/* 	sds1 = pot; */
	/* 	dds1 = pot1; */
	/* 	ddp1 = pot2; */
	/* 	ddd1 = pot3; */
	/* 	a1 = d1 - f6*(3.*dds1 + 2.*ddp1 + 4.*ddd1); */
	/* 	a2 = f6*(3.*dds1 - ddp1 - 2.*ddd1); */
	/* 	a3 = f6*sqrt(3.)*(ddp1 - ddd1); */
	/* 	a4 = f2*(ddp1 - ddd1); */
	/* 	a5 = d2 - f2*(2.*ddp1 + ddd1); */
	/* 	DD<< a1, -a2,  a2,   0, 2.*a3, */ 
	/* 	    -a2,  a1, -a2, -a4,    a3, */
	/* 	     a2, -a2,  a1, -a4,   -a3, */
	/* 	      0, -a4, -a4,  a5,     0, */
	/* 	  2.*a3,  a3, -a3,   0,    a5; */
	/* 	c1 = f2*sqrt(3.)*sds1; */
	/* 	c2 = i*f6*sqrt(3.)*(sqrt(3.)*pds1 + pdp1); */
	/* 	c3 = i*f6*sqrt(3.)*(sqrt(3.)*pds1 - 2.*pdp1); */
	/* 	c4 = i*f2*sqrt(3.)*pdp1; */
	/* 	c5 = i*f2*pdp1; */
	/* 	BBBB<<-c1, c1,-c1,  0,     0, */
	/* 	      -c2,-c3,-c2, c4,   -c5, */
	/* 	       c2,-c2,-c3,-c4,    c5, */
	/* 	      -c3,-c2, c2,  0,-2.*c5; */
	/* 	CCCC<<-c1,  c2, -c2,    c3, */
	/* 	       c1,  c3,  c2,    c2, */
	/* 	      -c1,  c2,  c3,   -c2, */
	/* 	        0, -c4,  c4,     0, */
	/* 		0,  c5, -c5, 2.*c5; */
	/* 	MASTER.topLeftCorner(4,4) = BB; */
	/* 	MASTER.bottomRightCorner(5,5) = DD; */
	/* 	MASTER.topRightCorner(4,5) = BBBB; */
	/* 	MASTER.bottomLeftCorner(5,4) = CCCC; */
	/* 	EB.compute(MASTER); */
	/* 	feigs2 = EB.eigenvalues().real(); */
	/* 	sort(feigs2.data(), feigs2.data()+feigs2.size()); */
	/* 	/1* cout<<eigs2.transpose()<<endl<<endl; *1/ */
	/* 	if ((abs(feigs(0)-feigs2(0)) < 1.23e-2) && (abs(feigs(1)-feigs2(1)) < 1.23e-2) */
	/* 		       	&& (abs(feigs(2)-feigs2(2)) < 1.23e-2) && (abs(feigs(3)-feigs2(3)) < 1.23e-2) */
	/* 		       	&& (abs(feigs(4)-feigs2(4)) < 1.23e-2) && (abs(feigs(5)-feigs2(5)) < 1.23e-2) */
	/* 		       	&& (abs(feigs(6)-feigs2(6)) < 1.23e-2) && (abs(feigs(7)-feigs2(7)) < 1.23e-2) */
	/* 		       	&& (abs(feigs(8)-feigs2(8)) < 1.23e-2) && (abs(feigs.sum()-feigs2.sum())<7e-3)){ */
	/* 		Myfile<<setprecision(9)<<"sds1 = "<<sds1<<"; "<< */
	/* 			"dds1 = "<<dds1<<"; "<<"ddp1 = "<<ddp1<<"; "<<"ddd1 = "<<ddd1<<";"<<endl; */
	/* 		Myfile2<<setprecision(9)<<sds1<<" " */
	/* 			<<dds1<<" "<<ddp1<<" "<<ddd1<<" "<<abs(feigs.sum()-feigs2.sum())<<endl; */
	/* 		Myfile4<<feigs2.transpose()<<endl; */
	/* 		cout<<feigs.transpose()<<endl<<feigs2.transpose()<<endl<<endl; */
	/* 	} */
	/* } */
	/* } */
	/* } */
	/* } */

	//This one may have chosen the wrong sign as d orbitals were incorrect
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
//obtained from conj.cpp
sss1 =  -0.118103; sps1 =   0.169403; pps1 =    0.269011; ppp1 =  -0.0186036;
sds1 = -0.0761861; pds1 =  -0.112363; pdp1 =   0.0305943; dds1 =  -0.0433613; ddp1 =   0.0207795; ddd1 = -0.00149068;
sss2 = -0.0228124; sps2 =  0.0616934; pps2 =    0.164078; ppp2 =   0.0300529;
sds2 = -0.0281175; pds2 = -0.0534488; pdp2 =  -0.0129463; dds2 =  -0.0244981; ddp2 = -0.00288585; ddd2 =  0.00130428;
//obtained from conj.cpp in one run
sss1 = -0.11806; sps1 = 0.169358; pps1 = 0.269424; ppp1 = -0.0188317; sds1 = -0.0756415;
pds1 = -0.110769; pdp1 = 0.0306304; dds1 = -0.0438971; ddp1 = 0.0200678; ddd1 = -0.000538349;
sss2 = -0.023016; sps2 = 0.061915; pps2 = 0.165029; ppp2 = 0.0301988; sds2 = -0.0318296;
pds2 =-0.0567376; pdp2 = -0.0105958; dds2 = -0.024202; ddp2 = -0.00369887; ddd2 = 0.00177699;
double sss3, sps3, pps3, ppp3, sds3, pds3, pdp3, dds3, ddp3, ddd3;
//obtained from conj.cpp keeping 1st and 2nd nn version fixed
sss3 = 0.0171521; sps3 = -0.0243898; pps3 = -0.0500761; ppp3 = 0.015545; sds3 = 0.00872122;
pds3 = 0.00898767; pdp3 = -0.00244515; dds3 = 0.00487065; ddp3 = -0.000943947; ddd3 = -0.000206366;
//obtained from conj.cpp in one run
sss1 = -0.117924; sps1 = 0.169252; pps1 = 0.268776; ppp1 = -0.0187581; sds1 = -0.0757635;
pds1 = -0.110995; pdp1 = 0.0287141; dds1 = -0.0447073; ddp1 = 0.0197826; ddd1 = 5.90525e-05;
sss2 = -0.0220203; sps2 = 0.0614591; pps2 = 0.1639; ppp2 = 0.0290812; sds2 = -0.0326461;
pds2 = -0.0571585; pdp2 = -0.0117408; dds2 = -0.0253075; ddp2 = -0.00318763; ddd2 = 0.0026379;
sss3 = 0.017236; sps3 = -0.0252936; pps3 = -0.0503531; ppp3 = 0.0160445; sds3 = 0.0088867; 
pds3 = 0.0103232; pdp3 = -0.00234456; dds3 = 0.00533451; ddp3 = -0.00119393; ddd3 = -0.000529859;
//obtained from conj.cpp in one run but over 35 k-points
/* sss1 = -0.117542; sps1 = 0.170921; pps1 = 0.269489; ppp1 = -0.0200408; sds1 = -0.0792361; */
/* pds1 = -0.115149; pdp1 = 0.0246214; dds1 = -0.0431099; ddp1 = 0.0194892; ddd1 = -0.000221514; */
/* sss2 = -0.0228803; sps2 = 0.0597177; pps2 = 0.161577; ppp2 = 0.0302199; sds2 = -0.0321991; */
/* pds2 = -0.0566746; pdp2 = -0.0116066; dds2 = -0.0252924; ddp2 = -0.00275971; ddd2 = 0.00362367; */
/* sss3 = 0.018064; sps3 = -0.0265867; pps3 = -0.0500519; ppp3 = 0.0163536; sds3 = 0.00738257; */ 
/* pds3 = 0.0102813; pdp3 = -0.00646097; dds3 = 0.00258293; ddp3 = -0.00116498; ddd3 = -0.000344613; */
//obtained from conj.cpp in one run with the usual 8 k-points, but with smaller step in dxargs
/* sss1 = -0.118215; sps1 = 0.170286; pps1 = 0.267537; ppp1 = -0.0181281; sds1 = -0.0774747; */
/* pds1 = -0.117903; pdp1 = 0.0210175; dds1 = -0.043845; ddp1 = 0.0205712; ddd1 = -0.000717575; */
/* sss2 = -0.0227063; sps2 = 0.05606; pps2 = 0.163461; ppp2 = 0.0288275; sds2 = -0.0335694; */
/* pds2 = -0.0568299; pdp2 = -0.0108152; dds2 = -0.024618; ddp2 = -0.00190861; ddd2 = 0.00202976; */
/* sss3 = 0.017559; sps3 = -0.0249802; pps3 = -0.0502737; ppp3 = 0.0164508; sds3 = 0.00599136; */ 
/* pds3 = 0.0087034; pdp3 = -0.000854381; dds3 = 0.00251543; ddp3 = -0.00209; ddd3 = 0.00113103; */
//obtained from conj.cpp in one run SPIN UP THIS TIME!
/* sss1 = -0.121162; sps1 = 0.173363; pps1 = 0.272863; ppp1 = -0.0178175; sds1 = -0.0822938; */
/* pds1 = -0.120081; pdp1 = 0.0308906; dds1 = -0.0505002; ddp1 = 0.0248364; ddd1 = -0.000290327; */
/* sss2 = -0.0225453; sps2 = 0.0654877; pps2 = 0.172982; ppp2 = 0.032321; sds2 = -0.035475; */
/* pds2 = -0.0673341; pdp2 = -0.0195645; dds2 = -0.0321466; ddp2 = -0.00496235; ddd2 = 0.00295384; */
/* sss3 = 0.0187611; sps3 = -0.0265531; pps3 = -0.0504351; ppp3 = 0.0150564; sds3 = 0.00808336; */ 
/* pds3 = 0.00602272; pdp3 = 0.00130759; dds3 = 0.00818363; ddp3 = -0.00161298; ddd3 = -0.000514812; */

//spin up obtained from conj.cpp in one run with 350 k-points over the same k-point path as bandstructure 
sss1 = -0.121541; sps1 = 0.172875; pps1 = 0.273266; ppp1 = -0.0182001; sds1 = -0.0852609;
pds1 = -0.124371; pdp1 = 0.0233553; dds1 = -0.049139; ddp1 = 0.0262377; ddd1 = -0.00143019;
sss2 = -0.0239523; sps2 = 0.0608069; pps2 = 0.170197; ppp2 = 0.0314321; sds2 = -0.0365085;
pds2 = -0.0703064; pdp2 = -0.0146721; dds2 = -0.032293; ddp2 = -0.00262184; ddd2 = 0.00332746;
sss3 = 0.0190502; sps3 = -0.0285496; pps3 = -0.0500113; ppp3 = 0.0167854; sds3 = 0.00457944; 
pds3 = 0.00712668; pdp3 = -0.00213332; dds3 = 0.00432888; ddp3 = -0.00220879; ddd3 = 0.000162835;
//spin down obtained from conj.cpp in one run with 350 k-points over the same k-point path as bandstructure 
/* sss1 = -0.118269; sps1 = 0.168832; pps1 = 0.267676; ppp1 = -0.0185253; sds1 = -0.0786542; */
/* pds1 = -0.117806; pdp1 = 0.0212018; dds1 = -0.0440309; ddp1 = 0.0211417; ddd1 = -0.000775478; */
/* sss2 = -0.02272; sps2 = 0.0567299; pps2 = 0.162659; ppp2 = 0.0295257; sds2 = -0.0334902; */
/* pds2 = -0.0567011; pdp2 = -0.0127372; dds2 = -0.0255078; ddp2 = -0.00205756; ddd2 = 0.00253642; */
/* sss3 = 0.0175224; sps3 = -0.0266653; pps3 = -0.0506264; ppp3 = 0.016701; sds3 = 0.00493575; */ 
/* pds3 = 0.0078857; pdp3 = -0.00184691; dds3 = 0.00297654; ddp3 = -0.00177099; ddd3 = 0.000399248; */
cout<<lambda<<endl<<endl;

//this is spin down from my conj-deriv.cpp
//1st nearest neighbour SK terms:
/* sss1 = -0.11827695; sps1 = +0.16913680; pps1 = +0.26953060; ppp1 = -0.01898993; sds1 = -0.07623039; */
/* pds1 = -0.11716609; pdp1 = +0.03441365; dds1 = -0.04355404; ddp1 = +0.02002808; ddd1 = -0.00025232; */

/* //2nd nearest neighbour SK terms: */
/* sss2 = -0.02253620; sps2 = +0.06249661; pps2 = +0.16390065; ppp2 = +0.03046760; sds2 = -0.02980746; */
/* pds2 = -0.05287675; pdp2 = -0.01406475; dds2 = -0.02403076; ddp2 = -0.00343147; ddd2 = +0.00115502; */

/* //3rd nearest neighbour SK terms: */
/* sss3 = +0.01755674; sps3 = -0.03033931; pps3 = -0.05395566; ppp3 = +0.01701694; sds3 = +0.00914059; */ 
/* pds3 = +0.01309417; pdp3 = -0.00075238; dds3 = +0.00384095; ddp3 = -0.00072702; ddd3 = +0.00031952; */

//this is spin down from my with loads of k-points conj-deriv.cpp
//1st nearest neighbour SK terms:
/* sss1 = -0.11857130; sps1 = +0.16926845; pps1 = +0.26868093; ppp1 = -0.01867662; sds1 = -0.07632015; */
/* pds1 = -0.11385718; pdp1 = +0.03100820; dds1 = -0.04270871; ddp1 = +0.01948932; ddd1 = +0.00025397; */

/* //2nd nearest neighbour SK terms: */
/* sss2 = -0.02273373; sps2 = +0.06211954; pps2 = +0.16409752; ppp2 = +0.03109120; sds2 = -0.03194757; */
/* pds2 = -0.05498245; pdp2 = -0.01386365; dds2 = -0.02345247; ddp2 = -0.00426026; ddd2 = +0.00109387; */

/* //3rd nearest neighbour SK terms: */
/* sss3 = +0.01769373; sps3 = -0.02979613; pps3 = -0.05439512; ppp3 = +0.01698135; sds3 = +0.00860351; */ 
/* pds3 = +0.01157403; pdp3 = -0.00217653; dds3 = +0.00474107; ddp3 = -0.00059868; ddd3 = +0.00017544; */

//1st nearest neighbour SK terms:
/* sss1 = -0.11815786; sps1 = +0.16934335; pps1 = +0.26867396; ppp1 = -0.01861928; sds1 = -0.07624969; */
/* pds1 = -0.11646945; pdp1 = +0.03343823; dds1 = -0.04301012; ddp1 = +0.01981637; ddd1 = -0.00023724; */

/* //2nd nearest neighbour SK terms: */
/* sss2 = -0.02276512; sps2 = +0.06188834; pps2 = +0.16346027; ppp2 = +0.03095968; sds2 = -0.03021270; */
/* pds2 = -0.05194449; pdp2 = -0.01348462; dds2 = -0.02370670; ddp2 = -0.00415149; ddd2 = +0.00185349; */

/* //3rd nearest neighbour SK terms: */
/* sss3 = +0.01743744; sps3 = -0.02983831; pps3 = -0.05287215; ppp3 = +0.01648297; sds3 = +0.00893271; */ 
/* pds3 = +0.01125425; pdp3 = -0.00114839; dds3 = +0.00459383; ddp3 = -0.00076054; ddd3 = -0.00028146; */

//1st nearest neighbour SK terms:
/* sss1 = -0.11819470; sps1 = +0.16932010; pps1 = +0.26893690; ppp1 = -0.01873386; sds1 = -0.07640028; */
/* pds1 = -0.11626894; pdp1 = +0.03329925; dds1 = -0.04294269; ddp1 = +0.01971803; ddd1 = -0.00019785; */

/* //2nd nearest neighbour SK terms: */
/* sss2 = -0.02274765; sps2 = +0.06205349; pps2 = +0.16357396; ppp2 = +0.03095619; sds2 = -0.03026521; */
/* pds2 = -0.05232678; pdp2 = -0.01359036; dds2 = -0.02364862; ddp2 = -0.00415891; ddd2 = +0.00170718; */

/* //3rd nearest neighbour SK terms: */
/* sss3 = +0.01739150; sps3 = -0.02996342; pps3 = -0.05305194; ppp3 = +0.01658685; sds3 = +0.00878227; */ 
/* pds3 = +0.01157414; pdp3 = -0.00128322; dds3 = +0.00456016; ddp3 = -0.00074586; ddd3 = -0.00027614; */

/* sss1 = -0.11784049; sps1 = +0.16914972; pps1 = +0.26684695; ppp1 = -0.01752626; sds1 = -0.07912834; */
/* pds1 = -0.11550098; pdp1 = +0.03719323; dds1 = -0.04080053; ddp1 = +0.01934707; ddd1 = -0.00063042; */

/* //2nd nearest neighbour SK terms: */
/* sss2 = -0.02261486; sps2 = +0.06190427; pps2 = +0.16294422; ppp2 = +0.03094901; sds2 = -0.02707729; */
/* pds2 = -0.04254137; pdp2 = -0.01730310; dds2 = -0.02269372; ddp2 = -0.00559271; ddd2 = +0.00263624; */

/* //3rd nearest neighbour SK terms: */
/* sss3 = +0.01767442; sps3 = -0.02705523; pps3 = -0.04979347; ppp3 = +0.01449936; sds3 = +0.01091518; */ 
/* pds3 = +0.01153676; pdp3 = -0.00146169; dds3 = +0.00501522; ddp3 = -0.00118016; ddd3 = -0.00033655; */

/* //1st nearest neighbour SK terms: */
/* sss1 = -0.11800888; sps1 = +0.16923521; pps1 = +0.26721163; ppp1 = -0.01777765; sds1 = -0.07963869; */
/* pds1 = -0.11419721; pdp1 = +0.03698884; dds1 = -0.04056494; ddp1 = +0.01927374; ddd1 = +0.00020437; */

/* //2nd nearest neighbour SK terms: */
/* sss2 = -0.02248562; sps2 = +0.06200947; pps2 = +0.16332625; ppp2 = +0.03114995; sds2 = -0.02716382; */
/* pds2 = -0.04630675; pdp2 = -0.01833265; dds2 = -0.02287608; ddp2 = -0.00650320; ddd2 = +0.00212757; */

/* //3rd nearest neighbour SK terms: */
/* sss3 = +0.01769323; sps3 = -0.02726719; pps3 = -0.04956013; ppp3 = +0.01432599; sds3 = +0.01110715; */ 
/* pds3 = +0.01173768; pdp3 = -0.00093273; dds3 = +0.00515386; ddp3 = -0.00116761; ddd3 = -0.00030219; */

//1st nearest neighbour SK terms:
/* sss1 = -0.11803129; sps1 = +0.16951493; pps1 = +0.26808825; ppp1 = -0.01796288; sds1 = -0.07687988; */
/* pds1 = -0.11862085; pdp1 = +0.03489328; dds1 = -0.04105560; ddp1 = +0.01914237; ddd1 = -0.00116003; */

/* //2nd nearest neighbour SK terms: */
/* sss2 = -0.02260240; sps2 = +0.06071330; pps2 = +0.16348580; ppp2 = +0.03062241; sds2 = -0.02915528; */
/* pds2 = -0.03976375; pdp2 = -0.01208336; dds2 = -0.02292477; ddp2 = -0.00524109; ddd2 = +0.00278232; */

/* //3rd nearest neighbour SK terms: */
/* sss3 = +0.01753364; sps3 = -0.02784025; pps3 = -0.04908586; ppp3 = +0.01465080; sds3 = +0.01060320; */ 
/* pds3 = +0.01140022; pdp3 = -0.00235438; dds3 = +0.00430353; ddp3 = -0.00073711; ddd3 = +0.00005169; */

//1st nearest neighbour SK terms:
sss1 = -0.11796930; sps1 = +0.16957449; pps1 = +0.26805300; ppp1 = -0.01785382; sds1 = -0.07686737;
pds1 = -0.11869704; pdp1 = +0.03510885; dds1 = -0.04102823; ddp1 = +0.01893526; ddd1 = -0.00095504;

//2nd nearest neighbour SK terms:
sss2 = -0.02256967; sps2 = +0.06095638; pps2 = +0.16354321; ppp2 = +0.03062320; sds2 = -0.02906798;
pds2 = -0.04110930; pdp2 = -0.01188413; dds2 = -0.02292941; ddp2 = -0.00540828; ddd2 = +0.00286941;

//3rd nearest neighbour SK terms:
sss3 = +0.01752226; sps3 = -0.02775716; pps3 = -0.04853709; ppp3 = +0.01430107; sds3 = +0.01095750; 
pds3 = +0.01172195; pdp3 = -0.00190808; dds3 = +0.00445770; ddp3 = -0.00078898; ddd3 = +0.00001258;
/* //this block from Papa */
/*       s =  1.14481; // on-site */
/*       p =  1.80769; */
/*       d1 =  0.5*(0.78456 + 0.75661); */
/*       d2 =  0.5*(0.78456 + 0.75661); */
/*       sss1 = -0.13243;   //  same atom hopping */
/*       sps1 =  0.17278; */
/*       pps1 =  0.25911; */
/*       ppp1 =  0.02653; */
/*       sds1 = -0.07145; */
/*       pds1 = -0.09702; */
/*       pdp1 =  0.02129; */
/*       dds1 = -0.05266; */
/*       ddp1 =  0.03276; */
/*       ddd1 = -0.00286; */
/*       sss2 = -0.03003; */
/*       sps2 =  0.07159; */
/*       pps2 =  0.18256; */
/*       ppp2 =  0.03703; */
/*       sds2 = -0.04075; */
/*       pds2 = -0.06522; */
/*       pdp2 = -0.00467; */
/*       dds2 = -0.03396; */
/*       ddp2 =  0.00581; */
/*       ddd2 =  0.00114; */
/*       sss3 =  0.01589; */
/*       sps3 = -0.02306; */
/*       pps3 = -0.04253; */
/*       ppp3 =  0.01538; */
/*       sds3 =  0.00016; */
/*       pds3 =  0.00222; */
/*       pdp3 = -0.00351; */
/*       dds3 =  0.00233; */
/*       ddp3 =  0.00013; */
/*       ddd3 = -0.00060; */

	VectorXd nn(10),nnn(20);
	nn<<sss1, sps1, pps1, ppp1, sds1, pds1, pdp1, dds1, ddp1, ddd1;
	nnn<<sss2, sps2, pps2, ppp2, sds2, pds2, pdp2, dds2, ddp2, ddd2, sss3, sps3, pps3, ppp3, sds3, pds3, pdp3, dds3, ddp3, ddd3;
	cout<<"full 1st nn SK potentials;"<<endl;
	cout<<"sss = "<<nn(0)<<endl;
	cout<<"sps = "<<nn(1)<<endl;
	cout<<"pps = "<<nn(2)<<endl;
	cout<<"ppp = "<<nn(3)<<endl;
	cout<<"sds = "<<nn(4)<<endl;
	cout<<"pds = "<<nn(5)<<endl;
	cout<<"pdp = "<<nn(6)<<endl;
	cout<<"dds = "<<nn(7)<<endl;
	cout<<"ddp = "<<nn(8)<<endl;
	cout<<"ddd = "<<nn(9)<<endl;
	cout<<endl;
	cout<<"full 2nd nn SK potentials;"<<endl;
	cout<<"sss = "<<nnn(0)<<endl;
	cout<<"sps = "<<nnn(1)<<endl;
	cout<<"pps = "<<nnn(2)<<endl;
	cout<<"ppp = "<<nnn(3)<<endl;
	cout<<"sds = "<<nnn(4)<<endl;
	cout<<"pds = "<<nnn(5)<<endl;
	cout<<"pdp = "<<nnn(6)<<endl;
	cout<<"dds = "<<nnn(7)<<endl;
	cout<<"ddp = "<<nnn(8)<<endl;
	cout<<"ddd = "<<nnn(9)<<endl;
	cout<<endl;
	cout<<"full 3rd nn SK potentials;"<<endl;
	cout<<"sss = "<<nnn(10)<<endl;
	cout<<"sps = "<<nnn(11)<<endl;
	cout<<"pps = "<<nnn(12)<<endl;
	cout<<"ppp = "<<nnn(13)<<endl;
	cout<<"sds = "<<nnn(14)<<endl;
	cout<<"pds = "<<nnn(15)<<endl;
	cout<<"pdp = "<<nnn(16)<<endl;
	cout<<"dds = "<<nnn(17)<<endl;
	cout<<"ddp = "<<nnn(18)<<endl;
	cout<<"ddd = "<<nnn(19)<<endl;
	/* cout<<endl<<t_13.real().bottomRightCorner(5,5)<<endl<<endl; */

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

	/* cout<<endl<<t_13.real().bottomRightCorner(5,5)<<endl<<endl; */

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

	int nt = 351;
	VectorXd vec1(nt), vec2(nt), vec3(nt), vec4(nt), vec5(nt), vec6(nt), vec7(nt), vec8(nt), vec9(nt);
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
		       	+ t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
			+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
				+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K))
				+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
				+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
				+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K))
				+ t_19*exp(i*d_19.dot(K)) + t_20*exp(i*d_20.dot(K))
				+ t_21*exp(i*d_21.dot(K)) + t_22*exp(i*d_22.dot(K))
				+ t_23*exp(i*d_23.dot(K)) + t_24*exp(i*d_24.dot(K))
				+ t_25*exp(i*d_25.dot(K)) + t_26*exp(i*d_26.dot(K))
				+ t_27*exp(i*d_27.dot(K)) + t_28*exp(i*d_28.dot(K))
				+ t_29*exp(i*d_29.dot(K)) + t_30*exp(i*d_30.dot(K))
				;

		/* E_mid = E.bottomRightCorner(5,5); */
		/* E_small = E.topLeftCorner(4,4); */

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

		vec1(k) = O(0);
		vec2(k) = O(1);
		vec3(k) = O(2);
		vec4(k) = O(3);
		vec5(k) = O(4);
		vec6(k) = O(5);
		vec7(k) = O(6);
		vec8(k) = O(7);
		vec9(k) = O(8);
		/* Myfile<<k<<" "<<O(0)<<endl; */
		/* Myfile<<k<<" "<<O(1)<<endl; */
		/* Myfile<<k<<" "<<O(2)<<endl; */
		/* Myfile<<k<<" "<<O(3)<<endl; */

		/* Myfile<<k<<" "<<O(4)<<endl; */
		/* Myfile<<k<<" "<<O(5)<<endl; */
		/* Myfile<<k<<" "<<O(6)<<endl; */
		/* Myfile<<k<<" "<<O(7)<<endl; */
		/* Myfile<<k<<" "<<O(8)<<endl; */

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
	for (int k = 0; k < vec1.size(); k++)
		Myfile<<k<<" "<<vec1(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec2.size(); k++)
		Myfile<<k<<" "<<vec2(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec3.size(); k++)
		Myfile<<k<<" "<<vec3(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec4.size(); k++)
		Myfile<<k<<" "<<vec4(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec5.size(); k++)
		Myfile<<k<<" "<<vec5(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec6.size(); k++)
		Myfile<<k<<" "<<vec6(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec7.size(); k++)
		Myfile<<k<<" "<<vec7(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec8.size(); k++)
		Myfile<<k<<" "<<vec8(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec9.size(); k++)
		Myfile<<k<<" "<<vec9(k)<<endl;

	Myfile.close();
	/* Myfile2.close(); */
	/* Myfile3.close(); */
	/* Myfile4.close(); */
	return 0;
}
