#include <iostream>
#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Dense>
/* #include "TB.h" */
/* #include "TBdynamic.h" */

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
//calculates the bandstructure of fcc Cu
//
Matrix<dcomp, 9, 9> read(int nn, Vector3d &dvec, int ispin){
      string input;
      if (ispin == +1)
        input = "iron_up_hr.dat";
      if (ispin == -1)
        input = "iron_dn_hr.dat";
      ifstream infile(input);
      Matrix<dcomp, 9, 9> rt;
      rt.fill(0.);
  
      dcomp i;
      i = -1;
      i = sqrt(i);
      string line;
      double a, b, c, d, e, f, g;
      double eV_Ry = 0.073498618;
      double iron_fermi = 13.4; // obtained from Wannier90 scf.out
      iron_fermi *= eV_Ry;
      double silver_fermi = 0.4635; //this is lazy - amend to carry this through
      double cshift = silver_fermi - iron_fermi;
      complex<double> im;
      if (nn == 0){
        while (!infile.eof()) 
        {
		getline(infile, line);
		istringstream iss(line);
		iss >> a >> b >> c >> d >> e >> f >> g;
		if ((a == 0) && (b == 0) && (c == 0))
		{
			if (d == e)
				rt(d-1,e-1) = (f + g*i)*eV_Ry + cshift;
		}
        }
      }
      else{
        while (!infile.eof()) 
        {
		getline(infile, line);
		istringstream iss(line);
		iss >> a >> b >> c >> d >> e >> f >> g;
		if ((a == dvec(0)) && (b == dvec(1)) && (c == dvec(2)))
			rt(d-1,e-1) = (f + g*i)*eV_Ry;
        }
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
	int ispin = +1;

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
	u = read(0, d_1, ispin);
	t_1 = read(1, d_1, ispin);
	t_2 = read(1, d_2, ispin);
	t_3 = read(1, d_3, ispin);
	t_4 = read(1, d_4, ispin);
	t_5 = read(1, d_5, ispin);
	t_6 = read(1, d_6, ispin);
	t_7 = read(1, d_7, ispin);
	t_8 = read(1, d_8, ispin);

	t_13 = read(2, d_13, ispin);
	t_14 = read(2, d_14, ispin);
	t_15 = read(2, d_15, ispin);
	t_16 = read(2, d_16, ispin);
	t_17 = read(2, d_17, ispin);
	t_18 = read(2, d_18, ispin);
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
