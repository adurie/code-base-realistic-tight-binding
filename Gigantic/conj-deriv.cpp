#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "cgls-derivs.h"
#include <vector>
#include "TB_sk.h"
#include <iomanip>
#include <fstream>
#include "derivs.h"

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

void eigs(const vector<m3> &T, const vector<Vector3d> &D, const m3 &U, bool is_E,
	       	bool is_eigval, bool is_eigvec,	vector<m3> &dummyE_or_vec, vector<VectorXd> &dummy_val){
  dummyE_or_vec.clear();
  dummy_val.clear();
  if ((is_E == true) && (is_eigval == true))
    cout<<"Error, is_E & is_eigval both set to true"<<endl;
  m3 E;
  dcomp i;
  i = -1.;
  i = sqrt(i);
  const int n = 20;
  Vector3d K;
  double k_x, k_y, k_z;
  const double A = M_PI;
  int it = D.size();
  for (int k = 0; k!=n+1; k++){
    if (k%2!=0){
      k_x = A*k/n;
      for (int l = 0; l!=k+1; l++){
        if (l%2!=0){
          k_y = A*l/n;
          for (int m = 0; m!=(k-l)/2 + 1; m++){
            if (m%2!=0){
              k_z = A*m/n;
	      K <<k_x, k_y, k_z;
	      E = U;
	      for (int iter = 0; iter < it; iter++)
	        E = E + T[iter]*exp(i*D[iter].dot(K));
	      if (is_E == true)
	        dummyE_or_vec.emplace_back(E);
	      if (is_eigval == true){
		ComplexEigenSolver<m3> CA;
		VectorXd eigs;
		if (is_eigvec == true){
		  CA.compute(E,true);
		  eigs = CA.eigenvalues().real();
		  dummy_val.emplace_back(eigs);
		  m3 O;
                  O = CA.eigenvectors();
		  dummyE_or_vec.emplace_back(O);
		}
		else{
		  CA.compute(E,false);
		  eigs = CA.eigenvalues().real();
		  /* sort(eigs.data(), eigs.data()+eigs.size()); */
		  dummy_val.emplace_back(eigs);
		}
	      }

	      
 	    }
          }
	}
      }
    }
  }
  return;
}

template <typename func1>
void dxargs(VectorXd& p, VectorXd& xi, func1&& func, const vector<VectorXd>& Weigs,
	       	const m3 &U, const vector<Vector3d> &D){
	m3 t;
	vector<m3> T;
	for (int it = 0; it < 8; it++){//!!! 8 will need changing as NN increases!
		t = TB(0,1,0,9,D[it],p,p);
		T.emplace_back(t);
	}
	vector<VectorXd> Teigs, dummy;
	vector<m3> Eig_vecs, dE;
	eigs(T, D, U, false, true, true, Eig_vecs, Teigs);
	double fret;
	double weight;
	double V;
	for (int j = 0; j < p.size(); j++){
		if (j == 0){
			T.clear();
			fret = 0.;
			for (int it = 0; it < 8; it++){//!!! 8 will need changing as NN increases!
				t = dsss();
				T.emplace_back(t);
			}
			eigs(T, D, U, true, false, false, dE, dummy);
			for (int it = 0; it < Teigs.size(); it ++){
				for (int k = 0; k < Teigs[0].size(); k++){
					if (k > 5)
						weight = 0.2;
					else
						weight = 1.;
					V = (Eig_vecs[it].col(k).adjoint()*dE[it]*Eig_vecs[it].col(k)).real()(0);
					fret += 2.*V*(Teigs[it](k) - Weigs[it](k))*weight;
				}
			}
		}
		if (j == 1){
			T.clear();
			fret = 0.;
			for (int it = 0; it < 8; it++){//!!! 8 will need changing as NN increases!
				t = dsps(D[it]);
				T.emplace_back(t);
			}
			eigs(T, D, U, true, false, false, dE, dummy);
			for (int it = 0; it < Teigs.size(); it ++){
				for (int k = 0; k < Teigs[0].size(); k++){
					if (k > 5)
						weight = 0.2;
					else
						weight = 1.;
					V = (Eig_vecs[it].col(k).adjoint()*dE[it]*Eig_vecs[it].col(k)).real()(0);
					fret += 2.*V*(Teigs[it](k) - Weigs[it](k))*weight;
				}
			}
		}
		if (j == 2){
			T.clear();
			fret = 0.;
			for (int it = 0; it < 8; it++){//!!! 8 will need changing as NN increases!
				t = dpps(D[it]);
				T.emplace_back(t);
			}
			eigs(T, D, U, true, false, false, dE, dummy);
			for (int it = 0; it < Teigs.size(); it ++){
				for (int k = 0; k < Teigs[0].size(); k++){
					if (k > 5)
						weight = 0.2;
					else
						weight = 1.;
					V = (Eig_vecs[it].col(k).adjoint()*dE[it]*Eig_vecs[it].col(k)).real()(0);
					fret += 2.*V*(Teigs[it](k) - Weigs[it](k))*weight;
				}
			}
		}
		if (j == 3){
			T.clear();
			fret = 0.;
			for (int it = 0; it < 8; it++){//!!! 8 will need changing as NN increases!
				t = dppp(D[it]);
				T.emplace_back(t);
			}
			eigs(T, D, U, true, false, false, dE, dummy);
			for (int it = 0; it < Teigs.size(); it ++){
				for (int k = 0; k < Teigs[0].size(); k++){
					if (k > 5)
						weight = 0.2;
					else
						weight = 1.;
					V = (Eig_vecs[it].col(k).adjoint()*dE[it]*Eig_vecs[it].col(k)).real()(0);
					fret += 2.*V*(Teigs[it](k) - Weigs[it](k))*weight;
				}
			}
		}
		if (j == 4){
			T.clear();
			fret = 0.;
			for (int it = 0; it < 8; it++){//!!! 8 will need changing as NN increases!
				t = dsds(D[it]);
				T.emplace_back(t);
			}
			eigs(T, D, U, true, false, false, dE, dummy);
			for (int it = 0; it < Teigs.size(); it ++){
				for (int k = 0; k < Teigs[0].size(); k++){
					if (k > 5)
						weight = 0.2;
					else
						weight = 1.;
					V = (Eig_vecs[it].col(k).adjoint()*dE[it]*Eig_vecs[it].col(k)).real()(0);
					fret += 2.*V*(Teigs[it](k) - Weigs[it](k))*weight;
				}
			}
		}
		if (j == 5){
			T.clear();
			fret = 0.;
			for (int it = 0; it < 8; it++){//!!! 8 will need changing as NN increases!
				t = dpds(D[it]);
				T.emplace_back(t);
			}
			eigs(T, D, U, true, false, false, dE, dummy);
			for (int it = 0; it < Teigs.size(); it ++){
				for (int k = 0; k < Teigs[0].size(); k++){
					if (k > 5)
						weight = 0.2;
					else
						weight = 1.;
					V = (Eig_vecs[it].col(k).adjoint()*dE[it]*Eig_vecs[it].col(k)).real()(0);
					fret += 2.*V*(Teigs[it](k) - Weigs[it](k))*weight;
				}
			}
		}
		if (j == 6){
			T.clear();
			fret = 0.;
			for (int it = 0; it < 8; it++){//!!! 8 will need changing as NN increases!
				t = dpdp(D[it]);
				T.emplace_back(t);
			}
			eigs(T, D, U, true, false, false, dE, dummy);
			for (int it = 0; it < Teigs.size(); it ++){
				for (int k = 0; k < Teigs[0].size(); k++){
					if (k > 5)
						weight = 0.2;
					else
						weight = 1.;
					V = (Eig_vecs[it].col(k).adjoint()*dE[it]*Eig_vecs[it].col(k)).real()(0);
					fret += 2.*V*(Teigs[it](k) - Weigs[it](k))*weight;
				}
			}
		}
		if (j == 7){
			T.clear();
			fret = 0.;
			for (int it = 0; it < 8; it++){//!!! 8 will need changing as NN increases!
				t = ddds(D[it]);
				T.emplace_back(t);
			}
			eigs(T, D, U, true, false, false, dE, dummy);
			for (int it = 0; it < Teigs.size(); it ++){
				for (int k = 0; k < Teigs[0].size(); k++){
					if (k > 5)
						weight = 0.2;
					else
						weight = 1.;
					V = (Eig_vecs[it].col(k).adjoint()*dE[it]*Eig_vecs[it].col(k)).real()(0);
					fret += 2.*V*(Teigs[it](k) - Weigs[it](k))*weight;
				}
			}
		}
		if (j == 8){
			T.clear();
			fret = 0.;
			for (int it = 0; it < 8; it++){//!!! 8 will need changing as NN increases!
				t = dddp(D[it]);
				T.emplace_back(t);
			}
			eigs(T, D, U, true, false, false, dE, dummy);
			for (int it = 0; it < Teigs.size(); it ++){
				for (int k = 0; k < Teigs[0].size(); k++){
					if (k > 5)
						weight = 0.2;
					else
						weight = 1.;
					V = (Eig_vecs[it].col(k).adjoint()*dE[it]*Eig_vecs[it].col(k)).real()(0);
					fret += 2.*V*(Teigs[it](k) - Weigs[it](k))*weight;
				}
			}
		}
		if (j == 9){
			T.clear();
			fret = 0.;
			for (int it = 0; it < 8; it++){//!!! 8 will need changing as NN increases!
				t = dddd(D[it]);
				T.emplace_back(t);
			}
			eigs(T, D, U, true, false, false, dE, dummy);
			for (int it = 0; it < Teigs.size(); it ++){
				for (int k = 0; k < Teigs[0].size(); k++){
					if (k > 5)
						weight = 0.2;
					else
						weight = 1.;
					V = (Eig_vecs[it].col(k).adjoint()*dE[it]*Eig_vecs[it].col(k)).real()(0);
					fret += 2.*V*(Teigs[it](k) - Weigs[it](k))*weight;
				}
			}
		}
		xi(j) = fret;
	}
}

double xargs(VectorXd& p, const vector<VectorXd>& Weigs, const m3 &U, const vector<Vector3d> &D){

	m3 t;
	vector<m3> T;
	for (int it = 0; it < 8; it++){
		t = TB(0,1,0,9,D[it],p,p);
		T.emplace_back(t);
	}

	vector<VectorXd> Teigs;
	vector<m3> dummyE;
	eigs(T, D, U, false, true, false, dummyE, Teigs);
	double fret = 0.;
	double weight;
	for (int it = 0; it < Teigs.size(); it ++){
		for (int k = 0; k < Teigs[0].size(); k++){
			if (k > 5)
				weight = 0.2;
			else
				weight = 1.;
			fret += (Teigs[it](k) - Weigs[it](k))*(Teigs[it](k) - Weigs[it](k))*weight;
		}
	}

	return fret;

}

int main(){
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
	vector<m3> T;
	T.emplace_back(t_1);
	T.emplace_back(t_2);
	T.emplace_back(t_3);
	T.emplace_back(t_4);
	T.emplace_back(t_5);
	T.emplace_back(t_6);
	T.emplace_back(t_7);
	T.emplace_back(t_8);
	vector<Vector3d> D;
	D.emplace_back(d_1);
	D.emplace_back(d_2);
	D.emplace_back(d_3);
	D.emplace_back(d_4);
	D.emplace_back(d_5);
	D.emplace_back(d_6);
	D.emplace_back(d_7);
	D.emplace_back(d_8);
	vector<m3> dummyE;
	vector<VectorXd> Weigs;
	eigs(T, D, lambda, false, true, false, dummyE, Weigs);

	double s, p, d1, d2, sss1, sss2, pps1, pps2, ppp1, ppp2, dds1, dds2, ddp1, ddp2, ddd1, ddd2, sps1, sps2, sds1, sds2, pds1, pds2, pdp1, pdp2;
	//intial guess
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
	m3 U;
	U(0,0) = lambda(0,0);
	U(1,1) = U(2,2) = U(3,3) = (lambda(1,1) + lambda(2,2) + lambda(3,3))/3.;
	U(4,4) = U(5,5) = U(6,6) = (lambda(4,4) + lambda(5,5) + lambda(6,6))/3.;
	U(7,7) = U(8,8) = (lambda(7,7) + lambda(8,8))/2.;

	double ftol = 1e-5;
	double fret;
	int iter;
	cout<<"Find the SK potentials of the 1st neighbour terms:"<<endl;
	frprmn(nn, ftol, iter, fret, xargs, Weigs, U, D);
	cout<<"sss1 = "<<nn(0)<<"; sps1 = "<<nn(1)<<"; pps1 = "<<nn(2)<<"; ppp1 = "
		<<nn(3)<<"; sds1 = "<<nn(4)<<";"<<endl<<"pds1 = "<<nn(5)<<"; pdp1 = "<<nn(6)<<
		"; dds1 = "<<nn(7)<<"; ddp1 = "<<nn(8)<<"; ddd1 = "<<nn(9)<<";"<<endl;
	cout<<endl;
	cout<<"number of iterations: "<<iter<<endl;
	cout<<"minimum of function: "<<fret<<endl;
	cout<<endl<<endl;

	return 0;
}
