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
  const int n = 8;
  Vector3d K;
  double k_x, k_y, k_z;
  const double A = M_PI;
  int it = T.size();
  for (int k = 0; k!=n+1; k++){
    /* if (k%2!=0){ */
      k_x = A*k/n;
      for (int l = 0; l!=k+1; l++){
        /* if (l%2!=0){ */
          k_y = A*l/n;
          for (int m = 0; m!=(k-l)/2 + 1; m++){
            /* if (m%2!=0){ */
              k_z = A*m/n;
	      K <<k_x, k_y, k_z;
	      E = U;
	      for (int iter = 0; iter < it; iter++)
	        E = E + T[iter]*exp(i*D[iter].dot(K)/2.);
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
 	    /* } */
          /* } */
	/* } */
      }
    }
  }
  return;
}

template <typename func1>
void dxargs(VectorXd& p, VectorXd& xi, func1&& func, const vector<VectorXd>& Weigs,
	       	const m3 &U, const vector<Vector3d> &D, 
		const VectorXd &nn, const int numnn, const int numnnn){
	m3 t;
	vector<m3> T;
	int sz, sz1, sz2, sz3, sz4, sz6, total;
	sz = nn.size();
	sz2 = p.size();
	sz3 = sz + sz2;
	sz1 = sz3 - 10;
	sz6 = D.size();
	sz4 = sz2 - 10;
	total = numnn + numnnn;
	VectorXd p1(10), p2(sz1);
	vector<Vector3d> Dderiv, Deriv2;
	int start1, end1, end2, end3;
	if (sz == 0){
		start1 = 0;
		end1 = numnn;
		end2 = total;
		end3 = sz6;
		Dderiv = D;
	       	p1 = p.head(10);
		if (sz2 > 10)
		        p2 = p.tail(sz4);
	}
	if (sz == 10){
		vector<Vector3d>::const_iterator first = D.begin() + numnn;
		vector<Vector3d>::const_iterator last = D.end();
		vector<Vector3d> Dderiv2(first,last);
		Dderiv = Dderiv2;
	       	p1 = nn;
		p2 = p;
		start1 = 0;
		end1 = numnnn;
		end2 = Dderiv.size();
		end3 = 0;
	}
	if (sz == 20){
		vector<Vector3d>::const_iterator first = D.begin() + total;
		vector<Vector3d>::const_iterator last = D.end();
		vector<Vector3d> Dderiv2(first,last);
		Dderiv = Dderiv2;
	       	p1 = nn.head(10);
		p2.head(10) = nn.tail(10);
		p2.tail(10) = p;
		start1 = 0;
		end1 = Dderiv.size();
		end2 = 0;
		end3 = 0;
	}
	for (int it = 0; it < numnn; it++){
		t = TB(0,1,0,9,D[it],p1,p1);
		T.emplace_back(t);
	}
	if (sz3 > 10){
		for (int it = numnn; it < total; it++){
			t = TB(0,1,1,9,D[it],p1,p2);
			T.emplace_back(t);
		}
		if (sz3 > 20){
			for (int it = total; it < sz6; it++){
				t = TB(0,1,2,9,D[it],p1,p2);
				T.emplace_back(t);
			}
		}
	}
	vector<VectorXd> Teigs, dummy;
	m3 dummy2;
	dummy2.fill(0.);
	vector<m3> Eig_vecs, dE;
	eigs(T, D, U, false, true, true, Eig_vecs, Teigs);
	double fret;
	double weight;
	double V;
	int start;
	int end;
	for (int j = 0; j < p.size(); j++){
		T.clear();
		if (j < 10){
			start = start1;
			end = end1;
		}
		if (j > 9){
			start = end1;
			end = end2;
		}
		if (j > 19){
			start = end2;
			end = end3;
		}
		if (j%10 == 0){
			for (int it = start; it < end; it++){
				t = dsss();
				T.emplace_back(t);
			}
		}
		if (j%10 == 1){
			for (int it = start; it < end; it++){
				t = dsps(Dderiv[it]);
				T.emplace_back(t);
			}
		}
		if (j%10 == 2){
			for (int it = start; it < end; it++){
				t = dpps(Dderiv[it]);
				T.emplace_back(t);
			}
		}
		if (j%10 == 3){
			for (int it = start; it < end; it++){
				t = dppp(Dderiv[it]);
				T.emplace_back(t);
			}
		}
		if (j%10 == 4){
			for (int it = start; it < end; it++){
				t = dsds(Dderiv[it]);
				T.emplace_back(t);
			}
		}
		if (j%10 == 5){
			for (int it = start; it < end; it++){
				t = dpds(Dderiv[it]);
				T.emplace_back(t);
			}
		}
		if (j%10 == 6){
			for (int it = start; it < end; it++){
				t = dpdp(Dderiv[it]);
				T.emplace_back(t);
			}
		}
		if (j%10 == 7){
			for (int it = start; it < end; it++){
				t = ddds(Dderiv[it]);
				T.emplace_back(t);
			}
		}
		if (j%10 == 8){
			for (int it = start; it < end; it++){
				t = dddp(Dderiv[it]);
				T.emplace_back(t);
			}
		}
		if (j%10 == 9){
			for (int it = start; it < end; it++){
				t = dddd(Dderiv[it]);
				T.emplace_back(t);
			}
		}
		eigs(T, Dderiv, dummy2, true, false, false, dE, dummy);
		fret = 0.;
		for (int it = 0; it < Teigs.size(); it ++){
			for (int k = 0; k < Teigs[0].size(); k++){
				if (k > 5)
					weight = 0.2;
				else
					weight = 2.;
				V = (Eig_vecs[it].col(k).adjoint()*dE[it]*Eig_vecs[it].col(k)).real()(0);
				fret += 2.*V*(Teigs[it](k) - Weigs[it](k))*weight;
			}
		}
		xi(j) = fret;
	}
}

double xargs(VectorXd& p, const vector<VectorXd>& Weigs, const m3 &U, const vector<Vector3d> &D,
		const VectorXd &nn, const int numnn, const int numnnn){

	m3 t;
	vector<m3> T;
	int sz, sz1, sz2, sz3, sz4, sz5, total;
	sz = nn.size();
	sz2 = p.size();
	sz3 = sz + sz2;
	sz1 = sz3 - 10;
	sz4 = sz2 - 10;
	total = numnn + numnnn;
	sz5 = D.size();
	VectorXd p1(10), p2(sz1);
	if (sz == 0){
	       	p1 = p.head(10);
		if (sz2 > 10)
		        p2 = p.tail(sz4);
	}
	if (sz == 10){
	       	p1 = nn;
		p2 = p;
	}
	if (sz == 20){
	       	p1 = nn.head(10);
		p2.head(10) = nn.tail(10);
		p2.tail(10) = p;
	}
	for (int it = 0; it < numnn; it++){
		t = TB(0,1,0,9,D[it],p1,p1);
		T.emplace_back(t);
	}
	if (sz3 > 10){
		for (int it = numnn; it < total; it++){
			t = TB(0,1,1,9,D[it],p1,p2);
			T.emplace_back(t);
		}
		if (sz3 > 20){
			for (int it = total; it < sz5; it++){
				t = TB(0,1,2,9,D[it],p1,p2);
				T.emplace_back(t);
			}
		}
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
				weight = 2.;
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
	sss2 = -0.0227063;
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
	VectorXd dummy;
	const int numnn = 8;
	const int numnnn = 6;
	cout<<"Number of k-points in k-point mesh: "<<Weigs.size()<<endl<<endl;
	cout<<"Find the SK potentials of the 1st neighbour terms:"<<endl;
	frprmn(nn, ftol, iter, fret, xargs, Weigs, U, D, dummy, numnn, numnnn);
	cout<<showpos<<fixed<<setprecision(8)<<"sss1 = "<<nn(0)<<"; sps1 = "<<nn(1)<<"; pps1 = "<<nn(2)<<"; ppp1 = "
		<<nn(3)<<"; sds1 = "<<nn(4)<<";"<<endl<<"pds1 = "<<nn(5)<<"; pdp1 = "<<nn(6)<<
		"; dds1 = "<<nn(7)<<"; ddp1 = "<<nn(8)<<"; ddd1 = "<<nn(9)<<";"<<endl;
	cout<<endl;
	cout<<noshowpos<<"number of iterations: "<<iter<<endl;
	cout<<"minimum of function: "<<fret<<endl;
	T.emplace_back(t_13);
	T.emplace_back(t_14);
	T.emplace_back(t_15);
	T.emplace_back(t_16);
	T.emplace_back(t_17);
	T.emplace_back(t_18);
	D.emplace_back(d_13);
	D.emplace_back(d_14);
	D.emplace_back(d_15);
	D.emplace_back(d_16);
	D.emplace_back(d_17);
	D.emplace_back(d_18);
	eigs(T, D, lambda, false, true, false, dummyE, Weigs);
	cout<<string(100, '*')<<endl;
	cout<<"Now find the SK potentials of the 2nd neighbour terms, keeping 1st neighbour terms fixed:"<<endl;
	frprmn(nnn, ftol, iter, fret, xargs, Weigs, U, D, nn, numnn, numnnn);
	cout<<showpos<<"sss2 = "<<nnn(0)<<"; sps2 = "<<nnn(1)<<"; pps2 = "<<nnn(2)<<"; ppp2 = "
		<<nnn(3)<<"; sds2 = "<<nnn(4)<<";"<<endl<<"pds2 = "<<nnn(5)<<"; pdp2 = "<<nnn(6)<<
		"; dds2 = "<<nnn(7)<<"; ddp2 = "<<nnn(8)<<"; ddd2 = "<<nnn(9)<<";"<<endl;
	cout<<endl;
	cout<<noshowpos<<"number of iterations: "<<iter<<endl;
	cout<<"minimum of function: "<<fret<<endl;
	cout<<string(100, '*')<<endl;
	VectorXd NNN(20);
	NNN.head(10) = nn;
	NNN.tail(10) = nnn;
	cout<<"Now find the SK potentials up to 2nd neighbour terms:"<<endl;
	frprmn(NNN, ftol, iter, fret, xargs, Weigs, U, D, dummy, numnn, numnnn);
	cout<<"1st nearest neighbour SK terms:"<<endl;
	cout<<showpos<<"sss1 = "<<NNN(0)<<"; sps1 = "<<NNN(1)<<"; pps1 = "<<NNN(2)<<"; ppp1 = "
		<<NNN(3)<<"; sds1 = "<<NNN(4)<<";"<<endl<<"pds1 = "<<NNN(5)<<"; pdp1 = "<<NNN(6)<<
		"; dds1 = "<<NNN(7)<<"; ddp1 = "<<NNN(8)<<"; ddd1 = "<<NNN(9)<<";"<<endl;
	cout<<endl;
	cout<<"2nd nearest neighbour SK terms:"<<endl;
	cout<<"sss2 = "<<NNN(10)<<"; sps2 = "<<NNN(11)<<"; pps2 = "<<NNN(12)<<"; ppp2 = "
		<<NNN(13)<<"; sds2 = "<<NNN(14)<<";"<<endl<<"pds2 = "<<NNN(15)<<"; pdp2 = "<<NNN(16)<<
		"; dds2 = "<<NNN(17)<<"; ddp2 = "<<NNN(18)<<"; ddd2 = "<<NNN(19)<<";"<<endl;
	cout<<endl;
	cout<<noshowpos<<"number of iterations: "<<iter<<endl;
	cout<<"minimum of function: "<<fret<<endl;
	cout<<string(100, '*')<<endl;
	T.emplace_back(t_19);
	T.emplace_back(t_20);
	T.emplace_back(t_21);
	T.emplace_back(t_22);
	T.emplace_back(t_23);
	T.emplace_back(t_24);
	T.emplace_back(t_25);
	T.emplace_back(t_26);
	T.emplace_back(t_27);
	T.emplace_back(t_28);
	T.emplace_back(t_29);
	T.emplace_back(t_30);
	D.emplace_back(d_19);
	D.emplace_back(d_20);
	D.emplace_back(d_21);
	D.emplace_back(d_22);
	D.emplace_back(d_23);
	D.emplace_back(d_24);
	D.emplace_back(d_25);
	D.emplace_back(d_26);
	D.emplace_back(d_27);
	D.emplace_back(d_28);
	D.emplace_back(d_29);
	D.emplace_back(d_30);
	eigs(T, D, lambda, false, true, false, dummyE, Weigs);
	cout<<"Now find 3rd neighbour terms, keeping 1st & 2nd neighbour terms fixed:"<<endl;
	VectorXd nnnn(10);
	nnnn.fill(0.);
	frprmn(nnnn, ftol, iter, fret, xargs, Weigs, U, D, NNN, numnn, numnnn);
	cout<<showpos<<"sss3 = "<<nnnn(0)<<"; sps3 = "<<nnnn(1)<<"; pps3 = "<<nnnn(2)<<"; ppp3 = "
		<<nnnn(3)<<"; sds3 = "<<nnnn(4)<<";"<<endl<<"pds3 = "<<nnnn(5)<<"; pdp3 = "<<nnnn(6)<<
		"; dds3 = "<<nnnn(7)<<"; ddp3 = "<<nnnn(8)<<"; ddd3 = "<<nnnn(9)<<";"<<endl;
	cout<<endl;
	cout<<noshowpos<<"number of iterations: "<<iter<<endl;
	cout<<"minimum of function: "<<fret<<endl;
	cout<<string(100, '*')<<endl;
	VectorXd NNNN(30);
	NNNN.head(20) = NNN;
	NNNN.tail(10) = nnnn;
	cout<<"Finally, find up to 3rd neighbour terms:"<<endl;
	frprmn(NNNN, ftol, iter, fret, xargs, Weigs, U, D, dummy, numnn, numnnn);
	cout<<"1st nearest neighbour SK terms:"<<endl;
	cout<<showpos<<"sss1 = "<<NNNN(0)<<"; sps1 = "<<NNNN(1)<<"; pps1 = "<<NNNN(2)<<"; ppp1 = "
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
	cout<<noshowpos<<"number of iterations: "<<iter<<endl;
	cout<<"minimum of function: "<<fret<<endl;
	cout<<string(100, '*')<<endl;

	return 0;
}
