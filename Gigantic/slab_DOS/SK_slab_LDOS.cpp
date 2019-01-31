#include <cmath>
#define EIGEN_USE_MKL_ALL
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <iomanip>
#include <fstream>
#include "TB_Co_Cu.h"

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
typedef Vector2d v2;
typedef Vector3d v3;
// TODO make Matrix sizes more general - 45x45 = (9x5)x(9x5) i.e 9 bands, 5 atoms, dependent on input file
typedef vector<Matrix<dcomp, 45, 45>, aligned_allocator<Matrix<dcomp, 45, 45>>> vm;
typedef vector<v2, aligned_allocator<v2>> vv;
typedef Matrix<dcomp, 45, 45> bigM;
typedef Matrix<dcomp, 9, 9> smallM;

//Purpose: To calculate the partial DOS. Current system, slab geometry with 5 atom unit cell
//but *hr.dat decides system.

bigM supercell(smallM &ii, smallM &iip1, smallM &iip2, smallM &ip1i, smallM &ip2i, smallM &Z){
	bigM H;
	/* H << ii, iip1, iip2, Z, Z, */
	/*      ip1i, ii, iip1, iip2, Z, */
	/*      ip2i, ip1i, ii, iip1, iip2, */
	/*      Z, ip2i, ip1i, ii, iip1, */
	/*      Z, Z, ip2i, ip1i, ii; */
	for (int i = 0; i < 5; i++)
		H.block(i*9, i*9, 9, 9) = ii;
	for (int i = 0; i < 4; i++){
		H.block(i*9, i*9 + 9, 9, 9) = iip1;
		H.block(i*9 + 9, i*9, 9, 9) = ip1i;
	}
	for (int i = 0; i < 3; i++){
		H.block(i*9, i*9 + 18, 9, 9) = iip2;
		H.block(i*9 + 18, i*9, 9, 9) = ip2i;
	}
	for (int i = 0; i < 2; i++){
		H.block(i*9, i*9 + 27, 9, 9) = Z;
		H.block(i*9 + 27, i*9, 9, 9) = Z;
	}
	H.block(0, 36, 9, 9) = Z;
	H.block(36, 0, 9, 9) = Z;
	double Ry_eV = 1./0.073498618; // convert from Rydbergs to eV 
	H = Ry_eV*H;
	return H;
}

Matrix<dcomp, 45, 45> H_k(double k_x, double k_y, vm &Ham, vv &pos){
	// diagonalise Hamiltonian in k-space representation, here H_k(k_x, k_y)
	dcomp im;
	im = -1.;
	im = sqrt(im);
	v2 K;
	K(0) = k_x;
	K(1) = k_y;
	Matrix<dcomp, 45, 45> Ham_k;
	Ham_k.fill(0.);
	for (int i = 0; i < pos.size(); i++)
		Ham_k = Ham_k + Ham[i]*exp(im*pos[i].dot(K));
	return Ham_k;
}

/* template <typename... Args> */
double Greens(double eps, vector<double> &ham){
	// Calculate G(eps, k||) = (eps - H(k||) + im*delta)^(-1)
	// TODO make this perform a trace over a subset of G, 
	// to calculate PDOS
	dcomp im;
	im = -1.;
	im = sqrt(im);
	const double delta = 1e-5;
	double rho = 0.; //DOS, rho = (-1/pi)*ImTrG
	dcomp G = 0.;
	for (int k = 0; k<ham.size(); k++)
		G += 1./(eps + im*delta - ham[k]);
	rho = imag(G);
	return rho;
}

int main(){
	const int ispin = +1; // TODO shoudn't this be determined at runtime?

	int index;
	if (ispin == +1)
		index = 0;
	else if (ispin == -1)
		index = 1;
	else{
		cout<<"Error spin polarisation not defined"<<endl;
		exit(EXIT_FAILURE);
	}

	// generate 9x9 SK hams, naming convention should be obvious;
	smallM U, tn100, tnm100, tn010, tnm010, tnhhh, tnmhmhh, tnmhhh, tnhmhh,
	       tnmhmhmh, tnhhmh, tnhmhmh, tnmhhmh, tnn110, tnnm1m10, tnn1m10,
	       tnnm110, tnn002, tnnm002;
	v3 v3tmp;

	v3tmp << 1, 0, 0;
	//onsite matrix (in 9x9)
	U = TB(index, 0, 0, 8, v3tmp);
	//first neighbours (in 9x9 - I say this as definition changes in supercell)
	tn100    = TB(index, 1, 0, 8, v3tmp);
	v3tmp << -1, 0, 0;
	tnm100   = TB(index, 1, 0, 8, v3tmp);
	v3tmp << 0, 1, 0;
	tn010    = TB(index, 1, 0, 8, v3tmp);
	v3tmp << 0, -1, 0;
	tnm010   = TB(index, 1, 0, 8, v3tmp);
	v3tmp << 0.5, 0.5, M_SQRT1_2;
	tnhhh    = TB(index, 1, 0, 8, v3tmp);
	v3tmp << -0.5, -0.5, M_SQRT1_2;
	tnmhmhh  = TB(index, 1, 0, 8, v3tmp);
	v3tmp << -0.5, 0.5, M_SQRT1_2;
	tnmhhh   = TB(index, 1, 0, 8, v3tmp);
	v3tmp << 0.5, -0.5, M_SQRT1_2;
	tnhmhh   = TB(index, 1, 0, 8, v3tmp);
	v3tmp << -0.5, -0.5, -M_SQRT1_2;
	tnmhmhmh = TB(index, 1, 0, 8, v3tmp);
	v3tmp << 0.5, 0.5, -M_SQRT1_2;
	tnhhmh   = TB(index, 1, 0, 8, v3tmp);
	v3tmp << 0.5, -0.5, -M_SQRT1_2;
	tnhmhmh  = TB(index, 1, 0, 8, v3tmp);
	v3tmp << -0.5, 0.5, -M_SQRT1_2;
	tnmhhmh  = TB(index, 1, 0, 8, v3tmp);
	//2nd neighbours in 9x9
	v3tmp << 1, 1, 0;
	tnn110   = TB(index, 1, 1, 8, v3tmp);
	v3tmp << -1, -1, 0;
	tnnm1m10 = TB(index, 1, 1, 8, v3tmp);
	v3tmp << 1, -1, 0;
	tnn1m10  = TB(index, 1, 1, 8, v3tmp);
	v3tmp << -1, 1, 0;
	tnnm110  = TB(index, 1, 1, 8, v3tmp);
	v3tmp << 0, 0, M_SQRT2;
	tnn002   = TB(index, 1, 1, 8, v3tmp);
	v3tmp << 0, 0, -M_SQRT2;
	tnnm002  = TB(index, 1, 1, 8, v3tmp);

	//now to generate 45x45 supercell hams onsite and offsite with position vecs
	bigM bigtmp;
	v2 in_plane_latt;
	vm Ham; // Dimensions based on number of bands & atoms in unit cell. TODO generalise
	vv pos; // pos to contain all in plane position vectors.
	smallM Z;
	Z = smallM::Zero();

	in_plane_latt << 0, 0;
	bigtmp = supercell(U, tnhhh, tnn002, tnmhmhmh, tnnm002, Z);
	pos.emplace_back(in_plane_latt);
	Ham.emplace_back(bigtmp);

	in_plane_latt << 1, 0;
	bigtmp = supercell(tn100, Z, Z, tnhmhmh, Z, Z);
	pos.emplace_back(in_plane_latt);
	Ham.emplace_back(bigtmp);

	in_plane_latt << -1, 0;
	bigtmp = supercell(tnm100, tnmhhh, Z, Z, Z, Z);
	pos.emplace_back(in_plane_latt);
	Ham.emplace_back(bigtmp);

	in_plane_latt << 0, 1;
	bigtmp = supercell(tn010, Z, Z, tnmhhmh, Z, Z);
	pos.emplace_back(in_plane_latt);
	Ham.emplace_back(bigtmp);

	in_plane_latt << 0, -1;
	bigtmp = supercell(tnm010, tnhmhh, Z, Z, Z, Z);
	pos.emplace_back(in_plane_latt);
	Ham.emplace_back(bigtmp);

	in_plane_latt << 1, 1;
	bigtmp = supercell(tnn110, Z, Z, tnhhmh, Z, Z);
	pos.emplace_back(in_plane_latt);
	Ham.emplace_back(bigtmp);

	in_plane_latt << -1, -1;
	bigtmp = supercell(tnnm1m10, tnmhmhh, Z, Z, Z, Z);
	pos.emplace_back(in_plane_latt);
	Ham.emplace_back(bigtmp);

	in_plane_latt << 1, -1;
	bigtmp = supercell(tnn1m10, Z, Z, Z, Z, Z);
	pos.emplace_back(in_plane_latt);
	Ham.emplace_back(bigtmp);

	in_plane_latt << -1, 1;
	bigtmp = supercell(tnnm110, Z, Z, Z, Z, Z);
	pos.emplace_back(in_plane_latt);
	Ham.emplace_back(bigtmp);

	double x,z;
	int n = 1600; // The number of k-points used per meridian (half)
	vector<int> KK, LL; // These keep track of l, n in the loop after
	vector<vector<double>> kHam;
	SelfAdjointEigenSolver<bigM> CA;
	bigM tmp;
	Matrix<double, 45, 1> eigs;
	vector<double> Eigs;

	/* string Mydata = "bands.dat"; */
	/* ofstream Myfile; */	
	/* Myfile.open( Mydata.c_str(),ios::trunc ); */
	/* double k_x, k_y, k_z, b, pi; */
	/* b = 1.; */
	/* for (int k = 0; k < 351; k++) */
	/* { */
	/* 	if (k < 101){ */
	/* 		pi = 2*M_PI*k/100.; */
	/* 		k_x = pi/b; */
	/* 		k_y = 0; */
	/* 		k_z = 0; */
	/* 	} */
	/* 	if ((k > 100) && (k < 201)){ */
	/* 		pi = M_PI*(k-100)/100.; */
	/* 		k_x = (2*M_PI-pi)/b; */
	/* 		k_y = pi/b; */
	/* 		k_z = pi/b; */
	/* 	} */	
	/* 	if ((k > 200) && (k < 251)){ */
	/* 		pi = M_PI*(k-200)/50.; */
	/* 		k_x = (M_PI - pi)/b; */
	/* 		k_y = M_PI/b; */
	/* 		k_z = M_PI/b; */
	/* 	} */
	/* 	if ((k > 250) && (k < 351)){ */
	/* 		pi = M_PI*(k-250)/100.; */
	/* 		k_x = 0; */
	/* 		k_y = (M_PI-pi)/b; */
	/* 		k_z = (M_PI-pi)/b; */
	/* 	} */
	/* 	tmp = H_k(k_x, k_y, Ham, pos); */
	/* 	CA.compute(tmp, false); */
	/* 	eigs = CA.eigenvalues(); */
	/* 	for (int ll = 0; ll < eigs.size(); ll++) */
	/* 		Myfile<<k<<" "<<eigs(ll)<<endl; */
	/* } */
	/* Myfile.close(); */

	cout<<endl;
	cout<<"Building a table of k-space Hamiltonians"<<endl;
	//more efficient to build k_hams out of integration loop
	//the loop utilises the suggestions of
	//https://math.stackexchange.com/questions/3086903/trace-of-a-the-inverse-of-a-complex-matrix
	//Tr(D_w)^{-1} where D_w is the eigenvalues of -H + w
	for (int k = 0; k!=n+1; k++){
		if (k%2!=0){
			x = M_PI*k/n;
			for (int l = 0; l!=k+1; l++){
				if (l%2!=0){
					z = M_PI*l/n;
					tmp = H_k(x, z, Ham, pos);
		  			CA.compute(tmp, false);
		  			eigs = CA.eigenvalues();
					for (int aa = 0; aa < eigs.size(); aa++)
						Eigs.emplace_back(eigs(aa));
					kHam.emplace_back(Eigs);
					Eigs.clear();
					KK.emplace_back(k);
					LL.emplace_back(l);
				}
			}
		}
		if (k%20 == 0)
			cout<<"     "<<k/(n*1.+1)*100.<<"% completed"<<endl;
	}
	cout<<"     100% completed"<<endl;
	cout<<"Done!"<<endl;

	// The following block checks whether the hamiltonian is hermitian
	/* Matrix<dcomp, 45, 45> mtmp; */
	/* mtmp = H_k(0.22, -2.7, Ham, pos); */
	/* cout<<mtmp-mtmp.adjoint()<<endl; */

	double rho;
	string Mydata = "SK_DOS_up_quick.dat";
	ofstream Myfile;	
	Myfile.open( Mydata.c_str(),ios::trunc );
	cout<<endl;
	cout<<"Calculating PDOS"<<endl; // perform integration over k_|| graph as a function of energy
	vector<double> vtmp;
	vtmp.reserve(kHam[0].size());
	int counter = 0; // This to give percentage output
	for (double eps = 0.; eps < 24.1; eps += 0.05){
		rho = 0.;
		for (int k = 0; k < KK.size(); k++){
			vtmp = kHam[k];
			if ((KK[k]==1) && (LL[k]==1))
				rho += 0.5*Greens(eps, vtmp);
			else{
				if (KK[k]==LL[k]){
					rho += 0.5*Greens(eps, vtmp);
				}
				else
					rho += Greens(eps, vtmp);
			}
		}
		rho *= (-8./M_PI)/(n*n);
		Myfile<<eps<<" "<<rho<<endl;
		counter++;
		if (counter%10 == 0)
			cout<<"     "<<((eps)/24.)*100<<"% complete"<<endl;
	}
	Myfile.close();
	cout<<"     100% completed"<<endl;
	cout<<"Done!"<<endl;

	return 0;
}
