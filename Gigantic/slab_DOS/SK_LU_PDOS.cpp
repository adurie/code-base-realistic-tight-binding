#include <cmath>
#define EIGEN_USE_MKL_ALL
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <iomanip>
#include <fstream>
#include "TB_Co_Cu.h"
// Orthogonal matrix: inverse is equal to it's conjugate transpose. The eigenvectors of a Hermitian matrix are orthogonal
// we will be exploiting this fact

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
typedef Vector2d v2;
typedef Vector3d v3;
// TODO make Matrix sizes more general - 45x45 = (9x5)x(9x5) i.e 9 bands, 5 atoms, dependent on input file
typedef vector<Matrix<dcomp, 45, 45>, aligned_allocator<Matrix<dcomp, 45, 45>>> vm;
typedef vector<v2, aligned_allocator<v2>> vv;
typedef vector<Matrix<double, 45, 1>, aligned_allocator<Matrix<double, 45, 1>>> VV;
typedef Matrix<dcomp, 45, 45> bigM;
typedef Matrix<dcomp, 9, 9> smallM;

//Purpose: To calculate the partial DOS. Current system, slab geometry with 5 atom unit cell
//but *hr.dat decides system.

bigM supercell(smallM &ii, smallM &iip1, smallM &iip2, smallM &ip1i, smallM &ip2i, smallM &Z){
	bigM H;
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

double Greens(double eps, bigM &tvecs, Matrix<double, 45, 1> &teigs, double &rho1, double &rho2, double &rho3, double &rho4, double &rho5){
	// Calculate G(eps, k||) = (eps - H(k||) + im*delta)^(-1)
	// TODO make this perform a trace over a subset of G, 
	// to calculate PDOS
	bigM adjvecs;
	adjvecs = tvecs.adjoint();
	dcomp im;
	im = -1.;
	im = sqrt(im);
	const double delta = 1e-5;
	Matrix<dcomp, 45, 1> lambda;
	lambda = teigs.cast<dcomp>();
	lambda = lambda.array() - (eps + delta*im);
	lambda = lambda.array().inverse();
	bigM G, Gtmp;
	Gtmp = lambda.asDiagonal();
	G = tvecs*Gtmp*adjvecs;
	double rho; //DOS, rho = (-1/pi)*ImTrG
	rho = static_cast<double>(G.imag().trace());
	rho1 = 0.; rho2 = 0.; rho3 = 0.; rho4 = 0.; rho5 = 0.;
	for (int k = 0; k < 9; k++){
		rho1 -= imag(G(k, k));
		rho2 -= imag(G(k + 9, k + 9));
		rho3 -= imag(G(k + 18, k + 18));
		rho4 -= imag(G(k + 27, k + 27));
		rho5 -= imag(G(k + 36, k + 36));
	}
	return -rho;
}

int main(){
	const int ispin = -1; // TODO shoudn't this be determined at runtime?

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
	int n = 500; // The number of k-points used per meridian (half)
	vector<int> KK, LL; // These keep track of l, n in the loop after
	vm vecs;
	VV lams;
	SelfAdjointEigenSolver<bigM> CA;
	bigM tmp;
	bigM eigv;
	Matrix<double, 45, 1> eigs;
	cout<<endl;
	cout<<"Building a table of k-space Hamiltonians"<<endl;
	//more efficient to build k_hams out of integration loop
	//the loop utilises the suggestions of
	//https://math.stackexchange.com/questions/3086903/trace-of-a-the-inverse-of-a-complex-matrix
	//Tr(D_w)^{-1} where D_w is the eigenvalues of -H + w

	/* //This for DOS at fixed k_|| */
	/* tmp = H_k(2.53, 0, Ham, pos); */
	/* CA.compute(tmp); */
	/* eigs = CA.eigenvalues(); */
	/* eigv = CA.eigenvectors(); */

	//This for integrated DOS in k||
	for (int k = 0; k!=n+1; k++){
		if (k%2!=0){
			x = M_PI*k/n;
			for (int l = 0; l!=k+1; l++){
				if (l%2!=0){
					z = M_PI*l/n;
					tmp = H_k(x, z, Ham, pos);
		  			CA.compute(tmp);
		  			eigs = CA.eigenvalues();
					eigv = CA.eigenvectors();
					vecs.emplace_back(eigv);
					lams.emplace_back(eigs);
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
	/* cout<<fixed<<mtmp<<endl; */

	double rho;
	double rhop1, rhop2, rhop3, rhop4, rhop5;
	double trhop1, trhop2, trhop3, trhop4, trhop5;
	string Mydata = "SK_DOS_dn_full.dat";
	ofstream Myfile, Myfile1, Myfile2, Myfile3, Myfile4, Myfile5;	
	Myfile.open( Mydata.c_str(),ios::trunc );
	Mydata = "SK_PDOS_l1_dn.dat";
	Myfile1.open( Mydata.c_str(),ios::trunc );
	Mydata = "SK_PDOS_l2_dn.dat";
	Myfile2.open( Mydata.c_str(),ios::trunc );
	Mydata = "SK_PDOS_l3_dn.dat";
	Myfile3.open( Mydata.c_str(),ios::trunc );
	Mydata = "SK_PDOS_l4_dn.dat";
	Myfile4.open( Mydata.c_str(),ios::trunc );
	Mydata = "SK_PDOS_l5_dn.dat";
	Myfile5.open( Mydata.c_str(),ios::trunc );
	cout<<endl;
	cout<<"Calculating PDOS"<<endl; // perform integration over k_|| graph as a function of energy
	bigM tvecs;
	Matrix<double, 45, 1> teigs;
	double factor = (-8./M_PI)/(n*n);
	int counter = 0; // This to give percentage output
	for (double eps = 0.; eps < 24.1; eps += 0.05){
		rho = 0.;
		rhop1 = 0.; rhop2 = 0.; rhop3 =0.; rhop4 = 0.; rhop5 = 0.;

		/* //This for DOS at fixed k_|| */
		/* rho += Greens(eps, eigv, eigs, trhop1, trhop2, trhop3, trhop4, trhop5); */
		/* rhop1 += trhop1; rhop2 += trhop2; rhop3 += trhop3; rhop4 += trhop4; rhop5 += trhop5; */

		//This for integrated DOS in k||
		for (int k = 0; k < KK.size(); k++){
			tvecs = vecs[k];
			teigs = lams[k];
			if ((KK[k]==1) && (LL[k]==1)){
				rho += 0.5*Greens(eps, tvecs, teigs, trhop1, trhop2, trhop3, trhop4, trhop5);
				rhop1 += 0.5*trhop1; rhop2 += 0.5*trhop2; rhop3 += 0.5*trhop3; rhop4 += 0.5*trhop4; rhop5 += 0.5*trhop5;
			}
			else{
				if (KK[k]==LL[k]){
					rho += 0.5*Greens(eps, tvecs, teigs, trhop1, trhop2, trhop3, trhop4, trhop5);
					rhop1 += 0.5*trhop1; rhop2 += 0.5*trhop2; rhop3 += 0.5*trhop3; rhop4 += 0.5*trhop4; rhop5 += 0.5*trhop5;
				}
				else{
					rho += Greens(eps, tvecs, teigs, trhop1, trhop2, trhop3, trhop4, trhop5);
					rhop1 += trhop1; rhop2 += trhop2; rhop3 += trhop3; rhop4 += trhop4; rhop5 += trhop5;
				}
			}
		}

		rho *= factor;
		rhop1 *= factor;
		rhop2 *= factor;
		rhop3 *= factor;
		rhop4 *= factor;
		rhop5 *= factor;
		Myfile<<eps<<" "<<rho<<endl;
		Myfile1<<eps<<" "<<rhop1<<endl;
		Myfile2<<eps<<" "<<rhop2<<endl;
		Myfile3<<eps<<" "<<rhop3<<endl;
		Myfile4<<eps<<" "<<rhop4<<endl;
		Myfile5<<eps<<" "<<rhop5<<endl;
		counter++;
		if (counter%10 == 0)
			cout<<"     "<<((eps)/24.)*100<<"% complete"<<endl;
	}
	Myfile.close();
	Myfile1.close(); Myfile2.close(); Myfile3.close(); Myfile4.close(); Myfile5.close();
	cout<<"     100% completed"<<endl;
	cout<<"Done!"<<endl;

	return 0;
}
