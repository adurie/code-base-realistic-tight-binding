#include <cmath>
#define EIGEN_USE_MKL_ALL
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <iomanip>
#include <fstream>
#include "int_full.h"

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
typedef Vector2d v2;
// TODO make Matrix sizes more general - 45x45 = (9x5)x(9x5) i.e 9 bands, 5 atoms, dependent on input file
typedef vector<Matrix<dcomp, 45, 45>, aligned_allocator<Matrix<dcomp, 45, 45>>> vm;
typedef vector<v2, aligned_allocator<v2>> vv;
typedef Matrix<dcomp, 45, 45> bigM;

//Purpose: To calculate the partial DOS. Current system, slab geometry with 5 atom unit cell
//but *hr.dat decides system.

Matrix<dcomp, 45, 45> read(Vector2d &dvec, int ispin){
      // read in Hamiltonians from *hr.dat
      string input;
      // TODO make more general
      if (ispin == +1)
        input = "cobalt_up_hr.dat";
      if (ispin == -1)
        input = "cobalt_dn_hr.dat";
      ifstream infile(input);
      Matrix<dcomp, 45, 45> rt;
      rt.fill(0.);
      Vector2d a_1, a_2;
      a_1 << 1, 0; // reinterpret *.hr.dat in terms of lattice vectors
      a_2 << 0, 1; // here 2D simple cubic
      Vector2d A;
  
      dcomp i;
      i = -1;
      i = sqrt(i);
      string line;
      double a, b, c, d, e, f, g;
      /* double eV_Ry = 0.073498618; // convert to Ryds */
      double eV_Ry = 1; // keep as eV for now
      for (int ii=0; ii<9; ii++) // ignore file header TODO make more general
      	infile.ignore(100,'\n');
      while (!infile.eof()) 
      {
	getline(infile, line);
	istringstream iss(line);
	iss >> a >> b >> c >> d >> e >> f >> g;
	A = a*a_1 + b*a_2;
	if ((A(0) == dvec(0)) && (A(1) == dvec(1))){
		rt(d-1,e-1) = (f + g*i)*eV_Ry;
		if ((d==45) && (e==45)) // no point in reading after NxNth element TODO make more general
			break;
	}
	else{
		for (int ii=0; ii<2024; ii++) // if (0,0) element doesn't match, skip by NxN TODO make more general
			infile.ignore(100,'\n');

	}
	
      }
      return rt;
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

	Vector2d lat_1, lat_2;
	lat_1 << 1, 0; // This and the below loop to generate 
	lat_2 << 0, 1; // position of each neighbour
	Vector2d result;
	vm Ham; // Dimensions based on number of bands & atoms in unit cell. TODO generalise
	vv pos; // by extracting information from *hr.dat. pos to contain all position vectors.
	int k = 0;
	cout<<"Reading in Wannier90 Hamiltonians"<<endl;
	for (int i = -4; i < 5; i++){ // i & j range taken from *hr.dat. TODO generalise
		for (int j = -4; j < 5; j++){
			result = i*lat_1 + j*lat_2;
			pos.emplace_back(result);
			Ham.emplace_back(read(result, ispin));
		}
		k++;
		cout<<"     "<<(k/9.)*100<<"% completed"<<endl; // TODO fix % so more general
	}
	cout<<"Done!"<<endl;
	double x,z;
	int n = 1600; // The number of k-points used per meridian (half)
	vector<int> KK, LL; // These keep track of l, n in the loop after
	vector<vector<double>> kHam;
	SelfAdjointEigenSolver<bigM> CA;
	bigM tmp;
	Matrix<double, 45, 1> eigs;
	vector<double> Eigs;
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
	/* cout<<fixed<<mtmp<<endl; */

	double rho;
	string Mydata = "DOS_quick.dat";
	ofstream Myfile;	
	Myfile.open( Mydata.c_str(),ios::trunc );
	cout<<endl;
	cout<<"Calculating PDOS"<<endl; // perform integration over k_|| graph as a function of energy
	vector<double> vtmp;
	vtmp.reserve(kHam[0].size());
	int counter = 0; // This to give percentage output
	for (double eps = -12.; eps < 12.1; eps += 0.05){
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
			cout<<"     "<<((eps+12)/24.)*100<<"% complete"<<endl;
	}
	Myfile.close();
	cout<<"     100% completed"<<endl;
	cout<<"Done!"<<endl;

	return 0;
}
