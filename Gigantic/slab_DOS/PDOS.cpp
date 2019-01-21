#include <cmath>
#define EIGEN_USE_MKL_ALL
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <iomanip>
#include <fstream>
#include "cunningham.h"

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
double Greens(double k_x, double k_y, double a, double eps, vm &Ham, vv &pos){ //Args&&... params){
	// 'a' still here for legacy reasons
	// Calculate G(eps, k||) = (eps - H(k||) + im*delta)^(-1)
	// TODO make this perform a trace over a subset of G, 
	// to calculate PDOS
	bigM I = bigM::Identity();
	dcomp im;
	im = -1.;
	im = sqrt(im);
	const double delta = 1e-5;
	bigM G, Gtmp, tmp;
	tmp = H_k(k_x, k_y, Ham, pos);//forward<Args>(params)...);
	Gtmp = I*(eps + delta*im) - tmp;
	G = Gtmp.inverse(); // This is the power of eigen!
	double rho; //DOS, rho = (-1/pi)*ImTrG
	rho = static_cast<double>(G.imag().trace());
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
		cout<<(k/9.)*100<<"% completed"<<endl; // TODO fix % so more general
	}
	cout<<"Done!"<<endl;
	double rho;
	//see integration header for what the following zero values are interpreted as
	int k_start = 0;
	int k_max = 0;
	double abs_error = 1e-1;
	const double a = 1; // for legacy reasons only
	string Mydata = "DOS.dat";
	ofstream Myfile;	
	Myfile.open( Mydata.c_str(),ios::trunc );
	cout<<"Calculating PDOS"<<endl; // perform integration over k_|| graph as a function of energy
	for (double eps = -12.; eps < 12.1; eps += 0.1){
		rho = kspace(&Greens, k_start, abs_error, k_max, a, eps, Ham, pos);
		rho *= (-1./M_PI);
		Myfile<<eps<<" "<<rho<<endl;
		cout<<((eps+12)/24.)*100<<"% complete"<<endl;
	}
	Myfile.close();
	/* rho = Greens(3., Ham, pos, M_PI, M_PI); */
	return 0;
}
