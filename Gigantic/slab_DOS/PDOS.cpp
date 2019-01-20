#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
typedef Vector2d v2;
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

Matrix<dcomp, 45, 45> H_k(vm &Ham, vv &pos, double k_x, double k_y){
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

int main(){
	const int ispin = -1; // TODO shoudn't this be determined at runtime?

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
	// Calculate G(eps, k||) = (eps - H(k||) + im*delta)^(-1)
	bigM I = bigM::Identity();
	dcomp im;
	im = -1.;
	im = sqrt(im);
	double eps = 3.;
	const double delta = 1e-5;
	bigM G, Gtmp, tmp;
	tmp = H_k(Ham, pos, M_PI, M_PI);
	Gtmp = I*(eps + delta*im) - tmp;
	G = Gtmp.inverse(); // This is the power of eigen!
	cout<<G.topLeftCorner(9,9).real();
	return 0;
}
