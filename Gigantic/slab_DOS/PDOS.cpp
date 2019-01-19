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

Matrix<dcomp, 45, 45> read(Vector2d &dvec, int ispin){
      string input;
      if (ispin == +1)
        input = "cobalt_up_hr.dat";
      if (ispin == -1)
        input = "cobalt_dn_hr.dat";
      ifstream infile(input);
      Matrix<dcomp, 45, 45> rt;
      rt.fill(0.);
      Vector2d a_1, a_2;
      a_1 << 1, 0;
      a_2 << 0, 1;
      Vector2d A;
  
      dcomp i;
      i = -1;
      i = sqrt(i);
      string line;
      double a, b, c, d, e, f, g;
      /* double eV_Ry = 0.073498618; */
      double eV_Ry = 1;
      for (int ii=0; ii<9; ii++)
      	infile.ignore(100,'\n');
      while (!infile.eof()) 
      {
	getline(infile, line);
	istringstream iss(line);
	iss >> a >> b >> c >> d >> e >> f >> g;
	A = a*a_1 + b*a_2;
	if ((A(0) == dvec(0)) && (A(1) == dvec(1))){
		rt(d-1,e-1) = (f + g*i)*eV_Ry;
		if ((d==45) && (e==45))
			break;
	}
	else{
		for (int ii=0; ii<2024; ii++)
			infile.ignore(100,'\n');

	}
	
      }
      return rt;
}

Matrix<dcomp, 45, 45> H_k(vm &Ham, vv &pos, double k_x, double k_y){
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
	int ispin = -1;

	Vector2d lat_1, lat_2;
	lat_1 << 1, 0;
	lat_2 << 0, 1;
	Vector2d result;
	vm Ham;
	vv pos;
	Matrix<dcomp, 45, 45> dummy;
	int k = 0;
	cout<<"Reading in Wannier90 Hamiltonians"<<endl;
	for (int i = -4; i < 5; i++){
		for (int j = -4; j < 5; j++){
			result = i*lat_1 + j*lat_2;
			pos.emplace_back(result);
			Ham.emplace_back(read(result, ispin));
		}
		k++;
		cout<<(k/9.)*100<<"% completed"<<endl;
	}
	bigM I = bigM::Identity();
	dcomp im;
	im = -1.;
	im = sqrt(im);
	double eps = 3.;
	double delta = 1e-5;
	bigM G, Gtmp, tmp;
	tmp = H_k(Ham, pos, M_PI, M_PI);
	Gtmp = I*(eps + delta*im) - tmp;
	G = Gtmp.inverse();
	cout<<G.topLeftCorner(9,9).real();
	return 0;
}
