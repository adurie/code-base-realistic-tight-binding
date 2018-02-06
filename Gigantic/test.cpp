#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <eigen3/Eigen/StdVector>
#include <gsl/gsl_math.h>

using namespace Eigen;
using namespace std;
// define a custom template binary functor

/* template<typename Scalar> struct MakeComplexOp { */
/*   EIGEN_EMPTY_STRUCT_CTOR(MakeComplexOp) */
/*   typedef complex<Scalar> result_type; */
/*   complex<Scalar> operator()(const Scalar& a, const Scalar& b) const { return complex<Scalar>(a,b); } */
/* }; */

typedef complex<double> dcomp;
bool WayToSort(double i, double j){ return abs(i) < abs(j);}

struct gmean {
	dcomp operator()(const dcomp& x, const dcomp& y) const {
		dcomp gmean;
		if (real(x*y) > 0){
			if (real(x) < 0)
				gmean=-sqrt(x*y);
			if (real(x) > 0)
			gmean = sqrt(x*y);
      		}
		else
       		gmean=(x+y)/2.;
     		return gmean;
      	}
};

Vector2d prestructij(int i, int j, int numn){
	Vector2d dist;
	if (i == j)
		dist << 3, 6;
	else
		dist << 2, 7;
	return dist;
}

int main()
{

	/* MatrixXcd m1(2,2), m2(2,2); */
	/* m1 << 0.5, -3.42, -0.883, 2.333; */
	/* m2 << -0.81, -3.42, 5.883, 2.833; */
	/* cout << m1.binaryExpr(m2, gmean()) << endl; */
	/* vector<MatrixXcd, aligned_allocator<MatrixXcd>> big; */
	/* big.reserve(2); */
	/* big.emplace_back(m1); */
	/* big.emplace_back(m2); */
	/* for (int it = 0; it <=1 ; it++) */
	/* 	cout<<big[it](1,1)<<endl; */

	/* int nlay = 7; */
	/* int numn = 2; */
      /* MatrixXd mdist(nlay-1, nlay-1); */
      /* mdist.fill(0); */
      /* Vector2d dist; */
      /* vector<MatrixXd, aligned_allocator<MatrixXd>> dnn; */
      /* dnn.reserve(numn); */
      /* for (int i=1; i <= numn; i++){ */
        /* for (int ilay=2; ilay <= nlay; ilay++){ */
          /* dist = prestructij(ilay,ilay,numn); */
          /* mdist(ilay-2,ilay-2)=dist(i-1); */
          /* dist = prestructij(ilay,ilay-1,numn); */
	/*   if (ilay > 2){ */
            /* mdist(ilay-2,ilay-3)=dist(i-1); */
            /* mdist(ilay-3,ilay-2)=dist(i-1); */
	/*   } */
	/* } */
	/* dnn.emplace_back(mdist); */
      /* } */
      /* cout<<dnn[0]<<endl<<endl<<dnn[1]<<endl; */
      /* cout<<gsl_hypot3(3,4,5)<<endl; */
      /* cout<<sqrt(50)<<endl; */

	/* Vector3d test; */
	/* test << 1, 3, 8; */
	/* cout<<test.hypotNorm()<<endl; */
	/* cout<<sqrt(1+9+64)<<endl; */
	/* test.normalize(); */
	/* cout<<test<<endl; */

/*       double xnorm = test.dot(test); */
/*         test = test/sqrt(xnorm); */
/* 	cout<<test<<endl; */


	/* vector<double> test; */
	/* test.reserve(4); */
	/* for (int i = 0; i < 4; i++) */
	/* 	test.emplace_back(3.14); */
	/* test[3] = 54.2; */
	/* for (int i = 0; i < 4; i++) */
	/* 	cout<<test[i]<<endl; */

	/* vector<double> test = {53.43, -52.28, -54.23, 0.32, -1.5}; */
      /* sort(test.begin(), test.end(), WayToSort); */
      /* for (double i:test) */
	/*       cout<<i<<endl; */
	
	Matrix3d test1, test2;
	test1 << 1,2,3,1,2,3,1,2,3;
	test2 << 4,5,6,4,5,6,4,5,6;
	/* cout<<test(0,1)<<endl; */
	/* Vector3d tmp, tmp2; */
	/* tmp << 12,13,14; */
	/* tmp2 << 32, 34, 35; */
	/* test.block<3,1>(0,1) = tmp; */
	/* test.rightCols(1) = tmp2; */
	cout<<test2*test1<<endl;


}
