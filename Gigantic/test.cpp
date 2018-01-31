#include <eigen3/Eigen/Core>
#include <iostream>
using namespace Eigen;
using namespace std;
// define a custom template binary functor

/* template<typename Scalar> struct MakeComplexOp { */
/*   EIGEN_EMPTY_STRUCT_CTOR(MakeComplexOp) */
/*   typedef complex<Scalar> result_type; */
/*   complex<Scalar> operator()(const Scalar& a, const Scalar& b) const { return complex<Scalar>(a,b); } */
/* }; */

typedef complex<double> dcomp;

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

int main()
{
	MatrixXcd m1(2,2), m2(2,2);
	m1 << 0.5, -3.42, -0.883, 2.333;
	m2 << -0.81, -3.42, 5.883, 2.833;
	cout << m1.binaryExpr(m2, gmean()) << endl;
	return 0;
}
