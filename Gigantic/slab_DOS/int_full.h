#ifndef INT_FULL_H
#define INT_FULL_H
//This has been edited as it wouldn't accept type double previously
//it will no longer take a general type, as a result
//routine verified with int_{-pi}^{pi} x^2 + y^2
//note integration is divided by 4pi^2

#include <cmath>
#include <eigen3/Eigen/Dense>
//this is the first TRUE cunningham points algorithm I have made
//with automatic convergence checks. Parameters include the max
//number of points in one meridian (iterations), the relative 
//error value, and the number of points in one meridian to start
//with. 27-06-17

using namespace std;

double Condition(std::complex<double> result, std::complex<double> part_result){
	double condition;
	if ((std::abs(part_result) == 0) && (std::abs(result) == 0))
		condition = 0;
	else
		condition = std::abs(std::abs(9.*part_result)/std::abs(result)-1.);
	return condition;
}
double Condition(double result, double part_result){
	double condition;
	if ((std::abs(part_result) == 0) && (std::abs(result) == 0))
		condition = 0;
	else
		condition = std::abs(std::abs(std::abs(9.*part_result)/(result))-1.);
	return condition;
}
double Condition(const Eigen::VectorXd &result, const Eigen::VectorXd &part_result){
	double condition;
	if ((std::abs(part_result.sum()) == 0) && (std::abs(result.sum()) == 0))
		condition = 0;
	else
		condition = std::abs(std::abs(9.*part_result.sum())/std::abs(result.sum())-1.);
	return condition;
}
double Condition(const Eigen::VectorXcd &result, const Eigen::VectorXcd &part_result){
	double condition;
	if ((std::abs(part_result.sum()) == 0) && (std::abs(result.sum()) == 0))
		condition = 0;
	else
		condition = std::abs(std::abs(9.*part_result.sum())/std::abs(result.sum())-1.);
	return condition;
}

template <typename func, typename... Args>
double kspace(func&& predicate, int iterations, double rel, int start, const double a, Args&&... params){

	double error;
	if (rel == 0)
		error = 0.002;
	else
		error = rel;
	int max_width;
	if (iterations == 0)
		max_width = 400;
	else
		max_width = 2*iterations; 	//2 here to cater for the doubling of N

	int N; 	//starting number of points in one meridian - editable
	if (start == 0)
		N = 1;
	else
		N = start;

	int n = 2*N;	//this caters for the special spacing but see conditions for odd numbers 
			//below. This ensures number of points in one meridian is indeed n
	
	/* string Mydata = "test.txt"; */
	/* ofstream Myfile; */	
	/* Myfile.open( Mydata.c_str(),ios::trunc ); */

	const double A = M_PI/a;
	double x,z;
	double integral;
	double tmp;
	double rel_error;
	integral = 0.;

	for (int k = -n; k!=n+1; k++){
		if (k%2!=0){
			x = A*k/n;
			for (int l = -n; l!=n+1; l++){
				if (l%2!=0){
					z = A*l/n;
					/* Myfile<<x<<" , "<<z<<endl; */
					/* if ((k==1) && (l==1)) */
					/* 	integral = 0.5*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...); */
					/* else{ */
					/* 	if (k==l){ */
							/* integral += 0.5*predicate(x,z,forward<T>(params)); */
							/* integral = integral + 0.5*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...); */
						/* } */
						/* else */
							integral = integral + std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);
					/* } */
				}
			}
		}
	}
	
	int Nnew;
	int condition = -2;
	double foo = 0;

	while (condition != 1){
		Nnew = 3*n;
		tmp = integral;
		for (int k = -Nnew; k!=Nnew+1; k++){
			if (k%2!=0){
				x = A*k/Nnew;
				for (int l = -Nnew; l!=Nnew+1; l++){
					if ((l%2!=0) && ((k%3!=0) || (l%3!=0))){
						z = A*l/Nnew;
						/* Myfile<<x<<" , "<<z<<endl; */
						/* if (k==l){ */
						/* 	/1* integral += 0.5*predicate(x,z,forward<T>(params)); *1/ */
						/* 	integral = integral + 0.5*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...); */
						/* } */
						/* else{ */
							/* integral += predicate(x,z,forward<T>(params)); */
							integral = integral + std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);
						/* } */
					}
				}
			}
		}
		n = Nnew;

		/* if (std::abs(integral) < 1e-25) */
		/* 	integral = 0.; */
		/* if (Nnew>64 && integral == 0. && tmp == 0.) */
		/* 	condition = 1; */

		//256 seems to be the magic number for the lower bound of sampling points
		//3% margin of error as below works well when benchmarked against fixed method 
		//of double Simpson's integral of rho-E/LDOS.cpp
		rel_error = Condition(integral, tmp);//condition overloaded for all data types this header is designed for
		if (Nnew > max_width){
			cout<<"error, maximum number of iteration points reached. Estimated error = "<<rel_error<<endl;
			condition = 1;
		}
		if (Nnew>=128 && rel_error <= error)
			condition = 1;
	}
	integral = integral*(1./(Nnew*Nnew));

	return integral;
}
#endif
