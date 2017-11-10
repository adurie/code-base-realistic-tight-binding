#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

double f(double x, double y, double z){
	return 1;
}

int main(){

	int N = 40;
	double x, y, z;
	double a = 1;
	double A = M_PI/a;
	int n = 2*N;
	double integral = 0;

	for (int k = 0; k!=n+1; k++){
		if (k%2!=0){
			x = A*k/n;
			for (int l = 0; l!=k+1; l++){
				if (l%2!=0){
					y = A*l/n;
					for (int m = 0; m!=l+1; m++){
						if (m%2!=0){
							z = A*m/n;
							if ((k==l) && (k==m) && (l==m)){
								/* integral = integral + 0.5*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...); */
								integral = integral + (1/6.)*f(x,y,z);
								/* f(x,y,z); */
							}
							else if ((k==l) || (k==m) || (l==m)){
								/* integral = integral + 0.5*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...); */
								integral = integral + 0.5*f(x,y,z);
								/* f(x,y,z); */
							}
							else
							{
								integral = integral + f(x,y,z);
								/* integral = integral + std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...); */
								/* f(x,y,z); */
							}
						}
					}
				}
			}
		}
	}
	integral = (48.*A*A*A/(N*N*N))*integral;
	cout<<integral<<endl;
	
	return 0;
}
