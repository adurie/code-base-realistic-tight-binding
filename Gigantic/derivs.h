#ifndef DTB_H
#define DTB_H

#include <cmath>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;

//these derivatives are actually dH/d***, and are calculated in the Maple
//worksheet sk_deriv.mw
Matrix<dcomp, 9, 9> dsss(){
	Matrix<dcomp, 9, 9> b;
	b.fill(0.);
	b(0,0) = 1.;
	return b;
}

Matrix<dcomp, 9, 9> dsps(const Vector3d &pos){
	Matrix<dcomp, 9, 9> b;
	b.fill(0.);
	double x, y, z;
	Vector3d X, Y, Z;
	X << 1, 0, 0;
	Y << 0, 1, 0;
	Z << 0, 0, 1;
	x = pos.dot(X)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	y = pos.dot(Y)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	z = pos.dot(Z)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	b(0,1) = x;
	b(0,2) = y;
	b(0,3) = z;
	b(1,0) = -x;
	b(2,0) = -y;
	b(3,0) = -z;
	return b;
}

Matrix<dcomp, 9, 9> dpps(const Vector3d &pos){
	Matrix<dcomp, 9, 9> b;
	b.fill(0.);
	double x, y, z;
	Vector3d X, Y, Z;
	X << 1, 0, 0;
	Y << 0, 1, 0;
	Z << 0, 0, 1;
	x = pos.dot(X)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	y = pos.dot(Y)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	z = pos.dot(Z)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
        double xx=x*x;
        double xy=x*y;
        double yy=y*y;
        double yz=y*z;
        double zz=z*z;
        double zx=z*x;
	b(1,1) = xx;
	b(2,2) = yy;
	b(3,3) = zz;
	b(1,2) = xy;
	b(1,3) = zx;
	b(2,1) = xy;
	b(2,3) = yz;
	b(3,1) = zx;
	b(3,2) = yz;
	return b;
}

Matrix<dcomp, 9, 9> dppp(const Vector3d &pos){
	Matrix<dcomp, 9, 9> b;
	b.fill(0.);
	double x, y, z;
	Vector3d X, Y, Z;
	X << 1, 0, 0;
	Y << 0, 1, 0;
	Z << 0, 0, 1;
	x = pos.dot(X)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	y = pos.dot(Y)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	z = pos.dot(Z)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
        double xx=x*x;
        double xy=x*y;
        double yy=y*y;
        double yz=y*z;
        double zz=z*z;
        double zx=z*x;
	b(1,1) = 1 - xx;
	b(2,2) = 1 - yy;
	b(3,3) = 1 - zz;
	b(1,2) = -xy;
	b(1,3) = -zx;
	b(2,1) = -xy;
	b(2,3) = -yz;
	b(3,1) = -zx;
	b(3,2) = -yz;
	return b;
}

Matrix<dcomp, 9, 9> dsds(const Vector3d &pos){
	double r3 = sqrt(3.);
	Matrix<dcomp, 9, 9> b;
	b.fill(0.);
	double x, y, z;
	Vector3d X, Y, Z;
	X << 1, 0, 0;
	Y << 0, 1, 0;
	Z << 0, 0, 1;
	x = pos.dot(X)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	y = pos.dot(Y)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	z = pos.dot(Z)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	b(0,4) = r3*x*y;
	b(0,5) = r3*y*z;
	b(0,6) = r3*z*x;
	b(0,7) = 0.5*r3*(x*x - y*y);
	b(0,8) = 0.5*(3*z*z - 1);
	b(4,0) = b(0,4);
	b(5,0) = b(0,5);
	b(6,0) = b(0,6);
	b(7,0) = b(0,7);
	b(8,0) = b(0,8);
	return b;
}

Matrix<dcomp, 9, 9> dpds(const Vector3d &pos){
	double r3 = sqrt(3.);
	Matrix<dcomp, 9, 9> b;
	b.fill(0.);
	double x, y, z;
	Vector3d X, Y, Z;
	X << 1, 0, 0;
	Y << 0, 1, 0;
	Z << 0, 0, 1;
	x = pos.dot(X)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	y = pos.dot(Y)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	z = pos.dot(Z)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
        double xx=x*x;
        double xy=x*y;
        double yy=y*y;
        double yz=y*z;
        double zz=z*z;
        double zx=z*x;
	b(1,4) = r3*xx*y;
	b(1,5) = r3*xy*z;
	b(1,6) = r3*zx*x;
	b(1,7) = 0.5*r3*(xx - yy)*x;
	b(1,8) = 0.5*(3*zz - 1)*x;
	b(4,1) = -b(1,4);
	b(5,1) = -b(1,5);
	b(6,1) = -b(1,6);
	b(7,1) = -b(1,7);
	b(8,1) = -b(1,8);
	b(2,4) = r3*x*yy;
	b(2,5) = r3*yy*z;
	b(2,6) = r3*zx*y;
	b(2,7) = 0.5*r3*(xx - yy)*y;
	b(2,8) = 0.5*(3*zz - 1)*y;
	b(4,2) = -b(2,4);
	b(5,2) = -b(2,5);
	b(6,2) = -b(2,6);
	b(7,2) = -b(2,7);
	b(8,2) = -b(2,8);
	b(3,4) = r3*xy*z;
	b(3,5) = r3*yz*z;
	b(3,6) = r3*zx*z;
	b(3,7) = 0.5*r3*(xx - yy)*z;
	b(3,8) = 0.5*(3*zz - 1)*z;
	b(4,3) = -b(3,4);
	b(5,3) = -b(3,5);
	b(6,3) = -b(3,6);
	b(7,3) = -b(3,7);
	b(8,3) = -b(3,8);
	return b;
}

Matrix<dcomp, 9, 9> dpdp(const Vector3d &pos){
	double r3 = sqrt(3.);
	Matrix<dcomp, 9, 9> b;
	b.fill(0.);
	double x, y, z;
	Vector3d X, Y, Z;
	X << 1, 0, 0;
	Y << 0, 1, 0;
	Z << 0, 0, 1;
	x = pos.dot(X)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	y = pos.dot(Y)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	z = pos.dot(Z)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
        double xx=x*x;
        double xy=x*y;
        double yy=y*y;
        double yz=y*z;
        double zz=z*z;
        double zx=z*x;
	b(1,4) = y*(1 - 2.*xx); 
	b(1,5) = -2*xy*z;
	b(1,6) = z*(1 - 2.*xx);
	b(1,7) = x*(1 - xx + yy);
	b(1,8) = -x*zz*r3;
	b(4,1) = -b(1,4);
	b(5,1) = -b(1,5);
	b(6,1) = -b(1,6);
	b(7,1) = -b(1,7);
	b(8,1) = -b(1,8);
	b(2,4) = x*(1 - 2.*yy); 
	b(2,5) = z*(1 - 2.*yy);
	b(2,6) = -2*xy*z;
	b(2,7) = y*(1 + xx - yy);
	b(2,8) = -yz*z*r3;
	b(4,2) = -b(2,4);
	b(5,2) = -b(2,5);
	b(6,2) = -b(2,6);
	b(7,2) = -b(2,7);
	b(8,2) = -b(2,8);
	b(3,4) = -2*xy*z;
	b(3,5) = y*(1 - 2.*zz);
	b(3,6) = x*(1 - 2.*zz); 
	b(3,7) = -z*(xx - yy);
	b(3,8) = z*r3*(xx + yy);
	b(4,3) = -b(3,4);
	b(5,3) = -b(3,5);
	b(6,3) = -b(3,6);
	b(7,3) = -b(3,7);
	b(8,3) = -b(3,8);
	return b;
}

Matrix<dcomp, 9, 9> ddds(const Vector3d &pos){
	double r3 = sqrt(3.);
	Matrix<dcomp, 9, 9> b;
	b.fill(0.);
	double x, y, z;
	Vector3d X, Y, Z;
	X << 1, 0, 0;
	Y << 0, 1, 0;
	Z << 0, 0, 1;
	x = pos.dot(X)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	y = pos.dot(Y)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	z = pos.dot(Z)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
        double xx=x*x;
        double xy=x*y;
        double yy=y*y;
        double yz=y*z;
        double zz=z*z;
        double zx=z*x;
        double xxyy=xx*yy;
        double yyzz=yy*zz;
        double zzxx=zz*xx;
	b(4,4) = 3*xxyy;
	b(4,5) = 3*yy*zx;
	b(4,6) = 3*xx*yz;
	b(4,7) = 0.5*xy*(3*xx - 3*yy);
	b(4,8) = 0.5*r3*xy*(2*zz - xx - yy);
	b(5,4) = b(4,5);
	b(6,4) = b(4,6);
	b(7,4) = b(4,7);
	b(8,4) = b(4,8);
	b(5,5) = 3*yyzz; 
	b(5,6) = 3*zz*xy;
	b(5,7) = 0.5*yz*(3*xx - 3*yy);
	b(5,8) = 0.5*r3*yz*(2*zz - xx - yy);
	b(6,5) = b(5,6);
	b(7,5) = b(5,7);
	b(8,5) = b(5,8);
	b(6,6) = 3*zzxx;
	b(6,7) = 0.5*zx*(3*xx - 3*yy);
	b(6,8) = 0.5*r3*zx*(2*zz - xx - yy);
	b(7,6) = b(6,7);
	b(8,6) = b(6,8);
	b(7,7) = 0.75*(xx - yy)*(xx - yy);
	b(7,8) = 0.25*r3*(xx - yy)*(2*zz - xx - yy);
	b(8,7) = b(7,8);
	b(8,8) = 0.25*(2*zz - xx - yy)*(2*zz - xx - yy);
	return b;
}

Matrix<dcomp, 9, 9> dddp(const Vector3d &pos){
	double r3 = sqrt(3.);
	Matrix<dcomp, 9, 9> b;
	b.fill(0.);
	double x, y, z;
	Vector3d X, Y, Z;
	X << 1, 0, 0;
	Y << 0, 1, 0;
	Z << 0, 0, 1;
	x = pos.dot(X)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	y = pos.dot(Y)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	z = pos.dot(Z)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
        double xx=x*x;
        double xy=x*y;
        double yy=y*y;
        double yz=y*z;
        double zz=z*z;
        double zx=z*x;
        double xxyy=xx*yy;
        double yyzz=yy*zz;
        double zzxx=zz*xx;
	b(4,4) = -4*xxyy + xx + yy;
	b(4,5) = zx*(1 - 4*yy);
	b(4,6) = yz*(1 - 4*xx);
	b(4,7) = 2*xy*(yy - xx);
	b(4,8) = -2*r3*zz*xy;
	b(5,4) = b(4,5);
	b(6,4) = b(4,6);
	b(7,4) = b(4,7);
	b(8,4) = b(4,8);
	b(5,5) = -4*yyzz + yy + zz;
	b(5,6) = xy*(1 - 4*zz);
	b(5,7) = yz*(2*yy - 2*xx - 1);
	b(5,8) = r3*yz*(xx + yy - zz);
	b(6,5) = b(5,6);
	b(7,5) = b(5,7);
	b(8,5) = b(5,8);
	b(6,6) = -4*zzxx + xx + zz;
	b(6,7) = zx*(2*yy - 2*xx + 1);
	b(6,8) = r3*zx*(xx + yy - zz);
	b(7,6) = b(6,7);
	b(8,6) = b(6,8);
	b(7,7) = xx + yy - (xx - yy)*(xx - yy);
	b(7,8) = -r3*zz*(xx - yy);
	b(8,7) = b(7,8);
	b(8,8) = 3*zz*(xx + yy);
	return b;
}

Matrix<dcomp, 9, 9> dddd(const Vector3d &pos){
	double r3 = sqrt(3.);
	Matrix<dcomp, 9, 9> b;
	b.fill(0.);
	double x, y, z;
	Vector3d X, Y, Z;
	X << 1, 0, 0;
	Y << 0, 1, 0;
	Z << 0, 0, 1;
	x = pos.dot(X)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	y = pos.dot(Y)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	z = pos.dot(Z)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
        double xx=x*x;
        double xy=x*y;
        double yy=y*y;
        double yz=y*z;
        double zz=z*z;
        double zx=z*x;
        double xxyy=xx*yy;
        double yyzz=yy*zz;
        double zzxx=zz*xx;
	b(4,4) = xxyy + zz;
	b(4,5) = zx*(yy - 1);
	b(4,6) = yz*(xx - 1);
	b(4,7) = 0.5*xy*(xx - yy);
	b(4,8) = 0.5*r3*xy*(zz + 1);
	b(5,4) = b(4,5);
	b(6,4) = b(4,6);
	b(7,4) = b(4,7);
	b(8,4) = b(4,8);
	b(5,5) = yyzz + xx;
	b(5,6) = xy*(zz - 1);
	b(5,7) = 0.5*yz*(2 + xx - yy);
	b(5,8) = -0.5*yz*r3*(xx + yy);
	b(6,5) = b(5,6);
	b(7,5) = b(5,7);
	b(8,5) = b(5,8);
	b(6,6) = zzxx + yy;
	b(6,7) = 0.5*zx*(xx - yy - 2);
	b(6,8) = -0.5*zx*r3*(xx + yy);
	b(7,6) = b(6,7);
	b(8,6) = b(6,8);
	b(7,7) = zz + 0.25*(xx - yy)*(xx - yy);
	b(7,8) = 0.25*r3*(zz + 1)*(xx - yy);
	b(8,7) = b(7,8);
	b(8,8) = 0.75*(xx + yy)*(xx + yy);
	return b;
}

#endif
