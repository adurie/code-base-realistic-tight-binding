#ifndef CUNNINGHAM_H
#define CUNNINGHAM_H

#include <cmath>
#include <eigen3/Eigen/Dense>
#include <string>
#include <fstream>
#include "calc_spawn.h"

using namespace std;
typedef complex<double> dcomp;

template <typename... Args>// a square can only spawn a square
void aux_square(int depth, double error, double n, double v, double w, VectorXcd &fu, VectorXcd &fd,
	       VectorXcd &fud, VectorXcd &fdu, VectorXcd &zresu, VectorXcd &zresd, VectorXcd &zresud,
	       VectorXcd &zresdu, double A, Vector3d &d1, Vector3d &d2, int nsub, int nsubat,
	       int ndiff, Vector3d &b1, Vector3d &b2, const vector<pair<int,int>> &ifold, int nfold, 
	       const Matrix3d &baib, double fact, int irecip, dcomp zener, Args&&... params)
{

	string Mydata;
	ofstream Myfile;	
	Mydata = "output2.txt";
	Myfile.open( Mydata.c_str(),ios::app );

        VectorXcd zconu(ndiff), zcond(ndiff), zconud(ndiff), zcondu(ndiff);
	Vector3d xk;
        vector<Vector3d, aligned_allocator<Vector3d>> xfold;
	VectorXd xcold, xc;
        xcold=fact*(fu + fd - fdu - fud).real()/nsub;
	VectorXcd f1u(ndiff), f1d(ndiff), f1du(ndiff), f1ud(ndiff);
	VectorXcd f2u(ndiff), f2d(ndiff), f2du(ndiff), f2ud(ndiff);
	VectorXcd f3u(ndiff), f3d(ndiff), f3du(ndiff), f3ud(ndiff);
	VectorXcd f4u(ndiff), f4d(ndiff), f4du(ndiff), f4ud(ndiff);
	double x, z;
	x = v - A;		z = w - A;	
        xk=x*d1+z*d2;
	/* Myfile<<x<<" "<<z<<endl; */
        int nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

        int ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f1u = zconu*(2/n);
        f1d = zcond*(2/n);
        f1ud = zconud*(2/n);
        f1du = zcondu*(2/n);

	x = v - A;	z = w + A;	
        xk=x*d1+z*d2;
	/* Myfile<<x<<" "<<z<<endl; */
        nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

        ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f2u = zconu*(2/n);
        f2d = zcond*(2/n);
        f2ud = zconud*(2/n);
        f2du = zcondu*(2/n);

	x = v + A;	z = w - A;	
        xk=x*d1+z*d2;
	/* Myfile<<x<<" "<<z<<endl; */
        nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

        ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f3u = zconu*(2/n);
        f3d = zcond*(2/n);
        f3ud = zconud*(2/n);
        f3du = zcondu*(2/n);

	x = v + A;	z = w + A;	
        xk=x*d1+z*d2;
	/* Myfile<<x<<" "<<z<<endl; */
        nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

        ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f4u = zconu*(2/n);
        f4d = zcond*(2/n);
        f4ud = zconud*(2/n);
        f4du = zcondu*(2/n);

	zresu = f1u + f2u + f3u + f4u;
	zresd = f1d + f2d + f3d + f4d;
	zresud = f1ud + f2ud + f3ud + f4ud;
	zresdu = f1du + f2du + f3du + f4du;
        xc=fact*(zresu+zresd-zresud-zresdu).real()/nsub;
	double rel_error;
	rel_error = conv(xcold, xc, ndiff);
	/* cout<<scientific<<xc.transpose()<<endl<<xcold.transpose()<<endl<<endl; */
	/* cout<<scientific<<rel_error<<endl; */

	if ((depth <= 0) && (rel_error >= error)){
	  /* cout<<scientific<<"Error, maximum depth exceeded, error: "<<rel_error<<endl; */
	  return;
	}
        if (rel_error < error){
	  /* cout<<"converged xcon = "<<error<<"; error = "<<rel_error<<endl; */
	  return;
	}

	if (rel_error >= error){
	VectorXcd g1u(ndiff), g1d(ndiff), g1du(ndiff), g1ud(ndiff);
	VectorXcd g2u(ndiff), g2d(ndiff), g2du(ndiff), g2ud(ndiff);
	VectorXcd g3u(ndiff), g3d(ndiff), g3du(ndiff), g3ud(ndiff);
	VectorXcd g4u(ndiff), g4d(ndiff), g4du(ndiff), g4ud(ndiff);

	aux_square(depth-1, error/n, n*4, v - A, w - A, f1u, f1d, f1ud, f1du, g1u, g1d, g1ud, g1du, A/2., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth-1, error/n, n*4, v - A, w + A, f2u, f2d, f2ud, f2du, g2u, g2d, g2ud, g2du, A/2., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth-1, error/n, n*4, v + A, w - A, f3u, f3d, f3ud, f3du, g3u, g3d, g3ud, g3du, A/2., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth-1, error/n, n*4, v + A, w + A, f4u, f4d, f4ud, f4du, g4u, g4d, g4ud, g4du, A/2., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);

	zresu = g1u + g2u + g3u + g4u;
	zresd = g1d + g2d + g3d + g4d;
	zresud = g1ud + g2ud + g3ud + g4ud;
	zresdu = g1du + g2du + g3du + g4du;

	return;
	}
}

template <typename... Args>// a triangle can spawn square and triangles, but only triangles when x = z
void aux_tri(int depth, double error, double n, double v, double w, VectorXcd &fu, VectorXcd &fd,
	       VectorXcd &fud, VectorXcd &fdu, VectorXcd &zresu, VectorXcd &zresd, VectorXcd &zresud,
	       VectorXcd &zresdu, double A, Vector3d &d1, Vector3d &d2, int nsub, int nsubat,
	       int ndiff, Vector3d &b1, Vector3d &b2, const vector<pair<int,int>> &ifold, int nfold, 
	       const Matrix3d &baib, double fact, int irecip, dcomp zener, Args&&... params)
{
	string Mydata;
	ofstream Myfile;	
	Mydata = "output2.txt";
	Myfile.open( Mydata.c_str(),ios::app );

        VectorXcd zconu(ndiff), zcond(ndiff), zconud(ndiff), zcondu(ndiff);
	Vector3d xk;
        vector<Vector3d, aligned_allocator<Vector3d>> xfold;
	VectorXd xcold, xc;
        xcold=fact*(fu + fd - fdu - fud).real()/nsub;
	VectorXcd f1u(ndiff), f1d(ndiff), f1du(ndiff), f1ud(ndiff);
	VectorXcd f2u(ndiff), f2d(ndiff), f2du(ndiff), f2ud(ndiff);
	VectorXcd f3u(ndiff), f3d(ndiff), f3du(ndiff), f3ud(ndiff);
	double x, z;
	x = v - A;	z = w - A;	
        xk=x*d1+z*d2;
	/* Myfile<<x<<" "<<z<<endl; */
	int nxfold;
        nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

	int ifail;
        ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f1u = zconu*(1/n);
        f1d = zcond*(1/n);
        f1ud = zconud*(1/n);
        f1du = zcondu*(1/n);

	x = v + A;	z = w - A;	
        xk=x*d1+z*d2;
	/* Myfile<<x<<" "<<z<<endl; */
        nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

        ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f2u = zconu*(2/n);
        f2d = zcond*(2/n);
        f2ud = zconud*(2/n);
        f2du = zcondu*(2/n);

	x = v + A;	z = w + A;	
        xk=x*d1+z*d2;
	/* Myfile<<x<<" "<<z<<endl; */
        nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

        ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f3u = zconu*(1/n);
        f3d = zcond*(1/n);
        f3ud = zconud*(1/n);
        f3du = zcondu*(1/n);

	zresu = f1u + f2u + f3u;
	zresd = f1d + f2d + f3d;
	zresud = f1ud + f2ud + f3ud;
	zresdu = f1du + f2du + f3du;
        xc=fact*(zresu+zresd-zresud-zresdu).real()/nsub;
	double rel_error;
	rel_error = conv(xcold, xc, ndiff);
	/* cout<<scientific<<xc.transpose()<<endl<<xcold.transpose()<<endl<<endl; */

	if ((depth <= 0) && (rel_error >= error)){
	  /* cout<<scientific<<"Error, maximum depth exceeded, error: "<<rel_error<<endl; */
	  return;
	}
        if (rel_error < error){
	  /* cout<<"converged xcon = "<<error<<"; error = "<<rel_error<<endl; */
	  return;
	}

	if (rel_error >= error){
	//begin spawning points around points already integrated over
	VectorXcd g1u(ndiff), g1d(ndiff), g1du(ndiff), g1ud(ndiff);
	VectorXcd g2u(ndiff), g2d(ndiff), g2du(ndiff), g2ud(ndiff);
	VectorXcd g3u(ndiff), g3d(ndiff), g3du(ndiff), g3ud(ndiff);

	aux_tri(depth-1, error/n, n*4, v - A, w - A, f1u, f1d, f1ud, f1du, g1u, g1d, g1ud, g1du, A/2., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth-1, error/n, n*4, v + A, w - A, f2u, f2d, f2ud, f2du, g2u, g2d, g2ud, g2du, A/2., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_tri(depth-1, error/n, n*4, v + A, w + A, f3u, f3d, f3ud, f3du, g3u, g3d, g3ud, g3du, A/2., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);

	zresu = g1u + g2u + g3u;
	zresd = g1d + g2d + g3d;
	zresud = g1ud + g2ud + g3ud;
	zresdu = g1du + g2du + g3du;

	return;
	}
}

template <typename... Args>
void kcon(int nsubat, const vector<pair<int,int>> &ifold, int nfold,const Matrix3d &baib, int nsub, int ndiff, double fact, 
	VectorXcd &zresu, VectorXcd &zresd, VectorXcd &zresud, VectorXcd &zresdu, int irecip, Vector3d &b1, Vector3d &b2,
	dcomp zener, Args&&... params){

	/* string Mydata; */
	/* ofstream Myfile; */	
	/* Mydata = "output2.txt"; */
	/* Myfile.open( Mydata.c_str(),ios::app ); */

        int ifail;
        int nxfold;
	int depth = 4;
        VectorXcd zconu(ndiff), zcond(ndiff), zconud(ndiff), zcondu(ndiff);
        zconu.fill(0);
        zcond.fill(0);
        zconud.fill(0);
        zcondu.fill(0);
        vector<Vector3d, aligned_allocator<Vector3d>> xfold;
        VectorXd xc(ndiff), xcold(ndiff);
	VectorXcd f1u(ndiff), f1d(ndiff), f1du(ndiff), f1ud(ndiff);
	VectorXcd f2u(ndiff), f2d(ndiff), f2du(ndiff), f2ud(ndiff);
	VectorXcd f3u(ndiff), f3d(ndiff), f3du(ndiff), f3ud(ndiff);
        Vector3d d1, d2, xk;
        d1=b1/2.;
        d2=b2/2.;
	double x, z;
	double n = 4.;

//       calculate the function

        double A = 1.;
	x = A/2.;
	z = A/2.;
        xk=x*d1+z*d2;
	//calculate the mean value point only	
        nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

        ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f1u = zconu;
        f1d = zcond;
        f1ud = zconud;
        f1du = zcondu;
        xcold=fact*(f1u + f1d - f1du - f1ud).real()/nsub;

	//calculate the 5 remaining points in the 6 point grid
	A /= 4.;
	x = A;		z = A;	
        xk=x*d1+z*d2;
	/* Myfile<<x<<" "<<z<<endl; */
        nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

        ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f1u = zconu*(1/n);
        f1d = zcond*(1/n);
        f1ud = zconud*(1/n);
        f1du = zcondu*(1/n);

	x = 3*A;	z = A;	
        xk=x*d1+z*d2;
	/* Myfile<<x<<" "<<z<<endl; */
        nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

        ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f2u = zconu*(2/n);
        f2d = zcond*(2/n);
        f2ud = zconud*(2/n);
        f2du = zcondu*(2/n);

	x = 3*A;	z = 3*A;	
        xk=x*d1+z*d2;
	/* Myfile<<x<<" "<<z<<endl; */
        nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

        ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f3u = zconu*(1/n);
        f3d = zcond*(1/n);
        f3ud = zconud*(1/n);
        f3du = zcondu*(1/n);

	zresu = f1u + f2u + f3u;
	zresd = f1d + f2d + f3d;
	zresud = f1ud + f2ud + f3ud;
	zresdu = f1du + f2du + f3du;
        xc=fact*(zresu+zresd-zresud-zresdu).real()/nsub;
	/* double rel_error; */
	/* rel_error = conv(xcold, xc, ndiff); */

        double error=5e-6;
        /* if (rel_error < error){ */
	  /* cout<<"converged xcon = "<<error<<"; error = "<<rel_error<<endl; */
	  /* return; */
	/* } */

	/* if (rel_error >= error){ */
	//begin spawning points around points already integrated over
	VectorXcd g1u(ndiff), g1d(ndiff), g1du(ndiff), g1ud(ndiff);
	VectorXcd g2u(ndiff), g2d(ndiff), g2du(ndiff), g2ud(ndiff);
	VectorXcd g3u(ndiff), g3d(ndiff), g3du(ndiff), g3ud(ndiff);

	aux_tri(depth, error/n, 4*n, A, A, f1u, f1d, f1ud, f1du, g1u, g1d, g1ud, g1du, A/2., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth, error/n, 4*n, 3*A, A, f2u, f2d, f2ud, f2du, g2u, g2d, g2ud, g2du, A/2.,
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_tri(depth, error/n, 4*n, 3*A, 3*A, f3u, f3d, f3ud, f3du, g3u, g3d, g3ud, g3du, A/2.,
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);

	zresu = g1u + g2u + g3u;
	zresd = g1d + g2d + g3d;
	zresud = g1ud + g2ud + g3ud;
	zresdu = g1du + g2du + g3du;

	return;
	/* } */
}
#endif
