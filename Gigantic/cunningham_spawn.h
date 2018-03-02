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
        VectorXcd zconu(ndiff), zcond(ndiff), zconud(ndiff), zcondu(ndiff);
	Vector3d xk;
        vector<Vector3d, aligned_allocator<Vector3d>> xfold;
	VectorXd xcold, xc;
        xcold=fact*(fu + fd - fdu - fud).real()/nsub;
	VectorXcd f1u(ndiff), f1d(ndiff), f1du(ndiff), f1ud(ndiff);
	VectorXcd f2u(ndiff), f2d(ndiff), f2du(ndiff), f2ud(ndiff);
	VectorXcd f3u(ndiff), f3d(ndiff), f3du(ndiff), f3ud(ndiff);
	VectorXcd f4u(ndiff), f4d(ndiff), f4du(ndiff), f4ud(ndiff);
	VectorXcd f5u(ndiff), f5d(ndiff), f5du(ndiff), f5ud(ndiff);
	VectorXcd f6u(ndiff), f6d(ndiff), f6du(ndiff), f6ud(ndiff);
	VectorXcd f7u(ndiff), f7d(ndiff), f7du(ndiff), f7ud(ndiff);
	VectorXcd f8u(ndiff), f8d(ndiff), f8du(ndiff), f8ud(ndiff);
	VectorXcd f9u(ndiff), f9d(ndiff), f9du(ndiff), f9ud(ndiff);
	double x, z;
	x = v;		z = w;	
	f1u = (1/9.)*fu;
	f1d = (1/9.)*fd;
	f1ud = (1/9.)*fud;
	f1du = (1/9.)*fdu;

	x = v;	z = w + 2*A;	
        xk=x*d1+z*d2;
        int nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

        int ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f2u = zconu*(2/n);
        f2d = zcond*(2/n);
        f2ud = zconud*(2/n);
        f2du = zcondu*(2/n);

	x = v;	z = w - 2*A;	
        xk=x*d1+z*d2;
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

	x = v + 2*A;	z = w;	
        xk=x*d1+z*d2;
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

	x = v + 2*A;	z = w + 2*A;	
        xk=x*d1+z*d2;
        nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

        ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f5u = zconu*(2/n);
        f5d = zcond*(2/n);
        f5ud = zconud*(2/n);
        f5du = zcondu*(2/n);

	x = v + 2*A;	z = w - 2*A;	
        xk=x*d1+z*d2;
        nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

        ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f6u = zconu*(2/n);
        f6d = zcond*(2/n);
        f6ud = zconud*(2/n);
        f6du = zcondu*(2/n);

	x = v - 2*A;	z = w;	
        xk=x*d1+z*d2;
        nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

        ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f7u = zconu*(2/n);
        f7d = zcond*(2/n);
        f7ud = zconud*(2/n);
        f7du = zcondu*(2/n);

	x = v - 2*A;	z = w + 2*A;	
        xk=x*d1+z*d2;
        nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

        ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f8u = zconu*(2/n);
        f8d = zcond*(2/n);
        f8ud = zconud*(2/n);
        f8du = zcondu*(2/n);

	x = v - 2*A;	z = w - 2*A;	
        xk=x*d1+z*d2;
        nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

        ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f9u = zconu*(2/n);
        f9d = zcond*(2/n);
        f9ud = zconud*(2/n);
        f9du = zcondu*(2/n);

	zresu = f1u + f2u + f3u + f4u + f5u + f6u + f7u + f8u + f9u;
	zresd = f1d + f2d + f3d + f4d + f5d + f6d + f7d + f8d + f9d;
	zresud = f1ud + f2ud + f3ud + f4ud + f5ud + f6ud + f7ud + f8ud + f9ud;
	zresdu = f1du + f2du + f3du + f4du + f5du + f6du + f7du + f8du + f9du;
        xc=fact*(zresu+zresd-zresud-zresdu).real()/nsub;
	double rel_error;
	rel_error = conv(xcold, xc, ndiff);
	/* cout<<scientific<<xc.transpose()<<endl<<xcold.transpose()<<endl<<endl; */
	/* cout<<scientific<<rel_error<<endl; */

	if ((depth <= 0) && (rel_error >= error)){
	  cout<<scientific<<"Error, maximum depth exceeded, error: "<<rel_error<<endl;
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
	VectorXcd g5u(ndiff), g5d(ndiff), g5du(ndiff), g5ud(ndiff);
	VectorXcd g6u(ndiff), g6d(ndiff), g6du(ndiff), g6ud(ndiff);
	VectorXcd g7u(ndiff), g7d(ndiff), g7du(ndiff), g7ud(ndiff);
	VectorXcd g8u(ndiff), g8d(ndiff), g8du(ndiff), g8ud(ndiff);
	VectorXcd g9u(ndiff), g9d(ndiff), g9du(ndiff), g9ud(ndiff);

	aux_square(depth-1, error/(0.8*depth), n*9, v, w, f1u, f1d, f1ud, f1du, g1u, g1d, g1ud, g1du, A/3., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth-1, error/(0.8*depth), n*9, v, w + 2*A, f2u, f2d, f2ud, f2du, g2u, g2d, g2ud, g2du, A/3., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth-1, error/(0.8*depth), n*9, v, w - 2*A, f3u, f3d, f3ud, f3du, g3u, g3d, g3ud, g3du, A/3., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth-1, error/(0.8*depth), n*9, v + 2*A, w, f4u, f4d, f4ud, f4du, g4u, g4d, g4ud, g4du, A/3., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth-1, error/(0.8*depth), n*9, v + 2*A, w + 2*A, f5u, f5d, f5ud, f5du, g5u, g5d, g5ud, g5du, A/3., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth-1, error/(0.8*depth), n*9, v + 2*A, w - 2*A, f6u, f6d, f6ud, f6du, g6u, g6d, g6ud, g6du, A/3., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth-1, error/(0.8*depth), n*9, v - 2*A, w, f7u, f7d, f7ud, f7du, g7u, g7d, g7ud, g7du, A/3., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth-1, error/(0.8*depth), n*9, v - 2*A, w + 2*A, f8u, f8d, f8ud, f8du, g8u, g8d, g8ud, g8du, A/3., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth-1, error/(0.8*depth), n*9, v - 2*A, w - 2*A, f9u, f9d, f9ud, f9du, g9u, g9d, g9ud, g9du, A/3., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);

	zresu = g1u + g2u + g3u + g4u + g5u + g6u + g7u + g8u + g9u;
	zresd = g1d + g2d + g3d + g4d + g5d + g6d + g7d + g8d + g9d;
	zresud = g1ud + g2ud + g3ud + g4ud + g5ud + g6ud + g7ud + g8ud + g9ud;
	zresdu = g1du + g2du + g3du + g4du + g5du + g6du + g7du + g8du + g9du;

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
        VectorXcd zconu(ndiff), zcond(ndiff), zconud(ndiff), zcondu(ndiff);
	Vector3d xk;
        vector<Vector3d, aligned_allocator<Vector3d>> xfold;
	VectorXd xcold, xc;
        xcold=fact*(fu + fd - fdu - fud).real()/nsub;
	VectorXcd f1u(ndiff), f1d(ndiff), f1du(ndiff), f1ud(ndiff);
	VectorXcd f2u(ndiff), f2d(ndiff), f2du(ndiff), f2ud(ndiff);
	VectorXcd f3u(ndiff), f3d(ndiff), f3du(ndiff), f3ud(ndiff);
	VectorXcd f4u(ndiff), f4d(ndiff), f4du(ndiff), f4ud(ndiff);
	VectorXcd f5u(ndiff), f5d(ndiff), f5du(ndiff), f5ud(ndiff);
	VectorXcd f6u(ndiff), f6d(ndiff), f6du(ndiff), f6ud(ndiff);
	double x, z;
	x = v;		z = w;	
	f1u = (1/9.)*fu;
	f1d = (1/9.)*fd;
	f1ud = (1/9.)*fud;
	f1du = (1/9.)*fdu;

	x = v + 2*A;	z = w;	
        xk=x*d1+z*d2;
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

        f2u = zconu*(2/n);
        f2d = zcond*(2/n);
        f2ud = zconud*(2/n);
        f2du = zcondu*(2/n);

	x = v + 2*A;	z = w + 2*A;	
        xk=x*d1+z*d2;
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

	x = v;	z = w - 2*A;	
        xk=x*d1+z*d2;
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

	x = v + 2*A;	z = w - 2*A;	
        xk=x*d1+z*d2;
        nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

        ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f5u = zconu*(2/n);
        f5d = zcond*(2/n);
        f5ud = zconud*(2/n);
        f5du = zcondu*(2/n);

	x = v - 2*A;	z = w - 2*A;	
        xk=x*d1+z*d2;
        nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

        ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f6u = zconu*(1/n);
        f6d = zcond*(1/n);
        f6ud = zconud*(1/n);
        f6du = zcondu*(1/n);

	zresu = f1u + f2u + f3u + f4u + f5u + f6u;
	zresd = f1d + f2d + f3d + f4d + f5d + f6d;
	zresud = f1ud + f2ud + f3ud + f4ud + f5ud + f6ud;
	zresdu = f1du + f2du + f3du + f4du + f5du + f6du;
        xc=fact*(zresu+zresd-zresud-zresdu).real()/nsub;
	double rel_error;
	rel_error = conv(xcold, xc, ndiff);
	/* cout<<scientific<<xc.transpose()<<endl<<xcold.transpose()<<endl<<endl; */

	if ((depth <= 0) && (rel_error >= error)){
	  cout<<scientific<<"Error, maximum depth exceeded, error: "<<rel_error<<endl;
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
	VectorXcd g4u(ndiff), g4d(ndiff), g4du(ndiff), g4ud(ndiff);
	VectorXcd g5u(ndiff), g5d(ndiff), g5du(ndiff), g5ud(ndiff);
	VectorXcd g6u(ndiff), g6d(ndiff), g6du(ndiff), g6ud(ndiff);

	aux_tri(depth-1, error/(0.8*depth), n*9, v, w, f1u, f1d, f1ud, f1du, g1u, g1d, g1ud, g1du, A/3., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth-1, error/(0.8*depth), n*9, v + 2*A, w, f2u, f2d, f2ud, f2du, g2u, g2d, g2ud, g2du, A/3., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_tri(depth-1, error/(0.8*depth), n*9, v + 2*A, w + 2*A, f3u, f3d, f3ud, f3du, g3u, g3d, g3ud, g3du, A/3., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth-1, error/(0.8*depth), n*9, v, w - 2*A, f4u, f4d, f4ud, f4du, g4u, g4d, g4ud, g4du, A/3., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth-1, error/(0.8*depth), n*9, v + 2*A, w - 2*A, f5u, f5d, f5ud, f5du, g5u, g5d, g5ud, g5du, A/3., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_tri(depth-1, error/(0.8*depth), n*9, v - 2*A, w - 2*A, f6u, f6d, f6ud, f6du, g6u, g6d, g6ud, g6du, A/3., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);

	zresu = g1u + g6u + g3u + g4u + g5u + g2u;
	zresd = g1d + g6d + g3d + g4d + g5d + g2d;
	zresud = g1ud + g6ud + g3ud + g4ud + g5ud + g2ud;
	zresdu = g1du + g6du + g3du + g4du + g5du + g2du;

	return;
	}
}

template <typename... Args>
void kcon(int nsubat, const vector<pair<int,int>> &ifold, int nfold,const Matrix3d &baib, int nsub, int ndiff, double fact, 
	VectorXcd &zresu, VectorXcd &zresd, VectorXcd &zresud, VectorXcd &zresdu, int irecip, Vector3d &b1, Vector3d &b2,
	dcomp zener, Args&&... params){

        int ifail;
        int nxfold;
	int depth = 3;
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
	VectorXcd f4u(ndiff), f4d(ndiff), f4du(ndiff), f4ud(ndiff);
	VectorXcd f5u(ndiff), f5d(ndiff), f5du(ndiff), f5ud(ndiff);
	VectorXcd f6u(ndiff), f6d(ndiff), f6du(ndiff), f6ud(ndiff);
        Vector3d d1, d2, xk;
        d1=b1/2.;
        d2=b2/2.;
	double x, z;
	double n = 9.;

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

        f5u = zconu;
        f5d = zcond;
        f5ud = zconud;
        f5du = zcondu;
        xcold=fact*(f5u + f5d - f5du - f5ud).real()/nsub;

	//calculate the 5 remaining points in the 6 point grid
	A /= 6.;
	f5u = (1/n)*f5u;
	f5d = (1/n)*f5d;
	f5ud = (1/n)*f5ud;
	f5du = (1/n)*f5du;

	x = A;		z = A;	
        xk=x*d1+z*d2;
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

	x = 5*A;	z = A;	
        xk=x*d1+z*d2;
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

	x = 5*A;	z = 3*A;	
        xk=x*d1+z*d2;
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

	x = 5*A;	z = 5*A;	
        xk=x*d1+z*d2;
        nxfold = testk(x,z,b1,b2,ifold,xfold,nfold,baib,irecip);
        if (nxfold != nsub/nsubat){
          cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	  exit(EXIT_FAILURE);
	}

        ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
        if (ifail != 0)
          cout<<x<<" "<<z<<endl;

        f6u = zconu*(1/n);
        f6d = zcond*(1/n);
        f6ud = zconud*(1/n);
        f6du = zcondu*(1/n);

	zresu = f1u + f2u + f3u + f4u + f5u + f6u;
	zresd = f1d + f2d + f3d + f4d + f5d + f6d;
	zresud = f1ud + f2ud + f3ud + f4ud + f5ud + f6ud;
	zresdu = f1du + f2du + f3du + f4du + f5du + f6du;
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
	VectorXcd g4u(ndiff), g4d(ndiff), g4du(ndiff), g4ud(ndiff);
	VectorXcd g5u(ndiff), g5d(ndiff), g5du(ndiff), g5ud(ndiff);
	VectorXcd g6u(ndiff), g6d(ndiff), g6du(ndiff), g6ud(ndiff);

	aux_tri(depth, error/(0.8*depth), 9*n, A, A, f1u, f1d, f1ud, f1du, g1u, g1d, g1ud, g1du, A/3., 
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth, error/(0.8*depth), 9*n, 3*A, A, f2u, f2d, f2ud, f2du, g2u, g2d, g2ud, g2du, A/3.,
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth, error/(0.8*depth), 9*n, 5*A, A, f3u, f3d, f3ud, f3du, g3u, g3d, g3ud, g3du, A/3.,
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth, error/(0.8*depth), 9*n, 5*A, 3*A, f4u, f4d, f4ud, f4du, g4u, g4d, g4ud, g4du, A/3.,
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_tri(depth, error/(0.8*depth), 9*n, 3*A, 3*A, f5u, f5d, f5ud, f5du, g5u, g5d, g5ud, g5du, A/3.,
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_tri(depth, error/(0.8*depth), 9*n, 5*A, 5*A, f6u, f6d, f6ud, f6du, g6u, g6d, g6ud, g6du, A/3.,
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);

	zresu = g1u + g5u + g6u + g2u + g3u + g4u;
	zresd = g1d + g5d + g6d + g2d + g3d + g4d;
	zresud = g1ud + g5ud + g6ud + g2ud + g3ud + g4ud;
	zresdu = g1du + g5du + g6du + g2du + g3du + g4du;

	return;
	/* } */
}
#endif
