#ifndef CUNNINGHAM_H
#define CUNNINGHAM_H

#include <cmath>
#include <eigen3/Eigen/Dense>
/* #include <string> */
/* #include <fstream> */
#include "calc_spawn.h"
#include <omp.h>

using namespace std;
typedef complex<double> dcomp;
typedef vector<VectorXcd, aligned_allocator<VectorXcd>> vVXcd;

template <typename... Args>// a square can only spawn a square
void aux_square(int depth, double error, int B, double quad, double dx, double dy, vVXcd &fu, vVXcd &fd,
	       vVXcd &fud, vVXcd &fdu, VectorXcd &zresu, VectorXcd &zresd, VectorXcd &zresud,
	       VectorXcd &zresdu, Vector3d &d1, Vector3d &d2, int nsub, int nsubat,
	       int ndiff, Vector3d &b1, Vector3d &b2, const vector<pair<int,int>> &ifold, int nfold, 
	       const Matrix3d &baib, double fact, int irecip, dcomp zener, Args&&... params)
{

	/* string Mydata; */
	/* ofstream Myfile; */	
	/* Mydata = "output2.txt"; */
	/* Myfile.open( Mydata.c_str(),ios::app ); */

	int nxfold, ifail;
        VectorXcd zconu(ndiff), zcond(ndiff), zconud(ndiff), zcondu(ndiff);
	Vector3d xk;
        zconu.fill(0);
        zcond.fill(0);
        zconud.fill(0);
        zcondu.fill(0);
        zresu.fill(0);
        zresd.fill(0);
        zresud.fill(0);
        zresdu.fill(0);
        vector<Vector3d, aligned_allocator<Vector3d>> xfold;
	vVXcd f1u, f1d, f1du, f1ud;
	vVXcd f2u, f2d, f2du, f2ud;
	vVXcd f3u, f3d, f3du, f3ud;
	vVXcd f4u, f4d, f4du, f4ud;
	double x, z;

//       calculate the function

	double nk = quad*quad;
        int A = 2*B;
	for (int i = 1; i<=A; i++){
	  for (int j = 1; j<=A; j++){
	    if ((i%2 != 0) && (j%2 != 0)){
	      x = i/(A*quad)+dx;
	      z = j/(A*quad)+dy;
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

	      if (i < B){
		if (j < B){
		  zconu = zconu*(2/(nk*B*B));
                  f1u.emplace_back(zconu);
		  zcond = zcond*(2/(nk*B*B));
                  f1d.emplace_back(zcond);
		  zconud = zconud*(2/(nk*B*B));
                  f1ud.emplace_back(zconud);
		  zcondu = zcondu*(2/(nk*B*B));
                  f1du.emplace_back(zcondu);
		}
		if (j > B){
		  zconu = zconu*(2/(nk*B*B));
                  f2u.emplace_back(zconu);
		  zcond = zcond*(2/(nk*B*B));
                  f2d.emplace_back(zcond);
		  zconud = zconud*(2/(nk*B*B));
                  f2ud.emplace_back(zconud);
		  zcondu = zcondu*(2/(nk*B*B));
                  f2du.emplace_back(zcondu);
		}
	      }

	      if (i > B){
	        if (j < B){
		  zconu = zconu*(2/(nk*B*B));
                  f3u.emplace_back(zconu);
		  zcond = zcond*(2/(nk*B*B));
                  f3d.emplace_back(zcond);
		  zconud = zconud*(2/(nk*B*B));
                  f3ud.emplace_back(zconud);
		  zcondu = zcondu*(2/(nk*B*B));
                  f3du.emplace_back(zcondu);
		}

		if (j > B){
		  zconu = zconu*(2/(nk*B*B));
                  f4u.emplace_back(zconu);
		  zcond = zcond*(2/(nk*B*B));
                  f4d.emplace_back(zcond);
		  zconud = zconud*(2/(nk*B*B));
                  f4ud.emplace_back(zconud);
		  zcondu = zcondu*(2/(nk*B*B));
                  f4du.emplace_back(zcondu);
		}
	      }
	    }
	  }
	}

	for (int it = 0; it < fu.size(); it++){
	  zresu = zresu + fu[it];
	  zresd = zresd + fd[it];
	  zresud = zresud + fud[it];
	  zresdu = zresdu + fdu[it];
	}
	VectorXd xcold, xc;
        xcold=fact*(zresu+zresd-zresud-zresdu).real()/nsub;
        zresu.fill(0);
        zresd.fill(0);
        zresud.fill(0);
        zresdu.fill(0);
	for (int it = 0; it < f1u.size(); it++){
	  zresu = zresu + f1u[it] + f2u[it] + f3u[it] + f4u[it]; 
	  zresd = zresd + f1d[it] + f2d[it] + f3d[it] + f4d[it]; 
	  zresud = zresud + f1ud[it] + f2ud[it] + f3ud[it] + f4ud[it]; 
	  zresdu = zresdu + f1du[it] + f2du[it] + f3du[it] + f4du[it]; 
	}
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

	quad*=2.;
	error /= 2.;
	aux_square(depth-1, error, B, quad, dx, dy, f1u, f1d, f1ud, f1du, g1u, g1d, g1ud, g1du,
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth-1, error, B, quad, dx + 1/quad, dy, f2u, f2d, f2ud, f2du, g2u, g2d, g2ud, g2du,
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth-1, error, B, quad, dx, dy + 1/quad, f3u, f3d, f3ud, f3du, g3u, g3d, g3ud, g3du,
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth-1, error, B, quad, dx + 1/quad, dy + 1/quad, f4u, f4d, f4ud, f4du, g4u, g4d, g4ud, g4du,
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);

        zresu.fill(0);
        zresd.fill(0);
        zresud.fill(0);
        zresdu.fill(0);

	zresu = g1u + g2u + g3u + g4u;
	zresd = g1d + g2d + g3d + g4d;
	zresud = g1ud + g2ud + g3ud + g4ud;
	zresdu = g1du + g2du + g3du + g4du;

	return;
	}
}

template <typename... Args>// a triangle can spawn square and triangles, but only triangles when x = z
void aux_tri(int depth, double error, int B, double quad, double dx, double dy, vVXcd &fu, vVXcd &fd,
	       vVXcd &fud, vVXcd &fdu, VectorXcd &zresu, VectorXcd &zresd, VectorXcd &zresud,
	       VectorXcd &zresdu, Vector3d &d1, Vector3d &d2, int nsub, int nsubat,
	       int ndiff, Vector3d &b1, Vector3d &b2, const vector<pair<int,int>> &ifold, int nfold, 
	       const Matrix3d &baib, double fact, int irecip, dcomp zener, Args&&... params)
{

	/* string Mydata; */
	/* ofstream Myfile; */	
	/* Mydata = "output2.txt"; */
	/* Myfile.open( Mydata.c_str(),ios::app ); */

	int nxfold, ifail;
        VectorXcd zconu(ndiff), zcond(ndiff), zconud(ndiff), zcondu(ndiff);
	Vector3d xk;
        zconu.fill(0);
        zcond.fill(0);
        zconud.fill(0);
        zcondu.fill(0);
        zresu.fill(0);
        zresd.fill(0);
        zresud.fill(0);
        zresdu.fill(0);
        vector<Vector3d, aligned_allocator<Vector3d>> xfold;
	vVXcd f1u, f1d, f1du, f1ud;
	vVXcd f2u, f2d, f2du, f2ud;
	vVXcd f3u, f3d, f3du, f3ud;
	double x, z;

//       calculate the function

	double nk = quad*quad;
        int A = 2*B;
	for (int i = 1; i<=A; i++){
	  for (int j = 1; j<=i; j++){
	    if ((i%2 != 0) && (j%2 != 0)){
	      x = i/(A*quad)+dx;
	      z = j/(A*quad)+dy;
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

	      if (i < B){
		if (i == j){
		  zconu = zconu*(1/(nk*B*B));
                  f1u.emplace_back(zconu);
		  zcond = zcond*(1/(nk*B*B));
                  f1d.emplace_back(zcond);
		  zconud = zconud*(1/(nk*B*B));
                  f1ud.emplace_back(zconud);
		  zcondu = zcondu*(1/(nk*B*B));
                  f1du.emplace_back(zcondu);
		}
		else{
		  zconu = zconu*(2/(nk*B*B));
                  f1u.emplace_back(zconu);
		  zcond = zcond*(2/(nk*B*B));
                  f1d.emplace_back(zcond);
		  zconud = zconud*(2/(nk*B*B));
                  f1ud.emplace_back(zconud);
		  zcondu = zcondu*(2/(nk*B*B));
                  f1du.emplace_back(zcondu);
		}
	      }

	      if (i > B){
	        if (j < B){
		  zconu = zconu*(2/(nk*B*B));
                  f2u.emplace_back(zconu);
		  zcond = zcond*(2/(nk*B*B));
                  f2d.emplace_back(zcond);
		  zconud = zconud*(2/(nk*B*B));
                  f2ud.emplace_back(zconud);
		  zcondu = zcondu*(2/(nk*B*B));
                  f2du.emplace_back(zcondu);
		}

		if (j > B){
		  if (i == j){
		    zconu = zconu*(1/(nk*B*B));
                    f3u.emplace_back(zconu);
		    zcond = zcond*(1/(nk*B*B));
                    f3d.emplace_back(zcond);
		    zconud = zconud*(1/(nk*B*B));
                    f3ud.emplace_back(zconud);
		    zcondu = zcondu*(1/(nk*B*B));
                    f3du.emplace_back(zcondu);
	  	  }
		  else{
		    zconu = zconu*(2/(nk*B*B));
                    f3u.emplace_back(zconu);
		    zcond = zcond*(2/(nk*B*B));
                    f3d.emplace_back(zcond);
		    zconud = zconud*(2/(nk*B*B));
                    f3ud.emplace_back(zconud);
		    zcondu = zcondu*(2/(nk*B*B));
                    f3du.emplace_back(zcondu);
		  }
		}
	      }
	    }
	  }
	}

	for (int it = 0; it < fu.size(); it++){
	  zresu = zresu + fu[it];
	  zresd = zresd + fd[it];
	  zresud = zresud + fud[it];
	  zresdu = zresdu + fdu[it];
	}
	VectorXd xcold, xc;
        xcold=fact*(zresu+zresd-zresud-zresdu).real()/nsub;
        zresu.fill(0);
        zresd.fill(0);
        zresud.fill(0);
        zresdu.fill(0);
	for (int it = 0; it < f1u.size(); it++){
	  zresu = zresu + f1u[it] + f3u[it]; 
	  zresd = zresd + f1d[it] + f3d[it]; 
	  zresud = zresud + f1ud[it] + f3ud[it]; 
	  zresdu = zresdu + f1du[it] + f3du[it]; 
	}
	for (int it = 0; it < f2u.size(); it++){
	  zresu = zresu + f2u[it];
	  zresd = zresd + f2d[it];
	  zresud = zresud + f2ud[it];
	  zresdu = zresdu + f2du[it];
	}
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

	quad*=2.;
	error /= 2.;
	aux_tri(depth-1, error, B, quad, dx, dy, f1u, f1d, f1ud, f1du, g1u, g1d, g1ud, g1du,
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_square(depth-1, error, B, quad, dx + 1/quad, dy, f2u, f2d, f2ud, f2du, g2u, g2d, g2ud, g2du,
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	aux_tri(depth-1, error, B, quad, dx + 1/quad, dy + 1/quad, f3u, f3d, f3ud, f3du, g3u, g3d, g3ud, g3du,
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);

        zresu.fill(0);
        zresd.fill(0);
        zresud.fill(0);
        zresdu.fill(0);
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
	int depth = 3;
        VectorXcd zconu(ndiff), zcond(ndiff), zconud(ndiff), zcondu(ndiff);
        zconu.fill(0);
        zcond.fill(0);
        zconud.fill(0);
        zcondu.fill(0);
        zresu.fill(0);
        zresd.fill(0);
        zresud.fill(0);
        zresdu.fill(0);
        vector<Vector3d, aligned_allocator<Vector3d>> xfold;
        VectorXd xc(ndiff), xcold(ndiff);
	vVXcd f1u, f1d, f1du, f1ud;
	vVXcd f2u, f2d, f2du, f2ud;
	vVXcd f3u, f3d, f3du, f3ud;
        Vector3d d1, d2, xk;
        d1=b1/2.;
        d2=b2/2.;
	double x, z;
	int tid;
        int maxt = omp_get_max_threads();

//       calculate the function

	double nk = 1.;
        int B = 4;//number of points per axis
        int A = 2*B;
	int size = B/2;
	int C = 0;
	for (int c = 1; c <= size; c++)
	  C = C + c;
	for (int it = 0; it < C; it++){
	  f1u.emplace_back(zresu);
	  f1d.emplace_back(zresu);
	  f1ud.emplace_back(zresu);
	  f1du.emplace_back(zresu);
	  f3u.emplace_back(zresu);
	  f3d.emplace_back(zresu);
	  f3ud.emplace_back(zresu);
	  f3du.emplace_back(zresu);
	}
	for (int it = 0; it < size*size; it++){
	  f2u.emplace_back(zresu);
	  f2d.emplace_back(zresu);
	  f2ud.emplace_back(zresu);
	  f2du.emplace_back(zresu);
	}
	int I;
	for (int i = 1; i<=A; i++){
	  if (i%2 != 0){
	    x = i/(A*1.);
	    I = (i-1)/2;
#pragma omp parallel private(tid, z, xk, nxfold, xfold, ifail)
		{
            tid = omp_get_thread_num();
	    for (int j = 1; j<=i; j++){
	      if (j%2 != 0){
		int J = (j-1)/2;
		int ind;
                if (J%maxt == tid){
		  int index;
		  ind = 0;
                  VectorXcd zconu(ndiff), zcond(ndiff), zconud(ndiff), zcondu(ndiff);
	          z = j/(A*1.);
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

	          if (i < B){
		    for (int c = 1; c<=I; c++)
	              ind += c;
	            index = J + ind;
		    if (i == j){
		      zconu = zconu*(1/(nk*B*B));
                      f1u[index] = zconu;
    		      zcond = zcond*(1/(nk*B*B));
                      f1d[index] = zcond;
    		      zconud = zconud*(1/(nk*B*B));
                      f1ud[index] = zconud;
    		      zcondu = zcondu*(1/(nk*B*B));
                      f1du[index] = zcondu;
		    }
		    else{
		      zconu = zconu*(2/(nk*B*B));
                      f1u[index] = zconu;
		      zcond = zcond*(2/(nk*B*B));
                      f1d[index] = zcond;
		      zconud = zconud*(2/(nk*B*B));
                      f1ud[index] = zconud;
		      zcondu = zcondu*(2/(nk*B*B));
                      f1du[index] = zcondu;
	  	    }
	          }

	          if (i > B){
	            if (j < B){
	              index = J + size*(I - size);
		      zconu = zconu*(2/(nk*B*B));
                      f2u[index] = zconu;
		      zcond = zcond*(2/(nk*B*B));
                      f2d[index] = zcond;
		      zconud = zconud*(2/(nk*B*B));
                      f2ud[index] = zconud;
		      zcondu = zcondu*(2/(nk*B*B));
                      f2du[index] = zcondu;
	  	    }
		    if (j > B){
		      ind = 0;
		      for (int c = 1; c<=I-size; c++)
	                ind += c;
	              index = J + ind - size;
		      if (i == j){
		        zconu = zconu*(1/(nk*B*B));
                        f3u[index] = zconu;
		        zcond = zcond*(1/(nk*B*B));
                        f3d[index] = zcond;
		        zconud = zconud*(1/(nk*B*B));
                        f3ud[index] = zconud;
		        zcondu = zcondu*(1/(nk*B*B));
                        f3du[index] = zcondu;
	  	      }
		      else{
		        zconu = zconu*(2/(nk*B*B));
                        f3u[index] = zconu;
		        zcond = zcond*(2/(nk*B*B));
                        f3d[index] = zcond;
		        zconud = zconud*(2/(nk*B*B));
                        f3ud[index] = zconud;
		        zcondu = zcondu*(2/(nk*B*B));
                        f3du[index] = zcondu;
		      }
		    }
	          }
	        }
	      }
	    }
		}
	  }
	}
        double error=5e-6;

	//begin spawning points around points already integrated over
	VectorXcd g1u(ndiff), g1d(ndiff), g1du(ndiff), g1ud(ndiff);
	VectorXcd g2u(ndiff), g2d(ndiff), g2du(ndiff), g2ud(ndiff);
	VectorXcd g3u(ndiff), g3d(ndiff), g3du(ndiff), g3ud(ndiff);

	double quad = 2.;
	error /= 2.;
#pragma omp parallel private(tid)
	{
        tid = omp_get_thread_num();
	if (tid == 0){
	aux_tri(depth, error, B, 2., 0, 0, f1u, f1d, f1ud, f1du, g1u, g1d, g1ud, g1du,
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	}
        if (tid == 1%maxt){ 
	aux_square(depth, error, B, 2., 0.5, 0, f2u, f2d, f2ud, f2du, g2u, g2d, g2ud, g2du,
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	}
	if (tid == 2%maxt){
	aux_tri(depth, error, B, 2., 0.5, 0.5, f3u, f3d, f3ud, f3du, g3u, g3d, g3ud, g3du,
		d1, d2, nsub, nsubat, ndiff, b1, b2, ifold, nfold, baib, fact, irecip, zener, params...);
	}
        }

	zresu = g1u + g2u + g3u; 
	zresd = g1d + g2d + g3d; 
	zresud = g1ud + g2ud + g3ud; 
	zresdu = g1du + g2du + g3du; 

	return;
}
#endif
