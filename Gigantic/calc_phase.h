#ifndef CALC_H 
#define CALC_H 

#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>
#include "greens_phase.h"
/* #include "greens_large_array.h" */

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
typedef vector<Vector3d, aligned_allocator<Vector3d>> vV3d;

int testk(double x, double y, const Vector3d &b1, const Vector3d &b2, const vector<pair<int,int>> &ifold,
	       	vV3d &xfold, int nfold, const Matrix3d &baib, int irecip){
      VectorXd xx1(nfold), xx2(nfold);
      xx1.fill(0);
      xx2.fill(0);
      int icnt = 0;
      int nxfold;
      int i, j;
      double xk, yk, x1, x2, dx, dy, ddx, ddy;
      int condition = 0;
      Vector3d xtmp;
      xfold.clear();
//     -----------------------------------------------------------------
      if (irecip >= 3){
        cout<<"reciprocal superlattice not cubic or rectangular"<<endl;
        cout<<"this case has not been coded"<<endl;
	exit(EXIT_FAILURE);
      }
      if (irecip == 0){
        cout<<"TESTK : irecip = 0 ---- NOT CODED"<<endl;
	exit(EXIT_FAILURE);
      }
      if ((irecip == 1) || (irecip == 2)){
//     -----------------------------------------------------------------
//     CUBIC AND RECTANGULAR ATOMIC RECIPROCAL LATTICE
//     first construct supercell reciprocal lattice cluster
//     for the vector J = i*b1 + j*b2
//     then           J = x1*ba1 + x2*ba2
        for (int iff=0; iff<nfold; iff++){
          i=ifold[iff].first;
          j=ifold[iff].second;
          xk=0.5*x+i;
          yk=0.5*y+j;
          x1=baib(0,0)*xk+baib(0,1)*yk;
          x2=baib(1,0)*xk+baib(1,1)*yk;
          if ((abs(x1) <= 0.50000001) && (abs(x2) <= 0.50000001)){
//       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//       first check if this point is not related to any of the others by
//       a reciprocal lattice vector
            condition = 0;
            for (int iii=0; iii<icnt; iii++){
              dx=x1-xx1(iii);
              dy=x2-xx2(iii);
              ddx=floor(abs(dx)+1e-10)-abs(dx);
              ddy=floor(abs(dy)+1e-10)-abs(dy);
              if ((abs(ddx) < 1e-10) && (abs(ddy) < 1e-10))
		condition = 1;
	    }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	    if (condition == 0){
              xx1(icnt)=x1;
              xx2(icnt)=x2;
              icnt++;
              xtmp=(0.5*x+i)*b1+(0.5*y+j)*b2;
      	      xfold.emplace_back(xtmp);
	    }
	  }
	}
        nxfold=icnt;
      }
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     HEXAGONAL ATOMIC RECIPROCAL LATTICE
      else if (irecip == 3){
        cout<<"HEXAGONAL ATOMIC RECIPROCAL LATTICE NOT CODED"<<endl;
	exit(EXIT_FAILURE);
      }
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//     CENTRED-RECTANGULAR ATOMIC RECIPROCAL LATTICE
      else if (irecip == 4){
        cout<<"CENTRED-RECTANGULAR ATOMIC RECIPROCAL LATTICE NOT CODED"<<endl;
	exit(EXIT_FAILURE);
      }
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return nxfold;
}

template <typename... Args>
void sumk(int nsub, int nsubat, const vector<pair<int,int>> &ifold, int nfold, const Matrix3d &baib, int ndiff,
	int nk, dcomp &zresu, dcomp &zresd, dcomp &zresud, dcomp &zresdu,
	int irecip, Vector3d &b1, Vector3d &b2, dcomp zener, Args&&... params){

      double dq;
      if (nk != 1)
        dq=1./(nk-1.);
      if (nk == 1)
	dq=1.;
      int max=nk-1;
      int nktot=0;
      Vector3d d1, d2, xk;
      double x, y;
      int iwght;
      int ifail;
      int nxfold;
      dcomp zconu, zcond, zconud, zcondu;
      zconu = 0.;
      zcond = 0.;
      zconud = 0.;
      zcondu = 0.;
      vector<Vector3d, aligned_allocator<Vector3d>> xfold;
//     -----------------------------------------------------------------
      if (irecip == 0){
        cout<<"SUMK : irecip = 0 ---- NOT CODED"<<endl;
	exit(EXIT_FAILURE);
      }
//     -----------------------------------------------------------------
//     CUBIC
      if (irecip == 1){
        d1=b1/2.;
        d2=b2/2.;

        for (int i=0; i<=max; i++){
          x=dq*i;
          for (int j=0; j<=i; j++){
            y=dq*j;
  
            xk=x*d1+y*d2;
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//         determine weights
            iwght=8;
            if (i == j)
	      iwght=4;
            if (j == 0)
	      iwght=4;
            if (i == max)
	      iwght=4;
            if (i == 0)
	      iwght=1;
            if (j == max)
	      iwght=1;
            if ((i == max) && (j == 0))
	      iwght=2;
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//         find the folding points

            nxfold = testk(x,y,b1,b2,ifold,xfold,nfold,baib,irecip);
            if (nxfold != nsub/nsubat){
              cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	      exit(EXIT_FAILURE);
	    }

//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            nktot+=iwght;
            ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
            if (ifail != 0)
              cout<<i<<" "<<j<<endl;

            zresu = zresu + zconu*static_cast<double>(iwght);
            zresd = zresd + zcond*static_cast<double>(iwght);
            zresud = zresud + zconud*static_cast<double>(iwght);
            zresdu = zresdu + zcondu*static_cast<double>(iwght);
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	  }
	}
      }
//     -----------------------------------------------------------------
//     PRIMITIVE-RECTANGULAR
      if (irecip == 2){
        d1=b1/2.;
        d2=b2/2.;
        for (int i=0; i<=max; i++){
          x=dq*i;
          for (int j=0; j<=max; j++){
            y=dq*j;
            xk=x*d1+y*d2;
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//         determine weights
            iwght=4;
            if (i == 0)
	      iwght=2;
            if (j == 0)
	      iwght=2;
            if (i == max)
	      iwght=2;
            if (j == max)
	      iwght=2;
            if ((i == 0) && (j == 0))
	      iwght=1;
            if ((i == 0) && (j == max))
  	      iwght=1;
            if ((i == max) && (j == 0)) 
	      iwght=1;
            if ((i == max) && (j == max))
	      iwght=1;
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//         find the folding points

            nxfold = testk(x,y,b1,b2,ifold,xfold,nfold,baib,irecip);
            if (nxfold != nsub/nsubat){
              cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	      exit(EXIT_FAILURE);
	    }

//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            nktot+=iwght;
            ifail=0;
            ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
            if (ifail != 0)
              cout<<i<<" "<<j<<endl;
//
            zresu = zresu + zconu*static_cast<double>(iwght);
            zresd = zresd + zcond*static_cast<double>(iwght);
            zresud = zresud + zconud*static_cast<double>(iwght);
            zresdu = zresdu + zcondu*static_cast<double>(iwght);
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	  }
	}
      }
//     -----------------------------------------------------------------
//     HEXAGONAL
      if (irecip == 3){
        d1=b1/2.;
        d2=(b1+b2)/3.;
        for (int i=0; i<=max; i++){
          x=dq*i;
          for (int j=0; j<=i; j++){
            y=dq*j;
            xk=x*d1+y*d2;
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//         determine weights
            iwght=12;
            if (i == j)
	      iwght=6;
            if (j == 0)
	      iwght=6;
            if (i == max)
	      iwght=4;
            if (i == 0)
	      iwght=1;
            if (j == max)
	      iwght=2;
            if ((i == max) && (j == 0))
	      iwght=3;
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//         find the folding points

            nxfold = testk(x,y,b1,b2,ifold,xfold,nfold,baib,irecip);
            if (nxfold != nsub/nsubat){
              cout<<"ERROR SUMK : nxfold is not equal to nsub/nsubat "<<nxfold<<" "<<nsub<<" "<<nsubat<<endl;
	      exit(EXIT_FAILURE);
	    }

//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            nktot+=iwght;
            ifail=0;
            ifail = cond(zener,xk,zconu,zcond,zconud,zcondu,nsub,nsubat,nxfold,xfold,forward<Args>(params)...);
            if (ifail != 0)
              cout<<i<<" "<<j;
// 
            zresu = zresu + zconu*static_cast<double>(iwght);
            zresd = zresd + zcond*static_cast<double>(iwght);
            zresud = zresud + zconud*static_cast<double>(iwght);
            zresdu = zresdu + zcondu*static_cast<double>(iwght);
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	  }
	}
      }
//     -----------------------------------------------------------------
//     CENTRED-RECTANGULAR
     if (irecip == 4){
//     ?????????????????????????????????????????????????????????????????
      cout<<"CENTRED-RECTANGULAR NOT CODED YET"<<endl;
      exit(EXIT_FAILURE);
     }
//     ?????????????????????????????????????????????????????????????????
      zresu = zresu/(1.*nktot);
      zresd = zresd/(1.*nktot);
      zresud = zresud/(1.*nktot);
      zresdu = zresdu/(1.*nktot);

      return;
}

template <typename... Args>
void kcon(int nsubat, const vector<pair<int,int>> &ifold, int nfold,const Matrix3d &baib, int nsub, int ndiff, double fact, 
	dcomp &zresu, dcomp &zresd, dcomp &zresud, dcomp &zresdu,
	Args&&... params){

      double xc, xcold;
      double xcon=1e-7;

      int nk;
      double diff, diffc;
      int iamcon;
      for (int kk=1; kk<=10; kk++){
        nk=1+pow(2,(kk+2));
        zresu = 0.;
        zresd = 0.;
        zresud = 0.;
        zresdu = 0.;

//       calculate the function
        sumk(nsub,nsubat,ifold,nfold,baib,ndiff,nk,zresu,zresd,zresud,zresdu,forward<Args>(params)...);
        xc=fact*real(zresu+zresd-zresud-zresdu)/nsub;

//       check for convergence
        diff=0.;
        iamcon=1;
        if (kk == 1)
	  iamcon=0;
	else{
          diffc=abs(xc-xcold);
          diff=max(diff,diffc);
          if (diff > xcon)
	    iamcon=0;
	}
        if (iamcon == 0){
          cout<<"NOT converged for nk = "<<nk<<"; xcon = "<<xcon<<"; error = "<<diff<<endl;
          if (nk > 150){
            cout<<"FORGETTING THIS ENERGY POINT !!!"<<endl;
            cout<<"K-POINT CONVERGENCE FAILED !! "<<diff<<endl;
	    break;
	  }

          xcold=xc;
	}
        if (iamcon == 1){
	  cout<<"converged for nk = "<<nk<<"; xcon = "<<xcon<<"; error = "<<diff<<endl;
	  break;
	}
      }
      return;
}
#endif
