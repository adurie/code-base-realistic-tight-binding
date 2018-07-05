#include <iostream>
/* #define EIGEN_USE_MKL_ALL */
/* #include "eigen3/Eigen/src/Core/util/MKL_support.h" */
#include <string>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <fstream>
#include <eigen3/Eigen/StdVector>
#include <utility>
#include <cstdlib>
/* #include "FeAgFe_new.h" */
#include "Fe_iron_dn_extract.h"
/* #include "cunningham_spawn.h" */
/* #include "cunningham_quad.h" */
/* #include "cunningham_multipoint.h" */
/* #include "cunningham_multipoint_par2.h" */
/* #include "cunningham_multipoint_par1.h" */
/* #include "cutest.h" */
/* #include "calc_spawn.h" */
/* #include "calc_par.h" */
#include "calc_dos.h"
/* #include "calc_cunningham.h" */
/* #include "dos.h" */
#include "dos_new.h"


//     This program calculates the coupling for a general multilayer,
//     with general supercell configuration, and general growth direction.
//     Each layer is allowed to have a different geometry.
//     However, the in-plane geometry (defined by (a1,a2)) of the lattice
//     must be common to all layers. This is so that the 2 
//     geometries can be simultaneously diagonalised in-plane.
//
//     There are nlay layers :
//     Layers 1 and 2 define the LH semi-infinite substrate.
//     Layers nlay-1 and nlay define the RH semi-infinite substrate.      
//     Layers 3..(nlay-2) define the intermediate multilayer. 
//     The ith multilayer has perpendicular position given by a3(i),
//     and basis vsub(i,1-->nsub). So the lattice for the ith layer is
//           R = \sum{n1,n2,m} [n1.a1 + n2.a2 + a3(i)  + vsub(i,m)]

//     There are significant changes to subroutine param
//     The hopping between different atoms is now given by the geometric
//     mean, unless this is overridden --- as in the case of Fe-MgO.
//     There are two spins per atom type.

//________                 _____________
//... x  x| o  o ... o    |  +     + ...
//________|               |_____________
//    1  2  3  4  nlay-2  nlay-1  nlay

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
typedef vector<Vector3d, aligned_allocator<Vector3d>> vV3d;
typedef vector<vector<Vector3d, aligned_allocator<Vector3d>>> vvV3d;
typedef vector<Matrix2d, aligned_allocator<Matrix2d>> vm2d;

bool WayToSort(double i, double j){ return abs(i) < abs(j);}

int recip(const Vector3d &a1, const Vector3d &a2, Vector3d &b1, Vector3d &b2){
//     construct the third (normalised) lattice vector
      Vector3d d;
      d = a1.cross(a2);
      d.normalize();
//     -----------------------------------------------------------------
//     now determine the reciprocal lattice vectors
      Vector3d an1, an2;
      an1=a1;
      an2=a2;
      b1 = an2.cross(d);
      b2 = an1.cross(d);
      double an1b1, an2b2;
      an1b1 = an1.dot(b1);
      an2b2 = an2.dot(b2);
      b1 = 2*M_PI*b1/an1b1;
      b2 = 2*M_PI*b2/an2b2;
//     now determine whether b2 or -b2 has the smallest angle with b1
      double b1b1, b2b2, b1b2;
      b1b1 = b1.dot(b1);
      b2b2 = b2.dot(b2);
      b1b2 = b1.dot(b2);
      if (b1b2 < 0)
        b2=-b2;
      b1b2=-b1b2;
      double cos12;
      cos12 = acos(b1b2/sqrt(b1b1*b2b2));
//     -----------------------------------------------------------------
//     now classify the reciprocal lattice
      cout<<endl<<endl;
      int irecip=0;
      if ((abs(b1b2) < 1e-8) && (abs(b1b1-b2b2) < 1e-8)){
        cout<<"reciprocal lattice is cubic"<<endl;
        irecip=1;
      }
      else if ((abs(b1b2) < 1e-8) && (abs(b1b1-b2b2) > 1e-8)){
        cout<<"recip lattice is primitive-rectangular"<<endl;
        irecip=2;
      }
      else if ((abs(cos12-M_PI/3.) < 1e-8) && (abs(b1b1-b2b2) < 1e-8)){
        cout<<"reciprocal lattice is hexagonal"<<endl;
        irecip=3;
      }
      else if ((abs(b1b2) > 1e-8) && (abs(b1b1-b2b2) < 1e-8)){
        cout<<"recip lattice is centred-rectangular"<<endl;
        irecip=4;
      }
      cout<<endl<<endl;
      return irecip;
}

vector<pair<int,int>> folding(Matrix3d &baib, const Vector3d &b1, const Vector3d &b2, const Vector3d &ba1, const Vector3d &ba2, int irecipa, int nsub, int &nfold){
//     check {b1,b2} and {ba1,ba2} in the same plane
      Vector3d b12, ba12;
      b12 = b1.cross(b2);
      ba12 = ba1.cross(ba2);
      if ((abs(b12.dot(ba1)) > 1e-10) || (abs(b12.dot(ba2)) > 1e-10)){
        cout<<"{b1,b2} not in same plane as {ba1,ba2}"<<endl;
	exit(EXIT_FAILURE);
      }
//     -----------------------------------------------------------------
//     construct the transformation matrix from {ba1,ba2} --> {b1,b2}
      Matrix3d bbatmp, bb, bba;
      bb.leftCols(1)        =b1;
      bb.block<3,1>(0,1)    =b2;
      bb.rightCols(1)       =b12;
      bbatmp.leftCols(1)    =ba1;
      bbatmp.block<3,1>(0,1)=ba2;
      bbatmp.rightCols(1)   =ba12;
      bba = bbatmp.inverse();
      baib = bba*bb;
//     -----------------------------------------------------------------
      vector<pair<int,int>> ifold;
      if (irecipa == 0){
        cout<<"TESTK : irecip = 0 ---- NOT CODED"<<endl;
	exit(EXIT_FAILURE);
      }
      if (irecipa == 3){
        cout<<"HEXAGONAL ATOMIC RECIPROCAL LATTICE NOT CODED"<<endl;
        exit(EXIT_FAILURE);
      }
      if (irecipa == 4){
        cout<<"CENTRED-RECTANGULAR ATOMIC RECIPROCAL LATTICE NOT CODED"<<endl;
        exit(EXIT_FAILURE);
      }
      if ((irecipa == 1) || (irecipa == 2)){
//     -----------------------------------------------------------------
//     CUBIC AND RECTANGULAR ATOMIC RECIPROCAL LATTICE
//     first construct supercell reciprocal lattice cluster
//     for the vector J = i*b1 + j*b2 + vsym(k)
//     then           J = x1*ba1 + x2*ba2
        int icnt=0;
        Vector3d tmp1, tmp2;
	double xsym, ysym, x1, x2;
	int cond;
	pair<int,int> foldtmp;
        for (int isym=1; isym <= 4; isym++){
          if (isym == 1){
            xsym=0.;
            ysym=0.;
	  }
          if (isym == 2){
            xsym=0.5;
            ysym=0.;
	  }
          if (isym == 3){
            xsym=0.;
            ysym=0.5;
  	  }
          if (isym == 4){
            xsym=0.5;
            ysym=0.5;
	  }
          for (int i=-nsub; i<=nsub; i++){
            for (int j=-nsub; j<=nsub; j++){
	      cond = 0;
              x1=baib(0,0)*(i+xsym)+baib(0,1)*(j+ysym);
              x2=baib(1,0)*(i+xsym)+baib(1,1)*(j+ysym);

//           check consistency
              tmp1=(i+xsym)*b1+(j+ysym)*b2;
              tmp2=x1*ba1+x2*ba2;
              for (int l=0; l<3; l++){
                if (abs(tmp1(l)-tmp2(l)) > 1e-10){
                  cout<<"ERROR FOLDING "<<tmp1<<" "<<tmp2<<endl;
		  exit(EXIT_FAILURE);
		}
	      }

//           check that this vector lies inside the atomic BZ
              if((abs(x1) <= 0.50000001) && (abs(x2) <= 0.50000001)){
//           check if this point the same as any of the others
                for (int iii=0; iii<icnt; iii++){
                  if((i == ifold[iii].first) && (j == ifold[iii].second)){
		    cond = 1;
		    break;
	          }
	        }

	        if (cond == false){
		  foldtmp = make_pair(i, j);
		  ifold.emplace_back(foldtmp);
                  icnt++;
		}
	      }
	    }
	  }
	}
        nfold=icnt;
      }

      cout<<endl<<"FOLDING VECTORS"<<endl;
      for ( int i=0; i < nfold; i++)
        cout<<i+1<<" "<<ifold[i].first<<" "<<ifold[i].second<<endl;
      cout<<endl<<endl;
      return ifold;
}

vector<double> prestructij(int ilay, int jlay, int nsub, const vvV3d &vsub, 
		int numnn, const Vector3d &a1, const Vector3d &a2, const vV3d &a3){
//     find the NN distances from layer ilay to the layer jlay.
//     -----------------------------------------------------------------
//     INPUT : {a1,a2,a3} - primitive lattice vectors
//             nsub       - number of sub-lattice vectors
//             vsub       - the set of nsub sublattice vectors
//                          {vsub(1),...,vsub(nsub)}
//             numnn      - upto order numnn.th N.N interactions
//
//     OUTPUT : dist(i)   - the distance from the origin to the ith N.N
//     -----------------------------------------------------------------
//     Construct vectors from (isub,ilay) --> (jsub,jlay)

      vector<double> disttmp;
      disttmp.reserve(numnn);
      for (int i = 0; i<numnn; i++)
        disttmp.emplace_back(1e10);
      Vector3d v;
      double dv;
      int kmax;
      int tempo;

      for (int i=-numnn; i<=numnn; i++){
        for (int j=-numnn; j<=numnn; j++){
          for (int isub=0; isub<nsub; isub++){
            for (int jsub=0; jsub<nsub; jsub++){
              v=(vsub[ilay-1][isub]+a3[ilay-1])-(vsub[jlay-1][jsub]+a3[jlay-1]+i*a1+j*a2);
	      dv = v.hypotNorm();

//   check dv not already in disttmp and replace the largest element it is smaller than
              if (dv > 1e-8){
		tempo = 1;
                for (int k=0; k<numnn; k++){
                  if (abs(dv-disttmp[k]) < 1e-10){
		    tempo = 0;
		    break;
		  }
		}
		if (tempo == 0)
		  continue;
                kmax=0;
                for (int k=1; k<numnn; k++){
                  if (disttmp[k] > disttmp[kmax])
		    kmax=k;
		}
                if (dv < disttmp[kmax])
		  disttmp[kmax]=dv;
	      }
	    }
	  }
	}
      }
//     -----------------------------------------------------------------
//     sort into distance from origin
      sort(disttmp.begin(), disttmp.end(), WayToSort);
//     -----------------------------------------------------------------
      return disttmp;
}

int main(){

//     DATA
      const int nspin=9;    // Number of energy bands
      const int numnn=2;    // No of nearest neighbours 

      const int nins=0;    // No of spacer principal layers
      const int mlay=0;     // No of substrate layers on each side of SGF
      const int numat=2;    // No of atom types: one for each element

      const int nlay=nins+2*mlay+4;  // total No of layers inc 4 lead layers

      const double ef=0.4635;  // fermi energy
      double wm = 1e-4;  // infinitesimal imaginary contribution to energy

      /* const int iwmax = 15;  // No of Matsubara frequencies */
      const double tfac     = 8.617e-5/13.6058;
      const double temp  = 300*tfac;
      /* const double temp  = 315.79*tfac; */

//     =================================================================
//     ATOMIC DATA FOR LEADS
//     =================================================================
//     in-plane atomic lattice vectors for substrate
//     in this code we assume that the LH and RH leads have the same 
//     lattice vectors and No of basis vectors

      const int nsubat=2; // No. of sublayer atoms in leads
      const int natom=nspin*nsubat;
      Vector3d aa1, aa2;

      aa1 << 1, 0, 0;
      aa2 << 0, 1, 0;

//     =================================================================
//     LH LEAD BASIS VECTORS
//     =================================================================
      vector<Vector3d, aligned_allocator<Vector3d>> aa3, a3;
      aa3.reserve(4);
      a3.reserve(nlay);
      vector<Vector3d, aligned_allocator<Vector3d>> vsubtmpFe, vsubtmpAg, vsubattmp;
      vsubattmp.reserve(nsubat);
      vector<vector<Vector3d, aligned_allocator<Vector3d>>> vsubat, vsub;
      vsubat.reserve(4);
      vsub.reserve(nlay);
      Vector3d tmp;
      VectorXd attype(nsubat);
      vector<VectorXd, aligned_allocator<VectorXd>> itypeat;

//       Sublattice
      tmp << 0, 0, 0;
      vsubattmp.emplace_back(tmp);
      tmp << 0.5, 0.5, 0.5;
      vsubattmp.emplace_back(tmp);
      for (int ilay=1; ilay<=2; ilay++){
//       Out of plane lattice vector
	tmp << 0, 0, ilay - 1;
	aa3.emplace_back(tmp);

	vsubat.emplace_back(vsubattmp);

//       Atom types
        for (int kk=1; kk <= nsubat; kk++)
	  attype(kk-1) = 1;
	itypeat.emplace_back(attype);
      }

//     =================================================================
//     SUPERCELL STRUCTURE
//     =================================================================

      const int nsub=2; //No. of sublayer atoms in spacer
      vsubtmpAg.reserve(nsub);
      vsubtmpFe.reserve(nsub);

//     2 in-plane lattice vectors
//     CUBIC
      Vector3d a1, a2;

      a1 << 1, 0, 0;
      a2 << 0, 1, 0;

//    ----------  Crystal structure ----------------------
      vector<vector<int>> itype;
      itype.reserve(nlay);
      vector<int> itmp;
      itmp.reserve(nsub);
//       Sublattice
      tmp << 0, 0, 0;
      vsubtmpAg.emplace_back(tmp);
      vsubtmpFe.emplace_back(tmp);
      tmp << 0.5, 0.5, 0.5;
      vsubtmpFe.emplace_back(tmp);
      tmp << 0.5, 0.5, 1./M_SQRT2;
      vsubtmpAg.emplace_back(tmp);
      for (int ilay=1; ilay <= nlay; ilay++){
//       Out of plane lattice vector

//       Atom types
        if((ilay == 3) && (nins > 0)){
	  vsub.emplace_back(vsubtmpAg);
	  tmp << 0, 0, (ilay - 1);
	  a3.emplace_back(tmp);
	  for (int isub = 0; isub < nsub; isub++)
	    itmp.emplace_back(2); // Ag 
	  itype.emplace_back(itmp);
	  itmp.clear();
	}
//this for odd layers
	/* if(ilay == 3){ */               
	/*   vsub.emplace_back(vsubtmpFe); */
	/*   tmp << 0, 0, (ilay - 1); */
	/*   a3.emplace_back(tmp); */
	/*   itmp.emplace_back(1); // Fe */
	/*   itmp.emplace_back(2); // Ag */
	/*   itype.emplace_back(itmp); */
	/*   itmp.clear(); */
	/* } */
	else if((ilay > 3) && (ilay <= nins+2)){
	  vsub.emplace_back(vsubtmpAg);
	  /* tmp << 0, 0, M_SQRT2*(ilay - 4) + 2 + (M_SQRT2 + 1)/2.; //This for odd layers */
	  /* tmp << 0, 0, M_SQRT2*(ilay - 3) + 2; //This for even layers */
	  tmp << 0, 0, ilay - 1;
	  a3.emplace_back(tmp);
	  for (int isub = 0; isub < nsub; isub++)
	    itmp.emplace_back(2); // Ag 
	  itype.emplace_back(itmp);
	  itmp.clear();
	}
	else if (ilay < 3){
	  vsub.emplace_back(vsubtmpFe);
	  tmp << 0, 0, ilay - 1;
	  a3.emplace_back(tmp);
	  for (int isub = 0; isub < nsub; isub++)
	    itmp.emplace_back(1); //  Fe
	  itype.emplace_back(itmp);
	  itmp.clear();
	}
	else {
	  vsub.emplace_back(vsubtmpFe);
	  /* tmp << 0, 0, (nins - 1)*M_SQRT2 - nins - 1 + (M_SQRT2 + 1)/2. + ilay; //This for even layers */
	  /* tmp << 0, 0, (nins - 1)*M_SQRT2 - nins + ilay; //This for odd layers */
	  tmp << 0, 0, ilay - 1;
	  a3.emplace_back(tmp);
	  for (int isub = 0; isub < nsub; isub++)
	    itmp.emplace_back(1); //  Fe
	  itype.emplace_back(itmp);
	  itmp.clear();
	}
      }

//     =================================================================
//     RH LEAD BASIS VECTORS
//     =================================================================

      for (int ilay=nlay-1; ilay <= nlay; ilay++){
//       Out of plane lattice vector
	/* tmp << 0, 0, (nins - 1)*M_SQRT2 - nins - 1 + (M_SQRT2 + 1)/2. + ilay; //This for even layers */
	tmp << 0, 0, ilay - 1;
	/* tmp << 0, 0, (nins - 1)*M_SQRT2 - nins + ilay; //This for odd layers */
	aa3.emplace_back(tmp);

	vsubat.emplace_back(vsubattmp);

//       Atom types
        for (int kk=1; kk <= nsubat; kk++)
	  attype(kk-1) = 1;
	itypeat.emplace_back(attype);
      }

//     =================================================================
//     The map between Supercell sublattice and LH atomic sublattice :
//     imap:   supercell --> atomic
//     Defined s.t.     vsub(k)=vsubat(imap(k)) + atomic lattice vector

      vector<int> imapl, imapr;
      for (int isub = 0; isub < nsub; isub++){
	if (isub == 0){
          imapl.emplace_back(1);
          imapr.emplace_back(1);
	}
	if (isub == 1){
          imapl.emplace_back(2);
          imapr.emplace_back(2);
	}
      }

//     =================================================================
//     THIS SECTION TO GET NN DISTANCES ONLY

//     In and out of plane distances:
      MatrixXd mdist(nlay-1, nlay-1);
      mdist.fill(0);
      vector<double> dist;
      dist.reserve(numnn);
      vector<MatrixXd, aligned_allocator<MatrixXd>> dnn;
      dnn.reserve(numnn);
      for (int i=1; i <= numnn; i++){

        for (int ilay=2; ilay <= nlay; ilay++){
          dist = prestructij(ilay,ilay, nsub, vsub, numnn, a1, a2, a3);
          mdist(ilay-2,ilay-2)=dist[i-1];
          dist = prestructij(ilay,ilay-1, nsub, vsub, numnn, a1, a2, a3);
	  if (ilay > 2){
            mdist(ilay-2,ilay-3)=dist[i-1];
            mdist(ilay-3,ilay-2)=dist[i-1];
	  }
	}

	/* cout<<mdist<<endl<<endl; */
	dnn.emplace_back(mdist);
      }

//     =================================================================
//     NN DISTANCES
//     1,2 Fe  ; 3 Mg ; 4 O 
//     ddnn(type1,type2,NN)

      // As in original code, ddnn is set up for an arbitrary
      // number of nearest neighbours, but the following code
      // only caters for up to two nearest neighbours
      vector<MatrixXd, aligned_allocator<MatrixXd>> ddnn;
      ddnn.reserve(numnn);
      MatrixXd ddnntmp(numat, numat);
      double M_SQRT3o2 = sqrt(3.)/2.;
      ddnntmp << M_SQRT3o2, M_SQRT3o2, M_SQRT3o2, 1.;
      /* ddnntmp << M_SQRT3o2, 1., 1., 1.; //this for odd layers? */
      ddnn.emplace_back(ddnntmp);
      ddnntmp << 1., 1., M_SQRT2/2. + 0.5, M_SQRT2; //original
      /* ddnntmp << 1., 1., 1., M_SQRT2; // this for even layers? */
      /* ddnntmp << 1., 1.20711, 1.20711, M_SQRT2; // this for odd layers? */
      ddnn.emplace_back(ddnntmp);

//     =================================================================
//     !!!!!!! OUTPUT ATOMIC POSITIONS FOR RASMOL VIEWER !!!!!!
//     load into rasmol with command :    > load xyz pos0.dat
      vector<string> atname;
      atname.reserve(2);
      atname.emplace_back("Co");
      atname.emplace_back("Cu");

//     whole cluster
      ofstream Myfile;	
      string Mydata = "pos0.dat";
      Myfile.open( Mydata.c_str(),ios::trunc );

      int idum0=0;

      Vector3d rr;
      for (int ilay = 0; ilay < nlay; ilay++){
	for (int iii = 0; iii < nsub; iii++){
          for (int i1=-nlay; i1 <= nlay; i1++){
            for (int i2=-nlay; i2 <= nlay; i2++){
              rr=a3[ilay]+vsub[ilay][iii]+i1*a1+i2*a2;
              if ((abs(rr(0)) < 1.3001) && (abs(rr(1)) < 1.3001))
                idum0++;
	    }
      	  }
	}
      }

      Myfile<<idum0<<endl<<"foo"<<endl;

      for (int ilay = 0; ilay < nlay; ilay++){
	for (int iii = 0; iii < nsub; iii++){
          for (int i1=-nlay; i1 <= nlay; i1++){
            for (int i2=-nlay; i2 <= nlay; i2++){
              rr=a3[ilay]+vsub[ilay][iii]+i1*a1+i2*a2;
              if ((abs(rr(0)) < 1.3001) && (abs(rr(1)) < 1.3001))
                Myfile<<atname[(itype[ilay][iii]) - 1]<<" "<<4*rr.transpose()<<endl;
	    }
	  }
	}
      }
      Myfile.close();

//     =================================================================
//     construct the growth direction:
      Vector3d xn;
      xn = a1.cross(a2);
      xn.normalize();

//     =================================================================

      cout<<"Supercell GMR"<<endl<<endl;
      cout<<"Substrate lattice vectors"<<endl;
      cout<<fixed<<a1.transpose()<<endl;
      cout<<a2.transpose()<<endl<<endl<<endl;
      cout<<"Growth direction"<<endl;
      cout<<xn.transpose()<<endl<<endl;
      cout<<"ef = "<<ef<<endl<<endl;
      cout<<"No of Cu at. layers = "<<nins<<endl<<endl;

//       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//       this prints out all interplanar distances
      for (int ii = 0; ii < nlay; ii++){
	cout<<string(72, '-')<<endl;
        cout<<"The "<<ii+1<<"'th Layer"<<endl<<endl;
        cout<<a3[ii].transpose()<<endl<<endl;
        cout<<"Sub-lattice vectors"<<endl;
	for (int iii = 0; iii < nsub; iii++)
	  cout<<"vsub("<<iii+1<<")= "<<(vsub[ii][iii]).transpose()<<endl;
        if (ii >= 1){
          cout<<endl<<"In-Plane NN distances"<<endl;
          for (int i=0; i < numnn; i++)
            cout<<ii+1<<" "<<dnn[i](ii-1,ii-1)<<endl;
          cout<<endl<<"Out of Plane NN distances"<<endl;
          for (int i=0; i < numnn; i++){
	    if (ii == 1)
              cout<<ii+1<<" "<<dnn[i](ii-1,ii)<<endl;
	    else
              cout<<ii+1<<" "<<dnn[i](ii-1,ii-2)<<endl;
	  }
	}
      }
      cout<<string(72, '-')<<endl;

//     -----------------------------------------------------------------
//     determine the reciprocal lattice structure
      Vector3d b1, b2;
      b1.fill(0);
      b2.fill(0);
      int irecip = recip(a1,a2,b1,b2);

      cout<<endl<<"perpr reciprocal lattice vectors"<<endl;
      cout<<b1.transpose()<<endl;
      cout<<b2.transpose()<<endl<<endl;
//     -----------------------------------------------------------------

//     determine the atomic reciprocal lattice vectors ba1,ba2
      Vector3d ba1, ba2;
      ba1.fill(0);
      ba2.fill(0);
      int irecipa = recip(aa1,aa2,ba1,ba2);
      cout<<string(52, '-')<<endl<<endl;
      cout<<"perpr atomic reciprocal lattice vectors"<<endl;
      cout<<ba1.transpose()<<endl;
      cout<<ba2.transpose()<<endl;
      cout<<endl<<string(52, '-')<<endl<<endl;

//     -----------------------------------------------------------------
//     now find the 'folding' vectors j[1:nfold]
      int nfold;
      Matrix3d baib;
      vector<pair<int,int>> ifold;
      ifold = folding(baib,b1,b2,ba1,ba2,irecipa,nsub,nfold);

//     =================================================================
//     CHECK THAT IMAP IS CORRECTLY DEFINED
//     ie.    vsub(k)=vsubat(imap(k)) + n1.aa1 + n2.aa2

//     LH lead
      Vector3d vtmp;
      double vba1, vba2;
      for (int ilay=1; ilay<=2; ilay++){
        for (int isub=1; isub<=nsub; isub++){
          vtmp=vsub[ilay - 1][isub - 1] - vsubat[ilay - 1][imapl[isub - 1]-1];
          vba1=(abs(vtmp.dot(ba1))+1e-10)/(2*M_PI);
          vba2=(abs(vtmp.dot(ba2))+1e-10)/(2*M_PI);
          if ((abs(floor(abs(vba1))-abs(vba1)) > 1e-10) || (abs(floor(abs(vba2))-abs(vba2)) > 1e-10)){
            cout<<"imapl is wrong "<<vba1<<" "<<vba2<<endl;
	    exit(EXIT_FAILURE);
	  }
          if (itype[ilay - 1][isub - 1] != (itypeat[ilay - 1])(imapl[isub - 1] - 1)){
            cout<<"imapl is wrong "<<itype[ilay - 1][isub - 1]<<" "<<(itypeat[ilay - 1])(imapl[isub - 1] - 1)<<endl;
 	    exit(EXIT_FAILURE);
	  }
	}
      }

//     RH lead
      for (int ilay=nlay-1; ilay<=nlay; ilay++){
        for (int isub=1; isub<=nsub; isub++){
          vtmp=vsub[ilay - 1][isub - 1] - vsubat[ilay - nlay + 3][imapr[isub - 1]-1];
          vba1=(abs(vtmp.dot(ba1))+1e-10)/(2*M_PI);
          vba2=(abs(vtmp.dot(ba2))+1e-10)/(2*M_PI);
          if ((abs(floor(abs(vba1))-abs(vba1)) > 1e-10) || (abs(floor(abs(vba2))-abs(vba2)) > 1e-10)){
            cout<<"imapr is wrong "<<vba1<<" "<<vba2<<endl;
	    exit(EXIT_FAILURE);
	  }
          if (itype[ilay - 1][isub - 1] != (itypeat[ilay - nlay + 3])(imapr[isub - 1] - 1)){
            cout<<"imapr is wrong "<<itype[ilay - 1][isub - 1]<<" "<<(itypeat[ilay - nlay + 3])(imapr[isub - 1] - 1)<<endl;
 	    exit(EXIT_FAILURE);
	  }
	}
      }

//     =================================================================
//     DO THE CALCULATION
      Matrix2d s0, p0, d0t, d0e;
      vm2d sssint, spsint, ppsint, pppint, sdsint, pdsint, pdpint, ddsint, ddpint, dddint;
      sssint.reserve(2); spsint.reserve(2); ppsint.reserve(2); pppint.reserve(2); sdsint.reserve(2); 
      pdsint.reserve(2); pdpint.reserve(2); ddsint.reserve(2); ddpint.reserve(2); dddint.reserve(2);
      param(numat, numnn, s0, p0, d0t, d0e, sssint, spsint, ppsint, pppint, sdsint, pdsint, pdpint, ddsint, ddpint, dddint);

      int ndiff = 0;

      /* double start = -1; */
      double start = 0.4196;
      double end = 1.;
      double step = 0.0026;

      string Mydata_up = "Fe_spin_up_new_extract_ignore.txt";
      string Mydata_down = "Fe_spin_down_sk.txt";
      string Mydata_total = "Fe_tdos_extract_ignore.txt";
      ofstream Myfile_up, Myfile_down, Myfile_total;	
      Myfile_up.open( Mydata_up.c_str(),ios::app );
      Myfile_down.open( Mydata_down.c_str(),ios::app );
      Myfile_total.open( Mydata_total.c_str(),ios::app );
      double up_result;
      double down_result;

      double zresu, zresd;

      double fact = (M_PI + M_PI)*temp;
      dcomp ec, i;
      i = -1.;
      i = sqrt(i);
      int nmat = nspin*nsub;
      for (double j = start; j<end + step; j=j+step){
        ec = j + i*wm;
        kcon(nsubat,ifold,nfold,baib,nsub,ndiff,fact,zresu,zresd,irecip,b1,b2,ec,nmat,mlay,nins,nlay,
  	  nspin,imapl,imapr,vsub,vsubat,numnn,a1,a2,a3,aa1,aa2,aa3,itype,itypeat,ddnn,s0, p0, d0t, d0e, sssint, spsint, ppsint, pppint, sdsint,
	  pdsint, pdpint, ddsint, ddpint, dddint);
	  up_result = zresu*(-1./(2.*M_PI));
	  down_result = zresd*(-1./(2.*M_PI));
	  Myfile_up<<j<<" "<<up_result<<endl;
  	  Myfile_down<<j<<" "<<down_result<<endl;
  	  Myfile_total<<j<<" "<<up_result + down_result<<endl;
      }
      Myfile_up.close();
      Myfile_down.close();
      Myfile_total.close();

      return 0;
}
