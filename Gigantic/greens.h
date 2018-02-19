#ifndef GREENS_H 
#define GREENS_H 

#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>
#include <string>

using namespace std;
using namespace Eigen;
typedef Matrix2d m2d;
typedef complex<double> dcomp;
typedef vector<Matrix2d, aligned_allocator<Matrix2d>> vm2d;
typedef vector<Vector3d, aligned_allocator<Vector3d>> vV3d;
typedef vector<vector<Vector3d, aligned_allocator<Vector3d>>> vvV3d;
typedef	vector<VectorXd, aligned_allocator<VectorXd>> vVXd;
typedef vector<MatrixXd, aligned_allocator<MatrixXd>> vMXd;

void adlayer1(MatrixXcd &zgl, MatrixXcd &zu, MatrixXcd &zt, dcomp zener, int n){
//     adlayer ontop of gl
      MatrixXcd zwrk(n, n), zgldag(n, n);
      zwrk = zgl*zt;
      zgl = zt.adjoint()*zwrk;

      for (int ir=0; ir<n; ir++){
        for (int is=0; is<=ir-1; is++){
          zgl(ir,is)=-zu(ir,is)-zgl(ir,is);
          zgl(is,ir)=-zu(is,ir)-zgl(is,ir);
	}
        zgl(ir,ir)=zener-zu(ir,ir)-zgl(ir,ir);
      }
//     find inverse to zginv
      zwrk = zgl.inverse();
      zgl = zwrk;
      return;
}

int surfacenew(MatrixXcd &zu, MatrixXcd &zt, dcomp zener, MatrixXcd &zsurfl, MatrixXcd &zsurfr, int natom){
//     this routine is written for the atomic cell SGF
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      int n=natom;                       
      int n2=2*n;

      dcomp zi=-1;
      zi = sqrt(zi);

      MatrixXcd zunit = MatrixXcd::Identity(n,n);

//     -----------------------------------------------------------------
//     -----------------------------------------------------------------
//     now calculate GF's from closed form
//     -----------------------------------------------------------------
//     -----------------------------------------------------------------
//     define zp
//     redefine gamma and delta
      MatrixXcd ztmp1(n,n), ztinv(n,n), zs(n,n), zsinv(n,n), zgamma(n,n);
      MatrixXcd zero = MatrixXcd::Zero(n,n);
      ztmp1 = zener*zunit - zu;
      ztinv = zt.inverse();
      zs = zt.adjoint();
      zsinv = zs.inverse();
      zgamma = ztmp1*ztinv;
      MatrixXcd zp(n2,n2), O(n2,n2);
      zp.topLeftCorner(n,n) = zero;
      zp.topRightCorner(n,n) = ztinv;
      zp.bottomLeftCorner(n,n) = -zs;
      zp.bottomRightCorner(n,n) = zgamma; 
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//     diagonalise zp
      ComplexEigenSolver<MatrixXcd> ces;
      ces.compute(zp);
      O = ces.eigenvectors();
      int ifail = ces.info();
      if (ifail != 0){
        cout<<"SURFACENEW : ifail = "<<ifail<<endl;
        ifail=1;
        return ifail;
      }
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//     calculate L.H. zsurf
      MatrixXcd b = O.topRightCorner(n, n);
      MatrixXcd d = O.bottomRightCorner(n, n);
      zsurfl = b*d.inverse();
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//     calculate R.H. zsurf
      MatrixXcd zevl, ztmp2, ztmp3;
      zevl = ces.eigenvalues().asDiagonal();
      ztmp1 = O.topLeftCorner(n,n);
      ztmp2 = ztmp1.inverse();
      ztmp1 = ztmp1*zevl.topLeftCorner(n,n);
      ztmp3 = ztmp1*ztmp2;
      zsurfr = ztmp3*zsinv;
//     -----------------------------------------------------------------
      ztmp1 = zsurfl;
      ztmp2 = zsurfr;
      adlayer1(ztmp1,zu,zt,zener,n);
      adlayer1(ztmp2,zu,zs,zener,n);
      /* double xminl=0.; */
      /* double xminr=0.; */
      /* double xmin1, xmin2; */
      /* for (int ir=0; ir<n; ir++){ */
      /*   for (int is=0; is<n; is++){ */
      /*     xmin1=abs(zsurfl(ir,is)-ztmp1(ir,is)); */
      /*     xmin2=abs(zsurfr(ir,is)-ztmp2(ir,is)); */
      /*     if (xmin1 > xminl) */
	    /* xminl=xmin1; */
      /*     if (xmin2 > xminr) */
      /*       xminr=xmin2; */
	/* } */
      /* } */
      zsurfl = ztmp1;
      zsurfr = ztmp2;
      /* double xmin = max(xminl,xminr); */
      /* if (xmin > 5e-5) */
        /* ifail=1; */

      return ifail;
}

Matrix<complex<double>, 9, 9> eint1(double sss, double sps, double pps, double ppp,
		     double ss, double ps, double pp, double ds, double dp,
		     double dd, double x, double y, double z)
{
//     THIS ROUTINE IS SPIN DEPENDENT :-
//     IT IS WRITTEN FOR s,p and d BANDS ONLY
      Matrix<complex<double>, 9, 9> b;
      double xx=x*x;
      double xy=x*y;
      double yy=y*y;
      double yz=y*z;
      double zz=z*z;
      double zx=z*x;
      double xxyy=xx*yy;
      double yyzz=yy*zz;
      double zzxx=zz*xx;
      double aux=pps-ppp;
      double r3=sqrt(3.);
      double aux1=r3*ss;
      double f8=3.*zz-1.;
      double f1=xx+yy;
      double f2=xx-yy;
      double f3=zz-.5*f1;
      double g1=1.5*f2*ds;
      double g2=r3*f3*ds;
      b(0,0)=sss;
      b(0,1)=x*sps;
      b(0,2)=y*sps;
      b(0,3)=z*sps;
      b(1,0)=-b(0,1);
      b(1,1)=xx*pps+(1.-xx)*ppp;
      b(1,2)=xy*aux;
      b(1,3)=zx*aux;
      b(2,0)=-b(0,2);
      b(2,1)=b(1,2);
      b(2,2)=yy*pps+(1.-yy)*ppp;
      b(2,3)=yz*aux;
      b(3,0)=-b(0,3);
      b(3,1)=b(1,3);
      b(3,2)=b(2,3);
      b(3,3)=zz*pps+(1.-zz)*ppp;
      b(0,4)=xy*aux1;
      b(0,5)=yz*aux1;
      b(0,6)=zx*aux1;
      b(0,7)=.5*f2*aux1;
      b(0,8)=.5*f8*ss;
      b(4,0)=b(0,4);
      b(5,0)=b(0,5);
      b(6,0)=b(0,6);
      b(7,0)=b(0,7);
      b(8,0)=b(0,8);
      double f4=.5*r3*f2*ps;
      double f5=.5*f8*ps;
      double aux2=r3*xx*ps+(1.-2.*xx)*pp;
      b(1,4)=aux2*y;
      b(1,5)=(r3*ps-2.*pp)*xy*z;
      b(1,6)=aux2*z;
      b(1,7)=(f4+(1.-f2)*pp)*x;
      b(1,8)=(f5-r3*zz*pp)*x;
      double aux3=(r3*yy*ps+(1.-2.*yy)*pp);
      b(2,4)=aux3*x;
      b(2,5)=aux3*z;
      b(2,6)=b(1,5);
      b(2,7)=(f4-(1.+f2)*pp)*y;
      b(2,8)=(f5-r3*zz*pp)*y;
      double aux4=r3*zz*ps+(1.-2.*zz)*pp;
      b(3,4)=b(1,5);
      b(3,5)=aux4*y;
      b(3,6)=aux4*x;
      b(3,7)=(f4-f2*pp)*z;
      b(3,8)=(f5+r3*f1*pp)*z;
      b(4,1)=-b(1,4);
      b(5,1)=-b(1,5);
      b(6,1)=-b(1,6);
      b(7,1)=-b(1,7);
      b(8,1)=-b(1,8);
      b(4,2)=-b(2,4);
      b(5,2)=-b(2,5);
      b(6,2)=-b(2,6);
      b(7,2)=-b(2,7);
      b(8,2)=-b(2,8);
      b(4,3)=-b(3,4);
      b(5,3)=-b(3,5);
      b(6,3)=-b(3,6);
      b(7,3)=-b(3,7);
      b(8,3)=-b(3,8);
      b(4,4)=3.*xxyy*ds+(f1-4.*xxyy)*dp+(zz+xxyy)*dd;
      b(4,5)=(3.*yy*ds+(1.-4.*yy)*dp+(yy-1.)*dd)*zx;
      b(4,6)=(3.*xx*ds+(1.-4.*xx)*dp+(xx-1.)*dd)*yz;
      b(4,7)=(g1-2.*f2*dp+.5*f2*dd)*xy;
      b(4,8)=(g2-2.0*r3*zz*dp+.5*r3*(1.+zz)*dd)*xy;
      b(5,4)=b(4,5);
      b(5,5)=3.*yyzz*ds+(yy+zz-4.0*yyzz)*dp+(xx+yyzz)*dd;
      b(5,6)=(3.0*zz*ds+(1.0-4.0*zz)*dp+(zz-1.0)*dd)*xy;
      b(5,7)=(g1-(1.0+2.0*f2)*dp+(1.0+.5*f2)*dd)*yz;
      b(5,8)=(g2+r3*(f1-zz)*dp-.5*r3*f1*dd)*yz;
      b(6,4)=b(4,6);
      b(6,5)=b(5,6);
      b(6,6)=3.*zzxx*ds+(zz+xx-4.*zzxx)*dp+(yy+zzxx)*dd;
      b(6,7)=(g1+(1.-2.*f2)*dp-(1.-.5*f2)*dd)*zx;
      b(6,8)=(g2+r3*(f1-zz)*dp-.5*r3*f1*dd)*zx;
      b(7,4)=b(4,7);
      b(7,5)=b(5,7);
      b(7,6)=b(6,7);
      b(7,7)=.75*f2*f2*ds+(f1-f2*f2)*dp+(zz+.25*f2*f2)*dd;
      b(7,8)=.5*f2*g2-r3*zz*f2*dp+.25*r3*(1.+zz)*f2*dd;
      b(8,4)=b(4,8);
      b(8,5)=b(5,8);
      b(8,6)=b(6,8);
      b(8,7)=b(7,8);
      b(8,8)=f3*f3*ds+3.*zz*f1*dp+.75*f1*f1*dd;
      return b;
}

MatrixXcd sk(int ind1, int ind2, int nn, Vector3d &d, double dd, const Vector3d &xk, Vector3d &dpar, int nspin, int ispin,
		const m2d &s0, const m2d &p0, const m2d &d0t, const m2d &d0e, const vm2d &sssint, 
		const vm2d &spsint, const vm2d &ppsint, const vm2d &pppint, const vm2d &sdsint, 
		const vm2d &pdsint, const vm2d &pdpint, const vm2d &ddsint, const vm2d &ddpint, const vm2d &dddint){
//         calculate the cosine angles :
      MatrixXcd rt(nspin,nspin);
      rt.fill(0.);
      Vector3d c;
      if (dd > 1e-8)
        c=d/dd;
  
      double dk;
      dcomp zexdk, i;
      dk=xk.dot(dpar); 
      i = -1;
      i = sqrt(i);
      zexdk=cos(dk) + i*sin(dk);
      int elem;
      if (ispin == +1)
	elem = 1;
      if (ispin == -1)
	elem = 0;
      double g1,g2,g3,g4,g5,g6,g7,g8,g9,g10;
  
  //         find the hopping for this nn.th N.N in direction d
  //         ensure that this routine gives U for d = 0
      if (nn == 0){
        rt(0,0)=s0(ind1,elem);
        rt(1,1)=p0(ind1,elem);
        rt(2,2)=p0(ind1,elem);
        rt(3,3)=p0(ind1,elem);
    
        rt(4,4)=d0e(ind1,elem);
        rt(5,5)=d0t(ind1,elem);
        rt(6,6)=d0t(ind1,elem);
        rt(7,7)=d0t(ind1,elem);
        rt(8,8)=d0e(ind1,elem);
      }
      else{
	nn--;
        g1=sssint[ind1](ind2,nn);
        g2=spsint[ind1](ind2,nn);
        g3=ppsint[ind1](ind2,nn);
        g4=pppint[ind1](ind2,nn);
        g5=sdsint[ind1](ind2,nn);
        g6=pdsint[ind1](ind2,nn);
        g7=pdpint[ind1](ind2,nn);
        g8=ddsint[ind1](ind2,nn);
        g9=ddpint[ind1](ind2,nn);
        g10=dddint[ind1](ind2,nn);
        rt = eint1(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,c(0),c(1),c(2));
      }
      return rt*zexdk;
}

template <typename... Args>
MatrixXcd helement(int n1, int n2, int isub, int jsub, const Vector3d &xk, int ispin, int ndim,
	       	int nspin, int natom, int nmat, const vvV3d &vsub, const vvV3d &vsubat, int numnn,
	       	const Vector3d &a1, const Vector3d &a2, const vV3d &a3,	const Vector3d &aa1,
	       	const Vector3d &aa2, const vV3d &aa3, const vector<vector<int>> &itype,
	       	const vVXd &itypeat, const vMXd &ddnn, Args&&... params){

      int ind1, ind2;
      if (ndim == 1){
        ind1=itypeat[n1](isub)-1;
        ind2=itypeat[n2](jsub)-1;
      }
      if (ndim == 0){
        ind1=itype[n1][isub]-1;
        ind2=itype[n2][jsub]-1;
      }

      MatrixXcd zh(nspin,nspin);
      zh.fill(0.);
      
//     Calculate 
//     H_{n1,n2}(isub,jsub) =
//          \sum{d_{//}} exp(i k_{//} d_{//}) <isub,n1|H|jsub,n2,d_{//}>
//
//     To do this we need to calculate 
//     D = (vsub(isub,n1)+a3(n1)) - (vsub(jsub,n2)+a3(n2)+d//)
//     We need to sum over d_{//}

      double ddmax=0.;
      Vector3d dpar, d;
      double dd;
      int nn;
      for (int inn=0; inn<numnn; inn++)
        ddmax=max(ddmax,ddnn[inn](ind1,ind2)); // sometimes ddnn(ind1,ind2,numnn)=0
      for (int i1=-numnn; i1<=numnn; i1++){
        for (int i2=-numnn; i2<=numnn; i2++){

//         First construct d_{//}
          if (ndim == 1){
            dpar=i1*aa1+i2*aa2+aa3[n2]-aa3[n1];
            d=vsubat[n2][jsub]-vsubat[n1][isub]+dpar;
	  }
          if(ndim == 0){
            dpar=i1*a1+i2*a2+a3[n2]-a3[n1];
            d=vsub[n2][jsub]-vsub[n1][isub]+dpar;
	  }
          dd=sqrt(d(0)*d(0) + d(1)*d(1) + d(2)*d(2));

          if (dd > (ddmax+1e-8))
		  continue;   // this is not a NN
          if (abs(dd) < 1e-8){
            nn=0;
	    zh = zh + sk(ind1, ind2, nn, d, dd, xk, dpar, nspin, ispin, forward<Args>(params)...);

	  }
          for (int inn=0; inn<numnn; inn++){
            if (abs(dd-ddnn[inn](ind1,ind2)) < 1e-8){
              nn=inn+1;
	      zh = zh + sk(ind1, ind2, nn, d, dd, xk, dpar, nspin, ispin, forward<Args>(params)...);
	    }
	  }
	}
      }
      return zh;
}

template <typename... Args>
MatrixXcd hamil(const Vector3d &xk, int n1, int n2, int ispin, int n, int ncell, int nsubat, int nmat, int nspin, vector<int> &imapl,
	vector<int> &imapr, const vvV3d &vsub, const vvV3d &vsubat, Args&&... params){
//     compute the hamiltonian matrices zu and zt

//     calculate for given layers n1,n2 ; sublattice points isub,jsub ;
//     orbital indices isp,jsp the Hamiltonian element
//     <n1 ; isub ; isp | H | n2 ; jsub ; jsp>

      int natom = nspin*nsubat;
      MatrixXcd zt(nmat,nmat);
      MatrixXcd zh(nspin,nspin);
      int jjj, iii;

//     compute U or T
      for (int isub=0; isub<ncell; isub++){
        for (int jsub=0; jsub<ncell; jsub++){
          zh = helement(n1,n2,isub,jsub,xk,ispin,n,nspin,natom,nmat,vsub,vsubat,forward<Args>(params)...);
          for (int isp=0; isp<nspin; isp++){       // orbital elements
            for (int jsp=0; jsp<nspin; jsp++){
              iii=(isub)*nspin+isp;
              jjj=(jsub)*nspin+jsp;
              zt(iii,jjj)=zh(isp,jsp);
	    }
	  }
	}
      }

      return zt;
}

dcomp coupl(MatrixXcd &zgsl, MatrixXcd &zgsr, MatrixXcd &zt){
      int rows = zgsr.rows();
      dcomp zcoupl;
      MatrixXcd ztdag(rows,rows), ztmp1(rows,rows), ztmp2(rows,rows);
      ztdag = zt.adjoint();
      ztmp1 = zgsl*zt;
      ztmp2 = ztdag*ztmp1;
      MatrixXcd I = MatrixXcd::Identity(rows, rows);

      ztmp1 = zgsr*ztmp2 - I;
      zcoupl=log(ztmp1.determinant())/M_PI;

      return zcoupl;
}

template <typename... Args>
int green(dcomp zener, int ispin, string &side, MatrixXcd &zgl, MatrixXcd &zgr, 
	int nsub, int nsubat, int nlay, int nmat, int nxfold, vV3d &xfold, vV3d &xshift, int nspin, vector<int> &imapl, 
	vector<int> &imapr, const vvV3d &vsub, const vvV3d &vsubat, Args&&... params){

//     Calculate the SGF at a given k// ,
//     in the supercell representation -- so that the SGF's are nmatx x nmatx

//     if side="LH" we calculate the LH SGF. if side="RH" the RH SGF
//     if ispin=-1 minority, +1 majority

//     calculate the LH & RH SGF in the atomic basis, and convert to the 
//     supercell basis.

      int natom = nspin*nsubat;
      MatrixXcd zgltmp(nmat,nmat), zgrtmp(nmat,nmat);
      MatrixXcd ztat(natom,natom), zuat(natom,natom);
      zgltmp.fill(0.);
      zgrtmp.fill(0.);
      Vector3d ds;
      vector<int> imap;
      imap.reserve(2);
      int ifail = 0;
      int ilay;
      int iii, jjj, ii, jj;
      dcomp i, zarg1;
      i = -1.;
      i = sqrt(i);
      Vector3d vtmp;
      if (side == "LH"){
        ilay=0;
        imap=imapl;
      }
      if (side == "RH"){
        ilay=nlay-2;
        imap=imapr;
      }

      Vector3d xkat;
      MatrixXcd zglat(natom,natom), zgrat(natom,natom);
      zglat.fill(0.); zgrat.fill(0.);
      for (int iff=0; iff<nxfold; iff++){
        xkat=xfold[iff];

//       define atomic hamiltonian elements in lead   -  u1, t1
//       '1' after ispin relates to this referring to the lead
        ztat = hamil(xkat,ilay,ilay+1,ispin,1,nsubat,nsubat,nmat,nspin,imapl,imapr,vsub,vsubat,forward<Args>(params)...);
        zuat = hamil(xkat,ilay+1,ilay+1,ispin,1,nsubat,nsubat,nmat,nspin,imapl,imapr,vsub,vsubat,forward<Args>(params)...);

        ifail=0;
        ifail = surfacenew(zuat,ztat,zener,zglat,zgrat,natom);
        if (ifail != 0)// zt has a near zero eigenvalue
	  cout<<"eigenvalues ill-conditioned. Consider coding to higher precision"<<endl;

//       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//       now load elements into zgl and zgr
        for (int isub=0; isub<nsub; isub++){           // sublattice elements
          for (int jsub=0; jsub<nsub; jsub++){
            ds=vsub[ilay][isub]-vsub[ilay][jsub] - (vsubat[ilay][imap[isub]-1]-vsubat[ilay][imap[jsub]-1]);
            zarg1=i*xkat.dot(ds);
            vtmp=xshift[iff];

            for (int isp=0; isp<nspin; isp++){       // orbital elements
              for (int jsp=0; jsp<nspin; jsp++){
                iii=isub*nspin+isp;
                jjj=jsub*nspin+jsp;

                ii=(imap[isub]-1)*nspin+isp;
                jj=(imap[jsub]-1)*nspin+jsp;

                zgltmp(iii,jjj)=zgltmp(iii,jjj)+zglat(ii,jj)*exp(zarg1)*(nsubat*1.)/(1.*nsub);
                zgrtmp(iii,jjj)=zgrtmp(iii,jjj)+zgrat(ii,jj)*exp(zarg1)*(nsubat*1.)/(1.*nsub);
	      }
	    }
	  }
	}
      }
      zgl=zgltmp;
      zgr=zgrtmp;

      return ifail;
}

template <typename... Args>
int cond(dcomp zener, const Vector3d &xk, VectorXcd &zconu, VectorXcd &zcond, VectorXcd &zconud, VectorXcd &zcondu,
	int nsub, int nsubat, int nxfold, vV3d &xfold, vV3d &xshift, int nmat, int mlay, int nins, int nlay, Args&&... params){

//     Calculate the coupling at a given k// , via the det formula,
//     in the supercell representation -- so that the SGF's are nmatx x nmatx
//     =================================================================
//     CALCULATE SGFs IN ATOMIC BASIS
//     =================================================================
//     DO THIS IF LH AND RH LEADS ARE THE SAME
      int ifail = 0;
      string st = "LH";
      MatrixXcd zglu(nmat, nmat), zgru(nmat, nmat), zgld(nmat, nmat), zgrd(nmat, nmat);
      MatrixXcd zt(nmat, nmat), zu(nmat, nmat), ztdag(nmat, nmat), zudag(nmat, nmat);
      ifail = green(zener,+1,st,zglu,zgru,nsub,nsubat,nlay,nmat,nxfold,xfold,xshift,forward<Args>(params)...);   // LH UP
      if (ifail != 0)
	return ifail;
      ifail = green(zener,-1,st,zgld,zgrd,nsub,nsubat,nlay,nmat,nxfold,xfold,xshift,forward<Args>(params)...);   // LH DOWN
      if (ifail != 0)
	return ifail;
//     -----------------------------------------------------------------
//     DO THIS IF LH AND RH LEADS DIFFER
//     call green(zener,xk,+1,"LH",zglu,zfoo,ifail)   // LH UP
//     if(ifail.ne.0)return
//     call green(zener,xk,-1,"LH",zgld,zfoo,ifail)   // LH DOWN
//     if(ifail.ne.0)return

//     call green(zener,xk,+1,"RH",zfoo,zgru,ifail)   // RH UP
//     if(ifail.ne.0)return
//     call green(zener,xk,-1,"RH",zfoo,zgrd,ifail)   // RH DOWN
//     if(ifail.ne.0)return
//     -----------------------------------------------------------------

//     =================================================================
//     ALTERNATIVELY DO EVERYTHING IN SUPERCELL
//     note we do not normally keep surfacenewcell as a subroutine as 
//     it adds to the size of the code
//     To create surfacenewcell make changes to surfacenew as indicated
//     ifail=0
//     call hamil(zt,xk,1,2,+1,nmat,nmatx)
//     call hamil(zu,xk,2,2,+1,nmat,nmatx)
//     call surfacenewcell(zu,zt,zener,zglu,zgru,ifail)
//     call hamil(zt,xk,1,2,-1,nmat,nmatx)
//     call hamil(zu,xk,2,2,-1,nmat,nmatx)
//     call surfacenewcell(zu,zt,zener,zgld,zgrd,ifail)
//     =================================================================
//
//     adlayer the LH & RH mlay substrate layers (ie interface layers)
//     '0' after '+1' indicates this is required for the spacer
      for (int ill=2; ill<2+mlay; ill++){
        zt = hamil(xk,ill-1,ill,+1,0,nsub,nsubat,nmat,forward<Args>(params)...);
        zu = hamil(xk,ill,ill,+1,0,nsub,nsubat,nmat,forward<Args>(params)...);
        adlayer1(zglu,zu,zt,zener,nmat);

        zt = hamil(xk,ill-1,ill,-1,0,nsub,nsubat,nmat,forward<Args>(params)...);
        zu = hamil(xk,ill,ill,-1,0,nsub,nsubat,nmat,forward<Args>(params)...);
        adlayer1(zgld,zu,zt,zener,nmat);
      }

      for (int ill=nlay-3; ill>nlay-mlay-1; ill--){
        zt = hamil(xk,ill,ill+1,+1,0,nsub,nsubat,nmat,forward<Args>(params)...);
        zu = hamil(xk,ill,ill,+1,0,nsub,nsubat,nmat,forward<Args>(params)...);
        ztdag=zt.adjoint();
        adlayer1(zgru,zu,ztdag,zener,nmat);

        zt = hamil(xk,ill,ill+1,-1,0,nsub,nsubat,nmat,forward<Args>(params)...);
        zu = hamil(xk,ill,ill,-1,0,nsub,nsubat,nmat,forward<Args>(params)...);
        ztdag=zt.adjoint();
        adlayer1(zgrd,zu,ztdag,zener,nmat);
      }

//     =================================================================
//     CALCULATE GFs IN SPACER AND THE CONDUCTANCE
//     =================================================================
      for (int ill=0; ill<nins; ill++){

//       adlayer LH  ----   zgl
        zt = hamil(xk,ill+mlay+1,ill+mlay+2,+1,0,nsub,nsubat,nmat,forward<Args>(params)...);
        zu = hamil(xk,ill+mlay+2,ill+mlay+2,+1,0,nsub,nsubat,nmat,forward<Args>(params)...);
        adlayer1(zglu,zu,zt,zener,nmat);
        zt = hamil(xk,ill+mlay+1,ill+mlay+2,-1,0,nsub,nsubat,nmat,forward<Args>(params)...);
        zu = hamil(xk,ill+mlay+2,ill+mlay+2,-1,0,nsub,nsubat,nmat,forward<Args>(params)...);
        adlayer1(zgld,zu,zt,zener,nmat);

//       SPIN UP
        zt = hamil(xk,1+nins+mlay,2+nins+mlay,+1,0,nsub,nsubat,nmat,forward<Args>(params)...);
        zconu(ill) = coupl(zglu,zgru,zt);
	
//       SPIN DOWN
        zt = hamil(xk,1+nins+mlay,2+nins+mlay,-1,0,nsub,nsubat,nmat,forward<Args>(params)...);
        zcond(ill) = coupl(zgld,zgrd,zt);

//       SPIN UP-DOWN
        zt = hamil(xk,1+nins+mlay,2+nins+mlay,-1,0,nsub,nsubat,nmat,forward<Args>(params)...);
        zconud(ill) = coupl(zglu,zgrd,zt);

//       SPIN DOWN-UP
        zt = hamil(xk,1+nins+mlay,2+nins+mlay,+1,0,nsub,nsubat,nmat,forward<Args>(params)...);
        zcondu(ill) = coupl(zgld,zgru,zt);
      }
      return ifail;
}
#endif
