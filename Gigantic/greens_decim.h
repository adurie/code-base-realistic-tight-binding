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

void remve(MatrixXcd &zh0, vector<int> &lst, int n, int nlst){
      MatrixXcd zh(n,n);
      vector<int> notlst;
//     remove rows and columns contained in lst from zh0

//     first create a list of rows not to eliminate
      int icount=0;
      int condition = 0;
      for (int i = 0; i < n; i++){
	for (int j = 0; j < nlst; j++){
	  if (lst[j] == i) 
            condition = 1;
	}
	if (condition == 0){
	  icount++;
	  notlst.emplace_back(i);
	}
	condition = 0;
      }

      int ii, jj;
      for (int i = 0; i < icount; i++){
	for (int j = 0; j < icount; j++){
          ii = notlst[i];
	  jj = notlst[j];
	  zh(i,j) = zh0(ii,jj);
	}
      }

      if (icount != n-nlst)
        cout<<"ERROR REMVE: "<<icount<<" "<<n-nlst<<endl;

      zh0 = zh;

      return;
}

void decim(MatrixXcd &zh0, int k, int n){
      MatrixXcd zh(n,n);
//     Decimate row/column k from zh0
      for (int i = 0; i < n; i++){
	for (int j = 0; j < n; j++){
          zh(i,j)=zh0(i,j)-zh0(i,k)*zh0(k,j)/zh0(k,k);
	}
      }
      zh0 = zh;
      return;
}

void adlayer1(MatrixXcd &zgl, MatrixXcd &zu, MatrixXcd &zt, dcomp zener, int n){
//     adlayer ontop of gl
      MatrixXcd zwrk(n, n);
      zwrk = zgl*zt;
      zgl = zt.adjoint()*zwrk;
      MatrixXcd I = MatrixXcd::Identity(n, n);

      zwrk = (zener*I-zu-zgl).inverse();
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
      double ener = real(zener);

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
      ztmp1 = ener*zunit - zu;
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
      VectorXcd rr;
      rr = ces.eigenvalues();
      int ifail = ces.info();
      if (ifail != 0){
        cout<<"SURFACENEW : ifail = "<<ifail<<endl;
        ifail=1;
        return ifail;
      }

//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//     sort the |evals| into increasing order
//     and check eval(n+1)>eval(n)

//     define zdp
      MatrixXcd zdp(n2,n2);
      zdp.fill(0.);
      zdp.bottomRightCorner(n,n) = ztinv;
      MatrixXcd zfoo1(n2,n2), zfoo2(n2,n2);
      zfoo1 = O;
      zfoo2 = O.inverse();

      double evlir;
      dcomp zevlir, zdevl, zdkde;
      VectorXd evlab(n2);
      evlab.fill(0.);
      for (int ir = 0; ir < n2; ir++){
	evlir = abs(rr(ir));
	zevlir = rr(ir);
	if (evlir > (1+1e-8))
	  evlab(ir) = evlir;
	else if (evlir < (1-1e-8))
	  evlab(ir) = evlir;
        else { // the eigenvalue lies on the unit circle .. check its derivative
	  zdevl = 0.;
	  for (int is = 0; is < n2; is++){
	    for (int it = 0; it < n2; it++){
	      zdevl = zdevl + zfoo2(ir,is)*zdp(is,it)*zfoo1(it,ir);
	    }
	  }
          zdkde=zdevl/(zevlir*zi);
          if (imag(zdkde) > 5.e-5)
	    cout<<"ERROR:dimag(zdkde)=/ 0 "<<imag(zdkde)<<endl;
          evlab(ir)=evlir*exp(-real(zdkde)*1.e-8);
	  }
	}
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      /* do ir=1,n */
      /*   do is=1,n */
      /*     ik=evlab(is+n) */
      /*     ztmp1(ir,is)=dcmplx(vr(ir,ik),vi(ir,ik)) */
      /*     ztmp2(ir,is)=dcmplx(vr(ir+n,ik),vi(ir+n,ik)) */
      /*   enddo */
      /* enddo */
//     calculate L.H. zsurf
      MatrixXcd b = O.topRightCorner(n, n);
      MatrixXcd d = O.bottomRightCorner(n, n);
      zsurfl = b*d.inverse();
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

/* //	** less efficient but ifail doesn't complain */
/*       zgamma = ztmp1*zsinv; */
/*       zp.topLeftCorner(n,n) = zero; */
/*       zp.topRightCorner(n,n) = zsinv; */
/*       zp.bottomLeftCorner(n,n) = -zt; */
/*       zp.bottomRightCorner(n,n) = zgamma; */ 
/* //     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
/* //     diagonalise zp */
/*       ces.compute(zp); */
/*       O = ces.eigenvectors(); */
/*       ifail = ces.info(); */
/*       if (ifail != 0){ */
/*         cout<<"SURFACENEW : ifail = "<<ifail<<endl; */
/*         ifail=1; */
/*         return ifail; */
/*       } */
/* //     calculate R.H. zsurf */
/*       b = O.topRightCorner(n, n); */
/*       d = O.bottomRightCorner(n, n); */
/*       zsurfr = b*d.inverse(); */
/*       MatrixXcd ztmp2; */
/*       ztmp1 = zsurfl; */
/*       ztmp2 = zsurfr; */
/*       adlayer1(ztmp1,zu,zt,zener,n); */
/*       adlayer1(ztmp2,zu,zs,zener,n); */
/*       ztmp1 = ztmp1 - zsurfl; */
/*       ztmp2 = ztmp2 - zsurfr; */
/*       if ((ztmp1.cwiseAbs().maxCoeff() > 5e-5) || (ztmp2.cwiseAbs().maxCoeff() > 5e-5)) */
/*         ifail=1; */
/* //	** */

//	more efficient but not fully functional (though final results don't seem impacted)
//	*UPDATE* since adding the below loop, ifail no longer complains
      /* do ir=1,n2 */
      /*   ik=isort(ir) */
      /*   zevl(ir)=dcmplx(rr(ik),ri(ik)) */
      /* enddo */
      /* do ir=1,n */
      /*   do is=1,n */
      /*     ik=isort(is) */
      /*     ztmp1(ir,is)=dcmplx(vr(ir,ik),vi(ir,ik)) */
      /*     ztmp2(ir,is)=ztmp1(ir,is) */
      /*   enddo */
      /* enddo */
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
      ztmp3 = ztmp1 - zsurfl;
      if (ztmp3.cwiseAbs().maxCoeff() > 5e-5)
        ifail=1;
      zsurfl = ztmp1;
      ztmp3 = ztmp2 - zsurfr;
      while (ztmp3.cwiseAbs().maxCoeff() > 5e-5){
	zsurfr = ztmp2;
        adlayer1(ztmp2,zu,zs,zener,n);
        ztmp3 = ztmp2 - zsurfr;
      }
      zsurfr = ztmp2;

      //This line shifts the bandstructure of Co at the surface if required
      /* adlayer1(zsurfr,zu,zs,zener-0.05,n); */
      /* adlayer1(zsurfr,zu,zs,zener-0.05,n); */

      return ifail;
}

double surfacedecim(MatrixXcd &zu0, MatrixXcd &zt0, dcomp zener, MatrixXcd &zsurfl, MatrixXcd &zsurfr, int n){
//     This program calculates the surface Green's function 
//     using decimation/reduction and Mobius transformation
//     The decimation part reduces the size of U and T is T is singular.
//     This part is described in Sanvito et. al. PRB 73, 085414 (2006)

//     A simpler version is given in (Andrey's)
//     Papers and Notes/Tight Binding/Surface Greens Functions/decimation*.mws

//     -----------------------------------------------------------------
//     The following is a key parameter: may need to change this for 
//     more accurate SGF

      double svmin=5e-7;      // if sv(i)<svmin then sv(i)=0
//     -----------------------------------------------------------------
      MatrixXcd zunit = MatrixXcd::Identity(n,n);
//     -----------------------------------------------------------------
//     Find SVD of zt : zt = Q.S.P     where S(1,1) >= S(2,2,) >= ... >= 0
      MatrixXcd zt(n,n), zu(n,n), zp(n,n), ztmp(n,n);
      zt = zt0;
      /* JacobiSVD<MatrixXcd> svd; */
      JacobiSVD<MatrixXcd,NoQRPreconditioner> svd(zt, ComputeFullV);
      VectorXd sv;
      sv = svd.singularValues();
      zp = svd.matrixV();
      dcomp i;
      i = -1;
      i = sqrt(i);
//     zp can have some very small numbers ... set these to 0
      for (int ii = 0; ii < n; ii++){
	for (int jj = 0; jj < n; jj++){
	  if ((abs(real(zp(ii,jj))) < 1e-40) && (abs(imag(zp(ii,jj))) < 1e-40))
	    zp(ii,jj) = 0.;
	  else if (abs(real(zp(ii,jj))) < 1e-40)
	    zp(ii,jj) = i*imag(zp(ii,jj));
	  else if (abs(imag(zp(ii,jj))) < 1e-40)
	    zp(ii,jj) = real(zp(ii,jj));
	}
      }

//     Now rearrange P and S so that 0 <= S(1,1) <= S(2,2) <= ...
      dcomp ztmp1;
      for (int ii = 0; ii < n/2; ii++){
	for (int jj = 0; jj < n; jj++){
          ztmp1 = zp(ii, jj);
	  zp(ii, jj) = zp(n-ii-1, jj);
	  zp(n-ii-1, jj) = ztmp1;
	}
      }
      sv.reverseInPlace();

//     Now transform all matrices  M -> P.M.P^h
      zu = zu0;
      zt = zt0;
      ztmp = zp*zu;
      zu = ztmp*zp.adjoint();
      ztmp = zp*zt;
      zt = ztmp*zp.adjoint();
//     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//     now set up lst and nlst:   the rows to be decimated
      int icnt=0;
      vector<int> lst;
      int nlst;
      for (int ii = 0; ii < n; ii++){
	if (abs(sv(ii)) < svmin){
	  lst.emplace_back(ii);
	  icnt++;
	}
      }
      int itmp;
      for (int ii = 0; ii < icnt; ii++){
	itmp = lst[ii] + n;
	lst.emplace_back(itmp);
      }
      nlst=2*icnt;

//     -----------------------------------------------------------------
//     Create a matrix to calculate the bulk decimated/reduced on and off 
//     site elements
      int dn = 2*n;
      MatrixXcd zh(dn,dn), zhbulk(dn,dn);
      ztmp = zener*zunit - zu;
      zh.topLeftCorner(n,n) = ztmp;
      zh.bottomRightCorner(n,n) = ztmp;
      zh.topRightCorner(n,n) = -zt;
      zh.bottomLeftCorner(n,n) = -zt.adjoint();

      zhbulk = zh;
      int iii;
      for (int ii = 0; ii<nlst; ii++){
	iii = lst[ii];
        decim(zhbulk,iii,dn);
      }
      remve(zhbulk,lst,dn,nlst);

//     zhbulk is (2.nred x 2.nred)
//     Now calculate the reduced/decimated SGF ie. the very bottom RH element
      int nred=n-nlst/2;
      MatrixXcd zunit2 = MatrixXcd::Identity(nred,nred);
      MatrixXcd zured1(nred,nred), zured2(nred,nred), ztred(nred,nred);
      zured1 = zener*zunit2 - zhbulk.topLeftCorner(nred,nred);
      zured2 = zener*zunit2 - zhbulk.block(nred,nred,nred,nred);
      ztred = -zhbulk.block(0,nred,nred,nred);
//     Create SGFs
      int ifail=0;
      MatrixXcd zglred(nred,nred), zgrred(nred,nred);
      ifail = surfacenew(zured2,ztred,zener,zglred,zgrred,nred);
      if (ifail != 0)// zt has a near zero eigenvalue
	cout<<"eigenvalues ill-conditioned. Consider coding to higher precision"<<endl;

//     =================================================================
//     Left Hand SGF:
//     =================================================================
//     Create the nxn non-reduced SGF by adlayering
//     First create u and t for the last layer
      MatrixXcd zhlast(dn,dn);
      MatrixXcd zulast(n,n), ztlast(n,n);
      zhlast = zh;
      for (int ii=0; ii < nlst/2; ii++){
        iii=lst[ii];
        decim(zhlast,iii,dn);
      }
      zulast = zener*zunit - zhlast.bottomRightCorner(n,n);
//     Bulk out ztlast to be (n x n)
      ztlast = -zhlast.topRightCorner(n,n);

//     Now bulk the SGF to (n x n)
//     (could use non-square matrices instead here!)
      zsurfl.fill(0.);
      zsurfl.block(n-nred,n-nred,nred,nred) = zglred;

//     Now adlayer on top of this
      adlayer1(zsurfl,zulast,ztlast,zener,n);

//     Now untransform this SGF
      zp.adjointInPlace();
      ztmp = zp*zsurfl;
      zsurfl = ztmp*zp.adjoint();

//     =================================================================
//     Right Hand SGF:
//     =================================================================
      MatrixXcd ztreddag(nred,nred);
      ztreddag = ztred.adjoint();
      adlayer1(zgrred,zured1,ztreddag,zener,nred);

//     Now bulk the SGF to (n x n)
      zsurfr.fill(0.);
      zsurfr.block(n-nred,n-nred,nred,nred) = zgrred;

//     Now adlayer on top of this
      MatrixXcd ztdag(n,n);
      ztdag = zt.adjoint();
      adlayer1(zsurfr,zu,ztdag,zener,n);

//     Now untransform this SGF
      ztmp = zp*zsurfr;
      zsurfr = ztmp*zp.adjoint();

//     =================================================================
//     =================================================================

//     Check these are the correct SGF
      MatrixXcd zsurf2(n,n);
      zsurf2 = zsurfl;
      adlayer1(zsurf2,zu0,zt0,zener,n);
      double xmaxerrl=0.;
      double errl;
      for (int ii = 0; ii < n; ii++){
	for (int jj = 0; jj < n; jj++){
          errl = abs(zsurfl(ii,jj) - zsurf2(ii,jj));
	  xmaxerrl = max(xmaxerrl, errl);
	}
      }

      MatrixXcd zt0dag(n,n);
      zsurf2 = zsurfr;
      zt0dag = zt0.adjoint();
      adlayer1(zsurf2,zu0,zt0dag,zener,n);
      double xmaxerrr=0.;
      double errr;
      for (int ii = 0; ii < n; ii++){
	for (int jj = 0; jj < n; jj++){
          errr = abs(zsurfr(ii,jj) - zsurf2(ii,jj));
	  xmaxerrr = max(xmaxerrr, errr);
	}
      }

      double errmax=max(xmaxerrl,xmaxerrr);

      return errmax;
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
      /* cout<<zh.real()<<endl<<endl; */
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
      /* cout<<zt.real()<<endl<<endl; */

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
	int nsub, int nsubat, int nlay, int nmat, int nxfold, vV3d &xfold, int nspin, vector<int> &imapl, 
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
	//this routine supercedes surfacenew below
	double errmax;
        errmax = surfacedecim(zuat,ztat,zener,zglat,zgrat,natom);
        if (errmax > 1e-3){   // zt has a near zero eigenvalue
          cout<<"ERROR surfacedecim:  errmax = "<<errmax<<endl;
          ifail=1;
	}
        /* ifail = surfacenew(zuat,ztat,zener,zglat,zgrat,natom); */
        /* if (ifail != 0)// zt has a near zero eigenvalue */
	  /* cout<<"eigenvalues ill-conditioned. Consider coding to higher precision"<<endl; */

//       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//       now load elements into zgl and zgr
        for (int isub=0; isub<nsub; isub++){           // sublattice elements
          for (int jsub=0; jsub<nsub; jsub++){
            ds=vsub[ilay][isub]-vsub[ilay][jsub] - (vsubat[ilay][imap[isub]-1]-vsubat[ilay][imap[jsub]-1]);
            zarg1=i*xkat.dot(ds);

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
	int nsub, int nsubat, int nxfold, vV3d &xfold, int nmat, int mlay, int nins, int nlay, Args&&... params){

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
      ifail = green(zener,+1,st,zglu,zgru,nsub,nsubat,nlay,nmat,nxfold,xfold,forward<Args>(params)...);   // LH UP
      if (ifail != 0)
	return ifail;
      ifail = green(zener,-1,st,zgld,zgrd,nsub,nsubat,nlay,nmat,nxfold,xfold,forward<Args>(params)...);   // LH DOWN
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
