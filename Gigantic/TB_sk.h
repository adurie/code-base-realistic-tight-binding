#ifndef TB_H
#define TB_H

#include <cmath>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

Matrix<complex<double>, 9, 9> eint1(double, double, double, double,
		     double, double, double, double, double,
		     double, double, double, double);

Matrix<complex<double>, 9, 9> TB(int ind, int dd, int nn, int nspin, const Vector3d &pos, VectorXd &NN, VectorXd &NNN){
//  lattice constant a set to 1 as is irrelevant, disappears after integration
//  ind = atom type,    0=Fe down, 1=Fe up,  2=Ag
//  dd = type of matrix, 0 = U, 1 = T
//  nn = 0 or 1 depending on 1st or 2nd NN
//  nspin is the number of orbitals to be modelled (1-9)
//  pos is the position vector of each atomic site

//    *****************************************************************;
//    THIS ROUTINE IS ATOM DEPENDENT :-;
//    it is written for Cu and Co-up and Co-down, with s,p,d bands.;

//    two centre integrals; slater-koster parameters for spin polarized;
//    bcc Iron was obtained starting from the paramagnetic parametrization;
//    of R H Victora "Magnetic and Electronic Properties of Transition;
//    Metals and Overlayers". Self consistency is performed shifting;
//    the centre of the d bands only (Udd=1.eV Usp=0);
//    bc//Cr bulk spin bands:;
//    Also from R H Victora's paper (Udd=Usp=0 eV);


//    -----------------------------------------------------------------;
//    The first index in the tight binding parameter arrays refers to;
//    1: Fe down;
//    2: Fe up;
//    3: Ag;

	Vector3d s0;
	Vector3d p0;
	Vector3d d0t;
	Vector3d d0e;

	Matrix<double, 3, 2> sss, sps, pps, ppp, sds, pds, pdp, dds, ddp, ddd;
//    -----------------------------------------------------------------;
      const double delta = 1.545054e-01;
      const double dels2 = 0.5*delta;
      const double cshift = .4635 + .02490135;
      /* const double cshift = .485 + .02490135; */
//
//     Fe down:
//      
//     on site
//      
      s0(0) =  0.514+cshift;
      p0(0) =  1.118+cshift;
      d0t(0) = -0.089 + dels2+cshift;
      d0e(0) = -0.068 + dels2+cshift;
      
//     Fe up:      
      
//     on site
      
      s0(1) =  0.514+cshift;
      p0(1) =  1.118+cshift;
      d0t(1) = -0.089 - dels2+cshift;
      d0e(1) = -0.068 - dels2+cshift;

      //Fe up and down nn
//
//     first n.n.
//     
      for (int ii = 0; ii < 2; ii++){
        sss(ii,0) = NN(0); 
        sps(ii,0) = NN(1);
        pps(ii,0) = NN(2);
        ppp(ii,0) = NN(3);
        sds(ii,0) = NN(4);
        pds(ii,0) = NN(5);
        pdp(ii,0) = NN(6);
        dds(ii,0) = NN(7);
        ddp(ii,0) = NN(8);
        ddd(ii,0) = NN(9);
      
//      
//     second n.n.
//
        sss(ii,1) = NNN(0);
        sps(ii,1) = NNN(1);
        pps(ii,1) = NNN(2);
        ppp(ii,1) = NNN(3);
        sds(ii,1) = NNN(4);
        pds(ii,1) = NNN(5);
        pdp(ii,1) = NNN(6);
        dds(ii,1) = NNN(7);
        ddp(ii,1) = NNN(8);
        ddd(ii,1) = NNN(9);
      }

//     -----------------------------------------------------------------
//     Ag:

      s0(2) =  0.68297;
      p0(2) =  1.13432;
      d0t(2) =  0.12249;
      d0e(2) =  0.12006;

//     first n.n.

      sss(2,0) = -0.06581;
      sps(2,0) =  0.09781;
      pps(2,0) =  0.15752;
      ppp(2,0) =  0.00649;
      sds(2,0) = -0.03110;
      pds(2,0) = -0.03905;
      pdp(2,0) =  0.01519;
      dds(2,0) = -0.03151;
      ddp(2,0) =  0.01757;
      ddd(2,0) = -0.00336;

//     second n.n.

      sss(2,1) =  0.00143;
      sps(2,1) =  0.00545;
      pps(2,1) =  0.03971;
      ppp(2,1) =  0.00434;
      sds(2,1) = -0.00462;
      pds(2,1) = -0.00065;
      pdp(2,1) =  0.00172;
      dds(2,1) = -0.00282;
      ddp(2,1) =  0.00171;
      ddd(2,1) = -0.00038;

      Matrix<complex<double>, 9, 9> MAT = Matrix<complex<double>, 9, 9>::Zero();

//  If I want to construct  <0|H|0> then the 9x9 matrix U is diagonal with elements:
        if(dd == 0)
	{
          MAT(0,0)=s0(ind);
          MAT(1,1)=p0(ind);
          MAT(2,2)=p0(ind);
          MAT(3,3)=p0(ind);
          for(int ir = 4; ir < 7; ir++)
            MAT(ir,ir)=d0t(ind);
          for(int ir = 7; ir < 9; ir++)
            MAT(ir,ir)=d0e(ind);
	}
        else
	{

//  If I want to construct <0|H|d> then
//  d(k)   !!! it has position vector d
//  c(k) = d.(i, j or k)/|d|   !!! direction cosines
	  double x, y, z;
	  Vector3d X, Y, Z;
	  X << 1, 0, 0;
	  Y << 0, 1, 0;
	  Z << 0, 0, 1;
	  x = pos.dot(X)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	  y = pos.dot(Y)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 
	  z = pos.dot(Z)/sqrt(pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2)); 

          double g1,g2,g3,g4,g5,g6,g7,g8,g9,g10;
//   the 9 x 9 matrix T is given by;
          g1=sss(ind,nn);
          g2=sps(ind,nn);
          g3=pps(ind,nn);
          g4=ppp(ind,nn);
          g5=sds(ind,nn);
          g6=pds(ind,nn);
          g7=pdp(ind,nn);
          g8=dds(ind,nn);
          g9=ddp(ind,nn);
          g10=ddd(ind,nn);
          MAT = eint1(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,x,y,z);
	}

	/* if (nspin > 8){ */
	/* 	Matrix<complex<double>, 9, 9> return_mat; */
	/* 	return_mat = MAT; */
	/* } */
	/* else{ */
	/* 	Matrix<complex<double>, nspin, nspin> return_mat; */
	/* 	return_mat = MAT.topLeftCorner(nspin, nspin); */
	/* } */
	/* return return_mat; */

	return MAT;

}

//    *****************************************************************;

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
#endif
