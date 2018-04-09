#ifndef COCUCO_H
#define COCUCO_H

#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>

using namespace std;
using namespace Eigen;
typedef Matrix2d m2d;
typedef vector<Matrix2d, aligned_allocator<Matrix2d>> vm2d;

double gmean(double x, double y){
	double gmean;
	if (x == y)
		gmean = x;
	else if (x*y > 0){
		if (x < 0)
			gmean=-sqrt(x*y);
		if (x > 0)
		gmean = sqrt(x*y);
      	}
	else
       		gmean=(x+y)/2.;
 	return gmean;
}

void param(int numat, int numnn, Matrix2d &s0, Matrix2d &p0, Matrix2d &d0t, Matrix2d &d0e,
	       	vm2d &sssint, vm2d &spsint, vm2d &ppsint, vm2d &pppint, vm2d &sdsint,
		vm2d &pdsint, vm2d &pdpint, vm2d &ddsint, vm2d &ddpint, vm2d &dddint){

      m2d sss, sps, pps, ppp, sds, pds, pdp, dds, ddp, ddd;
//     With these params the Fermi energy is (Alex you need to check this)
//      ef=-0.02490135;

      const double delta = 1.545054e-01;
      const double dels2 = 0.5*delta;
      const double cshift = .4635 + .02490135;
      /* const double cshift = .485 + .02490135; */
//
//     Fe down:
//      
//     on site
//      
      s0(0,0) =  0.514+cshift;
      p0(0,0) =  1.118+cshift;
      d0t(0,0) = -0.089 + dels2+cshift;
      d0e(0,0) = -0.068 + dels2+cshift;
      
//     Fe up:      
      
//     on site
      
      s0(0,1) =  0.514+cshift;
      p0(0,1) =  1.118+cshift;
      d0t(0,1) = -0.089 - dels2+cshift;
      d0e(0,1) = -0.068 - dels2+cshift;

      //Fe up and down nn
//
//     first n.n.
//      
      sss(0,0) = -0.10605;
      sps(0,0) = -0.19326;
      pps(0,0) =  0.26932;
      ppp(0,0) = -0.00672;
      sds(0,0) = -0.08018;
      pds(0,0) =  0.09915;
      pdp(0,0) = -0.02467;
      dds(0,0) = -0.05669;
      ddp(0,0) =  0.03738;
      ddd(0,0) = -0.00669;
//      
//     second n.n.
//
      sss(0,1) = -0.05276;
      sps(0,1) = -0.06725;
      pps(0,1) =  0.17199;
      ppp(0,1) =  0.03242;
      sds(0,1) = -0.03035;
      pds(0,1) =  0.05202;
      pdp(0,1) = -0.00092;
      dds(0,1) = -0.03367;
      ddp(0,1) =  0.01000;
      ddd(0,1) = -0.00075;
//     -----------------------------------------------------------------
//     Ag:

      for (int isp=0; isp<2; isp++){
        s0(1,isp) =  0.68297;
        p0(1,isp) =  1.13432;
        d0t(1,isp) =  0.5*(0.12249 + 0.12006);
        d0e(1,isp) =  0.5*(0.12249 + 0.12006);
      }

//     first n.n.

      sss(1,0) = -0.06581;
      sps(1,0) =  0.09781;
      pps(1,0) =  0.15752;
      ppp(1,0) =  0.00649;
      sds(1,0) = -0.03110;
      pds(1,0) = -0.03905;
      pdp(1,0) =  0.01519;
      dds(1,0) = -0.03151;
      ddp(1,0) =  0.01757;
      ddd(1,0) = -0.00336;

//     second n.n.

      sss(1,1) =  0.00143;
      sps(1,1) =  0.00545;
      pps(1,1) =  0.03971;
      ppp(1,1) =  0.00434;
      sds(1,1) = -0.00462;
      pds(1,1) = -0.00065;
      pdp(1,1) =  0.00172;
      dds(1,1) = -0.00282;
      ddp(1,1) =  0.00171;
      ddd(1,1) = -0.00038;

//     =================================================================
//     Now evaluate the inter-atomic hoppings 
      for (int iat=0; iat<numat; iat++){
        for (int jat=0; jat<numat; jat++){
          for (int inn=0; inn<numnn; inn++){
            sssint[iat](jat,inn)=gmean(sss(iat,inn),sss(jat,inn)); //different atom hopping
            ppsint[iat](jat,inn)=gmean(pps(iat,inn),pps(jat,inn));
            pppint[iat](jat,inn)=gmean(ppp(iat,inn),ppp(jat,inn));
            ddsint[iat](jat,inn)=gmean(dds(iat,inn),dds(jat,inn));
            ddpint[iat](jat,inn)=gmean(ddp(iat,inn),ddp(jat,inn));
            dddint[iat](jat,inn)=gmean(ddd(iat,inn),ddd(jat,inn));
            spsint[iat](jat,inn)=gmean(sps(iat,inn),sps(jat,inn));
            sdsint[iat](jat,inn)=gmean(sds(iat,inn),sds(jat,inn));
            pdsint[iat](jat,inn)=gmean(pds(iat,inn),pds(jat,inn));
            pdpint[iat](jat,inn)=gmean(pdp(iat,inn),pdp(jat,inn));
	  }
	}
      }

      return;
}
#endif
