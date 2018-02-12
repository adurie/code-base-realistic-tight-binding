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

void param(int numat, int numnn, m2d &s0, m2d &p0, m2d &d0t, m2d &d0e,
	       	vm2d &sssint, vm2d &spsint, vm2d &ppsint, vm2d &pppint, vm2d &sdsint,
		vm2d &pdsint, vm2d &pdpint, vm2d &ddsint, vm2d &ddpint, vm2d &dddint){
//     THIS ROUTINE IS ATOM DEPENDENT :-
//     -----------------------------------------------------------------
//     The first index in the tight binding parameter arrays refers to
//         1: Bulk Co
//         2: Bulk Cu

//         Interatomic hoppings at end
//     -----------------------------------------------------------------
//     Co up:

      m2d sss, sps, pps, ppp, sds, pds, pdp, dds, ddp, ddd;
      const double cshift = .575530 - .715751;
      const double delta = 1.113608939931278e-1;
      const double vex = delta*0.5;

      for (int isp=0; isp<2; isp++){    // isp=0 minority , =1 majority
        s0(0,isp) =  1.12946 + cshift; // on-site
        p0(0,isp) =  1.75262 + cshift;
        d0t(0,isp) =  0.5*(0.60547 + 0.60445) + cshift - (2*isp-1)*vex;
        d0e(0,isp) =  0.5*(0.60547 + 0.60445) + cshift - (2*isp-1)*vex;
      }

//     first n.n.

      sss(0,0) = -0.09043;   //  same atom hopping
      sps(0,0) =  0.13649;
      pps(0,0) =  0.23748;
      ppp(0,0) = -0.00142;
      sds(0,0) = -0.03806;
      pds(0,0) = -0.04069;
      pdp(0,0) =  0.02797;
      dds(0,0) = -0.04213;
      ddp(0,0) =  0.02976;
      ddd(0,0) = -0.00684;

//     second n.n.

      sss(0,1) = -0.00337;
      sps(0,1) =  0.00135;
      pps(0,1) =  0.02849;
      ppp(0,1) =  0.01099;
      sds(0,1) = -0.01119;
      pds(0,1) = -0.01061;
      pdp(0,1) =  0.01134;
      dds(0,1) = -0.00759;
      ddp(0,1) =  0.00495;
      ddd(0,1) = -0.00016;

//     -----------------------------------------------------------------
//     Cu:

      for (int isp=0; isp<2; isp++){
        s0(1,isp) =  0.79466;
        p0(1,isp) =  1.35351;
        d0t(1,isp) =  0.5*(0.37307 + 0.37180);
        d0e(1,isp) =  0.5*(0.37307 + 0.37180);
      }

//     first n.n.

      sss(1,0) = -0.07518;
      sps(1,0) =  0.11571;
      pps(1,0) =  0.19669;
      ppp(1,0) =  0.01940;
      sds(1,0) = -0.03107;
      pds(1,0) = -0.03289;
      pdp(1,0) =  0.01753;
      dds(1,0) = -0.02566;
      ddp(1,0) =  0.01800;
      ddd(1,0) = -0.00408;

//     second n.n.

      sss(1,1) = -0.00092;
      sps(1,1) =  0.01221;
      pps(1,1) =  0.05389;
      ppp(1,1) =  0.00846;
      sds(1,1) = -0.00852;
      pds(1,1) = -0.00536;
      pdp(1,1) =  0.00321;
      dds(1,1) = -0.00451;
      ddp(1,1) =  0.00241;
      ddd(1,1) = -0.00029;
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
