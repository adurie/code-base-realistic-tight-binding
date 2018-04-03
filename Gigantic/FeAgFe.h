#ifndef COCUCO_H
#define COCUCO_H

#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>

using namespace std;
using namespace Eigen;
typedef Matrix<double, 2, 3> m23;
typedef vector<m23, aligned_allocator<m23>> m2d;
typedef vector<vector<m23, aligned_allocator<m23>>> vm2d;

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
//     THIS ROUTINE IS ATOM DEPENDENT :-
//     -----------------------------------------------------------------
//     The first index in the tight binding parameter arrays refers to
//         1: Bulk Fe 
//         2: Bulk Ag 

//         Interatomic hoppings at end
//     -----------------------------------------------------------------
//     initialise m2d
      m2d sss, sps, pps, ppp, sds, pds, pdp, dds, ddp, ddd;
      Matrix<double, 2, 3> zero = Matrix<double, 2, 3>::Zero();
      for (int isp=0; isp<2; isp++){
       sss.emplace_back(zero); 
       sps.emplace_back(zero); 
       pps.emplace_back(zero); 
       ppp.emplace_back(zero); 
       sds.emplace_back(zero); 
       pds.emplace_back(zero); 
       pdp.emplace_back(zero); 
       dds.emplace_back(zero); 
       ddp.emplace_back(zero); 
       ddd.emplace_back(zero); 
      }
      for (int isp=0; isp<2; isp++){
       sssint.emplace_back(sss); 
       spsint.emplace_back(sss); 
       ppsint.emplace_back(sss); 
       pppint.emplace_back(sss); 
       sdsint.emplace_back(sss); 
       pdsint.emplace_back(sss); 
       pdpint.emplace_back(sss); 
       ddsint.emplace_back(sss); 
       ddpint.emplace_back(sss); 
       dddint.emplace_back(sss); 
      }

//     Fe down:

      /* const double cshift = .4635 - .7505; */
      //measured fermi level for Fe 0.841 for spin down, 0.7006 for spin up
      const double cshift = .485 - .7708;
      /* const double cshift = 0.; */

      s0(0,0) =  1.14481 + cshift; // on-site
      p0(0,0) =  1.80769 + cshift;
      d0t(0,0) =  0.5*(0.78456 + 0.75661) + cshift;
      d0e(0,0) =  0.5*(0.78456 + 0.75661) + cshift;

//     first n.n.

      sss[0](0,0) = -0.13243;   //  same atom hopping
      sps[0](0,0) =  0.17278;
      pps[0](0,0) =  0.25911;
      ppp[0](0,0) =  0.02653;
      sds[0](0,0) = -0.07145;
      pds[0](0,0) = -0.09702;
      pdp[0](0,0) =  0.02129;
      dds[0](0,0) = -0.05266;
      ddp[0](0,0) =  0.03276;
      ddd[0](0,0) = -0.00286;

//     second n.n.

      sss[0](0,1) = -0.03003;
      sps[0](0,1) =  0.07159;
      pps[0](0,1) =  0.18256;
      ppp[0](0,1) =  0.03703;
      sds[0](0,1) = -0.04075;
      pds[0](0,1) = -0.06522;
      pdp[0](0,1) = -0.00467;
      dds[0](0,1) = -0.03396;
      ddp[0](0,1) =  0.00581;
      ddd[0](0,1) =  0.00114;
      
//     third n.n.

      sss[0](0,2) =  0.01589;
      sps[0](0,2) = -0.02306;
      pps[0](0,2) = -0.04253;
      ppp[0](0,2) =  0.01538;
      sds[0](0,2) =  0.00016;
      pds[0](0,2) =  0.00222;
      pdp[0](0,2) = -0.00351;
      dds[0](0,2) =  0.00233;
      ddp[0](0,2) =  0.00013;
      ddd[0](0,2) = -0.00060;

//     -----------------------------------------------------------------
//      Fe up

      s0(0,1) =  1.13516 + cshift; // on-site
      p0(0,1) =  1.81739 + cshift;
      d0t(0,1) =  0.5*(0.64840 + 0.62960) + cshift;
      d0e(0,1) =  0.5*(0.64840 + 0.62960) + cshift;

//     first n.n.

      sss[1](0,0) = -0.12950;   //  same atom hopping
      sps[1](0,0) =  0.17363;
      pps[1](0,0) =  0.25741;
      ppp[1](0,0) =  0.02422;
      sds[1](0,0) = -0.06115;
      pds[1](0,0) = -0.08485;
      pdp[1](0,0) =  0.01778;
      dds[1](0,0) = -0.04541;
      ddp[1](0,0) =  0.02714;
      ddd[1](0,0) = -0.00260;

//     second n.n.

      sss[1](0,1) = -0.02915;
      sps[1](0,1) =  0.06571;
      pps[1](0,1) =  0.16827;
      ppp[1](0,1) =  0.04112;
      sds[1](0,1) = -0.03560;
      pds[1](0,1) = -0.05473;
      pdp[1](0,1) = -0.00280;
      dds[1](0,1) = -0.02713;
      ddp[1](0,1) =  0.00589;
      ddd[1](0,1) =  0.00060;
      
//     third n.n.

      sss[1](0,2) =  0.01595;
      sps[1](0,2) = -0.02477;
      pps[1](0,2) = -0.04985;
      ppp[1](0,2) =  0.01796;
      sds[1](0,2) = -0.00073;
      pds[1](0,2) = -0.00082;
      pdp[1](0,2) = -0.00241;
      dds[1](0,2) =  0.00112;
      ddp[1](0,2) =  0.00034;
      ddd[1](0,2) = -0.00056;
//     -----------------------------------------------------------------
//     Ag:

      for (int isp=0; isp<2; isp++){
        s0(1,isp) =  0.68297;
        p0(1,isp) =  1.13432;
        d0t(1,isp) =  0.5*(0.12249 + 0.12006);
        d0e(1,isp) =  0.5*(0.12249 + 0.12006);

//     first n.n.

        sss[isp](1,0) = -0.06581;
        sps[isp](1,0) =  0.09781;
        pps[isp](1,0) =  0.15752;
        ppp[isp](1,0) =  0.00649;
        sds[isp](1,0) = -0.03110;
        pds[isp](1,0) = -0.03905;
        pdp[isp](1,0) =  0.01519;
        dds[isp](1,0) = -0.03151;
        ddp[isp](1,0) =  0.01757;
        ddd[isp](1,0) = -0.00336;

//     second n.n.

        sss[isp](1,1) =  0.00143;
        sps[isp](1,1) =  0.00545;
        pps[isp](1,1) =  0.03971;
        ppp[isp](1,1) =  0.00434;
        sds[isp](1,1) = -0.00462;
        pds[isp](1,1) = -0.00065;
        pdp[isp](1,1) =  0.00172;
        dds[isp](1,1) = -0.00282;
        ddp[isp](1,1) =  0.00171;
        ddd[isp](1,1) = -0.00038;

//     third n.n.

        sss[isp](1,2) =  0.;
        sps[isp](1,2) =  0.;
        pps[isp](1,2) =  0.;
        ppp[isp](1,2) =  0.;
        sds[isp](1,2) =  0.;
        pds[isp](1,2) =  0.;
        pdp[isp](1,2) =  0.;
        dds[isp](1,2) =  0.;
        ddp[isp](1,2) =  0.;
        ddd[isp](1,2) =  0.;
      }
//     =================================================================
//     Now evaluate the inter-atomic hoppings 
      for (int isp=0; isp<2; isp++){
        for (int iat=0; iat<numat; iat++){
          for (int jat=0; jat<numat; jat++){
            for (int inn=0; inn<numnn; inn++){
              sssint[isp][iat](jat,inn)=gmean(sss[isp](iat,inn),sss[isp](jat,inn)); //different atom hopping
              ppsint[isp][iat](jat,inn)=gmean(pps[isp](iat,inn),pps[isp](jat,inn));
              pppint[isp][iat](jat,inn)=gmean(ppp[isp](iat,inn),ppp[isp](jat,inn));
              ddsint[isp][iat](jat,inn)=gmean(dds[isp](iat,inn),dds[isp](jat,inn));
              ddpint[isp][iat](jat,inn)=gmean(ddp[isp](iat,inn),ddp[isp](jat,inn));
              dddint[isp][iat](jat,inn)=gmean(ddd[isp](iat,inn),ddd[isp](jat,inn));
              spsint[isp][iat](jat,inn)=gmean(sps[isp](iat,inn),sps[isp](jat,inn));
              sdsint[isp][iat](jat,inn)=gmean(sds[isp](iat,inn),sds[isp](jat,inn));
              pdsint[isp][iat](jat,inn)=gmean(pds[isp](iat,inn),pds[isp](jat,inn));
              pdpint[isp][iat](jat,inn)=gmean(pdp[isp](iat,inn),pdp[isp](jat,inn));
	      //turn gmean off
	      /* if ((iat == 0) && (jat == 0)){ */
                /* sssint[isp][iat](jat,inn)=sss[isp](iat,inn); */
                /* ppsint[isp][iat](jat,inn)=pps[isp](iat,inn); */
                /* pppint[isp][iat](jat,inn)=ppp[isp](iat,inn); */
                /* ddsint[isp][iat](jat,inn)=dds[isp](iat,inn); */
                /* ddpint[isp][iat](jat,inn)=ddp[isp](iat,inn); */
                /* dddint[isp][iat](jat,inn)=ddd[isp](iat,inn); */
                /* spsint[isp][iat](jat,inn)=sps[isp](iat,inn); */
                /* sdsint[isp][iat](jat,inn)=sds[isp](iat,inn); */
                /* pdsint[isp][iat](jat,inn)=pds[isp](iat,inn); */
                /* pdpint[isp][iat](jat,inn)=pdp[isp](iat,inn); */
	      /* } */
	      /* else{ */
                /* sssint[isp][iat](jat,inn)=sss[isp](1,inn); */
                /* ppsint[isp][iat](jat,inn)=pps[isp](1,inn); */
                /* pppint[isp][iat](jat,inn)=ppp[isp](1,inn); */
                /* ddsint[isp][iat](jat,inn)=dds[isp](1,inn); */
                /* ddpint[isp][iat](jat,inn)=ddp[isp](1,inn); */
                /* dddint[isp][iat](jat,inn)=ddd[isp](1,inn); */
                /* spsint[isp][iat](jat,inn)=sps[isp](1,inn); */
                /* sdsint[isp][iat](jat,inn)=sds[isp](1,inn); */
                /* pdsint[isp][iat](jat,inn)=pds[isp](1,inn); */
                /* pdpint[isp][iat](jat,inn)=pdp[isp](1,inn); */
	      /* } */
	    }
	  }
	}
      }

      return;
}
#endif
