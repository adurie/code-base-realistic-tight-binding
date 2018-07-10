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

      const double cshift = .4635 - .9252;
      //measured fermi level for Fe 0.841 for spin down, 0.7006 for spin up
      /* const double cshift = .485 - .7708; */
      /* const double cshift = 0.; */

      double s, p, d1, d2, sss1, sss2, pps1, pps2, ppp1, ppp2, 
	     dds1, dds2, ddp1, ddp2, ddd1, ddd2, sps1, sps2, sds1, sds2, pds1, pds2, pdp1, pdp2;

      s = 1.33239;
      p = 1.94576;
      d1 = 0.846975;
      d2 = 0.79515;

      s0(0,0) = s + cshift; // on-site
      p0(0,0) = p + cshift;
      d0t(0,0) = d1 + cshift;
      d0e(0,0) = d2 + cshift;

//     first n.n.

      double sss3, sps3, pps3, ppp3, sds3, pds3, pdp3, dds3, ddp3, ddd3;

      sss1 = -0.117924; sps1 = 0.169252; pps1 = 0.268776; ppp1 = -0.0187581; sds1 = -0.0757635;
      pds1 = -0.110995; pdp1 = 0.0287141; dds1 = -0.0447073; ddp1 = 0.0197826; ddd1 = 5.90525e-05;
      sss2 = -0.0220203; sps2 = 0.0614591; pps2 = 0.1639; ppp2 = 0.0290812; sds2 = -0.0326461;
      pds2 = -0.0571585; pdp2 = -0.0117408; dds2 = -0.0253075; ddp2 = -0.00318763; ddd2 = 0.0026379;
      sss3 = 0.017236; sps3 = -0.0252936; pps3 = -0.0503531; ppp3 = 0.0160445; sds3 = 0.0088867; 
      pds3 = 0.0103232; pdp3 = -0.00234456; dds3 = 0.00533451; ddp3 = -0.00119393; ddd3 = -0.000529859;

      sss[0](0,0) = sss1;  //  same atom hopping
      sps[0](0,0) = sps1;
      pps[0](0,0) = pps1;
      ppp[0](0,0) = ppp1; 
      sds[0](0,0) = sds1;
      pds[0](0,0) = pds1;
      pdp[0](0,0) = pdp1;
      dds[0](0,0) = dds1;
      ddp[0](0,0) = ddp1;
      ddd[0](0,0) = ddd1;

//     second n.n.

      sss[0](0,1) = sss2;
      sps[0](0,1) = sps2;
      pps[0](0,1) = pps2;
      ppp[0](0,1) = ppp2; 
      sds[0](0,1) = sds2;
      pds[0](0,1) = pds2;
      pdp[0](0,1) = pdp2;
      dds[0](0,1) = dds2;
      ddp[0](0,1) = ddp2;
      ddd[0](0,1) = ddd2;
      
//     third n.n.

      sss[0](0,2) = sss3;
      sps[0](0,2) = sps3;
      pps[0](0,2) = pps3;
      ppp[0](0,2) = ppp3; 
      sds[0](0,2) = sds3;
      pds[0](0,2) = pds3;
      pdp[0](0,2) = pdp3;
      dds[0](0,2) = dds3;
      ddp[0](0,2) = ddp3;
      ddd[0](0,2) = ddd3;

//     -----------------------------------------------------------------
//      Fe up

      s = 1.37057;
      p = 1.97431;
      d1 = 1.00461;
      d2 = 0.996866;

      s0(0,1) = s + cshift; // on-site
      p0(0,1) = p + cshift;
      d0t(0,1) = d1 + cshift;
      d0e(0,1) = d2 + cshift;

//     first n.n.

      sss1 = -0.121541; sps1 = 0.172875; pps1 = 0.273266; ppp1 = -0.0182001; sds1 = -0.0852609;
      pds1 = -0.124371; pdp1 = 0.0233553; dds1 = -0.049139; ddp1 = 0.0262377; ddd1 = -0.00143019;
      sss2 = -0.0239523; sps2 = 0.0608069; pps2 = 0.170197; ppp2 = 0.0314321; sds2 = -0.0365085;
      pds2 = -0.0703064; pdp2 = -0.0146721; dds2 = -0.032293; ddp2 = -0.00262184; ddd2 = 0.00332746;
      sss3 = 0.0190502; sps3 = -0.0285496; pps3 = -0.0500113; ppp3 = 0.0167854; sds3 = 0.00457944; 
      pds3 = 0.00712668; pdp3 = -0.00213332; dds3 = 0.00432888; ddp3 = -0.00220879; ddd3 = 0.000162835;

      sss[1](0,0) = sss1;  //  same atom hopping
      sps[1](0,0) = sps1;
      pps[1](0,0) = pps1;
      ppp[1](0,0) = ppp1; 
      sds[1](0,0) = sds1;
      pds[1](0,0) = pds1;
      pdp[1](0,0) = pdp1;
      dds[1](0,0) = dds1;
      ddp[1](0,0) = ddp1;
      ddd[1](0,0) = ddd1;

//     second n.n.

      sss[1](0,1) = sss2;
      sps[1](0,1) = sps2;
      pps[1](0,1) = pps2;
      ppp[1](0,1) = ppp2; 
      sds[1](0,1) = sds2;
      pds[1](0,1) = pds2;
      pdp[1](0,1) = pdp2;
      dds[1](0,1) = dds2;
      ddp[1](0,1) = ddp2;
      ddd[1](0,1) = ddd2;
      
//     third n.n.

      sss[1](0,2) = sss3;
      sps[1](0,2) = sps3;
      pps[1](0,2) = pps3;
      ppp[1](0,2) = ppp3; 
      sds[1](0,2) = sds3;
      pds[1](0,2) = pds3;
      pdp[1](0,2) = pdp3;
      dds[1](0,2) = dds3;
      ddp[1](0,2) = ddp3;
      ddd[1](0,2) = ddd3;
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
