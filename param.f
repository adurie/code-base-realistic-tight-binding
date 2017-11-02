      subroutine param
c     THIS ROUTINE IS ATOM DEPENDENT :-
c     it is written for Cu and Co-up and Co-down, with s,p,d bands.
      implicit double precision (a-h,o-z)
      common/par1/ s0(5),p0(5),d0e(5),d0t(5)
      common/par2/ sss(5,2),sps(5,2),pps(5,2),ppp(5,2),sds(5,2)
      common/par3/ pds(5,2),pdp(5,2),dds(5,2),ddp(5,2),ddd(5,2)
c
c     two centre integrals; slater-koster parameters for spin polarized
c     bcc Iron was obtained starting from the paramagnetic parametrization
c     of R H Victora "Magnetic and Electronic Properties of Transition
c     Metals and Overlayers". Self consistency is performed shifting
c     the centre of the d bands only (Udd=1.eV Usp=0)
c     bcc Cr bulk spin bands:
c     Also from R H Victora's paper (Udd=Usp=0 eV)
c
c
c     -----------------------------------------------------------------
c     The first index in the tight binding parameter arrays refers to
c     1: Co up
c     2: Co down
c     3: Cu
c
c     -----------------------------------------------------------------
c     Co up:
c
      cshift = .575530d0 - .715751d0
      delta = 1.113608939931278d-001
      vex = delta*0.5d0
      s0(1) =  1.12946d0 + cshift
      p0(1) =  1.75262d0 + cshift
      d0t(1) =  0.5d0*(0.60547d0 + 0.60445d0) + cshift - vex
      d0e(1) =  0.5d0*(0.60547d0 + 0.60445d0) + cshift - vex
c
c     first n.n.
c
      sss(1,1) = -0.09043d0
      sps(1,1) =  0.13649d0
      pps(1,1) =  0.23748d0
      ppp(1,1) = -0.00142d0
      sds(1,1) = -0.03806d0
      pds(1,1) = -0.04069d0
      pdp(1,1) =  0.02797d0
      dds(1,1) = -0.04213d0
      ddp(1,1) =  0.02976d0
      ddd(1,1) = -0.00684d0
c
c     second n.n.
c
      sss(1,2) = -0.00337d0
      sps(1,2) =  0.00135d0
      pps(1,2) =  0.02849d0
      ppp(1,2) =  0.01099d0
      sds(1,2) = -0.01119d0
      pds(1,2) = -0.01061d0
      pdp(1,2) =  0.01134d0
      dds(1,2) = -0.00759d0
      ddp(1,2) =  0.00495d0
      ddd(1,2) = -0.00016d0
c
c     -----------------------------------------------------------------
c     Co down:
c
      cshift = .575530d0 - .715751d0
      delta = 1.113608939931278d-001
      vex = delta*0.5d0
      s0(2) =  1.12946d0 + cshift
      p0(2) =  1.75262d0 + cshift
      d0t(2) =  0.5d0*(0.60547d0 + 0.60445d0) + cshift + vex
      d0e(2) =  0.5d0*(0.60547d0 + 0.60445d0) + cshift + vex
c
c     first n.n.
c
      sss(2,1) = -0.09043d0
      sps(2,1) =  0.13649d0
      pps(2,1) =  0.23748d0
      ppp(2,1) = -0.00142d0
      sds(2,1) = -0.03806d0
      pds(2,1) = -0.04069d0
      pdp(2,1) =  0.02797d0
      dds(2,1) = -0.04213d0
      ddp(2,1) =  0.02976d0
      ddd(2,1) = -0.00684d0
c
c     second n.n.
c
      sss(2,2) = -0.00337d0
      sps(2,2) =  0.00135d0
      pps(2,2) =  0.02849d0
      ppp(2,2) =  0.01099d0
      sds(2,2) = -0.01119d0
      pds(2,2) = -0.01061d0
      pdp(2,2) =  0.01134d0
      dds(2,2) = -0.00759d0
      ddp(2,2) =  0.00495d0
      ddd(2,2) = -0.00016d0
c
c     -----------------------------------------------------------------
c     Cu:
c
      s0(3) =  0.79466d0
      p0(3) =  1.35351d0
      d0t(3) =  0.5d0*(0.37307d0 + 0.37180d0)
      d0e(3) =  0.5d0*(0.37307d0 + 0.37180d0)
c
c     first n.n.
c
      sss(3,1) = -0.07518d0
      sps(3,1) =  0.11571d0
      pps(3,1) =  0.19669d0
      ppp(3,1) =  0.01940d0
      sds(3,1) = -0.03107d0
      pds(3,1) = -0.03289d0
      pdp(3,1) =  0.01753d0
      dds(3,1) = -0.02566d0
      ddp(3,1) =  0.01800d0
      ddd(3,1) = -0.00408d0
c
c     second n.n.
c
      sss(3,2) = -0.00092d0
      sps(3,2) =  0.01221d0
      pps(3,2) =  0.05389d0
      ppp(3,2) =  0.00846d0
      sds(3,2) = -0.00852d0
      pds(3,2) = -0.00536d0
      pdp(3,2) =  0.00321d0
      dds(3,2) = -0.00451d0
      ddp(3,2) =  0.00241d0
      ddd(3,2) = -0.00029d0
c
c     -----------------------------------------------------------------
      return
      end
