c     *****************************************************************
c     PROGRAM jmbcocu.f
c
c     Last major modification: 27/2/12  inclusion of decimation for SGF
c     from ~/progs/supercell/new/surface/decimation/test25FeMgOFe_2.f
c
c
c     This program calculates the GMR for a general multilayer,
c     with general supercell configuration, and general growth direction.
c     Each layer is allowed to have a different geometry.
c     However, the in-plane geometry (defined by (a1,a2)) of the lattice
c     must be common to all layers. This is so that the 2 
c     geometries can be simulataneously diagonalised in-plane.
c
c     The in plane geometry of the substrate is defined by (a1,a2).
c     There are nlay layers :
c     Layers 1 and 2 define the LH semi-infinite substrate.
c     Layers nlay-1 and nlay define the RH semi-infinite substrate.      
c     Layers 3..(nlay-2) define the intermediate multilayer. 
c     The ith multilayer has perpendicular position given by a3(i),
c     and basis vsub(i,1-->nsub). So the lattice for the ith layer is
c           R = \sum{n1,n2,m} [n1.a1 + n2.a2 + a3(i)  + vsub(i,m)]

c     There are significant changes to subroutine param
c     The hopping between different atoms is now given by the geometric
c     mean, unless this is overridden --- as in the case of Fe-MgO.
c     There are two spins per atom type.

c     ******************************************************************
c     
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
c
CMPI      include 'mpif.h'
c
      parameter (nplx=75,nslx=1,nspinx=10,nsubx=4,nsubatx=4)
      parameter (nmatx=nslx*nsubx*nspinx,natomx=nslx*nspinx*nsubatx)
      parameter (nlayx=nplx*nslx)
      parameter (nnx=4,npmax=100) !! upto nnx.th N.N's
c                                    npmax - max No of points per plane
      parameter (numatx=8)       !! No of atom types
      parameter (nfoldx=100)
      complex*16 ec
      common/shift/xmgo_shift,ag_shift,xinas_shift
      common /fermi/ ef
      common/lbl/ pi, sq2, sq3
      common/data/nmat,nins,mlay,natom,nspin,nlay,nsub,numnn
      common/data1/frat,ifrac
      common/layer/ndiff
      common/ssstruct/a1(3),a2(3),a3(3,nlayx),vsub(3,nsubx,nlayx),
     $  dnn(nnx,nlayx,nlayx),ddnn(numatx,numatx,nnx),itype(nsubx,nlayx)
      dimension vtmp(3),dist(nnx)
      dimension b1(3),b2(3)
      dimension zresu(0:nlayx),zresd(0:nlayx)
      dimension zresud(0:nlayx),zresdu(0:nlayx)
      dimension rr(3)
      common/elements/iatfe(5),numat
      character*2 atname(nlayx)
      common/parallel/myid,numprocs
      dimension imp(0:nsubx)
      dimension xn(3)
      common/atmap/itypefe(3),itypeau(2)


      common/ssstructat/aa1(3),aa2(3),aa3(3,nlayx),vsubat(3,nsubx,nlayx)
     $,itypeat(nsubx,nlayx),imapl(nsubx),imapr(nsubx),nsubat
      common/fold/xfold(3,nfoldx),xshift(3,nfoldx),ifold(2,nfoldx),
     $  baib(3,3),nxfold,nfold,irecipa
      dimension ba1(3),ba2(3)
      dimension disttmp(numatx,numatx,nnx)

c     ------------------------------------------------------------
c     parallel processing initial setup commands
CMPI      call MPI_INIT(ierr)
CMPI      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
CMPI      call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
      

      myid=0
      numprocs=1


      if(myid.eq.0)then
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)' There are ',numprocs,' processors for this job '
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)
      endif
c     ------------------------------------------------------------
      pi    = 4.0d0*datan(1.0d0)
      sq2   = dsqrt(2.0d0)
      sq3   = dsqrt(3.0d0)
      ryd   = 13.6058d0
c     ------------------------------------------------------------
c     DATA
      nsl=1      !!!! Cannot change this in this code !!!!!
                 !!!! Can only layer principal layers

      nspin=10
      numnn=3

      mlay=4     !! No of Fe bi-layers
      nins=6     !! No of GaAs layers
      nau=15     !! Max No of Au bi-layers

      numat=8    !! No of atom types: one for each element

      nlay=nins+mlay+nau+6  !!! total No of layers inc 4 lead layers
      ndiff=nau

      ef=0.7505d0          !! Ef shifted into conduction band
      wm = 1.d-14

      nq=20

      dee=3.5         !!!   classic
      xmgo_shift=ef-dee/ryd
      ag_shift=ef-0.4635d0

      c=3.d8                  ! speed of light (m/s)
      hb=6.582119d-16         ! reduced Planks constant (eV s)
      alat=6.0583d-9          ! InAs lattice constant
      xxm=0.024d0*0.511d6/c**2 ! InAs effctv mass at G (eV) x electron mass (eV)
      xxn=1.d18                ! filling (No / cm^3 )
      write(*,*)"xxn?"
      read(*,*)xxn

      ef_shift=100**2*hb**2*(9*pi**4*xxn**2)**(1.d0/3.d0)/(2*xxm)!! InAs shift for filling xxn (eV)
      xinas_shift=ef-0.417717033707d0/ryd-ef_shift/ryd

      rmax=3.25d0*alat*(3*xxn*pi**2)**(1.d0/3.d0)  ! max radius in BZ where conduction occurs
      nqeff=int(dfloat(nq)/rmax) !!gives nq integration points inside region rmax
C     =================================================================
C     ATOMIC DATA FOR LEADS
C     =================================================================
c     in-plane atomic lattice vectors for substrate
c     in this code we assume that the LH and RH leads have the same 
C     lattice vectors and No of basis vectors

      nsubat=4
      natom=nsl*nspin*nsubat

      aa1(1)=0.5d0
      aa1(2)=0.5d0
      aa1(3)=0.d0
      aa2(1)=0.5d0
      aa2(2)=-0.5d0
      aa2(3)=0.d0

C     =================================================================
c     LH LEAD BASIS VECTORS
C     =================================================================

      do ilay=1,2
c       Out of plane lattice vector
        aa3(1,ilay)=0.d0
        aa3(2,ilay)=0.d0
        aa3(3,ilay)=(ilay-1.d0)*1/2.d0

c       Sublattice
        vsubat(1,1,ilay)=0.d0
        vsubat(2,1,ilay)=0.d0
        vsubat(3,1,ilay)=0.d0
        vsubat(1,2,ilay)=0.5d0
        vsubat(2,2,ilay)=0.d0
        vsubat(3,2,ilay)=0.d0
        vsubat(1,3,ilay)=0.25d0
        vsubat(2,3,ilay)=0.25d0
        vsubat(3,3,ilay)=-0.25d0
        vsubat(1,4,ilay)=-0.25d0
        vsubat(2,4,ilay)=0.25d0
        vsubat(3,4,ilay)=-0.25d0

c       Atom types
        itypeat(1,ilay)=1   !!!   Fe
        itypeat(2,ilay)=1   !!!   Fe
        itypeat(3,ilay)=1   !!!   Fe
        itypeat(4,ilay)=1   !!!   Fe
      enddo
c

C     =================================================================
C     SUPERCELL STRUCTURE
C     =================================================================

      nsub=4

c     2 in-plane lattice vectors
c     CUBIC
      a1(1)=0.5d0
      a1(2)=0.5d0
      a1(3)=0.d0
      a2(1)=0.5d0
      a2(2)=-0.5d0
      a2(3)=0.d0

c     ------------------- LH Fe lead --------------------
      do ilay=1,3
c       Out of plane lattice vector
        a3(1,ilay)=0.d0
        a3(2,ilay)=0.d0
        a3(3,ilay)=(ilay-1.d0)/2.d0

c       Sublattice
        vsub(1,1,ilay)=0.d0
        vsub(2,1,ilay)=0.d0
        vsub(3,1,ilay)=0.d0
        vsub(1,2,ilay)=0.5d0
        vsub(2,2,ilay)=0.d0
        vsub(3,2,ilay)=0.d0
        vsub(1,3,ilay)=0.25d0
        vsub(2,3,ilay)=0.25d0
        vsub(3,3,ilay)=-0.25d0
        vsub(1,4,ilay)=-0.25d0
        vsub(2,4,ilay)=0.25d0
        vsub(3,4,ilay)=-0.25d0


c       Atom types
        itype(1,ilay)=1   !!!   Fe
        itype(2,ilay)=1
        itype(3,ilay)=1
        itype(4,ilay)=1
      enddo

c     Atom types
      itype(1,3)=5   !!!   Fe interface with Au
      itype(2,3)=5   !!!   Fe interface with Au
      itype(3,3)=1   !!!   Fe
      itype(4,3)=1   !!!   Fe


c     ------------------- nau Au layers --------------------
      do ilay=4,nau+3
c       Out of plane lattice vector
        a3(1,ilay)=a3(1,ilay-1)+0.d0
        a3(2,ilay)=a3(2,ilay-1)+0.d0
        a3(3,ilay)=a3(3,ilay-1)+1.d0/sqrt(2.d0)

c       Sublattice
        vsub(1,1,ilay)=0.d0
        vsub(2,1,ilay)=0.d0
        vsub(3,1,ilay)=0.d0
        vsub(1,2,ilay)=0.5d0
        vsub(2,2,ilay)=0.d0
        vsub(3,2,ilay)=0.d0
        vsub(1,3,ilay)=0.25d0
        vsub(2,3,ilay)=0.25d0
        vsub(3,3,ilay)=-1/(2*sqrt(2.d0))
        vsub(1,4,ilay)=-0.25d0
        vsub(2,4,ilay)=0.25d0
        vsub(3,4,ilay)=-1/(2*sqrt(2.d0))

c       Atom types
        itype(1,ilay)=2   !!!   Au
        itype(2,ilay)=2   !!!   Au
        itype(3,ilay)=2   !!!   Au
        itype(4,ilay)=2   !!!   Au
      enddo

c     Atom types
      itype(1,4)=2   !!!   Au
      itype(2,4)=2   !!!   Au
      itype(3,4)=6   !!!   Au IF with Fe
      itype(4,4)=6   !!!   Au IF with Fe

      itype(1,nau+3)=6   !!!   Au IF with Fe
      itype(2,nau+3)=6   !!!   Au IF with Fe
      itype(3,nau+3)=2   !!!   Au
      itype(4,nau+3)=2   !!!   Au

c     ------------------- mlay Fe layers --------------------
      do ilay=nau+4,nau+4
c       Out of plane lattice vector
        a3(1,ilay)=0.d0
        a3(2,ilay)=0.d0
        a3(3,ilay)=a3(3,ilay-1)+0.25d0+0.5d0/sqrt(2.d0)

c       Sublattice
        vsub(1,1,ilay)=0.d0
        vsub(2,1,ilay)=0.d0
        vsub(3,1,ilay)=0.d0
        vsub(1,2,ilay)=0.5d0
        vsub(2,2,ilay)=0.d0
        vsub(3,2,ilay)=0.d0
        vsub(1,3,ilay)=0.25d0
        vsub(2,3,ilay)=0.25d0
        vsub(3,3,ilay)=-0.25d0
        vsub(1,4,ilay)=-0.25d0
        vsub(2,4,ilay)=0.25d0
        vsub(3,4,ilay)=-0.25d0

c       Atom types
        itype(1,ilay)=1   !!!   Fe
        itype(2,ilay)=1   !!!   Fe
        itype(3,ilay)=5   !!!   Fe interface with Au
        itype(4,ilay)=5   !!!   Fe interface with Au
      enddo

      do ilay=nau+5,nau+mlay+3
c       Out of plane lattice vector
        a3(1,ilay)=0.d0
        a3(2,ilay)=0.d0
        a3(3,ilay)=a3(3,ilay-1)+1.d0/2.d0

c       Sublattice
        vsub(1,1,ilay)=0.d0
        vsub(2,1,ilay)=0.d0
        vsub(3,1,ilay)=0.d0
        vsub(1,2,ilay)=0.5d0
        vsub(2,2,ilay)=0.d0
        vsub(3,2,ilay)=0.d0
        vsub(1,3,ilay)=0.25d0
        vsub(2,3,ilay)=0.25d0
        vsub(3,3,ilay)=-0.25d0
        vsub(1,4,ilay)=-0.25d0
        vsub(2,4,ilay)=0.25d0
        vsub(3,4,ilay)=-0.25d0


c       Atom types
        itype(1,ilay)=1   !!!   Fe
        itype(2,ilay)=1
        itype(3,ilay)=1
        itype(4,ilay)=1
      enddo

c     Atom types
      itype(1,nau+mlay+3)=7   !!!   Fe interface with Mg0
      itype(2,nau+mlay+3)=7   !!!   Fe interface with Mg0
      itype(3,nau+mlay+3)=1   !!!   Fe
      itype(4,nau+mlay+3)=1   !!!   Fe

c     ------------------- nins GaAs layers --------------------

      do ilay=nau+mlay+4,nlay
c       Out of plane lattice vector
        a3(1,ilay)=a3(1,ilay-1)+0.5d0
        a3(2,ilay)=a3(2,ilay-1)-0.d0
        a3(3,ilay)=a3(3,ilay-1)+0.5d0     !! approximate

c       Sublattice
        vsub(1,1,ilay)=0.d0
        vsub(2,1,ilay)=0.d0
        vsub(3,1,ilay)=0.d0
        vsub(1,2,ilay)=0.5d0
        vsub(2,2,ilay)=0.d0
        vsub(3,2,ilay)=0.d0
        vsub(1,3,ilay)=0.25d0
        vsub(2,3,ilay)=0.25d0
        vsub(3,3,ilay)=-0.25d0
        vsub(1,4,ilay)=-0.25d0
        vsub(2,4,ilay)=0.25d0
        vsub(3,4,ilay)=-0.25d0

c       Atom types
        itype(1,ilay)=3    !! Ga
        itype(2,ilay)=8    !! Vac
        itype(3,ilay)=4    !! As
        itype(4,ilay)=8    !! Vac
      enddo

C     =================================================================
c     RH LEAD BASIS VECTORS
C     =================================================================
      do ilay=nlay-1,nlay
c       Out of plane lattice vector
        aa3(1,ilay)=a3(1,ilay)
        aa3(2,ilay)=a3(2,ilay)
        aa3(3,ilay)=a3(3,ilay)

c       Sublattice
        vsubat(1,1,ilay)=0.d0
        vsubat(2,1,ilay)=0.d0
        vsubat(3,1,ilay)=0.d0
        vsubat(1,2,ilay)=0.5d0
        vsubat(2,2,ilay)=0.d0
        vsubat(3,2,ilay)=0.d0
        vsubat(1,3,ilay)=0.25d0
        vsubat(2,3,ilay)=0.25d0
        vsubat(3,3,ilay)=-0.25d0
        vsubat(1,4,ilay)=-0.25d0
        vsubat(2,4,ilay)=0.25d0
        vsubat(3,4,ilay)=-0.25d0


c       Atom types
        itypeat(1,ilay)=3
        itypeat(2,ilay)=8
        itypeat(3,ilay)=4
        itypeat(4,ilay)=8
      enddo

C     =================================================================
C     The map between Supercell sublattice and LH atomic sublattice :
c     imap:   supercell --> atomic
C     Defined s.t.     vsub(k)=vsubat(imap(k)) + atomic lattice vector

      imapl(1)=1
      imapl(2)=2
      imapl(3)=3
      imapl(4)=4

      imapr(1)=1
      imapr(2)=2
      imapr(3)=3
      imapr(4)=4

c
C     =================================================================
C     THIS SECTION TO GET NN DISTANCES ONLY

C     In and out of plane distances:
      do ilay=2,nlay
        call prestructij(ilay,ilay,dist)
        do i=1,numnn
          dnn(i,ilay,ilay)=dist(i)
        enddo

        call prestructij(ilay,ilay-1,dist)
        do i=1,numnn
          dnn(i,ilay,ilay-1)=dist(i)
          dnn(i,ilay-1,ilay)=dist(i)
        enddo
      end do

C     ddnn - atom-atom distances  ... this needs to be checked
C     Better to put NN distances in by hand as in next section.
      do iat=1,numat
        do jat=1,numat
          do inn=1,numnn
            disttmp(iat,jat,inn)=1.d10
          enddo
        enddo
      enddo
      do ilay=2,nlay
        call prestructijnew(ilay,ilay,disttmp)
        call prestructijnew(ilay,ilay-1,disttmp)
        call prestructijnew(ilay-1,ilay,disttmp)
      enddo
C     =================================================================
C     NN DISTANCES
C     1,2 Fe  ; 3 Mg ; 4 O 
C     ddnn(type1,type2,NN)

      itypefe(1)=1
      itypefe(2)=5
      itypefe(3)=7
      itypeau(1)=2
      itypeau(2)=6

      do i=1,numat
        do j=1,numat
          do k=1,numnn
            ddnn(i,j,k)=0.d0
          enddo
        enddo
      enddo

C     Fe-Fe
      do i=1,3
        do j=1,3
          ii=itypefe(i)
          jj=itypefe(j)
          ddnn(ii,jj,1)=0.4330127019d0
          ddnn(ii,jj,2)=0.5d0
          ddnn(ii,jj,3)=1.d0/sqrt(2.d0)
        enddo
      enddo

C     Au-Au
      do i=1,2
        do j=1,2
          ii=itypeau(i)
          jj=itypeau(j)
          ddnn(ii,jj,1)=0.5d0
          ddnn(ii,jj,2)=1.d0/sqrt(2.d0)
          ddnn(ii,jj,3)=0.8660254038d0
        enddo
      enddo

C     Ga-Ga
      ddnn(3,3,1)=1/sqrt(2.d0)
      ddnn(3,3,2)=1.d0
      ddnn(3,3,3)=1.222474d0

C     As-As
      ddnn(4,4,1)=1/sqrt(2.d0)
      ddnn(4,4,2)=1.d0
      ddnn(4,4,3)=1.222474d0

C     Vac-Vac
      ddnn(8,8,1)=sqrt(3.d0)/4.d0
      ddnn(8,8,2)=1/sqrt(2.d0)
      ddnn(8,8,3)=0.829156d0

C     Ga-As
      ddnn(3,4,1)=sqrt(3.d0)/4.d0
      ddnn(3,4,2)=0.829156d0
      ddnn(3,4,3)=1.08972d0
      ddnn(4,3,1)=sqrt(3.d0)/4.d0
      ddnn(4,3,2)=0.829156d0
      ddnn(4,3,3)=1.08972d0

C     Ga-Vac
      ddnn(3,8,1)=sqrt(3.d0)/4.d0
      ddnn(3,8,2)=0.5d0
      ddnn(3,8,3)=0.829156d0
      ddnn(8,3,1)=sqrt(3.d0)/4.d0
      ddnn(8,3,2)=0.5d0
      ddnn(8,3,3)=0.829156d0

C     Vac-As
      ddnn(4,8,1)=sqrt(3.d0)/4.d0
      ddnn(4,8,2)=0.5d0
      ddnn(4,8,3)=0.829156d0
      ddnn(8,4,1)=sqrt(3.d0)/4.d0
      ddnn(8,4,2)=0.5d0
      ddnn(8,4,3)=0.829156d0

C     Fe-GaAs
      do i=1,3
        ii=itypefe(i)
        ddnn(ii,4,1)=0.4330127019d0
        ddnn(ii,3,2)=0.5d0
        ddnn(ii,4,2)=0.5d0
        ddnn(ii,3,3)=0.7071067812d0
        ddnn(ii,4,3)=0.7071067812d0

        ddnn(4,ii,1)=0.4330127019d0
        ddnn(3,ii,2)=0.5d0
        ddnn(4,ii,2)=0.5d0
        ddnn(3,ii,3)=0.7071067812d0
        ddnn(4,ii,3)=0.7071067812d0
      enddo

C     Fe-Au
      do i=1,3
        do j=1,2
          ii=itypefe(i)
          jj=itypeau(j)
          ddnn(ii,jj,1)=0.5d0
          ddnn(ii,jj,2)=0.6035533906d0
          ddnn(ii,jj,3)=1.d0/sqrt(2.d0)
          ddnn(jj,ii,1)=0.5d0
          ddnn(jj,ii,2)=0.6035533906d0
          ddnn(jj,ii,3)=1.d0/sqrt(2.d0)
        enddo
      enddo

C     =================================================================
C     !!!!!!! OUTPUT ATOMIC POSITIONS FOR RASMOL VIEWER !!!!!!
C     load into rasmol with command :    > load xyz pos0.dat
      atname(1)="Fe"
      atname(2)="Au"
      atname(3)="Ga"
      atname(4)="As"
      atname(5)="Fe"
      atname(6)="Au"
      atname(7)="Fe"
      atname(8)="Vac"


c     whole cluster
      open(20,file='pos0.dat')
      idum0=0

      do ilay=1,nlay
        do isub=1,nsub
          do i1=-nlay,nlay
            do i2=-nlay,nlay
              do k=1,3
                rr(k)=a3(k,ilay)+vsub(k,isub,ilay)+i1*a1(k)+i2*a2(k)
              enddo
              if(abs(rr(1)).lt.0.5001.and.abs(rr(2)).lt.0.5001)then
                idum0=idum0+1
              endif
            enddo
          enddo
        enddo
      enddo


      if(myid.eq.0)then
        write(20,*)idum0
        write(20,*)"foo"
        do ilay=1,nlay
          do isub=1,nsub
            do i1=-nlay,nlay
              do i2=-nlay,nlay
                do k=1,3
                  rr(k)=a3(k,ilay)+vsub(k,isub,ilay)+i1*a1(k)+i2*a2(k)
                enddo
                if(abs(rr(1)).lt.0.5001.and.abs(rr(2)).lt.0.5001)then
                  write(20,'(a2,5x,3f12.8)')atname(itype(isub,ilay)),
     $           (4*rr(k),k=1,3)
                endif
              enddo
            enddo
          enddo
        enddo
      endif
C     =================================================================
c     Dimensionality Checks

      nmat=nsl*nspin*nsub

      if(myid.eq.0)then
        if(numnn.gt.nnx)then
          write(*,*)' ERROR : numnn > nnx'
          stop
        endif

        if(nsub.gt.nsubx)then
          write(*,*)' ERROR : nsub > nsubx'
          stop
        endif

        if(ndiff.gt.nlayx)then
          write(*,*)' ERROR : ndiff > nlayx ',ndiff,nlayx
          stop
        endif
      endif

      if(numat.gt.numatx)then
        write(*,*)"ERROR: numat>numatx",numat,numatx
        stop
      endif

c     construct the growth direction:
      call cross(a1,a2,xn,3)
      xnorm=dot(xn,xn,3)
      do i=1,3
        xn(i)=xn(i)/sqrt(xnorm)
      enddo
C     =================================================================
c
      if(myid.eq.0)then
        write(*,*)' Supercell GMR'
        write(*,*)
        write(*,*)' Substrate lattice vectors'
        write(*,*)(a1(i),i=1,3)
        write(*,*)(a2(i),i=1,3)
        write(*,*)
        write(*,*)
        write(*,*)' Growth direction'
        write(*,*)(xn(i),i=1,3)
        write(*,*)
        write(*,*)' nq,nqeff =',nq,nqeff
        write(*,*)" GaAs ef_shift, rmax = ",ef_shift/ryd,rmax
        write(*,*)' ef =',ef
        write(*,*)

C       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C       this prints out all interplanar distances
        do ii=1,nlay
          write(*,'(72("-"))')
          write(*,*)' The ',ii,'th Layer'
          write(*,*)
          write(*,'("a3 = ",3f12.8)')(a3(k,ii),k=1,3)
          write(*,*)
          write(*,*)' Sub-lattice vectors'
          do i=1,nsub
            write(*,'("vsub(",i2,")=",3f8.4)')i,vsub(1,i,ii),
     $vsub(2,i,ii),vsub(3,i,ii)
          enddo
          if(ii.ge.2)then
            write(*,*)
            write(*,*)"In-Plane NN distances"
            do i=1,numnn
              write(*,*)ii,dnn(i,ii,ii)
            enddo
            write(*,*)
            write(*,*)"Out of Plane NN distances"
            do i=1,numnn
              write(*,*)ii,dnn(i,ii,ii-1)
            enddo
          endif
        enddo
        write(*,'(72("-"))')

C       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        write(*,'(65("="))')
        write(*,*)"      CALCULATED NN DISTANCES --- check these"
        write(*,'(65("="))')
        write(*,'("  at1   at2  ",100("     dNN",i1,"        "))')
     $(k,k=1,numnn-1)
        do i=1,numat
          do j=1,numat
            write(*,'(2i5,100d18.10)')i,j,(disttmp(i,j,k),k=1,numnn)
          enddo
        enddo
        write(*,'(65("="))')
        write(*,'(65("="))')
      endif
c
c     -----------------------------------------------------------------

      call param

c     -----------------------------------------------------------------
c     determine the reciprocal lattice structure
      call recip(a1,a2,b1,b2,irecip)
      if(myid.eq.0)then
        write(*,*)
        write(*,*)' perpr reciprocal lattice vectors'
        write(*,'(3f12.4)')(b1(i),i=1,3)
        write(*,'(3f12.4)')(b2(i),i=1,3)
        write(*,*)
      endif
c     -----------------------------------------------------------------

c     determine the atomic reciprocal lattice vectors ba1,ba2
      call recip(aa1,aa2,ba1,ba2,irecipa)
      if(myid.eq.0)then
        write(*,*)
        write(*,*)'----------------------------------------------------'
        write(*,*)
        write(*,*)' perpr atomic reciprocal lattice vectors'
        write(*,'(3f12.4)')(ba1(i),i=1,3)
        write(*,'(3f12.4)')(ba2(i),i=1,3)
        write(*,*)
        write(*,*)'----------------------------------------------------'
        write(*,*)
      endif

c     -----------------------------------------------------------------
c     now find the 'folding' vectors j[1:nfold]
      call folding(b1,b2,ba1,ba2)

C     =================================================================
C     CHECK THAT IMAP IS CORRECTLY DEFINED
C     ie.    vsub(k)=vsubat(imap(k)) + n1.aa1 + n2.aa2

c     LH lead
      do ilay=1,2
        do isub=1,nsub
          do k=1,3
            vtmp(k)=vsub(k,isub,ilay)-vsubat(k,imapl(isub),ilay)
          enddo
          vba1=(abs(dot(vtmp,ba1,3))+1.d-10)/(2*pi)
          vba2=(abs(dot(vtmp,ba2,3))+1.d-10)/(2*pi)
          if(abs(int(vba1)-vba1).gt.1.d-10.or.
     $abs(int(vba2)-vba2).gt.1.d-10)then
            write(*,*)"imapl is wrong",vba1,vba2
            stop
          endif
          if(itype(isub,ilay).ne.itypeat(imapl(isub),ilay))then
            write(*,*)"imapl is wrong",itype(isub,ilay),
     $itypeat(imapl(isub),ilay)
            stop
          endif
        enddo
      enddo

c     RH lead
      do ilay=nlay-1,nlay
        do isub=1,nsub
          do k=1,3
            vtmp(k)=vsub(k,isub,ilay)-vsubat(k,imapr(isub),ilay)
          enddo
          vba1=(abs(dot(vtmp,ba1,3))+1.d-10)/(2*pi)
          vba2=(abs(dot(vtmp,ba2,3))+1.d-10)/(2*pi)
          if(abs(int(vba1)-vba1).gt.1.d-10.or.
     $abs(int(vba2)-vba2).gt.1.d-10)then
            write(*,*)"imapr is wrong",vba1,vba2
            stop
          endif
          if(itype(isub,ilay).ne.itypeat(imapr(isub),ilay))then
            write(*,*)"imapr is wrong",itype(isub,ilay),
     $itypeat(imapr(isub),ilay)
            stop
          endif
        enddo
      enddo

C     =================================================================
C     DO THE CALCULATION


      ec = dcmplx(ef,wm)
      if(myid.eq.0)write(*,*)' Complex Energy = ',ec
      do in=0,ndiff
        zresu(in)=0.d0
        zresd(in)=0.d0
        zresud(in)=0.d0
        zresdu(in)=0.d0
      enddo
      call sumk(irecip,ec,nqeff,rmax,b1,b2,zresu,zresd,zresud,zresdu)
c
      if(myid.eq.0)then
      write(*,*)
      write(*,*)
      write(*,*)'****************************************************'
      write(*,*)'****************************************************'
      write(*,*)' UP SPIN'
      do in=0,ndiff
        write(*,*)2*in,dreal(zresu(in))/nsub
      enddo
      write(*,*)
      write(*,*)'****************************************************'
      write(*,*)'****************************************************'
      write(*,*)' DOWN SPIN'
      do in=0,ndiff
        write(*,*)2*in,dreal(zresd(in))/nsub
      enddo
      write(*,*)
      write(*,*)'****************************************************'
      write(*,*)'****************************************************'
      write(*,*)' UP-DOWN SPIN'
      do in=0,ndiff
        write(*,*)2*in,dreal(zresud(in))/nsub
      enddo
      write(*,*)
      write(*,*)'****************************************************'
      write(*,*)'****************************************************'
      write(*,*)' DOWN-UP SPIN'
      do in=0,ndiff
        write(*,*)2*in,dreal(zresdu(in))/nsub
      enddo
      write(*,*)
      write(*,*)'****************************************************'
      write(*,*)'****************************************************'
      write(*,*)' (UP + DOWN - 2*AF)/(UD + DU)'
      do in=0,ndiff
        if(cdabs(zresud(in)*zresdu(in)).ne.0)then
          write(*,*)2*in,
     $      dreal(zresu(in)+zresd(in)-zresud(in)-zresdu(in))/
     $      dreal(zresud(in)+zresdu(in))
        else
          write(*,*)' ERROR : resaf = 0 '
        endif
      enddo
      endif
c
CMPI      call MPI_FINALIZE(iend)
      stop
      end
c
c     *****************************************************************
c
CMPI      subroutine sumk(irecip,zener,nk,b1,b2,zresu,zresd,zresud,zresdu)
CMPI      implicit double precision (a-h,o-y)
CMPI      implicit complex*16 (z)
CMPIc
CMPI      include 'mpif.h'
CMPIc
CMPI      parameter (nplx=75,nslx=1,nspinx=10,nsubx=4,nsubatx=4)
CMPI      parameter (nmatx=nslx*nsubx*nspinx,natomx=nslx*nspinx*nsubatx)
CMPI      parameter (nlayx=nplx*nslx)
CMPI      parameter (numatx=8)       !! No of atom types
CMPI      parameter (nnx=4,npmax=100) !! upto nnx.th N.N's
CMPI      parameter (numkx=1500*1500)
CMPI      parameter (nfoldx=100)
CMPI      common/lbl/ pi, sq2, sq3
CMPI      common/data/nmat,nins,mlay,natom,nspin,nlay,nsub,numnn
CMPI      common/data1/frat,ifrac
CMPI      common/layer/ndiff
CMPI      common/ssstructat/aa1(3),aa2(3),aa3(3,nlayx),vsubat(3,nsubx,nlayx)
CMPI     $,itypeat(nsubx,nlayx),imapl(nsubx),imapr(nsubx),nsubat
CMPI      dimension b1(3),b2(3),d1(3),d2(3),xk(3)
CMPI      dimension zresu(0:nlayx),zresd(0:nlayx)
CMPI      dimension zresud(0:nlayx),zresdu(0:nlayx)
CMPI      dimension zconu(0:nlayx),zcond(0:nlayx)
CMPI      dimension zconud(0:nlayx),zcondu(0:nlayx)
CMPI      common/parallel/myid,numprocs
CMPI      dimension ik(numkx),jk(numkx)
CMPI      dimension zans(4*nlayx+5),zres(4*nlayx+5)
CMPI      dimension ibuff(2)
CMPI      dimension numcalls(1000)
CMPI      common/fold/xfold(3,nfoldx),xshift(3,nfoldx),ifold(2,nfoldx),
CMPI     $  baib(3,3),nxfold,nfold,irecipa
CMPI      integer status(MPI_STATUS_SIZE)
CMPI
CMPIc     -----------------------------------------------------------------
CMPI      if(nk.ne.1)dq=1.d0/dfloat(nk-1)
CMPI      if(nk.eq.1)dq=1.d0
CMPI      max=nk-1
CMPIc     -----------------------------------------------------------------
CMPIC     Set up basis vectors and k-points
CMPI      if(irecip.eq.0)then
CMPI        write(*,*)' SUMK : irecip = 0 ---- NOT CODED '
CMPI        stop
CMPIc     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
CMPI      elseif(irecip.eq.1)then   !!    CUBIC
CMPI        do k=1,3
CMPI          d1(k)=b1(k)/2.d0
CMPI          d2(k)=b2(k)/2.d0
CMPI        enddo
CMPI
CMPI        ijsum=0
CMPI        do i=0,max
CMPI          do j=0,i
CMPI            ijsum=ijsum+1
CMPI            ik(ijsum)=i
CMPI            jk(ijsum)=j
CMPI          enddo
CMPI        enddo
CMPIc     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
CMPI      elseif(irecip.eq.2)then  !!     PRIMITIVE-RECTANGULAR
CMPI        do k=1,3
CMPI          d1(k)=b1(k)/2.d0
CMPI          d2(k)=b2(k)/2.d0
CMPI        enddo
CMPI
CMPI        ijsum=0
CMPI        do i=0,max
CMPI          do j=0,max
CMPI            ijsum=ijsum+1
CMPI            ik(ijsum)=i
CMPI            jk(ijsum)=j
CMPI          enddo
CMPI        enddo
CMPIc     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
CMPI      elseif(irecip.eq.3)then  !!     HEXAGONAL
CMPI        do k=1,3
CMPI          d1(k)=b1(k)/2.d0
CMPI          d2(k)=(b1(k)+b2(k))/3.d0 - b1(k)/2.d0
CMPI        enddo
CMPI
CMPI        ijsum=0
CMPI        do i=0,max
CMPI          do j=0,i
CMPI            ijsum=ijsum+1
CMPI            ik(ijsum)=i
CMPI            jk(ijsum)=j
CMPI          enddo
CMPI        enddo
CMPIc     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
CMPI      elseif(irecip.eq.4)then
CMPIc       CENTRED-RECTANGULAR
CMPI        write(*,*)'CENTRED-RECTANGULAR NOT CODED YET '
CMPI        stop
CMPI      endif
CMPIc     -----------------------------------------------------------------
CMPI
CMPI
CMPIC     Now begin set up of parallelisation
CMPI
CMPI      numk=ijsum
CMPI      if(numk.gt.numkx)then
CMPI        write(*,*)' ERROR SUMK : numk > numkx '
CMPI        stop
CMPI      endif
CMPI
CMPIC     Want on average 10 calls to each processor :
CMPIC     So give each processor ncols k-points at a time
CMPI      ncols=numk/(10*(numprocs-1))
CMPI      if(ncols.eq.0)ncols=1
CMPI      nrows=numk/ncols
CMPI      if(nrows*ncols.lt.numk)nrows=nrows+1
CMPI      if(myid.eq.0)then
CMPI        write(*,*)
CMPI        write(*,'(72("-"))')
CMPI        write(*,*)" there are ",numk,"k-points"
CMPI        write(*,*)" Each processor gets ",ncols," k-points ","per call"
CMPI        write(*,*)" There are ",nrows/(numprocs-1)," calls on average"
CMPI        write(*,'(72("-"))')
CMPI        write(*,*)
CMPI        numsent = 0
CMPI        do i=1,min(numprocs-1,nrows)
CMPI          ibuff(1)=numsent*ncols+1
CMPI          ibuff(2)=min((numsent+1)*ncols,numk)
CMPI          call MPI_SEND(ibuff,2,MPI_INTEGER,i,i, MPI_COMM_WORLD,ierr)
CMPI          numsent = numsent+1
CMPI        enddo
CMPI      endif
CMPI
CMPI      if(myid.gt.0.and.myid.le.nrows)then
CMPI90      call MPI_RECV(ibuff,2, MPI_INTEGER,0,
CMPI     &                 MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
CMPI        if (status(MPI_TAG).eq.0)goto 200
CMPI
CMPIc       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
CMPIC         NOW LET EACH NODE DO ITS SUM
CMPI          nktot=0
CMPI          do in=0,ndiff
CMPI            zresu(in) = 0.d0
CMPI            zresd(in) = 0.d0
CMPI            zresud(in) = 0.d0
CMPI            zresdu(in) = 0.d0
CMPI          enddo
CMPI
CMPI          do kk=ibuff(1),ibuff(2)
CMPI            i=ik(kk)
CMPI            j=jk(kk)
CMPI            x=dq*i
CMPI            y=dq*j
CMPI
CMPIc     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CMPIc           find the folding points
CMPI            call testk(x,y,b1,b2,irecip)
CMPI            if(nxfold.ne.nsub/nsubat)then
CMPI              write(*,*)' ERROR SUMK : nxfold ne nsub/nsubat',
CMPI     $          nxfold,nsub,nsubat
CMPI              stop
CMPI            endif
CMPIc     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CMPI
CMPI            do k=1,3
CMPI              xk(k)=x*d1(k)+y*d2(k)
CMPI            enddo
CMPIc
CMPI            ifail=0
CMPI            call cond(zener,xk,zconu,zcond,zconud,zcondu,ifail)
CMPI
CMPI            call weight(i,j,max,irecip,iwght)
CMPI            if(ifail.eq.0)then
CMPI              nktot=nktot+iwght
CMPI              do in=0,ndiff
CMPI                zresu(in) = zresu(in) + zconu(in)*iwght
CMPI                zresd(in) = zresd(in) + zcond(in)*iwght
CMPI                zresud(in) = zresud(in) + zconud(in)*iwght
CMPI                zresdu(in) = zresdu(in) + zcondu(in)*iwght
CMPI              enddo
CMPI            endif
CMPI          enddo
CMPIc       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
CMPI
CMPIC         Now send back all the results in a single array zans
CMPI          zans(1)=dfloat(nktot)
CMPI          do in=0,ndiff
CMPI            zans(in+2+0*ndiff)=zresu(in)
CMPI            zans(in+3+1*ndiff)=zresd(in)
CMPI            zans(in+4+2*ndiff)=zresud(in)
CMPI            zans(in+5+3*ndiff)=zresdu(in)
CMPI          enddo
CMPI
CMPI          irow = status(MPI_TAG)
CMPI          call MPI_SEND(zans,4*ndiff+5,MPI_DOUBLE_COMPLEX,0,
CMPI     $      irow,MPI_COMM_WORLD,ierr)
CMPI          go to 90
CMPI
CMPI200     continue
CMPI      endif
CMPI
CMPI
CMPIc     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CMPI      do i=1,numprocs
CMPI        numcalls(i)=0
CMPI      enddo
CMPI      if(myid.eq.0)then
CMPI        do k=1,4*ndiff+5
CMPI          zres(k)=0.d0
CMPI        enddo
CMPI        do i = 1,nrows
CMPI          call MPI_RECV(zans,4*ndiff+5, MPI_DOUBLE_COMPLEX,
CMPI     &                 MPI_ANY_SOURCE, MPI_ANY_TAG,
CMPI     &                 MPI_COMM_WORLD, status, ierr)
CMPI          nsender=status(MPI_SOURCE)
CMPI          numcalls(nsender+1)=numcalls(nsender+1)+1
CMPI          do k=1,4*ndiff+5
CMPI              zres(k)=zres(k)+zans(k)
CMPI          enddo
CMPI
CMPI          write(*,*)'node',nsender,' : ',numsent,'sets out of ',nrows,
CMPI     &      'received so far.'
CMPI
CMPI          if(numsent.lt.nrows)then        ! send another row
CMPI            ibuff(1)=numsent*ncols+1
CMPI            ibuff(2)=min((numsent+1)*ncols,numk)
CMPI
CMPI            call MPI_SEND(ibuff,2,MPI_INTEGER,
CMPI     &                   nsender,numsent+1, MPI_COMM_WORLD, ierr)
CMPI            numsent=numsent+1
CMPI          else      ! Tell nsender that there is no more work
CMPI            call MPI_SEND(MPI_BOTTOM, 0, MPI_DOUBLE_COMPLEX,
CMPI     &                    nsender, 0, MPI_COMM_WORLD, ierr)
CMPI          endif
CMPI        enddo
CMPI
CMPI
CMPIC       Extract the info received from the array
CMPI        xktot=dreal(zres(1))
CMPI        do in=0,ndiff
CMPI          zresu(in)=zres(in+2+0*ndiff)
CMPI          zresd(in)=zres(in+3+1*ndiff)
CMPI          zresud(in)=zres(in+4+2*ndiff)
CMPI          zresdu(in)=zres(in+5+3*ndiff)
CMPI        enddo
CMPI
CMPI        write(*,*)' processor                 number of calls'
CMPI        do i=1,numprocs
CMPI          write(*,*)i,numcalls(i)
CMPI        enddo
CMPI
CMPI        do in=0,ndiff
CMPI          zresu(in) = zresu(in)/xktot
CMPI          zresd(in) = zresd(in)/xktot
CMPI          zresud(in) = zresud(in)/xktot
CMPI          zresdu(in) = zresdu(in)/xktot
CMPI        enddo
CMPI
CMPIC       Send these results to all processors
CMPI        do i=1,numprocs-1
CMPI          call MPI_SEND(zresu,ndiff+1,MPI_DOUBLE_COMPLEX,
CMPI     &      i,1000,MPI_COMM_WORLD,ierr)
CMPI          call MPI_SEND(zresd,ndiff+1,MPI_DOUBLE_COMPLEX,
CMPI     &      i,1001,MPI_COMM_WORLD,ierr)
CMPI          call MPI_SEND(zresud,ndiff+1,MPI_DOUBLE_COMPLEX,
CMPI     &      i,1001,MPI_COMM_WORLD,ierr)
CMPI          call MPI_SEND(zresdu,ndiff+1,MPI_DOUBLE_COMPLEX,
CMPI     &      i,1001,MPI_COMM_WORLD,ierr)
CMPI        enddo
CMPI      endif
CMPI
CMPI      if(myid.gt.0)then
CMPI        call MPI_RECV(zresu,ndiff+1,MPI_DOUBLE_COMPLEX,0,
CMPI     &                 1000,MPI_COMM_WORLD,status,ierr)
CMPI        call MPI_RECV(zresd,ndiff+1,MPI_DOUBLE_COMPLEX,0,
CMPI     &                 1001,MPI_COMM_WORLD,status,ierr)
CMPI        call MPI_RECV(zresud,ndiff+1,MPI_DOUBLE_COMPLEX,0,
CMPI     &                 1001,MPI_COMM_WORLD,status,ierr)
CMPI        call MPI_RECV(zresdu,ndiff+1,MPI_DOUBLE_COMPLEX,0,
CMPI     &                 1001,MPI_COMM_WORLD,status,ierr)
CMPI
CMPI      endif
CMPI      return
CMPI      end
c
c     *****************************************************************
c
      subroutine sumk(irecip,zener,nk,rmax,b1,b2,zresu,zresd,zresud,
     $zresdu)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter (nplx=75,nslx=1,nspinx=10,nsubx=4,nsubatx=4)
      parameter (nmatx=nslx*nsubx*nspinx,natomx=nslx*nspinx*nsubatx)
      parameter (nlayx=nplx*nslx)
      parameter (numatx=8)       !! No of atom types
      parameter (nnx=4,npmax=100) !! upto nnx.th N.N's
      parameter (numkx=1500*1500)
      parameter (nfoldx=100)
      common/lbl/ pi, sq2, sq3
      common/data/nmat,nins,mlay,natom,nspin,nlay,nsub,numnn
      common/data1/frat,ifrac
      common/layer/ndiff
      common/ssstructat/aa1(3),aa2(3),aa3(3,nlayx),vsubat(3,nsubx,nlayx)
     $,itypeat(nsubx,nlayx),imapl(nsubx),imapr(nsubx),nsubat
      dimension b1(3),b2(3),d1(3),d2(3),xk(3)
      dimension zresu(0:nlayx),zresd(0:nlayx)
      dimension zresud(0:nlayx),zresdu(0:nlayx)
      dimension zconu(0:nlayx),zcond(0:nlayx)
      dimension zconud(0:nlayx),zcondu(0:nlayx)
      common/parallel/myid,numprocs
      dimension ik(numkx),jk(numkx)
      dimension zans(4*nlayx+5),zres(4*nlayx+5)
      dimension ibuff(2)
      dimension numcalls(1000)
      common/fold/xfold(3,nfoldx),xshift(3,nfoldx),ifold(2,nfoldx),
     $  baib(3,3),nxfold,nfold,irecipa


c     -----------------------------------------------------------------
      if(nk.ne.1)dq=1.d0/dfloat(nk-1)
      if(nk.eq.1)dq=1.d0
      max=nk-1
      nktot=0
c     -----------------------------------------------------------------
      if(irecip.eq.0)then
        write(*,*)' SUMK : irecip = 0 ---- NOT CODED '
        stop
      endif
      if(irecip.eq.1)goto 10
      if(irecip.eq.2)goto 20
      if(irecip.eq.3)goto 30
      if(irecip.eq.4)goto 40
c     -----------------------------------------------------------------
c     CUBIC
10    continue
      do k=1,3
        d1(k)=b1(k)/2.d0
        d2(k)=b2(k)/2.d0
      enddo

      do i=0,max
        x=dq*i
        do j=0,i
          y=dq*j

          do k=1,3
            xk(k)=x*d1(k)+y*d2(k)
          enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c         determine weights
          iwght=8
          if(i.eq.j)iwght=4
          if(j.eq.0)iwght=4
          if(i.eq.max)iwght=4
          if(i.eq.0)iwght=1
          if(j.eq.max)iwght=1
          if((i.eq.max).and.(j.eq.0))iwght=2
          nktot=nktot+iwght
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C         the next line drops the integration outside the region rmax
          if(x.le.rmax.and.y.le.rmax)then
c           find the folding points
            call testk(x,y,b1,b2,irecip)
C           write(*,*)
C           write(*,'(i3," folding k-points at k = ",3f10.6)')
C    $        nxfold,(xk(k),k=1,3)
C           do iff=1,nxfold
C             write(*,'(3f10.6)')(xfold(k,iff),k=1,3)
C           enddo
            if(nxfold.ne.nsub/nsubat)then
              write(*,*)' ERROR SUMK : nxfold ne nsub/nsubat',
     $          nxfold,nsub,nsubat
              stop
            endif
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            ifail=0
            call cond(zener,xk,zconu,zcond,zconud,zcondu,ifail)
            if(ifail.ne.0)then
              write(*,*)i,j
            endif
c
            do in=0,ndiff
              zresu(in) = zresu(in) + zconu(in)*iwght
              zresd(in) = zresd(in) + zcond(in)*iwght
              zresud(in) = zresud(in) + zconud(in)*iwght
              zresdu(in) = zresdu(in) + zcondu(in)*iwght
            enddo
C     ?????????????????????????????????????????????????????????????????
            write(30,*)x,y,dreal(zconu(ndiff))
            write(31,*)x,y,dreal(zcond(ndiff))
            write(32,*)x,y,dreal(zconud(ndiff))
            write(33,*)x,y,dreal(zcondu(ndiff))
C     ?????????????????????????????????????????????????????????????????
          endif
C         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        enddo
      enddo
      goto 5000
c     -----------------------------------------------------------------
c     PRIMITIVE-RECTANGULAR
20    continue
      do k=1,3
        d1(k)=b1(k)/2.d0
        d2(k)=b2(k)/2.d0
      enddo
      do i=0,max
        x=dq*i
        do j=0,max
          y=dq*j
          do k=1,3
            xk(k)=x*d1(k)+y*d2(k)
          enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c         determine weights
          iwght=4
          if(i.eq.0)iwght=2
          if(j.eq.0)iwght=2
          if(i.eq.max)iwght=2
          if(j.eq.max)iwght=2
          if((i.eq.0).and.(j.eq.0))iwght=1
          if((i.eq.0).and.(j.eq.max))iwght=1
          if((i.eq.max).and.(j.eq.0))iwght=1
          if((i.eq.max).and.(j.eq.max))iwght=1
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c         find the folding points
          call testk(x,y,b1,b2,irecip)
C         write(*,*)
C         write(*,'(i3," folding k-points at k = ",3f10.6)')
C    $      nxfold,(xk(k),k=1,3)
C         do iff=1,nxfold
C           write(*,'(3f10.6)')(xfold(k,iff),k=1,3)
C         enddo
          if(nxfold.ne.nsub/nsubat)then
            write(*,*)' ERROR SUMK : nxfold ne nsub/nsubat',
     $        nxfold,nsub,nsubat
            stop
          endif
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          nktot=nktot+iwght
          ifail=0
          call cond(zener,xk,zconu,zcond,zconud,zcondu,ifail)
          if(ifail.ne.0)then
            write(*,*)i,j
          endif
c
          do in=0,ndiff
            zresu(in) = zresu(in) + zconu(in)*iwght
            zresd(in) = zresd(in) + zcond(in)*iwght
            zresud(in) = zresud(in) + zconud(in)*iwght
            zresdu(in) = zresdu(in) + zcondu(in)*iwght
          enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        enddo
      enddo
      goto 5000
c     -----------------------------------------------------------------
c     HEXAGONAL
30    continue
      do k=1,3
        d1(k)=b1(k)/2.d0
        d2(k)=(b1(k)+b2(k))/3.d0
      enddo
      do i=0,max
        x=dq*i
        do j=0,i
          y=dq*j
          do k=1,3
            xk(k)=x*d1(k)+y*d2(k)
          enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c         determine weights
          iwght=12
          if(i.eq.j)iwght=6
          if(j.eq.0)iwght=6
          if(i.eq.max)iwght=4
          if(i.eq.0)iwght=1
          if(j.eq.max)iwght=2
          if((i.eq.max).and.(j.eq.0))iwght=3
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c         find the folding points
          call testk(x,y,b1,b2,irecip)
C         write(*,*)
C         write(*,'(i3," folding k-points at k = ",3f10.6)')
C    $      nxfold,(xk(k),k=1,3)
C         do iff=1,nxfold
C           write(*,'(3f10.6)')(xfold(k,iff),k=1,3)
C         enddo
          if(nxfold.ne.nsub/nsubat)then
            write(*,*)' ERROR SUMK : nxfold ne nsub/nsubat',
     $        nxfold,nsub,nsubat
            stop
          endif
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          nktot=nktot+iwght
          ifail=0
          call cond(zener,xk,zconu,zcond,zconud,zcondu,ifail)
          if(ifail.ne.0)then
            write(*,*)i,j
          endif
c
          do in=0,ndiff
            zresu(in) = zresu(in) + zconu(in)*iwght
            zresd(in) = zresd(in) + zcond(in)*iwght
            zresud(in) = zresud(in) + zconud(in)*iwght
            zresdu(in) = zresdu(in) + zcondu(in)*iwght
          enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        enddo
      enddo
      goto 5000
c     -----------------------------------------------------------------
c     CENTRED-RECTANGULAR
40    continue
c     ?????????????????????????????????????????????????????????????????
      write(*,*)'CENTRED-RECTANGULAR NOT CODED YET '
      stop
c     ?????????????????????????????????????????????????????????????????
C     b1b1=dot(b1,b1,3)
C     b1b2=dot(b1,b2,3)
C     alpha=b1b1/(2.d0*(b1b1+b1b2))
C     do k=1,3
C       d1(k)=(b1(k)-b2(k))/2.d0
C       s(k)=(0.5d0-alpha)*(b1(k)+b2(k))
C       t(k)=(2.d0*alpha-0.5d0)*(b1(k)+b2(k))
C     enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     sum over d1 x s  square
C     do i=0,max
C       x=dq*i
C       do j=0,max
C         y=dq*j
C         do k=1,3
C           xk(k)=x*d1(k)+y*s(k)
C         enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c         determine weights
C         iwght=4
C         if(i.eq.0)iwght=2
C         if(j.eq.0)iwght=2
C         if(i.eq.max)iwght=2
C         if(j.eq.max)iwght=2
C         if((i.eq.0).and.(j.eq.0))iwght=1
C         if((i.eq.0).and.(j.eq.max))iwght=1
C         if((i.eq.max).and.(j.eq.0))iwght=1
C         if((i.eq.max).and.(j.eq.max))iwght=1
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C         nktot=nktot+iwght
C         call cond(zener,xk,zconu,zcond,zconud,zcondu,ifail)
c
C         do in=0,ndiff
C           zresu(in) = zresu(in) + zconu(in)*iwght
C           zresd(in) = zresd(in) + zcond(in)*iwght
C           zresud(in) = zresud(in) + zconud(in)*iwght
C           zresdu(in) = zresdu(in) + zcondu(in)*iwght
C         enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C       enddo
C     enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     sum over d1 x t  triangle
C     do i=0,max
C       x=dq*i
C       do j=0,max
C         y=dq*j
C         do k=1,3
C           xk(k)=x*d1(k)+y*s(k)
C         enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c         determine weights
C         iwght=4
C         if(i.eq.0)iwght=2
C         if(j.eq.0)iwght=2
C         if(i.eq.max)iwght=2
C         if(j.eq.max)iwght=2
C         if((i.eq.0).and.(j.eq.0))iwght=1
C         if((i.eq.0).and.(j.eq.max))iwght=1
C         if((i.eq.max).and.(j.eq.0))iwght=1
C         if((i.eq.max).and.(j.eq.max))iwght=1
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C         nktot=nktot+iwght
C         call cond(zener,xk,zconu,zcond,zconud,zcondu,ifail)
C
C         do in=0,ndiff
C           zresu(in) = zresu(in) + zconu(in)*iwght
C           zresd(in) = zresd(in) + zcond(in)*iwght
C           zresud(in) = zresud(in) + zconud(in)*iwght
C           zresdu(in) = zresdu(in) + zcondu(in)*iwght
C         enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C     goto 5000
c     -----------------------------------------------------------------
c     =================================================================
5000  continue
      do in=0,ndiff
        zresu(in) = zresu(in)/dfloat(nktot)
        zresd(in) = zresd(in)/dfloat(nktot)
        zresud(in) = zresud(in)/dfloat(nktot)
        zresdu(in) = zresdu(in)/dfloat(nktot)
      enddo

      return
      end
c
c     *****************************************************************
c
      subroutine weight(i,j,max,irecip,iwght)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      if(irecip.eq.1)then
          iwght=8
          if(i.eq.j)iwght=4
          if(j.eq.0)iwght=4
          if(i.eq.max)iwght=4
          if(i.eq.0)iwght=1
          if(j.eq.max)iwght=1
          if((i.eq.max).and.(j.eq.0))iwght=2
      elseif(irecip.eq.2)then
          iwght=4
          if(i.eq.0)iwght=2
          if(j.eq.0)iwght=2
          if(i.eq.max)iwght=2
          if(j.eq.max)iwght=2
          if((i.eq.0).and.(j.eq.0))iwght=1
          if((i.eq.0).and.(j.eq.max))iwght=1
          if((i.eq.max).and.(j.eq.0))iwght=1
          if((i.eq.max).and.(j.eq.max))iwght=1
      elseif(irecip.eq.3)then
          iwght=12
          if(i.eq.j)iwght=6
          if(j.eq.0)iwght=6
          if(i.eq.max)iwght=6
          if(i.eq.0)iwght=1
          if(j.eq.max)iwght=2
          if((i.eq.max).and.(j.eq.0))iwght=3
      endif
      return
      end
c
c     *****************************************************************
c
      subroutine testk(x,y,b1,b2,irecip)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter (nplx=75,nslx=1,nspinx=10,nsubx=4,nsubatx=4)
      parameter (nmatx=nslx*nsubx*nspinx,natomx=nslx*nspinx*nsubatx)
      parameter (nlayx=nplx*nslx)
      parameter (nnx=4,npmax=100) !! upto nnx.th N.N's
      parameter (nfoldx=100)
      common/lbl/ pi, sq2, sq3
      common/data/nmat,nins,mlay,natom,nspin,nlay,nsub,numnn
      common/layer/ndiff
      common/fold/xfold(3,nfoldx),xshift(3,nfoldx),ifold(2,nfoldx),
     $  baib(3,3),nxfold,nfold,irecipa
      dimension b1(3),b2(3),xx1(nfoldx),xx2(nfoldx)
c     -----------------------------------------------------------------
      do i=1,nfold
        xx1(i)=0.d0
        xx2(i)=0.d0
      enddo
c     -----------------------------------------------------------------
      if(irecip.ge.3)then
        write(*,*)'reciprocal superlattice not cibic or rectangular'
        write(*,*)' this case has not been coded'
        stop
      endif
      if(irecipa.eq.0)then
        write(*,*)' TESTK : irecip = 0 ---- NOT CODED '
        stop
      endif
      if(irecipa.eq.1)goto 10
      if(irecipa.eq.2)goto 10
      if(irecipa.eq.3)goto 30
      if(irecipa.eq.4)goto 40
c     -----------------------------------------------------------------
C     CUBIC AND RECTANGULAR ATOMIC RECIPROCAL LATTICE
10    continue
c     first construct supercell reciprocal lattice cluster
c     for the vector J = i*b1 + j*b2
c     then           J = x1*ba1 + x2*ba2
      icnt=0
      do iff=1,nfold
        i=ifold(1,iff)
        j=ifold(2,iff)
        xk=0.5d0*x+i
        yk=0.5d0*y+j
        x1=baib(1,1)*xk+baib(1,2)*yk
        x2=baib(2,1)*xk+baib(2,2)*yk
        if((abs(x1).le.0.50000001d0).and.(abs(x2).le.0.50000001d0))then
c       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c       first check if this point is not related to any of the others by
c       a reciprocal lattice vector
          do iii=1,icnt
            dx=x1-xx1(iii)
            dy=x2-xx2(iii)
            ddx=int(abs(dx)+1.d-10)-abs(dx)
            ddy=int(abs(dy)+1.d-10)-abs(dy)
            if((abs(ddx).lt.1.d-10).and.(abs(ddy).lt.1.d-10))goto 55
          enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          icnt=icnt+1
          xx1(icnt)=x1
          xx2(icnt)=x2
          do k=1,3
            xfold(k,icnt)=(0.5d0*x+i)*b1(k)+(0.5d0*y+j)*b2(k)
            xshift(k,icnt)=i*b1(k)+j*b2(k)
          enddo
55        continue
        endif
      enddo
      nxfold=icnt
      if(nfold.gt.nfoldx)then
        write(*,*)' ERROR - FOLDING : nfold > nfoldx !! '
        stop
      endif
      goto 70
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C     HEXAGONAL ATOMIC RECIPROCAL LATTICE
30    continue
      write(*,*)' HEXAGONAL ATOMIC RECIPROCAL LATTICE NOT CODED'
      stop
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     CENTRED-RECTANGULAR ATOMIC RECIPROCAL LATTICE
40    continue
      write(*,*)' HEXAGONAL ATOMIC RECIPROCAL LATTICE NOT CODED'
      stop
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
70    continue

      return
      end
c
c     *****************************************************************
c
      subroutine folding(b1,b2,ba1,ba2)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter (nplx=75,nslx=1,nspinx=10,nsubx=4,nsubatx=4)
      parameter (nmatx=nslx*nsubx*nspinx,natomx=nslx*nspinx*nsubatx)
      parameter (nlayx=nplx*nslx)
      parameter (nnx=4,npmax=100) !! upto nnx.th N.N's
      parameter (nfoldx=100)
      common/lbl/ pi, sq2, sq3
      common/data/nmat,nins,mlay,natom,nspin,nlay,nsub,numnn
      common/layer/ndiff
      common/fold/xfold(3,nfoldx),xshift(3,nfoldx),ifold(2,nfoldx),
     $  baib(3,3),nxfold,nfold,irecipa
      dimension b1(3),b2(3),ba1(3),ba2(3),b12(3),ba12(3)
      dimension bba(3,3),bb(3,3),tmp1(3),tmp2(3)
      common/parallel/myid,numprocs

c     -----------------------------------------------------------------
c     check {b1,b2} and {ba1,ba2} in the same plane
      call cross(b1,b2,b12,3)
      call cross(ba1,ba2,ba12,3)
      if((abs(dot(b12,ba1,3)).gt.1.d-10).or.
     $  (abs(dot(b12,ba2,3)).gt.1.d-10))then
        write(*,*)'{b1,b2} not in same plane as {ba1,ba2}'
        stop
      endif
c     -----------------------------------------------------------------
c     construct the transformation matrix from {ba1,ba2} --> {b1,b2}
      do k=1,3
        bb(k,1)=b1(k)
        bb(k,2)=b2(k)
        bb(k,3)=b12(k)
        bba(k,1)=ba1(k)
        bba(k,2)=ba2(k)
        bba(k,3)=ba12(k)
      enddo
      call reinvers(bba,3,3)
      call remultiply(bba,bb,baib,3,3)
c     -----------------------------------------------------------------
      if(irecipa.eq.0)then
        write(*,*)' TESTK : irecip = 0 ---- NOT CODED '
        stop
      endif
      if(irecipa.eq.1)goto 10
      if(irecipa.eq.2)goto 10
      if(irecipa.eq.3)goto 30
      if(irecipa.eq.4)goto 40
c     -----------------------------------------------------------------
C     CUBIC AND RECTANGULAR ATOMIC RECIPROCAL LATTICE
10    continue
c     first construct supercell reciprocal lattice cluster
c     for the vector J = i*b1 + j*b2 + vsym(k)
c     then           J = x1*ba1 + x2*ba2
      icnt=0
      do isym=1,4
        if(isym.eq.1)then
          xsym=0.d0
          ysym=0.d0
        endif
        if(isym.eq.2)then
          xsym=0.5d0
          ysym=0.d0
        endif
        if(isym.eq.3)then
          xsym=0.d0
          ysym=0.5d0
        endif
        if(isym.eq.4)then
          xsym=0.5d0
          ysym=0.5d0
        endif
        do i=-nsub,nsub,1
          do j=-nsub,nsub,1
            x1=baib(1,1)*(i+xsym)+baib(1,2)*(j+ysym)
            x2=baib(2,1)*(i+xsym)+baib(2,2)*(j+ysym)

c           check consistency
            do l=1,3
              tmp1(l)=(i+xsym)*b1(l)+(j+ysym)*b2(l)
              tmp2(l)=x1*ba1(l)+x2*ba2(l)
              if(abs(tmp1(l)-tmp2(l)).gt.1.d-10)then
                write(*,*)' ERROR FOLDING ',tmp1,tmp2
                stop
              endif
            enddo

C           check that this vector lies inside the atomic BZ
            if((abs(x1).le.0.50000001d0).and.(abs(x2).le.0.50000001d0))
     $then

c           check if this point the same as any of the others
            do iii=1,icnt
              if((i.eq.ifold(1,iii)).and.(j.eq.ifold(2,iii)))goto 55
            enddo

              icnt=icnt+1
              ifold(1,icnt)=i
              ifold(2,icnt)=j
            endif
55        continue
          enddo
        enddo
      enddo
      nfold=icnt
      if(nfold.gt.nfoldx)then
        write(*,*)' ERROR - FOLDING : nfold > nfoldx !! '
        stop
      endif
      goto 70
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C     HEXAGONAL ATOMIC RECIPROCAL LATTICE
30    continue
      write(*,*)' HEXAGONAL ATOMIC RECIPROCAL LATTICE NOT CODED'
      stop
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     CENTRED-RECTANGULAR ATOMIC RECIPROCAL LATTICE
40    continue
      write(*,*)' HEXAGONAL ATOMIC RECIPROCAL LATTICE NOT CODED'
      stop
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
70    continue
      
      if(myid.eq.0)then
        write(*,*)
        write(*,*)'FOLDING VECTORS'
        do i=1,nfold
          write(*,'(i5,5x,2i5)')i,ifold(1,i),ifold(2,i)
        enddo
        write(*,*)
        write(*,*)
      endif
      return
      end
c
c     *****************************************************************
c
      subroutine cond(zener,xk,zconu,zcond,zconud,zcondu,ifail)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
C
C     Calculate the coupling at a given k// , via the det formula,
C     in the supercell representaion -- so that the SGF's are nmatx x nmatx
C
      parameter (nplx=75,nslx=1,nspinx=10,nsubx=4,nsubatx=4)
      parameter (nmatx=nslx*nsubx*nspinx,natomx=nslx*nspinx*nsubatx)
      parameter (nlayx=nplx*nslx)
      parameter (nnx=4,npmax=100) !! upto nnx.th N.N's
      parameter (nfoldx=100)
      common/lbl/ pi, sq2, sq3
      common/data/nmat,nins,mlay,natom,nspin,nlay,nsub,numnn
      common/data1/frat,ifrac
      common/layer/ndiff

      dimension zconu(0:nlayx),zcond(0:nlayx)
      dimension zconud(0:nlayx),zcondu(0:nlayx)
      dimension zglu(nmatx,nmatx),zgld(nmatx,nmatx)
      dimension zgru(nmatx,nmatx),zgrd(nmatx,nmatx)
      dimension zu(nmatx,nmatx),zt(nmatx,nmatx)
      dimension ztdag(nmatx,nmatx)
      dimension xk(3)
      dimension zfoo(nmatx,nmatx)

c     =================================================================
C     CALCULATE SGFs IN ATOMIC BASIS
c     =================================================================
C     DO THIS IF LH AND RH LEADS ARE THE SAME
C     call green(zener,xk,+1,"LH",zglu,zgru,ifail)   !! LH UP
C     if(ifail.ne.0)return
C     call green(zener,xk,-1,"LH",zgld,zgrd,ifail)   !! LH DOWN
C     if(ifail.ne.0)return

c     -----------------------------------------------------------------
C     DO THIS IF LH AND RH LEADS DIFFER
      call green(zener,xk,+1,"LH",zglu,zfoo,ifail)   !! LH UP
      if(ifail.ne.0)return
      call green(zener,xk,-1,"LH",zgld,zfoo,ifail)   !! LH DOWN
      if(ifail.ne.0)return

      call green(zener,xk,+1,"RH",zfoo,zgru,ifail)   !! RH UP
      if(ifail.ne.0)return
      call green(zener,xk,-1,"RH",zfoo,zgrd,ifail)   !! RH DOWN
      if(ifail.ne.0)return
c     -----------------------------------------------------------------

c     =================================================================
C     ALTERNATIVELY DO EVERYTHING IN SUPERCELL
c     note we do not normally keep surfacenewcell as a subroutine as 
c     it adds to the size of the code
c     To create surfacenewcell make changes to surfacenew as indicated
c     ifail=0
c     call hamil(zt,xk,1,2,+1,nmat,nmatx)
c     call hamil(zu,xk,2,2,+1,nmat,nmatx)
c     call surfacedecim(zu,zt,zener,zglu,zgru,errmax,nmat,nmatx)
c     call hamil(zt,xk,1,2,-1,nmat,nmatx)
c     call hamil(zu,xk,2,2,-1,nmat,nmatx)
c     call surfacedecim(zu,zt,zener,zgld,zgrd,errmax,nmat,nmatx)
c     =================================================================
c
C     adlayer the LH & RH mlay substrate layers (ie interface layers)
      mlayl=1
      mlayr=mlay+nins+1


      do ill=3,2+mlayl
        call hamil(zt,xk,ill-1,ill,+1,nmat,nmatx)
        call hamil(zu,xk,ill,ill,+1,nmat,nmatx)
        call adlayer1(zglu,zu,zt,zfoo,zener,nmat,nmatx)

        call hamil(zt,xk,ill-1,ill,-1,nmat,nmatx)
        call hamil(zu,xk,ill,ill,-1,nmat,nmatx)
        call adlayer1(zgld,zu,zt,zfoo,zener,nmat,nmatx)
      enddo

      do ill=nlay-2,nlay-mlayr-1,-1
        call hamil(zt,xk,ill,ill+1,+1,nmat,nmatx)
        call hamil(zu,xk,ill,ill,+1,nmat,nmatx)
        do ir=1,nmat
          do is=1,nmat
            ztdag(ir,is)=dconjg(zt(is,ir))
          enddo
        enddo
        call adlayer1(zgru,zu,ztdag,zfoo,zener,nmat,nmatx)

        call hamil(zt,xk,ill,ill+1,-1,nmat,nmatx)
        call hamil(zu,xk,ill,ill,-1,nmat,nmatx)
        do ir=1,nmat
          do is=1,nmat
            ztdag(ir,is)=dconjg(zt(is,ir))
          enddo
        enddo
        call adlayer1(zgrd,zu,ztdag,zfoo,zener,nmat,nmatx)
      enddo


c     =================================================================
C     CALCULATE GFs IN SPACER AND THE CONDUCTANCE
c     =================================================================
      icnt=0
      do ill=mlayl+3,nlay-mlayr-2
        icnt=icnt+1

c       adlayer LH  ----   zgl
        call hamil(zt,xk,ill-1,ill,+1,nmat,nmatx)
        call hamil(zu,xk,ill,ill,+1,nmat,nmatx)
        call adlayer1(zglu,zu,zt,zfoo,zener,nmat,nmatx)
        call hamil(zt,xk,ill-1,ill,-1,nmat,nmatx)
        call hamil(zu,xk,ill,ill,-1,nmat,nmatx)
        call adlayer1(zgld,zu,zt,zfoo,zener,nmat,nmatx)


c       SPIN UP
        call hamil(zt,xk,nlay-mlayr-2,nlay-mlayr-1,+1,nmat,nmatx)
        call kubo(zglu,zgru,zt,zkon)
        zconu(icnt)=zkon
c
c       SPIN DOWN
        call hamil(zt,xk,nlay-mlayr-2,nlay-mlayr-1,-1,nmat,nmatx)
        call kubo(zgld,zgrd,zt,zkon)
        zcond(icnt)=zkon
c
c       SPIN UP-DOWN
        call hamil(zt,xk,nlay-mlayr-2,nlay-mlayr-1,-1,nmat,nmatx)
        call kubo(zglu,zgrd,zt,zkon)
        zconud(icnt)=zkon
c
c       SPIN DOWN-UP
        call hamil(zt,xk,nlay-mlayr-2,nlay-mlayr-1,+1,nmat,nmatx)
        call kubo(zgld,zgru,zt,zkon)
        zcondu(icnt)=zkon

      enddo

      return
      end
c
c     *****************************************************************
c
      subroutine kubo(zgsl,zgsr,zt,zkon)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter (nplx=75,nslx=1,nspinx=10,nsubx=4,nsubatx=4)
      parameter (nmatx=nslx*nsubx*nspinx,natomx=nslx*nspinx*nsubatx)
      parameter (nlayx=nplx*nslx)
      parameter (nnx=4,npmax=100) !! upto nnx.th N.N's
      parameter (nfoldx=100)
      common/lbl/ pi, sq2, sq3
      common/data/nmat,nins,mlay,natom,nspin,nlay,nsub,numnn
      common/data1/frat,ifrac
      common/layer/ndiff
      dimension zt(nmatx,nmatx)
      dimension zgsl(nmatx,nmatx),zgsr(nmatx,nmatx)

      common/shareus/ztdag(nmatx,nmatx),zglt(nmatx,nmatx),
     $  zgrt(nmatx,nmatx),zltrt(nmatx,nmatx),zrtlt(nmatx,nmatx),
     $    ztltrt(nmatx,nmatx),ztrtlt(nmatx,nmatx)

c
c     this routine calculates the kubo formula for SGFs gl,gr
c     and hopping t
c
      zi=dcmplx(0.d0,1.d0)
      do ir=1,nmat
        do is=1,nmat
          ztdag(ir,is)=dconjg(zt(is,ir))
        enddo
      enddo
c
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     construct total greens functions

      call multiply(zgsl,zt,zglt,nmat,nmatx)
      call multiply(zgsr,ztdag,zgrt,nmat,nmatx)

      call multiply(zgrt,zglt,zrtlt,nmat,nmatx)

      do ir=1,nmat
        do is=1,nmat
          zrtlt(ir,is)=-zrtlt(ir,is)
        enddo
        zrtlt(ir,ir)=1.d0+zrtlt(ir,ir)
      enddo

      call invers(zrtlt,nmat,nmatx)  !!! (1 - gr t+ gl t)^{-1}

      call multiply(zt,zrtlt,ztrtlt,nmat,nmatx)

      do ir=1,nmat
        do is=1,nmat
          ztltrt(ir,is)=dconjg(ztrtlt(is,ir))
        enddo
      enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     construct the Imaginary parts
      do ir=1,nmat
        do is=1,nmat
          zgrt(ir,is)=zgsr(ir,is)-dconjg(zgsr(is,ir))
          zglt(ir,is)=zgsl(ir,is)-dconjg(zgsl(is,ir))
        enddo
      enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     construct kubo formula
      call multiply(ztrtlt,zgrt,zrtlt,nmat,nmatx)
      call multiply(ztltrt,zglt,zltrt,nmat,nmatx)

      call multiply(zrtlt,zltrt,zgrt,nmat,nmatx)

      zkon=0.d0
      do ir=1,nmat
        zkon=zkon-zgrt(ir,ir)
      enddo
c
      return
      end
c
c     *****************************************************************
c
      subroutine green(zener,xk,ispin,side,zgl,zgr,ifail)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
C
C     Calculate the SGF at a given k// ,
C     in the supercell representaion -- so that the SGF's are nmatx x nmatx
C
C     if side="LH" we calculate the LH SGF. if side="RH" the RH SGF
C     if ispin=-1 minority, +1 majority
C
      parameter (nplx=75,nslx=1,nspinx=10,nsubx=4,nsubatx=4)
      parameter (nmatx=nslx*nsubx*nspinx,natomx=nslx*nspinx*nsubatx)
      parameter (nlayx=nplx*nslx)
      parameter (nnx=4,npmax=100) !! upto nnx.th N.N's
      parameter (nfoldx=100)
      parameter (numatx=8)       !! No of atom types
      common/lbl/ pi, sq2, sq3
      common/data/nmat,nins,mlay,natom,nspin,nlay,nsub,numnn
      common/data1/frat,ifrac
      common/layer/ndiff

      dimension zgl(nmatx,nmatx)
      dimension zgr(nmatx,nmatx)

      common/shareus/ztmp1(nmatx,nmatx),ztmp2(nmatx,nmatx),
     $  zu(nmatx,nmatx),zt(nmatx,nmatx),ztdag(nmatx,nmatx),
     $    zgltmp(nmatx,nmatx),zgrtmp(nmatx,nmatx)


      dimension xk(3)

      dimension zuat(natomx,natomx),ztat(natomx,natomx)
      dimension zglat(natomx,natomx)
      dimension zgrat(natomx,natomx)

      common/ssstructat/aa1(3),aa2(3),aa3(3,nlayx),vsubat(3,nsubx,nlayx)
     $,itypeat(nsubx,nlayx),imapl(nsubx),imapr(nsubx),nsubat
      common/fold/xfold(3,nfoldx),xshift(3,nfoldx),ifold(2,nfoldx),
     $  baib(3,3),nxfold,nfold,irecipa
      common/ssstruct/a1(3),a2(3),a3(3,nlayx),vsub(3,nsubx,nlayx),
     $  dnn(nnx,nlayx,nlayx),ddnn(numatx,numatx,nnx),itype(nsubx,nlayx)
      dimension xkat(3),ds(3),vtmp(3)

      dimension imap(nsubx)
      character*2 side
c     -----------------------------------------------------------------
C     calculate the LH & RH SGF in the atomic basis, and convert to the 
C     supercell basis.


      if(side.eq."LH")then
        ilay=1
        do i=1,nsub
          imap(i)=imapl(i)
        enddo
      elseif(side.eq."RH")then
        ilay=nlay-1
        do i=1,nsub
          imap(i)=imapr(i)
        enddo
      endif


      do ir=1,nmat
        do is=1,nmat
          zgltmp(ir,is)=0.d0
          zgrtmp(ir,is)=0.d0
        enddo
      enddo
      do iff=1,nxfold
        do k=1,3
          xkat(k)=xfold(k,iff)
        enddo

c       define atomic hamiltonian elements in lead   -  u1, t1
        call hamil(ztat,xkat,ilay,ilay+1,ispin,natom,natomx)
        call hamil(zuat,xkat,ilay+1,ilay+1,ispin,natom,natomx)

        call surfacedecim(zuat,ztat,zener,zglat,zgrat,errmax,natom,
     $natomx)
        if(errmax.gt.1.d-3)then   !!! zt has a near zero eigenvalue
          write(*,*)"ERROR surfacedecim:  errmax =",errmax
          ifail=1
        endif

c       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c       now load elements into zgl and zgr
        do isub=1,nsub           !!! sublattice elements
          do jsub=1,nsub
            do l=1,3
              ds(l)=vsub(l,isub,ilay)-vsub(l,jsub,ilay) -
     $          (vsubat(l,imap(isub),ilay)-vsubat(l,imap(jsub),ilay))
            enddo
            arg1=dot(xkat,ds,3)
            zarg1=dcmplx(0.d0,arg1)
            do l=1,3
              vtmp(l)=xshift(l,iff)
            enddo

            do isp=1,nspin       !!! orbital elements
              do jsp=1,nspin
                iii=(isub-1)*nspin+isp
                jjj=(jsub-1)*nspin+jsp

                ii=(imap(isub)-1)*nspin+isp
                jj=(imap(jsub)-1)*nspin+jsp

                zgltmp(iii,jjj)=zgltmp(iii,jjj)+
     $        zglat(ii,jj)*cdexp(zarg1)*dfloat(nsubat)/dfloat(nsub)
                zgrtmp(iii,jjj)=zgrtmp(iii,jjj)+
     $        zgrat(ii,jj)*cdexp(zarg1)*dfloat(nsubat)/dfloat(nsub)

              enddo
            enddo
          enddo
        enddo
      enddo
      do ir=1,nmat
        do is=1,nmat
          zgl(ir,is)=zgltmp(ir,is)
          zgr(ir,is)=zgrtmp(ir,is)
        enddo
      enddo

      return
      end
c
c     *****************************************************************
c
      subroutine adlayer1(zgl,zu,zt,zwrk,zener,n,nx)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)

      dimension zgl(nx,nx),zu(nx,nx),zt(nx,nx)
      dimension zwrk(nx,nx)

c     -----------------------------------------------------------------
c     adlayer ontop of gl
      call multiply(zgl,zt,zwrk,n,nx)
      call hconjugate(zt,n,nx)
      call multiply(zt,zwrk,zgl,n,nx)
      call hconjugate(zt,n,nx)

      do ir=1,n
        do is=1,ir-1
          zgl(ir,is)=-zu(ir,is)-zgl(ir,is)
          zgl(is,ir)=-zu(is,ir)-zgl(is,ir)
        enddo
        zgl(ir,ir)=zener-zu(ir,ir)-zgl(ir,ir)
      enddo
c     find inverse to zginv
      call invers(zgl,n,nx)
c     -----------------------------------------------------------------

      return
      end
c
c     *****************************************************************
c
      subroutine hconjugate(zmat,n,nx)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      dimension zmat(nx,nx)
      do ir=1,n
        do is=1,ir-1
          ztmp=zmat(is,ir)
          zmat(is,ir)=dconjg(zmat(ir,is))
          zmat(ir,is)=dconjg(ztmp)
        enddo
        zmat(ir,ir)=dconjg(zmat(ir,ir))
      enddo

      return
      end
c
c     *****************************************************************
c     *****************************************************************
C     HAMILTONIAN ROUTINES -- GENERAL STRUCTURE
c     *****************************************************************
c     *****************************************************************
      subroutine hamil(zt,xk,n1tmp,n2tmp,ispin,ntmp,nx)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)

      parameter (nplx=75,nslx=1,nspinx=10,nsubx=4,nsubatx=4)
      parameter (nmatx=nslx*nsubx*nspinx,natomx=nslx*nspinx*nsubatx)
      parameter (nlayx=nplx*nslx)
      parameter (nnx=4,npmax=100) !! upto nnx.th N.N's
c                                    npmax - max No of points per plane
      parameter (numatx=8)       !! No of atom types
      common/lbl/ pi, sq2, sq3
      common/data/nmat,nins,mlay,natom,nspin,nlay,nsub,numnn
      common/data1/frat,ifrac
      common/layer/ndiff
      dimension xk(3)
      dimension zt(nx,nx)
      dimension zh(nspinx,nspinx)
      common/shift/xmgo_shift,ag_shift,xinas_shift

      common/ssstructat/aa1(3),aa2(3),aa3(3,nlayx),vsubat(3,nsubx,nlayx)
     $,itypeat(nsubx,nlayx),imapl(nsubx),imapr(nsubx),nsubat
c     -----------------------------------------------------------------
c     compute the hamiltonian matrices zu and zt

c     calculate for given layers n1,n2 ; sublattice points isub,jsub ;
c     orbital indices isp,jsp the Hamiltonian element
c     <n1 ; isub ; isp | H | n2 ; jsub ; jsp>

      n=ntmp              !!! these lines to stop ftnchek moaning
      n1=n1tmp
      n2=n2tmp

      if(n.eq.natom)ncell=nsubat
      if(n.eq.nmat)ncell=nsub

c     compute U or T
      do isub=1,ncell
        do jsub=1,ncell
          call helement(n1,n2,isub,jsub,xk,zh,ispin,n)
          do isp=1,nspin       !!! orbital elements
            do jsp=1,nspin
              iii=(isub-1)*nspin+isp
              jjj=(jsub-1)*nspin+jsp
              zt(iii,jjj)=zh(isp,jsp)
            enddo
          enddo
        enddo
      enddo

c     -----------------------------------------------------------------

      return
      end

c
c     *****************************************************************
c
      subroutine helement(n1,n2,isub,jsub,xk,zh,ispin,ndim)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter (nplx=75,nslx=1,nspinx=10,nsubx=4,nsubatx=4)
      parameter (nmatx=nslx*nsubx*nspinx,natomx=nslx*nspinx*nsubatx)
      parameter (nlayx=nplx*nslx)
      parameter (nnx=4,npmax=100) !! upto nnx.th N.N's
c                                    npmax - max No of points per plane
      parameter (numatx=8)       !! No of atom types
      common/lbl/ pi, sq2, sq3
      common/data/nmat,nins,mlay,natom,nspin,nlay,nsub,numnn
      common/data1/frat,ifrac
      common/layer/ndiff
      dimension d(3),c(3),xk(3),dpar(3)
      dimension zh(nspinx,nspinx),zt(nspinx,nspinx)
      dimension rt(nspinx,nspinx)
      common/par1/s0(numatx,-1:1),p0(numatx,-1:1),d0e(numatx,-1:1),
     $  d0t(numatx,-1:1)
      common/pard/d05(numatx,-1:1),d06(numatx,-1:1),d07(numatx,-1:1),
     $d08(numatx,-1:1),d09(numatx,-1:1)
      common/par2/sssint(numatx,numatx,nnx,-1:1),
     $  spsint(numatx,numatx,nnx,-1:1),
     $  ppsint(numatx,numatx,nnx,-1:1),pppint(numatx,numatx,nnx,-1:1),
     $  sdsint(numatx,numatx,nnx,-1:1),pdsint(numatx,numatx,nnx,-1:1),
     $  pdpint(numatx,numatx,nnx,-1:1),ddsint(numatx,numatx,nnx,-1:1),
     $  ddpint(numatx,numatx,nnx,-1:1),dddint(numatx,numatx,nnx,-1:1)
      common/par3/x0(numatx,-1:1), xxsint(numatx,numatx,nnx,-1:1),
     $  xssint(numatx,numatx,nnx,-1:1),xpsint(numatx,numatx,nnx,-1:1),
     $  xdsint(numatx,numatx,nnx,-1:1)
      common/shift/xmgo_shift,ag_shift,xinas_shift
      common/ssstruct/a1(3),a2(3),a3(3,nlayx),vsub(3,nsubx,nlayx),
     $  dnn(nnx,nlayx,nlayx),ddnn(numatx,numatx,nnx),itype(nsubx,nlayx)
      common/ssstructat/aa1(3),aa2(3),aa3(3,nlayx),vsubat(3,nsubx,nlayx)
     $,itypeat(nsubx,nlayx),imapl(nsubx),imapr(nsubx),nsubat


      if(ndim.eq.natom)then
        ind1=itypeat(isub,n1)
        ind2=itypeat(jsub,n2)
      endif
      if(ndim.eq.nmat)then
        ind1=itype(isub,n1)
        ind2=itype(jsub,n2)
      endif


      ijsub=jsub-isub
      do ir=1,nspin
        do is=1,nspin
          zh(ir,is)=0.d0
          zt(ir,is)=0.d0
        enddo
      enddo
      
C
C     Calculate 
C     H_{n1,n2}(isub,jsub) =
C          \sum{d_{//}} exp(i k_{//} d_{//}) <isub,n1|H|jsub,n2,d_{//}>
C
C     To do this we need to calculate 
C     D = (vsub(isub,n1)+a3(n1)) - (vsub(jsub,n2)+a3(n2)+d//)
C     We need to sum over d_{//}

      ddmax=0.d0
      do inn=1,numnn
        ddmax=max(ddmax,ddnn(ind1,ind2,inn)) ! sometimes ddnn(ind1,ind2,numnn)=0
      enddo
      do i1=-numnn-1,numnn+1,1
        do i2=-numnn-1,numnn+1,1

C         First construct d_{//}
          do k=1,3
            if(ndim.eq.natom)then
              dpar(k)=i1*aa1(k)+i2*aa2(k)+aa3(k,n2)-aa3(k,n1)
              d(k)=vsubat(k,jsub,n2)-vsubat(k,isub,n1)+dpar(k)
            endif
            if(ndim.eq.nmat)then
              dpar(k)=i1*a1(k)+i2*a2(k)+a3(k,n2)-a3(k,n1)
              d(k)=vsub(k,jsub,n2)-vsub(k,isub,n1)+dpar(k)
            endif
          enddo
          dd=sqrt(d(1)**2 + d(2)**2 + d(3)**2)


          if(dd.gt.(ddmax+1.d-8))goto 999   !! this is not a NN
          if(abs(dd).lt.1.d-8)then
            nn=0
            goto 888
          endif
          do inn=1,numnn
            if(abs(dd-ddnn(ind1,ind2,inn)).lt.1.d-8)then
              nn=inn
              goto 888
            endif
          enddo
CCC       write(*,*)' ERROR HELEMENT : dd NOT a NN'
CCC       stop
          go to 999
888       continue

C         calculate the cosine angles :
          if(dd.gt.1.d-8)then
            do k=1,3
              c(k)=d(k)/dd
            enddo
          endif

          do ir=1,nspin
            do is=1,nspin
              rt(ir,is)=0.d0
            enddo
          enddo

c       ?????????????????????????????????????????????????????????????????
c       ?????????????????????????????????????????????????????????????????
CCCCCCC   dk=dot(xk,d,3)
          dk=dot(xk,dpar,3)   !!! need this alteration with folding
c       ?????????????????????????????????????????????????????????????????
c       ?????????????????????????????????????????????????????????????????
          zexdk=dcmplx(cos(dk),sin(dk))

c         find the hopping for this nn.th N.N in direction d
c       ?????????????????????????????????????????????????????????????????
c       ?????????????????????????????????????????????????????????????????
c         ensure that this routine gives U for d = 0
CCCCC     call eint11(nn,ind1,zt)
c       ?????????????????????????????????????????????????????????????????
c       ?????????????????????????????????????????????????????????????????
          if(nn.eq.0)then
            rt(1,1)=s0(ind1,ispin)
            rt(2,2)=p0(ind1,ispin)
            rt(3,3)=p0(ind1,ispin)
            rt(4,4)=p0(ind1,ispin)
            rt(5,5)=d05(ind1,ispin)
            rt(6,6)=d06(ind1,ispin)
            rt(7,7)=d07(ind1,ispin)
            rt(8,8)=d08(ind1,ispin)
            rt(9,9)=d09(ind1,ispin)
            rt(10,10)=x0(ind1,ispin)
  
CCC         rt(5,5)=d0e(ind1,ispin)
CCC         rt(6,6)=d0t(ind1,ispin)
CCC         rt(7,7)=d0t(ind1,ispin)
CCC         rt(8,8)=d0t(ind1,ispin)
CCC         rt(9,9)=d0e(ind1,ispin)
          else
              s0s1s=sssint(ind1,ind2,nn,ispin)
              s0p1s=spsint(ind1,ind2,nn,ispin)
              p0p1s=ppsint(ind1,ind2,nn,ispin)
              p0p1p=pppint(ind1,ind2,nn,ispin)
              s0d1s=sdsint(ind1,ind2,nn,ispin)
              p0d1s=pdsint(ind1,ind2,nn,ispin)
              p0d1p=pdpint(ind1,ind2,nn,ispin)
              d0d1s=ddsint(ind1,ind2,nn,ispin)
              d0d1p=ddpint(ind1,ind2,nn,ispin)
              d0d1d=dddint(ind1,ind2,nn,ispin)
              x0x1s=xxsint(ind1,ind2,nn,ispin)
              x0s1s=xssint(ind1,ind2,nn,ispin)
              x0p1s=xpsint(ind1,ind2,nn,ispin)
              x0d1s=xdsint(ind1,ind2,nn,ispin)

              s1s0s=sssint(ind2,ind1,nn,ispin)
              s1p0s=spsint(ind2,ind1,nn,ispin)
              p1p0s=ppsint(ind2,ind1,nn,ispin)
              p1p0p=pppint(ind2,ind1,nn,ispin)
              s1d0s=sdsint(ind2,ind1,nn,ispin)
              p1d0s=pdsint(ind2,ind1,nn,ispin)
              p1d0p=pdpint(ind2,ind1,nn,ispin)
              d1d0s=ddsint(ind2,ind1,nn,ispin)
              d1d0p=ddpint(ind2,ind1,nn,ispin)
              d1d0d=dddint(ind2,ind1,nn,ispin)
              x1x0s=xxsint(ind2,ind1,nn,ispin)
              x1s0s=xssint(ind2,ind1,nn,ispin)
              x1p0s=xpsint(ind2,ind1,nn,ispin)
              x1d0s=xdsint(ind2,ind1,nn,ispin)

              call slater_koster(s0s1s,s0p1s,s0d1s,p0p1s,p0p1p,p0d1s,
     $p0d1p,d0d1s,d0d1p,d0d1d,x0x1s,x0s1s,x0p1s,x0d1s,s1s0s,s1p0s,
     $s1d0s,p1p0s,p1p0p,p1d0s,p1d0p,d1d0s,d1d0p,d1d0d,x1x0s,x1s0s,x1p0s,
     $x1d0s,c,rt)
          endif
          do ir=1,nspin
            do is=1,nspin
              zt(ir,is)=rt(ir,is)
            enddo
          enddo

          do ir=1,nspin
            do is=1,nspin
              zh(ir,is)=zh(ir,is)+zt(ir,is)*zexdk
            enddo
          enddo

999       continue
        enddo
      enddo

      return
      end
c
c     *****************************************************************
c
      subroutine prestructij(ilay,jlay,dist)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter (nplx=75,nslx=1,nspinx=10,nsubx=4,nsubatx=4)
      parameter (nnx=4,npmax=100)   !!! upto nnx th N.N
c                                    npmax - max No of points per plane
      parameter (nmatx=nslx*nsubx*nspinx,natomx=nslx*nspinx*nsubatx)
      parameter (nlayx=nplx*nslx)
      parameter (lattx=(2*nnx+1)*(2*nnx+1)*(2*nnx+1)*nsubx*nsubx)
      parameter (numatx=8)       !! No of atom types
c     find the NN distances from layer ilay to the layer jlay.
      dimension v(3)
      dimension dist(nnx)
      dimension disttmp(nnx),iorder(nnx)

      common/data/nmat,nins,mlay,natom,nspin,nlay,nsub,numnn

      common/ssstruct/a1(3),a2(3),a3(3,nlayx),vsub(3,nsubx,nlayx),
     $  dnn(nnx,nlayx,nlayx),ddnn(numatx,numatx,nnx),itype(nsubx,nlayx)
c     -----------------------------------------------------------------
c     INPUT : {a1,a2,a3} - primitive lattice vectors
c             nsub       - number of sub-lattice vectors
c             vsub       - the set of nsub sublattice vectors
c                          {vsub(1),...,vsub(nsub)}
c             numnn      - upto order numnn.th N.N interactions
c
c     OUTPUT : dist(i)   - the distance from the origin to the ith N.N
c     -----------------------------------------------------------------
C     Construct vectors from (isub,ilay) --> (jsub,jlay)

      do i=1,numnn
        disttmp(i)=1.d10
      enddo

      do i=-numnn,numnn,1
        do j=-numnn,numnn,1
          do isub=1,nsub
            do jsub=1,nsub
              do ll=1,3
                v(ll)=(vsub(ll,isub,ilay)+a3(ll,ilay))-
     $             (vsub(ll,jsub,jlay)+a3(ll,jlay)+i*a1(ll)+j*a2(ll))
              enddo
              dv=sqrt(v(1)**2 + v(2)**2 + v(3)**2)

C   check dv not already in disttmp and replace the largest element it is smaller than
              if(dv.gt.1.d-8)then
                do k=1,numnn
                  if(abs(dv-disttmp(k)).lt.1.d-10)goto 123
                enddo
                kmax=1
                do k=2,numnn
                  if(disttmp(k).gt.disttmp(kmax))kmax=k
                enddo
                if(dv.lt.disttmp(kmax))disttmp(kmax)=dv
              endif
123           continue
            enddo
          enddo
        enddo
      enddo
c     -----------------------------------------------------------------
c     sort into distance from origin
      call sort(disttmp,iorder,numnn)

      do k=1,numnn
        dist(k)=disttmp(iorder(k))
      enddo
c     -----------------------------------------------------------------
c
      return
      end
c
c     *****************************************************************
c
      subroutine prestructijnew(ilay,jlay,disttmp)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter (nplx=75,nslx=1,nspinx=10,nsubx=4,nsubatx=4)
      parameter (nnx=4,npmax=100)   !!! upto nnx th N.N
c                                    npmax - max No of points per plane
      parameter (nmatx=nslx*nsubx*nspinx,natomx=nslx*nspinx*nsubatx)
      parameter (nlayx=nplx*nslx)
      parameter (lattx=(2*nnx+1)*(2*nnx+1)*(2*nnx+1)*nsubx*nsubx)
      parameter (numatx=8)       !! No of atom types
c     find the NN distances from layer ilay to the layer jlay.
      dimension v(3)
      dimension dist(nnx)
      dimension disttmp(numatx,numatx,nnx),iorder(nnx)

      common/data/nmat,nins,mlay,natom,nspin,nlay,nsub,numnn

      common/ssstruct/a1(3),a2(3),a3(3,nlayx),vsub(3,nsubx,nlayx),
     $  dnn(nnx,nlayx,nlayx),ddnn(numatx,numatx,nnx),itype(nsubx,nlayx)
      common/elements/iatfe(5),numat
c     -----------------------------------------------------------------
c     INPUT : {a1,a2,a3} - primitive lattice vectors
c             nsub       - number of sub-lattice vectors
c             vsub       - the set of nsub sublattice vectors
c                          {vsub(1),...,vsub(nsub)}
c             numnn      - upto order numnn.th N.N interactions
c
c     OUTPUT : dist(i)   - the distance from the origin to the ith N.N
c     -----------------------------------------------------------------
C     Construct vectors from (isub,ilay) --> (jsub,jlay)

      do i=-numnn,numnn,1
        do j=-numnn,numnn,1
          do isub=1,nsub
            do jsub=1,nsub
              do ll=1,3
                v(ll)=(vsub(ll,isub,ilay)+a3(ll,ilay))-
     $             (vsub(ll,jsub,jlay)+a3(ll,jlay)+i*a1(ll)+j*a2(ll))
              enddo
              dv=sqrt(v(1)**2 + v(2)**2 + v(3)**2)

C   check dv not already in disttmp and replace the largest element it is smaller than
              iat=itype(isub,ilay)
              jat=itype(jsub,jlay)
              if(dv.gt.1.d-8)then
                do inn=1,numnn
                  if(abs(dv-disttmp(iat,jat,inn)).lt.1.d-10)goto 123
                enddo
                imax=1
                do inn=2,numnn
                  if(disttmp(iat,jat,inn).gt.disttmp(iat,jat,imax))
     $imax=inn
                enddo
                if(dv.lt.disttmp(iat,jat,imax))disttmp(iat,jat,imax)=dv
              endif
123           continue
            enddo
          enddo
        enddo
      enddo
c     -----------------------------------------------------------------
c     sort into distance from origin
      do iat=1,numat
        do jat=1,numat
          do inn=1,numnn
            dist(inn)=disttmp(iat,jat,inn)
          enddo
          call sort(dist,iorder,numnn)
          do inn=1,numnn
            disttmp(iat,jat,inn)=dist(iorder(inn))
          enddo
        enddo
      enddo
c     -----------------------------------------------------------------
c
      return
      end
c
c     *****************************************************************
c
      subroutine prestruct(vsub,a1,a2,a3,numnn,nsub,dist)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter (nplx=75,nslx=1,nspinx=10,nsubx=4,nsubatx=4)
      parameter (nnx=4,npmax=100)   !!! upto nnx th N.N
c                                    npmax - max No of points per plane
      parameter (nmatx=nslx*nsubx*nspinx,natomx=nslx*nspinx*nsubatx)
      parameter (nlayx=nplx*nslx)
      parameter (lattx=(2*nnx+1)*(2*nnx+1)*(2*nnx+1)*nsubx)
c     construct a lattice from the primitive lattice vectors {a1,a2,a3}
c     and order according to distance from origin
      dimension r(lattx,3),a1(3),a2(3),a3(3),v(3),d(lattx),iorder(lattx)
      dimension rtmp(lattx,3),dtmp(lattx),dist(0:nnx+1)
      dimension vsub(3,nsubx)
c     -----------------------------------------------------------------
c     INPUT : {a1,a2,a3} - primitive lattice vectors
c             nsub       - number of sub-lattice vectors
c             vsub       - the set of nsub sublattice vectors
c                          {vsub(1),...,vsub(nsub)}
c             numnn      - upto order numnn.th N.N interactions
c
c     OUTPUT : dist(i)   - the distance from the origin to the ith N.N
c     -----------------------------------------------------------------
c     now construct lattice
      icnt=0
      do i=-numnn,numnn,1
        do j=-numnn,numnn,1
          do k=-numnn,numnn,1
            do ms=1,nsub
              icnt=icnt+1
              iorder(icnt)=icnt
              do ll=1,3
                v(ll)=i*a1(ll)+j*a2(ll)+k*a3(ll)+vsub(ll,ms)
                r(icnt,ll)=v(ll)
              enddo
              d(icnt)=sqrt(v(1)**2 + v(2)**2 + v(3)**2)
            enddo
          enddo
        enddo
      enddo
      ntot=icnt
      if(ntot.gt.lattx)then
        write(*,*)' ERROR - STRUCT : ntot > lattx !! '
        stop
      endif
c     -----------------------------------------------------------------
c     sort into distance from origin
      call sort(d,iorder,ntot)
      do i=1,ntot
        do j=1,3
          rtmp(i,j)=r(iorder(i),j)
        enddo
        dtmp(i)=d(iorder(i))
      enddo
c     -----------------------------------------------------------------
c     count up number of points with same distance from origin
      d0=dtmp(1) !!! = 0
      id=1       !!! id th most distant set from origin
      dist(id-1)=dtmp(1)
      icnt=1
      do i=2,ntot
        if(abs(d0-dtmp(i)).gt.1.d-8)then
          d0=dtmp(i)
          id=id+1
          dist(id-1)=dtmp(i)
        endif
        if(id-1.gt.numnn)goto 100
        icnt=i
      enddo
100   continue
      return
      end
c
c     *****************************************************************
c
      subroutine recip(a1,a2,b1,b2,irecip)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter (nplx=75,nslx=1,nspinx=10,nsubx=4,nsubatx=4)
      parameter (nmatx=nslx*nsubx*nspinx,natomx=nslx*nspinx*nsubatx)
      parameter (nlayx=nplx*nslx)
      parameter (nnx=4,npmax=100) !! upto nnx.th N.N's
      parameter (lattx=(2*nnx+1)*(2*nnx+1)*(2*nnx+1))
      common/lbl/ pi, sq2, sq3
      common/data/nmat,nins,mlay,natom,nspin,nlay,nsub,numnn
      common/data1/frat,ifrac
      common/layer/ndiff
      dimension an1(3),an2(3),d(3),v(3),b1(3),b2(3),xn(3)
      dimension r(lattx,3),a1(3),a2(3),a3(3),dv(lattx),iorder(lattx)
      dimension rtmp(lattx,3),dtmp(lattx)
      common/parallel/myid,numprocs
 
c     initialise
      do l=1,3
        b1(l)=0.d0
        b2(l)=0.d0
      enddo
      irecip=0

c     now construct reciprocal lattice perp'r to xn
      do k=1,3
        an1(k)=a1(k)
        an2(k)=a2(k)
      enddo
c     -----------------------------------------------------------------
c     construct the third (normalised) lattice vector
      call cross(an1,an2,xn,3)
      xnorm=dot(xn,xn,3)
      do i=1,3
        d(i)=xn(i)/sqrt(xnorm)
      enddo
c     -----------------------------------------------------------------
c     now determine the reciprocal lattice vectors
      call cross(an2,d,b1,3)
      call cross(an1,d,b2,3)
      an1b1=dot(an1,b1,3)
      an2b2=dot(an2,b2,3)
      do k=1,3
        b1(k)=2*pi*b1(k)/an1b1
        b2(k)=2*pi*b2(k)/an2b2
      enddo
c     now determine whether b2 or -b2 has the smallest angle with b1
      b1b1=dot(b1,b1,3)
      b2b2=dot(b2,b2,3)
      b1b2=dot(b1,b2,3)
      if(b1b2.lt.0)then
        do k=1,3
          b2(k)=-b2(k)
        enddo
      endif
      b1b2=-b1b2
      cos12=acos(b1b2/sqrt(b1b1*b2b2))
c     -----------------------------------------------------------------
c     now classify the reciprocal lattice
      if(myid.eq.0)write(*,*)
      if(myid.eq.0)write(*,*)
      irecip=0
      if((abs(b1b2).lt.1.d-8).and.(abs(b1b1-b2b2).lt.1.d-8))then
        if(myid.eq.0)write(*,*)' reciprocal lattice is cubic'
        irecip=1
        goto 400
      elseif((abs(b1b2).lt.1.d-8).and.(abs(b1b1-b2b2).ge.1.d-8))then
        if(myid.eq.0)write(*,*)' recip lattice is primitive-rectangular'
        irecip=2
        goto 400
      elseif((abs(cos12-pi/3.d0).lt.1.d-8).and.
     $  (abs(b1b1-b2b2).lt.1.d-8))then
        if(myid.eq.0)write(*,*)' reciprocal lattice is hexagonal'
        irecip=3
        goto 400
      elseif((abs(b1b2).gt.1.d-8).and.(abs(b1b1-b2b2).lt.1.d-8))then
        if(myid.eq.0)write(*,*)' recip lattice is centred-rectangular'
        irecip=4
        goto 400
      endif
400   continue
      if(myid.eq.0)write(*,*)
      if(myid.eq.0)write(*,*)
      return
      end
c     *****************************************************************
c     *****************************************************************
c     *****************************************************************
C     MATHEMATICAL ROUTINES
c     *****************************************************************
c     *****************************************************************
c
      double precision function gmean(x,y)
      implicit double precision (a-h,o-y)
      isgn=1
      if(x*y.ge.0.d0)then
        if(x.lt.0.d0)isgn=-1
        gmean=isgn*sqrt(x*y)
      else
        gmean=(x+y)/2.d0
      endif
      return
      end
c
c     *****************************************************************
c
      double precision function dot(x,y,n)
      implicit double precision (a-h,o-z)
      dimension x(n),y(n)
      xy=0.d0
      do i=1,n
        xy=xy+x(i)*y(i)
      enddo
      dot=xy
      return
      end
c
c     *****************************************************************
c
      subroutine cross(x,y,xy,n)
      implicit double precision (a-h,o-z)
      dimension x(n),y(n),xy(n)
      xy(1)=x(2)*y(3)-x(3)*y(2)
      xy(2)=x(3)*y(1)-x(1)*y(3)
      xy(3)=x(1)*y(2)-x(2)*y(1)
      return
      end
c
c     *****************************************************************
c
      double precision function arg(z)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      pi=acos(-1.d0)
      eps=1.d-16
      xr=dreal(z)
      xi=dimag(z)

      if(cdabs(z).lt.eps)then
        write(*,*)' ERROR function arg : cdabs(z) = 0'
        stop
      endif

      if((xr.gt.0).and.(xi.ge.0))arg=atan(xi/xr)
      if((xr.gt.0).and.(xi.le.0))arg=atan(xi/xr)
      if((xr.lt.0).and.(xi.ge.0))arg=atan(xi/xr)+pi
      if((xr.lt.0).and.(xi.le.0))arg=atan(xi/xr)-pi

      if(xr.eq.0)then
        if(xi.gt.0)arg=pi/2
        if(xi.lt.0)arg=-pi/2
      endif
      return
      end
c
c     *****************************************************************
c
      subroutine remultiply(x,y,xy,nmat,nmatx)
      implicit double precision (a-h,o-z)
      dimension x(nmatx,nmatx),y(nmatx,nmatx),xy(nmatx,nmatx)
      do ir=1,nmat
        do is=1,nmat
          xy(ir,is)=0.d0
          do k=1,nmat
            xy(ir,is)=xy(ir,is)+x(ir,k)*y(k,is)
          enddo
        enddo
      enddo
      return
      end
c
c     *****************************************************************
c
      subroutine multiply(zx,zy,zxy,nmat,nmatx)
C     implicit double precision (a-h,o-y)
C     implicit complex*16 (z)
C     dimension zx(nmatx,nmatx),zy(nmatx,nmatx),zxy(nmatx,nmatx)
C     do ir=1,nmat
C       do is=1,nmat
C         zxy(ir,is)=0.d0
C         do k=1,nmat
C           zxy(ir,is)=zxy(ir,is)+zx(ir,k)*zy(k,is)
C         enddo
C       enddo
C     enddo
c
c     ----------------------------------
c     BLAS Level3 in DXML
c     ----------------------------------
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      dimension zx(nmatx,nmatx),zy(nmatx,nmatx),zxy(nmatx,nmatx)
      zalpha=1.0d0
      zbeta =0.0d0
      CALL ZGEMM('n','n',nmat,nmat,nmat,zalpha,zx,nmatx,zy,nmatx,
     &           zbeta,zxy,nmatx)

      return
      end
c
c     *****************************************************************
c
      subroutine invers(array,n,nx) 
C       complex*16 array(nx,nx)
C       complex*16 amax, save
C       parameter (iworkx=100)
C       integer ik(iworkx),jk(iworkx) 
C       if(iworkx.lt.n)then
C         write(*,*)' ERROR INVERS : iworkx < n'
C         stop
C       endif
C 11    do 100 k=1,n 
C         amax=0.d0 
C 21      do 30 i=k,n 
C           do 30 j=k,n 
C 23          if(cdabs(amax)-cdabs(array(i,j)))24,24,30 
C 24          amax=array(i,j) 
C             ik(k)=i 
C             jk(k)=j 
C 30      continue 
C 41      i=ik(k) 
C         if(i-k)21,51,43 
C 43      do 50 j=1,n 
C           save=array(k,j) 
C           array(k,j)=array(i,j) 
C 50        array(i,j)=-save 
C 51      j=jk(k) 
C         if(j-k)21,61,53 
C 53      do 60 i=1,n 
C           save=array(i,k) 
C           array(i,k)=array(i,j) 
C 60        array(i,j)=-save 
C 61      do 70 i=1,n 
C           if(i-k)63,70,63 
C 63        array(i,k)=-array(i,k)/amax 
C 70      continue 
C 71      do 80 i=1,n 
C           do 80 j=1,n 
C             if(i-k)74,80,74 
C 74          if(j-k)75,80,75 
C 75          array(i,j)=array(i,j)+array(i,k)*array(k,j) 
C 80      continue 
C 81      do 90 j=1,n 
C           if(j-k)83,90,83 
C 83        array(k,j)=array(k,j)/amax 
C 90      continue 
C         array(k,k)=1.d0/amax 
C 100   continue 
C 101   do 130 l=1,n 
C         k=n-l+1 
C         j=ik(k) 
C         if(j-k)111,111,105 
C 105     do 110 i=1,n 
C           save=array(i,k) 
C           array(i,k)=-array(i,j) 
C 110       array(i,j)=save 
C 111     i=jk(k) 
C         if(i-k)130,130,1013
C 1013    do 120 j=1,n 
C           save=array(k,j) 
C           array(k,j)=-array(i,j) 
C 120       array(i,j)=save 
C 130   continue 
C 140   return 
C       end

c
c     ----------------------------------
c     LAPACK routine in DXML
c     ----------------------------------
      parameter (iworkx=10000)
      parameter (npfm=1)                    !!!!! npfm may not be needed.
      parameter (jwork=npfm*iworkx)         !!!!!
      integer    ipiv(iworkx)
      complex*16 array(nx,nx),work(jwork)
      if(iworkx.lt.n)then
        write(*,*)' ERROR INVERS : iworkx < n'
        stop
      endif
*     LU factorization using partial pivoting with row interchanges
      CALL ZGETRF(n,n,array,nx,ipiv,info)
      if(info.ne.0) then
        write(*,*) 'ERROR occurs in ZGETRF.   info=',info
        stop
      endif
*     inverse of a matrix using the LU factorization computed by ZGETRF
      CALL ZGETRI(n,array,nx,ipiv,work,jwork,info)
ccc   write(*,*) 'work(1) =',work(1) !Check optimized npfm
      if(info.ne.0) then
        write(*,*) 'ERROR occurs in ZGETRI.   info=',info
        stop
        endif
      return
      end
c
c     *****************************************************************
c
c
c     *****************************************************************
c
      subroutine reinvers(array,n,nx) 
      double precision array(nx,nx)
      double precision amax, save
      parameter (iworkx=100)
      integer ik(iworkx),jk(iworkx) 
      if(iworkx.lt.n)then
        write(*,*)' ERROR INVERS : iworkx < n'
        stop
      endif
11    do 100 k=1,n 
        amax=0.d0 
21      do 30 i=k,n 
          do 30 j=k,n 
23          if(abs(amax)-abs(array(i,j)))24,24,30 
24          amax=array(i,j) 
            ik(k)=i 
            jk(k)=j 
30      continue 
41      i=ik(k) 
        if(i-k)21,51,43 
43      do 50 j=1,n 
          save=array(k,j) 
          array(k,j)=array(i,j) 
50        array(i,j)=-save 
51      j=jk(k) 
        if(j-k)21,61,53 
53      do 60 i=1,n 
          save=array(i,k) 
          array(i,k)=array(i,j) 
60        array(i,j)=-save 
61      do 70 i=1,n 
          if(i-k)63,70,63 
63        array(i,k)=-array(i,k)/amax 
70      continue 
71      do 80 i=1,n 
          do 80 j=1,n 
            if(i-k)74,80,74 
74          if(j-k)75,80,75 
75          array(i,j)=array(i,j)+array(i,k)*array(k,j) 
80      continue 
81      do 90 j=1,n 
          if(j-k)83,90,83 
83        array(k,j)=array(k,j)/amax 
90      continue 
        array(k,k)=1.d0/amax 
100   continue 
101   do 130 l=1,n 
        k=n-l+1 
        j=ik(k) 
        if(j-k)111,111,105 
105     do 110 i=1,n 
          save=array(i,k) 
          array(i,k)=-array(i,j) 
110       array(i,j)=save 
111     i=jk(k) 
        if(i-k)130,130,1013
1013    do 120 j=1,n 
          save=array(k,j) 
          array(k,j)=-array(i,j) 
120       array(i,j)=save 
130   continue 
140   return 
      end
c
c     *****************************************************************
c

c
c     *****************************************************************
c
      subroutine cdet(nn,a,ns,detr,deti,au,indx)
      implicit double precision (a-h,o-z)
      integer nn,ns,i,j
      integer indx(ns)
      double precision d,detr,deti
      complex*16 det
      complex*16 a(ns,ns),au(ns,ns)
c
c     calcula o determinante de uma matriz complexa
c
c     input:  a    - matriz de dimensao (ns,ns)
c             ns   - dimensao da matriz
c     output: indx - array de trabalho de dimensao nsize
c             detr - parte real do determinante
c             deti - parte imaginaria do determinante
c
      do 10 i=1,ns
      do 10 j=1,ns
         au(i,j) = a(i,j)
 10   continue
      call ludcmp(au,nn,ns,indx,d)
      det = dcmplx(1.d0,0.d0)
      do 11 j=1,nn
         det = det*dcmplx(d,0.d0)*au(j,j)
 11   continue
      detr = dreal(det)
      deti = dimag(det)
      return
      end
c
c     *****************************************************************
c
      subroutine ludcmp(a,n,np,indx,d)
      implicit double precision (a-h,o-z)
      integer  nmatx,n,np,i,j,k,imax
      double precision   tiny,d,aamax,dum
      parameter (nplx=75,nslx=1,nspinx=10,nsubx=4,nsubatx=4)
      parameter (nmatx=nslx*nsubx*nspinx,natomx=nslx*nspinx*nsubatx)
      integer    indx(np)
      double precision     vv(nmatx)
      complex*16 a(np,np)
      complex*16 czero,sum,zdum
c
      tiny=1.0d-20
      d = 1.d0
      czero = (0.d0,0.d0)
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (cdabs(a(i,j)).gt.aamax) aamax=cdabs(a(i,j))
11      continue
        if (aamax.eq.0.d0) pause 'Singular matrix.'
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        if (j.gt.1) then
          do 14 i=1,j-1
            sum=a(i,j)
            if (i.gt.1)then
              do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
13            continue
              a(i,j)=sum
            endif
14        continue
        endif
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          if (j.gt.1)then
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
15          continue
            a(i,j)=sum
          endif
          dum=vv(i)*cdabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            zdum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=zdum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(j.ne.n)then
          if(a(j,j).eq.czero)a(j,j)=dcmplx(tiny,tiny)
          zdum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*zdum
18        continue
        endif
19    continue
      if(a(n,n).eq.czero)a(n,n)=dcmplx(tiny,tiny)
      return
      end   
c
c     *****************************************************************
c
      subroutine sort(x,ka,l)
      implicit double precision (a-h,o-z)
c
c     subroutine sorts the array x
c     the indices (ascending order) are given in array ka
c
      parameter (nvecx=500000)
      dimension x(l),ka(l),kb(nvecx)
      if(l.gt.nvecx)then
        write(*,*)' ERROR SORT.F : l > nvecx ',l,nvecx
        stop
      endif
      do 13 j=1,l
        kb(j)=j
13      ka(j)=j
      l1=1
200   j=1
20    limi1=l1*2*(j-1)+1
      if(limi1.gt.l) go to 28
      limi2=limi1+l1
      lims1=min(limi2-1,l)
      if(limi2.gt.l) go to 28
      lims2=min(limi2+l1-1,l)
      ir=limi1
      i1=ir
      i2=limi2
21    k1=ka(i1)
      k2=ka(i2)
      if(x(k1).le.x(k2)) go to 22
      kb(ir)=k2
      i2=i2+1
      ir=ir+1
      if(i2.le.lims2) go to 21
      go to 23
22    kb(ir)=k1
      i1 = i1+1
      ir=ir+1
      if(i1.le.lims1) go to 21
24    k2=ka(i2)
      kb(ir) = k2
      i2 = i2+1
      ir = ir+1
      if(i2.le.lims2) go to 24
      go to 25
23    k1=ka(i1)
      kb(ir) = k1
      i1=i1+1
      ir=ir+1
      if(i1.le.lims1) go to 23
25    j=j+1
      go to 20
28    do 280 k=1,l
280     ka(k)=kb(k)
      l1 = l1*2
      if(l1.lt.2*l) go to 200
      return
      end
c
c     *****************************************************************
c     *****************************************************************
c     ATOM DEPENDENT ROUTINES
c     *****************************************************************
c     *****************************************************************
c
      subroutine param
c     THIS ROUTINE IS ATOM DEPENDENT :-
      implicit double precision (a-h,o-z)
      parameter (numatx=8)       !! No of atom types
      parameter (nnx=4,npmax=100) !! upto nnx.th N.N's
      dimension sss(numatx,nnx,-1:1),sps(numatx,nnx,-1:1)
      dimension pps(numatx,nnx,-1:1),pdp(numatx,nnx,-1:1)
      dimension ppp(numatx,nnx,-1:1),dds(numatx,nnx,-1:1)
      dimension sds(numatx,nnx,-1:1),pds(numatx,nnx,-1:1)
      dimension ddp(numatx,nnx,-1:1),ddd(numatx,nnx,-1:1)
      dimension xxs(numatx,nnx,-1:1),xss(numatx,nnx,-1:1)
      dimension xps(numatx,nnx,-1:1),xds(numatx,nnx,-1:1)

      common/par1/s0(numatx,-1:1),p0(numatx,-1:1),d0e(numatx,-1:1),
     $  d0t(numatx,-1:1)
      common/pard/d05(numatx,-1:1),d06(numatx,-1:1),d07(numatx,-1:1),
     $d08(numatx,-1:1),d09(numatx,-1:1)
      common/par2/sssint(numatx,numatx,nnx,-1:1),
     $  spsint(numatx,numatx,nnx,-1:1),
     $  ppsint(numatx,numatx,nnx,-1:1),pppint(numatx,numatx,nnx,-1:1),
     $  sdsint(numatx,numatx,nnx,-1:1),pdsint(numatx,numatx,nnx,-1:1),
     $  pdpint(numatx,numatx,nnx,-1:1),ddsint(numatx,numatx,nnx,-1:1),
     $  ddpint(numatx,numatx,nnx,-1:1),dddint(numatx,numatx,nnx,-1:1)
      common/par3/x0(numatx,-1:1), xxsint(numatx,numatx,nnx,-1:1),
     $  xssint(numatx,numatx,nnx,-1:1),xpsint(numatx,numatx,nnx,-1:1),
     $  xdsint(numatx,numatx,nnx,-1:1)
      common/shift/xmgo_shift,ag_shift,xinas_shift
      common/elements/iatfe(5),numat
      common/data/nmat,nins,mlay,natom,nspin,nlay,nsub,numnn
      common/atmap/itypefe(3),itypeau(2)
c
c     -----------------------------------------------------------------
c     The first index in the tight binding parameter arrays refers to
c         1: Bulk Fe
c         2: Interface Fe
c         3: Ga
c         4: As
c         5: Vacuum
c
c         Interatomic hoppings at end
c
c     -----------------------------------------------------------------
C     initialise everything
      do isp=-1,1,2
        do i=1,numat
          s0(i,isp)=0.d0                    !! on-site
          p0(i,isp)=0.d0
          d05(i,isp)=0.d0
          d06(i,isp)=0.d0
          d07(i,isp)=0.d0
          d08(i,isp)=0.d0
          d09(i,isp)=0.d0
          x0(i,isp)=0.d0                    !! on-site
          do k=1,numnn
            sss(i,k,isp)=0.d0               !! same atom hopping
            pps(i,k,isp)=0.d0
            ppp(i,k,isp)=0.d0
            dds(i,k,isp)=0.d0
            ddp(i,k,isp)=0.d0
            ddd(i,k,isp)=0.d0
            sps(i,k,isp)=0.d0
            sds(i,k,isp)=0.d0
            pds(i,k,isp)=0.d0
            pdp(i,k,isp)=0.d0

            xxs(i,k,isp)=0.d0
            xss(i,k,isp)=0.d0
            xps(i,k,isp)=0.d0
            xds(i,k,isp)=0.d0
            do j=1,numat
              sssint(i,j,k,isp)=0.d0
              spsint(i,j,k,isp)=0.d0
              ppsint(i,j,k,isp)=0.d0
              pppint(i,j,k,isp)=0.d0
              sdsint(i,j,k,isp)=0.d0
              pdsint(i,j,k,isp)=0.d0
              pdpint(i,j,k,isp)=0.d0
              ddsint(i,j,k,isp)=0.d0
              ddpint(i,j,k,isp)=0.d0
              dddint(i,j,k,isp)=0.d0
              xxsint(i,j,k,isp)=0.d0
              xssint(i,j,k,isp)=0.d0
              xpsint(i,j,k,isp)=0.d0
              xdsint(i,j,k,isp)=0.d0
            enddo
          enddo
        enddo
      enddo
c     -----------------------------------------------------------------
c
c     BULK Fe up:
c
c     on site
c
      s0(1,+1) =  1.13516d0
      p0(1,+1) =  1.81739d0
      d05(1,+1) = 0.64840d0
      d06(1,+1) = 0.64840d0
      d07(1,+1) = 0.64840d0
      d08(1,+1) = 0.62960d0
      d09(1,+1) = 0.62960d0
      x0(1,+1) =  100.d0
c
c     first n.n.
c
      sss(1,1,+1) = -0.12950d0
      pps(1,1,+1) =  0.25741d0
      ppp(1,1,+1) =  0.02422d0
      dds(1,1,+1) = -0.04541d0
      ddp(1,1,+1) =  0.02714d0
      ddd(1,1,+1) = -0.00260d0
      sps(1,1,+1) =  0.17363d0
      sds(1,1,+1) = -0.06115d0
      pds(1,1,+1) = -0.08485d0
      pdp(1,1,+1) =  0.01778d0
c
c     second n.n.
c
      sss(1,2,+1) = -0.02915d0
      pps(1,2,+1) =  0.16827d0
      ppp(1,2,+1) =  0.04112d0
      dds(1,2,+1) = -0.02713d0
      ddp(1,2,+1) =  0.00589d0
      ddd(1,2,+1) =  0.00060d0
      sps(1,2,+1) =  0.06571d0
      sds(1,2,+1) = -0.03560d0
      pds(1,2,+1) = -0.05473d0
      pdp(1,2,+1) = -0.00280d0
c
c     third n.n.
c
      sss(1,3,+1) =  0.01595d0
      pps(1,3,+1) = -0.04985d0
      ppp(1,3,+1) =  0.01796d0
      dds(1,3,+1) =  0.00112d0
      ddp(1,3,+1) =  0.00034d0
      ddd(1,3,+1) = -0.00056d0
      sps(1,3,+1) = -0.02477d0
      sds(1,3,+1) = -0.00073d0
      pds(1,3,+1) = -0.00082d0
      pdp(1,3,+1) = -0.00241d0
c
c     -----------------------------------------------------------------
c     BULK Fe down:
c
c     on site
c
      s0(1,-1) =  1.14481d0
      p0(1,-1) =  1.80769d0
      d05(1,-1) = 0.78456d0
      d06(1,-1) = 0.78456d0
      d07(1,-1) = 0.78456d0
      d08(1,-1) = 0.75661d0
      d09(1,-1) = 0.75661d0
      x0(1,-1) =  100.d0
c
c     first n.n.
c
      sss(1,1,-1) = -0.13243d0
      pps(1,1,-1) =  0.25911d0
      ppp(1,1,-1) =  0.02653d0
      dds(1,1,-1) = -0.05266d0
      ddp(1,1,-1) =  0.03276d0
      ddd(1,1,-1) = -0.00286d0
      sps(1,1,-1) =  0.17278d0
      sds(1,1,-1) = -0.07145d0
      pds(1,1,-1) = -0.09702d0
      pdp(1,1,-1) =  0.02129d0
c
c     second n.n.
c
      sss(1,2,-1) = -0.03003d0
      pps(1,2,-1) =  0.18256d0
      ppp(1,2,-1) =  0.03703d0
      dds(1,2,-1) = -0.03396d0
      ddp(1,2,-1) =  0.00581d0
      ddd(1,2,-1) =  0.00114d0
      sps(1,2,-1) =  0.07159d0
      sds(1,2,-1) = -0.04075d0
      pds(1,2,-1) = -0.06522d0
      pdp(1,2,-1) = -0.00467d0
c
c     third n.n.
c
      sss(1,3,-1) =  0.01589d0
      pps(1,3,-1) = -0.04253d0
      ppp(1,3,-1) =  0.01538d0
      dds(1,3,-1) =  0.00233d0
      ddp(1,3,-1) =  0.00013d0
      ddd(1,3,-1) = -0.00060d0
      sps(1,3,-1) = -0.02306d0
      sds(1,3,-1) =  0.00016d0
      pds(1,3,-1) =  0.00222d0
      pdp(1,3,-1) = -0.00351d0
c
c
c     -----------------------------------------------------------------
c     Ag
c
      do isp=-1,1,2

c       on site
        s0(2,isp) =  0.68297d0 + ag_shift
        p0(2,isp) =  1.13432d0 + ag_shift
        d06(2,isp) = 0.12249d0 + ag_shift   !! note rotated principal axes
        d07(2,isp) = 0.12249d0 + ag_shift
        d08(2,isp) = 0.12249d0 + ag_shift
        d05(2,isp) = 0.12006d0 + ag_shift
        d09(2,isp) = 0.12006d0 + ag_shift
        x0(2,isp) =  100.d0
c
c       first n.n.
c
        sss(2,1,isp) = -0.06581d0
        sps(2,1,isp) =  0.09781d0
        pps(2,1,isp) =  0.15752d0
        ppp(2,1,isp) =  0.00694d0
        sds(2,1,isp) = -0.03110d0
        pds(2,1,isp) = -0.03905d0
        pdp(2,1,isp) =  0.01519d0
        dds(2,1,isp) = -0.03151d0
        ddp(2,1,isp) =  0.01757d0
        ddd(2,1,isp) = -0.00336d0
c
c       second n.n.
c
        sss(2,2,isp) =  0.00143d0
        sps(2,2,isp) =  0.00545d0
        pps(2,2,isp) =  0.03971d0
        ppp(2,2,isp) =  0.00434d0
        sds(2,2,isp) = -0.00462d0
        pds(2,2,isp) = -0.00065d0
        pdp(2,2,isp) =  0.00172d0
        dds(2,2,isp) = -0.00282d0
        ddp(2,2,isp) =  0.00171d0
        ddd(2,2,isp) = -0.00038d0
c
c       third n.n.
c
        sss(2,3,isp) =  0.d0
        pps(2,3,isp) =  0.d0
        ppp(2,3,isp) =  0.d0
        dds(2,3,isp) =  0.d0
        ddp(2,3,isp) =  0.d0
        ddd(2,3,isp) =  0.d0
        sps(2,3,isp) =  0.d0
        sds(2,3,isp) =  0.d0
        pds(2,3,isp) =  0.d0
        pdp(2,3,isp) =  0.d0

      enddo
c
c     -----------------------------------------------------------------
c
c     In :
c     On-site InAs !!! from jancu:    Phys Rev B  57 6493 (1998)
c     Off site determined further down this routine
      ryd=13.6058d0

      do isp=-1,1,2
c
        s0(3,isp) =   0.3333d0/ryd + xinas_shift
        p0(3,isp) =   6.4939d0/ryd + xinas_shift
        d05(3,isp) = 12.1954d0/ryd + xinas_shift
        d06(3,isp) = 12.1954d0/ryd + xinas_shift
        d07(3,isp) = 12.1954d0/ryd + xinas_shift
        d08(3,isp) = 12.1954d0/ryd + xinas_shift
        d09(3,isp) = 12.1954d0/ryd + xinas_shift
        x0(3,isp) =  17.8411d0/ryd + xinas_shift
c
c       first n.n.
c
        sss(3,1,isp) = harrison(2.80d0,0.63d0,"sss")/ryd
        pps(3,1,isp) = harrison(2.80d0,0.63d0,"pps")/ryd
        ppp(3,1,isp) = harrison(2.80d0,0.63d0,"ppp")/ryd
        dds(3,1,isp) = harrison(2.80d0,0.63d0,"dds")/ryd
        ddp(3,1,isp) = harrison(2.80d0,0.63d0,"ddp")/ryd
        ddd(3,1,isp) = harrison(2.80d0,0.63d0,"ddp")/ryd
        sps(3,1,isp) = harrison(2.80d0,0.63d0,"sps")/ryd
        sds(3,1,isp) = harrison(2.80d0,0.63d0,"sds")/ryd
        pds(3,1,isp) = harrison(2.80d0,0.63d0,"pds")/ryd
        pdp(3,1,isp) = harrison(2.80d0,0.63d0,"pdp")/ryd
      enddo

c     -----------------------------------------------------------------
c     As :
c
      do isp=-1,1,2
        s0(4,isp) =  -5.9801d0/ryd + xinas_shift
        p0(4,isp) =   3.5813d0/ryd + xinas_shift
        d05(4,isp) = 12.1954d0/ryd + xinas_shift
        d06(4,isp) = 12.1954d0/ryd + xinas_shift
        d07(4,isp) = 12.1954d0/ryd + xinas_shift
        d08(4,isp) = 12.1954d0/ryd + xinas_shift
        d09(4,isp) = 12.1954d0/ryd + xinas_shift
        x0(4,isp) =  17.8411d0/ryd + xinas_shift

c
c       first n.n.
c
        sss(4,1,isp) = harrison(2.44d0,0.51d0,"sss")/ryd
        pps(4,1,isp) = harrison(2.44d0,0.51d0,"pps")/ryd
        ppp(4,1,isp) = harrison(2.44d0,0.51d0,"ppp")/ryd
        dds(4,1,isp) = harrison(2.44d0,0.51d0,"dds")/ryd
        ddp(4,1,isp) = harrison(2.44d0,0.51d0,"ddp")/ryd
        ddd(4,1,isp) = harrison(2.44d0,0.51d0,"ddp")/ryd
        sps(4,1,isp) = harrison(2.44d0,0.51d0,"sps")/ryd
        sds(4,1,isp) = harrison(2.44d0,0.51d0,"sds")/ryd
        pds(4,1,isp) = harrison(2.44d0,0.51d0,"pds")/ryd
        pdp(4,1,isp) = harrison(2.44d0,0.51d0,"pdp")/ryd
      enddo


c     -----------------------------------------------------------------
c     Fe IF with Ag up:
c
c     On site
      s0(5,+1) =  0.99333d0
      p0(5,+1) =  1.67556d0
      d05(5,+1) = 0.62020d0
      d06(5,+1) = 0.62020d0
      d07(5,+1) = 0.62020d0
      d08(5,+1) = 0.59682d0
      d09(5,+1) = 0.59682d0
      x0(5,+1) =  100.d0
      
c
c     first n.n.
c
      sss(5,1,+1) = -0.12950d0
      pps(5,1,+1) =  0.25741d0
      ppp(5,1,+1) =  0.02422d0
      dds(5,1,+1) = -0.04541d0
      ddp(5,1,+1) =  0.02714d0
      ddd(5,1,+1) = -0.00260d0
      sps(5,1,+1) =  0.17363d0
      sds(5,1,+1) = -0.06115d0
      pds(5,1,+1) = -0.08485d0
      pdp(5,1,+1) =  0.01778d0
c
c     second n.n.
c
      sss(5,2,+1) = -0.02915d0
      pps(5,2,+1) =  0.16827d0
      ppp(5,2,+1) =  0.04112d0
      dds(5,2,+1) = -0.02713d0
      ddp(5,2,+1) =  0.00589d0
      ddd(5,2,+1) =  0.00060d0
      sps(5,2,+1) =  0.06571d0
      sds(5,2,+1) = -0.03560d0
      pds(5,2,+1) = -0.05473d0
      pdp(5,2,+1) = -0.00280d0
c
c     third n.n.
c
      sss(5,3,+1) =  0.01595d0
      pps(5,3,+1) = -0.04985d0
      ppp(5,3,+1) =  0.01796d0
      dds(5,3,+1) =  0.00112d0
      ddp(5,3,+1) =  0.00034d0
      ddd(5,3,+1) = -0.00056d0
      sps(5,3,+1) = -0.02477d0
      sds(5,3,+1) = -0.00073d0
      pds(5,3,+1) = -0.00082d0
      pdp(5,3,+1) = -0.00241d0
c
c     -----------------------------------------------------------------
c     Fe IF with Ag down:
c
c     on site
c
      s0(5,-1) =  1.00298d0
      p0(5,-1) =  1.66586d0
      d05(5,-1) = 0.78840d0
      d06(5,-1) = 0.78840d0
      d07(5,-1) = 0.78840d0
      d08(5,-1) = 0.76502d0
      d09(5,-1) = 0.76502d0
      x0(5,-1) =  100.d0
c
c     first n.n.
c
      sss(5,1,-1) = -0.13243d0
      pps(5,1,-1) =  0.25911d0
      ppp(5,1,-1) =  0.02653d0
      dds(5,1,-1) = -0.05266d0
      ddp(5,1,-1) =  0.03276d0
      ddd(5,1,-1) = -0.00286d0
      sps(5,1,-1) =  0.17278d0
      sds(5,1,-1) = -0.07145d0
      pds(5,1,-1) = -0.09702d0
      pdp(5,1,-1) =  0.02129d0
c
c     second n.n.
c
      sss(5,2,-1) = -0.03003d0
      pps(5,2,-1) =  0.18256d0
      ppp(5,2,-1) =  0.03703d0
      dds(5,2,-1) = -0.03396d0
      ddp(5,2,-1) =  0.00581d0
      ddd(5,2,-1) =  0.00114d0
      sps(5,2,-1) =  0.07159d0
      sds(5,2,-1) = -0.04075d0
      pds(5,2,-1) = -0.06522d0
      pdp(5,2,-1) = -0.00467d0
c
c     third n.n.
c
      sss(5,3,-1) =  0.01589d0
      pps(5,3,-1) = -0.04253d0
      ppp(5,3,-1) =  0.01538d0
      dds(5,3,-1) =  0.00233d0
      ddp(5,3,-1) =  0.00013d0
      ddd(5,3,-1) = -0.00060d0
      sps(5,3,-1) = -0.02306d0
      sds(5,3,-1) =  0.00016d0
      pds(5,3,-1) =  0.00222d0
      pdp(5,3,-1) = -0.00351d0
c
c
c     -----------------------------------------------------------------
c     Ag IF with Fe
c     -----------------------------------------------------------------
CCC   csp6=0.0852229d0
CCC   cd6=-0.07448307d0
      csp8=0.06889933d0
      cd8=-0.164046d0

      do isp=-1,1,2

c       on site
        s0(6,isp) = 0.68297d0 + ag_shift + csp8
        p0(6,isp) = 1.13432d0 + ag_shift + csp8
        d06(6,isp) =0.12249d0 + ag_shift + cd8   !! note rotated principal axes
        d07(6,isp) =0.12249d0 + ag_shift + cd8
        d08(6,isp) =0.12249d0 + ag_shift + cd8
        d05(6,isp) =0.12006d0 + ag_shift + cd8
        d09(6,isp) =0.12006d0 + ag_shift + cd8
        x0(6,isp)=100.d0
c
c       first n.n.
c
        sss(6,1,isp) = -0.06581d0
        sps(6,1,isp) =  0.09781d0
        pps(6,1,isp) =  0.15752d0
        ppp(6,1,isp) =  0.00694d0
        sds(6,1,isp) = -0.03110d0
        pds(6,1,isp) = -0.03905d0
        pdp(6,1,isp) =  0.01519d0
        dds(6,1,isp) = -0.03151d0
        ddp(6,1,isp) =  0.01757d0
        ddd(6,1,isp) = -0.00336d0
c
c       second n.n.
c
        sss(6,2,isp) =  0.00143d0
        sps(6,2,isp) =  0.00545d0
        pps(6,2,isp) =  0.03971d0
        ppp(6,2,isp) =  0.00434d0
        sds(6,2,isp) = -0.00462d0
        pds(6,2,isp) = -0.00065d0
        pdp(6,2,isp) =  0.00172d0
        dds(6,2,isp) = -0.00282d0
        ddp(6,2,isp) =  0.00171d0
        ddd(6,2,isp) = -0.00038d0
c
c       third n.n.
c
        sss(6,3,isp) =  0.d0
        pps(6,3,isp) =  0.d0
        ppp(6,3,isp) =  0.d0
        dds(6,3,isp) =  0.d0
        ddp(6,3,isp) =  0.d0
        ddd(6,3,isp) =  0.d0
        sps(6,3,isp) =  0.d0
        sds(6,3,isp) =  0.d0
        pds(6,3,isp) =  0.d0
        pdp(6,3,isp) =  0.d0

      enddo
c     -----------------------------------------------------------------
c     Fe IF with MgO up:
c     -----------------------------------------------------------------
      csp7=-0.310269d0
      cd7=4.9467782d-003
c
c     on site
c
      s0(7,+1) =  1.13516d0 + csp7
      p0(7,+1) =  1.81739d0 + csp7
      d05(7,+1) = 0.64840d0 - cd7
      d06(7,+1) = 0.64840d0 - cd7
      d07(7,+1) = 0.64840d0 - cd7
      d08(7,+1) = 0.62960d0 - cd7
      d09(7,+1) = 0.62960d0 - cd7
      x0(7,+1) =  100.d0
c
c     first n.n.
c
      sss(7,1,+1) = -0.12950d0
      pps(7,1,+1) =  0.25741d0
      ppp(7,1,+1) =  0.02422d0
      dds(7,1,+1) = -0.04541d0
      ddp(7,1,+1) =  0.02714d0
      ddd(7,1,+1) = -0.00260d0
      sps(7,1,+1) =  0.17363d0
      sds(7,1,+1) = -0.06115d0
      pds(7,1,+1) = -0.08485d0
      pdp(7,1,+1) =  0.01778d0
c
c     second n.n.
c
      sss(7,2,+1) = -0.02915d0
      pps(7,2,+1) =  0.16827d0
      ppp(7,2,+1) =  0.04112d0
      dds(7,2,+1) = -0.02713d0
      ddp(7,2,+1) =  0.00589d0
      ddd(7,2,+1) =  0.00060d0
      sps(7,2,+1) =  0.06571d0
      sds(7,2,+1) = -0.03560d0
      pds(7,2,+1) = -0.05473d0
      pdp(7,2,+1) = -0.00280d0
c
c     third n.n.
c
      sss(7,3,+1) =  0.01595d0
      pps(7,3,+1) = -0.04985d0
      ppp(7,3,+1) =  0.01796d0
      dds(7,3,+1) =  0.00112d0
      ddp(7,3,+1) =  0.00034d0
      ddd(7,3,+1) = -0.00056d0
      sps(7,3,+1) = -0.02477d0
      sds(7,3,+1) = -0.00073d0
      pds(7,3,+1) = -0.00082d0
      pdp(7,3,+1) = -0.00241d0
c
c     -----------------------------------------------------------------
c     Fe IF with MgO down:
c
c     on site
c
      s0(7,-1) =  1.14481d0 + csp7
      p0(7,-1) =  1.80769d0 + csp7
      d05(7,-1) = 0.78456d0 + cd7
      d06(7,-1) = 0.78456d0 + cd7
      d07(7,-1) = 0.78456d0 + cd7
      d08(7,-1) = 0.75661d0 + cd7
      d09(7,-1) = 0.75661d0 + cd7
      x0(7,-1) =  100.d0
c
c     first n.n.
c
      sss(7,1,-1) = -0.13243d0
      pps(7,1,-1) =  0.25911d0
      ppp(7,1,-1) =  0.02653d0
      dds(7,1,-1) = -0.05266d0
      ddp(7,1,-1) =  0.03276d0
      ddd(7,1,-1) = -0.00286d0
      sps(7,1,-1) =  0.17278d0
      sds(7,1,-1) = -0.07145d0
      pds(7,1,-1) = -0.09702d0
      pdp(7,1,-1) =  0.02129d0
c
c     second n.n.
c
      sss(7,2,-1) = -0.03003d0
      pps(7,2,-1) =  0.18256d0
      ppp(7,2,-1) =  0.03703d0
      dds(7,2,-1) = -0.03396d0
      ddp(7,2,-1) =  0.00581d0
      ddd(7,2,-1) =  0.00114d0
      sps(7,2,-1) =  0.07159d0
      sds(7,2,-1) = -0.04075d0
      pds(7,2,-1) = -0.06522d0
      pdp(7,2,-1) = -0.00467d0
c
c     third n.n.
c
      sss(7,3,-1) =  0.01589d0
      pps(7,3,-1) = -0.04253d0
      ppp(7,3,-1) =  0.01538d0
      dds(7,3,-1) =  0.00233d0
      ddp(7,3,-1) =  0.00013d0
      ddd(7,3,-1) = -0.00060d0
      sps(7,3,-1) = -0.02306d0
      sds(7,3,-1) =  0.00016d0
      pds(7,3,-1) =  0.00222d0
      pdp(7,3,-1) = -0.00351d0
c
c
c     -----------------------------------------------------------------
c     Vacuum
c 
      do isp=-1,1,2

c       on site
        s0(8,isp) =  100.d0
        p0(8,isp) =  100.d0
        d05(8,isp) = 100.d0
        d06(8,isp) = 100.d0
        d07(8,isp) = 100.d0
        d08(8,isp) = 100.d0
        d09(8,isp) = 100.d0
        x0(8,isp) =  100.d0

c       ALL NN = 0
      enddo
c
c     -----------------------------------------------------------------
c
c     -----------------------------------------------------------------
c
c     =================================================================
c     Now evaluate the inter-atomic hoppings 
      do isp=-1,1,2
        do iat=1,numat
          do jat=1,numat
            do inn=1,numnn
              sssint(iat,jat,inn,isp)=
     $gmean(sss(iat,inn,isp),sss(jat,inn,isp))
              ppsint(iat,jat,inn,isp)=
     $gmean(pps(iat,inn,isp),pps(jat,inn,isp))
              pppint(iat,jat,inn,isp)=
     $gmean(ppp(iat,inn,isp),ppp(jat,inn,isp))
              ddsint(iat,jat,inn,isp)=
     $gmean(dds(iat,inn,isp),dds(jat,inn,isp))
              ddpint(iat,jat,inn,isp)=
     $gmean(ddp(iat,inn,isp),ddp(jat,inn,isp))
              dddint(iat,jat,inn,isp)=
     $gmean(ddd(iat,inn,isp),ddd(jat,inn,isp))
              spsint(iat,jat,inn,isp)=
     $gmean(sps(iat,inn,isp),sps(jat,inn,isp))
              sdsint(iat,jat,inn,isp)=
     $gmean(sds(iat,inn,isp),sds(jat,inn,isp))
              pdsint(iat,jat,inn,isp)=
     $gmean(pds(iat,inn,isp),pds(jat,inn,isp))
              pdpint(iat,jat,inn,isp)=
     $gmean(pdp(iat,inn,isp),pdp(jat,inn,isp))
            enddo
          enddo
        enddo
      enddo
c
c     -----------------------------------------------------------------
c     InAs !!! from jancu:    Phys Rev B  57 6493 (1998)
c     Note the asymmetry in Ga-As   w.r.t.   As-Ga   
c     ... need Gabriel's Slater_Koster routine to evaluate these correctly
c     ... eint2 fails!
c
      do isp=-1,1,2
        do iat=3,4
          do jat=3,4
            do inn=1,numnn
              sssint(iat,jat,inn,isp)=0.d0
              ppsint(iat,jat,inn,isp)=0.d0
              pppint(iat,jat,inn,isp)=0.d0
              ddsint(iat,jat,inn,isp)=0.d0
              ddpint(iat,jat,inn,isp)=0.d0
              dddint(iat,jat,inn,isp)=0.d0
              spsint(iat,jat,inn,isp)=0.d0
              sdsint(iat,jat,inn,isp)=0.d0
              pdsint(iat,jat,inn,isp)=0.d0
              pdpint(iat,jat,inn,isp)=0.d0
            enddo
          enddo
        enddo
      enddo

      do isp=-1,1,2
c
c       first n.n.
c
        sssint(3,4,1,isp) =  -1.4789d0/ryd
        spsint(3,4,1,isp) =   2.8006d0/ryd
        sdsint(3,4,1,isp) =  -2.4499d0/ryd
        ppsint(3,4,1,isp) =   4.1188d0/ryd
        pppint(3,4,1,isp) =  -1.3687d0/ryd
        pdsint(3,4,1,isp) =  -2.0584d0/ryd
        pdpint(3,4,1,isp) =   1.7106d0/ryd
        ddsint(3,4,1,isp) =  -1.2009d0/ryd
        ddpint(3,4,1,isp) =   2.1820d0/ryd
        dddint(3,4,1,isp) =  -1.7788d0/ryd
        xxsint(3,4,1,isp) =  -3.8514d0/ryd
        xssint(3,4,1,isp) =  -2.1320d0/ryd
        xpsint(3,4,1,isp) =   1.9012d0/ryd
        xdsint(3,4,1,isp) =  -0.8371d0/ryd
c
        sssint(4,3,1,isp) =  -1.4789d0/ryd
        spsint(4,3,1,isp) =   2.3159d0/ryd
        sdsint(4,3,1,isp) =  -2.5828d0/ryd
        ppsint(4,3,1,isp) =   4.1188d0/ryd
        pppint(4,3,1,isp) =  -1.3687d0/ryd
        pdsint(4,3,1,isp) =  -2.1222d0/ryd
        pdpint(4,3,1,isp) =   1.5462d0/ryd
        ddsint(4,3,1,isp) =  -1.2009d0/ryd
        ddpint(4,3,1,isp) =   2.1820d0/ryd
        dddint(4,3,1,isp) =  -1.7788d0/ryd
        xxsint(4,3,1,isp) =  -3.8514d0/ryd
        xssint(4,3,1,isp) =  -1.2219d0/ryd
        xpsint(4,3,1,isp) =   2.6467d0/ryd
        xdsint(4,3,1,isp) =  -0.8497d0/ryd
c
      enddo
c
c     =================================================================
      return
      end


c
c     *****************************************************************
c
      double precision function harrison(d,rd,vxxx)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      character*3 vxxx

      c=7.62d0

      if(vxxx.eq."sss")harrison=-1.40d0*c/d**2
      if(vxxx.eq."sps")harrison=1.84d0*c/d**2
      if(vxxx.eq."pps")harrison=3.24d0*c/d**2
      if(vxxx.eq."ppp")harrison=-0.81d0*c/d**2

      if(vxxx.eq."sds")harrison=-3.16d0*c*rd**(1.5)/d**(3.5)
      if(vxxx.eq."pds")harrison=-2.95d0*c*rd**(1.5)/d**(3.5)
      if(vxxx.eq."pdp")harrison=1.36d0*c*rd**(1.5)/d**(3.5)

      if(vxxx.eq."dds")harrison=-16.2d0*c*rd**3/d**5
      if(vxxx.eq."ddp")harrison=8.75d0*c*rd**3/d**5
      if(vxxx.eq."ddd")harrison=0.0d0*c*rd**3/d**5

      return
      end
c
c     *****************************************************************
c
      subroutine eint1(sss,sps,pps,ppp,ss,ps,pp,ds,dp,dd,w,b)
C     THIS ROUTINE IS SPIN DEPENDENT :-
C     IT IS WRITTEN FOR s,p and d BANDS ONLY
      implicit double precision (a-h,o-z)      
      dimension b(9,9)
      dimension w(3)
      x=w(1)
      y=w(2)
      z=w(3)
      xx=x*x
      xy=x*y
      yy=y*y
      yz=y*z
      zz=z*z
      zx=z*x
      xxyy=xx*yy
      yyzz=yy*zz
      zzxx=zz*xx
      aux=pps-ppp
      r3=dsqrt(3.d0)
      aux1=r3*ss
      f8=3.d0*zz-1.d0
      f1=xx+yy
      f2=xx-yy
      f3=zz-.5d0*f1
      g1=1.5d0*f2*ds
      g2=r3*f3*ds
      b(1,1)=sss
      b(1,2)=x*sps
      b(1,3)=y*sps
      b(1,4)=z*sps
      b(2,1)=-b(1,2)
      b(2,2)=xx*pps+(1.d0-xx)*ppp
      b(2,3)=xy*aux
      b(2,4)=zx*aux
      b(3,1)=-b(1,3)
      b(3,2)=b(2,3)
      b(3,3)=yy*pps+(1.d0-yy)*ppp
      b(3,4)=yz*aux
      b(4,1)=-b(1,4)
      b(4,2)=b(2,4)
      b(4,3)=b(3,4)
      b(4,4)=zz*pps+(1.d0-zz)*ppp
      b(1,5)=xy*aux1
      b(1,6)=yz*aux1
      b(1,7)=zx*aux1
      b(1,8)=.5d0*f2*aux1
      b(1,9)=.5d0*f8*ss
      b(5,1)=b(1,5)
      b(6,1)=b(1,6)
      b(7,1)=b(1,7)
      b(8,1)=b(1,8)
      b(9,1)=b(1,9)
      f4=.5d0*r3*f2*ps
      f5=.5d0*f8*ps
      aux2=r3*xx*ps+(1.d0-2.d0*xx)*pp
      b(2,5)=aux2*y
      b(2,6)=(r3*ps-2.d0*pp)*xy*z
      b(2,7)=aux2*z
      b(2,8)=(f4+(1.d0-f2)*pp)*x
      b(2,9)=(f5-r3*zz*pp)*x
      aux3=(r3*yy*ps+(1.d0-2.d0*yy)*pp)
      b(3,5)=aux3*x
      b(3,6)=aux3*z
      b(3,7)=b(2,6)
      b(3,8)=(f4-(1.d0+f2)*pp)*y
      b(3,9)=(f5-r3*zz*pp)*y
      aux4=r3*zz*ps+(1.d0-2.d0*zz)*pp
      b(4,5)=b(2,6)
      b(4,6)=aux4*y
      b(4,7)=aux4*x
      b(4,8)=(f4-f2*pp)*z
      b(4,9)=(f5+r3*f1*pp)*z
      b(5,2)=-b(2,5)
      b(6,2)=-b(2,6)
      b(7,2)=-b(2,7)
      b(8,2)=-b(2,8)
      b(9,2)=-b(2,9)
      b(5,3)=-b(3,5)
      b(6,3)=-b(3,6)
      b(7,3)=-b(3,7)
      b(8,3)=-b(3,8)
      b(9,3)=-b(3,9)
      b(5,4)=-b(4,5)
      b(6,4)=-b(4,6)
      b(7,4)=-b(4,7)
      b(8,4)=-b(4,8)
      b(9,4)=-b(4,9)
      b(5,5)=3.d0*xxyy*ds+(f1-4.d0*xxyy)*dp+(zz+xxyy)*dd
      b(5,6)=(3.d0*yy*ds+(1.d0-4.d0*yy)*dp+(yy-1.d0)*dd)*zx
      b(5,7)=(3.d0*xx*ds+(1.d0-4.d0*xx)*dp+(xx-1.d0)*dd)*yz
      b(5,8)=(g1-2.d0*f2*dp+.5d0*f2*dd)*xy
      b(5,9)=(g2-2.d0*r3*zz*dp+.5d0*r3*(1.d0+zz)*dd)*xy
      b(6,5)=b(5,6)
      b(6,6)=3.d0*yyzz*ds+(yy+zz-4.d0*yyzz)*dp+(xx+yyzz)*dd
      b(6,7)=(3.d0*zz*ds+(1.d0-4.d0*zz)*dp+(zz-1.d0)*dd)*xy
      b(6,8)=(g1-(1.d0+2.d0*f2)*dp+(1.d0+.5d0*f2)*dd)*yz
      b(6,9)=(g2+r3*(f1-zz)*dp-.5d0*r3*f1*dd)*yz
      b(7,5)=b(5,7)
      b(7,6)=b(6,7)
      b(7,7)=3.d0*zzxx*ds+(zz+xx-4.d0*zzxx)*dp+(yy+zzxx)*dd
      b(7,8)=(g1+(1.d0-2.d0*f2)*dp-(1.d0-.5d0*f2)*dd)*zx
      b(7,9)=(g2+r3*(f1-zz)*dp-.5d0*r3*f1*dd)*zx
      b(8,5)=b(5,8)
      b(8,6)=b(6,8)
      b(8,7)=b(7,8)
      b(8,8)=.75d0*f2*f2*ds+(f1-f2*f2)*dp+(zz+.25d0*f2*f2)*dd
      b(8,9)=.5d0*f2*g2-r3*zz*f2*dp+.25d0*r3*(1.d0+zz)*f2*dd
      b(9,5)=b(5,9)
      b(9,6)=b(6,9)
      b(9,7)=b(7,9)
      b(9,8)=b(8,9)
      b(9,9)=f3*f3*ds+3.d0*zz*f1*dp+.75d0*f1*f1*dd
      return
      end
c
c     *****************************************************************
c
      subroutine slater_koster(s0s1s,s0p1s,s0d1s,p0p1s,p0p1p,p0d1s,
     $p0d1p,d0d1s,d0d1p,d0d1d,x0x1s,x0s1s,x0p1s,x0d1s,s1s0s,s1p0s,
     $s1d0s,p1p0s,p1p0p,p1d0s,p1d0p,d1d0s,d1d0p,d1d0d,x1x0s,x1s0s,x1p0s,
     $x1d0s,w,b)
C     THIS ROUTINE IS SPIN DEPENDENT :-
C     IT IS WRITTEN FOR s,p,d and s* BANDS
      implicit double precision (a-h,o-z)      
      dimension b(10,10)
      dimension w(3)

      dd=d0d1d
      dp=d0d1p
      ds=d0d1s
      pps=p0p1s
      ppp=p0p1p       
      x=w(1)
      y=w(2)
      z=w(3)
      xx=x*x
      xy=x*y
      yy=y*y
      yz=y*z
      zz=z*z
      zx=z*x
      xxyy=xx*yy
      yyzz=yy*zz
      zzxx=zz*xx
      r3=dsqrt(3.d0)
      f8=3.d0*zz-1.d0
      f1=xx+yy
      f2=xx-yy
      f3=zz-.5d0*f1
! sss
      b(1,1)=s0s1s
! sp
      b(1,2)=x*s0p1s
      b(1,3)=y*s0p1s
      b(1,4)=z*s0p1s
      b(2,1)=-x*s1p0s
      b(3,1)=-y*s1p0s
      b(4,1)=-z*s1p0s
! pp 
      aux=pps-ppp 
      b(2,2)=xx*pps+(1.d0-xx)*ppp
      b(2,3)=xy*aux
      b(2,4)=zx*aux
      b(3,2)=b(2,3)
      b(3,3)=yy*pps+(1.d0-yy)*ppp
      b(3,4)=yz*aux
      b(4,2)=b(2,4)
      b(4,3)=b(3,4)
      b(4,4)=zz*pps+(1.d0-zz)*ppp
    
! sd
      b(1,5)=xy*r3*s0d1s
      b(1,6)=yz*r3*s0d1s
      b(1,7)=zx*r3*s0d1s
      b(1,8)=.5d0*f2*r3*s0d1s
      b(1,9)=.5d0*f8*s0d1s

      b(5,1)=xy*r3*s1d0s
      b(6,1)=yz*r3*s1d0s
      b(7,1)=zx*r3*s1d0s
      b(8,1)=.5d0*f2*r3*s1d0s
      b(9,1)=.5d0*f8*s1d0s
! pd
      f4=.5d0*r3*f2*p0d1s
      f5=.5d0*f8*p0d1s
      aux2=r3*xx*p0d1s+(1.d0-2.d0*xx)*p0d1p
      b(2,5)=aux2*y
      b(2,6)=(r3*p0d1s-2.d0*p0d1p)*xy*z
      b(2,7)=aux2*z
      b(2,8)=(f4+(1.d0-f2)*p0d1p)*x
      b(2,9)=(f5-r3*zz*p0d1p)*x
      aux3=(r3*yy*p0d1s+(1.d0-2.d0*yy)*p0d1p)
      b(3,5)=aux3*x
      b(3,6)=aux3*z
      b(3,7)=b(2,6)
      b(3,8)=(f4-(1.d0+f2)*p0d1p)*y
      b(3,9)=(f5-r3*zz*p0d1p)*y
      aux4=r3*zz*p0d1s+(1.d0-2.d0*zz)*p0d1p
      b(4,5)=b(2,6)
      b(4,6)=aux4*y
      b(4,7)=aux4*x
      b(4,8)=(f4-f2*p0d1p)*z
      b(4,9)=(f5+r3*f1*p0d1p)*z

      f4=.5d0*r3*f2*p1d0s
      f5=.5d0*f8*p1d0s
      aux2=r3*xx*p1d0s+(1.d0-2.d0*xx)*p1d0p
      b(5,2)=-aux2*y
      b(6,2)=-(r3*p1d0s-2.d0*p1d0p)*xy*z
      b(7,2)=-aux2*z
      b(8,2)=-(f4+(1.d0-f2)*p1d0p)*x
      b(9,2)=-(f5-r3*zz*p1d0p)*x
      aux3=(r3*yy*p1d0s+(1.d0-2.d0*yy)*p1d0p)
      b(5,3)=-aux3*x
      b(6,3)=-aux3*z
      b(7,3)=b(6,2)
      b(8,3)=-(f4-(1.d0+f2)*p1d0p)*y
      b(9,3)=-(f5-r3*zz*p1d0p)*y
      aux4=r3*zz*p1d0s+(1.d0-2.d0*zz)*p1d0p
      b(5,4)=b(6,2)
      b(6,4)=-aux4*y
      b(7,4)=-aux4*x
      b(8,4)=-(f4-f2*p1d0p)*z
      b(9,4)=-(f5+r3*f1*p1d0p)*z
! dd
      g1=1.5d0*f2*ds
      g2=r3*f3*ds
      b(5,5)=3.d0*xxyy*ds+(f1-4.d0*xxyy)*dp+(zz+xxyy)*dd
      b(5,6)=(3.d0*yy*ds+(1.d0-4.d0*yy)*dp+(yy-1.d0)*dd)*zx
      b(5,7)=(3.d0*xx*ds+(1.d0-4.d0*xx)*dp+(xx-1.d0)*dd)*yz
      b(5,8)=(g1-2.d0*f2*dp+.5d0*f2*dd)*xy
      b(5,9)=(g2-2.d0*r3*zz*dp+.5d0*r3*(1.d0+zz)*dd)*xy
      b(6,5)=b(5,6)
      b(6,6)=3.d0*yyzz*ds+(yy+zz-4.d0*yyzz)*dp+(xx+yyzz)*dd
      b(6,7)=(3.d0*zz*ds+(1.d0-4.d0*zz)*dp+(zz-1.d0)*dd)*xy
      b(6,8)=(g1-(1.d0+2.d0*f2)*dp+(1.d0+.5d0*f2)*dd)*yz
      b(6,9)=(g2+r3*(f1-zz)*dp-.5d0*r3*f1*dd)*yz
      b(7,5)=b(5,7)
      b(7,6)=b(6,7)
      b(7,7)=3.d0*zzxx*ds+(zz+xx-4.d0*zzxx)*dp+(yy+zzxx)*dd
      b(7,8)=(g1+(1.d0-2.d0*f2)*dp-(1.d0-.5d0*f2)*dd)*zx
      b(7,9)=(g2+r3*(f1-zz)*dp-.5d0*r3*f1*dd)*zx
      b(8,5)=b(5,8)
      b(8,6)=b(6,8)
      b(8,7)=b(7,8)
      b(8,8)=.75d0*f2*f2*ds+(f1-f2*f2)*dp+(zz+.25d0*f2*f2)*dd
      b(8,9)=.5d0*f2*g2-r3*zz*f2*dp+.25d0*r3*(1.d0+zz)*f2*dd
      b(9,5)=b(5,9)
      b(9,6)=b(6,9)
      b(9,7)=b(7,9)
      b(9,8)=b(8,9)
      b(9,9)=f3*f3*ds+3.d0*zz*f1*dp+.75d0*f1*f1*dd

! sxsx
      b(10,10)=x0x1s
! sxs
      b(1,10)=x0s1s
      b(10,1)=x1s0s
! sxp
      b(10,2)=x*x0p1s
      b(10,3)=y*x0p1s
      b(10,4)=z*x0p1s
      b(2,10)=-x*x1p0s
      b(3,10)=-y*x1p0s
      b(4,10)=-z*x1p0s
! sxd
      b(10,5)=xy*r3*x0d1s
      b(10,6)=yz*r3*x0d1s
      b(10,7)=zx*r3*x0d1s
      b(10,8)=.5d0*f2*r3*x0d1s
      b(10,9)=.5d0*f8*x0d1s
      b(5,10)=xy*r3*x1d0s
      b(6,10)=yz*r3*x1d0s
      b(7,10)=zx*r3*x1d0s
      b(8,10)=.5d0*f2*r3*x1d0s
      b(9,10)=.5d0*f8*x1d0s

      return
      end
c
c     *****************************************************************
      subroutine eint11(nn,ind,ispin,zt)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter (nplx=75,nslx=1,nspinx=10,nsubx=4,nsubatx=4)
      common/data/nmat,nins,mlay,natom,nspin,nlay,nsub,numnn
      common/data1/frat,ifrac
      common/layer/ndiff
      dimension zt(nspinx,nspinx)
      dimension d(3),c(3)
c     ?????????????????????????????????????????????????????????????????
c     ensure that this routine gives U for d = 0
c     ?????????????????????????????????????????????????????????????????
c     ?????????????????????????????????????????????????????????????????
      do ir=1,nspin
        do is=1,nspin
          zt(ir,is)=0.0d0
        enddo
      enddo
      if(nn.eq.0)then
        if(ind.eq.1.and.ispin.eq.1)zui=0.9d0
        if(ind.eq.1.and.ispin.eq.-1)zui=0.55d0
        if(ind.eq.2)zui=0.90d0
      elseif(nn.eq.1)then
        zui=0.5d0
      elseif(nn.eq.2)then
        zui=0.25d0
      elseif(nn.eq.3)then
        zui=0.1d0
      endif
      do ir=1,nspin
        zt(ir,ir)=zui
      enddo
      return
      end
c
c     *****************************************************************
c
      subroutine surfacedecim(zu0,zt0,zener,zsurfl,zsurfr,errmax,n,nx)
C     This program calculates the surface Green's function 
C     using decimation/reduction and Mobius transformation
C     The decimation part reduces the size of U and T is T is singular.
C     This part is described in Sanvito et. al. PRB 73, 085414 (2006)

C     A simpler version is given in 
C     Papers and Notes/Tight Binding/Surface Greens Functions/decimation*.mws

C     
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      dimension zu0(nx,nx),zt0(nx,nx)
      dimension zu(nx,nx),zt(nx,nx)
      dimension zured1(nx,nx),zured2(nx,nx),ztred(nx,nx)
      dimension zglred(nx,nx),zgrred(nx,nx)
      dimension zh(2*nx,2*nx),zhbulk(2*nx,2*nx),zhlast(2*nx,2*nx)
      dimension zulast(nx,nx),ztlast(nx,nx)
      dimension zsurfl(nx,nx),zsurfr(nx,nx),zsurf2(nx,nx)
      dimension zwrk(nx,nx),zunit(nx,nx)
      dimension ztdag(nx,nx),zt0dag(nx,nx),ztreddag(nx,nx)
      dimension lst(2*nx)

      dimension zp(nx,nx),zq(nx,nx),sv(nx)
      dimension rwork(5*nx),zb(1),zwork(nx+nx*nx)
      logical wantp,wantq

C     -----------------------------------------------------------------
C     The following is a key parameter: may need to change this for 
C     more accurate SGF

      svmin=5.d-7      !! if sv(i)<svmin then sv(i)=0
C     -----------------------------------------------------------------
      do ir=1,n
        do is=1,n
          zunit(ir,is)=0.d0
        enddo
        zunit(ir,ir)=1.d0
      enddo
C     -----------------------------------------------------------------
C     Find SVD of zt : zt = Q.S.P     where S(1,1) >= S(2,2,) >= ... >= 0
      wantp=.true.
      wantq=.true.
      ifail=0
      call equals(zt0,zt,n,nx)   !! note zt will be written over
      call f02xef(n,n,zt,nx,0,zb,1,wantq,zq,nx,sv,wantp,zp,nx,
     $rwork,zwork,ifail)
      if(ifail.ne.0)then
        write(*,*)" ERROR surfacedecim: f02xef ifail= ",ifail
      endif
C     zp can have some very small numbers ... set these to 0
      do i=1,n
        do j=1,n
          xr=dreal(zp(i,j))
          xi=dimag(zp(i,j))
          if(abs(xr).lt.1.d-40)xr=0.d0
          if(abs(xi).lt.1.d-40)xi=0.d0
          zp(i,j)=dcmplx(xr,xi)
        enddo
      enddo


C     Now rearrange P and S so that 0 <= S(1,1) <= S(2,2) <= ...
      do i=1,n/2
        do j=1,n
          zptmp=zp(i,j)
          zp(i,j)=zp(n-i+1,j)
          zp(n-i+1,j)=zptmp
C         zqtmp=zq(j,i)
C         zq(j,i)=zq(j,n-i+1)
C         zq(j,n-i+1)=zqtmp
        enddo
        svtmp=sv(i)
        sv(i)=sv(n-i+1)
        sv(n-i+1)=svtmp
      enddo

C     Now transform all matrices  M -> P.M.P^h
      call equals(zu0,zu,n,nx)
      call equals(zt0,zt,n,nx)
      call simtran(zu,zp,n,nx)
      call simtran(zt,zp,n,nx)
C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C     now set up lst and nlst:   the rows to be decimated
      icnt=0
      do i=1,n
        if(abs(sv(i)).lt.svmin)then
          icnt=icnt+1
          lst(icnt)=i
        endif
      enddo
      do i=1,icnt
          lst(icnt+i)=lst(i)+n
      enddo
      nlst=2*icnt

CCC   write(*,*)" n, nlst/2 =",n,nlst/2
C     -----------------------------------------------------------------
C     Create a matrix to calculate the bulk decimated/reduced on and off 
C     site elements
      do ir=1,n
        do is=1,n
          zh(ir,is)=zener*zunit(ir,is)-zu(ir,is)
          zh(ir,is+n)=-zt(ir,is)
          zh(ir+n,is)=-dconjg(zt(is,ir))
          zh(ir+n,is+n)=zener*zunit(ir,is)-zu(ir,is)
        enddo
      enddo

      call equals(zh,zhbulk,2*n,2*nx)
      do i=1,nlst
        ii=lst(i)
        call decim(zhbulk,ii,2*n,2*nx)
      enddo
      call remve(zhbulk,lst,2*n,nlst,2*nx)

C     zhbulk is (2.nred x 2.nred)
C     Now calculate the reduced/decimated SGF ie. the very bottom RH element
      nred=n-nlst/2
      do i=1,nred
        do j=1,nred
          zured1(i,j)=zener*zunit(i,j)-zhbulk(i,j)
          zured2(i,j)=zener*zunit(i,j)-zhbulk(nred+i,nred+j)
          ztred(i,j)=-zhbulk(i,nred+j)
        enddo
      enddo
C     Create SGFs
      ifail=0
      call surfacenew(zured2,ztred,zener,zglred,zgrred,nred,nx,ifail)
      if(ifail.ne.0)then
        ifail=0
        call surfacenew16(zured2,ztred,zener,zglred,zgrred,nred,nx,
     $ifail)
        if(ifail.ne.0)write(*,*)'surfacedecim  re*16 ifail =',ifail
        if(ifail.ne.0)return
      endif

C     =================================================================
C     Left Hand SGF:
C     =================================================================
C     Create the nxn non-reduced SGF by adlayering
C     First create u and t for the last layer
      call equals(zh,zhlast,2*n,2*nx)
      do i=1,nlst/2
        ii=lst(i)
        call decim(zhlast,ii,2*n,2*nx)
      enddo
      do i=1,n
        do j=1,n
          zulast(i,j)=zener*zunit(i,j)-zhlast(i+n,j+n)
        enddo
      enddo
C     Bulk out ztlast to be (n x n)
      do i=1,n
        do j=1,n
          ztlast(i,j)=-zhlast(i,j+n)
        enddo
      enddo

C     Now bulk the SGF to (n x n)
C     (could use non-square matrices instead here!)
      do i=1,n
        do j=1,n
          zsurfl(i,j)=0.d0
        enddo
      enddo
      do i=1,nred
        do j=1,nred
          zsurfl(i+n-nred,j+n-nred)=zglred(i,j)
        enddo
      enddo

C     Now adlayer on top of this
      call adlayer1(zsurfl,zulast,ztlast,zwrk,zener,n,nx)

C     Now untransform this SGF
      call hconjugate(zp,n,nx)
      call simtran(zsurfl,zp,n,nx)

C     =================================================================
C     Right Hand SGF:
C     =================================================================
      call equals(ztred,ztreddag,nred,nx)
      call hconjugate(ztreddag,nred,nx)
      call adlayer1(zgrred,zured1,ztreddag,zwrk,zener,nred,nx)

C     Now bulk the SGF to (n x n)
      do i=1,n
        do j=1,n
          zsurfr(i,j)=0.d0
        enddo
      enddo
      do i=1,nred
        do j=1,nred
          zsurfr(i+n-nred,j+n-nred)=zgrred(i,j)
        enddo
      enddo

C     Now adlayer on top of this
      call equals(zt,ztdag,n,nx)
      call hconjugate(ztdag,n,nx)
      call adlayer1(zsurfr,zu,ztdag,zwrk,zener,n,nx)

C     Now untransform this SGF
      call simtran(zsurfr,zp,n,nx)

C     =================================================================
C     =================================================================

C     Check these are the correct SGF
      call equals(zsurfl,zsurf2,n,nx)
      call adlayer1(zsurf2,zu0,zt0,zwrk,zener,n,nx)
      xmaxerrl=0.d0
      do i=1,n
        do j=1,n
          errl=abs(zsurfl(i,j)-zsurf2(i,j))
CCC       if(errl.gt.1.d-6)write(*,*)"ERROR decim zsurfl =/ zsurf2",errl,
CCC  $abs(zsurfl(i,j)),abs(zsurf2(i,j)),errl
          xmaxerrl=max(xmaxerrl,errl)
        enddo
      enddo

      call equals(zsurfr,zsurf2,n,nx)
      call equals(zt0,zt0dag,n,nx)
      call hconjugate(zt0dag,n,nx)
      call adlayer1(zsurf2,zu0,zt0dag,zwrk,zener,n,nx)
      xmaxerrr=0.d0
      do i=1,n
        do j=1,n
          errr=abs(zsurfr(i,j)-zsurf2(i,j))
CCC       if(errr.gt.1.d-6)write(*,*)"ERROR decim zsurfl =/ zsurf2",errr,
CCC  $abs(zsurfl(i,j)),abs(zsurf2(i,j)),errl
          xmaxerrr=max(xmaxerrr,errr)
        enddo
      enddo

      errmax=max(xmaxerrl,xmaxerrr)

      return
      end
C
C     *****************************************************************
C
      subroutine simtran(zt,zp,n,nx)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      dimension zt(nx,nx),zp(nx,nx),ztmp(nx,nx)
C     zt -> zp.zt.zp^h

      call multiply(zp,zt,ztmp,n,nx)
      call hconjugate(zp,n,nx)
      call multiply(ztmp,zp,zt,n,nx)
      call hconjugate(zp,n,nx)

      return
      end
C
C     *****************************************************************
C
      subroutine equals(zh0,zh1,n,nx)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      dimension zh1(nx,nx),zh0(nx,nx)
      do i=1,n
        do j=1,n
          zh1(i,j)=zh0(i,j)
        enddo
      enddo
      return
      end
C
C     *****************************************************************
C
      subroutine decim(zh0,k,n,nx)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      dimension zh(nx,nx),zh0(nx,nx)
C     Decimate row/column k from zh0
      do i=1,n
        do j=1,n
          zh(i,j)=zh0(i,j)-zh0(i,k)*zh0(k,j)/zh0(k,k)
        enddo
      enddo
      call equals(zh,zh0,n,nx)

      return
      end
C
C     *****************************************************************
C
      subroutine remve(zh0,lst,n,nlst,nx)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      dimension zh(nx,nx),zh0(nx,nx)
      dimension lst(nx),notlst(nx)
C     remove rows and columns contained in lst from zh0

C     first create a list of rows not to eliminate
      icount=0
      do i=1,n
        do j=1,nlst
          if(lst(j).eq.i)goto 100
        enddo
        icount=icount+1
        notlst(icount)=i
100     continue
      enddo

      do i=1,icount
        do j=1,icount
          ii=notlst(i)
          jj=notlst(j)
          zh(i,j)=zh0(ii,jj)
        enddo
      enddo

      if(icount.ne.n-nlst)then
        write(*,*)"ERROR REMVE: ",icount,n-nlst 
      endif

      call equals(zh,zh0,n,nx)

      return
      end
C
C     *****************************************************************
C
c
c     *****************************************************************
c
      subroutine surfacenew(zu,zt,zener,zsurfl,zsurfr,n,nx,ifail)
      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      dimension zu(nx,nx),zt(nx,nx)
      dimension zsurfl(nx,nx),zsurfr(nx,nx)
      dimension zgamma(nx,nx)
      dimension zunit(nx,nx)
      dimension ztmp1(nx,nx),ztmp2(nx,nx)
      dimension ztmp3(nx,nx)
      dimension ztinv(nx,nx),zsinv(nx,nx)
      dimension zs(nx,nx)
      dimension zp(2*nx,2*nx)
      dimension isort(2*nx)
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      dimension wk1(2*nx),wk2(2*nx),wk3(2*nx)
      dimension rr(2*nx),ri(2*nx),evlab(2*nx),zevl(2*nx)
      dimension ar(2*nx,2*nx),ai(2*nx,2*nx)
      dimension vr(2*nx,2*nx),vi(2*nx,2*nx)
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      dimension zfoo1(2*nx,2*nx),zfoo2(2*nx,2*nx)
      dimension zv(2*nx),zw(2*nx)
      dimension zdp(2*nx,2*nx)

      dimension itran1(nx),itran2(nx)
c     -----------------------------------------------------------------

      n2=2*n

      ener=dreal(zener)
      zi=dcmplx(0.d0,1.d0)

      do i=1,n
        do j=1,n
          zunit(i,j)=0.d0
        enddo
        zunit(i,i)=1.d0
      enddo

c     -----------------------------------------------------------------
c     -----------------------------------------------------------------
c     now calculate GF's from closed form
c     -----------------------------------------------------------------
c     -----------------------------------------------------------------
c     define zp
c     redefine gamma and delta
      do i=1,n
        do j=1,n
          ztinv(i,j)=zt(i,j)
          zs(i,j)=dconjg(zt(j,i))
          zsinv(i,j)=zs(i,j)
          ztmp1(i,j)=ener*zunit(i,j)-zu(i,j)
        enddo
      enddo
      call invers(ztinv,n,nx)
      call invers(zsinv,n,nx)
      call multiply(ztmp1,ztinv,zgamma,n,nx)
      do ir=1,n
        do is=1,n
          zp(ir,is)=0.d0
          zp(ir,is+n)=ztinv(ir,is)
          zp(ir+n,is)=-zs(ir,is)
          zp(ir+n,is+n)=zgamma(ir,is)
        enddo
      enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     diagonalise zp
      do ir=1,n2
        do is=1,n2
          ar(ir,is)=dreal(zp(ir,is))
          ai(ir,is)=dimag(zp(ir,is))
        enddo
      enddo 
      idoevec=1
      ifail=0
      call cg(2*nx,n2,ar,ai,rr,ri,idoevec,vr,vi,wk1,wk2,wk3,ifail)
      if(ifail.ne.0)then
        write(*,*)'SURFACENEW : ifail =',ifail
        ifail=1
        return
      endif
c
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     sort the |evals| into increasing order
c     and check eval(n+1)>eval(n)
c
c     define zdp
      do ir=1,n
        do is=1,n
          zdp(ir,is)=0.d0
          zdp(ir,is+n)=0.d0
          zdp(ir+n,is)=0.d0
          zdp(ir+n,is+n)=ztinv(ir,is)
        enddo
      enddo
      do ir=1,n2
        do is=1,n2
          zfoo1(ir,is)=dcmplx(vr(ir,is),vi(ir,is))
          zfoo2(ir,is)=dcmplx(vr(ir,is),vi(ir,is))
        enddo
      enddo
      call invers(zfoo2,n2,2*nx)

      nless=0
      nmore=0
      do ir=1,n2
        evlir=sqrt((rr(ir)**2)+(ri(ir)**2))
        zevlir=dcmplx(rr(ir),ri(ir))
        if(evlir.gt.1.d0+1.d-8)then
          evlab(ir)=evlir
        elseif(evlir.lt.1.d0-1.d-8)then
          evlab(ir)=evlir
        else !!! the eigenvalue lies on the unit circle .. check its derivative
          zdevl=0.d0
          do is=1,n2
            do it=1,n2
              zdevl=zdevl+zfoo2(ir,is)*zdp(is,it)*zfoo1(it,ir)
            enddo
          enddo
          zdkde=zdevl/(zevlir*zi)
          if(dimag(zdkde).gt.5.d-5)write(*,*)"ERROR:dimag(zdkde)=/ 0",
     $dimag(zdkde)
          evlab(ir)=evlir*exp(-dreal(zdkde)*1.d-8)
        endif
      enddo
      call sort(evlab,isort,n2)


CCC   evln1=evlab(isort(n+1))
CCC   evln=evlab(isort(n))
CCC   if(abs(evln1-evln).lt.1.d-16)then
CCC     write(*,*)' ERROR SURFACE : degenerate eigenstates '
CCC     stop
CCC   endif
c
c     -----------------------------------------------------------------
c     load the highest n eigenvectors
      do ir=1,n
        do is=1,n
          ik=isort(is+n)
          ztmp1(ir,is)=dcmplx(vr(ir,ik),vi(ir,ik))
          ztmp2(ir,is)=dcmplx(vr(ir+n,ik),vi(ir+n,ik))
        enddo
      enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     calculate L.H. zsurf
      call invers(ztmp2,n,nx)
      call multiply(ztmp1,ztmp2,zsurfl,n,nx)
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     calculate R.H. zsurf
      do ir=1,n2
        ik=isort(ir)
        zevl(ir)=dcmplx(rr(ik),ri(ik))
      enddo
      do ir=1,n
        do is=1,n
          ik=isort(is)
          ztmp1(ir,is)=dcmplx(vr(ir,ik),vi(ir,ik))
          ztmp2(ir,is)=ztmp1(ir,is)
        enddo
      enddo
      call invers(ztmp2,n,nx)
      do ir=1,n
        do is=1,n
          ztmp1(ir,is)=ztmp1(ir,is)*zevl(is)
        enddo
      enddo
      call multiply(ztmp1,ztmp2,ztmp3,n,nx)
      call multiply(ztmp3,zsinv,zsurfr,n,nx)
c
c     -----------------------------------------------------------------
c     calculate the RH SGF a little more accurately
C     do ir=1,9
C       do is=1,9
C         ztmp1(ir,is)=zsurfl(ir+9,is+9)
C         ztmp1(ir+9,is+9)=zsurfl(ir,is)
C         ztmp1(ir,is+9)=zsurfl(is,ir+9)
C         ztmp1(ir+9,is)=zsurfl(is+9,ir)
C       enddo
C     enddo
C     itran1(1)=1
C     itran1(2)=1
C     itran1(3)=1
C     itran1(4)=-1
C     itran1(5)=1
C     itran1(6)=-1
C     itran1(7)=-1
C     itran1(8)=1
C     itran1(9)=1
C     itran2(1)=1
C     itran2(2)=-1
C     itran2(3)=-1
C     itran2(4)=-1
C     itran2(5)=1
C     itran2(6)=1
C     itran2(7)=1
C     itran2(8)=1
C     itran2(9)=1
C     do ir=1,9
C       do is=1,9
C         ztmp1(ir,is)=itran1(ir)*itran1(is)*ztmp1(ir,is)
C         ztmp1(ir+9,is+9)=itran1(ir)*itran1(is)*ztmp1(ir+9,is+9)
C         ztmp1(ir,is+9)=itran2(ir)*itran2(is)*ztmp1(ir,is+9)
C         ztmp1(ir+9,is)=itran2(ir)*itran2(is)*ztmp1(ir+9,is)
C       enddo
C     enddo
C     do ir=1,n
C       do is=1,n
C         zsurfr(ir,is)=ztmp1(ir,is)
C       enddo
C     enddo
c     -----------------------------------------------------------------
      do ir=1,n
        do is=1,n
          ztmp1(ir,is)=zsurfl(ir,is)
          ztmp2(ir,is)=zsurfr(ir,is)
        enddo
      enddo
      call adlayer1(ztmp1,zu,zt,ztmp3,zener,n,nx)
      call adlayer1(ztmp2,zu,zs,ztmp3,zener,n,nx)
      xminl=0.d0
      xminr=0.d0
      do ir=1,n
        do is=1,n
          xmin1=cdabs(zsurfl(ir,is)-ztmp1(ir,is))
          xmin2=cdabs(zsurfr(ir,is)-ztmp2(ir,is))
          if(xmin1.gt.xminl)xminl=xmin1
          if(xmin2.gt.xminr)xminr=xmin2
        enddo
      enddo
      do ir=1,n
        do is=1,n
          zsurfl(ir,is)=ztmp1(ir,is)
          zsurfr(ir,is)=ztmp2(ir,is)
        enddo
      enddo
      xmin=max(xminl,xminr)
      if(xmin.gt.5.d-5)then
        ifail=1
      endif

      return
      end
c
c     *****************************************************************
c
      subroutine cg(nm,n,ar,ai,wr,wi,matz,vr,vi,fv1,fv2,fv3,ierr)
c
      integer n,nm,is1,is2,ierr,matz
      double precision ar(nm,n),ai(nm,n),wr(n),wi(n),vr(nm,n),vi(nm,n),
     x       fv1(n),fv2(n),fv3(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a complex general matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a=(ar,ai).
c
c        ar  and  ai  contain the real and imaginary parts,
c        respectively, of the complex general matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        wr  and  wi  contain the real and imaginary parts,
c        respectively, of the eigenvalues.
c
c        vr  and  vi  contain the real and imaginary parts,
c        respectively, of the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for comqr
c           and comqr2.  the normal completion code is zero.
c
c        fv1, fv2, and  fv3  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  cbal(nm,n,ar,ai,is1,is2,fv1)
      call  corth(nm,n,is1,is2,ar,ai,fv2,fv3)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  comqr(nm,n,is1,is2,ar,ai,wr,wi,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  comqr2(nm,n,is1,is2,fv2,fv3,ar,ai,wr,wi,vr,vi,ierr)
      if (ierr .ne. 0) go to 50
      call  cbabk2(nm,n,is1,is2,fv1,n,vr,vi)
   50 return
      end

c     -----------------------------------------------------------------
c     ----------------------------------
c     LAPACK routine in DXML
c     ----------------------------------

C       implicit double precision (a-h,o-y)
C       implicit complex*16 (z)
C       parameter (nplx=75,nslx=1,nspinx=10,nsubx=4,nsubatx=4)
C       parameter (nmatx=nslx*nsubx*nspinx,natomx=nslx*nspinx*nsubatx)
C       parameter (n2x=2*nmatx)
C 
C       parameter (nworkx=4*n2x)
C       dimension zwork(nworkx),rwork(nworkx)
C       dimension za(n2x,n2x),zv(n2x),zvl(n2x,n2x),zvr(n2x,n2x)
C 
C c
C       integer n,nm,is1,is2,ierr,matz
C       double precision ar(nm,n),ai(nm,n),wr(n),wi(n),vr(nm,n),vi(nm,n),
C      $       fv1(n),fv2(n),fv3(n)
C 
C 
C       if(nm.ne.n2x)then
C         write(*,*)'ERROR cg : nm /= n2x'
C         stop
C       endif
C 
C       do ir=1,n
C         do is=1,n
C           za(ir,is)=dcmplx(ar(ir,is),ai(ir,is))
C         enddo
C       enddo
C 
C       call zgeev('n','v',n,za,n2x,zv,zvl,n2x,zvr,n2x,zwork,
C      $  nworkx,rwork,ifail)
C 
C       do ir=1,n
C         do is=1,n
C           vr(ir,is)=dreal(zvr(ir,is))
C           vi(ir,is)=dimag(zvr(ir,is))
C         enddo
C         wr(ir)=dreal(zv(ir))
C         wi(ir)=dimag(zv(ir))
C       enddo
C       return
C       end
c
c     *****************************************************************
c
      subroutine cbabk2(nm,n,low,igh,scale,m,zr,zi)
c
      integer i,j,k,m,n,ii,nm,igh,low
      double precision scale(n),zr(nm,m),zi(nm,m)
      double precision s
c
c     this subroutine is a translation of the algol procedure
c     cbabk2, which is a complex version of balbak,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine forms the eigenvectors of a complex general
c     matrix by back transforming those of the corresponding
c     balanced matrix determined by  cbal.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by  cbal.
c
c        scale contains information determining the permutations
c          and scaling factors used by  cbal.
c
c        m is the number of eigenvectors to be back transformed.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors to be
c          back transformed in their first m columns.
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
c
      do 110 i = low, igh
         s = scale(i)
c     .......... left hand eigenvectors are back transformed
c                if the foregoing statement is replaced by
c                s=1.0d0/scale(i). ..........
         do 100 j = 1, m
            zr(i,j) = zr(i,j) * s
            zi(i,j) = zi(i,j) * s
  100    continue
c
  110 continue
c     .......... for i=low-1 step -1 until 1,
c                igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
c
         do 130 j = 1, m
            s = zr(i,j)
            zr(i,j) = zr(k,j)
            zr(k,j) = s
            s = zi(i,j)
            zi(i,j) = zi(k,j)
            zi(k,j) = s
  130    continue
c
  140 continue
c
  200 return
      end
c
c     *****************************************************************
c
      subroutine cbal(nm,n,ar,ai,low,igh,scale)
c
      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      double precision ar(nm,n),ai(nm,n),scale(n)
      double precision c,f,g,r,s,b2,radix
      logical noconv
c
c     this subroutine is a translation of the algol procedure
c     cbalance, which is a complex version of balance,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine balances a complex matrix and isolates
c     eigenvalues whenever possible.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex matrix to be balanced.
c
c     on output
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the balanced matrix.
c
c        low and igh are two integers such that ar(i,j) and ai(i,j)
c          are equal to zero if
c           (1) i is greater than j and
c           (2) j=1,...,low-1 or i=igh+1,...,n.
c
c        scale contains information determining the
c           permutations and scaling factors used.
c
c     suppose that the principal submatrix in rows low through igh
c     has been balanced, that p(j) denotes the index interchanged
c     with j during the permutation step, and that the elements
c     of the diagonal matrix used are denoted by d(i,j).  then
c        scale(j) = p(j),    for j = 1,...,low-1
c                 = d(j,j)       j = low,...,igh
c                 = p(j)         j = igh+1,...,n.
c     the order in which the interchanges are made is n to igh+1,
c     then 1 to low-1.
c
c     note that 1 is returned for igh if igh is zero formally.
c
c     the algol procedure exc contained in cbalance appears in
c     cbal  in line.  (note that the algol roles of identifiers
c     k,l have been reversed.)
c
c     arithmetic is real throughout.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      radix = 16.0d0
c
      b2 = radix * radix
      k = 1
      l = n
      go to 100
c     .......... in-line procedure for row and
c                column exchange ..........
   20 scale(m) = j
      if (j .eq. m) go to 50
c
      do 30 i = 1, l
         f = ar(i,j)
         ar(i,j) = ar(i,m)
         ar(i,m) = f
         f = ai(i,j)
         ai(i,j) = ai(i,m)
         ai(i,m) = f
   30 continue
c
      do 40 i = k, n
         f = ar(j,i)
         ar(j,i) = ar(m,i)
         ar(m,i) = f
         f = ai(j,i)
         ai(j,i) = ai(m,i)
         ai(m,i) = f
   40 continue
c
   50 go to (80,130), iexc
c     .......... search for rows isolating an eigenvalue
c                and push them down ..........
   80 if (l .eq. 1) go to 280
      l = l - 1
c     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l
         j = l + 1 - jj
c
         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (ar(j,i) .ne. 0.0d0 .or. ai(j,i) .ne. 0.0d0) go to 120
  110    continue
c
         m = l
         iexc = 1
         go to 20
  120 continue
c
      go to 140
c     .......... search for columns isolating an eigenvalue
c                and push them left ..........
  130 k = k + 1
c
  140 do 170 j = k, l
c
         do 150 i = k, l
            if (i .eq. j) go to 150
            if (ar(i,j) .ne. 0.0d0 .or. ai(i,j) .ne. 0.0d0) go to 170
  150    continue
c
         m = k
         iexc = 2
         go to 20
  170 continue
c     .......... now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 scale(i) = 1.0d0
c     .......... iterative loop for norm reduction ..........
  190 noconv = .false.
c
      do 270 i = k, l
         c = 0.0d0
         r = 0.0d0
c
         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + abs(ar(j,i)) + abs(ai(j,i))
            r = r + abs(ar(i,j)) + abs(ai(i,j))
  200    continue
c     .......... guard against zero c or r due to underflow ..........
         if (c .eq. 0.0d0 .or. r .eq. 0.0d0) go to 270
         g = r / radix
         f = 1.0d0
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
c     .......... now balance ..........
  240    if ((c + r) / f .ge. 0.95d0 * s) go to 270
         g = 1.0d0 / f
         scale(i) = scale(i) * f
         noconv = .true.
c
         do 250 j = k, n
            ar(i,j) = ar(i,j) * g
            ai(i,j) = ai(i,j) * g
  250    continue
c
         do 260 j = 1, l
            ar(j,i) = ar(j,i) * f
            ai(j,i) = ai(j,i) * f
  260    continue
c
  270 continue
c
      if (noconv) go to 190
c
  280 low = k
      igh = l
      return
      end
c
c     *****************************************************************
c
      subroutine cdiv(ar,ai,br,bi,cr,ci)
      double precision ar,ai,br,bi,cr,ci
c
c     complex division, (cr,ci) = (ar,ai)/(br,bi)
c
      double precision s,ars,ais,brs,bis
      s = abs(br) + abs(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2 + bis**2
      cr = (ars*brs + ais*bis)/s
      ci = (ais*brs - ars*bis)/s
      return
      end
c
c     *****************************************************************
c
      subroutine comqr(nm,n,low,igh,hr,hi,wr,wi,ierr)
c
      integer i,j,l,n,en,ll,nm,igh,itn,its,low,lp1,enm1,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag
c
c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr, num. math. 12, 369-376(1968) by martin
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
c
c     this subroutine finds the eigenvalues of a complex
c     upper hessenberg matrix by the qr method.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain
c          information about the unitary transformations used in
c          the reduction by  corth, if performed.
c
c     on output
c
c        the upper hessenberg portions of hr and hi have been
c          destroyed.  therefore, they must be saved before
c          calling  comqr  if subsequent calculation of
c          eigenvectors is to be performed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (low .eq. igh) go to 180
c     .......... create real subdiagonal elements ..........
      l = low + 1
c
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0d0) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0d0
c
         do 155 j = i, igh
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
c
         do 160 j = low, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
c
  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0d0
      ti = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 1001
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low d0 -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1))
     x            + abs(hr(l,l)) + abs(hi(l,l))
         tst2 = tst1 + abs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0d0
      yi = (hi(enm1,enm1) - si) / 2.0d0
      call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
      si = 0.0d0
c
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
c
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0d0
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0d0
         hi(i,i-1) = sr / norm
c
         do 490 j = i, en
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
c
  500 continue
c
      si = hi(en,en)
      if (si .eq. 0.0d0) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0d0
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
c
         do 580 i = l, j
            yr = hr(i,j-1)
            yi = 0.0d0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
c
  600 continue
c
      if (si .eq. 0.0d0) go to 240
c
      do 630 i = l, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
c
      go to 240
c     .......... a root found ..........
  660 wr(en) = hr(en,en) + tr
      wi(en) = hi(en,en) + ti
      en = enm1
      go to 220
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
c
c     *****************************************************************
c
      subroutine comqr2(nm,n,low,igh,ortr,orti,hr,hi,wr,wi,zr,zi,ierr)
C  MESHED overflow control WITH vectors of isolated roots (10/19/89 BSG)
C  MESHED overflow control WITH triangular multiply (10/30/89 BSG)
c
      integer i,j,k,l,m,n,en,ii,jj,ll,nm,nn,igh,ip1,
     x        itn,its,low,lp1,enm1,iend,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),
     x       ortr(igh),orti(igh)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag
c
c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr2, num. math. 16, 181-204(1970) by peters
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a complex upper hessenberg matrix by the qr
c     method.  the eigenvectors of a complex general matrix
c     can also be found if  corth  has been used to reduce
c     this general matrix to hessenberg form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        ortr and orti contain information about the unitary trans-
c          formations used in the reduction by  corth, if performed.
c          only elements low through igh are used.  if the eigenvectors
c          of the hessenberg matrix are desired, set ortr(j) and
c          orti(j) to 0.0d0 for these elements.
c
c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain further
c          information about the transformations which were used in the
c          reduction by  corth, if performed.  if the eigenvectors of
c          the hessenberg matrix are desired, these elements may be
c          arbitrary.
c
c     on output
c
c        ortr, orti, and the upper hessenberg portions of hr and hi
c          have been destroyed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors.  the eigenvectors
c          are unnormalized.  if an error exit is made, none of
c          the eigenvectors has been found.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated october 1989.
c
c     ------------------------------------------------------------------
c
      ierr = 0
c     .......... initialize eigenvector matrix ..........
      do 101 j = 1, n
c
         do 100 i = 1, n
            zr(i,j) = 0.0d0
            zi(i,j) = 0.0d0
  100    continue
         zr(j,j) = 1.0d0
  101 continue
c     .......... form the matrix of accumulated transformations
c                from the information left by corth ..........
      iend = igh - low - 1
      if (iend) 180, 150, 105
c     .......... for i=igh-1 step -1 until low+1 do -- ..........
  105 do 140 ii = 1, iend
         i = igh - ii
         if (ortr(i) .eq. 0.0d0 .and. orti(i) .eq. 0.0d0) go to 140
         if (hr(i,i-1) .eq. 0.0d0 .and. hi(i,i-1) .eq. 0.0d0) go to 140
c     .......... norm below is negative of h formed in corth ..........
         norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)
         ip1 = i + 1
c
         do 110 k = ip1, igh
            ortr(k) = hr(k,i-1)
            orti(k) = hi(k,i-1)
  110    continue
c
         do 130 j = i, igh
            sr = 0.0d0
            si = 0.0d0
c
            do 115 k = i, igh
               sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
               si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
  115       continue
c
            sr = sr / norm
            si = si / norm
c
            do 120 k = i, igh
               zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
               zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
  120       continue
c
  130    continue
c
  140 continue
c     .......... create real subdiagonal elements ..........
  150 l = low + 1
c
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0d0) go to 170
         norm = pythag(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0d0
c
         do 155 j = i, n
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
c
         do 160 j = 1, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
c
         do 165 j = low, igh
            si = yr * zi(j,i) + yi * zr(j,i)
            zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
            zi(j,i) = si
  165    continue
c
  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0d0
      ti = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 680
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1))
     x            + abs(hr(l,l)) + abs(hi(l,l))
         tst2 = tst1 + abs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0d0
      yi = (hi(enm1,enm1) - si) / 2.0d0
      call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
      si = 0.0d0
c
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
c
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0d0
         norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0d0
         hi(i,i-1) = sr / norm
c
         do 490 j = i, n
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
c
  500 continue
c
      si = hi(en,en)
      if (si .eq. 0.0d0) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0d0
      if (en .eq. n) go to 540
      ip1 = en + 1
c
      do 520 j = ip1, n
         yr = hr(en,j)
         yi = hi(en,j)
         hr(en,j) = sr * yr + si * yi
         hi(en,j) = sr * yi - si * yr
  520 continue
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
c
         do 580 i = 1, j
            yr = hr(i,j-1)
            yi = 0.0d0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
c
         do 590 i = low, igh
            yr = zr(i,j-1)
            yi = zi(i,j-1)
            zzr = zr(i,j)
            zzi = zi(i,j)
            zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
            zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  590    continue
c
  600 continue
c
      if (si .eq. 0.0d0) go to 240
c
      do 630 i = 1, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
c
      do 640 i = low, igh
         yr = zr(i,en)
         yi = zi(i,en)
         zr(i,en) = sr * yr - si * yi
         zi(i,en) = sr * yi + si * yr
  640 continue
c
      go to 240
c     .......... a root found ..........
  660 hr(en,en) = hr(en,en) + tr
      wr(en) = hr(en,en)
      hi(en,en) = hi(en,en) + ti
      wi(en) = hi(en,en)
      en = enm1
      go to 220
c     .......... all roots found.  backsubstitute to find
c                vectors of upper triangular form ..........
  680 norm = 0.0d0
c
      do 720 i = 1, n
c
         do 720 j = i, n
            tr = abs(hr(i,j)) + abs(hi(i,j))
            if (tr .gt. norm) norm = tr
  720 continue
c
      if (n .eq. 1 .or. norm .eq. 0.0d0) go to 1001
c     .......... for en=n step -1 until 2 do -- ..........
      do 800 nn = 2, n
         en = n + 2 - nn
         xr = wr(en)
         xi = wi(en)
         hr(en,en) = 1.0d0
         hi(en,en) = 0.0d0
         enm1 = en - 1
c     .......... for i=en-1 step -1 until 1 do -- ..........
         do 780 ii = 1, enm1
            i = en - ii
            zzr = 0.0d0
            zzi = 0.0d0
            ip1 = i + 1
c
            do 740 j = ip1, en
               zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
               zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
  740       continue
c
            yr = xr - wr(i)
            yi = xi - wi(i)
            if (yr .ne. 0.0d0 .or. yi .ne. 0.0d0) go to 765
               tst1 = norm
               yr = tst1
  760          yr = 0.01d0 * yr
               tst2 = norm + yr
               if (tst2 .gt. tst1) go to 760
  765       continue
            call cdiv(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
c     .......... overflow control ..........
            tr = abs(hr(i,en)) + abs(hi(i,en))
            if (tr .eq. 0.0d0) go to 780
            tst1 = tr
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 780
            do 770 j = i, en
               hr(j,en) = hr(j,en)/tr
               hi(j,en) = hi(j,en)/tr
  770       continue
c
  780    continue
c
  800 continue
c     .......... end backsubstitution ..........
c     .......... vectors of isolated roots ..........
      do  840 i = 1, N
         if (i .ge. low .and. i .le. igh) go to 840
c
         do 820 j = I, n
            zr(i,j) = hr(i,j)
            zi(i,j) = hi(i,j)
  820    continue
c
  840 continue
c     .......... multiply by transformation matrix to give
c                vectors of original full matrix.
c                for j=n step -1 until low do -- ..........
      do 880 jj = low, N
         j = n + low - jj
         m = min0(j,igh)
c
         do 880 i = low, igh
            zzr = 0.0d0
            zzi = 0.0d0
c
            do 860 k = low, m
               zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
               zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
  860       continue
c
            zr(i,j) = zzr
            zi(i,j) = zzi
  880 continue
c
      go to 1001
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
c
c     *****************************************************************
c
      subroutine corth(nm,n,low,igh,ar,ai,ortr,orti)
c
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      double precision ar(nm,n),ai(nm,n),ortr(igh),orti(igh)
      double precision f,g,h,fi,fr,scale,pythag
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure orthes, num. math. 12, 349-368(1968)
c     by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a complex general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     unitary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex input matrix.
c
c     on output
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the hessenberg matrix.  information
c          about the unitary transformations used in the reduction
c          is stored in the remaining triangles under the
c          hessenberg matrix.
c
c        ortr and orti contain further information about the
c          transformations.  only elements low through igh are used.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c
      do 180 m = kp1, la
         h = 0.0d0
         ortr(m) = 0.0d0
         orti(m) = 0.0d0
         scale = 0.0d0
c     .......... scale column (algol tol then not needed) ..........
         do 90 i = m, igh
   90    scale = scale + abs(ar(i,m-1)) + abs(ai(i,m-1))
c
         if (scale .eq. 0.0d0) go to 180
         mp = m + igh
c     .......... for i=igh step -1 until m do -- ..........
         do 100 ii = m, igh
            i = mp - ii
            ortr(i) = ar(i,m-1) / scale
            orti(i) = ai(i,m-1) / scale
            h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
  100    continue
c
         g = dsqrt(h)
         f = pythag(ortr(m),orti(m))
         if (f .eq. 0.0d0) go to 103
         h = h + f * g
         g = g / f
         ortr(m) = (1.0d0 + g) * ortr(m)
         orti(m) = (1.0d0 + g) * orti(m)
         go to 105
c
  103    ortr(m) = g
         ar(m,m-1) = scale
c     .......... form (i-(u*ut)/h) * a ..........
  105    do 130 j = m, n
            fr = 0.0d0
            fi = 0.0d0
c     .......... for i=igh step -1 until m do -- ..........
            do 110 ii = m, igh
               i = mp - ii
               fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
               fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
  110       continue
c
            fr = fr / h
            fi = fi / h
c
            do 120 i = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)
               ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)
  120       continue
c
  130    continue
c     .......... form (i-(u*ut)/h)*a*(i-(u*ut)/h) ..........
         do 160 i = 1, igh
            fr = 0.0d0
            fi = 0.0d0
c     .......... for j=igh step -1 until m do -- ..........
            do 140 jj = m, igh
               j = mp - jj
               fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
               fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
  140       continue
c
            fr = fr / h
            fi = fi / h
c
            do 150 j = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)
               ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)
  150       continue
c
  160    continue
c
         ortr(m) = scale * ortr(m)
         orti(m) = scale * orti(m)
         ar(m,m-1) = -g * ar(m,m-1)
         ai(m,m-1) = -g * ai(m,m-1)
  180 continue
c
  200 return
      end
c
c     *****************************************************************
c
      subroutine csroot(xr,xi,yr,yi)
      double precision xr,xi,yr,yi
c
c     (yr,yi) = complex dsqrt(xr,xi) 
c     branch chosen so that yr .ge. 0.0 and sign(yi) .eq. sign(xi)
c
      double precision s,tr,ti,pythag
      tr = xr
      ti = xi
      s = dsqrt(0.5d0*(pythag(tr,ti) + abs(tr)))
      if (tr .ge. 0.0d0) yr = s
      if (ti .lt. 0.0d0) s = -s
      if (tr .le. 0.0d0) yi = s
      if (tr .lt. 0.0d0) yr = 0.5d0*(ti/yi)
      if (tr .gt. 0.0d0) yi = 0.5d0*(ti/yr)
      return
      end
c
c     *****************************************************************
c
      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = max(abs(a),abs(b))
      if (p .eq. 0.0d0) go to 20
      r = (min(abs(a),abs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
c
c     *****************************************************************
c
      subroutine surfacenew16(zu,zt,zener,zsurfl,zsurfr,n,nx,ifail)
c
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      complex*16 zu(nx,nx),zt(nx,nx),zener
      complex*16 zsurfl(nx,nx),zsurfr(nx,nx)
      real*16 xur(nx,nx),xtr(nx,nx)
      real*16 xui(nx,nx),xti(nx,nx)
      real*16 xsurflr(nx,nx),xsurfli(nx,nx)
      real*16 xsurfrr(nx,nx),xsurfri(nx,nx)
      real*16 xenerr,xeneri
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      n2=2*n

      do ir=1,n
        do is=1,n
          xur(ir,is)=dreal(zu(ir,is))
          xui(ir,is)=dimag(zu(ir,is))
          xtr(ir,is)=dreal(zt(ir,is))
          xti(ir,is)=dimag(zt(ir,is))
        enddo
      enddo
      xenerr=dreal(zener)
      xeneri=dimag(zener)

      call surfsurf16at(xur,xui,xtr,xti,xenerr,xeneri,xsurflr,xsurfli,
     $  xsurfrr,xsurfri,n,nx,ifail)

      do ir=1,n
        do is=1,n
          zsurfl(ir,is)=dcmplx(xsurflr(ir,is),xsurfli(ir,is))
          zsurfr(ir,is)=dcmplx(xsurfrr(ir,is),xsurfri(ir,is))
        enddo
      enddo
      return
      end
c     *****************************************************************
c     *****************************************************************
c     *****************************************************************
c     *****************************************************************
c     *****************************************************************
c     *****************************************************************
C     Quadruple Precision Routines
C     QUADRUPLE PRECISION ROUTINES
c     *****************************************************************
c     *****************************************************************
c     *****************************************************************
c     *****************************************************************
c     *****************************************************************
c     *****************************************************************
c     *****************************************************************
c     *****************************************************************

      subroutine surfsurf16at(xur,xui,xtr,xti,xenerr,xeneri,xsurflr,
     $  xsurfli,xsurfrr,xsurfri,n,nx,ifail)
      implicit real*16 (a-h,o-y)
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      dimension xgammar(nx,nx),xgammai(nx,nx)
      dimension xunit(nx,nx)
      dimension xtmp1r(nx,nx),xtmp2r(nx,nx)
      dimension xtmp1i(nx,nx),xtmp2i(nx,nx)
      dimension xtmp3r(nx,nx)
      dimension xtmp3i(nx,nx)
      dimension xur(nx,nx),xtr(nx,nx)
      dimension xui(nx,nx),xti(nx,nx)
      dimension xtinvr(nx,nx),xsinvr(nx,nx)
      dimension xtinvi(nx,nx),xsinvi(nx,nx)
      dimension xsr(nx,nx),xsi(nx,nx)
      dimension xpr(2*nx,2*nx),xpi(2*nx,2*nx)
      dimension xdpr(2*nx,2*nx),xdpi(2*nx,2*nx)
      dimension xsurflr(nx,nx),xsurfrr(nx,nx)
      dimension xsurfli(nx,nx),xsurfri(nx,nx)
      dimension isort(2*nx)
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      dimension wk1(2*nx),wk2(2*nx),wk3(2*nx)
      dimension rr(2*nx),ri(2*nx),evlab(2*nx)
      dimension xevlr(2*nx),xevli(2*nx),evlphi(2*nx)
      dimension vr(2*nx,2*nx),vi(2*nx,2*nx)
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      dimension xfoo1r(2*nx,2*nx),xfoo2r(2*nx,2*nx)
      dimension xfoo1i(2*nx,2*nx),xfoo2i(2*nx,2*nx)
      dimension itran1(nx),itran2(nx)
      common/data/nmat,nins,mlay,natom,nspin,nlay,nsub,numnn
c     -----------------------------------------------------------------
      n2=2*n

      do i=1,n
        do j=1,n
          xunit(i,j)=0.d0
        enddo
        xunit(i,i)=1.d0
      enddo

c     -----------------------------------------------------------------
c     -----------------------------------------------------------------
c     now calculate GF's from closed form
c     -----------------------------------------------------------------
c     -----------------------------------------------------------------
c     define xp
c     redefine gamma and delta
      do i=1,n
        do j=1,n
          xtinvr(i,j)=xtr(i,j)
          xtinvi(i,j)=xti(i,j)
          xsr(i,j)=xtr(j,i)
          xsi(i,j)=-xti(j,i)
          xsinvr(i,j)=xsr(i,j)
          xsinvi(i,j)=xsi(i,j)
          xtmp1r(i,j)=xenerr*xunit(i,j)-xur(i,j)
          xtmp1i(i,j)=0.d0*xunit(i,j)-xui(i,j)
        enddo
      enddo
      call invers16(xtinvr,xtinvi,n,nx)
      call invers16(xsinvr,xsinvi,n,nx)
      call multiply16(xtmp1r,xtmp1i,xtinvr,xtinvi,xgammar,xgammai,n,
     $  nx)
      do ir=1,n
        do is=1,n
          xpr(ir,is)=0.d0
          xpi(ir,is)=0.d0
          xpr(ir,is+n)=xtinvr(ir,is)
          xpi(ir,is+n)=xtinvi(ir,is)
          xpr(ir+n,is)=-xsr(ir,is)
          xpi(ir+n,is)=-xsi(ir,is)
          xpr(ir+n,is+n)=xgammar(ir,is)
          xpi(ir+n,is+n)=xgammai(ir,is)
        enddo
      enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     diagonalise xp
      idoevec=1
      ifail=0
      call cg16(2*nx,n2,xpr,xpi,rr,ri,idoevec,vr,vi,wk1,wk2,wk3,
     $ifail)
      if(ifail.ne.0)then
        write(*,*)'SURFACENEW : ifail =',ifail
        ifail=1
        return
      endif
c
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     sort the |evals| into increasing order
c     and check eval(n+1)>eval(n)
c
c     define zdp
      do ir=1,n
        do is=1,n
          xdpr(ir,is)=0.d0
          xdpr(ir,is+n)=0.d0
          xdpr(ir+n,is)=0.d0
          xdpr(ir+n,is+n)=xtinvr(ir,is)
          xdpi(ir,is)=0.d0
          xdpi(ir,is+n)=0.d0
          xdpi(ir+n,is)=0.d0
          xdpi(ir+n,is+n)=xtinvi(ir,is)
        enddo
      enddo
      do ir=1,n2
        do is=1,n2
          xfoo1r(ir,is)=vr(ir,is)
          xfoo1i(ir,is)=vi(ir,is)
          xfoo2r(ir,is)=vr(ir,is)
          xfoo2i(ir,is)=vi(ir,is)
        enddo
      enddo
      call invers16(xfoo2r,xfoo2i,n2,2*nx)

      nless=0
      nmore=0
      do ir=1,n2
        evlir=sqrt((rr(ir)**2)+(ri(ir)**2))

CCC     zevlir=dcmplx(rr(ir),ri(ir))
        xevlirr=rr(ir)
        xevliri=ri(ir)

        if(evlir.gt.1.d0+1.d-8)then
          evlab(ir)=evlir
        elseif(evlir.lt.1.d0-1.d-8)then
          evlab(ir)=evlir
        else !!! the eigenvalue lies on the unit circle .. check its derivative
          xdevlr=0.d0
          xdevli=0.d0
          do is=1,n2
            do it=1,n2
CCCC          zdevl=zdevl+zfoo2(ir,is)*zdp(is,it)*zfoo1(it,ir)
              xdevlr=xdevlr+xfoo2r(ir,is)*
     $(xdpr(is,it)*xfoo1r(it,ir)-xdpi(is,it)*xfoo1i(it,ir))-
     $xfoo2i(ir,is)*
     $(xdpi(is,it)*xfoo1r(it,ir)+xdpr(is,it)*xfoo1i(it,ir))

              xdevli=xdevli+xfoo2i(ir,is)*
     $(xdpr(is,it)*xfoo1r(it,ir)-xdpi(is,it)*xfoo1i(it,ir))+
     $xfoo2r(ir,is)*
     $(xdpi(is,it)*xfoo1r(it,ir)+xdpr(is,it)*xfoo1i(it,ir))
            enddo
          enddo
CCC       zdkde=zdevl/(zevlir*zi)
          xdkder=(xdevli*xevlirr-xdevlr*xevliri)/(evlir**2)
          xdkdei=-(xdevlr*xevlirr+xdevli*xevliri)/(evlir**2)


          if(xdkdei.gt.5.d-5)write(*,*)"ERROR:dimag(zdkde)=/ 0",
     $xdkdei
          evlab(ir)=evlir*exp(-xdkder*1.d-8)
        endif
      enddo
      call sort16(evlab,isort,n2)

CCC   do ir=1,n2
CCC     evlab(ir)=sqrt((rr(ir)**2)+(ri(ir)**2))
CCC     isort(ir)=ir
CCC   enddo
CCC   call sort16(evlab,isort,n2)
CCC   do ir=1,n2
CCC     evlphi(ir)=arg16(xevlr(ir),xevli(ir))
CCC     evlab(ir)=sqrt(xevlr(ir)*xevlr(ir)+xevli(ir)*xevli(ir))
CCC   enddo
c
c     -----------------------------------------------------------------
c     load the highest n eigenvectors
      do ir=1,n
        do is=1,n
          ik=isort(is+n)
          xtmp1r(ir,is)=vr(ir,ik)
          xtmp1i(ir,is)=vi(ir,ik)
          xtmp2r(ir,is)=vr(ir+n,ik)
          xtmp2i(ir,is)=vi(ir+n,ik)
        enddo
      enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     calculate L.H. xsurf
      call invers16(xtmp2r,xtmp2i,n,nx)
      call multiply16(xtmp1r,xtmp1i,xtmp2r,xtmp2i,xsurflr,xsurfli,n,
     $  nx)
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c     calculate R.H. xsurf
      do ir=1,n2
        ik=isort(ir)
        xevlr(ir)=rr(ik)
        xevli(ir)=ri(ik)
      enddo
      do ir=1,n
        do is=1,n
          ik=isort(is)
          xtmp3r(ir,is)=vr(ir,ik)
          xtmp3i(ir,is)=vi(ir,ik)
          xtmp2r(ir,is)=vr(ir,ik)
          xtmp2i(ir,is)=vi(ir,ik)
        enddo
      enddo
      call invers16(xtmp2r,xtmp2i,n,nx)
      do ir=1,n
        do is=1,n
          xtmp1r(ir,is)=xtmp3r(ir,is)*xevlr(is)-xtmp3i(ir,is)*xevli(is)
          xtmp1i(ir,is)=xtmp3r(ir,is)*xevli(is)+xtmp3i(ir,is)*xevlr(is)
        enddo
      enddo
      call multiply16(xtmp1r,xtmp1i,xtmp2r,xtmp2i,xtmp3r,xtmp3i,
     $  n,nx)
      call multiply16(xtmp3r,xtmp3i,xsinvr,xsinvi,xsurfrr,xsurfri,
     $  n,nx)
c
c     -----------------------------------------------------------------
c     calculate the RH SGF a little more accurately
C     do ir=1,9
C       do is=1,9
C         xtmp1r(ir,is)=xsurflr(ir+9,is+9)
C         xtmp1i(ir,is)=xsurfli(ir+9,is+9)
C         xtmp1r(ir+9,is+9)=xsurflr(ir,is)
C         xtmp1i(ir+9,is+9)=xsurfli(ir,is)
C         xtmp1r(ir,is+9)=xsurflr(is,ir+9)
C         xtmp1i(ir,is+9)=xsurfli(is,ir+9)
C         xtmp1r(ir+9,is)=xsurflr(is+9,ir)
C         xtmp1i(ir+9,is)=xsurfli(is+9,ir)
C       enddo
C     enddo
C     itran1(1)=1
C     itran1(2)=1
C     itran1(3)=1
C     itran1(4)=-1
C     itran1(5)=1
C     itran1(6)=-1
C     itran1(7)=-1
C     itran1(8)=1
C     itran1(9)=1
C     itran2(1)=1
C     itran2(2)=-1
C     itran2(3)=-1
C     itran2(4)=-1
C     itran2(5)=1
C     itran2(6)=1
C     itran2(7)=1
C     itran2(8)=1
C     itran2(9)=1
C     do ir=1,9
C       do is=1,9
C         xtmp1r(ir,is)=itran1(ir)*itran1(is)*xtmp1r(ir,is)
C         xtmp1i(ir,is)=itran1(ir)*itran1(is)*xtmp1i(ir,is)
C         xtmp1r(ir+9,is+9)=itran1(ir)*itran1(is)*xtmp1r(ir+9,is+9)
C         xtmp1i(ir+9,is+9)=itran1(ir)*itran1(is)*xtmp1i(ir+9,is+9)
C         xtmp1r(ir,is+9)=itran2(ir)*itran2(is)*xtmp1r(ir,is+9)
C         xtmp1i(ir,is+9)=itran2(ir)*itran2(is)*xtmp1i(ir,is+9)
C         xtmp1r(ir+9,is)=itran2(ir)*itran2(is)*xtmp1r(ir+9,is)
C         xtmp1i(ir+9,is)=itran2(ir)*itran2(is)*xtmp1i(ir+9,is)
C       enddo
C     enddo
C     do ir=1,n
C       do is=1,n
C         xsurfrr(ir,is)=xtmp1r(ir,is)
C         xsurfri(ir,is)=xtmp1i(ir,is)
C       enddo
C     enddo
c     -----------------------------------------------------------------
      do ir=1,n
        do is=1,n
          xtmp1r(ir,is)=xsurflr(ir,is)
          xtmp1i(ir,is)=xsurfli(ir,is)
          xtmp2r(ir,is)=xsurfrr(ir,is)
          xtmp2i(ir,is)=xsurfri(ir,is)
        enddo
      enddo
      call adlayer116(xtmp1r,xtmp1i,xur,xui,xtr,xti,xtmp3r,xtmp3i,
     $xenerr,xeneri,n,nx)
      call adlayer116(xtmp2r,xtmp2i,xur,xui,xsr,xsi,xtmp3r,xtmp3i,
     $xenerr,xeneri,n,nx)
      xminl=0.d0
      xminr=0.d0
      do ir=1,n
        do is=1,n
          xmin1=abs((xsurflr(ir,is)-xtmp1r(ir,is))**2+
     $      (xsurfli(ir,is)-xtmp1i(ir,is))**2)
          xmin2=abs((xsurfrr(ir,is)-xtmp2r(ir,is))**2+
     $      (xsurfri(ir,is)-xtmp2i(ir,is))**2)
          if(xmin1.gt.xminl)xminl=xmin1
          if(xmin2.gt.xminr)xminr=xmin2
        enddo
      enddo
      do ir=1,n
        do is=1,n
          xsurflr(ir,is)=xtmp1r(ir,is)
          xsurfli(ir,is)=xtmp1i(ir,is)
          xsurfrr(ir,is)=xtmp2r(ir,is)
          xsurfri(ir,is)=xtmp2i(ir,is)
        enddo
      enddo
      xmin=max(xminl,xminr)
      if(xmin.gt.5.d-5)then
        write(*,*)"surfsurf16 :",xmin
        ifail=1
      endif

      return
      end
c
c     *****************************************************************
c
      real*16 function arg16(xr,xi)
      implicit real*16 (a-h,o-y)
      pi=acos(-1.d0)
      eps=1.d-24
      xmod=sqrt(xr*xr+xi*xi)

      if(xmod.lt.eps)then
        write(*,*)' ERROR function arg16 : cdabs(z) = 0'
        stop
      endif

      if((xr.gt.0).and.(xi.ge.0))arg16=atan(xi/xr)
      if((xr.gt.0).and.(xi.le.0))arg16=atan(xi/xr)
      if((xr.lt.0).and.(xi.ge.0))arg16=atan(xi/xr)+pi
      if((xr.lt.0).and.(xi.le.0))arg16=atan(xi/xr)-pi

      if(xr.eq.0)then
        if(xi.gt.0)arg16=pi/2
        if(xi.lt.0)arg16=-pi/2
      endif
      return
      end
c
c     *****************************************************************
c
      subroutine adlayer116(xglr,xgli,xur,xui,xtr,xti,xwrkr,xwrki,
     $xenerr,xeneri,n,nx)
      implicit real*16 (a-h,o-y)

      dimension xglr(nx,nx),xur(nx,nx),xtr(nx,nx),xwrkr(nx,nx)
      dimension xgli(nx,nx),xui(nx,nx),xti(nx,nx),xwrki(nx,nx)

c     -----------------------------------------------------------------
c     adlayer ontop of gl
      call multiply16(xglr,xgli,xtr,xti,xwrkr,xwrki,n,nx)
      call hconjugate16(xtr,xti,n,nx)
      call multiply16(xtr,xti,xwrkr,xwrki,xglr,xgli,n,nx)
      call hconjugate16(xtr,xti,n,nx)

      do ir=1,n
        do is=1,ir-1
          xglr(ir,is)=-xur(ir,is)-xglr(ir,is)
          xgli(ir,is)=-xui(ir,is)-xgli(ir,is)
          xglr(is,ir)=-xur(is,ir)-xglr(is,ir)
          xgli(is,ir)=-xui(is,ir)-xgli(is,ir)
        enddo
        xglr(ir,ir)=xenerr-xur(ir,ir)-xglr(ir,ir)
        xgli(ir,ir)=xeneri-xui(ir,ir)-xgli(ir,ir)
      enddo
      call invers16(xglr,xgli,n,nx)
c     -----------------------------------------------------------------

      return
      end
c
c     *****************************************************************
c
      subroutine hconjugate16(xmatr,xmati,n,nx)
      implicit real*16 (a-h,o-y)
      dimension xmatr(nx,nx),xmati(nx,nx)
      do ir=1,n
        do is=1,ir-1
          xtmpr=xmatr(is,ir)
          xtmpi=xmati(is,ir)
          xmatr(is,ir)=xmatr(ir,is)
          xmati(is,ir)=-xmati(ir,is)
          xmatr(ir,is)=xtmpr
          xmati(ir,is)=-xtmpi
        enddo
        xmatr(ir,ir)=xmatr(ir,ir)
        xmati(ir,ir)=-xmati(ir,ir)
      enddo

      return
      end
c
c     *****************************************************************
c
      subroutine multiply16(xr,xi,yr,yi,xyr,xyi,nmat,nmatx)
      implicit real*16 (a-h,o-y)
      dimension xr(nmatx,nmatx),yr(nmatx,nmatx),xyr(nmatx,nmatx)
      dimension xi(nmatx,nmatx),yi(nmatx,nmatx),xyi(nmatx,nmatx)
      do ir=1,nmat
        do is=1,nmat
          xyr(ir,is)=0.d0
          xyi(ir,is)=0.d0
          do k=1,nmat
            xyr(ir,is)=xyr(ir,is)+xr(ir,k)*yr(k,is)-xi(ir,k)*yi(k,is)
            xyi(ir,is)=xyi(ir,is)+xr(ir,k)*yi(k,is)+xi(ir,k)*yr(k,is)
          enddo
        enddo
      enddo
c
      return
      end
c
c     *****************************************************************
c
      subroutine remultiply16(xr,yr,xyr,nmat,nmatx)
      implicit real*16 (a-h,o-y)
      dimension xr(nmatx,nmatx),yr(nmatx,nmatx),xyr(nmatx,nmatx)
      do ir=1,nmat
        do is=1,nmat
          xyr(ir,is)=0.d0
          do k=1,nmat
            xyr(ir,is)=xyr(ir,is)+xr(ir,k)*yr(k,is)
          enddo
        enddo
      enddo
c
      return
      end
c
c     *****************************************************************
c 
      subroutine invers16(ar,ai,n,nx) 
      real*16 ar(nx,nx),ai(nx,nx)
      real*16 amaxr,amaxi,saver,savei,cr,ci,xtmp,ytmp
      parameter (iworkx=10000)
      integer ik(iworkx),jk(iworkx)
      if(iworkx.lt.n)then
        write(*,*)' ERROR INVERS16 : iworkx < n'
        stop
      endif
      do 100 k=1,n 
        amaxr=0.d0 
        amaxi=0.d0 
21      do 30 i=k,n 
          do 30 j=k,n 
23          if((amaxr**2+amaxi**2)-(ar(i,j)**2+ai(i,j)**2))24,24,30 
24          amaxr=ar(i,j) 
            amaxi=ai(i,j) 
            ik(k)=i 
            jk(k)=j 
30      continue 
41      i=ik(k) 
        if(i-k)21,51,43 
43      do 50 j=1,n 
          saver=ar(k,j) 
          savei=ai(k,j) 
          ar(k,j)=ar(i,j) 
          ai(k,j)=ai(i,j) 
          ar(i,j)=-saver
50        ai(i,j)=-savei
51      j=jk(k) 
        if(j-k)21,61,53 
53      do 60 i=1,n 
          saver=ar(i,k) 
          savei=ai(i,k) 
          ar(i,k)=ar(i,j) 
          ai(i,k)=ai(i,j) 
          ar(i,j)=-saver
60        ai(i,j)=-savei
61      do 70 i=1,n 
          if(i-k)63,70,63 
63        call cdiv16(ar(i,k),ai(i,k),amaxr,amaxi,cr,ci)
          ar(i,k)=-cr
          ai(i,k)=-ci
70      continue 
71      do 80 i=1,n 
          do 80 j=1,n 
            if(i-k)74,80,74 
74          if(j-k)75,80,75 
75          ar(i,j)=ar(i,j)+(ar(i,k)*ar(k,j)-ai(i,k)*ai(k,j)) 
            ai(i,j)=ai(i,j)+(ai(i,k)*ar(k,j)+ar(i,k)*ai(k,j)) 
80      continue 
81      do 90 j=1,n 
          if(j-k)83,90,83 
83        call cdiv16(ar(k,j),ai(k,j),amaxr,amaxi,cr,ci)
          ar(k,j)=cr
          ai(k,j)=ci
90      continue 
        xtmp=1.d0
        ytmp=0.d0
        call cdiv16(xtmp,ytmp,amaxr,amaxi,cr,ci)
        ar(k,k)=cr
        ai(k,k)=ci
100   continue 
101   do 130 l=1,n 
        k=n-l+1 
        j=ik(k) 
        if(j-k)111,111,105 
105     do 110 i=1,n 
          saver=ar(i,k) 
          savei=ai(i,k) 
          ar(i,k)=-ar(i,j) 
          ai(i,k)=-ai(i,j) 
          ar(i,j)=saver
110       ai(i,j)=savei
111     i=jk(k) 
        if(i-k)130,130,1013
1013    do 120 j=1,n 
          saver=ar(k,j) 
          savei=ai(k,j) 
          ar(k,j)=-ar(i,j) 
          ai(k,j)=-ai(i,j) 
          ar(i,j)=saver
120       ai(i,j)=savei
130   continue 
140   return 
      end

c

c
c
c     *****************************************************************
c
      subroutine sort16(x,ka,l)
      implicit real*16 (a-h,o-z)
c
c     subroutine sorts the array x
c     the indices (ascending order) are given in array ka
c
      parameter (nvecx=100000)
      dimension x(l),ka(l),kb(nvecx)
      if(l.gt.nvecx)then
        write(*,*)' ERROR SORT.F : l > nvecx ',l,nvecx
        stop
      endif
      do 13 j=1,l
        kb(j)=j
13      ka(j)=j
      l1=1
200   j=1
20    limi1=l1*2*(j-1)+1
      if(limi1.gt.l) go to 28
      limi2=limi1+l1
      lims1=min(limi2-1,l)
      if(limi2.gt.l) go to 28
      lims2=min(limi2+l1-1,l)
      ir=limi1
      i1=ir
      i2=limi2
21    k1=ka(i1)
      k2=ka(i2)
      if(x(k1).le.x(k2)) go to 22
      kb(ir)=k2
      i2=i2+1
      ir=ir+1
      if(i2.le.lims2) go to 21
      go to 23
22    kb(ir)=k1
      i1 = i1+1
      ir=ir+1
      if(i1.le.lims1) go to 21
24    k2=ka(i2)
      kb(ir) = k2
      i2 = i2+1
      ir = ir+1
      if(i2.le.lims2) go to 24
      go to 25
23    k1=ka(i1)
      kb(ir) = k1
      i1=i1+1
      ir=ir+1
      if(i1.le.lims1) go to 23
25    j=j+1
      go to 20
28    do 280 k=1,l
280     ka(k)=kb(k)
      l1 = l1*2
      if(l1.lt.2*l) go to 200
      return
      end
c
c     *****************************************************************
c
      subroutine cg16(nm,n,ar,ai,wr,wi,matz,vr,vi,fv1,fv2,fv3,ierr)
c
      integer n,nm,is1,is2,ierr,matz
      real*16 ar(nm,n),ai(nm,n),wr(n),wi(n),vr(nm,n),vi(nm,n),
     x       fv1(n),fv2(n),fv3(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a complex general matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a=(ar,ai).
c
c        ar  and  ai  contain the real and imaginary parts,
c        respectively, of the complex general matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        wr  and  wi  contain the real and imaginary parts,
c        respectively, of the eigenvalues.
c
c        vr  and  vi  contain the real and imaginary parts,
c        respectively, of the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for comqr
c           and comqr2.  the normal completion code is zero.
c
c        fv1, fv2, and  fv3  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  cbal16(nm,n,ar,ai,is1,is2,fv1)
      call  corth16(nm,n,is1,is2,ar,ai,fv2,fv3)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  comqr16(nm,n,is1,is2,ar,ai,wr,wi,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  comqr216(nm,n,is1,is2,fv2,fv3,ar,ai,wr,wi,vr,vi,ierr)
      if (ierr .ne. 0) go to 50
      call  cbabk216(nm,n,is1,is2,fv1,n,vr,vi)
   50 return
      end

c
c     *****************************************************************
c
      subroutine cbabk216(nm,n,low,igh,scale,m,zr,zi)
c
      integer i,j,k,m,n,ii,nm,igh,low
      real*16 scale(n),zr(nm,m),zi(nm,m)
      real*16 s
c
c     this subroutine is a translation of the algol procedure
c     cbabk2, which is a complex version of balbak,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine forms the eigenvectors of a complex general
c     matrix by back transforming those of the corresponding
c     balanced matrix determined by  cbal.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by  cbal.
c
c        scale contains information determining the permutations
c          and scaling factors used by  cbal.
c
c        m is the number of eigenvectors to be back transformed.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors to be
c          back transformed in their first m columns.
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120
c
      do 110 i = low, igh
         s = scale(i)
c     .......... left hand eigenvectors are back transformed
c                if the foregoing statement is replaced by
c                s=1.0d0/scale(i). ..........
         do 100 j = 1, m
            zr(i,j) = zr(i,j) * s
            zi(i,j) = zi(i,j) * s
  100    continue
c
  110 continue
c     .......... for i=low-1 step -1 until 1,
c                igh+1 step 1 until n do -- ..........
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140
c
         do 130 j = 1, m
            s = zr(i,j)
            zr(i,j) = zr(k,j)
            zr(k,j) = s
            s = zi(i,j)
            zi(i,j) = zi(k,j)
            zi(k,j) = s
  130    continue
c
  140 continue
c
  200 return
      end
c
c     *****************************************************************
c
      subroutine cbal16(nm,n,ar,ai,low,igh,scale)
c
      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      real*16 ar(nm,n),ai(nm,n),scale(n)
      real*16 c,f,g,r,s,b2,radix
      logical noconv
c
c     this subroutine is a translation of the algol procedure
c     cbalance, which is a complex version of balance,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine balances a complex matrix and isolates
c     eigenvalues whenever possible.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex matrix to be balanced.
c
c     on output
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the balanced matrix.
c
c        low and igh are two integers such that ar(i,j) and ai(i,j)
c          are equal to zero if
c           (1) i is greater than j and
c           (2) j=1,...,low-1 or i=igh+1,...,n.
c
c        scale contains information determining the
c           permutations and scaling factors used.
c
c     suppose that the principal submatrix in rows low through igh
c     has been balanced, that p(j) denotes the index interchanged
c     with j during the permutation step, and that the elements
c     of the diagonal matrix used are denoted by d(i,j).  then
c        scale(j) = p(j),    for j = 1,...,low-1
c                 = d(j,j)       j = low,...,igh
c                 = p(j)         j = igh+1,...,n.
c     the order in which the interchanges are made is n to igh+1,
c     then 1 to low-1.
c
c     note that 1 is returned for igh if igh is zero formally.
c
c     the algol procedure exc contained in cbalance appears in
c     cbal  in line.  (note that the algol roles of identifiers
c     k,l have been reversed.)
c
c     arithmetic is real throughout.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      radix = 16.0d0
c
      b2 = radix * radix
      k = 1
      l = n
      go to 100
c     .......... in-line procedure for row and
c                column exchange ..........
   20 scale(m) = j
      if (j .eq. m) go to 50
c
      do 30 i = 1, l
         f = ar(i,j)
         ar(i,j) = ar(i,m)
         ar(i,m) = f
         f = ai(i,j)
         ai(i,j) = ai(i,m)
         ai(i,m) = f
   30 continue
c
      do 40 i = k, n
         f = ar(j,i)
         ar(j,i) = ar(m,i)
         ar(m,i) = f
         f = ai(j,i)
         ai(j,i) = ai(m,i)
         ai(m,i) = f
   40 continue
c
   50 go to (80,130), iexc
c     .......... search for rows isolating an eigenvalue
c                and push them down ..........
   80 if (l .eq. 1) go to 280
      l = l - 1
c     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l
         j = l + 1 - jj
c
         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (ar(j,i) .ne. 0.0d0 .or. ai(j,i) .ne. 0.0d0) go to 120
  110    continue
c
         m = l
         iexc = 1
         go to 20
  120 continue
c
      go to 140
c     .......... search for columns isolating an eigenvalue
c                and push them left ..........
  130 k = k + 1
c
  140 do 170 j = k, l
c
         do 150 i = k, l
            if (i .eq. j) go to 150
            if (ar(i,j) .ne. 0.0d0 .or. ai(i,j) .ne. 0.0d0) go to 170
  150    continue
c
         m = k
         iexc = 2
         go to 20
  170 continue
c     .......... now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 scale(i) = 1.0d0
c     .......... iterative loop for norm reduction ..........
  190 noconv = .false.
c
      do 270 i = k, l
         c = 0.0d0
         r = 0.0d0
c
         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + abs(ar(j,i)) + abs(ai(j,i))
            r = r + abs(ar(i,j)) + abs(ai(i,j))
  200    continue
c     .......... guard against zero c or r due to underflow ..........
         if (c .eq. 0.0d0 .or. r .eq. 0.0d0) go to 270
         g = r / radix
         f = 1.0d0
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
c     .......... now balance ..........
  240    if ((c + r) / f .ge. 0.95d0 * s) go to 270
         g = 1.0d0 / f
         scale(i) = scale(i) * f
         noconv = .true.
c
         do 250 j = k, n
            ar(i,j) = ar(i,j) * g
            ai(i,j) = ai(i,j) * g
  250    continue
c
         do 260 j = 1, l
            ar(j,i) = ar(j,i) * f
            ai(j,i) = ai(j,i) * f
  260    continue
c
  270 continue
c
      if (noconv) go to 190
c
  280 low = k
      igh = l
      return
      end
c
c     *****************************************************************
c
      subroutine cdiv16(ar,ai,br,bi,cr,ci)
      real*16 ar,ai,br,bi,cr,ci
c
c     complex division, (cr,ci) = (ar,ai)/(br,bi)
c
      real*16 s,ars,ais,brs,bis
      s = abs(br) + abs(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2 + bis**2
      cr = (ars*brs + ais*bis)/s
      ci = (ais*brs - ars*bis)/s
      return
      end
c
c     *****************************************************************
c
      subroutine comqr16(nm,n,low,igh,hr,hi,wr,wi,ierr)
c
      integer i,j,l,n,en,ll,nm,igh,itn,its,low,lp1,enm1,ierr
      real*16 hr(nm,n),hi(nm,n),wr(n),wi(n)
      real*16 si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag16
c
c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr, num. math. 12, 369-376(1968) by martin
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
c
c     this subroutine finds the eigenvalues of a complex
c     upper hessenberg matrix by the qr method.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain
c          information about the unitary transformations used in
c          the reduction by  corth, if performed.
c
c     on output
c
c        the upper hessenberg portions of hr and hi have been
c          destroyed.  therefore, they must be saved before
c          calling  comqr  if subsequent calculation of
c          eigenvectors is to be performed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag16 for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (low .eq. igh) go to 180
c     .......... create real subdiagonal elements ..........
      l = low + 1
c
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0d0) go to 170
         norm = pythag16(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0d0
c
         do 155 j = i, igh
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
c
         do 160 j = low, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
c
  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0d0
      ti = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 1001
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low d0 -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1))
     x            + abs(hr(l,l)) + abs(hi(l,l))
         tst2 = tst1 + abs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0d0
      yi = (hi(enm1,enm1) - si) / 2.0d0
      call csroot16(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv16(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
      si = 0.0d0
c
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
c
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0d0
         norm = pythag16(pythag16(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0d0
         hi(i,i-1) = sr / norm
c
         do 490 j = i, en
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
c
  500 continue
c
      si = hi(en,en)
      if (si .eq. 0.0d0) go to 540
      norm = pythag16(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0d0
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
c
         do 580 i = l, j
            yr = hr(i,j-1)
            yi = 0.0d0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
c
  600 continue
c
      if (si .eq. 0.0d0) go to 240
c
      do 630 i = l, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
c
      go to 240
c     .......... a root found ..........
  660 wr(en) = hr(en,en) + tr
      wi(en) = hi(en,en) + ti
      en = enm1
      go to 220
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
c
c     *****************************************************************
c
      subroutine comqr216(nm,n,low,igh,ortr,orti,hr,hi,wr,wi,zr,zi,ierr)
C  MESHED overflow control WITH vectors of isolated roots (10/19/89 BSG)
C  MESHED overflow control WITH triangular multiply (10/30/89 BSG)
c
      integer i,j,k,l,m,n,en,ii,jj,ll,nm,nn,igh,ip1,
     x        itn,its,low,lp1,enm1,iend,ierr
      real*16 hr(nm,n),hi(nm,n),wr(n),wi(n),zr(nm,n),zi(nm,n),
     x       ortr(igh),orti(igh)
      real*16 si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     x       pythag16
c
c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr2, num. math. 16, 181-204(1970) by peters
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a complex upper hessenberg matrix by the qr
c     method.  the eigenvectors of a complex general matrix
c     can also be found if  corth  has been used to reduce
c     this general matrix to hessenberg form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        ortr and orti contain information about the unitary trans-
c          formations used in the reduction by  corth, if performed.
c          only elements low through igh are used.  if the eigenvectors
c          of the hessenberg matrix are desired, set ortr(j) and
c          orti(j) to 0.0d0 for these elements.
c
c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain further
c          information about the transformations which were used in the
c          reduction by  corth, if performed.  if the eigenvectors of
c          the hessenberg matrix are desired, these elements may be
c          arbitrary.
c
c     on output
c
c        ortr, orti, and the upper hessenberg portions of hr and hi
c          have been destroyed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the eigenvectors.  the eigenvectors
c          are unnormalized.  if an error exit is made, none of
c          the eigenvectors has been found.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag16 for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated october 1989.
c
c     ------------------------------------------------------------------
c
      ierr = 0
c     .......... initialize eigenvector matrix ..........
      do 101 j = 1, n
c
         do 100 i = 1, n
            zr(i,j) = 0.0d0
            zi(i,j) = 0.0d0
  100    continue
         zr(j,j) = 1.0d0
  101 continue
c     .......... form the matrix of accumulated transformations
c                from the information left by corth ..........
      iend = igh - low - 1
      if (iend) 180, 150, 105
c     .......... for i=igh-1 step -1 until low+1 do -- ..........
  105 do 140 ii = 1, iend
         i = igh - ii
         if (ortr(i) .eq. 0.0d0 .and. orti(i) .eq. 0.0d0) go to 140
         if (hr(i,i-1) .eq. 0.0d0 .and. hi(i,i-1) .eq. 0.0d0) go to 140
c     .......... norm below is negative of h formed in corth ..........
         norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)
         ip1 = i + 1
c
         do 110 k = ip1, igh
            ortr(k) = hr(k,i-1)
            orti(k) = hi(k,i-1)
  110    continue
c
         do 130 j = i, igh
            sr = 0.0d0
            si = 0.0d0
c
            do 115 k = i, igh
               sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
               si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
  115       continue
c
            sr = sr / norm
            si = si / norm
c
            do 120 k = i, igh
               zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
               zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
  120       continue
c
  130    continue
c
  140 continue
c     .......... create real subdiagonal elements ..........
  150 l = low + 1
c
      do 170 i = l, igh
         ll = min0(i+1,igh)
         if (hi(i,i-1) .eq. 0.0d0) go to 170
         norm = pythag16(hr(i,i-1),hi(i,i-1))
         yr = hr(i,i-1) / norm
         yi = hi(i,i-1) / norm
         hr(i,i-1) = norm
         hi(i,i-1) = 0.0d0
c
         do 155 j = i, n
            si = yr * hi(i,j) - yi * hr(i,j)
            hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
            hi(i,j) = si
  155    continue
c
         do 160 j = 1, ll
            si = yr * hi(j,i) + yi * hr(j,i)
            hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
            hi(j,i) = si
  160    continue
c
         do 165 j = low, igh
            si = yr * zi(j,i) + yi * zr(j,i)
            zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
            zi(j,i) = si
  165    continue
c
  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
         if (i .ge. low .and. i .le. igh) go to 200
         wr(i) = hr(i,i)
         wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0d0
      ti = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 680
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
  240 do 260 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 300
         tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1))
     x            + abs(hr(l,l)) + abs(hi(l,l))
         tst2 = tst1 + abs(hr(l,l-1))
         if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0d0
      yi = (hi(enm1,enm1) - si) / 2.0d0
      call csroot16(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv16(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
      si = 0.0d0
c
  340 do 360 i = low, en
         hr(i,i) = hr(i,i) - sr
         hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
c
      do 500 i = lp1, en
         sr = hr(i,i-1)
         hr(i,i-1) = 0.0d0
         norm = pythag16(pythag16(hr(i-1,i-1),hi(i-1,i-1)),sr)
         xr = hr(i-1,i-1) / norm
         wr(i-1) = xr
         xi = hi(i-1,i-1) / norm
         wi(i-1) = xi
         hr(i-1,i-1) = norm
         hi(i-1,i-1) = 0.0d0
         hi(i,i-1) = sr / norm
c
         do 490 j = i, n
            yr = hr(i-1,j)
            yi = hi(i-1,j)
            zzr = hr(i,j)
            zzi = hi(i,j)
            hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
            hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
            hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
            hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490    continue
c
  500 continue
c
      si = hi(en,en)
      if (si .eq. 0.0d0) go to 540
      norm = pythag16(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0d0
      if (en .eq. n) go to 540
      ip1 = en + 1
c
      do 520 j = ip1, n
         yr = hr(en,j)
         yi = hi(en,j)
         hr(en,j) = sr * yr + si * yi
         hi(en,j) = sr * yi - si * yr
  520 continue
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
         xr = wr(j-1)
         xi = wi(j-1)
c
         do 580 i = 1, j
            yr = hr(i,j-1)
            yi = 0.0d0
            zzr = hr(i,j)
            zzi = hi(i,j)
            if (i .eq. j) go to 560
            yi = hi(i,j-1)
            hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580    continue
c
         do 590 i = low, igh
            yr = zr(i,j-1)
            yi = zi(i,j-1)
            zzr = zr(i,j)
            zzi = zi(i,j)
            zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
            zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
            zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
            zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  590    continue
c
  600 continue
c
      if (si .eq. 0.0d0) go to 240
c
      do 630 i = 1, en
         yr = hr(i,en)
         yi = hi(i,en)
         hr(i,en) = sr * yr - si * yi
         hi(i,en) = sr * yi + si * yr
  630 continue
c
      do 640 i = low, igh
         yr = zr(i,en)
         yi = zi(i,en)
         zr(i,en) = sr * yr - si * yi
         zi(i,en) = sr * yi + si * yr
  640 continue
c
      go to 240
c     .......... a root found ..........
  660 hr(en,en) = hr(en,en) + tr
      wr(en) = hr(en,en)
      hi(en,en) = hi(en,en) + ti
      wi(en) = hi(en,en)
      en = enm1
      go to 220
c     .......... all roots found.  backsubstitute to find
c                vectors of upper triangular form ..........
  680 norm = 0.0d0
c
      do 720 i = 1, n
c
         do 720 j = i, n
            tr = abs(hr(i,j)) + abs(hi(i,j))
            if (tr .gt. norm) norm = tr
  720 continue
c
      if (n .eq. 1 .or. norm .eq. 0.0d0) go to 1001
c     .......... for en=n step -1 until 2 do -- ..........
      do 800 nn = 2, n
         en = n + 2 - nn
         xr = wr(en)
         xi = wi(en)
         hr(en,en) = 1.0d0
         hi(en,en) = 0.0d0
         enm1 = en - 1
c     .......... for i=en-1 step -1 until 1 do -- ..........
         do 780 ii = 1, enm1
            i = en - ii
            zzr = 0.0d0
            zzi = 0.0d0
            ip1 = i + 1
c
            do 740 j = ip1, en
               zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
               zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
  740       continue
c
            yr = xr - wr(i)
            yi = xi - wi(i)
            if (yr .ne. 0.0d0 .or. yi .ne. 0.0d0) go to 765
               tst1 = norm
               yr = tst1
  760          yr = 0.01d0 * yr
               tst2 = norm + yr
               if (tst2 .gt. tst1) go to 760
  765       continue
            call cdiv16(zzr,zzi,yr,yi,hr(i,en),hi(i,en))
c     .......... overflow control ..........
            tr = abs(hr(i,en)) + abs(hi(i,en))
            if (tr .eq. 0.0d0) go to 780
            tst1 = tr
            tst2 = tst1 + 1.0d0/tst1
            if (tst2 .gt. tst1) go to 780
            do 770 j = i, en
               hr(j,en) = hr(j,en)/tr
               hi(j,en) = hi(j,en)/tr
  770       continue
c
  780    continue
c
  800 continue
c     .......... end backsubstitution ..........
c     .......... vectors of isolated roots ..........
      do  840 i = 1, N
         if (i .ge. low .and. i .le. igh) go to 840
c
         do 820 j = I, n
            zr(i,j) = hr(i,j)
            zi(i,j) = hi(i,j)
  820    continue
c
  840 continue
c     .......... multiply by transformation matrix to give
c                vectors of original full matrix.
c                for j=n step -1 until low do -- ..........
      do 880 jj = low, N
         j = n + low - jj
         m = min0(j,igh)
c
         do 880 i = low, igh
            zzr = 0.0d0
            zzi = 0.0d0
c
            do 860 k = low, m
               zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
               zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
  860       continue
c
            zr(i,j) = zzr
            zi(i,j) = zzi
  880 continue
c
      go to 1001
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
c
c     *****************************************************************
c
      subroutine corth16(nm,n,low,igh,ar,ai,ortr,orti)
c
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      real*16 ar(nm,n),ai(nm,n),ortr(igh),orti(igh)
      real*16 f,g,h,fi,fr,scale,pythag16
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure orthes, num. math. 12, 349-368(1968)
c     by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a complex general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     unitary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex input matrix.
c
c     on output
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the hessenberg matrix.  information
c          about the unitary transformations used in the reduction
c          is stored in the remaining triangles under the
c          hessenberg matrix.
c
c        ortr and orti contain further information about the
c          transformations.  only elements low through igh are used.
c
c     calls pythag16 for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c
      do 180 m = kp1, la
         h = 0.0d0
         ortr(m) = 0.0d0
         orti(m) = 0.0d0
         scale = 0.0d0
c     .......... scale column (algol tol then not needed) ..........
         do 90 i = m, igh
   90    scale = scale + abs(ar(i,m-1)) + abs(ai(i,m-1))
c
         if (scale .eq. 0.0d0) go to 180
         mp = m + igh
c     .......... for i=igh step -1 until m do -- ..........
         do 100 ii = m, igh
            i = mp - ii
            ortr(i) = ar(i,m-1) / scale
            orti(i) = ai(i,m-1) / scale
            h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
  100    continue
c
         g = sqrt(h)
         f = pythag16(ortr(m),orti(m))
         if (f .eq. 0.0d0) go to 103
         h = h + f * g
         g = g / f
         ortr(m) = (1.0d0 + g) * ortr(m)
         orti(m) = (1.0d0 + g) * orti(m)
         go to 105
c
  103    ortr(m) = g
         ar(m,m-1) = scale
c     .......... form (i-(u*ut)/h) * a ..........
  105    do 130 j = m, n
            fr = 0.0d0
            fi = 0.0d0
c     .......... for i=igh step -1 until m do -- ..........
            do 110 ii = m, igh
               i = mp - ii
               fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
               fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
  110       continue
c
            fr = fr / h
            fi = fi / h
c
            do 120 i = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)
               ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)
  120       continue
c
  130    continue
c     .......... form (i-(u*ut)/h)*a*(i-(u*ut)/h) ..........
         do 160 i = 1, igh
            fr = 0.0d0
            fi = 0.0d0
c     .......... for j=igh step -1 until m do -- ..........
            do 140 jj = m, igh
               j = mp - jj
               fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
               fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
  140       continue
c
            fr = fr / h
            fi = fi / h
c
            do 150 j = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)
               ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)
  150       continue
c
  160    continue
c
         ortr(m) = scale * ortr(m)
         orti(m) = scale * orti(m)
         ar(m,m-1) = -g * ar(m,m-1)
         ai(m,m-1) = -g * ai(m,m-1)
  180 continue
c
  200 return
      end
c
c     *****************************************************************
c
      subroutine csroot16(xr,xi,yr,yi)
      real*16 xr,xi,yr,yi
c
c     (yr,yi) = complex dsqrt(xr,xi) 
c     branch chosen so that yr .ge. 0.0 and sign(yi) .eq. sign(xi)
c
      real*16 s,tr,ti,pythag16
      tr = xr
      ti = xi
      s = sqrt(0.5d0*(pythag16(tr,ti) + abs(tr)))
      if (tr .ge. 0.0d0) yr = s
      if (ti .lt. 0.0d0) s = -s
      if (tr .le. 0.0d0) yi = s
      if (tr .lt. 0.0d0) yr = 0.5d0*(ti/yi)
      if (tr .gt. 0.0d0) yi = 0.5d0*(ti/yr)
      return
      end
c
c     *****************************************************************
c
      real*16 function pythag16(a,b)
      real*16 a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      real*16 p,r,s,t,u
      p = max(abs(a),abs(b))
      if (p .eq. 0.0d0) go to 20
      r = (min(abs(a),abs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag16 = p
      return
      end
c
c     *****************************************************************
c
c
