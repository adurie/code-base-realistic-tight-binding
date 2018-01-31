



      implicit double precision (a-h,o-y)
      implicit complex*16 (z)
      parameter (nplx=75,nslx=1,nspinx=9,nsubx=2,nsubatx=2)
      parameter (nmatx=nslx*nsubx*nspinx,natomx=nslx*nspinx*nsubatx)
      parameter (nlayx=nplx*nslx)
      parameter (nnx=4,npmax=100) // upto nnx.th N.N's
//                                    npmax - max No of points per plane
      parameter (numatx=7)       // No of atom types
      parameter (nfoldx=100)
      complex*16 ec
      common/shift/xmgo_shift,au_shift
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
      dimension vcuu(0:nlayx),vcdd(0:nlayx),vcdu(0:nlayx),vcud(0:nlayx)
      dimension rr(3)
      common/elements/iatfe(5),numat
      character*2 atname(nlayx)
      common/parallel/myid,numprocs
      dimension imp(0:nsubx)
      dimension xn(3)


      common/ssstructat/aa1(3),aa2(3),aa3(3,nlayx),vsubat(3,nsubx,nlayx)
     $,itypeat(nsubx,nlayx),imapl(nsubx),imapr(nsubx),nsubat
      common/fold/xfold(3,nfoldx),xshift(3,nfoldx),ifold(2,nfoldx),
     $  baib(3,3),nxfold,nfold,irecipa
      dimension ba1(3),ba2(3)
      dimension disttmp(numatx,numatx,nnx)

//     ------------------------------------------------------------
      double pi    = M_PI; 
      double sq2   = sqrt(2.0);
      double sq3   = sqrt(3.0);
//     ------------------------------------------------------------
//     DATA

      int nspin=9;
      int numnn=2;

      int nins=10;
      int mlay=0;     // No of substrate layers on each side of SGF
      int numat=2;    // No of atom types: one for each element

      int nlay=nins+2*mlay+4;  // total No of layers inc 4 lead layers
      int ndiff=nins;

      double ef=0.57553;
      double wm = 1e-14;

      int iwmax = 15;
      double tfac     = 8.617e-5/13.6058;
      double temp  = 315.79*tfac;


//     =================================================================
//     ATOMIC DATA FOR LEADS
//     =================================================================
//     in-plane atomic lattice vectors for substrate
//     in this code we assume that the LH and RH leads have the same 
//     lattice vectors and No of basis vectors

      nsubat=2
      natom=nsl*nspin*nsubat

      aa1(1)=0.5d0
      aa1(2)=0.5d0
      aa1(3)=0.d0
      aa2(1)=0.5d0
      aa2(2)=-0.5d
      aa2(3)=0.d0

//     =================================================================
//     LH LEAD BASIS VECTORS
//     =================================================================

      do ilay=1,2
//       Out of plane lattice vector
        aa3(1,ilay)=0.d0
        aa3(2,ilay)=0.d0
        aa3(3,ilay)=ilay-1.d0

//       Sublattice
        vsubat(1,1,ilay)=0.0d0
        vsubat(2,1,ilay)=0.0d0
        vsubat(3,1,ilay)=0.0d0
        vsubat(1,2,ilay)=0.5d0
        vsubat(2,2,ilay)=0.0d0
        vsubat(3,2,ilay)=-0.5d0


//       Atom types
        do kk=1,nsubat
          itypeat(kk,ilay)=1
        enddo
      enddo

//     =================================================================
//     SUPERCELL STRUCTURE
//     =================================================================

      nsub=2

//     2 in-plane lattice vectors
//     CUBIC
      a1(1)=0.5d0
      a1(2)=0.5d0
      a1(3)=0.d0
      a2(1)=0.5d0
      a2(2)=-0.5d0
      a2(3)=0.d0

//    ----------  Crystal structure ----------------------
      do ilay=1,nlay
//       Out of plane lattice vector
        a3(1,ilay)=0.d0
        a3(2,ilay)=0.d0
        a3(3,ilay)=ilay-1.d0

//       Sublattice
        vsub(1,1,ilay)=0.d0
        vsub(2,1,ilay)=0.d0
        vsub(3,1,ilay)=0.d0
        vsub(1,2,ilay)=0.5d0
        vsub(2,2,ilay)=0.d0
        vsub(3,2,ilay)=-0.5d0

//       Atom types
        itype(1,ilay)=1   !!!   Co
        itype(2,ilay)=1
        if((ilay >= 3) && (ilay <= nins+2))then
          itype(1,ilay)=2   !!!   Cu
          itype(2,ilay)=2
        endif
      enddo

//     =================================================================
//     RH LEAD BASIS VECTORS
//     =================================================================

      do ilay=nlay-1,nlay
//       Out of plane lattice vector
        aa3(1,ilay)=0.d0
        aa3(2,ilay)=0.d0
        aa3(3,ilay)=ilay-1.d0

//       Sublattice
        vsubat(1,1,ilay)=0.d0
        vsubat(2,1,ilay)=0.d0
        vsubat(3,1,ilay)=0.d0
        vsubat(1,2,ilay)=0.5d0
        vsubat(2,2,ilay)=0.d0
        vsubat(3,2,ilay)=-0.5d0

//       Atom types
        do kk=1,nsubat
          itypeat(kk,ilay)=1
        enddo
      enddo
//

//     =================================================================
//     The map between Supercell sublattice and LH atomic sublattice :
//     imap:   supercell --> atomic
//     Defined s.t.     vsub(k)=vsubat(imap(k)) + atomic lattice vector

      imapl(1)=1
      imapl(2)=2

      imapr(1)=1
      imapr(2)=2

//
//     =================================================================
//     THIS SECTION TO GET NN DISTANCES ONLY

//     In and out of plane distances:
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

//     ddnn - atom-atom distances  ... this needs to be checked
//     Better to put NN distances in by hand as in next section.
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
//     =================================================================
//     NN DISTANCES
//     1,2 Fe  ; 3 Mg ; 4 O 
//     ddnn(type1,type2,NN)

      do i=1,numat
        do j=1,numat
          do k=1,numnn
            ddnn(i,j,k)=0.d0
          enddo 
        enddo
      enddo

      do ii=1,numat
        do jj=1,numat
          ddnn(ii,jj,1)=1.d0/sqrt(2.d0)
          ddnn(ii,jj,2)=1.d0
        enddo
      enddo


//     =================================================================
//     !!!!!!! OUTPUT ATOMIC POSITIONS FOR RASMOL VIEWER !!!!!!
//     load into rasmol with command :    > load xyz pos0.dat
      atname(1)="Co"
      atname(2)="Cu"

//     whole cluster
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

//     =================================================================
//     Dimensionality Checks

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

//     construct the growth direction:
      call cross(a1,a2,xn,3)
      xnorm=dot(xn,xn,3)
      do i=1,3
        xn(i)=xn(i)/sqrt(xnorm)
      enddo
//     =================================================================

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
        write(*,*)' nq =',nq
        write(*,*)' ef =',ef
        write(*,*)
        write(*,*)' No of MgO at. layers =',nins
        write(*,*)

//       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//       this prints out all interplanar distances
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

//       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        write(*,'(65("="))')
        write(*,*)"      CALCULATED NN DISTANCES --- check these"
        write(*,'(65("="))')
        write(*,'("  at1   at2  ",100("    dNN",i1,"       "))')
     $(k,k=1,numnn-1)
        do i=1,numat
          do j=1,numat
            write(*,'(2i5,100d15.6)')i,j,(disttmp(i,j,k),k=1,numnn)
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


C     ec = dcmplx(ef,wm)
C     if(myid.eq.0)write(*,*)' Complex Energy = ',ec
C     do in=0,ndiff
C       zresu(in)=0.d0
C       zresd(in)=0.d0
C       zresud(in)=0.d0
C       zresdu(in)=0.d0
C     enddo
C     call sumk(irecip,ec,nq,b1,b2,zresu,zresd,zresud,zresdu)



      do in=0,ndiff
        vcuu(in) = 0.d0
        vcud(in) = 0.d0
        vcdu(in) = 0.d0
        vcdd(in) = 0.d0
        zresu(in) = 0.d0
        zresd(in) = 0.d0
        zresud(in) = 0.d0
        zresdu(in) = 0.d0
      enddo

      fact = (pi + pi)*temp
      do iw=1,iwmax
        wm = pi*dfloat(2*iw-1)*temp
        ec = dcmplx(ef,wm)
        call kcon(irecip,ec,fact,b1,b2,zresu,zresd,zresud,zresdu)
//     call sumk(irecip,ec,nq,b1,b2,zresu,zresd,zresud,zresdu)
        do in=0,ndiff
          vcuu(in) = vcuu(in) - fact*dreal(zresu(in))
          vcud(in) = vcud(in) - fact*dreal(zresud(in))
          vcdu(in) = vcdu(in) - fact*dreal(zresdu(in))
          vcdd(in) = vcdd(in) - fact*dreal(zresd(in))
        enddo
        if(myid.eq.0)write(*,*)' iw =',iw,'   complete'
        if(myid.eq.0)write(*,*)
        call flush(6)
      enddo

      if (myid == 0){
	      cout<<endl;
	      cout<<"****************************************************"<<endl;
	      cout<<"****************************************************"<<endl;
	      cout<<"UP SPIN"<<endl;
	      for (int in=0; in <=ndiff; in++){
	        cout<<2*in<<" "<<vcuu(in)/nsub<<endl;
	      }
	      cout<<endl;
	      cout<<"****************************************************"<<endl;
	      cout<<"****************************************************"<<endl;
	      cout<<"DOWN SPIN"<<endl;
	      for (int in=0; in <=ndiff; in++){
	        cout<<2*in<<" "<<vcdd(in)/nsub<<endl;
	      }
	      cout<<endl;
	      cout<<"****************************************************"<<endl;
	      cout<<"****************************************************"<<endl;
	      cout<<"UP-DOWN SPIN"<<endl;
	      for (int in=0; in <=ndiff; in++){
	        cout<<2*in<<" "<<vcud(in)/nsub<<endl;
	      }
	      cout<<endl;
	      cout<<"****************************************************"<<endl;
	      cout<<"****************************************************"<<endl;
	      cout<<"DOWN-UP SPIN"<<endl;
	      for (int in=0; in <=ndiff; in++){
	        cout<<2*in<<" "<<vcdu(in)/nsub<<endl;
	      }
	      cout<<endl;
	      cout<<"****************************************************"<<endl;
	      cout<<"****************************************************"<<endl;
	      cout<<"(UP + DOWN - 2*AF)"<<endl;
	      for (int in=0; in <=ndiff; in++){
	        cout<<2*in<<" "<<-(vcuu(in)+vcdd(in)-vcdu(in)-vcud(in))/nsub<<endl;
	      }
      }
}
