
      dimension bba(3,3),bb(3,3),baib(3,3)
      do k=1,3
        bb(k,1)=1
        bb(k,2)=2
        bb(k,3)=3
        bba(k,1)=4
        bba(k,2)=5
        bba(k,3)=6
      enddo
      do k = 1,3
        write(*,*)(bb(k,i),i=1,3)
        write(*,*)
      enddo
      write(*,*)
      do k = 1,3
        write(*,*)(bba(k,i),i=1,3)
        write(*,*)
      enddo
      call remultiply(bba,bb,baib,3,3)
      write(*,*)
      do k = 1,3
        write(*,*)(baib(k,i),i=1,3)
        write(*,*)
      enddo
      write(*,*)0.50000001d0
      end

      subroutine remultiply(x,y,xy,nmat,nmatx)
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
