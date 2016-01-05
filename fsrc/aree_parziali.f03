!***********************************************************************
      subroutine aree_parziali(myid,nproc,lett,i_rest,bodyforce,area1,area2,area3,area4,area5,area6,kparasta,kparaend)
!***********************************************************************
      use myarrays_metri3
      use myarrays2
      use scala3
      use tipologia
      use orl
      !
      use mpi

      implicit none

!-----------------------------------------------------------------------
! array declaration
      integer i_rest,lett
      integer i,j,k,isc,myid,ierr,nproc
      integer l,numero_celle_bloccate
      integer bodyforce,kparasta,kparaend
      integer cont1,cont2,cont3,cont4,cont5,cont6

      integer inf1,inf2,inf3,inf4,inf5,inf6  

      real area1,area2,area3,area4,area5,area6
      real parziale1,parziale2
      real parziale3,parziale4
      real parziale5,parziale6
!
!-----------------------------------------------------------------------
!    OUTFLOW
!-----------------------------------------------------------------------
! in case of partial outflow
! if no outflow index = 0
! if    outflow index = 1

!-----------------------------------------------------------------------
!     INDEXES
      inf1 = 0
      inf2 = 0
      inf3 = 0
      inf4 = 0
      inf5 = 0
      inf6 = 0
      
      !case inflow outflow
      if(lett==1 .and. infout1==1)inf1=1
      if(lett==1 .and. infout2==1)inf2=1
      if(lett==1 .and. infout3==1)inf3=1
      if(lett==1 .and. infout4==1)inf4=1
      if(lett==1 .and. infout5==1)inf5=1
      if(lett==1 .and. infout6==1)inf6=1
      
      !case nesting
      if(i_rest==3 .and. infout1/=0)inf1=1
      if(i_rest==3 .and. infout2/=0)inf2=1
      if(i_rest==3 .and. infout5/=0)inf5=1
      if(i_rest==3 .and. infout6/=0)inf6=1      
      
      
            
!-----------------------------------------------------------------------
!     initialization

      area_bagnata1 = area1
      area_bagnata2 = area2
      area_bagnata3 = area3
      area_bagnata4 = area4
      area_bagnata5 = area5
      area_bagnata6 = area6

      parziale1 = 0.
      parziale2 = 0.
      parziale3 = 0.
      parziale4 = 0.
      parziale5 = 0.
      parziale6 = 0.

      do k=1,jz
      do j=1,jy
        index_out1(j,k)=1
        index_out2(j,k)=1

        index_rho1(j,k)=1
        index_rho2(j,k)=1
      end do
      end do

      do k=1,jz
      do i=1,jx
        index_out3(i,k)=1
        index_out4(i,k)=1

        index_rho3(i,k)=1
        index_rho4(i,k)=1
      end do
      end do

      do j=1,jy
      do i=1,jx
        index_out5(i,j)=1
        index_out6(i,j)=1

        index_rho5(i,j)=1
        index_rho6(i,j)=1
      end do
      end do

      cont1=0
      cont2=0
      cont3=0
      cont4=0
      cont5=0
      cont6=0

!-----------------------------------------------------------------------
!     IBM settings for outflow condition


      if( (lett==1 .or. i_rest==3) .and. bodyforce .ge.1)then
!      stop orlansky on solid cells
       open (22,file='Celle_Bloccate_Indici.inp',status='old')
       read (22,*) numero_celle_bloccate
       do l=1,numero_celle_bloccate
          read(22,*)i,j,k
!         face 1
          if(i.eq.1)then
             index_out1(j,k)=0
             cont1=cont1+1
          end if
!         face 2
          if(i.eq.jx)then
             index_out2(j,k)=0
             cont2=cont2+1
          end if
!         face 3
          if(j.eq.1)then
             index_out3(i,k)=0
             cont3=cont3+1
          end if
!         face 4
          if(j.eq.jy)then
             index_out4(i,k)=0
             cont4=cont4+1
          end if
!         face 5
          if(k.eq.1)then
             index_out5(i,j)=0
             cont5=cont5+1
          end if
!         face 6
          if(k.eq.jz)then
             index_out6(i,j)=0
             cont6=cont6+1
          end if
       end do
       close(22)
 
!      stop orlansky on ib cells
       open (22,file='Celle_IB_indici.inp',status='old')
       read (22,*) numero_celle_bloccate
       do l=1,numero_celle_bloccate
          read(22,*)i,j,k
!         face 1
          if(i.eq.1)then
             index_out1(j,k)=1 !0
             cont1=cont1+1
          end if
!         face 2
          if(i.eq.jx)then
             index_out2(j,k)=1 !0
             cont2=cont2+1
          end if
!         face 3
          if(j.eq.1)then
             index_out3(i,k)=1 !0
             cont3=cont3+1
          end if
!         face 4
          if(j.eq.jy)then
             index_out4(i,k)=1 !0
             cont4=cont4+1
          end if
!         face 5
          if(k.eq.1)then
             index_out5(i,j)=1 !0
             cont5=cont5+1
          end if
!         face 6
          if(k.eq.jz)then
             index_out6(i,j)=1 !0
             cont6=cont6+1
          end if
       end do
       close(22)      
      end if

!-----------------------------------------------------------------------
! compute partial face surface (it will be reduced if IBM is on)

!     faces 1 and 2
      if(inf1==1)then
        parziale1=0.
        area_bagnata1=0.
        do k=kparasta,kparaend  !1,jz
        do j=1,jy
           parziale1 = parziale1+areola1(j,k)*real(index_out1(j,k))
        end do
        end do

        call MPI_ALLREDUCE(parziale1,area_bagnata1,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
      end if
!
      if(inf2==1)then
        parziale2=0.
        area_bagnata2=0.
        do k=kparasta,kparaend  !1,jz
        do j=1,jy
           parziale2 = parziale2+areola2(j,k)*real(index_out2(j,k))
        end do
        end do

        call MPI_ALLREDUCE(parziale2,area_bagnata2,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
      end if
!
!     faces 3 and 4
      if(inf3==1)then
        parziale3=0.
        area_bagnata3=0.
        do k=kparasta,kparaend  !1,jz
        do i=1,jx
           parziale3 = parziale3+areola3(i,k)*real(index_out3(i,k))
        end do
        end do

        call MPI_ALLREDUCE(parziale3,area_bagnata3,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
      end if
!
      if(inf4==1)then
        parziale4=0.
        area_bagnata4=0.
        do k=kparasta,kparaend  !1,jz
        do i=1,jx
           parziale4 = parziale4+areola4(i,k)*real(index_out4(i,k))
        end do
        end do

        call MPI_ALLREDUCE(parziale4,area_bagnata4,1,MPI_REAL_SD,MPI_SUM,MPI_COMM_WORLD,ierr)
      end if
!
!     faces 5 and 6
      if(inf5==1)then
      if(myid.eq.0)then
        parziale5=0.
        area_bagnata5=0.
        do j=1,jy
        do i=1,jx
           parziale5 = parziale5+areola5(i,j)*real(index_out5(i,j))
        end do
        end do
        area_bagnata5=parziale5
      end if
      call MPI_BCAST(area_bagnata5,1,MPI_REAL_SD,0,MPI_COMM_WORLD,ierr)
      end if

      if(inf6==1)then
      if(myid.eq.nproc-1)then
        parziale6=0.
        area_bagnata6=0.
        do j=1,jy
        do i=1,jx
           parziale6 = parziale6+areola6(i,j)*real(index_out6(i,j))
        end do
        end do
        area_bagnata6=parziale6
      end if
      call MPI_BCAST(area_bagnata6,1,MPI_REAL_SD,nproc-1,MPI_COMM_WORLD,ierr)
      end if
      
      
      return
      end
      
