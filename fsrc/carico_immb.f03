!***********************************************************************
subroutine carico_immb(tipo)
   !***********************************************************************
   !     read ibm input file
   use myarrays_ibm
   use mysending
   !
   use scala3
   use period
   !
   use mpi

   implicit none


   !-----------------------------------------------------------------------
   !     array declaration
   integer l,i,j,k,in,jn,kn
   integer iloop,jloop,kloop
      
   integer ireq1,ireq2,ireq3,ireq4
   integer status(MPI_STATUS_SIZE),ierror,ierr

   integer itipo
   integer contatore,num_solide_real
   integer ib_totali,solide_totali

   real a1,a2,a3
   real dist_x,dist_y,dist_z,dist1,dist2
      
   real rot11,rot12,rot13
   real rot21,rot22,rot23
   real rot31,rot32,rot33
       
   real irot11,irot12,irot13
   real irot21,irot22,irot23
   real irot31,irot32,irot33
      
   real tri1,tri2,tri3,tri4
   integer i_tri1,i_tri2,i_tri3,i_tri4
   integer j_tri1,j_tri2,j_tri3,j_tri4
   integer k_tri1,k_tri2,k_tri3,k_tri4
      
   !      integer tipo(0:n1+1,0:n2+1,0:n3+1)
   integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

   logical icheck
   integer np
   !-----------------------------------------------------------------------
   icheck = .false.
   allocate(tipo_spedito(0:n1+1,0:n2+1,kparasta-1:kparaend+1))
   !-----------------------------------------------------------------------
   if(myid.eq.0)then
      WRITE(*,*)' '
      write(*,*)'*****************************************'
      write(*,*)'     LOAD IBM'
   end if
   !-----------------------------------------------------------------------
   !     initialization
   do k=kparasta-deepl,kparaend+deepr !0,jz+1
      do j=0,jy+1
         do i=0,jx+1
            tipo(i,j,k)=2
         end do
      end do
   end do

   do i=1,MN
      do j=1,6
         indici_CELLE_IB(i,j)=0
      end do
   end do

   do i=1,MN
      do j=1,3
         distanze_CELLE_IB(i,j)=0
      end do
   end do

   do i=1,MP
      do j=1,3
         indici_celle_bloccate(i,j)=0
      end do
   end do

   if(myid.eq.0)then
      write(*,*)'matrix ---> OK'
   end if

   !-----------------------------------------------------------------------
   !     read input files for ibm


   !     ib point
   open (20,file='Celle_IB_indici.inp',status='old')
   open (21,file='Celle_IB_distanze.inp',status='old')
   open (23,file='distanze_interpolazioni.inp',status='old')
   open (24,file='rotazione.inp',status='old')
      
   open (25,file='trilinear_ibm.inp',status='old')
   open (26,file='trilinear_i.inp',status='old')
   open (27,file='trilinear_j.inp',status='old')
   open (28,file='trilinear_k.inp',status='old')

   read (20,*) numero_celle_IB

   num_ib=0

   do l=1,numero_celle_IB

      read(20,*)i,j,k,in,jn,kn

      read(21,310)dist_x,dist_y,dist_z


      read(23,320)dist1,dist2,a1,a2,a3
	
      read(24,350)rot11,rot12,rot13, &
         rot21,rot22,rot23, &
         rot31,rot32,rot33

      read(24,350)irot11,irot12,irot13, &
         irot21,irot22,irot23, &
         irot31,irot32,irot33
     
      read(25,360)tri1,tri2,tri3,tri4

      read(26,370)i_tri1,i_tri2,i_tri3,i_tri4
     
      read(27,370)j_tri1,j_tri2,j_tri3,j_tri4
     
      read(28,370)k_tri1,k_tri2,k_tri3,k_tri4

      if((k.ge.kparasta-deepl).and.(k.le.kparaend+deepr))tipo(i,j,k)=1
	
      if((k.ge.kparasta).and.(k.le.kparaend))then

         tipo(i,j,k)=1

         num_ib=num_ib+1

         indici_CELLE_IB(num_ib,1)=i
         indici_CELLE_IB(num_ib,2)=j
         indici_CELLE_IB(num_ib,3)=k
         indici_CELLE_IB(num_ib,4)=in
         indici_CELLE_IB(num_ib,5)=jn
         indici_CELLE_IB(num_ib,6)=kn

         distanze_CELLE_IB(num_ib,1)=dist_x
         distanze_CELLE_IB(num_ib,2)=dist_y
         distanze_CELLE_IB(num_ib,3)=dist_z

         dist_pp_ib(num_ib)=dist1
         dist_ib_parete(num_ib)=dist2

         proiezioni(num_ib,1)=a1
         proiezioni(num_ib,2)=a2
         proiezioni(num_ib,3)=a3
	  
	  
         !         rotation matrix (to construct tangential and normal velocity)
         rot(num_ib,1,1)=rot11
         rot(num_ib,1,2)=rot12
         rot(num_ib,1,3)=rot13

         rot(num_ib,2,1)=rot21
         rot(num_ib,2,2)=rot22
         rot(num_ib,2,3)=rot23
	  
         rot(num_ib,3,1)=rot31
         rot(num_ib,3,2)=rot32
         rot(num_ib,3,3)=rot33


         !         inverse rotation matrix
         rot_inverse(num_ib,1,1)=irot11
         rot_inverse(num_ib,1,2)=irot12
         rot_inverse(num_ib,1,3)=irot13

         rot_inverse(num_ib,2,1)=irot21
         rot_inverse(num_ib,2,2)=irot22
         rot_inverse(num_ib,2,3)=irot23
	  
         rot_inverse(num_ib,3,1)=irot31
         rot_inverse(num_ib,3,2)=irot32
         rot_inverse(num_ib,3,3)=irot33
	  
         !         trilinear coefficent
         tricoef(num_ib,1) = tri1
         tricoef(num_ib,2) = tri2
         tricoef(num_ib,3) = tri3
         tricoef(num_ib,4) = tri4

	  
         !         trilinear index
         trind(num_ib,1,1) = i_tri1
         trind(num_ib,2,1) = i_tri2
         trind(num_ib,3,1) = i_tri3
         trind(num_ib,4,1) = i_tri4
	   

         trind(num_ib,1,2) = j_tri1
         trind(num_ib,2,2) = j_tri2
         trind(num_ib,3,2) = j_tri3
         trind(num_ib,4,2) = j_tri4


         trind(num_ib,1,3) = k_tri1
         trind(num_ib,2,3) = k_tri2
         trind(num_ib,3,3) = k_tri3
         trind(num_ib,4,3) = k_tri4
				   
				  
      end if
   end do

   close(20)
   close(21)
   close(23)
   close(24)

   close(25)
   close(26)
   close(27)
   close(28)
      
      

300 format(6(i4,1x))
310 format(3e15.8)
320 format (5e15.8)
350 format(9e15.8)
360 format(4e15.8)     ! tri_ibm
370 format(4(i8,1x))

   !      write(2000+myid,*)myid,'number of ib points: ',num_ib
   !.......................................................................


   !     solid points
   open (22,file='Celle_Bloccate_Indici.inp',status='old')

   read (22,*) numero_celle_bloccate

   num_solide=0
   contatore=0
   do l=1,numero_celle_bloccate

      read(22,*)i,j,k

      if((k.ge.kparasta-deepl).and.(k.le.kparaend+deepr))then

         tipo(i,j,k)=0

         num_solide=num_solide+1

         indici_celle_bloccate(num_solide,1)=i
         indici_celle_bloccate(num_solide,2)=j
         indici_celle_bloccate(num_solide,3)=k

         if(k.gt.kparaend .or. k.lt.kparasta)then
            contatore=contatore+1
         end if
      end if
   end do

   close (22)
   !      if(myid.eq.0)then
   !      write(*,*)'read solid cells ---> OK'
   !      end if
   write(200,*)'border solid',contatore

   num_solide_real=num_solide-contatore !without border cells

   write(200,*)myid,'number of solid cells: ',num_solide
      
   !-----------------------------------------------------------------------
   !     check number of ibm for each proc

   call MPI_REDUCE(num_ib,ib_totali,1,MPI_INTEGER,MPI_SUM, &
      0,MPI_COMM_WORLD,ierr)

   call MPI_REDUCE(num_solide_real,solide_totali,1,MPI_INTEGER, &
      MPI_SUM,0,MPI_COMM_WORLD,ierr)

   if(myid.eq.0)then
      write(*,*)'--------------------------------------------'
      if(ib_totali.eq.numero_celle_IB)then
         write(*,*)'check IB for each proc --> OK'
      else
         write(*,*)'check IB for each proc --> NO'
      end if

      if(solide_totali.eq.numero_celle_bloccate)then
         write(*,*)'check solid cells for each proc --> OK'
      else
         write(*,*)'check solid cells for each proc --> NO'
      end if

      write(*,*)'--------------------------------------------'
      write(*,*)'read input file for Immersed Boundaries'
      write(*,*)'total number of IB',numero_celle_IB
      write(*,*)'total number of solid',numero_celle_bloccate
      write(*,*)'--------------------------------------------'
   end if

   write(200,*)'number IB proc --->',num_ib
   write(200,*)'number solid cells proc --->',num_solide
   write(200,*)'--------------------------------------------'

   if(myid.eq.0)then
      write(*,*)'       LOAD IBM finished'
      write(*,*)'*****************************************'
      write (*,*)' '
   end if

   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------
   !     PREPARE COMUNICATION BETWEEN PROCS FOR IB
   allocate( stencil_left_snd(2*jx*jy,3))
   allocate(stencil_right_snd(2*jx*jy,3))
   allocate( stencil_left_rcv(2*jx*jy,3))
   allocate(stencil_right_rcv(2*jx*jy,3))

   stencil_left_rcv  = 0
   stencil_right_rcv = 0
   stencil_left_snd  = 0
   stencil_right_snd = 0

   ! my stencil requires node of the close proc
   ! therefore kparasta-1 and kparaend+1
       
   num_left_snd = 0
   num_right_snd = 0
   num_right_rcv = 0
   num_left_rcv = 0
      
   tipo_spedito = 0
      
   do l=1,num_ib
      do np = 1,4
         i=trind(l,np,1)
         j=trind(l,np,2)
         k=trind(l,np,3)


         if(tipo_spedito(i,j,k)==0)then
	 
            if(k==kparasta-1)then
               num_left_rcv = num_left_rcv + 1
               stencil_left_rcv(num_left_rcv,1) = i
               stencil_left_rcv(num_left_rcv,2) = j
               stencil_left_rcv(num_left_rcv,3) = k
	
               tipo_spedito(i,j,k) = 1
            end if


            if(k==kparaend+1)then
               num_right_rcv = num_right_rcv + 1
               stencil_right_rcv(num_right_rcv,1) = i
               stencil_right_rcv(num_right_rcv,2) = j
               stencil_right_rcv(num_right_rcv,3) = k

               tipo_spedito(i,j,k) = 1
            end if

         end if ! tipo_spedito
	 
      end do
   end do


   !      write(*,*)myid,'num left recv',num_left_rcv
   !      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
   !      write(*,*)myid,'num right recv',num_right_rcv
   !      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   !     COMUNICATION OF INDEX AND COUNTER
   !
   !     what I send in stencil_left_send will be recived in stencil_right_recive
   !     in closer left proc
   !
   !     and the same for the right
   !
   !     what I send in stencil_right_send will be recived in stencil_left_recive
   !     in closer right proc
   !
   !     step 1, comunication of counter, how much I send/recive









   !     comunicate left the point left has to send right
   if(myid.ne.0)then
      call MPI_SSEND(num_left_rcv,1,MPI_INTEGER, &
         leftpe,tagls,MPI_COMM_WORLD,ierror)
   !      call MPI_WAIT (ireq1,status,ierror)
   end if
   if(myid.ne.nproc-1)then
      call MPI_RECV(num_right_snd,1,MPI_INTEGER, &
         rightpe,tagrr,MPI_COMM_WORLD,status,ierror)
   !      call MPI_WAIT (ireq2,status,ierror)
   end if

   !     comunicate right the point right has to send left
   if(myid.ne.nproc-1)then
      call MPI_SSEND(num_right_rcv,1,MPI_INTEGER, &
         rightpe,tagrs,MPI_COMM_WORLD,ierror)
   !      call MPI_WAIT (ireq3,status,ierror)
   end if
   if(myid.ne.0)then
      call MPI_RECV(num_left_snd,1,MPI_INTEGER, &
         leftpe,taglr,MPI_COMM_WORLD,status,ierror)
   !      call MPI_WAIT (ireq4,status,ierror)
   end if

   !     if periodic
   if(kp==0)then
      if(myid.eq.0)then
         call MPI_SSEND(num_left_rcv,1,MPI_INTEGER, &
            nproc-1,tagls,MPI_COMM_WORLD,ierror)
      !          call MPI_WAIT (ireq1,status,ierror)
      end if
      if(myid.eq.nproc-1)then
         call MPI_RECV(num_right_snd,1,MPI_INTEGER, &
            0,tagrr,MPI_COMM_WORLD,status,ierror)
      !          call MPI_WAIT (ireq2,status,ierror)
      end if

      if(myid.eq.nproc-1)then
         call MPI_SSEND(num_right_rcv,1,MPI_INTEGER, &
            0,tagrs,MPI_COMM_WORLD,ierror)
      !          call MPI_WAIT (ireq3,status,ierror)
      end if
      if(myid.eq.0)then
         call MPI_RECV(num_left_snd,1,MPI_INTEGER, &
            nproc-1,taglr,MPI_COMM_WORLD,status,ierror)
      !          call MPI_WAIT (ireq4,status,ierror)
      end if
   else
      if(myid.eq.0)then
         num_left_snd = 0
         num_left_rcv = 0
      end if
      if(myid.eq.nproc-1)then
         num_right_rcv = 0
         num_right_snd = 0
      end if
   end if

   !     I know the counter! now I send the values

   !     send left the indeces left has to send right
   if(myid.ne.0)then
      call MPI_SSEND(stencil_left_rcv(1,1),3*2*jx*jy,MPI_INTEGER, &
         leftpe,tagls,MPI_COMM_WORLD,ierror)
   !      call MPI_WAIT (ireq1,status,ierror)
   end if
   if(myid.ne.nproc-1)then
      call MPI_RECV(stencil_right_snd(1,1),3*2*jx*jy,MPI_INTEGER, &
         rightpe,tagrr,MPI_COMM_WORLD,status,ierror)
   !      call MPI_WAIT (ireq2,status,ierror)
   end if

   !     send right the indeces right has to send left
   if(myid.ne.nproc-1)then
      call MPI_SSEND(stencil_right_rcv(1,1),3*2*jx*jy,MPI_INTEGER, &
         rightpe,tagrs,MPI_COMM_WORLD,ierror)
   !      call MPI_WAIT (ireq3,status,ierror)
   end if
   if(myid.ne.0)then
      call MPI_RECV(stencil_left_snd(1,1),3*2*jx*jy,MPI_INTEGER, &
         leftpe,taglr,MPI_COMM_WORLD,status,ierror)
   !      call MPI_WAIT (ireq4,status,ierror)
   end if

   !     if periodic
   if(kp==0)then
      
      !       left to 0 is nproc-1
      if(myid.eq.0)then
         call MPI_SSEND(stencil_left_rcv(1,1),3*2*jx*jy,MPI_INTEGER,   &
            nproc-1,tagls,MPI_COMM_WORLD,ierror)
      !          call MPI_WAIT (ireq1,status,ierror)
      end if
      if(myid.eq.nproc-1)then
         call MPI_RECV(stencil_right_snd(1,1),3*2*jx*jy,MPI_INTEGER, &
            0,tagrr,MPI_COMM_WORLD,status,ierror)
      !          call MPI_WAIT (ireq2,status,ierror)
      end if
	
      !	right to nproc-1 is 0
      if(myid.eq.nproc-1)then
         call MPI_SSEND(stencil_right_rcv(1,1),3*2*jx*jy,MPI_INTEGER, &
            0,tagrs,MPI_COMM_WORLD,ierror)
      !          call MPI_WAIT (ireq3,status,ierror)
      end if
      if(myid.eq.0)then
         call MPI_RECV(stencil_left_snd(1,1),3*2*jx*jy,MPI_INTEGER, &
            nproc-1,taglr,MPI_COMM_WORLD,status,ierror)
      !          call MPI_WAIT (ireq4,status,ierror)
      end if
	
      !       now myid=0 and myid=nproc-1 know the index they will recive,
      !       but they need to save
      !       on the border, so I change the k

      if(myid == nproc-1)stencil_right_snd(:,3)=jz
      if(myid == 0)      stencil_left_snd(:,3)=1
	   
   end if
      
   !-----------------------------------------------------------------------
   !-----------------------------------------------------------------------

      

   !--------------------------------------------------------------------------
   if(icheck)then
      !     check
      write(*,*)myid,'send left', num_left_snd
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
             
      write(*,*)myid,'send right', num_right_snd 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)      
      
      write(*,*)myid,'recive left', num_left_rcv 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            
      write(*,*)myid,'recive right', num_right_rcv             
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)     
      !      write(*,*)myid,'send recv solid',
      !     > numsolid_left_snd  ,
      !     > numsolid_right_snd ,
      !     > numsolid_right_rcv ,
      !     > numsolid_left_rcv
      

             
      !      do l=1,numsolid_left_snd
      !	 write(4100+myid,*)solid_left_snd(l,1),
      !     >  	           solid_left_snd(l,2),
      !     >  	           solid_left_snd(l,3)
      !      end do

      !      do l=1,numsolid_right_snd
      !	 write(4200+myid,*)solid_right_snd(l,1),
      !     >  	           solid_right_snd(l,2),
      !     >  	           solid_right_snd(l,3)
      !      end do

      !      do l=1,numsolid_left_rcv
      !	 write(4300+myid,*)solid_left_rcv(l,1),
      !     >  	           solid_left_rcv(l,2),
      !     >  	           solid_left_rcv(l,3)
      !      end do

      !      do l=1,numsolid_right_rcv
      !	 write(4400+myid,*)solid_right_rcv(l,1),
      !     >  	           solid_right_rcv(l,2),
      !     >  	           solid_right_rcv(l,3)
      !      end do



      do l=1,num_left_snd
         write(5100+myid,*)stencil_left_snd(l,1), &
            stencil_left_snd(l,2), &
            stencil_left_snd(l,3)
      end do 

      do l=1,num_right_snd
         write(5200+myid,*)stencil_right_snd(l,1), &
            stencil_right_snd(l,2), &
            stencil_right_snd(l,3)
      end do 

      do l=1,num_left_rcv
         write(5300+myid,*)stencil_left_rcv(l,1), &
            stencil_left_rcv(l,2), &
            stencil_left_rcv(l,3)
      end do 

      do l=1,num_right_rcv
         write(5400+myid,*)stencil_right_rcv(l,1), &
            stencil_right_rcv(l,2), &
            stencil_right_rcv(l,3)
      end do 

   end if

   return
end
