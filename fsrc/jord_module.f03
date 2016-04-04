module jord_module

   use myarrays_metri3
   use myarrays_velo3
   !
   use scala3
   use tipologia
   use velpar
   !
   use mpi

   implicit none

   private

   public :: jord1,jord2,jord3

contains
   !***********************************************************************
   subroutine jord1(in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1, &
      kparasta,kparaend,rightpe,leftpe, &
      tagls,taglr,tagrs,tagrr,rightpem,leftpem)
      !***********************************************************************
      ! compute the gradient dudx, dvdx, dwdx
      !

      !-----------------------------------------------------------------------
      !     array declaration
      integer ierr,myid,nproc,status(MPI_STATUS_SIZE)
      integer ncolperproc,kparasta,kparaend,m
      !
      integer tagls,taglr,tagrs,tagrr
      integer leftpe,rightpe
      integer leftpem,rightpem
      integer req1,req2,req3,req4,req5,req6
      integer istatus(MPI_STATUS_SIZE)
      !
      integer i,j,k
      integer in_dx1(n1,n2,kparasta:kparaend) ! n3)
      integer in_sn1(n1,n2,kparasta:kparaend) !n3)
      integer in_sp1(n1,n2,kparasta:kparaend) !n3)
      integer in_st1(n1,n2,kparasta:kparaend) !n3)
      integer in_av1(n1,n2,kparasta:kparaend) !n3)
      integer in_in1(n1,n2,kparasta:kparaend) !n3)
      !
      real an
      real u1,u2,u3,u4,u5,u6
      real v1,v2,v3,v4,v5,v6
      real w1,w2,w3,w4,w5,w6
      !-----------------------------------------------------------------------
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      !-----------------------------------------------------------------------
      an=1./re
      !
      do k=kparasta,kparaend
         do j=1,jy
            do i=1,jx
               !
               u2=in_dx1(i,j,k)*u(i+1,j,k)+(1-in_dx1(i,j,k))*udx(j,k)
               u1=in_sn1(i,j,k)*u(i-1,j,k)+(1-in_sn1(i,j,k))*usn(j,k)
               !
               u4=in_sp1(i,j,k)*u(i,j+1,k)+(1-in_sp1(i,j,k))*usp(i,k)
               u3=in_st1(i,j,k)*u(i,j-1,k)+(1-in_st1(i,j,k))*ust(i,k)
               !
               u6=in_av1(i,j,k)*u(i,j,k+1)+(1-in_av1(i,j,k))*uav(i,j)
               u5=in_in1(i,j,k)*u(i,j,k-1)+(1-in_in1(i,j,k))*uin(i,j)
               !
               v2=in_dx1(i,j,k)*v(i+1,j,k)+(1-in_dx1(i,j,k))*vdx(j,k)
               v1=in_sn1(i,j,k)*v(i-1,j,k)+(1-in_sn1(i,j,k))*vsn(j,k)
               !
               v4=in_sp1(i,j,k)*v(i,j+1,k)+(1-in_sp1(i,j,k))*vsp(i,k)
               v3=in_st1(i,j,k)*v(i,j-1,k)+(1-in_st1(i,j,k))*vst(i,k)
               !
               v6=in_av1(i,j,k)*v(i,j,k+1)+(1-in_av1(i,j,k))*vav(i,j)
               v5=in_in1(i,j,k)*v(i,j,k-1)+(1-in_in1(i,j,k))*vin(i,j)
               !
               w2=in_dx1(i,j,k)*w(i+1,j,k)+(1-in_dx1(i,j,k))*wdx(j,k)
               w1=in_sn1(i,j,k)*w(i-1,j,k)+(1-in_sn1(i,j,k))*wsn(j,k)
               !
               w4=in_sp1(i,j,k)*w(i,j+1,k)+(1-in_sp1(i,j,k))*wsp(i,k)
               w3=in_st1(i,j,k)*w(i,j-1,k)+(1-in_st1(i,j,k))*wst(i,k)
               !
               w6=in_av1(i,j,k)*w(i,j,k+1)+(1-in_av1(i,j,k))*wav(i,j)
               w5=in_in1(i,j,k)*w(i,j,k-1)+(1-in_in1(i,j,k))*win(i,j)


               !
               ! compute nu_T*dudx
               !
               gra1(i,j,k)=(csx(i,j,k)*.5*(u(i,j,k)+u2)- &
                  csx(i-1,j,k)*.5*(u(i,j,k)+u1)+ &
                  etx(i,j,k)*.5*(u(i,j,k)+u4)- &
                  etx(i,j-1,k)*.5*(u(i,j,k)+u3)+ &
                  ztx(i,j,k)*.5*(u(i,j,k)+u6)- &
                  ztx(i,j,k-1)*.5*(u(i,j,k)+u5)) &
                  *(annit(i,j,k)-an)/giac(i,j,k)
               !
               ! compute nu_T*dvdx
               !
               gra2(i,j,k)=(csx(i,j,k)*.5*(v(i,j,k)+v2)- &
                  csx(i-1,j,k)*.5*(v(i,j,k)+v1)+ &
                  etx(i,j,k)*.5*(v(i,j,k)+v4)- &
                  etx(i,j-1,k)*.5*(v(i,j,k)+v3)+ &
                  ztx(i,j,k)*.5*(v(i,j,k)+v6)- &
                  ztx(i,j,k-1)*.5*(v(i,j,k)+v5)) &
                  *(annit(i,j,k)-an)/giac(i,j,k)
               !
               ! compute nu_T*dwdx
               !
               gra3(i,j,k)=(csx(i,j,k)*.5*(w(i,j,k)+w2)- &
                  csx(i-1,j,k)*.5*(w(i,j,k)+w1)+ &
                  etx(i,j,k)*.5*(w(i,j,k)+w4)- &
                  etx(i,j-1,k)*.5*(w(i,j,k)+w3)+ &
                  ztx(i,j,k)*.5*(w(i,j,k)+w6)- &
                  ztx(i,j,k-1)*.5*(w(i,j,k)+w5)) &
                  *(annit(i,j,k)-an)/giac(i,j,k)
            !

            enddo
         enddo
      enddo
      !
      ! make gra1, gra2 and gra3 at the border known between the procs
      !
      if(leftpem /= MPI_PROC_NULL) then
         call MPI_SSEND(gra1(1,1,kparasta),jx*jy,MPI_REAL_SD,leftpem ,tagls,MPI_COMM_WORLD,ierr)
         call MPI_SSEND(gra2(1,1,kparasta),jx*jy,MPI_REAL_SD,leftpem ,tagls,MPI_COMM_WORLD,ierr)
         call MPI_SSEND(gra3(1,1,kparasta),jx*jy,MPI_REAL_SD, leftpem ,tagls,MPI_COMM_WORLD,ierr)
      endif
      if(rightpem /= MPI_PROC_NULL) then
         call MPI_RECV(gra1(1,1,kparaend+1),jx*jy,MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
         call MPI_RECV(gra2(1,1,kparaend+1),jx*jy,MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
         call MPI_RECV(gra3(1,1,kparaend+1),jx*jy,MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
      endif

      if(leftpem /= MPI_PROC_NULL) then
      !      call MPI_WAIT(req1,istatus,ierr)
      !      call MPI_WAIT(req3,istatus,ierr)
      !      call MPI_WAIT(req5,istatus,ierr)
      endif
      if(rightpem /= MPI_PROC_NULL) then
      !      call MPI_WAIT(req2,istatus,ierr)
      !      call MPI_WAIT(req4,istatus,ierr)
      !      call MPI_WAIT(req6,istatus,ierr)
      endif

      ! send k=jz to P0 and k=1 to Pn-1 to compute at the border
      ! cgra3 in flu_turbo

      if (myid.eq.nproc-1) then

         call MPI_SENDRECV(gra1(1,1,jz),jx*jy,MPI_REAL_SD,0,101,gra1_appoggio(1,1,1),jx*jy, &
         MPI_REAL_SD,0,201,MPI_COMM_WORLD,status,ierr)

         call MPI_SENDRECV(gra2(1,1,jz),jx*jy,MPI_REAL_SD,0,102,gra2_appoggio(1,1,1),jx*jy, &
         MPI_REAL_SD,0,202,MPI_COMM_WORLD,status,ierr)

         call MPI_SENDRECV(gra3(1,1,jz),jx*jy,MPI_REAL_SD,0,103,gra3_appoggio(1,1,1),jx*jy, &
         MPI_REAL_SD,0,203,MPI_COMM_WORLD,status,ierr)

      endif

      if (myid.eq.0) then

         call MPI_SENDRECV(gra1(1,1,1),jx*jy,MPI_REAL_SD,nproc-1,201,gra1_appoggio(1,1,jz),jx*jy, &
         MPI_REAL_SD,nproc-1,101,MPI_COMM_WORLD,status,ierr)

         call MPI_SENDRECV(gra2(1,1,1),jx*jy,MPI_REAL_SD,nproc-1,202,gra2_appoggio(1,1,jz),jx*jy, &
         MPI_REAL_SD,nproc-1,102,MPI_COMM_WORLD,status,ierr)

         call MPI_SENDRECV(gra3(1,1,1),jx*jy,MPI_REAL_SD,nproc-1,203,gra3_appoggio(1,1,jz),jx*jy, &
         MPI_REAL_SD,nproc-1,103,MPI_COMM_WORLD,status,ierr)

      endif

      return
   end

   !***********************************************************************
   subroutine jord2(in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1, &
      kparasta,kparaend,rightpe,leftpe, &
      tagls,taglr,tagrs,tagrr,rightpem,leftpem)
      !***********************************************************************
      ! compute the gradient dudy, dvdy, dwdy
      !

      !-----------------------------------------------------------------------
      !     array declaration
      integer ierr,myid,nproc,status(MPI_STATUS_SIZE)
      integer ncolperproc,kparasta,kparaend,m
      !
      integer tagls,taglr,tagrs,tagrr
      integer leftpe,rightpe
      integer leftpem,rightpem
      integer req1,req2,req3,req4,req5,req6
      integer istatus(MPI_STATUS_SIZE)
      !
      integer i,j,k
      integer in_dx1(n1,n2,kparasta:kparaend) ! n3)
      integer in_sn1(n1,n2,kparasta:kparaend) !n3)
      integer in_sp1(n1,n2,kparasta:kparaend) !n3)
      integer in_st1(n1,n2,kparasta:kparaend) !n3)
      integer in_av1(n1,n2,kparasta:kparaend) !n3)
      integer in_in1(n1,n2,kparasta:kparaend) !n3)
      !
      real an
      real u1,u2,u3,u4,u5,u6
      real v1,v2,v3,v4,v5,v6
      real w1,w2,w3,w4,w5,w6
      !-----------------------------------------------------------------------
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      !-----------------------------------------------------------------------
      an=1./re
      !
      do k=kparasta,kparaend
         do j=1,jy
            do i=1,jx
               !
               u2=in_dx1(i,j,k)*u(i+1,j,k)+(1-in_dx1(i,j,k))*udx(j,k)
               u1=in_sn1(i,j,k)*u(i-1,j,k)+(1-in_sn1(i,j,k))*usn(j,k)
               !
               u4=in_sp1(i,j,k)*u(i,j+1,k)+(1-in_sp1(i,j,k))*usp(i,k)
               u3=in_st1(i,j,k)*u(i,j-1,k)+(1-in_st1(i,j,k))*ust(i,k)
               !
               u6=in_av1(i,j,k)*u(i,j,k+1)+(1-in_av1(i,j,k))*uav(i,j)
               u5=in_in1(i,j,k)*u(i,j,k-1)+(1-in_in1(i,j,k))*uin(i,j)
               !
               v2=in_dx1(i,j,k)*v(i+1,j,k)+(1-in_dx1(i,j,k))*vdx(j,k)
               v1=in_sn1(i,j,k)*v(i-1,j,k)+(1-in_sn1(i,j,k))*vsn(j,k)
               !
               v4=in_sp1(i,j,k)*v(i,j+1,k)+(1-in_sp1(i,j,k))*vsp(i,k)
               v3=in_st1(i,j,k)*v(i,j-1,k)+(1-in_st1(i,j,k))*vst(i,k)
               !
               v6=in_av1(i,j,k)*v(i,j,k+1)+(1-in_av1(i,j,k))*vav(i,j)
               v5=in_in1(i,j,k)*v(i,j,k-1)+(1-in_in1(i,j,k))*vin(i,j)
               !
               w2=in_dx1(i,j,k)*w(i+1,j,k)+(1-in_dx1(i,j,k))*wdx(j,k)
               w1=in_sn1(i,j,k)*w(i-1,j,k)+(1-in_sn1(i,j,k))*wsn(j,k)
               !
               w4=in_sp1(i,j,k)*w(i,j+1,k)+(1-in_sp1(i,j,k))*wsp(i,k)
               w3=in_st1(i,j,k)*w(i,j-1,k)+(1-in_st1(i,j,k))*wst(i,k)
               !
               w6=in_av1(i,j,k)*w(i,j,k+1)+(1-in_av1(i,j,k))*wav(i,j)
               w5=in_in1(i,j,k)*w(i,j,k-1)+(1-in_in1(i,j,k))*win(i,j)
               !
               ! compute nu_T*dudy
               !
               gra1(i,j,k)=(csy(i,j,k)*.5*(u(i,j,k)+u2)- &
                  csy(i-1,j,k)*.5*(u(i,j,k)+u1)+ &
                  ety(i,j,k)*.5*(u(i,j,k)+u4)- &
                  ety(i,j-1,k)*.5*(u(i,j,k)+u3)+ &
                  zty(i,j,k)*.5*(u(i,j,k)+u6)- &
                  zty(i,j,k-1)*.5*(u(i,j,k)+u5)) &
                  *(annit(i,j,k)-an)/giac(i,j,k)
               !
               ! compute nu_T*dvdy
               !
               gra2(i,j,k)=(csy(i,j,k)*.5*(v(i,j,k)+v2)- &
                  csy(i-1,j,k)*.5*(v(i,j,k)+v1)+ &
                  ety(i,j,k)*.5*(v(i,j,k)+v4)- &
                  ety(i,j-1,k)*.5*(v(i,j,k)+v3)+ &
                  zty(i,j,k)*.5*(v(i,j,k)+v6)- &
                  zty(i,j,k-1)*.5*(v(i,j,k)+v5)) &
                  *(annit(i,j,k)-an)/giac(i,j,k)
               !
               ! compute nu_T*dwdy
               !
               gra3(i,j,k)=(csy(i,j,k)*.5*(w(i,j,k)+w2)- &
                  csy(i-1,j,k)*.5*(w(i,j,k)+w1)+ &
                  ety(i,j,k)*.5*(w(i,j,k)+w4)- &
                  ety(i,j-1,k)*.5*(w(i,j,k)+w3)+ &
                  zty(i,j,k)*.5*(w(i,j,k)+w6)- &
                  zty(i,j,k-1)*.5*(w(i,j,k)+w5)) &
                  *(annit(i,j,k)-an)/giac(i,j,k)

            !
            enddo
         enddo
      enddo
      !
      ! make gra1, gra2 and gra3 at the border known between the procs
      !

      if(leftpem /= MPI_PROC_NULL) then
         call MPI_SSEND(gra1(1,1,kparasta),jx*jy,MPI_REAL_SD,leftpem ,tagls,MPI_COMM_WORLD,ierr)
         call MPI_SSEND(gra2(1,1,kparasta),jx*jy,MPI_REAL_SD,leftpem ,tagls,MPI_COMM_WORLD,ierr)
         call MPI_SSEND(gra3(1,1,kparasta),jx*jy,MPI_REAL_SD,leftpem ,tagls,MPI_COMM_WORLD,ierr)
      endif
      if(rightpem /= MPI_PROC_NULL) then
         call MPI_RECV(gra1(1,1,kparaend+1),jx*jy,MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
         call MPI_RECV(gra2(1,1,kparaend+1),jx*jy,MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
         call MPI_RECV(gra3(1,1,kparaend+1),jx*jy,MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
      endif

      if(leftpem /= MPI_PROC_NULL) then
      !      call MPI_WAIT(req1,istatus,ierr)
      !      call MPI_WAIT(req3,istatus,ierr)
      !      call MPI_WAIT(req5,istatus,ierr)
      endif
      if(rightpem /= MPI_PROC_NULL) then
      !      call MPI_WAIT(req2,istatus,ierr)
      !      call MPI_WAIT(req4,istatus,ierr)
      !      call MPI_WAIT(req6,istatus,ierr)
      endif

      ! send k=jz to P0 and k=1 to Pn-1 to compute at the border
      ! cgra3 in flu_turbo


      if (myid.eq.nproc-1) then

         call MPI_SENDRECV(gra1(1,1,jz),jx*jy,MPI_REAL_SD,0,101,gra1_appoggio(1,1,1),jx*jy, &
         MPI_REAL_SD,0,201,MPI_COMM_WORLD,status,ierr)

         call MPI_SENDRECV(gra2(1,1,jz),jx*jy,MPI_REAL_SD,0,102,gra2_appoggio(1,1,1),jx*jy, &
         MPI_REAL_SD,0,202,MPI_COMM_WORLD,status,ierr)

         call MPI_SENDRECV(gra3(1,1,jz),jx*jy,MPI_REAL_SD,0,103,gra3_appoggio(1,1,1),jx*jy, &
         MPI_REAL_SD,0,203,MPI_COMM_WORLD,status,ierr)

      endif

      if (myid.eq.0) then

         call MPI_SENDRECV(gra1(1,1,1),jx*jy,MPI_REAL_SD,nproc-1,201,gra1_appoggio(1,1,jz),jx*jy, &
         MPI_REAL_SD,nproc-1,101,MPI_COMM_WORLD,status,ierr)

         call MPI_SENDRECV(gra2(1,1,1),jx*jy,MPI_REAL_SD,nproc-1,202,gra2_appoggio(1,1,jz),jx*jy, &
         MPI_REAL_SD,nproc-1,102,MPI_COMM_WORLD,status,ierr)

         call MPI_SENDRECV(gra3(1,1,1),jx*jy,MPI_REAL_SD,nproc-1,203,gra3_appoggio(1,1,jz),jx*jy, &
         MPI_REAL_SD,nproc-1,103,MPI_COMM_WORLD,status,ierr)

      endif

      return
   end

   !***********************************************************************
   subroutine jord3(in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1, &
      kparasta,kparaend,rightpe,leftpe, &
      tagls,taglr,tagrs,tagrr,rightpem,leftpem)
      !***********************************************************************
      ! compute the gradient dudz, dvdz, dwdz
      !

      !----------------------------------------------------------------------
      !     array declaration
      integer ierr,myid,nproc,status(MPI_STATUS_SIZE)
      integer ncolperproc,kparasta,kparaend,m
      !
      integer tagls,taglr,tagrs,tagrr
      integer leftpe,rightpe
      integer leftpem,rightpem
      integer req1,req2,req3,req4,req5,req6
      integer istatus(MPI_STATUS_SIZE)
      !
      integer i,j,k
      integer in_dx1(n1,n2,kparasta:kparaend) ! n3)
      integer in_sn1(n1,n2,kparasta:kparaend) !n3)
      integer in_sp1(n1,n2,kparasta:kparaend) !n3)
      integer in_st1(n1,n2,kparasta:kparaend) !n3)
      integer in_av1(n1,n2,kparasta:kparaend) !n3)
      integer in_in1(n1,n2,kparasta:kparaend) !n3)
      !
      real an
      real u1,u2,u3,u4,u5,u6
      real v1,v2,v3,v4,v5,v6
      real w1,w2,w3,w4,w5,w6
      !-----------------------------------------------------------------------
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      !-----------------------------------------------------------------------
      an=1./re
      !
      do k=kparasta,kparaend
         do j=1,jy
            do i=1,jx
               !
               u2=in_dx1(i,j,k)*u(i+1,j,k)+(1-in_dx1(i,j,k))*udx(j,k)
               u1=in_sn1(i,j,k)*u(i-1,j,k)+(1-in_sn1(i,j,k))*usn(j,k)
               !
               u4=in_sp1(i,j,k)*u(i,j+1,k)+(1-in_sp1(i,j,k))*usp(i,k)
               u3=in_st1(i,j,k)*u(i,j-1,k)+(1-in_st1(i,j,k))*ust(i,k)
               !
               u6=in_av1(i,j,k)*u(i,j,k+1)+(1-in_av1(i,j,k))*uav(i,j)
               u5=in_in1(i,j,k)*u(i,j,k-1)+(1-in_in1(i,j,k))*uin(i,j)
               !
               v2=in_dx1(i,j,k)*v(i+1,j,k)+(1-in_dx1(i,j,k))*vdx(j,k)
               v1=in_sn1(i,j,k)*v(i-1,j,k)+(1-in_sn1(i,j,k))*vsn(j,k)
               !
               v4=in_sp1(i,j,k)*v(i,j+1,k)+(1-in_sp1(i,j,k))*vsp(i,k)
               v3=in_st1(i,j,k)*v(i,j-1,k)+(1-in_st1(i,j,k))*vst(i,k)
               !
               v6=in_av1(i,j,k)*v(i,j,k+1)+(1-in_av1(i,j,k))*vav(i,j)
               v5=in_in1(i,j,k)*v(i,j,k-1)+(1-in_in1(i,j,k))*vin(i,j)
               !
               w2=in_dx1(i,j,k)*w(i+1,j,k)+(1-in_dx1(i,j,k))*wdx(j,k)
               w1=in_sn1(i,j,k)*w(i-1,j,k)+(1-in_sn1(i,j,k))*wsn(j,k)
               !
               w4=in_sp1(i,j,k)*w(i,j+1,k)+(1-in_sp1(i,j,k))*wsp(i,k)
               w3=in_st1(i,j,k)*w(i,j-1,k)+(1-in_st1(i,j,k))*wst(i,k)
               !
               w6=in_av1(i,j,k)*w(i,j,k+1)+(1-in_av1(i,j,k))*wav(i,j)
               w5=in_in1(i,j,k)*w(i,j,k-1)+(1-in_in1(i,j,k))*win(i,j)
               !
               ! compute nu_T*dudz
               !
               gra1(i,j,k)=(csz(i,j,k)*.5*(u(i,j,k)+u2)- &
                  csz(i-1,j,k)*.5*(u(i,j,k)+u1)+ &
                  etz(i,j,k)*.5*(u(i,j,k)+u4)- &
                  etz(i,j-1,k)*.5*(u(i,j,k)+u3)+ &
                  ztz(i,j,k)*.5*(u(i,j,k)+u6)- &
                  ztz(i,j,k-1)*.5*(u(i,j,k)+u5)) &
                  *(annit(i,j,k)-an)/giac(i,j,k)
               !
               ! compute nu_T*dvdz
               !
               gra2(i,j,k)=(csz(i,j,k)*.5*(v(i,j,k)+v2)- &
                  csz(i-1,j,k)*.5*(v(i,j,k)+v1)+ &
                  etz(i,j,k)*.5*(v(i,j,k)+v4)- &
                  etz(i,j-1,k)*.5*(v(i,j,k)+v3)+ &
                  ztz(i,j,k)*.5*(v(i,j,k)+v6)- &
                  ztz(i,j,k-1)*.5*(v(i,j,k)+v5)) &
                  *(annit(i,j,k)-an)/giac(i,j,k)
               !
               ! compute nu_T*dwdz
               !
               gra3(i,j,k)=(csz(i,j,k)*.5*(w(i,j,k)+w2)- &
                  csz(i-1,j,k)*.5*(w(i,j,k)+w1)+ &
                  etz(i,j,k)*.5*(w(i,j,k)+w4)- &
                  etz(i,j-1,k)*.5*(w(i,j,k)+w3)+ &
                  ztz(i,j,k)*.5*(w(i,j,k)+w6)- &
                  ztz(i,j,k-1)*.5*(w(i,j,k)+w5)) &
                  *(annit(i,j,k)-an)/giac(i,j,k)
            !
            enddo
         enddo
      enddo
      !
      ! make gra1, gra2 and gra3 at the border known between the procs
      !

      if(leftpem /= MPI_PROC_NULL) then
         call MPI_SSEND(gra1(1,1,kparasta),jx*jy,MPI_REAL_SD,leftpem ,tagls,MPI_COMM_WORLD,ierr)
         call MPI_SSEND(gra2(1,1,kparasta),jx*jy,MPI_REAL_SD,leftpem ,tagls,MPI_COMM_WORLD,ierr)
         call MPI_SSEND(gra3(1,1,kparasta),jx*jy,MPI_REAL_SD,leftpem ,tagls,MPI_COMM_WORLD,ierr)
      endif
      if(rightpem /= MPI_PROC_NULL) then
         call MPI_RECV(gra1(1,1,kparaend+1),jx*jy,MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
         call MPI_RECV(gra2(1,1,kparaend+1),jx*jy,MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
         call MPI_RECV(gra3(1,1,kparaend+1),jx*jy,MPI_REAL_SD,rightpem,tagrr,MPI_COMM_WORLD,status,ierr)
      endif

      if(leftpem /= MPI_PROC_NULL) then
      !      call MPI_WAIT(req1,istatus,ierr)
      !      call MPI_WAIT(req3,istatus,ierr)
      !      call MPI_WAIT(req5,istatus,ierr)
      endif
      if(rightpem /= MPI_PROC_NULL) then
      !      call MPI_WAIT(req2,istatus,ierr)
      !      call MPI_WAIT(req4,istatus,ierr)
      !      call MPI_WAIT(req6,istatus,ierr)
      endif

      ! send k=jz to P0 and k=1 to Pn-1 to compute at the border
      ! cgra3 in flu_turbo

      if (myid.eq.nproc-1) then

         call MPI_SENDRECV(gra1(1,1,jz),jx*jy,MPI_REAL_SD,0,101,gra1_appoggio(1,1,1),jx*jy, &
         MPI_REAL_SD,0,201,MPI_COMM_WORLD,status,ierr)

         call MPI_SENDRECV(gra2(1,1,jz),jx*jy,MPI_REAL_SD,0,102,gra2_appoggio(1,1,1),jx*jy, &
         MPI_REAL_SD,0,202,MPI_COMM_WORLD,status,ierr)

         call MPI_SENDRECV(gra3(1,1,jz),jx*jy,MPI_REAL_SD,0,103,gra3_appoggio(1,1,1),jx*jy, &
         MPI_REAL_SD,0,203,MPI_COMM_WORLD,status,ierr)

      endif

      if (myid.eq.0) then

         call MPI_SENDRECV(gra1(1,1,1),jx*jy,MPI_REAL_SD,nproc-1,201,gra1_appoggio(1,1,jz),jx*jy, &
         MPI_REAL_SD,nproc-1,101,MPI_COMM_WORLD,status,ierr)

         call MPI_SENDRECV(gra2(1,1,1),jx*jy,MPI_REAL_SD,nproc-1,202,gra2_appoggio(1,1,jz),jx*jy, &
         MPI_REAL_SD,nproc-1,102,MPI_COMM_WORLD,status,ierr)

         call MPI_SENDRECV(gra3(1,1,1),jx*jy,MPI_REAL_SD,nproc-1,203,gra3_appoggio(1,1,jz),jx*jy, &
         MPI_REAL_SD,nproc-1,103,MPI_COMM_WORLD,status,ierr)

      endif


      return
   end

end module jord_module
