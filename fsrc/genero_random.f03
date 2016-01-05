!***********************************************************************
subroutine genero_random(corr_factor,nreal,nstep,disturbance,istampo,myid)
   !***********************************************************************
   !
   implicit none
   ! program to illustrate the colored Gaussian Noise generator cgaussA
   ! The routine must be initialized with cgaus0A and calls a flat distribution
   ! random number generator available with most compilers or you can write your
   ! own. Here we used the routine RAN1 from Numerical Recipes 2nd Edition, by
   ! Press, Teukolsky, Vetterling, and Flannery.

   ! It now uses the F90 intrinsic subroutine RANDOM_NUMBER.

   ! The White Guassian noise generator GASDEV from Numerical Recipes was
   ! adapted to produce Colored Gaussian noise. The basic equations for this
   ! computation are presented in the article by
   ! Fox et al., Physical Review A vol.38(1988) page 5938.
   ! This code was [originally] compiled and tested with Microsoft Powerstation.

   ! It was modified by Walt Brainerd to be standard Fortran and
   ! compiled on NAGWare F90.

   integer i,j
   integer istampo,myid
   integer nreal,nstep,npts,idly
   real, allocatable :: eps(:,:),sum(:)
   real smean
   real deviation,variance
   real std,mean,cgaussA,dt
   real cortim
   real corr_factor
   real disturbance(nreal,-1:nstep*2)
   real start_time,end_time
   ! get input parameters (typical values shown)
   !c        open(1,file='cnoise.dat')
   !c        read(1,*)nreal             !number of realizations=1000
   !c        read(1,*)nstep             !max delay in corr. func=10
   !c        read(1,*)dt                !time step size=.5
   !c        read(1,*)cortim            !corr. time in the same units as DT=5

   !         nreal = 10
   !	 nstep = 1000
   dt = 1.
   !	 corr_factor = 0.
   cortim = corr_factor*dt
 
   allocate(eps(nreal,-1:nstep*2))
   allocate(sum(nreal))
   eps = 0.
   sum = 0.
	         
   ! initialize
   !         call cpu_time(start_time)
   call cgaus0A(dt,cortim)
   !	 call cpu_time(end_time)
   !	 write(*,*)'TIME'
   !	 write(*,'(1e18.10)')end_time-start_time
     
   ! store several series of Gaussian noise values in array EPS.
   !        call cpu_time(start_time)

   do i=1,nreal
      mean = 0.
      do j=1,nstep*2 !0,nstep*2
         eps(i,j) = cgaussA()
         mean = mean + eps(i,j)
      end do
      mean = mean/float(nstep*2)  !float(nstep*2+1)
	   
      do j=1,nstep*2 !0,nstep*2
         eps(i,j) = eps(i,j) - mean
      end do
   end do
   mean = 0.
        

   do i=1,nreal             !realizations
      sum(i) = 0.
      do j=1,nstep*2  !0,nstep*2          !time delays
         !cccc          eps(i,j)=cgaussA()
         !          write(100+i,*)j,eps(i,j)
         sum(i) = sum(i) + eps(i,j)
      enddo
      sum(i)=sum(i)/float(nstep*2)   !float(nstep*2+1)
	 
      deviation = 0.
      do j=1,nstep*2  !0,nstep*2          !time delays
         deviation = deviation + (eps(i,j)-sum(i))*(eps(i,j)-sum(i))
      enddo
      deviation = deviation/float(nstep*2) !float(nstep*2+1)
      deviation = sqrt(deviation)

      variance = 1
      variance = variance/deviation
      do j=1,nstep*2 !0,nstep*2
         eps(i,j)=eps(i,j)*variance
      end do
   enddo
	

   ! calculate the autocorrelation function in variable MEAN.
   !        call cpu_time(start_time)
   npts=nstep*nreal
   do idly=1,nstep !0,nstep
      mean=0.
      std=0.
      do i=1,nreal
         do j=1,nstep !0,nstep
            mean=mean+real(eps(i,j)*eps(i,j+idly))
         enddo
      enddo
      mean=mean/real(npts)
      smean=sngl(mean)          !single precision speeds up calculations

      ! calculate the error in autocorrelation function in variable STD.
      do i=1,nreal
         do j=1,nstep !0,nstep
            std=std+real((eps(i,j)*eps(i,j+idly)-smean)**2.)
         enddo
      enddo
      std=sqrt(std)/real(npts) !dble(npts-1.)
      write(200+myid,*)idly,mean,std            !output results
   enddo
   !	call cpu_time(end_time)
   !	write(*,*)'TIME'
   !	write(*,'(1e18.10)')end_time-start_time


   !       storage
   do i=1,nreal
      do j=0,nstep*2
         disturbance(i,j) = eps(i,j)
      end do
   end do


   deallocate(eps)
   deallocate(sum)
   !        end
   return
end
!==========================================================================
! initialize the RNG's
! and set the color of gaussian noise
! DT is the time step used in whatever process the colored Gaussian noise
!   is used.
! CORTIM is correlation time in the same units as time step DT.
! WHITE=.true. means generate white gaussian noise which happens when
!   CORTIM=0. This flag is used in cgaussA.
! Here we use the flat distribution RAN1 also taken from Numerical Recipe
! but any other good flat distribution random number generator will do.

subroutine cgaus0A(dt,cortim)
   !        double precision ran1,cape,dt,l1me2,cgaussA
   real cape,dt,l1me2,cgaussA
   real cortim,x
   logical white
   common /color/l1me2,cape,white
   if(cortim.eq.0.)then
      white=.true.
      l1me2=-2.000                        !white noise
      cape=0.0
   else
      white=.false.
      cape=exp(-dt/real(cortim))
      !parameter needed in cgaussA
      l1me2=-(real(1.)-cape*cape)*real(2./cortim)
   endif
   !        idum=-1
   !        x=ran1(idum)            !initialize flat rng
   x=cgaussA()            !initialize cgaussA value
   return
END

!==========================================================================
! Program to produce exponentially correlated colored (Gaussian) noise.
! based on Fox et al Physical Review A vol.38(1988)5938 and
! modification of GASDEV from Numerical Recipes for Fortran(2nd ed.pg279)

! CAPE is capital E in the article by Fox et. al.
! PREV is the previous value of cgaussA used in the next iteration
! L1ME2 is the main parameters causing colored noise in Fox et al
!       and represents (lamda*(1-E)^2). Ditto for H in that article.

! routine is illustrated in Double Precision in case it is needed in this
! mode, otherwise all Double Precision variables maybe changed to REAL
! but the corresponding changes must be made to cgaus0A and the calling
! programs.


real FUNCTION cgaussA()

   Implicit none

   !      INTEGER idum,iset
   INTEGER iset
   logical white
   !      double precision  fac,gset,rsq,v1,v2,ran1,l1me2,h,cape
   real  fac,gset,rsq,v1,v2,l1me2,h,cape,prev
   common /color/l1me2,cape,white
      
   SAVE iset,gset,prev
   DATA iset/0/
   DATA prev/0.0d0/

   if (iset.eq.0) then
      !1       v1=2.*ran1(idum)-1.
1     call random_number(v1)
      v1=2.*v1-1
      !        v2=2.*ran1(idum)-1.
      call random_number(v2)
      v2=2.*v2-1
      rsq=v1**2.+v2**2.
      if(rsq>=1. .or. rsq==0.)goto 1
      !took out sqrt(2) vs eq(28) Fox etal
      fac=sqrt(l1me2*log(rsq)/rsq)
      gset=v1*fac
      h=v2*fac
      iset=1
   else
      h=gset
      iset=0
   endif
      
   if(white)then  !please note that the time step vs its sqrt
      cgaussA=h      !in integration is previously set in PARAM
   else
      cgaussA=prev*cape+h
      prev=cgaussA
   endif
      
   return
end
