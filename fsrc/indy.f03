!***********************************************************************
subroutine indy(nlevel,jxc,jyc,jzc)
   !***********************************************************************
   ! find the multigrid level and number of cells for each level
   ! it is necessary to solve poisson equation with multigrid
   !
   use mysending, only: myid, nproc
   use scala3, only: n1,n2,n3
       
   implicit none

   !-----------------------------------------------------------------------
   !      variables declaration
   integer i,j,k,kx,ky,kz
   integer mezx,mezy,mez1,mez2,mez3
   integer,intent(out) :: nlevel
   integer,intent(out) :: jxc(0:4),jyc(0:4),jzc(0:4)

   !-----------------------------------------------------------------------
   ! find number of level in x
   !
   mez1=n1
   mez3=mez1
   mez2=mez1
   !
   k=0
   do while (mez3.eq.mez1.and.k.lt.20)
      !
      mez1=mez2
      k=k+1
      mez2=mez1/2
      mez3=2*mez2
   !
   end do
   mezx=mez1
   kx=k
   if (mezx.eq.1) kx=kx-1
   !
   ! find number of level in y
   !
   mez1=n2
   mez3=mez1
   mez2=mez1
   !
   k=0
   do while (mez3.eq.mez1.and.k.lt.4)
      !
      mez1=mez2
      k=k+1
      mez2=mez1/2
      mez3=2*mez2
   !
   end do
   mezy=mez1
   ky=k
   if (mezy.eq.1) ky=ky-1
   !
   ! find number of level in z
   !
   mez1=n3/nproc
   mez3=mez1
   mez2=mez1
   !
   k=0
   do while (mez3.eq.mez1.and.k.lt.4)
      !
      mez1=mez2
      k=k+1
      mez2=mez1/2
      mez3=2*mez2
   !
   end do
   mezy=mez1
   kz=k
   if (mezy.eq.1) kz=kz-1
   !
   nlevel=min(kx,ky,kz)

   if(myid==0)write(*,*)'NLEVEL:',nlevel,kx,ky,kz
   !
   ! find number of cells for each level
   !
   jxc(0)=0
   jyc(0)=0
   jzc(0)=0
   !
   jzc(1)=n3
   do j=2,nlevel
      jzc(j)=jzc(j-1)/2
   end do
   !
   jyc(1)=n2
   do j=2,nlevel
      jyc(j)=jyc(j-1)/2
   end do
   !
   jxc(1)=n1
   do i=2,nlevel
      jxc(i)=jxc(i-1)/2
   end do
   !
   if(myid==0)write(*,*)'check1:'
   return

end
