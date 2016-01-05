!***********************************************************************
subroutine indy(jx,jy,jz,nlevel,jxc,jyc,jzc)
   !***********************************************************************
   ! find the multigrid level and number of cells for each level
   ! it is necessary to solve poisson equation with multigrid
   !
   use mysending
       
   implicit none
   !
   !-----------------------------------------------------------------------
   !      variables declaration
   integer i,j,k,kx,ky,kz,jx,jy,jz,nlevel
   integer mezx,mezy,mez1,mez2,mez3
   integer jxc(0:4),jyc(0:4),jzc(0:4)
   !-----------------------------------------------------------------------
   ! find number of level in x
   !
   mez1=jx
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
   mez1=jy
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
   mez1=jz/nproc
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

   if(myid==0)write(*,*)'NL:',nlevel,kx,ky,kz
   !
   ! find number of cells for each level
   !
   jxc(0)=0
   jyc(0)=0
   jzc(0)=0
   !
   jzc(1)=jz
   do j=2,nlevel
      jzc(j)=jzc(j-1)/2
   end do
   !
   jyc(1)=jy
   do j=2,nlevel
      jyc(j)=jyc(j-1)/2
   end do
   !
   jxc(1)=jx
   do i=2,nlevel
      jxc(i)=jxc(i-1)/2
   end do
   !
   return
end
