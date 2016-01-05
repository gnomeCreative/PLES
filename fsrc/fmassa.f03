!***********************************************************************
subroutine fmassa(dbbx,dbby,dbbz,rich, &
   bbx,bby,bbz, &
   bcsi,beta,bzet, &
   myid,nproc,kparasta,kparaend, &
   ti,langyes,wavebk,ktime)
   !***********************************************************************
   ! compute covariant component of bodyforce
   !
   use myarrays_WB
   use myarrays_LC
   use myarrays_cor
   use myarrays_metri3
   use myarrays_velo3
   use myarrays_moisture
   !
   use scala3
   use tipologia
   !
   use mpi

   implicit none

   !-----------------------------------------------------------------------
   !     array declaration
   !

   integer ierr,myid,nproc
   integer kparasta,kparaend
   !
   integer i,j,k
   integer langyes,wavebk,ktime
   real dbbx,dbby,dbbz,ti
   real bbx,bby,bbz,rich
   real bcsi(n1,n2,kparasta:kparaend)
   real beta(n1,n2,kparasta:kparaend)
   real bzet(n1,n2,kparasta:kparaend)
      
   real tm_loc,tm_tot
   real qm_loc,qm_tot
      
   !real omega1,omega2,omega3
   real earth_rotation,ang,pi
   real beta_exp
      
   !-------ANDREA
   real smooth
   !-----------------------------------------------------------------------
   beta_exp = 1.
   earth_rotation = 0.000073
   pi = acos(-1.)
   ang = LATITUDE*pi/180.
   if(abs(LATITUDE).lt. 1.d-10)then
      omega1 = 0.
      omega2 = 0.
      omega3 = 0.
   else
      omega1 = 0.
      omega2 = 2.*earth_rotation*cos(ang)
      omega3 = 2.*earth_rotation*sin(ang)
   end if
   !-----------------------------------------------------------------------
   ! omega1, omega2, omega3 are the three component of the rotation vector
   ! must be defined in Agenerale.in
   smooth=1.0
   if(ktime.lt.50)smooth=(real(ktime)/50.)
   if(myid.eq.0)write(*,*)'smooth',smooth
      
   if(imoist==0)then
      do k=kparasta,kparaend
         do j=1,jy
            do i=1,jx

               bcsi(i,j,k) = bbx + dbbx &
                  +omega2* w(i,j,k)  &
                  -omega3* v(i,j,k)

               beta(i,j,k) = bby + dbby &
                  +omega3* u(i,j,k) &
                  -omega1* w(i,j,k) &
               !     >              -rich*beta_exp*rhov(1,i,j,k) &
               +smooth*rich*(  &
                  1.5*1E-04*(rhov(1,i,j,k)-6.7)) !temperatura minima MONFY
               !     >              - 7.5*1E-04*(rhov(2,i,j,k)-0.0)) !salinità di riferimento nulla

               bzet(i,j,k) = bbz + dbbz &
                  +omega1* v(i,j,k) &
                  -omega2* u(i,j,k)

            end do
         end do
      end do
   end if
   !-----------------------------------------------------------------------
   if(imoist==1)then
      ! Compute the mean potential T at each plane
      tpotm=0.
      do j=1,jy
         tm_loc = 0.
         tm_tot = 0.
         do k=kparasta,kparaend
            do i=1,jx
               tm_loc = tm_loc + rhov(1,i,j,k)
            end do
         end do
         !        sum on Tm_loc
         call MPI_ALLREDUCE(tm_loc,tm_tot,1,MPI_REAL_SD, &
            MPI_SUM,MPI_COMM_WORLD,ierr)
         tpotm(j) = tm_tot/real(jx*jz)
      enddo
      
      ! Compute the mean humidity at each plane
      qm=0.
      do j=1,jy
         qm_loc = 0.
         qm_tot = 0.
         do k=kparasta,kparaend
            do i=1,jx
               qm_loc = qm_loc + rhov(2,i,j,k)
            end do
         end do
         !        sum on qm_loc
         call MPI_ALLREDUCE(qm_loc,qm_tot,1,MPI_REAL_SD, &
            MPI_SUM,MPI_COMM_WORLD,ierr)
         qm(j) = qm_tot/real(jx*jz)
      enddo

      !      if(myid.eq.0)then
      !         write(*,*)'Tp mean ', tpotm
      !      end if
      
      do k=kparasta,kparaend
         do j=1,jy
            do i=1,jx

               bcsi(i,j,k) = bbx + dbbx &
                  +2.*omega2* w(i,j,k)  &
                  -2.*omega3* v(i,j,k)

               beta(i,j,k) = bby + dbby &
                  +2.*omega3* u(i,j,k) &
                  -2.*omega1* w(i,j,k) &
                  + rich*(  &
                  (rhov(1,i,j,k)-tpotm(j))/(tpotm(j)+Tref) &
                  + ( (Ma-Mv)/( Ma*qm(j)+Mv*(1-qm(j))) ) &
                  *(rhov(2,i,j,k)-qm(j)))

               bzet(i,j,k) = bbz + dbbz &
                  +2.*omega1* v(i,j,k) &
                  -2.*omega2* u(i,j,k)

            end do
         end do
      end do

   end if

   !-----------------------------------------------------------------------
   ! forcing with imposed pressure gradient (bbx,bby,bbz), that change at
   ! every iteration (dbbx, dbby, dbbz).
   ! Forcing likes tide is imposed like A*cos(omega*t)
   ! Coriolis rotation term + interaction with tide
   ! term for vertical stratification

   !      do k=kparasta,kparaend
   !      do j=1,jy
   !      do i=1,jx

   !      bcsi(i,j,k) = bbx + dbbx + A1*cos(omegaM2*ti)+
   !     >            ( 2.*omega2*(V0*sin(omegaM2*ti)-v(i,j,k))-
   !     >              2.*omega3*(W0*sin(omegaM2*ti)-w(i,j,k))  )
   !      beta(i,j,k) = bby + dbby + A3*cos(omegaM2*ti)+
   !     >            ( 2.*omega1*(W0*sin(omegaM2*ti)-w(i,j,k))-
   !     >              2.*omega2*(U0*sin(omegaM2*ti)-u(i,j,k))  )-
   !     >              rich*rho(i,j,k)
   !      bzet(i,j,k) = bbz + dbbz + A2*cos(omegaM2*ti)+
   !     >            ( 2.*omega3*(U0*sin(omegaM2*ti)-u(i,j,k))-
   !     >              2.*omega1*(V0*sin(omegaM2*ti)-v(i,j,k))  )

   !      end do
   !      end do
   !      end do

   !
   !-----------------------------------------------------------------------
   ! wavebreaking at the surface

   if (wavebk.eq.1) then

      call gauss_random(ktime,n1,n2,n3, &
         myid,nproc,kparasta,kparaend)

      j=jy

      do k=kparasta,kparaend
         do i=1,jx
            bcsi(i,j,k) = bcsi(i,j,k) + alpha*u_att(i,k)*alpha*u_att(i,k)* &
               Fx(i,k)/(0.1*l_0)
            beta(i,j,k) = beta(i,j,k)
            bzet(i,j,k) = bzet(i,j,k) + alpha*w_att(i,k)*alpha*w_att(i,k)* &
               Fz(i,k)/(0.1*l_0)
         end do
      end do

   end if
   !
   !-----------------------------------------------------------------------
   ! Langmuir term
   if (langyes.eq.1) then

      do k=kparasta,kparaend
         do j=1,jy
            do i=1,jx
               bcsi(i,j,k) = bcsi(i,j,k) + w_drift(i,j,k)*vorty(i,j,k)
               beta(i,j,k) = beta(i,j,k) + u_drift(i,j,k)*vortz(i,j,k) - &
                  w_drift(i,j,k)*vortx(i,j,k)
               bzet(i,j,k) = bzet(i,j,k) - u_drift(i,j,k)*vorty(i,j,k)
            end do
         end do
      end do

   end if
   !
   !-----------------------------------------------------------------------
   !
   return
end
