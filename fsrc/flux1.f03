!***********************************************************************
subroutine flux1(rc,cgra1,r,insc,tipo,body)
   !***********************************************************************
   ! compute explict flux on csi component:
   !
   ! convective  uc*(u,v,w)  with cenetered scheme or quick
   ! diffusive   nni*g12*d(u,v,w)/d(eta)
   ! diffusive   nni*g13*d(u,v,w)/d(zita)
   !
   use mysending
   use myarrays_LC
   use myarrays_metri3
   !
   use scala3
   use period
   use convex
   !
   use mpi

   implicit none

   !-----------------------------------------------------------------------
   !     array declaration
   integer ierr,insc
   integer ncolperproc,m
   integer kparastak,kparaendk
   double precision bulk_loc
   !
   integer i,j,k,jj,kk,is,iss
   real xf,yf,zf,xc1,yc1,zc1,xc2,yc2,zc2,ddi,dsi
   real fg,r0,r1,r2,an
   real rc(0:n1,1:n2,kparasta  :kparaend)
   real cgra1(0:n1,n2,kparasta-1:kparaend+1)
   real r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
   real ravanti,rindietro1,rindietro2
   !      integer tipo2(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
   !      integer tipo2(0:n1+1,0:n2+1,0:n3+1) c7

   integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
      
   real phic(n1,n2,kparasta:kparaend)
   real den, num, myv
   integer body
   !-----------------------------------------------------------------------
   ! CONVECTIVE TERM Uc*u
   ! implemented with central scheme or with quick depending on setting
   ! in point closer to wall it uses uf as ghost cell
   !
   !
   !     side left and right
   do i=1,ip
      !
      do k=kparasta,kparaend
         do j=1,jy
            !
            f1(0,j,k) = (rc(0 ,j,k)+ucs(0,j,k))*r(0   ,j,k)- &
               cgra1(0 ,j,k)
            f1(jx,j,k)= (rc(jx,j,k)+ucs(jx,j,k))*r(jx+1,j,k)- &
               cgra1(jx,j,k)
         !
         enddo
      enddo

   enddo
   !
   !     into the field
   do k=kparasta,kparaend
      do j=1,jy
         do i=ip,jx-ip
            !
            f1(i,j,k)=(rc(i,j,k)+ucs(i,j,k)) &
               *.5*(r(i,j,k)+r(i+1,j,k))-cgra1(i,j,k)
         !
         end do
      end do
   end do
   !
   !
   !
   ! quick
   !
   if(insc.eq.1)then
      !     sides 1 and 2
      do k=kparasta,kparaend     
         do j=1,jy
      
            i=1
            if(rc(i,j,k).gt.0.)then
               ravanti    =        r(2,j,k)
               rindietro1 =        r(1,j,k)
               rindietro2 = ip*(2.*r(0,j,k) - r(1,j,k)) &
                  +(1.-ip)*r(0,j,k)
            else
               ravanti    =    r(1,j,k)
               rindietro1 =    r(2,j,k)
               rindietro2 =    r(3,j,k)
            end if
        
            if(tipo(i,j,k)==2)then
               f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                  *.125*( -ravanti +2.*rindietro1 - rindietro2 )
            end if
 
            i=jx-1
            if(rc(i,j,k).gt.0.)then
               ravanti    = r(jx  ,j,k)
               rindietro1 = r(jx-1,j,k)
               rindietro2 = r(jx-2,j,k)
            else
               ravanti    =        r(jx-1,j,k)
               rindietro1 =        r(jx  ,j,k)
               rindietro2 = ip*(2.*r(jx+1,j,k) - r(jx,j,k)) &
                  +(1.-ip)*r(jx+1,j,k)
            end if

            if(tipo(i,j,k)==2)then
               f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                  *.125*( -ravanti +2.*rindietro1 - rindietro2 )
            end if
         end do
      end do    
       
      !     into the field
      do k=kparasta,kparaend
         do j=1,jy
            do i=2,jx-2
               !Giulia non comprende tutti i casi, da un lato è quick dall'altro dc
               !      if(tipo(i,j,k)==2)then
               !     if (rc(i,j,k).gt.0.) then
               !      f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*
               !     > .125*(-r(i+1,j,k)+2*r(i,j,k)-r(i-1,j,k))
               !      else
               !      f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*
               !     > .125*(-r(i,j,k)+2*r(i+1,j,k)-r(i+2,j,k))
               !      end if
               !         end if
               !
               if (rc(i,j,k).gt.0.) then

                  if(tipo(i,j,k)==1)then !ib cell
                     if(tipo(i+1,j,k)==0)then !solido davanti
                        f1(i,j,k)=0.
                     elseif(tipo(i-1,j,k)==0)then !uso diff centrate
                        f1(i,j,k)=f1(i,j,k)
                     else ! solido sta nell'altra direzione quick
                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                           *.125*(-r(i+1,j,k)+2*r(i,j,k)-r(i-1,j,k))

                     endif !solido avanti/indietro/lato
                  else !i è fluido allora quick normale

                     f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                        *.125*(-r(i+1,j,k)+2*r(i,j,k)-r(i-1,j,k))
                  endif !tipo
    
               else !rc(i,j,k).le.0.
   
                  if(tipo(i+1,j,k)==1)then
                     if(tipo(i,j,k)==0)then !solido indietro
                        f1(i,j,k)=0.
                     elseif(tipo(i+2,j,k)==0)then !uso diff centrate
                        f1(i,j,k)=f1(i,j,k)
                     else !solido di lato
                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                           *.125*(-r(i,j,k)+2*r(i+1,j,k)-r(i+2,j,k))
                     endif !solido avanti/indietro/lato
                  else !i+1 è fluido allora quick normale
                     f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                        *.125*(-r(i,j,k)+2*r(i+1,j,k)-r(i+2,j,k))
                  endif !tipo
	
               endif !rc

            end do
         end do
      end do
   ! SMART/QUICK modificato upwind
   elseif(insc.eq.2)then

      ! GIULIA UCS - CGRA
      !     side left and right
      do i=1,ip
         !
         do k=kparasta,kparaend
            do j=1,jy
               !
               f1(0,j,k) = (rc(0 ,j,k)+ucs(0,j,k))*r(0   ,j,k)- &
                  cgra1(0 ,j,k)
               f1(jx,j,k)= (rc(jx,j,k)+ucs(jx,j,k))*r(jx+1,j,k)- &
                  cgra1(jx,j,k)
            !
            enddo
         enddo

      enddo
      !
      !     into the field
      do k=kparasta,kparaend
         do j=1,jy
            do i=ip,jx-ip
               !
               f1(i,j,k)=ucs(i,j,k) &
                  *.5*(r(i,j,k)+r(i+1,j,k))-cgra1(i,j,k)
            !
            end do
         end do
      end do


      !     sides 1 and 2  sia per bodyforce 0 o 1
      do k=kparasta,kparaend
         do j=1,jy
            i=1
            if(rc(i,j,k).gt.0.)then
               ravanti    =        r(2,j,k)
               rindietro1 =        r(1,j,k)
               rindietro2 = ip*(2.*r(0,j,k) - r(1,j,k)) &
                  +(1.-ip)*r(0,j,k)
            else
               ravanti    =    r(1,j,k)
               rindietro1 =    r(2,j,k)
               rindietro2 =    r(3,j,k)
            end if

            f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
               *.125*( 3.*ravanti +6.*rindietro1 - rindietro2 )


            i=jx-1
            if(rc(i,j,k).gt.0.)then
               ravanti    =    r(jx   ,j,k)
               rindietro1 =    r(jx-1,j, k)
               rindietro2 =    r(jx-2,j, k)
            else
               ravanti    =        r(jx-1,j,k)
               rindietro1 =        r(jx  ,j,k)
               rindietro2 = ip*(2.*r(jx+1,j,k) - r(jx,j,k)) &
                  +(1.-ip)*r(jx+1,j,k)
            end if


            f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
               *.125*( 3.*ravanti +6.*rindietro1 - rindietro2 )

         end do
      end do
    
      
      !     into the field without bodyforce
      if(body.eq.0)then
         do k=kparasta,kparaend
            do j=1,jy
               do i=2,jx-2
                  if(rc(i,j,k).gt.0.)then
 
                     f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                        *.125*(3.*r(i+1,j,k)+6*r(i,j,k)-r(i-1,j,k))



                  else  !(rc(i,j,k).le.0.)

                     f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                        *.125*(3.*r(i,j,k)+6*r(i+1,j,k)-r(i+2,j,k)) !my index +1

                  endif !rc</>0
               end do
            end do
         end do

      !
      else !bodyforce attivo



         do k=kparasta,kparaend
            do j=1,jy
               do i=2,jx-2
                  if(rc(i,j,k).gt.0.)then
                     if(tipo(i,j,k)==1)then !ib cell
                        if(tipo(i+1,j,k)==0)then !solido davanti

                           f1(i,j,k)=0.
                        ! dopo diventerà nullo

                        elseif(tipo(i-1,j,k)==0)then

                           f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)


                        else ! solido sta nell'altra direzione

                           f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                              *.125*(3.*r(i+1,j,k)+6*r(i,j,k)-r(i-1,j,k))

                        endif !solido avanti/indietro/lato


                     else !i non è IB procedo come se non ci fossero gli ib
                        ! se i-1 è ib faccio la correzione per evitare la scia
                        !	      if(tipo(i-1,j,k)==1)then !la cella prima è un ib
                        !		if(tipo(i-2,j,k)==0)then !perché considero anche questa?
                        !   		   rindietro2 =r(i,j,k)

                        !		endif
                        !	      endif

                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                           *.125*(3.*r(i+1,j,k)+6*r(i,j,k)-r(i-1,j,k))

                     endif !tipo

                  else !rc(i,j,k).le.0.
   
                     if(tipo(i+1,j,k)==1)then
                        if(tipo(i,j,k)==0)then !solido indietro
                           f1(i,j,k)=0.
                        elseif(tipo(i+2,j,k)==0)then !solido avanti
                           f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)
                        else !solido di lato

                           f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                              *.125*(3.*r(i,j,k)+6*r(i+1,j,k)-r(i+2,j,k))

                        endif !solido avanti/indietro/lato

                     else !i non è IB procedo come se non ci fossero gli ib

                        ! se i+1 è ib faccio la correzione per evitare la scia
                        !	        if(tipo(i+1,j,k)==1)then
                        !		if(tipo(i+2,j,k)==0)then
                        !	          rindietro2=r(i+1,j,k)
                        !		endif
                        !		endif


                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                           *.125*(3.*r(i,j,k)+6*r(i+1,j,k)-r(i+2,j,k))
                     endif !tipo
	
                  endif !rc
      
               end do
            end do
         end do

      endif !body yes/no
   !	endif !SMART/quick modificato upwind
   ! SMART/QUICK modificato upwind
   elseif(insc.eq.3)then

      ! GIULIA UCS - CGRA
      !     side left and right
      do i=1,ip
         !
         do k=kparasta,kparaend
            do j=1,jy
               !
               f1(0,j,k) = (rc(0 ,j,k)+ucs(0,j,k))*r(0   ,j,k)- &
                  cgra1(0 ,j,k)
               f1(jx,j,k)= (rc(jx,j,k)+ucs(jx,j,k))*r(jx+1,j,k)- &
                  cgra1(jx,j,k)
            !
            enddo
         enddo

      enddo
      !
      !     into the field
      do k=kparasta,kparaend
         do j=1,jy
            do i=ip,jx-ip
               !
               f1(i,j,k)=ucs(i,j,k) &
                  *.5*(r(i,j,k)+r(i+1,j,k))-cgra1(i,j,k)
            !
            end do
         end do
      end do


      !     sides 1 and 2  sia per bodyforce 0 o 1
      do k=kparasta,kparaend
         do j=1,jy
            i=1
            if(rc(i,j,k).gt.0.)then
               rindietro1 =        r(1,j,k)
 		 
            else
               rindietro1 =    r(2,j,k)
            end if

            f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*rindietro1


            i=jx-1
            if(rc(i,j,k).gt.0.)then
               rindietro1 =    r(jx-1,j, k)
 
            else
 
               rindietro1 =        r(jx  ,j,k)
 	 	 
            end if

            f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*rindietro1

         end do
      end do
    
      
      !     into the field without bodyforce
      if(body.eq.0)then
         do k=kparasta,kparaend
            do j=1,jy
               do i=2,jx-2
                  if(rc(i,j,k).gt.0.)then
 
                     f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)

                  else  !(rc(i,j,k).le.0.)

                     f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)

                  endif !rc</>0
               end do
            end do
         end do

      !
      else !bodyforce attivo
         do k=kparasta,kparaend
            do j=1,jy
               do i=2,jx-2
                  if(rc(i,j,k).gt.0.)then
                     if(tipo(i,j,k)==1)then !ib cell
                        if(tipo(i+1,j,k)==0)then !solido davanti

                           f1(i,j,k)=0.
                        ! dopo diventerà nullo
    
                        else ! solido sta nell'altra direzione

                           f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)


                        endif !solido avanti/indietro/lato


                     else !i non è IB procedo come se non ci fossero gli ib
                        ! se i-1 è ib faccio la correzione per evitare la scia
                        !	      if(tipo(i-1,j,k)==1)then !la cella prima è un ib
                        !		if(tipo(i-2,j,k)==0)then !perché considero anche questa?
                        !   		   rindietro2 =r(i,j,k)

                        !		endif
                        !	      endif

                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)

                     endif !tipo

                  else !rc(i,j,k).le.0.
   
                     if(tipo(i+1,j,k)==1)then
                        if(tipo(i,j,k)==0)then !solido indietro
                           f1(i,j,k)=0.
                        !		elseif(tipo(i+2,j,k)==0)then !solido avanti
                        !	          f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)
                        else !solido di lato

                           f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)

                        endif !solido avanti/indietro/lato

                     else !i non è IB procedo come se non ci fossero gli ib

                        ! se i+1 è ib faccio la correzione per evitare la scia
                        !	        if(tipo(i+1,j,k)==1)then
                        !		if(tipo(i+2,j,k)==0)then
                        !	          rindietro2=r(i+1,j,k)
                        !		endif
                        !		endif


                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)
                     endif !tipo
	
                  endif !rc
      
               end do
            end do
         end do

      endif !body yes/no
   endif !SMART/quick modificato upwind

   !
   !
   !-----------------------------------------------------------------------
   ! DIFFUSIVE TERM NNI*G12*DU/D(ETA)
   !
   !     sides left and right
   is=0
   iss=0

   do i=1,2*ip
      !
      do k=kparasta,kparaend
         !
         do j=1+jp,jy-jp
            !
            fg=.5*(r(iss,j+1,k)-r(iss,j-1,k))
            f1(is,j,k)=-f1(is,j,k) &
               +annit(iss,j,k)*g12(is,j,k)*fg
         !
         end do
         !
         do jj=1,jp
            !
            !        check derivative at wall bottom and upper
            !
            j=1
            fg=.5*(-3.*r(iss,j,k)+4.*r(iss,j+1,k)-r(iss,j+2,k))
            f1(is,j,k)=-f1(is,j,k) &
               +annit(iss,j,k)*g12(is,j,k)*fg
            !
            j=jy
            fg=.5*(3.*r(iss,jy,k)-4.*r(iss,jy-1,k)+r(iss,jy-2,k))
            f1(is,j,k)=-f1(is,j,k) &
               +annit(iss,j,k)*g12(is,j,k)*fg
         !
         end do
      !
      enddo
      !
      is=jx
      iss=jx+1
   !
   enddo
   !
   !
   !     into the field
   do k=kparasta,kparaend
      do i=ip,jx-ip
         !
         do j=1+jp,jy-jp
            !
            r1=.5*(r(i,j-1,k)+r(i+1,j-1,k))
            r2=.5*(r(i,j+1,k)+r(i+1,j+1,k))
            fg=.5*(r2-r1)
            an=.5*(annit(i,j,k)+annit(i+1,j,k))
            f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*an*fg
         !
         end do
         !
         !     check derivative on sides bottom and upper
         !
         do jj=1,jp
            !
            j=1
            r0=.5*(r(i,j  ,k)+r(i+1,j  ,k))
            r1=.5*(r(i,j+1,k)+r(i+1,j+1,k))
            r2=.5*(r(i,j+2,k)+r(i+1,j+2,k))
            fg=.5*(-3.*r0+4.*r1-r2)
            an=.5*(annit(i,j,k)+annit(i+1,j,k))
            f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*an*fg
            !
            j=jy
            r0=.5*(r(i,jy  ,k)+r(i+1,jy  ,k))
            r1=.5*(r(i,jy-1,k)+r(i+1,jy-1,k))
            r2=.5*(r(i,jy-2,k)+r(i+1,jy-2,k))
            fg=.5*(3.*r0-4.*r1+r2)
            an=.5*(annit(i,j,k)+annit(i+1,j,k))
            f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*an*fg
         !
         end do
      !
      end do
   end do
   !
   !-----------------------------------------------------------------------
   ! DIFFUSIVE TERM NNI*G13*DU/D(ZITA)
   !
   ! sides left and right
   !
   ! define the limit depending on periodicity in z
   !
   if (myid.eq.0) then
      kparastak=kparasta+kp
      kparaendk=kparaend
   else if (myid.eq.nproc-1) then
      kparastak=kparasta
      kparaendk=kparaend-kp
   else
      kparastak=kparasta
      kparaendk=kparaend
   endif

   is=0
   iss=0

   do i=1,2*ip
      !
      do j=1,jy
         !
         do k=kparastak,kparaendk
            !
            fg=.5*(r(iss,j,k+1)-r(iss,j,k-1))
            f1(is,j,k)=f1(is,j,k) &
               +annit(iss,j,k)*g13(is,j,k)*fg
         !
         end do
         !
         do kk=1,kp
            !
            ! check derivative on sides back and front
            !
            if (myid.eq.0) then

               k=1
               fg=.5*(-3.*r(iss,j,k)+4.*r(iss,j,k+1)-r(iss,j,k+2))
               f1(is,j,k)=f1(is,j,k) &
                  +annit(iss,j,k)*g13(is,j,k)*fg

            endif
            !
            if (myid.eq.nproc-1) then

               k=jz
               fg=.5*(3.*r(iss,j,jz)-4.*r(iss,j,jz-1)+r(iss,j,jz-2))
               f1(is,j,k)=f1(is,j,k) &
                  +annit(iss,j,k)*g13(is,j,k)*fg

            endif
         !
         end do
      !
      enddo
      !
      is=jx
      iss=jx+1
   !
   enddo


   !
   ! into the field
   !
   do j=1,jy
      do i=ip,jx-ip
         do k=kparastak,kparaendk
            !
            !
            an=.5*(annit(i,j,k)+annit(i+1,j,k))
            r1=.5*(r(i,j,k-1)+r(i+1,j,k-1))
            r2=.5*(r(i,j,k+1)+r(i+1,j,k+1))
            fg=.5*(r2-r1)
            f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*an*fg
         !

         end do
         !

         do kk=1,kp
            !
            ! check derivative on sides back and front
            !
            if (myid.eq.0) then

               k=1
               r0=.5*(r(i,j,k  )+r(i+1,j,k  ))
               r1=.5*(r(i,j,k+1)+r(i+1,j,k+1))
               r2=.5*(r(i,j,k+2)+r(i+1,j,k+2))
               fg=.5*(-3.*r0+4.*r1-r2)
               an=.5*(annit(i,j,k)+annit(i+1,j,k))
               f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*an*fg

            endif
            !
            if (myid.eq.nproc-1) then

               k=jz
               r0=.5*(r(i,j,jz  )+r(i+1,j,jz  ))
               r1=.5*(r(i,j,jz-1)+r(i+1,j,jz-1))
               r2=.5*(r(i,j,jz-2)+r(i+1,j,jz-2))
               fg=.5*(3.*r0-4.*r1+r2)
               an=.5*(annit(i,j,k)+annit(i+1,j,k))
               f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*an*fg

            endif
         !
         end do

      end do

   end do


   do k=kparasta,kparaend
      do j=1,jy
         do i=ip,jx-ip
            !

            if(body.eq.1)then
               if(tipo(i,j,k).eq.0)f1(i,j,k)=0.
               if(i.lt.jx)then
                  if(tipo(i+1,j,k).eq.0)f1(i,j,k)=0.
               endif
            endif
         enddo
      enddo
   enddo




   !
   !
   ! integral of convective plus diffusive off diagonal
   !
   bulk_loc=0.
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,jx
            bulk_loc=bulk_loc+f1(i,j,k)-f1(i-1,j,k)
         end do
      end do
   end do
   !
   ! make the value known to all procs to compute the total
   !
   call MPI_ALLREDUCE(bulk_loc,bulk,1,MPI_DOUBLE_PRECISION, &
      MPI_SUM, &
      MPI_COMM_WORLD,ierr)

   ! now bulk is known to all procs
   !
   return
end
