!***********************************************************************
subroutine flux2(rc,cgra2,r,insc,tipo,body)
   !***********************************************************************
   ! compute explicit flux on eta component
   !
   ! convective  vc*(u,v,w)   with central schema or quick
   ! diffusive   nni*g21*d(u,v,w)/d(csi)
   ! diffusive   nni*g23*d(u,v,w)/d(zita)
   !
   use mysending
   use myarrays_wallmodel
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
   double precision bulk_loc,bulk2

   integer i,j,k,ii,kk,js,jss
   real xf,yf,zf,xc1,yc1,zc1,xc2,yc2,zc2,ddj,dsj
   real fg,r0,r1,r2,an
   real rc(n1,0:n2,kparasta:kparaend)
   real r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)    !0:n3+1)
   real cgra2(n1,0:n2,kparasta-1:kparaend+1)
   real ravanti,rindietro1,rindietro2
   integer iwall,iwfp
   integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
   !
   real phic(n1,n2,kparasta:kparaend)
   real den, num, myv
   integer body
   !-----------------------------------------------------------------------
   ! CONVECTIVE TERM Vc*u
   ! implemented with central scheme or with quick depending on settings
   ! in points closer to wall it uses uf as ghost cell

   !     sides bottom and upper
   do j=1,jp
      !
      do k=kparasta,kparaend
         do i=1,jx
            !
            f2(i,0,k)  = rc(i,0 ,k)*r(i,0   ,k)-cgra2(i,0 ,k)
            f2(i,jy,k) = rc(i,jy,k)*r(i,jy+1,k)-cgra2(i,jy,k)
         !
         enddo
      enddo

   enddo

   !     into the field
   do k=kparasta,kparaend
      do j=jp,jy-jp
         do i=1,jx
            !
            f2(i,j,k)=rc(i,j,k)*.5*(r(i,j,k)+r(i,j+1,k))-cgra2(i,j,k)
         !
         end do
      end do
   end do
   !
   ! quick
   !
   if(insc.eq.1)then
      
      !     sides 3 and 4
      do k=kparasta,kparaend
         do i=1,jx
          
            j=1
            if(rc(i,j,k).gt.0.)then
               ravanti    =        r(i,2,k)
               rindietro1 =        r(i,1,k)
               rindietro2 = jp*(2.*r(i,0,k) - r(i,1,k)) &
                  +(1.-jp)*r(i,0,k)
            else
               ravanti    = r(i,1,k)
               rindietro1 = r(i,2,k)
               rindietro2 = r(i,3,k)
            end if

            if(tipo(i,j,k)==2)then
               f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                  *.125*( -ravanti +2.*rindietro1 - rindietro2 )
            end if
     
            j=jy-1
            if(rc(i,j,k).gt.0.)then
               ravanti    = r(i,jy  ,k)
               rindietro1 = r(i,jy-1,k)
               rindietro2 = r(i,jy-2,k)
            else
               ravanti    =        r(i,jy-1,k)
               rindietro1 =        r(i,jy  ,k)
               rindietro2 = jp*(2.*r(i,jy+1,k) - r(i,jy,k)) &
                  +(1.-jp)*r(i,jy+1,k)
            end if

            if(tipo(i,j,k)==2)then
               f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                  *.125*( -ravanti +2.*rindietro1 - rindietro2 )
            end if
         end do
      end do         
          
      !     into the field
      do k=kparasta,kparaend
         do j=2,jy-2
            do i=1,jx
               !
               !         if(tipo(i,j,k)==2)then
               !      if (rc(i,j,k).gt.0.) then
               !      f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*
               !     > .125*(-r(i,j+1,k)+2*r(i,j,k)-r(i,j-1,k))
               !      else
               !      f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*
               !     > .125*(-r(i,j,k)+2*r(i,j+1,k)-r(i,j+2,k))
               !      end if
               !      end if
               !

               if (rc(i,j,k).gt.0.) then

                  if(tipo(i,j,k)==1)then !ib cell
                     if(tipo(i,j+1,k)==0)then !solido davanti
                        f2(i,j,k)=0.
                     elseif(tipo(i,j-1,k)==0)then !uso diff centrate
                        f2(i,j,k)=f2(i,j,k)
                     else ! solido sta nell'altra direzione quick
                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                           .125*(-r(i,j+1,k)+2*r(i,j,k)-r(i,j-1,k))

                     endif !solido avanti/indietro/lato
                  else !i è fluido allora quick normale
                     f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                        .125*(-r(i,j+1,k)+2*r(i,j,k)-r(i,j-1,k))
                  endif !tipo
    
               else !rc(i,j,k).le.0.
                  if(tipo(i,j+1,k)==1)then
                     if(tipo(i,j,k)==0)then !solido indietro
                        f2(i,j,k)=0.
                     elseif(tipo(i,j+2,k)==0)then !uso diff centrate
                        f2(i,j,k)=f2(i,j,k)

                     else !solido di lato
                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                           .125*(-r(i,j,k)+2*r(i,j+1,k)-r(i,j+2,k))


                     endif !solido avanti/indietro/lato
                  else !i+1 è fluido allora quick normale
                     f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                        .125*(-r(i,j,k)+2*r(i,j+1,k)-r(i,j+2,k))
                  endif !tipo
	
               endif !rc
            end do
         end do
      end do

   ! quick upwind
   elseif(insc.eq.2)then

      ! GIULIA UCS - CGRA
      !     side left and right
      do j=1,jp
         !
         do k=kparasta,kparaend
            do i=1,jx
               !
               f2(i,0,k)  = rc(i,0 ,k)*r(i,0   ,k)-cgra2(i,0 ,k)
               f2(i,jy,k) = rc(i,jy,k)*r(i,jy+1,k)-cgra2(i,jy,k)
            !
            enddo
         enddo

      enddo
      !
      !     into the field
      do k=kparasta,kparaend
         do j=jp,jy-jp
            do i=1,jx
               !
               f2(i,j,k)=-cgra2(i,j,k)
            !
            end do
         end do
      end do


      !     sides 1 and 2  sia per bodyforce 0 o 1
      do k=kparasta,kparaend
         do i=1,jx
            j=1
            if(rc(i,j,k).gt.0.)then
               ravanti    =        r(i,2,k)
               rindietro1 =        r(i,1,k)
               rindietro2 = jp*(2.*r(i,0,k) - r(i,1,k)) &
                  +(1.-jp)*r(i,0,k)
            else
               ravanti    =    r(i,1,k)
               rindietro1 =    r(i,2,k)
               rindietro2 =    r(i,3,k)
            end if

            f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
               *.125*( 3.*ravanti +6.*rindietro1 - rindietro2 )

            j=jy-1
            if(rc(i,j,k).gt.0.)then
               ravanti    =    r(i,jy  , k)
               rindietro1 =    r(i,jy-1, k)
               rindietro2 =    r(i,jy-2, k)
            else
               ravanti    =        r(i,jy-1,k)
               rindietro1 =        r(i,jy  ,k)
               rindietro2 = jp*(2.*r(i,jy+1,k) - r(i,jy,k)) &
                  +(1.-jp)*r(i,jy+1,k)
            end if


            f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
               *.125*( 3.*ravanti +6.*rindietro1 - rindietro2 )

         end do
      end do
    
      
      !     into the field without bodyforce
      if(body.eq.0)then
         do k=kparasta,kparaend
            do j=2,jy-2
               do i=1,jx
                  if(rc(i,j,k).gt.0.)then

                     f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                        *.125*(3.*r(i,j+1,k)+6*r(i,j,k)-r(i,j-1,k))


                  else  !(rc(i,j,k).le.0.)

                     f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                        *.125*(3.*r(i,j,k)+6*r(i,j+1,k)-r(i,j+2,k)) !my index +1


                  endif !rc</>0
               end do
            end do
         end do

      !
      else !bodyforce attivo

 
         do k=kparasta,kparaend
            do j=2,jy-2
               do i=1,jx
                  if(rc(i,j,k).gt.0.)then
                     if(tipo(i,j,k)==1)then !ib cell
                        if(tipo(i,j+1,k)==0)then !solido davanti

                           f2(i,j,k)=0.
                        ! dopo diventerà nullo
                        elseif(tipo(i,j-1,k)==0)then

                           f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)
                        else ! solido sta nell'altra direzione

                           f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                              *.125*(3.*r(i,j+1,k)+6*r(i,j,k)-r(i,j-1,k))

                        endif !solido avanti/indietro/lato


                     else !i non è IB procedo come se non ci fossero gli ib

                        ! se i-1 è ib faccio la correzione per evitare la scia
                        !	      if(tipo(i-1,j,k)==1)then !la cella prima è un ib
                        !		if(tipo(i-2,j,k)==0)then !perché considero anche questa?
                        !   		   rindietro2 =r(i,j,k)

                        !		endif
                        !	      endif


                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                           *.125*(3.*r(i,j+1,k)+6*r(i,j,k)-r(i,j-1,k))


                     endif !tipo

                  else !rc(i,j,k).le.0.
   
                     if(tipo(i,j+1,k)==1)then
                        if(tipo(i,j,k)==0)then !solido indietro
                           f2(i,j,k)=0.
                        elseif(tipo(i,j+2,k)==0)then !solido avanti
                           f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)
                        else !solido di lato

                           f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                              *.125*(3.*r(i,j,k)+6*r(i,j+1,k)-r(i,j+2,k))

                        endif !solido avanti/indietro/lato

                     else !i non è IB procedo come se non ci fossero gli ib

                        ! se i+1 è ib faccio la correzione per evitare la scia
                        !	        if(tipo(i+1,j,k)==1)then
                        !		if(tipo(i+2,j,k)==0)then
                        !	          rindietro2=r(i+1,j,k)
                        !		endif
                        !		endif


                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                           *.125*(3.*r(i,j,k)+6*r(i,j+1,k)-r(i,j+2,k))

                     endif !tipo
	
                  endif !rc
      
               end do
            end do
         end do

      endif !body yes/no
   !	endif !SMART/	quick modificato upwind

   ! quick upwind
   elseif(insc.eq.3)then

      ! GIULIA UCS - CGRA
      !     side left and right
      do j=1,jp
         !
         do k=kparasta,kparaend
            do i=1,jx
               !
               f2(i,0,k)  = rc(i,0 ,k)*r(i,0   ,k)-cgra2(i,0 ,k)
               f2(i,jy,k) = rc(i,jy,k)*r(i,jy+1,k)-cgra2(i,jy,k)
            !
            enddo
         enddo

      enddo
      !
      !     into the field
      do k=kparasta,kparaend
         do j=jp,jy-jp
            do i=1,jx
               !
               f2(i,j,k)=-cgra2(i,j,k)
            !
            end do
         end do
      end do


      !     sides 1 and 2  sia per bodyforce 0 o 1
      do k=kparasta,kparaend
         do i=1,jx
            j=1
            if(rc(i,j,k).gt.0.)then
 	 
               rindietro1 =        r(i,1,k)
 		 
            else
 
               rindietro1 =    r(i,2,k)
	 	 
            end if

            f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*rindietro1

            j=jy-1
            if(rc(i,j,k).gt.0.)then
 
               rindietro1 =    r(i,jy-1, k)

            else

               rindietro1 =        r(i,jy  ,k)
	 	 	 
            end if


            f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*rindietro1

         end do
      end do
    
      
      !     into the field without bodyforce
      if(body.eq.0)then
         do k=kparasta,kparaend
            do j=2,jy-2
               do i=1,jx
                  if(rc(i,j,k).gt.0.)then

                     f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)


                  else  !(rc(i,j,k).le.0.)

                     f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)


                  endif !rc</>0
               end do
            end do
         end do

      !
      else !bodyforce attivo
         do k=kparasta,kparaend
            do j=2,jy-2
               do i=1,jx
                  if(rc(i,j,k).gt.0.)then
                     if(tipo(i,j,k)==1)then !ib cell
                        if(tipo(i,j+1,k)==0)then !solido davanti

                           f2(i,j,k)=0.
                        ! dopo diventerà nullo

                        else ! solido sta nell'altra direzione

                           f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)

                        endif !solido avanti/indietro/lato


                     else !i non è IB procedo come se non ci fossero gli ib

                        ! se i-1 è ib faccio la correzione per evitare la scia
                        !	      if(tipo(i-1,j,k)==1)then !la cella prima è un ib
                        !		if(tipo(i-2,j,k)==0)then !perché considero anche questa?
                        !   		   rindietro2 =r(i,j,k)

                        !		endif
                        !	      endif

                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)

                     endif !tipo

                  else !rc(i,j,k).le.0.
   
                     if(tipo(i,j+1,k)==1)then
                        if(tipo(i,j,k)==0)then !solido indietro
                           f2(i,j,k)=0.
                        elseif(tipo(i,j+2,k)==0)then !solido avanti
                           f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)
                        else !solido di lato

                           f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)

                        endif !solido avanti/indietro/lato

                     else !i non è IB procedo come se non ci fossero gli ib

                        ! se i+1 è ib faccio la correzione per evitare la scia
                        !	        if(tipo(i+1,j,k)==1)then
                        !		if(tipo(i+2,j,k)==0)then
                        !	          rindietro2=r(i+1,j,k)
                        !		endif
                        !		endif

                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)

                     endif !tipo
	
                  endif !rc
      
               end do
            end do
         end do

      endif !body yes/no
   endif !SMART/	quick modificato upwind
   !
   !-----------------------------------------------------------------------
   ! DIFFUSIVE TERM NNI*G21*DU/D(CSI)
   !
   !     sides bottom and upper
   !
   js=0
   jss=0
   iwfp = wfp3
      
   do j=1,2*jp
      !
      do k=kparasta,kparaend
         !


         do  i=1+ip,jx-ip
            fg=.5*(r(i+1,jss,k)-r(i-1,jss,k))
            !         WALL MODEL : all in flucn
            do iwall = 1,iwfp
               f2(i,js,k)=-f2(i,js,k)
            end do
            !         NO WALL MODEL
            do iwall = 1,1-iwfp
               f2(i,js,k)=-f2(i,js,k) &
                  +annit(i,jss,k)*g21(i,js,k)*fg
            end do
         end do
      
         !
         do ii=1,ip
            !     check derivative on sides left and right
            i=1
            fg=.5*(-3.*r(i,jss,k)+4.*r(i+1,jss,k)-r(i+2,jss,k))
            f2(i,js,k)=-f2(i,js,k) &
               +annit(i,jss,k)*g21(i,js,k)*fg
            i=jx
            fg=.5*(3.*r(jx,jss,k)-4.*r(jx-1,jss,k)+r(jx-2,jss,k))
            f2(i,js,k)=-f2(i,js,k) &
               +annit(i,jss,k)*g21(i,js,k)*fg
         end do
      
      
      !
      enddo
      !
      js=jy
      jss=jy+1
      iwfp = wfp4
   !
   enddo
   !
   ! into the field
   !
   do k=kparasta,kparaend
      do j=jp,jy-jp
         !
         do i=1+ip,jx-ip
            !
            r1=.5*(r(i-1,j,k)+r(i-1,j+1,k))
            r2=.5*(r(i+1,j,k)+r(i+1,j+1,k))
            fg=.5*(r2-r1)
            an=.5*(annit(i,j,k)+annit(i,j+1,k))
            f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*an*fg
         !
         end do
         !
         do ii=1,ip
            !
            ! check derivative on sides left and right
            !
            i=1
            r0=.5*(r(i  ,j,k)+r(i  ,j+1,k))
            r1=.5*(r(i+1,j,k)+r(i+1,j+1,k))
            r2=.5*(r(i+2,j,k)+r(i+2,j+1,k))
            fg=.5*(-3.*r0+4.*r1-r2)
            an=.5*(annit(i,j,k)+annit(i,j+1,k))
            f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*an*fg
            !
            i=jx
            r0=.5*(r(jx  ,j,k)+r(jx  ,j+1,k))
            r1=.5*(r(jx-1,j,k)+r(jx-1,j+1,k))
            r2=.5*(r(jx-2,j,k)+r(jx-2,j+1,k))
            fg=.5*(3.*r0-4.*r1+r2)
            an=.5*(annit(i,j,k)+annit(i,j+1,k))
            f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*an*fg
         !
         end do
      !
      enddo
   enddo
   !
   !-----------------------------------------------------------------------
   ! DIFFUSIVE TERM NNI*G23*DU/D(ZITA)
   !
   ! sides bottom and upper
   !
   ! define computation limit depending on periodicity in z
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

   js=0
   jss=0
   iwfp = wfp3
   do j=1,2*jp
      !
      do i=1,jx
         !
         do k=kparastak,kparaendk
            !
            fg=.5*(r(i,jss,k+1)-r(i,jss,k-1))

            !           WALL MODEL: all in flucn
            do iwall=1,iwfp
               f2(i,js,k)=f2(i,js,k)	    
            end do
            !           NO WALL MODEL
            do iwall=1,1-iwfp	    
               f2(i,js,k)=f2(i,js,k) &
                  +annit(i,jss,k)*g23(i,js,k)*fg
            end do
         !
         end do
         !

         do kk=1,kp
            !
            ! check derivative on sides back and front
            !
            if (myid.eq.0) then

               k=1
               fg=.5*(-3.*r(i,jss,k)+4.*r(i,jss,k+1)-r(i,jss,k+2))
               f2(i,js,k)=f2(i,js,k) &
                  +annit(i,jss,k)*g23(i,js,k)*fg

            endif
            !
            if (myid.eq.nproc-1) then

               k=jz
               fg=.5*(3.*r(i,jss,jz)-4.*r(i,jss,jz-1)+r(i,jss,jz-2))
               f2(i,js,k)=f2(i,js,k) &
                  +annit(i,jss,k)*g23(i,js,k)*fg

            endif
         !
         end do
      !
      enddo
      !
      js=jy
      jss=jy+1
      iwfp = wfp4
   !
   enddo
   !
   ! into the field
   !
   do j=jp,jy-jp
      do i=1,jx
         !
         do k=kparastak,kparaendk
            !
            r1=.5*(r(i,j,k-1)+r(i,j+1,k-1))
            r2=.5*(r(i,j,k+1)+r(i,j+1,k+1))
            fg=.5*(r2-r1)
            an=.5*(annit(i,j,k)+annit(i,j+1,k))
            f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*an*fg
         !
         end do
         !
         do kk=1,kp
            !
            ! check derivative on sides back and front
            !
            if (myid.eq.0) then

               k=1
               r0=.5*(r(i,j,k  )+r(i,j+1,k  ))
               r1=.5*(r(i,j,k+1)+r(i,j+1,k+1))
               r2=.5*(r(i,j,k+2)+r(i,j+1,k+2))
               fg=.5*(-3.*r0+4.*r1-r2)
               an=.5*(annit(i,j,k)+annit(i,j+1,k))
               f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*an*fg

            endif
            !
            if (myid.eq.nproc-1) then

               k=jz
               r0=.5*(r(i,j,jz  )+r(i,j+1,jz  ))
               r1=.5*(r(i,j,jz-1)+r(i,j+1,jz-1))
               r2=.5*(r(i,j,jz-2)+r(i,j+1,jz-2))
               fg=.5*(3.*r0-4.*r1+r2)
               an=.5*(annit(i,j,k)+annit(i,j+1,k))
               f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*an*fg

            endif
         !
         end do
      !
      enddo
   enddo


   do k=kparasta,kparaend
      do j=jp,jy-jp
         do i=1,jx
            !

            if(body.eq.1)then
               if(tipo(i,j,k).eq.0)f2(i,j,k)=0.
               if(j.lt.jy)then
                  if(tipo(i,j+1,k).eq.0)f2(i,j,k)=0.
               endif
            endif
         enddo
      enddo
   enddo


   !
   ! integral on f2
   bulk_loc=0.
   !
   do k=kparasta,kparaend
      do j=1,jy
         do i=1,jx
            bulk_loc=bulk_loc+f2(i,j,k)-f2(i,j-1,k)
         end do
      end do
   end do

   ! make the value known to all procs
   call MPI_ALLREDUCE(bulk_loc,bulk2,1,MPI_DOUBLE_PRECISION, &
      MPI_SUM, &
      MPI_COMM_WORLD,ierr)

   ! now bulk is known by all procs

   bulk=bulk+bulk2

   return
end
