!***********************************************************************
subroutine flud1(rc,cgra1,r,akapt,insc,isc,tipo,body)
   !***********************************************************************
   ! compute explicit flux on csi component for scalar equation
   !
   ! convective uc*(rho) with centered scheme or quick
   ! diffusive akapt*g12*d(rho)/d(eta)
   ! diffusive akapt*g13*d(rho)/d(zita)
   !
   use mysending
   use myarrays_metri3
   !
   use scala3
   use period
   !
   use mpi

   implicit none
   !-----------------------------------------------------------------------
   !     array declaration
   integer ierr,insc,isc
   integer ncolperproc,m
   integer kparastak,kparaendk
   integer i,j,k,is,iss,jj,kk,ii,iii,jjj,kkk
   !
   real ak,r0,r1,r2,fg
   real rc(0:n1,n2,kparasta:kparaend) !n3)
   real  r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
   real akapt(nscal,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
   real cgra1(0:n1,n2,kparasta-1:kparaend+1)
   real ravanti,rindietro1,rindietro2
   integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
   !      integer tipo(0:n1+1,0:n2+1,0:n3+1) c7

   real phic(n1,n2,kparasta:kparaend)
   real den, num, myv
   integer body
   logical found_solid
   !-----------------------------------------------------------------------
   ! periodicity not in eta
   !
   !-----------------------------------------------------------------------
   ! CONVECTIVE TERM Uc*rho:
   ! implemented with central scheme or with quick depending on settings
   ! in points closer to wall it uses uf as ghost cell

   !     sides left and right
   do i=1,ip

      do k=kparasta,kparaend
         do j=1,jy
            !
            f1(0,j,k) = rc(0 ,j,k)*r(0   ,j,k)-cgra1(0 ,j,k)
            f1(jx,j,k)= rc(jx,j,k)*r(jx+1,j,k)-cgra1(jx,j,k)
         !
         enddo
      enddo

   enddo
   !
   ! into the field
   !
   do k=kparasta,kparaend
      do j=1,jy
         do i=ip,jx-ip
            !
            f1(i,j,k)=rc(i,j,k)*.5*(r(i,j,k)+r(i+1,j,k))-cgra1(i,j,k)
         !
         end do
      end do
   end do
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

               if(tipo(i,j,k)==1)then !ib cell
                  if(tipo(i+1,j,k)==0)then !solido davanti
                     f1(i,j,k)=0.
                  elseif(tipo(i-1,j,k)==0)then !uso diff centrate
                     f1(i,j,k)=f1(i,j,k)
                  else ! solido sta nell'altra direzione quick
                     f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                        *.125*(-ravanti +2.*rindietro1 - rindietro2)
                  endif !solido avanti/indietro/lato
               else !i è fluido allora quick normale
                  f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                     *.125*(-ravanti +2.*rindietro1 - rindietro2)
               endif !tipo
		 
            else
               ravanti    =    r(1,j,k)
               rindietro1 =    r(2,j,k)
               rindietro2 =    r(3,j,k)

               if(tipo(i+1,j,k)==1)then
                  if(tipo(i,j,k)==0)then !solido indietro
                     f1(i,j,k)=0.
                  elseif(tipo(i+2,j,k)==0)then !uso diff centrate
                     f1(i,j,k)=f1(i,j,k)
                  else !solido di lato
                     f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                        *.125*( -ravanti +2.*rindietro1 - rindietro2)
                  endif !solido avanti/indietro/lato
               else !i+1 è fluido allora quick normale
                  f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                     *.125*(-ravanti +2.*rindietro1 - rindietro2)
               endif !tipo
            end if

            !         if(tipo(i,j,k)==2)then
            !	 f1(i,j,k)=f1(i,j,k)+rc(i,j,k)
            !     >	          *.125*( -ravanti +2.*rindietro1 - rindietro2 )
            !         end if
     
            i=jx-1
            if(rc(i,j,k).gt.0.)then
               ravanti    =    r(jx  ,j,k)
               rindietro1 =    r(jx-1,j,k)
               rindietro2 =    r(jx-2,j,k)

               if(tipo(i,j,k)==1)then !ib cell
                  if(tipo(i+1,j,k)==0)then !solido davanti
                     f1(i,j,k)=0.
                  elseif(tipo(i-1,j,k)==0)then !uso diff centrate
                     f1(i,j,k)=f1(i,j,k)
                  else ! solido sta nell'altra direzione quick
                     f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                        *.125*(-ravanti +2.*rindietro1 - rindietro2)
                  endif !solido avanti/indietro/lato
               else !i è fluido allora quick normale
                  f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                     *.125*(-ravanti +2.*rindietro1 - rindietro2)
               endif !tipo
            else
               ravanti    =        r(jx-1,j,k)
               rindietro1 =        r(jx  ,j,k)
               rindietro2 = ip*(2.*r(jx+1,j,k) - r(jx,j,k)) &
                  +(1.-ip)*r(jx+1,j,k)

               if(tipo(i+1,j,k)==1)then
                  if(tipo(i,j,k)==0)then !solido indietro
                     f1(i,j,k)=0.
                  elseif(tipo(i+2,j,k)==0)then !uso diff centrate
                     f1(i,j,k)=f1(i,j,k)
                  else !solido di lato
                     f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                        *.125*( -ravanti +2.*rindietro1 - rindietro2)
                  endif !solido avanti/indietro/lato
               else !i+1 è fluido allora quick normale
                  f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                     *.125*(-ravanti +2.*rindietro1 - rindietro2)
               endif !tipo
            end if

         !         if(tipo(i,j,k)==2)then
         !	 f1(i,j,k)=f1(i,j,k)+rc(i,j,k)
         !     >	          *.125*( -ravanti +2.*rindietro1 - rindietro2 )
         !         end if
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


   ! QUICK MODIFICATO /SMART
   elseif(insc.eq.2)then

      ! GIULIA UCS - CGRA
      !     side left and right
      do i=1,ip
         !
         do k=kparasta,kparaend
            do j=1,jy
               !
               f1(0,j,k) = rc(0 ,j,k)*r(0   ,j,k)-cgra1(0 ,j,k)
               f1(jx,j,k)= rc(jx,j,k)*r(jx+1,j,k)-cgra1(jx,j,k)
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
               f1(i,j,k)=-cgra1(i,j,k)
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
            !	den=(ravanti-rindietro2)
            !	num=(rindietro1-rindietro2)
            !if (den.lt.0.0000000001.and.den.ge.0) den=0.0000000001
            !if (den.gt.-0.0000000001.and.den.lt.0) den=-0.000000001
            !        if(den.eq.0)then
            !        phic(i,j,k)= 1.
            !	else
            !	phic(i,j,k)= (num /den)
            !        end if
            ! 	!if(phic(i,j,k).gt.1.or.phic(i,j,k).lt.0)phic(i,j,k)=-0.1
            !
            !	 if(phic(i,j,k).ge.0.17.and.phic(i,j,k).le.0.83)then !1/6<phic<4/5
            !	    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)
            !     >	          *.125*( 3*ravanti +6*rindietro1 - rindietro2 )
            !	else if(phic(i,j,k).gt.0.and.phic(i,j,k).lt.0.17)then !0<phic<1/6
            !	     f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*(3*rindietro1-2*rindietro2)
            !	else if(phic(i,j,k).gt.0.83.and.phic(i,j,k).lt.1)then !4/5<phic<1
            !	     f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*.125*(7*ravanti+rindietro1)
            !	else
            !	     f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*rindietro1
            !	endif!phic
            f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
               *.125*(3*ravanti+6*rindietro1-rindietro2)
	     
	     
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

            !	den=(ravanti-rindietro2)
            !	num=(rindietro1-rindietro2)
            !	!if (den.lt.0.0000000001.and.den.ge.0) den=0.0000000001
            !	!if (den.gt.-0.0000000001.and.den.lt.0) den=-0.000000001
            !        if(den.eq.0)then
            !        phic(i,j,k)= 1.
            !	else
            !	phic(i,j,k)= (num /den)
            !        end if
            ! 	!if(phic(i,j,k).gt.1.or.phic(i,j,k).lt.0)phic(i,j,k)=-0.1

            !	if(phic(i,j,k).ge.0.17.and.phic(i,j,k).le.0.83)then
            !	    f1(i,j,k)=f1(i,j,k)+rc(i,j,k)
            !     >	          *.125*( 3.*ravanti +6.*rindietro1 - rindietro2 )
            !	else if(phic(i,j,k).gt.0.and.phic(i,j,k).lt.0.17)then
            !	     f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*(3*rindietro1-2*rindietro2)
            !	else if(phic(i,j,k).gt.0.83.and.phic(i,j,k).lt.1)then
            !	     f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*.125*(7*ravanti+rindietro1)
            !	else
            !	     f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*rindietro1
            !	endif!phic
	
            f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
               *.125*(3*ravanti+6*rindietro1-rindietro2)
         end do
      end do
    
      
      !     into the field without bodyforce
      if(body.eq.0)then
         do k=kparasta,kparaend
            do j=1,jy
               do i=2,jx-2
                  if(rc(i,j,k).gt.0.)then
 
                     f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                        *.125*(3*r(i+1,j,k)+6*r(i,j,k)-r(i-1,j,k))

                  else  !(rc(i,j,k).le.0.)

                     f1(i,j,k)=f1(i,j,k)+rc(i,j,k)* &
                        .125*(3*r(i,j,k)+6*r(i+1,j,k)-r(i+2,j,k))
                  endif
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
                              *.125*(3*r(i+1,j,k)+6*r(i,j,k)-r(i-1,j,k))
     
                        endif !solido avanti/indietro/lato


                     else !i non è IB procedo come se non ci fossero gli ib

                        ! se i-1 è ib faccio la correzione per evitare la scia
                        !	      if(tipo(i-1,j,k)==1)then !la cella prima è un ib
                        !		if(tipo(i-2,j,k)==0)then !perché considero anche questa?
                        !   		   rindietro2 =r(i,j,k)

                        !		endif
                        !	      endif

                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k) &
                           *.125*(3*r(i+1,j,k)+6*r(i,j,k)-r(i-1,j,k))
      
                     endif !tipo

                  else !rc(i,j,k).le.0.
   
                     if(tipo(i+1,j,k)==1)then
                        if(tipo(i,j,k)==0)then !solido indietro
                           f1(i,j,k)=0.
                        elseif(tipo(i+2,j,k)==0)then !solido avanti
                           f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i+1,j,k)
                        else !solido di lato
                           f1(i,j,k)=f1(i,j,k)+rc(i,j,k)* &
                              .125*(3*r(i,j,k)+6*r(i+1,j,k)-r(i+2,j,k))
     
                        endif !solido avanti/indietro/lato

                     else !i non è IB procedo come se non ci fossero gli ib



                        ! se i+1 è ib faccio la correzione per evitare la scia
                        !	        if(tipo(i+1,j,k)==1)then
                        !		if(tipo(i+2,j,k)==0)then
                        !	          rindietro2=r(i+1,j,k)
                        !		endif
                        !		endif

                        f1(i,j,k)=f1(i,j,k)+rc(i,j,k)* &
                           .125*(3*r(i,j,k)+6*r(i+1,j,k)-r(i+2,j,k))
                     endif !tipo
	
                  endif !rc
      
               end do
            end do
         end do

      endif !body yes/no
   !	endif !SMART /quick modificato upwind

   elseif(insc.eq.3)then !upwind

      ! GIULIA UCS - CGRA
      !     side left and right
      do i=1,ip
         !
         do k=kparasta,kparaend
            do j=1,jy
               !
               f1(0,j,k) = rc(0 ,j,k)*r(0   ,j,k)-cgra1(0 ,j,k)
               f1(jx,j,k)= rc(jx,j,k)*r(jx+1,j,k)-cgra1(jx,j,k)
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
               f1(i,j,k)=-cgra1(i,j,k)
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
                  endif
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
    
                        !		elseif(tipo(i-1,j,k)==0)then
                        !
                        !	            f1(i,j,k)=f1(i,j,k)+rc(i,j,k)*r(i,j,k)
                        else ! solido sta nell'altra direzione o dietro

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
                           !               else !solido di lato
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
   endif !quick





   !      end if
   !-----------------------------------------------------------------------
   ! DIFFUSIVE TERM NNI*G12*Drho/D(ETA)
   !
   !     sides left and right
   is=0
   iss=0
   do i=1,2*ip
      do k=kparasta,kparaend
         !
         do j=1+jp,jy-jp
            !
            fg=.5*(r(iss,j+1,k)-r(iss,j-1,k))
            !             f1(is,j,k)=-f1(is,j,k)
            !     >	                  +  akapt(isc,iss,j,k)*g12(is,j,k)*fg
            !
            if(tipo(iss,j-1,k).eq.0.or.tipo(iss,j+1,k).eq.0)then
               f1(is,j,k)=-f1(is,j,k)
            !         elseif(tipo(iss,j+1,k).eq.0)then
            !           f1(is,j,k)=-f1(is,j,k)
            else
               f1(is,j,k)=-f1(is,j,k) &
                  +  akapt(isc,iss,j,k)*g12(is,j,k)*fg

            end if

         end do

         !
         do jj=1,jp
            !
            !        check derivative at wall bottom and upper
            !
            j=1
            fg=.5*(-3.*r(iss,j,k)+4.*r(iss,j+1,k)-r(iss,j+2,k))
            !         f1(is,j,k)=-f1(is,j,k)
            !     >	               +akapt(isc,iss,j,k)*g12(is,j,k)*fg
            !
            if(tipo(iss,j+1,k).eq.0.or.tipo(iss,j+2,k).eq.0)then
               f1(is,j,k)=-f1(is,j,k)
            !         elseif(tipo(iss,j+2,k).eq.0)then
            !           f1(is,j,k)=-f1(is,j,k)
            else
               f1(is,j,k)=-f1(is,j,k) &
                  +akapt(isc,iss,j,k)*g12(is,j,k)*fg
            end if


            j=jy
            fg=.5*(3.*r(iss,jy,k)-4.*r(iss,jy-1,k)+r(iss,jy-2,k))
            !         f1(is,j,k)=-f1(is,j,k)
            !     >	               +akapt(isc,iss,j,k)*g12(is,j,k)*fg

            if(tipo(iss,j-1,k).eq.0.or.tipo(iss,j-2,k).eq.0)then
               f1(is,j,k)=-f1(is,j,k)
            !         elseif(tipo(iss,j-2,k).eq.0)then
            !           f1(is,j,k)=-f1(is,j,k)
            else
               f1(is,j,k)=-f1(is,j,k) &
                  +akapt(isc,iss,j,k)*g12(is,j,k)*fg
            end if
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
   ! inside the field
   !
   do k=kparasta,kparaend
      do i=ip,jx-ip
         !
         do j=1+jp,jy-jp
            !
            r1=.5*(r(i,j-1,k)+r(i+1,j-1,k))
            r2=.5*(r(i,j+1,k)+r(i+1,j+1,k))
            fg=.5*(r2-r1)
            ak=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))
 
	 
            found_solid=.false.
            do kk=k,k
               do jj=j-1,j+1
                  do ii=i,i+1
                     if(tipo(ii,jj,kk)==0)found_solid=.true.
                  end do
               end do
            end do
	 
            if(found_solid)then
               f1(i,j,k)=-f1(i,j,k)
            else
               f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*ak*fg
            end if
 
         !
         end do
         !
         !     check derivative on sides bottom and upper
         !
         do jjj=1,jp
            !
            j=1
            r0=.5*(r(i,j  ,k)+r(i+1,j  ,k))
            r1=.5*(r(i,j+1,k)+r(i+1,j+1,k))
            r2=.5*(r(i,j+2,k)+r(i+1,j+2,k))
            fg=.5*(-3.*r0+4.*r1-r2)
            ak=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))
            !      f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*ak*fg

            found_solid=.false.
            do kk=k,k
               do jj=j,j+2
                  do ii=i,i+1
                     if(tipo(ii,jj,kk)==0)found_solid=.true.
                  end do
               end do
            end do

            if(found_solid)then
               f1(i,j,k)=-f1(i,j,k)
            else
               f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*ak*fg
            end if

            !
            j=jy
            r0=.5*(r(i,jy  ,k)+r(i+1,jy  ,k))
            r1=.5*(r(i,jy-1,k)+r(i+1,jy-1,k))
            r2=.5*(r(i,jy-2,k)+r(i+1,jy-2,k))
            fg=.5*(3.*r0-4.*r1+r2)
            ak=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))
            !      f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*ak*fg
            !
            found_solid=.false.
            do kk=k,k
               do jj=j-2,j
                  do ii=i,i+1
                     if(tipo(ii,jj,kk)==0)found_solid=.true.
                  end do
               end do
            end do

            if(found_solid)then
               f1(i,j,k)=-f1(i,j,k)
            else
               f1(i,j,k)=-f1(i,j,k)+g12(i,j,k)*ak*fg
            end if
         end do
      !
      enddo
   enddo
   !
   !-----------------------------------------------------------------------
   ! DIFFUSIVE TERM NNI*G13*Drho/D(ZITA)
   !
   ! sides left and right
   !
   ! define the limit dependeing on periodicity in z
   !
   if (myid.eq.0) then
      kparastak=kparasta+kp
      kparaendk=kparaend
   else if (myid.eq.nproc-1) then
      kparastak=kparasta
      kparaendk=kparaend-kp
   else if ((myid.ne.0).and.(myid.ne.nproc-1)) then
      kparastak=kparasta
      kparaendk=kparaend
   endif
   !
   is=0
   iss=0
   do i=1,2*ip
      !
      do j=1,jy
         !
         do k=kparastak,kparaendk
            !
            fg=.5*(r(iss,j,k+1)-r(iss,j,k-1))
            !       f1(is,j,k)=f1(is,j,k)+akapt(isc,iss,j,k)*g13(is,j,k)*fg
            !

            if(tipo(iss,j,k-1).eq.0.or.tipo(iss,j,k+1).eq.0)then
               f1(is,j,k)=f1(is,j,k)
            !         elseif(tipo(iss,j,k+1).eq.0)then
            !           f1(is,j,k)=f1(is,j,k)
            else
               f1(is,j,k)=f1(is,j,k)+akapt(isc,iss,j,k)*g13(is,j,k)*fg
            end if
         end do
         !
         do kk=1,kp
            !
            ! check derivative on sides back and front
            !
            if (myid.eq.0) then

               k=1
               fg=.5*(-3.*r(iss,j,k)+4.*r(iss,j,k+1)-r(iss,j,k+2))
               !      f1(is,j,k)=f1(is,j,k)+akapt(isc,iss,j,k)*g13(is,j,k)*fg

               if(tipo(iss,j,k+1).eq.0.or.tipo(iss,j,k+2).eq.0)then
                  f1(is,j,k)=f1(is,j,k)
               !         elseif(tipo(iss,j,k+2).eq.0)then
               !           f1(is,j,k)=f1(is,j,k)
               else
                  f1(is,j,k)=f1(is,j,k)+akapt(isc,iss,j,k)*g13(is,j,k)*fg
               end if

            endif
            !
            if (myid.eq.nproc-1) then
               k=jz
               fg=.5*(3.*r(iss,j,jz)-4.*r(iss,j,jz-1)+r(iss,j,jz-2))
               !      f1(is,j,k)=f1(is,j,k)+akapt(isc,iss,j,k)*g13(is,j,k)*fg

               if(tipo(iss,j,k-1).eq.0.or.tipo(iss,j,k-2).eq.0)then
                  f1(is,j,k)=f1(is,j,k)
               !         elseif(tipo(iss,j,k-2).eq.0)then
               !           f1(is,j,k)=f1(is,j,k)
               else
                  f1(is,j,k)=f1(is,j,k)+akapt(isc,iss,j,k)*g13(is,j,k)*fg
               end if
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
         !
         do k=kparastak,kparaendk
            !
            ak=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))
            r1=.5*(r(i,j,k-1)+r(i+1,j,k-1))
            r2=.5*(r(i,j,k+1)+r(i+1,j,k+1))
            fg=.5*(r2-r1)

            !
            found_solid=.false.
            do kk=k-1,k+1
               do jj=j,j
                  do ii=i,i+1
                     if(tipo(ii,jj,kk)==0)found_solid=.true.
                  end do
               end do
            end do
	 
            if(found_solid)then
               f1(i,j,k)=f1(i,j,k)
            else
               f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*ak*fg
            end if




         end do
         !
         do kkk=1,kp
            !
            ! check derivative on sides back and front
            !
            if (myid.eq.0) then

               k=1
               r0=.5*(r(i,j,k  )+r(i+1,j,k  ))
               r1=.5*(r(i,j,k+1)+r(i+1,j,k+1))
               r2=.5*(r(i,j,k+2)+r(i+1,j,k+2))
               fg=.5*(-3.*r0+4.*r1-r2)
               ak=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))
               !      f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*ak*fg


               found_solid=.false.
               do kk=k,k+2
                  do jj=j,j
                     do ii=i,i+1
                        if(tipo(ii,jj,kk)==0)found_solid=.true.
                     end do
                  end do
               end do

               if(found_solid)then
                  f1(i,j,k)=f1(i,j,k)
               else
                  f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*ak*fg
               end if



            endif
            !
            if (myid.eq.nproc-1) then

               k=jz
               r0=.5*(r(i,j,jz  )+r(i+1,j,jz  ))
               r1=.5*(r(i,j,jz-1)+r(i+1,j,jz-1))
               r2=.5*(r(i,j,jz-2)+r(i+1,j,jz-2))
               fg=.5*(3.*r0-4.*r1+r2)
               ak=.5*(akapt(isc,i,j,k)+akapt(isc,i+1,j,k))
               !      f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*ak*fg

               found_solid=.false.
               do kk=k-2,k
                  do jj=j,j
                     do ii=i,i+1
                        if(tipo(ii,jj,kk)==0)found_solid=.true.
                     end do
                  end do
               end do

               if(found_solid)then
                  f1(i,j,k)=f1(i,j,k)
               else
                  f1(i,j,k)=f1(i,j,k)+g13(i,j,k)*ak*fg
               end if

            endif
         !
         end do
      !
      enddo
   enddo
	

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
   return
end
