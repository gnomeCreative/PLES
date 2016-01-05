!***********************************************************************
subroutine flud2(rc,cgra2,r,akapt,insc,isc,tipo,body)
   !***********************************************************************
   ! compute explicit flux on eta component for scalar equation
   !
   ! convective vc*(rho) with centered scheme or quick
   ! diffusive akapt*g21*d(rho)/d(csi)
   ! diffusive akapt*g23*d(rho)/d(zita)
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
   integer i,j,k,js,jss,ii,kk,jj,iii,jjj,kkk

   real ak,r0,r1,r2,fg
   real rc(n1,0:n2,kparasta:kparaend) !n3)
   real r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
   real akapt(nscal,0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
   real cgra2(n1,0:n2,kparasta-1:kparaend+1)
   real ravanti,rindietro1,rindietro2
   integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)
   !      integer tipo(0:n1+1,0:n2+1,0:n3+1) c7
   !
   real phic(n1,n2,kparasta:kparaend)
   real den, num, myv
   integer body
   logical found_solid
   !-----------------------------------------------------------------------
   ! CONVECTIVE TERM Vc*rho:
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
   !
   ! into the field
   !
   do k=kparasta,kparaend
      do i=1,jx
         do j=jp,jy-jp
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

               if(tipo(i,j,k)==1)then !ib cell
                  if(tipo(i,j+1,k)==0)then !solido davanti
                     f2(i,j,k)=0.
                  elseif(tipo(i,j-1,k)==0)then !uso diff centrate
                     f2(i,j,k)=f2(i,j,k)
                  else ! solido sta nell'altra direzione quick
                     f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                        .125*(-ravanti +2.*rindietro1 - rindietro2)

                  endif !solido avanti/indietro/lato
               else !i è fluido allora quick normale
                  f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                     .125*(-ravanti +2.*rindietro1 - rindietro2)
               endif !tipo
            else
               ravanti    = r(i,1,k)
               rindietro1 = r(i,2,k)
               rindietro2 = r(i,3,k)

               if(tipo(i,j+1,k)==1)then
                  if(tipo(i,j,k)==0)then !solido indietro
                     f2(i,j,k)=0.
                  elseif(tipo(i,j+2,k)==0)then !uso diff centrate
                     f2(i,j,k)=f2(i,j,k)

                  else !solido di lato
                     f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                        .125*(-ravanti +2.*rindietro1 - rindietro2)


                  endif !solido avanti/indietro/lato
               else !i+1 è fluido allora quick normale
                  f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                     .125*(-ravanti +2.*rindietro1 - rindietro2)
               endif !tipo
            end if

            !         if(tipo(i,j,k)==2)then
            !	 f2(i,j,k)=f2(i,j,k)+rc(i,j,k)
            !     >	          *.125*( -ravanti +2.*rindietro1 - rindietro2 )
            !         end if
     
            j=jy-1
            if(rc(i,j,k).gt.0.)then
               ravanti    = r(i,jy  ,k)
               rindietro1 = r(i,jy-1,k)
               rindietro2 = r(i,jy-2,k)

               if(tipo(i,j,k)==1)then !ib cell
                  if(tipo(i,j+1,k)==0)then !solido davanti
                     f2(i,j,k)=0.
                  elseif(tipo(i,j-1,k)==0)then !uso diff centrate
                     f2(i,j,k)=f2(i,j,k)
                  else ! solido sta nell'altra direzione quick
                     f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                        .125*(-ravanti +2.*rindietro1 - rindietro2)

                  endif !solido avanti/indietro/lato
               else !i è fluido allora quick normale
                  f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                     .125*(-ravanti +2.*rindietro1 - rindietro2)
               endif !tipo
            else
               ravanti    =        r(i,jy-1,k)
               rindietro1 =        r(i,jy  ,k)
               rindietro2 = jp*(2.*r(i,jy+1,k) - r(i,jy,k)) &
                  +(1.-jp)*r(i,jy+1,k)

               if(tipo(i,j+1,k)==1)then
                  if(tipo(i,j,k)==0)then !solido indietro
                     f2(i,j,k)=0.
                  elseif(tipo(i,j+2,k)==0)then !uso diff centrate
                     f2(i,j,k)=f2(i,j,k)

                  else !solido di lato
                     f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                        .125*(-ravanti +2.*rindietro1 - rindietro2)


                  endif !solido avanti/indietro/lato
               else !i+1 è fluido allora quick normale
                  f2(i,j,k)=f2(i,j,k)+rc(i,j,k)* &
                     .125*(-ravanti +2.*rindietro1 - rindietro2)
               endif !tipo
            end if
         !         if(tipo(i,j,k)==2)then
         !	 f2(i,j,k)=f2(i,j,k)+rc(i,j,k)
         !     >	          *.125*( -ravanti +2.*rindietro1 - rindietro2 )
         !         end if
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
            !
            end do
         end do
      end do

   ! smart/quick modificato
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

            !	 if(phic(i,j,k).ge.0.17.and.phic(i,j,k).le.0.83)then !1/6<phic<4/5
            f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
               *.125*( 3.*ravanti +6.*rindietro1 - rindietro2 )
            !	else if(phic(i,j,k).gt.0.and.phic(i,j,k).lt.0.17)then !0<phic<1/6
            !	     f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*(3*rindietro1-2*rindietro2)
            !	else if(phic(i,j,k).gt.0.83.and.phic(i,j,k).lt.1)then !4/5<phic<1
            !	     f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*.125*(7*ravanti+rindietro1)
            !	else
            !	     f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*rindietro1
            !	endif!phic

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

            !	den=(ravanti-rindietro2)
            !	num=(rindietro1-rindietro2)
 	
            !        !if (den.lt.0.0000000001.and.den.ge.0) den=0.0000000001
            !	!if (den.gt.-0.0000000001.and.den.lt.0) den=-0.000000001
            !        if(den.eq.0)then
            !        phic(i,j,k)= 1.
            !	else
            !	phic(i,j,k)= (num /den)
            !        end if
            ! 	!if(phic(i,j,k).gt.1.or.phic(i,j,k).lt.0)phic(i,j,k)=-0.1

            !	if(phic(i,j,k).ge.0.17.and.phic(i,j,k).le.0.83)then
            f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
               *.125*( 3.*ravanti +6.*rindietro1 - rindietro2 )
         !	else if(phic(i,j,k).gt.0.and.phic(i,j,k).lt.0.17)then
         !	     f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*(3*rindietro1-2*rindietro2)
         !	else if(phic(i,j,k).gt.0.83.and.phic(i,j,k).lt.1)then
         !	     f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*.125*(7*ravanti+rindietro1)
         !	else
         !	     f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*rindietro1
         !	endif!phic
         end do
      end do
    
      
      !     into the field without bodyforce
      if(body.eq.0)then
         do k=kparasta,kparaend
            do j=2,jy-2
               do i=1,jx
                  if(rc(i,j,k).gt.0.)then

                     f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                        *.125*(3*r(i,j+1,k)+6*r(i,j,k)-r(i,j-1,k))


                  else  !(rc(i,j,k).le.0.)

                     f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                        *.125*(3*r(i,j,k)+6*r(i,j+1,k)-r(i,j+2,k)) !my index +1

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

                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                           *.125*(3*r(i,j+1,k)+6*r(i,j,k)-r(i,j-1,k))

                     endif   !tipo
             
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

                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k) &
                           *.125*(3.*r(i,j,k)+6*r(i,j+1,k)-r(i,j+2,k))

                     endif !tipo
	
                  endif !rc
      
               end do
            end do
         end do

      endif !body yes/no
   !	endif !SMART/QUICK modificato upwind

   ! smart/quick modificato
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
                        !		elseif(tipo(i,j-1,k)==0)then

                        !	            f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)
                        else ! solido sta nell'altra direzione


                           f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)

                        endif !solido avanti/indietro/lato


                     else !i non è IB procedo come se non ci fossero gli ib

                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j,k)
                     endif   !tipo
             
                  else !rc(i,j,k).le.0.
   
                     if(tipo(i,j+1,k)==1)then
                        if(tipo(i,j,k)==0)then !solido indietro
                           f2(i,j,k)=0.
                        else
                           f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)
 
                        endif !solido avanti/indietro/lato

                     else !i non è IB procedo come se non ci fossero gli ib

                        f2(i,j,k)=f2(i,j,k)+rc(i,j,k)*r(i,j+1,k)
                     endif !tipo
	
                  endif !rc
      
               end do
            end do
         end do

      endif !body yes/no
   endif !quick
   !-----------------------------------------------------------------------
   ! DIFFUSIVE TERM NNI*G21*Drho/D(CSI)
   !
   !     sides bottom and upper
   js=0
   jss=0
   do j=1,2*jp
      !
      do k=kparasta,kparaend
         !
         do  i=1+ip,jx-ip
            !
            fg=.5*(r(i+1,jss,k)-r(i-1,jss,k))
            !      f2(i,js,k)=-f2(i,js,k)+akapt(isc,i,jss,k)*g21(i,js,k)*fg
            !

            if(tipo(i-1,jss,k).eq.0.or.tipo(i+1,jss,k).eq.0)then
               f2(i,js,k)=-f2(i,js,k)
            !         elseif(tipo(i+1,jss,k).eq.0)then
            !           f2(i,js,k)=-f2(i,js,k)
            else
               f2(i,js,k)=-f2(i,js,k)+akapt(isc,i,jss,k)*g21(i,js,k)*fg
            end if

         end do
         !
         do ii=1,ip
            !
            !     check derivative on sides left and right
            !
            i=1
            fg=.5*(-3.*r(i,jss,k)+4.*r(i+1,jss,k)-r(i+2,jss,k))
            !      f2(i,js,k)=-f2(i,js,k)+akapt(isc,i,jss,k)*g21(i,js,k)*fg

            if(tipo(i+1,jss,k).eq.0.or.tipo(i+2,jss,k).eq.0)then
               f2(i,js,k)=-f2(i,js,k)
            !         elseif(tipo(i+2,jss,k).eq.0)then
            !           f2(i,js,k)=-f2(i,js,k)
            else
               f2(i,js,k)=-f2(i,js,k)+akapt(isc,i,jss,k)*g21(i,js,k)*fg
            end if
            !
            i=jx
            fg=.5*(3.*r(jx,jss,k)-4.*r(jx-1,jss,k)+r(jx-2,jss,k))
            !      f2(i,js,k)=-f2(i,js,k)+akapt(isc,i,jss,k)*g21(i,js,k)*fg

            if(tipo(i-1,jss,k).eq.0.or.tipo(i-2,jss,k).eq.0)then
               f2(i,js,k)=-f2(i,js,k)
            !         elseif(tipo(i-2,jss,k).eq.0)then
            !           f2(i,js,k)=-f2(i,js,k)
            else
               f2(i,js,k)=-f2(i,js,k)+akapt(isc,i,jss,k)*g21(i,js,k)*fg
            end if
         !
         end do
      !
      enddo
      !
      js=jy
      jss=jy+1
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
            ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j+1,k))
	 
	 
            found_solid=.false.
            do kk=k,k
               do jj=j,j+1
                  do ii=i-1,i+1
                     if(tipo(ii,jj,kk)==0)found_solid=.true.
                  end do
               end do
            end do
	 
            if(found_solid)then
               f2(i,j,k)=-f2(i,j,k)
            else
               f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*ak*fg
            end if


         end do
         !
         do iii=1,ip
            !
            ! check derivative on sides left and right
            !
            i=1
            r0=.5*(r(i  ,j,k)+r(i  ,j+1,k))
            r1=.5*(r(i+1,j,k)+r(i+1,j+1,k))
            r2=.5*(r(i+2,j,k)+r(i+2,j+1,k))
            fg=.5*(-3.*r0+4.*r1-r2)
            ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j+1,k))
            !      f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*ak*fg

            found_solid=.false.
            do kk=k,k
               do jj=j,j+1
                  do ii=i,i+2
                     if(tipo(ii,jj,kk)==0)found_solid=.true.
                  end do
               end do
            end do
	 
            if(found_solid)then
               f2(i,j,k)=-f2(i,j,k)
            else
               f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*ak*fg
            end if
            !
            i=jx
            r0=.5*(r(jx  ,j,k)+r(jx  ,j+1,k))
            r1=.5*(r(jx-1,j,k)+r(jx-1,j+1,k))
            r2=.5*(r(jx-2,j,k)+r(jx-2,j+1,k))
            fg=.5*(3.*r0-4.*r1+r2)
            ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j+1,k))
            !      f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*ak*fg
            !
            found_solid=.false.
            do kk=k,k
               do jj=j,j+1
                  do ii=i-2,i
                     if(tipo(ii,jj,kk)==0)found_solid=.true.
                  end do
               end do
            end do
	 
            if(found_solid)then
               f2(i,j,k)=-f2(i,j,k)
            else
               f2(i,j,k)=-f2(i,j,k)+g21(i,j,k)*ak*fg
            end if

         end do
      !
      enddo
   enddo
   !
   !-----------------------------------------------------------------------
   ! DIFFUSIVE TERM NNI*G23*Drho/D(ZITA)
   !
   ! sides bottom and upper
   !
   if (myid.eq.0) then
      kparastak=kparasta+kp
      kparaendk=kparaend
   else if (myid.eq.nproc-1) then
      kparastak=kparasta
      kparaendk=kparaend-kp
   else if ((myid.gt.0).and.(myid.lt.nproc-1)) then
      kparastak=kparasta
      kparaendk=kparaend
   endif

   js=0
   jss=0
   do j=1,2*jp
      !
      do i=1,jx
         !
         do k=kparastak,kparaendk
            !
            fg=.5*(r(i,jss,k+1)-r(i,jss,k-1))
            !         f2(i,js,k)=f2(i,js,k)+akapt(isc,i,jss,k)*g23(i,js,k)*fg

            if(tipo(i,jss,k-1).eq.0.or.tipo(i+1,jss,k+1).eq.0)then
               f2(i,js,k)=f2(i,js,k)
            !         elseif(tipo(i+1,jss,k+1).eq.0)then
            !           f2(i,js,k)=f2(i,js,k)
            else
               f2(i,js,k)=f2(i,js,k)+akapt(isc,i,jss,k)*g23(i,js,k)*fg
            end if
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
               !      f2(i,js,k)=f2(i,js,k)+akapt(isc,i,jss,k)*g23(i,js,k)*fg

               if(tipo(i,jss,k+1).eq.0.or.tipo(i+1,jss,k+2).eq.0)then
                  f2(i,js,k)=f2(i,js,k)
               !         elseif(tipo(i+1,jss,k+2).eq.0)then
               !           f2(i,js,k)=f2(i,js,k)
               else
                  f2(i,js,k)=f2(i,js,k)+akapt(isc,i,jss,k)*g23(i,js,k)*fg
               end if

            endif
            !
            if (myid.eq.nproc-1) then

               k=jz
               fg=.5*(3.*r(i,jss,jz)-4.*r(i,jss,jz-1)+r(i,jss,jz-2))
               !      f2(i,js,k)=f2(i,js,k)+akapt(isc,i,jss,k)*g23(i,js,k)*fg

               if(tipo(i,jss,k-1).eq.0.or.tipo(i+1,jss,k-2).eq.0)then
                  f2(i,js,k)=f2(i,js,k)
               !         elseif(tipo(i+1,jss,k-2).eq.0)then
               !           f2(i,js,k)=f2(i,js,k)
               else
                  f2(i,js,k)=f2(i,js,k)+akapt(isc,i,jss,k)*g23(i,js,k)*fg
               end if

            endif
         !
         end do
      !
      enddo
      !
      js=jy
      jss=jy+1
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
            ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j+1,k))
            !
            found_solid=.false.
            do kk=k-1,k+1
               do jj=j,j+1
                  do ii=i,i
                     if(tipo(ii,jj,kk)==0)found_solid=.true.
                  end do
               end do
            end do
	 
            if(found_solid)then
               f2(i,j,k)=f2(i,j,k)
            else
               f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*ak*fg
            end if



         end do
         !
         do kkk=1,kp
            !
            ! check derivative on sides back and front
            !
            if (myid.eq.0) then

               k=1
               r0=.5*(r(i,j,k  )+r(i,j+1,k  ))
               r1=.5*(r(i,j,k+1)+r(i,j+1,k+1))
               r2=.5*(r(i,j,k+2)+r(i,j+1,k+2))
               fg=.5*(-3.*r0+4.*r1-r2)
               ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j+1,k))
               !      f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*ak*fg

               found_solid=.false.
               do kk=k,k+2
                  do jj=j,j+1
                     do ii=i,i
                        if(tipo(ii,jj,kk)==0)found_solid=.true.
                     end do
                  end do
               end do
	 
               if(found_solid)then
                  f2(i,j,k)=f2(i,j,k)
               else
                  f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*ak*fg
               end if


            endif
            !
            if (myid.eq.nproc-1) then

               k=jz
               r0=.5*(r(i,j,jz  )+r(i,j+1,jz  ))
               r1=.5*(r(i,j,jz-1)+r(i,j+1,jz-1))
               r2=.5*(r(i,j,jz-2)+r(i,j+1,jz-2))
               fg=.5*(3.*r0-4.*r1+r2)
               ak=.5*(akapt(isc,i,j,k)+akapt(isc,i,j+1,k))
               !      f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*ak*fg

               found_solid=.false.
               do kk=k-2,k
                  do jj=j,j+1
                     do ii=i,i
                        if(tipo(ii,jj,kk)==0)found_solid=.true.
                     end do
                  end do
               end do
	 
               if(found_solid)then
                  f2(i,j,k)=f2(i,j,k)
               else
                  f2(i,j,k)=f2(i,j,k)+g23(i,j,k)*ak*fg
               end if
            endif
         !
         end do
      !
      enddo
   enddo
   !


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

   return
end
