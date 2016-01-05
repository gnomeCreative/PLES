!***********************************************************************
subroutine flux3(rc,cgra3,r,insc,tipo,body)
   !***********************************************************************
   ! compute explicit flux on component zita
   !
   ! convective  wc*(u,v,w)   with central schema or quick
   ! diffusive   nni*g31*d(u,v,w)/d(csi)
   ! diffusive   nni*g32*d(u,v,w)/d(zita)
   !
   use mysending
   use myarrays_LC
   use myarrays_metri3
   !
   use scala3
   use period
   use convex
   use tipologia
   !
   use mpi

   implicit none

   !-----------------------------------------------------------------------
   !     array declaration
   integer ierr,status(MPI_STATUS_SIZE)
   integer ncolperproc,m,insc
   integer kparastal,kparaendl
   integer kparastall,kparaendll
   double precision bulk_loc,bulk3
   integer i,j,k,ii,jj,ks,kss
   real xf,yf,zf,xc1,yc1,zc1,xc2,yc2,zc2,ddk,dsk
   real fg,r0,r1,r2,an
   real rc(n1,n2,kparasta-1:kparaend) !0:n3)
   real r(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr) !0:n3+1)
   real cgra3(n1,n2,kparasta-1:kparaend+1)
   real ravanti,rindietro1,rindietro2
   integer req1,req2
   integer istatus(MPI_STATUS_SIZE)
   !      integer tipo2(0:n1+1,0:n2+1,kparasta-1:kparaend+1)
   !      integer tipo2(0:n1+1,0:n2+1,0:n3+1) c7
   integer tipo(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr)

   real phic(n1,n2,kparasta:kparaend)
   real den, num, myv
   integer body
   !-----------------------------------------------------------------------
   ! CONVECTIVE TERM  Wc*u:
   ! implemented with central scheme or with quick depending on settings
   ! in points closer to wall it uses uf as ghost cell
   !
   !     sides back and front
   do k=1,kp

      do j=1,jy
         do i=1,jx
            !
            if (myid.eq.0) then
               f3(i,j,0)  = (rc(i,j,0 )+wcs(i,j,0))*r(i,j,0   )- &
                  cgra3(i,j,0 )
            endif

            if (myid.eq.nproc-1) then
               f3(i,j,jz) = (rc(i,j,jz)+wcs(i,j,jz))*r(i,j,jz+1)- &
                  cgra3(i,j,jz)
            endif
         !
         enddo
      enddo
   !
   enddo
   !
   !     into the field
   !
   if (myid.eq.0) then
      kparastal=kp
      kparaendl=kparaend
   else if (myid.eq.nproc-1) then
      kparastal=kparasta
      kparaendl=kparaend-kp
   else
      kparastal=kparasta
      kparaendl=kparaend
   endif


   do k=kparastal,kparaendl
      do j=1,jy
         do i=1,jx
            !
            f3(i,j,k)=(rc(i,j,k)+wcs(i,j,k)) &
               *.5*(r(i,j,k)+r(i,j,k+1))-cgra3(i,j,k)
         !
         end do
      end do
   end do

   ! quick
   if(insc.eq.1)then

      if (myid.eq.0) then
         kparastall=kparasta +2
         kparaendll=kparaend
      else if (myid.eq.nproc-1) then
         kparastall=kparasta
         kparaendll=kparaend-2
      else
         kparastall=kparasta
         kparaendll=kparaend
      endif

      !     sides 5 and 6

      if(myid .eq. 0)then
         do j=1,jy
            do i=1,jx
            
               k=1
               if(rc(i,j,k).gt.0.)then
                  ravanti    =        r(i,j,2)
                  rindietro1 =        r(i,j,1)
                  rindietro2 = kp*(2.*r(i,j,0) - r(i,j,1)) &
                     +(1.-kp)*r(i,j,0)
               else
                  ravanti    =    r(i,j,1)
                  rindietro1 =    r(i,j,2)
                  rindietro2 =    r(i,j,3)
               end if

               if(tipo(i,j,k)==2)then
                  f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                     *.125*( -ravanti +2.*rindietro1 - rindietro2 )
               end if

            end do
         end do
      end if
      
      if(myid .eq. nproc-1)then
         do j=1,jy
            do i=1,jx
               
               k=jz-1
               if(rc(i,j,k).gt.0.)then
                  ravanti    = r(i,j,jz  )
                  rindietro1 = r(i,j,jz-1)
                  rindietro2 = r(i,j,jz-2)
               else
                  ravanti    =        r(i,j,jz-1)
                  rindietro1 =        r(i,j,jz  )
                  rindietro2 = kp*(2.*r(i,j,jz+1) - r(i,j,jz)) &
                     +(1.-kp)*r(i,j,jz+1)
               end if

               if(tipo(i,j,k)==2)then
                  f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                     *.125*( -ravanti +2.*rindietro1 - rindietro2 )
               end if
            end do
         end do
      end if    
       
      !     into the field
      do k=kparastall,kparaendll
         do j=1,jy
            do i=1,jx
               !
               !      if(tipo(i,j,k)==2)then
               !      if (rc(i,j,k).gt.0.) then
               !      f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*
               !     > .125*(-r(i,j,k+1)+2.*r(i,j,k)-r(i,j,k-1))
               !      else
               !      f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*
               !     > .125*(-r(i,j,k)+2.*r(i,j,k+1)-r(i,j,k+2))
               !      end if
               !      end if

               if (rc(i,j,k).gt.0.) then

                  if(tipo(i,j,k)==1)then !ib cell
                     if(tipo(i,j,k+1)==0)then !solido davanti
                        f3(i,j,k)=0.
                     elseif(tipo(i,j,k-1)==0)then !uso diff centrate
                        f3(i,j,k)=f3(i,j,k)
                     else ! solido sta nell'altra direzione quick
                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                           .125*(-r(i,j,k+1)+2.*r(i,j,k)-r(i,j,k-1))

                     endif !solido avanti/indietro/lato
                  else !i è fluido allora quick normale
                     f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                        .125*(-r(i,j,k+1)+2.*r(i,j,k)-r(i,j,k-1))
                  endif !tipo
    
               else !rc(i,j,k).le.0.
   
                  if(tipo(i,j,k+1)==1)then
                     if(tipo(i,j,k)==0)then !solido indietro
                        f3(i,j,k)=0.
                     elseif(tipo(i,j,k+2)==0)then !uso diff centrate
                        f3(i,j,k)=f3(i,j,k)
                     else !solido di lato
                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                           .125*(-r(i,j,k)+2.*r(i,j,k+1)-r(i,j,k+2))
                     endif !solido avanti/indietro/lato
                  else !i+1 è fluido allora quick normale
                     f3(i,j,k)=f3(i,j,k)+rc(i,j,k)* &
                        .125*(-r(i,j,k)+2.*r(i,j,k+1)-r(i,j,k+2))
                  endif !tipo
	
               endif !rc
            end do
         end do
      end do
      
   ! SMART/quick modificato upwind
   elseif(insc.eq.2)then

      ! GIULIA UCS - CGRA
      !     side left and right
      !     sides back and front
      do k=1,kp

         do j=1,jy
            do i=1,jx
               !
               if (myid.eq.0) then
                  f3(i,j,0)  = (rc(i,j,0 )+wcs(i,j,0))*r(i,j,0   )- &
                     cgra3(i,j,0 )
               endif

               if (myid.eq.nproc-1) then
                  f3(i,j,jz) = (rc(i,j,jz)+wcs(i,j,jz))*r(i,j,jz+1)- &
                     cgra3(i,j,jz)
               endif
            !
            enddo
         enddo
      !
      enddo
      
      if (myid.eq.0) then
         kparastal=kp
         kparaendl=kparaend
      else if (myid.eq.nproc-1) then
         kparastal=kparasta
         kparaendl=kparaend-kp
      else
         kparastal=kparasta
         kparaendl=kparaend
      endif


      do k=kparastal,kparaendl
         do j=1,jy
            do i=1,jx
               !
               f3(i,j,k)=wcs(i,j,k) &
                  *.5*(r(i,j,k)+r(i,j,k+1))-cgra3(i,j,k)
            !
            end do
         end do
      end do

      if (myid.eq.0) then
         kparastall=2
         kparaendll=kparaend
      else if (myid.eq.nproc-1) then
         kparastall=kparasta
         kparaendll=kparaend-2
      else
         kparastall=kparasta
         kparaendll=kparaend
      endif

      !     sides 5 and 6 sia per bodyforce 0 o 1
      if(myid .eq. 0)then
         do j=1,jy
            do i=1,jx
               k=1
               if(rc(i,j,k).gt.0.)then
                  ravanti    =        r(i,j,2)
                  rindietro1 =        r(i,j,1)
                  rindietro2 = kp*(2.*r(i,j,0) - r(i,j,1)) &
                     +(1.-kp)*r(i,j,0)
               else
                  ravanti    =    r(i,j,1)
                  rindietro1 =    r(i,j,2)
                  rindietro2 =    r(i,j,3)
               end if

               f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                  *.125*( 3.*ravanti +6.*rindietro1 - rindietro2 )

      
            end do
         end do

      end if !myid==0
      
      if(myid .eq. nproc-1)then
         do j=1,jy
            do i=1,jx
               
               k=jz-1
               if(rc(i,j,k).gt.0.)then
                  ravanti    = r(i,j,jz  )
                  rindietro1 = r(i,j,jz-1)
                  rindietro2 = r(i,j,jz-2)
               else
                  ravanti    =        r(i,j,jz-1)
                  rindietro1 =        r(i,j,jz  )
                  rindietro2 = kp*(2.*r(i,j,jz+1) - r(i,j,jz)) &
                     +(1.-kp)*r(i,j,jz+1)
               end if


               f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                  *.125*( 3.*ravanti +6.*rindietro1 - rindietro2 )

            end do
         end do
      end if    
  
   
      
      !     into the field without bodyforce
      if(body.eq.0)then
         do k=kparastall,kparaendll
            do j=1,jy
               do i=1,jx
                  if(rc(i,j,k).gt.0.)then

                     f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                        *.125*(3.*r(i,j,k+1)+6*r(i,j,k)-r(i,j,k-1))

                  else  !(rc(i,j,k).le.0.)

                     f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                        *.125*(3.*r(i,j,k)+6*r(i,j,k+1)-r(i,j,k+2)) !my index +1


                  endif !rc</>0
               end do
            end do
         end do

      !
      else !bodyforce attivo
         do k=kparastall,kparaendll
            do j=1,jy
               do i=1,jx
                  if(rc(i,j,k).gt.0.)then
                     if(tipo(i,j,k)==1)then !ib cell
                        if(tipo(i,j,k+1)==0)then !solido davanti

                           f3(i,j,k)=0.
	    
                        elseif(tipo(i,j,k-1)==0)then !solido dietro

                           f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)
                        else !solido lati

                           f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                              *.125*(3.*r(i,j,k+1)+6*r(i,j,k)-r(i,j,k-1))

                        endif !solido avanti/indietro/lati

                     else !k non è IB procedo come se non ci fossero gli ib

                        ! se k-1 è ib faccio la correzione per evitare la scia
                        !	      if(tipo(i,j,k-1)==1)then !la cella prima è un ib
                        !		if(tipo(i,j,k-2)==0)then !perché considero anche questa?
                        !		  rindietro2 =r(i,j,k)
                        !		endif
                        !	      endif !fine k-1 IB


                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                           *.125*(3.*r(i,j,k+1)+6*r(i,j,k)-r(i,j,k-1))


                     endif !tip02

                  else !rc(i,j,k).le.0.
   
                     if(tipo(i,j,k+1)==1)then !ib
                        if(tipo(i,j,k)==0)then !solido avanti
                           f3(i,j,k)=0.
                        elseif(tipo(i,j,k+2)==0)then !solido indietro
                           f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)
                        else !solido lati

                           f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                              *.125*(3.*r(i,j,k)+6*r(i,j,k+1)-r(i,j,k+2))

            
                        endif !solido avanti indietro e lati


                     else !k non è IB procedo come se non ci fossero gli ib

                        ! se k-1 è ib faccio la correzione per evitare la scia
                        !	      if(tipo(i,j,k+1)==1)then
                        !		if(tipo(i,j,k+2)==0)then
                        !		rindietro2=r(i,j,k+1)
                        !		endif
                        !	      endif

                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k) &
                           *.125*(3.*r(i,j,k)+6*r(i,j,k+1)-r(i,j,k+2))

                     endif !tipo
	
                  endif !rc
      
               end do
            end do
         end do

      endif !body yes/no
   !	endif !SMART
   !-----SMART /quick modificato upwind
   !  upwind
   elseif(insc.eq.3)then

      ! GIULIA UCS - CGRA
      !     side left and right
      !     sides back and front
      do k=1,kp

         do j=1,jy
            do i=1,jx
               !
               if (myid.eq.0) then
                  f3(i,j,0)  = (rc(i,j,0 )+wcs(i,j,0))*r(i,j,0   )- &
                     cgra3(i,j,0 )
               endif

               if (myid.eq.nproc-1) then
                  f3(i,j,jz) = (rc(i,j,jz)+wcs(i,j,jz))*r(i,j,jz+1)- &
                     cgra3(i,j,jz)
               endif
            !
            enddo
         enddo
      !
      enddo
      
      if (myid.eq.0) then
         kparastal=kp
         kparaendl=kparaend
      else if (myid.eq.nproc-1) then
         kparastal=kparasta
         kparaendl=kparaend-kp
      else
         kparastal=kparasta
         kparaendl=kparaend
      endif


      do k=kparastal,kparaendl
         do j=1,jy
            do i=1,jx
               !
               f3(i,j,k)=wcs(i,j,k) &
                  *.5*(r(i,j,k)+r(i,j,k+1))-cgra3(i,j,k)
            !
            end do
         end do
      end do

      if (myid.eq.0) then
         kparastall=2
         kparaendll=kparaend
      else if (myid.eq.nproc-1) then
         kparastall=kparasta
         kparaendll=kparaend-2
      else
         kparastall=kparasta
         kparaendll=kparaend
      endif

      !     sides 5 and 6 sia per bodyforce 0 o 1
      if(myid .eq. 0)then
         do j=1,jy
            do i=1,jx
               k=1
               if(rc(i,j,k).gt.0.)then
                  ravanti    =        r(i,j,2)
                  rindietro1 =        r(i,j,1)
                  rindietro2 = kp*(2.*r(i,j,0) - r(i,j,1)) &
                     +(1.-kp)*r(i,j,0)
               else
                  ravanti    =    r(i,j,1)
                  rindietro1 =    r(i,j,2)
                  rindietro2 =    r(i,j,3)
               end if

               f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*rindietro1

      
            end do
         end do

      end if !myid==0
      
      if(myid .eq. nproc-1)then
         do j=1,jy
            do i=1,jx
               
               k=jz-1
               if(rc(i,j,k).gt.0.)then
                  ravanti    = r(i,j,jz  )
                  rindietro1 = r(i,j,jz-1)
                  rindietro2 = r(i,j,jz-2)
               else
                  ravanti    =        r(i,j,jz-1)
                  rindietro1 =        r(i,j,jz  )
                  rindietro2 = kp*(2.*r(i,j,jz+1) - r(i,j,jz)) &
                     +(1.-kp)*r(i,j,jz+1)
               end if


               f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*rindietro1

            end do
         end do
      end if    
  
   
      
      !     into the field without bodyforce
      if(body.eq.0)then
         do k=kparastall,kparaendll
            do j=1,jy
               do i=1,jx
                  if(rc(i,j,k).gt.0.)then

                     f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)

                  else  !(rc(i,j,k).le.0.)

                     f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)


                  endif !rc</>0
               end do
            end do
         end do

      !
      else !bodyforce attivo
         do k=kparastall,kparaendll
            do j=1,jy
               do i=1,jx
                  if(rc(i,j,k).gt.0.)then
                     if(tipo(i,j,k)==1)then !ib cell
                        if(tipo(i,j,k+1)==0)then !solido davanti

                           f3(i,j,k)=0.
	    

                        else !solido lati

                           f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)

                        endif !solido avanti/indietro/lati

                     else !k non è IB procedo come se non ci fossero gli ib

                        ! se k-1 è ib faccio la correzione per evitare la scia
                        !	      if(tipo(i,j,k-1)==1)then !la cella prima è un ib
                        !		if(tipo(i,j,k-2)==0)then !perché considero anche questa?
                        !		  rindietro2 =r(i,j,k)
                        !		endif
                        !	      endif !fine k-1 IB


                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k)


                     endif !tip02

                  else !rc(i,j,k).le.0.
   
                     if(tipo(i,j,k+1)==1)then !ib
                        if(tipo(i,j,k)==0)then !solido avanti
                           f3(i,j,k)=0.

                        else !solido lati

                           f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)

            
                        endif !solido avanti indietro e lati


                     else !k non è IB procedo come se non ci fossero gli ib

                        ! se k-1 è ib faccio la correzione per evitare la scia
                        !	      if(tipo(i,j,k+1)==1)then
                        !		if(tipo(i,j,k+2)==0)then
                        !		rindietro2=r(i,j,k+1)
                        !		endif
                        !	      endif

                        f3(i,j,k)=f3(i,j,k)+rc(i,j,k)*r(i,j,k+1)

                     endif !tipo
	
                  endif !rc
      
               end do
            end do
         end do

      endif !body yes/no
   endif !SMART
   !-----SMART /quick modificato upwind

   !----------------------------------------------------------------------
   ! DIFFUSIVE TERM NNI*G31*DU/D(CSI)
   !
   !     sides back and front
   ks=0
   kss=0
   do k=1,2*kp

      if (myid.eq.(k-1)*(nproc-1)) then

         do j=1,jy
            !
            do i=1+ip,jx-ip
               !
               fg=.5*(r(i+1,j,kss)-r(i-1,j,kss))
               f3(i,j,ks)=-f3(i,j,ks) &
                  +annit(i,j,kss)*g31(i,j,ks)*fg
            !
            end do
            !
            do ii=1,ip
               !
               ! check derivative on sides left and right
               !
               i=1
               fg=.5*(-3.*r(i,j,kss)+4.*r(i+1,j,kss)-r(i+2,j,kss))
               f3(i,j,ks)=-f3(i,j,ks) &
                  +annit(i,j,kss)*g31(i,j,ks)*fg
               !
               i=jx
               fg=.5*(3.*r(jx,j,kss)-4.*r(jx-1,j,kss)+r(jx-2,j,kss))
               f3(i,j,ks)=-f3(i,j,ks) &
                  +annit(i,j,kss)*g31(i,j,ks)*fg
            !
            end do
         !
         enddo
      !
      endif

      ks=jz
      kss=jz+1
   !
   enddo
   !
   !     into the field
   do k=kparastal,kparaendl
      do j=1,jy
         !
         do i=1+ip,jx-ip
            !
            r1=.5*(r(i-1,j,k)+r(i-1,j,k+1))
            r2=.5*(r(i+1,j,k)+r(i+1,j,k+1))
            fg=.5*(r2-r1)
            an=.5*(annit(i,j,k)+annit(i,j,k+1))
            f3(i,j,k)=-f3(i,j,k)+g31(i,j,k)*an*fg
         !
         end do
         !
         do ii=1,ip
            !
            ! check derivative on sides left and right
            !
            i=1
            r0=.5*(r(i  ,j,k)+r(i  ,j,k+1))
            r1=.5*(r(i+1,j,k)+r(i+1,j,k+1))
            r2=.5*(r(i+2,j,k)+r(i+2,j,k+1))
            fg=.5*(-3.*r0+4.*r1-r2)
            an=.5*(annit(i,j,k)+annit(i,j,k+1))
            f3(i,j,k)=-f3(i,j,k)+g31(i,j,k)*an*fg
            !
            i=jx
            r0=.5*(r(jx  ,j,k)+r(jx  ,j,k+1))
            r1=.5*(r(jx-1,j,k)+r(jx-1,j,k+1))
            r2=.5*(r(jx-2,j,k)+r(jx-2,j,k+1))
            fg=.5*(3.*r0-4.*r1+r2)
            an=.5*(annit(i,j,k)+annit(i,j,k+1))
            f3(i,j,k)=-f3(i,j,k)+g31(i,j,k)*an*fg
         !
         end do
      !
      enddo
   enddo
   !
   !-----------------------------------------------------------------------
   ! DIFFUSIVE TERM NNI*G32*DU/D(ZITA)
   !
   !     sides back and front
   !
   ks=0
   kss=0

   do k=1,2*kp

      if (myid.eq.(k-1)*(nproc-1)) then

         do i=1,jx
            !
            do j=1+jp,jy-jp
               !
               fg=.5*(r(i,j+1,kss)-r(i,j-1,kss))
               f3(i,j,ks)=f3(i,j,ks) &
                  +annit(i,j,kss)*g32(i,j,ks)*fg
            !
            end do
            !
            do jj=1,jp
               !
               ! check derivative on sides back and front
               !
               j=1
               fg=.5*(-3.*r(i,j,kss)+4.*r(i,j+1,kss)-r(i,j+2,kss))
               f3(i,j,ks)=f3(i,j,ks) &
                  +annit(i,j,kss)*g32(i,j,ks)*fg
               !
               j=jy
               fg=.5*(3.*r(i,jy,kss)-4.*r(i,jy-1,kss)+r(i,jy-2,kss))
               f3(i,j,ks)=f3(i,j,ks) &
                  +annit(i,j,kss)*g32(i,j,ks)*fg
            !
            end do
         !
         enddo
      !
      endif

      ks=jz
      kss=jz+1
   !
   enddo

   !
   !     into the field
   do k=kparastal,kparaendl
      do i=1,jx
         !
         do j=1+jp,jy-jp
            !
            r1=.5*(r(i,j-1,k)+r(i,j-1,k+1))
            r2=.5*(r(i,j+1,k)+r(i,j+1,k+1))
            fg=.5*(r2-r1)
            an=.5*(annit(i,j,k)+annit(i,j,k+1))
            f3(i,j,k)=f3(i,j,k)+g32(i,j,k)*an*fg
         !
         end do
         !
         do jj=1,jp
            !
            ! check derivative on sides back and front
            !
            j=1
            r0=.5*(r(i,j  ,k)+r(i,j  ,k+1))
            r1=.5*(r(i,j+1,k)+r(i,j+1,k+1))
            r2=.5*(r(i,j+2,k)+r(i,j+2,k+1))
            fg=.5*(-3.*r0+4.*r1-r2)
            an=.5*(annit(i,j,k)+annit(i,j,k+1))
            f3(i,j,k)=f3(i,j,k)+g32(i,j,k)*an*fg
            !
            j=jy
            r0=.5*(r(i,jy  ,k)+r(i,jy  ,k+1))
            r1=.5*(r(i,jy-1,k)+r(i,jy-1,k+1))
            r2=.5*(r(i,jy-2,k)+r(i,jy-2,k+1))
            fg=.5*(3.*r0-4.*r1+r2)
            an=.5*(annit(i,j,k)+annit(i,j,k+1))
            f3(i,j,k)=f3(i,j,k)+g32(i,j,k)*an*fg
         !
         end do
      !
      enddo
   enddo
   !

   do j=1,jy
      do i=1,jx
         do k=kparastal,kparaendl
            !
            if(body.eq.1)then
               if(tipo(i,j,k).eq.0)f3(i,j,k)=0.
               if(k.lt.jz)then
                  if(tipo(i,j,k+1).eq.0)f3(i,j,k)=0.
               endif
            endif
         enddo
      enddo
   enddo
   !     pass f3 at the border between procs
   if(myid.eq.0)then
      leftpem=MPI_PROC_NULL
      rightpem=rightpe
   else if(myid.eq.nproc-1)then
      leftpem=leftpe
      rightpem=MPI_PROC_NULL
   else
      leftpem=leftpe
      rightpem=rightpe
   endif

   if(rightpem /= MPI_PROC_NULL) then
      call MPI_ISEND(f3(1,1,kparaend),jx*jy,MPI_REAL_SD, &
         rightpem ,tagrs,MPI_COMM_WORLD,req1,ierr)
   endif
   if(leftpem /= MPI_PROC_NULL) then
      call MPI_IRECV(f3(1,1,kparasta-1),jx*jy,MPI_REAL_SD, &
         leftpem  ,taglr,MPI_COMM_WORLD,req2,ierr)
   endif

   if(rightpem /= MPI_PROC_NULL) then
      call MPI_WAIT(req1,istatus,ierr)
   endif
   if(leftpem /= MPI_PROC_NULL) then
      call MPI_WAIT(req2,istatus,ierr)
   endif

   !     integral for f3
   bulk_loc=0.

   do k=kparasta,kparaend
      do j=1,jy
         do i=1,jx
            bulk_loc=bulk_loc+f3(i,j,k)-f3(i,j,k-1)
         end do
      end do
   end do

   !     make the value known to all procs

   call MPI_ALLREDUCE(bulk_loc,bulk3,1,MPI_DOUBLE_PRECISION, &
      MPI_SUM, &
      MPI_COMM_WORLD,ierr)

   !    now bulk is known by all procs

   bulk=bulk+bulk3
   !
   return
end
