!***********************************************************************
subroutine iniz_metrica(ncolperproc,nlevel,nlevmultimax)
   !***********************************************************************
   use myarrays_metri3
   use mysending
   use scala3

   implicit none
   !-----------------------------------------------------------------------
   !     arrays declaration
   integer i,j,k
   integer kpsta,kpend
   integer ncolperproc,nlevel,nlevmultimax
   integer var_piani,livello
   !-----------------------------------------------------------------------
   ! griglia

   !    qui decido la profondita' di allocazione della griglia
   !    per ogni processore in funzione della subroutine mul_met
   !    questo perche' in relazione al nlevel di multigrid la
   !    allocazione locale di griglia deve essere maggiore
   var_piani=ncolperproc
   livello=1
   do while(var_piani.ne.1.and.livello.ne.5)
      var_piani=int(var_piani/2)
      livello=livello+1
   end do
   livello=livello-1
   !      write(*,*)myid,'multigrid level for each proc',livello
   if(livello.ne.1)then
      deep_mul=2*2**(livello-2) ! allocation deep
   !        write(*,*)'grid allocation deep',deep_mul
   else
      if(myid==0)then
         write(*,*)'PROBLEM on multigrid level for proc'
      end if
      stop
   end if
      
   if(livello .lt. nlevel)then
      if(myid==0)then
         write(*,*)'PROBLEM on multigrid level for proc'
         write(*,*)'level reduced from ',nlevel,' to ',livello
      end if
      nlevel = livello
      nlevmultimax = nlevel
   end if

   !-----------------------------------------------------------------------
   !     coordinates allocation
   if(myid .eq. 0)then
      allocate (x(-8:n1+8,-8:n2+8,-8:kparaend+deep_mul+1))
      allocate (y(-8:n1+8,-8:n2+8,-8:kparaend+deep_mul+1))
      allocate (z(-8:n1+8,-8:n2+8,-8:kparaend+deep_mul+1))
   elseif(myid.eq.nproc-1)then
      allocate (x(-8:n1+8,-8:n2+8,kparasta-deep_mul-1:kparaend+8))
      allocate (y(-8:n1+8,-8:n2+8,kparasta-deep_mul-1:kparaend+8))
      allocate (z(-8:n1+8,-8:n2+8,kparasta-deep_mul-1:kparaend+8))
   else
      allocate(x(-8:n1+8,-8:n2+8, &
         kparasta-deep_mul-1:kparaend+deep_mul+1))
      allocate(y(-8:n1+8,-8:n2+8, &
         kparasta-deep_mul-1:kparaend+deep_mul+1))
      allocate(z(-8:n1+8,-8:n2+8, &
         kparasta-deep_mul-1:kparaend+deep_mul+1))
   end if
      
   x = 0.
   y = 0.
   z = 0.
         
   if(myid.eq.0)write(*,*)myid,'coordinates allocation done'
   !-----------------------------------------------------------------------
   !     metric terms allocation

   if(myid.eq.0)then
      kpsta = kparasta
      kpend = kparaend !+1
   elseif(myid.eq.nproc-1)then
      kpsta = kparasta !-1
      kpend = kparaend
   else
      kpsta = kparasta !-1
      kpend = kparaend !+1
   end if

   allocate (csx(0:n1,1:n2,kpsta:kpend))
   allocate (csy(0:n1,1:n2,kpsta:kpend))
   allocate (csz(0:n1,1:n2,kpsta:kpend))
   allocate (g11(0:n1,1:n2,kpsta:kpend))
   allocate (g12(0:n1,1:n2,kpsta:kpend))
   allocate (g13(0:n1,1:n2,kpsta:kpend))

   allocate (etx(1:n1,0:n2,kpsta:kpend))
   allocate (ety(1:n1,0:n2,kpsta:kpend))
   allocate (etz(1:n1,0:n2,kpsta:kpend))
   allocate (g21(1:n1,0:n2,kpsta:kpend))
   allocate (g22(1:n1,0:n2,kpsta:kpend))
   allocate (g23(1:n1,0:n2,kpsta:kpend))

   allocate (ztx(1:n1,1:n2,kparasta-1:kpend)) !NOTA: ksta per myid=0 vale 0
   allocate (zty(1:n1,1:n2,kparasta-1:kpend))
   allocate (ztz(1:n1,1:n2,kparasta-1:kpend))
   allocate (g31(1:n1,1:n2,kparasta-1:kpend))
   allocate (g32(1:n1,1:n2,kparasta-1:kpend))
   allocate (g33(1:n1,1:n2,kparasta-1:kpend))


   allocate (giac(1:n1,1:n2,kparasta:kparaend))
      
   csx = 0.
   csy = 0.
   csz = 0.

   etx = 0.
   ety = 0.
   etz = 0.

   ztx = 0.
   zty = 0.
   ztz = 0.

   g11 = 0.
   g12 = 0.
   g13 = 0.

   g21 = 0.
   g22 = 0.
   g23 = 0.

   g31 = 0.
   g32 = 0.
   g33 = 0.

   giac = 0.
      
   !CORIOLISOIL
   allocate (g_co11(0:n1,kpsta:kpend))
   allocate (g_co12(0:n1,kpsta:kpend))
   allocate (g_co13(0:n1,kpsta:kpend))

   allocate (g_co31(1:n1,kparasta-1:kpend))
   allocate (g_co32(1:n1,kparasta-1:kpend))
   allocate (g_co33(1:n1,kparasta-1:kpend))

   g_co11 = 0.
   g_co12 = 0.
   g_co13 = 0.

   g_co31 = 0.
   g_co32 = 0.
   g_co33 = 0.
      
   if(myid.eq.0)write(*,*)myid,'metric allocation done'
   !-----------------------------------------------------------------------
   !     fluxes allocation

   allocate (f1(0:n1,1:n2,kparasta  :kparaend))
   allocate (f2(1:n1,0:n2,kparasta  :kparaend))
   allocate (f3(1:n1,1:n2,kparasta-1:kparaend))
      
   f1 = 0.
   f2 = 0.
   f3 = 0.
      
   if(myid.eq.0)write(*,*)myid,'fluxes allocation done'
   !-----------------------------------------------------------------------

   ! viscosita' turbolenta

   allocate ( annit(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
   allocate (annitV(0:n1+1,0:n2+1,kparasta-deepl:kparaend+deepr))
      
   annit  = 0.
   annitV = 0.
      
   if(myid.eq.0)then
      allocate( annit_piano(0:n1+1,0:n2+1,n3:n3))
      allocate(annitV_piano(0:n1+1,0:n2+1,n3:n3))
      annit_piano  = 0.
      annitV_piano = 0.
   elseif(myid.eq. nproc-1)then
      allocate( annit_piano(0:n1+1,0:n2+1,1:1))
      allocate(annitV_piano(0:n1+1,0:n2+1,1:1))
      annit_piano  = 0.
      annitV_piano = 0.
   end if

   if(myid.eq.0)write(*,*)myid,'eddy viscosity allocation done'
      
   return
end


