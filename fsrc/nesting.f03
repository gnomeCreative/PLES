!***********************************************************************
subroutine nesting(bodyforce,ti,area1,area2,area5,area6,&
   kparasta,kparaend,myid,nproc,freesurface)
   !***********************************************************************
   ! subroutine for nesting procedure
   use myarrays_nesting
   !
   use scala3
   use tipologia
   use orl
   !
   use mpi


   implicit none

   !-----------------------------------------------------------------------
   ! array declaration
   integer bodyforce
   integer tii,myid,nproc,kparasta,kparaend
   real area1,area2,area3,area4,area5,area6
   real ti
   integer freesurface
   !-----------------------------------------------------------------------
   ! check that dt not larger than t_new_pom

   if(ti<=ti_pom_new)then

      call interpola_pareti_pom(ti)
      call contrin_lat

      if(freesurface.eq.0)then !if freesurface off
         if(myid.eq.0)then
            write(*,*)'freesurface is off and entering redistribuzione.'
         end if
         call redistribuzione(bodyforce,area1,area2,area5,area6,kparasta,kparaend,myid,nproc)
      elseif(freesurface.eq.1)then !if freesurface is on
         if(myid.eq.0)then
            write(*,*)'free surface is on. redistribuzione skipped.'
         end if
      end if !if freesurface on/off


   elseif(ti>ti_pom_new.and.ti<ti_pom_fin)then

      if(myid==0)write(*,*)'NESTING READ NEW PLANE'

      ti=ti_pom_new
      ti_pom_old=ti_pom_new
      
      if(infout1 /=0)then
         uo1=un1
         vo1=vn1
         wo1=wn1
         rhovo1=rhovn1
      end if
      
      if(infout2 /=0)then
         uo2=un2
         vo2=vn2
         wo2=wn2
         rhovo2=rhovn2
      end if
                  
      if(infout5 /=0)then
         uo5=un5
         vo5=vn5
         wo5=wn5
         rhovo5=rhovn5
      end if
                  
      if(infout6 /=0)then
         uo6=un6
         vo6=vn6
         wo6=wn6
         rhovo6=rhovn6
      end if

      call leggi_pareti_pom(ti)
      call interpola_pareti_pom(ti)
      call contrin_lat

      if(freesurface.eq.0)then !if freesurface off
         if(myid.eq.0)then
            write(*,*)'freesurface is off and entering redistribuzione.'
         end if
         call redistribuzione(bodyforce,area1,area2,area5,area6,kparasta,kparaend,myid,nproc)
      elseif(freesurface.eq.1)then !if freesurface is on
         if(myid.eq.0)then
            write(*,*)'free surface is on. redistribuzione skipped.'
         end if
      end if !if freesurface on/off

      
   end if
   !
   ! if ti greater than t_pom_fin the computation end
   if(ti>=ti_pom_fin)then
      ti=ti_pom_fin
      termina=.true.
   !cc      close(1313)
   end if

   return
end
