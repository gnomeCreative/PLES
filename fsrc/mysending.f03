!***********************************************************************
      module mysending
!***********************************************************************
! contains info for the message passaging
!-----------------------------------------------------------------------
      use iso_c_binding
      integer tagls,taglr,tagrs,tagrr
      integer leftpe,rightpe
      integer leftpem,rightpem
      integer,bind(C) :: myid,nproc
      integer kparasta,kparaend
      integer iparasta,iparaend
      
      integer deepl,deepr
      integer deepgl,deepgr
      integer deep_mul


!***********************************************************************
      end module mysending
!***********************************************************************
