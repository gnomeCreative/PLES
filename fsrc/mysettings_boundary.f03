!c***********************************************************************
      module mysettings_boundary
!***********************************************************************
! contains settings for boundary conditions at wall
!-----------------------------------------------------------------------
use iso_c_binding

      integer,bind(C) :: iboun1,iboun2
      integer,bind(C) :: iboun3,iboun4
      integer,bind(C) :: iboun5,iboun6

      integer,bind(C) :: ibodybuffer1,ibodybuffer2
      integer,bind(C) :: ibodybuffer3,ibodybuffer4
      integer,bind(C) :: ibodybuffer5,ibodybuffer6
      
      end module 
