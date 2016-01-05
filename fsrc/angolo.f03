!***********************************************************************
       subroutine angolo ( x1, y1, z1, x2, y2, z2, x3, &
                          y3, z3, x4, y4, z4, angle_deg, cor )
!***********************************************************************
! THIS SUBROUTINE IS TAKEN FROM LIBRARY GEOMETRY PACK
!
!c LINES_EXP_ANGLE_3D finds the angle between two explicit lines in 3D.
!
!
!  Formula:
!
!    The explicit form of a line in 3D is:
!
!      (X1,Y1,Z1), (X2,Y2,Z2).
!
!  Modified:
!
!    24 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, two distince points on the first line.
!
!    Input, real X3, Y3, Z3, X4, Y4, Z4, two distinct points on the second line.
!
!    Output, real ANGLE, the angle in radians between the two lines.
!    ANGLE is computed using the ACOS function, and so lies between 0 and PI.
!    But if one of the lines is degenerate, ANGLE is returned as -1.0.
!
       implicit none
!
       real:: angle
       real:: pi
       real:: angle_deg
       real:: arc_cosine
       real:: ctheta
       real:: enorm0_3d
       real:: pdotq
       real:: pnorm
       real:: qnorm
       real:: x1,x2,x3,x4
       real:: y1,y2,y3,y4
       real:: z1,z2,z3,z4

       integer :: cor
       real :: tol

       tol=.000001

       pi=acos(-1.)

       pnorm = sqrt ( ( x1 - x2 )**2 + ( y1 - y2 )**2 + ( z1 - z2 )**2 )
       qnorm = sqrt ( ( x3 - x4 )**2 + ( y3 - y4 )**2 + ( z3 - z4 )**2 )

       pdotq =    ( x2 - x1 ) * ( x4 - x3 ) &
           + ( y2 - y1 ) * ( y4 - y3 ) &
           + ( z2 - z1 ) * ( z4 - z3 )

       if ( pnorm.eq.0.0 .or. qnorm.eq.0.0 ) then
         write (*,*) ' '
         write (*,*) 'LINES_EXP_ANGLE_3D - Fatal error!'
         write (*,*) '  One of the lines is degenerate!'
         angle = - 1.0
       else

!      write(*,*)'pi: ',pi
         ctheta = pdotq / ( pnorm * qnorm )

         if(abs(ctheta).ge.1.)then
            write(1000,*)'TROVATO ANGOLO DEGENERE',ctheta
            cor=1
            angle_deg = 0.0
         else
            angle =abs(acos(ctheta)) ! arc_cosine ( ctheta )
!         write(*,*)'angolo',angle
            angle_deg=angle*180./pi
            angle_deg=angle
!         write(*,*)'angolo deg: ',angle_deg
         end if
       end if


       return
       end
