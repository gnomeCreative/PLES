module trilinear

    use myarrays_metri3, only: xcd,ycd,zcd
    use mysending, only: myid, nproc
    use geometricRoutines

    implicit none

    ! for trilinear computation
    integer :: index_triangle(8,3)
    integer :: quadrilatero(8,4)
    integer :: index_points(9,2)

contains

    subroutine find_trilinear_coef2(l,trind,trilinear_coef, &
        pp_position,icheck_print,i_ib,j_ib,k_ib,num_ib)

        implicit none
        !-----------------------------------------------------------------------
        real, parameter :: tolE=1.10-6
        integer num_ib
        integer icheck_print
        integer i_ib,j_ib,k_ib
        integer ival
        integer i,np
        integer il(4),jl(4),kl(4)
        integer nx,ny,nz
        integer l
        integer sum_degenere,idegenere
        integer trind(num_ib,4,3)

        real trilinear_coef(4)
        real pp_position(3)
        real v(4,3),v_pro(4,2),PP(3),PP_pro(2)
        real dir(3),direction,val,dist,t
        real p1(2),p2(2),p3(2),p4(2),M(2),N(2),pn(2),delta1,delta2
        real PA(2),PB(2)
        real dA1,dA2,dB3,dB4,dPA,dPB
        real dist12,dist14
        real sum
        real area(3) !area1,area2,area3

        !-----------------------------------------------------------------------
        ! point definition
        PP(:) = pp_position(:)

        do np = 1,4
            il(np) = trind(l,np,1)
            jl(np) = trind(l,np,2)
            kl(np) = trind(l,np,3)
            !write(*,*)'np',il(np),jl(np),kl(np)
            ! point of quadrilater
            v(np,1) = xcd(il(np),jl(np),kl(np))
            v(np,2) = ycd(il(np),jl(np),kl(np))
            v(np,3) = zcd(il(np),jl(np),kl(np))
        end do
        !-----------------------------------------------------------------------
        ! decide the direction for the projection 3D-2D
        ! dir(1) = abs( v(1,1) - v(3,1) )
        ! dir(2) = abs( v(1,2) - v(3,2) )
        ! dir(3) = abs( v(1,3) - v(3,3) )

        ! val = 1.d10
        ! do i=1,3
        !    if (dir(i) < val) then

        !       direction = i
        !       val = dir(i)
        !    end if
        ! end do


        ! direction= 1
        ! dir 1
        !check line 14
        ! if ( abs(v(1,2) - v(4,2))==tolE .and. abs(v(1,3) - v(4,3))==tolE ) then
        !    direction=2
        ! else
          !check line 12
        !    if ( abs(v(1,2) - v(2,2))==tolE .and. abs(v(1,3) - v(2,3))==tolE  ) then
        !      direction=2
        !    end if
        ! end if


        ! dir 2
        ! if (direction==2) then
        !   if ( abs(v(1,1) - v(4,1))==tolE .and. abs(v(1,3) - v(4,3))==tolE ) then
        !      direction=3
        !   else
        !      if ( abs(v(1,1) - v(2,1))==tolE .and. abs(v(1,3) - v(2,3))==tolE ) then
        !        direction=3
        !      end if
        !   end if
        ! end if


        ! dir 3
        ! if ( (v(1,1) - v(4,1))==0 .and. (v(1,2) - v(4,2))==0 )


        !  dir 1
        dist12 = sqrt((v(1,2) - v(2,2))**2. + (v(1,3) - v(2,3))**2.)
        dist14 = sqrt((v(1,2) - v(4,2))**2. + (v(1,3) - v(4,3))**2.)
        area(1) = dist12*dist14
        !  dir 2
        dist12 = sqrt((v(1,1) - v(2,1))**2. + (v(1,3) - v(2,3))**2.)
        dist14 = sqrt((v(1,1) - v(4,1))**2. + (v(1,3) - v(4,3))**2.)
        area(2) = dist12*dist14
        !  dir 3
        dist12 = sqrt((v(1,2) - v(2,2))**2. + (v(1,1) - v(2,1))**2.)
        dist14 = sqrt((v(1,2) - v(4,2))**2. + (v(1,1) - v(4,1))**2.)
        area(3) = dist12*dist14

        !  if (area1 /= 0.)direction=1
        !  if (area2 /= 0.)direction=2
        !  if (area3 /= 0.)direction=3

        val = 0.
        do i=1,3
            if (area(i) > val) then
                direction = i
                val = area(i)
            end if
        end do

        !-----------------------------------------------------------------------
        ! 3D to 2D
        ! projection for the point
        if (direction==1) then
            do np=1,4
                ! j,k
                v_pro(np,1)=v(np,2)
                v_pro(np,2)=v(np,3)
            end do

            PP_pro(1) = PP(2)
            PP_pro(2) = PP(3)
        elseif (direction==2) then
            do np=1,4
                ! i,k
                v_pro(np,1)=v(np,1)
                v_pro(np,2)=v(np,3)
            end do

            PP_pro(1) = PP(1)
            PP_pro(2) = PP(3)
        elseif (direction==3) then
            do np=1,4
                ! i,j
                v_pro(np,1)=v(np,1)
                v_pro(np,2)=v(np,2)
            end do

            PP_pro(1) = PP(1)
            PP_pro(2) = PP(2)
        end if
        !-----------------------------------------------------------------------
        ! find on line 14 the closer point to PP
        do i=1,2
            p1(i)=v_pro(1,i)
            p2(i)=v_pro(2,i)
            p3(i)=v_pro(3,i)
            p4(i)=v_pro(4,i)
        end do




        if (p1(1)==p4(1) .and. p1(2)==p4(2)) then
            write(*,*)'line 14'
            write(*,*)'direction',direction,dir(1),dir(2),dir(3)
            write(*,*)v(1,1),v(1,2),v(1,3)
            write(*,*)v(2,1),v(2,2),v(2,3)
            write(*,*)v(3,1),v(3,2),v(3,3)
            write(*,*)v(4,1),v(4,2),v(4,3)
            write(*,*)myid,p1(1),p1(2),p2(1),p2(2),p3(1),p3(2),p4(1),p4(2)
        end if


        if (p1(1)==p2(1) .and. p1(2)==p2(2)) then
            write(*,*)'line 12'
            write(*,*)'direction',direction,dir(1),dir(2),dir(3)
            write(*,*)v(1,1),v(1,2),v(1,3)
            write(*,*)v(2,1),v(2,2),v(2,3)
            write(*,*)v(3,1),v(3,2),v(3,3)
            write(*,*)v(4,1),v(4,2),v(4,3)
            write(*,*)myid,p1(1),p1(2),p2(1),p2(2),p3(1),p3(2),p4(1),p4(2)
        end if




        call line_exp_point_near_2d ( p1, p4, PP_pro, pn, dist, t )

        delta1 = PN(1)-PP_pro(1)
        delta2 = PN(2)-PP_pro(2)

        ! M and N represent two points on a parallel line to line 14 passing through PP
        M(1) = p1(1)-delta1
        M(2) = p1(2)-delta2

        N(1) = p4(1)-delta1
        N(2) = p4(2)-delta2

        if (M(1)==N(1) .and. M(2)==N(2)) then
            write(*,*)'line 12',delta1,delta2
            write(*,*)'direction',direction,dir(1),dir(2),dir(3)
            write(*,*)'M and N',M(1),M(2),N(1),N(2)
            write(*,*)v(1,1),v(1,2),v(1,3)
            write(*,*)v(2,1),v(2,2),v(2,3)
            write(*,*)v(3,1),v(3,2),v(3,3)
            write(*,*)v(4,1),v(4,2),v(4,3)
            write(*,*)myid,p1(1),p1(2),p2(1),p2(2),p3(1),p3(2),p4(1),p4(2)
        end if




        !intersection with  line12 of MN line
        call lines_exp_int_2d(p1,p2,M,N,ival,PA)
        if (ival==0) then
            write(*,*)'INTER: no intersection',i_ib,j_ib,k_ib
            write(*,*)'line 12',delta1,delta2
            write(*,*)'direction',direction,dir(1),dir(2),dir(3)
            write(*,*)'M and N',M(1),M(2),N(1),N(2)
            write(*,*)'PN, PP_pro',PN(1),PN(2),PP_pro(1),PP_pro(2)
            write(*,*)v(1,1),v(1,2),v(1,3)
            write(*,*)v(2,1),v(2,2),v(2,3)
            write(*,*)v(3,1),v(3,2),v(3,3)
            write(*,*)v(4,1),v(4,2),v(4,3)
            write(*,*)myid,p1(1),p1(2),p2(1),p2(2),p3(1),p3(2),p4(1),p4(2)

        end if

        if (ival==2)write(*,*)'INTER: identical'

        !intersection with  line34 of MN line
        call lines_exp_int_2d ( p3, p4, M, N, ival, PB )
        if (ival==0)write(*,*)'INTER: no intersection'
        if (ival==2)write(*,*)'INTER: identical'

        !find the coefficent
        dA1 = sqrt((PA(1)-p1(1))**2. + (PA(2)-p1(2))**2.)
        dA2 = sqrt((PA(1)-p2(1))**2. + (PA(2)-p2(2))**2.)

        dB3 = sqrt((PB(1)-p3(1))**2. + (PB(2)-p3(2))**2.)
        dB4 = sqrt((PB(1)-p4(1))**2. + (PB(2)-p4(2))**2.)

        dPA = sqrt((PP_pro(1)-PA(1))**2. + (PP_pro(2)-PA(2))**2.)
        dPB = sqrt((PP_pro(1)-PB(1))**2. + (PP_pro(2)-PB(2))**2.)

        trilinear_coef(1) = ( dA2/(dA1+dA2) )*( dPB/(dPA+dPB) )
        trilinear_coef(2) = ( dA1/(dA1+dA2) )*( dPB/(dPA+dPB) )
        trilinear_coef(3) = ( dB4/(dB3+dB4) )*( dPA/(dPA+dPB) )
        trilinear_coef(4) = ( dB3/(dB3+dB4) )*( dPA/(dPA+dPB) )

        if (icheck_print==1) then
            write(3700,*)(p1(i),i=1,2)
            write(3700,*)(p2(i),i=1,2)
            write(3700,*)(p3(i),i=1,2)
            write(3700,*)(p4(i),i=1,2)

            write(3701,*)(PA(i),i=1,2)
            write(3701,*)(PB(i),i=1,2)

            write(3702,*)(PP_pro(i),i=1,2)

            write(3703,*)(M(i),i=1,2)
            write(3703,*)(N(i),i=1,2)


            write(3704,*)(v(1,i),i=1,3)
            write(3704,*)(v(2,i),i=1,3)
            write(3704,*)(v(3,i),i=1,3)
            write(3704,*)(v(4,i),i=1,3)
            write(3705,*)(PP(i),i=1,3)

            write(3705,*)(dir(i),i=1,3),'val',val
        end if

        !-----------------------------------------------------------------------
        ! check that the sum must be equal to 1
        sum = 0.
        do i=1,4
            sum = sum+trilinear_coef(i)
        end do

        if ( abs(1.-sum) > 1.d-10)write(*,*)'problem sum coef trilinear',l

        return
    end subroutine find_trilinear_coef2

    subroutine pp_on_meshbox(direction,ibase,jbase,kbase,ib_position,ip_position, &
        pp_position,l,trovato,trind,F,G,H,x0,y0,z0,num_ib)

        implicit none
        !-----------------------------------------------------------------------------
        integer,intent(in) :: direction
        integer,intent(in) :: ibase,jbase,kbase
        integer,intent(in) :: l,num_ib
        real,intent(in) :: ib_position(3)
        real,intent(in) :: ip_position(3)

        real,intent(out) :: F,G,H,x0,y0,z0
        logical,intent(out) :: trovato
        real,intent(out) :: pp_position(3)
        integer,intent(out) :: trind(num_ib,4,3)

        real,parameter :: tolC=1.d-12

        integer :: nt,np
        integer :: il(4),jl(4),kl(4)
        integer :: point,inout
        integer :: coincide,intersezione_faccio
        integer :: intersezioni,intersezioni_lati
        integer :: triangle_ok(8)
        real,dimension(3) :: p1,p2,p3

        logical :: intersect

        ! pp_position (local)
        real :: xpp,ypp,zpp
        real :: v(3,3),nrml(3),vert(3,3),pn(3),dist,tt
        real :: A,B,C,E

        !-----------------------------------------------------------------------------
        trovato=.false.
        intersezioni=0
        intersezioni_lati=0
        triangle_ok=0
        nt=1

        do while (nt<=8 .and. intersezioni==0 .and. intersezioni_lati==0)
            ! do nt = 1,8
            do np = 1,3
                point = index_triangle(nt,np)

                if (direction==1) then
                    il(np)=ibase
                    jl(np)=jbase+index_points(point,2)
                    kl(np)=kbase+index_points(point,1)
                else if (direction ==2) then
                    il(np)=ibase+index_points(point,1)
                    jl(np)=jbase
                    kl(np)=kbase+index_points(point,2)
                else if (direction==3) then
                    il(np)=ibase+index_points(point,1)
                    jl(np)=jbase+index_points(point,2)
                    kl(np)=kbase
                end if

                v(np,1)=xcd(il(np),jl(np),kl(np))
                v(np,2)=ycd(il(np),jl(np),kl(np))
                v(np,3)=zcd(il(np),jl(np),kl(np))

                vert(1,np)=xcd(il(np),jl(np),kl(np))
                vert(2,np)=ycd(il(np),jl(np),kl(np))
                vert(3,np)=zcd(il(np),jl(np),kl(np))

            end do

            p1(:)=v(1,:)
            p2(:)=v(2,:)
            p3(:)=v(3,:)

            ! plane for the 3 vertex
            call plane_3points(p1,p2,p3,A,B,C,E)
!            call piano_per3punti(v(1,1),v(1,2),v(1,3), &
!                v(2,1),v(2,2),v(2,3), &
!                v(3,1),v(3,2),v(3,3), &
!                A,B,C,E )

            coincide=0
            intersezione_faccio=0
            xpp=0.0
            ypp=0.0
            zpp=0.0

            call plane_imp_line_par_int_3d(A,B,C,E,x0,y0,z0,F,G,H,intersect, &
                xpp,ypp,zpp,coincide,intersezione_faccio)

            if (intersezione_faccio==1) then

                call line_exp_point_near_3d_bis(ip_position(1),ip_position(2),ip_position(3), &
                    ib_position(1),ib_position(2),ib_position(3),xpp,ypp,zpp,pn,dist,tt)

                if (tt>1.0) then

                    ! normal to the plane
                    call plane_exp_normal_3d(vert,nrml)

                    ! determine if q is inside the triangle
                    call ptpolg(3,vert,xpp,ypp,zpp,nrml,tolC,inout)

                    ! 1 inside the triangle, 0 on the border
                    if (inout>=0) then

                        if (inout==1) intersezioni=intersezioni+1
                        if (inout==0) intersezioni_lati=intersezioni_lati+1

                        triangle_ok(nt)=1
                        pp_position(1)=xpp
                        pp_position(2)=ypp
                        pp_position(3)=zpp

                    end if
                end if

            end if
            nt = nt + 1
        end do

        if (intersezioni==1 .or. intersezioni_lati==1) then
            trovato=.true.
        end if

        if (trovato) then
            do nt=1,8
                if (triangle_ok(nt)==1) then
                    do np=1,4

                        point=quadrilatero(nt,np)

                        if (direction==1) then
                            il(np)=ibase
                            jl(np)=jbase+index_points(point,2)
                            kl(np)=kbase+index_points(point,1)
                        else if (direction==2) then
                            il(np)=ibase+index_points(point,1)
                            jl(np)=jbase
                            kl(np)=kbase+index_points(point,2)
                        else if (direction==3) then
                            il(np)=ibase+index_points(point,1)
                            jl(np)=jbase+index_points(point,2)
                            kl(np)=kbase
                        end if

                        trind(l,np,1)=il(np)
                        trind(l,np,2)=jl(np)
                        trind(l,np,3)=kl(np)

                    end do
                end if
            end do
        end if

        return
    end subroutine pp_on_meshbox

    subroutine line_exp_point_near_3d_bis(xip,yip,zip,sx,sy,sz,xpp,ypp,zpp,pn,dist,t)

        !*****************************************************************************80
        !
        !! LINE_EXP_POINT_NEAR_3D: nearest point on explicit line to point in 3D.
        !
        !  Discussion:
        !
        !    The explicit form of a line in 3D is:
        !
        !      the line through the points P1 and P2.
        !
        !    The nearest point PN has the form:
        !
        !      PN = ( 1 - T ) * P1 + T * P2.
        !
        !    If T is less than 0, PN is furthest away from P2.
        !    If T is between 0 and 1, PN is between P1 and P2.
        !    If T is greater than 1, PN is furthest away from P1.
        !
        !  Modified:
        !
        !    06 May 2005
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, real :: P1(3), P2(3), two points on the line.
        !
        !    Input, real :: P(3), the point whose nearest neighbor on
        !    the line is to be determined.
        !
        !    Output, real :: PN(3), the point which is the nearest
        !    point on the line to P.
        !
        !    Output, real :: DIST, the distance from the point to the
        !    nearest point on the line.
        !
        !    Output, real :: T, the relative position of the point
        !    PN to P1 and P2.
        !
        implicit none

        integer, parameter :: dim_num = 3

        real :: bot
        real :: dist
        !logical line_exp_is_degenerate_nd
        real :: p(dim_num)
        real :: p1(dim_num)
        real :: p2(dim_num)
        real :: pn(dim_num)
        real :: t
        real :: xip,yip,zip
        real :: sx,sy,sz
        real :: xpp,ypp,zpp

        p1(1) = xip
        p1(2) = yip
        p1(3) = zip

        p2(1) = sx
        p2(2) = sy
        p2(3) = sz

        p(1) = xpp
        p(2) = ypp
        p(3) = zpp



        if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'LINE_EXP_POINT_NEAR_3D - Fatal error!'
            write ( *, '(a)' ) '  The line is degenerate.'
            stop
        end if
        !
        !  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
        !
        !  (P-P1) dot (P2-P1) / Norm(P-P1)**2 = normalized coordinate T
        !  of the projection of (P-P1) onto (P2-P1).
        !
        bot = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )

        t = sum ( ( p(1:dim_num) - p1(1:dim_num) ) &
            * ( p2(1:dim_num) - p1(1:dim_num) ) ) / bot
        !
        !  Now compute the location of the projection point.
        !
        pn(1:dim_num) = p1(1:dim_num) + t * ( p2(1:dim_num) - p1(1:dim_num) )
        !
        !  Now compute the distance between the projection point and P.
        !
        dist = sqrt ( sum ( ( pn(1:dim_num) - p(1:dim_num) )**2 ) )

        return
    end

    subroutine triangles_mesh()
        !***********************************************************************

        implicit none
        !-----------------------------------------------------------------------
        index_triangle(1,1) = 1
        index_triangle(1,2) = 2
        index_triangle(1,3) = 5

        index_triangle(2,1) = 1
        index_triangle(2,2) = 5
        index_triangle(2,3) = 4

        index_triangle(3,1) = 2
        index_triangle(3,2) = 3
        index_triangle(3,3) = 6

        index_triangle(4,1) = 2
        index_triangle(4,2) = 6
        index_triangle(4,3) = 5

        index_triangle(5,1) = 4
        index_triangle(5,2) = 5
        index_triangle(5,3) = 8

        index_triangle(6,1) = 4
        index_triangle(6,2) = 8
        index_triangle(6,3) = 7

        index_triangle(7,1) = 5
        index_triangle(7,2) = 6
        index_triangle(7,3) = 9

        index_triangle(8,1) = 5
        index_triangle(8,2) = 9
        index_triangle(8,3) = 8

        index_points(1,1) = -1
        index_points(1,2) = -1

        index_points(2,1) =  0
        index_points(2,2) = -1

        index_points(3,1) =  1
        index_points(3,2) = -1

        index_points(4,1) = -1
        index_points(4,2) =  0

        index_points(5,1) =  0
        index_points(5,2) =  0

        index_points(6,1) =  1
        index_points(6,2) =  0

        index_points(7,1) = -1
        index_points(7,2) =  1

        index_points(8,1) =  0
        index_points(8,2) =  1

        index_points(9,1) =  1
        index_points(9,2) =  1


        quadrilatero(1,1) = 1
        quadrilatero(1,2) = 2
        quadrilatero(1,3) = 5
        quadrilatero(1,4) = 4

        quadrilatero(2,1) = 1
        quadrilatero(2,2) = 2
        quadrilatero(2,3) = 5
        quadrilatero(2,4) = 4

        quadrilatero(3,1) = 2
        quadrilatero(3,2) = 3
        quadrilatero(3,3) = 6
        quadrilatero(3,4) = 5

        quadrilatero(4,1) = 2
        quadrilatero(4,2) = 3
        quadrilatero(4,3) = 6
        quadrilatero(4,4) = 5

        quadrilatero(5,1) = 4
        quadrilatero(5,2) = 5
        quadrilatero(5,3) = 8
        quadrilatero(5,4) = 7

        quadrilatero(6,1) = 4
        quadrilatero(6,2) = 5
        quadrilatero(6,3) = 8
        quadrilatero(6,4) = 7

        quadrilatero(7,1) = 5
        quadrilatero(7,2) = 6
        quadrilatero(7,3) = 9
        quadrilatero(7,4) = 8

        quadrilatero(8,1) = 5
        quadrilatero(8,2) = 6
        quadrilatero(8,3) = 9
        quadrilatero(8,4) = 8
        !-----------------------------------------------------------------------
        return
    end

    subroutine line_exp_point_near_2d(p1,p2,p,pn,dist,t)

        !*****************************************************************************80
        !
        !! LINE_EXP_POINT_NEAR_2D computes the point on an explicit line nearest a point in 2D.
        !
        !  Discussion:
        !
        !    The explicit form of a line in 2D is:
        !
        !      the line through the points P1 and P2.
        !
        !    The nearest point PN = (XN,YN) has the form:
        !
        !      PN = (1-T) * P1 + T * P2.
        !
        !    If T is less than 0, PN is furthest from P2.
        !    If T is between 0 and 1, PN is between P1 and P2.
        !    If T is greater than 1, PN is furthest from P1.
        !
        !  Modified:
        !
        !    06 May 2005
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, real :: P1(2), P2(2), two points on the line.
        !
        !    Input, real :: P(2), the point whose nearest neighbor on the
        !    line is to be determined.
        !
        !    Output, real :: PN(2), the nearest point on the line to P.
        !
        !    Output, real :: DIST, the distance from the point to the line.
        !
        !    Output, real :: T, the relative position of the point
        !    PN to the points P1 and P2.
        !
        implicit none

        integer, parameter :: dim_num = 2

        real :: bot
        real :: dist
        !logical line_exp_is_degenerate_nd
        real :: p(dim_num)
        real :: p1(dim_num)
        real :: p2(dim_num)
        real :: pn(dim_num)
        real :: t

        if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'LINE_POINT_NEAR_2D - Fatal error!'
            write ( *, '(a)' ) '  The line is degenerate.'
            stop
        end if
        !
        !  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
        !
        !  (P-P1) dot (P2-P1) / Norm(P-P1)**2 = normalized coordinate T
        !  of the projection of (P-P1) onto (P2-P1).
        !
        bot = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )

        t = sum ( ( p1(1:dim_num) - p(1:dim_num) ) &
            * ( p1(1:dim_num) - p2(1:dim_num) ) ) / bot

        pn(1:dim_num) = p1(1:dim_num) + t * ( p2(1:dim_num) - p1(1:dim_num) )

        dist = sqrt ( sum ( ( pn(1:dim_num) - p(1:dim_num) )**2 ) )

        return
    end subroutine line_exp_point_near_2d

    subroutine lines_exp_int_2d(p1,p2,q1,q2,ival,p)

        !*****************************************************************************80
        !
        !! LINES_EXP_INT_2D determines where two explicit lines intersect in 2D.
        !
        !  Discussion:
        !
        !    The explicit form of a line in 2D is:
        !
        !      the line through the points P1 and P2.
        !
        !  Modified:
        !
        !    02 January 2005
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, real :: P1(2), P2(2), two points on the first line.
        !
        !    Input, real :: Q1(2), Q2(2), two points on the second line.
        !
        !    Output, integer IVAL, reports on the intersection:
        !    0, no intersection, the lines may be parallel or degenerate.
        !    1, one intersection point, returned in P.
        !    2, infinitely many intersections, the lines are identical.
        !
        !    Output, real :: P(2), if IVAl = 1, P is
        !    the intersection point.  Otherwise, P = 0.
        !
        implicit none

        integer, parameter :: dim_num = 2

        real :: a1
        real :: a2
        real :: b1
        real :: b2
        real :: c1
        real :: c2
        integer ival
        logical point_1
        logical point_2
        real :: p(dim_num)
        real :: p1(dim_num)
        real :: p2(dim_num)
        real :: q1(dim_num)
        real :: q2(dim_num)

        ival = 0
        p(1:dim_num) = 0.0D+00
        !
        !  Check whether either line is a point.
        !
        if ( all ( p1(1:dim_num)==p2(1:dim_num) ) ) then
            point_1 = .true.
        else
            point_1 = .false.
        end if

        if ( all ( q1(1:dim_num)==q2(1:dim_num) ) ) then
            point_2 = .true.
        else
            point_2 = .false.
        end if
        !
        !  Convert the lines to ABC format.
        !
        if ( .not. point_1 ) then
            call line_exp2imp_2d ( p1, p2, a1, b1, c1 )
        end if

        if ( .not. point_2 ) then
            call line_exp2imp_2d ( q1, q2, a2, b2, c2 )
        end if
        !
        !  Search for intersection of the lines.
        !
        if ( point_1 .and. point_2 ) then
            if ( all ( p1(1:dim_num)==q1(1:dim_num) ) ) then
                ival = 1
                p(1:dim_num) = p1(1:dim_num)
            end if
        else if ( point_1 ) then
            if ( a2 * p1(1) + b2 * p1(2)==c2 ) then
                ival = 1
                p(1:dim_num) = p1(1:dim_num)
            end if
        else if ( point_2 ) then
            if ( a1 * q1(1) + b1 * q1(2)==c1 ) then
                ival = 1
                p(1:dim_num) = q1(1:dim_num)
            end if
        else
            call lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, p )
        end if

        return
    end subroutine lines_exp_int_2d

    subroutine line_exp2imp_2d(p1,p2,a,b,c)

        !*****************************************************************************80
        !
        !! LINE_EXP2IMP_2D converts an explicit line to implicit form in 2D.
        !
        !  Discussion:
        !
        !    The explicit form of a line in 2D is:
        !
        !      the line through the points P1 and P2.
        !
        !    The implicit form of a line in 2D is:
        !
        !      A * X + B * Y + C = 0
        !
        !  Modified:
        !
        !    06 May 2005
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, real :: P1(2), P2(2), two points on the line.
        !
        !    Output, real :: A, B, C, the implicit form of the line.
        !
        implicit none

        integer, parameter :: dim_num = 2

        real :: a
        real :: b
        real :: c
        !logical line_exp_is_degenerate_nd
        real :: norm
        real :: p1(dim_num)
        real :: p2(dim_num)
        !
        !  Take care of degenerate cases.
        !
        if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'LINE_EXP2IMP_2D - Warning!'
            write ( *, '(a)' ) '  The line is degenerate.'
        end if

        a = p2(2) - p1(2)
        b = p1(1) - p2(1)
        c = p2(1) * p1(2) - p1(1) * p2(2)

        norm = a * a + b * b + c * c

        if ( 0.0D+00 < norm ) then
            a = a / norm
            b = b / norm
            c = c / norm
        end if

        if ( a < 0.0D+00 ) then
            a = -a
            b = -b
            c = -c
        end if

        return
    end subroutine line_exp2imp_2d

    subroutine lines_imp_int_2d(a1,b1,c1,a2,b2,c2,ival,p)

        !*****************************************************************************80
        !
        !! LINES_IMP_INT_2D determines where two implicit lines intersect in 2D.
        !
        !  Discussion:
        !
        !    The implicit form of a line in 2D is:
        !
        !      A * X + B * Y + C = 0
        !
        !  Modified:
        !
        !    25 February 2005
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, real :: A1, B1, C1, define the first line.
        !    At least one of A1 and B1 must be nonzero.
        !
        !    Input, real :: A2, B2, C2, define the second line.
        !    At least one of A2 and B2 must be nonzero.
        !
        !    Output, integer IVAL, reports on the intersection.
        !
        !    -1, both A1 and B1 were zero.
        !    -2, both A2 and B2 were zero.
        !     0, no intersection, the lines are parallel.
        !     1, one intersection point, returned in P.
        !     2, infinitely many intersections, the lines are identical.
        !
        !    Output, real :: P(2), if IVAL = 1, then P is
        !    the intersection point.  Otherwise, P = 0.
        !
        implicit none

        integer, parameter :: dim_num = 2

        real :: a(dim_num,dim_num+1)
        real :: a1
        real :: a2
        real :: b1
        real :: b2
        real :: c1
        real :: c2
        integer info
        integer ival
        !logical line_imp_is_degenerate_2d
        real :: p(dim_num)

        p(1:dim_num) = 0.0D+00
        !
        !  Refuse to handle degenerate lines.
        !
        if ( line_imp_is_degenerate_2d ( a1, b1, c1 ) ) then
            ival = -1
            return
        end if

        if ( line_imp_is_degenerate_2d ( a2, b2, c2 ) ) then
            ival = -2
            return
        end if
        !
        !  Set up and solve a linear system.
        !
        a(1,1) = a1
        a(1,2) = b1
        a(1,3) = -c1

        a(2,1) = a2
        a(2,2) = b2
        a(2,3) = -c2

        call r8mat_solve ( 2, 1, a, info )
        !
        !  If the inverse exists, then the lines intersect at the solution point.
        !
        if ( info==0 ) then

            ival = 1
            p(1:dim_num) = a(1:dim_num,3)
        !
        !  If the inverse does not exist, then the lines are parallel
        !  or coincident.  Check for parallelism by seeing if the
        !  C entries are in the same ratio as the A or B entries.
        !
        else

            ival = 0

            if ( a1==0.0D+00 ) then
                if ( b2 * c1==c2 * b1 ) then
                    ival = 2
                end if
            else
                if ( a2 * c1==c2 * a1 ) then
                    ival = 2
                end if
            end if

        end if

        return
    end subroutine lines_imp_int_2d

    subroutine r8mat_solve(n,rhs_num,a,info)

        !*****************************************************************************80
        !
        !! R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
        !
        !  Modified:
        !
        !    29 August 2003
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, integer N, the order of the matrix.
        !
        !    Input, integer rhs_num, the number of right hand sides.  rhs_num
        !    must be at least 0.
        !
        !    Input/output, real :: A(N,N+rhs_num), contains in rows and
        !    columns 1 to N the coefficient matrix, and in columns N+1 through
        !    N+rhs_num, the right hand sides.  On output, the coefficient matrix
        !    area has been destroyed, while the right hand sides have
        !    been overwritten with the corresponding solutions.
        !
        !    Output, integer INFO, singularity flag.
        !    0, the matrix was not singular, the solutions were computed;
        !    J, factorization failed on step J, and the solutions could not
        !    be computed.
        !
        implicit none

        integer n
        integer rhs_num

        real :: a(n,n+rhs_num)
        real :: apivot
        real :: factor
        integer i
        integer info
        integer ipivot
        integer j

        info = 0

        do j = 1, n
            !
            !  Choose a pivot row.
            !
            ipivot = j
            apivot = a(j,j)

            do i = j+1, n
                if ( abs ( apivot ) < abs ( a(i,j) ) ) then
                    apivot = a(i,j)
                    ipivot = i
                end if
            end do

            if ( apivot==0.0D+00 ) then
                info = j
                return
            end if
            !
            !  Interchange.
            !
            do i = 1, n + rhs_num
                call r8_swap ( a(ipivot,i), a(j,i) )
            end do
            !
            !  A(J,J) becomes 1.
            !
            a(j,j) = 1.0D+00
            a(j,j+1:n+rhs_num) = a(j,j+1:n+rhs_num) / apivot
            !
            !  A(I,J) becomes 0.
            !
            do i = 1, n

                if ( i /= j ) then

                    factor = a(i,j)
                    a(i,j) = 0.0D+00
                    a(i,j+1:n+rhs_num) = a(i,j+1:n+rhs_num) - factor * a(j,j+1:n+rhs_num)

                end if

            end do

        end do

        return
    end subroutine r8mat_solve

    subroutine r8_swap(x,y)

        !*****************************************************************************80
        !
        !! R8_SWAP switches two R8's.
        !
        !  Modified:
        !
        !    01 May 2000
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input/output, real :: X, Y.  On output, the values of X and
        !    Y have been interchanged.
        !
        implicit none

        real :: x
        real :: y
        real :: z

        z = x
        x = y
        y = z

        return
    end subroutine r8_swap

    function line_imp_is_degenerate_2d ( a, b, c )

        !*****************************************************************************80
        !
        !! LINE_IMP_IS_DEGENERATE_2D finds if an implicit point is degenerate in 2D.
        !
        !  Discussion:
        !
        !    The implicit form of a line in 2D is:
        !
        !      A * X + B * Y + C = 0
        !
        !  Modified:
        !
        !    06 May 2005
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, real :: A, B, C, the implicit line parameters.
        !
        !    Output, logical LINE_IMP_IS_DEGENERATE_2D, is true if the
        !    line is degenerate.
        !
        implicit none

        integer, parameter :: dim_num = 2

        real :: a
        real :: b
        real :: c
        logical line_imp_is_degenerate_2d

        line_imp_is_degenerate_2d=(a*a+b*b==0.0D+00)

        return
    end function line_imp_is_degenerate_2d

end module trilinear
