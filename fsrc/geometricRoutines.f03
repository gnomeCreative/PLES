module geometricRoutines

    implicit none

contains

    subroutine plane_imp_line_par_int_3d(a,b,c,d,x0,y0,z0,f,g,h,intersect,x,y,z,coincide,intersezione)
        !

        !
        !! PLANE_IMP_LINE_PAR_INT_3D: intersection ( implicit plane, parametric line ) in 3D.
        !
        !
        !  Definition:
        !
        !    The implicit form of a plane in 3D is:
        !
        !      A * X + B * Y + C * Z + D = 0
        !
        !  Reference:
        !
        !    Adrian Bowyer and John Woodwark,
        !    A Programmer's Geometry,
        !    Butterworths, 1983, page 111.
        !
        !  Modified:
        !
        !    16 April 1999
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, real A, B, C, D, parameters that define the implicit plane.
        !
        !    Input, real X0, Y0, Z0, F, G, H, parameters that define the
        !    parametric line.
        !
        !    Output, logical INTERSECT, is TRUE if the line and the plane
        !    intersect, and false otherwise.
        !
        !    Output, real X, Y, Z, is a point of intersection of the line
        !    and the plane, if INTERSECT is TRUE.
        !
        implicit none
        !
        real :: a
        real :: b
        real :: c
        real :: d
        real :: denom
        real :: f
        real :: g
        real :: h
        logical :: intersect
        real :: norm1
        real :: norm2
        real :: t
        real, parameter :: tol = 0.00001E+00
        real :: x
        real :: x0
        real :: y
        real :: y0
        real :: z
        real :: z0
        integer :: coincide,intersezione
        !
        !  Check.
        !
        norm1 = sqrt ( a * a + b * b + c * c )

        if ( norm1 == 0.0E+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PLANE_IMP_LINE_PAR_INT_3D - Fatal error!'
            write ( *, '(a)' ) '  The plane normal vector is null.'
            stop
        end if

        norm2 = sqrt ( f * f + g * g + h * h )

        if ( norm2 == 0.0E+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PLANE_IMP_LINE_PAR_INT_3D - Fatal error!'
            write ( *, '(a)' ) '  The line direction vector is null.'
            stop
        end if

        denom = a * f + b * g + c * h
        !
        !  The line and the plane may be parallel.
        !
        if ( abs ( denom ) < TOL * norm1 * norm2 ) then

            if ( a * x0 + b * y0 + c * z0 + d == 0.0E+00 ) then
                intersect = .TRUE.
                intersezione = 1
                x = x0
                y = y0
                z = z0
                coincide = 1
            else
                intersect = .FALSE.
                intersezione = 0
                x = 0.0E+00
                y = 0.0E+00
                z = 0.0E+00
            end if
        !
        !  If they are not parallel, they must intersect.
        !
        else

            intersect = .TRUE.
            intersezione = 1
            t = - ( a * x0 + b * y0 + c * z0 + d ) / denom
            x = x0 + t * f
            y = y0 + t * g
            z = z0 + t * h

        end if

        return
    end subroutine plane_imp_line_par_int_3d

    subroutine plane_exp_normal_3d (vertice,normal)
        !! PLANE_EXP_NORMAL_3D finds the normal to an explicit plane in 3D.
        !
        !  Discussion:
        !    The explicit form of a plane in 3D is
        !      the plane through P1, P2 and P3.
        !  Modified:
        !    10 February 2005
        !  Author:
        !    John Burkardt
        !  Parameters:
        !
        !    Input, real ( kind = 8 ) P1(3), P2(3), P3(3), three points on the plane.
        !
        !    Output, real ( kind = 8 ) NORMAL(3), the coordinates of the unit normal
        !    vector to the plane containing the three points.
        !
        implicit none

        integer, parameter :: dim_num = 3
        integer :: i

        real ( kind = 8 ) normal(dim_num)
        real ( kind = 8 ) normal_norm
        real ( kind = 8 ) p1(dim_num)
        real ( kind = 8 ) p2(dim_num)
        real ( kind = 8 ) p3(dim_num)
        real ( kind = 8 ) vertice(3,3)

        do i=1,3
            p1(i)=vertice(i,1)
            p2(i)=vertice(i,2)
            p3(i)=vertice(i,3)
        end do
        !
        !  The cross product (P2-P1) x (P3-P1) is normal to (P2-P1) and (P3-P1).
        !
        normal(1) = ( p2(2) - p1(2) ) * ( p3(3) - p1(3) ) &
            - ( p2(3) - p1(3) ) * ( p3(2) - p1(2) )

        normal(2) = ( p2(3) - p1(3) ) * ( p3(1) - p1(1) ) &
            - ( p2(1) - p1(1) ) * ( p3(3) - p1(3) )

        normal(3) = ( p2(1) - p1(1) ) * ( p3(2) - p1(2) ) &
            - ( p2(2) - p1(2) ) * ( p3(1) - p1(1) )

        normal_norm = sqrt ( sum ( normal(1:dim_num)**2 ) )

        if ( normal_norm == 0.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PLANE_EXP_NORMAL_3D - Fatal error!'
            write ( *, '(a)' ) '  The plane is poorly defined.'
            stop
        end if

        normal(1:dim_num) = normal(1:dim_num) / normal_norm

        return
    end subroutine plane_exp_normal_3d

    subroutine ptpolg(dim,vcl,qx,qy,qz,nrml,dtol,inout)

        !*****************************************************************************80
        !
        !! PTPOLG determines if a point is in, on or outside a polygon.
        !
        !  Purpose:
        !
        !    Determine whether a point lies inside, outside, or on
        !    boundary of a planar polygon in 2 or 3 dimensional space.
        !    It is assumed that point lies in plane of polygon.
        !
        !  Modified:
        !
        !    12 July 1999
        !
        !  Author:
        !
        !    Barry Joe,
        !    Department of Computing Science,
        !    University of Alberta,
        !    Edmonton, Alberta, Canada  T6G 2H1
        !
        !  Reference:
        !
        !    Barry Joe,
        !    GEOMPACK - a software package for the generation of meshes
        !    using geometric algorithms,
        !    Advances in Engineering Software,
        !    Volume 13, pages 325-331, 1991.
        !
        !  Parameters:
        !
        !    Input, integer DIM, the dimension of the polygon (2 or 3).
        !
        !    Input, integer LDV, the leading dimension of VCL array in calling routine.
        !
        !    Input, integer NV, the number of vertices in polygon.
        !
        !    Input, integer INC, the increment for PGIND indicating indices of polygon.
        !
        !    Input, integer PGIND(0:NV*INC), indices in VCL of polygon vertices are in
        !    PGIND(0), PGIND(INC), ..., PGIND(NV*INC) with first and
        !    last vertices identical.
        !
        !    Input, real ( kind = 8 ) VCL(1:DIM,1:*), the vertex coordinate list.
        !
        !    Input, real ( kind = 8 ) PT(1:DIM), the point for which in/out test is
        !    applied.
        !
        !    Input, real ( kind = 8 ) NRML(1:3), the unit normal vector of plane
        !    containing polygon, with vertices oriented counter clockwise with
        !    respect to the normal (used iff DIM = 3);
        !    The normal is assumed to be (0,0,1) if DIM = 2.
        !
        !    Input, real ( kind = 8 ) DTOL, an absolute tolerance to determine
        !    whether a point is on a line or plane.
        !
        !    Output, integer INOUT, point PT is:
        !    +1, inside the polygon,
        !     0, on boundary of polygon,
        !    -1, outside polygon;
        !    -2 if error in input parameters.
        !
        integer dim
        integer ldv

        real ( kind = 8 ) cp(3)
        real ( kind = 8 ) de(3)
        real ( kind = 8 ) dir(3)
        real ( kind = 8 ) dist
        real ( kind = 8 ) dotp
        real ( kind = 8 ) dtol
        integer i
        integer inc
        integer inout
        integer j
        integer k
        integer l
        integer la
        integer lb
        real ( kind = 8 ) len1
        real ( kind = 8 ) len2
        integer m
        integer n
        real ( kind = 8 ) nr(4)
        real ( kind = 8 ) nrml(3)
        integer nv
        integer pgind(0:3)
        real ( kind = 8 ) pt(dim)
        real ( kind = 8 ) rhs(3)
        integer s
        integer sa
        integer sb
        real ( kind = 8 ) t
        real ( kind = 8 ) ta
        real ( kind = 8 ) tol
        real ( kind = 8 ) vcl(3,3)
        real ( kind = 8 ) qx,qy,qz

        !---------------------------------------------
        pgind(0)=1
        pgind(1)=2
        pgind(2)=3
        pgind(3)=1
        nv = 3    !umero di vertici per poligono
        inc = 1  !increment for pgind
        ldv  = 3  !leading dimension of vertices

        pt(1)=qx  !punto da testare
        pt(2)=qy
        pt(3)=qz
        !---------------------------------------------

        tol = 100.0D+00 * epsilon ( tol )

        inout = -2

        if ( dim < 2 .or. 3 < dim ) then
            return
        end if
        !
        !  Direction of ray is from PT through midpoint of first edge
        !  such that PT is not collinear with edge. NR is normal of plane
        !  containing ray, which is also orthogonal to NRML.
        !
        i = 0
        lb = pgind(0)

10  continue

    i = i + 1

    if ( nv <= i ) then
        return
    end if

    la = lb
    lb = pgind(i*inc)

    do j = 1, dim
        de(j) = vcl(j,lb) - vcl(j,la)
        dir(j) = pt(j) - vcl(j,la)
    end do

    if ( dim == 2 ) then
        len1 = de(1)**2 + de(2)**2
        len2 = dir(1)**2 + dir(2)**2
    else
        len1 = de(1)**2 + de(2)**2 + de(3)**2
        len2 = dir(1)**2 + dir(2)**2 + dir(3)**2
    end if

    if ( len1 == 0.0D+00 ) then
        goto 10
    else if ( len2 == 0.0D+00 ) then
        inout = 0
        return
    end if

    if ( dim == 2 ) then
        dotp = abs ( de(1) * dir(1) + de(2) * dir(2)) / sqrt(len1*len2)
    else if ( dim == 3 ) then
        dotp = abs ( de(1) * dir(1) + de(2) * dir(2) + de(3) * dir(3) ) &
            / sqrt(len1*len2)
    end if

    if ( 1.0D+00 - tol <= dotp ) then
        goto 10
    end if

    if ( dim == 2 ) then
        dir(1) = 0.5E+00 * ( vcl(1,la) + vcl(1,lb) ) - pt(1)
        dir(2) = 0.5E+00 * ( vcl(2,la) + vcl(2,lb) ) - pt(2)
        dist = sqrt ( dir(1)**2 + dir(2)**2 )
        dir(1) = dir(1) / dist
        dir(2) = dir(2) / dist
        dir(3) = 0.0D+00
        nr(1) = -dir(2)
        nr(2) = dir(1)
        nr(3) = 0.0D+00
        nr(4) = nr(1) * pt(1) + nr(2) * pt(2)
        dist = nr(1) * vcl(1,lb) + nr(2) * vcl(2,lb) - nr(4)
    else if ( dim == 3 ) then
        dir(1) = 0.5E+00 * ( vcl(1,la) + vcl(1,lb) ) - pt(1)
        dir(2) = 0.5E+00 * ( vcl(2,la) + vcl(2,lb) ) - pt(2)
        dir(3) = 0.5E+00 * ( vcl(3,la) + vcl(3,lb) ) - pt(3)
        dist = sqrt ( dir(1)**2 + dir(2)**2 + dir(3)**2 )

        dir(1) = dir(1) / dist
        dir(2) = dir(2) / dist
        dir(3) = dir(3) / dist

        nr(1) = nrml(2)*dir(3) - nrml(3)*dir(2)
        nr(2) = nrml(3)*dir(1) - nrml(1)*dir(3)
        nr(3) = nrml(1)*dir(2) - nrml(2)*dir(1)
        nr(4) = nr(1)*pt(1) + nr(2)*pt(2) + nr(3)*pt(3)
        dist = nr(1)*vcl(1,lb)+nr(2)*vcl(2,lb)+nr(3)*vcl(3,lb) - nr(4)
    end if

    if ( 0.0D+00 < dist ) then
        sb = 1
    else
        sb = -1
    end if

    m = 1
    if ( abs ( dir(1) ) < abs ( dir(2) ) ) then
        m = 2
    end if

    if ( abs ( dir(m) ) < abs ( dir(3) ) ) then
        m = 3
    end if

    k = 1
    !
    !  For remaining edges of polygon, check whether ray intersects edge.
    !  Vertices or edges lying on ray are handled by looking at preceding
    !  and succeeding edges not lying on ray.
    !
    n = i
    i = i + 1

30 continue

   la = lb
   lb = pgind(i*inc)
   sa = sb

   if ( dim == 2 ) then
       dist = nr(1) * vcl(1,lb) + nr(2)*vcl(2,lb) - nr(4)
   else if ( dim == 3 ) then
       dist = nr(1) * vcl(1,lb) + nr(2)*vcl(2,lb) + nr(3)*vcl(3,lb)- nr(4)
   end if

   if ( abs ( dist) <= dtol ) then
       sb = 0
   else if ( 0.0D+00 < dist ) then
       sb = 1
   else
       sb = -1
   end if

   s = sa * sb

   if ( s < 0 ) then

       if ( dim == 2 ) then
           de(1) = vcl(1,la) - vcl(1,lb)
           de(2) = vcl(2,la) - vcl(2,lb)
           rhs(1) = vcl(1,la) - pt(1)
           rhs(2) = vcl(2,la) - pt(2)
           t = ( rhs(1) * de(2) - rhs(2) * de(1) ) &
               / ( dir(1) * de(2) - dir(2) * de(1) )
       else if ( dim == 3 ) then
           de(1) = vcl(1,la) - vcl(1,lb)
           de(2) = vcl(2,la) - vcl(2,lb)
           de(3) = vcl(3,la) - vcl(3,lb)
           rhs(1) = vcl(1,la) - pt(1)
           rhs(2) = vcl(2,la) - pt(2)
           rhs(3) = vcl(3,la) - pt(3)
           cp(1) = dir(2) * de(3) - dir(3) * de(2)
           cp(2) = dir(3) * de(1) - dir(1) * de(3)
           cp(3) = dir(1) * de(2) - dir(2) * de(1)

           l = 1
           if ( abs ( cp(1) ) < abs ( cp(2) ) ) then
               l = 2
           end if
           if ( abs ( cp(l) ) < abs ( cp(3) ) ) then
               l = 3
           end if

           if ( l == 1 ) then
               t = ( rhs(2) * de(3) - rhs(3) * de(2) ) / cp(1)
           else if ( l == 2 ) then
               t = ( rhs(3) * de(1) - rhs(1) * de(3) ) / cp(2)
           else
               t = ( rhs(1) * de(2) - rhs(2) * de(1) ) / cp(3)
           end if

       end if

       if ( dtol < t ) then
           k = k + 1
       else if ( -dtol <= t ) then
           inout = 0
           return
       end if

   else if ( s == 0 ) then

       l = lb
40 continue
   i = i + 1

   if ( nv < i ) then
       i = 1
   end if

   if ( i == n) then
       return
   end if

   la = lb
   lb = pgind(i*inc)

   if ( dim == 2 ) then
       dist = nr(1) * vcl(1,lb) + nr(2) * vcl(2,lb) - nr(4)
   else
       dist = nr(1) * vcl(1,lb) + nr(2) * vcl(2,lb) + nr(3) * vcl(3,lb) - nr(4)
   end if

   if ( abs ( dist) <= dtol ) then
       go to 40
   else if ( 0.0D+00 < dist ) then
       sb = 1
   else
       sb = -1
   end if

   t = ( vcl(m,l) - pt(m) ) / dir(m)

   if ( abs ( t) <= dtol ) then
       inout = 0
       return
   end if

   if ( la /= l ) then
       ta = ( vcl(m,la) - pt(m) ) / dir(m)
       if ( abs ( ta ) <= dtol .or. t * ta < 0.0D+00 ) then
           inout = 0
           return
       end if
   end if

   if ( sa * sb < 0 .and. 0.0D+00 < t ) then
       k = k + 1
   end if

   end if

   i = i + 1

   if ( nv < i ) then
       i = 1
   end if

   if ( i /= n ) then
       go to 30
   end if
   !
   !  Point lies inside polygon if number of intersections K is odd.
   !
   if ( mod ( k, 2 ) == 1 ) then
       inout = 1
   else
       inout = -1
   end if

   return
   end subroutine ptpolg

   subroutine line_exp2par_3d(x1,y1,z1,x2,y2,z2,f,g,h,x0,y0,z0)
       !

       !
       !! LINE_EXP2PAR_3D converts a line from explicit to parametric form in 3D.
       !
       !
       !  Formula:
       !
       !    The explicit form of a line in 3D is:
       !
       !      (X1,Y1,Z1), (X2,Y2,Z2).
       !
       !    The parametric form of a line in 3D is:
       !
       !      X = X0 + F * T
       !      Y = Y0 + G * T
       !      Z = Z0 + H * T
       !
       !  Modified:
       !
       !    30 January 1999
       !
       !  Author:
       !
       !    John Burkardt
       !
       !  Parameters:
       !
       !    Input, real X1, Y1, Z1, X2, Y2, Z3, two points on the line.
       !
       !    Output, real F, G, H, X0, Y0, Z0, the parametric parameters of the line.
       !
       implicit none
       !
       real :: f
       real :: g
       real :: h
       real :: norm
       real :: x0
       real :: x1
       real :: x2
       real :: y0
       real :: y1
       real :: y2
       real :: z0
       real :: z1
       real :: z2
       !
       x0 = x1
       y0 = y1
       z0 = z1

       f = x2 - x1
       g = y2 - y1
       h = z2 - z1

       norm = sqrt ( f * f + g * g + h * h )

       if ( norm /= 0.0E+00 ) then
           f = f / norm
           g = g / norm
           h = h / norm
       end if

       return
   end subroutine line_exp2par_3d

   function line_exp_is_degenerate_nd(dim_num,p1,p2)

       !*****************************************************************************80
       !
       !! LINE_EXP_IS_DEGENERATE_ND finds if an explicit line is degenerate in ND.
       !
       !  Discussion:
       !
       !    The explicit form of a line in ND is:
       !
       !      the line through the points P1 and P2.
       !
       !    An explicit line is degenerate if the two defining points are equal.
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
       !    Input, integer DIM_NUM, the spatial dimension.
       !
       !    Input, real :: P1(DIM_NUM), P2(DIM_NUM), two points on the line.
       !
       !    Output, logical LINE_EXP_IS_DEGENERATE_ND, is TRUE if the line
       !    is degenerate.
       !
       implicit none

       integer dim_num

       logical line_exp_is_degenerate_nd
       real :: p1(dim_num)
       real :: p2(dim_num)

       line_exp_is_degenerate_nd = ( all ( p1(1:dim_num) == p2(1:dim_num) ) )

       return
   end function line_exp_is_degenerate_nd

   function parallelepiped_contains_point_3d(p1,p2,p3,p4,p)

       !*********************************************************************
       !
       !! PARALLELEPIPED_CONTAINS_POINT_3D: point inside parallelepiped in 3D.
       !
       !  Discussion:
       !
       !    A parallelepiped is a "slanted box", that is, opposite
       !    sides are parallel planes.
       !
       !         *------------------*
       !        / .                / \
       !       /   .              /   \
       !      /     .            /     \
       !    P4------------------*       \
       !      \        .         \       \
       !       \        .         \       \
       !        \        .         \       \
       !         \       P2.........\.......\
       !          \     .            \     /
       !           \   .              \   /
       !            \ .                \ /
       !            P1-----------------P3
       !
       !  Modified:
       !
       !    05 March 2005
       !
       !  Author:
       !
       !    John Burkardt
       !
       !  Parameters:
       !
       !    Input, real :: P1(3), P2(3), P3(3), P4(3), four corners
       !    of the parallelepiped.  It is assumed that P2, P3 and P4 are
       !    immediate neighbors of P1.
       !
       !    Input, real :: P(3), the point to be checked.
       !
       !    Output, logical PARALLELEPIPED_CONTAINS_POINT_3D, is true if P
       !    is inside the parallelepiped, or on its boundary.
       !
       implicit none

       integer, parameter :: dim_num = 3

       real :: dot
       logical parallelepiped_contains_point_3d
       real :: p(dim_num)
       real :: p1(dim_num)
       real :: p2(dim_num)
       real :: p3(dim_num)
       real :: p4(dim_num)

       parallelepiped_contains_point_3d = .false.

       dot = dot_product ( p(1:dim_num)  - p1(1:dim_num), &
           p2(1:dim_num) - p1(1:dim_num) )

       if ( dot < 0.0D+00 ) then
           return
       end if

       if ( sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 ) < dot ) then
           return
       end if

       dot = dot_product ( p(1:dim_num)  - p1(1:dim_num), &
           p3(1:dim_num) - p1(1:dim_num) )

       if ( dot < 0.0D+00 ) then
           return
       end if

       if ( sum ( ( p3(1:dim_num) - p1(1:dim_num) )**2 ) < dot ) then
           return
       end if

       dot = dot_product ( p(1:dim_num)  - p1(1:dim_num), &
           p4(1:dim_num) - p1(1:dim_num) )

       if ( dot < 0.0D+00 ) then
           return
       end if

       if ( sum ( ( p4(1:dim_num) - p1(1:dim_num) )**2 ) < dot ) then
           return
       end if

       parallelepiped_contains_point_3d = .true.

       return
   end function parallelepiped_contains_point_3d

   subroutine r8vec_cross_3d(v1,v2,v3)

       !*****************************************************************************80
       !
       !! R8VEC_CROSS_3D computes the cross product of two vectors in 3D.
       !
       !  Discussion:
       !
       !    The cross product in 3D can be regarded as the determinant of the
       !    symbolic matrix:
       !
       !          |  i  j  k |
       !      det | x1 y1 z1 |
       !          | x2 y2 z2 |
       !
       !      = ( y1 * z2 - z1 * y2 ) * i
       !      + ( z1 * x2 - x1 * z2 ) * j
       !      + ( x1 * y2 - y1 * x2 ) * k
       !
       !  Modified:
       !
       !    07 August 2005
       !
       !  Author:
       !
       !    John Burkardt
       !
       !  Parameters:
       !
       !    Input, real :: V1(3), V2(3), the two vectors.
       !
       !    Output, real :: V3(3), the cross product vector.
       !
       implicit none

       integer, parameter :: dim_num = 3

       real :: v1(dim_num)
       real :: v2(dim_num)
       real :: v3(dim_num)

       v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
       v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
       v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

       return
   end subroutine r8vec_cross_3d

   function r8vec_length(dim_num,x)

       !*****************************************************************************80
       !
       !! R8VEC_LENGTH returns the Euclidean length of a vector.
       !
       !  Modified:
       !
       !    08 August 2005
       !
       !  Author:
       !
       !    John Burkardt
       !
       !  Parameters:
       !
       !    Input, integer DIM_NUM, the spatial dimension.
       !
       !    Input, real :: X(DIM_NUM), the vector.
       !
       !    Output, real :: R8VEC_LENGTH, the Euclidean length of the vector.
       !
       implicit none

       integer dim_num

       real :: r8vec_length
       real :: x(dim_num)

       r8vec_length = sqrt ( sum ( ( x(1:dim_num) )**2 ) )

       return
   end function r8vec_length

   subroutine rotation2(ib_position,ip_position,rot,rot_inverse)
       !-----------------------------------------------------------------------
       ! this subroutine compute a rotation matrix with eulerian angles between
       ! two catersian system 123 and 1'2'3' with angles phi,theta,psi.
       ! It uses x-convention and the right hand rules.
       !
       ! phi is the rotation around 3 that bring axis 1 to N the line of Nodes,
       ! where the line N is the intersection between plane 12 and 1'2'
       !
       ! theta is the rotation from 3 to 3' around N
       !
       ! psi is the rotation from N to 1' around 3
       !-----------------------------------------------------------------------
       implicit none
       !-----------------------------------------------------------------------
       ! declaration
       real :: p0(3),p1(3),p2(3),p3(3)
       real :: p0p(3),p1p(3),p3p(3)
       real :: p1_tan(3),p2_tan(3),p3_tan(3)
       real :: pN0(3),pNa(3)
       real :: p_base(3),F,G,H
       real :: normal_12(3),normal_1p2p(3),node_line(3)

       real :: check1,check2,check3
       real :: traslo_1,traslo_2,traslo_3
       real :: ip1,ip2,ip3
       real :: ip_position(3)
       real :: s1,s2,s3
       real :: ib_position(3)
       !real :: sx,sy,sz

       real :: rot(3,3),rot_inverse(3,3),det

       real :: phi,theta,psi,pi

       real :: c_phi,c_theta,c_psi
       real :: s_phi,s_theta,s_psi

       real :: p1_p1p,p2_p1p,p3_p1p
       real :: p1_p3p,p2_p3p,p3_p3p
       real :: pp1_p1p,pp2_p1p,pp3_p1p
       !-----------------------------------------------------------------------
       pi = acos(-1.)
       !-----------------------------------------------------------------------
       ! points on the general cartesian system 123

       p0(1)= 0.
       p0(2)= 0.
       p0(3)= 0.

       p1(1)= 1.
       p1(2)= 0.
       p1(3)= 0.

       p2(1)= 0.
       p2(2)= 1.
       p2(3)= 0.

       p3(1)= 0.
       p3(2)= 0.
       p3(3)= 1.

       !find the normal to plane 012
       call plane_exp_normal_3dbis( p0,p1,p2,normal_12)
       !-----------------------------------------------------------------------
       ! I need to work with right hand rule therefore I change y to z
       ip1 = ip_position(1)
       ip2 = ip_position(3) !ipz
       ip3 = ip_position(2) !ipy

       s1 = ib_position(1)
       s2 = ib_position(3) !sz
       s3 = ib_position(2) !sy
       !-----------------------------------------------------------------------
       ! points on the local cartesian axys 1'2'3' but written in the general
       ! coordinates system
       traslo_1 = ip1
       traslo_2 = ip2
       traslo_3 = ip3

       p0p(1)= ip1 - traslo_1
       p0p(2)= ip2 - traslo_2
       p0p(3)= ip3 - traslo_3

       p3p(1)= s1 - traslo_1
       p3p(2)= s2 - traslo_2
       p3p(3)= s3 - traslo_3

       call line_exp2par_3d (p0p(1),p0p(2),p0p(3),p3p(1),p3p(2),p3p(3), &
           F, G, H, p_base(1), p_base(2), p_base(3) )

       normal_1p2p(1)=F
       normal_1p2p(2)=G
       normal_1p2p(3)=H

       call plane_normal2exp_3d ( p0p, normal_1p2p,p1_tan,p2_tan,p3_tan)

       ! fix a point on axys 1'
       check1 = abs(p0p(1)-p1_tan(1)) &
           + abs(p0p(2)-p1_tan(2)) &
           + abs(p0p(3)-p1_tan(3))

       check2 = abs(p0p(1)-p2_tan(1)) &
           + abs(p0p(2)-p2_tan(2)) &
           + abs(p0p(3)-p2_tan(3))

       check3 = abs(p0p(1)-p3_tan(1)) &
           + abs(p0p(2)-p3_tan(2)) &
           + abs(p0p(3)-p3_tan(3))


       if(check1 /= 0)then
           p1p(1) = p1_tan(1)
           p1p(2) = p1_tan(2)
           p1p(3) = p1_tan(3)
       else
           if(check2 /= 0)then
               p1p(1) = p2_tan(1)
               p1p(2) = p2_tan(2)
               p1p(3) = p2_tan(3)
           else
               if(check2 /= 0)then
                   p1p(1) = p3_tan(1)
                   p1p(2) = p3_tan(2)
                   p1p(3) = p3_tan(3)
               end if
           end if
       end if

       !-----------------------------------------------------------------------
       ! construct N the line of node, intersection between 12 and 1'2'

       !vector product between normal to the plane
       call r8vec_cross_3d ( normal_1p2p, normal_12, node_line )

       pN0 = p0
       pNa(1) = node_line(1)*3.
       pNa(2) = node_line(2)*3.
       pNa(3) = node_line(3)*3.

       !-----------------------------------------------------------------------
       ! phi rotation from 1 to N
       ! phi0 = angle_rad_3d ( p1, p0, pNa )
       !-----------------------------------------------------------------------
       ! theta rotation from 3 to 3' [0,pi]
       ! theta0 = angle_rad_3d ( p3, p0, p3p )
       !-----------------------------------------------------------------------
       ! psi rotation from N to 1'
       ! psi0 = angle_rad_3d ( pNa, p0, p1p )
       !-----------------------------------------------------------------------

       ! phi rotate 1 to N around 3
       phi = atan2(pNa(2),pNa(1))
       c_phi = cos(phi)
       s_phi = sin(phi)

       !rotate points p1p and p3p using phi
       p1_p3p = c_phi*p3p(1) + s_phi*p3p(2)
       p2_p3p =-s_phi*p3p(1) + c_phi*p3p(2)
       p3_p3p = p3p(3)

       p1_p1p = c_phi*p1p(1) + s_phi*p1p(2)
       p2_p1p =-s_phi*p1p(1) + c_phi*p1p(2)
       p3_p1p = p1p(3)


       !------------------------------------------------------------------------
       ! theta rotate 3 to 3' around N [0,pi] or [-pi/2,pi/2]
       !the angle between 3 and 3', using atan2 I need to change the sign
       theta = -atan2(p2_p3p,p3_p3p)
       c_theta = cos(theta)
       s_theta = sin(theta)

       ! rotate p_p1p with theta to obtain pp_p1p
       pp1_p1p = p1_p1p
       pp2_p1p = c_theta*p2_p1p + s_theta*p3_p1p
       pp3_p1p =-s_theta*p2_p1p + c_theta*p3_p1p

       !psi rotate N to 1' around 3
       psi = atan2(pp2_p1p,pp1_p1p)

       c_psi = cos(psi)
       s_psi = sin(psi)

       ! rotation matrix construction
       rot(1,1) =  c_psi * c_phi  -  s_psi * c_theta * s_phi
       rot(1,2) =  c_psi * s_phi  +  s_psi * c_theta * c_phi
       rot(1,3) =  s_theta * s_psi

       rot(2,1) = -s_psi * c_phi  -  c_psi * c_theta * s_phi
       rot(2,2) = -s_psi * s_phi  +  c_psi * c_theta * c_phi
       rot(2,3) =  s_theta * c_psi

       rot(3,1) =  s_theta * s_phi
       rot(3,2) = -s_theta * c_phi
       rot(3,3) =  c_theta

       call r8mat_inverse_3d ( rot, rot_inverse, det )

       return
   end subroutine rotation2

   subroutine plane_exp_normal_3dbis(p1,p2,p3,normal)
       !! PLANE_EXP_NORMAL_3D finds the normal to an explicit plane in 3D.
       !
       !  Discussion:
       !    The explicit form of a plane in 3D is
       !      the plane through P1, P2 and P3.
       !  Modified:
       !    10 February 2005
       !  Author:
       !    John Burkardt
       !  Parameters:
       !
       !    Input, real :: P1(3), P2(3), P3(3), three points on the plane.
       !
       !    Output, real :: NORMAL(3), the coordinates of the unit normal
       !    vector to the plane containing the three points.
       !
       implicit none

       integer, parameter :: dim_num = 3

       real :: normal(dim_num)
       real :: normal_norm
       real :: p1(dim_num)
       real :: p2(dim_num)
       real :: p3(dim_num)

       !  do i=1,3
       !   p1(i)=vertice(i,1)
       !   p2(i)=vertice(i,2)
       !   p3(i)=vertice(i,3)
       !  end do
       !
       !  The cross product (P2-P1) x (P3-P1) is normal to (P2-P1) and (P3-P1).
       !
       normal(1) = ( p2(2) - p1(2) ) * ( p3(3) - p1(3) ) &
           - ( p2(3) - p1(3) ) * ( p3(2) - p1(2) )

       normal(2) = ( p2(3) - p1(3) ) * ( p3(1) - p1(1) ) &
           - ( p2(1) - p1(1) ) * ( p3(3) - p1(3) )

       normal(3) = ( p2(1) - p1(1) ) * ( p3(2) - p1(2) ) &
           - ( p2(2) - p1(2) ) * ( p3(1) - p1(1) )

       normal_norm = sqrt ( sum ( normal(1:dim_num)**2 ) )

       if ( normal_norm == 0.0D+00 ) then
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) 'PLANE_EXP_NORMAL_3D - Fatal error!'
           write ( *, '(a)' ) '  The plane is poorly defined.'
           stop
       end if

       normal(1:dim_num) = normal(1:dim_num) / normal_norm

       return
   end

   subroutine plane_normal2exp_3d(pp,normal,p1,p2,p3)

       !*****************************************************************************80
       !
       !! PLANE_NORMAL2EXP_3D converts a normal plane to explicit form in 3D.
       !
       !  Discussion:
       !
       !    The normal form of a plane in 3D is
       !
       !      PP, a point on the plane, and
       !      N, the unit normal to the plane.
       !
       !    The explicit form of a plane in 3D is
       !
       !      the plane through P1, P2 and P3.
       !
       !  Modified:
       !
       !    26 February 2005
       !
       !  Author:
       !
       !    John Burkardt
       !
       !  Parameters:
       !
       !    Input, real :: PP(3), a point on the plane.
       !
       !    Input, real :: NORMAL(3), a normal vector N to the plane.  The
       !    vector must not have zero length, but it is not necessary for N
       !    to have unit length.
       !
       !    Output, real :: P1(3), P2(3), P3(3), three points on the plane.
       !
       implicit none

       integer, parameter :: dim_num = 3

       real :: normal(dim_num)
       real :: p1(dim_num)
       real :: p2(dim_num)
       real :: p3(dim_num)
       real :: pp(dim_num)
       real :: pq(dim_num)
       real :: pr(dim_num)

       call plane_normal_basis_3d (normal,pq,pr)

       p1(1:dim_num) = pp(1:dim_num)
       p2(1:dim_num) = pp(1:dim_num) + pq(1:dim_num)
       p3(1:dim_num) = pp(1:dim_num) + pr(1:dim_num)

       return
   end

   subroutine plane_normal_basis_3d (normal,pq,pr)

       !*****************************************************************************80
       !
       !! PLANE_NORMAL_BASIS_3D finds two perpendicular vectors in a plane in 3D.
       !
       !  Discussion:
       !
       !    The normal form of a plane in 3D is:
       !
       !      PP is a point on the plane,
       !      N is a normal vector to the plane.
       !
       !    The two vectors to be computed, PQ and PR, can be regarded as
       !    the basis of a Cartesian coordinate system for points in the plane.
       !    Any point in the plane can be described in terms of the "origin"
       !    point PP plus a weighted sum of the two vectors PQ and PR:
       !
       !      P = PP + a * PQ + b * PR.
       !
       !    The vectors PQ and PR have unit length, and are perpendicular to N
       !    and to each other.
       !
       !  Modified:
       !
       !    27 August 2005
       !
       !  Author:
       !
       !    John Burkardt
       !
       !  Parameters:
       !
       !    Input, real :: PP(3), a point on the plane.  (Actually,
       !    we never need to know these values to do the calculation!)
       !
       !    Input, real :: NORMAL(3), a normal vector N to the plane.  The
       !    vector must not have zero length, but it is not necessary for N
       !    to have unit length.
       !
       !    Output, real :: PQ(3), a vector of unit length,
       !    perpendicular to the vector N and the vector PR.
       !
       !    Output, real :: PR(3), a vector of unit length,
       !    perpendicular to the vector N and the vector PQ.
       !
       implicit none

       integer, parameter :: dim_num = 3

       !real :: r8vec_length
       real :: normal(dim_num)
       real :: normal_norm
       real :: pq(dim_num)
       real :: pr(dim_num)
       real :: pr_norm
       !
       !  Compute the length of NORMAL.
       !
       normal_norm = r8vec_length ( dim_num, normal )

       if ( normal_norm == 0.0D+00 ) then
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) 'PLANE_NORMAL_BASIS_3D - Fatal error!'
           write ( *, '(a)' ) '  The normal vector is 0.'
           stop
       end if
       !
       !  Find a vector PQ that is normal to NORMAL and has unit length.
       !
       call r8vec_any_normal ( dim_num, normal, pq )
       !
       !  Now just take the cross product NORMAL x PQ to get the PR vector.
       !
       call r8vec_cross_3d ( normal, pq, pr )

       pr_norm = r8vec_length ( dim_num, pr )

       pr(1:dim_num) = pr(1:dim_num) / pr_norm

       return
   end

   subroutine r8vec_any_normal (dim_num,v1,v2)

       !*****************************************************************************
       !
       !! R8VEC_ANY_NORMAL returns some normal vector to V1.
       !
       !  Discussion:
       !
       !    If DIM_NUM < 2, then no normal vector can be returned.
       !
       !    If V1 is the zero vector, then any unit vector will do.
       !
       !    No doubt, there are better, more robust algorithms.  But I will take
       !    just about ANY reasonable unit vector that is normal to V1.
       !
       !  Modified:
       !
       !    23 August 2005
       !
       !  Author:
       !
       !    John Burkardt
       !
       !  Parameters:
       !
       !    Input, integer DIM_NUM, the spatial dimension.
       !
       !    Input, real :: V1(DIM_NUM), the vector.
       !
       !    Output, real :: V2(DIM_NUM), a vector that is
       !    normal to V2, and has unit Euclidean length.
       !
       implicit none

       integer dim_num

       !real :: r8vec_length
       integer i
       integer j
       integer k
       real :: v1(dim_num)
       real :: v2(dim_num)
       real :: vj
       real :: vk

       if ( dim_num < 2 ) then
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) 'R8VEC_ANY_NORMAL - Fatal error!'
           write ( *, '(a)' ) '  Called with DIM_NUM < 2.'
           stop
       end if

       if ( r8vec_length ( dim_num, v1 ) == 0.0D+00 ) then
           v2(1) = 1.0D+00
           v2(2:dim_num) = 0.0D+00
           return
       end if
       !
       !  Seek the largest entry in V1, VJ = V1(J), and the
       !  second largest, VK = V1(K).
       !
       !  Since V1 does not have zero norm, we are guaranteed that
       !  VJ, at least, is not zero.
       !
       j = -1
       vj = 0.0D+00

       k = -1
       vk = 0.0D+00

       do i = 1, dim_num

           if ( abs ( vk ) < abs ( v1(i) ) .or. k < 1 ) then

               if ( abs ( vj ) < abs ( v1(i) ) .or. j < 1 ) then
                   k = j
                   vk = vj
                   j = i
                   vj = v1(i)
               else
                   k = i
                   vk = v1(i)
               end if

           end if

       end do
       !
       !  Setting V2 to zero, except that V2(J) = -VK, and V2(K) = VJ,
       !  will just about do the trick.
       !
       v2(1:dim_num) = 0.0D+00

       v2(j) = -vk / sqrt ( vk * vk + vj * vj )
       v2(k) =  vj / sqrt ( vk * vk + vj * vj )

       return
   end

   subroutine r8mat_inverse_3d (a,b,det)

       !*****************************************************************************80
       !
       !! R8MAT_INVERSE_3D inverts a 3 by 3 real matrix using Cramer's rule.
       !
       !  Discussion:
       !
       !    If DET is zero, then A is singular, and does not have an
       !    inverse.  In that case, B is simply set to zero, and a
       !    message is printed.
       !
       !    If DET is nonzero, then its value is roughly an estimate
       !    of how nonsingular the matrix A is.
       !
       !  Modified:
       !
       !    16 April 1999
       !
       !  Author:
       !
       !    John Burkardt
       !
       !  Parameters:
       !
       !    Input, real :: A(3,3), the matrix to be inverted.
       !
       !    Output, real :: B(3,3), the inverse of the matrix A.
       !
       !    Output, real :: DET, the determinant of the matrix A.
       !
       implicit none

       real :: a(3,3)
       real :: b(3,3)
       real :: det
       !
       !  Compute the determinant of A
       !
       det =   a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
           + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
           + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )
       !
       !  If the determinant is zero, bail out.
       !
       if ( det == 0.0D+00 ) then

           b(1:3,1:3) = 0.0D+00

           return
       end if
       !
       !  Compute the entries of the inverse matrix using an explicit
       !  formula.
       !
       b(1,1) = + ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) / det
       b(1,2) = - ( a(1,2) * a(3,3) - a(1,3) * a(3,2) ) / det
       b(1,3) = + ( a(1,2) * a(2,3) - a(1,3) * a(2,2) ) / det

       b(2,1) = - ( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) / det
       b(2,2) = + ( a(1,1) * a(3,3) - a(1,3) * a(3,1) ) / det
       b(2,3) = - ( a(1,1) * a(2,3) - a(1,3) * a(2,1) ) / det

       b(3,1) = + ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) / det
       b(3,2) = - ( a(1,1) * a(3,2) - a(1,2) * a(3,1) ) / det
       b(3,3) = + ( a(1,1) * a(2,2) - a(1,2) * a(2,1) ) / det

       return
   end

   subroutine rotationAle(sx,sy,sz,ipx,ipy,ipz,rot,rot_inverse)
       !-----------------------------------------------------------------------
       ! This is specific for the IBM, the complexity of rotation2 is reduced

       implicit none

       real,intent(in) :: sx,sy,sz
       real,intent(in) :: ipx,ipy,ipz

       real,intent(out) :: rot(:,:),rot_inverse(:,:)

       ! declaration
       ! normal vector
       real :: nx,ny,nz
       ! tangent vector 1
       real :: tx,ty,tz
       ! tangent vector 2 (binormal)
       real :: bx,by,bz

       real :: norm_t,norm_n,norm_b

       real,parameter :: pi = acos(-1.)
       real,parameter :: toll=1.0e-5

       !-----------------------------------------------------------------------
       ! normal vector n
       nx=sx-ipx
       ny=sy-ipy
       nz=sz-ipz

       ! norm of n
       norm_n=sqrt(nx*nx+ny*ny+nz*nz)

       ! normalization of n
       nx=nx/norm_n
       ny=ny/norm_n
       nz=nz/norm_n

       ! recompute norm
       norm_n=sqrt(nx*nx+ny*ny+nz*nz)

       if (norm_n<1.0-toll .or. norm_n>1.0+toll) then
           write(*,*) 'Very bad error, norm fucked up'
       end if

       ! we try to align t1 or t2 to the x axes

       if (abs(nx)>toll) then
           tz=1/sqrt(1+(nz/nx)**2)
           tx=-1.0*tz*nz/nx
           ty=0.0
       else if (abs(nz)>toll) then
           ty=1/sqrt(1+(ny/nz)**2)
           tz=-1.0*ty*ny/nz
           tx=0.0
       else if (abs(ny)>toll) then
           tx=1/sqrt(1+(nx/ny)**2)
           ty=-1.0*tx*nx/ny
           tz=0.0
       else
           write(*,*) 'Problem in rotation routine, t computation. n_x=',nx,' ny=',ny,' nz=',nz
       end if

       ! norm of t
       norm_t=sqrt(tx*tx+ty*ty+tz*tz)

       ! compute secont tangent vector with cross product
       bx=ny*tz-nz*ty
       by=nz*tx-nx*tz
       bz=nx*ty-ny*tx

       ! norm of t
       norm_b=sqrt(bx*bx+by*by+bz*bz)

       if (norm_b<1.0-toll .or. norm_b>1.0+toll .or. &
           norm_t<1.0-toll .or. norm_t>1.0+toll) then
           write(*,*) 'Problem in rotation routine, norm_b=',norm_b,' norm_t=',norm_t
       end if

       ! rotation matrix construction
       rot(1,1)=tx
       rot(1,2)=ty
       rot(1,3)=tz

       rot(2,1)=bx
       rot(2,2)=by
       rot(2,3)=bz

       rot(3,1)=nx
       rot(3,2)=ny
       rot(3,3)=nz

       ! inverse is just the transpose
       rot_inverse(1,1)=rot(1,1)
       rot_inverse(1,2)=rot(2,1)
       rot_inverse(1,3)=rot(3,1)

       rot_inverse(2,1)=rot(1,2)
       rot_inverse(2,2)=rot(2,2)
       rot_inverse(2,3)=rot(3,2)

       rot_inverse(3,1)=rot(1,3)
       rot_inverse(3,2)=rot(2,3)
       rot_inverse(3,3)=rot(3,3)

       return

   end subroutine rotationAle

    subroutine spline_cubic_set(n,t,y,ibcbeg,ybcbeg,ibcend,ybcend,ypp)

        !c SPLINE_CUBIC_SET computes the second derivatives of a cubic spline.
        !
        !
        !  Discussion:
        !
        !    For data interpolation, the user must call SPLINE_CUBIC_SET to
        !    determine the second derivative data, passing in the data to be
        !    interpolated, and the desired boundary conditions.
        !
        !    The data to be interpolated, plus the SPLINE_CUBIC_SET output,
        !    defines the spline.  The user may then call SPLINE_CUBIC_VAL to
        !    evaluate the spline at any point.
        !
        !    The cubic spline is a piecewise cubic polynomial.  The intervals
        !    are determined by the "knots" or abscissas of the data to be
        !    interpolated.  The cubic spline has continous first and second
        !    derivatives over the entire interval of interpolation.
        !
        !    For any point T in the interval T(IVAL), T(IVAL+1), the form of
        !    the spline is
        !
        !      SPL(T) = A(IVAL)
        !             + B(IVAL) * ( T - T(IVAL) )
        !             + C(IVAL) * ( T - T(IVAL) )**2
        !             + D(IVAL) * ( T - T(IVAL) )**3
        !
        !    If we assume that we know the values Y(*) and YPP(*), which represent
        !    the values and second derivatives of the spline at each knot, then
        !    the coefficients can be computed as:
        !
        !      A(IVAL) = Y(IVAL)
        !      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
        !        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
        !      C(IVAL) = YPP(IVAL) / 2
        !      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
        !
        !    Since the first derivative of the spline is
        !
        !      SPL'(T) =     B(IVAL)
        !              + 2 * C(IVAL) * ( T - T(IVAL) )
        !              + 3 * D(IVAL) * ( T - T(IVAL) )**2,
        !
        !    the requirement that the first derivative be continuous at interior
        !    knot I results in a total of N-2 equations, of the form:
        !
        !      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1))
        !      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))**2 = B(IVAL)
        !
        !    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
        !
        !      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
        !      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
        !      + YPP(IVAL-1) * H(IVAL-1)
        !      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
        !      =
        !      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
        !      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
        !
        !    or
        !
        !      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
        !      + YPP(IVAL) * H(IVAL)
        !      =
        !      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
        !      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
        !
        !    Boundary conditions must be applied at the first and last knots.
        !    The resulting tridiagonal system can be solved for the YPP values.
        !
        !  Modified:
        !
        !    20 November 2000
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, integer N, the number of data points; N must be at least 2.
        !
        !    Input, real T(N), the points where data is specified.
        !    The values should be distinct, and increasing.
        !
        !    Input, real Y(N), the data values to be interpolated.
        !
        !    Input, integer IBCBEG, the left boundary condition flag:
        !
        !      0: the spline should be a quadratic over the first interval;
        !      1: the first derivative at the left endpoint should be YBCBEG;
        !      2: the second derivative at the left endpoint should be YBCBEG.
        !
        !    Input, real YBCBEG, the left boundary value, if needed.
        !
        !    Input, integer IBCEND, the right boundary condition flag:
        !
        !      0: the spline should be a quadratic over the last interval;
        !      1: the first derivative at the right endpoint should be YBCEND;
        !      2: the second derivative at the right endpoint should be YBCEND.
        !
        !    Input, real YBCEND, the right boundary value, if needed.
        !
        !    Output, real YPP(N), the second derivatives of the cubic spline.
        !
        implicit none
        !
        integer n
        !
        real diag(n)
        integer i
        integer ibcbeg
        integer ibcend
        real sub(2:n)
        real sup(1:n-1)
        real t(n)
        real y(n)
        real ybcbeg
        real ybcend
        real ypp(n)
        !
        !  Check.
        !
        if ( n <= 1 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal errorc'
            write ( *, '(a)' ) '  The number of knots must be at least 2.'
            write ( *, '(a,i6)' ) '  The input value of N = ', n
            stop
        end if

        do i = 1, n-1
            if ( t(i) >= t(i+1) ) then
                write(*,'(a)')' '
                write(*,'(a)')'SPLINE_CUBIC_SET - Fatal errorc'
                write(*,'(a)')'  The knots must be strictly increasing, but'
                write(*,'(a,i6,a,g14.6)' ) '  T(',  i,') = ', t(i)
                write(*,'(a,i6,a,g14.6)' ) '  T(',i+1,') = ', t(i+1)
                stop
            end if
        end do
        !
        !  Set the first equation.
        !
        if ( ibcbeg == 0 ) then
            ypp(1) = 0.0E+00
            diag(1) = 1.0E+00
            sup(1) = -1.0E+00
        else if ( ibcbeg == 1 ) then
            ypp(1) = ( y(2) - y(1) ) / ( t(2) - t(1) ) - ybcbeg
            diag(1) = ( t(2) - t(1) ) / 3.0E+00
            sup(1) = ( t(2) - t(1) ) / 6.0E+00
        else if ( ibcbeg == 2 ) then
            ypp(1) = ybcbeg
            diag(1) = 1.0E+00
            sup(1) = 0.0E+00
        else
            write(*,'(a)') ' '
            write(*,'(a)') 'SPLINE_CUBIC_SET - Fatal errorc'
            write(*,'(a)') '  The boundary flag IBCBEG must be 0, 1 or 2.'
            write(*,'(a,i6)' ) '  The input value is IBCBEG = ', ibcbeg
            stop
        end if
        !
        !  Set the intermediate equations.
        !
        do i = 2, n-1
            ypp(i) = ( y(i+1) - y(i) ) / ( t(i+1) - t(i) )  &
                - ( y(i) - y(i-1) ) / ( t(i) - t(i-1) )
            sub(i) = ( t(i) - t(i-1) ) / 6.0E+00
            diag(i) = ( t(i+1) - t(i-1) ) / 3.0E+00
            sup(i) = ( t(i+1) - t(i) ) / 6.0E+00
        end do
        !
        !  Set the last equation.
        !
        if ( ibcend == 0 ) then
            ypp(n) = 0.0E+00
            sub(n) = -1.0E+00
            diag(n) = 1.0E+00
        else if ( ibcend == 1 ) then
            ypp(n) = ybcend - ( y(n) - y(n-1) ) / ( t(n) - t(n-1) )
            sub(n) = ( t(n) - t(n-1) ) / 6.0E+00
            diag(n) = ( t(n) - t(n-1) ) / 3.0E+00
        else if ( ibcend == 2 ) then
            ypp(n) = ybcend
            sub(n) = 0.0E+00
            diag(n) = 1.0E+00
        else
            write (*,'(a)') ' '
            write (*,'(a)') 'SPLINE_CUBIC_SET - Fatal errorc'
            write (*,'(a)') '  The boundary flag IBCEND must be 0, 1 or 2.'
            write (*,'(a,i6)') '  The input value is IBCEND = ', ibcend
            stop
        end if
        !
        !  Special case:
        !    N = 2, IBCBEG = IBCEND = 0.
        !
        if ( n == 2 .and. ibcbeg == 0 .and. ibcend == 0 ) then

            ypp(1) = 0.0E+00
            ypp(2) = 0.0E+00
        !
        !  Solve the linear system.
        !
        else

            call s3_fs ( sub, diag, sup, n, ypp, ypp )

        end if

        return
    end

    subroutine spline_cubic_val(n,t,y,ypp,tval,yval,ypval,yppval)

        !c SPLINE_CUBIC_VAL evaluates a cubic spline at a specific point.
        !
        !
        !  Discussion:
        !
        !    SPLINE_CUBIC_SET must have already been called to define the
        !    values of YPP.
        !
        !    For any point T in the interval T(IVAL), T(IVAL+1), the form of
        !    the spline is
        !
        !      SPL(T) = A
        !             + B * ( T - T(IVAL) )
        !             + C * ( T - T(IVAL) )**2
        !             + D * ( T - T(IVAL) )**3
        !
        !    Here:
        !      A = Y(IVAL)
        !      B = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
        !        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
        !      C = YPP(IVAL) / 2
        !      D = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
        !
        !  Modified:
        !
        !    20 November 2000
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, integer N, the number of data values.
        !
        !    Input, real T(N), the knot values.
        !
        !    Input, real Y(N), the data values at the knots.
        !
        !    Input, real YPP(N), the second derivatives of the spline at the knots.
        !
        !    Input, real TVAL, a point, typically between T(1) and T(N), at
        !    which the spline is to be evalulated.  If TVAL lies outside
        !    this range, extrapolation is used.
        !
        !    Output, real YVAL, YPVAL, YPPVAL, the value of the spline, and
        !    its first two derivatives at TVAL.
        !
        implicit none
        !
        integer n
        !
        real dt
        real h
        integer left
        integer right
        real t(n)
        real tval
        real y(n)
        real ypp(n)
        real yppval
        real ypval
        real yval
        !
        !  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
        !  Values below T(1) or above T(N) use extrapolation.
        !
        call rvec_bracket ( n, t, tval, left, right )
        !
        !  Evaluate the polynomial.
        !
        dt = tval - t(left)
        h = t(right) - t(left)

        yval = y(left)   &
            + dt * ( ( y(right) - y(left) ) / h  &
            - ( ypp(right) / 6.0E+00 + ypp(left) / 3.0E+00 ) * h  &
            + dt * ( 0.5E+00 * ypp(left)  &
            + dt * ( ( ypp(right) - ypp(left) ) / ( 6.0E+00 * h ) ) ) )

        ypval = ( y(right) - y(left) ) / h  &
            - ( ypp(right) / 6.0E+00 + ypp(left) / 3.0E+00 ) * h  &
            + dt * ( ypp(left)  &
            + dt * ( 0.5E+00 * ( ypp(right) - ypp(left) ) / h ) )

        yppval = ypp(left) + dt * ( ypp(right) - ypp(left) ) / h

        return
    end

    subroutine s3_fs(a1,a2,a3,n,b,x)

        !c S3_FS factors and solves a tridiagonal linear system.
        !
        !
        !  Note:
        !
        !    This algorithm requires that each diagonal entry be nonzero.
        !
        !  Modified:
        !
        !    05 December 1998
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input/output, real A1(2:N), A2(1:N), A3(1:N-1).
        !    On input, the nonzero diagonals of the linear system.
        !    On output, the data in these vectors has been overwritten
        !    by factorization information.
        !
        !    Input, integer N, the order of the linear system.
        !
        !    Input/output, real B(N).
        !    On input, B contains the right hand side of the linear system.
        !    On output, B has been overwritten by factorization information.
        !
        !    Output, real X(N), the solution of the linear system.
        !
        implicit none
        !
        integer n
        !
        real a1(2:n)
        real a2(1:n)
        real a3(1:n-1)
        real b(n)
        integer i
        real x(n)
        real xmult
        !
        !  The diagonal entries can't be zero.
        !
        do i = 1, n
            if ( a2(i) == 0.0E+00 ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'S3_FS - Fatal errorc'
                write ( *, '(a,i6,a)' ) '  A2(', i, ') = 0.'
                return
            end if
        end do

        do i = 2, n-1

            xmult = a1(i) / a2(i-1)
            a2(i) = a2(i) - xmult * a3(i-1)

            b(i) = b(i) - xmult * b(i-1)

        end do

        xmult = a1(n) / a2(n-1)
        a2(n) = a2(n) - xmult * a3(n-1)

        x(n) = ( b(n) - xmult * b(n-1) ) / a2(n)
        do i = n-1, 1, -1
            x(i) = ( b(i) - a3(i) * x(i+1) ) / a2(i)
        end do

        return
    end

    subroutine rvec_bracket(n,x,xval,left,right)

        !c RVEC_BRACKET searches a sorted array for successive brackets of a value.
        !
        !
        !  Discussion:
        !
        !    If the values in the vector are thought of as defining intervals
        !    on the real line, then this routine searches for the interval
        !    nearest to or containing the given value.
        !
        !  Modified:
        !
        !    06 April 1999
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, integer N, length of input array.
        !
        !    Input, real X(N), an array sorted into ascending order.
        !
        !    Input, real XVAL, a value to be bracketed.
        !
        !    Output, integer LEFT, RIGHT, the results of the search.
        !    Either:
        !      XVAL < X(1), when LEFT = 1, RIGHT = 2;
        !      XVAL > X(N), when LEFT = N-1, RIGHT = N;
        !    or
        !      X(LEFT) <= XVAL <= X(RIGHT).
        !
        implicit none
        !
        integer n
        !
        integer i
        integer left
        integer right
        real x(n)
        real xval
        !
        do i = 2, n - 1

            if ( xval < x(i) ) then
                left = i - 1
                right = i
                return
            end if

        end do

        left = n - 1
        right = n

        return
    end

    subroutine line_angle(p1,p2,q1,q2,angle_deg)

        ! THIS SUBROUTINE IS TAKEN FROM LIBRARY GEOMETRY PACK
        !
        !c LINES_EXP_ANGLE_3D finds the angle between two explicit lines in 3D.
        !
        !  Formula:
        !
        !    The explicit form of a line in 3D is:
        !
        !      P1(X1,Y1,Z1), P2(X2,Y2,Z2).
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
        !    Input, real P1(X1,Y1,Z1) P2(X2,Y2,Z2), two distince points on the first line.
        !
        !    Input, real Q1(X3,Y3,Z3),Q2(X4,Y4,Z4), two distinct points on the second line.
        !
        !    Output, real ANGLE, the angle in radians between the two lines.
        !    ANGLE is computed using the ACOS function, and so lies between 0 and PI.
        !    But if one of the lines is degenerate, ANGLE is returned as -1.0.
        !
        implicit none
        !
        !real,intent(in) :: x1,x2,x3,x4
        !real,intent(in) :: y1,y2,y3,y4
        !real,intent(in) :: z1,z2,z3,z4
        !real,intent(in) :: z1,z2,z3,z4
        real,dimension(3),intent(in) :: p1,p2,q1,q2
        real,intent(out) :: angle_deg

        ! line direction
        real,dimension(3) :: p,q
        real :: angle,arc_cosine
        real :: ctheta,enorm0_3d,pdotq,pnorm,qnorm
        real,parameter :: tol=0.000001
        real,parameter :: pi=acos(-1.0)

        ! directions of two lines
        p(:)=p2(:)-p1(:)
        q(:)=q2(:)-q1(:)

        ! distance between points in the lines
        pnorm=norm2(p(:))
        qnorm=norm2(q(:))

        ! dot product of two directions
        pdotq=dot_product(p(:),q(:))


        if (pnorm==0.0 .or. qnorm==0.0) then
            write (*,*) ' '
            write (*,*) 'LINES_EXP_ANGLE_3D - Fatal error!'
            write (*,*) '  One of the lines is degenerate!'
            angle = - 1.0
        else

            ctheta = pdotq / ( pnorm * qnorm )

            if (abs(ctheta)>=1.0) then
                write(1000,*)'TROVATO ANGOLO DEGENERE',ctheta
                angle_deg = 0.0
            else
                ! arc_cosine ( ctheta )
                angle =abs(acos(ctheta))
                angle_deg=angle*180./pi
                angle_deg=angle
            end if
        end if


        return
    end

    subroutine plane_3points(p1,p2,p3,a,b,c,d)
        !-------------------------------------------------------------------------------
        ! dati tre punti di coordinate p1(x1,y1,z1),p2(x2,y2,z2),p3(x3,y3,z3) trovo il piano
        ! che vi passa di eq: ax+by+cz+d=0
        !-------------------------------------------------------------------------------
        implicit none
        !
        real,intent(out) :: a,b,c,d
        real,dimension(3),intent(in) :: p1,p2,p3

        a=(p2(2)-p1(2))*(p3(3)-p1(3))-(p2(3)-p1(3))*(p3(2)-p1(2))
        b=(p2(3)-p1(3))*(p3(1)-p1(1))-(p2(1)-p1(1))*(p3(3)-p1(3))
        c=(p2(1)-p1(1))*(p3(2)-p1(2))-(p2(2)-p1(2))*(p3(1)-p1(1))
        d=-p2(1)*a-p2(2)*b-p2(3)*c

        return
    end subroutine plane_3points

    subroutine find_plane_projection(a,b,c,d,p,pn)
        !------------------------------------------------------------
        ! dato un piano di equazione ax+by+cz+d=0 ed un punto p(x,y,z)
        ! trovo il punto appartenente al piano pn(xn,yn,zn) piu' vicino
        ! al punto sopra
        !------------------------------------------------------------
        ! COME PROCEDO:
        ! normale N al piano (a,b,c)
        ! la linea definita da (xn-x)/a = (yn-y)/b = (zn-z)/c = t
        ! passa per il punto (x,y,z) ed e' parallela a N
        !
        ! risolvendo per il punto (xn,yn,zn) ottengo:
        !
        !  xn = a * t + x
        !  yn = b * t + y
        !  zn = c * t + z
        !
        ! pongo questi valori nell'equazione del piano e risolvo per t
        !
        ! a*(a * t + x) + b*(b * t + y) + c*(c * t + z) + d = 0
        !
        ! t=(-a*x -b*y -c*z -d) / (a*a + b*b + c*c)
        !------------------------------------------------------------

        implicit none

        real,intent(in):: a,b,c,d
        real,dimension(3),intent(in) :: p
        real,dimension(3),intent(out) :: pn

        real:: t

        if ( a == 0.0D+00 .and. b == 0.0D+00 .and. c == 0.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PLANE_IMP_POINT_NEAR_3D - Fatal error!'
            write ( *, '(a)' ) '  A = B = C = 0.'
            stop
        else

            t=-(a*p(1)+b*p(2)+c*p(3)+d)/(a*a+b*b+c*c)

            pn(1)=p(1)+a*t
            pn(2)=p(2)+b*t
            pn(3)=p(3)+c*t

        end if

        return

    end subroutine find_plane_projection

   end module geometricRoutines
