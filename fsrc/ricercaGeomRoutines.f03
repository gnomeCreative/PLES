!---------------------------------------------------------------------------------------------------------------------------------------------
!*********************************************************************************************************************************************
!---------------------------------------------------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
!********************************************************************************
!--------------------------------------------------------------------------------
subroutine piano_per3punti ( x1, y1, z1, x2, y2, z2, x3, y3, z3, a, b, c, d )
    !-------------------------------------------------------------------------------
    ! dati tre punti di coordinate (x1,y1,z1),(x2,y2,z2),(x3,y3,z3) trovo il piano
    ! che vi passa di eq: ax+by+cz+d=0
    !-------------------------------------------------------------------------------
    implicit none
    !
    real(8):: a,b,c,d
    real(8):: x1,y1,z1
    real(8):: x2,y2,z2
    real(8):: x3,y3,z3

    a = ( y2 - y1 ) * ( z3 - z1 ) - ( z2 - z1 ) * ( y3 - y1 )
    b = ( z2 - z1 ) * ( x3 - x1 ) - ( x2 - x1 ) * ( z3 - z1 )
    c = ( x2 - x1 ) * ( y3 - y1 ) - ( y2 - y1 ) * ( x3 - x1 )
    d = - x2 * a - y2 * b - z2 * c

    return
end

!--------------------------------------------------------------------------------
!********************************************************************************
!--------------------------------------------------------------------------------
subroutine punto_proiezione_piano(a,b,c,d,x,y,z,xn,yn,zn)
    !------------------------------------------------------------
    ! dato un piano di equazione ax+by+cz+d=0 ed un punto (x,y,z)
    ! trovo il punto appartenente al piano (xn,yn,zn) piu' vicino
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

    real(8):: a,b,c,d
    real(8):: x,y,z
    real(8):: xn,yn,zn
    real(8):: t

    if ( a == 0.0D+00 .and. b == 0.0D+00 .and. c == 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PLANE_IMP_POINT_NEAR_3D - Fatal error!'
        write ( *, '(a)' ) '  A = B = C = 0.'
        stop
    else

        t = - ( a * x + b * y + c * z + d ) / ( a * a + b * b + c * c )

        xn = x + a * t
        yn = y + b * t
        zn = z + c * t

    end if

    return
end

!---------------------------------------------------------------------
!********************************************************************************
!--------------------------------------------------------------------------------
subroutine plane_imp_line_par_int_3d (a,b,c,d,x0,y0,z0,f,g,h,intersect,x,y,z,coincide,intersezione)
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
    real(8) :: a
    real(8) :: b
    real(8) :: c
    real(8) :: d
    real(8) :: denom
    real(8) :: f
    real(8) :: g
    real(8) :: h
    logical :: intersect
    real(8) :: norm1
    real(8) :: norm2
    real(8) :: t
    real(8), parameter :: tol = 0.00001E+00
    real(8) :: x
    real(8) :: x0
    real(8) :: y
    real(8) :: y0
    real(8) :: z
    real(8) :: z0
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
end

!--------------------------------------------------------------------------------
!********************************************************************************
!--------------------------------------------------------------------------------
subroutine line_exp2par_3d ( x1, y1, z1, x2, y2, z2, f, g, h, x0, y0, z0 )
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
    real(8) :: f
    real(8) :: g
    real(8) :: h
    real(8) :: norm
    real(8) :: x0
    real(8) :: x1
    real(8) :: x2
    real(8) :: y0
    real(8) :: y1
    real(8) :: y2
    real(8) :: z0
    real(8) :: z1
    real(8) :: z2
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
end

!--------------------------------------------------------------------------------
!********************************************************************************
!--------------------------------------------------------------------------------
subroutine ptpolg ( dim, vcl, qx, qy, qz, nrml, dtol, inout )

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

10 continue

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
       go to 10
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
       go to 10
   end if

   if ( dim == 2 ) then
       dir(1) = 0.5D+00 * ( vcl(1,la) + vcl(1,lb) ) - pt(1)
       dir(2) = 0.5D+00 * ( vcl(2,la) + vcl(2,lb) ) - pt(2)
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
       dir(1) = 0.5D+00 * ( vcl(1,la) + vcl(1,lb) ) - pt(1)
       dir(2) = 0.5D+00 * ( vcl(2,la) + vcl(2,lb) ) - pt(2)
       dir(3) = 0.5D+00 * ( vcl(3,la) + vcl(3,lb) ) - pt(3)
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
   end

   !--------------------------------------------------------------------------------
   !********************************************************************************
   !--------------------------------------------------------------------------------
   subroutine plane_exp_normal_3d ( vertice, normal )
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
   end


   !--------------------------------------------------------------------------------
   !********************************************************************************
   !--------------------------------------------------------------------------------
   subroutine line_exp_point_near_3d ( p1, p2, p, pn, dist, t )

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
       !    Input, real ( kind = 8 ) P1(3), P2(3), two points on the line.
       !
       !    Input, real ( kind = 8 ) P(3), the point whose nearest neighbor on
       !    the line is to be determined.
       !
       !    Output, real ( kind = 8 ) PN(3), the point which is the nearest
       !    point on the line to P.
       !
       !    Output, real ( kind = 8 ) DIST, the distance from the point to the
       !    nearest point on the line.
       !
       !    Output, real ( kind = 8 ) T, the relative position of the point
       !    PN to P1 and P2.
       !
       implicit none

       integer, parameter :: dim_num = 3

       real ( kind = 8 ) bot
       real ( kind = 8 ) dist
       logical line_exp_is_degenerate_nd
       real ( kind = 8 ) p(dim_num)
       real ( kind = 8 ) p1(dim_num)
       real ( kind = 8 ) p2(dim_num)
       real ( kind = 8 ) pn(dim_num)
       real ( kind = 8 ) t

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

   !--------------------------------------------------------------------------------
   !********************************************************************************
   !--------------------------------------------------------------------------------
   function line_exp_is_degenerate_nd ( dim_num, p1, p2 )

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
       !    Input, real ( kind = 8 ) P1(DIM_NUM), P2(DIM_NUM), two points on the line.
       !
       !    Output, logical LINE_EXP_IS_DEGENERATE_ND, is TRUE if the line
       !    is degenerate.
       !
       implicit none

       integer dim_num

       logical line_exp_is_degenerate_nd
       real ( kind = 8 ) p1(dim_num)
       real ( kind = 8 ) p2(dim_num)

       line_exp_is_degenerate_nd = ( all ( p1(1:dim_num) == p2(1:dim_num) ) )

       return
   end


   !--------------------------------------------------------------------------------
   !********************************************************************************
   !--------------------------------------------------------------------------------
   function parallelepiped_contains_point_3d ( p1, p2, p3, p4, p )

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
       !    Input, real ( kind = 8 ) P1(3), P2(3), P3(3), P4(3), four corners
       !    of the parallelepiped.  It is assumed that P2, P3 and P4 are
       !    immediate neighbors of P1.
       !
       !    Input, real ( kind = 8 ) P(3), the point to be checked.
       !
       !    Output, logical PARALLELEPIPED_CONTAINS_POINT_3D, is true if P
       !    is inside the parallelepiped, or on its boundary.
       !
       implicit none

       integer, parameter :: dim_num = 3

       real ( kind = 8 ) dot
       logical parallelepiped_contains_point_3d
       real ( kind = 8 ) p(dim_num)
       real ( kind = 8 ) p1(dim_num)
       real ( kind = 8 ) p2(dim_num)
       real ( kind = 8 ) p3(dim_num)
       real ( kind = 8 ) p4(dim_num)

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
   end

   !--------------------------------------------------------------------------------
   !********************************************************************************
   !--------------------------------------------------------------------------------
   subroutine poliedro(ii,jj,kk,jx,jy,jz,xc,yc,zc,dentro,tx,ty,tz)
       implicit none
       integer :: ii,jj,kk
       integer :: jx,jy,jz
       real(8) :: xc(0:jx+1,0:jy+1,0:jz+1),yc(0:jx+1,0:jy+1,0:jz+1),zc(0:jx+1,0:jy+1,0:jz+1)
       real(8) :: v(3,8),p(3),p1(3),p2(3)
       integer :: face_order(6),i,j
       integer :: face_point(4,6)
       logical :: inside,inside_line,dentro
       real(8) :: dist,dist1,dist2,differenza,tx,ty,tz

       ! poliedro vertices
       v(1,1)=xc(ii-1,jj-1,kk-1)
       v(2,1)=yc(ii-1,jj-1,kk-1)
       v(3,1)=zc(ii-1,jj-1,kk-1)

       v(1,2)=xc(ii+1,jj-1,kk-1)
       v(2,2)=yc(ii+1,jj-1,kk-1)
       v(3,2)=zc(ii+1,jj-1,kk-1)

       v(1,3)=xc(ii+1,jj+1,kk-1)
       v(2,3)=yc(ii+1,jj+1,kk-1)
       v(3,3)=zc(ii+1,jj+1,kk-1)

       v(1,4)=xc(ii-1,jj+1,kk-1)
       v(2,4)=yc(ii-1,jj+1,kk-1)
       v(3,4)=zc(ii-1,jj+1,kk-1)
       !---------------
       v(1,5)=xc(ii-1,jj-1,kk+1)
       v(2,5)=yc(ii-1,jj-1,kk+1)
       v(3,5)=zc(ii-1,jj-1,kk+1)

       v(1,6)=xc(ii+1,jj-1,kk+1)
       v(2,6)=yc(ii+1,jj-1,kk+1)
       v(3,6)=zc(ii+1,jj-1,kk+1)

       v(1,7)=xc(ii+1,jj+1,kk+1)
       v(2,7)=yc(ii+1,jj+1,kk+1)
       v(3,7)=zc(ii+1,jj+1,kk+1)

       v(1,8)=xc(ii-1,jj+1,kk+1)
       v(2,8)=yc(ii-1,jj+1,kk+1)
       v(3,8)=zc(ii-1,jj+1,kk+1)

       do i=1,6
           face_order(i)=4
       end do

       ! surface defined by points
       face_point(1,1)=1
       face_point(2,1)=2
       face_point(3,1)=3
       face_point(4,1)=4


       face_point(1,2)=2
       face_point(2,2)=6
       face_point(3,2)=7
       face_point(4,2)=3


       face_point(1,3)=6
       face_point(2,3)=5
       face_point(3,3)=8
       face_point(4,3)=7


       face_point(1,4)=5
       face_point(2,4)=1
       face_point(3,4)=4
       face_point(4,4)=8


       face_point(1,5)=2
       face_point(2,5)=1
       face_point(3,5)=5
       face_point(4,5)=6


       face_point(1,6)=4
       face_point(2,6)=3
       face_point(3,6)=7
       face_point(4,6)=8


       p(1)=tx
       p(2)=ty
       p(3)=tz


       inside = .false.

       ! check if the point stay on the surface
       inside_line = .false.
       do i=1,8
           p1(1)=v(1,i)
           p1(2)=v(2,i)
           p1(3)=v(3,i)

           do j=1,8
               if(j==i)cycle
               p2(1)=v(1,j)
               p2(2)=v(2,j)
               p2(3)=v(3,j)

               dist1=sqrt( sum((p(1:3)-p1(1:3))**2.))
               dist2=sqrt( sum((p(1:3)-p2(1:3))**2.))

               dist=sqrt( sum((p1(1:3)-p2(1:3))**2.))

               !call line_exp_point_dist_3d ( p1, p2, p, dist )

               differenza = abs(dist - (dist1 + dist2))

               if(differenza.lt.1.0D-20)then
                   inside_line=.true.
                   goto 10
               end if
           end do
       end do
10 continue

   call polyhedron_contains_point_3d( 8, 6, 4, v, face_order, face_point, p, inside )

   dentro=.false.

   if(inside)then
       dentro=.true.
   else
       if(inside_line)then
           dentro=.true.
       else
           dentro=.false.
       end if
   end if

   return
   end

   !--------------------------------------------------------------------------------
   !********************************************************************************
   !--------------------------------------------------------------------------------
   subroutine polyhedron_contains_point_3d  ( node_num, face_num, &
       face_order_max, v, face_order, face_point, p, inside )

       !*********************************************************************
       !
       !! POLYHEDRON_CONTAINS_POINT_3D determines if a point is inside a polyhedron.
       !
       !  Discussion:
       !
       !    The reference states that the polyhedron should be simple (that
       !    is, the faces should form a single connected surface), and that
       !    the individual faces should be consistently oriented.
       !
       !    However, the polyhedron does not, apparently, need to be convex.
       !
       !  Modified:
       !
       !    30 August 2005
       !
       !  Author:
       !
       !    John Burkardt
       !
       !  Reference:
       !
       !    Paulo Cezar Pinto Carvalho, Paulo Roma Cavalcanti,
       !    Point in Polyhedron Testing Using Spherical Polygons,
       !    in Graphics Gems V,
       !    edited by Alan Paeth,
       !    Academic Press, 1995, T385.G6975.
       !
       !  Parameters:
       !
       !    Input, integer NODE_NUM, the number of vertices.
       !
       !    Input, integer FACE_NUM, the number of faces.
       !
       !    Input, integer FACE_ORDER_MAX, the maximum order of any face.
       !
       !    Input, real ( kind = 8 ) V(3,NODE_NUM), the coordinates of the vertices.
       !
       !    Input, integer FACE_ORDER(FACE_NUM), the order of each face.
       !
       !    Input, integer FACE_POINT(FACE_ORDER_MAX,FACE_NUM), the indices of the
       !    nodes that make up each face.
       !
       !    Input, real ( kind = 8 ) P(3), the point to be tested.
       !
       !    Output, logical INSIDE, is true if the point
       !    is inside the polyhedron.
       !
       implicit none

       integer, parameter :: dim_num = 3
       integer face_num
       integer face_order_max
       integer node_num

       real ( kind = 8 ) area
       integer face
       integer face_order(face_num)
       integer face_point(face_order_max,face_num)
       logical inside
       integer k
       integer node
       integer node_num_face
       real ( kind = 8 ) p(dim_num)
       real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
       real ( kind = 8 ) solid_angle
       real ( kind = 8 ) v(dim_num,node_num)
       real ( kind = 8 ) v_face(dim_num,face_order_max)
       !-------------------------------------------------------------





       !--------------------------------------------------------------
       area = 0.0D+00

       do face = 1, face_num

           node_num_face = face_order(face)

           do k = 1, node_num_face

               node = face_point(k,face)

               v_face(1:dim_num,k) = v(1:dim_num,node)

           end do

           call polygon_solid_angle_3d ( node_num_face, v_face, p, solid_angle )

           area = area + solid_angle

       end do
       !
       !  AREA should be -4*PI, 0, or 4*PI.
       !  So this test should be quite safe!
       !
       if ( area < -2.0D+00 * pi .or. 2.0D+00 * pi < area ) then
           inside = .true.
       else
           inside = .false.
       end if

       return
   end

   !--------------------------------------------------------------------------------
   !********************************************************************************
   !--------------------------------------------------------------------------------
   subroutine polygon_solid_angle_3d ( n, v, p, solid_angle )

       !*********************************************************************
       !
       !! POLYGON_SOLID_ANGLE_3D: projected solid angle of a 3D plane polygon.
       !
       !  Discussion:
       !
       !    A point P is at the center of the unit sphere.  A planar polygon
       !    is to be projected onto the surface of the unit sphere, by drawing
       !    the ray from P to each polygonal vertex, and noting where this ray
       !    intersects the unit sphere.  The area of the projected polygon is
       !    equal to the solid angle, since we are considering the unit sphere.
       !
       !  Modified:
       !
       !    30 August 2005
       !
       !  Author:
       !
       !    John Burkardt
       !
       !  Reference:
       !
       !    Paulo Cezar Pinto Carvalho, Paulo Roma Cavalcanti,
       !    Point in Polyhedron Testing Using Spherical Polygons,
       !    in Graphics Gems V,
       !    edited by Alan Paeth,
       !    Academic Press, 1995, T385.G6975.
       !
       !  Parameters:
       !
       !    Input, integer N, the number of vertices.
       !
       !    Input, real ( kind = 8 ) V(3,N), the coordinates of the vertices.
       !
       !    Input, real ( kind = 8 ) P(3), the point at the center of the unit sphere.
       !
       !    Output, double SOLID_ANGLE, the solid angle subtended
       !    by the polygon, as projected onto the unit sphere around the point P.
       !
       implicit none

       integer, parameter :: dim_num = 3
       integer n

       real ( kind = 8 ) a(dim_num)
       real ( kind = 8 ) angle
       real ( kind = 8 ) arc_cosine
       real ( kind = 8 ) area
       real ( kind = 8 ) b(dim_num)
       real ( kind = 8 ) r8vec_length
       real ( kind = 8 ) r8vec_triple_product
       integer i4_wrap
       integer j
       integer jp1
       real ( kind = 8 ) normal1(dim_num)
       real ( kind = 8 ) normal1_norm
       real ( kind = 8 ) normal2(dim_num)
       real ( kind = 8 ) normal2_norm
       real ( kind = 8 ) p(dim_num)
       real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
       real ( kind = 8 ) plane(dim_num)
       real ( kind = 8 ) r1(dim_num)
       real ( kind = 8 ) s
       real ( kind = 8 ) solid_angle
       real ( kind = 8 ) v(dim_num,n)

       if ( n < 3 ) then
           solid_angle = 0.0D+00
           return
       end if

       call polygon_normal_3d ( n, v, plane )

       a(1:dim_num) = v(1:dim_num,n) - v(1:dim_num,1)

       area = 0.0D+00

       do j = 1, n

           r1(1:dim_num) = v(1:dim_num,j) - p(1:dim_num)

           jp1 = i4_wrap ( j + 1, 1, n )

           b(1:dim_num) = v(1:dim_num,jp1) - v(1:dim_num,j)

           call r8vec_cross_3d ( a, r1, normal1 )

           normal1_norm = r8vec_length ( dim_num, normal1 )

           call r8vec_cross_3d ( r1, b, normal2 )

           normal2_norm = r8vec_length ( dim_num, normal2 )

           s = dot_product ( normal1(1:dim_num), normal2(1:dim_num) ) &
               / ( normal1_norm * normal2_norm )

           angle = arc_cosine ( s )

           s = r8vec_triple_product ( b, a, plane )

           if ( 0.0D+00 < s ) then
               area = area + pi - angle
           else
               area = area + pi + angle
           end if

           a(1:dim_num) = -b(1:dim_num)

       end do

       area = area - pi * real ( n - 2, kind = 8 )

       if ( 0.0D+00 < dot_product ( plane(1:dim_num), r1(1:dim_num) ) ) then
           solid_angle = -area
       else
           solid_angle = area
       end if

       return
   end

   !--------------------------------------------------------------------------------
   !********************************************************************************
   !--------------------------------------------------------------------------------
   subroutine polygon_normal_3d ( n, v, normal )

       !*********************************************************************
       !
       !! POLYGON_NORMAL_3D computes the normal vector to a polygon in 3D.
       !
       !  Discussion:
       !
       !    If the polygon is planar, then this calculation is correct.
       !
       !    Otherwise, the normal vector calculated is the simple average
       !    of the normals defined by the planes of successive triples
       !    of vertices.
       !
       !    If the polygon is "almost" planar, this is still acceptable.
       !    But as the polygon is less and less planar, so this averaged normal
       !    vector becomes more and more meaningless.
       !
       !  Modified:
       !
       !    12 August 2005
       !
       !  Author:
       !
       !    John Burkardt
       !
       !  Reference:
       !
       !    Paulo Cezar Pinto Carvalho, Paulo Roma Cavalcanti,
       !    Point in Polyhedron Testing Using Spherical Polygons,
       !    in Graphics Gems V,
       !    edited by Alan Paeth,
       !    Academic Press, 1995, T385.G6975.
       !
       !  Parameters:
       !
       !    Input, integer N, the number of vertices.
       !
       !    Input, real ( kind = 8 ) V(3,N), the coordinates of the vertices.
       !
       !    Output, real ( kind = 8 ) NORMAL(3), the averaged normal vector
       !    to the polygon.
       !
       implicit none

       integer, parameter :: dim_num = 3
       integer n

       real ( kind = 8 ) r8vec_length
       integer j
       real ( kind = 8 ) normal(dim_num)
       real ( kind = 8 ) normal_norm
       real ( kind = 8 ) p(dim_num)
       real ( kind = 8 ) v(dim_num,n)
       real ( kind = 8 ) v1(dim_num)
       real ( kind = 8 ) v2(dim_num)

       normal(1:dim_num) = 0.0D+00

       v1(1:dim_num) = v(1:dim_num,2) - v(1:dim_num,1)

       do j = 3, n

           v2(1:dim_num) = v(1:dim_num,j) - v(1:dim_num,1)

           call r8vec_cross_3d ( v1, v2, p )

           normal(1:dim_num) = normal(1:dim_num) + p(1:dim_num)

           v1(1:dim_num) = v2(1:dim_num)

       end do
       !
       !  Normalize.
       !
       normal_norm = r8vec_length ( dim_num, normal )

       if ( normal_norm == 0.0D+00 ) then
           return
       end if

       normal(1:dim_num) = normal(1:dim_num) / normal_norm

       return
   end

   !--------------------------------------------------------------------------------
   !********************************************************************************
   !--------------------------------------------------------------------------------
   subroutine r8vec_cross_3d ( v1, v2, v3 )

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
       !    Input, real ( kind = 8 ) V1(3), V2(3), the two vectors.
       !
       !    Output, real ( kind = 8 ) V3(3), the cross product vector.
       !
       implicit none

       integer, parameter :: dim_num = 3

       real ( kind = 8 ) v1(dim_num)
       real ( kind = 8 ) v2(dim_num)
       real ( kind = 8 ) v3(dim_num)

       v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
       v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
       v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

       return
   end

   !--------------------------------------------------------------------------------
   !********************************************************************************
   !--------------------------------------------------------------------------------
   function r8vec_length ( dim_num, x )

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
       !    Input, real ( kind = 8 ) X(DIM_NUM), the vector.
       !
       !    Output, real ( kind = 8 ) R8VEC_LENGTH, the Euclidean length of the vector.
       !
       implicit none

       integer dim_num

       real ( kind = 8 ) r8vec_length
       real ( kind = 8 ) x(dim_num)

       r8vec_length = sqrt ( sum ( ( x(1:dim_num) )**2 ) )

       return
   end

   !--------------------------------------------------------------------------------
   !********************************************************************************
   !--------------------------------------------------------------------------------
   function r8vec_triple_product ( v1, v2, v3 )

       !*****************************************************************************80
       !
       !! R8VEC_TRIPLE_PRODUCT finds the triple product in 3D.
       !
       !  Discussion:
       !
       !    [A,B,C] = A dot ( B cross C )
       !            = B dot ( C cross A )
       !            = C dot ( A cross B )
       !
       !    The volume of a parallelepiped, whose sides are given by
       !    vectors A, B, and C, is abs ( A dot ( B cross C ) ).
       !
       !    Three vectors are coplanar if and only if their triple product vanishes.
       !
       !  Modified:
       !
       !    07 August 2005
       !
       !  Author:
       !
       !    John Burkardt
       !
       !  Reference:
       !
       !    Eric Weisstein,
       !    "Scalar Triple Product",
       !    CRC Concise Encyclopedia of Mathematics,
       !    CRC, 1999
       !
       !  Parameters:
       !
       !    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the vectors.
       !
       !    Output, real ( kind = 8 ) R8VEC_TRIPLE_PRODUCT, the triple product.
       !
       implicit none

       integer, parameter :: dim_num = 3

       real ( kind = 8 ) r8vec_triple_product
       real ( kind = 8 ) v1(dim_num)
       real ( kind = 8 ) v2(dim_num)
       real ( kind = 8 ) v3(dim_num)
       real ( kind = 8 ) v4(dim_num)

       call r8vec_cross_3d ( v2, v3, v4 )

       r8vec_triple_product = dot_product ( v1(1:dim_num), v4(1:dim_num) )

       return
   end

   !--------------------------------------------------------------------------------
   !********************************************************************************
   !--------------------------------------------------------------------------------
   function i4_wrap ( ival, ilo, ihi )

       !*****************************************************************************80
       !
       !! I4_WRAP forces an I4 to lie between given limits by wrapping.
       !
       !  Example:
       !
       !    ILO = 4, IHI = 8
       !
       !    I  I4_WRAP
       !
       !    -2     8
       !    -1     4
       !     0     5
       !     1     6
       !     2     7
       !     3     8
       !     4     4
       !     5     5
       !     6     6
       !     7     7
       !     8     8
       !     9     4
       !    10     5
       !    11     6
       !    12     7
       !    13     8
       !    14     4
       !
       !  Modified:
       !
       !    19 August 2003
       !
       !  Author:
       !
       !    John Burkardt
       !
       !  Parameters:
       !
       !    Input, integer IVAL, an integer value.
       !
       !    Input, integer ILO, IHI, the desired bounds for the integer value.
       !
       !    Output, integer I4_WRAP, a "wrapped" version of IVAL.
       !
       implicit none

       integer i4_modp
       integer i4_wrap
       integer ihi
       integer ilo
       integer ival
       integer jhi
       integer jlo
       integer wide

       jlo = min ( ilo, ihi )
       jhi = max ( ilo, ihi )

       wide = jhi - jlo + 1

       if ( wide == 1 ) then
           i4_wrap = jlo
       else
           i4_wrap = jlo + i4_modp ( ival - jlo, wide )
       end if

       return
   end

   !--------------------------------------------------------------------------------
   !********************************************************************************
   !--------------------------------------------------------------------------------
   function i4_modp ( i, j )

       !*****************************************************************************80
       !
       !! I4_MODP returns the nonnegative remainder of integer division.
       !
       !  Discussion:
       !
       !    If
       !      NREM = I4_MODP ( I, J )
       !      NMULT = ( I - NREM ) / J
       !    then
       !      I = J * NMULT + NREM
       !    where NREM is always nonnegative.
       !
       !    The MOD function computes a result with the same sign as the
       !    quantity being divided.  Thus, suppose you had an angle A,
       !    and you wanted to ensure that it was between 0 and 360.
       !    Then mod(A,360) would do, if A was positive, but if A
       !    was negative, your result would be between -360 and 0.
       !
       !    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
       !
       !  Examples:
       !
       !        I     J     MOD  I4_MODP    Factorization
       !
       !      107    50       7       7    107 =  2 *  50 + 7
       !      107   -50       7       7    107 = -2 * -50 + 7
       !     -107    50      -7      43   -107 = -3 *  50 + 43
       !     -107   -50      -7      43   -107 =  3 * -50 + 43
       !
       !  Modified:
       !
       !    02 March 1999
       !
       !  Author:
       !
       !    John Burkardt
       !
       !  Parameters:
       !
       !    Input, integer I, the number to be divided.
       !
       !    Input, integer J, the number that divides I.
       !
       !    Output, integer I4_MODP, the nonnegative remainder when I is
       !    divided by J.
       !
       implicit none

       integer i
       integer i4_modp
       integer j

       if ( j == 0 ) then
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) 'I4_MODP - Fatal error!'
           write ( *, '(a,i8)' ) '  I4_MODP ( I, J ) called with J = ', j
           stop
       end if

       i4_modp = mod ( i, j )

       if ( i4_modp < 0 ) then
           i4_modp = i4_modp + abs ( j )
       end if

       return
   end

   !--------------------------------------------------------------------------------
   !********************************************************************************
   !--------------------------------------------------------------------------------
   function arc_cosine ( c )

       !*****************************************************************************80
       !
       !! ARC_COSINE computes the arc cosine function, with argument truncation.
       !
       !  Discussion:
       !
       !    If you call your system ACOS routine with an input argument that is
       !    even slightly outside the range [-1.0, 1.0 ], you may get an unpleasant
       !    surprise (I did).
       !
       !    This routine simply truncates arguments outside the range.
       !
       !  Modified:
       !
       !    02 December 2000
       !
       !  Author:
       !
       !    John Burkardt
       !
       !  Parameters:
       !
       !    Input, real ( kind = 8 ) C, the argument.
       !
       !    Output, real ( kind = 8 ) ARC_COSINE, an angle whose cosine is C.
       !
       implicit none

       real ( kind = 8 ) arc_cosine
       real ( kind = 8 ) c
       real ( kind = 8 ) c2

       c2 = c
       c2 = max ( c2, -1.0D+00 )
       c2 = min ( c2, +1.0D+00 )

       arc_cosine = acos ( c2 )

       return
   end

   subroutine rotation2(sx,sy,sz,ipx,ipy,ipz,rot,rot_inverse)
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
       real*8 p0(3),p1(3),p2(3),p3(3)
       real*8 p0p(3),p1p(3),p2p(3),p3p(3)
       real*8 p1_tan(3),p2_tan(3),p3_tan(3)
       real*8 pN0(3),pNa(3)
       real*8 p_base(3),F,G,H
       real*8 normal_12(3),normal_1p2p(3),node_line(3)

       real*8 check1,check2,check3
       real*8 traslo_1,traslo_2,traslo_3
       real*8 ip1,ip2,ip3
       real*8 ipx,ipy,ipz
       real*8 s1,s2,s3
       real*8 sx,sy,sz

       real*8 rot(3,3),rot_inverse(3,3),det

       real*8 phi,theta,psi,pi,alfa

       real*8 c_phi,c_theta,c_psi
       real*8 s_phi,s_theta,s_psi

       real*8 angle_rad_3d

       real*8 p1_p1p,p2_p1p,p3_p1p
       real*8 p1_p3p,p2_p3p,p3_p3p
       real*8 pp1_p1p,pp2_p1p,pp3_p1p
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
       ip1 = ipx
       ip2 = ipz !ipz
       ip3 = ipy !ipy

       s1 = sx
       s2 = sz !sz
       s3 = sy !sy
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



   subroutine plane_exp_normal_3dbis ( p1,p2,p3, normal )
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




   function angle_rad_3d ( p1, p2, p3 )

       !*****************************************************************************80
       !
       !! ANGLE_RAD_3D returns the angle in radians between two rays in 3D.
       !
       !  Discussion:
       !
       !    The routine always computes the SMALLER of the two angles between
       !    two rays.  Thus, if the rays make an (exterior) angle of
       !    1.5 pi radians, the (interior) angle of 0.5 pi radians will be reported.
       !
       !    X dot Y = Norm(X) * Norm(Y) * Cos ( Angle(X,Y) )
       !
       !  Modified:
       !
       !    21 January 2005
       !
       !  Author:
       !
       !    John Burkardt
       !
       !  Parameters:
       !
       !    Input, real ( kind = 8 ) P1(3), P2(3), P3(3), points defining an angle.
       !    The rays are P1 - P2 and P3 - P2.
       !
       !    Output, real ( kind = 8 ) ANGLE_RAD_3D, the angle between the two rays,
       !    in radians.  This value will always be between 0 and PI.  If either ray has
       !    zero length, then the angle is returned as zero.
       !
       implicit none

       integer, parameter :: dim_num = 3

       real ( kind = 8 ) angle_rad_3d
       real ( kind = 8 ) arc_cosine
       real ( kind = 8 ) dot
       real ( kind = 8 ) p1(dim_num)
       real ( kind = 8 ) p2(dim_num)
       real ( kind = 8 ) p3(dim_num)
       real ( kind = 8 ) v1norm
       real ( kind = 8 ) v2norm

       v1norm = sqrt ( sum ( ( p1(1:dim_num) - p2(1:dim_num) )**2 ) )

       if ( v1norm == 0.0D+00 ) then
           angle_rad_3d = 0.0D+00
           return
       end if

       v2norm = sqrt ( sum ( ( p3(1:dim_num) - p2(1:dim_num) )**2 ) )

       if ( v2norm == 0.0D+00 ) then
           angle_rad_3d = 0.0D+00
           return
       end if

       dot = sum ( ( p1(1:dim_num) - p2(1:dim_num) ) &
           * ( p3(1:dim_num) - p2(1:dim_num) ) )

       angle_rad_3d = arc_cosine ( dot / ( v1norm * v2norm ) )

       return
   end







   !***********************************************************************
   !***********************************************************************
   !***********************************************************************
   subroutine plane_normal2imp_3d ( pp, normal, a, b, c, d )

       !*****************************************************************************80
       !
       !! PLANE_NORMAL2IMP_3D converts a normal form plane to implicit form in 3D.
       !
       !  Discussion:
       !
       !    The normal form of a plane in 3D is
       !
       !      PP, a point on the plane, and
       !      N, the unit normal to the plane.
       !
       !    The implicit form of a plane in 3D is
       !
       !      A * X + B * Y + C * Z + D = 0.
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
       !    Input, real ( kind = 8 ) PP(3), a point on the plane.
       !
       !    Input, real ( kind = 8 ) NORMAL(3), the unit normal vector to the plane.
       !
       !    Output, real ( kind = 8 ) A, B, C, D, the implicit plane parameters.
       !
       implicit none

       integer, parameter :: dim_num = 3

       real ( kind = 8 ) a
       real ( kind = 8 ) b
       real ( kind = 8 ) c
       real ( kind = 8 ) d
       real ( kind = 8 ) normal(dim_num)
       real ( kind = 8 ) pp(dim_num)

       a = normal(1)
       b = normal(2)
       c = normal(3)
       d = - a * pp(1) - b * pp(2) - c * pp(3)

       return
   end




   !***********************************************************************
   !***********************************************************************
   !***********************************************************************
   subroutine plane_exp2normal_3d ( p1, p2, p3, pp, normal )

       !*****************************************************************************80
       !
       !! PLANE_EXP2NORMAL_3D converts an explicit plane to normal form in 3D.
       !
       !  Discussion:
       !
       !    The explicit form of a plane in 3D is
       !
       !      the plane through P1, P2 and P3.
       !
       !    The normal form of a plane in 3D is
       !
       !      PP, a point on the plane, and
       !      N, the unit normal to the plane.
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
       !    Input, real ( kind = 8 ) P1(3), P2(3), P3(3), three points on the plane.
       !
       !    Output, real ( kind = 8 ) PP(3), a point on the plane.
       !
       !    Output, real ( kind = 8 ) NORMAL(3), a unit normal vector to the plane.
       !
       implicit none

       integer, parameter :: dim_num = 3

       real ( kind = 8 ) norm
       real ( kind = 8 ) normal(dim_num)
       real ( kind = 8 ) p1(dim_num)
       real ( kind = 8 ) p2(dim_num)
       real ( kind = 8 ) p3(dim_num)
       real ( kind = 8 ) pp(dim_num)

       pp(1:dim_num) = p1(1:dim_num)

       normal(1) = ( p2(2) - p1(2) ) * ( p3(3) - p1(3) ) &
           - ( p2(3) - p1(3) ) * ( p3(2) - p1(2) )

       normal(2) = ( p2(3) - p1(3) ) * ( p3(1) - p1(1) ) &
           - ( p2(1) - p1(1) ) * ( p3(3) - p1(3) )

       normal(3) = ( p2(1) - p1(1) ) * ( p3(2) - p1(2) ) &
           - ( p2(2) - p1(2) ) * ( p3(1) - p1(1) )

       norm = sqrt ( sum ( normal(1:dim_num)**2 ) )

       if ( norm == 0.0D+00 ) then
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) 'PLANE_EXP2NORMAL_3D - Fatal error!'
           write ( *, '(a)' ) '  The normal vector is null.'
           write ( *, '(a)' ) '  Two points coincide, or nearly so.'
           stop
       end if

       normal(1:dim_num) = normal(1:dim_num) / norm

       return
   end






   !***********************************************************************
   !***********************************************************************
   !***********************************************************************
   subroutine plane_normal2exp_3d ( pp, normal, p1, p2, p3 )

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
       !    Input, real ( kind = 8 ) PP(3), a point on the plane.
       !
       !    Input, real ( kind = 8 ) NORMAL(3), a normal vector N to the plane.  The
       !    vector must not have zero length, but it is not necessary for N
       !    to have unit length.
       !
       !    Output, real ( kind = 8 ) P1(3), P2(3), P3(3), three points on the plane.
       !
       implicit none

       integer, parameter :: dim_num = 3

       real ( kind = 8 ) normal(dim_num)
       real ( kind = 8 ) p1(dim_num)
       real ( kind = 8 ) p2(dim_num)
       real ( kind = 8 ) p3(dim_num)
       real ( kind = 8 ) pp(dim_num)
       real ( kind = 8 ) pq(dim_num)
       real ( kind = 8 ) pr(dim_num)

       call plane_normal_basis_3d ( pp, normal, pq, pr )

       p1(1:dim_num) = pp(1:dim_num)
       p2(1:dim_num) = pp(1:dim_num) + pq(1:dim_num)
       p3(1:dim_num) = pp(1:dim_num) + pr(1:dim_num)

       return
   end




   !***********************************************************************
   !***********************************************************************
   !***********************************************************************
   subroutine plane_normal_basis_3d ( pp, normal, pq, pr )

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
       !    Input, real ( kind = 8 ) PP(3), a point on the plane.  (Actually,
       !    we never need to know these values to do the calculation!)
       !
       !    Input, real ( kind = 8 ) NORMAL(3), a normal vector N to the plane.  The
       !    vector must not have zero length, but it is not necessary for N
       !    to have unit length.
       !
       !    Output, real ( kind = 8 ) PQ(3), a vector of unit length,
       !    perpendicular to the vector N and the vector PR.
       !
       !    Output, real ( kind = 8 ) PR(3), a vector of unit length,
       !    perpendicular to the vector N and the vector PQ.
       !
       implicit none

       integer, parameter :: dim_num = 3

       real ( kind = 8 ) r8vec_length
       real ( kind = 8 ) normal(dim_num)
       real ( kind = 8 ) normal_norm
       real ( kind = 8 ) pp(dim_num)
       real ( kind = 8 ) pq(dim_num)
       real ( kind = 8 ) pr(dim_num)
       real ( kind = 8 ) pr_norm
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




   !***********************************************************************
   !***********************************************************************
   !***********************************************************************
   subroutine r8vec_any_normal ( dim_num, v1, v2 )

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
       !    Input, real ( kind = 8 ) V1(DIM_NUM), the vector.
       !
       !    Output, real ( kind = 8 ) V2(DIM_NUM), a vector that is
       !    normal to V2, and has unit Euclidean length.
       !
       implicit none

       integer dim_num

       real ( kind = 8 ) r8vec_length
       integer i
       integer j
       integer k
       real ( kind = 8 ) v1(dim_num)
       real ( kind = 8 ) v2(dim_num)
       real ( kind = 8 ) vj
       real ( kind = 8 ) vk

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




   !***********************************************************************
   !***********************************************************************
   !***********************************************************************
   subroutine lines_exp_angle_3d ( p1, p2, q1, q2, angle )

       !*****************************************************************************80
       !
       !! LINES_EXP_ANGLE_3D finds the angle between two explicit lines in 3D.
       !
       !  Discussion:
       !
       !    The explicit form of a line in 3D is:
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
       !    Input, real ( kind = 8 ) P1(3), P2(3), two points on the first line.
       !
       !    Input, real ( kind = 8 ) Q1(3), Q2(3), two points on the second line.
       !
       !    Output, real ( kind = 8 ) ANGLE, the angle in radians between the two
       !    lines.  The angle is computed using the ACOS function, and so lies between
       !    0 and PI.  But if one of the lines is degenerate, the angle is
       !    returned as -1.0.
       !
       implicit none

       integer, parameter :: dim_num = 3

       real ( kind = 8 ) angle
       real ( kind = 8 ) arc_cosine
       real ( kind = 8 ) ctheta
       logical line_exp_is_degenerate_nd
       real ( kind = 8 ) p1(dim_num)
       real ( kind = 8 ) p2(dim_num)
       real ( kind = 8 ) pdotq
       real ( kind = 8 ) pnorm
       real ( kind = 8 ) q1(dim_num)
       real ( kind = 8 ) q2(dim_num)
       real ( kind = 8 ) qnorm

       if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
           !   write ( *, '(a)' ) ' '
           !   write ( *, '(a)' ) 'LINES_EXP_ANGLE_3D - Fatal error!'
           !   write ( *, '(a)' ) '  The line (P1,P2) is degenerate!'
           angle = -1.0D+00
           return
       end if

       if ( line_exp_is_degenerate_nd ( dim_num, q1, q2 ) ) then
           !   write ( *, '(a)' ) ' '
           !   write ( *, '(a)' ) 'LINES_EXP_ANGLE_3D - Fatal error!'
           !   write ( *, '(a)' ) '  The line (Q1,Q2) is degenerate!'
           angle = -1.0D+00
           return
       end if

       pnorm = sqrt ( sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 ) )

       qnorm = sqrt ( sum ( ( q2(1:dim_num) - q1(1:dim_num) )**2 ) )

       pdotq = sum ( ( p2(1:dim_num) - p1(1:dim_num) ) * ( q2(1:dim_num) - q1(1:dim_num) ) )

       ctheta = pdotq / ( pnorm * qnorm )

       angle = arc_cosine ( ctheta )

       return
   end




   !***********************************************************************
   !***********************************************************************
   !***********************************************************************
   subroutine r8mat_inverse_3d ( a, b, det )

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
       !    Input, real ( kind = 8 ) A(3,3), the matrix to be inverted.
       !
       !    Output, real ( kind = 8 ) B(3,3), the inverse of the matrix A.
       !
       !    Output, real ( kind = 8 ) DET, the determinant of the matrix A.
       !
       implicit none

       real ( kind = 8 ) a(3,3)
       real ( kind = 8 ) b(3,3)
       real ( kind = 8 ) det
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




   subroutine rotation_mat_vector_3d ( a, v, w )

       !*****************************************************************************80
       !
       !! ROTATION_MAT_VECTOR_3D applies a marix rotation to a vector in 3d.
       !
       !  Modified:
       !
       !    30 July 1999
       !
       !  Author:
       !
       !    John Burkardt
       !
       !  Parameters:
       !
       !    Input, real ( kind = 8 ) A(3,3), the matrix defining the rotation.
       !
       !    Input, real ( kind = 8 ) V(3), the vector to be rotated.
       !
       !    Output, real ( kind = 8 ) W(3), the rotated vector.
       !
       implicit none

       integer, parameter :: dim_num = 3

       real ( kind = 8 ) a(dim_num,dim_num)
       real ( kind = 8 ) v(dim_num)
       real ( kind = 8 ) w(dim_num)

       w(1:dim_num) = matmul ( a(1:dim_num,1:dim_num), v(1:dim_num) )

       return
   end
