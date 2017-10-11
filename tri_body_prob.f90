module tri_body_prob
    
    implicit none
    
    !-------------------------------------------------------
    ! definition of universal constants
    !-------------------------------------------------------
    
    real*8, parameter   ::  G  = 1.0d0  ! gravitation constant
    real*8, parameter   ::  pi = 3.1415926535897932384626433d0
    real*8, parameter   ::  forbidden_dist = 1.0e-14
    real*8, parameter   ::  metric(3) = [1.0d0,1.0d0,1.0d0]
    
    !-------------------------------------------------------
    ! definition of the basic type characterizing the problem
    !-------------------------------------------------------
    
    type planet
        real*8                ::  m
        real*8, dimension(3)  ::  x, v
      contains
        procedure   ::  Ek
        procedure   ::  x_val
        procedure   ::  v_val
        procedure   ::  show_planet
    end type planet
    
    type triBody
        type(planet)    ::  a, b, c
      contains
        procedure   ::  energy
        procedure   ::  crush
        procedure   ::  modulus
        procedure   ::  show_triBody
    end type triBody
    
    interface operator(+)
        procedure   ::  triBody_plus, planet_plus
    end interface
    
    interface operator(-)
        procedure   ::  triBody_subtract, planet_subtract
    end interface
    
    interface operator(*)
        procedure   ::  triBody_mul, planet_mul
    end interface
    
  contains
     
     !-------------------------------------------------------
     ! procedures associated with class planet
     !-------------------------------------------------------
     
    function Ek( obj )
        class(planet)   ::  obj
        real*8          ::  Ek
        Ek = ( obj%m/2.0d0 ) * dot_product( obj%v, obj%v )
    end function Ek
    
    function x_val( obj )
        class(planet)   ::  obj
        real*8          ::  x_val
        x_val = sqrt( dot_product( obj%x, obj%x ) ) 
    end function x_val
     
    function v_val( obj )
        class(planet)   ::  obj
        real*8          ::  v_val
        v_val = sqrt( dot_product( obj%v, obj%v ) ) 
    end function v_val
    
    subroutine show_planet( obj )
        class(planet)   ::  obj
        print *, 'Planet Information:'
        print '(a,f20.16)',  '  mass      : ', obj%m
        print '(a,3f20.16)', '  coordinate: ', obj%x
        print '(a,3f20.16)', '  velocity  : ', obj%v
    end subroutine show_planet
    
    !-------------------------------------------------------
    ! functions characterizing two body interaction
    !-------------------------------------------------------
    
    function r_12( obj1, obj2 )
        type(planet)    ::  obj1, obj2
        real*8          ::  r_12
        r_12 = sqrt( dot_product( obj1%x-obj2%x, obj1%x-obj2%x ) )
    end function r_12
    
    function U_12( obj1, obj2 )
        type(planet)    ::  obj1, obj2
        real*8          ::  U_12
        U_12 = - G*obj1%m*obj2%m/r_12( obj1,obj2 )
    end function U_12
    
    function scaled_f_12( obj1, obj2 )
        ! vector originating from obj1 to obj2, hence renders 2-1
        ! the real force is f_12 = m1*m2*scaled_f_12
        real*8, dimension(3)    :: scaled_f_12
        type(planet)            :: obj1, obj2
        real*8                  :: r             
        r = r_12( obj1,obj2 )
        scaled_f_12 = ( G/r**3 )*( obj2%x - obj1%x )
    end function scaled_f_12
    
    function collision( obj1, obj2 )
        type(planet)    ::  obj1, obj2
        logical         ::  collision
        real*8          ::  r
        collision = .false.
        r = r_12( obj1,obj2 )
        if ( r <= forbidden_dist ) then
            collision = .true.
        end if
    end function collision
    
    !--------------------------------------------------------
    ! set elements of the class instances
    !--------------------------------------------------------
    
    subroutine set_planet( obj, m, x, v )
        type(planet)            ::  obj
        real*8                  ::  m
        real*8, dimension(3)    ::  x, v
        obj%m = m
        obj%x = x
        obj%v = v
    end subroutine set_planet
    
    subroutine set_triBody( triObj, obj1, obj2, obj3 )
        type(triBody)  ::  triObj
        type(planet)    ::  obj1, obj2, obj3
        call set_planet( triObj%a, obj1%m, obj1%x, obj1%v )
        call set_planet( triObj%b, obj2%m, obj2%x, obj2%v )
        call set_planet( triObj%c, obj3%m, obj3%x, obj3%v )
     end subroutine set_triBody
    
     !-------------------------------------------------------
     ! procedures associated with class triBody
     !-------------------------------------------------------

    function energy( triObj )
        class(triBody)  ::  triObj
        real*8          ::  energy
        real*8          ::  potential, kinetic
        potential = 0.50d0*( U_12( triObj%a, triObj%b )&
                    + U_12( triObj%b, triObj%c )&
                    + U_12( triObj%a, triObj%c ) )
        kinetic   = triObj%a%Ek() + triObj%b%Ek() + triObj%c%Ek()
        energy    = potential + kinetic
    end function energy
    
    function crush( triObj )
        class(triBody)   ::  triObj
        logical         ::  crush
        crush = .false.
        if ( collision(triObj%a, triObj%b ) .eqv. .true. ) then
            crush = .true.
        else
        if ( collision(triObj%a, triObj%c ) .eqv. .true. ) then
            crush = .true.
        else
        if ( collision(triObj%b, triObj%c ) .eqv. .true. ) then
            crush = .true.
        end if
        end if
        end if
    end function crush
    
    function modulus( triObj )
        class(triBody)  ::  triObj
        real*8          ::  modulus
        real*8          ::  m_m_square, m_x_square, m_v_square
        m_m_square = metric(1)*( (triObj%a%m)**2 &
                                + (triObj%b%m)**2 &
                                + (triObj%c%m)**2 )
        m_x_square = metric(2)*( dot_product( triObj%a%x, triObj%a%x )&
                        + dot_product( triObj%b%x, triObj%b%x )&
                        + dot_product( triObj%c%x, triObj%c%x ) )
        m_v_square = metric(3)*( dot_product( triObj%a%v, triObj%a%v )&
                        + dot_product( triObj%b%v, triObj%b%v )&
                        + dot_product( triObj%c%v, triObj%c%v ) )
        modulus = sqrt( m_m_square + m_x_square + m_v_square )
    end function modulus
    
    subroutine show_triBody( triObj )
        class(triBody)  ::  triObj
        print *, "Tri-body System Information:"
        print *, 'Planet a:'
        print '(a,f20.16)',  '  mass      : ', triObj%a%m
        print '(a,3f20.16)', '  coordinate: ', triObj%a%x
        print '(a,3f20.16)', '  velocity  : ', triObj%a%v
        print *, ' Planet b:'
        print '(a,f20.16)',  '  mass      : ', triObj%b%m
        print '(a,3f20.16)', '  coordinate: ', triObj%b%x
        print '(a,3f20.16)', '  velocity  : ', triObj%b%v
        print *, 'Planet c:'
        print '(a,f20.16)',  '  mass      : ', triObj%c%m
        print '(a,3f20.16)', '  coordinate: ', triObj%c%x
        print '(a,3f20.16)', '  velocity  : ', triObj%c%v
        if ( triObj%crush() .eqv. .false. ) then
            print '(a,f20.16)', '  Total energy:', triObj%energy()
        else
            print *, 'Tri-Body system may crush. Energy is not well defined.'
        end if
    end subroutine show_triBody
    
    !------------------------------------------------------------
    ! overload of operator + and -, namely plus and subtract
    !------------------------------------------------------------
    
    function planet_plus( obj1, obj2 )
        type(planet), intent(in)    ::  obj1, obj2
        type(planet)                ::  planet_plus
        planet_plus%m = obj1%m + obj2%m
        planet_plus%x = obj1%x + obj2%x
        planet_plus%v = obj1%v + obj2%v
    end function planet_plus
    
    function planet_subtract( obj1, obj2 )
        type(planet), intent(in)    ::  obj1, obj2
        type(planet)                ::  planet_subtract
        planet_subtract%m = obj1%m - obj2%m
        planet_subtract%x = obj1%x - obj2%x
        planet_subtract%v = obj1%v - obj2%v
    end function planet_subtract
    
    function triBody_plus( triObj1, triObj2 )
        type(triBody), intent(in)   ::  triObj1, triObj2
        type(triBody)               ::  triBody_plus
        triBody_plus%a = planet_plus( triObj1%a, triObj2%a )
        triBody_plus%b = planet_plus( triObj1%b, triObj2%b )
        triBody_plus%c = planet_plus( triObj1%c, triObj2%c )
    end function triBody_plus
    
    function triBody_subtract( triObj1, triObj2 )
        type(triBody), intent(in)   ::  triObj1, triObj2
        type(triBody)               ::  triBody_subtract
        triBody_subtract%a = planet_subtract( triObj1%a, triObj2%a )
        triBody_subtract%b = planet_subtract( triObj1%b, triObj2%b )
        triBody_subtract%c = planet_subtract( triObj1%c, triObj2%c )
    end function triBody_subtract
    
    function planet_mul( coeff, obj )
        type(planet), intent(in)    ::  obj
        real*8, intent(in)          ::  coeff
        type(planet)                ::  planet_mul
        planet_mul%m = coeff*obj%m
        planet_mul%x = coeff*obj%x
        planet_mul%v = coeff*obj%v
    end function planet_mul
    
    function triBody_mul( coeff, triObj )
        type(triBody), intent(in)   ::  triObj
        real*8, intent(in)          ::  coeff
        type(triBody)               ::  triBody_mul
        triBody_mul%a = planet_mul( coeff, triObj%a )
        triBody_mul%b = planet_mul( coeff, triObj%b )
        triBody_mul%c = planet_mul( coeff, triObj%c )
    end function triBody_mul
    
    !------------------------------------------------------------
    ! main procedure characterizing evolution of a triBody system
    ! the whole lot of work is to construct this function properly
    !------------------------------------------------------------
    
    function fdt(t,state)
        type(triBody)    ::  state, fdt
        real*8           ::  t
        fdt%a%m = 0.0d0
        fdt%a%x = state%a%v
        fdt%a%v = -state%b%m * scaled_f_12( state%b, state%a ) &
                  -state%c%m * scaled_f_12( state%c, state%a )
        fdt%b%m = 0.0d0
        fdt%b%x = state%b%v
        fdt%b%v = -state%a%m * scaled_f_12( state%a, state%b ) &
                  -state%c%m * scaled_f_12( state%c, state%b )
        fdt%c%m = 0.0d0
        fdt%c%x = state%c%v
        fdt%c%v = -state%a%m * scaled_f_12( state%a, state%c ) &
                  -state%b%m * scaled_f_12( state%b, state%c )
    end function fdt
    
    !------------------------------------------------------------
    ! End of the module
    !------------------------------------------------------------
    
end module tri_body_prob