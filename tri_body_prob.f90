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
!        procedure   ::  show_planet
    end type planet
    
    type triBody
        integer :: p_num
        type(planet), allocatable    ::  p(:)
      contains
        procedure   ::  energy
        procedure   ::  crush
        procedure   ::  modulus
!        procedure   ::  show_triBody
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

    interface triBody
        module procedure init_triBody
    end interface triBody

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
    
!    subroutine show_planet( obj )
!        class(planet)   ::  obj
!        print *, 'Planet Information:'
!        print '(a,f20.16)',  '  mass      : ', obj%m
!        print '(a,3f20.16)', '  coordinate: ', obj%x
!        print '(a,3f20.16)', '  velocity  : ', obj%v
!    end subroutine show_planet
    
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
    
    function init_triBody(obj_arr)
        type(triBody)  init_triBody
        type(planet)    ::  obj_arr(:)
        integer :: n, i
        n = size(obj_arr, 1)
        allocate(init_triBody%p(n))
        init_triBody%p_num = n
        do i=1, n
            call set_planet( init_triBody%p(i), obj_arr(i)%m, obj_arr(i)%x, obj_arr(i)%v )
        end do
     end function init_triBody
    
     !-------------------------------------------------------
     ! procedures associated with class triBody
     !-------------------------------------------------------

    function energy( triObj )
        class(triBody)  ::  triObj
        real*8          ::  energy
        real*8          ::  potential, kinetic
        integer :: i, j, n
        n = triObj%p_num
        potential = 0.d0
        do i=1, n
            do j=1, n
                if (i /= j) potential = potential + 0.5d0 * U_12(triObj%p(i), triObj%p(j))
            end do
        end do
        kinetic   = 0.d0
        do i=1, n
            kinetic = kinetic + triObj%p(i)%Ek()
        end do
        energy    = potential + kinetic
    end function energy
    
    function crush( triObj )
        class(triBody)   ::  triObj
        logical         ::  crush
        integer :: i,j,n
        n = triObj%p_num
        crush = .false.
        i = 1
        do while(i < n .and. crush .eqv. .false.)
            do j=i+1, n
                if ( collision(triObj%p(i), triObj%p(j) ) .eqv. .true. ) crush = .true.
                end do
        end do
    end function crush
    
    function modulus( triObj )
        class(triBody)  ::  triObj
        real*8          ::  modulus
        real*8          ::  m_m_square, m_x_square, m_v_square
        integer :: i, n
        n = triObj%p_num
        m_m_square = 0.d0
        m_x_square = 0.d0
        m_v_square = 0.d0
        do i=1, n
            m_m_square = m_m_square + metric(1)*((triObj%p(i)%m)**2)
            m_x_square = m_x_square + metric(2)*(dot_product( triObj%p(i)%x, triObj%p(i)%x ))
            m_v_square = m_v_square + metric(3)*(dot_product( triObj%p(i)%v, triObj%p(i)%v ))
        end do
        modulus = sqrt( m_m_square + m_x_square + m_v_square )
    end function modulus
    
!    subroutine show_triBody( triObj )
!        class(triBody)  ::  triObj
!        print *, "Tri-body System Information:"
!        print *, 'Planet a:'
!        print '(a,f20.16)',  '  mass      : ', triObj%a%m
!        print '(a,3f20.16)', '  coordinate: ', triObj%a%x
!        print '(a,3f20.16)', '  velocity  : ', triObj%a%v
!        print *, ' Planet b:'
!        print '(a,f20.16)',  '  mass      : ', triObj%b%m
!        print '(a,3f20.16)', '  coordinate: ', triObj%b%x
!        print '(a,3f20.16)', '  velocity  : ', triObj%b%v
!        print *, 'Planet c:'
!        print '(a,f20.16)',  '  mass      : ', triObj%c%m
!        print '(a,3f20.16)', '  coordinate: ', triObj%c%x
!        print '(a,3f20.16)', '  velocity  : ', triObj%c%v
!        if ( triObj%crush() .eqv. .false. ) then
!            print '(a,f20.16)', '  Total energy:', triObj%energy()
!        else
!            print *, 'Tri-Body system may crush. Energy is not well defined.'
!        end if
!    end subroutine show_triBody
    
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
        type(planet), allocatable :: temp_plist(:)
        integer :: n,i,j
        if (triObj1%p_num /= triObj2%p_num) print *, "p_num mismatch!"
        n = triObj1%p_num
        allocate(temp_plist(n))
        do i=1, n
            temp_plist = planet_plus( triObj1%p(i), triObj2%p(i) )
        end do
        triBody_plus = triBody(temp_plist)
    end function triBody_plus
    
    function triBody_subtract( triObj1, triObj2 )
        type(triBody), intent(in)   ::  triObj1, triObj2
        type(triBody)               ::  triBody_subtract
        type(planet), allocatable :: temp_plist(:)
        integer :: n,i,j
        if (triObj1%p_num /= triObj2%p_num) print *, "p_num mismatch!"
        n = triObj1%p_num
        allocate(temp_plist(n))
        do i=1, n
            temp_plist = planet_subtract( triObj1%p(i), triObj2%p(i) )
        end do
        triBody_subtract = triBody(temp_plist)
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
        type(planet), allocatable :: temp_plist(:)
        integer :: n,i,j
        n = triObj%p_num
        allocate(temp_plist(n))
        do i=1, n
            temp_plist(i) = planet_mul( coeff, triObj%p(i))
        end do
        triBody_mul = triBody(temp_plist)
    end function triBody_mul
    
    !------------------------------------------------------------
    ! main procedure characterizing evolution of a triBody system
    ! the whole lot of work is to construct this function properly
    !------------------------------------------------------------
    
    function fdt(t,state)
        type(triBody)    ::  state, fdt
        real*8           ::  t
        type(planet), allocatable :: temp_plist(:)
        integer :: n,i,j
        n = state%p_num
        allocate(temp_plist(n))
        do i=1, n
            temp_plist(i)%m = 0.d0
            temp_plist(i)%x = state%p(i)%v
            temp_plist(i)%v = 0.d0
            do j=1, n
                if (j /= i) temp_plist(i)%v = temp_plist(i)%v - state%p(j)%m * scaled_f_12( state%p(j), state%p(i))
            end do
        end do
        fdt = triBody(temp_plist)
    end function fdt
    
    !------------------------------------------------------------
    ! End of the module
    !------------------------------------------------------------
    
end module tri_body_prob