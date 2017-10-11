program main

    use tri_body_prob
    use solver
    use ploter
    implicit none

    !-----------------------------------------------------
    !-----------------------------------------------------

    integer                 ::  N, skip, frame_digs
    integer                 ::  i, k
    integer                 ::  data_unit, plot_unit
    character(50)           ::  data_file, plot_file
    character(255)          ::  plot_title
    real*8, dimension(2)    ::  xrange, yrange
    real*8                  ::  fps, movie_sec
    real*8                  ::  t0, t, h, delta_t
    real*8                  ::  err, tol
    type(triBody)           ::  state, iState
    type(planet), allocatable     ::  temp_planet(:)

    !------------------------------------------------------------
    !------------------------------------------------------------

    fps = 30
    movie_sec = 4


    t0  = 0.0d0        ! 0.0d0
    h   = 0.001d0      ! 1.0d0/1000.0d0
    tol = 1.0e-8       ! relative tolerance
    delta_t = 0.01

    skip    = int(1.d0 / (fps * delta_t))     ! say, 10, namely make a frame per skip iterations
    N = int(fps * movie_sec * skip)

!     iState%a%m = 1.0d0              ! mass1
!     iState%a%x = [ 0.0d0, 1.0d0, 0.0d0 ]     ! coordinate1
!     iState%a%v = [ -1.0d0/2.0d0, 0.0d0, 0.0d0 ]     ! velocity1
!
!     iState%b%m = 1.0d0              ! mass2
!     iState%b%x = [ 0.0d0, -1.0d0, 0.0d0 ]     ! coordinate2
!     iState%b%v = [ 1.0d0/2.0d0, 0.0d0, 0.0d0 ]     ! velocity2
!
!     iState%c%m = 0.0d0              ! mass3
!     iState%c%x = [ 0.0d0, 0.0d0, 0.0d0 ]     ! coordinate3
!     iState%c%v = [ 0.0d0, 0.0d0, 0.0d0 ]     ! velocity3
    allocate(temp_planet(3))
    temp_planet(1)%m = 1.0d0              ! mass1
    temp_planet(1)%x = [ 0.0d0, sqrt(3.0d0), 0.0d0 ]     ! coordinate1
    temp_planet(1)%v = [ -sqrt(1.d0/2.0d0), 0.0d0, 0.0d0 ]     ! velocity1

    temp_planet(2)%m = 1.0d0              ! mass2
    temp_planet(2)%x = [ -1.d0, 0.d0, 0.0d0 ]     ! coordinate2
    temp_planet(2)%v = [ sqrt(1.d0/2.0d0) * 0.5d0, -sqrt(1.d0/2.0d0) * sqrt(3.d0) / 2.d0, 0.0d0 ]     ! velocity2

    temp_planet(3)%m = 1.0d0              ! mass3
    temp_planet(3)%x = [ 1.d0, 0.d0, 0.0d0 ]     ! coordinate3
    temp_planet(3)%v = [ sqrt(1.d0/2.0d0) * 0.5d0, sqrt(1.d0/2.0d0) * sqrt(3.d0) / 2.d0, 0.0d0 ]     ! velocity3
    iState = triBody(temp_planet)

    xrange  =   [ -5, 5 ]        ! [ -10, 10 ]
    yrange  =   [ -5, 5 ]        ! [ -10, 10 ]

    data_unit = 35
    plot_unit = 6

    data_file = 'data.dat'
    plot_file = 'plot.gnu'
    plot_title = 'Animation of three planets'

    !--------------------------------------------------------------
    !--------------------------------------------------------------

    frame_digs = int( log10( dble(N)/dble(skip) ) ) + 1

    open(unit=data_unit,file=data_file,status='replace')
!     open(unit=plot_unit,file=plot_file,status='replace')

    t = t0
    state = iState
    k = 0
    err = 0
    call write_state(data_unit,t,err,state)
!    call plot_state(plot_unit,state,k,frame_digs,plot_title,xrange,yrange)
    do i=1,N
        call rkf45(h,t,delta_t,state,err,tol)
        call write_state(data_unit,t,err,state)
        if ( mod(i,skip) == 0 ) then
            k = k+1
!            call plot_state(plot_unit,state,k,frame_digs,plot_title,xrange,yrange)
        end if
    end do

!     close(plot_unit)
    close(data_unit)

end program