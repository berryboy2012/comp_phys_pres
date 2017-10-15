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
    real*16, dimension(2)    ::  xrange, yrange
    real*16                  ::  fps, movie_sec
    real*16                  ::  t0, t, h, delta_t
    real*16                  ::  err, tol
    type(triBody)           ::  state, iState
    real*16                  ::  cpu_cons(2)
    type(planet), allocatable     ::  temp_planet(:)

    integer                 ::  init_i, init_j, init_n
    !------------------------------------------------------------
    !------------------------------------------------------------



!     iState%a%m = 1.0q0              ! mass1
!     iState%a%x = [ 0.0q0, 1.0q0, 0.0q0 ]     ! coordinate1
!     iState%a%v = [ -1.0q0/2.0q0, 0.0q0, 0.0q0 ]     ! velocity1
!
!     iState%b%m = 1.0q0              ! mass2
!     iState%b%x = [ 0.0q0, -1.0q0, 0.0q0 ]     ! coordinate2
!     iState%b%v = [ 1.0q0/2.0q0, 0.0q0, 0.0q0 ]     ! velocity2
!
!     iState%c%m = 0.0q0              ! mass3
!     iState%c%x = [ 0.0q0, 0.0q0, 0.0q0 ]     ! coordinate3
!     iState%c%v = [ 0.0q0, 0.0q0, 0.0q0 ]     ! velocity3

    print *, 'Please specify the number of planets(input 0 to use samples instead):'
    read (*, *) init_n
    if (init_n /= 0) then
        allocate(temp_planet(init_n))
        print *, 'Input m_i, x_i, y_i, z_i, vx_i, vy_i, vz_i:'
        do init_i=1, init_n
            print *, 'Planet', init_i, ':'
            read (*, *) temp_planet(init_i)%m, temp_planet(init_i)%x, temp_planet(init_i)%v
        end do
    else
        print *, 'tell me which sample to use:'
        read (*, *) init_i
        call sample_temp(temp_planet, init_i)
    end if
    iState = triBody(temp_planet)

    print *, 'Please tell me how long you want to emulate(in seconds):'
    read (*, *) movie_sec
    fps = 60


    t0  = 0.0q0        ! 0.0q0
    h   = 1q0 / fps      ! 1.0q0/1000.0q0
    tol = 1.0e-16       ! relative tolerance
    delta_t = 1.q0 / fps

    skip    = int(1.q0 / (fps * h))     ! say, 10, namely make a frame per skip iterations
    N = int(fps * movie_sec * skip)

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
!    N = N * 100
    call cpu_time(cpu_cons(1))
    do i=1,N
        call rkf45(h,t,delta_t,state,err,tol)

        !h = 0.0001q0
        !call rkf45_step(t,h,state,err)
        !t = t + h

        !call euler(h,t,delta_t,state,err,tol)
        if ( mod(i,skip) == 0 ) then
            k = k+1
            !print *, k
            call write_state(data_unit,t,err,state)
!            call plot_state(plot_unit,state,k,frame_digs,plot_title,xrange,yrange)
        end if
    end do
    call cpu_time(cpu_cons(2))
!    print *, k
!     close(plot_unit)
    print *, 'CPU time elapsed: ', cpu_cons(2) - cpu_cons(1), ' sec.'
    close(data_unit)

    contains
        subroutine sample_temp(temp_planet, choice)
            type(planet), allocatable :: temp_planet(:)
            integer :: choice
            if (choice == 2) then
                allocate(temp_planet(3))
                temp_planet(1)%m = 2.0q0              ! mass1
                temp_planet(1)%x = [ 0.0q0, sqrt(3.0q0), 0.0q0 ]     ! coordinate1
                temp_planet(1)%v = [ -1.q0, 0.0q0, 0.0q0 ]     ! velocity1

                temp_planet(2)%m = 2.0q0              ! mass2
                temp_planet(2)%x = [ -1.q0, 0.q0, 0.0q0 ]     ! coordinate2
                temp_planet(2)%v = [ 0.5q0, -0.5q0 * sqrt(3.q0) , 0.0q0 ]
!        temp_planet(2)%v = [ 0.q0, 0.4q0, 0.0q0 ]     ! velocity2

                temp_planet(3)%m = 2.0q0              ! mass3
                temp_planet(3)%x = [ 1.q0, 0.q0, 0.0q0 ]     ! coordinate3
                temp_planet(3)%v = [ 0.5q0, 0.5q0 * sqrt(3.q0) , 0.0q0 ]     ! velocity3
!        temp_planet(3)%v = [ 0.q0, -0.4q0, 0.0q0 ]
            end if
            if (choice == 1) then
                allocate(temp_planet(3))
                temp_planet(1)%m = 1.0q0              ! mass1
                temp_planet(1)%x = [ 0.0q0, sqrt(3.0q0), 0.0q0 ]     ! coordinate1
                temp_planet(1)%v = [ -sqrt(1.q0/2.0q0), 0.0q0, 0.0q0 ]     ! velocity1

                temp_planet(2)%m = 1.0q0              ! mass2
                temp_planet(2)%x = [ -1.q0, 0.q0, 0.0q0 ]     ! coordinate2
                temp_planet(2)%v = [ sqrt(1.q0/2.0q0) * 0.5q0, -sqrt(1.q0/2.0q0) * sqrt(3.q0) / 2.q0, 0.0q0 ]     ! velocity2

                temp_planet(3)%m = 1.0q0              ! mass3
                temp_planet(3)%x = [ 1.q0, 0.q0, 0.0q0 ]     ! coordinate3
                temp_planet(3)%v = [ sqrt(1.q0/2.0q0) * 0.5q0, sqrt(1.q0/2.0q0) * sqrt(3.q0) / 2.q0, 0.0q0 ]     ! velocity3
            end if
            if (choice == 3) then
                allocate(temp_planet(3))
                temp_planet(1)%m = 1.0q0              ! mass1
                temp_planet(1)%x = [ 0.0q0, sqrt(3.0q0), 0.0q0 ]     ! coordinate1
                temp_planet(1)%v = [ -sqrt(1.q0/2.0q0), 0.0q0, 0.0q0 ]     ! velocity1

                temp_planet(2)%m = 1.0q0              ! mass2
                temp_planet(2)%x = [ -1.q0, 0.q0, 0.0q0 ]     ! coordinate2
                temp_planet(2)%v = [ sqrt(1.q0/2.0q0) * 0.5q0, -sqrt(1.q0/2.0q0) * sqrt(3.q0) / 2.q0, 0.0q0 ]     ! velocity2

                temp_planet(3)%m = 1.0q0              ! mass3
                temp_planet(3)%x = [ 1.q0, 1.q-32, 0.0q0 ]     ! coordinate3
                temp_planet(3)%v = [ sqrt(1.q0/2.0q0) * 0.5q0, sqrt(1.q0/2.0q0) * sqrt(3.q0) / 2.q0, 0.0q0 ]     ! velocity3
            end if
            if (choice == 4) then
                allocate(temp_planet(3))
                temp_planet(1)%m = 1.0q0              ! mass1
                temp_planet(1)%x = [ 0.0q0, sqrt(3.0q0), 0.0q0 ]     ! coordinate1
                temp_planet(1)%v = [ -sqrt(1.q0/2.0q0), 0.0q0, 0.0q0 ]     ! velocity1

                temp_planet(2)%m = 1.0q0              ! mass2
                temp_planet(2)%x = [ -1.q0, 0.q0, 0.0q0 ]     ! coordinate2
                temp_planet(2)%v = [ sqrt(1.q0/2.0q0) * 0.5q0, -sqrt(1.q0/2.0q0) * sqrt(3.q0) / 2.q0, 0.0q0 ]     ! velocity2

                temp_planet(3)%m = 1.0q0              ! mass3
                temp_planet(3)%x = [ 1.q0, 1.q-33, 0.0q0 ]     ! coordinate3
                temp_planet(3)%v = [ sqrt(1.q0/2.0q0) * 0.5q0, sqrt(1.q0/2.0q0) * sqrt(3.q0) / 2.q0, 0.0q0 ]     ! velocity3
            end if


        end subroutine sample_temp

end program