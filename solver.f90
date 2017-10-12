module solver

    use tri_body_prob
    implicit none

    !------------------------------------------------------
    ! constants used in rkf45 method
    !------------------------------------------------------

    real*8, parameter, dimension(5)     ::  rkf_ct = &
        [1.0d0/4.0d0, 3.0d0/8.0d0, 12.0d0/13.0d0,  1.0d0, 1.0d0/2.0d0]

    real*8, parameter, dimension(5,5)   ::  rkf_csu = transpose( reshape( [&
        &      1.0d0,       0.0d0,      0.0d0,     0.0d0,      0.0d0, &
        &      3.0d0,       9.0d0,      0.0d0,     0.0d0,      0.0d0, &
        &   1932.0d0,   -7200.0d0,   7296.0d0,     0.0d0,      0.0d0, &
        &   8341.0d0,  -32832.0d0,  29440.0d0,  -845.0d0,      0.0d0, &
        &  -6080.0d0,   41040.0d0, -28352.0d0,  9295.0d0,  -5643.0d0  &
        &],[5,5] ) )

    real*8, parameter, dimension(5)     ::  rkf_csl = &
        [ 4.0d0, 32.0d0, 2197.0d0, 4104.0d0, 20520.0d0 ]

    real*8, parameter, dimension(6)     ::  rkf_cw = &
        [    25.0d0/216.0d0,          0.0d0, 1408.0d0/2565.0d0, &
          2197.0d0/4104.0d0,   -1.0d0/5.0d0,             0.0d0 ]

    real*8, parameter, dimension(6)     ::  rkf_cz = &
        [    16.0d0/135.0d0,          0.0d0, 6656.0d0/12825.0d0, &
        28561.0d0/56430.0d0,  -9.0d0/50.0d0,       2.0d0/55.0d0 ]

    real*8, parameter, dimension(6)     ::  rkf_cerr = &
        [     1.0d0/360.0d0,          0.0d0, -128.0d0/4275.0d0, &
        -2197.0d0/75240.0d0,   1.0d0/50.0d0,      2.0d0/55.0d0 ]

  contains

    !************************************************************
    !the main subprocedures in rkf45 method
    !************************************************************

    !------------------------------------------------------------
    !the core iteration part of rkf45 method
    !------------------------------------------------------------

    subroutine rkf45_step(t,h,w,err)

        real*8, intent(in)          ::  t, h
        type(triBody)               ::  w, w1, z1
        real*8                      ::  err
        type(triBody), dimension(6) ::  s
        integer                     ::  i, j

        s(1) = fdt( t,w )
        ! w1 currently used as a temporary variable
        do i=1,5
            w1 = w
            do j=1,i
                w1 = w1 + ( rkf_csu(i,j)/rkf_csl(j) )*h*s(j)
            end do
            s(i+1) = fdt( t,w1 )
        end do

        w1 = w
        do j=1,6
            w1 = w1 + ( rkf_cw(j) )*h*s(j)
        end do

        z1 = w
        do j=1,6
            z1 = z1 + ( rkf_cz(j) )*h*s(j)
        end do

        err = modulus(z1-w1)/modulus(z1)
        w = z1

    end subroutine rkf45_step

    !--------------------------------------------------------
    ! integration over a interval with multi-steps
    !--------------------------------------------------------

    subroutine rkf45(h,t,delta_t,state,err,tol)

        real*8          ::  h, dt, mindt
        real*8          ::  t, ct, tt, delta_t
        type(triBody)   ::  state, cstate, tstate
        real*8          ::  err, cerr, terr, tol
        integer         ::  its, max_devisions, p

        max_devisions = 1000
        mindt = 2*epsilon(h)

        ct = t
        cstate = state
        cerr = 0.0d0

        dt = h
        tt = ct
        tstate = cstate
        call rkf45_step(tt,dt,tstate,terr)

        do      ! main loop
            dt = 0.8*( (tol/terr)**(1/dble(p+1)) )*dt
            its = 0
            do  ! subloop
                tt = ct
                tstate = cstate
                call rkf45_step(tt,dt,tstate,terr)
                if ( terr > tol .and. dt > 2.0d0*mindt .and. its < max_devisions) then
                    its = its+1
                    if(its==1) then
                        dt = 0.8*( (tol/terr)**(1/dble(p+1)) )*dt
                    else
                        dt = dt/2.0d0
                    end if
                else
                    exit
                end if
            end do
            ct = ct+dt
            cstate = tstate
            if(terr>cerr) then
                cerr = terr
            end if
            !print *, its
            !print *, ct
            if ( ct > t+delta_t ) then !!!!!!!!!!!
                dt = t+delta_t-tt
                call rkf45_step(tt,dt,tstate,terr)
                cstate = tstate
                if(terr>cerr) then
                    cerr = terr
                end if
                exit
            end if
        end do

        t = t+delta_t
        !print *, t
        state = cstate
        err = cerr

    end subroutine rkf45

end module