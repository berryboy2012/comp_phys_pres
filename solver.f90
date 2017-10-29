module solver

    use tri_body_prob
    implicit none

    !------------------------------------------------------
    ! constants used in rkf45 method
    !------------------------------------------------------

    real*16, parameter, dimension(6)     ::  rkf_ct = &
        [0.q0, 1.0q0/4.0q0, 3.0q0/8.0q0, 12.0q0/13.0q0,  1.0q0, 1.0q0/2.0q0]

    real*16, parameter, dimension(5,5)   ::  rkf_csu = transpose( reshape( [&
        &      1.0q0,       0.0q0,      0.0q0,     0.0q0,      0.0q0, &
        &      3.0q0,       9.0q0,      0.0q0,     0.0q0,      0.0q0, &
        &   1932.0q0,   -7200.0q0,   7296.0q0,     0.0q0,      0.0q0, &
        &   8341.0q0,  -32832.0q0,  29440.0q0,  -845.0q0,      0.0q0, &
        &  -6080.0q0,   41040.0q0, -28352.0q0,  9295.0q0,  -5643.0q0  &
        &],[5,5] ) )

    real*16, parameter, dimension(5)     ::  rkf_csl = &
        [ 4.0q0, 32.0q0, 2197.0q0, 4104.0q0, 20520.0q0 ]

    real*16, parameter, dimension(6)     ::  rkf_cw = &
        [    25.0q0/216.0q0,          0.0q0, 1408.0q0/2565.0q0, &
          2197.0q0/4104.0q0,   -1.0q0/5.0q0,             0.0q0 ]

    real*16, parameter, dimension(6)     ::  rkf_cz = &
        [    16.0q0/135.0q0,          0.0q0, 6656.0q0/12825.0q0, &
        28561.0q0/56430.0q0,  -9.0q0/50.0q0,       2.0q0/55.0q0 ]

    real*16, parameter, dimension(6)     ::  rkf_cerr = &
        [     1.0q0/360.0q0,          0.0q0, -128.0q0/4275.0q0, &
        -2197.0q0/75240.0q0,   1.0q0/50.0q0,      2.0q0/55.0q0 ]

  contains

    !************************************************************
    !the main subprocedures in rkf45 method
    !************************************************************

    !------------------------------------------------------------
    !the core iteration part of rkf45 method
    !------------------------------------------------------------

    subroutine rkf45_step(t,h,w,err)
        real*16, intent(in)          ::  t, h
        type(triBody)               ::  w, w1, z1
        real*16                      ::  err
        type(triBody), dimension(6) ::  s
        integer                     ::  i, j
        do i=1,6
            w1 = w
            do j=1,i - 1
                w1 = w1 + ( rkf_csu(i - 1,j)/rkf_csl(i - 1) )*h*s(j)
            end do
            s(i) = fdt(w1)
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

        real*16          ::  h, dt, mindt
        real*16          ::  t, ct, tt, delta_t
        type(triBody)   ::  state, cstate, tstate
        real*16          ::  err, cerr, terr, tol
        integer         ::  its, max_devisions, p

        max_devisions = 10
        mindt = 2.0Q-9
        p = 4

        ct = t
        cstate = state
        cerr = 0.0q0

        dt = h
        tt = ct
        tstate = cstate
        call rkf45_step(tt,dt,tstate,terr)

        do      ! main loop
            dt = 0.84q0*( (tol/terr)**(1/dble(p)) )*dt
            its = 0
            do  ! subloop
                tt = ct
                tstate = cstate
                call rkf45_step(tt,dt,tstate,terr)
                if ( terr > tol .and. dt > mindt .and. its < max_devisions) then
                    its = its+1
                    if(its==1) then
                        dt = 0.8*( (tol/terr)**(1/dble(p+1)) )*dt
                    else
                        dt = dt/2.0q0
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
        state = cstate
        err = cerr

    end subroutine rkf45

    subroutine euler(h,t,delta_t,state,err,tol)
        ! h: inital step length
        ! t: time
        ! delta_t: desired interval
        ! state: state of systems
        ! err: error
        ! tol: tolerance of relative error
        real*16          ::  h, dt, mindt
        real*16          ::  t, ct, tt, delta_t
        type(triBody)   ::  state, cstate, tstate
        real*16          ::  err, cerr, terr, tol
        integer         ::  its, max_its, p
        dt = h
        dt = 0.0000001q0
        max_its = int(delta_t / dt)
        do its=1, max_its
            do p=1, state%p_num
                state%p(p)%x = state%p(p)%x + dt * state%p(p)%v + 0.5q0 * dt ** 2.q0 * state_acc(state, p)
                state%p(p)%v = state%p(p)%v +  dt * state_acc(state, p)
            end do
            t = t + dt
        end do
    end subroutine euler

    function state_acc(state, p)
        type(triBody) :: state
        real*16 :: state_acc(3)
        integer :: i,p
        state_acc = 0.q0
        do i=1, state%p_num
            if (i /= p) state_acc = state_acc - state%p(i)%m * scaled_f_12(state%p(i), state%p(p))
        end do
    end function state_acc

end module