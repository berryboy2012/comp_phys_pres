module solver_new
    use tri_body_prob
    implicit none
    real*8, parameter :: min_adj_dt = 1.d-6
    real*8, parameter, dimension(6) :: rkf_tki = &
            [0.d0, 1.d0 / 4.d0, 3.d0 / 8.d0, 12.d0 / 13.d0, 1.d0, 1.d0 / 2.d0]
    real*8, parameter, dimension(6, 6) :: rkf_yki = transpose( reshape( [&
            0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
            1.d0 / 4.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
            3.d0 / 32.d0, 9.d0 / 32.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
            1932.d0 / 2197.d0, -7200.d0 / 2197.d0, 7296.d0 / 2197.d0, 0.d0, 0.d0, 0.d0, &
            439.d0 / 216.d0, -8.d0, 3680.d0 / 513.d0, -845.d0 / 4104.d0, 0.d0, 0.d0, &
            -8.d0 / 27.d0, 2.d0, -3544.d0 / 2565.d0, 1859.d0 / 4104.d0, -11.d0 / -40.d0, 0.d0 &
    ], [6, 6]))
    real*8, parameter, dimension(6) :: rkf_yk = &
            [25.d0 / 216.d0, 0.d0, 1408.d0 / 2565.d0, 2197.d0 / 4101.d0, -1.d0 / 5.d0, 0.d0]
    real*8, parameter, dimension(6) :: rkf_zk = &
            [16.d0 / 135.d0, 0.d0, 6656.d0 / 12825.d0, 28561.d0 / 56430.d0, -9.d0 / 50.d0, 2.d0 / 55.d0]
    integer :: loop_depth=1
    contains
    subroutine rk5(h, w, err)
        real*8, intent(in) :: h
        type(triBody) :: w, yl
        real*8 :: err
        type(triBody), dimension(6) :: k
        integer :: i, j

        do i=1, 6
            k(i) = w
            do j=1, i-1
                k(i) = k(i) + rkf_yki(i, j) * h * k(j)
            end do
            k(i) = fdt(k(i))
            !call k(i)%show_triBody()
        end do
        yl = w
        do i=1, 5
            yl = yl + rkf_yk(i) * h * k(i)
        end do
        do i=1, 6
            w = w + rkf_zk(i) * h * k(i)
        end do
        err = modulus(w - yl) / modulus(w)
        !call w%show_triBody()
        !print *, 'err in rk5',  err
    end
    recursive subroutine rkf45(t, delta_t, state, err, tol)
        real*8 :: adjusted_dt, err_tmp
        real*8, intent(inout) :: t, err
        real*8, intent(in) :: delta_t, tol
        type(triBody) :: state, test_state
        test_state = state
        call rk5(delta_t, test_state, err_tmp)
        if ((isnan(err_tmp) .eqv. .true.) .or. (err_tmp < tol)) then
            !print *, 'passed on first try'
            state = test_state
            t = t + delta_t
            err = err + err_tmp
        else
            !print *, 'F'
            adjusted_dt = delta_t * 0.84d0 * (tol / err_tmp)**(0.25d0)
            !print *, adjusted_dt / delta_t * 100.d0, '%'
            if (adjusted_dt < min_adj_dt) then
                !print *, 'adjusted length too short'
                call rk5(adjusted_dt, state, err_tmp)
                t = t + adjusted_dt
                err = err + err_tmp
            else
                !print *, 'L'
                call rkf45(t, adjusted_dt, state, err, tol)
            end if
            if (min_adj_dt > delta_t - adjusted_dt) then
                !print *, 'length left too short'
                call rk5(delta_t - adjusted_dt, state, err_tmp)
                t = t + delta_t - adjusted_dt
                err = err + err_tmp
            else
                !print *, 'R'
                call rkf45(t, delta_t - adjusted_dt, state, err, tol)
            end if
        end if
        !print *, 'EOF'
    end
    subroutine rkf45loop(t, delta_t, state, err, tol)
        implicit none
        real*8 :: adjusted_dt, err_tmp
        real*8, intent(inout) :: t, err
        real*8, intent(in) :: delta_t, tol
        type(triBody) :: state, test_state
        integer :: i
        adjusted_dt = delta_t / loop_depth
        err = 2.d0 * tol
        do while ((isnan(err) .eqv. .false.) .and. (err > tol) .and. adjusted_dt > min_adj_dt)
            err = 0.d0
            test_state = state
            do i=1, loop_depth
                call rk5(adjusted_dt, test_state, err_tmp)
                err = err + err_tmp
            end do
            adjusted_dt = adjusted_dt / 2.d0
            loop_depth = loop_depth * 2
            !print *, err
        end do
        if (tol > 8.d0 * err .and. loop_depth > 4) loop_depth = loop_depth / 4
        if (adjusted_dt < min_adj_dt) loop_depth = int(2.d0 * log(delta_t / min_adj_dt) / log(2.d0))
        state = test_state
        t = t + delta_t
    end
    subroutine rkf45_iter(t, delta_t, state, err, tol)
        real*8 :: adjusted_dt
        real*8, intent(inout) :: t, err
        real*8, intent(in) :: delta_t, tol
        type(triBody) :: state, test_state
        test_state = state
        call rk5(delta_t, test_state, err)
        adjusted_dt = delta_t
        do while ((isnan(err) .eqv. .false.) .and. (err > tol) .and. (adjusted_dt > min_adj_dt))
            !print *, 'passed on first try'
            call rk5(adjusted_dt, test_state, err)
            adjusted_dt = adjusted_dt * 0.84d0 * (tol / err)**(0.25d0)
        end do
        state = test_state
        t = t + adjusted_dt
    end
    subroutine rkf45_interpolation(t, delta_t, state, err, tol)
        real*8, intent(inout) :: t, err
        real*8, intent(in) :: delta_t, tol
        real*8 :: t_r, t_l, err_tmp
        type(triBody) :: state, state_r
        t_l = t
        t_r = t_l
        t = t + delta_t
        state_r = state
        call rkf45_iter(t_r, delta_t, state_r, err, tol)
        do while(t_r < t)
            t_l = t_r
            state = state_r
            call rkf45_iter(t_r, t - t_l, state_r, err_tmp, tol)
            err = err + err_tmp
        end do
        state = state + ((t - t_l) / (t_r - t_l)) * (state_r - state)
    end
end module solver_new