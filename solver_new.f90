module solver_new
    use tri_body_prob
    implicit none
    real*16, parameter :: min_adj_dt = 1.Q-9
    real*16, parameter, dimension(6) :: rkf_tki = &
            [0.q0, 1.q0 / 4.q0, 3.q0 / 8.q0, 12.q0 / 13.q0, 1.q0, 1.q0 / 2.q0]
    real*16, parameter, dimension(6, 6) :: rkf_yki = transpose( reshape( [&
            0.q0, 0.q0, 0.q0, 0.q0, 0.q0, 0.q0, &
            1.q0 / 4.q0, 0.q0, 0.q0, 0.q0, 0.q0, 0.q0, &
            3.q0 / 32.q0, 9.q0 / 32.q0, 0.q0, 0.q0, 0.q0, 0.q0, &
            1932.q0 / 2197.q0, -7200.q0 / 2197.q0, 7296.q0 / 2197.q0, 0.q0, 0.q0, 0.q0, &
            439.q0 / 216.q0, -8.q0, 3680.q0 / 513.q0, -845.q0 / 4104.q0, 0.q0, 0.q0, &
            -8.q0 / 27.q0, 2.q0, -3544.q0 / 2565.q0, 1859.q0 / 4104.q0, -11.q0 / -40.q0, 0.q0 &
    ], [6, 6]))
    real*16, parameter, dimension(6) :: rkf_yk = &
            [25.q0 / 216.q0, 0.q0, 1408.q0 / 2565.q0, 2197.q0 / 4101.q0, -1.q0 / 5.q0, 0.q0]
    real*16, parameter, dimension(6) :: rkf_zk = &
            [16.q0 / 135.q0, 0.q0, 6656.q0 / 12825.q0, 28561.q0 / 56430.q0, -9.q0 / 50.q0, 2.q0 / 55.q0]
    contains
    subroutine rk5(h, w, err)
        real*16, intent(in) :: h
        type(triBody) :: w, yl
        real*16 :: err
        type(triBody), dimension(6) :: k
        integer :: i, j

        do i=1, 6
            k(i) = w
            do j=1, i-1
                k(i) = k(i) + rkf_yki(i, j) * h * k(j)
            end do
            k(i) = fdt(rkf_tki(i) * h, k(i))
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
        real*16 :: adjusted_dt, err_tmp
        real*16, intent(inout) :: t, err
        real*16, intent(in) :: delta_t, tol
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
            adjusted_dt = delta_t * 0.84q0 * (tol / err_tmp)**(0.25q0)
            !print *, adjusted_dt / delta_t * 100.q0, '%'
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
    subroutine rkf45_loop(t, delta_t, state, err, tol)
        real*16 :: adjusted_dt, err_tmp, t_f
        real*16, intent(inout) :: t, err
        real*16, intent(in) :: delta_t, tol
        type(triBody) :: state, test_state
        test_state = state
        t_f = t + delta_t
        call rk5(delta_t, test_state, err_tmp)
        if ((isnan(err_tmp) .eqv. .true.) .or. (err_tmp < tol)) then
            !print *, 'passed on first try'
            state = test_state
            t = t + delta_t
            err = err + err_tmp
            return
        end if
        adjusted_dt = delta_t * 0.84q0 * (tol / err_tmp)**(0.25q0)
        do while (t + adjusted_dt < t_f .and. adjusted_dt > min_adj_dt)
            call rk5(adjusted_dt, test_state, err_tmp)
            t = t + adjusted_dt
            err = err + err_tmp
            adjusted_dt = adjusted_dt * 0.84q0 * (tol / err_tmp)**(0.25q0)
        end do
        call rk5(t_f - t, test_state, err_tmp)
        state = test_state
        err = err + err_tmp
        t = t_f
    end
    subroutine rkf45_iter(t, delta_t, state, err, tol)
        real*16 :: adjusted_dt
        real*16, intent(inout) :: t, err
        real*16, intent(in) :: delta_t, tol
        type(triBody) :: state, test_state
        test_state = state
        call rk5(delta_t, test_state, err)
        adjusted_dt = delta_t
        do while ((isnan(err) .eqv. .false.) .and. (err > tol) .and. (adjusted_dt > min_adj_dt))
            !print *, 'passed on first try'
            call rk5(adjusted_dt, test_state, err)
            adjusted_dt = adjusted_dt * 0.84q0 * (tol / err)**(0.25q0)
        end do
        state = test_state
        t = t + adjusted_dt
    end
    subroutine rkf45_interpolation(t, delta_t, state, err, tol)
        real*16, intent(inout) :: t, err
        real*16, intent(in) :: delta_t, tol
        real*16 :: t_r, t_l, err_tmp
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