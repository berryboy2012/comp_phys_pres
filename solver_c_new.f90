module new_glue
    use tri_body_prob
    use solver_new
    implicit none
    type(triBody) :: gstate
    real*16 :: gt, gdt, gerr, gtol
    contains
    subroutine init(t, delta_t, tol, p_num, p_m, p_x, p_v) bind(c, name='init')
        use, intrinsic :: ISO_C_BINDING
        integer(c_int), intent(in) :: p_num
        real(c_double), intent(in) :: t, delta_t, tol, p_m(p_num), p_x(p_num, 3), p_v(p_num, 3)
        type(planet) :: plist(p_num)
        integer :: i
        !print *, t, delta_t, tol, p_num
        gt = t
        gdt = delta_t
        gtol = tol
        do i=1, p_num
            !print *, 'No.', i
            plist(i)%m = p_m(i)
            !print *, p_m(i)
            plist(i)%x(:) = p_x(i,:)
            !print *, p_x(i,:)
            plist(i)%v(:) = p_v(i,:)
            !print *, p_v(i,:)
        end do
        gstate = triBody(plist)
    end
    subroutine iter_step(t, err, energy, p_num, p_x, p_v) bind(c, name='iter_step')
        use, intrinsic :: ISO_C_BINDING
        integer(c_int), intent(in) :: p_num
        real(c_double), intent(out) :: t, err, energy, p_x(p_num, 3), p_v(p_num, 3)
        type(planet) :: plist(p_num)
        integer :: i
        gerr = 0.q0
        !print *, gt, gdt, gerr, gtol
        call rkf45(gt,gdt,gstate,gerr,gtol)

        do i=1, p_num
            !p_m(i) = state%p(i)%m
            p_x(i,:) = gstate%p(i)%x(:)
            p_v(i,:) = gstate%p(i)%v(:)
        end do
        t = gt
        err = gerr
        energy = gstate%energy()
    end subroutine iter_step
    subroutine iter_step_new(t, err, energy, p_num, p_x, p_v) bind(c, name='iter_step_new')
        use, intrinsic :: ISO_C_BINDING
        integer(c_int), intent(in) :: p_num
        real(c_double), intent(out) :: t, err, energy, p_x(p_num, 3), p_v(p_num, 3)
        type(planet) :: plist(p_num)
        integer :: i
        real*16 :: t_r, t_l, err_tmp
        type(triBody) :: state, state_r
        gerr = 0.q0
        !print *, gt, gdt, gerr, gtol
        t_l = gt
        t_r = t_l
        gt = gt + gdt
        state_r = gstate
        call rkf45_iter(t_r, gdt, state_r, gerr, gtol)
        do while(t_r < gt)
            t_l = t_r
            gstate = state_r
            call rkf45_iter(t_r, gt - t_l, state_r, err_tmp, gtol)
            gerr = gerr + err_tmp
        end do
        state = gstate + ((gt - t_l) / (t_r - t_l)) * (state_r - gstate)
        do i=1, p_num
            !p_m(i) = state%p(i)%m
            p_x(i,:) = state%p(i)%x(:)
            p_v(i,:) = state%p(i)%v(:)
        end do
        t = gt
        err = gerr
        energy = state%energy()
    end subroutine
end module new_glue