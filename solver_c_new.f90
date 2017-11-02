module new_glue
    use tri_body_prob
    use solver_new
    implicit none
    type(triBody) :: gstate
    real*8 :: gt, gdt, gerr, gtol
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
        gerr = 0.d0
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
    subroutine iter_step_alter(t, err, energy, p_num, p_x, p_v) bind(c, name='iter_step_alter')
        use, intrinsic :: ISO_C_BINDING
        implicit none
        integer(c_int), intent(in) :: p_num
        real(c_double), intent(out) :: t, err, energy, p_x(p_num, 3), p_v(p_num, 3)
        type(planet) :: plist(p_num)
        integer :: i
        gerr = 0.d0
        !print *, 'ping'
        !call gstate%show_triBody()
        call rkf45loop(gt,gdt,gstate,gerr,gtol)
        !print *, gt, gdt, gerr, gtol
        !call gstate%show_triBody()
        do i=1, p_num
            !print *, 'o'
            !p_m(i) = state%p(i)%m
            p_x(i,:) = gstate%p(i)%x(:)
            p_v(i,:) = gstate%p(i)%v(:)
        end do
        !print *, 'pong'
        t = gt
        err = gerr
        energy = gstate%energy()
        !print *, 'pong'
        !energy = 0.d0
    end
end module new_glue