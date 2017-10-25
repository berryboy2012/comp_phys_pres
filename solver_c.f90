module glue
    use tri_body_prob
    use solver
    type(triBody) :: state
    real*16 :: gh, gt, gdt, gerr, gtol
contains
    subroutine init(h, t, delta_t, err, tol, p_num, p_m, p_x, p_v) bind(c, name='init')
        use, intrinsic :: ISO_C_BINDING
        integer(c_int), intent(in) :: p_num
        real(c_double), intent(in) :: h, t, delta_t, err, tol, p_m(p_num), p_x(p_num, 3), p_v(p_num, 3)
        type(planet) :: plist(p_num)
        integer :: i,j
        gh = h
        gt = t
        gdt = delta_t
        gerr = err
        gtol = tol
        do i=1, p_num
            !print *, i
            plist(i)%m = p_m(i)
            !print *, p_m(i)
            plist(i)%x(:) = p_x(i,:)
            !print *, p_x(i,:)
            plist(i)%v(:) = p_v(i,:)
            !print *, p_v(i,:)
        end do
        state = triBody(plist)
    end
    subroutine iter_step(h, t, delta_t, err, tol, p_num, p_m, p_x, p_v) bind(c, name='iter_step')
        use, intrinsic :: ISO_C_BINDING
        integer(c_int), intent(in) :: p_num
        real(c_double), intent(inout) :: h, t, delta_t, err, tol, p_m(p_num), p_x(p_num, 3), p_v(p_num, 3)
        type(planet) :: plist(p_num)
        integer :: i,j
        call rkf45(gh,gt,gdt,state,gerr,gtol)
        do i=1, p_num
            !p_m(i) = state%p(i)%m
            p_x(i,:) = state%p(i)%x(:)
            p_v(i,:) = state%p(i)%v(:)
        end do
        h = gh
        t = gt
        delta_t = gdt
        err = gerr
        tol = gtol
    end subroutine iter_step
end module