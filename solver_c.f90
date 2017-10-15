subroutine iter_step(h, t, delta_t, err, tol, p_num, p_m, p_x, p_v) bind(c, name='iter_step')
    use, intrinsic :: ISO_C_BINDING
    use tri_body_prob
    use solver
    integer(c_int), intent(in) :: p_num
    real(c_double), intent(inout) :: h, t, delta_t, err, tol, p_m(p_num), p_x(p_num, 3), p_v(p_num, 3)
    real*16 :: hh, tt, ddtt, eerrrr, ttooll
    type(planet) :: plist(p_num)
    type(triBody) :: state
    integer :: i,j
!    print *, h, t, delta_t, err, tol
!    print *, p_num
!    print *, p_m
!    print *, p_x
!    print *, p_v
    hh = h
    tt = t
    ddtt = delta_t
    eerrrr = err
    ttooll = tol
    do i=1, p_num
        plist(i)%m = p_m(i)
        plist(i)%x(:) = p_x(i,:)
        plist(i)%v(:) = p_v(i,:)
    end do
    state = triBody(plist)
    call rkf45(hh,tt,ddtt,state,eerrrr,ttooll)
    do i=1, p_num
        !p_m(i) = state%p(i)%m
        p_x(i,:) = state%p(i)%x(:)
        p_v(i,:) = state%p(i)%v(:)
    end do
    h = hh
    t = tt
    delta_t = ddtt
    err = eerrrr
    tol = ttooll
end subroutine iter_step