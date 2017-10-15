module ploter

    use tri_body_prob
    implicit none
    contains
    
    subroutine write_state(data_unit,t,err,state)
        
        integer         ::  data_unit
        integer         ::  i
        real*16          ::  t, Etotal, err
        type(triBody)   ::  state
        
        Etotal = energy(state)
        write( data_unit, '(3f, $)' ) t, Etotal , log10(err)
        do i=1, state%p_num
            write( data_unit, '(3f, $)') state%p(i)%x
        end do
        write( data_unit, *) ' '
    end subroutine write_state
 
!    subroutine plot_state(plot_unit,state,k,frame_digs,plot_title,xrange,yrange)
!
!        integer         ::  plot_unit
!        integer         ::  k, frame_digs, zeros, i
!        type(triBody)   ::  state
!        character(255)  ::  plot_title
!        real*16          ::  xrange(2), yrange(2)
!        character(255)  ::  string, cxrange, cyrange, cxy1, cxy2, cxy3
!
!        if (k == 0) then
!            zeros = frame_digs - 1
!        else
!            zeros = frame_digs - int( log10( dble(k) ) ) - 1
!        end if
!        write(string, *) k
!        string = adjustl( string )
!        do i=1,zeros
!            string = trim('0')//trim(string)
!        end do
!
!        write(cxrange,*) xrange(1),':',xrange(2)
!        write(cyrange,*) yrange(1),':',yrange(2)
!
!        write(cxy1,*) state%a%x(1),' ', state%a%x(2)
!        write(cxy2,*) state%b%x(1),' ', state%b%x(2)
!        write(cxy3,*) state%c%x(1),' ', state%c%x(2)
!
!        write(plot_unit,*) '# Start a new plot: plot'//trim(string)
!        write(plot_unit,*) 'reset'
!        write(plot_unit,*)
!        write(plot_unit,*) "set terminal pngcairo size 2000,2000 enhanced font 'Verdana,25'"
!        write(plot_unit,*) 'set output "frames/'//trim(trim('plot')//string)//'.png"'
!        write(plot_unit,*)
!        write(plot_unit,*) '# Styling'
!        write(plot_unit,*) 'set border linewidth 1.5'
!        write(plot_unit,*) 'set pointsize 1.5'
!        write(plot_unit,*) "set style line 1 lc rgb '#0060ad' pt 5   # square"
!        write(plot_unit,*) "set style line 2 lc rgb '#0060ad' pt 7   # circle"
!        write(plot_unit,*) "set style line 3 lc rgb '#0060ad' pt 9   # triangle"
!        write(plot_unit,*)
!        write(plot_unit,*) 'unset key'
!        write(plot_unit,*)
!        write(plot_unit,*) 'set xrange['//trim(cxrange)//']'
!        write(plot_unit,*) 'set yrange['//trim(cyrange)//']'
!        write(plot_unit,*) "set xlabel 'x'"
!        write(plot_unit,*) "set ylabel 'y'"
!        write(plot_unit,*) 'set title "'//trim(plot_title)//'"'
!        write(plot_unit,*)
!        write(plot_unit,*) "plot '-' w p ls 1, '-' w p ls 2, '-' w p ls 3"
!        write(plot_unit,*) trim(cxy1)
!        write(plot_unit,*) 'e'
!        write(plot_unit,*) trim(cxy2)
!        write(plot_unit,*) 'e'
!        write(plot_unit,*) trim(cxy3)
!        write(plot_unit,*) 'e'
!        write(plot_unit,*)
!        write(plot_unit,*)
!
!    end subroutine plot_state
    
end module