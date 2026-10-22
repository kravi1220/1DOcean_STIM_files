! Persistent Icepack column driven through stdin/stdout (Block D: MPAS-ADC x Icepack).
!
! Uses the same STIM wrapper as GOTM (stim_icepack: Icepack's own ice_step every
! sub-step, Icepack clock convention, frazil/frzmlt as in GOTM, icepack_in in cwd).
!
! stdin, first line : lat lon dt start_sec            (start_sec = s since Jan 1 of year_init)
! stdin, then per sub-step one line:
!   Tair_C Qa u10 v10 rhoa flwd fswd prec_mps snow_mps sst_C sss_psu h_top rho0 cp
! stdout, after every 'nsub' sub-steps (nsub from the first line's 5th value):
!   aice vice vsno Tsfc Sice  <fhocn_ai> <fswthru_ai> <fresh_ai> <fsalt_ai> <frazil_heat>
!   (state: end of the block; fluxes: means over the block, grid-box means, positive into ocean)
! A line starting with 'END' stops the program.
program icepack_pipe
   use stim_icepack
   implicit none
   integer, parameter :: rk = kind(1.d0)
   real(rk) :: lat, lon, dt, start_sec, f(14), acc(5)
   integer :: nsub, k, ios
   character(len=1024) :: line

   read(*,*) lat, lon, dt, start_sec, nsub
   call init_stim_icepack(lat, lon, dt, start_sec)
   write(*,'(a)') 'READY'
   flush(6)
   do
      acc = 0._rk
      do k = 1, nsub
         read(*,'(a)', iostat=ios) line
         if (ios /= 0) stop 'icepack_pipe: unexpected end of input'
         if (line(1:3) == 'END') then
            call clean_stim_icepack()
            stop
         end if
         read(line,*) f
         call do_stim_icepack(dt, f(1), f(2), f(3), f(4), f(5), f(6), f(7), f(8), f(9), &
                              f(10), f(11), f(12), f(13), f(14))
         acc = acc + (/ ip_fhocn, ip_fswthru, ip_fresh, ip_fsalt, ip_frazil_heat /)
      end do
      write(*,'(10es24.15)') ip_aice, ip_vice, ip_vsno, ip_Tsfc, ip_Sice, acc/nsub
      flush(6)
   end do
end program icepack_pipe
