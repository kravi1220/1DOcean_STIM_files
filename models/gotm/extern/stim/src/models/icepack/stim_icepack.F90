   MODULE stim_icepack
     !! Icepack (CICE-Consortium) column physics as a STIM ice model for GOTM.
     !!
     !! One Icepack column is advanced with Icepack's own driver step (icedrv
     !! ice_step: prep_radiation, step_therm1, step_therm2, update_state, step_snow,
     !! step_radiation, coupling_prep) once per GOTM time step, exactly as Icepack's
     !! standalone run loop does, except that the forcing arrays are filled from GOTM
     !! instead of from Icepack's forcing files and Icepack's slab ocean is not used
     !! (oceanmixed_ice must be .false.).
     !!
     !! Atmosphere -> ice part: ERA5 meteorology as used by GOTM (air T, dew-point
     !!   humidity, 10 m wind, pressure via air density, downward SW and LW, rain and
     !!   snow); Icepack computes its own fluxes over ice (icepack_atm_boundary).
     !! Ocean -> ice: top-layer in-situ T and practical S, freezing/melting potential
     !!   frzmlt = (Tf - sst) rho0 cp h_top / dt, limited to +-frzmlt_max (CICE
     !!   convention); no ice motion, so ice-ocean stress = 0 (ustar = ustar_min).
     !! Ice -> ocean (grid-box means, positive into the ocean): fhocn_ai, fswthru_ai,
     !!   fresh_ai, fsalt_ai, plus the frazil heat the ocean gives up when frzmlt > 0.
     !! The open-water fraction (1 - aice) is handled by GOTM itself (Kondo fluxes).
     !!
     !! Settings are read by Icepack from 'icepack_in' in the run directory.
     !! Author: block C of the KG6 experiment design (2026-10-05), AI-assisted (Claude).

   use stim_variables, only: rk, Hice
   use icedrv_kinds, only: dbl_kind, int_kind
   IMPLICIT NONE

   private
   public init_stim_icepack, do_stim_icepack, clean_stim_icepack

   real(rk), public :: ip_aice=0._rk, ip_vice=0._rk, ip_vsno=0._rk, ip_Tsfc=0._rk, ip_Sice=0._rk
   real(rk), public :: ip_fhocn=0._rk, ip_fswthru=0._rk, ip_fresh=0._rk, ip_fsalt=0._rk
   real(rk), public :: ip_frazil_heat=0._rk, ip_frzmlt=0._rk, ip_Tf=0._rk

   real(rk), parameter :: frzmlt_max = 1000._rk   ! W/m2, as CICE
   integer :: diag_unit = -1
   real(rk) :: acc(8) = 0._rk, acc_t = 0._rk
   integer  :: nt_Tsfc, nt_sice, nilyr

   contains

!-----------------------------------------------------------------------

   SUBROUTINE init_stim_icepack(lat_deg, lon_deg, dt_gotm, start_sec)
     !! lat/lon of the column [deg], GOTM time step [s], GOTM start time in seconds
     !! since 1 January of Icepack's year_init (used only to verify the calendar).
   use icedrv_InitMod, only: icedrv_initialize
   use icedrv_init, only: TLAT, TLON
   use icedrv_calendar, only: dt, time
   use icedrv_forcing, only: oceanmixed_ice
   use icepack_intfc, only: icepack_query_tracer_indices, icepack_query_parameters, &
                            icepack_query_tracer_sizes
   use icedrv_domain_size, only: nilyr_d => nilyr
   real(rk), intent(in) :: lat_deg, lon_deg, dt_gotm, start_sec
   real(dbl_kind) :: pi
!-----------------------------------------------------------------------
   write(*,*) 'init_stim_icepack()'
   call icedrv_initialize()
   call icepack_query_parameters(pi_out=pi)
   TLAT(:) = lat_deg*pi/180._dbl_kind
   TLON(:) = lon_deg*pi/180._dbl_kind
   if (abs(dt - dt_gotm) > 1.e-9_rk) then
      write(*,*) 'icepack_in dt =', dt, ' differs from GOTM dt =', dt_gotm
      stop 'init_stim_icepack(): dt mismatch'
   end if
   ! icedrv_initialize advances the clock by one step: during ice_step, time is the END
   ! of the step being computed (Icepack's own convention), i.e. start + dt here
   if (abs(time - dt - start_sec) > 0.5_rk) then
      write(*,*) 'icepack_in start (istep0*dt) =', time-dt, ' differs from GOTM start =', start_sec
      stop 'init_stim_icepack(): start time mismatch (set istep0 = start_sec/dt)'
   end if
   if (oceanmixed_ice) stop 'init_stim_icepack(): set oceanmixed_ice = .false. in icepack_in'
   call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_sice_out=nt_sice)
   nilyr = nilyr_d
   open(newunit=diag_unit, file='icepack_column.dat', status='replace', action='write')
   write(diag_unit,'(a)') '# hours_since_Jan1_year_init aice vice_m vsno_m hi_m Tsfc_C Sice_bulk_ppt' // &
        ' fhocn_ai_Wm2 fswthru_ai_Wm2 fresh_ai_kgm2s fsalt_ai_kgm2s frazil_heat_Wm2 frzmlt_Wm2 Tf_C' // &
        ' (fluxes: means over the preceding hour; state: instantaneous)'
   write(*,*) 'Icepack column at lat,lon =', lat_deg, lon_deg, ' dt =', dt, ' start time =', time-dt
   write(*,*) 'done.'
   END SUBROUTINE init_stim_icepack

!-----------------------------------------------------------------------

   SUBROUTINE do_stim_icepack(dtg, Tair_C, Qa_in, u10, v10, rhoa_in, flwd, fswd, &
                              prec_mps, snow_mps, sst_C, sss_psu, h_top, rho0, cp)
     !! Advance Icepack by one GOTM step. Inputs: air T [C], specific humidity
     !! [kg/kg], 10 m wind [m/s], air density [kg/m3], downward LW and SW [W/m2],
     !! total precipitation and snowfall [m/s water], top-layer in-situ T [C] and
     !! practical S [psu], top-layer thickness [m], rho0 [kg/m3], cp [J/kg/K].
   use icedrv_flux, only: Tair, Qa, uatm, vatm, wind, potT, rhoa, zlvl, flw, &
                          swvdr, swvdf, swidr, swidf, fsw, frain, fsnow, strax, stray, &
                          sst, sss, Tf, frzmlt, uocn, vocn, qdp, hmix, &
                          fhocn_ai, fswthru_ai, fresh_ai, fsalt_ai, init_flux_atm_ocn
   use icedrv_state, only: aice, vice, vsno, trcr
   use icedrv_RunMod, only: ice_step
   use icedrv_calendar, only: istep, istep1, time, dt, calendar
   use icepack_intfc, only: icepack_sea_freezing_temperature, icepack_query_parameters
   real(rk), intent(in) :: dtg, Tair_C, Qa_in, u10, v10, rhoa_in, flwd, fswd
   real(rk), intent(in) :: prec_mps, snow_mps, sst_C, sss_psu, h_top, rho0, cp
   real(dbl_kind) :: Tffresh, rhofresh
   real(rk) :: s, t_end
   integer :: k
!-----------------------------------------------------------------------
   call icepack_query_parameters(Tffresh_out=Tffresh, rhofresh_out=rhofresh)

   ! atmosphere, as icedrv prepare_forcing does it (frc* = Icepack band fractions)
   Tair(:) = Tair_C + Tffresh
   potT(:) = Tair(:)
   Qa(:)   = max(Qa_in, 0._rk)
   uatm(:) = u10;  vatm(:) = v10
   wind(:) = sqrt(u10*u10 + v10*v10)
   strax(:) = 0._rk;  stray(:) = 0._rk          ! calc_strair = .true.: Icepack computes stress
   rhoa(:) = max(rhoa_in, 0._rk)
   zlvl(:) = 10._rk
   flw(:)  = flwd
   fsw(:)  = max(fswd, 0._rk)
   swvdr(:) = fsw(:)*0.28_rk; swvdf(:) = fsw(:)*0.24_rk
   swidr(:) = fsw(:)*0.31_rk; swidf(:) = fsw(:)*0.17_rk
   fsnow(:) = max(snow_mps, 0._rk)*rhofresh
   frain(:) = max(prec_mps - snow_mps, 0._rk)*rhofresh

   ! ocean
   sss(:) = max(sss_psu, 0._rk)
   sst(:) = sst_C
   Tf(:)  = icepack_sea_freezing_temperature(sss(1))
   frzmlt(:) = max(min((Tf(1) - sst_C)*rho0*cp*h_top/dtg, frzmlt_max), -frzmlt_max)
   uocn(:) = 0._rk; vocn(:) = 0._rk; qdp(:) = 0._rk; hmix(:) = h_top

   call init_flux_atm_ocn()
   t_end = time                      ! end of the step being computed (Icepack convention)
   call ice_step()
   istep  = istep  + 1
   istep1 = istep1 + 1
   time   = time + dt
   call calendar(time)

   ip_aice = aice(1); ip_vice = vice(1); ip_vsno = vsno(1)
   ip_Tsfc = trcr(1,nt_Tsfc)
   s = 0._rk
   do k = 1, nilyr
      s = s + trcr(1,nt_sice+k-1)
   end do
   ip_Sice = s/nilyr
   ip_fhocn = fhocn_ai(1); ip_fswthru = fswthru_ai(1)
   ip_fresh = fresh_ai(1); ip_fsalt = fsalt_ai(1)
   ip_frzmlt = frzmlt(1); ip_Tf = Tf(1)
   ip_frazil_heat = max(frzmlt(1), 0._rk)   ! heat the ocean releases to reach Tf (W/m2)
   Hice = ip_vice                            ! grid-box mean thickness, for GOTM output

   acc = acc + dtg*(/ ip_fhocn, ip_fswthru, ip_fresh, ip_fsalt, ip_frazil_heat, ip_frzmlt, 0._rk, 0._rk /)
   acc_t = acc_t + dtg
   if (mod(t_end, 3600._rk) < 0.5_rk*dtg .or. 3600._rk - mod(t_end, 3600._rk) < 0.5_rk*dtg) then
      write(diag_unit,'(f12.4,3es14.6,es14.6,2f10.4,6es14.6,f9.4)') t_end/3600._rk, ip_aice, ip_vice, ip_vsno, &
           ip_vice/max(ip_aice,1.e-10_rk), ip_Tsfc, ip_Sice, acc(1:6)/acc_t, ip_Tf
      acc = 0._rk; acc_t = 0._rk
   end if
   END SUBROUTINE do_stim_icepack

!-----------------------------------------------------------------------

   SUBROUTINE clean_stim_icepack()
   if (diag_unit > 0) close(diag_unit)
   END SUBROUTINE clean_stim_icepack

!-----------------------------------------------------------------------

   END MODULE stim_icepack
