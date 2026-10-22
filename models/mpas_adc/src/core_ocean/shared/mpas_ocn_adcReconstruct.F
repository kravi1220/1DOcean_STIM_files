module ocn_adc_mixing_fused

  use mpas_kind_types
  use ocn_constants
  use ocn_turbulence
  use mpas_log
  use ocn_config
  use mpas_constants

  implicit none

  logical :: defineFirst, stopflag

  integer :: i1_f, i2_f, i3_f, i1, i2
  real(kind=RKIND) :: Cw1, Cw2, Cw3

  contains

  subroutine get_array_pointers

      i1 = mod(iterCount - 1, 2) + 1
      i2 = mod(iterCount, 2) + 1
      
      i3_f = mod(iterCount+2, 3) + 1 ! {3,1,2}
      i2_f = mod(iterCount+1, 3) + 1 ! {2,3,1}
      i1_f = mod(iterCount  , 3) + 1 ! {1,2,3}
      
  end subroutine get_array_pointers

  subroutine get_weights

      real (kind=RKIND) :: w3_ab3 = 1.91666666667_RKIND, &
        w2_ab3 = -1.3333333333_RKIND, &
        w1_ab3 = 0.41666666667_RKIND

      if (iterCount <= 1) then
         Cw3 = 1.0_RKIND
         Cw2 = 0.0_RKIND
         Cw1 = 0.0_RKIND
      elseif (iterCount == 2) then
         Cw3 = 1.5_RKIND
         Cw2 = -0.5_RKIND
         Cw1 = 0.0_RKIND
      else
         Cw3 = w3_ab3
         Cw2 = w2_ab3
         Cw1 = w1_ab3
      endif

  end subroutine get_weights

  subroutine dissipation_lengths2(nCells,nVertLevels,activeTracers,alphaT,betaS)
      integer,intent(in) :: nVertLevels, nCells
      real (kind=RKIND), dimension(2,nVertLevels,nCells), intent(in) :: activeTracers
      real (kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: alphaT, betaS
      integer :: iCell, k, ij, i

      real (kind=RKIND), dimension(nVertLevels) :: B, Bup, Bdn
      real (kind=RKIND), dimension(nVertLevels+1) :: tedge, sedge, KEsps, tke, BupEdge, BdnEdge, BEdge
      real (kind=RKIND), dimension(nVertLevels+1) :: alphaEdge, betaEdge
      real (kind=RKIND) :: sigav, tumdav, sumdav, Tup, Tdn, Sup, Sdn
      real (kind=RKIND) :: sumv, sumv0
      real (kind=RKIND), parameter :: refT = 15.0_RKIND, refS = 35.0_RKIND, minlen = 0.55_RKIND

!     NOTE: will need to convert to some form of displaced density in the mpas framework soon
!     possibly go back to a more traditional length scale formulation
      do iCell = 1,nCells
         do k = 2,nVertLevels
            tedge(k) = 0.5_RKIND*(activeTracers(1,k,iCell) + activeTracers(1,k-1,iCell))
            sedge(k) = 0.5_RKIND*(activeTracers(2,k,iCell) + activeTracers(2,k-1,iCell))
            alphaEdge(k) = 0.5_RKIND*(alphaT(k,iCell) + alphaT(k-1,iCell))
            betaEdge(k) = 0.5_RKIND*(betaS(k,iCell) + betaS(k-1,iCell))
         end do

         betaEdge(1) = betaEdge(2)
         alphaEdge(1) = alphaEdge(2)
         betaEdge(nVertLevels+1) = betaEdge(nVertLevels)
         alphaEdge(nVertLevels+1) = alphaEdge(nVertLevels)
         tedge(1) = tedge(2)
         sedge(1) = sedge(2)
         tedge(nVertLevels+1) = tedge(nVertLevels)
         sedge(nVertLevels+1) = sedge(nVertLevels)

         do k = 1,nVertLevels+1
            KEsps(k) = areaFraction(k,iCell) * KspsD(i1,k,iCell) + &
              (1.0_RKIND - areaFraction(k,iCell)) * KspsU(i1,k,iCell)
            tke(k) = 0.5_RKIND*(u2(i1,k,iCell) + v2(i1,k,iCell) + w2(i1,k,iCell))
         enddo

         do k = 1,nVertLevels+1
            Bedge(k) = gravity * (alphaEdge(k) * (tedge(k) - refT) - &
              betaEdge(k)  * (sedge(k) - refS))

            sigav = areaFraction(k,iCell)
            tumdav = tumd(k,iCell)
            sumdav = sumd(k,iCell)

            Tup = tedge(k) + sigav*tumdav
            Tdn = tedge(k) - (1.0_RKIND - sigav)*tumdav
            Sup = sedge(k) + sigav*sumdav
            Sdn = sedge(k) - (1.0_RKIND - sigav)*sumdav
            
            Bupedge(k) = gravity * (alphaEdge(k) * (Tup - refT) - betaEdge(k) * (Sup - refS))
            Bdnedge(k) = gravity * (alphaEdge(k) * (Tdn - refT) - betaEdge(k) * (Sdn - refS))
         enddo

         do k = 2,nVertLevels
!           updraft length scale
            sumv  = 0.0_RKIND
            sumv0 = 0.0_RKIND
            ij=k
            lenup(k,iCell) = 0.0_RKIND
            do while(sumv0 <= tke(k) .and. ij > 1)
               sumv = sumv0 + (BEdge(ij) - BupEdge(k)) * (ze(ij-1,iCell)-ze(ij,iCell))
               if (sumv > tke(k)) THEN
                  lenup(k,iCell) = max(minlen, lenup(k,iCell) + &
                    abs(ze(ij-1,iCell)-ze(ij,iCell))*(tke(k)-sumv0)/(sumv-sumv0))
                    exit
               else
                  lenup(k,iCell) = lenup(k,iCell) + abs(ze(ij-1,iCell)-ze(ij,iCell))
               endif
               sumv0 = sumv
               ij = ij - 1
            end do

!           downdraft length scale
            sumv  = 0.0_RKIND
            sumv0 = 0.0_RKIND
            ij=k
            lendn(k,iCell) = 0.0_RKIND
            do while(sumv0 <= tke(k) .and. ij < nVertLevels+1)
               sumv = sumv0 + (BdnEdge(k) - BEdge(ij)) * (ze(ij-1,iCell)-ze(ij,iCell))
               if (sumv > tke(k)) THEN
                  lendn(k,iCell) = max(minlen, lendn(k,iCell) + &
                    abs(ze(ij-1,iCell)-ze(ij,iCell))*(tke(k)-sumv0)/(sumv-sumv0))
                    exit
               else
                  lendn(k,iCell) = lendn(k,iCell) + abs(ze(ij-1,iCell)-ze(ij,iCell))
               endif
               sumv0 = sumv
               ij = ij + 1
            end do

            length(k,iCell) = 2.0_RKIND * lenup(k,iCell) * lendn(k,iCell)  &
              / (lenup(k,iCell) + lendn(k,iCell))
         enddo
      enddo

      length(1,iCell) = 0.55_RKIND
      length(nVertLevels+1,iCell) = 0.55_RKIND

  end subroutine dissipation_lengths2

  subroutine compute_ADC_tends(nCells,nVertLevels, nTracers, dt,activeTracers, uvel, vvel, BVF,  &
                  uwsfc, vwsfc, wtsfc, wssfc, alphaT, betaS, fCell, boundaryLayerDepth)

      integer,intent(in) :: nCells, nVertLevels, nTracers
      real(kind=RKIND),intent(in) :: dt

      real(kind=RKIND),dimension(nTracers,nVertLevels,nCells),intent(inout) :: activeTracers
      real(kind=RKIND),dimension(nVertLevels,nCells),intent(inout) :: uvel, vvel, alphaT, betaS
      real(kind=RKIND),dimension(nCells),intent(in) :: uwsfc, vwsfc, wtsfc, wssfc, fCell
      real(kind=RKIND),dimension(nCells),intent(in) :: boundaryLayerDepth
      real(kind=RKIND),dimension(nVertLevels,nCells),intent(inout) :: BVF
      integer :: niter, iIter,iCell,k, supercycle

      real(kind=RKIND) :: dt_small

      real(kind=RKIND) :: Sw, St, Ss, Eav, Dav, sigav, sigavp1, wumdAv, tumdAv, sumdAv, wumdAvp1, tumdAvp1, sumdAvp1
      real(kind=RKIND) :: Swup, KspsUav, KspsDav, KspsUavp1, KspsDavp1, KE, Mcav, lenav,u2av,v2av,w2av
      real(kind=RKIND) :: w2tTemp, w2tCheck, w2sTemp, w2sCheck, w3temp, w3check2, w3check, mval, KEsps
      real(kind=RKIND) :: Uz, Vz, Ustz, Vstz, dz

      real(kind=RKIND) :: invLen, l, len1, len2, lenmax, integrandTop, integrandBot
      real(kind=RKIND) :: len0, len2_1, len2_2, sfcBuoy, lengthT, bvfT
      real(kind=RKIND) :: KEm1, KEp1, tauUP, tauDN, tomUP, tomDN
      real(kind=RKIND) :: tauM1, tau, tauP1, tauAV, utemp, vtemp
      real(kind=RKIND) :: B, Cval, diff, wtav, dzmid, Ksps, Sz, Tz, w4k, w4kp1, w2k, w2kp1
      real(kind=RKIND) :: lareaFraction, wstar, Q, w3av, tempMoment, frictionVelocity
      real(kind=RKIND) :: sfcFrictionVelocitySquared, wtSumUp, wtSumDn, wsSumUp, wsSumDn
      real(kind=RKIND) :: uStokes_surf, vStokes_surf

      real(kind=RKIND),dimension(nVertLevels,nCells) :: Swumd
      real(kind=RKIND),dimension(nVertLevels,nCells) :: tauw3, tau_tracer, tauVel, tauvVel
      real(kind=RKIND),dimension(nVertLevels,nCells) :: uStokes, vStokes, uStokesm, vStokesm
      real(kind=RKIND),dimension(nVertLevels,nCells) :: areaFractionMid, tumdMid, McMid, wumdMid, sumdMid

      real(kind=RKIND) :: SwU, SwD, tau_sfc, d_sqrt_wp2_dz, tp2, sp2, min_wp2_sfc_val, wp2_splat_sfc_correction
      real(kind=RKIND) :: min_wps_sfc_val

      min_wps_sfc_val = 1.0E-10_RKIND
      dt_small = config_adc_timestep
      niter = dt / dt_small

      call get_array_pointers
      call get_weights

      stopflag = .false.

      do iIter=1,niter
!     on further examination build_diagnostics array can live outside the iter loop
         do iCell=1,nCells
            Q = grav*(alphaT(1,iCell)*wtsfc(iCell) + betaS(1,iCell)*wssfc(iCell))* &
              boundaryLayerDepth(iCell)
              if(Q > 0) then
                 wstar = abs(Q)**(1.0_RKIND/3.0_RKIND)
              else
                 wstar = 0.0_RKIND
              endif

              tumd(1,iCell) = 0.0_RKIND
              wumd(1,iCell) = 0.0_RKIND
              areaFraction(1,iCell) = 0.5_RKIND
              Mc(1,iCell) = 0.0_RKIND
              w2t(1,iCell) = -0.3_RKIND*wstar * wtsfc(iCell)
              w2s(1,iCell) = 0.3_RKIND*wstar * wssfc(iCell)

              sfcFrictionVelocitySquared = uwsfc(iCell) + vwsfc(iCell)
              frictionVelocity = sqrt(sfcFrictionVelocitySquared + config_adc_bc_wstar * wstar * wstar)
              frictionVelocity = max( config_adc_frictionVelocityMin, frictionVelocity )

              ! Initialize Stokes drift surface values as zero, overwrite these if config_adc_langmuir = .True.
              uStokes_surf = 0.0_RKIND
              vStokes_surf = 0.0_RKIND
              if (config_adc_langmuir) then
                  uStokes_surf = COS(stokes_angle)*frictionVelocity / (La_t**2.0_RKIND)
                  vStokes_surf = SIN(stokes_angle)*frictionVelocity / (La_t**2.0_RKIND)
              endif

!             Calculate exponential Stokes drift profiles from surface values & decay depth
!             [stokes_decay_depth]. This is temporary; eventually MPAS will pass profiles
!             For now this is calculated using the turbulent Langmuir number La_t
!             La_t = (friction velocity / surface stokes drift magnitude)^(0.5) and the
!             stokes_angle (angle between surface stress and Stokes drift vectors)

              do k=1,nVertLevels
!                Calculate Stokes drift at vertical edge of cells
                 uStokes(k,iCell) = uStokes_surf*EXP(ze(k,iCell)/stokes_depth)
                 vStokes(k,iCell) = vStokes_surf*EXP(ze(k,iCell)/stokes_depth)
!                Calculate Stokes drift at vertical midpoint of cells
                 uStokesm(k,iCell) = uStokes_surf*EXP(zm(k,iCell)/stokes_depth)
                 vStokesm(k,iCell) = vStokes_surf*EXP(zm(k,iCell)/stokes_depth)
              enddo

              do k=1,3
                 u2(k,1,iCell) = config_adc_up2_vp2_factor*config_adc_bc_const*frictionVelocity**2.0_RKIND
                 v2(k,1,iCell) = config_adc_up2_vp2_factor*config_adc_bc_const*frictionVelocity**2.0_RKIND
                 w2(k,1,iCell) = config_adc_bc_const_wp2*frictionVelocity**2.0_RKIND
                 uw(k,1,iCell) = -uwsfc(iCell)
                 vw(k,1,iCell) = -vwsfc(iCell)
                 wt(k,1,iCell) = wtsfc(iCell)
                 ws(k,1,iCell) = -wssfc(iCell)
                 KE = 0.5_RKIND*(u2(i1,1,iCell) + v2(i1,1,iCell))
                 eps(k,1,iCell) = KE**1.5_RKIND/(0.5_RKIND*(ze(1,iCell) - ze(2,iCell) + 1.0E-10_RKIND))
              enddo
           enddo

           if(config_adc_use_splat_parameterization) then

!   COMPUTE SPLATTING EFFECT ON W2, U2, V2 AND W3
!   This splatting effect was adopted from the CLUBB code.
!   It computes the (negative) tendencies of w2 and w3 due to "splatting" of eddies as
!   their vertical motion is inhibited, e.g., near the surface or mixed layer base.
!   Near these 'splatting layers', wtend6 removes w2 and redistributes it evenly
!   between u2 and v2 [via u2tend6 and v2tend6].
!   The functional forms of splatting terms (w2tend6 and w2tend6) are
!
!   w2tend6 \propto - w'2 * (turbulent time scale) * ( d/dz( sqrt(w'2) ) )^2
!
!   The functional form of wp3_splat is
!
!   wp3_splat \propto - w'3 * (turbulent time scale) * ( d/dz( sqrt(w'2) ) )^2
!
!   If the coefficient on wp3_splat is at least 1.5 times greater than the
!   coefficient on wp2_splat, then skewness will be damped, promoting
!   more stratiform layers.
!
!   The calculation for k=1 is separate for performance reasons (one-sided derivatives)
!   Note that the splat_factor forces the term to be smaller than factor * dt

              k=1
              do iCell=1,nCells
                 d_sqrt_wp2_dz = (sqrt(w2(i1,1,iCell)) - sqrt(0.5_RKIND*(w2(i1,1,iCell)+ &
                   w2(i1,2,iCell)))) / (ze(1,iCell) - zm(1,iCell))
                 tau_sfc = length(1,iCell) / sqrt(0.5_RKIND*(u2(i1,k,iCell) + v2(i1,k,iCell)))
                 w2tend6(1,iCell) = min(max(-config_adc_splat_tend_max, -w2(i1,1,iCell)*      &
                   config_adc_splat_wp2_val*tau_sfc*d_sqrt_wp2_dz**2),      &
                   config_adc_splat_tend_max)
                 tau_sfc = 0.5_RKIND*(length(1,iCell) + length(2,iCell)) / sqrt(0.5_RKIND* &
                   (0.5_RKIND*(u2(i1,1,iCell) + u2(i1,2,iCell)) + 0.5_RKIND*(          &
                   v2(i1,1,iCell) + v2(i1,2,iCell)) + 0.5_RKIND*(w2(i1,1,iCell) +      &
                   w2(i1,2,iCell))))
                 d_sqrt_wp2_dz = (sqrt(w2(i1,1,iCell)) - sqrt(w2(i1,2,iCell))) / &
                   (ze(1,iCell) - ze(2,iCell))
                 w3tend6(1,iCell) = min(max(-config_adc_splat_tend_max, -w3(i1,1,iCell)*      &
                   config_adc_splat_wp2_val*tau_sfc*d_sqrt_wp2_dz**2),   &
                   config_adc_splat_tend_max)
              end do

              do iCell=1,nCells
                 do k=2,nVertLevels
                    d_sqrt_wp2_dz = (sqrt(w2(i1,k-1,iCell)) - sqrt(w2(i1,k+1,iCell))) / &
                      (ze(k-1,iCell) - ze(k+1,iCell))
                    tau_sfc = length(k,iCell) / sqrt(0.5_RKIND*(u2(i1,k,iCell) + v2(i1,k,iCell) + &
                      w2(i1,k,iCell)))
                    w2tend6(k,iCell) = min(max(-config_adc_splat_tend_max, -w2(i1,k,iCell)*      &
                      config_adc_splat_wp2_val*tau_sfc*d_sqrt_wp2_dz**2),      &
                      config_adc_splat_tend_max)
                    tau_sfc = 0.5_RKIND*(length(k,iCell) + length(k+1,iCell)) / sqrt(0.5_RKIND* &
                      (0.5_RKIND*(u2(i1,k,iCell) + u2(i1,k+1,iCell)) + 0.5_RKIND*(          &
                      v2(i1,k,iCell) + v2(i1,k+1,iCell)) + 0.5_RKIND*(w2(i1,k,iCell) +      &
                      w2(i1,k+1,iCell))))
                    d_sqrt_wp2_dz = (sqrt(w2(i1,k,iCell)) - sqrt(w2(i1,k+1,iCell))) / &
                      (ze(k,iCell) - ze(k+1,iCell))
                    w3tend6(1,iCell) = min(max(-config_adc_splat_tend_max, -w3(i1,k,iCell)*         &
                      config_adc_splat_wp2_val*tau_sfc*d_sqrt_wp2_dz**2), &
                      config_adc_splat_tend_max)
                 end do
              end do

              do iCell=1,nCells
                 tp2 = 0.4_RKIND * config_adc_bc_const * (wt(i1,1,iCell) / frictionVelocity)**2
                 sp2 = 0.4_RKIND * config_adc_bc_const * (ws(i1,1,iCell) / frictionVelocity)**2
                 min_wp2_sfc_val = max(1.0E-10_RKIND, wt(i1,1,iCell)**2 / (tp2 * 0.99_RKIND**2 + 1.0E-15_RKIND), &
                   ws(i1,1,iCell)**2 / (sp2 * 0.99_RKIND**2 + 1.0E-15_RKIND))
                 tau_sfc = length(1,iCell) / sqrt(0.5_RKIND*(u2(i1,1,iCell) + v2(i1,1,iCell)))

                 if(w2(i1,1,iCell) + tau_sfc * w2tend6(1,iCell) < min_wps_sfc_val) then
                    wp2_splat_sfc_correction = -w2(i1,1,iCell) + min_wp2_sfc_val
                    w2(i1,1,iCell) = min_wp2_sfc_val
                 else
                    wp2_splat_sfc_correction = tau_sfc * w2tend6(1,iCell)
                    w2(i1,1,iCell) = w2(i1,1,iCell) + wp2_splat_sfc_correction
                 end if
                 u2(i1,1,iCell) = u2(i1,1,iCell) - 0.5_RKIND * wp2_splat_sfc_correction
                 v2(i1,1,iCell) = v2(i1,1,iCell) - 0.5_RKIND * wp2_splat_sfc_correction
              end do

!     build up splat tend for u2 and v2
              do iCell=1,nCells
                 do k=2,nVertLevels
                    u2tend6(k,iCell) = -0.5_RKIND*w2tend6(k,iCell)
                    v2tend6(k,iCell) = -0.5_RKIND*w2tend6(k,iCell)
                 end do
              end do
           else
              do iCell=1,nCells
                 do k=1,nVertLevels
                    w3tend6(k,iCell) = 0.0_RKIND
                    w2tend6(k,iCell) = 0.0_RKIND
                    u2tend6(k,iCell) = 0.0_RKIND
                    v2tend6(k,iCell) = 0.0_RKIND
                 end do
              end do
           end if               ! end use splat correction


!     Kernel 1 inlined versions of the base arrays, needed for later to make them collapsible
           do iCell=1,nCells
              do k=2,nVertLevels
                 w3av = 0.5_RKIND*(w3(i1,k-1,iCell) + w3(i1,k,iCell))

                 Sw = w3av / (w2(i1,k,iCell)**1.5_RKIND + 1.0E-15_RKIND)
                 lareaFraction = 0.5_RKIND + 0.5_RKIND*Sw / sqrt(4.0_RKIND + Sw**2)

                 if(lareaFraction < 0.01_RKIND) then
                    lareaFraction = 0.01_RKIND
                 end if

                 if(lareaFraction > 0.99_RKIND) then
                    lareaFraction = 0.99_RKIND
                 end if

                 areaFraction(k,iCell) = lareaFraction

                 wumd(k,iCell) = sqrt(w2(i1,k,iCell) / (areaFraction(k,iCell) * &
                   (1.0_RKIND - areaFraction(k,iCell))))
                 Mc(k,iCell) = areaFraction(k,iCell)*(1.0_RKIND - &
                   areaFraction(k,iCell)) * wumd(k,iCell)

                 tumd(k,iCell) = wt(i1,k,iCell) / (1.0E-12_RKIND + Mc(k,iCell))
                 sumd(k,iCell) = ws(i1,k,iCell) / (1.0E-12_RKIND + Mc(k,iCell))

                 if(w2(i1,k,iCell) <= epsilon + 1.0e-9_RKIND) then
                    tumd(k,iCell) = 0.0_RKIND
                    sumd(k,iCell) = 0.0_RKIND
                 endif
!     redo mid quantities
                 areaFractionMid(k-1,iCell) = 0.5_RKIND*(areaFraction(k-1,iCell) + areaFraction(k,iCell))
                 wumdMid(k-1,iCell) = 0.5_RKIND*(wumd(k-1,iCell) + wumd(k,iCell))
                 tumdMid(k-1,iCell) = 0.5_RKIND*(tumd(k-1,iCell) + tumd(k,iCell))
                 sumdMid(k-1,iCell) = 0.5_RKIND*(sumd(k-1,iCell) + sumd(k,iCell))
                 McMid(k-1,iCell) = 0.5_RKIND*(Mc(k-1,iCell) + Mc(k,iCell))
              enddo
              areaFractionMid(nVertLevels,iCell) = 0.5_RKIND
              wumdMid(nVertLevels,iCell) = 0.0_RKIND
              tumdMid(nVertLevels,iCell) = 0.0_RKIND
              sumdMid(nVertLevels,iCell) = 0.0_RKIND
              McMid(nVertLevels,iCell) = 0.0_RKIND
           enddo

           do iCell = 1,nCells
              do k=2,nVertLevels
                 sigav = 0.5_RKIND*(areaFraction(k,iCell) + areaFraction(k+1,iCell))
                 tumdav = 0.5_RKIND*(tumd(k,iCell) + tumd(k+1,iCell))
                 sumdav = 0.5_RKIND*(sumd(k,iCell) + sumd(k+1,iCell))
                 wumdav = 0.5_RKIND*(wumd(k,iCell) + wumd(k+1,iCell))
                 sigav = areaFractionMid(k,iCell)
                 tumdav = tumdMid(k,iCell)
                 wumdav = wumdMid(k,iCell)
                 sumdav = sumdMid(k,iCell)
                 w2tTemp = -sigav*(1.0_RKIND - sigav)*(1.0_RKIND - 2.0_RKIND*sigav)*wumdav**2.0_RKIND*tumdav
                 w2tCheck = sigav*(1.0_RKIND - sigav)*sqrt(1.0_RKIND - 3.0_RKIND*sigav + 3.0_RKIND*sigav**2)* &
                   wumdav**2.0_RKIND*tumdav
                 w2t(k,iCell) = w2tTemp

                 w2sTemp = -sigav*(1.0_RKIND - sigav)*(1.0_RKIND - 2.0_RKIND*sigav)*wumdav**2.0_RKIND*sumdav
                 w2sCheck = sigav*(1.0_RKIND - sigav)*sqrt(1.0_RKIND - 3.0_RKIND*sigav + 3.0_RKIND*sigav**2)* &
                   wumdav**2.0_RKIND*sumdav
                 w2s(k,iCell) = w2sTemp

!     also use this loop to reset the cliptends for the step
                 u2cliptend(k,iCell) = 0.0_RKIND
                 v2cliptend(k,iCell) = 0.0_RKIND
                 w2cliptend(k,iCell) = 0.0_RKIND
              enddo
           enddo
           if(config_adc_use_old_length_scale)then
              call dissipation_lengths2(nCells,nVertLevels,activeTracers,alphaT,betaS)
           endif
           do iCell = 1,nCells
              do k=1,nVertLevels
!     change length scale to lenup and lendown
                 Entrainment(k,iCell) = Cww_E*areaFraction(k,iCell)*(1.0_RKIND- &
                   areaFraction(k,iCell))*Mc(k,iCell) / ( lenup(k,iCell) + 1.0E-15_RKIND )
                 Detrainment(k,iCell) = Cww_D*areaFraction(k,iCell)*(1.0_RKIND- &
                   areaFraction(k,iCell))*Mc(k,iCell) / ( lendn(k,iCell) +1.0E-15_RKIND )
              enddo
           enddo

!     This cell loop computes w3Tend and TOMs for later tendencies.  Should be okay to collapse loops here
           do iCell = 1,nCells
              do k=1,nVertLevels
                 Eav = 0.5_RKIND*(Entrainment(k+1,iCell) + Entrainment(k,iCell))
                 Dav = 0.5_RKIND*(Detrainment(k+1,iCell) + Detrainment(k,iCell))
                 u2av = 0.5_RKIND*(u2(i1,k,iCell) + u2(i1,k+1,iCell))
                 v2av = 0.5_RKIND*(v2(i1,k,iCell) + v2(i1,k+1,iCell))
                 w2av = 0.5_RKIND*(w2(i1,k,iCell) + w2(i1,k+1,iCell))

                 sigav = 0.5_RKIND*(areaFraction(k,iCell) + areaFraction(k+1,iCell))
                 wumdav = 0.5_RKIND*(wumd(k,iCell) + wumd(k+1,iCell))
                 tumdav = 0.5_RKIND*(tumd(k,iCell) + tumd(k+1,iCell))
                 sumdav = 0.5_RKIND*(sumd(k,iCell) + sumd(k+1,iCell))
                 Mcav = 0.5_RKIND*(Mc(k,iCell) + Mc(k+1,iCell))
                 lenav = 0.5_RKIND*(length(k,iCell) + length(k+1,iCell))
                 KE = sqrt(0.5_RKIND*(u2av+v2av+w2av))

                 dz = ze(k,iCell) - ze(k+1,iCell)
                 Swumd(k,iCell) = - 2.0_RKIND/3.0_RKIND*(1.0_RKIND/sigav*(areaFraction(k,iCell)*KspsD(i1,k,iCell) - &
                   areaFraction(k+1,iCell)*KspsD(i1,k+1,iCell)) / dz - 1.0_RKIND / (1.0_RKIND - &
                   sigav)*((1.0_RKIND - areaFraction(k,iCell))*KspsU(i1,k,iCell) - (1.0_RKIND - &
                   areaFraction(k+1,iCell))*KspsU(i1,k+1,iCell)) / dz)

                 tauw3(k,iCell) = C_mom_w3*KE / (1.0E-15_RKIND + sqrt(2.0_RKIND)*lenAv)

                 wumdav = wumdMid(k,ICell)
                 sigav = areaFractionMid(k,iCell)
                 Mcav = McMid(k,iCell)

                SwU = -(1.0_RKIND-2.0_RKIND*areaFraction(k,iCell)) &
                   /(areaFraction(k,iCell)*(1.0_RKIND-areaFraction(k,iCell)))**0.5_RKIND
                SwD = -(1.0_RKIND-2.0_RKIND*areaFraction(k+1,iCell)) &
                   /(areaFraction(k+1,iCell)*(1.0_RKIND-areaFraction(k+1,iCell)))**0.5_RKIND

!    Calculate Stokes drift derivatives at vertical cell centers (where w3 is located)
                Ustz = (uStokes(k,iCell) - uStokes(k+1,iCell)) / dz
                Vstz = (vStokes(k,iCell) - vStokes(k+1,iCell)) / dz

!      CALCULATE W3 TENDENCY TERMS
!      Calculate the tendency (rate-of-change) terms for the w3 budget.
!      The w3 tendency is the sum of the following terms (saved individually)
!                   (1) entrainment and detrainment between plumes
!                   (2) turbulent transport of w3 [ -d(w^4)/dz ] where w^4 is estimated using ADC
!                   (3) gradient production of w3 [ (3/2)*d( [w2]^2 )/dz ]
!                   (4) parameterized pressure terms (Rotta, buoyancy)
!                   (5) buoyancy production and sub-plume scale effects
!                   (6) splatting (calculated earlier)
!                   (7) Stokes drift production terms -3*[uw2*d(ustokes)/dz + vw2*d(ustokes)/dz]
!      The tendency term is truncated if config_adc_truncate_tend flag it switched on
                 w3tend1(k,iCell) = -wumdav**3.0_RKIND*(Eav*(3.0_RKIND*sigav - 2.0_RKIND) &
                   + Dav*(3.0_RKIND*sigav - 1.0_RKIND))

                 w3tend2(k,iCell) = -((3.0_RKIND + SwU**2.0_RKIND)*(w2(i1,k,iCell)**2.0_RKIND) - &
                   (3.0_RKIND + SwD**2.0_RKIND)*(w2(i1,k+1,iCell)**2.0_RKIND) ) / dz

                 w3tend3(k,iCell) = 1.5_RKIND*(w2(i1,k,iCell)**2.0_RKIND - w2(i1,k+1,iCell)**2.0_RKIND) / dz

                 w3tend4(k,iCell) = - tauw3(k,iCell)*w3(i1,k,iCell) - &
                   3.0_RKIND*c11*grav*(alphaT(k,iCell)*w2t(k,iCell) - betaS(k,iCell)*w2s(k,iCell))

                 w3tend5(k,iCell) = 3.0_RKIND*grav*(alphaT(k,iCell)*w2t(k,iCell) - betaS(k,iCell)* &
                   w2s(k,iCell)) - 3.0_RKIND*(1.0_RKIND - 2.0_RKIND*sigav)*Mcav*wumdav*Swumd(k,iCell)

                 w3tend7(k,iCell) = - 3.0_RKIND*(uw2(k,iCell)*Ustz + vw2(k,iCell)*Vstz)

                 w3tend(i3_f,k,iCell) = w3tend1(k,iCell) + w3tend2(k,iCell) + w3tend3(k,iCell) + &
                   w3tend4(k,iCell) + w3tend5(k,iCell) + w3tend6(k,iCell) + w3tend7(k,iCell)

                 if(k>1 .and. k < nVertLevels .and. kappa_w3 > 0.0_RKIND) then
                    w3tend(i3_f,k,iCell) = w3tend(i3_f,k,iCell) + kappa_w3*(w3(i1,k-1,iCell) &
                      - w3(i1,k+1,iCell)) / (zm(k-1,iCell) - zm(k+1,iCell))**2.0_RKIND
                 endif

                 if (config_adc_truncate_tend) then
                    w3tend(i3_f,k,iCell) = FLOAT (INT(w3tend(i3_f,k,iCell) * adcRound &
                      + 0.5_RKIND)) / adcRound
                 endif

!     now get all the downgradient third-order moments
                 KE = 0.5_RKIND*(u2(i1,k,iCell) + v2(i1,k,iCell) + w2(i1,k,iCell))
                 KEp1 = 0.5_RKIND*(u2(i1,k+1,iCell) + v2(i1,k+1,iCell) + w2(i1,k+1,iCell))
                 lenav = 0.5_RKIND*(length(k,iCell) + length(k+1,iCell))
                 diff = C_mom * sqrt(0.5_RKIND*(KE + KEp1)) * lenav
                 dz = ze(k,iCell) - ze(k+1,iCell)
                 uw2(k,iCell) = -diff*(uw(i1,k,iCell) - uw(i1,k+1,iCell)) / dz
                 vw2(k,iCell) = -diff*(vw(i1,k,iCell) - vw(i1,k+1,iCell)) / dz
                 u2w(k,iCell) = -diff*(u2(i1,k,iCell) - u2(i1,k+1,iCell)) / dz
                 v2w(k,iCell) = -diff*(v2(i1,k,iCell) - v2(i1,k+1,iCell)) / dz
                 uvw(k,iCell) = -diff*(uv(i1,k,iCell) - uv(i1,k+1,iCell)) / dz

                 diff = C_therm*sqrt(0.5_RKIND*(KE + KEp1)) * lenav
                 uwt(k,iCell) = -diff*(ut(i1,k,iCell) - ut(i1,k+1,iCell)) / dz
                 vwt(k,iCell) = -diff*(vt(i1,k,iCell) - vt(i1,k+1,iCell)) / dz
                 uws(k,iCell) = -diff*(us(i1,k,iCell) - us(i1,k+1,iCell)) / dz
                 vws(k,iCell) = -diff*(vs(i1,k,iCell) - vs(i1,k+1,iCell)) / dz

              enddo             !nVertLevels
           enddo                !nCells

!     next all second order moment tendencies also apply tendencies
           do iCell=1,nCells
              do k=2,nVertLevels
                 dzmid = (zm(k-1,iCell) - zm(k,iCell))
                 dz = ze(k-1,iCell) - ze(k,iCell)
                 B = grav*(alphaT(k,iCell)*areaFraction(k,iCell)*(1.0_RKIND - &
                   areaFraction(k,iCell))*wumd(k,iCell)*tumd(k,iCell) -      &
                   betaS(k,iCell)*areaFraction(k,iCell)*(1.0_RKIND -          &
                   areaFraction(k,iCell))*wumd(k,iCell)*sumd(k,iCell))

                 Uz = (uvel(k-1,iCell) - uvel(k,iCell)) / dzmid
                 Vz = (vvel(k-1,iCell) - vvel(k,iCell)) / dzmid
!    Calculate Stokes drift derivatives at vertical cell edges (where w2 etc. are located)
                 Ustz = (uStokesm(k-1,iCell) - uStokesm(k,iCell)) / dzmid
                 Vstz = (vStokesm(k-1,iCell) - vStokesm(k,iCell)) / dzmid
                 Tz = (activeTracers(1,k-1,iCell) - activeTracers(1,k,iCell)) / dzmid
                 Sz = (activeTracers(2,k-1,iCell) - activeTracers(2,k,iCell)) / dzmid
                 
                 KE = sqrt(0.5_RKIND*(u2(i1,k,iCell) + v2(i1,k,iCell) + areaFraction(k,iCell)* &
                   (1.0_RKIND - areaFraction(k,iCell))*wumd(k,iCell)**2.0_RKIND))

!      Define turbulence timescales/rates for return-to-isotropy tendency terms (Rotta, 1951)
!      Note that tau_tracer, tauVel and tauvVel are INVERSE timescales (i.e. rates)
                 tau_tracer(k,iCell) = c_slow_tracer*KE / (1.0E-15_RKIND + length(k,iCell))
                 tauVel(k,iCell) = C_slow*KE / (1.0E-15_RKIND + length(k,iCell))
                 tauvVel(k,iCell) = slow_w_factor*C_slow*KE / (1.0E-15_RKIND + length(k,iCell))


!      CALCULATE W2 TENDENCY TERMS
!      Calculate the tendency (rate-of-change) terms for the w2 budget. These terms
!      are separated based on the physical processes that cause each tendency term.
!      The w2 tendency is the sum of the following terms (which are saved individually)
!                   (1) entrainment and detrainment between plumes
!                   (2) turbulent transport of w2 [ -d(www)/dz ]
!                   (3) parameterized pressure-strain terms (Rotta, rapid, buoyancy, Stokes)
!                   (4) production of w2 by buoyancy fluxes (conversion to/from potential energy)
!                   (5) sub-plume scale effects (based on Lappen & Randall 2001)
!                   (6) splatting (calculated earlier; converts w2 to u2/v2 near surface)
!                   (7) Stokes drift production terms -2*[uw*d(ustokes)/dz + vw*d(vstokes)/dz]
!      The tendency term is truncated if config_adc_truncate_tend flag it switched on
!      The parameterization consists of a slow (Rotta) term [Rotta, 1951],
!      a rapid term that depends on gradients in mean velocity [Speziale et al., 1991],
!      a buoyancy term [Umlauf & Burchard, 2005], and
!      a Stokes drift term based on Pearson et al., 2019.
!      Note that the slow pressure w2 term is included below so that w2tend3 an be compared directly to
!      pressure-strain terms, but it is removed  as it is included implicitly in the timestepping
                 w2tend1(k,iCell) = -wumd(k,iCell)**2.0_RKIND*(         &
                   Entrainment(k,iCell) + Detrainment(k,iCell))

                 w2tend2(k,iCell) = (McMid(k-1,iCell)*(1.0_RKIND - 2.0_RKIND* &
                   areaFractionMid(k-1,iCell))*wumdMid(k-1,iCell)**2.0_RKIND - McMid(k,iCell)* &
                   (1.0_RKIND - 2.0_RKIND*areaFractionMid(k,iCell))*wumdMid(k,iCell)**2.0_RKIND) / dzmid

                 w2tend3(k,iCell) = tauvVel(k,iCell)*(u2(i1,k,iCell) + v2(i1,k,iCell) &
                 -2.0_RKIND*w2(i1,k,iCell))/3.0_RKIND + (1.0_RKIND/3.0_RKIND*alpha_1 - &
                   alpha_2)*(uw(i1,k,iCell)*Uz + vw(i1,k,iCell)*Vz) - 4.0_RKIND/3.0_RKIND* &
                   C_b*Mc(k,iCell)*(grav*alphaT(k,iCell)*tumd(k,iCell) - &
                   grav*betaS(k,iCell)*sumd(k,iCell)) + &
                   (alpha_st_1 / 3.0_RKIND - alpha_st_2)*(uw(i1,k,iCell)*Ustz + vw(i1,k,iCell)*Vstz)

                 w2tend4(k,iCell) = 2.0_RKIND*Mc(k,iCell)* &
                   (grav*alphaT(k,iCell)*tumd(k,iCell) - grav*betaS(k,iCell)*sumd(k,iCell))

                 w2tend5(k,iCell) = Mc(k,iCell)*(Swumd(k-1,iCell) + Swumd(k,iCell))

                 w2tend7(k,iCell) = -2.0_RKIND*(uw(i1,k,iCell)*Ustz + vw(i1,k,iCell)*Vstz)
!      Sum the above tendency terms to find the total w2 tendency.
!      Note this includes a splatting term calculated earlier (w2tend6).
!      Also, the w2 slow pressure term (within w2tend3) is removed here as it is
!      included implicitly in the timestepping
                 w2tend(i3_f,k,iCell) = w2tend1(k,iCell) + w2tend2(k,iCell) + &
                   w2tend3(k,iCell) + w2tend4(k,iCell) + w2tend5(k,iCell) + &
                   w2tend6(k,iCell) + w2tend7(k,iCell) +  &
                   2.0_RKIND*tauvVel(k,iCell)*w2(i1,k,iCell)/3.0_RKIND

                 if (config_adc_truncate_tend) then
                    w2tend(i3_f,k,iCell) = FLOAT (INT(w2tend(i3_f,k,iCell) * adcRound &
                      + 0.5_RKIND)) / adcRound
                 endif


!      CALCULATE WT TENDENCY TERMS
!      The wt tendency is the sum of the following terms (which are saved individually)
!                   (1) entrainment and detrainment between plumes
!                   (2) turbulent transport of wt [ -d(wtw)/dz ]
!                   (3) parameterized pressure-scalar terms (Rotta, rapid, buoyancy, Stokes)
!                   (4) gradient and buoyancy production terms of wt [ -w2*dT/dz + t*buoyancy]
!                   (5) dissipation of wt [ -kappa_FL*[d(wt)/dz]^2 ]
!                   (6) sub-plume scale effects (based on Lappen & Randall 2001)
!                   (7) Stokes drift production terms -[ut*d(ustokes)/dz + vt*d(vstokes)/dz]
!      The tendency term is truncated if config_adc_truncate_tend flag it switched on
                 wttend1(k,iCell) = -1.0_RKIND*(Entrainment(k,iCell) + Detrainment(k,iCell)) * &
                   wumd(k,iCell)*tumd(k,iCell)

                 wttend2(k,iCell) = -(w2t(k-1,iCell) - w2t(k,iCell)) / (zm(k-1,iCell) - zm(k,iCell))

                 wttend3(k,iCell) = - tau_tracer(k,iCell)*wt(i1,k,iCell) &
                    - C_b_tracer*areaFraction(k,iCell)*  &
                   (1.0_RKIND - areaFraction(k,iCell))*grav*(alphaT(k,iCell) &
                   *tumd(k,iCell)**2.0_RKIND - betaS(k,iCell)*tumd(k,iCell)*sumd(k,iCell)) + &
                   0.5_RKIND*(alpha_tracer1 - alpha_tracer2)*(ut(i1,k,iCell)*Uz + &
                   vt(i1,k,iCell)*Vz) &
                   + 0.5_RKIND*(alpha_st_tracer1 - alpha_st_tracer2)* &
                   (ut(i1,k,iCell)*Ustz + vt(i1,k,iCell)*Vstz)

                 wttend4(k,iCell) = - Mc(k,iCell)*wumd(k,iCell)*Tz + &
                   areaFraction(k,iCell)*(1.0_RKIND - areaFraction(k,iCell))*grav*(alphaT(k,iCell) &
                   *tumd(k,iCell)**2.0_RKIND - betaS(k,iCell)*tumd(k,iCell)*sumd(k,iCell))

                 wttend5(k,iCell) = kappa_FL*(wt(i1,k-1,iCell) -       &
                   wt(i1,k+1,iCell)) / (ze(k-1,iCell) - ze(k+1,iCell))**2.0_RKIND

                 wttend6(k,iCell) = areaFraction(k,iCell)*(1.0_RKIND - areaFraction(k,iCell))* &
                   tumd(k,iCell)*0.5_RKIND*(Swumd(k-1,iCell) + Swumd(k,iCell)) - Mc(k,iCell)* &
                   (1.0_RKIND / areaFraction(k,iCell) * (areaFraction(k-1,iCell)*wt_spsU(k-1,iCell) - &
                   areaFraction(k+1,iCell)*wt_spsU(k+1,iCell)) / (ze(k-1,iCell) - ze(k+1,iCell))  &
                   -1.0_RKIND/(1.0_RKIND - areaFraction(k,iCell))*((1.0_RKIND - areaFraction(k-1,iCell))* &
                   wt_spsD(k-1,iCell) - (1.0_RKIND - areaFraction(k+1,iCell))*wt_spsD(k+1,iCell)) / &
                   (ze(k-1,iCell) - ze(k+1,iCell)))

                 wttend7(k,iCell) = -ut(i1,k,iCell)*Ustz - vt(i1,k,iCell)*Vstz
!      Sum the above tendency terms to find the total wt tendency.
!      Also, the wt slow pressure term (within wttend3) is removed here as it is
!      included implicitly in the timestepping
                 wttend(i3_f,k,iCell) = wttend1(k,iCell) + wttend2(k,iCell) +        &
                   wttend3(k,iCell) + wttend4(k,iCell) + wttend5(k,iCell) + wttend6(k,iCell) + &
                   wttend7(k,iCell) + tau_tracer(k,iCell)*wt(i1,k,iCell)

                 if (config_adc_truncate_tend) then
                    wttend(i3_f,k,iCell) = FLOAT (INT(wttend(i3_f,k,iCell) * adcRound &
                      + 0.5_RKIND)) / adcRound
                 endif


!      CALCULATE WS TENDENCY TERMS
!      The ws tendency is the sum of the following terms (which are saved individually)
!                   (1) entrainment and detrainment between plumes
!                   (2) turbulent transport of ws [ -d(wsw)/dz ]
!                   (3) parameterized pressure-scalar terms (Rotta, rapid, buoyancy, Stokes)
!                   (4) gradient and buoyancy production terms of ws [ -w2*dS/dz + s*buoyancy]
!                   (5) dissipation of ws [ -kappa_FL*[d(ws)/dz]^2 ]
!                   (6) sub-plume scale effects (based on Lappen & Randall 2001)
!                   (7) Stokes drift production terms -[us*d(ustokes)/dz + vs*d(vstokes)/dz]
!      The tendency term is truncated if config_adc_truncate_tend flag it switched on
                 wstend1(k,iCell) = -(Entrainment(k,iCell) + Detrainment(k,iCell)) * &
                   wumd(k,iCell)*sumd(k,iCell)

                 wstend2(k,iCell) = ((1.0_RKIND - 2.0_RKIND*areaFraction(k-1,iCell))* &
                   wumd(k-1,iCell)*sumd(k-1,iCell)*Mc(k-1,iCell) - (1.0_RKIND - 2.0_RKIND* &
                   areaFraction(k+1,iCell))*wumd(k+1,iCell)*sumd(k+1,iCell)*Mc(k+1,iCell)) &
                   / (ze(k-1,iCell) - ze(k+1,iCell))

                 wstend3(k,iCell) = - tau_tracer(k,iCell)*ws(i1,k,iCell) - &
                   C_b_tracer*grav*areaFraction(k,iCell)* &
                   (1.0_RKIND - areaFraction(k,iCell))*(alphaT(k,iCell) &
                   *tumd(k,iCell)*sumd(k,iCell) - betaS(k,iCell)*sumd(k,iCell)*sumd(k,iCell)) + &
                   0.5_RKIND*(alpha_tracer1 - alpha_tracer2)*(us(i1,k,iCell)*Uz + vs(i1,k,iCell)*Vz) &
                   + 0.5_RKIND*(alpha_st_tracer1 - alpha_st_tracer2)* &
                   (us(i1,k,iCell)*Ustz + vs(i1,k,iCell)*Vstz)

                 wstend4(k,iCell) = - Mc(k,iCell)*wumd(k,iCell)*Sz + &
                    grav*areaFraction(k,iCell)* &
                   (1.0_RKIND - areaFraction(k,iCell))*(alphaT(k,iCell) &
                   *tumd(k,iCell)*sumd(k,iCell) - betaS(k,iCell)*sumd(k,iCell)*sumd(k,iCell))

                 wstend5(k,iCell) = kappa_FL*(ws(i1,k-1,iCell) -       &
                   ws(i1,k+1,iCell)) / (ze(k-1,iCell) - ze(k+1,iCell))**2.0_RKIND

                 wstend6(k,iCell) = areaFraction(k,iCell)*(1.0_RKIND - areaFraction(k,iCell))* &
                   sumd(k,iCell)*0.5_RKIND*(Swumd(k-1,iCell) + Swumd(k,iCell)) - Mc(k,iCell)* &
                   (1.0_RKIND / areaFraction(k,iCell) * (areaFraction(k-1,iCell)*ws_spsU(k-1,iCell) - &
                   areaFraction(k+1,iCell)*ws_spsU(k+1,iCell)) / (ze(k-1,iCell) - ze(k+1,iCell))  &
                   -1.0_RKIND/(1.0_RKIND - areaFraction(k,iCell))*((1.0_RKIND - areaFraction(k-1,iCell))* &
                   ws_spsD(k-1,iCell) - (1.0_RKIND - areaFraction(k+1,iCell))*ws_spsD(k+1,iCell)) / &
                   (ze(k-1,iCell) - ze(k+1,iCell)))

                 wstend7(k,iCell) = -us(i1,k,iCell)*Ustz - vs(i1,k,iCell)*Vstz
!      Sum the above tendency terms to find the total ws tendency.
!      Also, the ws slow pressure term (within wstend3) is removed here as it is
!      included implicitly in the timestepping
                 wstend(i3_f,k,iCell) = wstend1(k,iCell) + wstend2(k,iCell) + &
                   wstend3(k,iCell) + wstend4(k,iCell) + wstend5(k,iCell) + wstend6(k,iCell) + &
                   wstend7(k,iCell) + tau_tracer(k,iCell)*ws(i1,k,iCell)

                 if (config_adc_truncate_tend) then
                    wstend(i3_f,k,iCell) = FLOAT (INT(wstend(i3_f,k,iCell) * adcRound &
                      + 0.5_RKIND)) / adcRound
                 endif


!      CALCULATE UW TENDENCY TERMS
!      The uw tendency is the sum of the following terms (which are saved individually)
!                   (1) turbulent transport of uw [ -d(uww)/dz ]
!                   (2) parameterized pressure-strain terms (Rotta, rapid, buoyancy, Stokes)
!                   (3) shear production of uw [-w2*dU/dz]
!                   (4) buoyancy production of uw [ u*buoyancy]
!                   (5) dissipation of uw [ -kappa_FL*[d(uw)/dz]^2 ]
!                   (6) Stokes production term -[u2*d(ustokes)/dz + uv*d(vstokes)/dz]
!      The tendency term is truncated if config_adc_truncate_tend flag it switched on
                 uwtend1(k,iCell) = -(uw2(k-1,iCell) - uw2(k,iCell)) / dzmid

                 uwtend2(k,iCell) = -tauVel(k,iCell)*uw(i1,k,iCell) + &
                   0.5_RKIND*((alpha_0 - 4.0_RKIND*alpha_1/3.0_RKIND)*KE**2.0_RKIND +  &
                   (alpha_1 - alpha_2)*u2(i1,k,iCell) + (alpha_1 +  &
                   alpha_2)*w2(i1,k,iCell))*Uz + 0.5_RKIND*(alpha_1 - alpha_2)*    &
                   uv(i1,k,iCell)*Vz - C_b*grav*(alphaT(k,iCell)*   &
                   ut(i1,k,iCell) - betaS(k,iCell)*us(i1,k,iCell)) + &
                   0.5_RKIND*((alpha_st_0 - 4.0_RKIND*alpha_st_1/3.0_RKIND)*KE**2.0_RKIND +  &
                   (alpha_st_1 - alpha_st_2)*u2(i1,k,iCell) + (alpha_st_1 +  &
                   alpha_st_2)*w2(i1,k,iCell))*Ustz + 0.5_RKIND*(alpha_st_1 - alpha_st_2)*    &
                   uv(i1,k,iCell)*Vstz

                 uwtend3(k,iCell) = - w2(i1,k,iCell)*Uz

                 uwtend4(k,iCell) = grav*(alphaT(k,iCell)*   &
                   ut(i1,k,iCell) - betaS(k,iCell)*us(i1,k,iCell))

                 uwtend5(k,iCell) = kappa_FL*(uw(i1,k-1,iCell) - uw(i1,k+1,iCell)) / &
                   (ze(k-1,iCell) - ze(k+1,iCell))**2.0_RKIND

                 uwtend6(k,iCell) = -u2(i1,k,iCell)*Ustz - uv(i1,k,iCell)*Vstz
!      Sum the above tendency terms to find the total uw tendency.
                 uwtend(i3_f,k,iCell) = uwtend1(k,iCell) + uwtend2(k,iCell) + &
                   uwtend3(k,iCell) + uwtend4(k,iCell) + uwtend5(k,iCell) + uwtend6(k,iCell)

                 if (config_adc_truncate_tend) then
                    uwtend(i3_f,k,iCell) = FLOAT (INT(uwtend(i3_f,k,iCell) * adcRound &
                      + 0.5_RKIND)) / adcRound
                 endif


!      CALCULATE VW TENDENCY TERMS
!      The vw tendency is the sum of the following terms (which are saved individually)
!                   (1) turbulent transport of vw [ -d(vww)/dz ]
!                   (2) parameterized pressure-strain terms (Rotta, rapid, buoyancy, Stokes)
!                   (3) shear production of vw [-w2*dV/dz]
!                   (4) buoyancy production of vw [ v*buoyancy]
!                   (5) dissipation of vw [ -kappa_FL*[d(vw)/dz]^2 ]
!                   (6) Stokes production term -[uv*d(ustokes)/dz + v2*d(vstokes)/dz]
!      The tendency term is truncated if config_adc_truncate_tend flag it switched on
                 vwtend1(k,iCell) = -(vw2(k-1,iCell) - vw2(k,iCell)) / dzmid

                 vwtend2(k,iCell) = -tauVel(k,iCell)*vw(i1,k,iCell) + &
                   0.5_RKIND*((alpha_0 - 4.0_RKIND*alpha_1/3.0_RKIND)*KE**2.0_RKIND +  &
                   (alpha_1 - alpha_2)*v2(i1,k,iCell) + (alpha_1 +  &
                   alpha_2)*w2(i1,k,iCell))*Vz + 0.5_RKIND*(alpha_1 - alpha_2)*    &
                   uv(i1,k,iCell)*Uz - C_b*grav*(alphaT(k,iCell)*   &
                   vt(i1,k,iCell) - betaS(k,iCell)*vs(i1,k,iCell))+ &
                   0.5_RKIND*((alpha_st_0 - 4.0_RKIND*alpha_st_1/3.0_RKIND)*KE**2.0_RKIND +  &
                   (alpha_st_1 - alpha_st_2)*v2(i1,k,iCell) + (alpha_st_1 +  &
                   alpha_st_2)*w2(i1,k,iCell))*Vstz + 0.5_RKIND*(alpha_st_1 - alpha_st_2)*    &
                   uv(i1,k,iCell)*Ustz

                 vwtend3(k,iCell) = - w2(i1,k,iCell)*Vz

                 vwtend4(k,iCell) = grav*(alphaT(k,iCell)*   &
                   vt(i1,k,iCell) - betaS(k,iCell)*vs(i1,k,iCell))

                 vwtend5(k,iCell) = kappa_FL*(vw(i1,k-1,iCell) - vw(i1,k+1,iCell)) / &
                   (ze(k-1,iCell) - ze(k+1,iCell))**2.0_RKIND

                 vwtend6(k,iCell) = -uv(i1,k,iCell)*Ustz - v2(i1,k,iCell)*Vstz
!      Sum the above tendency terms to find the total vw tendency.
                 vwtend(i3_f,k,iCell) = vwtend1(k,iCell) + vwtend2(k,iCell) + &
                   vwtend3(k,iCell) + vwtend4(k,iCell) + vwtend5(k,iCell) + vwtend6(k,iCell)

                 if (config_adc_truncate_tend) then
                    vwtend(i3_f,k,iCell) = FLOAT (INT(vwtend(i3_f,k,iCell) * adcRound &
                      + 0.5_RKIND)) / adcRound
                 endif


!      CALCULATE UV TENDENCY TERMS (these are not saved individually)
!      The uv tendency is the sum of
!                   (1) turbulent transport of vw [ -d(vww)/dz ]
!                   (2) parameterized pressure-strain terms (Rotta, rapid, Stokes)
!                   (3) shear production of uv [uw*dV/dz + vw*dU/dz]
!                   (4) dissipation of uv [ -kappa_FL*[d(uv)/dz]^2 ]
!      The tendency term is truncated if config_adc_truncate_tend flag it switched on
                uvtend(i3_f,k,iCell) = -(uvw(k-1,iCell) - uvw(k,iCell)) / dz & ! transport
                   - tauVel(k,iCell)*uv(i1,k,iCell)  & ! Rotta
                   + 0.5_RKIND*(alpha_1+alpha_2)* &
                   (uw(i1,k,iCell)*Vz + vw(i1,k,iCell)*Uz) & ! rapid
                   + 0.5_RKIND*(alpha_st_1 + alpha_st_2)* &
                   (uw(i1,k,iCell)*Vstz + vw(i1,k,iCell)*Ustz) & ! Stokes
                   - (uw(i1,k,iCell)*Vz + vw(i1,k,iCell)*Uz) & ! shear production
                   + kappa_VAR*(uv(i1,k-1,iCell) - uv(i1,k+1,iCell)) / &
                   (ze(k-1,iCell) - ze(k+1,iCell))**2.0_RKIND ! dissipation

!      Truncate the uv tendency term if the truncation flag is switched on
                 if (config_adc_truncate_tend) then
                    uvtend(i3_f,k,iCell) = FLOAT (INT(uvtend(i3_f,k,iCell) * adcRound &
                      + 0.5_RKIND)) / adcRound
                 endif


!      CALCULATE U2 TENDENCY TERMS
!      The u2 tendency is the sum of the following terms (which are saved individually)
!                   (1) turbulent transport of u2 [ -d(uuw)/dz ]
!                   (2) parameterized pressure-strain terms (Rotta, rapid, buoyancy, Stokes)
!                   (3) shear production of u2 [-2*uw*dU/dz]
!                   (4) dummy variable (should be zero)
!                   (5) dissipation of u2 [-(2/3)*dissipation]
!                   (6) splatting (calculated earlier)
!      The tendency term is truncated if config_adc_truncate_tend flag it switched on
                 u2tend1(k,iCell) = -(u2w(k-1,iCell) - u2w(k,iCell)) / dzmid

                 u2tend2(k,iCell) = tauVel(k,iCell)*(v2(i1,k,iCell) - 2.0_RKIND* &
                   u2(i1,k,iCell) + areaFraction(k,iCell)*(1.0_RKIND - &
                   areaFraction(k,iCell))*wumd(k,iCell)**2)/3.0_RKIND + &
                   (1.0_RKIND/3.0_RKIND*alpha_1 + alpha_2) &
                   *uw(i1,k,iCell)*Uz - 2.0_RKIND/3.0_RKIND*alpha_1*vw(i1,k,iCell)*Vz &
                   + 2.0_RKIND/3.0_RKIND*C_b*B + (alpha_st_1/3.0_RKIND + alpha_st_2) &
                   *uw(i1,k,iCell)*Ustz &
                   - 2.0_RKIND*alpha_st_1/3.0_RKIND*vw(i1,k,iCell)*Vstz

                 u2tend3(k,iCell) = -2.0_RKIND*uw(i1,k,iCell)*Uz

                 u2tend4(k,iCell) = 0.0_RKIND

                 u2tend5(k,iCell) = -2.0_RKIND/3.0_RKIND*eps(i1,k,iCell)
!      Sum the above tendency terms to find the total u2 tendency.
!      Note, the u2 slow pressure term (within u2tend3) is removed here as it is
!      included implicitly in the timestepping
                 u2tend(i3_f,k,iCell) = u2tend1(k,iCell) + u2tend2(k,iCell) + &
                   u2tend3(k,iCell) + u2tend4(k,iCell) + u2tend5(k,iCell) + &
                   u2tend6(k,iCell) + 2.0_RKIND*tauVel(k,iCell)*u2(i1,k,iCell)/3.0_RKIND

                 if (config_adc_truncate_tend) then
                    u2tend(i3_f,k,iCell) = FLOAT (INT(u2tend(i3_f,k,iCell) * adcRound &
                      + 0.5_RKIND)) / adcRound
                 endif


!      CALCULATE V2 TENDENCY TERMS
!      The v2 tendency is the sum of the following terms (which are saved individually)
!                   (1) turbulent transport of v2 [ -d(vvw)/dz ]
!                   (2) parameterized pressure-strain terms (Rotta, rapid, buoyancy, Stokes)
!                   (3) shear production of v2 [-2*vw*dV/dz]
!                   (4) dummy variable (should be zero)
!                   (5) dissipation of v2 [-(2/3)*dissipation]
!                   (6) splatting (calculated earlier)
!      The tendency term is truncated if config_adc_truncate_tend flag it switched on
                 v2tend1(k,iCell) = -(v2w(k-1,iCell) - v2w(k,iCell)) / dzmid

                 v2tend2(k,iCell) = tauVel(k,iCell)*(u2(i1,k,iCell) - 2.0_RKIND* &
                   v2(i1,k,iCell) + areaFraction(k,iCell)*(1.0_RKIND - &
                   areaFraction(k,iCell))*wumd(k,iCell)**2)/3.0_RKIND + &
                   (1.0_RKIND/3.0_RKIND*alpha_1 + alpha_2) &
                   *vw(i1,k,iCell)*Vz - 2.0_RKIND/3.0_RKIND*alpha_1*uw(i1,k,iCell)*Uz &
                   + 2.0_RKIND/3.0_RKIND*C_b*B + (alpha_st_1/3.0_RKIND + alpha_st_2) &
                   *vw(i1,k,iCell)*Vstz &
                   - 2.0_RKIND*alpha_st_1/3.0_RKIND*uw(i1,k,iCell)*Ustz

                 v2tend3(k,iCell) = -2.0_RKIND*vw(i1,k,iCell)*Vz

                 v2tend4(k,iCell) = 0.0_RKIND

                 v2tend5(k,iCell) = -2.0_RKIND/3.0_RKIND*eps(i1,k,iCell)
!      Sum the above tendency terms to find the total v2 tendency.
!      Note, the v2 slow pressure term (within v2tend3) is removed here as it is
!      included implicitly in the timestepping
                 v2tend(i3_f,k,iCell) = v2tend1(k,iCell) + v2tend2(k,iCell) + &
                   v2tend3(k,iCell) + v2tend4(k,iCell) + v2tend5(k,iCell) + &
                   v2tend6(k,iCell) + 2.0_RKIND*tauVel(k,iCell)*v2(i1,k,iCell)/3.0_RKIND

                 if (config_adc_truncate_tend) then
                    v2tend(i3_f,k,iCell) = FLOAT (INT(v2tend(i3_f,k,iCell) * adcRound &
                      + 0.5_RKIND)) / adcRound
                 endif
                

!      CALCULATE UT TENDENCY TERMS (these are not saved individually)
!      The ut tendency is the sum of
!                   (1) turbulent transport of ut [ -d(uwt)/dz ]
!                   (2) parameterized pressure-strain terms (Rotta, rapid, Stokes)
!                   (3) shear production of ut [-wt*dU/dz]
!                   (4) gradient production of ut [-uw*dT/dz]
!      The tendency term is truncated if config_adc_truncate_tend flag it switched on
                 uttend(i3_f,k,iCell) = -(uwt(k-1,iCell) - uwt(k,iCell))/dz  & ! Transport
                   - ut(i1,k,iCell)*tau_tracer(k,iCell) & ! Rotta
                   + 0.5_RKIND*(alpha_tracer1 + alpha_tracer2)*wt(i1,k,iCell)*Uz & ! Rapid
                   + 0.5_RKIND*(alpha_st_tracer1 + alpha_st_tracer2) &
                    *wt(i1,k,iCell)*Ustz & ! Stokes
                   - wt(i1,k,iCell)*Uz & ! Shear production
                   - uw(i1,k,iCell)*Tz ! Gradient production

                 if (config_adc_truncate_tend) then
                    uttend(i3_f,k,iCell) = FLOAT (INT(uttend(i3_f,k,iCell) * adcRound &
                      + 0.5_RKIND)) / adcRound
                 endif


!      CALCULATE VT TENDENCY TERMS (these are not saved individually)
!      The vt tendency is the sum of
!                   (1) turbulent transport of vt [ -d(vwt)/dz ]
!                   (2) parameterized pressure-strain terms (Rotta, rapid, Stokes)
!                   (3) shear production of vt [-wt*dV/dz]
!                   (4) gradient production of vt [-vw*dT/dz]
!      The tendency term is truncated if config_adc_truncate_tend flag it switched on
                 vttend(i3_f,k,iCell) = -(vwt(k-1,iCell) - vwt(k,iCell))/dz  & ! Transport
                   - vt(i1,k,iCell)*tau_tracer(k,iCell) & ! Rotta
                   + 0.5_RKIND*(alpha_tracer1 + alpha_tracer2)*wt(i1,k,iCell)*Vz & ! Rapid
                   + 0.5_RKIND*(alpha_st_tracer1 + alpha_st_tracer2) &
                    *wt(i1,k,iCell)*Vstz & ! Stokes
                   - wt(i1,k,iCell)*Vz & ! Shear production
                   - vw(i1,k,iCell)*Tz ! Gradient production

                 if (config_adc_truncate_tend) then
                    vttend(i3_f,k,iCell) = FLOAT (INT(vttend(i3_f,k,iCell) * adcRound &
                      + 0.5_RKIND)) / adcRound
                 endif

!      CALCULATE US TENDENCY TERMS (these are not saved individually)
!      The us tendency is the sum of
!                   (1) turbulent transport of us [ -d(uws)/dz ]
!                   (2) parameterized pressure-strain terms (Rotta, rapid, Stokes)
!                   (3) shear production of us [-ws*dU/dz]
!                   (4) gradient production of us [-uw*dS/dz]
!      The tendency term is truncated if config_adc_truncate_tend flag it switched on
                 ustend(i3_f,k,iCell) = -(uws(k-1,iCell) - uws(k,iCell))/dz  & ! Transport
                   - us(i1,k,iCell)*tau_tracer(k,iCell) & ! Rotta
                   + 0.5_RKIND*(alpha_tracer1 + alpha_tracer2)*ws(i1,k,iCell)*Uz & ! Rapid
                   + 0.5_RKIND*(alpha_st_tracer1 + alpha_st_tracer2) &
                    *ws(i1,k,iCell)*Ustz & ! Stokes
                   - ws(i1,k,iCell)*Uz & ! Shear production
                   - uw(i1,k,iCell)*Sz ! Gradient production

                 if (config_adc_truncate_tend) then
                    ustend(i3_f,k,iCell) = FLOAT (INT(ustend(i3_f,k,iCell) * adcRound &
                      + 0.5_RKIND)) / adcRound
                 endif


!      CALCULATE VS TENDENCY TERMS
!      The vs tendency is the sum of the following terms (which are NOT saved individually)
!                   (1) turbulent transport of vs [ -d(vws)/dz ]
!                   (2) parameterized pressure-strain terms (Rotta, rapid, Stokes)
!                   (3) shear production of vs [-ws*dV/dz]
!                   (4) gradient production of vs [-vw*dS/dz]
!      The tendency term is truncated if config_adc_truncate_tend flag it switched on
                 vstend(i3_f,k,iCell) = -(vws(k-1,iCell) - vws(k,iCell))/dz  & ! Transport
                   - vs(i1,k,iCell)*tau_tracer(k,iCell) & ! Rotta
                   + 0.5_RKIND*(alpha_tracer1 + alpha_tracer2)*ws(i1,k,iCell)*Vz & ! Rapid
                   + 0.5_RKIND*(alpha_st_tracer1 + alpha_st_tracer2) &
                    *ws(i1,k,iCell)*Vstz & ! Stokes
                   - ws(i1,k,iCell)*Vz & ! Shear production
                   - vw(i1,k,iCell)*Sz ! Gradient production

                 if (config_adc_truncate_tend) then
                    vstend(i3_f,k,iCell) = FLOAT (INT(vstend(i3_f,k,iCell) * adcRound &
                      + 0.5_RKIND)) / adcRound
                 endif



                 t2(i2,k,iCell) = tumd(k,iCell)**2.0_RKIND*areaFraction(k,iCell)*&
                   (1.0_RKIND-areaFraction(k,iCell))
                 s2(i2,k,iCell) = sumd(k,iCell)**2.0_RKIND*areaFraction(k,iCell)*&
                   (1.0_RKIND-areaFraction(k,iCell))
                 ts(i2,k,iCell) = tumd(k,iCell)*sumd(k,iCell)* &
                   areaFraction(k,iCell)*(1.0_RKIND-areaFraction(k,iCell))

                 if(.not. config_adc_use_old_length_scale)then
!     epsilon tendency is next
                    KEm1 = 0.5_RKIND*(u2(i1,k-1,iCell) + v2(i1,k-1,iCell) + w2(i1,k-1,iCell))
                    KE = 0.5_RKIND*(u2(i1,k,iCell) + v2(i1,k,iCell) + w2(i1,k,iCell))
                    KEp1 = 0.5_RKIND*(u2(i1,k+1,iCell) + v2(i1,k+1,iCell) + w2(i1,k+1,iCell))
                    tauM1 = 2.0_RKIND*KEm1 / (eps(i1,k-1,iCell) + 1.0E-15_RKIND)
                    tau = 2.0_RKIND*KE / (eps(i1,k,iCell) + 1.0E-15_RKIND)
                    tauP1 = 2.0_RKIND*KEp1 / (eps(i1,k+1,iCell) + 1.0E-15_RKIND)
                    tauAV = 0.5_RKIND*(tauUP + tauDN)
                    
                    tomUP = 0.5577_RKIND/1.3_RKIND*(0.5_RKIND*(KEm1+KE))**2/(0.5_RKIND*(eps(i1,k-1,iCell) + &
                      eps(i1,k,iCell)) + 1.0E-15_RKIND)*(eps(i1,k-1,iCell) - &
                      eps(i1,k,iCell)) / (ze(k-1,iCell) - ze(k,iCell) + 1.0E-10_RKIND)
                    tomDN = 0.5577_RKIND/1.3_RKIND*(0.5_RKIND*(KE+KEp1))**2/(0.5_RKIND*(eps(i1,k,iCell) + &
                      eps(i1,k+1,iCell)) + 1.0E-15_RKIND)*(eps(i1,k,iCell) - &
                      eps(i1,k+1,iCell)) / (ze(k,iCell) - ze(k+1,iCell) + 1.0E-10_RKIND)
                    
                    epstend(i3_f,k,iCell) = (tomUP - tomDN) / (zm(k-1,iCell) - zm(k,iCell) + 1.0E-10_RKIND) &
                      - 2.88_RKIND/(tau + 1.0E-10_RKIND)*(uw(i1,k,iCell)*Uz + vw(i1,k,iCell)*Vz) &
                      - 2.88_RKIND/(tau + 1.0E-10_RKIND)*(uw(i1,k,iCell)*Ustz + vw(i1,k,iCell)*Vstz) &
                      + 2.88_RKIND/(tau + 1.0E-10_RKIND)*grav*(alphaT(k,iCell)*wt(i1,k,iCell) - betaS(k,iCell)* &
                      ws(i1,k,iCell)) - 3.84_RKIND*eps(i1,k,iCell)/(tau + 1.0E-10_RKIND)
                 endif

!     finally update all subplume fluxes
                 if(BVF(k,iCell) <= 0.0_RKIND) then
                    lenspsU(k,iCell) = dzmid
                    lenspsD(k,iCell) = dzmid
                 else
                    lenspsU(k,iCell) = min(dz,0.76_RKIND*sqrt(KspsU(i1,k,iCell)/BVF(k,iCell)))
                    lenspsD(k,iCell) = min(dz,0.76_RKIND*sqrt(KspsD(i1,k,iCell)/BVF(k,iCell)))
                 endif

                 KmU(k,iCell) = 0.1_RKIND*lenspsU(k,iCell)*sqrt( KspsU(i1,k,iCell) )
                 KhU(k,iCell) = ( 1.0_RKIND+2.0_RKIND*lenspsU(k,iCell)/dz)*KmU(k,iCell)
                 wt_spsU(k,iCell) =  -KhU(k,iCell)*Tz
                 ws_spsU(k,iCell) =  -KhU(k,iCell)*Sz
                 
                 KmD(k,iCell) = 0.1_RKIND*lenspsD(k,iCell)*sqrt( KspsD(i1,k,iCell) )
                 KhD(k,iCell) = ( 1.0_RKIND+2.0_RKIND*lenspsD(k,iCell)/dz)*KmD(k,iCell)
                 wt_spsD(k,iCell) = -KhD(k,iCell)*Tz
                 ws_spsD(k,iCell) = -KhD(k,iCell)*Sz

                 if(k==2) then
                    Cval = 3.96_RKIND
                 else
                    Cval = (0.19_RKIND+0.51_RKIND*lenspsU(k,iCell)/dzmid)
                 endif

                 KspsUtend(i3_f,k,iCell) = grav*(alphaT(k,iCell)*wt_spsU(k,iCell) - &
                   betaS(k,iCell)*ws_spsU(k,ICell)) + ((KmU(k-1,iCell) +                 &
                   KmU(k,iCell))*(KspsU(i1,k-1,iCell) - KspsU(i1,k,iCell)) /         &
                   (ze(k-1,iCell) - ze(k,iCell)) - (KmU(k,iCell) +             &
                   KmU(k+1,iCell)) * (KspsU(i1,k,iCell) - KspsU(i1,k+1,iCell)) /     &
                   (ze(k,iCell) - ze(k+1,iCell))) / dzmid -                                 &
                   Cval*KspsU(i1,k,iCell)**1.5_RKIND/(1.0E-15_RKIND + lenspsU(k,iCell)) +                     &
                   eps(i1,k,iCell) / (2.0_RKIND*(1.0_RKIND - areaFraction(k,iCell))) + KmU(k,iCell)* &
                   (Uz**2.0_RKIND + Vz**2.0_RKIND + Uz*Ustz + Vz*Vstz)

                 if (config_adc_truncate_tend) then
                    KspsUtend(i3_f,k,iCell) = FLOAT (INT(KspsUtend(i3_f,k,iCell) * adcRound &
                      + 0.5_RKIND)) / adcRound
                 endif

                 if(k==2) then
                    Cval = 3.96_RKIND
                 else
                    Cval = (0.19_RKIND+0.51_RKIND*lenspsD(k,iCell)/dz)
                 endif

                 KspsDtend(i3_f,k,iCell) = grav*(alphaT(k,iCell)*wt_spsD(k,iCell) - &
                   betaS(k,iCell)*ws_spsD(k,iCell)) + ((KmD(k-1,iCell) +                 &
                   KmD(k,iCell))*(KspsD(i1,k-1,iCell) - KspsD(i1,k,iCell)) /         &
                   (ze(k-1,iCell) - ze(k,iCell)) - (KmD(k,iCell) +             &
                   KmD(k+1,iCell)) * (KspsD(i1,k,iCell) - KspsD(i1,k+1,iCell)) /     &
                   (ze(k,iCell) - ze(k+1,iCell))) / dz -                                 &
                   Cval*KspsD(i1,k,iCell)**1.5_RKIND /(1.0E-15_RKIND + lenspsD(k,iCell)) +                   &
                   eps(i1,k,iCell) / (2.0_RKIND*(areaFraction(k,iCell))) + &
                   KmD(k,iCell)*(Uz**2.0_RKIND + Vz**2.0_RKIND + Uz*Ustz + Vz*Vstz)

                 if (config_adc_truncate_tend) then
                    KspsDtend(i3_f,k,iCell) = FLOAT (INT(KspsDtend(i3_f,k,iCell) * adcRound &
                      + 0.5_RKIND)) / adcRound
                 endif

              enddo           ! nVertLevels
           enddo                ! nCells

!     In this step we update second moments except w3 which needs updated w2
           do iCell = 1,nCells
              do k=2,nVertLevels

!      STEP THE SECOND-MOMENTS FORWARD IN TIME
!      The tendency terms [w2tend etc.] are used to step the moments forward.
!      These tendencies are calculated for each moment (arbitrarily XY) using,
!
!      XYtend = d(XY)/dt \approx Sources and Sinks of XY
!
!      Then, XY_future = XY_past + f(XYtend at different times in past, timestep)
!
!      Most of the sources and sinks are included explicitly using the above method.
!      For the variances (w2, u2, v2) and vertical tracer fluxes (wt, ws), the
!      'return-to-isotropy' mechanism of pressure terms is included implicitly
!      instead. This term is parameterized as a sink of XY with
!               tendency \propto - A* XY/(turbulence timescale)
!      where A is a constant that can vary between each moment [Rotta, 1951].
!      The implicit inclusion moves this term from the right- to left-hand side of the
!      above equation. This results in a modified step-forward,
!
!      XY_future = XY_past + f(XYtendencies, timestep)/(1+A*timestep/turbulence_timescale)
!
!      Note that the tauvVel, tauVel, tau_tracer terms below are INVERSE timescales.

                 w2(i2,k,iCell) = (w2(i1,k,iCell) + dt_small*(Cw3 * w2tend(i3_f,k,iCell)  &
                   + Cw2 * w2tend(i2_f,k,iCell) + Cw1 * w2tend(i1_f,k,iCell))) /           &
                   (1.0_RKIND + dt_small*tauvVel(k,iCell)*2.0_RKIND/3.0_RKIND)
                 if(w2(i2,k,iCell) < epsilon) then
                    w2cliptend(k,iCell) = epsilon - w2(i2,k,iCell)
                    w2(i2,k,iCell) = epsilon
                 endif

                 if(abs(w2(i2,k,iCell)) > 1.0_RKIND) then
                    call mpas_log_write("ERROR: w2 out of range, w2 = $r, location = $i, $i", &
                      MPAS_LOG_CRIT,realArgs=(/w2(i2,k,iCell)/),intArgs=(/k,iCell/))
                 endif

                 u2(i2,k,iCell) = (u2(i1,k,iCell) + dt_small*(Cw3 * u2tend(i3_f,k,iCell) +  &
                   Cw2 * u2tend(i2_f,k,iCell) + Cw1 * u2tend(i1_f,k,iCell))) /          &
                   (1.0_RKIND + dt_small*tauVel(k,iCell)*2.0_RKIND/3.0_RKIND)
                 if(u2(i2,k,iCell) < epsilon) then
                    u2cliptend(k,iCell) = epsilon - u2(i2,k,iCell)
                    u2(i2,k,iCell) = epsilon
                 endif

                 v2(i2,k,iCell) = (v2(i1,k,iCell) + dt_small*(Cw3 * v2tend(i3_f,k,iCell) + &
                   Cw2 * v2tend(i2_f,k,iCell) + Cw1 * v2tend(i1_f,k,iCell))) /              &
                   (1.0_RKIND + dt_small*tauVel(k,iCell)*2.0_RKIND/3.0_RKIND)
                 if(v2(i2,k,iCell) < epsilon) then
                    v2cliptend(k,iCell) = epsilon - v2(i2,k,iCell)
                    v2(i2,k,iCell) = epsilon
                 endif

                 uw(i2,k,iCell) = uw(i1,k,iCell) + dt_small*(Cw3 * uwtend(i3_f,k,iCell) + &
                   Cw2 * uwtend(i2_f,k,iCell) + Cw1 * uwtend(i1_f,k,iCell))
                 vw(i2,k,iCell) = vw(i1,k,iCell) + dt_small*(Cw3 * vwtend(i3_f,k,iCell) + &
                   Cw2 * vwtend(i2_f,k,iCell) + Cw1 * vwtend(i1_f,k,iCell))
                 uv(i2,k,iCell) = uv(i1,k,iCell) + dt_small*(Cw3 * uvtend(i3_f,k,iCell) + &
                   Cw2 * uvtend(i2_f,k,iCell) + Cw1 * uvtend(i1_f,k,iCell))
                 ut(i2,k,iCell) = ut(i1,k,iCell) + dt_small*(Cw3 * uttend(i3_f,k,iCell) + &
                   Cw2 * uttend(i2_f,k,iCell) + Cw1 * uttend(i1_f,k,iCell))
                 wt(i2,k,iCell) = (wt(i1,k,iCell) + dt_small*(Cw3 * wttend(i3_f,k,iCell) + &
                   Cw2 * wttend(i2_f,k,iCell) + Cw1 * wttend(i1_f,k,iCell))) /              &
                   (1.0_RKIND + dt_small*tau_tracer(k,iCell))
                 vt(i2,k,iCell) = vt(i1,k,iCell) + dt_small*(Cw3 * vttend(i3_f,k,iCell) + &
                   Cw2 * vttend(i2_f,k,iCell) + Cw1 * vttend(i1_f,k,iCell))
                 us(i2,k,iCell) = us(i1,k,iCell) + dt_small*(Cw3 * ustend(i3_f,k,iCell) + &
                   Cw2 * ustend(i2_f,k,iCell) + Cw1 * ustend(i1_f,k,iCell))
                 vs(i2,k,iCell) = vs(i1,k,iCell) + dt_small*(Cw3 * vstend(i3_f,k,iCell) + &
                   Cw2 * vstend(i2_f,k,iCell) + Cw1 * vstend(i1_f,k,iCell))
                 ws(i2,k,iCell) = (ws(i1,k,iCell) + dt_small*(Cw3 * wstend(i3_f,k,iCell) + &
                   Cw2 * wstend(i2_f,k,iCell) + Cw1 * wstend(i1_f,k,iCell))) /              &
                   (1.0_RKIND + dt_small*tau_tracer(k,iCell))
                 if(abs(wt(i2,k,iCell)) > 1.0_RKIND) then
                    call mpas_log_write("ERROR: wt out of range, wt = $r, location = $i, $i", &
                      MPAS_LOG_CRIT,realArgs=(/wt(i2,k,iCell)/),intArgs=(/k,iCell/))
                 endif
                 
                 if(abs(ws(i2,k,iCell)) > 1.0_RKIND) then
                    call mpas_log_write("ERROR: ws out of range, ws = $r, location = $i, $i", &
                      MPAS_LOG_CRIT,realArgs=(/ws(i2,k,iCell)/),intArgs=(/k,iCell/))
                 endif
                 
                 if(abs(u2(i2,k,iCell)) > 1.0_RKIND) then
                    call mpas_log_write("ERROR: u2 out of range, u2 = $r, location = $i, $i", &
                      MPAS_LOG_CRIT,realArgs=(/u2(i2,k,iCell)/),intArgs=(/k,iCell/))
                 endif

                 if(abs(v2(i2,k,iCell)) > 1.0_RKIND) then
                    call mpas_log_write("ERROR: v2 out of range, v2 = $r, location = $i, $i", &
                      MPAS_LOG_CRIT,realArgs=(/v2(i2,k,iCell)/),intArgs=(/k,iCell/))
                 endif

                 if(.not. config_adc_use_old_length_scale)then
                    eps(i2,k,iCell) = eps(i1,k,iCell) + dt_small*(Cw3 * epstend(i3_f,k,iCell) + &
                      Cw2 * epstend(i2_f,k,iCell) + Cw1 * epstend(i1_f,k,iCell))
                    if (eps(i2,k,iCell) < 1.0E-10_RKIND) then
                       eps(i2,k,iCell) = 1.0E-10_RKIND
                    endif
!     Amrapalli: 16.6 here can be tuned if needed, would show up in u2 and v2 primarily
!     also in KEspsU/D
                    KE = 0.5_RKIND*(u2(i2,k,iCell) + v2(i2,k,iCell) + w2(i2,k,iCell))
                    length(k,iCell) = KE**1.5_RKIND / (c_epsilon * eps(i2,k,iCell))
                 else
                    KE = 0.5_RKIND*(u2(i2,k,iCell) + v2(i2,k,iCell) + w2(i2,k,iCell))
                    eps(:,k,iCell) = KE**1.5_RKIND / (c_epsilon * length(k,iCell))
                 endif
                 KspsU(i2,k,iCell) = max(epsilon,KspsU(i1,k,iCell) + dt_small*(Cw3 * &
                   KspsUtend(i3_f,k,iCell) + Cw2 * KspsUtend(i2_f,k,iCell) + Cw1 *    &
                   KspsUtend(i1_f,k,iCell)))
                 KspsD(i2,k,iCell) = max(epsilon,KspsD(i1,k,iCell) + dt_small*(Cw3 * &
                   KspsDtend(i3_f,k,iCell) + Cw2 * KspsDtend(i2_f,k,iCell) + Cw1 *    &
                   KspsDtend(i1_f,k,iCell)))

              enddo           !nVertLevels
           enddo                !nCells for second order moment tendencies

!     update the third order w3 and mean fields
           do iCell = 1,nCells
              do k=1,nVertLevels
                 w3check = (w2(i2,k,iCell) + w2(i2,k+1,iCell))**1.5_RKIND
                 if(w3check < 3e-14_RKIND) w3check = 1.0_RKIND
                 w3(i2,k,iCell) = max(-w3check,min((w3(i1,k,iCell) + dt_small*(Cw3 * w3tend(i3_f,k,iCell) + &
                   Cw2 * w3tend(i2_f,k,iCell) + Cw1 * w3tend(i1_f,k,iCell))) / &
                   (1.0_RKIND + 0.0_RKIND*dt_small*tauw3(k,iCell)),w3check))


                 if(abs(w3(i2,k,iCell)) > 1.0_RKIND) then
                    call mpas_log_write("ERROR: w3 out of range, w3 = $r, location = $i, $i", &
                      MPAS_LOG_CRIT,realArgs=(/w3(i2,k,iCell)/),intArgs=(/k,iCell/))
                 endif
              enddo
           enddo
           iterCount = iterCount + 1
           call get_array_pointers
           call get_weights
!     enddo !end iteration loop -- substepping is done.
!     now that substepping is done, apply computed fluxes to update mean fields.
!     you can collapse this loop too.
           do iCell = 1, nCells
              do k=1,nVertLevels
                 utemp = uvel(k,iCell)
                 vtemp = vvel(k,iCell)
                 uvel(k,iCell) = uvel(k,iCell) - dt_small*(uw(i1,k,iCell) - uw(i1,k+1,iCell)) /  &
                   (ze(k,iCell) - ze(k+1,iCell)) + dt_small*fCell(iCell)*vStokes(k,iCell) !+ dt_small*fCell(iCell)*vtemp

                 vvel(k,iCell) = vvel(k,iCell) - dt_small*(vw(i1,k,iCell) - vw(i1,k+1,iCell)) /  &
                   (ze(k,iCell) - ze(k+1,iCell)) - dt_small*fCell(iCell)*uStokes(k,iCell) !- dt_small*fCell(iCell)*utemp
                 
                 wtSumUp = wt(i1,k,iCell) + areaFraction(k,iCell)*wt_spsD(k,iCell) + &
                   (1.0_RKIND - areaFraction(k,iCell))*wt_spsU(k,iCell)
                 wtSumDn = wt(i1,k+1,iCell) + areaFraction(k+1,iCell)*wt_spsD(k+1,iCell) + &
                   (1.0_RKIND - areaFraction(k+1,iCell))*wt_spsU(k+1,iCell)
                 wsSumUp = ws(i1,k,iCell) + areaFraction(k,iCell)*ws_spsD(k,iCell) + &
                   (1.0_RKIND - areaFraction(k,iCell))*ws_spsU(k,iCell)
                 ! saltfix: spsD/spsU were swapped here (not as in wsSumUp/wtSumDn), so the
                 ! interface salt flux did not telescope; same fix as later upstream (MPAS-Model a0d295fb)
                 wsSumDn = ws(i1,k+1,iCell) + areaFraction(k+1,iCell)*ws_spsD(k+1,iCell) + &
                   (1.0_RKIND - areaFraction(k+1,iCell))*ws_spsU(k+1,iCell)
                 
                 activeTracers(1,k,iCell) = activeTracers(1,k,iCell) - dt_small*(wtSumUp - &
                   wtSumDn) / (ze(k,iCell) - ze(k+1,iCell))
                 activeTracers(2,k,iCell) = activeTracers(2,k,iCell) - dt_small*(wsSumUp - &
                   wsSumDn) / (ze(k,iCell) - ze(k+1,iCell))
              enddo
           enddo
        enddo
 end subroutine compute_ADC_tends

end module ocn_adc_mixing_fused
