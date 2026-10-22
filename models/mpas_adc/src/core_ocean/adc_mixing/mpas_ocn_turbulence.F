!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! \file mpas_ocn_turbulence.F
!
! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_turbulence
!
!>  \brief MPAS ocean turbulence arrays for ADC mixing structure
!!  \author Luke Van Roekel
!!  \date   14 Jan 2020
!!  \details
!!  This module creates and maintains a primary ocean mesh structure
!!  for the ADC mixing scheme.  Eventually it will send arrays over to GPU
!!  but for now just loads on to cpu
!
!-------------------------------------------------------------------------------

module ocn_turbulence

   ! module dependencies
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_log

   implicit none
   private

   !----------------------------------------------------------------------------
   !
   ! Public parameters
   !
   !----------------------------------------------------------------------------

   type, public :: adcMixingType

   real(kind=RKIND) :: epsilon, &! small value below which w2, v2, u2 cannot go
                       sigmat, Ko, gamma1, beta5, c1, c2, c8, c10, c11, &
                       alpha1, alpha2, alpha3, B1, Kt, grav, c_1, c_2, c_mom, &
                       c_therm, c_mom_w3, c_mom_w3, c_pt, c_ps, c_pv, &
                       kappa_FL, kappa_w3, kappa_VAR, Cww_D, Cww_E

   real(kind=RKIND), dimension(:,:), pointer :: KspsU, KspsD, ze, zm, &
                       eps, length, lenspsD, lenspsU, KhU, KhD, KmU,  &
                       KmD, wt_spsU, wt_spsD, ws_spsU, ws_spsD, uw2,  &
                       vw2, u2w, v2w, w2t, w2s, wts, uvw, uwt, vwt,   &
                       uws, vws, ws2, wt2, sigma, Entrainment,        &
                       Detrainment, tumd, sumd, wumd, Mc

   real(kind=RKIND), dimension(:,:), pointer :: w2tend1, w2tend2,     &
                       w2tend3, w2tend4, w2tend5, w2tend6, wttend1, wttend2,   &
                       wttend3, wttend4, wttend5, w3tend1, w3tend2,   &
                       w3tend3, w3tend5, w3tend4, wstend1, wstend2,   &
                       wstend3, wstend4, wstend5, uwtend1, uwtend2,   &
                       uwtend3, uwtend4, uwtend5, vwtend1, vwtend2,   &
                       vwtend3, vwtend4, vwtend5, u2tend1, u2tend2,   &
                       u2tend3, u2tend4, u2tend5, v2tend1, v2tend2,   &
                       v2tend3, v2tend4, v2tend5, u2cliptend,         &
                       v2cliptend, w2cliptend

   real(kind=RKIND), dimension(:,:,:), pointer :: u2, v2, w2, t2, s2, &
                       uw, vw, wt, ws, w3, uv, ut, vt, us, vs, ts

   end type adcMixingType

   type (adcMixingType) :: adcMixing

   !----------------------------------------------------------------------------
   !
   ! Public member functions
   !
   !----------------------------------------------------------------------------

   public ::                      &
      ocn_turbulenceCreate,       &
      ocn_turbulenceUpdateFields, &
      ocn_turbulenceDestroy

   subroutine ocn_turbulenceCreate(domain)

   type(domain_type) :: &
         domain                    !< [in] MPAS type to describe domain

   integer :: &
         blockCount               ! counter for number of blocks

   type(block_type), pointer :: &
         block                    ! variables in current subblock

   type (mpas_pool_type), pointer :: &
         adcDiagnosticsPool, adcPrognosticsPool, adcTendPool  ! mesh variables in MPAS pool structure

   real(kind=RKIND), dimension(:,:), pointer :: KspsUTmp, KspsDTmp, zeTmp, zmTmp, &
                       epsTmp, lengthTmp, lenspsDTmp, lenspsUTmp, KhUTmp, KhDTmp, KmUTmp,  &
                       KmDTmp, wt_spsUTmp, wt_spsDTmp, ws_spsUTmp, ws_spsDTmp, uw2Tmp,  &
                       vw2Tmp, u2wTmp, v2wTmp, w2tTmp, w2sTmp, wtsTmp, uvwTmp, uwtTmp, vwtTMP,   &
                       uwsTmp, vwsTmp, ws2Tmp, wt2Tmp, sigmaTmp, EntrainmentTmp,        &
                       DetrainmentTmp, tumdTmp, sumdTmp, wumdTmp, McTmp

   real(kind=RKIND), dimension(:,:), pointer :: w2tend1Tmp, w2tend2Tmp,     &
                       w2tend3Tmp, w2tend4Tmp, w2tend5Tmp, w2tend6Tmp, &
                       wttend1Tmp, wttend2Tmp,   &
                       wttend3Tmp, wttend4Tmp, wttend5Tmp, w3tend1Tmp, w3tend2Tmp,   &
                       w3tend3Tmp, w3tend5Tmp, w3tend4Tmp, wstend1Tmp, wstend2Tmp,   &
                       wstend3Tmp, wstend4Tmp, wstend5Tmp, uwtend1Tmp, uwtend2Tmp,   &
                       uwtend3Tmp, uwtend4Tmp, uwtend5Tmp, vwtend1Tmp, vwtend2Tmp,   &
                       vwtend3Tmp, vwtend4Tmp, vwtend5Tmp, u2tend1Tmp, u2tend2Tmp,   &
                       u2tend3Tmp, u2tend4Tmp, u2tend5Tmp, v2tend1Tmp, v2tend2Tmp,   &
                       v2tend3Tmp, v2tend4Tmp, v2tend5Tmp, u2cliptendTmp,         &
                       v2cliptendTmp, w2cliptendTmp

   real(kind=RKIND), dimension(:,:,:), pointer :: u2Tmp, v2Tmp, w2Tmp, t2Tmp, s2Tmp, &
                       uwTmp, vwTmp, wtTmp, wsTmp, w3Tmp, uvTmp, utTmp, vtTmp, usTmp,&
                       vsTmp, tsTmp

   block => domain % blocklist
   do while (associated(block))
      call mpas_pool_get_subpool(block % structs, 'adcDiagnosticArrays', adcDiagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'adcTendArrays', adcTendPool)
      call mpas_pool_get_subpool(block % structs, 'adcPrognosticArrays', adcPrognosticPool)
      adcMixing%sigmat = config_adc_sigmat
      adcMixing%epsilon = config_adc_epsilon
      adcMixing%Ko = config_adc_Ko
      adcMixing%gamma1 = config_adc_gamma1
      adcMixing%beta5 = config_adc_beta5
      adcMixing%c1 = 2.88_RKIND
      adcMixing%c2 = 3.84_RKIND
      adcMixing%c8 = 3.0_RKIND
      adcMixing%c10 = 5.0_RKIND
      adcMixing%c11 = 0.1_RKIND
      adcMixing%alpha1 = 0.984_RKIND
      adcMixing%alpha2 = 0.568_RKIND
      adcMixing%alpha3 = 0.29_RKIND
      adcMixing%B1 = 16.6_RKIND
      adcMixing%Kt = 0.4_RKIND
      adcMixing%c_1 = config_adc_c1
      adcMixing%c_2 = config_adc_c2
      adcMixing%c_mom = config_adc_Cmom
      adcMixing%c_therm = config_adc_Ctherm
      adcMixing%c_mom_w3 = config_adc_Cmom_w3
      adcMixing%c_pt = config_adc_cpt
      adcMixing%c_ps = config_adc_cps
      adcMixing%c_pv = config_adc_cpv
      adcMixing%Cww_D = config_adc_CwwD
      adcMixing%Cww_E = config_adc_CwwE
      adcMixing%kappa_FL = config_adc_kappaFL
      adcMixing%kappa_VAR = config_adc_kappaVAR
      adcMixing%kappa_w3 = config_adc_kappaW3

      call mpas_pool_get_array(adcDiagnosticsPool, 'KspsD', KspsDTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'KspsU', KspsUTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'ze', zeTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'zm', zmTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'eps', epsTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'length', lengthTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'lenspsD', lenspsDTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'lenspsU', lenspsUTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'KhU', KhUTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'KhD', KhDTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'KmU', KmUTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'KmD', KmDTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'wt_spsU', wt_spsUTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'wt_spsD', wt_spsDTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'ws_spsU', ws_spsUTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'ws_spsD', ws_spsDTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'uw2', uw2Tmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'vw2', vw2Tmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'u2w', u2wTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'v2w', v2wTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'w2t', w2tTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'w2s', w2sTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'wts', wtsTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'uvw', uvwTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'uwt', uwtTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'vwt', vwtTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'uws', uwsTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'vws', vwsTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'ws2', ws2Tmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'wt2', wt2Tmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'sigma', sigmaTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'Entrainment', EntrainmentTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'Detrainment', DetrainmentTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'tumd', tumdTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'sumd', sumdTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'wumd', wumdTmp)

      adcMixing%KspsD = KspsDTmp
      adcMixing%KspsU = KspsUTmp
      adcMixing%ze = zeTmp
      adcMixing%zm = zmTmp
      adcMixing%eps = epsTmp
      adcMixing%length = lengthTmp
      adcMixing%lenspsD = lenspsDTmp
      adcMixing%lenspsU = lenspsUTmp
      adcMixing%KhU = KhUTmp
      adcMixing%KhD = KhDTmp
      adcMixing%KmU = KmUTmp
      adcMixing%KmD = KmDTmp
      adcMixing%wt_spsU = wt_spsUTmp
      adcMixing%wt_spsD = wt_spsDTmp
      adcMixing%ws_spsU = ws_spsUTmp
      adcMixing%ws_spsD = ws_spsDTmp
      adcMixing%uw2 = uw2Tmp
      adcMixing%vw2 = vw2Tmp
      adcMixing%u2w = u2wTmp
      adcMixing%v2w = v2wTmp
      adcMixing%w2t = w2tTmp
      adcMixing%w2s = w2sTmp
      adcMixing%wts = wtsTmp
      adcMixing%uvw = uvwTmp
      adcMixing%uwt = uwtTmp
      adcMixing%vwt = vwtTmp
      adcMixing%uws = uwsTmp
      adcMixing%vws = vwsTmp
      adcMixing%ws2 = ws2Tmp
      adcMixing%wt2 = wt2Tmp
      adcMixing%sigma = simgmaTmp
      adcMixing%Entrainment = EntrainmentTmp
      adcMixing%Detrainment = DetrainmentTmp
      adcMixing%tumd = tumdTmp
      adcMixing%sumd = sumdTmp
      adcMixing%wumd = wumdTmp

      call mpas_pool_get_array(adcTendPool, 'w2tend1', w2tend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend2', w2tend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend3', w2tend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend4', w2tend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend5', w2tend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend6', w2tend6Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend1', w3tend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend2', w3tend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend3', w3tend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend4', w3tend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend5', w3tend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend1', wttend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend2', wttend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend3', wttend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend4', wttend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend5', wttend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend1', wstend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend2', wstend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend3', wstend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend4', wstend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend5', wstend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend1', uwtend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend2', uwtend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend3', uwtend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend4', uwtend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend5', uwtend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend1', vwtend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend2', vwtend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend3', vwtend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend4', vwtend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend5', vwtend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend1', u2tend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend2', u2tend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend3', u2tend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend4', u2tend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend5', u2tend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend1', v2tend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend2', v2tend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend3', v2tend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend4', v2tend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend5', v2tend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2cliptend', u2cliptendTmp)
      call mpas_pool_get_array(adcTendPool, 'v2cliptend', v2cliptendTmp)
      call mpas_pool_get_array(adcTendPool, 'w2cliptend', w2cliptendTmp)

      adcMixing%w2tend1 = w2tend1Tmp
      adcMixing%w2tend2 = w2tend2Tmp
      adcMixing%w2tend3 = w2tend3Tmp
      adcMixing%w2tend4 = w2tend4Tmp
      adcMixing%w2tend5 = w2tend5Tmp
      adcMixing%w2tend6 = w2tend6Tmp
      adcMixing%w3tend1 = w3tend1Tmp
      adcMixing%w3tend2 = w3tend2Tmp
      adcMixing%w3tend3 = w3tend3Tmp
      adcMixing%w3tend4 = w3tend4Tmp
      adcMixing%w3tend5 = w3tend5Tmp
      adcMixing%wttend1 = wttend1Tmp
      adcMixing%wttend2 = wttend2Tmp
      adcMixing%wttend3 = wttend3Tmp
      adcMixing%wttend4 = wttend4Tmp
      adcMixing%wttend5 = wttend5Tmp
      adcMixing%wstend1 = wstend1Tmp
      adcMixing%wstend2 = wstend2Tmp
      adcMixing%wstend3 = wstend3Tmp
      adcMixing%wstend4 = wstend4Tmp
      adcMixing%wstend5 = wstend5Tmp
      adcMixing%uwtend1 = uwtend1Tmp
      adcMixing%uwtend2 = uwtend2Tmp
      adcMixing%uwtend3 = uwtend3Tmp
      adcMixing%uwtend4 = uwtend4Tmp
      adcMixing%uwtend5 = uwtend5Tmp
      adcMixing%vwtend1 = vwtend1Tmp
      adcMixing%vwtend2 = vwtend2Tmp
      adcMixing%vwtend3 = vwtend3Tmp
      adcMixing%vwtend4 = vwtend4Tmp
      adcMixing%vwtend5 = vwtend5Tmp
      adcMixing%u2tend1 = u2tend1Tmp
      adcMixing%u2tend2 = u2tend2Tmp
      adcMixing%u2tend3 = u2tend3Tmp
      adcMixing%u2tend4 = u2tend4Tmp
      adcMixing%u2tend5 = u2tend5Tmp
      adcMixing%v2tend1 = v2tend1Tmp
      adcMixing%v2tend2 = v2tend2Tmp
      adcMixing%v2tend3 = v2tend3Tmp
      adcMixing%v2tend4 = v2tend4Tmp
      adcMixing%v2tend5 = v2tend5Tmp
      adcMixing%u2cliptend = u2cliptendTmp
      adcMixing%v2cliptend = v2cliptendTmp
      adcMixing%w2cliptend = w2cliptendTmp

      call mpas_pool_get_array(adcPrognosticPool, 'u2', u2Tmp)
      call mpas_pool_get_array(adcPrognosticPool, 'v2', v2Tmp)
      call mpas_pool_get_array(adcPrognosticPool, 'w2', w2Tmp)
      call mpas_pool_get_array(adcPrognosticPool, 'wt', wtTmp)
      call mpas_pool_get_array(adcPrognosticPool, 'ws', wsTmp)
      call mpas_pool_get_array(adcPrognosticPool, 'uw', uwTmp)
      call mpas_pool_get_array(adcPrognosticPool, 'vw', vwTmp)
      call mpas_pool_get_array(adcPrognosticPool, 'w3', w3Tmp)
      call mpas_pool_get_array(adcPrognosticPool, 'uv', uvTmp)
      call mpas_pool_get_array(adcPrognosticPool, 't2', t2Tmp)
      call mpas_pool_get_array(adcPrognosticPool, 'ut', utTmp)
      call mpas_pool_get_array(adcPrognosticPool, 'vt', vtTmp)
      call mpas_pool_get_array(adcPrognosticPool, 'us', usTmp)
      call mpas_pool_get_array(adcPrognosticPool, 'vs', vsTmp)
      call mpas_pool_get_array(adcPrognosticPool, 'ts', tsTmp)
      call mpas_pool_get_array(adcPrognosticPool, 's2', s2Tmp)

      adcMixing%u2 = u2Tmp
      adcMixing%v2 = v2Tmp
      adcMixing%w2 = w2Tmp
      adcMixing%wt = wtTmp
      adcMixing%ws = wsTmp
      adcMixing%uw = uwTmp
      adcMixing%vw = vwTmp
      adcMixing%w3 = w3Tmp
      adcMixing%uv = uvTmp
      adcMixing%t2 = t2Tmp
      adcMixing%ut = utTmp
      adcMixing%vt = vtTmp
      adcMixing%us = usTmp
      adcMixing%vs = vsTmp
      adcMixing%ts = tsTmp
      adcMixing%s2 = s2Tmp

!   at end could do a copy to accelerator

   end subroutine ocn_turbulenceCreate

   subroutine ocn_turbulenceUpdateFields

!   would only be needed to do device updates

   end subroutine ocn_turbulenceUpdateFields

   subroutine ocn_turbulenceDestroy(ierr)

!   deallocate stuff, set integers and constants to zero, nullify stuff on gpu

   end subroutine ocn_turbulenceDestroy

   end module ocn_turbulence
