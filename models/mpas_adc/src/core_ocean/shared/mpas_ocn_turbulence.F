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
   use mpas_constants
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_log

   use ocn_config

   implicit none
   private

   !----------------------------------------------------------------------------
   !
   ! Public parameters
   !
   !----------------------------------------------------------------------------

   !type, public :: adcMixingType

   real(kind=RKIND),public :: epsilon, &! small value below which w2, v2, u2 cannot go
                       sigmat, Ko, c_b, c_b_tracer, alpha_tracer1, alpha_tracer2, c11, &
                       alpha_st_tracer1, alpha_st_tracer2, alpha_st_0, alpha_st_1, &
                       alpha_st_2, alpha_0, alpha_1, alpha_2, B1, Kt, grav, c_mom, &
                       c_therm, c_mom_w3, c_epsilon, c_slow, c_slow_tracer, slow_w_factor, &
                       kappa_FL, kappa_w3, kappa_VAR, Cww_D, Cww_E, adcRound, &
                       La_t, stokes_depth, stokes_angle

   integer,public :: iterCount, nCells, nVertLevels

   real(kind=RKIND),public, dimension(:,:), allocatable :: ze, zm, &
                       epsSPS, lendn, lenup, length, lenspsD, lenspsU, KhU, KhD, KmU,  &
                       KmD, wt_spsU, wt_spsD, ws_spsU, ws_spsD, uw2,  &
                       vw2, u2w, v2w, w2t, w2s, wts, uvw, uwt, vwt,   &
                       uws, vws, ws2, wt2, areaFraction, Entrainment,        &
                       Detrainment, tumd, sumd, wumd, Mc

   real(kind=RKIND), public, dimension(:,:), allocatable :: w2tend1, w2tend2, w2tend3, &
                       w2tend4, w2tend5, w2tend6, w2tend7, wttend1, wttend2,   &
                       wttend3, wttend4, wttend5, wttend6, wttend7, w3tend1,   &
                       w3tend2, w3tend3, w3tend4, w3tend5, w3tend6, w3tend7,   &
                       wstend1, wstend2, wstend3, wstend4, wstend5, wstend6,   &
                       wstend7, uwtend1, uwtend2, uwtend3, uwtend4, uwtend5,   &
                       uwtend6, vwtend1, vwtend2, vwtend3, vwtend4, vwtend5,   &
                       vwtend6, u2tend1, u2tend2, u2tend3, u2tend4, u2tend5,   &
                       u2tend6, v2tend1, v2tend2, v2tend3, v2tend4, v2tend5,   &
                       v2tend6, u2cliptend, v2cliptend, w2cliptend

   real(kind=RKIND), public, dimension(:,:,:), allocatable :: u2, v2, w2, t2, s2, &
                       uw, vw, wt, ws, w3, uv, ut, vt, us, vs, ts, eps, KspsU, KspsD, &
                       u2tend, v2tend, w2tend, uwtend, vwtend, wttend, wstend, w3tend, &
                       uvtend, uttend, vttend, ustend, vstend, tstend, epstend, KspsUtend, &
                       KspsDtend

!   end type adcMixingType

   !ntype (adcMixingType), public :: adcMixing

   !----------------------------------------------------------------------------
   !
   ! Public member functions
   !
   !----------------------------------------------------------------------------

   public ::                      &
      ocn_turbulenceCreate,       &
      ocn_turbulenceUpdateFields, &
      ocn_turbulenceDestroy

   contains

   subroutine ocn_turbulenceCreate(domain)

   type(domain_type) :: &
         domain                    !< [in] MPAS type to describe domain

   integer :: &
         blockCount               ! counter for number of blocks

   type(block_type), pointer :: &
         block                    ! variables in current subblock

   type (mpas_pool_type), pointer :: &
         adcDiagnosticsPool, adcPrognosticsPool, adcTendPool, meshPool  ! mesh variables in MPAS pool structure

   real(kind=RKIND), dimension(:,:), pointer :: zeTmp, zmTmp, &
                       epsSPSTmp, lengthTmp, lenspsDTmp, lenspsUTmp, KhUTmp, KhDTmp, KmUTmp,  &
                       KmDTmp, wt_spsUTmp, wt_spsDTmp, ws_spsUTmp, ws_spsDTmp, uw2Tmp,  &
                       vw2Tmp, u2wTmp, v2wTmp, w2tTmp, w2sTmp, wtsTmp, uvwTmp, uwtTmp, vwtTMP,   &
                       uwsTmp, vwsTmp, ws2Tmp, wt2Tmp, areaFractionTmp, EntrainmentTmp,        &
                       DetrainmentTmp, tumdTmp, sumdTmp, wumdTmp, McTmp, lendnTmp, lenupTmp

   real(kind=RKIND), dimension(:,:), pointer :: w2tend1Tmp, w2tend2Tmp, w2tend3Tmp,  &
                       w2tend4Tmp, w2tend5Tmp, w2tend6Tmp, w2tend7Tmp, wttend1Tmp,   &
                       wttend2Tmp, wttend3Tmp, wttend4Tmp, wttend5Tmp, wttend6Tmp,   &
                       wttend7Tmp, w3tend1Tmp, w3tend2Tmp, w3tend3Tmp,               &
                       w3tend4Tmp, w3tend5Tmp, w3tend6Tmp, w3tend7Tmp, wstend1Tmp,   &
                       wstend2Tmp, wstend3Tmp, wstend4Tmp, wstend5Tmp, wstend6Tmp,   &
                       wstend7Tmp, uwtend1Tmp, uwtend2Tmp, uwtend3Tmp, uwtend4Tmp,   &
                       uwtend5Tmp, uwtend6Tmp, vwtend1Tmp, vwtend2Tmp, vwtend3Tmp,   &
                       vwtend4Tmp, vwtend5Tmp, vwtend6Tmp, u2tend1Tmp, u2tend2Tmp,   &
                       u2tend3Tmp, u2tend4Tmp, u2tend5Tmp, u2tend6Tmp, v2tend1Tmp,   &
                       v2tend2Tmp, v2tend3Tmp, v2tend4Tmp, v2tend5Tmp, v2tend6Tmp,   &
                       u2cliptendTmp, v2cliptendTmp, w2cliptendTmp

   real(kind=RKIND), dimension(:,:), pointer :: u2Tmp, v2Tmp, w2Tmp, t2Tmp, s2Tmp, &
                       uwTmp, vwTmp, wtTmp, wsTmp, w3Tmp, uvTmp, utTmp, vtTmp, usTmp,&
                       vsTmp, tsTmp, epsTmp, KspsUTmp, KspsDTmp, u2tendTmp, v2tendTmp, &
                       w2tendTmp, uwtendTmp, vwtendTmp, wttendTmp, wstendTmp, w3tendTmp, &
                       uvtendTmp, uttendTmp, vttendTmp, ustendTmp, vstendTmp, tstendTmp, &
                       epstendTmp, KspsUtendTmp, KspsDtendTmp

   integer, pointer :: iterCountTmp, nVertLevelsTmp, nVertLevelsP1
   integer :: i,nCellsTmp
   integer, dimension(:), pointer :: nCellsArray

   block => domain % blocklist
   do while (associated(block))
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevelsTmp)
      call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)

      nCellsTmp = nCellsArray( size(nCellsArray) )
      call mpas_pool_get_subpool(block % structs, 'adcDiagnosticArrays', adcDiagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'adcTendArrays', adcTendPool)
      call mpas_pool_get_subpool(block % structs, 'adcPrognosticArrays', adcPrognosticsPool)
      nCells = nCellsTmp
      nVertLevels = nVertLevelsTmp
      sigmat = config_adc_sigmat
      grav = gravity
      epsilon = config_adc_epsilon
      Ko = config_adc_Ko
      C_b_tracer = config_adc_c_b_tracer
      alpha_0 = config_adc_alpha_0
      alpha_1 = config_adc_alpha_1
      alpha_2 = config_adc_alpha_2
      alpha_tracer1 = config_adc_alpha_tracer1
      alpha_tracer2 = config_adc_alpha_tracer2
      alpha_st_0 = config_adc_alpha_st_0
      alpha_st_1 = config_adc_alpha_st_1
      alpha_st_2 = config_adc_alpha_st_2
      alpha_st_tracer1 = config_adc_alpha_st_tracer1
      alpha_st_tracer2 = config_adc_alpha_st_tracer2
      B1 = 16.6_RKIND
      Kt = 0.4_RKIND
      c_b = config_adc_c_b
      c11 = config_adc_c11
      c_mom = config_adc_Cmom
      c_therm = config_adc_Ctherm
      c_mom_w3 = config_adc_Cmom_w3
      c_slow_tracer = config_adc_c_slow_tracer
      c_slow = config_adc_c_slow
      slow_w_factor = config_adc_slow_w_factor
      c_epsilon = config_adc_dissipation_constant
      Cww_D = config_adc_CwwD
      Cww_E = config_adc_CwwE
      kappa_FL = config_adc_kappaFL
      kappa_VAR = config_adc_kappaVAR
      kappa_w3 = config_adc_kappaW3
      La_t = config_adc_La_t
      stokes_depth = config_adc_stokes_depth
      stokes_angle = config_adc_stokes_angle
      adcRound = 10.0_RKIND**config_adc_decimals_to_keep
      iterCount = 1

      call mpas_pool_get_array(adcDiagnosticsPool, 'ze', zeTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'zm', zmTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'epsSPS', epsSPSTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'length', lengthTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'lenup', lenupTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'lendn', lendnTmp)
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
      call mpas_pool_get_array(adcDiagnosticsPool, 'areaFraction', areaFractionTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'Entrainment', EntrainmentTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'Detrainment', DetrainmentTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'tumd', tumdTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'sumd', sumdTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'wumd', wumdTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'Mc', McTmp)

      allocate( ze( nVertLevelsP1, nCells ), zm( nVertLevels, nCells ))
      allocate( epsSPS( nVertLevelsP1, nCells ), length( nVertLevelsP1, nCells ))
      allocate( lendn( nVertLevelsP1, nCells ), lenup( nVertLevelsP1, nCells ))
      allocate( lenspsU( nVertLevelsP1, nCells), lenspsD( nVertLevelsP1, nCells))
      allocate( KhU( nVertLevelsP1, nCells), KhD( nVertLevelsP1, nCells))
      allocate( KmU( nVertLevelsP1, nCells), KmD( nVertLevelsP1, nCells))
      allocate( wt_spsU( nVertLevelsP1, nCells), wt_spsD( nVertLevelsP1, nCells))
      allocate( ws_spsU( nVertLevelsP1, nCells), ws_spsD( nVertLevelsP1, nCells))
      allocate( uw2( nVertLevels, nCells ), vw2( nVertLevels, nCells))
      allocate( u2w( nVertLevels, nCells ), v2w( nVertLevels, nCells ))
      allocate( w2t( nVertLevels, nCells ), w2s( nVertLevels, nCells ))
      allocate( wts( nVertLevels, nCells ), uvw( nVertLevels, nCells ))
      allocate( uwt( nVertLevels, nCells ), vwt( nVertLevels, nCells ))
      allocate( uws( nVertLevels, nCells ), vws( nVertLevels, nCells ))
      allocate( ws2( nVertLevels, nCells ), wt2( nVertLevels, nCells ))
      allocate( areaFraction( nVertLevelsP1, nCells ))
      allocate( Entrainment( nVertLevelsP1, nCells ), Detrainment( nVertLevelsP1, nCells ))
      allocate( tumd( nVertLevelsP1, nCells ), sumd( nVertLevelsP1, nCells ))
      allocate( wumd( nVertLevelsP1, nCells ), Mc( nVertLevelsP1, nCells ))

      ze = zeTmp
      zm = zmTmp
      epsSPS = epsSPSTmp
      length = lengthTmp
      lenup = lenupTmp
      lendn = lendnTmp
      lenspsD = lenspsDTmp
      lenspsU = lenspsUTmp
      KhU = KhUTmp
      KhD = KhDTmp
      KmU = KmUTmp
      KmD = KmDTmp
      wt_spsU = wt_spsUTmp
      wt_spsD = wt_spsDTmp
      ws_spsU = ws_spsUTmp
      ws_spsD = ws_spsDTmp
      uw2 = uw2Tmp
      vw2 = vw2Tmp
      u2w = u2wTmp
      v2w = v2wTmp
      w2t = w2tTmp
      w2s = w2sTmp
      wts = wtsTmp
      uvw = uvwTmp
      uwt = uwtTmp
      vwt = vwtTmp
      uws = uwsTmp
      vws = vwsTmp
      ws2 = ws2Tmp
      wt2 = wt2Tmp
      areaFraction = areaFractionTmp
      Entrainment = EntrainmentTmp
      Detrainment = DetrainmentTmp
      tumd = tumdTmp
      sumd = sumdTmp
      wumd = wumdTmp
      Mc = McTmp

      call mpas_pool_get_array(adcTendPool, 'w2tend1', w2tend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend2', w2tend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend3', w2tend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend4', w2tend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend5', w2tend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend6', w2tend6Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend7', w2tend7Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend1', w3tend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend2', w3tend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend3', w3tend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend4', w3tend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend5', w3tend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend6', w3tend6Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend7', w3tend7Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend1', wttend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend2', wttend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend3', wttend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend4', wttend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend5', wttend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend6', wttend6Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend7', wttend7Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend1', wstend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend2', wstend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend3', wstend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend4', wstend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend5', wstend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend6', wstend6Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend7', wstend7Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend1', uwtend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend2', uwtend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend3', uwtend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend4', uwtend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend5', uwtend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend6', uwtend6Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend1', vwtend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend2', vwtend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend3', vwtend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend4', vwtend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend5', vwtend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend6', vwtend6Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend1', u2tend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend2', u2tend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend3', u2tend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend4', u2tend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend5', u2tend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend6', u2tend6Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend1', v2tend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend2', v2tend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend3', v2tend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend4', v2tend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend5', v2tend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend6', v2tend6Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2cliptend', u2cliptendTmp)
      call mpas_pool_get_array(adcTendPool, 'v2cliptend', v2cliptendTmp)
      call mpas_pool_get_array(adcTendPool, 'w2cliptend', w2cliptendTmp)

      allocate( w2tend1( nVertLevelsP1, nCells ))
      allocate( w2tend2( nVertLevelsP1, nCells ))
      allocate( w2tend3( nVertLevelsP1, nCells ))
      allocate( w2tend4( nVertLevelsP1, nCells ))
      allocate( w2tend5( nVertLevelsP1, nCells ))
      allocate( w2tend6( nVertLevelsP1, nCells ))
      allocate( w2tend7( nVertLevelsP1, nCells ))
      allocate( wttend1( nVertLevelsP1, nCells ))
      allocate( wttend2( nVertLevelsP1, nCells ))
      allocate( wttend3( nVertLevelsP1, nCells ))
      allocate( wttend4( nVertLevelsP1, nCells ))
      allocate( wttend5( nVertLevelsP1, nCells ))
      allocate( wttend6( nVertLevelsP1, nCells ))
      allocate( wttend7( nVertLevelsP1, nCells ))
      allocate( wstend1( nVertLevelsP1, nCells ))
      allocate( wstend2( nVertLevelsP1, nCells ))
      allocate( wstend3( nVertLevelsP1, nCells ))
      allocate( wstend4( nVertLevelsP1, nCells ))
      allocate( wstend5( nVertLevelsP1, nCells ))
      allocate( wstend6( nVertLevelsP1, nCells ))
      allocate( wstend7( nVertLevelsP1, nCells ))
      allocate( uwtend1( nVertLevelsP1, nCells ))
      allocate( uwtend2( nVertLevelsP1, nCells ))
      allocate( uwtend3( nVertLevelsP1, nCells ))
      allocate( uwtend4( nVertLevelsP1, nCells ))
      allocate( uwtend5( nVertLevelsP1, nCells ))
      allocate( uwtend6( nVertLevelsP1, nCells ))
      allocate( vwtend1( nVertLevelsP1, nCells ))
      allocate( vwtend2( nVertLevelsP1, nCells ))
      allocate( vwtend3( nVertLevelsP1, nCells ))
      allocate( vwtend4( nVertLevelsP1, nCells ))
      allocate( vwtend5( nVertLevelsP1, nCells ))
      allocate( vwtend6( nVertLevelsP1, nCells ))
      allocate( v2tend1( nVertLevelsP1, nCells ))
      allocate( v2tend2( nVertLevelsP1, nCells ))
      allocate( v2tend3( nVertLevelsP1, nCells ))
      allocate( v2tend4( nVertLevelsP1, nCells ))
      allocate( v2tend5( nVertLevelsP1, nCells ))
      allocate( v2tend6( nVertLevelsP1, nCells ))
      allocate( u2tend1( nVertLevelsP1, nCells ))
      allocate( u2tend2( nVertLevelsP1, nCells ))
      allocate( u2tend3( nVertLevelsP1, nCells ))
      allocate( u2tend4( nVertLevelsP1, nCells ))
      allocate( u2tend5( nVertLevelsP1, nCells ))
      allocate( u2tend6( nVertLevelsP1, nCells ))
      allocate( w3tend1( nVertLevels, nCells ))
      allocate( w3tend2( nVertLevels, nCells ))
      allocate( w3tend3( nVertLevels, nCells ))
      allocate( w3tend4( nVertLevels, nCells ))
      allocate( w3tend5( nVertLevels, nCells ))
      allocate( w3tend6( nVertLevels, nCells ))
      allocate( w3tend7( nVertLevels, nCells ))
      allocate( u2cliptend( nVertLevelsP1, nCells ))
      allocate( v2cliptend( nVertLevelsP1, nCells ))
      allocate( w2cliptend( nVertLevelsP1, nCells ))

      w2tend1 = w2tend1Tmp
      w2tend2 = w2tend2Tmp
      w2tend3 = w2tend3Tmp
      w2tend4 = w2tend4Tmp
      w2tend5 = w2tend5Tmp
      w2tend6 = w2tend6Tmp
      w2tend7 = w2tend7Tmp
      w3tend1 = w3tend1Tmp
      w3tend2 = w3tend2Tmp
      w3tend3 = w3tend3Tmp
      w3tend4 = w3tend4Tmp
      w3tend5 = w3tend5Tmp
      w3tend6 = w3tend6Tmp
      w3tend7 = w3tend7Tmp
      wttend1 = wttend1Tmp
      wttend2 = wttend2Tmp
      wttend3 = wttend3Tmp
      wttend4 = wttend4Tmp
      wttend5 = wttend5Tmp
      wttend6 = wttend6Tmp
      wttend7 = wttend7Tmp
      wstend1 = wstend1Tmp
      wstend2 = wstend2Tmp
      wstend3 = wstend3Tmp
      wstend4 = wstend4Tmp
      wstend5 = wstend5Tmp
      wstend6 = wstend6Tmp
      wstend7 = wstend7Tmp
      uwtend1 = uwtend1Tmp
      uwtend2 = uwtend2Tmp
      uwtend3 = uwtend3Tmp
      uwtend4 = uwtend4Tmp
      uwtend5 = uwtend5Tmp
      uwtend6 = uwtend6Tmp
      vwtend1 = vwtend1Tmp
      vwtend2 = vwtend2Tmp
      vwtend3 = vwtend3Tmp
      vwtend4 = vwtend4Tmp
      vwtend5 = vwtend5Tmp
      vwtend6 = vwtend6Tmp
      u2tend1 = u2tend1Tmp
      u2tend2 = u2tend2Tmp
      u2tend3 = u2tend3Tmp
      u2tend4 = u2tend4Tmp
      u2tend5 = u2tend5Tmp
      u2tend6 = u2tend6Tmp
      v2tend1 = v2tend1Tmp
      v2tend2 = v2tend2Tmp
      v2tend3 = v2tend3Tmp
      v2tend4 = v2tend4Tmp
      v2tend5 = v2tend5Tmp
      v2tend6 = v2tend6Tmp
      u2cliptend = u2cliptendTmp
      v2cliptend = v2cliptendTmp
      w2cliptend = w2cliptendTmp

      allocate( KspsD( 3, nVertLevelsP1, nCells ) )
      allocate( KspsU( 3, nVertLevelsP1, nCells) )
      allocate( eps( 3, nVertLevelsP1, nCells ) )
      allocate( u2( 3, nVertLevelsP1, nCells ) )
      allocate( v2( 3, nVertLevelsP1, nCells ) )
      allocate( w2( 3, nVertLevelsP1, nCells ) )
      allocate( wt( 3, nVertLevelsP1, nCells ) )
      allocate( ws( 3, nVertLevelsP1, nCells ) )
      allocate( uw( 3, nVertLevelsP1, nCells ) )
      allocate( vw( 3, nVertLevelsP1, nCells ) )
      allocate( w3( 3, nVertLevels, nCells ) )
      allocate( uv( 3, nVertLevelsP1, nCells ) )
      allocate( t2( 3, nVertLevelsP1, nCells ) )
      allocate( ut( 3, nVertLevelsP1, nCells ) )
      allocate( vt( 3, nVertLevelsP1, nCells ) )
      allocate( us( 3, nVertLevelsP1, nCells ) )
      allocate( vs( 3, nVertLevelsP1, nCells ) )
      allocate( ts( 3, nVertLevelsP1, nCells ) )
      allocate( s2( 3, nVertLevelsP1, nCells ) )
      allocate( KspsDtend( 3, nVertLevelsP1, nCells ) )
      allocate( KspsUtend( 3, nVertLevelsP1, nCells) )
      allocate( epstend( 3, nVertLevelsP1, nCells ) )
      allocate( u2tend( 3, nVertLevelsP1, nCells ) )
      allocate( v2tend( 3, nVertLevelsP1, nCells ) )
      allocate( w2tend( 3, nVertLevelsP1, nCells ) )
      allocate( wttend( 3, nVertLevelsP1, nCells ) )
      allocate( wstend( 3, nVertLevelsP1, nCells ) )
      allocate( uwtend( 3, nVertLevelsP1, nCells ) )
      allocate( vwtend( 3, nVertLevelsP1, nCells ) )
      allocate( w3tend( 3, nVertLevels, nCells ) )
      allocate( uvtend( 3, nVertLevelsP1, nCells ) )
      allocate( uttend( 3, nVertLevelsP1, nCells ) )
      allocate( vttend( 3, nVertLevelsP1, nCells ) )
      allocate( ustend( 3, nVertLevelsP1, nCells ) )
      allocate( vstend( 3, nVertLevelsP1, nCells ) )
      allocate( tstend( 3, nVertLevelsP1, nCells ) )

      do i = 1,3
      call mpas_pool_get_array(adcPrognosticsPool, 'KspsD', KspsDTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'KspsU', KspsUTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'u2', u2Tmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'v2', v2Tmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'w2', w2Tmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'wt', wtTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'ws', wsTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'uw', uwTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'vw', vwTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'w3', w3Tmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'uv', uvTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 't2', t2Tmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'ut', utTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'vt', vtTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'us', usTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'vs', vsTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'ts', tsTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 's2', s2Tmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'eps', epsTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'KspsDtend', KspsDtendTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'KspsUtend', KspsUtendTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'u2tend', u2tendTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'v2tend', v2tendTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'w2tend', w2tendTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'wttend', wttendTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'wstend', wstendTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'uwtend', uwtendTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'vwtend', vwtendTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'w3tend', w3tendTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'uvtend', uvtendTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'uttend', uttendTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'vttend', vttendTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'ustend', ustendTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'vstend', vstendTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'tstend', tstendTmp, i)
      call mpas_pool_get_array(adcPrognosticsPool, 'epstend', epstendTmp, i)
      u2(i,:,:) = u2Tmp(:,:nCells)
      v2(i,:,:) = v2Tmp(:,:nCells)
      w2(i,:,:) = w2Tmp(:,:nCells)
      wt(i,:,:) = wtTmp(:,:nCells)
      ws(i,:,:) = wsTmp(:,:nCells)
      uw(i,:,:) = uwTmp(:,:nCells)
      vw(i,:,:) = vwTmp(:,:nCells)
      w3(i,:,:) = w3Tmp(:,:nCells)
      uv(i,:,:) = uvTmp(:,:nCells)
      t2(i,:,:) = t2Tmp(:,:nCells)
      ut(i,:,:) = utTmp(:,:nCells)
      vt(i,:,:) = vtTmp(:,:nCells)
      us(i,:,:) = usTmp(:,:nCells)
      vs(i,:,:) = vsTmp(:,:nCells)
      ts(i,:,:) = tsTmp(:,:nCells)
      s2(i,:,:) = s2Tmp(:,:nCells)
      eps(i,:,:) = epsTmp(:,:nCells)
      KspsU(i,:,:) = KspsUTmp(:,:nCells)
      KspsD(i,:,:) = KspsDTmp(:,:nCells)
      u2tend(i,:,:) = u2tendTmp(:,:nCells)
      v2tend(i,:,:) = v2tendTmp(:,:nCells)
      w2tend(i,:,:) = w2tendTmp(:,:nCells)
      wttend(i,:,:) = wttendTmp(:,:nCells)
      wstend(i,:,:) = wstendTmp(:,:nCells)
      uwtend(i,:,:) = uwtendTmp(:,:nCells)
      vwtend(i,:,:) = vwtendTmp(:,:nCells)
      w3tend(i,:,:) = w3tendTmp(:,:nCells)
      uvtend(i,:,:) = uvtendTmp(:,:nCells)
      uttend(i,:,:) = uttendTmp(:,:nCells)
      vttend(i,:,:) = vttendTmp(:,:nCells)
      ustend(i,:,:) = ustendTmp(:,:nCells)
      vstend(i,:,:) = vstendTmp(:,:nCells)
      tstend(i,:,:) = tstendTmp(:,:nCells)
      epstend(i,:,:) = epstendTmp(:,:nCells)
      KspsUtend(i,:,:) = KspsUtendTmp(:,:nCells)
      KspsDtend(i,:,:) = KspsDtendTmp(:,:nCells)
      end do
      block => block % next
   end do
!   at end could do a copy to accelerator
! need to copyin the parameters I guess...

      !$acc enter data copyin(KspsD,    &
      !$acc                   nCells,   &
      !$acc                   nVertLevels,  &
      !$acc                   sigmat,       &
      !$acc                   grav,         &
      !$acc                   epsilon,      &
      !$acc                   Ko,           &
      !$acc                   gamma1,       &
      !$acc                   beta5,        &
      !$acc                   c1,           &
      !$acc                   c2,           &
      !$acc                   c8,           &
      !$acc                   c10,          &
      !$acc                   c11,          &
      !$acc                   alpha1,       &
      !$acc                   alpha2,       &
      !$acc                   alpha3,       &
      !$acc                   B1,           &
      !$acc                   Kt,           &
      !$acc                   c_1,          &
      !$acc                   c_2,          &
      !$acc                   c_mom,        &
      !$acc                   c_therm,      &
      !$acc                   c_mom_w3,     &
      !$acc                   c_Temp,       &
      !$acc                   c_Salt,       &
      !$acc                   c_Vel,        &
      !$acc                   c_vVel,       &
      !$acc                   c_epsilon,    &
      !$acc                   Cww_D,        &
      !$acc                   Cww_E,        &
      !$acc                   kappa_FL,     &
      !$acc                   kappa_VAR,    &
      !$acc                   kappa_w3,     &
      !$acc                   iterCount,    &
      !$acc                  KspsU,    &
      !$acc                  ze,       &
      !$acc                  zm,       &
      !$acc                  epsSPS,      &
      !$acc                  length,   &
      !$acc                  lenup,    &
      !$acc                  lendn,    &
      !$acc                  lenspsD,  &
      !$acc                  lenspsU,  &
      !$acc                  KhU,      &
      !$acc                  KhD,      &
      !$acc                  KmU,      &
      !$acc                  KmD,      &
      !$acc                  wt_spsU,  &
      !$acc                  wt_spsD,  &
      !$acc                  ws_spsU,  &
      !$acc                  ws_spsD,  &
      !$acc                  uw2,      &
      !$acc                  vw2,      &
      !$acc                  u2w,      &
      !$acc                  v2w,      &
      !$acc                  w2t,      &
      !$acc                  w2s,      &
      !$acc                  wts,      &
      !$acc                  uvw,      &
      !$acc                  uwt,      &
      !$acc                  vwt,      &
      !$acc                  uws,      &
      !$acc                  vws,      &
      !$acc                  ws2,      &
      !$acc                  wt2,      &
      !$acc                  areaFraction,    &
      !$acc                  Entrainment, &
      !$acc                  Detrainment, &
      !$acc                  tumd,        &
      !$acc                  sumd,        &
      !$acc                  wumd,        &
      !$acc                  w2tend1,     &
      !$acc                  w2tend2,     &
      !$acc                  w2tend3,     &
      !$acc                  w2tend4,     &
      !$acc                  w2tend5,     &
      !$acc                  w2tend6,     &
      !$acc                  w3tend1,     &
      !$acc                  w3tend2,     &
      !$acc                  w3tend3,     &
      !$acc                  w3tend4,     &
      !$acc                  w3tend5,     &
      !$acc                  w3tend6,     &
      !$acc                  wttend1,     &
      !$acc                  wttend2,     &
      !$acc                  wttend3,     &
      !$acc                  wttend4,     &
      !$acc                  wttend5,     &
      !$acc                  wttend6,     &
      !$acc                  wstend1,     &
      !$acc                  wstend2,     &
      !$acc                  wstend3,     &
      !$acc                  wstend4,     &
      !$acc                  wstend5,     &
      !$acc                  wstend6,     &
      !$acc                  uwtend1,     &
      !$acc                  uwtend2,     &
      !$acc                  uwtend3,     &
      !$acc                  uwtend4,     &
      !$acc                  uwtend5,     &
      !$acc                  vwtend1,     &
      !$acc                  vwtend2,     &
      !$acc                  vwtend3,     &
      !$acc                  vwtend4,     &
      !$acc                  vwtend5,     &
      !$acc                  u2tend1,     &
      !$acc                  u2tend2,     &
      !$acc                  u2tend3,     &
      !$acc                  u2tend4,     &
      !$acc                  u2tend5,     &
      !$acc                  u2tend6,     &
      !$acc                  v2tend1,     &
      !$acc                  v2tend2,     &
      !$acc                  v2tend3,     &
      !$acc                  v2tend4,     &
      !$acc                  v2tend5,     &
      !$acc                  v2tend6,     &
      !$acc                  u2cliptend,  &
      !$acc                  v2cliptend,  &
      !$acc                  w2cliptend,  &
      !$acc                  u2tend,      &
      !$acc                  v2tend,      &
      !$acc                  w2tend,      &
      !$acc                  wttend,      &
      !$acc                  wstend,      &
      !$acc                  uwtend,      &
      !$acc                  vwtend,      &
      !$acc                  w3tend,      &
      !$acc                  uvtend,      &
      !$acc                  epstend,     &
      !$acc                  uttend,      &
      !$acc                  vttend,      &
      !$acc                  ustend,      &
      !$acc                  vstend,      &
      !$acc                  tstend,      &
      !$acc                  u2,          &
      !$acc                  v2,          &
      !$acc                  w2,          &
      !$acc                  wt,          &
      !$acc                  ws,          &
      !$acc                  uw,          &
      !$acc                  vw,          &
      !$acc                  w3,          &
      !$acc                  uv,          &
      !$acc                  t2,          &
      !$acc                  ut,          &
      !$acc                  vt,          &
      !$acc                  us,          &
      !$acc                  vs,          &
      !$acc                  ts,          &
      !$acc                  s2,          &
      !$acc                  eps)


   end subroutine ocn_turbulenceCreate

   subroutine ocn_turbulenceUpdateFields(domain)

   type(domain_type) :: &
         domain                    !< [in] MPAS type to describe domain

   integer :: &
         blockCount               ! counter for number of blocks

   type(block_type), pointer :: &
         block                    ! variables in current subblock

   type (mpas_pool_type), pointer :: &
         adcDiagnosticsPool, adcPrognosticsPool, adcTendPool, meshPool  ! mesh variables in MPAS pool structure

   real(kind=RKIND), dimension(:,:), pointer :: zeTmp, zmTmp, &
                       epsSPSTmp, lengthTmp, lenspsDTmp, lenspsUTmp, KhUTmp, KhDTmp, KmUTmp,  &
                       KmDTmp, wt_spsUTmp, wt_spsDTmp, ws_spsUTmp, ws_spsDTmp, uw2Tmp,  &
                       vw2Tmp, u2wTmp, v2wTmp, w2tTmp, w2sTmp, wtsTmp, uvwTmp, uwtTmp, vwtTMP,   &
                       uwsTmp, vwsTmp, ws2Tmp, wt2Tmp, areaFractionTmp, EntrainmentTmp,        &
                       DetrainmentTmp, tumdTmp, sumdTmp, wumdTmp, McTmp, lendnTmp, lenupTmp

   real(kind=RKIND), dimension(:,:), pointer :: w2tend1Tmp, w2tend2Tmp, w2tend3Tmp,  &
                       w2tend4Tmp, w2tend5Tmp, w2tend6Tmp, w2tend7Tmp, wttend1Tmp,   &
                       wttend2Tmp, wttend3Tmp, wttend4Tmp, wttend5Tmp, wttend6Tmp,   &
                       wttend7Tmp, w3tend1Tmp, w3tend2Tmp, w3tend3Tmp,               &
                       w3tend4Tmp, w3tend5Tmp, w3tend6Tmp, w3tend7Tmp, wstend1Tmp,   &
                       wstend2Tmp, wstend3Tmp, wstend4Tmp, wstend5Tmp, wstend6Tmp,   &
                       wstend7Tmp, uwtend1Tmp, uwtend2Tmp, uwtend3Tmp, uwtend4Tmp,   &
                       uwtend5Tmp, uwtend6Tmp, vwtend1Tmp, vwtend2Tmp, vwtend3Tmp,   &
                       vwtend4Tmp, vwtend5Tmp, vwtend6Tmp, u2tend1Tmp, u2tend2Tmp,   &
                       u2tend3Tmp, u2tend4Tmp, u2tend5Tmp, u2tend6Tmp, v2tend1Tmp,   &
                       v2tend2Tmp, v2tend3Tmp, v2tend4Tmp, v2tend5Tmp, v2tend6Tmp,   &
                       u2cliptendTmp, v2cliptendTmp, w2cliptendTmp

   real(kind=RKIND), dimension(:,:), pointer :: u2Tmp, v2Tmp, w2Tmp, t2Tmp, s2Tmp, &
                       uwTmp, vwTmp, wtTmp, wsTmp, w3Tmp, uvTmp, utTmp, vtTmp, usTmp,&
                       vsTmp, tsTmp, epsTmp, KspsUTmp, KspsDTmp, u2tendTmp, v2tendTmp, &
                       w2tendTmp, uwtendTmp, vwtendTmp, wttendTmp, wstendTmp, w3tendTmp, &
                       uvtendTmp, uttendTmp, vttendTmp, ustendTmp, vstendTmp, tstendTmp, &
                       epstendTmp, KspsUtendTmp, KspsDtendTmp

   integer, pointer :: nVertLevels, nVertLevelsP1
   integer :: nCells, i
   integer, dimension(:), pointer :: nCellsArray

   block => domain % blocklist
   do while (associated(block))
      !For runs on device, must update host before write out, could make this list smaller
      !at some point for production runs
      !should not ever need an update device no quantities get calculated on host

      !$acc update host(KspsD,    &
      !$acc                   nCells,   &
      !$acc                   nVertLevels,  &
      !$acc                   sigmat,       &
      !$acc                   grav,         &
      !$acc                   epsilon,      &
      !$acc                   Ko,           &
      !$acc                   gamma1,       &
      !$acc                   beta5,        &
      !$acc                   c1,           &
      !$acc                   c2,           &
      !$acc                   c8,           &
      !$acc                   c10,          &
      !$acc                   c11,          &
      !$acc                   alpha1,       &
      !$acc                   alpha2,       &
      !$acc                   alpha3,       &
      !$acc                   B1,           &
      !$acc                   Kt,           &
      !$acc                   c_1,          &
      !$acc                   c_2,          &
      !$acc                   c_mom,        &
      !$acc                   c_therm,      &
      !$acc                   c_mom_w3,     &
      !$acc                   c_Temp,       &
      !$acc                   c_Salt,       &
      !$acc                   c_Vel,        &
      !$acc                   c_vVel,       &
      !$acc                   c_epsilon,    &
      !$acc                   Cww_D,        &
      !$acc                   Cww_E,        &
      !$acc                   kappa_FL,     &
      !$acc                   kappa_VAR,    &
      !$acc                   kappa_w3,     &
      !$acc                   iterCount,    &
      !$acc                  KspsU,    &
      !$acc                  ze,       &
      !$acc                  zm,       &
      !$acc                  epsSPS,      &
      !$acc                  length,   &
      !$acc                  lenup,    &
      !$acc                  lendn,    &
      !$acc                  lenspsD,  &
      !$acc                  lenspsU,  &
      !$acc                  KhU,      &
      !$acc                  KhD,      &
      !$acc                  KmU,      &
      !$acc                  KmD,      &
      !$acc                  wt_spsU,  &
      !$acc                  wt_spsD,  &
      !$acc                  ws_spsU,  &
      !$acc                  ws_spsD,  &
      !$acc                  uw2,      &
      !$acc                  vw2,      &
      !$acc                  u2w,      &
      !$acc                  v2w,      &
      !$acc                  w2t,      &
      !$acc                  w2s,      &
      !$acc                  wts,      &
      !$acc                  uvw,      &
      !$acc                  uwt,      &
      !$acc                  vwt,      &
      !$acc                  uws,      &
      !$acc                  vws,      &
      !$acc                  ws2,      &
      !$acc                  wt2,      &
      !$acc                  areaFraction,    &
      !$acc                  Entrainment, &
      !$acc                  Detrainment, &
      !$acc                  tumd,        &
      !$acc                  sumd,        &
      !$acc                  wumd,        &
      !$acc                  w2tend1,     &
      !$acc                  w2tend2,     &
      !$acc                  w2tend3,     &
      !$acc                  w2tend4,     &
      !$acc                  w2tend5,     &
      !$acc                  w2tend6,     &
      !$acc                  w3tend1,     &
      !$acc                  w3tend2,     &
      !$acc                  w3tend3,     &
      !$acc                  w3tend4,     &
      !$acc                  w3tend5,     &
      !$acc                  w3tend6,     &
      !$acc                  wttend1,     &
      !$acc                  wttend2,     &
      !$acc                  wttend3,     &
      !$acc                  wttend4,     &
      !$acc                  wttend5,     &
      !$acc                  wttend6,     &
      !$acc                  wstend1,     &
      !$acc                  wstend2,     &
      !$acc                  wstend3,     &
      !$acc                  wstend4,     &
      !$acc                  wstend5,     &
      !$acc                  wstend6,     &
      !$acc                  uwtend1,     &
      !$acc                  uwtend2,     &
      !$acc                  uwtend3,     &
      !$acc                  uwtend4,     &
      !$acc                  uwtend5,     &
      !$acc                  vwtend1,     &
      !$acc                  vwtend2,     &
      !$acc                  vwtend3,     &
      !$acc                  vwtend4,     &
      !$acc                  vwtend5,     &
      !$acc                  u2tend1,     &
      !$acc                  u2tend2,     &
      !$acc                  u2tend3,     &
      !$acc                  u2tend4,     &
      !$acc                  u2tend5,     &
      !$acc                  u2tend6,     &
      !$acc                  v2tend1,     &
      !$acc                  v2tend2,     &
      !$acc                  v2tend3,     &
      !$acc                  v2tend4,     &
      !$acc                  v2tend5,     &
      !$acc                  v2tend6,     &
      !$acc                  u2cliptend,  &
      !$acc                  v2cliptend,  &
      !$acc                  w2cliptend,  &
      !$acc                  u2tend,      &
      !$acc                  v2tend,      &
      !$acc                  w2tend,      &
      !$acc                  wttend,      &
      !$acc                  wstend,      &
      !$acc                  uwtend,      &
      !$acc                  vwtend,      &
      !$acc                  w3tend,      &
      !$acc                  uvtend,      &
      !$acc                  epstend,     &
      !$acc                  uttend,      &
      !$acc                  vttend,      &
      !$acc                  ustend,      &
      !$acc                  vstend,      &
      !$acc                  tstend,      &
      !$acc                  u2,          &
      !$acc                  v2,          &
      !$acc                  w2,          &
      !$acc                  wt,          &
      !$acc                  ws,          &
      !$acc                  uw,          &
      !$acc                  vw,          &
      !$acc                  w3,          &
      !$acc                  uv,          &
      !$acc                  t2,          &
      !$acc                  ut,          &
      !$acc                  vt,          &
      !$acc                  us,          &
      !$acc                  vs,          &
      !$acc                  ts,          &
      !$acc                  s2,          &
      !$acc                  eps)

     call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)

      nCells = nCellsArray( size(nCellsArray) )

      call mpas_pool_get_subpool(block % structs, 'adcDiagnosticArrays', adcDiagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'adcTendArrays', adcTendPool)
      call mpas_pool_get_subpool(block % structs, 'adcPrognosticArrays', adcPrognosticsPool)
      call mpas_pool_get_array(adcDiagnosticsPool, 'ze', zeTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'zm', zmTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'epsSPS', epsSPSTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'length', lengthTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'lenup', lenupTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'lendn', lendnTmp)
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
      call mpas_pool_get_array(adcDiagnosticsPool, 'areaFraction', areaFractionTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'Entrainment', EntrainmentTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'Detrainment', DetrainmentTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'tumd', tumdTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'sumd', sumdTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'wumd', wumdTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'Mc', McTmp)

      McTmp = Mc
      zeTmp = ze
      zmTmp = zm
      epsSPSTmp = epsSPS
      lengthTmp = length
      lenupTmp = lenup
      lendnTmp = lendn
      lenspsDTmp = lenspsD
      lenspsUTmp = lenspsU
      KhUTmp = KhU
      KhDTmp = KhD
      KmUTmp = KmU
      KmDTmp = KmD
      wt_spsUTmp = wt_spsU
      wt_spsDTmp = wt_spsD
      ws_spsUTmp = ws_spsU
      ws_spsDTmp = ws_spsD
      uw2Tmp = uw2
      vw2Tmp = vw2
      u2wTmp = u2w
      v2wTmp = v2w
      w2tTmp = w2t
      w2sTmp = w2s
      wtsTmp = wts
      uvwTmp = uvw
      uwtTmp = uwt
      vwtTmp = vwt
      uwsTmp = uws
      vwsTmp = vws
      ws2Tmp = ws2
      wt2Tmp = wt2
      areaFractionTmp = areaFraction
      EntrainmentTmp = Entrainment
      DetrainmentTmp = Detrainment
      tumdTmp = tumd
      sumdTmp = sumd
      wumdTmp = wumd

      call mpas_pool_get_array(adcTendPool, 'w2tend1', w2tend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend2', w2tend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend3', w2tend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend4', w2tend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend5', w2tend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend6', w2tend6Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend7', w2tend7Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend1', w3tend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend2', w3tend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend3', w3tend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend4', w3tend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend5', w3tend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend6', w3tend6Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend7', w3tend7Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend1', wttend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend2', wttend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend3', wttend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend4', wttend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend5', wttend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend6', wttend6Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend7', wttend7Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend1', wstend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend2', wstend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend3', wstend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend4', wstend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend5', wstend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend6', wstend6Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend7', wstend7Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend1', uwtend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend2', uwtend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend3', uwtend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend4', uwtend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend5', uwtend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend6', uwtend6Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend1', vwtend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend2', vwtend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend3', vwtend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend4', vwtend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend5', vwtend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend6', vwtend6Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend1', u2tend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend2', u2tend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend3', u2tend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend4', u2tend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend5', u2tend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend6', u2tend6Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend1', v2tend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend2', v2tend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend3', v2tend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend4', v2tend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend5', v2tend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend6', v2tend6Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2cliptend', u2cliptendTmp)
      call mpas_pool_get_array(adcTendPool, 'v2cliptend', v2cliptendTmp)
      call mpas_pool_get_array(adcTendPool, 'w2cliptend', w2cliptendTmp)

      w2tend1Tmp = w2tend1
      w2tend2Tmp = w2tend2
      w2tend3Tmp = w2tend3
      w2tend4Tmp = w2tend4
      w2tend5Tmp = w2tend5
      w2tend6Tmp = w2tend6
      w2tend7Tmp = w2tend7
      w3tend1Tmp = w3tend1
      w3tend2Tmp = w3tend2
      w3tend3Tmp = w3tend3
      w3tend4Tmp = w3tend4
      w3tend5Tmp = w3tend5
      w3tend6Tmp = w3tend6
      w3tend7Tmp = w3tend7
      wttend1Tmp = wttend1
      wttend2Tmp = wttend2
      wttend3Tmp = wttend3
      wttend4Tmp = wttend4
      wttend5Tmp = wttend5
      wttend6Tmp = wttend6
      wttend7Tmp = wttend7
      wstend1Tmp = wstend1
      wstend2Tmp = wstend2
      wstend3Tmp = wstend3
      wstend4Tmp = wstend4
      wstend5Tmp = wstend5
      wstend6Tmp = wstend6
      wstend7Tmp = wstend7
      uwtend1Tmp = uwtend1
      uwtend2Tmp = uwtend2
      uwtend3Tmp = uwtend3
      uwtend4Tmp = uwtend4
      uwtend5Tmp = uwtend5
      uwtend6Tmp = uwtend6
      vwtend1Tmp = vwtend1
      vwtend2Tmp = vwtend2
      vwtend3Tmp = vwtend3
      vwtend4Tmp = vwtend4
      vwtend5Tmp = vwtend5
      vwtend6Tmp = vwtend6
      u2tend1Tmp = u2tend1
      u2tend2Tmp = u2tend2
      u2tend3Tmp = u2tend3
      u2tend4Tmp = u2tend4
      u2tend5Tmp = u2tend5
      u2tend6Tmp = u2tend6
      v2tend1Tmp = v2tend1
      v2tend2Tmp = v2tend2
      v2tend3Tmp = v2tend3
      v2tend4Tmp = v2tend4
      v2tend5Tmp = v2tend5
      v2tend6Tmp = v2tend6
      u2cliptendTmp = u2cliptend
      v2cliptendTmp = v2cliptend
      w2cliptendTmp = w2cliptend

      do i=1,3
        call mpas_pool_get_array(adcPrognosticsPool, 'KspsD', KspsDTmp,i)
        call mpas_pool_get_array(adcPrognosticsPool, 'KspsU', KspsUTmp,i)
        call mpas_pool_get_array(adcPrognosticsPool, 'u2', u2Tmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'v2', v2Tmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'w2', w2Tmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'wt', wtTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'ws', wsTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'uw', uwTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'vw', vwTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'w3', w3Tmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'uv', uvTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 't2', t2Tmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'ut', utTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'vt', vtTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'us', usTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'vs', vsTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'ts', tsTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 's2', s2Tmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'eps', epsTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'KspsDtend', KspsDtendTmp,i)
        call mpas_pool_get_array(adcPrognosticsPool, 'KspsUtend', KspsUtendTmp,i)
        call mpas_pool_get_array(adcPrognosticsPool, 'u2tend', u2tendTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'v2tend', v2tendTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'w2tend', w2tendTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'wttend', wttendTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'wstend', wstendTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'uwtend', uwtendTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'vwtend', vwtendTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'w3tend', w3tendTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'uvtend', uvtendTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'uttend', uttendTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'vttend', vttendTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'ustend', ustendTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'vstend', vstendTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'tstend', tstendTmp, i)
        call mpas_pool_get_array(adcPrognosticsPool, 'epstend', epstendTmp, i)

        epsTmp(:,:nCells) = eps(i,:,:nCells)
        u2Tmp(:,:nCells) = u2(i,:,:nCells)
        v2Tmp(:,:nCells) = v2(i,:,:nCells)
        w2Tmp(:,:nCells) = w2(i,:,:nCells)
        wtTmp(:,:nCells) = wt(i,:,:nCells)
        wsTmp(:,:nCells) = ws(i,:,:nCells)
        uwTmp(:,:nCells) = uw(i,:,:nCells)
        vwTmp(:,:nCells) = vw(i,:,:nCells)
        w3Tmp(:,:nCells) = w3(i,:,:nCells)
        uvTmp(:,:nCells) = uv(i,:,:nCells)
        t2Tmp(:,:nCells) = t2(i,:,:nCells)
        utTmp(:,:nCells) = ut(i,:,:nCells)
        vtTmp(:,:nCells) = vt(i,:,:nCells)
        usTmp(:,:nCells) = us(i,:,:nCells)
        vsTmp(:,:nCells) = vs(i,:,:nCells)
        tsTmp(:,:nCells) = ts(i,:,:nCells)
        s2Tmp(:,:nCells) = s2(i,:,:nCells)
        KspsUTmp(:,:nCells) = KspsU(i,:,:nCells)
        KspsDTmp(:,:nCells) = KspsD(i,:,:nCells)
        epstendTmp(:,:nCells) = epstend(i,:,:nCells)
        u2tendTmp(:,:nCells) = u2tend(i,:,:nCells)
        v2tendTmp(:,:nCells) = v2tend(i,:,:nCells)
        w2tendTmp(:,:nCells) = w2tend(i,:,:nCells)
        wttendTmp(:,:nCells) = wttend(i,:,:nCells)
        wstendTmp(:,:nCells) = wstend(i,:,:nCells)
        uwtendTmp(:,:nCells) = uwtend(i,:,:nCells)
        vwtendTmp(:,:nCells) = vwtend(i,:,:nCells)
        w3tendTmp(:,:nCells) = w3tend(i,:,:nCells)
        uvtendTmp(:,:nCells) = uvtend(i,:,:nCells)
        uttendTmp(:,:nCells) = uttend(i,:,:nCells)
        vttendTmp(:,:nCells) = vttend(i,:,:nCells)
        ustendTmp(:,:nCells) = ustend(i,:,:nCells)
        vstendTmp(:,:nCells) = vstend(i,:,:nCells)
        tstendTmp(:,:nCells) = tstend(i,:,:nCells)
        KspsUtendTmp(:,:nCells) = KspsUtend(i,:,:nCells)
        KspsDtendTmp(:,:nCells) = KspsDtend(i,:,:nCells)

      end do
      block => block % next
   end do
!   at end could do a copy to accelerator


   end subroutine ocn_turbulenceUpdateFields

   subroutine ocn_turbulenceDestroy(ierr)

     integer, intent(inout) :: ierr
!   deallocate stuff, set integers and constants to zero, nullify stuff on gpu

      sigmat = 0.0_RKIND
      epsilon = 0.0_RKIND
      Ko = 0.0_RKIND
      c_b = 0.0_RKIND
      c_b_tracer = 0.0_RKIND
      c11 = 0.0_RKIND
      alpha_0 = 0.0_RKIND
      alpha_1 = 0.0_RKIND
      alpha_2 = 0.0_RKIND
      alpha_tracer1 = 0.0_RKIND
      alpha_tracer2 = 0.0_RKIND
      alpha_st_0 = 0.0_RKIND
      alpha_st_1 = 0.0_RKIND
      alpha_st_2 = 0.0_RKIND
      alpha_st_tracer1 = 0.0_RKIND
      alpha_st_tracer2 = 0.0_RKIND
      B1 = 0.0_RKIND
      Kt = 0.0_RKIND
      c_mom = 0.0_RKIND
      c_therm = 0.0_RKIND
      c_mom_w3 = 0.0_RKIND
      c_slow = 0.0_RKIND
      c_slow_tracer = 0.0_RKIND
      slow_w_factor = 0.0_RKIND
      c_epsilon = 0.0_RKIND
      Cww_D = 0.0_RKIND
      Cww_E = 0.0_RKIND
      kappa_FL = 0.0_RKIND
      kappa_VAR = 0.0_RKIND
      kappa_w3 = 0.0_RKIND
      La_t = 0.0_RKIND
      stokes_depth = 0.0_RKIND
      stokes_angle = 0.0_RKIND

      !$acc exit data delete(KspsD,    &
      !$acc                   nCells,   &
      !$acc                   nVertLevels,  &
      !$acc                   sigmat,       &
      !$acc                   grav,         &
      !$acc                   epsilon,      &
      !$acc                   Ko,           &
      !$acc                   gamma1,       &
      !$acc                   beta5,        &
      !$acc                   c1,           &
      !$acc                   c2,           &
      !$acc                   c8,           &
      !$acc                   c10,          &
      !$acc                   c11,          &
      !$acc                   alpha1,       &
      !$acc                   alpha2,       &
      !$acc                   alpha3,       &
      !$acc                   B1,           &
      !$acc                   Kt,           &
      !$acc                   c_1,          &
      !$acc                   c_2,          &
      !$acc                   c_mom,        &
      !$acc                   c_therm,      &
      !$acc                   c_mom_w3,     &
      !$acc                   c_Temp,       &
      !$acc                   c_Salt,       &
      !$acc                   c_Vel,        &
      !$acc                   c_vVel,       &
      !$acc                   c_epsilon,    &
      !$acc                   Cww_D,        &
      !$acc                   Cww_E,        &
      !$acc                   kappa_FL,     &
      !$acc                   kappa_VAR,    &
      !$acc                   kappa_w3,     &
      !$acc                   iterCount,    &
    !$acc                  KspsU,    &
      !$acc                  ze,       &
      !$acc                  zm,       &
      !$acc                  epsSPS,      &
      !$acc                  length,   &
      !$acc                  lenup,    &
      !$acc                  lendn,    &
      !$acc                  lenspsD,  &
      !$acc                  lenspsU,  &
      !$acc                  KhU,      &
      !$acc                  KhD,      &
      !$acc                  KmU,      &
      !$acc                  KmD,      &
      !$acc                  wt_spsU,  &
      !$acc                  wt_spsD,  &
      !$acc                  ws_spsU,  &
      !$acc                  ws_spsD,  &
      !$acc                  uw2,      &
      !$acc                  vw2,      &
      !$acc                  u2w,      &
      !$acc                  v2w,      &
      !$acc                  w2t,      &
      !$acc                  w2s,      &
      !$acc                  wts,      &
      !$acc                  uvw,      &
      !$acc                  uwt,      &
      !$acc                  vwt,      &
      !$acc                  uws,      &
      !$acc                  vws,      &
      !$acc                  ws2,      &
      !$acc                  wt2,      &
      !$acc                  areaFraction,    &
      !$acc                  Entrainment, &
      !$acc                  Detrainment, &
      !$acc                  tumd,        &
      !$acc                  sumd,        &
      !$acc                  wumd,        &
      !$acc                  w2tend1,     &
      !$acc                  w2tend2,     &
      !$acc                  w2tend3,     &
      !$acc                  w2tend4,     &
      !$acc                  w2tend5,     &
      !$acc                  w2tend6,     &
      !$acc                  w3tend1,     &
      !$acc                  w3tend2,     &
      !$acc                  w3tend3,     &
      !$acc                  w3tend4,     &
      !$acc                  w3tend5,     &
      !$acc                  w3tend6,     &
      !$acc                  wttend1,     &
      !$acc                  wttend2,     &
      !$acc                  wttend3,     &
      !$acc                  wttend4,     &
      !$acc                  wttend5,     &
      !$acc                  wttend6,     &
      !$acc                  wstend1,     &
      !$acc                  wstend2,     &
      !$acc                  wstend3,     &
      !$acc                  wstend4,     &
      !$acc                  wstend5,     &
      !$acc                  wstend6,     &
      !$acc                  uwtend1,     &
      !$acc                  uwtend2,     &
      !$acc                  uwtend3,     &
      !$acc                  uwtend4,     &
      !$acc                  uwtend5,     &
      !$acc                  vwtend1,     &
      !$acc                  vwtend2,     &
      !$acc                  vwtend3,     &
      !$acc                  vwtend4,     &
      !$acc                  vwtend5,     &
      !$acc                  u2tend1,     &
      !$acc                  u2tend2,     &
      !$acc                  u2tend3,     &
      !$acc                  u2tend4,     &
      !$acc                  u2tend5,     &
      !$acc                  u2tend6,     &
      !$acc                  v2tend1,     &
      !$acc                  v2tend2,     &
      !$acc                  v2tend3,     &
      !$acc                  v2tend4,     &
      !$acc                  v2tend5,     &
      !$acc                  v2tend6,     &
      !$acc                  u2cliptend,  &
      !$acc                  v2cliptend,  &
      !$acc                  w2cliptend,  &
      !$acc                  u2tend,      &
      !$acc                  v2tend,      &
      !$acc                  w2tend,      &
      !$acc                  wttend,      &
      !$acc                  wstend,      &
      !$acc                  uwtend,      &
      !$acc                  vwtend,      &
      !$acc                  w3tend,      &
      !$acc                  uvtend,      &
      !$acc                  epstend,     &
      !$acc                  uttend,      &
      !$acc                  vttend,      &
      !$acc                  ustend,      &
      !$acc                  vstend,      &
      !$acc                  tstend,      &
      !$acc                  u2,          &
      !$acc                  v2,          &
      !$acc                  w2,          &
      !$acc                  wt,          &
      !$acc                  ws,          &
      !$acc                  uw,          &
      !$acc                  vw,          &
      !$acc                  w3,          &
      !$acc                  uv,          &
      !$acc                  t2,          &
      !$acc                  ut,          &
      !$acc                  vt,          &
      !$acc                  us,          &
      !$acc                  vs,          &
      !$acc                  ts,          &
      !$acc                  s2,          &
      !$acc                  eps)

      deallocate(KspsD,    &
                 KspsU,    &
                 ze,       &
                 zm,       &
                 epsSPS,      &
                 length,   &
                 lenup,    &
                 lendn,    &
                 lenspsD,  &
                 lenspsU,  &
                 KhU,      &
                 KhD,      &
                 KmU,      &
                 KmD,      &
                 wt_spsU,  &
                 wt_spsD,  &
                 ws_spsU,  &
                 ws_spsD,  &
                 uw2,      &
                 vw2,      &
                 u2w,      &
                 v2w,      &
                 w2t,      &
                 w2s,      &
                 wts,      &
                 uvw,      &
                 uwt,      &
                 vwt,      &
                 uws,      &
                 vws,      &
                 ws2,      &
                 wt2,      &
                 areaFraction,    &
                 Entrainment, &
                 Detrainment, &
                 tumd,        &
                 sumd,        &
                 wumd,        &
                 w2tend1,     &
                 w2tend2,     &
                 w2tend3,     &
                 w2tend4,     &
                 w2tend5,     &
                 w2tend6,     &
                 w2tend7,     &
                 w3tend1,     &
                 w3tend2,     &
                 w3tend3,     &
                 w3tend4,     &
                 w3tend5,     &
                 w3tend6,     &
                 w3tend7,     &
                 wttend1,     &
                 wttend2,     &
                 wttend3,     &
                 wttend4,     &
                 wttend5,     &
                 wttend6,     &
                 wttend7,     &
                 wstend1,     &
                 wstend2,     &
                 wstend3,     &
                 wstend4,     &
                 wstend5,     &
                 wstend6,     &
                 wstend7,     &
                 uwtend1,     &
                 uwtend2,     &
                 uwtend3,     &
                 uwtend4,     &
                 uwtend5,     &
                 uwtend6,     &
                 vwtend1,     &
                 vwtend2,     &
                 vwtend3,     &
                 vwtend4,     &
                 vwtend5,     &
                 vwtend6,     &
                 u2tend1,     &
                 u2tend2,     &
                 u2tend3,     &
                 u2tend4,     &
                 u2tend5,     &
                 u2tend6,     &
                 v2tend1,     &
                 v2tend2,     &
                 v2tend3,     &
                 v2tend4,     &
                 v2tend5,     &
                 v2tend6,     &
                 u2cliptend,  &
                 v2cliptend,  &
                 w2cliptend,  &
                 u2tend,      &
                 v2tend,      &
                 w2tend,      &
                 wttend,      &
                 wstend,      &
                 uwtend,      &
                 vwtend,      &
                 w3tend,      &
                 uvtend,      &
                 epstend,     &
                 uttend,      &
                 vttend,      &
                 ustend,      &
                 vstend,      &
                 tstend,      &
                 u2,          &
                 v2,          &
                 w2,          &
                 wt,          &
                 ws,          &
                 uw,          &
                 vw,          &
                 w3,          &
                 uv,          &
                 t2,          &
                 ut,          &
                 vt,          &
                 us,          &
                 vs,          &
                 ts,          &
                 s2,          &
                 eps,         &
                 Mc,          &
                 KspsDtend,   &
                 KspsUtend)

   end subroutine ocn_turbulenceDestroy

   end module ocn_turbulence
