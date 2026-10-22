module ocn_adc_mixing_constants

   use mpas_derived_types
   use mpas_pool_routines

   type :: adc_mixing_constants
     real :: grav,sigmat,Ko,gamma1,beta5,c1,c2,ce,alpha1,&
        alpha2,alpha3,c8,c10,c11,B1,Kt,cp,rho0,c_1,c_2,   &
        c_mom,c_therm,c_mom_w3,c_ps,c_pt,c_pv,kappa_FL,   &
        kappa_w3,kappa_VAR,Cww_E, Cww_D
   end type adc_mixing_constants

   real :: grav,sigmat,Ko,gamma1,beta5,c1,c2,ce,alpha1,&
     alpha2,alpha3,c8,c10,c11,B1,Kt,cp,rho0,c_1,c_2,   &
     c_mom,c_therm,c_pt,c_ps,c_pv,c_mom_w3,kappa_FL,   &
     kappa_w3,kappa_VAR, Cww_E, Cww_D

contains

subroutine init_mixing_constants
!Move all this to registry
      Ko = 1.66_RKIND**3.0_RKIND
      gamma1 = 1.0_RKIND/3.0_RKIND
      beta5 = 0.5_RKIND
      c1 = 2.88_RKIND
    c2 = 3.84_RKIND
    ce = 0.8_RKIND
    alpha1 = 0.984_RKIND
    alpha2 = 0.568_RKIND
    alpha3 = 0.29_RKIND
    c8 = 3.0_RKIND
    c10 = 5.0_RKIND
    c11 = 0.1_RKIND
    B1 = 16.6_RKIND
    Kt = 0.4_RKIND
  sigmat = 0.72_RKIND
  kappa_FL = 0.0_RKIND
  kappa_w3 = 0.0_RKIND
  kappa_VAR = 0.0_RKIND
  Cww_D = 1.5_RKIND
  Cww_E = 1.0_RKIND
end subroutine init_mixing_constants

end module mixing_constants
