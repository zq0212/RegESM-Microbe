module mod_clm_cndecompcascadebgc
#ifdef CN

#ifndef CENTURY_DECOMP

  !
  ! Module that sets the coeffiecients used in the decomposition cascade
  ! submodel.  This uses the BGC parameters as in CLMCN 4.0
  !
  use mod_intkinds
  use mod_realkinds
  use mod_runparams , only : dtsrf
  use mod_clm_varcon , only : tfrz , zsoi
  use mod_clm_varpar , only : nlevsoi, nlevgrnd, nlevdecomp, &
          ndecomp_cascade_transitions, ndecomp_pools, nsompools
#ifdef MICROBE
   use mod_clm_varpar   , only: i_met_lit, i_cel_lit, i_lig_lit, i_cwd, i_bacteria, i_fungi, i_dom, cn_dom, CUEmax
   use mod_clm_microbevarcon
#else
   use mod_clm_varpar   , only: i_met_lit, i_cel_lit, i_lig_lit, i_cwd
#endif
  use mod_clm_varctl , only : spinup_state
#if (defined LCH4) || (defined MICROBE)
   use mod_clm_varctl, only: anoxia
#endif
#ifdef LCH4
  use mod_clm_ch4varcon , only : mino2lim
#endif

  implicit none

  save

  private

  public:: init_decompcascade, decomp_rate_constants
!!! zhang add for Microbe
#if (defined VERTSOILC) || (defined MICROBE)
  ! (meters) e-folding depth for reduction in decomposition
  ! [set to large number for depth-independance]
  !real(rk8), public :: decomp_depth_efolding = 0.5_rk8
#endif
  ! separate q10 for frozen soil respiration rates.
  !  default to same as above zero rates
  real(rk8), public :: froz_q10 = 1.5_rk8
#if (defined LCH4) || (defined MICROBE)
  ! true ==> weight anoxia by inundated fraction
  logical,  public :: anoxia_wtsat = .false.
#endif
  integer(ik4), public :: nlev_soildecomp_standard   ! used here and in ch4Mod

  !! parameters for AD spinup
  ! multipliers for soil decomp during accelerated spinup
  real(rk8), public, parameter :: spinup_vector(nsompools) = &
          (/ 1.0_rk8, 1.0_rk8, 5.0_rk8, 70.0_rk8 /)

  contains

  subroutine init_decompcascade(begc, endc)
    !
    ! initialize rate constants and decomposition pathways for the BGC
    ! model originally implemented in CLM-CN
    !  written by C. Koven based on original CLM4 decomposition cascade
    ! by P. Thornton
    use mod_clm_type
    use mod_clm_varpar      , only: max_pft_per_col, mach_eps
    use mod_clm_varcon, only : istsoil, istcrop
    implicit none
    ! column level
    integer(ik4) :: begc, endc ! per-proc beginning and ending column indices

    !-- properties of each pathway along decomposition cascade
    character(len=8), pointer :: cascade_step_name(:) ! name of transition
    ! respired fraction in decomposition step (frac)
    real(rk8), pointer :: rf_decomp_cascade(:,:,:)
    ! which pool is C taken from for a given decomposition step
    integer(ik4),  pointer :: cascade_donor_pool(:)
    ! which pool is C added to for a given decomposition step
    integer(ik4),  pointer :: cascade_receiver_pool(:)
    ! what fraction of C leaving a given pool passes through a given transition
    ! (frac)
    real(rk8), pointer :: pathfrac_decomp_cascade(:,:,:)
    !-- properties of each decomposing pool
    ! TRUE => pool has fixed C:N ratio
    logical,  pointer :: floating_cn_ratio_decomp_pools(:)
    ! name of pool for restart files
    character(len=8), pointer :: decomp_pool_name_restart(:)
    ! name of pool for history files
    character(len=8), pointer :: decomp_pool_name_history(:)
    ! name of pool for netcdf long names
    character(len=20), pointer :: decomp_pool_name_long(:)
    ! name of pool for netcdf short names
    character(len=8), pointer :: decomp_pool_name_short(:)
    logical, pointer :: is_litter(:)     ! TRUE => pool is a litter pool
    logical, pointer :: is_soil(:)       ! TRUE => pool is a soil pool
    logical, pointer :: is_cwd(:)        ! TRUE => pool is a cwd pool
!!! zhangq add for Microbe
#ifdef MICROBE
   logical, pointer :: is_microbe(:)     ! TRUE => pool is a microbe pool
   integer(ik4) :: c, fc, g, l
   real(rk8), pointer :: wtcol(:)         ! pft weight relative to column (0-1)
   integer(ik4), allocatable :: pft_index(:)
   !integer(ik4) :: pft_index(:)
   integer(ik4), pointer :: pfti(:)         ! pft index array
   integer(ik4), pointer :: pftf(:)         ! pft index array
#endif
    ! c:n ratio for initialization of pools
    real(rk8), pointer :: initial_cn_ratio(:)
    ! initial concentration for seeding at spinup
    real(rk8), pointer :: initial_stock(:)
    logical, pointer :: is_metabolic(:)   ! TRUE => pool is metabolic material
    logical, pointer :: is_cellulose(:)   ! TRUE => pool is cellulose
    logical, pointer :: is_lignin(:)      ! TRUE => pool is lignin
    ! factor for AD spinup associated with each pool
    real(rk8), pointer :: spinup_factor(:)
    real(rk8):: rf_l1s1      !respiration fraction litter 1 -> SOM 1
    real(rk8):: rf_l2s2      !respiration fraction litter 2 -> SOM 2
    real(rk8):: rf_l3s3      !respiration fraction litter 3 -> SOM 3
    real(rk8):: rf_s1s2      !respiration fraction SOM 1 -> SOM 2
    real(rk8):: rf_s2s3      !respiration fraction SOM 2 -> SOM 3
    real(rk8):: rf_s3s4      !respiration fraction SOM 3 -> SOM 4
    real(rk8):: cwd_fcel
    real(rk8):: cwd_flig
    real(rk8) :: cn_s1
    real(rk8) :: cn_s2
    real(rk8) :: cn_s3
    real(rk8) :: cn_s4
!!! zhangq add for Microbe
#ifdef MICROBE
   real(rk8) :: rf_l1m      ! fraction carbon release to atmosphere for l1 to microbe
   real(rk8) :: rf_l2m      ! fraction of carbon release to atmosphere for l2 to microbe
   real(rk8) :: rf_l3m      ! fraction of carbon release to atmosphere for l3 to microbe
   real(rk8),allocatable :: rf_s1m(:)      ! fraction of carbon release to atmosphere for s1 to microbe
   real(rk8),allocatable :: rf_s2m(:)      ! fraction of carbon release to atmosphere for s2 to microbe
   real(rk8),allocatable :: rf_s3m(:)      ! fraction of carbon release to atmosphere for s3 to microbe
   real(rk8),allocatable :: rf_s4m(:)      ! fraction of carbon release to atmosphere for s4 to microbe
   real(rk8),allocatable :: l1m_fb(:)      ! fraction of carbon to bacteria from l1
   real(rk8),allocatable :: l1m_ff(:)      ! fraction of carbon to fungi from l1
   real(rk8),allocatable :: l2m_fb(:)      ! fraction of carbon to bacteria from l2
   real(rk8),allocatable :: l2m_ff(:)      ! fraction of carbon to fungi from l2
   real(rk8),allocatable :: l3m_fb(:)      ! fraction of carbon to bacteria from l3
   real(rk8),allocatable :: l3m_ff(:)      ! fraction of carbon to fungi from l3
   real(rk8),allocatable :: s1m_fb(:)      ! fraction of carbon to bacteria from s1
   real(rk8),allocatable :: s1m_ff(:)      ! fraction of carbon to fungi from s1
   real(rk8),allocatable :: s2m_fb(:)      ! fraction of carbon to bacteria from s2
   real(rk8),allocatable :: s2m_ff(:)      ! fraction of carbon to fungi from s2
   real(rk8),allocatable :: s3m_fb(:)      ! fraction of carbon to bacteria from s3
   real(rk8),allocatable :: s3m_ff(:)      ! fraction of carbon to fungi from s3
   real(rk8),allocatable :: s4m_fb(:)      ! fraction of carbon to bacteria from s4
   real(rk8),allocatable :: s4m_ff(:)      ! fraction of carbon to fungi from s4

   real(rk8) :: l1dom_f   ! fraction of carbon from l1 to DOM
   real(rk8) :: l2dom_f   ! fraction of carbon from l2 to DOM
   real(rk8) :: l3dom_f   ! fraction of carbon from l3 to DOM
   real(rk8) :: s1dom_f   ! fraction of carbon from s1 to DOM
   real(rk8) :: s2dom_f   ! fraction of carbon from s2 to DOM
   real(rk8) :: s3dom_f   ! fraction of carbon from s3 to DOM
   real(rk8) :: s4dom_f   ! fraction of carbon from s4 to DOM

   real(rk8) :: l1s1_f      ! fraction of carbon from l1 to S1
   real(rk8) :: l2s2_f      ! fraction of carbon from l2 to S2
   real(rk8) :: l3s3_f      ! fraction of carbon from l3 to S3
   real(rk8) :: s1s2_f    ! fraction of carbon from s1 to S2
   real(rk8) :: s2s3_f    ! fraction of carbon from s2 to S3
   real(rk8) :: s3s4_f    ! fraction of carbon from s3 to S4
   
   real(rk8),allocatable :: batm_f(:)       ! fraction of carbon from bacteria to atmosphere
   real(rk8),allocatable :: bdom_f(:)    ! fraction of carbon from bacteria to dissolved organic matter
   real(rk8),allocatable :: bs1_f(:)     ! fraction of carbon from bacteria to s1
   real(rk8),allocatable :: bs2_f(:)     ! fraction of carbon from bacteria to s2
   real(rk8),allocatable :: bs3_f(:)     ! fraction of carbon from bacteria to s3
   real(rk8),allocatable :: bs4_f(:)     ! fraction of carbon from bacteria to s4
   real(rk8),allocatable :: fatm_f(:)       ! fraction of carbon from fungi to atmosphere
   real(rk8),allocatable :: fdom_f(:)       ! fraction of carbon from fungi to dissolved organic matter
   real(rk8),allocatable :: fs1_f(:)     ! fraction of carbon from fungi to s1
   real(rk8),allocatable :: fs2_f(:)     ! fraction of carbon from fungi to s2
   real(rk8),allocatable :: fs3_f(:)     ! fraction of carbon from fungi to s3
   real(rk8),allocatable :: fs4_f(:)     ! fraction of carbon from fungi to s4

   real(rk8),allocatable :: domb_f(:)    ! fraction of carbon to bacteria from dissolved organic matter
   real(rk8),allocatable :: domf_f(:)       ! fraction of carbon to fungi from dissolved organic matter
   real(rk8),allocatable :: doms1_f(:)   ! fraction of carbon to bacteria from s1
   real(rk8),allocatable :: doms2_f(:)   ! fraction of carbon to fungi from s2
   real(rk8),allocatable :: doms3_f(:)   ! fraction of carbon to bacteria from s3
   real(rk8),allocatable :: doms4_f(:)   ! fraction of carbon to fungi from s4

   real(rk8),allocatable :: cn_fungi_in(:)  ! C:N ratio of fungi
   real(rk8),allocatable :: cn_bacteria_in(:)  ! C:N ratio of bacteria

   integer(ik4) , pointer :: cgridcell(:)       ! gridcell index of column

   real(rk8), pointer :: m_bdom_f(:)     ! bacterial C -> DOM
   real(rk8), pointer :: m_bs1_f(:)     ! bacterial C -> SOM1
   real(rk8), pointer :: m_bs2_f(:)     ! bacterial C -> SOM2
   real(rk8), pointer :: m_bs3_f(:)     ! bacterial C -> SOM3
   real(rk8), pointer :: m_fdom_f(:)     ! fungal C -> DOM
   real(rk8), pointer :: m_fs1_f(:)     ! fungal C -> SOM1
   real(rk8), pointer :: m_fs2_f(:)     ! fungal C -> SOM2
   real(rk8), pointer :: m_fs3_f(:)     ! fungal C -> SOM3
   real(rk8), pointer :: m_rf_s1m(:)     ! SOM1 -> microbes
   real(rk8), pointer :: m_rf_s2m(:)     ! SOM2 -> microbes
   real(rk8), pointer :: m_rf_s3m(:)     ! SOM3 -> microbes
   real(rk8), pointer :: m_rf_s4m(:)     ! SOM4 -> microbes
   real(rk8), pointer :: m_batm_f(:)     ! bacteria -> atmosphere
   real(rk8), pointer :: m_fatm_f(:)     ! fungi -> atmosphere
   real(rk8), pointer :: m_domb_f(:)     ! DOM -> bacteria
   real(rk8), pointer :: m_domf_f(:)     ! DOM -> fungi
   real(rk8), pointer :: m_doms1_f(:)     ! DOM -> SOM1
   real(rk8), pointer :: m_doms2_f(:)     ! DOM -> SOM2
   real(rk8), pointer :: m_doms3_f(:)     ! DOM -> SOM3
   real(rk8), pointer :: cn_bacteria(:)     ! C:N ratio of bacteria
   real(rk8), pointer :: cn_fungi(:)     ! C:N ratio of fungi

#endif

    integer(ik4) :: i_litr1
    integer(ik4) :: i_litr2
    integer(ik4) :: i_litr3
    integer(ik4) :: i_soil1
    integer(ik4) :: i_soil2
    integer(ik4) :: i_soil3
    integer(ik4) :: i_soil4
    integer(ik4) :: i_atm

!!! zhangq add for Microbe
#ifdef MICROBE
   integer(ik4) :: i_cwdl2
   integer(ik4) :: i_cwdl3
   integer(ik4) :: i_l1b
   integer(ik4) :: i_l1f
   integer(ik4) :: i_l1s1
   integer(ik4) :: i_l2b
   integer(ik4) :: i_l2f
   integer(ik4) :: i_l2s2
   integer(ik4) :: i_l3b
   integer(ik4) :: i_l3f
   integer(ik4) :: i_l3s3   
   integer(ik4) :: i_s1b
   integer(ik4) :: i_s1f
   integer(ik4) :: i_s1s2      
   integer(ik4) :: i_s2b
   integer(ik4) :: i_s2f
   integer(ik4) :: i_s2s3         
   integer(ik4) :: i_s3b
   integer(ik4) :: i_s3f
   integer(ik4) :: i_s3s4         
   integer(ik4) :: i_s4b
   integer(ik4) :: i_s4f        !22
   integer(ik4) :: i_bs1
   integer(ik4) :: i_fs1
   integer(ik4) :: i_bs2
   integer(ik4) :: i_fs2
   integer(ik4) :: i_bs3
   integer(ik4) :: i_fs3
   integer(ik4) :: i_bs4
   integer(ik4) :: i_fs4        !30
   integer(ik4) :: i_doms1
   integer(ik4) :: i_doms2
   integer(ik4) :: i_doms3
   integer(ik4) :: i_doms4
   integer(ik4) :: i_bdom
   integer(ik4) :: i_fdom
   integer(ik4) :: i_domb
   integer(ik4) :: i_domf       !
   integer(ik4) :: i_batm
   integer(ik4) :: i_fatm    !34
   integer(ik4) :: i_l1dom
   integer(ik4) :: i_l2dom
   integer(ik4) :: i_l3dom
   integer(ik4) :: i_s1dom      !38
   integer(ik4) :: i_s2dom
   integer(ik4) :: i_s3dom      !46
   integer(ik4) :: i_s4dom      !47
#else
   integer(ik4) :: i_l1s1
   integer(ik4) :: i_l2s2
   integer(ik4) :: i_l3s3
   integer(ik4) :: i_s1s2
   integer(ik4) :: i_s2s3
   integer(ik4) :: i_s3s4
   integer(ik4) :: i_s4atm
   integer(ik4) :: i_cwdl2
   integer(ik4) :: i_cwdl3
#endif

    cascade_step_name       => decomp_cascade_con%cascade_step_name
    rf_decomp_cascade       => clm3%g%l%c%cps%rf_decomp_cascade
    cascade_donor_pool      => decomp_cascade_con%cascade_donor_pool
    cascade_receiver_pool   => decomp_cascade_con%cascade_receiver_pool
    pathfrac_decomp_cascade => clm3%g%l%c%cps%pathfrac_decomp_cascade
    floating_cn_ratio_decomp_pools => &
            decomp_cascade_con%floating_cn_ratio_decomp_pools
    decomp_pool_name_restart => decomp_cascade_con%decomp_pool_name_restart
    decomp_pool_name_history => decomp_cascade_con%decomp_pool_name_history
    decomp_pool_name_long   => decomp_cascade_con%decomp_pool_name_long
    decomp_pool_name_short  => decomp_cascade_con%decomp_pool_name_short
    is_litter          => decomp_cascade_con%is_litter
    is_soil            => decomp_cascade_con%is_soil
    is_cwd             => decomp_cascade_con%is_cwd

!!! zhangq modified for Micorbe
#ifdef MICROBE
   allocate(rf_s1m(1:endc-begc+1))      ! fraction of carbon release toatmosphere for s1 to microbe
   allocate(rf_s2m(1:endc-begc+1))      ! fraction of carbon release toatmosphere for s2 to microbe
   allocate(rf_s3m(1:endc-begc+1))      ! fraction of carbon release toatmosphere for s3 to microbe
   allocate(rf_s4m(1:endc-begc+1))      ! fraction of carbon release toatmosphere for s4 to microbe
   allocate(l1m_fb(1:endc-begc+1))      ! fraction of carbon to bacteria from l1
   allocate(l1m_ff(1:endc-begc+1))      ! fraction of carbon to fungi from l1
   allocate(l2m_fb(1:endc-begc+1))      ! fraction of carbon to bacteria from l2
   allocate(l2m_ff(1:endc-begc+1))      ! fraction of carbon to fungi from l2
   allocate(l3m_fb(1:endc-begc+1))      ! fraction of carbon to bacteria from l3
   allocate(l3m_ff(1:endc-begc+1))      ! fraction of carbon to fungi from l3
   allocate(s1m_fb(1:endc-begc+1))      ! fraction of carbon to bacteria from s1
   allocate(s1m_ff(1:endc-begc+1))      ! fraction of carbon to fungi from s1
   allocate(s2m_fb(1:endc-begc+1))      ! fraction of carbon to bacteria from s2
   allocate(s2m_ff(1:endc-begc+1))      ! fraction of carbon to fungi from s2
   allocate(s3m_fb(1:endc-begc+1))      ! fraction of carbon to bacteria from s3
   allocate(s3m_ff(1:endc-begc+1))      ! fraction of carbon to fungi from s3
   allocate(s4m_fb(1:endc-begc+1))      ! fraction of carbon to bacteria from s4
   allocate(s4m_ff(1:endc-begc+1))      ! fraction of carbon to fungi from s4

   allocate(batm_f(1:endc-begc+1))       ! fraction of carbon from bacteria toatmosphere
   allocate(bdom_f(1:endc-begc+1))    ! fraction of carbon from bacteria todissolved organic matter
   allocate(bs1_f(1:endc-begc+1))     ! fraction of carbon from bacteria to s1
   allocate(bs2_f(1:endc-begc+1))     ! fraction of carbon from bacteria to s2
   allocate(bs3_f(1:endc-begc+1))     ! fraction of carbon from bacteria to s3
   allocate(bs4_f(1:endc-begc+1))     ! fraction of carbon from bacteria to s4
   allocate(fatm_f(1:endc-begc+1))       ! fraction of carbon from fungi toatmosphere
   allocate(fdom_f(1:endc-begc+1))       ! fraction of carbon from fungi todissolved organic matter
   allocate(fs1_f(1:endc-begc+1))     ! fraction of carbon from fungi to s1
   allocate(fs2_f(1:endc-begc+1))     ! fraction of carbon from fungi to s2
   allocate(fs3_f(1:endc-begc+1))     ! fraction of carbon from fungi to s3
   allocate(fs4_f(1:endc-begc+1))     ! fraction of carbon from fungi to s4

   allocate(domb_f(1:endc-begc+1))    ! fraction of carbon to bacteria fromdissolved organic matter
   allocate(domf_f(1:endc-begc+1))       ! fraction of carbon to fungi fromdissolved organic matter
   allocate(doms1_f(1:endc-begc+1))   ! fraction of carbon to bacteria from s1
   allocate(doms2_f(1:endc-begc+1))   ! fraction of carbon to fungi from s2
   allocate(doms3_f(1:endc-begc+1))   ! fraction of carbon to bacteria from s3
   allocate(doms4_f(1:endc-begc+1))   ! fraction of carbon to fungi from s4

   allocate(cn_fungi_in(1:endc-begc+1))  ! C:N ratio of fungi
   allocate(cn_bacteria_in(1:endc-begc+1))  ! C:N ratio of bacteria
#endif

!!! zhangq add for Microbe
#ifdef MICROBE
    is_microbe    => decomp_cascade_con%is_microbe

    m_bdom_f      => pftcon%m_bdom_f
    m_bs1_f       => pftcon%m_bs1_f
    m_bs2_f       => pftcon%m_bs2_f
    m_bs3_f       => pftcon%m_bs3_f
    m_fdom_f      => pftcon%m_fdom_f
    m_fs1_f       => pftcon%m_fs1_f
    m_fs2_f       => pftcon%m_fs2_f
    m_fs3_f       => pftcon%m_fs3_f
    m_rf_s1m      => pftcon%m_rf_s1m
    m_rf_s2m      => pftcon%m_rf_s2m
    m_rf_s3m      => pftcon%m_rf_s3m
    m_rf_s4m      => pftcon%m_rf_s4m
    m_batm_f      => pftcon%m_batm_f
    m_fatm_f      => pftcon%m_fatm_f
    m_domb_f      => pftcon%m_domb_f
    m_domf_f      => pftcon%m_domf_f
    m_doms1_f     => pftcon%m_doms1_f
    m_doms2_f     => pftcon%m_doms2_f
    m_doms3_f     => pftcon%m_doms3_f
    cn_bacteria   => pftcon%cn_bacteria
    cn_fungi      => pftcon%cn_fungi

    cgridcell     => clm3%g%l%c%gridcell

    pfti          => clm3%g%l%c%pfti
    pftf          => clm3%g%l%c%pftf
    wtcol         => clm3%g%l%c%p%wtcol

!!! zhangq modified for Microbe
    allocate(pft_index(1:endc-begc+1))

    do c = begc,endc

    if (wtcol(clm3%g%l%c%pfti(c)) .ge. 0.96_rk8 .or. SUM(wtcol((clm3%g%l%c%pfti(c) + 1):clm3%g%l%c%pftf(c))) .lt. 0.001_rk8) then
       !pft_index = [pft_index, 0]
       pft_index(c-begc+1) = 0
    else
       !pft_index = [pft_index, MAXLOC(wtcol((clm3%g%l%c%pfti(c) + 1):clm3%g%l%c%pftf(c)), DIM=1)]
       pft_index(c-begc+1) = MAXLOC(wtcol((clm3%g%l%c%pfti(c) + 1):clm3%g%l%c%pftf(c)),DIM=1)
    end if

     !write (*, *) "col%pfti(c) 1", clm3%g%l%c%pfti(c)

     !write (*, *) "SUM(wtcol((col%pfti(c) + 1):col%pftf(c)))", SUM(wtcol((clm3%g%l%c%pfti(c) + 1):clm3%g%l%c%pftf(c)))

    end do

     !write (*, *) "pft_index 1", pft_index


    cn_bacteria_in = cn_bacteria(pft_index(:))
    cn_fungi_in = cn_fungi(pft_index(:))

     !write (*, *) "cn_bacteria_in 1", cn_bacteria_in
     !write (*, *) "cn_fungi_in 1", cn_fungi_in

#endif

    initial_cn_ratio   => decomp_cascade_con%initial_cn_ratio
    initial_stock      => decomp_cascade_con%initial_stock
    is_metabolic       => decomp_cascade_con%is_metabolic
    is_cellulose       => decomp_cascade_con%is_cellulose
    is_lignin          => decomp_cascade_con%is_lignin
    spinup_factor      => decomp_cascade_con%spinup_factor

    !------- time-constant coefficients ---------- !
    ! set soil organic matter compartment C:N ratios (from Biome-BGC v4.2.0)
    cn_s1 = 12.0_rk8
    cn_s2 = 12.0_rk8
    cn_s3 = 10.0_rk8
    cn_s4 = 10.0_rk8

    ! set respiration fractions for fluxes between compartments
    ! (from Biome-BGC v4.2.0)
    rf_l1s1 = 0.39_rk8
    rf_l2s2 = 0.55_rk8
    rf_l3s3 = 0.29_rk8
    rf_s1s2 = 0.28_rk8
    rf_s2s3 = 0.46_rk8
    rf_s3s4 = 0.55_rk8

    ! set the cellulose and lignin fractions for coarse woody debris
    cwd_fcel = 0.76_rk8
    cwd_flig = 0.24_rk8

    ! list of pools and their attributes

!!! zhangq add for Microbe
#ifdef MICROBE
   i_litr1 = i_met_lit
   floating_cn_ratio_decomp_pools(i_litr1) = .true.
   decomp_pool_name_restart(i_litr1) = 'litr1'
   decomp_pool_name_history(i_litr1) = 'LITR1'
   decomp_pool_name_long(i_litr1) = 'litter 1'
   decomp_pool_name_short(i_litr1) = 'L1'
   is_litter(i_litr1) = .true.
   is_soil(i_litr1) = .false.
   is_cwd(i_litr1) = .false.
   is_microbe(i_litr1) = .false.
   initial_cn_ratio(i_litr1) = 90._rk8
   initial_stock(i_litr1) = 0._rk8
   is_metabolic(i_litr1) = .true.
   is_cellulose(i_litr1) = .false.
   is_lignin(i_litr1) = .false.
   
   i_litr2 = i_cel_lit
   floating_cn_ratio_decomp_pools(i_litr2) = .true.
   decomp_pool_name_restart(i_litr2) = 'litr2'
   decomp_pool_name_history(i_litr2) = 'LITR2'
   decomp_pool_name_long(i_litr2) = 'litter 2'
   decomp_pool_name_short(i_litr2) = 'L2'
   is_litter(i_litr2) = .true.
   is_soil(i_litr2) = .false.
   is_cwd(i_litr2) = .false.
   is_microbe(i_litr2) = .false.
   initial_cn_ratio(i_litr2) = 90._rk8
   initial_stock(i_litr2) = 0._rk8
   is_metabolic(i_litr2) = .false.
   is_cellulose(i_litr2) = .true.
   is_lignin(i_litr2) = .false.

   i_litr3 = i_lig_lit
   floating_cn_ratio_decomp_pools(i_litr3) = .true.
   decomp_pool_name_restart(i_litr3) = 'litr3'
   decomp_pool_name_history(i_litr3) = 'LITR3'
   decomp_pool_name_long(i_litr3) = 'litter 3'
   decomp_pool_name_short(i_litr3) = 'L3'
   is_litter(i_litr3) = .true.
   is_soil(i_litr3) = .false.
   is_cwd(i_litr3) = .false.
   is_microbe(i_litr3) = .false.
   initial_cn_ratio(i_litr3) = 90._rk8
   initial_stock(i_litr3) = 0._rk8
   is_metabolic(i_litr3) = .false.
   is_cellulose(i_litr3) = .false.
   is_lignin(i_litr3) = .true.
   
! i_cwd
   floating_cn_ratio_decomp_pools(i_cwd) = .true.
   decomp_pool_name_restart(i_cwd) = 'cwd'
   decomp_pool_name_history(i_cwd) = 'CWD'
   decomp_pool_name_long(i_cwd) = 'coarse woody debris'
   decomp_pool_name_short(i_cwd) = 'CWD'
   is_litter(i_cwd) = .false.
   is_soil(i_cwd) = .false.
   is_cwd(i_cwd) = .true.
   is_microbe(i_cwd) = .true.
   initial_cn_ratio(i_cwd) = 500._rk8
   initial_stock(i_cwd) = 0._rk8
   is_metabolic(i_cwd) = .false.
   is_cellulose(i_cwd) = .false.
   is_lignin(i_cwd) = .false.

   i_soil1 = 5
   floating_cn_ratio_decomp_pools(i_soil1) = .false.
   decomp_pool_name_restart(i_soil1) = 'soil1'
   decomp_pool_name_history(i_soil1) = 'SOIL1'
   decomp_pool_name_long(i_soil1) = 'soil 1'
   decomp_pool_name_short(i_soil1) = 'S1'
   is_litter(i_soil1) = .false.
   is_soil(i_soil1) = .true.
   is_cwd(i_soil1) = .false.
   is_microbe(i_soil1) = .false.
   initial_cn_ratio(i_soil1) = cn_s1
   initial_stock(i_soil1) = 0._rk8
   is_metabolic(i_soil1) = .false.
   is_cellulose(i_soil1) = .false.
   is_lignin(i_soil1) = .false.

   i_soil2 = 6
   floating_cn_ratio_decomp_pools(i_soil2) = .false.
   decomp_pool_name_restart(i_soil2) = 'soil2'
   decomp_pool_name_history(i_soil2) = 'SOIL2'
   decomp_pool_name_long(i_soil2) = 'soil 2'
   decomp_pool_name_short(i_soil2) = 'S2'
   is_litter(i_soil2) = .false.
   is_soil(i_soil2) = .true.
   is_cwd(i_soil2) = .false.
   is_microbe(i_soil2) = .false.
   initial_cn_ratio(i_soil2) = cn_s2
   initial_stock(i_soil2) = 0._rk8
   is_metabolic(i_soil2) = .false.
   is_cellulose(i_soil2) = .false.
   is_lignin(i_soil2) = .false.

   i_soil3 = 7
   floating_cn_ratio_decomp_pools(i_soil3) = .false.
   decomp_pool_name_restart(i_soil3) = 'soil3'
   decomp_pool_name_history(i_soil3) = 'SOIL3'
   decomp_pool_name_long(i_soil3) = 'soil 3'
   decomp_pool_name_short(i_soil3) = 'S3'
   is_litter(i_soil3) = .false.
   is_soil(i_soil3) = .true.
   is_cwd(i_soil3) = .false.
   is_microbe(i_soil3) = .false.
   initial_cn_ratio(i_soil3) = cn_s3
   initial_stock(i_soil3) = 0._rk8
   is_metabolic(i_soil3) = .false.
   is_cellulose(i_soil3) = .false.
   is_lignin(i_soil3) = .false.

   i_soil4 = 8
   floating_cn_ratio_decomp_pools(i_soil4) = .false.
   decomp_pool_name_restart(i_soil4) = 'soil4'
   decomp_pool_name_history(i_soil4) = 'SOIL4'
   decomp_pool_name_long(i_soil4) = 'soil 4'
   decomp_pool_name_short(i_soil4) = 'S4'
   is_litter(i_soil4) = .false.
   is_soil(i_soil4) = .true.
   is_cwd(i_soil4) = .false.
   is_microbe(i_soil4) = .false.
   initial_cn_ratio(i_soil4) = cn_s4
   initial_stock(i_soil4) = 10._rk8
   is_metabolic(i_soil4) = .false.
   is_cellulose(i_soil4) = .false.
   is_lignin(i_soil4) = .false.

   i_atm = 0  !! for terminal pools (i.e. 100% respiration)
   floating_cn_ratio_decomp_pools(i_atm) = .false.
   decomp_pool_name_restart(i_atm) = 'atmosphere'
   decomp_pool_name_history(i_atm) = 'atmosphere'
   decomp_pool_name_long(i_atm) = 'atmosphere'
   decomp_pool_name_short(i_atm) = 'ATM'
   is_litter(i_atm) = .true.
   is_soil(i_atm) = .false.
   is_cwd(i_atm) = .false.
   is_microbe(i_atm) = .false.
   initial_cn_ratio(i_atm) = 0._rk8
   initial_stock(i_atm) = 0._rk8
   is_metabolic(i_atm) = .false.
   is_cellulose(i_atm) = .false.
   is_lignin(i_atm) = .false.

!   i_bacteria = 9  
   floating_cn_ratio_decomp_pools(i_bacteria) = .false.
   decomp_pool_name_restart(i_bacteria) = 'bacteria'
   decomp_pool_name_history(i_bacteria) = 'BACTERIA'
   decomp_pool_name_long(i_bacteria) = 'bacteria'
   decomp_pool_name_short(i_bacteria) = 'B'
   is_litter(i_bacteria) = .false.
   is_soil(i_bacteria) = .true.
   is_cwd(i_bacteria) = .false.
   is_microbe(i_bacteria) = .true.
   initial_cn_ratio(i_bacteria) = 5
   initial_stock(i_bacteria) = 1.e-5
   is_metabolic(i_bacteria) = .false.
   is_cellulose(i_bacteria) = .false.
   is_lignin(i_bacteria) = .false.

!   i_fungi = 10  
   floating_cn_ratio_decomp_pools(i_fungi) = .false.
   decomp_pool_name_restart(i_fungi) = 'fungi'
   decomp_pool_name_history(i_fungi) = 'FUNGI'
   decomp_pool_name_long(i_fungi) = 'fungi'
   decomp_pool_name_short(i_fungi) = 'F'
   is_litter(i_fungi) = .false.
   is_soil(i_fungi) = .true.
   is_cwd(i_fungi) = .false.
   is_microbe(i_fungi) = .true.
   initial_cn_ratio(i_fungi) = 15
   initial_stock(i_fungi) = 1.e-5
   is_metabolic(i_fungi) = .false.
   is_cellulose(i_fungi) = .false.
   is_lignin(i_fungi) = .false.

!   i_dom = 11
   floating_cn_ratio_decomp_pools(i_dom) = .false.
   decomp_pool_name_restart(i_dom) = 'dom'
   decomp_pool_name_history(i_dom) = 'DOM'
   decomp_pool_name_long(i_dom) = 'dissolved organic matter'
   decomp_pool_name_short(i_dom) = 'DOM'
   is_litter(i_dom) = .false.
   is_soil(i_dom) = .true.
   is_cwd(i_dom) = .false.
   is_microbe(i_dom) = .false.
   initial_cn_ratio(i_dom) = cn_dom
   initial_stock(i_dom) = 0._rk8
   is_metabolic(i_dom) = .false.
   is_cellulose(i_dom) = .false.
   is_lignin(i_dom) = .false.

#else
   i_litr1 = i_met_lit
   floating_cn_ratio_decomp_pools(i_litr1) = .true.
   decomp_pool_name_restart(i_litr1) = 'litr1'
   decomp_pool_name_history(i_litr1) = 'LITR1'
   decomp_pool_name_long(i_litr1) = 'litter 1'
   decomp_pool_name_short(i_litr1) = 'L1'
   is_litter(i_litr1) = .true.
   is_soil(i_litr1) = .false.
   is_cwd(i_litr1) = .false.
   initial_cn_ratio(i_litr1) = 90._rk8
   initial_stock(i_litr1) = 0._rk8
   is_metabolic(i_litr1) = .true.
   is_cellulose(i_litr1) = .false.
   is_lignin(i_litr1) = .false.
   
   i_litr2 = i_cel_lit
   floating_cn_ratio_decomp_pools(i_litr2) = .true.
   decomp_pool_name_restart(i_litr2) = 'litr2'
   decomp_pool_name_history(i_litr2) = 'LITR2'
   decomp_pool_name_long(i_litr2) = 'litter 2'
   decomp_pool_name_short(i_litr2) = 'L2'
   is_litter(i_litr2) = .true.
   is_soil(i_litr2) = .false.
   is_cwd(i_litr2) = .false.
   initial_cn_ratio(i_litr2) = 90._rk8
   initial_stock(i_litr2) = 0._rk8
   is_metabolic(i_litr2) = .false.
   is_cellulose(i_litr2) = .true.
   is_lignin(i_litr2) = .false.

   i_litr3 = i_lig_lit
   floating_cn_ratio_decomp_pools(i_litr3) = .true.
   decomp_pool_name_restart(i_litr3) = 'litr3'
   decomp_pool_name_history(i_litr3) = 'LITR3'
   decomp_pool_name_long(i_litr3) = 'litter 3'
   decomp_pool_name_short(i_litr3) = 'L3'
   is_litter(i_litr3) = .true.
   is_soil(i_litr3) = .false.
   is_cwd(i_litr3) = .false.
   initial_cn_ratio(i_litr3) = 90._rk8
   initial_stock(i_litr3) = 0._rk8
   is_metabolic(i_litr3) = .false.
   is_cellulose(i_litr3) = .false.
   is_lignin(i_litr3) = .true.

   floating_cn_ratio_decomp_pools(i_cwd) = .true.
   decomp_pool_name_restart(i_cwd) = 'cwd'
   decomp_pool_name_history(i_cwd) = 'CWD'
   decomp_pool_name_long(i_cwd) = 'coarse woody debris'
   decomp_pool_name_short(i_cwd) = 'CWD'
   is_litter(i_cwd) = .false.
   is_soil(i_cwd) = .false.
   is_cwd(i_cwd) = .true.
   initial_cn_ratio(i_cwd) = 500._rk8
   initial_stock(i_cwd) = 0._rk8
   is_metabolic(i_cwd) = .false.
   is_cellulose(i_cwd) = .false.
   is_lignin(i_cwd) = .false.

   i_soil1 = 5
   floating_cn_ratio_decomp_pools(i_soil1) = .false.
   decomp_pool_name_restart(i_soil1) = 'soil1'
   decomp_pool_name_history(i_soil1) = 'SOIL1'
   decomp_pool_name_long(i_soil1) = 'soil 1'
   decomp_pool_name_short(i_soil1) = 'S1'
   is_litter(i_soil1) = .false.
   is_soil(i_soil1) = .true.
   is_cwd(i_soil1) = .false.
   initial_cn_ratio(i_soil1) = cn_s1
   initial_stock(i_soil1) = 0._rk8
   is_metabolic(i_soil1) = .false.
   is_cellulose(i_soil1) = .false.
   is_lignin(i_soil1) = .false.

   i_soil2 = 6
   floating_cn_ratio_decomp_pools(i_soil2) = .false.
   decomp_pool_name_restart(i_soil2) = 'soil2'
   decomp_pool_name_history(i_soil2) = 'SOIL2'
   decomp_pool_name_long(i_soil2) = 'soil 2'
   decomp_pool_name_short(i_soil2) = 'S2'
   is_litter(i_soil2) = .false.
   is_soil(i_soil2) = .true.
   is_cwd(i_soil2) = .false.
   initial_cn_ratio(i_soil2) = cn_s2
   initial_stock(i_soil2) = 0._rk8
   is_metabolic(i_soil2) = .false.
   is_cellulose(i_soil2) = .false.
   is_lignin(i_soil2) = .false.

   i_soil3 = 7
   floating_cn_ratio_decomp_pools(i_soil3) = .false.
   decomp_pool_name_restart(i_soil3) = 'soil3'
   decomp_pool_name_history(i_soil3) = 'SOIL3'
   decomp_pool_name_long(i_soil3) = 'soil 3'
   decomp_pool_name_short(i_soil3) = 'S3'
   is_litter(i_soil3) = .false.
   is_soil(i_soil3) = .true.
   is_cwd(i_soil3) = .false.
   initial_cn_ratio(i_soil3) = cn_s3
   initial_stock(i_soil3) = 0._rk8
   is_metabolic(i_soil3) = .false.
   is_cellulose(i_soil3) = .false.
   is_lignin(i_soil3) = .false.

   i_soil4 = 8
   floating_cn_ratio_decomp_pools(i_soil4) = .false.
   decomp_pool_name_restart(i_soil4) = 'soil4'
   decomp_pool_name_history(i_soil4) = 'SOIL4'
   decomp_pool_name_long(i_soil4) = 'soil 4'
   decomp_pool_name_short(i_soil4) = 'S4'
   is_litter(i_soil4) = .false.
   is_soil(i_soil4) = .true.
   is_cwd(i_soil4) = .false.
   initial_cn_ratio(i_soil4) = cn_s4
   initial_stock(i_soil4) = 10._rk8
   is_metabolic(i_soil4) = .false.
   is_cellulose(i_soil4) = .false.
   is_lignin(i_soil4) = .false.

   i_atm = 0  !! for terminal pools (i.e. 100% respiration)
   floating_cn_ratio_decomp_pools(i_atm) = .false.
   decomp_pool_name_restart(i_atm) = 'atmosphere'
   decomp_pool_name_history(i_atm) = 'atmosphere'
   decomp_pool_name_long(i_atm) = 'atmosphere'
   decomp_pool_name_short(i_atm) = ''
   is_litter(i_atm) = .true.
   is_soil(i_atm) = .false.
   is_cwd(i_atm) = .false.
   initial_cn_ratio(i_atm) = 0._rk8
   initial_stock(i_atm) = 0._rk8
   is_metabolic(i_atm) = .false.
   is_cellulose(i_atm) = .false.
   is_lignin(i_atm) = .false.
#endif

    spinup_factor(i_litr1) = 1._rk8
    spinup_factor(i_litr2) = 1._rk8
    spinup_factor(i_litr3) = 1._rk8
    spinup_factor(i_cwd) = 1._rk8
    spinup_factor(i_soil1) = spinup_vector(1)
    spinup_factor(i_soil2) = spinup_vector(2)
    spinup_factor(i_soil3) = spinup_vector(3)
    spinup_factor(i_soil4) = spinup_vector(4)

!!! zhangq add for Microbe
#ifdef MICROBE
    spinup_factor(i_bacteria) = 1._rk8
    spinup_factor(i_fungi) = 1._rk8
    spinup_factor(i_dom) = 1._rk8
#endif

!!! zhangq add for Microbe
#ifdef MICROBE
! the following parameters will be changed upon finalizing
! keep this as backup

   !~ rf_l1m = 0.39
   !~ rf_l2m = 0.55
   !~ rf_l3m = 0.29
   !~ rf_s1m = 0.28
   !~ rf_s2m = 0.46
   !~ rf_s3m = 0.55
   !~ rf_s4m = 0.55
   ! Sinsabaugh et al., 2013, Ecology Letters, 16, 930-939
   rf_l1m = CUEmax * min(1._rk8, 8.0 / (initial_cn_ratio(i_litr1) * CUEmax))
   rf_l2m = CUEmax * min(1._rk8, 8.0 / (initial_cn_ratio(i_litr2) * CUEmax))
   rf_l3m = CUEmax * min(1._rk8, 8.0 / (initial_cn_ratio(i_litr3) * CUEmax))
   !rf_s1m = 0.28
   !rf_s2m = 0.46
   !rf_s3m = 0.55
   !rf_s4m = 0.75
   
   rf_s1m = m_rf_s1m(pft_index(:))
   rf_s2m = m_rf_s2m(pft_index(:))
   rf_s3m = m_rf_s3m(pft_index(:))
   rf_s4m = m_rf_s4m(pft_index(:))

   l1m_ff = ((cn_fungi_in / initial_cn_ratio(i_litr1))**0.6 + (1 / 30)**0.6)/ &
     ((cn_bacteria_in / initial_cn_ratio(i_litr1))**0.6 + (cn_fungi_in / &
     initial_cn_ratio(i_litr1))**0.6 + (1 / 30)**0.6)
   l1m_fb = 1.0 - l1m_ff 
   l2m_ff = ((cn_fungi_in / initial_cn_ratio(i_litr2))**0.6 + (1 / 30)**0.6)/ &
     ((cn_bacteria_in / initial_cn_ratio(i_litr2))**0.6 + (cn_fungi_in / &
     initial_cn_ratio(i_litr2))**0.6 + (1 / 30)**0.6)
   l2m_fb = 1.0 - l2m_ff
   l3m_ff = ((cn_fungi_in / initial_cn_ratio(i_litr3))**0.6 + (1 / 30)**0.6)/ &
     ((cn_bacteria_in / initial_cn_ratio(i_litr3))**0.6 + (cn_fungi_in / &
     initial_cn_ratio(i_litr3))**0.6 + (1 / 30)**0.6)
   l3m_fb = 1.0 - l3m_ff
   s1m_ff = ((cn_fungi_in / initial_cn_ratio(i_soil1))**0.6 + (1 / 30)**0.6)/ &
     ((cn_bacteria_in / initial_cn_ratio(i_soil1))**0.6 + (cn_fungi_in / &
     initial_cn_ratio(i_soil1))**0.6 + (1 / 30)**0.6)
   s1m_fb = 1.0 - s1m_ff
   s2m_ff = ((cn_fungi_in / initial_cn_ratio(i_soil2))**0.6 + (1 / 30)**0.6)/ &
     ((cn_bacteria_in / initial_cn_ratio(i_soil2))**0.6 + (cn_fungi_in / &
     initial_cn_ratio(i_soil2))**0.6 + (1 / 30)**0.6)
   s2m_fb = 1.0 - s2m_ff
   s3m_ff = ((cn_fungi_in / initial_cn_ratio(i_soil3))**0.6 + (1 / 30)**0.6)/ &
     ((cn_bacteria_in / initial_cn_ratio(i_soil3))**0.6 + (cn_fungi_in / &
     initial_cn_ratio(i_soil3))**0.6 + (1 / 30)**0.6)
   s3m_fb = 1.0 - s3m_ff
   s4m_ff = ((cn_fungi_in / initial_cn_ratio(i_soil4))**0.6 + (1 / 30)**0.6)/ &
     ((cn_bacteria_in / initial_cn_ratio(i_soil4))**0.6 + (cn_fungi_in / &
     initial_cn_ratio(i_soil4))**0.6 + (1 / 30)**0.6)
   s4m_fb = 1.0 - s4m_ff
   
!   l1dom_f  = 0.2
!   l2dom_f  = 0.16
!   l3dom_f  = 0.12
!   s1dom_f   = 0.36
!   s2dom_f   = 0.28
!   s3dom_f   = 0.20
!   s4dom_f   = 0.12

! Jan 27, 2022, Xiaofeng replaced the above code with the below ones to add the
! C flow to SOM, supruss microbes

   l1dom_f  = 0.1   ! half of previous version
   l2dom_f  = 0.08    ! half of previous version
   l3dom_f  = 0.06    ! half of previous version
   s1dom_f   = 0.18   ! half of previous version
   s2dom_f   = 0.14   ! half of previous version
   s3dom_f   = 0.10   ! half of previous version
   s4dom_f   = 0.06   ! half of previous version

! Zheng et la. 2021; Carbon and nitrogen transfer from litter to soil is higher
! in slow than rapid decompsoing plant litter: a symthesis of stable isotope
! studies, SBB, 156, 108196.
   l1s1_f =  0.19     ! half of 0.39 in CLM-CN
   l2s2_f = 0.21    ! higher than l1s1
   l3s3_f = 0.23    ! higher than l2s2, both is less than L2 - S2 in CLM-CN
   s1s2_f = 0.14    ! half of 0.39 in CLM-CN
   s2s3_f = 0.23    ! half of 0.39 in CLM-CN
   s3s4_f = 0.27      ! half of 0.39 in CLM-CN

! update the below calucation given the above C flow form L to SOM
   l1m_fb = l1m_fb * (1.0 - l1dom_f - l1s1_f)
   l1m_ff = l1m_ff * (1.0 - l1dom_f - l1s1_f)
   l2m_fb = l2m_fb * (1.0 - l2dom_f - l2s2_f)
   l2m_ff = l2m_ff * (1.0 - l2dom_f - l2s2_f)
   l3m_fb = l3m_fb * (1.0 - l3dom_f - l3s3_f)
   l3m_ff = l3m_ff * (1.0 - l3dom_f - l3s3_f)
   s1m_fb = s1m_fb * (1.0 - s1dom_f - s1s2_f)
   s1m_ff = s1m_ff * (1.0 - s1dom_f - s1s2_f)
   s2m_fb = s2m_fb * (1.0 - s2dom_f - s2s3_f)
   s2m_ff = s2m_ff * (1.0 - s2dom_f - s2s3_f)
   s3m_fb = s3m_fb * (1.0 - s3dom_f - s3s4_f)
   s3m_ff = s3m_ff * (1.0 - s3dom_f - s3s4_f)
   s4m_fb = s4m_fb * (1.0 - s4dom_f)
   s4m_ff = s4m_ff * (1.0 - s4dom_f)

   !write(*,*) "debug:l1m_ff:",l1m_ff

   !batm_f = 0.05
   !bdom_f = 0.25
   !bs1_f = 0.1
   !bs2_f = 0.12
   !bs3_f = 0.18
   !bs4_f = 0.30
   !fatm_f = 0.05
   !fdom_f = 0.25
   !fs1_f = 0.1
   !fs2_f = 0.12
   !fs3_f = 0.18
   !fs4_f = 0.30
   !domb_f = 0.3
   !domf_f = 0.3
   !doms1_f = 0.2
   !doms2_f = 0.15
   !doms3_f = 0.05
   !doms4_f = 0.0
   
   batm_f = m_batm_f(pft_index(:))
   bdom_f = m_bdom_f(pft_index(:))
   bs1_f = m_bs1_f(pft_index(:))
   bs2_f = m_bs2_f(pft_index(:))
   bs3_f = m_bs3_f(pft_index(:))
   bs4_f = 1.0 - batm_f - bdom_f - bs1_f - bs2_f - bs3_f
   fatm_f = m_fatm_f(pft_index(:))
   fdom_f = m_fdom_f(pft_index(:))
   fs1_f = m_fs1_f(pft_index(:))
   fs2_f = m_fs2_f(pft_index(:))
   fs3_f = m_fs3_f(pft_index(:))
   fs4_f = 1.0 - fatm_f - fdom_f - fs1_f - fs2_f - fs3_f
   domb_f = m_domb_f(pft_index(:))
   domf_f = m_domf_f(pft_index(:))
   doms1_f = m_doms1_f(pft_index(:))
   doms2_f = m_doms2_f(pft_index(:))
   doms3_f = m_doms3_f(pft_index(:))
   doms4_f = 1.0 - domb_f - domf_f - doms1_f - doms2_f - doms3_f
   bs4_f = max(0._rk8, bs4_f)
   fs4_f = max(0._rk8, fs4_f)
   doms4_f = max(0._rk8, doms4_f)

#endif

    ! list of transitions and their time-independent coefficients
    !!! zhangq add for Microbe
#ifdef MICROBE

   do c = begc,endc

       l =clm3%g%l%c%landunit(c)
       if (clm3%g%l%itype(l) == istsoil .or. clm3%g%l%itype(l) == istcrop) then

   i_cwdl2 = 1
   cascade_step_name(i_cwdl2) = 'CWDL2'
   rf_decomp_cascade(c,1:nlevdecomp,i_cwdl2) = 0._rk8
   cascade_donor_pool(i_cwdl2) = i_cwd
   cascade_receiver_pool(i_cwdl2) = i_litr2
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_cwdl2) = cwd_fcel

   i_cwdl3 = 2
   cascade_step_name(i_cwdl3) = 'CWDL3'
   rf_decomp_cascade(c,1:nlevdecomp,i_cwdl3) = 0._rk8
   cascade_donor_pool(i_cwdl3) = i_cwd
   cascade_receiver_pool(i_cwdl3) = i_litr3
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_cwdl3) = cwd_flig
   
   i_l1b = 3
   cascade_step_name(i_l1b) = 'L1B'
   rf_decomp_cascade(c,1:nlevdecomp,i_l1b) = 1.0 - rf_l1m
   cascade_donor_pool(i_l1b) = i_litr1
   cascade_receiver_pool(i_l1b) = i_bacteria
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_l1b) = l1m_fb(c-begc+1) !* rf_l1m

   i_l1f = 4
   cascade_step_name(i_l1f) = 'L1F'
   rf_decomp_cascade(c,1:nlevdecomp,i_l1f) = 1.0 - rf_l1m
   cascade_donor_pool(i_l1f) = i_litr1
   cascade_receiver_pool(i_l1f) = i_fungi
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_l1f) = l1m_ff(c-begc+1) !* rf_l1m

   i_l1s1 = 5
   cascade_step_name(i_l1s1) = 'L1S1'
   rf_decomp_cascade(c,1:nlevdecomp,i_l1s1) = 0.0
   cascade_donor_pool(i_l1s1) = i_litr1
   cascade_receiver_pool(i_l1s1) = i_soil1
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_l1s1) = l1s1_f
   
   i_l2b = 6
   cascade_step_name(i_l2b) = 'L2B'
   rf_decomp_cascade(c,1:nlevdecomp,i_l2b) = 1.0 - rf_l2m
   cascade_donor_pool(i_l2b) = i_litr2
   cascade_receiver_pool(i_l2b) = i_bacteria
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_l2b) = l2m_fb(c-begc+1) !* rf_l2m

   i_l2f = 7
   cascade_step_name(i_l2f) = 'L2F'
   rf_decomp_cascade(c,1:nlevdecomp,i_l2f) = 1.0 - rf_l2m
   cascade_donor_pool(i_l2f) = i_litr2
   cascade_receiver_pool(i_l2f) = i_fungi
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_l2f) = l2m_ff(c-begc+1) !* rf_l2m

   i_l2s2 = 8
   cascade_step_name(i_l2s2) = 'L2S2'
   rf_decomp_cascade(c,1:nlevdecomp,i_l2s2) = 0.0
   cascade_donor_pool(i_l2s2) = i_litr2
   cascade_receiver_pool(i_l2s2) = i_soil2
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_l2s2) = l2s2_f

   i_l3b = 9
   cascade_step_name(i_l3b) = 'L3B'
   rf_decomp_cascade(c,1:nlevdecomp,i_l3b) = 1.0 - rf_l3m
   cascade_donor_pool(i_l3b) = i_litr3
   cascade_receiver_pool(i_l3b) = i_bacteria
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_l3b) = l3m_fb(c-begc+1) !* rf_l3m

   i_l3f = 10
   cascade_step_name(i_l3f) = 'L3F'
   rf_decomp_cascade(c,1:nlevdecomp,i_l3f) = 1.0 - rf_l3m
   cascade_donor_pool(i_l3f) = i_litr3
   cascade_receiver_pool(i_l3f) = i_fungi
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_l3f) = l3m_ff(c-begc+1) !* rf_l3m

   i_l3s3 = 11
   cascade_step_name(i_l3s3) = 'L3S3'
   rf_decomp_cascade(c,1:nlevdecomp,i_l3s3) = 0.0
   cascade_donor_pool(i_l3s3) = i_litr3
   cascade_receiver_pool(i_l3s3) = i_soil3
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_l3s3) = l3s3_f

   i_s1b = 12
   cascade_step_name(i_s1b) = 'S1B'
   rf_decomp_cascade(c,1:nlevdecomp,i_s1b) = 1.0 - rf_s1m(c-begc+1)
   cascade_donor_pool(i_s1b) = i_soil1
   cascade_receiver_pool(i_s1b) = i_bacteria
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_s1b) = s1m_fb(c-begc+1) !* rf_s1m(c-begc+1)

   i_s1f = 13
   cascade_step_name(i_s1f) = 'S1F'
   rf_decomp_cascade(c,1:nlevdecomp,i_s1f) = 1.0 - rf_s1m(c-begc+1)
   cascade_donor_pool(i_s1f) = i_soil1
   cascade_receiver_pool(i_s1f) = i_fungi
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_s1f) = s1m_ff(c-begc+1) !* rf_s1m(c-begc+1)

   i_s1s2 = 14
   cascade_step_name(i_s1s2) = 'S1S2'
   rf_decomp_cascade(c,1:nlevdecomp,i_s1s2) = 0.0
   cascade_donor_pool(i_s1s2) = i_soil1
   cascade_receiver_pool(i_s1s2) = i_soil2
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_s1s2) = s1s2_f

   i_s2b = 15
   cascade_step_name(i_s2b) = 'S2B'
   rf_decomp_cascade(c,1:nlevdecomp,i_s2b) = 1.0 - rf_s2m(c-begc+1)
   cascade_donor_pool(i_s2b) = i_soil2
   cascade_receiver_pool(i_s2b) = i_bacteria
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_s2b) = s2m_fb(c-begc+1) !* rf_s2m(c-begc+1)

   i_s2f = 16
   cascade_step_name(i_s2f) = 'S2F'
   rf_decomp_cascade(c,1:nlevdecomp,i_s2f) = 1.0 - rf_s2m(c-begc+1)
   cascade_donor_pool(i_s2f) = i_soil2
   cascade_receiver_pool(i_s2f) = i_fungi
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_s2f) = s2m_ff(c-begc+1) !* rf_s2m(c-begc+1)

   i_s2s3 = 17
   cascade_step_name(i_s2s3) = 'S2S3'
   rf_decomp_cascade(c,1:nlevdecomp,i_s2s3) = 0.0
   cascade_donor_pool(i_s2s3) = i_soil2
   cascade_receiver_pool(i_s2s3) = i_soil3
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_s2s3) = s2s3_f

   i_s3b = 18
   cascade_step_name(i_s3b) = 'S3B'
   rf_decomp_cascade(c,1:nlevdecomp,i_s3b) = 1.0 - rf_s3m(c-begc+1)
   cascade_donor_pool(i_s3b) = i_soil3
   cascade_receiver_pool(i_s3b) = i_bacteria
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_s3b) = s3m_fb(c-begc+1) !* rf_s3m(c-begc+1)

   i_s3f = 19
   cascade_step_name(i_s3f) = 'S3F'
   rf_decomp_cascade(c,1:nlevdecomp,i_s3f) = 1.0 - rf_s3m(c-begc+1)
   cascade_donor_pool(i_s3f) = i_soil3
   cascade_receiver_pool(i_s3f) = i_fungi
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_s3f) = s3m_ff(c-begc+1) !* rf_s3m(c-begc+1)

   i_s3s4 = 20
   cascade_step_name(i_s3s4) = 'S3S4'
   rf_decomp_cascade(c,1:nlevdecomp,i_s3s4) = 0.0
   cascade_donor_pool(i_s3s4) = i_soil3
   cascade_receiver_pool(i_s3s4) = i_soil4
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_s3s4) = s3s4_f

   i_s4b = 21
   cascade_step_name(i_s4b) = 'S4B'
   rf_decomp_cascade(c,1:nlevdecomp,i_s4b) = 1.0 - rf_s4m(c-begc+1)
   cascade_donor_pool(i_s4b) = i_soil4
   cascade_receiver_pool(i_s4b) = i_bacteria
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_s4b) = s4m_fb(c-begc+1) !* rf_s4m(c-begc+1)

   i_s4f = 22
   cascade_step_name(i_s4f) = 'S4F'
   rf_decomp_cascade(c,1:nlevdecomp,i_s4f) = 1.0 - rf_s4m(c-begc+1)
   cascade_donor_pool(i_s4f) = i_soil4
   cascade_receiver_pool(i_s4f) = i_fungi
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_s4f) = s4m_ff(c-begc+1) !* rf_s4m(c-begc+1)
   
   i_bs1 = 23
   cascade_step_name(i_bs1) = 'BS1'
   rf_decomp_cascade(c,1:nlevdecomp,i_bs1) = 0._rk8
   cascade_donor_pool(i_bs1) = i_bacteria
   cascade_receiver_pool(i_bs1) = i_soil1
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_bs1) = bs1_f(c-begc+1)

   i_fs1 = 24
   cascade_step_name(i_fs1) = 'FS1'
   rf_decomp_cascade(c,1:nlevdecomp,i_fs1) = 0._rk8
   cascade_donor_pool(i_fs1) = i_fungi
   cascade_receiver_pool(i_fs1) = i_soil1
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_fs1) = fs1_f(c-begc+1)
   
   i_bs2 = 25
   cascade_step_name(i_bs2) = 'BS2'
   rf_decomp_cascade(c,1:nlevdecomp,i_bs2) = 0._rk8
   cascade_donor_pool(i_bs2) = i_bacteria
   cascade_receiver_pool(i_bs2) = i_soil2
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_bs2) = bs2_f(c-begc+1)

   i_fs2 = 26
   cascade_step_name(i_fs2) = 'FS2'
   rf_decomp_cascade(c,1:nlevdecomp,i_fs2) = 0._rk8
   cascade_donor_pool(i_fs2) = i_fungi
   cascade_receiver_pool(i_fs2) = i_soil2
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_fs2) = fs2_f(c-begc+1)

   i_bs3 = 27
   cascade_step_name(i_bs3) = 'BS3'
   rf_decomp_cascade(c,1:nlevdecomp,i_bs3) = 0._rk8
   cascade_donor_pool(i_bs3) = i_bacteria
   cascade_receiver_pool(i_bs3) = i_soil3
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_bs3) =bs3_f(c-begc+1)

   i_fs3 = 28
   cascade_step_name(i_fs3) = 'FS3'
   rf_decomp_cascade(c,1:nlevdecomp,i_fs3) = 0._rk8
   cascade_donor_pool(i_fs3) = i_fungi
   cascade_receiver_pool(i_fs3) = i_soil3
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_fs3) = fs3_f(c-begc+1)

   i_bs4 = 29
   cascade_step_name(i_bs4) = 'BS4'
   rf_decomp_cascade(c,1:nlevdecomp,i_bs4) = 0._rk8
   cascade_donor_pool(i_bs4) = i_bacteria
   cascade_receiver_pool(i_bs4) = i_soil4
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_bs4) = bs4_f(c-begc+1)

   i_fs4 = 30
   cascade_step_name(i_fs4) = 'FS4'
   rf_decomp_cascade(c,1:nlevdecomp,i_fs4) = 0._rk8
   cascade_donor_pool(i_fs4) = i_fungi
   cascade_receiver_pool(i_fs4) = i_soil4
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_fs4) = fs4_f(c-begc+1)

   i_doms1 = 31
   cascade_step_name(i_doms1) = 'DOMS1'
   rf_decomp_cascade(c,1:nlevdecomp,i_doms1) = 0._rk8
   cascade_donor_pool(i_doms1) = i_dom
   cascade_receiver_pool(i_doms1) = i_soil1
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_doms1) = doms1_f(c-begc+1)
   
   i_doms2 = 32
   cascade_step_name(i_doms2) = 'DOMS2'
   rf_decomp_cascade(c,1:nlevdecomp,i_doms2) = 0._rk8
   cascade_donor_pool(i_doms2) = i_dom
   cascade_receiver_pool(i_doms2) = i_soil2
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_doms2) = doms2_f(c-begc+1)

   i_doms3 = 33
   cascade_step_name(i_doms3) = 'DOMS3'
   rf_decomp_cascade(c,1:nlevdecomp,i_doms3) = 0._rk8
   cascade_donor_pool(i_doms3) = i_dom
   cascade_receiver_pool(i_doms3) = i_soil3
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_doms3) = doms3_f(c-begc+1)

   i_doms4 = 34
   cascade_step_name(i_doms4) = 'DOMS4'
   rf_decomp_cascade(c,1:nlevdecomp,i_doms4) = 0._rk8
   cascade_donor_pool(i_doms4) = i_dom
   cascade_receiver_pool(i_doms4) = i_soil4
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_doms4) = doms4_f(c-begc+1)

   i_bdom = 35
   cascade_step_name(i_bdom) = 'BDOM'
   rf_decomp_cascade(c,1:nlevdecomp,i_bdom) = 0._rk8
   cascade_donor_pool(i_bdom) = i_bacteria
   cascade_receiver_pool(i_bdom) = i_dom
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_bdom) = bdom_f(c-begc+1)

   i_fdom = 36
   cascade_step_name(i_fdom) = 'FDOM'
   rf_decomp_cascade(c,1:nlevdecomp,i_fdom) = 0._rk8
   cascade_donor_pool(i_fdom) = i_fungi
   cascade_receiver_pool(i_fdom) = i_dom
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_fdom) = fdom_f(c-begc+1)

   i_domb = 37
   cascade_step_name(i_domb) = 'DOMB'
   rf_decomp_cascade(c,1:nlevdecomp,i_domb) = 0._rk8
   cascade_donor_pool(i_domb) = i_dom
   cascade_receiver_pool(i_domb) = i_bacteria
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_domb) = domb_f(c-begc+1)
  
   i_domf = 38
   cascade_step_name(i_domf) = 'DOMF'
   rf_decomp_cascade(c,1:nlevdecomp,i_domf) = 0._rk8
   cascade_donor_pool(i_domf) = i_dom
   cascade_receiver_pool(i_domf) = i_fungi
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_domf) = domf_f(c-begc+1)
   
   i_batm = 39
   cascade_step_name(i_batm) = 'BATM'
   rf_decomp_cascade(c,1:nlevdecomp,i_batm) = 1.0
   cascade_donor_pool(i_batm) = i_bacteria
   cascade_receiver_pool(i_batm) = i_atm
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_batm) = batm_f(c-begc+1)
   
   i_fatm = 40
   cascade_step_name(i_fatm) = 'FATM'
   rf_decomp_cascade(c,1:nlevdecomp,i_fatm) = 1.0
   cascade_donor_pool(i_fatm) = i_fungi
   cascade_receiver_pool(i_fatm) = i_atm
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_fatm) = fatm_f(c-begc+1)
    
   i_l1dom = 41
   cascade_step_name(i_l1dom) = 'L1DOM'
   rf_decomp_cascade(c,1:nlevdecomp,i_l1dom) = 0._rk8
   cascade_donor_pool(i_l1dom) = i_litr1
   cascade_receiver_pool(i_l1dom) = i_dom
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_l1dom) = l1dom_f
   
   i_l2dom = 42
   cascade_step_name(i_l2dom) = 'L2DOM'
   rf_decomp_cascade(c,1:nlevdecomp,i_l2dom) = 0._rk8
   cascade_donor_pool(i_l2dom) = i_litr2
   cascade_receiver_pool(i_l2dom) = i_dom
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_l2dom) = l2dom_f

   i_l3dom = 43
   cascade_step_name(i_l3dom) = 'L3DOM'
   rf_decomp_cascade(c,1:nlevdecomp,i_l3dom) = 0._rk8
   cascade_donor_pool(i_l3dom) = i_litr3
   cascade_receiver_pool(i_l3dom) = i_dom
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_l3dom) = l3dom_f
   
   i_s1dom = 44
   cascade_step_name(i_s1dom) = 'S1DOM'
   rf_decomp_cascade(c,1:nlevdecomp,i_s1dom) = 0._rk8
   cascade_donor_pool(i_s1dom) = i_soil1
   cascade_receiver_pool(i_s1dom) = i_dom
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_s1dom) = s1dom_f
   
   i_s2dom = 45
   cascade_step_name(i_s2dom) = 'S2DOM'
   rf_decomp_cascade(c,1:nlevdecomp,i_s2dom) = 0._rk8
   cascade_donor_pool(i_s2dom) = i_soil2
   cascade_receiver_pool(i_s2dom) = i_dom
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_s2dom) = s2dom_f
   
   i_s3dom = 46
   cascade_step_name(i_s3dom) = 'S3DOM'
   rf_decomp_cascade(c,1:nlevdecomp,i_s3dom) = 0._rk8
   cascade_donor_pool(i_s3dom) = i_soil3
   cascade_receiver_pool(i_s3dom) = i_dom
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_s3dom) = s3dom_f 
   
   i_s4dom = 47
   cascade_step_name(i_s4dom) = 'S4DOM'
   rf_decomp_cascade(c,1:nlevdecomp,i_s4dom) = 0._rk8
   cascade_donor_pool(i_s4dom) = i_soil4
   cascade_receiver_pool(i_s4dom) = i_dom
   pathfrac_decomp_cascade(c,1:nlevdecomp,i_s4dom) = s4dom_f

   end if

   end do

#else
   i_l1s1 = 1
   cascade_step_name(i_l1s1) = 'L1S1'
   rf_decomp_cascade(:,1:nlevdecomp,i_l1s1) = rf_l1s1
   cascade_donor_pool(i_l1s1) = i_litr1
   cascade_receiver_pool(i_l1s1) = i_soil1
   pathfrac_decomp_cascade(:,1:nlevdecomp,i_l1s1) = 1.0_rk8

   i_l2s2 = 2
   cascade_step_name(i_l2s2) = 'L2S2'
   rf_decomp_cascade(:,1:nlevdecomp,i_l2s2) = rf_l2s2
   cascade_donor_pool(i_l2s2) = i_litr2
   cascade_receiver_pool(i_l2s2) = i_soil2
   pathfrac_decomp_cascade(:,1:nlevdecomp,i_l2s2) = 1.0_rk8

   i_l3s3 = 3
   cascade_step_name(i_l3s3) = 'L3S3'
   rf_decomp_cascade(:,1:nlevdecomp,i_l3s3) = rf_l3s3
   cascade_donor_pool(i_l3s3) = i_litr3
   cascade_receiver_pool(i_l3s3) = i_soil3
   pathfrac_decomp_cascade(:,1:nlevdecomp,i_l3s3) = 1.0_rk8

   i_s1s2 = 4
   cascade_step_name(i_s1s2) = 'S1S2'
   rf_decomp_cascade(:,1:nlevdecomp,i_s1s2) = rf_s1s2
   cascade_donor_pool(i_s1s2) = i_soil1
   cascade_receiver_pool(i_s1s2) = i_soil2
   pathfrac_decomp_cascade(:,1:nlevdecomp,i_s1s2) = 1.0_rk8

   i_s2s3 = 5
   cascade_step_name(i_s2s3) = 'S2S3'
   rf_decomp_cascade(:,1:nlevdecomp,i_s2s3) = rf_s2s3
   cascade_donor_pool(i_s2s3) = i_soil2
   cascade_receiver_pool(i_s2s3) = i_soil3
   pathfrac_decomp_cascade(:,1:nlevdecomp,i_s2s3) = 1.0_rk8

   i_s3s4 = 6 
   cascade_step_name(i_s3s4) = 'S3S4'
   rf_decomp_cascade(:,1:nlevdecomp,i_s3s4) = rf_s3s4
   cascade_donor_pool(i_s3s4) = i_soil3
   cascade_receiver_pool(i_s3s4) = i_soil4
   pathfrac_decomp_cascade(:,1:nlevdecomp,i_s3s4) = 1.0_rk8

   i_s4atm = 7
   cascade_step_name(i_s4atm) = 'S4'
   rf_decomp_cascade(:,1:nlevdecomp,i_s4atm) = 1.
   cascade_donor_pool(i_s4atm) = i_soil4
   cascade_receiver_pool(i_s4atm) = i_atm
   pathfrac_decomp_cascade(:,1:nlevdecomp,i_s4atm) = 1.0_rk8

   i_cwdl2 = 8
   cascade_step_name(i_cwdl2) = 'CWDL2'
   rf_decomp_cascade(:,1:nlevdecomp,i_cwdl2) = 0._rk8
   cascade_donor_pool(i_cwdl2) = i_cwd
   cascade_receiver_pool(i_cwdl2) = i_litr2
   pathfrac_decomp_cascade(:,1:nlevdecomp,i_cwdl2) = cwd_fcel

   i_cwdl3 = 9
   cascade_step_name(i_cwdl3) = 'CWDL3'
   rf_decomp_cascade(:,1:nlevdecomp,i_cwdl3) = 0._rk8
   cascade_donor_pool(i_cwdl3) = i_cwd
   cascade_receiver_pool(i_cwdl3) = i_litr3
   pathfrac_decomp_cascade(:,1:nlevdecomp,i_cwdl3) = cwd_flig
!!! zhangq modified for Microbe
#ifdef MICROBE
   !deallocate(rf_s1m)
   !deallocate(rf_s2m)
   !deallocate(rf_s3m)
   !deallocate(rf_s4m)
   !deallocate(l1m_fb)
   !deallocate(l1m_ff)
   !deallocate(l2m_fb)
   !deallocate(l2m_ff)
   !deallocate(l3m_fb)
   !deallocate(l3m_ff)
   !deallocate(s1m_fb)
   !deallocate(s1m_ff)
   !deallocate(s2m_fb)
   !deallocate(s2m_ff)
   !deallocate(s3m_fb)
   !deallocate(s3m_ff)
   !deallocate(s4m_fb)
   !deallocate(s4m_ff)
   !deallocate(batm_f)
   !deallocate(bdom_f)
   !deallocate(bs1_f)
   !deallocate(bs2_f)
   !deallocate(bs3_f)
   !deallocate(bs4_f)
   !deallocate(fatm_f)
   !deallocate(fdom_f)
   !deallocate(fs1_f)
   !deallocate(fs2_f)
   !deallocate(fs3_f)
   !deallocate(fs4_f)
   !deallocate(domb_f)
   !deallocate(domf_f)
   !deallocate(doms1_f)
   !deallocate(doms2_f)
   !deallocate(doms3_f)
   !deallocate(doms4_f)
   !deallocate(cn_fungi_in)
   !deallocate(cn_bacteria_in)
#endif

#endif

  end subroutine init_decompcascade

  subroutine decomp_rate_constants(lbc, ubc, num_soilc, filter_soilc)
    !
    ! calculate rate constants and decomposition pathways for the BGC
    ! model originally implemented in CLM-CN
    ! written by C. Koven based on original CLM4 decomposition cascade
    ! by P. Thornton
    !
    use mod_clm_type
    use mod_clm_varcon, only: secspday
    implicit none
    integer(ik4), intent(in) :: lbc, ubc   ! column bounds
    integer(ik4), intent(in) :: num_soilc  ! number of soil columns in filter
    integer(ik4), intent(in) :: filter_soilc(:) ! filter for soil columns
    ! column level
    ! rate constant for decomposition (1./sec)
    real(rk8), pointer :: decomp_k(:,:,:)
    real(rk8), pointer :: t_scalar(:,:) ! soil temperature scalar for decomp
    real(rk8), pointer :: w_scalar(:,:) ! soil water scalar for decomp
    ! fraction by which decomposition is limited by anoxia
    real(rk8), pointer :: o_scalar(:,:)

    real(rk8), pointer :: dz(:,:)   ! soil layer thickness (m)
    ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(rk8), pointer :: t_soisno(:,:)
    real(rk8), pointer :: sucsat(:,:)   ! minimum soil suction (mm)
    ! soil water potential in each soil layer (MPa)
    real(rk8), pointer :: soilpsi(:,:)
!!! zhangq add for Microbe
#if (defined LCH4) || (defined MICROBE)
    ! Ratio of oxygen available to that demanded by roots,
    ! aerobes, & methanotrophs (nlevsoi)
    real(rk8), pointer :: o2stress_unsat(:,:)
    ! Ratio of oxygen available to that demanded by roots,
    ! aerobes, & methanotrophs (nlevsoi)
    real(rk8), pointer :: o2stress_sat(:,:)
    ! fractional inundated area (excluding dedicated wetland columns)
    real(rk8), pointer :: finundated(:)
#endif
    integer(ik4), pointer :: alt_indx(:)  ! current depth of thaw

    real(rk8) :: dt  ! decomp timestep (seconds)
    real(rk8):: dtd  ! decomp timestep (days)

    real(rk8):: frw(lbc:ubc)      ! rooting fraction weight
    real(rk8), pointer:: fr(:,:)  ! column-level rooting fraction by soil depth
    real(rk8):: minpsi, maxpsi    ! limits for soil water scalar for decomp
    real(rk8):: psi               ! temporary soilpsi for water scalar
    ! real(rk8):: w_scalar(lbc:ubc,1:nlevdecomp) !soil water scalar for decomp
    real(rk8):: rate_scalar  ! combined rate scalar for decomp

    real(rk8):: k_l1         ! decomposition rate constant litter 1
    real(rk8):: k_l2         ! decomposition rate constant litter 2
    real(rk8):: k_l3         ! decomposition rate constant litter 3
    real(rk8):: k_s1         ! decomposition rate constant SOM 1
    real(rk8):: k_s2         ! decomposition rate constant SOM 2
    real(rk8):: k_s3         ! decomposition rate constant SOM 3
    real(rk8):: k_s4         ! decomposition rate constant SOM 3
    real(rk8):: k_frag       ! fragmentation rate constant CWD
    real(rk8):: ck_l1        ! corrected decomposition rate constant litter 1
    real(rk8):: ck_l2        ! corrected decomposition rate constant litter 2
    real(rk8):: ck_l3        ! corrected decomposition rate constant litter 3
    real(rk8):: ck_s1        ! corrected decomposition rate constant SOM 1
    real(rk8):: ck_s2        ! corrected decomposition rate constant SOM 2
    real(rk8):: ck_s3        ! corrected decomposition rate constant SOM 3
    real(rk8):: ck_s4        ! corrected decomposition rate constant SOM 3
    real(rk8):: ck_frag      ! corrected fragmentation rate constant CWD
    real(rk8):: cwd_fcel     ! cellulose fraction of coarse woody debris
    real(rk8):: cwd_flig     ! lignin fraction of coarse woody debris
    real(rk8):: cwdc_loss    ! fragmentation rate for CWD carbon (gC/m2/s)
    real(rk8):: cwdn_loss    ! fragmentation rate for CWD nitrogen (gN/m2/s)
!!! zhangq add for Microbe
#ifdef MICROBE

    real(rk8),allocatable :: decomp_depth_efolding_in(:) ! define the active decomposing soil depth
    real(rk8),allocatable :: k_dom_in(:)         ! decomposition rate constant dissolved organic matter, revised to read in from parameter file on April 16, 2016
    real(rk8),allocatable :: k_bacteria_in(:)        ! decomposition rate constant biomass of bacteria, revised to read in from parameter file on April 16, 2016
    real(rk8),allocatable :: k_fungi_in(:)         ! decomposition rate constant fungi biomass, revised to read in from parameter file on April 16, 2016

    real(rk8), pointer :: decomp_depth_efolding(:) ! define the active decomposing soil depth
    real(rk8), pointer :: k_dom(:)     ! turnover rate of DOM
    real(rk8), pointer :: k_bacteria(:)     ! turnover rate of bacteria
    real(rk8), pointer :: k_fungi(:)     ! turnover rate of fungi

    real(rk8),allocatable :: ck_dom(:)        ! corrected decomposition rate constant dissolved organic matter
    real(rk8),allocatable :: ck_bacteria(:)        ! corrected decomposition rate constant bacteria biomass
    real(rk8),allocatable :: ck_fungi(:)        ! corrected decomposition rate constant fungi biomass

    real(rk8), pointer :: wtcol(:)                ! pft weight relative to column (0-1)
    integer(ik4), allocatable :: pft_index(:)
    !integer(ik4) :: pft_index(:)
    integer(ik4), pointer :: pfti(:)        ! pft index array
    integer(ik4), pointer :: pftf(:)        ! pft index array

#endif

    integer(ik4) :: i_litr1
    integer(ik4) :: i_litr2
    integer(ik4) :: i_litr3
    integer(ik4) :: i_soil1
    integer(ik4) :: i_soil2
    integer(ik4) :: i_soil3
    integer(ik4) :: i_soil4
    integer(ik4) :: c, fc, j

!!! zhangq add for Microbe
#if (defined VERTSOILC) || (defined MICROBE)
    real(rk8) :: depth_scalar(lbc:ubc,1:nlevdecomp)
#endif

!!! zhangq modified for Microbe
#ifdef MICROBE
    allocate(decomp_depth_efolding_in(1:num_soilc)) ! define the active decomposing soil depth
    allocate(k_dom_in(1:num_soilc))         ! decomposition rate constant dissolved organic matter, revised to read in from parameter file on April 16, 2016
    allocate(k_bacteria_in(1:num_soilc))        ! decomposition rate constant biomass of bacteria, revised to read in from parameter file on April 16, 2016
    allocate(k_fungi_in(1:num_soilc))         ! decomposition rate constant fungi biomass, revised to read in from parameter file on April 16, 2016

    allocate(ck_dom(1:num_soilc))        ! corrected decomposition rate constant dissolved organic matter
    allocate(ck_bacteria(1:num_soilc))        ! corrected decomposition rate constant bacteria biomass
    allocate(ck_fungi(1:num_soilc))        ! corrected decomposition rate constant fungi biomass
#endif

    ! Assign local pointers to derived type arrays
    t_soisno              => clm3%g%l%c%ces%t_soisno
    sucsat                => clm3%g%l%c%cps%sucsat
    soilpsi               => clm3%g%l%c%cps%soilpsi
    dz                    => clm3%g%l%c%cps%dz
    t_scalar              => clm3%g%l%c%ccf%t_scalar
    w_scalar              => clm3%g%l%c%ccf%w_scalar
    o_scalar              => clm3%g%l%c%ccf%o_scalar
    decomp_k              => clm3%g%l%c%ccf%decomp_k
#ifdef LCH4
    o2stress_sat          => clm3%g%l%c%cch4%o2stress_sat
    o2stress_unsat        => clm3%g%l%c%cch4%o2stress_unsat
    finundated            => clm3%g%l%c%cws%finundated
#endif

!!! zhangq add for Microbe
#ifdef MICROBE
   !   o2stress_sat          => cch4%o2stress_sat
   !   o2stress_unsat        => cch4%o2stress_unsat
   !   finundated            => cws%finundated
   !depth_scalar  => ccf%depth_scalar

   decomp_depth_efolding =>pftcon%decomp_depth_efolding

   k_bacteria            => pftcon%k_bacteria
   k_fungi               => pftcon%k_fungi
   k_dom                 => pftcon%k_dom

   wtcol                 =>clm3%g%l%c%p%wtcol
   pfti                  =>clm3%g%l%c%pfti
   pftf                  =>clm3%g%l%c%pftf

!!! zhangq modified for Microbe
   allocate(pft_index(1:num_soilc))

  do fc = 1,num_soilc
  c = filter_soilc(fc)

    if (wtcol(clm3%g%l%c%pfti(c)) .ge. 0.96_rk8 .or. SUM(wtcol((clm3%g%l%c%pfti(c) + 1):clm3%g%l%c%pftf(c))) .lt. 0.001_rk8) then
       !pft_index = [pft_index, 0]
       pft_index(fc) = 0
    else
       pft_index(fc) = MAXLOC(wtcol((clm3%g%l%c%pfti(c) + 1):clm3%g%l%c%pftf(c)),DIM=1)
       !pft_index = [pft_index, MAXLOC(wtcol((clm3%g%l%c%pfti(c) + 1):clm3%g%l%c%pftf(c)), DIM=1)]
    end if

  ! write (*, *) "SUM(wtcol((col%pfti(c) + 1):col%pftf(c)))",
  ! SUM(wtcol((col%pfti(c) + 1):col%pftf(c)))
  ! write (*, *) "pft_index 2", pft_index

  end do

  k_dom_in = k_dom(pft_index(:))
  k_bacteria_in = k_bacteria(pft_index(:))
  k_fungi_in = k_fungi(pft_index(:))
  decomp_depth_efolding_in = decomp_depth_efolding(pft_index(:))

  ! write (*, *) "k_dom_in 2", k_dom_in
  ! write (*, *) "k_bacteria_in 2", k_bacteria_in
  ! write (*, *) "k_fungi_in 2", k_fungi_in
  ! write (*, *) "decomp_depth_efolding_in 2", decomp_depth_efolding_in

#endif

    alt_indx              => clm3%g%l%c%cps%alt_indx

    ! set time steps
    dt = dtsrf
    dtd = dt/secspday

    ! set initial base rates for decomposition mass loss (1/day)
    ! (from Biome-BGC v4.2.0, using three SOM pools)
    ! Value inside log function is the discrete-time values for a
    ! daily time step model, and the result of the log function is
    ! the corresponding continuous-time decay rate (1/day), following
    ! Olson, 1963.
    k_l1 = -log(1.0_rk8-0.7_rk8)
    k_l2 = -log(1.0_rk8-0.07_rk8)
    k_l3 = -log(1.0_rk8-0.014_rk8)
    k_s1 = -log(1.0_rk8-0.07_rk8)
    k_s2 = -log(1.0_rk8-0.014_rk8)
    k_s3 = -log(1.0_rk8-0.0014_rk8)
    k_s4 = -log(1.0_rk8-0.0001_rk8)
    k_frag = -log(1.0_rk8-0.001_rk8)

!!! zhangq add for Microbe
#ifdef MICROBE
   ck_dom = -log(1.0_rk8 - k_dom_in)
   ck_bacteria = -log(1.0_rk8 - k_bacteria_in)
   ck_fungi = -log(1.0_rk8 - k_fungi_in)
#endif

    ! calculate the new discrete-time decay rate for model timestep
    k_l1 = 1.0_rk8-exp(-k_l1*dtd)
    k_l2 = 1.0_rk8-exp(-k_l2*dtd)
    k_l3 = 1.0_rk8-exp(-k_l3*dtd)
    k_s1 = 1.0_rk8-exp(-k_s1*dtd)
    k_s2 = 1.0_rk8-exp(-k_s2*dtd)
    k_s3 = 1.0_rk8-exp(-k_s3*dtd)
    k_s4 = 1.0_rk8-exp(-k_s4*dtd)
    k_frag = 1.0_rk8-exp(-k_frag*dtd)

#ifdef MICROBE
    ck_dom = 1.0_rk8-exp(-ck_dom*dtd)
    ck_bacteria = 1.0_rk8-exp(-ck_bacteria*dtd)
    ck_fungi = 1.0_rk8-exp(-ck_fungi*dtd)
#endif

    ! The following code implements the acceleration part of the AD spinup
    ! algorithm, by multiplying all of the SOM decomposition base rates by 10.0.

    if ( spinup_state == 1 ) then
      k_s1 = k_s1 * spinup_vector(1)
      k_s2 = k_s2 * spinup_vector(2)
      k_s3 = k_s3 * spinup_vector(3)
      k_s4 = k_s4 * spinup_vector(4)
    end if

    i_litr1 = 1
    i_litr2 = 2
    i_litr3 = 3
    i_soil1 = 5
    i_soil2 = 6
    i_soil3 = 7
    i_soil4 = 8

    ! ! CWD fragmentation -> litter pools
    !  thse have now been put into the regular decomposition cascade
    ! cwdc_loss = cwdc_vr(c,j) * ck_frag / dt
    ! cwdc_to_litr2c_vr(c,j) = cwdc_loss * cwd_fcel
    ! cwdc_to_litr3c_vr(c,j) = cwdc_loss * cwd_flig
    ! cwdn_loss = cwdn_vr(c,j) * ck_frag / dt
    ! cwdn_to_litr2n_vr(c,j) = cwdn_loss * cwd_fcel
    ! cwdn_to_litr3n_vr(c,j) = cwdn_loss * cwd_flig

    !--- time dependent coefficients-----!
    if ( nlevdecomp == 1 ) then

      ! calculate function to weight the temperature and water potential scalars
      ! for decomposition control.

      ! the following normalizes values in fr so that they
      ! sum to 1.0 across top nlevdecomp levels on a column
      frw(lbc:ubc) = 0._rk8
      nlev_soildecomp_standard = 5
      allocate(fr(lbc:ubc,nlev_soildecomp_standard))
      do j = 1 , nlev_soildecomp_standard
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
          frw(c) = frw(c) + dz(c,j)
        end do
      end do
      do j = 1 , nlev_soildecomp_standard
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
          if ( frw(c) /= 0._rk8 ) then
            fr(c,j) = dz(c,j) / frw(c)
          else
            fr(c,j) = 0._rk8
          end if
        end do
      end do

      ! calculate rate constant scalar for soil temperature
      ! assuming that the base rate constants are assigned for non-moisture
      ! limiting conditions at 25 C.
      ! Peter Thornton: 3/13/09
      ! Replaced the Lloyd and Taylor function with a Q10 formula,
      ! with Q10 = 1.5 as part of the modifications made to improve the
      ! seasonal cycle of atmospheric CO2 concentration in global simulations.
      ! This does not impact the base rates at 25 C, which are calibrated
      ! from microcosm studies.
      do j = 1 , nlev_soildecomp_standard
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
          if ( j == 1 ) t_scalar(c,:) = 0._rk8
          ! use separate (possibly equal) t funcs above and below
          ! freezing point
          if ( t_soisno(c,j) >= tfrz ) then
            t_scalar(c,1) = t_scalar(c,1) + &
                    (1.5_rk8**((t_soisno(c,j)-(tfrz+25._rk8))/10._rk8))*fr(c,j)
          else
            t_scalar(c,1) = t_scalar(c,1) + &
                    (1.5_rk8**(-25._rk8/10._rk8))* &
                    (froz_q10**((t_soisno(c,j)-tfrz)/10._rk8))*fr(c,j)
          endif
        end do
      end do

      ! calculate the rate constant scalar for soil water content.
      ! Uses the log relationship with water potential given in
      ! Andren, O., and K. Paustian, 1987.
      ! Barley straw decomposition in the field: a comparison of models.
      ! Ecology, 68(5):1190-1200.
      ! and supported by data in
      ! Orchard, V.A., and F.J. Cook, 1983.
      ! Relationship between soil respiration
      ! and soil moisture. Soil Biol. Biochem., 15(4):447-453.

      minpsi = -10.0_rk8;

      do j = 1 , nlev_soildecomp_standard
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
          if ( j == 1 ) w_scalar(c,:) = 0._rk8
          maxpsi = sucsat(c,j) * (-9.8e-6_rk8)
          psi = min(soilpsi(c,j),maxpsi)
          ! decomp only if soilpsi is higher than minpsi
          if ( psi > minpsi ) then
            w_scalar(c,1) = w_scalar(c,1) + &
                    (log(minpsi/psi)/log(minpsi/maxpsi))*fr(c,j)
          end if
        end do
      end do

#if (defined LCH4) || (defined MICROBE)
      if ( anoxia_wtsat ) then ! Adjust for saturated fraction if unfrozen.
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
          if ( alt_indx(c) >= nlev_soildecomp_standard .and. &
               t_soisno(c,1) > tfrz ) then
            w_scalar(c,1) = w_scalar(c,1)*(1._rk8 - finundated(c)) + finundated(c)
          end if
        end do
      end if
#endif

#if (defined LCH4) || (defined MICROBE)
      ! Calculate ANOXIA
      if (anoxia) then
        ! Check for anoxia w/o LCH4 now done in controlMod.
        do j = 1 , nlev_soildecomp_standard
          do fc = 1 , num_soilc
            c = filter_soilc(fc)
            if ( j == 1 ) o_scalar(c,:) = 0._rk8
            if ( .not. anoxia_wtsat ) then
              o_scalar(c,1) = o_scalar(c,1) + &
                      fr(c,j) * max(o2stress_unsat(c,j), mino2lim)
            else
              o_scalar(c,1) = o_scalar(c,1) + fr(c,j) * &
                              (max(o2stress_unsat(c,j), mino2lim) * &
                              (1._rk8 - finundated(c)) + &
                              max(o2stress_sat(c,j), mino2lim)*finundated(c) )
            end if
          end do
        end do
      else
        o_scalar(lbc:ubc,1:nlevdecomp) = 1._rk8
      end if
#else
      o_scalar(lbc:ubc,1:nlevdecomp) = 1._rk8
#endif
      deallocate(fr)
    else

      ! calculate rate constant scalar for soil temperature
      ! assuming that the base rate constants are assigned for non-moisture
      ! limiting conditions at 25 C.
      ! Peter Thornton: 3/13/09
      ! Replaced the Lloyd and Taylor function with a Q10 formula,
      ! with Q10 = 1.5 as part of the modifications made to improve
      ! the seasonal cycle of atmospheric CO2 concentration in global
      ! simulations. This does not impact the base rates at 25 C, which
      ! are calibrated from microcosm studies.
      do j = 1 , nlevdecomp
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
          !! use separate (possibly equal) t funcs above and below
          !! freezing point
          !! t_scalar(c,j)= (1.5**((t_soisno(c,j)-(tfrz+25._rk8))/10._rk8))
          if (t_soisno(c,j) >= tfrz) then
            t_scalar(c,j)= (1.5_rk8**((t_soisno(c,j)-(tfrz+25._rk8))/10._rk8))
          else
            t_scalar(c,j)= (1.5_rk8**(-25._rk8/10._rk8)) * &
                    (froz_q10**((t_soisno(c,j)-tfrz)/10._rk8))
          endif
        end do
      end do

      ! calculate the rate constant scalar for soil water content.
      ! Uses the log relationship with water potential given in
      ! Andren, O., and K. Paustian, 1987.
      ! Barley straw decomposition in the field: a comparison of models.
      ! Ecology, 68(5):1190-1200.
      ! and supported by data in
      ! Orchard, V.A., and F.J. Cook, 1983.
      ! Relationship between soil respiration
      ! and soil moisture. Soil Biol. Biochem., 15(4):447-453.
      minpsi = -10.0_rk8;
      do j = 1 , nlevdecomp
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
          maxpsi = sucsat(c,j) * (-9.8e-6_rk8)
          psi = min(soilpsi(c,j),maxpsi)
          ! decomp only if soilpsi is higher than minpsi
          if ( psi > minpsi ) then
            w_scalar(c,j) = (log(minpsi/psi)/log(minpsi/maxpsi))
          else
            w_scalar(c,j) = 0._rk8
          end if
!!! zhangq add for Microbe
!#if (defined LCH4) || (defined MICROBE)
#ifdef LCH4
          ! wet area will have w_scalar of 1 if unfrozen
          if ( anoxia_wtsat .and. t_soisno(c,j) > tfrz ) then
            w_scalar(c,j) = w_scalar(c,j)*(1._rk8 - finundated(c)) + finundated(c)
          end if
#endif
!!! zhangq midified for Microbe
          w_scalar(c,j) = min(w_scalar(c,j), 1._rk8)
          w_scalar(c,j) = max(w_scalar(c,j), 0.3_rk8)
          w_scalar(c,j) = w_scalar(c,j)**0.5
        end do
      end do
    end if

#ifdef LCH4
    ! Calculate ANOXIA
    ! Check for anoxia w/o LCH4 now done in controlMod.
    if ( anoxia ) then
      do j = 1 , nlevdecomp
        do fc = 1 , num_soilc
          c = filter_soilc(fc)
          if ( .not. anoxia_wtsat ) then
            o_scalar(c,j) = max(o2stress_unsat(c,j), mino2lim)
          else
            o_scalar(c,j) = max(o2stress_unsat(c,j), mino2lim) * &
                            (1._rk8 - finundated(c)) + &
                            max(o2stress_sat(c,j), mino2lim) * finundated(c)
          end if
        end do
      end do
    else
      o_scalar(lbc:ubc,1:nlevdecomp) = 1._rk8
    end if
#else
    o_scalar(lbc:ubc,1:nlevdecomp) = 1._rk8
#endif

!write(*, *) "num_soilc=", num_soilc

#if (defined VERTSOILC) || (defined MICROBE)
    ! add a term to reduce decomposition rate at depth
    ! for now used a fixed e-folding depth
    do j = 1 , nlevdecomp
      do fc = 1 , num_soilc
        c = filter_soilc(fc)
        depth_scalar(c,j) = exp(-zsoi(j)/decomp_depth_efolding_in(fc))
      end do
    end do
#endif

#if (defined VERTSOILC) || (defined MICROBE)
    do j = 1 , nlevdecomp
      do fc = 1 , num_soilc
        c = filter_soilc(fc)
#ifdef MICROBE
         decomp_k(c,j,i_litr1) = k_l1 * t_scalar(c,j) * w_scalar(c,j) * &
           depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_litr2) = k_l2 * t_scalar(c,j) * w_scalar(c,j) * &
           depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_litr3) = k_l3 * t_scalar(c,j) * w_scalar(c,j) * &
           depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_cwd) = k_frag * t_scalar(c,j) * w_scalar(c,j) * &
           depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_soil1) = k_s1 * t_scalar(c,j) * w_scalar(c,j) * &
           depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_soil2) = k_s2 * (1.5 ** ((t_soisno(c,j)-(tfrz+25.0))/ &
           10._rk8)) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt 
         decomp_k(c,j,i_soil3) = k_s3 * (2.0 ** ((t_soisno(c,j)-(tfrz+25.0))/ &
           10._rk8)) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_soil4) = k_s4 * (2.5 ** ((t_soisno(c,j)-(tfrz+25.0))/ &
           10._rk8)) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_bacteria) = ck_bacteria(fc) * t_scalar(c,j) * w_scalar(c,j) * &
           depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_fungi) = ck_fungi(fc) * t_scalar(c,j) * w_scalar(c,j) * &
           depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_dom) = ck_dom(fc) * (1.25 ** &
           ((t_soisno(c,j)-(tfrz+25.0))/10._rk8))* w_scalar(c,j) * &
           depth_scalar(c,j) * o_scalar(c,j) / dt
!write(*,*) "scarlar:",o_scalar(fc,j),w_scalar(c,j),t_scalar(c,j),depth_scalar(c,j),c,j
#else
         decomp_k(c,j,i_litr1) = k_l1 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_litr2) = k_l2 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_litr3) = k_l3 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_cwd) = k_frag * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_soil1) = k_s1 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_soil2) = k_s2 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_soil3) = k_s3 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_soil4) = k_s4 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
#endif
      end do
   end do
#else
   do j = 1,nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         decomp_k(c,j,i_litr1) = k_l1 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_litr2) = k_l2 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_litr3) = k_l3 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_cwd) = k_frag * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_soil1) = k_s1 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_soil2) = k_s2 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_soil3) = k_s3 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_soil4) = k_s4 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
      end do
   end do

!!! zhangq modified for Microbe
#ifdef MICROBE
   !deallocate(decomp_depth_efolding_in)
   !deallocate(k_dom_in)
   !deallocate(k_bacteria_in)
   !deallocate(k_fungi_in)
   deallocate(ck_dom)
   deallocate(ck_bacteria)
   deallocate(ck_fungi)
#endif

#endif
  end subroutine decomp_rate_constants
#endif

#endif

end module mod_clm_cndecompcascadebgc
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
