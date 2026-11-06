module mod_clm_cndecomp
#ifdef CN
  !
  ! Module holding routines used in litter and soil decomposition model
  ! for coupled carbon-nitrogen code.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_clm_varcon , only : dzsoi_decomp, zsoi
#ifdef MICROBE
    use mod_clm_varcon, only: MBC_k
#endif
#ifndef CENTURY_DECOMP
  use mod_clm_cndecompcascadebgc , only : decomp_rate_constants
#else
  use mod_clm_cndecompcascadecentury , only : decomp_rate_constants
#endif
#ifdef NITRIF_DENITRIF
  use mod_clm_cnnitrifdenitrif , only : nitrif_denitrif
#endif
  use mod_clm_cnverticalprofile , only : decomp_vertprofiles

  implicit none

  save

  private

  public :: CNDecompAlloc

!!! zhangq add for Microbe
!#ifdef MICROBE
!   real(rk8), public :: decomp_depth_efolding = 0.17_rk8 !0.5    ! (meters) e-folding depth for reduction in decomposition [set to large number for
!   e-folding depth for reduction in decomposition [set to large number for depth-independance]
!#endif

  contains

  subroutine CNDecompAlloc(lbp, ubp, lbc, ubc, num_soilc, filter_soilc, &
                           num_soilp, filter_soilp)
    use mod_clm_type
    use mod_clm_cnallocation , only : CNAllocation
    use mod_clm_varpar , only : nlevdecomp , &
      ndecomp_cascade_transitions , ndecomp_pools
    use mod_clm_subgridave , only : p2c

!!! zhangq add for Microbe
#ifdef MICROBE
   use mod_runparams , only : dtsrf
   use mod_clm_microbevarcon, only: plant2doc, som2doc, micbiocn, micbioMR, CUEref, CUEt, Tcueref 
   use mod_clm_microbevarcon, only: Tsref, Tmref, Tsmin, Tmmin, Msmin, Msmax, Mmmin, dock 
   use mod_clm_cndecompcascadebgc
   use mod_clm_varcon, only: secspday   
   use mod_clm_varpar, only: i_bacteria, i_fungi, i_dom, cn_dom, numpft
   use mod_clm_varpar, only: max_pft_per_col
#endif

    implicit none
    integer(ik4), intent(in) :: lbp, ubp    ! pft-index bounds
    integer(ik4), intent(in) :: lbc, ubc    ! column-index bounds
    integer(ik4), intent(in) :: num_soilc   ! number of soil columns in filter
    ! filter for soil columns
    integer(ik4), intent(in) :: filter_soilc(ubc-lbc+1)
    integer(ik4), intent(in) :: num_soilp   ! number of soil pfts in filter
    integer(ik4), intent(in) :: filter_soilp(ubp-lbp+1) ! filter for soil pfts

    ! all c pools involved in decomposition
    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    real(rk8), pointer :: decomp_cpools_vr(:,:,:)
    ! all n pools involved in decomposition
    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
    real(rk8), pointer :: decomp_npools_vr(:,:,:)

    ! index into landunit level quantities
    integer(ik4), pointer :: clandunit(:)
    ! landunit type
    integer(ik4) , pointer :: itypelun(:)
    ! pft level
    ! fraction of roots in each soil layer  (nlevgrnd)
    real(rk8), pointer :: rootfr(:,:)

    ! fraction of potential immobilization (no units)
    real(rk8), pointer :: fpi_vr(:,:)
    ! vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
    real(rk8), pointer :: decomp_cascade_hr_vr(:,:,:)
    ! vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
    real(rk8), pointer :: decomp_cascade_ctransfer_vr(:,:,:)
    ! vert-res transfer of N from donor to receiver pool along
    ! decomp. cascade (gN/m3/s)
    real(rk8), pointer :: decomp_cascade_ntransfer_vr(:,:,:)
    ! vert-res mineral N flux for transition along
    ! decomposition cascade (gN/m3/s)
    real(rk8), pointer :: decomp_cascade_sminn_flux_vr(:,:,:)

    real(rk8), pointer :: potential_immob_vr(:,:)
#ifndef NITRIF_DENITRIF
    real(rk8), pointer :: sminn_to_denit_decomp_cascade_vr(:,:,:)
    real(rk8), pointer :: sminn_to_denit_excess_vr(:,:)
#endif
    real(rk8), pointer :: gross_nmin_vr(:,:)
    real(rk8), pointer :: net_nmin_vr(:,:)
    ! gross rate of N mineralization (gN/m2/s)
    real(rk8), pointer :: gross_nmin(:)
    ! net rate of N mineralization (gN/m2/s)
    real(rk8), pointer :: net_nmin(:)
    ! For methane code
#ifdef LCH4
    ! fraction of potential SOM + LITTER heterotrophic respiration
    real(rk8), pointer :: fphr(:,:)
    ! fraction by which decomposition is limited by moisture availability
    real(rk8), pointer :: w_scalar(:,:)
#endif

!!! zhangq add for Microbe
#ifdef MICROBE
!  real(rk8), pointer :: fphr(:,:)               ! fraction of potential SOM + LITTER heterotrophic respiration
   real(rk8), pointer :: w_scalar(:,:)            ! fraction by which decomposition is limited by moisture availability
   real(rk8), pointer :: t_scalar(:,:)            ! fraction by which decomposition is limited by temperature
   real(rk8), pointer :: o_scalar(:,:)            ! fraction by which decomposition is limited by oxygen
!  real(rk8), pointer :: depth_scalar(:,:)       ! fraction by which decomposition is reduced along soil profile
   real(rk8), pointer :: cdocs_pre(:,:)           ! the available carbon for methane processes
   real(rk8), pointer :: cdocs(:,:)               ! the available carbon for methane processes
   real(rk8), pointer :: cdocs_unsat(:,:)         ! the available carbon for methane processes in unsaturated fraction
   real(rk8), pointer :: cdocs_sat(:,:)           ! the available carbon for methane processes in saturated fraction
   real(rk8), pointer :: cmicbiocs(:,:)           ! microbial biomass carbon at columne-level
   real(rk8), pointer :: cdons(:,:)               ! availabe nitrogen
   real(rk8), pointer :: cdons_unsat(:,:)         ! nitrogen
   real(rk8), pointer :: cdons_sat(:,:)           ! nitrogen
   real(rk8), pointer :: cdons_min(:,:)           ! nitrogen
   real(rk8), pointer :: cmicbions(:,:)           ! microbial biomass nitrogen at columne-level
   real(rk8), pointer :: h2osoi_vol(:,:)          ! water volume for each column
   real(rk8), pointer :: watsat(:,:)              ! the available carbon for methane processes
   real(rk8), pointer :: t_soisno(:,:)            ! the available carbon for methane processes   
   real(rk8), pointer :: dochr_vr(:,:)            ! vertically resolved decomposition of dissolved organic carbon 
   real(rk8), pointer :: micbio_hr_vr(:,:,:)      ! vertically resolved decomposition of dissolved organic carbon 
   real(rk8), pointer :: micbio_hr(:,:)           ! vertically resolved decomposition of dissolved organic carbon 
   real(rk8), pointer :: finundated(:)            !inundated fractin of soil column
   real(rk8), Pointer :: froot_r(:,:)
   real(rk8), Pointer :: root2doc(:,:)            ! part of fine root exudate carbon to soil to develop dissolved organic carbon
                                                 ! we set it as part of fine root respiration
   real(rk8), Pointer :: cn_microbe(:,:)          ! C:N ration of the microbe as a combination of bacteria and fungi
   logical, pointer :: is_microbe(:)             ! TRUE => pool is a microbe pool

   real(rk8),allocatable :: decomp_depth_efolding_in(:) ! define the active decomposing soil depth
   real(rk8), pointer :: decomp_depth_efolding(:)       ! define the active decomposing soil depth

   integer(ik4), pointer :: ivt(:)                   ! pft vegetation type
   integer(ik4):: pi, p

   real(rk8), pointer :: wtcol(:)                ! pft weight relative to column (0-1)
   integer(ik4), allocatable :: pft_index(:)
   integer(ik4), pointer :: pfti(:)                 ! pft index array
   integer(ik4), pointer :: npfts(:)                ! number of pfts on the column
#endif
!!! zhangq midified for Microbe
   logical, pointer :: is_litter(:)                         ! TRUE => pool is a litter pool
   logical, pointer :: is_soil(:)                           ! TRUE => pool is a soil pool

    ! rate constant for decomposition (1./sec)
    real(rk8), pointer :: decomp_k(:,:,:)
    ! respired fraction in decomposition step (frac)
    real(rk8), pointer :: rf_decomp_cascade(:,:,:)
    ! which pool is C taken from for a given decomposition step
    integer(ik4),  pointer :: cascade_donor_pool(:)
    ! which pool is C added to for a given decomposition step
    integer(ik4),  pointer :: cascade_receiver_pool(:)
    ! what fraction of C leaving a given pool passes through
    ! a given transition (frac)
    real(rk8), pointer :: pathfrac_decomp_cascade(:,:,:)
    ! TRUE => pool has fixed C:N ratio
    logical,  pointer :: floating_cn_ratio_decomp_pools(:)

    integer(ik4) :: c,j,k,l    !indices
    integer(ik4) :: fc           !lake filter column index
    !potential C loss from one pool to another
    real(rk8):: p_decomp_cpool_loss(lbc:ubc,1:nlevdecomp, &
                                            1:ndecomp_cascade_transitions)
    !potential mineral N flux, from one pool to another
    real(rk8):: pmnf_decomp_cascade(lbc:ubc,1:nlevdecomp, &
                                            1:ndecomp_cascade_transitions)
    !potential N immobilization
    real(rk8):: immob(lbc:ubc,1:nlevdecomp)
    real(rk8):: ratio        !temporary variable
    real(rk8):: dnp          !denitrification proportion
    real(rk8):: cn_decomp_pools(lbc:ubc,1:nlevdecomp,1:ndecomp_pools)
    ! c:n ratio for initialization of pools
    real(rk8), pointer :: initial_cn_ratio(:)
    integer(ik4), parameter :: i_atm = 0
    ! maximum annual depth of thaw
    integer(ik4), pointer :: altmax_indx(:)
    ! prior year maximum annual depth of thaw
    integer(ik4), pointer :: altmax_lastyear_indx(:)

!!! zhangq add for Microbe
#ifdef MICROBE
   real(rk8) :: CUE, fm_t,  fm_m
   real(rk8) :: depth_scalar(lbc:ubc,1:nlevdecomp) 
   real(rk8) :: microbeMR
   real(rk8) :: doc_k
   integer :: dt, dtd
#endif

    ! For methane code
#ifndef NITRIF_DENITRIF
    real(rk8):: phr_vr(lbc:ubc,1:nlevdecomp) !potential HR (gC/m3/s)
#else
    real(rk8), pointer :: phr_vr(:,:)        !potential HR (gC/m3/s)
#endif
#ifdef LCH4
    real(rk8):: hrsum(lbc:ubc,1:nlevdecomp)  !sum of HR (gC/m2/s)
#endif

    decomp_cpools_vr             => clm3%g%l%c%ccs%decomp_cpools_vr
    decomp_cascade_hr_vr         => clm3%g%l%c%ccf%decomp_cascade_hr_vr
    decomp_cascade_ctransfer_vr  => clm3%g%l%c%ccf%decomp_cascade_ctransfer_vr
    decomp_npools_vr             => clm3%g%l%c%cns%decomp_npools_vr
    decomp_cascade_ntransfer_vr  => clm3%g%l%c%cnf%decomp_cascade_ntransfer_vr
    decomp_cascade_sminn_flux_vr => clm3%g%l%c%cnf%decomp_cascade_sminn_flux_vr
    fpi_vr                => clm3%g%l%c%cps%fpi_vr
    potential_immob_vr    => clm3%g%l%c%cnf%potential_immob_vr

    decomp_k                   => clm3%g%l%c%ccf%decomp_k
    rf_decomp_cascade          => clm3%g%l%c%cps%rf_decomp_cascade
    cascade_donor_pool         => decomp_cascade_con%cascade_donor_pool
    cascade_receiver_pool      => decomp_cascade_con%cascade_receiver_pool
    pathfrac_decomp_cascade    => clm3%g%l%c%cps%pathfrac_decomp_cascade
    floating_cn_ratio_decomp_pools => &
      decomp_cascade_con%floating_cn_ratio_decomp_pools
    initial_cn_ratio       => decomp_cascade_con%initial_cn_ratio
    altmax_indx            => clm3%g%l%c%cps%altmax_indx
    altmax_lastyear_indx   => clm3%g%l%c%cps%altmax_lastyear_indx

#ifndef NITRIF_DENITRIF
    sminn_to_denit_decomp_cascade_vr => &
      clm3%g%l%c%cnf%sminn_to_denit_decomp_cascade_vr
    sminn_to_denit_excess_vr => clm3%g%l%c%cnf%sminn_to_denit_excess_vr
#else
    phr_vr           => clm3%g%l%c%ccf%phr_vr
#endif
    gross_nmin_vr    => clm3%g%l%c%cnf%gross_nmin_vr
    net_nmin_vr      => clm3%g%l%c%cnf%net_nmin_vr
    gross_nmin       => clm3%g%l%c%cnf%gross_nmin
    net_nmin         => clm3%g%l%c%cnf%net_nmin
    ! For methane code
#ifdef LCH4
    fphr             => clm3%g%l%c%cch4%fphr
    w_scalar         => clm3%g%l%c%ccf%w_scalar
#endif

    rootfr           => clm3%g%l%c%p%pps%rootfr
    clandunit        => clm3%g%l%c%landunit
    itypelun         => clm3%g%l%itype

!!! zhangq add for Microbe
#ifdef MICROBE
   w_scalar                   => clm3%g%l%c%ccf%w_scalar
   t_scalar                   => clm3%g%l%c%ccf%t_scalar
   o_scalar                   => clm3%g%l%c%ccf%o_scalar
!   depth_scalar              => clm3%g%l%c%ccf%depth_scalar
   cdocs                      => clm3%g%l%c%cmic%cdocs
   cdocs_unsat                => clm3%g%l%c%cmic%cdocs_unsat
   cdocs_sat                  => clm3%g%l%c%cmic%cdocs_sat
   cmicbiocs                  => clm3%g%l%c%cmic%cmicbiocs
   cdons                      => clm3%g%l%c%cmic%cdons
   cdons_unsat                => clm3%g%l%c%cmic%cdons_unsat
   cdons_sat                  => clm3%g%l%c%cmic%cdons_sat
   cdons_min                  => clm3%g%l%c%cmic%cdons_min
   cmicbions                  => clm3%g%l%c%cmic%cmicbions
   h2osoi_vol                 => clm3%g%l%c%cws%h2osoi_vol
   watsat                     => clm3%g%l%c%cps%watsat
   t_soisno                   => clm3%g%l%c%ces%t_soisno
   dochr_vr                   => clm3%g%l%c%cmic%dochr_vr
   micbio_hr_vr               => clm3%g%l%c%cmic%micbio_hr_vr
   micbio_hr                  => clm3%g%l%c%cmic%micbio_hr
   root2doc                   => clm3%g%l%c%cmic%root2doc
   froot_r                    => clm3%g%l%c%cmic%froot_r
   cn_microbe                 => clm3%g%l%c%cmic%cn_microbe
   finundated                 => clm3%g%l%c%cws%finundated   
   is_microbe                 => decomp_cascade_con%is_microbe

   decomp_depth_efolding      => pftcon%decomp_depth_efolding 
   ivt                        => clm3%g%l%c%p%itype
  
   wtcol                      => clm3%g%l%c%p%wtcol
   pfti                       => clm3%g%l%c%pfti
   npfts                      => clm3%g%l%c%npfts

  if (num_soilc .gt. 0) then
  !!! zhangq modified for Microbe
  !allocate(pft_index(0))
  allocate(pft_index(1:num_soilc))
  allocate(decomp_depth_efolding_in(1:num_soilc))  

  do fc = 1,num_soilc
  c = filter_soilc(fc)

  do pi = 1,max_pft_per_col

        if (pi <=  npfts(c)) then
           p = pfti(c) + pi - 1
        end if
 
  end do


    if (wtcol(pfti(c)) .ge. 0.96_rk8 .or. SUM(wtcol((pfti(c) + 1):p)) .lt. 0.001_rk8) then
       !pft_index = [pft_index, 0]
       pft_index(fc) = 0
    else
       !pft_index = [pft_index, MAXLOC(wtcol((pfti(c) + 1):p), DIM=1)]
       pft_index(fc) = MAXLOC(wtcol((pfti(c) + 1):p),DIM=1)
    end if

! pft_index = [pft_index, (MAXLOC(wtcol(pfti(c):p), DIM=1, mask =
! wtcol(pfti(c):p) .gt. 0) - 1)]
! write (*, *) "pfti(c) 1", pfti(c)
! write (*, *) "p 1", p

! write (*, *) "SUM(wtcol((pfti(c) + 1):p)", SUM(wtcol((pfti(c) + 1):p))
! write (*, *) "pft_index 3", pft_index

  end do
  
  end if

decomp_depth_efolding_in = decomp_depth_efolding(pft_index(:))

! write (*, *) "decomp_depth_efolding_in 3", decomp_depth_efolding_in

#endif  
!!! zhangq midified for Microbe
   is_litter                  => decomp_cascade_con%is_litter
   is_soil                    => decomp_cascade_con%is_soil

    call decomp_rate_constants(lbc, ubc, num_soilc, filter_soilc)

    ! set initial values for potential C and N fluxes
    p_decomp_cpool_loss(:,:,:) = 0._rk8
    pmnf_decomp_cascade(:,:,:) = 0._rk8

    ! column loop to calculate potential decomp rates and total immobilization
    ! demand.

    ! calculate c:n ratios of applicable pools
    do l = 1 , ndecomp_pools
      if ( floating_cn_ratio_decomp_pools(l) ) then
        do j = 1 , nlevdecomp
          do fc = 1 , num_soilc
            c = filter_soilc(fc)
            if ( decomp_npools_vr(c,j,l) > 0._rk8 ) then
              cn_decomp_pools(c,j,l) = decomp_cpools_vr(c,j,l) / &
                decomp_npools_vr(c,j,l)
            end if
          end do
        end do
      else
        do j = 1,nlevdecomp
          do fc = 1,num_soilc
            c = filter_soilc(fc)
            cn_decomp_pools(c,j,l) = initial_cn_ratio(l)
          end do
        end do
      end if
    end do

!!! zhangq add for Microbe
#ifdef MICROBE
      ! add a term to reduce decomposition rate at depth
      ! for now used a fixed e-folding depth
      do j = 1, nlevdecomp
        do fc = 1,num_soilc
          c = filter_soilc(fc)
      if((decomp_npools_vr(c,j,i_bacteria) + decomp_npools_vr(c,j,i_fungi)) < 1e-15) then
      cn_microbe(c,j) = 10._rk8
      else
      cn_microbe(c,j) = (decomp_cpools_vr(c,j,i_bacteria) + &
                        decomp_cpools_vr(c,j,i_fungi)) / (decomp_npools_vr(c,j,i_bacteria) + &
                        decomp_npools_vr(c,j,i_fungi))
      end if
            depth_scalar(c,j) = exp(-zsoi(j)/decomp_depth_efolding_in(fc))
            cdocs(c,j) = max(cdocs(c,j),1e-20)
            !write(*,*) "here debugging1: ", c,j, cdocs(c,j), cdons(c,j),
            !decomp_cpools_vr(c,j,i_dom),decomp_npools_vr(c,j,i_dom)      
            !decomp_cpools_vr(c,j,i_dom) = cdocs(c,j)
            !decomp_npools_vr(c,j,i_dom) = cdocs(c,j) / cn_dom
            !write(*,*) "here debugging: ", c,j, cdocs(c,j), cdons(c,j)
        end do
      end do

   dt = dtsrf
   dtd = dt/secspday
   
   !~ microbeMR = -log(1.0_rk8-micbioMR)   
   !~ doc_k = -log(1.0_rk8-dock)  
   !~ ! calculate the new discrete-time decay rate for model timestep
   !~ microbeMR = 1.0_rk8-exp(-microbeMR*dtd)   
   !~ doc_k = 1.0_rk8-exp(-doc_k*dtd)  
#endif   

    ! calculate the non-nitrogen-limited fluxes
    ! these fluxes include the  "/ dt" term to put them on a
    ! per second basis, since the rate constants have been
    ! calculated on a per timestep basis.

    do k = 1 , ndecomp_cascade_transitions
      do j = 1 , nlevdecomp
        do fc = 1 , num_soilc
          c = filter_soilc(fc)

          if ( decomp_cpools_vr(c,j,cascade_donor_pool(k)) > 0._rk8 .and. &
               decomp_k(c,j,cascade_donor_pool(k)) > 0._rk8 ) then
!!! zhangq add for Microbe
#ifdef MICROBE
            if ((is_litter(cascade_donor_pool(k)) .or. is_soil(cascade_donor_pool(k))) .and. &
               (is_microbe(cascade_donor_pool(k)) .eqv. .false.)) then
              p_decomp_cpool_loss(c,j,k) = decomp_cpools_vr(c,j,cascade_donor_pool(k)) * &
                       decomp_k(c,j,cascade_donor_pool(k))  * pathfrac_decomp_cascade(c,j,k) !*(cmicbiocs(c,j) + 1e-5) / (cmicbiocs(c,j) + MBC_k(cascade_donor_pool(k)))
            else
              p_decomp_cpool_loss(c,j,k) = decomp_cpools_vr(c,j,cascade_donor_pool(k)) * &
                       decomp_k(c,j,cascade_donor_pool(k))  * pathfrac_decomp_cascade(c,j,k)
            end if
#else
              p_decomp_cpool_loss(c,j,k) = decomp_cpools_vr(c,j,cascade_donor_pool(k)) * &
                       decomp_k(c,j,cascade_donor_pool(k))  * pathfrac_decomp_cascade(c,j,k)
#endif

           if ( .not. &
             floating_cn_ratio_decomp_pools(cascade_receiver_pool(k)) ) then
             !! not transition of cwd to litter

             if (cascade_receiver_pool(k) /= i_atm ) then
               ! not 100% respiration
               ratio = 0._rk8
               if ( decomp_npools_vr(c,j,cascade_donor_pool(k)) > 0._rk8 ) then
                 ratio = cn_decomp_pools(c,j,cascade_receiver_pool(k)) / &
                   cn_decomp_pools(c,j,cascade_donor_pool(k))
               endif
               pmnf_decomp_cascade(c,j,k) = (p_decomp_cpool_loss(c,j,k) * &
                 (1.0_rk8 - rf_decomp_cascade(c,j,k) - ratio) / &
                  cn_decomp_pools(c,j,cascade_receiver_pool(k)) )
             else   ! 100% respiration
               pmnf_decomp_cascade(c,j,k) = - p_decomp_cpool_loss(c,j,k) / &
                 cn_decomp_pools(c,j,cascade_donor_pool(k))
             end if
           else   ! CWD -> litter
             pmnf_decomp_cascade(c,j,k) = 0._rk8
           end if
         end if
       end do
     end do
   end do

   ! Sum up all the potential immobilization fluxes (positive pmnf flux)
   ! and all the mineralization fluxes (negative pmnf flux)
   do j = 1 , nlevdecomp
     do fc = 1 , num_soilc
       c = filter_soilc(fc)
       immob(c,j) = 0._rk8
     end do
   end do
   do k = 1 , ndecomp_cascade_transitions
     do j = 1 , nlevdecomp
       do fc = 1 , num_soilc
         c = filter_soilc(fc)
         if (pmnf_decomp_cascade(c,j,k) > 0._rk8) then
           immob(c,j) = immob(c,j) + pmnf_decomp_cascade(c,j,k)
         else
           gross_nmin_vr(c,j) = gross_nmin_vr(c,j) - pmnf_decomp_cascade(c,j,k)
         end if
       end do
     end do
   end do

   do j = 1 , nlevdecomp
     do fc = 1 , num_soilc
       c = filter_soilc(fc)
       potential_immob_vr(c,j) = immob(c,j)
     end do
   end do

   ! Add up potential hr for methane calculations
   do j = 1 , nlevdecomp
     do fc = 1 , num_soilc
       c = filter_soilc(fc)
       phr_vr(c,j) = 0._rk8
     end do
   end do
   do k = 1 , ndecomp_cascade_transitions
     do j = 1 , nlevdecomp
       do fc = 1 , num_soilc
         c = filter_soilc(fc)
         phr_vr(c,j) = phr_vr(c,j) + rf_decomp_cascade(c,j,k) * &
           p_decomp_cpool_loss(c,j,k)
       end do
     end do
   end do

   call decomp_vertprofiles(lbp, ubp, lbc,ubc,num_soilc, &
     filter_soilc,num_soilp,filter_soilp)

#ifdef NITRIF_DENITRIF
   ! calculate nitrification and denitrification rates
   call nitrif_denitrif(lbc, ubc, num_soilc, filter_soilc)
#endif

   ! now that potential N immobilization is known, call allocation
   ! to resolve the competition between plants and soil heterotrophs
   ! for available soil mineral N resource.

   call CNAllocation(lbp, ubp, lbc,ubc,num_soilc,filter_soilc,num_soilp, &
        filter_soilp)

   ! column loop to calculate actual immobilization and decomp rates, following
   ! resolution of plant/heterotroph  competition for mineral N

   dnp = 0.01_rk8

   ! calculate c:n ratios of applicable pools
   do l = 1 , ndecomp_pools
     if ( floating_cn_ratio_decomp_pools(l) ) then
       do j = 1 , nlevdecomp
         do fc = 1 , num_soilc
           c = filter_soilc(fc)
           if ( decomp_npools_vr(c,j,l) > 0._rk8 ) then
             cn_decomp_pools(c,j,l) = decomp_cpools_vr(c,j,l) / &
               decomp_npools_vr(c,j,l)
           end if
         end do
       end do
     else
       do j = 1 , nlevdecomp
         do fc = 1 , num_soilc
           c = filter_soilc(fc)
           cn_decomp_pools(c,j,l) = initial_cn_ratio(l)
         end do
       end do
     end if
   end do

   ! upon return from CNAllocation, the fraction of potential immobilization
   ! has been set (cps%fpi_vr). now finish the decomp calculations.
   ! Only the immobilization steps are limited by fpi_vr (pmnf > 0)
   ! Also calculate denitrification losses as a simple proportion
   ! of mineralization flux.

   do k = 1 , ndecomp_cascade_transitions
     do j = 1 , nlevdecomp
       do fc = 1 , num_soilc
         c = filter_soilc(fc)

         if ( decomp_cpools_vr(c,j,cascade_donor_pool(k)) > 0._rk8 ) then
           if ( pmnf_decomp_cascade(c,j,k) > 0._rk8 ) then
             p_decomp_cpool_loss(c,j,k) = &
               p_decomp_cpool_loss(c,j,k) * fpi_vr(c,j)
             pmnf_decomp_cascade(c,j,k) = &
               pmnf_decomp_cascade(c,j,k) * fpi_vr(c,j)
#ifndef NITRIF_DENITRIF
             sminn_to_denit_decomp_cascade_vr(c,j,k) = 0._rk8
           else
             sminn_to_denit_decomp_cascade_vr(c,j,k) = -dnp * &
               pmnf_decomp_cascade(c,j,k)
#endif
           end if
           decomp_cascade_hr_vr(c,j,k) = rf_decomp_cascade(c,j,k) * &
             p_decomp_cpool_loss(c,j,k)
           decomp_cascade_ctransfer_vr(c,j,k) = (1._rk8 - &
             rf_decomp_cascade(c,j,k)) * p_decomp_cpool_loss(c,j,k)
           if ( decomp_npools_vr(c,j,cascade_donor_pool(k)) > 0._rk8 .and. &
                cascade_receiver_pool(k) /= i_atm ) then
             decomp_cascade_ntransfer_vr(c,j,k) = p_decomp_cpool_loss(c,j,k) / &
               cn_decomp_pools(c,j,cascade_donor_pool(k))
           else
             decomp_cascade_ntransfer_vr(c,j,k) = 0._rk8
           endif
           if ( cascade_receiver_pool(k) /= 0 ) then
             decomp_cascade_sminn_flux_vr(c,j,k) = pmnf_decomp_cascade(c,j,k)
           else  ! keep sign convention negative for terminal pools
             decomp_cascade_sminn_flux_vr(c,j,k) = - pmnf_decomp_cascade(c,j,k)
           endif
           net_nmin_vr(c,j) = net_nmin_vr(c,j) - pmnf_decomp_cascade(c,j,k)
         else
           decomp_cascade_ntransfer_vr(c,j,k) = 0._rk8
#ifndef NITRIF_DENITRIF
           sminn_to_denit_decomp_cascade_vr(c,j,k) = 0._rk8
#endif
           decomp_cascade_sminn_flux_vr(c,j,k) = 0._rk8
         end if
!!! zhangq add for Microbe
!~ #ifdef MICROBE    
  !~ if(decomp_cpools_vr(c,j,cascade_donor_pool(k)) .gt. 0._rk8)
  !~ decomp_cascade_hr_vr(c,j,k) = rf_decomp_cascade(c,j,k) *
  !p_decomp_cpool_loss(c,j,k)
!~ #endif
       end do
     end do
   end do

!!! zhangq add for Microbe
#ifdef MICROBE
!   do l = 1, ndecomp_pools
       do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
      !~ decomp_cpools_vr(c,j,i_dom) = decomp_cpools_vr(c,j,i_dom) !- cdons_min(c,j) * cn_dom
!write(iulog, *) "xiaofeng here1 ", cdocs(c,j), decomp_cpools_vr(c,j,i_dom)  
!write(*,*) "here debugging2: ", c,j, cdocs(c,j), cdons(c,j), decomp_cpools_vr(c,j,i_dom),decomp_npools_vr(c,j,i_dom)      
      !~ cdocs(c,j) = decomp_cpools_vr(c,j,i_dom)   ! gC/m3
!write(iulog, *) "xiaofeng here2 ", cdocs(c,j), decomp_cpools_vr(c,j,i_dom)     
      !~ cdocs_unsat(c,j) = cdocs(c,j) * (1. - finundated(c))
      !~ cdocs_sat(c,j) = cdocs(c,j) * finundated(c)
      
      !~ decomp_npools_vr(c,j,i_dom) = decomp_npools_vr(c,j,i_dom) !- cdons_min(c,j)
      !~ cdons(c,j) = decomp_npools_vr(c,j,i_dom)
      !~ cdons_unsat(c,j) = cdons(c,j) * (1. - finundated(c))
      !~ cdons_sat(c,j) = cdons(c,j) * finundated(c)
      
      cmicbiocs(c,j) = decomp_cpools_vr(c,j,i_bacteria) + decomp_cpools_vr(c,j,i_fungi)
      cmicbions(c,j) = decomp_npools_vr(c,j,i_bacteria) + decomp_npools_vr(c,j,i_fungi)
           end do
      end do
!end do
   !~ do j = 1,nlevsoi
      !~ do fc = 1,num_soilc
         !~ c = filter_soilc(fc)
!~ !write(iulog, *) "here2", cdocs(c,j), " doc_k ", doc_k
  !~ cdocs(c,j) = cdocs_unsat(c,j) * (1. - finundated(c)) + cdocs(c,j) * finundated(c)
        !~ dochr_vr(c,j) = cdocs(c,j) * doc_k * 1.5**((t_soisno(c,j) - 273.15 - 25.) / 10) ! * t_scalar(c,j) * w_scalar(c,j)
  !~ cdocs_unsat(c,j) = cdocs_unsat(c,j) * (1. - finundated(c)) * (1. - doc_k * 1.5**((t_soisno(c,j) - 273.15 - 25.) / 10)) ! * t_scalar(c,j) * w_scalar(c,j)
  !~ cdocs_sat(c,j) = cdocs_sat(c,j) * finundated(c) * (1.0 - doc_k * 1.5**((t_soisno(c,j) - 273.15 - 25.) / 10)) ! * t_scalar(c,j) * w_scalar(c,j)
  !~ cdocs_unsat(c,j) = max(0._rk8, cdocs_unsat(c,j))
  !~ cdocs_sat(c,j) = max(0._rk8, cdocs_sat(c,j)) cdocs(c,j) = cdocs_unsat(c,j) * (1. - finundated(c)) + cdocs(c,j) * finundated(c)
  !~ cdocs(c,j) = cdocs_unsat(c,j) * (1. - finundated(c)) + cdocs(c,j) * finundated(c)
!~ !write(iulog, *) dochr_vr(c,j), " dochr doc ", cdocs(c,j)
!~ !  calculation of death of the microbes

!~ if(t_soisno(c,j) < Tmmin) then
!~ fm_t = (2. ** ((Tmmin - Tmref)/10)) ** 2
!~ else
!~ fm_t = 2. ** ((t_soisno(c,j) - Tmref)/10)
!~ endif

!~ if((h2osoi_vol(c,j) / watsat(c,j)) < Mmmin) then
!~ fm_m = 0._rk8
!~ else
!~ fm_m = 1._rk8
!~ endif

!~ micbiohr_vr(c,j) = cmicbios(c,j) * microbeMR / dt * fm_t * fm_m *
!depth_scalar(c,j) * o_scalar(c,j)
!~ if(cmicbios(c,j) - micbiohr_vr(c,j) > 1e-15) then
!~ cmicbios(c,j) =  cmicbios(c,j) - micbiohr_vr(c,j)
!~ else
!~ micbiohr_vr(c,j) = cmicbios(c,j) - 1e-15
!~ cmicbios(c,j) =  15e-15
!~ endif
      !~ end do
   !~ end do
#endif   

#ifdef LCH4
   ! Calculate total fraction of potential HR, for methane code
   do j = 1,nlevdecomp
     do fc = 1,num_soilc
       c = filter_soilc(fc)
       hrsum(c,j) = 0._rk8
     end do
   end do
   do k = 1 , ndecomp_cascade_transitions
     do j = 1 , nlevdecomp
       do fc = 1 , num_soilc
         c = filter_soilc(fc)
         hrsum(c,j) = hrsum(c,j) + &
           rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(c,j,k)
       end do
     end do
   end do

   ! Nitrogen limitation / (low)-moisture limitation
   do j = 1 , nlevdecomp
     do fc = 1 , num_soilc
       c = filter_soilc(fc)
       if (phr_vr(c,j) > 0._rk8) then
         fphr(c,j) = hrsum(c,j) / phr_vr(c,j) * w_scalar(c,j)
         ! Prevent overflow errors for 0 respiration
         fphr(c,j) = max(fphr(c,j), 0.01_rk8)
       else
         fphr(c,j) = 1._rk8
       end if
     end do
   end do
#endif

   ! vertically integrate net and gross mineralization fluxes
   ! for diagnostic output
   do j = 1 , nlevdecomp
     do fc = 1 , num_soilc
       c = filter_soilc(fc)
       net_nmin(c) = net_nmin(c) + net_nmin_vr(c,j) * dzsoi_decomp(j)
       gross_nmin(c) = gross_nmin(c) + gross_nmin_vr(c,j) * dzsoi_decomp(j)
     end do
   end do
 end subroutine CNDecompAlloc

#endif

end module mod_clm_cndecomp
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
