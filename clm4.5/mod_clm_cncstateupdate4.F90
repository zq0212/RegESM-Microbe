module mod_clm_cncstateupdate4
#ifdef MICROBE

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CStateUpdate3Mod
!
! !DESCRIPTION:
! Module for carbon state variable update, mortality fluxes.
!
! !USES:
    use mod_intkinds
    use mod_realkinds
    use mod_runparams , only : dtsrf
    use mod_mpmessage
    use mod_clm_varctl, only : use_c13, use_c14
    use mod_clm_varcon  , only: istcrop, c13ratio, c14ratio
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public:: CStateUpdate4, vert_diffusion
!
! !REVISION HISTORY:
! Xiaofeng Xu (4/1/2020) created for methane module
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CStateUpdate4
!
! !INTERFACE:
subroutine CStateUpdate4(num_soilc, filter_soilc, num_soilp, filter_soilp, isotope)
!
! !DESCRIPTION:
! On the radiation time step, update all the prognostic carbon state
! variables affected by fire fluxes
!
! !USES:

! !USES:
        use mod_clm_type
        use mod_clm_subgridave , only : p2c, c2g
        use mod_clm_varpar, only : nlevgrnd, nlevdecomp, i_bacteria, i_fungi, i_dom, cn_dom
        use mod_clm_pftvarcon,  only : noveg
        use mod_clm_varcon, only : secspday, istwet, istsoil, istdlak, spval, istcrop, tfrz, rgas
        use mod_clm_microbevarcon 
        use mod_clm_filter,  only : filter

! !ARGUMENTS:
   implicit none
   integer(ik4), intent(in) :: num_soilc       ! number of soil columns in filter
   integer(ik4), intent(in) :: filter_soilc(:) ! filter for soil columns
   integer(ik4), intent(in) :: num_soilp       ! number of soil pfts in filter
   integer(ik4), intent(in) :: filter_soilp(:) ! filter for soil pfts
   character(len=*), intent(in) :: isotope         ! 'bulk', 'c13' or 'c14'
!
! !LOCAL VARIABLES:
! local pointers to implicit in variables
	real(rk8), pointer :: decomp_cpools_vr(:,:,:) 		! bring DOMC to methane module
	real(rk8), pointer :: cdocs(:,:)					! column-level concentration of DOC molC/m3 
	
	real(rk8), pointer :: cmicbiocs(:,:)				! column-level biomass of all microbes molC/m3
	real(rk8), pointer :: caces(:,:)					! column-level concentration of acetate molC/m3
	real(rk8), pointer :: cacebios(:,:)				! column-level biomass of methanogen based on acetate molC/m3
	real(rk8), pointer :: cco2bios(:,:)				! column-level biomass of methanogen based on CO2/H2 molC/m3
	real(rk8), pointer :: caerch4bios(:,:)			! column-level biomass of aerobix methanotrophy molC/m3
	real(rk8), pointer :: canaerch4bios(:,:)			! column-level biomass of anaerobic methanotrophy molC/m3

	real(rk8), pointer :: caces_prod(:,:)				! acetate acid production
	real(rk8), pointer :: caces_prod_h2(:,:)			! acetogenesis
	
	real(rk8), pointer :: ccon_ch4s(:,:)				! concentration of CH4 mol C/m3 at each layer
	real(rk8), pointer :: ccon_co2s(:,:)				! column-level concentration of CO2 mol C/m3
      
	real(rk8), pointer :: ch4_prod_ace_depth(:,:)		! column-level aceteclastic methanogenesis
	real(rk8), pointer :: ch4_prod_co2_depth(:,:)		! column-level hydrogenotrophic methanogenesis
	real(rk8), pointer :: ch4_oxid_o2_depth(:,:)		! column-level oxic methanotrophy
	real(rk8), pointer :: ch4_oxid_aom_depth(:,:)		! column-level anaerobic oxidation of methane
	
	real(rk8), pointer :: ch4_aere_depth(:,:)			! column-level 
	real(rk8), pointer :: ch4_dif_depth(:,:)
	real(rk8), pointer :: ch4_ebul_depth(:,:)
	
	real(rk8), pointer :: co2_prod_ace_depth(:,:)
	real(rk8), pointer :: co2_decomp_depth(:,:)
	real(rk8), pointer :: co2_cons_depth(:,:)
	
	real(rk8), pointer :: co2_ebul_depth(:,:)
	real(rk8), pointer :: co2_aere_depth(:,:)
	real(rk8), pointer :: co2_dif_depth(:,:)
	
	real(rk8), pointer :: ch4_surf_aere(:)
	real(rk8), pointer :: ch4_surf_ebul(:)
	real(rk8), pointer :: ch4_surf_dif(:)
	real(rk8), pointer :: ch4_surf_netflux(:)
	
	real(rk8), pointer :: co2_surf_aere(:)
	real(rk8), pointer :: co2_surf_ebul(:)
	real(rk8), pointer :: co2_surf_dif(:)
	real(rk8), pointer :: co2_surf_netflux(:)

! !OTHER LOCAL VARIABLES:
   type(column_microbe_type), pointer :: cmiciso
    type(column_cstate_type), pointer :: ccisos
   integer(ik4) :: c,p,j,l,k,fc      ! indices
   real(rk8):: dt       ! radiation time step (seconds)

!EOP
!-----------------------------------------------------------------------
   ! select which isotope
   select case (isotope)
   case ('bulk')
      ccisos =>  clm3%g%l%c%ccs
      cmiciso =>  clm3%g%l%c%cmic
   case ('c14')
      ccisos =>  clm3%g%l%c%cc14s
      cmiciso =>  clm3%g%l%c%cmicc14
   case ('c13')
      ccisos =>  clm3%g%l%c%cc13s
      cmiciso =>  clm3%g%l%c%cmicc13
   case default
      call fatal(__FILE__,__LINE__,'CNCIsoStateUpdate3Mod: iso must be bulk, c13 or c14')
   end select

    
! !ARGUMENTS:

   ! Column level pointers
	decomp_cpools_vr			=> ccisos%decomp_cpools_vr
	cdocs 					=> cmiciso%cdocs
	
	cmicbiocs 					=> cmiciso%cmicbiocs
	caces 					=> cmiciso%caces
	caces_prod 				=> cmiciso%caces_prod
	caces_prod_h2 				=> cmiciso%caces_prod_h2

	cacebios 					=> cmiciso%cacebios
	cco2bios 					=> cmiciso%cco2bios
	caerch4bios 				=> cmiciso%caerch4bios
	canaerch4bios 				=> cmiciso%canaerch4bios
     
	ccon_ch4s        				=> cmiciso%ccon_ch4s
	ccon_co2s         			=> cmiciso%ccon_co2s
 
	ch4_prod_ace_depth 			=> cmiciso%ch4_prod_ace_depth
	ch4_prod_co2_depth 			=> cmiciso%ch4_prod_co2_depth
	ch4_oxid_o2_depth 			=> cmiciso%ch4_oxid_o2_depth
	ch4_oxid_aom_depth 		=> cmiciso%ch4_oxid_aom_depth
	
	ch4_aere_depth 			=> cmiciso%ch4_aere_depth
	ch4_dif_depth 				=> cmiciso%ch4_dif_depth
	ch4_ebul_depth 			=> cmiciso%ch4_ebul_depth
	
	co2_prod_ace_depth 			=> cmiciso%co2_prod_ace_depth
	co2_decomp_depth 			=> cmiciso%co2_decomp_depth
	co2_cons_depth 			=> cmiciso%co2_cons_depth
	
	co2_ebul_depth 			=> cmiciso%co2_ebul_depth
	co2_aere_depth 			=> cmiciso%co2_aere_depth
	co2_dif_depth 				=> cmiciso%co2_dif_depth
	
	ch4_surf_aere				=> cmiciso%ch4_surf_aere
	ch4_surf_ebul 				=> cmiciso%ch4_surf_ebul
	ch4_surf_dif 				=> cmiciso%ch4_surf_dif
	ch4_surf_netflux 			=> cmiciso%ch4_surf_netflux
	
	co2_surf_aere				=> cmiciso%co2_surf_aere
	co2_surf_ebul 				=> cmiciso%co2_surf_ebul
	co2_surf_dif 				=> cmiciso%co2_surf_dif
	co2_surf_netflux 			=> cmiciso%co2_surf_netflux


    ! set time steps
    dt = dtsrf
    
	do j = 1,nlevdecomp
	do fc = 1,num_soilc
		c = filter_soilc(fc)
		cdocs(c,j) = decomp_cpools_vr(c,j,i_dom) * 0.0112372
		caces(c,j) = caces(c,j) * 0.0112372
		cacebios(c,j) = cacebios(c,j) * 0.0112372
		cco2bios(c,j) = cco2bios(c,j) * 0.0112372
		caerch4bios(c,j) = caerch4bios(c,j) * 0.0112372
		canaerch4bios(c,j) = canaerch4bios(c,j) * 0.0112372
		ccon_ch4s(c,j) = ccon_ch4s(c,j) * 0.0112372
		ccon_co2s(c,j) = ccon_co2s(c,j) * 0.0112372


		
		!call vert_diffusion(1, num_soilc, num_soilc, filter_soilc, isotope)
		cdocs(c,j) = max(0._rk8, cdocs(c,j) - caces_prod(c,j))

!		xiaofeng commented this line 5/9/2021
!		decomp_cpools_vr(c,j,i_dom) = cdocs(c,j)

		caces(c,j) = max(0._rk8, (caces(c,j) + caces_prod(c,j)*dt + caces_prod_h2(c,j)*dt - ch4_prod_ace_depth(c,j)*dt))
!		cacebios(c,j) = cacebios(c,j)
!		cco2bios(c,j) = cco2bios(c,j)
!		caerch4bios(c,j) = caerch4bios(c,j)
!		canaerch4bios(c,j) = canaerch4bios(c,j)
		ccon_ch4s(c,j)=max(0._rk8, (ccon_ch4s(c,j) + ch4_prod_ace_depth(c,j)*dt + ch4_prod_co2_depth(c,j)*dt - &
		ch4_oxid_o2_depth(c,j)*dt - ch4_oxid_aom_depth(c,j)*dt - &
		ch4_aere_depth(c,j)*dt - ch4_dif_depth(c,j)*dt - ch4_ebul_depth(c,j)*dt))
	end do
	end do

end subroutine CStateUpdate4
!-----------------------------------------------------------------------

subroutine vert_diffusion(lbc, ubc,num_micbioc, filter_micbioc, isotope)
! The gas diffusion was simulated based on Fick's law; the concentration gradient is the key driver for gas diffusion along soil profile
! The gas diffusion occur in saturated layers, the unsaturated layers was not considered
! !USES:
	use mod_clm_type
	use mod_clm_varcon, only : secspday, istwet, istsoil, istdlak, spval, istcrop, Fick_D_w
	use mod_clm_microbevarcon
	use mod_clm_varpar, only : nlevgrnd, nlevdecomp
!
! !ARGUMENTS:
implicit none
	integer(ik4), intent(in) :: lbc, ubc                		! column bounds
	integer(ik4), intent(in) :: num_micbioc               		! number of soil columns in filter
	integer(ik4), intent(in) :: filter_micbioc(ubc-lbc+1) 	! filter for soil columns
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
	real(rk8), pointer :: z(:,:)         				! layer depth (m) (-nlevsno+1:nlevsoi)
	real(rk8), pointer :: dz(:,:)         				! layer thickness (m) (-nlevsno+1:nlevgrnd)

	real(rk8), pointer :: ccon_ch4s(:,:)			! column-level concentration of CH4 in saturated fraction
	real(rk8), pointer :: ccon_o2s(:,:)			! column-level concentration of O2 in saturated fraction
	real(rk8), pointer :: ccon_co2s(:,:)			! column-level concentration of CO2 in saturated fraction
	real(rk8), pointer :: ccon_h2s(:,:)			! column-level concentration of H2 in saturated fraction
	real(rk8), pointer :: t_soisno(:,:)
		
	real(rk8) :: ccon_ch4s_temp(lbc:ubc,1:nlevgrnd)		! temporary array 
	real(rk8) :: ccon_o2s_temp(lbc:ubc,1:nlevgrnd) 		! temporary array 
	real(rk8) :: ccon_co2s_temp(lbc:ubc,1:nlevgrnd)	! temporary array 
	real(rk8) :: ccon_h2s_temp(lbc:ubc,1:nlevgrnd)	! temporary array 
	real(rk8) :: temp
	integer(ik4) :: jwaterhead_unsat(lbc:ubc)			! layer of the water head in unsaturated fraction

	integer(ik4), pointer :: ltype(:)
	integer(ik4), pointer :: clandunit(:)
!
! local pointers to implicit in/out scalars
! !OTHER LOCAL VARIABLES:
	integer(ik4) :: c, j, l      ! indices
	integer(ik4) :: fc        ! soil column filter indices
	real(rk8) :: dlt
	integer(ik4) :: nlevsoi
	
	character(len=*), intent(in) :: isotope  ! 'bulk', 'c13' or 'c14'
!EOP
!-----------------------------------------------------------------------
! assign local pointers to derived type arrays
type(column_microbe_type), pointer :: ccisos
   select case (isotope)
   case ('c14')
!      pcisof =>  pc14f
!      pcisos =>  pc14s
!      ccisof =>  cc14f
!      ccisos =>  cc14s
     ccisos => clm3%g%l%c%cmicc14
   case ('c13')
!      pcisof =>  pc13f
!      pcisos =>  pc13s
!      ccisof =>  cc13f
!      ccisos =>  cc13s
      ccisos => clm3%g%l%c%cmicc13
    case ('bulk')
!      pcisof =>  pc13f
!      pcisos =>  pc13s
!      ccisof =>  cc13f
!      ccisos =>  cc13s
      ccisos => clm3%g%l%c%cmic
         case default
      call fatal(__FILE__,__LINE__,'CNCIsoStateUpdate4Mod: iso must be bulk, c13 or c14')
   end select
   
	z                      				=> clm3%g%l%c%cps%z
	dz                      			=> clm3%g%l%c%cps%dz

	ccon_ch4s        			=> ccisos%ccon_ch4s
	ccon_co2s        			=> ccisos%ccon_co2s
	ccon_o2s          			=> ccisos%ccon_o2s
	ccon_h2s          			=> ccisos%ccon_h2s
	
	ltype               			=> clm3%g%l%c%p%itype
	clandunit       			=> clm3%g%l%c%landunit
	t_soisno					=> clm3%g%l%c%ces%t_soisno
	
	nlevsoi = 10
	ccon_ch4s_temp 			= 0._rk8
	ccon_o2s_temp 			= 0._rk8
	ccon_co2s_temp 			= 0._rk8
	ccon_h2s_temp 			= 0._rk8

!	call get_waterhead(lbc, ubc, num_micbioc, filter_micbioc,jwaterhead_unsat)
	
	do fc = 1,num_micbioc
        c = filter_micbioc(fc)
	l = clandunit(c)     
	if (ltype(l) == istsoil .or. ltype(l) == istcrop) then 
	do j = 2,nlevsoi
	if(j > jwaterhead_unsat(c) .and. j < nlevsoi) then
	!write(iulog,*) "above water table"
	end if
	! for the saturation portion
	!write(iulog,*) "before ", j, ccon_ch4s_sat(c,j-1), ccon_ch4s_sat(c,j)
		ccon_ch4s_temp(c,j) = (ccon_ch4s(c,j-1) - ccon_ch4s(c,j)) * Fick_D_w(1) * m_Fick_ad * 1.0e-3 * (t_soisno(c,j)/298)**1.87 * dtsrf / (z(c,j) - z(c,j-1)) !CH4_dif
		ccon_o2s_temp(c,j) = (ccon_o2s(c,j-1) - ccon_o2s(c,j)) * Fick_D_w(2) * m_Fick_ad * 1.0e-3 * (t_soisno(c,j)/298)**1.87 * dtsrf / (z(c,j) - z(c,j-1)) !O2_dif
		ccon_co2s_temp(c,j) = (ccon_co2s(c,j-1) - ccon_co2s(c,j)) * Fick_D_w(3) * m_Fick_ad * 1.0e-3 *  (t_soisno(c,j)/298)**1.87 * dtsrf / (z(c,j) - z(c,j-1)) !CO2_dif
		ccon_h2s_temp(c,j) = (ccon_h2s(c,j-1) - ccon_h2s(c,j)) * Fick_D_w(4) * m_Fick_ad * 1.0e-3 * (t_soisno(c,j)/298)**1.87 * dtsrf / (z(c,j) - z(c,j-1)) !H2_dif

		ccon_ch4s(c,j-1) = (ccon_ch4s(c,j-1) * dz(c,j-1) - ccon_ch4s_temp(c,j)) / dz(c,j-1)
		ccon_ch4s(c,j) = (ccon_ch4s(c,j) * dz(c,j) + ccon_ch4s_temp(c,j)) / dz(c,j)
		
		ccon_o2s(c,j-1) = (ccon_o2s(c,j-1) * dz(c,j-1) - ccon_o2s_temp(c,j)) / dz(c,j-1)
		ccon_o2s(c,j) = (ccon_o2s(c,j) * dz(c,j) + ccon_o2s_temp(c,j)) / dz(c,j)

		ccon_co2s(c,j-1) = (ccon_co2s(c,j-1) * dz(c,j-1) - ccon_co2s_temp(c,j)) / dz(c,j-1)
		ccon_co2s(c,j) = (ccon_co2s(c,j) * dz(c,j) + ccon_co2s_temp(c,j)) / dz(c,j)

		ccon_h2s(c,j-1) = (ccon_h2s(c,j-1) * dz(c,j-1) - ccon_h2s_temp(c,j)) / dz(c,j-1)
		ccon_h2s(c,j) = (ccon_h2s(c,j) * dz(c,j) + ccon_h2s_temp(c,j)) / dz(c,j)
	!write(iulog,*) "after ", ccon_ch4s(c,j-1), ccon_ch4s_sat(c,j)
	end do
        end if
	end do

end subroutine vert_diffusion

#endif

end module mod_clm_cncstateupdate4
