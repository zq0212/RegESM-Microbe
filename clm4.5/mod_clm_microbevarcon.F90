module mod_clm_microbevarcon

#ifdef MICROBE
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: microbevarcon
!
! !DESCRIPTION:
! Module containing microbe parameters and logical switches and routine to read constants from CLM namelist.
!
! !USES:
  use mod_stdio
  use mod_realkinds
  use mod_mpmessage
  use mod_dynparam
  use mod_mppparam
  use mod_runparams , only : namelistfile
!
! !PUBLIC TYPES:
  implicit none
  save
!
! Microbe Model Parameters
! 
!!! zhangq modified for Microbe (define rgas)
  ! public
  private

!!! zhangq add public for Microbe (define rgas)
  real(rk8), public, PARAMETER :: micbiocn = 8  
  real(rk8), public, PARAMETER :: micbioMR = 0.123
  real(rk8), public, PARAMETER :: dock = 0.031
  real(rk8), public, PARAMETER :: CUEref = 0.43 
  real(rk8), public, PARAMETER :: CUEt = -0.012  
  real(rk8), public, PARAMETER :: Tcueref = 288.15  
  real(rk8), public, PARAMETER :: Tsref = 283.15  
  real(rk8), public, PARAMETER :: Tmref = 285.15  
  real(rk8), public, PARAMETER :: Tsmin = 272.15 
  real(rk8), public, PARAMETER :: Tmmin = 271.15 
  real(rk8), public, PARAMETER :: Msmin = 0.01  
  real(rk8), public, PARAMETER :: Msmax = 1.0  
  real(rk8), public, PARAMETER :: Mmmin = 0.1
  real(rk8), public, PARAMETER :: MFGbiomin = 1e-15

  integer, public, parameter :: nummicrobepar = 88
  real(rk8), public :: q10ch4base = 295._rk8  ! Rough estimate from comparison between Walter and previous CLM-CH4 data
  ! Uses Michigan bog data from Shannon & White
  ! This is the temperature at which the effective f_ch4 actually equals the constant f_ch4.

  real(rk8), public :: q10ch4 = 1.33_rk8 ! Production Q10
  ! Note that this is the additional Q10 for methane production ABOVE the soil decomposition temperature relationship.
  ! Corresponds to a methanogenesis Q10 of 2 when SOM HR has Q10 of 1.5.
  ! (No assumption is made in CH4 code about SOM HR temperature relationship.)
  ! Note that this formulation should be improved by making anaerobic decomposition in general a stronger function of
  ! temperature, not just methane production.

  real(rk8), public :: vmax_ch4_oxid = 45.e-6_rk8 * 1000._rk8 / 3600._rk8 ! [mol/m3-w/s];
  ! 45 uM/h from Walter and Heimann for the Mich. site (2000)
  ! oxidation rate constant (Walter and Heimann 2000)

  real(rk8), public :: k_m = 5.e-6_rk8 * 1000._rk8 ! [mol/m3-w]
  ! Michaelis-Menten oxidation rate constant for CH4 concentration (Walter and Heimann 2000)

  real(rk8), public :: q10_ch4oxid = 1.9_rk8 ! Segers, 1998
  ! Q10 oxidation constant (Walter and Heimann 2000)

  real(rk8), public :: smp_crit = -2.4e5_rk8 ! mm. From Schnell & King, 1996.
  ! Critical soil moisture potential to reduce oxidation (mm) due to dessication of methanotrophs above the water table.
  ! To turn off limitation, set to very large negative value.

  real(rk8), public :: aereoxid = -1._rk8 ! fraction of methane flux entering aerenchyma rhizosphere that will be
  ! oxidized rather than emitted.  From Wania.
  ! Note, this has been replaced by prognostic O2 diffusion into aerenchyma and is set to -1 by default.
  ! Set to value between 0 & 1 (inclusive) for sensitivity tests.  In particular, to calculate the 
  ! prognostic fraction analogous to the Wania parameter, set to 0. and compare to a normal run.

  real(rk8), public :: mino2lim = 0.2_rk8 ! minimum anaerobic decomposition rate as a fraction of potential aerobic rate
  ! for soil decomposition (or diagnostic O2-limitation / seasonal inundation factor)

  real(rk8), public :: rootlitfrac = 0.50_rk8 ! Fraction of soil organic matter associated with roots
  ! Used to partition the production between rootfr and the top 5 layers (ifndef VERTSOILC)

  real(rk8), public :: scale_factor_aere = 1.0_rk8 ! scale factor on the aerenchyma area for sensitivity tests

  real(rk8), public :: vgc_max = 0.15_rk8 ! ratio of saturation pressure triggering ebullition

  real(rk8), public :: organic_max  = 130._rk8 ! organic matter content (kg/m3) where soil is assumed to act like peat
  ! for diffusion. Very large values will lead to all soil being treated as mineral. Negative values will lead   ! to all soil being treated as peat.

  real(rk8), public :: satpow = 2._rk8 ! exponent on watsat for saturated soil solute diffusion
  ! (2 = Buckingham / Moldrup; 4/3 = Millington-Quirk)

  real(rk8), public :: cnscalefactor = 1._rk8 ! scale factor on CN decomposition for assigning methane flux
                                    ! This should equal 1 except for sensitivity studies.

  real(rk8), public :: f_ch4 = 0.2_rk8
                          ! originally 25% / (100% + 25%) from Wania. This is the ratio of CH4 production to total C
                          ! mineralization.
                          ! fraction of total decomposition that comes off as CH4 rather than CO2
                          ! Effective value will depend on temperature, redox, & pH but cannot exceed 50%
                          ! based on stoichiometry.
                          ! Note this is a crude parameterization: values in the field vary widely.

  logical, public :: transpirationloss = .true. ! switch for activating CH4 loss from transpiration
                                      ! Transpiration loss assumes that the methane concentration in dissolved soil
                                      ! water remains constant through the plant and is released when the water evaporates
                                      ! from the stomata.
                                      ! Currently hard-wired to true; impact is < 1 Tg CH4/yr

  real(rk8), public :: k_m_o2 = 20.e-6_rk8 * 1000._rk8 ! [mol/m3-w]
  ! Michaelis-Menten oxidation rate constant for O2 concentration (Segers 1998, Lidstrom and Somers 1984)

  real(rk8), public :: nongrassporosratio = 0.33_rk8 ! Ratio of root porosity in non-grass to grass, used for aerenchyma transport
                                           ! Some values in Colmer 2003

  logical, public :: allowlakeprod = .false. ! Switch to allow production under lakes based on soil carbon dataset
                                     ! (Methane can be produced, and CO2 produced from methane oxidation,
                                     ! which will slowly reduce the available carbon stock, if ! replenishlakec, but no other biogeochem is done.)
                                     ! Note: switching this off turns off ALL lake methane biogeochem. However, 0 values
                                     ! will still be averaged into the concentration _sat history fields.

  ! Base decomposition rate (1/s) at 25C
  ! Equates to about a 200 year lifetime.
  real(rk8), public :: lake_decomp_fact = 9.e-11_rk8
 
 
  ! Switch for keeping carbon storage under lakes constant
  ! so that lakes do not affect the carbon balance
  ! Good for long term rather than transient warming experiments
  ! NOTE SWITCHING THIS OFF ASSUMES TRANSIENT CARBON SUPPLY FROM LAKES;
  ! COUPLED MODEL WILL NOT CONSERVE CARBON IN THIS MODE.
  logical, public :: replenishlakec = .true.

  logical, public :: usephfact = .false. ! Switch to use pH factor in methane production

  real(rk8), public :: k_m_unsat = 5.e-6_rk8 * 1000._rk8 / 10._rk8 ! [mol/m3-w]
  ! Michaelis-Menten oxidation rate constant for CH4 concentration: literature suggests that methanotrophs
  ! in upland areas have higher affinity for methane in order to access the low ambient concentrations above the water
  ! table. (See Bender & Conrad, 1992, etc.)

  real(rk8), public :: vmax_oxid_unsat = 45.e-6_rk8 * 1000._rk8 / 3600._rk8 / 10._rk8 ! [mol/m3-w/s]
  ! Literature suggests that while k_m is lower, vmax is also lower in upland areas.

  real(rk8), public :: scale_factor_gasdiff = 1.0_rk8 ! For sensitivity tests; convection would allow this to be > 1

  real(rk8), public :: scale_factor_liqdiff = 1.0_rk8 ! For sensitivity tests; convection would allow this to be > 1

  real(rk8), public :: redoxlag = 30.0_rk8 ! Number of days to lag in the calculation of finundated_lag, which will
                                 ! be used to assess the availability of alternative electron acceptors in recently
                                 ! inundated area, reducing production.
                                 ! Set to 0 to turn off this feature.

  ! New namelists added 6/12/11

  logical, public :: fin_use_fsat = .false. ! Use fsat rather than the inversion to Prigent satellite inundation obs. (applied to
                                    ! CLM water table depth and surface runoff) to calculated finundated which is
                                    ! used in methane code and potentially soil code
                                    !!!! Attn EK: Set this to true when Sean Swenson's prognostic, tested
                                       ! fsat is integrated. (CLM4 fsat is bad for these purposes.)

  real(rk8), public :: unsat_aere_ratio = 0.05_rk8 / 0.3_rk8 ! Ratio to multiply upland vegetation aerenchyma porosity by compared to
                                        ! inundated systems. Note: porosity will be kept at above a minimum residual value
                                        ! porosmin set in subroutine ch4_aere.

  logical, public :: usefrootc = .false.    ! Use CLMCN fine root C rather than ann NPP & LAI based parameterization to
                                    ! calculate tiller C for aerenchyma area calculation.
                                    ! The NPP & LAI param. was based on Wania for Arctic sedges and may not be
                                    ! appropriate for woody PFTs, although nongrassporosratio above partly adjusts
                                    ! for this.  However, using fine root C reduces the aerenchyma area by a large
                                    ! factor.

  logical, public :: ch4offline = .false.    ! true --> Methane is not passed between the land & atmosphere.
                                    ! NEM is not added to NEE flux to atm. to correct for methane production,
                                    ! and ambient CH4 is set to constant 2009 value.

  logical, public :: ch4rmcnlim = .false.   ! Remove the N and low moisture limitations on SOM HR when calculating
                                    ! methanogenesis.
                                    ! Note: this option has not been extensively tested.
                                    ! Currently hardwired off.

  logical, public :: anoxicmicrosites = .false. ! Use Arah & Stephen 1998 expression to allow production above the water table
                                        ! Currently hardwired off; expression is crude.

  logical, public :: ch4frzout = .false.    ! Exclude CH4 from frozen fraction of soil pore H2O, to simulate "freeze-out" pulse
                                    ! as in Mastepanov 2008.
                                    ! Causes slight increase in emissions in the fall and decrease in the spring.
                                    ! Currently hardwired off; small impact.

  real(rk8), public :: redoxlag_vertical = 0._rk8   ! time lag (days) to inhibit production for newly unsaturated layers
                                          ! when decreasing WT depth for unsat. zone. See the description for redoxlag.

  real(rk8), public :: atmch4 = 1.79e-6_rk8    ! Atmospheric CH4 mixing ratio to prescribe if not provided by the atmospheric model  ! (mol/mol)
  real(rk8), public :: atmo2 = 209460e-6_rk8    ! Atmospheric CH4 mixing ratio to prescribe if not provided by the atmospheric model
  real(rk8), public :: atmco2 = 397e-6_rk8    ! Atmospheric CH4 mixing ratio to prescribe if not provided by the atmospheric model
  real(rk8), public :: atmh2 = 0.55e-6_rk8    ! Atmospheric CH4 mixing ratio to prescribe if not provided by the atmospheric model

  real(rk8), public :: plant2doc = 0.2             ! fraction of litter decomposition to available carbon
  real(rk8), public :: som2doc = 0.05          ! fraction of soil organic matter decomposition to available carbon

!xiaofeng xu creared new mechanisms and  added the following new parameters
  	real(rk8), public :: m_dKAce = 16 ! //3000	//36*1000/12/50/4    /50 for the microbe c in soc, // mMol C m-3 : convert from 36 gC m-3
	real(rk8), public :: m_dAceProdACmax = 0.0000024 !;//14.4; //14.4E3; // mMol/m3/d 0.999995; // 1-pow(0.6, 23)  0.4h-1  // Grant R.F, 1998
	real(rk8), public :: m_dKAceProdO2 = 0.004 !;  // mMol/m3
	real(rk8), public :: m_dH2ProdAcemax = 0.00000005 !;  // mMol/m3   1.35; // 0.05625 * 24mMol / mMol Acetate-C h-1 : convert from 6.75e-3 mMol acetate g-1 h-1
	real(rk8), public :: m_dKH2ProdAce = 1.65e-3 !;	 // mMol/m3   //16.5e-3;
	real(rk8), public :: m_dKCO2ProdAce = 8.25e-2 !;		// mMol/m3   8.25e-3;
	real(rk8), public :: m_dGrowRH2Methanogens = 0.01 !;		// d-1  0.023*24=0.552, (1-0.023)^24 = 0.5589
	real(rk8), public :: m_dDeadRH2Methanogens = 0.001 !;		//0.0735 d-1   0.016*24  //0.384;  0.007 * 24 = 0.168; 0.002 - 0.012 h-1
	real(rk8), public :: m_dYH2Methanogens = 0.015 !;
	real(rk8), public :: m_dGrowRAceMethanogens = 0.008 !;		// d-1
	real(rk8), public :: m_dDeadRAceMethanogens = 0.002 !;		 //0.01 d-1
	real(rk8), public :: m_dYAceMethanogens = 0.2 !;   // 0.04; molC/(molAcetate -C) : convert from 0.08 mol C (mol acetate)-1
	real(rk8), public :: m_dGrowRMethanotrophs = 0.04 !;		// d-1
	real(rk8), public :: m_dDeadRMethanotrophs = 0.005 !;		//d-1
	real(rk8), public :: m_dYMethanotrophs = 0.40 !;	// molC/(mol CH4 -C) : convert from 0.4 mol C (mol CH4)-1
	real(rk8), public :: m_dGrowRAOMMethanotrophs = 0.004 !;		// d-1
	real(rk8), public :: m_dDeadRAOMMethanotrophs = 0.002 !;		//d-1
	real(rk8), public :: m_dYAOMMethanotrophs = 0.15 !;	// molC/(mol CH4 -C) : convert from 0.4 mol C (mol CH4)-1
	real(rk8), public :: m_dAceProdQ10 = 2.0 !;
	real(rk8), public :: m_dACProdQ10 = 2.0 !;
	real(rk8), public :: m_dACMinQ10 = 3.0 !;
	real(rk8), public :: m_dAceH2min = 0.48e-8 !;		 //0.48e-3;   mMol/m3
	real(rk8), public :: m_dCH4H2min = 0.41e-9 !;	 // 0.27e-4;   mMol/m3
	real(rk8), public :: m_dKH2ProdCH4 = 7.75e-5 !;	// mMol/m3
	real(rk8), public :: m_dKCO2ProdCH4 = 3.1e-4 !;	 // mMol/m3
	real(rk8), public :: m_dH2CH4ProdQ10 = 2.0 !;
	real(rk8), public :: m_dH2AceProdQ10 = 2.0 !;
	real(rk8), public :: m_dKCH4ProdAce = 0.05 !; //0.005 !;		 // mMol/m3
	real(rk8), public :: m_dKCH4ProdO2 = 0.015 !; // 0.01;	// mMol/m3
	real(rk8), public :: m_dCH4ProdQ10 = 2.0 !;
	real(rk8), public :: m_drCH4Prod = 0.5 !;		// MolCH4(Mol acetate -C)-1 : convert from 2 MolCH4(Mol acetate)-1
	real(rk8), public :: m_dKCH4OxidCH4 = 0.1 !; //0.005;		// mMol/m3
	real(rk8), public :: m_dKAOMCH4OxidCH4 = 1.5
	real(rk8), public :: m_dKCH4OxidO2 = 1.0 !; // 0.01;		// mMol/m3
	real(rk8), public :: m_dCH4OxidQ10 = 3.0 !;
	real(rk8), public :: m_dAOMCH4OxidQ10 = 1.2 !
	real(rk8), public :: m_drAer = 0.002 !;					// MolO2(mol C)-1
	real(rk8), public :: m_dKAerO2 = 0.0010 !;				 //0.01 mMol/m3
	real(rk8), public :: m_dAerDecomQ10 = 2.0 !;
	real(rk8), public :: m_drCH4Oxid = 2 !;			// MolO2(molCH4)-1
	real(rk8), public :: m_dKe = 0.03 !;			 // 0.05 h-1
        !!! zhangq modified for Microbe !!! 0.05 >> 0.0005
	real(rk8), public :: m_dCH4min = 0.0005 !;			 //0.5 mMol/m3 //程序中对此参数的判断条件单位为 mol C/m3
	real(rk8), public :: m_dAirCH4 = 0.0893 !;  // mM/m3
	real(rk8), public :: m_dAirH2 = 0.0257 !; // mM/m3
	real(rk8), public :: m_dAirO2 = 9.375e3 !; // mM/m3
	real(rk8), public :: m_dAirCO2 = 16.295 !; // mM/m3
	
	! not in the parameter file
	real(rk8), public :: m_dPlantTrans = 0.007
	real(rk8), public :: g_dMaxH2inWater = 4.73e-4
	
	!Penning, H., P. Claus, P. Casper, and R. Conrad. 2006. Carbon isotope fractionation during acetoclastic methanogenesis by Methanosaeta concilii in culture and a lake sediment. 
	! Applied and Environmental Microbiology 72:5648-5652.
	real(rk8), public :: frac_doc = 1.0			! fractionation factor for DOC production
	real(rk8), public :: frac_ace = 1.0			! fractionation factor for acetate production
	real(rk8), public :: frac_acch4 = 1.0		! fractionation factor for aceclatic methanogenesis
	real(rk8), public :: frac_hych4 = 1.0		! fractionation factor for hydrogenotrophic methanogenesis
 	real(rk8), public :: frac_acetogenesis = 1.0	! fractionation factor for acetogenesis
	real(rk8), public :: frac_ch4ox = 1.0		! fractionation factor for methanotrophy
 	real(rk8), public :: frac_ch4aom = 1.0		! fractionation factor for anaerobic oxidation of methane
  
	real(rk8), public :: pHmin = 4.0
	real(rk8), public :: pHmax = 10.
	real(rk8), public :: pHopt = 7.
  
	real(rk8), public :: dom_diffus = 0.00000018
	real(rk8), public :: m_Fick_ad = 2.0
	  
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: ch4conrd_microbe ! Read and initialize CH4 constants
!
! !REVISION HISTORY:
! Created by Zack Subin
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4conrd
!
! !INTERFACE:
  subroutine ch4conrd_microbe ()
!
! !DESCRIPTION:
! Read and initialize CH4 constants
!
! !USES:
    !!! zhangq modified
    !use fileutils , only : opnfil, getfil, relavu, getavu
    !use spmdMod   , only : masterproc, mpicom, MPI_REAL8, MPI_LOGICAL
    !use controlMod  , only: NLFilename

!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! routine initch4
!
! !REVISION HISTORY:
! Created by Zack Subin
! Modified 6/12/2011 to use namelist rather than ASCII file.
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,n                ! loop indices
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=32) :: subname = 'ch4conrd'
    !character(len=40) :: ch4parname(nummicrobepar)  ! subroutine name
    !character(LEN=256)::locfn='./microbepar_in'
!-----------------------------------------------------------------------

    ! Driver
    namelist /ch4par_in/ &
        ch4offline, fin_use_fsat, replenishlakec, allowlakeprod, atmch4

    ! Production
    namelist /ch4par_in/ &
        q10ch4base, q10ch4, rootlitfrac, f_ch4, cnscalefactor, &
        redoxlag, &
        mino2lim, lake_decomp_fact, ch4rmcnlim,anoxicmicrosites,redoxlag_vertical

    ! Oxidation
    namelist /ch4par_in/ &
       vmax_ch4_oxid, k_m, q10_ch4oxid, smp_crit, k_m_o2, k_m_unsat, &
       vmax_oxid_unsat

    ! Aerenchyma
    namelist /ch4par_in/ &
         aereoxid, scale_factor_aere,  &
         nongrassporosratio, unsat_aere_ratio, transpirationloss, &
         usefrootc

    ! Ebullition
    namelist /ch4par_in/ &
         vgc_max

    ! Transport
    namelist /ch4par_in/ &
         organic_max, satpow, scale_factor_gasdiff, &
         scale_factor_liqdiff, ch4frzout    
    
    ! Microbe
    namelist /ch4par_in/ &
         plant2doc, som2doc,m_dKAce, &
         m_dAceProdACmax,m_dKAceProdO2,m_dH2ProdAcemax, m_dKH2ProdAce, &
         m_dKCO2ProdAce, m_dGrowRH2Methanogens,m_dDeadRH2Methanogens, &
         m_dYH2Methanogens, m_dGrowRAceMethanogens, m_dDeadRAceMethanogens, &
         m_dYAceMethanogens,m_dGrowRMethanotrophs,m_dDeadRMethanotrophs, &
         m_dYMethanotrophs, m_dGrowRAOMMethanotrophs, &
         m_dDeadRAOMMethanotrophs, m_dYAOMMethanotrophs, m_dAceProdQ10, &
         m_dACProdQ10,m_dACMinQ10,m_dAceH2min,m_dCH4H2min,m_dKH2ProdCH4, &
         m_dKCO2ProdCH4, m_dH2CH4ProdQ10, m_dH2AceProdQ10, &
         m_dKCH4ProdAce, m_dKCH4ProdO2,m_dCH4ProdQ10,m_drCH4Prod, &
         m_dKCH4OxidCH4, m_dKAOMCH4OxidCH4, m_dKCH4OxidO2,m_dCH4OxidQ10, &
         m_dAOMCH4OxidQ10, m_drAer, m_dKAerO2, m_dAerDecomQ10, &
         m_drCH4Oxid, m_dKe,m_dCH4min, m_dAirCH4, m_dAirH2, &
         m_dAirO2, m_dAirCO2,frac_doc,frac_ace, &
         frac_acch4,frac_hych4, frac_acetogenesis, frac_ch4ox, &
         frac_ch4aom,dom_diffus, m_Fick_ad, m_dPlantTrans

    call bcast(q10ch4base)         
    call bcast(q10ch4)             
    call bcast(vmax_ch4_oxid)     
    call bcast(k_m)               
    call bcast(q10_ch4oxid)       
    call bcast(smp_crit)          
    call bcast(aereoxid)          
    call bcast(mino2lim)          
    call bcast(rootlitfrac)       
    call bcast(scale_factor_aere) 
    call bcast(vgc_max)           
    call bcast(organic_max)       
    call bcast(satpow)            
    call bcast(cnscalefactor)            
    call bcast(f_ch4)            
    call bcast(transpirationloss)            
    call bcast(k_m_o2)            
    call bcast(nongrassporosratio)            
    call bcast(allowlakeprod)            
    call bcast(lake_decomp_fact)            
    call bcast(usephfact)            
    call bcast(k_m_unsat)            
    call bcast(vmax_oxid_unsat)            
    call bcast(replenishlakec)            
    call bcast(scale_factor_gasdiff)            
    call bcast(scale_factor_liqdiff)            
    call bcast(redoxlag)            
    call bcast(fin_use_fsat)            
    call bcast(unsat_aere_ratio)            
    call bcast(usefrootc)            
    call bcast(ch4offline)            
    call bcast(ch4rmcnlim)            
    call bcast(anoxicmicrosites)            
    call bcast(ch4frzout)            
    call bcast(redoxlag_vertical)            
    call bcast(atmch4)
    call bcast(plant2doc)
    call bcast(som2doc)
   
    call bcast(m_dKAce)
    call bcast(m_dAceProdACmax)
    call bcast(m_dKAceProdO2)
    call bcast(m_dH2ProdAcemax)
    call bcast(m_dKH2ProdAce)
    call bcast(m_dKCO2ProdAce)
    call bcast(m_dGrowRH2Methanogens)
    call bcast(m_dDeadRH2Methanogens)
    call bcast(m_dYH2Methanogens)
    call bcast(m_dGrowRAceMethanogens)
    call bcast(m_dDeadRAceMethanogens)    
    call bcast(m_dYAceMethanogens)
    call bcast(m_dGrowRMethanotrophs)
    call bcast(m_dDeadRMethanotrophs)
    call bcast(m_dYMethanotrophs)
    call bcast(m_dGrowRAOMMethanotrophs)    
    call bcast(m_dDeadRAOMMethanotrophs)
    call bcast(m_dYAOMMethanotrophs)
    call bcast(m_dAceProdQ10)
    call bcast(m_dACProdQ10)
    call bcast(m_dACMinQ10)    
    call bcast(m_dAceH2min)
    call bcast(m_dCH4H2min)
    call bcast(m_dKH2ProdCH4)
    call bcast(m_dKCO2ProdCH4)
    call bcast(m_dH2CH4ProdQ10)    
    call bcast(m_dH2AceProdQ10)
    call bcast(m_dKCH4ProdAce)
    call bcast(m_dKCH4ProdO2)
    call bcast(m_dCH4ProdQ10)
    call bcast(m_drCH4Prod)    
    call bcast(m_dKCH4OxidCH4)
    call bcast(m_dKAOMCH4OxidCH4)
    call bcast(m_dKCH4OxidO2)
    call bcast(m_dCH4OxidQ10)
    call bcast(m_dAOMCH4OxidQ10)    
    call bcast(m_drAer)
    call bcast(m_dKAerO2)
    call bcast(m_dAerDecomQ10)
    call bcast(m_drCH4Oxid)
    call bcast(m_dKe)    
    call bcast(m_dCH4min)
    call bcast(m_dAirCH4)
    call bcast(m_dAirH2)
    call bcast(m_dAirO2)
    call bcast(m_dAirCO2)
 
    call bcast(frac_doc)
    call bcast(frac_ace)
    call bcast(frac_acch4)
    call bcast(frac_hych4)
    call bcast(frac_acetogenesis)
    call bcast(frac_ch4ox)
    call bcast(frac_ch4aom)
        
    call bcast(dom_diffus)
    call bcast(m_Fick_ad)
    call bcast(m_dPlantTrans)   
    
    if (myid == iocpu) then
      
      write(stdout,*) 'Successfully read CH4 namelist'
      write(stdout,*)' '
      write(stdout,*)'q10ch4base = ', q10ch4base
      write(stdout,*)'q10ch4 = ', q10ch4
      write(stdout,*)'vmax_ch4_oxid = ', vmax_ch4_oxid
      write(stdout,*)'k_m = ', k_m
      write(stdout,*)'q10_ch4oxid = ', q10_ch4oxid
      write(stdout,*)'smp_crit = ', smp_crit
      write(stdout,*)'aereoxid = ', aereoxid
      write(stdout,*)'mino2lim = ', mino2lim
      write(stdout,*)'rootlitfrac = ', rootlitfrac
      write(stdout,*)'scale_factor_aere = ', scale_factor_aere
      write(stdout,*)'vgc_max = ', vgc_max
      write(stdout,*)'organic_max = ', organic_max
      write(stdout,*)'satpow = ', satpow
      write(stdout,*)'cnscalefactor = ', cnscalefactor
      write(stdout,*)'f_ch4 = ', f_ch4
      write(stdout,*)'transpirationloss = ', transpirationloss
      write(stdout,*)'k_m_o2 = ', k_m_o2
      write(stdout,*)'nongrassporosratio = ', nongrassporosratio
      write(stdout,*)'allowlakeprod = ', allowlakeprod
      write(stdout,*)'lake_decomp_fact = ', lake_decomp_fact
      write(stdout,*)'usephfact = ', usephfact
      write(stdout,*)'k_m_unsat = ', k_m_unsat
      write(stdout,*)'vmax_oxid_unsat = ', vmax_oxid_unsat
      write(stdout,*)'replenishlakec = ', replenishlakec
      write(stdout,*)'scale_factor_gasdiff = ', scale_factor_gasdiff
      write(stdout,*)'scale_factor_liqdiff = ', scale_factor_liqdiff
      write(stdout,*)'redoxlag = ', redoxlag
      write(stdout,*)'fin_use_fsat = ', fin_use_fsat
      write(stdout,*)'unsat_aere_ratio = ', unsat_aere_ratio
      write(stdout,*)'usefrootc = ', usefrootc
      write(stdout,*)'ch4offline = ', ch4offline
      write(stdout,*)'ch4rmcnlim = ', ch4rmcnlim
      write(stdout,*)'anoxicmicrosites = ', anoxicmicrosites
      write(stdout,*)'ch4frzout = ', ch4frzout
      write(stdout,*)'redoxlag_vertical = ', redoxlag_vertical
      write(stdout,*)'atmch4 = ', atmch4
      write(stdout,*)'plant2doc = ', plant2doc
      write(stdout,*)'som2doc = ', som2doc

      write(stdout,*)'m_dKAce = ', m_dKAce
      write(stdout,*)'m_dAceProdACmax = ', m_dAceProdACmax
      write(stdout,*)'m_dKAceProdO2 = ', m_dKAceProdO2
      write(stdout,*)'m_dH2ProdAcemax = ', m_dH2ProdAcemax
      write(stdout,*)'m_dKH2ProdAce = ', m_dKH2ProdAce
      write(stdout,*)'m_dKCO2ProdAce = ', m_dKCO2ProdAce
      write(stdout,*)'m_dGrowRH2Methanogens = ', m_dGrowRH2Methanogens
      write(stdout,*)'m_dDeadRH2Methanogens = ', m_dDeadRH2Methanogens
      write(stdout,*)'m_dYH2Methanogens = ', m_dYH2Methanogens
      write(stdout,*)'m_dGrowRAceMethanogens = ', m_dGrowRAceMethanogens
      write(stdout,*)'m_dDeadRAceMethanogens = ', m_dDeadRAceMethanogens
      write(stdout,*)'m_dYAceMethanogens = ', m_dYAceMethanogens
      write(stdout,*)'m_dGrowRMethanotrophs = ', m_dGrowRMethanotrophs
      write(stdout,*)'m_dDeadRMethanotrophs = ', m_dDeadRMethanotrophs
      write(stdout,*)'m_dYMethanotrophs = ', m_dYMethanotrophs
      write(stdout,*)'m_dGrowRAOMMethanotrophs = ', m_dGrowRAOMMethanotrophs
      write(stdout,*)'m_dDeadRAOMMethanotrophs = ', m_dDeadRAOMMethanotrophs
      write(stdout,*)'m_dYAOMMethanotrophs = ', m_dYAOMMethanotrophs
      write(stdout,*)'m_dAceProdQ10 = ', m_dAceProdQ10
      write(stdout,*)'m_dACProdQ10 = ', m_dACProdQ10
      write(stdout,*)'m_dACMinQ10 = ', m_dACMinQ10
      write(stdout,*)'m_dAceH2min = ', m_dAceH2min
      write(stdout,*)'m_dCH4H2min = ', m_dCH4H2min
      write(stdout,*)'m_dKH2ProdCH4 = ', m_dKH2ProdCH4
      write(stdout,*)'m_dKCO2ProdCH4 = ', m_dKCO2ProdCH4
      write(stdout,*)'m_dH2CH4ProdQ10 = ', m_dH2CH4ProdQ10
      write(stdout,*)'m_dH2AceProdQ10 = ', m_dH2AceProdQ10
      write(stdout,*)'m_dKCH4ProdAce = ', m_dKCH4ProdAce
      write(stdout,*)'m_dKCH4ProdO2 = ', m_dKCH4ProdO2
      write(stdout,*)'m_dCH4ProdQ10 = ', m_dCH4ProdQ10
      write(stdout,*)'m_drCH4Prod = ', m_drCH4Prod
      write(stdout,*)'m_dKCH4OxidCH4 = ', m_dKCH4OxidCH4
      write(stdout,*)'m_dKAOMCH4OxidCH4 = ', m_dKAOMCH4OxidCH4
      write(stdout,*)'m_dKCH4OxidO2 = ', m_dKCH4OxidO2
      write(stdout,*)'m_dCH4OxidQ10 = ', m_dCH4OxidQ10
      write(stdout,*)'m_dAOMCH4OxidQ10 = ', m_dAOMCH4OxidQ10
      write(stdout,*)'m_drAer = ', m_drAer
      write(stdout,*)'m_dKAerO2 = ', m_dKAerO2
      write(stdout,*)'m_dAerDecomQ10 = ', m_dAerDecomQ10
      write(stdout,*)'m_drCH4Oxid = ', m_drCH4Oxid
      write(stdout,*)'m_dKe = ', m_dKe
      write(stdout,*)'m_dCH4min = ', m_dCH4min
      write(stdout,*)'m_dAirCH4 = ', m_dAirCH4
      write(stdout,*)'m_dAirH2 = ', m_dAirH2
      write(stdout,*)'m_dAirO2 = ', m_dAirO2        
      write(stdout,*)'m_dAirCO2 = ', m_dAirCO2
        
      write(stdout,*)'frac_doc = ', frac_doc
      write(stdout,*)'frac_ace = ', frac_ace
      write(stdout,*)'frac_acch4 = ', frac_acch4
      write(stdout,*)'frac_hych4 = ', frac_hych4
      write(stdout,*)'frac_acetogenesis = ', frac_acetogenesis
      write(stdout,*)'frac_ch4ox = ', frac_ch4ox
      write(stdout,*)'frac_ch4aom = ', frac_ch4aom
        
      write(stdout,*)'dom_diffus = ', dom_diffus 
      write(stdout,*)'m_dPlantTrans = ', m_dPlantTrans

      if (ch4offline) then
        write(stdout,*)'CH4 Model will be running offline and not affect fluxes to atmosphere.'
      end if
      if (aereoxid >= 0._rk8) then
        write(stdout,*) 'Fixed aerenchyma oxidation has been selected.'
      end if
      if (.not. allowlakeprod) then
        write(stdout,*) 'Lake production has been disabled. Lakes will not factor into CH4 BGC.'//&
                       '"Sat" history fields will not average over lakes except for concentrations,'//&
                       'which will average zero from lakes.'
      end if
        write(stdout,*)'Successfully initialized CH4 parameters'
        write(stdout,*)

    end if

  end subroutine ch4conrd_microbe

#endif
! defined MICROBE
end module mod_clm_microbevarcon

