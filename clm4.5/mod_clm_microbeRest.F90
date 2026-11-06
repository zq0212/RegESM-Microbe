module mod_clm_microbeRest
#ifdef MICROBE

!
! Reads from or writes restart data
!
  use mod_intkinds
  use mod_realkinds
  use mod_mpmessage
  use mod_stdio
  use mod_runparams
  use mod_clm_type
  use mod_clm_decomp
  use mod_clm_nchelper
  use mod_clm_varctl , only : nsrest, use_c13, use_c14
  use mod_clm_time_manager , only : is_restart 
  
!
! !PUBLIC TYPES:
  implicit none
  private
! save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: microbeRest
!
! !REVISION HISTORY:
! Created by Xiaofeng Xu in 2013
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: micrbeRest
!
! !INTERFACE:
  subroutine microbeRest( ncid, flag )
!
! !DESCRIPTION:
! Read/Write biogeophysics information to/from restart file.
!
! !ARGUMENTS:
    implicit none
    type(clm_filetype), intent(inout) :: ncid ! netcdf id
    character(len=*), intent(in) :: flag     ! 'read' or 'write'
!
! !CALLED FROM:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: c,l,g,j      ! indices
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    logical :: readvar      ! determine if variable is on initial file
    character(len=128) :: varname         ! temporary
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
!-----------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p
    !finundated => cptr%cws%finundated

    ! column microbial state variable - cdocs

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'pH', &
              cdims=(/'column '/), &
              long_name='pH', units='unitless')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'pH') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'pH', &
                cptr%cps%pH,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'pH', &
              cptr%cps%pH,gcomm_column)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'soilpH_unsat', &
              cdims=(/'column ','levgrnd'/), &
              long_name='soilpH in unsat', units='unitless', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'soilpH_unsat') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'soilpH_unsat', &
                cptr%cps%soilpH_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'soilpH_unsat', &
              cptr%cps%soilpH_unsat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'soilpH_sat', &
              cdims=(/'column ','levgrnd'/), &
              long_name='soilpH in sat', units='unitless', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'soilpH_sat') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'soilpH_sat', &
                cptr%cps%soilpH_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'soilpH_sat', &
              cptr%cps%soilpH_unsat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'FINUNDATED', &
              cdims=(/'column '/), &
              long_name='fraction of inundated area', units='m^2/m^2')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'FINUNDATED') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'FINUNDATED', &
                cptr%cws%finundated,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'FINUNDATED', &
              cptr%cws%finundated,gcomm_column)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'fsat_pre', &
              cdims=(/'column '/), &
              long_name='previous time step fraction of inundated area', units='m^2/m^2')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'fsat_pre') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'fsat_pre', &
                cptr%cmic%fsat_pre,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'fsat_pre', &
              cptr%cmic%fsat_pre,gcomm_column)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'fsat', &
              cdims=(/'column '/), &
              long_name='fraction of inundated area',units='m^2/m^2')
    else if ( flag == 'read' ) then      
      if ( .not. clm_check_var(ncid,'fsat') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping') 
        end if
      else
        call clm_readvar(ncid,'fsat', &
                cptr%cws%fsat,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'fsat', &
              cptr%cws%fsat,gcomm_column)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'waterhead_unsat', &
              cdims=(/'column '/), &
              long_name='fraction of inundated area',units='m^2/m^2')
    else if ( flag == 'read' ) then      
      if ( .not. clm_check_var(ncid,'waterhead_unsat') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'waterhead_unsat', &
                cptr%cmic%waterhead_unsat,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'waterhead_unsat', &
              cptr%cmic%waterhead_unsat,gcomm_column)
    end if

!20150831
   !~ if (flag == 'define') then
       !~ call ncd_defvar(ncid=ncid, varname='rr', xtype=ncd_double, &
            !~ dim1name='column', long_name='respiration', units='mol/s')
    !~ else if (flag == 'read' .or. flag == 'write') then
       !~ call ncd_io(varname='rr', data=pcf%rr, &
            !~ dim1name='column', ncid=ncid, flag=flag, readvar=readvar)
       !~ if (flag=='read' .and. .not. readvar) then
          !~ if (is_restart()) call endrun()
       !~ end if
    !~ end if
    
   !~ if (flag == 'define') then
       !~ call ncd_defvar(ncid=ncid, varname='froot_mr', xtype=ncd_double, &
            !~ dim1name='column', long_name='fine root maintenance respiration', units='mol/s')
    !~ else if (flag == 'read' .or. flag == 'write') then
       !~ call ncd_io(varname='froot_mr', data=pcf%froot_mr, &
            !~ dim1name='column', ncid=ncid, flag=flag, readvar=readvar)
       !~ if (flag=='read' .and. .not. readvar) then
          !~ if (is_restart()) call endrun()
       !~ end if
    !~ end if

    !~ if (flag == 'define') then
       !~ call ncd_defvar(ncid=ncid, varname='rootfr', xtype=ncd_double, &
            !~ dim1name='column', dim2name='levgrnd', switchdim=.true., &
            !~ long_name='root rerspiration fraction', units='%')
    !~ else if (flag == 'read' .or. flag == 'write') then
       !~ call ncd_io(varname='rootfr', data=pps%rootfr, &
            !~ dim1name='column', switchdim=.true., &
            !~ ncid=ncid, flag=flag, readvar=readvar)
       !~ if (flag=='read' .and. .not. readvar) then
          !~ if (is_restart()) call endrun()
       !~ end if
    !~ end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'hr_vr', &
           cdims=(/'column ','levgrnd'/), &
           long_name='hetrotrophic respiration fraction', units='mol/m^2', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'hr_vr') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'hr_vr', &
                cptr%ccf%hr_vr,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'hr_vr', &
              cptr%ccf%hr_vr,gcomm_column, switchdim=.true.)
    end if
    
 !20150831   
    !~ if (flag == 'define') then
       !~ call ncd_defvar(ncid=ncid, varname='sminn_vr', xtype=ncd_double, &
            !~ dim1name='column', dim2name='levgrnd', switchdim=.true., &
            !~ long_name='soil mineral nitrogen', units='mol/m^2')
    !~ else if (flag == 'read' .or. flag == 'write') then
       !~ call ncd_io(varname='sminn_vr', data=cns%sminn_vr, &
            !~ dim1name='column', switchdim=.true., &
            !~ ncid=ncid, flag=flag, readvar=readvar)
       !~ if (flag=='read' .and. .not. readvar) then
          !~ if (is_restart()) call endrun()
       !~ end if
    !~ end if

!20150830

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'froot_r', &
              cdims=(/'column ','levgrnd'/), &
              long_name='fine root respiration', units='mol/s', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'froot_r') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'froot_r', &
                cptr%cmic%froot_r,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'froot_r', &
              cptr%cmic%froot_r,gcomm_column, switchdim=.true.)
    end if

!~ #if (defined HUM_HOL)	
   !~ if (flag == 'define') then
       !~ call ncd_defvar(ncid=ncid, varname='qflx_lat_aqu_layer', xtype=ncd_double, &
            !~ dim1name='column', dim2name='levgrnd', switchdim=.true., &
            !~ long_name='lateral water flux', units='mol/s')
    !~ else if (flag == 'read' .or. flag == 'write') then
       !~ call ncd_io(varname='qflx_lat_aqu_layer', data=cwf%qflx_lat_aqu_layer, &
            !~ dim1name='column', switchdim=.true., &
            !~ ncid=ncid, flag=flag, readvar=readvar)
       !~ if (flag=='read' .and. .not. readvar) then
          !~ if (is_restart()) call endrun()
       !~ end if
    !~ end if
!~ #endif
!20150830

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CACES_PROD', &
           cdims=(/'column ','levgrnd'/), &              
           long_name='acetate production', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CACES_PROD') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CACES_PROD', &
                cptr%cmic%caces_prod,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CACES_PROD', &
              cptr%cmic%caces_prod,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CACES_UNSAT_PROD', &
           cdims=(/'column ','levgrnd'/), &
           long_name='acetate production in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CACES_UNSAT_PROD') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CACES_UNSAT_PROD', &
             cptr%cmic%caces_unsat_prod,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CACES_UNSAT_PROD', &
              cptr%cmic%caces_unsat_prod,gcomm_column, switchdim=.true.)
    end if
 
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CACES_SAT_PROD', &
           cdims=(/'column ','levgrnd'/), &
           long_name='acetate production in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CACES_SAT_PROD') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'CACES_SAT_PROD', &
                cptr%cmic%caces_sat_prod,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CACES_SAT_PROD', &
              cptr%cmic%caces_sat_prod,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CACES_PROD_H2', &
           cdims=(/'column ','levgrnd'/), &
           long_name='acetogenesis', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CACES_PROD_H2') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CACES_PROD_H2', &
                cptr%cmic%caces_prod_h2,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CACES_PROD_H2', &
              cptr%cmic%caces_prod_h2,gcomm_column, switchdim=.true.)
    end if
  
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CACES_UNSAT_PROD_H2', &
           cdims=(/'column ','levgrnd'/), & 
           long_name='acetogenesis in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CACES_UNSAT_PROD_H2') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CACES_UNSAT_PROD_H2', &
                cptr%cmic%caces_unsat_prod_h2,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CACES_UNSAT_PROD_H2', &
              cptr%cmic%caces_unsat_prod_h2,gcomm_column, switchdim=.true.)
    end if
 
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CACES_SAT_PROD_H2', &
           cdims=(/'column ','levgrnd'/), & 
           long_name='acetogenesis in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CACES_SAT_PROD_H2') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'CACES_SAT_PROD_H2', &
                cptr%cmic%caces_sat_prod_h2,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CACES_SAT_PROD_H2', &
              cptr%cmic%caces_sat_prod_h2,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CDOCS_PRE', &
           cdims=(/'column ','levgrnd'/), & 
           long_name='previous time step DOC concentration', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CDOCS_PRE') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'CDOCS_PRE', &
                cptr%cmic%cdocs_pre,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CDOCS_PRE', &
              cptr%cmic%cdocs_pre,gcomm_column, switchdim=.true.)
    end if
  
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CDONS_MIN', &
           cdims=(/'column ','levgrnd'/), & 
           long_name='dissolved organic nitrogen minterlization ', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CDONS_MIN') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'CDONS_MIN', &
                cptr%cmic%cdons_min,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CDONS_MIN', &
              cptr%cmic%cdons_min,gcomm_column, switchdim=.true.)
    end if
   
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CDOCS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='dissolved organic carbon',units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CDOCS') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CDOCS', &
                cptr%cmic%cdocs,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CDOCS', &
              cptr%cmic%cdocs,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CDOCS_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='dissolved organic carbon in unsaturated fraction',units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CDOCS_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CDOCS_UNSAT', &
                cptr%cmic%cdocs_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CDOCS_UNSAT', &
              cptr%cmic%cdocs_unsat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CDOCS_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='dissolved organic carbon in saturated fraction',units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CDOCS_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CDOCS_SAT', &
                cptr%cmic%cdocs_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CDOCS_SAT', &
              cptr%cmic%cdocs_sat,gcomm_column, switchdim=.true.)
    end if
   
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CDONS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='dissolved organic nitrogen',units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CDONS') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CDONS', &
                cptr%cmic%cdons,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CDONS', &
              cptr%cmic%cdons,gcomm_column, switchdim=.true.)
    end if
   
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CDONS_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='dissolved organic nitrogen in unsaturated fraction',units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CDONS_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CDONS_UNSAT', &
                cptr%cmic%cdons_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CDONS_UNSAT', &
              cptr%cmic%cdons_unsat,gcomm_column, switchdim=.true.)
    end if
 
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CDONS_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='dissolved organic nitrogen in saturated fraction',units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CDONS_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CDONS_SAT', &
                cptr%cmic%cdons_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CDONS_SAT', &
              cptr%cmic%cdons_sat,gcomm_column, switchdim=.true.)
    end if
    
    ! column microbial state variable - cmicbiocs

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CMICBIOCS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='microbial biomass carbon',units='molC/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CMICBIOCS') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'CMICBIOCS', &
                cptr%cmic%cmicbiocs,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CMICBIOCS', &
              cptr%cmic%cmicbiocs,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CMICBIONS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='microbial biomass nitrogen',units='molN/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CMICBIONS') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'CMICBIONS', &
                cptr%cmic%cmicbions,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CMICBIONS', &
              cptr%cmic%cmicbions,gcomm_column, switchdim=.true.)
    end if
    
    ! column microbial state variable - caces
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CACES', &
           cdims=(/'column ','levgrnd'/), &
           long_name='acetic acid',units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CACES') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'CACES', &
                cptr%cmic%caces,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CACES', &
              cptr%cmic%caces,gcomm_column, switchdim=.true.)
     end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CACES_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='acetic acid in unsaturated fraction',units='mol/m^3',switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CACES_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CACES_UNSAT', &
                cptr%cmic%caces_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CACES_UNSAT', &
              cptr%cmic%caces_unsat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CACES_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='acetic acid in saturated fraction',units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CACES_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CACES_SAT', &
                cptr%cmic%caces_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CACES_SAT', &
              cptr%cmic%caces_sat,gcomm_column, switchdim=.true.)
    end if
    
    ! column microbial state variable - cacebios
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CACEBIOS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='methanogenesis based on acetic acid',units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CACEBIOS') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'CACEBIOS', &
                cptr%cmic%cacebios,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CACEBIOS', &
              cptr%cmic%cacebios,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CACEBIOS_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='methanogenesis based on acetic acid in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CACEBIOS_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping') 
        end if
      else  
        call clm_readvar(ncid,'CACEBIOS_UNSAT', &
                cptr%cmic%cacebios_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CACEBIOS_UNSAT', &
              cptr%cmic%cacebios_unsat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CACEBIOS_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='methanogenesis based on acetic acid in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CACEBIOS_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'CACEBIOS_SAT', &
                cptr%cmic%cacebios_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CACEBIOS_SAT', &
              cptr%cmic%cacebios_sat,gcomm_column, switchdim=.true.)
    end if
   
    ! column microbial state variable - cco2bios
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CCO2BIOS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='methanogensis based on co2', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CCO2BIOS') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'CCO2BIOS', &
                cptr%cmic%cco2bios,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CCO2BIOS', &
              cptr%cmic%cco2bios,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CCO2BIOS_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='methanogensis based on co2 in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CCO2BIOS_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping') 
        end if
      else  
        call clm_readvar(ncid,'CCO2BIOS_UNSAT', &
                cptr%cmic%cco2bios_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CCO2BIOS_UNSAT', &
              cptr%cmic%cco2bios_unsat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CCO2BIOS_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='methanogensis based on co2 in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CCO2BIOS_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CCO2BIOS_SAT', &
                cptr%cmic%cco2bios_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CCO2BIOS_SAT', &
              cptr%cmic%cco2bios_sat,gcomm_column, switchdim=.true.)
    end if

    !~ ! column microbial state variable - caerch4bios
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CAERCH4BIOS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='aerobic oxidationo of ch4', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CAERCH4BIOS') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CAERCH4BIOS', &
                cptr%cmic%caerch4bios,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CAERCH4BIOS', &
              cptr%cmic%caerch4bios,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CAERCH4BIOS_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='aerobic oxidationo of ch4 in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CAERCH4BIOS_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'CAERCH4BIOS_UNSAT', &
                cptr%cmic%caerch4bios_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CAERCH4BIOS_UNSAT', &
              cptr%cmic%caerch4bios_unsat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CAERCH4BIOS_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='aerobic oxidationo of ch4 in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CAERCH4BIOS_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'CAERCH4BIOS_SAT', &
                cptr%cmic%caerch4bios_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CAERCH4BIOS_SAT', &
              cptr%cmic%caerch4bios_sat,gcomm_column, switchdim=.true.)
    end if

    !~ ! column microbial state variable - canaerch4bios   
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CANAERCH4BIOS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='biomass for microbes for anaerobic methane oxidation', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CANAERCH4BIOS') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'CANAERCH4BIOS', &
                cptr%cmic%canaerch4bios,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CANAERCH4BIOS', &
              cptr%cmic%canaerch4bios,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CANAERCH4BIOS_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='biomass for microbes for anaerobic methane oxidation in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CANAERCH4BIOS_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CANAERCH4BIOS_UNSAT', &
                cptr%cmic%canaerch4bios_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CANAERCH4BIOS_UNSAT', &
              cptr%cmic%canaerch4bios_unsat,gcomm_column, switchdim=.true.)
    end if
 
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CANAERCH4BIOS_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='biomass for microbes for anaerobic methane oxidation in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CANAERCH4BIOS_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CANAERCH4BIOS_SAT', &
                cptr%cmic%canaerch4bios_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CANAERCH4BIOS_SAT', &
              cptr%cmic%canaerch4bios_sat,gcomm_column, switchdim=.true.)
    end if
   
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CCON_CH4S', &
           cdims=(/'column ','levgrnd'/), &
           long_name='column ch4 concentration', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CCON_CH4S') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'CCON_CH4S', &
                cptr%cmic%ccon_ch4s,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CCON_CH4S', &
              cptr%cmic%ccon_ch4s,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CCON_CH4S_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='column ch4 concentration in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CCON_CH4S_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CCON_CH4S_UNSAT', &
                cptr%cmic%ccon_ch4s_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CCON_CH4S_UNSAT', &
              cptr%cmic%ccon_ch4s_unsat,gcomm_column, switchdim=.true.)
    end if
 
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CCON_CH4S_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='column ch4 concentration in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CCON_CH4S_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CCON_CH4S_SAT', &
                cptr%cmic%ccon_ch4s_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CCON_CH4S_SAT', &
              cptr%cmic%ccon_ch4s_sat,gcomm_column, switchdim=.true.)   
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CCON_CO2S', &
           cdims=(/'column ','levgrnd'/), &
           long_name='column co2 concentration', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CCON_CO2S') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'CCON_CO2S', &
                cptr%cmic%ccon_co2s,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CCON_CO2S', &
              cptr%cmic%ccon_co2s,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CCON_CO2S_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='column co2 concentration in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CCON_CO2S_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CCON_CO2S_UNSAT', &
                cptr%cmic%ccon_co2s_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CCON_CO2S_UNSAT', &
              cptr%cmic%ccon_co2s_unsat,gcomm_column, switchdim=.true.)
    end if
 
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CCON_CO2S_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='column co2 concentration in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CCON_CO2S_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CCON_CO2S_SAT', &
                cptr%cmic%ccon_co2s_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CCON_CO2S_SAT', &
              cptr%cmic%ccon_co2s_sat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CCON_H2S', &
           cdims=(/'column ','levgrnd'/), &
           long_name='column h2 concentration', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CCON_H2S') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'CCON_H2S', &
                cptr%cmic%ccon_h2s,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CCON_H2S', &
              cptr%cmic%ccon_h2s,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CCON_H2S_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='column h2 concentration in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CCON_H2S_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CCON_H2S_UNSAT', &
                cptr%cmic%ccon_h2s_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CCON_H2S_UNSAT', &
              cptr%cmic%ccon_h2s_unsat,gcomm_column, switchdim=.true.)
    end if
 
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CCON_H2S_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='column h2 concentration in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CCON_H2S_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CCON_H2S_SAT', &
                cptr%cmic%ccon_h2s_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CCON_H2S_SAT', &
              cptr%cmic%ccon_h2s_sat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CCON_O2S', &
           cdims=(/'column ','levgrnd'/), &
           long_name='column o2 concentration', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CCON_O2S') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'CCON_O2S', &
                cptr%cmic%ccon_o2s,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CCON_O2S', &
              cptr%cmic%ccon_o2s,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CCON_O2S_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='column o2 concentration in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CCON_O2S_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CCON_O2S_UNSAT', &
                cptr%cmic%ccon_o2s_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CCON_O2S_UNSAT', &
              cptr%cmic%ccon_o2s_unsat,gcomm_column, switchdim=.true.)
    end if
 
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CCON_O2S_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='column o2 concentration in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CCON_O2S_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CCON_O2S_SAT', &
                cptr%cmic%ccon_o2s_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CCON_O2S_SAT', &
              cptr%cmic%ccon_o2s_sat,gcomm_column, switchdim=.true.)   
    end if
    
! below are for C13
 if ( use_c13 ) then
   if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CACES_PROD', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 acetate production', units='mol/m^3',switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CACES_PROD') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'C13_CACES_PROD', &
                cptr%cmicc13%caces_prod,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CACES_PROD', &
              cptr%cmicc13%caces_prod,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CACES_UNSAT_PROD', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 acetate production in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CACES_UNSAT_PROD') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CACES_UNSAT_PROD', &
                cptr%cmicc13%caces_unsat_prod,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CACES_UNSAT_PROD', &
              cptr%cmicc13%caces_unsat_prod,gcomm_column, switchdim=.true.)
    end if
 
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CACES_SAT_PROD', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 acetate production in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CACES_SAT_PROD') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CACES_SAT_PROD', &
                cptr%cmicc13%caces_sat_prod,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CACES_SAT_PROD', &
              cptr%cmicc13%caces_sat_prod,gcomm_column, switchdim=.true.) 
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CDOCS_PRE', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 previous time step DOC concentration', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CDOCS_PRE') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CDOCS_PRE', &
                cptr%cmicc13%cdocs_pre,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CDOCS_PRE', &
              cptr%cmicc13%cdocs_pre,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CDOCS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 dissolved organic carbon', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CDOCS') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'C13_CDOCS', &
                cptr%cmicc13%cdocs,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CDOCS', &
              cptr%cmicc13%cdocs,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CACES_UNSAT_PROD', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 dissolved organic carbon in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CACES_UNSAT_PROD') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CACES_UNSAT_PROD', &
                cptr%cmicc13%caces_unsat_prod,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CACES_UNSAT_PROD', &
              cptr%cmicc13%caces_unsat_prod,gcomm_column, switchdim=.true.)
    end if
 
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CACES_SAT_PROD', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 dissolved organic carbon in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CACES_SAT_PROD') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CACES_SAT_PROD', &
                cptr%cmicc13%caces_sat_prod,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CACES_SAT_PROD', &
              cptr%cmicc13%caces_sat_prod,gcomm_column, switchdim=.true.)     
    end if
       
    ! column microbial state variable - cmicbiocs
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CMICBIOCS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 microbial biomass carbon', units='molC/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CMICBIOCS') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CMICBIOCS', &
                cptr%cmicc13%cmicbiocs,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CMICBIOCS', &
              cptr%cmicc13%cmicbiocs,gcomm_column, switchdim=.true.)
    end if

    ! column microbial state variable - caces
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CACES', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 acetic acid', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CACES') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'C13_CACES', &
                cptr%cmicc13%caces,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CACES', &
              cptr%cmicc13%caces,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CACES_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 acetic acid in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CACES_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CACES_UNSAT', &
                cptr%cmicc13%caces_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CACES_UNSAT', &
              cptr%cmicc13%caces_unsat,gcomm_column, switchdim=.true.)
    end if
 
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CACES_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 acetic acid in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CACES_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CACES_SAT', &
                cptr%cmicc13%caces_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CACES_SAT', &
              cptr%cmicc13%caces_sat,gcomm_column, switchdim=.true.)
    end if
    
    ! column microbial state variable - cacebios
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CACEBIOS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 methanogenesis based on acetic acid', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CACEBIOS') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'C13_CACEBIOS', &
                cptr%cmicc13%cacebios,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CACEBIOS', &
              cptr%cmicc13%cacebios,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CACEBIOS_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 methanogenesis based on acetic acid in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CACEBIOS_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CACEBIOS_UNSAT', &
                cptr%cmicc13%cacebios_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CACEBIOS_UNSAT', &
              cptr%cmicc13%cacebios_unsat,gcomm_column, switchdim=.true.)
    end if
 
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CACEBIOS_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 methanogenesis based on acetic acid in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CACEBIOS_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CACEBIOS_SAT', &
                cptr%cmicc13%cacebios_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CACEBIOS_SAT', &
              cptr%cmicc13%cacebios_sat,gcomm_column, switchdim=.true.)
    end if
    
    ! column microbial state variable - cco2bios
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CCO2BIOS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 methanogensis based on co2', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CCO2BIOS') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'C13_CCO2BIOS', &
                cptr%cmicc13%cco2bios,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CCO2BIOS', &
              cptr%cmicc13%cco2bios,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CCO2BIOS_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 methanogensis based on co2 in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CCO2BIOS_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CCO2BIOS_UNSAT', &
                cptr%cmicc13%cco2bios_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CCO2BIOS_UNSAT', &
              cptr%cmicc13%cco2bios_unsat,gcomm_column, switchdim=.true.)
    end if
 
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CCO2BIOS_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 methanogensis based on co2 in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CCO2BIOS_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CCO2BIOS_SAT', &
                cptr%cmicc13%cco2bios_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CCO2BIOS_SAT', &
              cptr%cmicc13%cco2bios_sat,gcomm_column, switchdim=.true.)
    end if

    !~ ! column microbial state variable - caerch4bios
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CAERCH4BIOS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 aerobic oxidationo of ch4', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CAERCH4BIOS') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'C13_CAERCH4BIOS', &
                cptr%cmicc13%caerch4bios,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CAERCH4BIOS', &
              cptr%cmicc13%caerch4bios,gcomm_column, switchdim=.true.)
    end if 
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CAERCH4BIOS_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 aerobic oxidationo of ch4 in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CAERCH4BIOS_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CAERCH4BIOS_UNSAT', &
                cptr%cmicc13%caerch4bios_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CAERCH4BIOS_UNSAT', &
              cptr%cmicc13%caerch4bios_unsat,gcomm_column, switchdim=.true.)
    end if
 
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CAERCH4BIOS_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 aerobic oxidationo of ch4 in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CAERCH4BIOS_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CAERCH4BIOS_SAT', &
                cptr%cmicc13%caerch4bios_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CAERCH4BIOS_SAT', &
              cptr%cmicc13%caerch4bios_sat,gcomm_column, switchdim=.true.)
    end if
    
    !~ ! column microbial state variable - canaerch4bios   
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CANAERCH4BIOS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 biomass for microbes for anaerobic methane oxidation', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CANAERCH4BIOS') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'C13_CANAERCH4BIOS', &
                cptr%cmicc13%canaerch4bios,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CANAERCH4BIOS', &
              cptr%cmicc13%canaerch4bios,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CANAERCH4BIOS_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 biomass for microbes for anaerobic methane oxidation in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CANAERCH4BIOS_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CANAERCH4BIOS_UNSAT', &
                cptr%cmicc13%canaerch4bios_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CANAERCH4BIOS_UNSAT', &
              cptr%cmicc13%canaerch4bios_unsat,gcomm_column, switchdim=.true.)
    end if
 
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CANAERCH4BIOS_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 biomass for microbes for anaerobic methane oxidation in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CANAERCH4BIOS_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CANAERCH4BIOS_SAT', &
                cptr%cmicc13%canaerch4bios_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CANAERCH4BIOS_SAT', &
              cptr%cmicc13%canaerch4bios_sat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CCON_CH4S', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 column ch4 concentration', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CCON_CH4S') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'C13_CCON_CH4S', &
                cptr%cmicc13%ccon_ch4s,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CCON_CH4S', &
              cptr%cmicc13%ccon_ch4s,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CCON_CH4S_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 column ch4 concentration in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CCON_CH4S_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CCON_CH4S_UNSAT', &
                cptr%cmicc13%ccon_ch4s_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CCON_CH4S_UNSAT', &
              cptr%cmicc13%ccon_ch4s_unsat,gcomm_column, switchdim=.true.)
    end if
 
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CCON_CH4S_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 column ch4 concentration in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CCON_CH4S_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CCON_CH4S_SAT', &
                cptr%cmicc13%ccon_ch4s_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CCON_CH4S_SAT', &
              cptr%cmicc13%ccon_ch4s_sat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CCON_CO2S', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 column co2 concentration', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CCON_CO2S') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'C13_CCON_CO2S', &
                cptr%cmicc13%ccon_co2s,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CCON_CO2S', &
              cptr%cmicc13%ccon_co2s,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CCON_CO2S_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 column co2 concentration in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CCON_CO2S_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CCON_CO2S_UNSAT', &
                cptr%cmicc13%ccon_co2s_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CCON_CO2S_UNSAT', &
              cptr%cmicc13%ccon_co2s_unsat,gcomm_column, switchdim=.true.)
    end if
 
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C13_CCON_CO2S_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C13 column co2 concentration in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C13_CCON_CO2S_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C13_CCON_CO2S_SAT', &
                cptr%cmicc13%ccon_co2s_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C13_CCON_CO2S_SAT', &
              cptr%cmicc13%ccon_co2s_sat,gcomm_column, switchdim=.true.)
    end if

end if
! end of C13 

! below are for C14

  if ( use_c14 ) then
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CACES_PROD', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 acetate production', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CACES_PROD') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'C14_CACES_PROD', &
                cptr%cmicc14%caces_prod,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CACES_PROD', &
              cptr%cmicc14%caces_prod,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CACES_UNSAT_PROD', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 acetate production in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CACES_UNSAT_PROD') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CACES_UNSAT_PROD', &
                cptr%cmicc14%caces_unsat_prod,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CACES_UNSAT_PROD', &
              cptr%cmicc14%caces_unsat_prod,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CACES_SAT_PROD', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 acetate production in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CACES_SAT_PROD') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CACES_SAT_PROD', &
                cptr%cmicc14%caces_sat_prod,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CACES_SAT_PROD', &
              cptr%cmicc14%caces_sat_prod,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CDOCS_PRE', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 previous time step DOC concentration', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CDOCS_PRE') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CDOCS_PRE', &
                cptr%cmicc14%cdocs_pre,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CDOCS_PRE', & 
              cptr%cmicc14%cdocs_pre,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CDOCS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 dissolved organic carbon', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CDOCS') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'C14_CDOCS', &
                cptr%cmicc14%cdocs,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CDOCS', &
              cptr%cmicc14%cdocs,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CDOCS_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 dissolved organic carbon in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CDOCS_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CDOCS_UNSAT', &
                cptr%cmicc14%cdocs_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CDOCS_UNSAT', &
              cptr%cmicc14%cdocs_unsat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CDOCS_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 dissolved organic carbon in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CDOCS_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CDOCS_SAT', &
                cptr%cmicc14%cdocs_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CDOCS_SAT', &
              cptr%cmicc14%cdocs_sat,gcomm_column, switchdim=.true.)
    end if

    ! column microbial state variable - cmicbiocs
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CMICBIOCS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 microbial biomass carbon', units='molC/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CDOCS_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CMICBIOCS', &
                cptr%cmicc14%cmicbiocs,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CMICBIOCS', &
              cptr%cmicc14%cmicbiocs,gcomm_column, switchdim=.true.)
    end if
    
    ! column microbial state variable - caces
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CACES', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 acetic acid', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CACES') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'C14_CACES', &
                cptr%cmicc14%caces,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CACES', &
              cptr%cmicc14%caces,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CACES_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 acetic acid in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CACES_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CACES_UNSAT', &
                cptr%cmicc14%caces_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CACES_UNSAT', &
              cptr%cmicc14%caces_unsat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CACES_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 acetic acid in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CACES_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CACES_SAT', &
                cptr%cmicc14%caces_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CACES_SAT', &
              cptr%cmicc14%caces_sat,gcomm_column, switchdim=.true.)
    end if
    
    ! column microbial state variable - cacebios
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CACEBIOS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 methanogenesis based on acetic acid', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CACEBIOS') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'C14_CACEBIOS', &
                cptr%cmicc14%cacebios,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CACEBIOS', &
              cptr%cmicc14%cacebios,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CACEBIOS_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 methanogenesis based on acetic acid in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CACEBIOS_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CACEBIOS_UNSAT', &
                cptr%cmicc14%cacebios_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CACEBIOS_UNSAT', &
              cptr%cmicc14%cacebios_unsat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CACEBIOS_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 methanogenesis based on acetic acid in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CACEBIOS_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CACEBIOS_SAT', &
                cptr%cmicc14%cacebios_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CACEBIOS_SAT', &
              cptr%cmicc14%cacebios_sat,gcomm_column, switchdim=.true.)
    end if

    ! column microbial state variable - cco2bios
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CCO2BIOS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 methanogensis based on co2', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CCO2BIOS') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'C14_CCO2BIOS', &
                cptr%cmicc14%cco2bios,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CCO2BIOS', &
              cptr%cmicc14%cco2bios,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CCO2BIOS_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 methanogensis based on co2 in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CCO2BIOS_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CCO2BIOS_UNSAT', &
                cptr%cmicc14%cco2bios_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CCO2BIOS_UNSAT', &
              cptr%cmicc14%cco2bios_unsat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CCO2BIOS_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 methanogensis based on co2 in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CCO2BIOS_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CCO2BIOS_SAT', &
                cptr%cmicc14%cco2bios_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CCO2BIOS_SAT', &
              cptr%cmicc14%cco2bios_sat,gcomm_column, switchdim=.true.)
    end if

    !~ ! column microbial state variable - caerch4bios
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CAERCH4BIOS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 aerobic oxidationo of ch4', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CAERCH4BIOS') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'C14_CAERCH4BIOS', &
                cptr%cmicc14%caerch4bios,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CAERCH4BIOS', &
              cptr%cmicc14%caerch4bios,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CAERCH4BIOS_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 aerobic oxidationo of ch4 in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CAERCH4BIOS_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CAERCH4BIOS_UNSAT', &
                cptr%cmicc14%caerch4bios_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CAERCH4BIOS_UNSAT', &
              cptr%cmicc14%caerch4bios_unsat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CAERCH4BIOS_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 aerobic oxidationo of ch4 in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CAERCH4BIOS_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CAERCH4BIOS_SAT', &
                cptr%cmicc14%caerch4bios_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CAERCH4BIOS_SAT', &
              cptr%cmicc14%caerch4bios_sat,gcomm_column, switchdim=.true.)
    end if

    !~ ! column microbial state variable - canaerch4bios    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CANAERCH4BIOS', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 biomass for microbes for anaerobic methane oxidation', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CANAERCH4BIOS') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'C14_CANAERCH4BIOS', &
                cptr%cmicc14%canaerch4bios,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CANAERCH4BIOS', &
              cptr%cmicc14%canaerch4bios,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CANAERCH4BIOS_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 biomass for microbes for anaerobic methane oxidation in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CANAERCH4BIOS_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CANAERCH4BIOS_UNSAT', &
                cptr%cmicc14%canaerch4bios_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CANAERCH4BIOS_UNSAT', &
              cptr%cmicc14%canaerch4bios_unsat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CANAERCH4BIOS_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 biomass for microbes for anaerobic methane oxidation in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CANAERCH4BIOS_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CANAERCH4BIOS_SAT', &
                cptr%cmicc14%canaerch4bios_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CANAERCH4BIOS_SAT', &
              cptr%cmicc14%canaerch4bios_sat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CCON_CH4S', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 column ch4 concentration', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CCON_CH4S') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'C14_CCON_CH4S', &
                cptr%cmicc14%ccon_ch4s,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CCON_CH4S', &
              cptr%cmicc14%ccon_ch4s,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CCON_CH4S_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 column ch4 concentration in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CCON_CH4S_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CCON_CH4S_UNSAT', &
                cptr%cmicc14%ccon_ch4s_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CCON_CH4S_UNSAT', &
              cptr%cmicc14%ccon_ch4s_unsat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CCON_CH4S_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 column ch4 concentration in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CCON_CH4S_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CCON_CH4S_SAT', &
                cptr%cmicc14%ccon_ch4s_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CCON_CH4S_SAT', &
              cptr%cmicc14%ccon_ch4s_sat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CCON_CO2S', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 column co2 concentration', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CCON_CO2S') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else  
        call clm_readvar(ncid,'C14_CCON_CO2S', &
                cptr%cmicc14%ccon_co2s,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CCON_CO2S', &
              cptr%cmicc14%ccon_co2s,gcomm_column, switchdim=.true.)
    end if
    
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CCON_CO2S_UNSAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 column co2 concentration in unsaturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CCON_CO2S_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CCON_CO2S_UNSAT', &
                cptr%cmicc14%ccon_co2s_unsat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CCON_CO2S_UNSAT', &
              cptr%cmicc14%ccon_co2s_unsat,gcomm_column, switchdim=.true.)
    end if

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'C14_CCON_CO2S_SAT', &
           cdims=(/'column ','levgrnd'/), &
           long_name='C14 column co2 concentration in saturated fraction', units='mol/m^3', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'C14_CCON_CO2S_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'C14_CCON_CO2S_SAT', &
                cptr%cmicc14%ccon_co2s_sat,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'C14_CCON_CO2S_SAT', &
              cptr%cmicc14%ccon_co2s_sat,gcomm_column, switchdim=.true.)
    end if


end if
! end of C14

    ! pft ch4 state variable - tempavg_agnpp
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'tempavg_agnpp', &
              cdims=(/'pft '/), &
              long_name='Temp. Average AGNPP',units='gC/m^2/s')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'tempavg_agnpp') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'tempavg_agnpp',pptr%pcf%tempavg_agnpp,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'tempavg_agnpp',pptr%pcf%tempavg_agnpp,gcomm_pft)
    end if    

    ! pft ch4 state variable - tempavg_bgnpp
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'tempavg_bgnpp', &
              cdims=(/'pft '/), &
              long_name='Temp. Average BGNPP',units='gC/m^2/s')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'tempavg_bgnpp') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'tempavg_bgnpp',pptr%pcf%tempavg_bgnpp,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'tempavg_bgnpp',pptr%pcf%tempavg_bgnpp,gcomm_pft)
    end if

    ! pft ch4 state variable - annavg_agnpp
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'annavg_agnpp', &
              cdims=(/'pft '/), &
              long_name='Ann. Average AGNPP',units='gC/m^2/s')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'annavg_agnpp') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'annavg_agnpp',pptr%pcf%annavg_agnpp,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'annavg_agnpp',pptr%pcf%annavg_agnpp,gcomm_pft)
    end if

    ! pft ch4 state variable - annavg_bgnpp
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'annavg_bgnpp', &
              cdims=(/'pft '/), &
              long_name='Ann. Average BGNPP',units='gC/m^2/s')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'annavg_bgnpp') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'annavg_bgnpp',pptr%pcf%annavg_bgnpp,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'annavg_bgnpp',pptr%pcf%annavg_bgnpp,gcomm_pft)
    end if

!!! zhangq midified
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'QFLX_SURF_LAG', &
              cdims=(/'column '/), &
              long_name='time-lagged surface runoff',units='mm/s')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'QFLX_SURF_LAG') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'QFLX_SURF_LAG',cptr%cmic%qflx_surf_lag,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'QFLX_SURF_LAG',cptr%cmic%qflx_surf_lag,gcomm_column)
    end if
    
  end subroutine microbeRest

#endif

end module mod_clm_microbeRest
