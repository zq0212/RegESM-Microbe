!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_che_emission
  !
  ! Chemical and aerosol surface emission
  !
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_mpmessage
  use mod_service
  use mod_dynparam
  use mod_che_common
  use mod_che_param
  use mod_che_mppio
  use mod_che_dust
  use mod_che_ncio
  use mod_che_indices
  use netcdf
  !
  implicit none

  private
  !
  public :: chem_emission , emis_tend ,chem_emission_co2, chem_emission_ch4_clm45 !!! xiaoxie 2017 !!! zhangq 2024
  !
  contains
  !
  ! SURFACE EMIOSSION flux
  !
  subroutine chem_emission(lyear,lmonth,lday,lhour)
    implicit none
    integer(ik4) :: i, j !!! zhangq add 2024
    integer(ik4) , intent(in) :: lyear , lmonth , lday , lhour
    integer(ik4) , save :: curry = -1
    integer(ik4) , save :: currm = -1
    integer(ik4) , save :: currd = -1
    integer(ik4) , save :: ifreq = -1
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'chem_emission'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! read the aerosol emission files
    !
    if ( chemsimtype(1:2) == 'DU' .or. &
         chemsimtype(1:4) == 'SSLT' .or. &
         chemsimtype(1:4) == 'DUSS' .or. &
         chemsimtype(1:4) == 'MINE' ) return

    if ( ifreq == ifrqmon ) then
      if ( curry == lyear .and. currm == lmonth ) then
        if ( myid == italk ) then
          write(stdout,*) &
            'EMISSION for  ',lyear*1000000+lmonth*10000+100*lday,' ready', &
            ' from ',curry*1000000+currm*10000+100*currd
        end if
        return
      end if
    else if ( ifreq == ifrqday ) then
      if ( curry == lyear .and. currm == lmonth .and. currd == lday ) return
    end if
    curry = lyear
    currm = lmonth
    currd = lday
    if ( myid == italk ) then
      write(stdout,*)'READ CHEM EMISSION for ',lyear*1000000+lmonth*10000+lday
    end if
    ! Also lmonth is not really necessary here, but KEEP THIS DIMENSION
    ! FOR HIGHER TEMPORAL RESOLUTION INVENTORIES
    call read_emission(ifreq,lyear,lmonth,lday,lhour,chemsrc)
    
    !!! zhangq add
    !if ( ich4 > 0 ) then
    !  do i = ici1, ici2
    !    do j = jci1 , jci2
    !      chemsrc(j,i,ich4) = chemsrc(j,i,ich4) - tempch4_em_clm_microbe(j,i)  + ch4_em_clm_microbe(j,i)/1000._rkx
    !      tempch4_em_clm_microbe(j,i) = ch4_em_clm_microbe(j,i)/1000._rkx
    !    end do
    !  end do
    !end if

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine chem_emission
  
!!! xiaoxie 2017 for reading co2 flux ---start
  subroutine chem_emission_co2(lyear,lmonth,lday,lhour,idate)
    implicit none
    !!! zhangq add
    integer(ik4) :: i, j
    integer(ik4) , intent(in) :: lyear , lmonth , lday , lhour
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) , save :: curry_x = -1
    integer(ik4) , save :: currm_x = -1
    integer(ik4) , save :: currd_x = -1
    integer(ik4) , save :: currh_x = 1
    integer(ik4) , save :: ifreq_x = -1
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'chem_emission'
    integer(ik4) , save :: idindx_x = 0
    call time_begin(subroutine_name,idindx_x)
#endif

    if ( ifreq_x == ifrqmon ) then
      if ( curry_x == lyear .and. currm_x == lmonth ) then
        if ( myid == italk ) then
          write(stdout,*) &
            'EMISSION for  ',lyear*1000000+lmonth*10000+100*lday,' ready', &
            ' from ',curry_x*1000000+currm_x*10000+100*currd_x
        end if
        return
      end if
    else if ( ifreq_x == ifrqday ) then
      if ( curry_x == lyear .and. currm_x == lmonth .and. currd_x == lday ) return
    else if ( ifreq_x == ifrqhrs) then
      if (mod(currh_x,3).ne.0) return
    end if
    curry_x = lyear
    currm_x = lmonth
    currd_x = lday
    currh_x = lhour

    if ( myid == italk ) then
      write(stdout,*) 'currh_x : ',currh_x,' lhour: ',lhour

      write(stdout,*)'READ CHEM EMISSION for ',lyear*1000000+lmonth*10000+lday
    end if
    ! Also lmonth is not really necessary here, but KEEP THIS DIMENSION
    ! FOR HIGHER TEMPORAL RESOLUTION INVENTORIES
    
    !!! zhangq add 
    if ( ich4 > 0 ) then
      do i = ici1, ici2
        do j = jci1 , jci2
          chemsrc(j,i,ich4) = chemsrc(j,i,ich4) - tempch4_em_clm_microbe(j,i) + ch4_em_clm_microbe(j,i)/1000._rkx
          tempch4_em_clm_microbe(j,i) = ch4_em_clm_microbe(j,i)/1000._rkx
        end do
      end do
    end if
 
    call read_co2flux_as_emission(ifreq_x,lyear,lmonth,lday,lhour,chemsrc,idate)

#ifdef DEBUG
    call time_end(subroutine_name,idindx_x)
#endif
  end subroutine chem_emission_co2 
!!! xiaoxie 2017 for reading co2 flux ---end  

!!! zhangq add for reading CH4 flux from CLM45Miceobe!!!CLM45Me !!! hourly !!! start
    subroutine chem_emission_ch4_clm45(lyear,lmonth,lday,lhour,idate)
    implicit none
    integer(ik4) :: i, j
    integer(ik4) , intent(in) :: lyear , lmonth , lday , lhour
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) , save :: curry_x = -1
    integer(ik4) , save :: currm_x = -1
    integer(ik4) , save :: currd_x = -1
    integer(ik4) , save :: currh_x = 1
    integer(ik4) , save :: ifreq_x = -1
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'chem_emission_ch4_clm45'
    integer(ik4) , save :: idindx_x = 0
    call time_begin(subroutine_name,idindx_x)
#endif

    if ( ifreq_x == ifrqmon ) then
      if ( curry_x == lyear .and. currm_x == lmonth ) then
        if ( myid == italk ) then
          write(stdout,*) &
            'EMISSION for  ',lyear*1000000+lmonth*10000+100*lday,' ready', &
            ' from ',curry_x*1000000+currm_x*10000+100*currd_x
        end if
        return
      end if
    else if ( ifreq_x == ifrqday ) then
      if ( curry_x == lyear .and. currm_x == lmonth .and. currd_x == lday ) return
    else if ( ifreq_x == ifrqhrs) then
      if (mod(currh_x,3).ne.0) return
    end if
    curry_x = lyear
    currm_x = lmonth
    currd_x = lday
    currh_x = lhour

    if ( myid == italk ) then
      write(stdout,*) 'currh_x : ',currh_x,' lhour: ',lhour

      write(stdout,*)'READ CHEM EMISSION for ',lyear*1000000+lmonth*10000+lday
    end if
    ! Also lmonth is not really necessary here, but KEEP THIS DIMENSION
    ! FOR HIGHER TEMPORAL RESOLUTION INVENTORIES

    !!! zhangq add 
    !if ( ich4 > 0 ) then
    !  do i = ici1, ici2
    !    do j = jci1 , jci2
    !    chemsrc(j,i,ich4) = chemsrc(j,i,ich4) + ch4_em_clm_microbe(j,i)/1000._rkx
    !    !print*,"chemsrc=>",chemsrc(j,i,ich4),"cvoc_em_clm=>",cvoc_em_clm(j,i,ich4)
    !    end do
    !  end do
    !end if

#ifdef DEBUG
    call time_end(subroutine_name,idindx_x)
#endif
  end subroutine chem_emission_ch4_clm45
!!! zhangq add for reading CH4 flux from CLM45Microbe!!!CLM45Me !!! hourly !!! start

  !
  ! Calculation of emission tendency
  !
#ifndef CLM45
  subroutine emis_tend(j,declin)
    implicit none
    integer(ik4) , intent(in) :: j
    real(rkx) , intent(in) :: declin
#else
  subroutine emis_tend(j)
    implicit none
    integer(ik4) , intent(in) :: j
#endif

    integer(ik4)  :: i , itr
    real(rkx) ::  amp   !!! xiaoxie change
    integer(ik4) :: lyear , lmonth , lday , lhour  !!! xiaoxie add
#ifndef CLM45
    real(rkx) :: daylen , fact , maxelev , dayhr
#endif
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'emis_tend'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    ! gather anthropogenic and biogenic emissions and
    ! calculate the emission tendency from emission fluxes

    !call split_idate(idate1,lyear,lmonth,lday,lhour)    !!! xiaoxie add 2018
    !if ( ino > 0 ) then
    !  amp = 1.22632*(exp(-(((lhour+8-15.041288)*3600.)**2)/(43200.**2))+0.1)
    !  do i = ici1 , ici2
    !      tmpsrc(j,i,ino)  = chemsrc(j,i,ino)
    !      chemsrc(j,i,ino) = amp*chemsrc(j,i,ino)
    !  end do
    !end if

#ifdef CLM45
    ! overwrite or add chemsrc for biogenic in case of CLM45/BVOC option
    ! test if CLM BVOC is activated and add bio component to chemsrc.

    ! isoprene ( overwrite biogenic emission from inventory  / rq
    ! this way  we  negelect the anthropog isoprene emissions )
    if ( iisop > 0 ) then
      do i = ici1, ici2
        chemsrc(j,i,iisop) = chemsrc(j,i,iisop)+ cvoc_em_clm(j,i,iisop) !dyma
      end do
    end if

    ! other BVOC : to do !!

    ! methane : add the biogenic part from CLM to anthropogenic
    ! nb the CN/LCH4 CLM flag have to be activated
    if ( ich4 > 0 ) then
      do i = ici1, ici2
        chemsrc(j,i,ich4) = chemsrc(j,i,ich4) + cvoc_em_clm(j,i,ich4)
      end do
    end if
#else
    ! Modify chemsrc for isoprene from inventory to account for
    ! simplified dirunal cycle
    if ( iisop > 0 ) then
      do i = ici1 , ici2
        dayhr = -tan(declin)*tan(cxlat(j,i)*degrad)
        if ( dayhr < -1 .or. dayhr > 1 ) then
          tmpsrc(j,i,iisop)  = chemsrc(j,i,iisop)
          chemsrc(j,i,iisop) = 0.0
        else
          daylen = d_two*acos(-tan(declin)*tan(cxlat(j,i)*degrad))*raddeg
          daylen = d_two*acos(dayhr)*raddeg
          daylen = daylen*24.0_rkx/360.0_rkx
          ! Maximum sun elevation
          maxelev = halfpi - ((cxlat(j,i)*degrad)-declin)
          fact = (halfpi-acos(czen(j,i)))/(d_two*maxelev)
          amp = 12.0_rkx*mathpi/daylen
          tmpsrc(j,i,iisop)  = chemsrc(j,i,iisop)
          chemsrc(j,i,iisop) = amp*chemsrc(j,i,iisop) * &
                              sin(mathpi*fact)*egrav/(dsigma(kz)*1.0e3_rkx)
        end if
      end do
    end if
#endif
    !
    ! add the source term to tracer tendency
    !
    if ( ichdrdepo /= 2 ) then
      do itr = 1 , ntr
         if ( chtrname(itr)(1:2) == 'DU' .or. &
              chtrname(itr)(1:4) == 'SSLT' .or. &
              chtrname(itr)(1:6) == 'POLLEN'.or. &
              any ( mine_name == chtrname(itr)(1:4) )) cycle
          do i = ici1 , ici2
            chiten(j,i,kz,itr) = chiten(j,i,kz,itr) + &
                chemsrc(j,i,itr)*egrav/(dsigma(kz)*1.0e3_rkx)
            ! diagnostic for source, cumul
            cemtrac(j,i,itr) = cemtrac(j,i,itr) + chemsrc(j,i,itr)*cfdout
            if ( ichdiag == 1 ) then
              cemisdiag(j,i,itr) = cemisdiag(j,i,itr) + &
                  chemsrc(j,i,itr)/ ( cdzq(j,i,kz)*crhob3d(j,i,kz)) * cfdout
            end if
          end do
      end do
    else if ( ichdrdepo == 2 ) then
      do itr = 1 , ntr
        if ( chtrname(itr)(1:2) == 'DU' .or. &
             chtrname(itr)(1:4) == 'SSLT' .or. &
             chtrname(itr)(1:6) == 'POLLEN'.or. &
             any ( mine_name == chtrname(itr)(1:4) )) cycle
        do i = ici1 , ici2
          ! if PBL scheme is not UW then calculate emission tendency
          if ( ibltyp /= 2 ) then
            chiten(j,i,kz,itr) = chiten(j,i,kz,itr) + &
                 chemsrc(j,i,itr)*egrav/(dsigma(kz)*1.0e3_rkx)
          end if
          ! otherwise emission is injected in the PBL scheme ( together
          ! with dry deposition) for tend calculation
          chifxuw(j,i,itr) = chifxuw(j,i,itr) +  chemsrc(j,i,itr)
          ! diagnostic for source, cumul
          cemtrac(j,i,itr) = cemtrac(j,i,itr) + chemsrc(j,i,itr)*cfdout
          if ( ichdiag == 1 ) then
            ! in this case we calculate emission tendency diagnostic, but
            ! this term will also be included in BL tendency diagnostic
            ! if UW scheme is used.
            if ( ibltyp /= 2 ) then
              cemisdiag(j,i,itr) = cemisdiag(j,i,itr) + &
                 chemsrc(j,i,itr)/ ( cdzq(j,i,kz)*crhob3d(j,i,kz)) * cfdout
            end if
          end if
        end do
      end do
    end if

    !if ( ino > 0 ) then   !!! xiaoxie 2018
    ! chemsrc(j,:,ino) = tmpsrc(j,:,ino)   
    !end if

#ifndef CLM45
! put back chemsrc to its mean value after diurnal cycle
    if ( iisop > 0 ) then
      chemsrc(j,:,iisop) = tmpsrc(j,:,iisop)
    end if
#endif

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine emis_tend
!
end module mod_che_emission
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
