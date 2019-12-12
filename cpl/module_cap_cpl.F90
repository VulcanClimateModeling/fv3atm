module module_cap_cpl
!
!*** this module contains the debug subroutines for fv3 coupled run
!
! revision history
!  12 Mar 2018: J. Wang       Pull coupled subroutines from fv3_cap.F90 to this module
!
  implicit none
  private
  public clock_cplIntval
  public realizeConnectedInternCplField
  public realizeConnectedCplFields
  public Dump_cplFields
!
  contains

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

    subroutine clock_cplIntval(gcomp, CF)

      integer :: gcomp, CF

    end subroutine clock_cplIntval

  !-----------------------------------------------------------------------------

    subroutine realizeConnectedInternCplField(state, field, standardName, grid, rc)

      integer :: state, field, grid
      character(len=*), optional      :: standardName
      integer, intent(out), optional  :: rc

    end subroutine realizeConnectedInternCplField

  !-----------------------------------------------------------------------------

    subroutine realizeConnectedCplFields(state, grid,                                      &
                                         numLevels, numSoilLayers, numTracers,             &
                                         num_diag_sfc_emis_flux, num_diag_down_flux,       &
                                         num_diag_type_down_flux, num_diag_burn_emis_flux, &
                                         num_diag_cmass, fieldNames, fieldTypes, fieldList, rc)

      integer,            intent(inout)  :: state
      integer,                intent(in)  :: grid
      integer,                        intent(in)  :: numLevels
      integer,                        intent(in)  :: numSoilLayers
      integer,                        intent(in)  :: numTracers
      integer,                        intent(in)  :: num_diag_sfc_emis_flux
      integer,                        intent(in)  :: num_diag_down_flux
      integer,                        intent(in)  :: num_diag_type_down_flux
      integer,                        intent(in)  :: num_diag_burn_emis_flux
      integer,                        intent(in)  :: num_diag_cmass
      character(len=*), dimension(:), intent(in)  :: fieldNames
      character(len=*), dimension(:), intent(in)  :: fieldTypes
      integer, dimension(:), intent(out) :: fieldList
      integer,                        intent(out) :: rc

    end subroutine realizeConnectedCplFields

  !-----------------------------------------------------------------------------

    subroutine Dump_cplFields(gcomp, importState, exportstate, clock_fv3,    &
         statewrite_flag, timeslice)

      integer, intent(in)       :: gcomp, clock_fv3
      integer                      :: importState, exportstate
      logical, intent(in)                   :: statewrite_flag
      integer                               :: timeslice

    end subroutine Dump_cplFields

  !-----------------------------------------------------------------------------

    subroutine ESMFPP_RegridWriteState(state, fileName, timeslice, rc)

      character(len=*), intent(in)          :: fileName
      integer, intent(in)                   :: state, timeslice
      integer, intent(out)                  :: rc

    end subroutine ESMFPP_RegridWriteState

    subroutine ESMFPP_RegridWrite(inField, outGrid, regridMethod, fileName, fieldName, timeslice, rc)

      ! input arguments
      integer, intent(in)             :: inField, outGrid, regridMethod
      character(len=*), intent(in)             :: filename
      character(len=*), intent(in)             :: fieldName
      integer,          intent(in)             :: timeslice
      integer,          intent(inout)          :: rc

    end subroutine ESMFPP_RegridWrite


  !-----------------------------------------------------------------------------

end module module_cap_cpl
