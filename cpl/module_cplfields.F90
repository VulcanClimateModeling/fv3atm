module module_cplfields

  !-----------------------------------------------------------------------------
  ! This module contains the fv3 Coupling Fields: export and import
  !
  !-----------------------------------------------------------------------------

  implicit none

  private

! Export Fields ----------------------------------------
  integer,          public, parameter :: NexportFields = 71
  character(len=*), public, parameter :: exportFieldsList(NexportFields) = (/ &
       "inst_pres_interface                      ", &
       "inst_pres_levels                         ", &
       "inst_geop_interface                      ", &
       "inst_geop_levels                         ", &
       "inst_temp_levels                         ", &
       "inst_zonal_wind_levels                   ", &
       "inst_merid_wind_levels                   ", &
       "inst_omega_levels                        ", &
       "inst_tracer_mass_frac                    ", &
       "soil_type                                ", &
       "inst_pbl_height                          ", &
       "surface_cell_area                        ", &
       "inst_convective_rainfall_amount          ", &
       "inst_exchange_coefficient_heat_levels    ", &
       "inst_spec_humid_conv_tendency_levels     ", &
       "inst_friction_velocity                   ", &
       "inst_rainfall_amount                     ", &
       "inst_soil_moisture_content               ", &
       "inst_up_sensi_heat_flx                   ", &
       "inst_lwe_snow_thickness                  ", &
       "vegetation_type                          ", &
       "inst_vegetation_area_frac                ", &
       "inst_surface_roughness                   ", &
       "mean_zonal_moment_flx                    ", &
       "mean_merid_moment_flx                    ", &
       "mean_sensi_heat_flx                      ", &
       "mean_laten_heat_flx                      ", &
       "mean_down_lw_flx                         ", &
       "mean_down_sw_flx                         ", &
       "mean_prec_rate                           ", &
       "inst_zonal_moment_flx                    ", &
       "inst_merid_moment_flx                    ", &
       "inst_sensi_heat_flx                      ", &
       "inst_laten_heat_flx                      ", &
       "inst_down_lw_flx                         ", &
       "inst_down_sw_flx                         ", &
       "inst_temp_height2m                       ", &
       "inst_spec_humid_height2m                 ", &
       "inst_zonal_wind_height10m                ", &
       "inst_merid_wind_height10m                ", &
       "inst_temp_height_surface                 ", &
       "inst_pres_height_surface                 ", &
       "inst_surface_height                      ", &
       "mean_net_lw_flx                          ", &
       "mean_net_sw_flx                          ", &
       "inst_net_lw_flx                          ", &
       "inst_net_sw_flx                          ", &
       "mean_down_sw_ir_dir_flx                  ", &
       "mean_down_sw_ir_dif_flx                  ", &
       "mean_down_sw_vis_dir_flx                 ", &
       "mean_down_sw_vis_dif_flx                 ", &
       "inst_down_sw_ir_dir_flx                  ", &
       "inst_down_sw_ir_dif_flx                  ", &
       "inst_down_sw_vis_dir_flx                 ", &
       "inst_down_sw_vis_dif_flx                 ", &
       "mean_net_sw_ir_dir_flx                   ", &
       "mean_net_sw_ir_dif_flx                   ", &
       "mean_net_sw_vis_dir_flx                  ", &
       "mean_net_sw_vis_dif_flx                  ", &
       "inst_net_sw_ir_dir_flx                   ", &
       "inst_net_sw_ir_dif_flx                   ", &
       "inst_net_sw_vis_dir_flx                  ", &
       "inst_net_sw_vis_dif_flx                  ", &
       "inst_land_sea_mask                       ", &
       "inst_temp_height_lowest                  ", &
       "inst_spec_humid_height_lowest            ", &
       "inst_zonal_wind_height_lowest            ", &
       "inst_merid_wind_height_lowest            ", &
       "inst_pres_height_lowest                  ", &
       "inst_height_lowest                       ", &
       "mean_fprec_rate                          " &
  /)
  ! Field types should be provided for proper handling
  ! according to the table below:
  !  g : soil levels (3D)
  !  i : interface (3D)
  !  l : model levels (3D)
  !  s : surface (2D)
  !  t : tracers (4D)
  character(len=*), public, parameter :: exportFieldTypes(NexportFields) = (/ &
       "i","l","i","l","l","l","l","l","t", &
       "s","s","s","s","l","l","s","s","g", &
       "s","s","s","s","s","s","s","s",     &
       "s","s","s","s","s","s","s","s",     &
       "s","s","s","s","s","s","s","s",     &
       "s","s","s","s","s","s","s","s",     &
       "s","s","s","s","s","s","s","s",     &
       "s","s","s","s","s","s","s","s",     &
       "s","s","s","s","s"                  &
!      "l","l","l","l","l","l","l","s",     &
  /)
  ! Set exportFieldShare to .true. if field is provided as memory reference
  ! to coupled components
  logical, public, parameter :: exportFieldShare(NexportFields) = (/ &
       .true. ,.true. ,.true. ,.true. ,.true. , &
       .true. ,.true. ,.true. ,.true. ,.true. , &
       .true. ,.true. ,.true. ,.true. ,.true. , &
       .true. ,.true. ,.true. ,.true. ,.true. , &
       .true. ,.true. ,.true. ,.false.,.false., &
       .false.,.false.,.false.,.false.,.false., &
       .false.,.false.,.false.,.false.,.false. , &
       .true. ,.false.,.false.,.false.,.false. , &
       .true. ,.false.,.false.,.false.,.false., &
       .false.,.false.,.false.,.false.,.false., &
       .false.,.false.,.false.,.false.,.false., &
       .false.,.false.,.false.,.false.,.false., &
       .false.,.false.,.false.,.true. ,.false., &
       .false.,.false.,.false.,.false.,.false., &
       .false.                                  &
  /)
  real(kind=8), allocatable, public :: exportData(:,:,:)

! Import Fields ----------------------------------------
  integer,          public, parameter :: NimportFields = 16
  logical,          public            :: importFieldsValid(NimportFields)
  character(len=*), public, parameter :: importFieldsList(NimportFields) = (/ &
       "inst_tracer_mass_frac                  ", &
       "land_mask                              ", &
       "sea_ice_surface_temperature            ", &
       "sea_surface_temperature                ", &
       "ice_fraction                           ", &
       "mean_up_lw_flx                         ", &
       "mean_laten_heat_flx                    ", &
       "mean_sensi_heat_flx                    ", &
       "mean_zonal_moment_flx                  ", &
       "mean_merid_moment_flx                  ", &
       "mean_ice_volume                        ", &
       "mean_snow_volume                       ", &
       "inst_tracer_up_surface_flx             ", &
       "inst_tracer_down_surface_flx           ", &
       "inst_tracer_clmn_mass_dens             ", &
       "inst_tracer_anth_biom_flx              "  &
  /)
  character(len=*), public, parameter :: importFieldTypes(NimportFields) = (/ &
       "t",                                 &
       "s","s","s","s","s",                 &
       "s","s","s","s","s",                 &
       "s","u","d","c","b"                  &
  /)
  ! Set importFieldShare to .true. if field is provided as memory reference
  ! from coupled components
  logical, public, parameter :: importFieldShare(NimportFields) = (/ &
       .true. ,                                 &
       .false.,.false.,.false.,.false.,.false., &
       .false.,.false.,.false.,.false.,.false., &
       .false.,.true. ,.true. ,.true. ,.true.   &
  /)

  ! Methods
  public fillExportFields
  public queryFieldList
  public cplFieldGet

!-----------------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------------

  subroutine fillExportFields(data_a2oi, rc)
    ! Fill updated data into the export Fields.
    real(kind=8), target, intent(in)            :: data_a2oi(:,:,:)
    integer, intent(out), optional              :: rc

  end subroutine fillExportFields
!
!------------------------------------------------------------------------------
!
  integer function queryFieldList(fieldlist, fieldname, abortflag, rc)
    ! returns integer index of first found fieldname in fieldlist
    ! by default, will abort if field not found, set abortflag to false
    ! to turn off the abort.
    ! return value of < 1 means the field was not found

    character(len=*),intent(in) :: fieldlist(:)
    character(len=*),intent(in) :: fieldname
    logical, optional           :: abortflag
    integer, optional           :: rc

  end function queryFieldList
!
!------------------------------------------------------------------------------
!
  subroutine cplStateGet(state, fieldList, fieldCount, rc)

    character(len=*), intent(in)            :: state
    integer, pointer,     optional :: fieldList(:)
    integer,          intent(out), optional :: fieldCount
    integer,          intent(out), optional :: rc

  end subroutine cplStateGet


  subroutine cplFieldGet(state, name, localDe, &
                         farrayPtr2d, farrayPtr3d, farrayPtr4d, rc)

    character(len=*),   intent(in)            :: state
    character(len=*),   intent(in)            :: name
    integer,            intent(in),  optional :: localDe
    real(kind=8), pointer,     optional :: farrayPtr2d(:,:)
    real(kind=8), pointer,     optional :: farrayPtr3d(:,:,:)
    real(kind=8), pointer,     optional :: farrayPtr4d(:,:,:,:)
    integer,            intent(out), optional :: rc

  end subroutine cplFieldGet
!
!------------------------------------------------------------------------------
!
end module module_cplfields
