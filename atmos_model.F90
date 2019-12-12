!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of fvGFS.                                       *
!*                                                                     *
!* fvGFS is free software; you can redistribute it and/or modify it    *
!* and are expected to follow the terms of the GNU General Public      *
!* License as published by the Free Software Foundation; either        *
!* version 2 of the License, or (at your option) any later version.    *
!*                                                                     *
!* fvGFS is distributed in the hope that it will be useful, but        *
!* WITHOUT ANY WARRANTY; without even the implied warranty of          *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU   *
!* General Public License for more details.                            *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************
module atmos_model_mod
!-----------------------------------------------------------------------
!<OVERVIEW>
!  Driver for the atmospheric model, contains routines to advance the
!  atmospheric model state by one time step.
!</OVERVIEW>

!<DESCRIPTION>
!     This version of atmos_model_mod has been designed around the implicit
!     version diffusion scheme of the GCM. It requires two routines to advance
!     the atmospheric model one time step into the future. These two routines
!     correspond to the down and up sweeps of the standard tridiagonal solver.
!     Most atmospheric processes (dynamics,radiation,etc.) are performed
!     in the down routine. The up routine finishes the vertical diffusion
!     and computes moisture related terms (convection,large-scale condensation,
!     and precipitation).

!     The boundary variables needed by other component models for coupling
!     are contained in a derived data type. A variable of this derived type
!     is returned when initializing the atmospheric model. It is used by other
!     routines in this module and by coupling routines. The contents of
!     this derived type should only be modified by the atmospheric model.

!</DESCRIPTION>

use mpp_mod,            only: mpp_pe, mpp_root_pe, mpp_clock_id, mpp_clock_begin
use mpp_mod,            only: mpp_clock_end, CLOCK_COMPONENT, MPP_CLOCK_SYNC
use mpp_mod,            only: FATAL, mpp_min, mpp_max, mpp_error, mpp_chksum
use mpp_domains_mod,    only: domain2d
use mpp_mod,            only: mpp_get_current_pelist_name
#ifdef INTERNAL_FILE_NML
use mpp_mod,            only: input_nml_file
#else
use fms_mod,            only: open_namelist_file
#endif
use fms_mod,            only: file_exist, error_mesg
use fms_mod,            only: close_file, write_version_number, stdlog, stdout
use fms_mod,            only: clock_flag_default
use fms_mod,            only: check_nml_error
use diag_manager_mod,   only: diag_send_complete_instant
use time_manager_mod,   only: time_type, get_time, get_date, &
                              operator(+), operator(-),real_to_time_type
use field_manager_mod,  only: MODEL_ATMOS
use tracer_manager_mod, only: get_number_tracers, get_tracer_names, &
                              get_tracer_index, NO_TRACER
use xgrid_mod,          only: grid_box_type
use atmosphere_mod,     only: atmosphere_init
use atmosphere_mod,     only: atmosphere_restart
use atmosphere_mod,     only: atmosphere_end
use atmosphere_mod,     only: atmosphere_state_update
use atmosphere_mod,     only: atmos_phys_driver_statein
use atmosphere_mod,     only: atmosphere_control_data
use atmosphere_mod,     only: atmosphere_resolution, atmosphere_domain
use atmosphere_mod,     only: atmosphere_grid_bdry, atmosphere_grid_ctr
use atmosphere_mod,     only: atmosphere_dynamics, atmosphere_diag_axes
use atmosphere_mod,     only: atmosphere_etalvls, atmosphere_hgt
!rab use atmosphere_mod,     only: atmosphere_tracer_postinit
use atmosphere_mod,     only: atmosphere_diss_est, atmosphere_nggps_diag
use atmosphere_mod,     only: atmosphere_scalar_field_halo
use atmosphere_mod,     only: atmosphere_get_bottom_layer
use atmosphere_mod,     only: set_atmosphere_pelist
use atmosphere_mod,     only: Atm, mytile
use block_control_mod,  only: block_control_type, define_blocks_packed
use DYCORE_typedefs,    only: DYCORE_data_type, DYCORE_diag_type
#ifdef CCPP
use IPD_typedefs,       only: IPD_init_type, IPD_diag_type,    &
                              IPD_restart_type, IPD_kind_phys, &
                              IPD_func0d_proc, IPD_func1d_proc
#else
use IPD_typedefs,       only: IPD_init_type, IPD_control_type, &
                              IPD_data_type, IPD_diag_type,    &
                              IPD_restart_type, IPD_kind_phys, &
                              IPD_func0d_proc, IPD_func1d_proc
#endif

#ifdef CCPP
use CCPP_data,          only: ccpp_suite,                      &
                              IPD_control => GFS_control,      &
                              IPD_data => GFS_data,            &
                              IPD_interstitial => GFS_interstitial
use IPD_driver,         only: IPD_initialize, IPD_initialize_rst
use CCPP_driver,        only: CCPP_step, non_uniform_blocks
#else
use IPD_driver,         only: IPD_initialize, IPD_initialize_rst, IPD_step
use physics_abstraction_layer, only: time_vary_step, radiation_step1, physics_step1, physics_step2
#endif

use stochastic_physics, only: init_stochastic_physics,         &
                              run_stochastic_physics
use stochastic_physics_sfc, only: run_stochastic_physics_sfc

use FV3GFS_io_mod,      only: FV3GFS_restart_read, FV3GFS_restart_write, &
                              FV3GFS_IPD_checksum,                       &
                              FV3GFS_diag_register, FV3GFS_diag_output,  &
                              DIAG_SIZE
use fv_iau_mod,         only: iau_external_data_type,getiauforcing,iau_initialize
use module_fv3_config,  only: output_1st_tstep_rst, first_kdt, nsout

!-----------------------------------------------------------------------

implicit none
private

public update_atmos_radiation_physics
public update_atmos_model_state
public update_atmos_model_dynamics
public atmos_model_init, atmos_model_end, atmos_data_type
public atmos_model_exchange_phase_1, atmos_model_exchange_phase_2
public atmos_model_restart
public get_atmos_model_ungridded_dim
public addLsmask2grid
!-----------------------------------------------------------------------

!<PUBLICTYPE >
 type atmos_data_type
     integer                       :: axes(4)            ! axis indices (returned by diag_manager) for the atmospheric grid 
                                                         ! (they correspond to the x, y, pfull, phalf axes)
     integer, pointer              :: pelist(:) =>null() ! pelist where atmosphere is running.
     integer                       :: layout(2)          ! computer task laytout
     logical                       :: regional           ! true if domain is regional
     logical                       :: nested             ! true if there is a nest
     integer                       :: mlon, mlat
     integer                       :: iau_offset         ! iau running window length
     logical                       :: pe                 ! current pe.
     real(kind=8),             pointer, dimension(:)     :: ak, bk
     real,                     pointer, dimension(:,:)   :: lon_bnd  => null() ! local longitude axis grid box corners in radians.
     real,                     pointer, dimension(:,:)   :: lat_bnd  => null() ! local latitude axis grid box corners in radians.
     real(kind=IPD_kind_phys), pointer, dimension(:,:)   :: lon      => null() ! local longitude axis grid box centers in radians.
     real(kind=IPD_kind_phys), pointer, dimension(:,:)   :: lat      => null() ! local latitude axis grid box centers in radians.
     real(kind=IPD_kind_phys), pointer, dimension(:,:)   :: dx, dy
     real(kind=8),             pointer, dimension(:,:)   :: area
     real(kind=8),             pointer, dimension(:,:,:) :: layer_hgt, level_hgt
     type(domain2d)                :: domain             ! domain decomposition
     type(time_type)               :: Time               ! current time
     type(time_type)               :: Time_step          ! atmospheric time step.
     type(time_type)               :: Time_init          ! reference time.
     type(grid_box_type)           :: grid               ! hold grid information needed for 2nd order conservative flux exchange 
     type(IPD_diag_type), pointer, dimension(:) :: Diag
 end type atmos_data_type
                                                         ! to calculate gradient on cubic sphere grid.
!</PUBLICTYPE >

integer :: fv3Clock, getClock, updClock, setupClock, radClock, physClock

!-----------------------------------------------------------------------
integer :: blocksize    = 1
logical :: chksum_debug = .false.
logical :: dycore_only  = .false.
logical :: debug        = .false.
!logical :: debug        = .true.
logical :: sync         = .false.
integer, parameter     :: maxhr = 4096
real, dimension(maxhr) :: fdiag = 0.
real                   :: fhmax=384.0, fhmaxhf=120.0, fhout=3.0, fhouthf=1.0,avg_max_length=3600.
#ifdef CCPP
namelist /atmos_model_nml/ blocksize, chksum_debug, dycore_only, debug, sync, fdiag, fhmax, fhmaxhf, fhout, fhouthf, ccpp_suite, avg_max_length
#else
namelist /atmos_model_nml/ blocksize, chksum_debug, dycore_only, debug, sync, fdiag, fhmax, fhmaxhf, fhout, fhouthf, avg_max_length
#endif

type (time_type) :: diag_time, diag_time_fhzero

!--- concurrent and decoupled radiation and physics variables
!-------------------
!  DYCORE containers
!-------------------
type(DYCORE_data_type),    allocatable :: DYCORE_Data(:)  ! number of blocks
type(DYCORE_diag_type)                 :: DYCORE_Diag(25)

!----------------
!  IPD containers
!----------------
#ifndef CCPP
type(IPD_control_type)              :: IPD_Control
type(IPD_data_type),    allocatable :: IPD_Data(:)  ! number of blocks
type(IPD_diag_type),    target      :: IPD_Diag(DIAG_SIZE)
type(IPD_restart_type)              :: IPD_Restart
#else
! IPD_Control and IPD_Data are coming from CCPP_data
type(IPD_diag_type),    target      :: IPD_Diag(DIAG_SIZE)
type(IPD_restart_type)              :: IPD_Restart
#endif

!--------------
! IAU container
!--------------
type(iau_external_data_type)        :: IAU_Data ! number of blocks

!-----------------
!  Block container
!-----------------
type (block_control_type), target   :: Atm_block

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

#ifdef NAM_phys
  logical,parameter :: flip_vc = .false.
#else
  logical,parameter :: flip_vc = .true.
#endif

  real(kind=IPD_kind_phys), parameter :: zero=0.0, one=1.0

contains

!#######################################################################
! <SUBROUTINE NAME="update_radiation_physics">
!
!<DESCRIPTION>
!   Called every time step as the atmospheric driver to compute the
!   atmospheric tendencies for dynamics, radiation, vertical diffusion of
!   momentum, tracers, and heat/moisture.  For heat/moisture only the
!   downward sweep of the tridiagonal elimination is performed, hence
!   the name "_down". 
!</DESCRIPTION>

!   <TEMPLATE>
!     call  update_atmos_radiation_physics (Atmos)
!   </TEMPLATE>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
!   These fields describe the atmospheric grid and are needed to
!   compute/exchange fluxes with other component models.  All fields in this
!   variable type are allocated for the global grid (without halo regions).
! </INOUT>

subroutine update_atmos_radiation_physics (Atmos)
#ifdef OPENMP
    use omp_lib
#endif
!-----------------------------------------------------------------------
  type (atmos_data_type), intent(in) :: Atmos
!--- local variables---
    integer :: nb, jdat(8), rc
    procedure(IPD_func0d_proc), pointer :: Func0d => NULL()
    procedure(IPD_func1d_proc), pointer :: Func1d => NULL()
    integer :: nthrds
#ifdef CCPP
    integer :: ierr
#endif

#ifdef OPENMP
    nthrds = omp_get_max_threads()
#else
    nthrds = 1
#endif

    if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "statein driver"
!--- get atmospheric state from the dynamic core
    call set_atmosphere_pelist()
    call mpp_clock_begin(getClock)
    if (IPD_control%do_skeb) call atmosphere_diss_est (IPD_control%skeb_npass) !  do smoothing for SKEB
    call atmos_phys_driver_statein (IPD_data, Atm_block, flip_vc)
    call mpp_clock_end(getClock)

!--- if dycore only run, set up the dummy physics output state as the input state
    if (dycore_only) then
      do nb = 1,Atm_block%nblks
        IPD_Data(nb)%Stateout%gu0 = IPD_Data(nb)%Statein%ugrs
        IPD_Data(nb)%Stateout%gv0 = IPD_Data(nb)%Statein%vgrs
        IPD_Data(nb)%Stateout%gt0 = IPD_Data(nb)%Statein%tgrs
        IPD_Data(nb)%Stateout%gq0 = IPD_Data(nb)%Statein%qgrs
      enddo
    else
      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "setup step"

!--- update IPD_Control%jdat(8)
      jdat(:) = 0
      call get_date (Atmos%Time, jdat(1), jdat(2), jdat(3),  &
                                 jdat(5), jdat(6), jdat(7))
      IPD_Control%jdat(:) = jdat(:)

!--- execute the IPD atmospheric setup step
      call mpp_clock_begin(setupClock)
#ifdef CCPP
      call CCPP_step (step="time_vary", nblks=Atm_block%nblks, ierr=ierr)
      if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP time_vary step failed')
#else
      Func1d => time_vary_step
      call IPD_step (IPD_Control, IPD_Data(:), IPD_Diag, IPD_Restart, IPD_func1d=Func1d)
#endif

!--- call stochastic physics pattern generation / cellular automata
    if (IPD_Control%do_sppt .OR. IPD_Control%do_shum .OR. IPD_Control%do_skeb .OR. IPD_Control%do_sfcperts) then
       call run_stochastic_physics(IPD_Control, IPD_Data(:)%Grid, IPD_Data(:)%Coupling, nthrds)
    end if

    if(IPD_Control%do_ca)then
       ! DH* The current implementation of cellular_automata assumes that all blocksizes are the
       ! same, this is tested in the initialization call to cellular_automata, no need to redo *DH
       call cellular_automata(IPD_Control%kdt, IPD_Data(:)%Statein, IPD_Data(:)%Coupling, IPD_Data(:)%Intdiag, &
                              Atm_block%nblks, IPD_Control%levs, IPD_Control%nca, IPD_Control%ncells,          &
                              IPD_Control%nlives, IPD_Control%nfracseed, IPD_Control%nseed,                    &
                              IPD_Control%nthresh, IPD_Control%ca_global, IPD_Control%ca_sgs,                  &
                              IPD_Control%iseed_ca, IPD_Control%ca_smooth, IPD_Control%nspinup,                &
                              Atm_block%blksz(1))
    endif

!--- if coupled, assign coupled fields
      if( IPD_Control%cplflx .or. IPD_Control%cplwav ) then
        call assign_importdata(rc)
      endif

      call mpp_clock_end(setupClock)

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "radiation driver"

!--- execute the IPD atmospheric radiation subcomponent (RRTM)

      call mpp_clock_begin(radClock)
#ifdef CCPP
      ! Performance improvement. Only enter if it is time to call the radiation physics.
      if (IPD_Control%lsswr .or. IPD_Control%lslwr) then
        call CCPP_step (step="radiation", nblks=Atm_block%nblks, ierr=ierr)
        if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP radiation step failed')
      endif
#else
      Func0d => radiation_step1
!$OMP parallel do default (none)       &
!$OMP            schedule (dynamic,1), &
!$OMP            shared   (Atm_block, IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, Func0d) &
!$OMP            private  (nb)
      do nb = 1,Atm_block%nblks
        call IPD_step (IPD_Control, IPD_Data(nb:nb), IPD_Diag, IPD_Restart, IPD_func0d=Func0d)
      enddo
#endif
      call mpp_clock_end(radClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'RADIATION STEP  ', IPD_Control%kdt, IPD_Control%fhour
        call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
      endif

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "physics driver"

!--- execute the IPD atmospheric physics step1 subcomponent (main physics driver)

      call mpp_clock_begin(physClock)
#ifdef CCPP
      call CCPP_step (step="physics", nblks=Atm_block%nblks, ierr=ierr)
      if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP physics step failed')
#else
      Func0d => physics_step1
!$OMP parallel do default (none) &
!$OMP            schedule (dynamic,1), &
!$OMP            shared   (Atm_block, IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, Func0d) &
!$OMP            private  (nb)
      do nb = 1,Atm_block%nblks
        call IPD_step (IPD_Control, IPD_Data(nb:nb), IPD_Diag, IPD_Restart, IPD_func0d=Func0d)
      enddo
#endif
      call mpp_clock_end(physClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'PHYSICS STEP1   ', IPD_Control%kdt, IPD_Control%fhour
        call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
      endif

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "stochastic physics driver"

!--- execute the IPD atmospheric physics step2 subcomponent (stochastic physics driver)

      call mpp_clock_begin(physClock)
#ifdef CCPP
      call CCPP_step (step="stochastics", nblks=Atm_block%nblks, ierr=ierr)
      if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP stochastics step failed')
#else
      Func0d => physics_step2
!$OMP parallel do default (none) &
!$OMP            schedule (dynamic,1), &
!$OMP            shared   (Atm_block, IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, Func0d) &
!$OMP            private  (nb)
      do nb = 1,Atm_block%nblks
        call IPD_step (IPD_Control, IPD_Data(nb:nb), IPD_Diag, IPD_Restart, IPD_func0d=Func0d)
      enddo
#endif
      call mpp_clock_end(physClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'PHYSICS STEP2   ', IPD_Control%kdt, IPD_Control%fhour
        call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
      endif
      call getiauforcing(IPD_Control,IAU_data)
      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "end of radiation and physics step"
    endif

#ifdef CCPP
    ! Update flag for first time step of time integration
    IPD_Control%first_time_step = .false.
#endif
!-----------------------------------------------------------------------
 end subroutine update_atmos_radiation_physics
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="atmos_model_init">
!
! <OVERVIEW>
! Routine to initialize the atmospheric model
! </OVERVIEW>

subroutine atmos_model_init (Atmos, Time_init, Time, Time_step)

#ifdef OPENMP
  use omp_lib
#endif
#ifdef CCPP
  use fv_mp_mod, only: commglobal
#endif
  use mpp_mod, only: mpp_npes

  type (atmos_data_type), intent(inout) :: Atmos
  type (time_type), intent(in) :: Time_init, Time, Time_step
!--- local variables ---
  integer :: unit, ntdiag, ntfamily, i, j, k
  integer :: mlon, mlat, nlon, nlat, nlev, sec, dt
  integer :: ierr, io, logunit
  integer :: idx, tile_num
  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: blk, ibs, ibe, jbs, jbe
  real(kind=IPD_kind_phys) :: dt_phys
  real, allocatable    :: q(:,:,:,:), p_half(:,:,:)
  character(len=80)    :: control
  character(len=64)    :: filename, filename2, pelist_name
  character(len=132)   :: text
  logical              :: p_hydro, hydro, fexist
  logical, save        :: block_message = .true.
  type(IPD_init_type)  :: Init_parm
  integer              :: bdat(8), cdat(8)
  integer              :: ntracers, maxhf, maxh
  character(len=32), allocatable, target :: tracer_names(:)
  integer :: nthrds

!-----------------------------------------------------------------------

!---- set the atmospheric model time ------

   Atmos % Time_init = Time_init
   Atmos % Time      = Time
   Atmos % Time_step = Time_step
   call get_time (Atmos % Time_step, sec)
   dt_phys = real(sec)      ! integer seconds

   logunit = stdlog()

!-----------------------------------------------------------------------
! initialize atmospheric model -----

#ifndef CCPP
!---------- initialize atmospheric dynamics -------
   call atmosphere_init (Atmos%Time_init, Atmos%Time, Atmos%Time_step,&
                         Atmos%grid, Atmos%area)
#endif

   IF ( file_exist('input.nml')) THEN
#ifdef INTERNAL_FILE_NML
      read(input_nml_file, nml=atmos_model_nml, iostat=io)
      ierr = check_nml_error(io, 'atmos_model_nml')
#else
      unit = open_namelist_file ( )
      ierr=1
      do while (ierr /= 0)
         read  (unit, nml=atmos_model_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'atmos_model_nml')
      enddo
 10     call close_file (unit)
#endif
   endif

#ifdef CCPP
!---------- initialize atmospheric dynamics after reading the namelist -------
!---------- (need name of CCPP suite definition file from input.nml) ---------
   call atmosphere_init (Atmos%Time_init, Atmos%Time, Atmos%Time_step,&
                         Atmos%grid, Atmos%area)
#endif

!-----------------------------------------------------------------------
   call atmosphere_resolution (nlon, nlat, global=.false.)
   call atmosphere_resolution (mlon, mlat, global=.true.)
   call alloc_atmos_data_type (nlon, nlat, Atmos)
   call atmosphere_domain (Atmos%domain, Atmos%layout, Atmos%regional, Atmos%nested, Atmos%pelist)
   call atmosphere_diag_axes (Atmos%axes)
   call atmosphere_etalvls (Atmos%ak, Atmos%bk, flip=flip_vc)
   call atmosphere_grid_bdry (Atmos%lon_bnd, Atmos%lat_bnd, global=.false.)
   call atmosphere_grid_ctr (Atmos%lon, Atmos%lat)
   call atmosphere_hgt (Atmos%layer_hgt, 'layer', relative=.false., flip=flip_vc)
   call atmosphere_hgt (Atmos%level_hgt, 'level', relative=.false., flip=flip_vc)

   Atmos%mlon = mlon
   Atmos%mlat = mlat
!-----------------------------------------------------------------------
!--- before going any further check definitions for 'blocks'
!-----------------------------------------------------------------------
   call atmosphere_control_data (isc, iec, jsc, jec, nlev, p_hydro, hydro, tile_num)
   call define_blocks_packed ('atmos_model', Atm_block, isc, iec, jsc, jec, nlev, &
                              blocksize, block_message)
   
   allocate(DYCORE_Data(Atm_block%nblks))
   allocate(IPD_Data(Atm_block%nblks))

#ifdef OPENMP
   nthrds = omp_get_max_threads()
#else
   nthrds = 1
#endif

#ifdef CCPP
   ! This logic deals with non-uniform block sizes for CCPP.
   ! When non-uniform block sizes are used, it is required
   ! that only the last block has a different (smaller)
   ! size than all other blocks. This is the standard in
   ! FV3. If this is the case, set non_uniform_blocks (a
   ! variable imported from CCPP_driver) to .true. and
   ! allocate nthreads+1 elements of the interstitial array.
   ! The extra element will be used by the thread that
   ! runs over the last, smaller block.
   if (minval(Atm_block%blksz)==maxval(Atm_block%blksz)) then
      non_uniform_blocks = .false.
      allocate(IPD_Interstitial(nthrds))
   else if (all(minloc(Atm_block%blksz)==(/size(Atm_block%blksz)/))) then
      non_uniform_blocks = .true.
      allocate(IPD_Interstitial(nthrds+1))
   else
      call mpp_error(FATAL, 'For non-uniform blocksizes, only the last element ' // &
                            'in Atm_block%blksz can be different from the others')
   end if

#endif

!--- update IPD_Control%jdat(8)
   bdat(:) = 0
   call get_date (Time_init, bdat(1), bdat(2), bdat(3),  &
                             bdat(5), bdat(6), bdat(7))
   cdat(:) = 0
   call get_date (Time,      cdat(1), cdat(2), cdat(3),  &
                             cdat(5), cdat(6), cdat(7))
   call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
   allocate (tracer_names(ntracers))
   do i = 1, ntracers
     call get_tracer_names(MODEL_ATMOS, i, tracer_names(i))
   enddo
!--- setup IPD Init_parm
   Init_parm%me              =  mpp_pe()
   Init_parm%master          =  mpp_root_pe()
   Init_parm%tile_num        =  tile_num
   Init_parm%isc             =  isc
   Init_parm%jsc             =  jsc
   Init_parm%nx              =  nlon
   Init_parm%ny              =  nlat
   Init_parm%levs            =  nlev
   Init_parm%cnx             =  mlon
   Init_parm%cny             =  mlat
   Init_parm%gnx             =  Init_parm%cnx*4
   Init_parm%gny             =  Init_parm%cny*2
   Init_parm%nlunit          =  9999
   Init_parm%logunit         =  logunit
   Init_parm%bdat(:)         =  bdat(:)
   Init_parm%cdat(:)         =  cdat(:)
   Init_parm%dt_dycore       =  dt_phys
   Init_parm%dt_phys         =  dt_phys
   Init_parm%iau_offset      =  Atmos%iau_offset
   Init_parm%blksz           => Atm_block%blksz
   Init_parm%ak              => Atmos%ak
   Init_parm%bk              => Atmos%bk
   Init_parm%xlon            => Atmos%lon
   Init_parm%xlat            => Atmos%lat
   Init_parm%area            => Atmos%area
   Init_parm%tracer_names    => tracer_names
#ifdef CCPP
   Init_parm%restart         = Atm(mytile)%flagstruct%warm_start
   Init_parm%hydrostatic     = Atm(mytile)%flagstruct%hydrostatic
#endif

#ifdef INTERNAL_FILE_NML
   Init_parm%input_nml_file  => input_nml_file
   Init_parm%fn_nml='using internal file'
#else
   pelist_name=mpp_get_current_pelist_name()
   Init_parm%fn_nml='input_'//trim(pelist_name)//'.nml'
   inquire(FILE=Init_parm%fn_nml, EXIST=fexist)
   if (.not. fexist ) then
      Init_parm%fn_nml='input.nml'
   endif
#endif

#ifdef CCPP
   call IPD_initialize (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, &
                        IPD_Interstitial, commglobal, mpp_npes(), Init_parm)
#else
   call IPD_initialize (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, Init_parm)
#endif

   if (IPD_Control%do_sppt .OR. IPD_Control%do_shum .OR. IPD_Control%do_skeb .OR. IPD_Control%do_sfcperts) then
      ! Initialize stochastic physics
      call init_stochastic_physics(IPD_Control, Init_parm, mpp_npes(), nthrds)
      if(IPD_Control%me == IPD_Control%master) print *,'do_skeb=',IPD_Control%do_skeb
   end if

#ifdef CCPP
   ! Initialize the CCPP framework
   call CCPP_step (step="init", nblks=Atm_block%nblks, ierr=ierr)
   if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP init step failed')
   ! Doing the init here requires logic in thompson aerosol init if no aerosol
   ! profiles are specified and internal profiles are calculated, because these
   ! require temperature/geopotential etc which are not yet set. Sim. for RUC LSM.
   call CCPP_step (step="physics_init", nblks=Atm_block%nblks, ierr=ierr)
   if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP physics_init step failed')
#endif

   Atmos%Diag => IPD_Diag

   if (IPD_Control%do_sfcperts) then
      ! Get land surface perturbations here (move to GFS_time_vary
      ! step if wanting to update each time-step)
      call run_stochastic_physics_sfc(IPD_Control, IPD_Data(:)%Grid, IPD_Data(:)%Coupling)
   end if

   ! Initialize cellular automata
   if(IPD_Control%do_ca)then
      ! DH* The current implementation of cellular_automata assumes that all blocksizes are the
      ! same - abort if this is not the case, otherwise proceed with Atm_block%blksz(1) below
      if (.not. minval(Atm_block%blksz)==maxval(Atm_block%blksz)) then
         call mpp_error(FATAL, 'Logic errror: cellular_automata not compatible with non-uniform blocksizes')
      end if
      ! *DH
      call cellular_automata(IPD_Control%kdt, IPD_Data(:)%Statein, IPD_Data(:)%Coupling, IPD_Data(:)%Intdiag, &
                             Atm_block%nblks, IPD_Control%levs, IPD_Control%nca, IPD_Control%ncells,          &
                             IPD_Control%nlives, IPD_Control%nfracseed, IPD_Control%nseed,                    &
                             IPD_Control%nthresh, IPD_Control%ca_global, IPD_Control%ca_sgs,                  &
                             IPD_Control%iseed_ca, IPD_Control%ca_smooth, IPD_Control%nspinup,                &
                             Atm_block%blksz(1))
   endif

   Atm(mytile)%flagstruct%do_skeb = IPD_Control%do_skeb

!  initialize the IAU module
   call iau_initialize (IPD_Control,IAU_data,Init_parm)

   Init_parm%blksz           => null()
   Init_parm%ak              => null()
   Init_parm%bk              => null()
   Init_parm%xlon            => null()
   Init_parm%xlat            => null()
   Init_parm%area            => null()
   Init_parm%tracer_names    => null()
   deallocate (tracer_names)

   !--- update tracers in FV3 with any initialized during the physics/radiation init phase
!rab   call atmosphere_tracer_postinit (IPD_Data, Atm_block)

   call atmosphere_nggps_diag (Time, init=.true.)
   call FV3GFS_diag_register (IPD_Diag, Time, Atm_block, IPD_Control, Atmos%lon, Atmos%lat, Atmos%axes)
   call IPD_initialize_rst (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, Init_parm)
#ifdef CCPP
   call FV3GFS_restart_read (IPD_Data, IPD_Restart, Atm_block, IPD_Control, Atmos%domain, Atm(mytile)%flagstruct%warm_start)
#else
   call FV3GFS_restart_read (IPD_Data, IPD_Restart, Atm_block, IPD_Control, Atmos%domain)
#endif

   !--- set the initial diagnostic timestamp
   diag_time = Time 
   if (output_1st_tstep_rst) then
     diag_time = Time - real_to_time_type(mod(int((first_kdt - 1)*dt_phys/3600.),6)*3600.0)
   endif
   if (Atmos%iau_offset > zero) then
     diag_time = Atmos%Time_init
     diag_time_fhzero = Atmos%Time
   endif

   !---- print version number to logfile ----

   call write_version_number ( version, tagname )
   !--- write the namelist to a log file
   if (mpp_pe() == mpp_root_pe()) then
      unit = stdlog( )
      write (unit, nml=atmos_model_nml)
      call close_file (unit)
   endif

   !--- get fdiag
#ifdef GFS_PHYS
!--- check fdiag to see if it is an interval or a list
   if (nint(fdiag(2)) == 0) then
     if (fhmaxhf > 0) then
       maxhf = fhmaxhf / fhouthf
       maxh  = maxhf + (fhmax-fhmaxhf) / fhout
       fdiag(1) = fhouthf
       do i=2,maxhf
        fdiag(i) = fdiag(i-1) + fhouthf
       enddo
       do i=maxhf+1,maxh
         fdiag(i) = fdiag(i-1) + fhout
       enddo
     else
       maxh  = fhmax / fhout
       do i = 2, maxh
         fdiag(i) = fdiag(i-1) + fhout
       enddo
     endif
   endif
   if (mpp_pe() == mpp_root_pe()) write(6,*) "---fdiag",fdiag(1:40)
#endif

   setupClock = mpp_clock_id( 'GFS Step Setup        ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   radClock   = mpp_clock_id( 'GFS Radiation         ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   physClock  = mpp_clock_id( 'GFS Physics           ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   getClock   = mpp_clock_id( 'Dynamics get state    ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   updClock   = mpp_clock_id( 'Dynamics update state ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   if (sync) then
     fv3Clock = mpp_clock_id( 'FV3 Dycore            ', flags=clock_flag_default+MPP_CLOCK_SYNC, grain=CLOCK_COMPONENT )
   else
     fv3Clock = mpp_clock_id( 'FV3 Dycore            ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   endif

#ifdef CCPP
   ! Set flag for first time step of time integration
   IPD_Control%first_time_step = .true.
#endif
!-----------------------------------------------------------------------
end subroutine atmos_model_init
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_dynamics"
!
! <OVERVIEW>
subroutine update_atmos_model_dynamics (Atmos)
! run the atmospheric dynamics to advect the properties
  type (atmos_data_type), intent(in) :: Atmos

    call set_atmosphere_pelist()
    call mpp_clock_begin(fv3Clock)
    call atmosphere_dynamics (Atmos%Time)
    call mpp_clock_end(fv3Clock)

end subroutine update_atmos_model_dynamics
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="atmos_model_exchange_phase_1"
!
! <OVERVIEW>
!   Perform data exchange with coupled components in run phase 1
! </OVERVIEW>
!
! <DESCRIPTION>
!  This subroutine currently exports atmospheric fields and tracers
!  to the chemistry component during the model's run phase 1, i.e.
!  before chemistry is run.
! </DESCRIPTION>

subroutine atmos_model_exchange_phase_1 (Atmos, rc)

  type (atmos_data_type), intent(inout) :: Atmos
  integer, optional,      intent(out)   :: rc
!--- local variables
  integer :: localrc

    !--- begin

    !--- if coupled, exchange coupled fields
    if( IPD_Control%cplchm ) then
      ! -- export fields to chemistry
      call update_atmos_chemistry('export', rc=localrc)
    endif

 end subroutine atmos_model_exchange_phase_1
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="atmos_model_exchange_phase_2"
!
! <OVERVIEW>
!   Perform data exchange with coupled components in run phase 2
! </OVERVIEW>
!
! <DESCRIPTION>
!  This subroutine currently imports fields updated by the coupled
!  chemistry component back into the atmospheric model during run
!  phase 2.
! </DESCRIPTION>

subroutine atmos_model_exchange_phase_2 (Atmos, rc)

  type (atmos_data_type), intent(inout) :: Atmos
  integer, optional,      intent(out)   :: rc
!--- local variables
  integer :: localrc

    !--- if coupled, exchange coupled fields
    if( IPD_Control%cplchm ) then
      ! -- import fields from chemistry
      call update_atmos_chemistry('import', rc=localrc)
    endif

 end subroutine atmos_model_exchange_phase_2
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_state"
!
! <OVERVIEW>
subroutine update_atmos_model_state (Atmos)
! to update the model state after all concurrency is completed
  type (atmos_data_type), intent(inout) :: Atmos
!--- local variables
  integer :: isec, seconds, isec_fhzero
  integer :: rc
  real(kind=IPD_kind_phys) :: time_int, time_intfull
!
    call set_atmosphere_pelist()
    call mpp_clock_begin(fv3Clock)
    call mpp_clock_begin(updClock)
    call atmosphere_state_update (Atmos%Time, IPD_Data, IAU_Data, Atm_block, flip_vc)
    call mpp_clock_end(updClock)
    call mpp_clock_end(fv3Clock)

    if (chksum_debug) then
      if (mpp_pe() == mpp_root_pe()) print *,'UPDATE STATE    ', IPD_Control%kdt, IPD_Control%fhour
      if (mpp_pe() == mpp_root_pe()) print *,'in UPDATE STATE    ', size(IPD_Data(1)%SfcProp%tsfc),'nblks=',Atm_block%nblks
      call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
    endif

    !--- advance time ---
    Atmos % Time = Atmos % Time + Atmos % Time_step

    call get_time (Atmos%Time - diag_time, isec)
    call get_time (Atmos%Time - Atmos%Time_init, seconds)
    call atmosphere_nggps_diag(Atmos%Time,ltavg=.true.,avg_max_length=avg_max_length)
    if (ANY(nint(fdiag(:)*3600.0) == seconds) .or. (IPD_Control%kdt == first_kdt) .or. nsout > 0) then
      if (mpp_pe() == mpp_root_pe()) write(6,*) "---isec,seconds",isec,seconds
      time_int = real(isec)
      if(Atmos%iau_offset > zero) then
        if( time_int - Atmos%iau_offset*3600. > zero ) then
          time_int = time_int - Atmos%iau_offset*3600.
        else if(seconds == Atmos%iau_offset*3600) then
          call get_time (Atmos%Time - diag_time_fhzero, isec_fhzero)
          time_int = real(isec_fhzero)
          if (mpp_pe() == mpp_root_pe()) write(6,*) "---iseczero",isec_fhzero
        endif
      endif
      time_intfull = real(seconds)
      if(Atmos%iau_offset > zero) then
        if( time_intfull - Atmos%iau_offset*3600. > zero) then
          time_intfull = time_intfull - Atmos%iau_offset*3600.
        endif
      endif
      if (mpp_pe() == mpp_root_pe()) write(6,*) ' gfs diags time since last bucket empty: ',time_int/3600.,'hrs'
      call atmosphere_nggps_diag(Atmos%Time)
      call FV3GFS_diag_output(Atmos%Time, IPD_DIag, Atm_block, IPD_Control%nx, IPD_Control%ny, &
                            IPD_Control%levs, 1, 1, 1.d0, time_int, time_intfull,              &
                            IPD_Control%fhswr, IPD_Control%fhlwr)
      if (nint(IPD_Control%fhzero) > 0) then 
        if (mod(isec,3600*nint(IPD_Control%fhzero)) == 0) diag_time = Atmos%Time
      else
        if (mod(isec,nint(3600*IPD_Control%fhzero)) == 0) diag_time = Atmos%Time
      endif
      call diag_send_complete_instant (Atmos%Time)
    endif

    !--- this may not be necessary once write_component is fully implemented
    !!!call diag_send_complete_extra (Atmos%Time)

    !--- get bottom layer data from dynamical core for coupling
    call atmosphere_get_bottom_layer (Atm_block, DYCORE_Data) 

    !if in coupled mode, set up coupled fields
    if (IPD_Control%cplflx .or. IPD_Control%cplwav) then
      call setup_exportdata(rc)
    endif

 end subroutine update_atmos_model_state
! </SUBROUTINE>



!#######################################################################
! <SUBROUTINE NAME="atmos_model_end">
!
! <OVERVIEW>
!  termination routine for atmospheric model
! </OVERVIEW>

! <DESCRIPTION>
!  Call once to terminate this module and any other modules used.
!  This routine writes a restart file and deallocates storage
!  used by the derived-type variable atmos_boundary_data_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_model_end (Atmos)
! </TEMPLATE>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
! </INOUT>

subroutine atmos_model_end (Atmos)
  type (atmos_data_type), intent(inout) :: Atmos
!---local variables
  integer :: idx
#ifdef CCPP
  integer :: ierr
#endif

!-----------------------------------------------------------------------
!---- termination routine for atmospheric model ----
                                              
    call atmosphere_end (Atmos % Time, Atmos%grid)
    call FV3GFS_restart_write (IPD_Data, IPD_Restart, Atm_block, &
                               IPD_Control, Atmos%domain)

#ifdef CCPP
!   Fast physics (from dynamics) are finalized in atmosphere_end above;
!   standard/slow physics (from IPD) are finalized in CCPP_step 'finalize'.
!   The CCPP framework for all cdata structures is finalized in CCPP_step 'finalize'.
    call CCPP_step (step="finalize", nblks=Atm_block%nblks, ierr=ierr)
    if (ierr/=0)  call mpp_error(FATAL, 'Call to CCPP finalize step failed')
#endif

end subroutine atmos_model_end

! </SUBROUTINE>
!#######################################################################
! <SUBROUTINE NAME="atmos_model_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine atmos_model_restart(Atmos, timestamp)
  type (atmos_data_type),   intent(inout) :: Atmos
  character(len=*),  intent(in)           :: timestamp

    call atmosphere_restart(timestamp)
    call FV3GFS_restart_write (IPD_Data, IPD_Restart, Atm_block, &
                               IPD_Control, Atmos%domain, timestamp)

end subroutine atmos_model_restart
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="get_atmos_model_ungridded_dim">
!
! <DESCRIPTION>
!  Retrieve ungridded dimensions of atmospheric model arrays
! </DESCRIPTION>

subroutine get_atmos_model_ungridded_dim(nlev, nsoillev, ntracers,     &
  num_diag_sfc_emis_flux, num_diag_down_flux, num_diag_type_down_flux, &
  num_diag_burn_emis_flux, num_diag_cmass)

  integer, optional, intent(out) :: nlev, nsoillev, ntracers,            &
    num_diag_sfc_emis_flux, num_diag_down_flux, num_diag_type_down_flux, &
    num_diag_burn_emis_flux, num_diag_cmass

  !--- number of atmospheric vertical levels
  if (present(nlev)) nlev = Atm_block%npz

  !--- number of soil levels
  if (present(nsoillev)) then
    nsoillev = 0
    if (allocated(IPD_Data)) then
      if (associated(IPD_Data(1)%Sfcprop%slc)) &
        nsoillev = size(IPD_Data(1)%Sfcprop%slc, dim=2)
    end if
  end if

  !--- total number of atmospheric tracers
  if (present(ntracers)) call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)

  !--- number of tracers used in chemistry diagnostic output
  if (present(num_diag_down_flux)) then
    num_diag_down_flux = 0
    if (associated(IPD_Data(1)%IntDiag%sedim)) &
      num_diag_down_flux = size(IPD_Data(1)%IntDiag%sedim, dim=2)
    if (present(num_diag_type_down_flux)) then
      num_diag_type_down_flux = 0
      if (associated(IPD_Data(1)%IntDiag%sedim))  &
        num_diag_type_down_flux = num_diag_type_down_flux + 1
      if (associated(IPD_Data(1)%IntDiag%drydep)) &
        num_diag_type_down_flux = num_diag_type_down_flux + 1
      if (associated(IPD_Data(1)%IntDiag%wetdpl)) &
        num_diag_type_down_flux = num_diag_type_down_flux + 1
      if (associated(IPD_Data(1)%IntDiag%wetdpc)) &
        num_diag_type_down_flux = num_diag_type_down_flux + 1
    end if
  end if

  !--- number of bins for chemistry diagnostic output
  if (present(num_diag_sfc_emis_flux)) then
    num_diag_sfc_emis_flux = 0
    if (associated(IPD_Data(1)%IntDiag%duem)) &
      num_diag_sfc_emis_flux = size(IPD_Data(1)%IntDiag%duem, dim=2)
    if (associated(IPD_Data(1)%IntDiag%ssem)) &
      num_diag_sfc_emis_flux = &
        num_diag_sfc_emis_flux + size(IPD_Data(1)%IntDiag%ssem, dim=2)
  end if

  !--- number of tracers used in emission diagnostic output
  if (present(num_diag_burn_emis_flux)) then
    num_diag_burn_emis_flux = 0
    if (associated(IPD_Data(1)%IntDiag%abem)) &
      num_diag_burn_emis_flux = size(IPD_Data(1)%IntDiag%abem, dim=2)
  end if

  !--- number of tracers used in column mass density diagnostics
  if (present(num_diag_cmass)) then
    num_diag_cmass = 0
    if (associated(IPD_Data(1)%IntDiag%aecm)) &
      num_diag_cmass = size(IPD_Data(1)%IntDiag%aecm, dim=2)
  end if

end subroutine get_atmos_model_ungridded_dim
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="update_atmos_chemistry">
! <DESCRIPTION>
!  Populate exported chemistry fields with current atmospheric state
!  data (state='export'). Update tracer concentrations for atmospheric
!  chemistry with values from chemistry component (state='import').
!  Fields should be exported/imported from/to the atmospheric state
!  after physics calculations.
!
!  NOTE: It is assumed that all the chemical tracers follow the standard
!  atmospheric tracers, which end with ozone. The order of the chemical
!  tracers must match their order in the chemistry component.
!
!  Requires:
!         IPD_Data
!         Atm_block
! </DESCRIPTION>
subroutine update_atmos_chemistry(state, rc)

  character(len=*),  intent(in)  :: state
  integer, optional, intent(out) :: rc

end subroutine update_atmos_chemistry
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="atmos_data_type_chksum">
!
! <OVERVIEW>
!  Print checksums of the various fields in the atmos_data_type.
! </OVERVIEW>

! <DESCRIPTION>
!  Routine to print checksums of the various fields in the atmos_data_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_data_type_chksum(id, timestep, atm)
! </TEMPLATE>

! <IN NAME="Atm" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields in the atmos_data_type.
! </INOUT>
!
! <IN NAME="id" TYPE="character">
!   Label to differentiate where this routine in being called from.
! </IN>
!
! <IN NAME="timestep" TYPE="integer">
!   An integer to indicate which timestep this routine is being called for.
! </IN>
!
subroutine atmos_data_type_chksum(id, timestep, atm)
type(atmos_data_type), intent(in) :: atm 
    character(len=*),  intent(in) :: id
    integer         ,  intent(in) :: timestep
    integer :: n, outunit

100 format("CHECKSUM::",A32," = ",Z20)
101 format("CHECKSUM::",A16,a,'%',a," = ",Z20)

  outunit = stdout()
  write(outunit,*) 'BEGIN CHECKSUM(Atmos_data_type):: ', id, timestep
  write(outunit,100) ' atm%lon_bnd                ', mpp_chksum(atm%lon_bnd)
  write(outunit,100) ' atm%lat_bnd                ', mpp_chksum(atm%lat_bnd)
  write(outunit,100) ' atm%lon                    ', mpp_chksum(atm%lon)
  write(outunit,100) ' atm%lat                    ', mpp_chksum(atm%lat)

end subroutine atmos_data_type_chksum

! </SUBROUTINE>

  subroutine alloc_atmos_data_type (nlon, nlat, Atmos)
   integer, intent(in) :: nlon, nlat
   type(atmos_data_type), intent(inout) :: Atmos
    allocate ( Atmos % lon_bnd  (nlon+1,nlat+1), &
               Atmos % lat_bnd  (nlon+1,nlat+1), &
               Atmos % lon      (nlon,nlat),     &
               Atmos % lat      (nlon,nlat)      )

  end subroutine alloc_atmos_data_type

  subroutine dealloc_atmos_data_type (Atmos)
   type(atmos_data_type), intent(inout) :: Atmos
    deallocate (Atmos%lon_bnd, &
                Atmos%lat_bnd, &
                Atmos%lon,     &
                Atmos%lat      )
  end subroutine dealloc_atmos_data_type

  subroutine assign_importdata(rc)

    implicit none
    integer, intent(out) :: rc

  end subroutine assign_importdata

!
  subroutine setup_exportdata (rc)

    implicit none
    integer, intent(out) :: rc

  end subroutine setup_exportdata

  subroutine addLsmask2grid(fcstgrid, rc)

    implicit none
    integer :: fcstgrid
    integer, optional, intent(out) :: rc

  end subroutine addLsmask2grid
!------------------------------------------------------------------------------

end module atmos_model_mod
