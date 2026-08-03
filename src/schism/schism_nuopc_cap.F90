! This code is part of the SCHISM-ESMF interface.  It defines
! the schism component for a NUOPC coupled system
!
! @copyright (C) 2021-2023 Helmholtz-Zentrum Hereon
! @copyright (C) 2022-2023 Virginia Institute of Marine Science
! @copyright (C) 2020-2021 Helmholtz-Zentrum Geesthacht
!
! @author Carsten Lemmen <carsten.lemmen@hereon.de>
! @author Joseph Y. Zhang >jzhang@vims.edu>
!
! @license Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
! 		http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!

#define ESMF_CONTEXT  line=__LINE__,file=ESMF_FILENAME,method=ESMF_METHOD
#define ESMF_ERR_PASSTHRU msg="SCHISM subroutine call returned error"
#undef ESMF_FILENAME
#define ESMF_FILENAME "schism_nuopc_cap.F90"

#define _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(X) if (ESMF_LogFoundError(rcToCheck=localrc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)
#define _SCHISM_LOG_AND_FINALIZE_ON_ERRORS_(X) if (ESMF_LogFoundError(rcToCheck=localRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X) .or. ESMF_LogFoundError(rcToCheck=userRc, ESMF_ERR_PASSTHRU, ESMF_CONTEXT, rcToReturn=X)) call ESMF_Finalize(rc=localrc, endflag=ESMF_END_ABORT)

module schism_nuopc_cap

  use esmf
  use nuopc
  use NUOPC_Model, &
    model_routine_SS           => SetServices, &
    model_label_DataInitialize => label_DataInitialize, &
    model_label_SetClock       => label_SetClock, &
    model_label_CheckImport    => label_CheckImport, &
    model_label_Advance        => label_Advance


  use schism_bmi
  use schism_esmf_util
  use schism_nuopc_util

  implicit none

  private
  public SetServices

#ifdef USE_NUOPC_RIVER
  !> One-way NWM river forcing stub: injects a uniform constant discharge
  !> [m^3/s] into every SCHISM source element until the per-element field is
  !> connected. Overridable via the 'river_stub_q' component attribute.
  real(ESMF_KIND_R8), save :: river_stub_q = 100.0_ESMF_KIND_R8
  !> Set from 'nwm_coupling': exempts 'river_volume_flux' from
  !> SCHISM_RemoveUnconnectedFields, since the CDEPS connector finishes
  !> connecting after this cap's IPDv00 realize phase would otherwise drop it.
  logical, save :: keep_river_field = .false.
  !> Set from 'river_stub': opt-in to substitute river_stub_q when
  !> river_volume_flux is unconnected; otherwise the cap aborts rather than
  !> silently fabricating inflow.
  logical, save :: river_stub_mode = .false.
  !> Frozen-provider tripwire: tracks the connected river field's NUOPC
  !> TimeStamp across coupling windows and warns if it stops advancing.
  logical, save :: river_stamp_seen = .false.
  integer(ESMF_KIND_I4), allocatable, save :: river_prev_stamp(:)
  type(ESMF_Time), save :: river_prev_time
  integer, save :: river_stall_count = 0
  !> Consecutive non-advancing windows tolerated before warning.
  integer, parameter :: RIVER_STALL_LIMIT = 2

  !> SCHISM flux-boundary coupling: drives qthcon for every ifltype=1
  !> open-boundary segment (flux.th replacement) from the N-point import
  !> field 'river_flux_segment', one entry per segment in bctides.in scan
  !> order. dnwm carries POSITIVE = inflow; the cap negates into qthcon
  !> (qthcon<0 = inflow). PET 0 reads the field and ESMF_VMBroadcasts it so
  !> every rank writes qthcon identically (global-replicated, unlike the
  !> element-based river_volume_flux path above). Independent of
  !> river_volume_flux/nwm_coupling; v1 is inflow-only, ifltype=2 untouched.
  !>
  !> 'nwm_flux_coupling' (component attribute): activation gate, read at
  !> init. true + unconnected field = hard abort unless 'flux_stub'=true;
  !> absent = legacy flux.th behavior.
  logical, save :: keep_flux_field = .false.
  !> 'flux_stub': opt-in to the constant-discharge stub (flux_stub_q) with
  !> no NWM provider, for standalone SCHISM-side validation.
  logical, save :: flux_stub_mode = .false.
  !> 'flux_stub_q' (m3/s, POSITIVE inflow): uniform constant discharge
  !> applied to every ifltype=1 segment in flux_stub_mode.
  real(ESMF_KIND_R8), save :: flux_stub_q = 0.0_ESMF_KIND_R8
  !> Frozen-provider tripwire for the flux-segment field, tracked
  !> independently from the river fields above.
  logical, save :: flux_stamp_seen = .false.
  integer(ESMF_KIND_I4), allocatable, save :: flux_prev_stamp(:)
  type(ESMF_Time), save :: flux_prev_time
  integer, save :: flux_stall_count = 0

  !> Tiny ESMF_Mesh of nfltype dummy triangles, built once by
  !> SCHISM_FluxMeshCreate to realize 'river_flux_segment' element-located
  !> (mesh-based, unlike the gridded 1-point + ungridded-N-axis layout used
  !> elsewhere). Module-scope, never destroyed.
  type(ESMF_Mesh), save :: fluxMesh
#endif

contains

#undef ESMF_METHOD
#define ESMF_METHOD "SetServices"
subroutine SetServices(comp, rc)

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  integer(ESMF_KIND_I4)             :: localrc

  character(len=ESMF_MAXSTR), parameter :: label_InternalState = 'InternalState'
  type(type_InternalStateWrapper) :: internalState
  type(type_InternalState), pointer :: isDataPtr => null()

  rc = ESMF_SUCCESS

  call NUOPC_CompDerive(comp, model_routine_SS, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(internalState%wrap, stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_UserCompSetInternalState(comp, label_InternalState, &
    internalState, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! NUOPC automatically has an entry point for InitializeP0, so do not?
  call NUOPC_CompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
    phaseLabelList=(/"IPDv00p0"/), userRoutine=InitializeP0, rc=localrc)

  call NUOPC_CompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
    phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeAdvertise, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
    phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeRealize, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSpecialize(comp, specLabel=model_label_SetClock, &
    specRoutine=SetClock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSpecialize(comp, specLabel=model_label_DataInitialize, &
    specRoutine=DataInitialize, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSpecialize(comp, specLabel=model_label_Advance, &
    specRoutine=ModelAdvance, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

#ifdef USE_NUOPC_RIVER
  ! Override the default run-phase import time check: the direct NWM->OCN
  ! connector is not mediator-time-brokered (see CheckImportRiver).
  call NUOPC_CompSpecialize(comp, specLabel=model_label_CheckImport, &
    specRoutine=CheckImportRiver, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
#endif

  !> Do we need a specialization of Finalize, by adding a label?
  !call NUOPC_CompSpecialize(comp, specLabel=model_label_Finalize, &
  !  specRoutine=Finalize, rc=localrc)
  ! Yes, we should release the haloHandle
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine SetServices

#undef ESMF_METHOD
#define ESMF_METHOD "InitializeP0"
subroutine InitializeP0(comp, importState, exportState, parentClock, rc)

!  Phase 0: sets the NUOPC IPDv00 InitializePhaseMap. p1 Advertise / p2 Realize
!  fields (implementor-provided); p3 check Connected status / p4 data-init +
!  timestamp exports (NUOPC-provided).

  type(ESMF_GridComp)   :: comp
  type(ESMF_State)      :: importState
  type(ESMF_State)      :: exportState
  type(ESMF_Clock)      :: parentClock
  integer, intent(out)  :: rc

  character(len=10)           :: InitializePhaseMap(4)
  integer                     :: localrc
  logical                     :: isPresent
  character(len=ESMF_MAXSTR)  :: configFileName, compName
  character(len=ESMF_MAXSTR)  :: message
  type(ESMF_Config)           :: config
  type(ESMF_Info)             :: info

  rc=ESMF_SUCCESS

  call NUOPC_CompGet(comp, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_CompSetClock(comp, externalClock=parentClock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A)') trim(compName)//' initializing (p=0) component ...'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  !> @todo add other IPD versions mappings, i.e. IPDv00p1, IPDv01p1, IPDv02p1,
  !> IPDv03p1, IPDv04p1, IPDv05p1; all map to 1 ?
  InitializePhaseMap = (/"IPDv00p1=1","IPDv00p2=2", &
    "IPDv00p3=3","IPDv00p4=4"/)

  !call ESMF_AttributeAdd(comp, convention="NUOPC", &
  !  purpose="General", &
  !  attrList=(/"InitializePhaseMap"/), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_InfoGetFromHost(comp, info=info, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_InfoSet(info, key="NUOPC/General/InitializePhaseMap", &
    values=InitializePhaseMap, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Read the configuration for this component from file if not
  ! already present in the component
  call ESMF_GridCompGet(comp, configIsPresent=isPresent, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (isPresent) then
    call ESMF_GridCompGet(comp, config=config, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    write(message, '(A)') trim(compName)//' uses internal configuration'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  else
    configfilename=trim(compName)//'.cfg'
    inquire(file=trim(configfilename), exist=isPresent)

    config = ESMF_ConfigCreate(rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    if (isPresent) then
      call ESMF_ConfigLoadFile(config, trim(configfilename), rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

      write(message,'(A)')  trim(compName)//' read configuration from '// trim(configFileName)
      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
    else
      write(message,'(A)')  trim(compName)//' has no configuration; use global config'
      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
    endif

    call ESMF_GridCompSet(comp, config=config, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

end subroutine InitializeP0

#undef ESMF_METHOD
#define ESMF_METHOD "InitializeAdvertise"
!> @description Advertises import/export fields by standard name only
!> (no memory allocated yet), for NUOPC to match against other components.
subroutine InitializeAdvertise(comp, importState, exportState, clock, rc)

#ifdef USE_NUOPC_RIVER
  use schism_glbl, only: nuopc_flux_active, nfltype, ifltype, nope_global
#endif

  implicit none

  type(ESMF_GridComp)  :: comp
  type(ESMF_State)     :: importState, exportState
  type(ESMF_Clock)     :: clock
  integer, intent(out) :: rc

  integer(ESMF_KIND_I4)       :: localrc, mpiCommunicator, mpiCommDuplicate
  integer(ESMF_KIND_I4)       :: ntime=0, iths=0, i, j, k
  character(len=ESMF_MAXSTR)  :: message, compName, cvalue
  character(len=ESMF_MAXSTR), allocatable :: itemNameList(:)
  logical                     :: isPresent, isSet

  type(ESMF_VM)                     :: vm
  type(type_InternalStateWrapper)   :: internalState
  type(type_InternalState), pointer :: isDataPtr => null()

  rc = ESMF_SUCCESS
  localrc = ESMF_SUCCESS

  !> more possikeywordds lities for interface: verbosity, profiling, diagnostic
  call NUOPC_CompGet(comp, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A)') trim(compName)//' initialize (p=1) component ...'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  allocate(internalState%wrap, stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompSetInternalState(comp, internalState, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_GridCompGetInternalState(comp, internalState, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  isDataPtr => internalState%wrap
  isDataPtr%numOwnedNodes = 0
  isDataPtr%numForeignNodes = 0

  call SCHISM_InitializePtrMap(comp, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (.not.ESMF_StateIsCreated(importState)) then
    importState=ESMF_StateCreate(name=trim(compName)//'Import', stateintent= &
      ESMF_STATEINTENT_IMPORT, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    write(message,'(A)') trim(compName)//' created state "'//trim(compName)// &
      'Import" for import'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  endif

  if (.not.ESMF_StateIsCreated(exportState)) then
    exportState=ESMF_StateCreate(name=trim(compName)//'Export', stateintent= &
      ESMF_STATEINTENT_EXPORT, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    write(message,'(A)') trim(compName)//' created state "'//trim(compName)// &
      'Export" for export'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  endif

  ! Get VM for this component
  call ESMF_GridCompGet(comp, vm=vm, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_VMGet(vm, mpiCommunicator=mpiCommunicator, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

#ifndef ESMF_MPIUNI
  call MPI_Comm_dup(mpiCommunicator, mpiCommDuplicate, rc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  write(message, '(A)') trim(compName)//' initializing parallel environment ...'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  call schism_parallel_init(communicator=mpiCommDuplicate)

  write(message, '(A)') trim(compName)//' initialized parallel environment.'
#endif

  call ESMF_UtilIOMkDir ('./outputs',  relaxedFlag=.true., rc=localrc)
  write(message, '(A)') trim(compName)//' writes results to directory "./outputs".'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  write(message, '(A)') trim(compName)//' initializing science model ...'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  inquire(file='param.nml', exist=isPresent)
  if (.not.isPresent) then
    write(message, '(A)') trim(compName)//' cannot start without required file "param.nml".'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_ERROR)
    localrc = ESMF_RC_FILE_OPEN
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  ! query attributes
  call NUOPC_CompAttributeGet(comp, name='meshloc', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  if (isPresent .and. isSet) then
    if (trim(cvalue) == 'node') then
      meshloc = ESMF_MESHLOC_NODE
    else
      meshloc = ESMF_MESHLOC_ELEMENT
    end if
  else
    cvalue = 'node'
    meshloc = ESMF_MESHLOC_NODE
  end if
  write(message, '(A)') trim(compName)//' meshloc is set to '//trim(cvalue)
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  ! debug option
  call NUOPC_CompAttributeGet(comp, name='debug_level', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  debug_level = 0
  if (isPresent .and. isSet) then
     read(cvalue,*) debug_level
  end if
  write(message, '(A,I1)') trim(compName)//' debug_level is set to ', debug_level
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

#ifdef USE_NUOPC_RIVER
  ! constant-discharge stub [m^3/s] (see river_stub_q declaration above)
  call NUOPC_CompAttributeGet(comp, name='river_stub_q', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  if (isPresent .and. isSet) then
     read(cvalue,*) river_stub_q
  end if
  write(message, '(A,F12.3)') trim(compName)//' river_stub_q [m^3/s] = ', river_stub_q
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  ! see keep_river_field declaration above
  call NUOPC_CompAttributeGet(comp, name='nwm_coupling', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  if (isPresent .and. isSet) then
     if (trim(cvalue)=='true' .or. trim(cvalue)=='.true.' .or. trim(cvalue)=='T' .or. trim(cvalue)=='TRUE') then
        keep_river_field = .true.
     end if
  end if
  write(message, '(A,L1)') trim(compName)//' nwm_coupling (keep river_volume_flux) = ', keep_river_field
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  ! see river_stub_mode declaration above
  call NUOPC_CompAttributeGet(comp, name='river_stub', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  if (isPresent .and. isSet) then
     if (trim(cvalue)=='true' .or. trim(cvalue)=='.true.' .or. trim(cvalue)=='T' .or. trim(cvalue)=='TRUE') then
        river_stub_mode = .true.
     end if
  end if
  write(message, '(A,L1)') trim(compName)//' river_stub (explicit stub mode) = ', river_stub_mode
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
#endif

#ifdef USE_NUOPC_RIVER
  ! Decide nuopc_flux_active before schism_init below, which needs it set to
  ! gate the flux.th cold-start read; a promise backed by a hard abort if unmet.
  call NUOPC_CompAttributeGet(comp, name='flux_stub_q', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  if (isPresent .and. isSet) then
     read(cvalue,*) flux_stub_q
  end if
  write(message, '(A,F12.3)') trim(compName)//' flux_stub_q [m^3/s] = ', flux_stub_q
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  ! see keep_flux_field declaration above
  call NUOPC_CompAttributeGet(comp, name='nwm_flux_coupling', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  if (isPresent .and. isSet) then
     if (trim(cvalue)=='true' .or. trim(cvalue)=='.true.' .or. trim(cvalue)=='T' .or. trim(cvalue)=='TRUE') then
        keep_flux_field = .true.
     end if
  end if
  write(message, '(A,L1)') trim(compName)//' nwm_flux_coupling (keep river_flux_segment) = ', keep_flux_field
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  ! see flux_stub_mode declaration above
  call NUOPC_CompAttributeGet(comp, name='flux_stub', value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  if (isPresent .and. isSet) then
     if (trim(cvalue)=='true' .or. trim(cvalue)=='.true.' .or. trim(cvalue)=='T' .or. trim(cvalue)=='TRUE') then
        flux_stub_mode = .true.
     end if
  end if
  write(message, '(A,L1)') trim(compName)//' flux_stub (explicit stub mode) = ', flux_stub_mode
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  ! legacy flux.th behavior unless a provider or stub was requested above
  nuopc_flux_active = keep_flux_field .or. flux_stub_mode
  write(message, '(A,L1)') trim(compName)//' nuopc_flux_active (skip native flux.th) = ', nuopc_flux_active
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
#endif

  ! init schism
  call schism_init(0, './', iths, ntime)
  write(message, '(A)') trim(compName)//' initialized science model'
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

#ifdef USE_NUOPC_RIVER
  ! Fail-fast: activating flux-boundary coupling with no ifltype=1 segments
  ! (nfltype known only now that bctides.in has been parsed) is a config error.
  if (nuopc_flux_active .and. nfltype==0) then
    write(message, '(A)') trim(compName)//': flux-boundary coupling requested '// &
      '(nwm_flux_coupling/flux_stub) but bctides.in has no ifltype=1 segments (nfltype=0).'
    call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=trim(message), ESMF_CONTEXT, rcToReturn=rc)
    return
  end if

  ! Loud logging of the segment j -> global open-boundary k mapping (EXACTLY the
  ! flux.th column order) and the sign convention, once at init.
  if (nuopc_flux_active) then
    write(message, '(A)') trim(compName)//' river_flux_segment sign convention: '// &
      'positive stream inflow -> qthcon negative (SCHISM outward-normal convention).'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
    j = 0
    do k = 1, nope_global
      if (ifltype(k) == 1) then
        j = j + 1
        write(message, '(A,I0,A,I0)') trim(compName)//' river_flux_segment segment j=', j, &
          ' -> global open boundary k=', k
        call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
      end if
    end do
  end if
#endif

  ! Set up the field dictionary and advertise the variables.
  ! @todo customize via a field dictionary-like configuration file

  call NUOPC_FieldDictionaryAddIfNeeded("air_pressure_at_sea_level", "N m-2", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("surface_downwelling_photosynthetic_radiative_flux", "W m-2 s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("sea_surface_temperature", "K", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("inst_temp_height2m", "K", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("sea_surface_salinity", "PSU", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("inst_merid_wind_height10m", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("inst_zonal_wind_height10m", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("eastward_wave_radiation_stress", "N m-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("eastward_northward_wave_radiation_stress", "N m-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("northward_wave_radiation_stress", "N m-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("depth-averaged_x-velocity", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("depth-averaged_y-velocity", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("ocn_current_zonal", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("ocn_current_merid", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("sea_surface_height_above_sea_level", "m", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("ocean_mask", "1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("sea_surface_slope_zonal", "1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("sea_surface_slope_merid", "1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_FieldDictionaryAddIfNeeded("mixed_layer_depth", "1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !  Adding CICE fields to ufs_fd
  !  Zonal ice velocity [m/s] --------------------------------------
  call NUOPC_FieldDictionaryAddIfNeeded("Si_uvel", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !  Merid ice velocity [m/s] --------------------------------------
  call NUOPC_FieldDictionaryAddIfNeeded("Si_vvel", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !  Zonal ice-to-ocean stress [N/m/m] -------------------------------
  call NUOPC_FieldDictionaryAddIfNeeded("Fioi_taux", "N m-2", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !  Merid ice-to-ocean stress [N/m/m] -------------------------------
  call NUOPC_FieldDictionaryAddIfNeeded("Fioi_tauy", "N m-2", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !  Snow volume (per unit area) [m] ----------------------------
  call NUOPC_FieldDictionaryAddIfNeeded("Si_vsno", "m", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !  Ice volume (per unit area) [m] ---------- ------------------
  call NUOPC_FieldDictionaryAddIfNeeded("Si_vice", "m", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !  Fresh water flux due to ice melting [kg/m/m/s] ------------------------
  call NUOPC_FieldDictionaryAddIfNeeded("Fioi_meltw", "kg m-2 s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
 
  !  Salt flux due to ice formation/melt [kg/m/m/s] -----------------------
  call NUOPC_FieldDictionaryAddIfNeeded("Fioi_salt", "kg m-2 s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !  Heat flux at the base of the ice [W/m/m] --------------------------
  call NUOPC_FieldDictionaryAddIfNeeded("Fioi_melth", "W m-2", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !  Shortwave radiation through sea ice coverage [W/m/m] --------------
  call NUOPC_FieldDictionaryAddIfNeeded("Fioi_swpen", "W m-2", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  
  !  Freezing melting potential (energy for ice formation) [W/m/m] ---
  call NUOPC_FieldDictionaryAddIfNeeded("Si_frzmlt", "W m-2", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !  Ice-to-ocean drage coef [Unitless] -------------------------
  call NUOPC_FieldDictionaryAddIfNeeded("Si_CdnIO", "1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Advertizing CICE states
  ! Sea surface height-------------------------------------------------------------- 
  call NUOPC_Advertise(importState, "sea_surface_height_above_sea_level", rc=localrc)
  
  !  Zonal/Merid ice velocity -----------------------------
  call NUOPC_Advertise(importState, "Si_uvel", rc=localrc)
  call NUOPC_Advertise(importState, "Si_vvel", rc=localrc)

  !  Zonal/Merid ice-to-ocn stress ------------------------
  call NUOPC_Advertise(importState, "Fioi_taux", rc=localrc)
  call NUOPC_Advertise(importState, "Fioi_tauy", rc=localrc)
 
  !  Snow volume ----------------------------------------
  call NUOPC_Advertise(importState, "Si_vsno", rc=localrc)

  !  Ice volume -----------------------------------------
  call NUOPC_Advertise(importState, "Si_vice", rc=localrc)

  !  Ice fraction ----------------------------------------
  call NUOPC_Advertise(importState, "Si_ifrac", rc=localrc)

  !  Melt water flux ---------------------------------------
  call NUOPC_Advertise(importState, "Fioi_meltw", rc=localrc)

  !  Salt flux --------------------------------------------
  call NUOPC_Advertise(importState, "Fioi_salt", rc=localrc)

  !  Heat flux ice to ocn ----------------------------------
  call NUOPC_Advertise(importState, "Fioi_melth", rc=localrc)

  !  Pen shortwave rad through ice -------------------------
  call NUOPC_Advertise(importState, "Fioi_swpen", rc=localrc)
  
  !  freezing melting potential  --------------------------
  call NUOPC_Advertise(importState, "Si_frzmlt", rc=localrc)

  !  freezing melting potential  -------------------------
  call NUOPC_Advertise(importState, "Si_CdnIO", rc=localrc)

  !  ice thickness potential  ------------------------- 
  call NUOPC_Advertise(importState, "Si_hi", rc=localrc)
  
  ! for coupling to ATM/DATM
  call NUOPC_Advertise(importState, "air_pressure_at_sea_level", rc=localrc)
  call NUOPC_Advertise(importState, "inst_zonal_wind_height10m", rc=localrc)
  call NUOPC_Advertise(importState, "inst_merid_wind_height10m", rc=localrc)
  call NUOPC_Advertise(importState, "inst_temp_height2m", rc=localrc)
  call NUOPC_Advertise(importState, "inst_spec_humid_height2m", rc=localrc)
  call NUOPC_Advertise(importState, "inst_net_sw_flx", rc=localrc)
  call NUOPC_Advertise(importState, "inst_down_lw_flx", rc=localrc)
  call NUOPC_Advertise(importState, "inst_prec_rate", rc=localrc)

  ! for coupling to WW3/WDAT
  call NUOPC_Advertise(importState, "sea_surface_wave_significant_height", rc=localrc)
  call NUOPC_Advertise(importState, "sea_water_waves_effect_on_currents_bernoulli_head_adjustment", rc=localrc)
  call NUOPC_Advertise(importState, "sea_surface_x_stress_due_to_waves", rc=localrc)
  call NUOPC_Advertise(importState, "sea_surface_y_stress_due_to_waves", rc=localrc)
  call NUOPC_Advertise(importState, "sea_bottom_upward_x_stress_due_to_waves", rc=localrc)
  call NUOPC_Advertise(importState, "sea_bottom_upward_y_stress_due_to_waves", rc=localrc)
  call NUOPC_Advertise(importState, "sea_bed_orbital_x_velocity_due_to_waves", rc=localrc)
  call NUOPC_Advertise(importState, "sea_bed_orbital_y_velocity_due_to_waves", rc=localrc)
  call NUOPC_Advertise(importState, "sea_surface_wave_mean_direction", rc=localrc)
  call NUOPC_Advertise(importState, "sea_surface_wave_mean_period", rc=localrc)
  call NUOPC_Advertise(importState, "sea_surface_wave_mean_number", rc=localrc)
  call NUOPC_Advertise(importState, "eastward_surface_stokes_drift_current", rc=localrc)
  call NUOPC_Advertise(importState, "northward_surface_stokes_drift_current", rc=localrc)
  call NUOPC_Advertise(importState, "eastward_wave_radiation_stress", rc=localrc)
  call NUOPC_Advertise(importState, "eastward_northward_wave_radiation_stress", rc=localrc)
  call NUOPC_Advertise(importState, "northward_wave_radiation_stress", rc=localrc)

#ifdef USE_NUOPC_RIVER
  ! for one-way NWM river forcing (from the dnwm CDEPS component): volumetric
  ! discharge [m3 s-1] per source element. Not a CF/standard field, so register
  ! it in the dictionary first (same helper used for the mesh fields below).
  call NUOPC_FieldDictionaryAddIfNeeded("river_volume_flux", "m3 s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_Advertise(importState, "river_volume_flux", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! see the flux-boundary module doc above; advertised unconditionally,
  ! InitializeRealize only gives it real geometry when nfltype>0.
  call NUOPC_FieldDictionaryAddIfNeeded("river_flux_segment", "m3 s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call NUOPC_Advertise(importState, "river_flux_segment", rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
#endif

  !> The mesh information is usually not in CF standard and therefore needs
  !> to be added to the FieldDictionary before advertising
  allocate(itemNameList(4))
  itemNameList=(/ 'mesh_topology                 ', &
                  'mesh_global_node_id           ', &
                  'mesh_global_element_id        ', &
                  'mesh_element_node_connectivity'/)

  do i=1,4
    call NUOPC_FieldDictionaryAddIfNeeded(trim(itemNameList(i)), "1", localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call NUOPC_Advertise(exportState, StandardName=trim(itemNameList(i)), &
      SharePolicyField='share', SharePolicyGeomObject='share', &
      TransferOfferGeomObject='will provide', rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  enddo

  ! for coupling to WW3
  call NUOPC_FieldAdvertise(exportState, "ocn_current_zonal", "m s-1", localrc)
  call NUOPC_FieldAdvertise(exportState, "ocn_current_merid", "m s-1", localrc)

  call NUOPC_FieldAdvertise(exportState, "sea_surface_temperature", "K", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FieldAdvertise(exportState, "temperature", "K", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FieldAdvertise(exportState, "sea_surface_salinity", "PSU", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FieldAdvertise(exportState, "depth-averaged_x-velocity", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FieldAdvertise(exportState, "depth-averaged_y-velocity", "m s-1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FieldAdvertise(exportState, "sea_surface_height_above_sea_level", "m", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FieldAdvertise(exportState, "sea_surface_slope_zonal", "1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FieldAdvertise(exportState, "sea_surface_slope_merid", "1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call NUOPC_FieldAdvertise(exportState, "mixed_layer_depth", "m", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! required for coupling through CMEPS mediator
  call NUOPC_FieldAdvertise(exportState, "ocean_mask", "1", localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine

#undef ESMF_METHOD
#define ESMF_METHOD "InitializeRealize"
!> @description Realizes the advertised fields as ESMF_Fields on their
!> Grid/Mesh and adds them to the import/export ESMF_States.
subroutine InitializeRealize(comp, importState, exportState, clock, rc)

  use schism_esmf_util, only : SCHISM_MeshCreateNode
  use schism_esmf_util, only : SCHISM_MeshCreateElement
  
  !> @todo move all use statements of schism into schism_bmi
  use schism_glbl, only: np, pr2, windx2, windy2, srad, nws, rkind, npa, nfltype
  use schism_esmf_util, only: SCHISM_StateFieldCreateRealize
  implicit none

  type(ESMF_GridComp)  :: comp
  type(ESMF_State)     :: importState, exportState
  type(ESMF_Clock)     :: clock
  integer, intent(out) :: rc

  type(ESMF_TimeInterval) :: stabilityTimeStep
  type(ESMF_Field)        :: field
  integer(ESMF_KIND_I4)   :: localrc, i

  type(ESMF_CoordSys_Flag) :: coordsys
  type(ESMF_Mesh)          :: mesh2d

  character(len=ESMF_MAXSTR)              :: message, compName
  character(len=ESMF_MAXSTR), allocatable :: itemNameList(:)
  type(ESMF_StateItem_Flag), allocatable  :: itemTypeList(:)
  integer(ESMF_KIND_I4)                   :: itemCount

  type(type_InternalStateWrapper)    :: internalState
  type(type_InternalState), pointer  :: isDataPtr => null()
  type(ESMF_DistGrid)                :: nodalDistgrid
  type(ESMF_Array)                   :: array

  real(ESMF_KIND_R8), pointer :: farrayPtr1(:) => null()

  ! Geometry for the 'river_flux_segment' N-point import (a tiny element-located
  ! mesh, NOT the SCHISM mesh above -- see the realize block below and
  ! SCHISM_FluxMeshCreate).
  real(ESMF_KIND_R8), pointer :: fluxFarrayPtr(:) => null()
  integer(ESMF_KIND_I4)       :: fluxLocalPet

  !> @todo move to internal state
  !> Maybe think more generally about how to handle these intermediate varialbes, 
  !> best of course dealt with in a mediator
  real(ESMF_KIND_R8), target, allocatable :: eastward_wave_radiation_stress(:)
  real(ESMF_KIND_R8), target, allocatable :: eastward_northward_wave_radiation_stress(:)
  real(ESMF_KIND_R8), target, allocatable :: northward_wave_radiation_stress(:)

  rc = ESMF_SUCCESS
  localrc= ESMF_SUCCESS

  if (meshloc == ESMF_MESHLOC_NODE) then
    call SCHISM_MeshCreateNode(comp, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  else
    call SCHISM_MeshCreateElement(comp, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  end if

  call ESMF_GridCompGet(comp, mesh=mesh2d, name=compName, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_MeshGet(mesh2d, nodalDistgrid=nodalDistgrid, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="inst_zonal_wind_height10m", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="inst_merid_wind_height10m", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="air_pressure_at_sea_level", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="downwelling_short_photosynthetic_radiation_at_water_surface", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_surface_wave_significant_height", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_water_waves_effect_on_currents_bernoulli_head_adjustment", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_surface_x_stress_due_to_waves", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_surface_y_stress_due_to_waves", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_bottom_upward_x_stress_due_to_waves", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_bottom_upward_y_stress_due_to_waves", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_bed_orbital_x_velocity_due_to_waves", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_bed_orbital_y_velocity_due_to_waves", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_surface_wave_mean_direction", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_surface_wave_mean_period", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="sea_surface_wave_mean_number", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="eastward_surface_stokes_drift_current", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="northward_surface_stokes_drift_current", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> @todo add more atmospheric fields (like humidity)
  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="inst_temp_height2m", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="inst_spec_humid_height2m", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="inst_net_sw_flx", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="inst_down_lw_flx", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="inst_prec_rate", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !  Adding CICE vars to import state 
  !> Ice_velocity ---------------------------------------------
  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="Si_uvel", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="Si_vvel", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Ice_to_ocean_stress --------------------------------------
  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="Fioi_taux", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="Fioi_tauy", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Volume of snow -------------------------------------------
  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="Si_vsno", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Volume of ice --------------------------------------------
  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="Si_vice", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> ---- ice_fraction ----------------------------------------
  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="Si_ifrac", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
 
  !> Melt_water_flux ------------------------------------------
  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="Fioi_meltw", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Salinity_flux --------------------------------------------
  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="Fioi_salt", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Heat_flux_ice_to_ocn -------------------------------------
  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="Fioi_melth", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Pen_shortwave_rad_through_ice ----------------------------
  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="Fioi_swpen", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> freezing melting potential  ------------------------------
  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="Si_frzmlt", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  
   
  !> Ice ocean drag coeff -------------------------------------
  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="Si_CdnIO", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)



  !> Wave parameters, for now we only have those from the WW3DATA cap in 
  !> NOAA's CoastalApp.  
  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="eastward_wave_radiation_stress", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="eastward_northward_wave_radiation_stress", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="northward_wave_radiation_stress", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

#ifdef USE_NUOPC_RIVER
  ! Discharge is injected at element-based sources, so fail loudly if the
  ! component-wide meshloc is not element-based, rather than mis-index later.
  if (meshloc /= ESMF_MESHLOC_ELEMENT) then
    call ESMF_LogWrite('USE_NUOPC_RIVER requires the OCN (SCHISM) component meshloc=element', &
      ESMF_LOGMSG_ERROR)
    localrc = ESMF_RC_NOT_VALID
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  end if
  ! If no provider connects, SCHISM_RemoveUnconnectedFields drops this and
  ! SCHISM_ImportRiver falls back to the constant river_stub_q.
  call SCHISM_StateFieldCreateRealize(comp, state=importState, &
    name="river_volume_flux", field=field, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  ! NOTE: the provider (dnwm cap) stamps its export with the component clock,
  ! so the "not at current time" run-phase check is satisfied without an
  ! import-side toggle here.

  ! Realize 'river_flux_segment' on a dedicated dummy-triangle mesh, not the
  ! SCHISM mesh above (see SCHISM_FluxMeshCreate). Realized regardless of
  ! activation so a late connector can still complete; skipped if nfltype==0.
  if (nfltype > 0) then
    call SCHISM_FluxMeshCreate(comp, nfltype, fluxMesh, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    field = ESMF_FieldCreate(fluxMesh, name="river_flux_segment", &
      meshloc=ESMF_MESHLOC_ELEMENT, typekind=ESMF_TYPEKIND_R8, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    ! zero-init; guard with root-PET test since only PET 0 owns elements
    ! (an empty-DE pointer dereference elsewhere would be unsafe).
    call ESMF_GridCompGet(comp, localPet=fluxLocalPet, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    if (fluxLocalPet == 0) then
      call ESMF_FieldGet(field, farrayPtr=fluxFarrayPtr, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
      fluxFarrayPtr(:) = 0.0_ESMF_KIND_R8
    end if

    call NUOPC_Realize(importState, field=field, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    write(message,'(A,I0,A)') trim(compName)//' realized field river_flux_segment on ', &
      nfltype, ' element(s) (PET 0 only)'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  end if
#endif

  !> The list of export states is declared in InitializeAdvertise
  call ESMF_StateGet(exportState, itemCount=itemCount, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(itemNameList(itemCount), itemTypeList(itemCount), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateGet(exportState, itemTypeList=itemTypeList,  &
    itemNameList=itemNameList, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  do i=1, itemCount

    if (itemTypeList(i) /= ESMF_STATEITEM_FIELD) cycle

    call SCHISM_FieldRealize(exportState, itemNameList(i), &
      mesh=mesh2d, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    write(message,'(A)') trim(compName)//' realized field '//trim(itemNameList(i))// &
      ' in its export state'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  enddo

  if (allocated(itemNameList)) deallocate(itemNameList)
  if (allocated(itemTypeList)) deallocate(itemTypeList)

  call ESMF_StateGet(importState, itemCount=itemCount, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (itemCount > 0) then 
    allocate(itemNameList(itemCount), itemTypeList(itemCount), stat=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

    call ESMF_StateGet(importState, itemTypeList=itemTypeList,  &
      itemNameList=itemNameList, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  endif

  do i=1, itemCount

    if (itemTypeList(i) /= ESMF_STATEITEM_FIELD) cycle

    write(message,'(A)') trim(compName)//' realized field '//trim(itemNameList(i))// &
      ' in its import state'
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  enddo

  
  !> As IPDv01p3 (NUOPC Provided) fails when fields are not connected, we here
  !> remove all unconnected fields from the import and export stabilityTimeStep
  call SCHISM_RemoveUnconnectedFields(importState, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call SCHISM_RemoveUnconnectedFields(exportState, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  
  !@TODO: should we destroy field & array vars?
end subroutine

#ifdef USE_NUOPC_RIVER
#undef ESMF_METHOD
#define ESMF_METHOD "SCHISM_FluxMeshCreate"
!> @description Builds a tiny ESMF_Mesh of n dummy triangles to realize
!> 'river_flux_segment'; only element count/ids matter, geometry is arbitrary.
subroutine SCHISM_FluxMeshCreate(comp, n, mesh, rc)

  implicit none

  type(ESMF_GridComp), intent(in)    :: comp
  integer(ESMF_KIND_I4), intent(in)  :: n
  type(ESMF_Mesh), intent(out)       :: mesh
  integer(ESMF_KIND_I4), intent(out) :: rc

  integer(ESMF_KIND_I4) :: localrc, fluxLocalPet, i
  integer(ESMF_KIND_I4) :: numFluxNodes, numFluxElements

  type(ESMF_DistGrid) :: fluxNodeDistgrid, fluxElementDistgrid

  integer, allocatable            :: fluxNodeIds(:), fluxNodeOwners(:)
  integer, allocatable            :: fluxElementIds(:), fluxElementTypes(:)
  integer, allocatable            :: fluxElementConn(:)
  real(ESMF_KIND_R8), allocatable :: fluxNodeCoords(:)

  rc = ESMF_SUCCESS
  localrc = ESMF_SUCCESS

  call ESMF_GridCompGet(comp, localPet=fluxLocalPet, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! PET 0 owns all n elements/3n nodes; every other PET owns none.
  if (fluxLocalPet == 0) then
    numFluxNodes = 3*n
    numFluxElements = n
  else
    numFluxNodes = 0
    numFluxElements = 0
  end if

  allocate(fluxNodeIds(numFluxNodes), fluxNodeOwners(numFluxNodes), &
    fluxNodeCoords(2*numFluxNodes), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  allocate(fluxElementIds(numFluxElements), fluxElementTypes(numFluxElements), &
    fluxElementConn(3*numFluxElements), stat=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! PET 0 only (loop no-op elsewhere): fill n disjoint unit triangles, spaced
  ! along x so every node has a distinct coordinate. Element global id = j = i
  ! (flux.th column order).
  do i = 1, numFluxElements
    fluxNodeIds(3*i-2) = 3*i-2
    fluxNodeIds(3*i-1) = 3*i-1
    fluxNodeIds(3*i)   = 3*i
    fluxNodeOwners(3*i-2:3*i) = 0

    fluxNodeCoords(2*(3*i-2)-1) = 2.0_ESMF_KIND_R8*real(i-1,ESMF_KIND_R8)
    fluxNodeCoords(2*(3*i-2))   = 0.0_ESMF_KIND_R8
    fluxNodeCoords(2*(3*i-1)-1) = 2.0_ESMF_KIND_R8*real(i-1,ESMF_KIND_R8) + 1.0_ESMF_KIND_R8
    fluxNodeCoords(2*(3*i-1))   = 0.0_ESMF_KIND_R8
    fluxNodeCoords(2*(3*i)-1)   = 2.0_ESMF_KIND_R8*real(i-1,ESMF_KIND_R8) + 0.5_ESMF_KIND_R8
    fluxNodeCoords(2*(3*i))     = 1.0_ESMF_KIND_R8

    fluxElementIds(i) = i
    fluxElementTypes(i) = ESMF_MESHELEMTYPE_TRI
    fluxElementConn(3*i-2) = 3*i-2
    fluxElementConn(3*i-1) = 3*i-1
    fluxElementConn(3*i)   = 3*i
  end do

  ! Collective calls: every PET (including those with zero-size local arrays
  ! above) participates with its own local piece.
  fluxNodeDistgrid = ESMF_DistGridCreate(fluxNodeIds, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  fluxElementDistgrid = ESMF_DistGridCreate(fluxElementIds, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  mesh = ESMF_MeshCreate(parametricDim=2, spatialDim=2, coordSys=ESMF_COORDSYS_CART, &
    nodeIds=fluxNodeIds, nodeCoords=fluxNodeCoords, nodeOwners=fluxNodeOwners, &
    nodalDistgrid=fluxNodeDistgrid, &
    elementIds=fluxElementIds, elementTypes=fluxElementTypes, elementConn=fluxElementConn, &
    elementDistgrid=fluxElementDistgrid, &
    rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  deallocate(fluxNodeIds, fluxNodeOwners, fluxNodeCoords)
  deallocate(fluxElementIds, fluxElementTypes, fluxElementConn)

end subroutine SCHISM_FluxMeshCreate
#endif

#undef ESMF_METHOD
#define ESMF_METHOD "SetClock"
!> @description A model's clock copies start/stop time and timestep from its
!> parent's clock.  If a model has a timestep that is different (smaller) than
!> the parent's it needs to be set here.
subroutine SetClock(comp, rc)

  use schism_glbl, only : start_year, start_month, start_day, start_hour, rnday
  use schism_glbl, only : wtiminc, wtime2

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  type(ESMF_Clock)           :: driverClock, modelClock
  type(ESMF_TimeInterval)    :: runDur, timeStep
  type(ESMF_Time)            :: startTime, stopTime
  integer                    :: localrc, d, h, m
  real(ESMF_KIND_R8)         :: seconds
  character(len=ESMF_MAXSTR) :: message

  rc = ESMF_SUCCESS

  call NUOPC_ModelGet(comp, driverClock=driverClock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Set start time
  call ESMF_TimeSet(startTime, yy=start_year, mm=start_month, dd=start_day, h=int(start_hour), rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Set stop time
  d = int(rnday)
  h = int((rnday-d)*24)
  m = int((rnday-d)*24*60-h*60)
  call ESMF_TimeIntervalSet(runDur, d=d, h=h, m=m, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  
  stopTime = startTime+runDur

  !> Time step must be same with the driver
  call ESMF_ClockGet(driverClock, timeStep=timeStep, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Create component clock
  modelClock = ESMF_ClockCreate(timeStep, startTime, stopTime=stopTime, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Update component clock
  call ESMF_GridCompSet(comp, clock=modelClock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Check that wtiminc, i.e. the time between two new atmospheric inputs
  !> corresponds to the parent (coupling) time step
  call ESMF_TimeIntervalGet(timeStep, s_r8=seconds, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  if (abs(wtiminc - seconds) > 1e-6) then 
    write(message, '(A,I7,A,I7)') 'Check setting of wtiminc = ', int(wtiminc), &
      ' in param.nml! Auto-resetting to wtiminc = ', int(seconds)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
  endif

  wtiminc = seconds
  wtime2 = wtiminc !Also need to reset wtime2 to make it consistent

end subroutine

#undef ESMF_METHOD
#define ESMF_METHOD "SetRunClock"
!> @description If the timestep of the parent is dynamic, then there might
!> be mismatches between intergral timesteps of a model and the parent
!> timestep.  This is reported as an error unless dealt with in this label
!> @todo (not registered and fully implemented yet, so not used)
subroutine SetRunClock(comp, rc)

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  type(ESMF_Clock) :: driverClock, modelClock
  type(ESMF_Time)  :: driverCurrTime
  integer          :: localrc

  rc = ESMF_SUCCESS

  !> query driver and the component clocks
  call NUOPC_ModelGet(comp, driverClock=driverClock, modelClock=modelClock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_ClockGet(driverClock, currTime=driverCurrTime, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> set model clock to have the current start time as the driver clock
  call ESMF_ClockSet(modelClock, currTime=driverCurrTime, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> check the component clock against the driver clock
  call NUOPC_CompCheckSetClock(comp, driverClock, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

end subroutine

#undef ESMF_METHOD
#define ESMF_METHOD "DataInitialize"
subroutine DataInitialize(comp, rc)

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  ! local variables
  type(ESMF_Time) :: currTime
  type(ESMF_Clock) :: clock
  type(ESMF_State) :: exportState
  integer(ESMF_KIND_I4) :: localrc
  character(len=*), parameter :: subname = '(DataInitialize): '
  !--------------------------------

  rc = ESMF_SUCCESS
  call ESMF_LogWrite(trim(subname)//' called', ESMF_LOGMSG_INFO)

  !> Query component for its clock, import and export states
  call NUOPC_ModelGet(comp, modelClock=clock, exportState=exportState, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Update fields on export state
  call SCHISM_Export(comp, exportState, clock, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_LogWrite(trim(subname)//' done', ESMF_LOGMSG_INFO)

end subroutine DataInitialize

#undef ESMF_METHOD
#define ESMF_METHOD "ModelAdvance"
!> @description Requests SCHISM to take one timestep forward; NUOPC advances
!> the clock automatically.
subroutine ModelAdvance(comp, rc)

  use schism_glbl, only: dt

  implicit none

  !> Input/output variables
  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  !> Local variables 
  type(ESMF_Clock)        :: clock
  type(ESMF_State)        :: importState, exportState
  type(ESMF_Time)         :: currTime
  type(ESMF_TimeInterval) :: timeStep
  character(len=160)      :: message
  integer(ESMF_KIND_I4)   :: localrc
  integer, save           :: it = 1
  integer                 :: i, num_schism_steps
  real(ESMF_KIND_R8)      :: seconds
  character(len=*), parameter :: subname = '(ModelAdvance): '
  !--------------------------------

  rc = ESMF_SUCCESS
  call ESMF_LogWrite(trim(subname)//' called', ESMF_LOGMSG_INFO)

  !> Query component for its clock, import and export states
  call NUOPC_ModelGet(comp, modelClock=clock, importState=importState, &
    exportState=exportState, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Update fields on import state
  call SCHISM_Import(comp, importState, clock, rc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

#ifdef USE_NUOPC_RIVER
  !> Fill ath3 before the timestep loop below, so schism_step interpolates
  !> the new level (see SCHISM_ImportRiver).
  call SCHISM_ImportRiver(comp, importState, rc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Fill qthcon before the timestep loop; no-op unless nuopc_flux_active
  !> (see SCHISM_ImportFlux).
  call SCHISM_ImportFlux(comp, importState, rc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
#endif

  !> Write log about advance
  call ESMF_ClockPrint(clock, options="currTime", &
      preString="--- advancing schism from ", unit=message, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_TimePrint(currTime + timeStep, &
      preString=trim(message)//" to ", unit=message, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_LogWrite(message, ESMF_LOGMSG_INFO, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Run SCHISM
  call ESMF_TimeIntervalGet(timeStep, s_r8=seconds, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  num_schism_steps=int(seconds/dt) 
  if (mod(seconds,dt) /= 0) then
          write(message, '(A)') 'Coupling step cannot be divided by dt, please adjust to avoid lack of steps! '
     call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
  end if

  do i = it, it+num_schism_steps-1
     call schism_step(i)
     it = it + 1
  end do

  !> Update fields on export state
  call SCHISM_Export(comp, exportState, clock, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_LogWrite(trim(subname)//' done', ESMF_LOGMSG_INFO)

end subroutine ModelAdvance

#undef ESMF_METHOD
#define ESMF_METHOD "SCHISM_RemoveUnconnectedFields"
subroutine SCHISM_RemoveUnconnectedFields(state, rc)

  implicit none

  type(ESMF_State), intent(inout)      :: state
  integer(kind=ESMF_KIND_I4), optional :: rc

  integer(kind=ESMF_KIND_I4)              :: rc_, localrc, itemCount, i
  type(ESMF_StateItem_Flag), allocatable  :: itemTypeList(:)
  character(len=ESMF_MAXSTR), allocatable :: itemNameList(:)
  character(len=ESMF_MAXSTR)              :: message
  type(ESMF_Field)                        :: field
  logical                                 :: isPresent

  if (present(rc)) rc = ESMF_SUCCESS

  call ESMF_StateGet(state, itemCount=itemCount, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  allocate(itemTypeList(itemCount))
  allocate(itemNameList(itemCount))

  call ESMF_StateGet(state, itemTypeList=itemTypeList,  &
    itemNameList=itemNameList, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

  do i=1, itemCount

    if (itemTypeList(i) /= ESMF_STATEITEM_FIELD) cycle

#ifdef USE_NUOPC_RIVER
    ! Exempt the river import from removal (nwm_coupling=true): this IPDv00
    ! pass runs before the direct NWM->OCN connector has connected the field.
    if (keep_river_field .and. trim(itemNameList(i))=='river_volume_flux') then
      write(message,'(A)') '--- keeping river_volume_flux (nwm_coupling=true); connection completes later'
      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
      cycle
    end if

    ! Same exemption for the flux-boundary N-point field: nwm_flux_coupling=true
    ! promises a dnwm provider that has not connected yet.
    if (keep_flux_field .and. trim(itemNameList(i))=='river_flux_segment') then
      write(message,'(A)') '--- keeping river_flux_segment (nwm_flux_coupling=true); connection completes later'
      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
      cycle
    end if
#endif

    call ESMF_StateGet(state, trim(itemNameList(i)), field=field, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    isPresent=.true.
!    call ESMF_AttributeGet(field, name='Connected', isPresent=isPresent, rc=localrc)
!    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    write(message,'(A)') '--- checking connection state of '//trim(itemNameList(i))
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
    !if (isPresent) isPresent = NUOPC_IsConnected(state,trim(itemNameList(i)), rc=localrc)
    if (isPresent) isPresent = NUOPC_IsConnected(field, rc=localrc)
    if(localrc/=ESMF_SUCCESS) then
      isPresent=.false.
      localrc=ESMF_SUCCESS
    endif
    
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

    if (.not.isPresent) then
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

      call ESMF_StateRemove(state, itemNameList(i:i), rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc_)

      write(message,'(A)') '--- removed unconnected field '//trim(itemNameList(i))
      call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
    endif

  enddo
end subroutine SCHISM_RemoveUnconnectedFields

#undef ESMF_METHOD
#define ESMF_METHOD "SCHISM_Export"
subroutine SCHISM_Export(comp, exportState, clock, rc)

  use schism_glbl,      only: nvrt, eta2, dav, uu2, vv2, tr_nd, idry_e, npa, deta1_dxy_elem, dp, znl, nvrt
  use schism_esmf_util, only: SCHISM_StateUpdate

  implicit none

  !> Input/output variables
  type(ESMF_GridComp), intent(in)    :: comp
  type(ESMF_State)   , intent(inout) :: exportState
  type(ESMF_Clock)   , intent(in)    :: clock
  integer            , intent(inout) :: rc

  !> Local variables
  type(ESMF_Time) :: currTime
  type(type_InternalStateWrapper) :: internalState
  type(type_InternalState), pointer :: isDataPtr => null()
  integer(ESMF_KIND_I4) :: localrc
  real(ESMF_KIND_R8), allocatable, save, target :: idry_r8(:)
  real(ESMF_KIND_R8), allocatable, save, target :: sst_K(:)
  real(ESMF_KIND_R8), allocatable, save, target :: hmix(:)
  character(len=ESMF_MAXSTR) :: timeStr
  character(len=*), parameter :: subname = '(SCHISM_Export): '
  !--------------------------------

  rc = ESMF_SUCCESS
  call ESMF_LogWrite(trim(subname)//' called', ESMF_LOGMSG_INFO)

  !> Query internal state
  call ESMF_GridCompGetInternalState(comp, internalState, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  isDataPtr => internalState%wrap

  if(.not.associated(isDataPtr)) localrc = ESMF_RC_PTR_NULL
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Update fields on export state
  !> sea surface height
  call SCHISM_StateUpdate(exportState, 'sea_surface_height_above_sea_level', eta2, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Sea surface graidend (for cice coupling) 
  call SCHISM_StateUpdate(exportState, 'sea_surface_slope_zonal', deta1_dxy_elem(:,1), &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call SCHISM_StateUpdate(exportState, 'sea_surface_slope_merid', deta1_dxy_elem(:,2), &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  
  !> depth average current in x direction
  call SCHISM_StateUpdate(exportState, 'depth-averaged_x-velocity', dav(1,:), &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> depth average current in y direction
  call SCHISM_StateUpdate(exportState, 'depth-averaged_y-velocity', dav(2,:), &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> surface current in x direction
  call SCHISM_StateUpdate(exportState, 'ocn_current_zonal', uu2(nvrt,:), &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> surface current in y direction
  call SCHISM_StateUpdate(exportState, 'ocn_current_merid', vv2(nvrt,:), &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)


  !> mixedlayer depth (CICE coupling) 
  if (.not. allocated(hmix)) then
     allocate(hmix(npa))
     hmix(:) = min( abs(znl(nvrt,:) - znl(nvrt-1,:)) , max( 0.d0, dp(:) + eta2(:) ))
  end if

  call SCHISM_StateUpdate(exportState, 'mixed_layer_depth', hmix(:), &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> surface temperature
  if (.not. allocated(sst_K)) then
     allocate(sst_K(npa))
     sst_K(:) =0.0d0+273.15d0
  end if
  sst_K(:) = tr_nd(1,nvrt,:)+273.15d0 !Change unit to K
  call SCHISM_StateUpdate(exportState, 'sea_surface_temperature', sst_K(:), &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> surface salinity
  call SCHISM_StateUpdate(exportState, 'sea_surface_salinity', tr_nd(2,nvrt,:), &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> ocean mask
  !> mediator expects ocean mask in double rather then integer
  if (.not. allocated(idry_r8)) then
     allocate(idry_r8(size(idry_e)))
     idry_r8(:) =0.0d0
  end if
  idry_r8(:) = dble(idry_e(:))

  call SCHISM_StateUpdate(exportState, 'ocean_mask', idry_r8, &
    isPtr=isDataPtr, onElement=.true., rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Write fields on export state for debugging
  if (debug_level > 0) then
     call ESMF_ClockGet(clock, currTime=currTime, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

     call ESMF_TimeGet(currTime, timeStringISOFrac=timeStr , rc=localrc)
     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

     call SCHISM_StateWriteVTK(exportState, 'export_'//trim(timeStr), rc)
     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  end if

  call ESMF_LogWrite(trim(subname)//' done', ESMF_LOGMSG_INFO)

end subroutine SCHISM_Export

#undef ESMF_METHOD
#define ESMF_METHOD "SCHISM_Import"
subroutine SCHISM_Import(comp, importState, clock, rc)

  use schism_glbl     , only: RADFLAG, windx2, windy2, pr2
  use schism_glbl     , only: airt2,shum2,srad,hradd,fluxprc,npa
  use schism_glbl     , only: uvice,vvice,taux,tauy,vsno,vice, &
                              aice,ifresh_flux,isalt_flux,iheat_flux, &
                              isw_pen,frzmlt,tau_oi,fresh_wa_flux, &
                              salinity_flux,net_heat_flux,srad_th_ice,CdnIO,znl, nvrt
  use schism_esmf_util, only: SCHISM_StateImportWaveTensor
  use schism_esmf_util, only: SCHISM_StateImportWave3dVortex
  use schism_esmf_util, only: SCHISM_StateUpdate

  implicit none

  !> Input/output variables
  type(ESMF_GridComp), intent(in)    :: comp
  type(ESMF_State)   , intent(inout) :: importState
  type(ESMF_Clock)   , intent(in)    :: clock
  integer            , intent(inout) :: rc
  
  !> Local variables
  type(ESMF_Time) :: currTime
  integer         :: i
  type(type_InternalStateWrapper) :: internalState
  type(type_InternalState), pointer :: isDataPtr => null()
  integer(ESMF_KIND_I4) :: localrc
  character(len=ESMF_MAXSTR) :: timeStr
  character(len=*), parameter :: subname = '(SCHISM_Import): '
  !--------------------------------

  rc = ESMF_SUCCESS
  call ESMF_LogWrite(trim(subname)//' called', ESMF_LOGMSG_INFO)

  !> Query internal state
  call ESMF_GridCompGetInternalState(comp, internalState, localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  isDataPtr => internalState%wrap

  if(.not.associated(isDataPtr)) localrc = ESMF_RC_PTR_NULL
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Update fields on import state
  if (RADFLAG == 'VOR') then
     !> Obtain required variables from wave component to do coupling with vortex formulation
     call SCHISM_StateImportWave3dVortex(importState, isDataPtr, rc=localrc)
     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  else
     !> Obtain radiation tensor from wave component and calculate the wave stress
     call SCHISM_StateImportWaveTensor(importState, isDataPtr, rc=localrc)
     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  end if

  !> surface wind component in x direction
  call SCHISM_StateUpdate(importState, 'inst_zonal_wind_height10m', windx2, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Surface wind component in y direction
  call SCHISM_StateUpdate(importState, 'inst_merid_wind_height10m', windy2, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Surface air pressure
  call SCHISM_StateUpdate(importState, 'air_pressure_at_sea_level', pr2, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Air temperature at sea level
  call SCHISM_StateUpdate(importState, 'inst_temp_height2m', airt2, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Air humidity at sea level
  call SCHISM_StateUpdate(importState, 'inst_spec_humid_height2m', shum2, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Downward shortwave at sea level
  call SCHISM_StateUpdate(importState, 'inst_net_sw_flx', srad, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Downward longwave at sea level
  call SCHISM_StateUpdate(importState, 'inst_down_lw_flx', hradd, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Precipatation rate
  call SCHISM_StateUpdate(importState, 'inst_prec_rate', fluxprc, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    
  !>-------------------------------------------------------
  !>              Allocated vars to dump import to
  !>-------------------------------------------------------
  
  if (.not. allocated(uvice)) then
      allocate(uvice(npa))
      uvice(:) =0.0d0
  end if
  if (.not. allocated(vvice)) then
      allocate(vvice(npa))
      vvice(:) =0.0d0
  end if
  if (.not. allocated(taux)) then
      allocate(taux(npa))
      taux(:) =0.0d0
  end if
  if (.not. allocated(tauy)) then
      allocate(tauy(npa))
      tauy(:) =0.0d0
  end if
  if (.not. allocated(vsno)) then
      allocate(vsno(npa))
      vsno(:) =0.0d0
  end if
  if (.not. allocated(vice)) then
      allocate(vice(npa))
      vice(:) =0.0d0
  end if
  if (.not. allocated(aice)) then
      allocate(aice(npa))
      aice(:) =0.0d0
  end if
  if (.not. allocated(ifresh_flux)) then
      allocate(ifresh_flux(npa))
      ifresh_flux(:) =0.0d0
  end if
  if (.not. allocated(isalt_flux)) then
      allocate(isalt_flux(npa))
      isalt_flux(:) =0.0d0
  end if
  if (.not. allocated(iheat_flux)) then
      allocate(iheat_flux(npa))
      iheat_flux(:) =0.0d0
  end if
  if (.not. allocated(isw_pen)) then
      allocate(isw_pen(npa))
      isw_pen(:) =0.0d0
  end if
  if (.not. allocated(frzmlt)) then
      allocate(frzmlt(npa))
      frzmlt(:) =0.0d0
  end if
  if (.not. allocated(CdnIO)) then
      allocate(CdnIO(npa))
      CdnIO(:) =0.0d0
  end if
  
  !  Importing CICE vars into schism
  !> Zonal-direction ------------------------------------
  call SCHISM_StateUpdate(importState, 'Si_uvel', uvice, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  !> Merid-direction ------------------------------------
  call SCHISM_StateUpdate(importState, 'Si_vvel', vvice, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  !> ice-stress Zonal-direction -------------------------
  call SCHISM_StateUpdate(importState, 'Fioi_taux', taux, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  !> ice-stress Merid-direction -------------------------
  call SCHISM_StateUpdate(importState, 'Fioi_tauy', tauy, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  !> Volume of snow ------------------------------------
  call SCHISM_StateUpdate(importState, 'Si_vsno', vsno, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Volume of ice -------------------------------------
  call SCHISM_StateUpdate(importState, 'Si_vice', vice, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> ice_fraction ---------------------------------------
  call SCHISM_StateUpdate(importState, 'Si_ifrac', aice, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Melt_water_flux --------------------------------------------
  call SCHISM_StateUpdate(importState, 'Fioi_meltw', ifresh_flux, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Melt_water_flux -------------------------------------------
  call SCHISM_StateUpdate(importState, 'Fioi_salt', isalt_flux, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  !> Heat_flux_ice_to_ocn ---------------------------------------
  call SCHISM_StateUpdate(importState, 'Fioi_melth', iheat_flux, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
 
  !> Heat_flux_ice_to_ocn ------------------------------------
  call SCHISM_StateUpdate(importState, 'Fioi_swpen', isw_pen, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
 
  !> freezing melting potential  ---------------------------
  call SCHISM_StateUpdate(importState, 'Si_frzmlt', frzmlt, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Ice-to-ocean drag coeff  ------------------------------
  call SCHISM_StateUpdate(importState, 'Si_CdnIO', CdnIO, &
    isPtr=isDataPtr, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)


  ! CICE Export to hydro
  i=0
  do i = 1,npa
    if (aice(i) > real(0.0)) then
      !>Ice ocean stress -------------------------------
      !>Taux,Tauy are in units of [N/m/m]

      tau_oi(1,i)=taux(i)
      tau_oi(2,i)=tauy(i)


      !> Salinity flux ---------------------------------
      !> This is the slainity flux to the ocean from ice 
      !> formaiton and melt. isalt has units [kg/s/m/m] 
      
      salinity_flux(i) = aice(i)*(isalt_flux(i))/real(1000)

      !>Fresh water flux -------------------------------
      !> ifresh_flux is in units of [kg/s/m/m]
      !> Water is distributed across whole element

      fresh_wa_flux(i) = ifresh_flux(i)

      !>Heat flux ice to ocean  ------------------------
      !> iheat_flux is in units of [W/m/m]
      !> Energy is distributed across whole element

      net_heat_flux(i) = iheat_flux(i)

      !>Short-wave pen. flux ---------------------------
      !>isw_pen is in units of [W/m/m]
      !>Weighted by ice area
      srad_th_ice(i)  = isw_pen(i)

    else

      !> No ice so all these values are zero
      tau_oi(1,i)      = real(0)
      tau_oi(2,i)      = real(0)
      salinity_flux(i) = real(0)
      fresh_wa_flux(i) = real(0) !+ max(real(0.0),frzmlt(i))
      net_heat_flux(i) = real(0)
      srad_th_ice(i)  = real(0)
    
    endif
  enddo

  !> Write fields on import state for debugging
  if (debug_level > 0) then
     call ESMF_ClockGet(clock, currTime=currTime, rc=localrc)
     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

     call ESMF_TimeGet(currTime, timeStringISOFrac=timeStr , rc=localrc)
     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

     call SCHISM_StateWriteVTK(importState, 'import_'//trim(timeStr), rc)
     _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  end if

  call ESMF_LogWrite(trim(subname)//' done', ESMF_LOGMSG_INFO)

end subroutine SCHISM_Import

#ifdef USE_NUOPC_RIVER
#undef ESMF_METHOD
#define ESMF_METHOD "SCHISM_ImportRiver"
!> @description Fills SCHISM's rank-replicated ath3 (zero-order hold) before
!> the timestep, from 'river_volume_flux' or the river_stub_q fallback.
subroutine SCHISM_ImportRiver(comp, importState, rc)

  use NUOPC,       only: NUOPC_IsConnected
  use schism_glbl, only: nsources, nsinks, ath3, if_source, ieg_source, iegl, ne
  use schism_msgp, only: myrank

  type(ESMF_GridComp)  :: comp
  type(ESMF_State)     :: importState
  integer, intent(out) :: rc

  integer(ESMF_KIND_I4)           :: localrc
  type(ESMF_Field)                :: field
  type(ESMF_VM)                   :: vm
  real(ESMF_KIND_R8), pointer     :: farrayPtr1(:) => null()
  real(ESMF_KIND_R8), allocatable :: q_local(:), q_global(:)
  integer                         :: i, ie
  logical                         :: connected
  type(ESMF_StateItem_Flag)       :: itemtype_river
  character(len=ESMF_MAXSTR)      :: message

  rc = ESMF_SUCCESS

  ! USE_NUOPC_RIVER requires if_source=1 (sparse source_sink.in read in
  ! schism_init). Nothing to do if sources are off or none were paired.
  if (if_source == 0 .or. nsources <= 0) return

  ! Sinks are not supported by this one-way prototype (only vsource is filled);
  ! fail loudly rather than silently ignore a configured sink.
  if (nsinks > 0) then
    write(message,'(A,I0,A)') 'SCHISM_ImportRiver: nsinks=', nsinks, &
      ' but sinks (vsink) are NOT supported by USE_NUOPC_RIVER; '// &
      'use a source-only source_sink.in (nsinks=0).'
    call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=trim(message), ESMF_CONTEXT, rcToReturn=rc)
    return
  end if

  ! Check presence before connection (stub runs have it removed already) to
  ! avoid an ESMF_StateGet abort on a missing item; absence/disconnection is
  ! handled explicitly below.
  connected = .false.
  call ESMF_StateGet(importState, itemName="river_volume_flux", itemType=itemtype_river, rc=localrc)
  if (localrc == ESMF_SUCCESS .and. itemtype_river == ESMF_STATEITEM_FIELD) then
    call ESMF_StateGet(importState, itemName="river_volume_flux", field=field, rc=localrc)
    if (localrc == ESMF_SUCCESS) then
      connected = NUOPC_IsConnected(field, rc=localrc)
      if (localrc /= ESMF_SUCCESS) connected = .false.
    end if
  end if

  if (connected) then
    ! Connected: field is realized on the SCHISM element mesh (owned elements
    ! only), but ath3 must be rank-replicated. Each rank reads only the
    ! sources it owns, then VMAllReduce(SUM) replicates the full list.
    call ESMF_FieldGet(field, farrayPtr=farrayPtr1, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    call ESMF_GridCompGet(comp, vm=vm, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    allocate(q_local(nsources), q_global(nsources))
    q_local = 0.0_ESMF_KIND_R8
    do i = 1, nsources
      if (iegl(ieg_source(i))%rank == myrank .and. iegl(ieg_source(i))%id <= ne) then
        ie = iegl(ieg_source(i))%id
        q_local(i) = max(0.0_ESMF_KIND_R8, farrayPtr1(ie))   ! owner reads its element [m^3/s]
      end if
    end do
    call ESMF_VMAllReduce(vm, q_local, q_global, nsources, ESMF_REDUCE_SUM, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    ath3(1:nsources,1,2,1) = real(q_global, 4)        ! new vsource, replicated on all ranks
    ath3(1:nsources,1,1,1) = ath3(1:nsources,1,2,1)   ! zero-order hold (old=new)
    write(message,'(A,I0,A,F0.3)') 'SCHISM_ImportRiver: CONNECTED to NWM provider; nsources=', &
      nsources, ', vsource(1) [m^3/s] = ', ath3(1,1,2,1)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
    deallocate(q_local, q_global)
  else if (keep_river_field) then
    ! nwm_coupling=true but the field is absent/unconnected: provider not
    ! wired. Abort rather than silently substitute a fabricated discharge.
    write(message,'(A)') 'SCHISM_ImportRiver: nwm_coupling=true but river_volume_flux '// &
      'is absent/unconnected -- NWM provider not wired. Refusing to fabricate flow. '// &
      'Wire the NWM component+connector, or set river_stub=true for the constant stub.'
    call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=trim(message), ESMF_CONTEXT, rcToReturn=rc)
    return
  else if (river_stub_mode) then
    ! Stub mode: inject a uniform constant discharge into every source element
    ! (old level set equal to new, i.e. a zero-order hold).
    ath3(1:nsources,1,2,1) = real(max(0.0_ESMF_KIND_R8, river_stub_q), 4)
    ath3(1:nsources,1,1,1) = ath3(1:nsources,1,2,1)
    write(message,'(A,F0.3)') 'SCHISM_ImportRiver: STUB MODE (river_stub=true, no NWM provider); '// &
      'vsource [m^3/s] = ', real(max(0.0_ESMF_KIND_R8, river_stub_q), 4)
    call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)
  else
    ! USE_NUOPC_RIVER is active with sources, but neither a provider (nwm_coupling)
    ! nor the explicit stub (river_stub) was selected. Ambiguous -- abort rather
    ! than guess. (Previously this branch silently fabricated river_stub_q.)
    write(message,'(A)') 'SCHISM_ImportRiver: river_volume_flux absent/unconnected and '// &
      'no mode selected. Set nwm_coupling=true (connect an NWM provider) or '// &
      'river_stub=true (constant river_stub_q stub).'
    call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=trim(message), ESMF_CONTEXT, rcToReturn=rc)
    return
  end if

end subroutine SCHISM_ImportRiver

#undef ESMF_METHOD
#define ESMF_METHOD "SCHISM_ImportFlux"
!> @description Drives qthcon for every ifltype=1 segment from
!> 'river_flux_segment', replacing flux.th (see module-level flux doc).
subroutine SCHISM_ImportFlux(comp, importState, rc)

  use NUOPC,       only: NUOPC_IsConnected
  use schism_glbl, only: nuopc_flux_active, nfltype, nope_global, ifltype, qthcon

  type(ESMF_GridComp)  :: comp
  type(ESMF_State)     :: importState
  integer, intent(out) :: rc

  integer(ESMF_KIND_I4)           :: localrc, localPet
  type(ESMF_Field)                :: field
  type(ESMF_VM)                   :: vm
  real(ESMF_KIND_R8), pointer     :: farrayPtr1(:) => null()
  real(ESMF_KIND_R8), allocatable :: q_flux(:)
  integer                         :: j, k
  logical                         :: connected
  type(ESMF_StateItem_Flag)       :: itemtype_flux
  character(len=ESMF_MAXSTR)      :: message

  rc = ESMF_SUCCESS

  ! Not activated at init: honor the legacy flux.th path untouched. Gated on
  ! the init-time attribute, not a live NUOPC_IsConnected() check.
  if (.not. nuopc_flux_active) return

  ! Defensive only: nfltype==0 with nuopc_flux_active is already a hard abort
  ! in InitializeAdvertise, so not reachable via normal startup.
  if (nfltype <= 0) return

  call ESMF_GridCompGet(comp, vm=vm, localPet=localPet, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  ! Is 'river_flux_segment' present in the import state AND connected to a
  ! provider? Mirrors SCHISM_ImportRiver's presence-then-connection query.
  connected = .false.
  call ESMF_StateGet(importState, itemName="river_flux_segment", itemType=itemtype_flux, rc=localrc)
  if (localrc == ESMF_SUCCESS .and. itemtype_flux == ESMF_STATEITEM_FIELD) then
    call ESMF_StateGet(importState, itemName="river_flux_segment", field=field, rc=localrc)
    if (localrc == ESMF_SUCCESS) then
      connected = NUOPC_IsConnected(field, rc=localrc)
      if (localrc /= ESMF_SUCCESS) connected = .false.
    end if
  end if

  allocate(q_flux(nfltype))

  if (connected) then
    ! PET 0 reads the field (realized entirely on PET 0), then broadcasts the
    ! N-vector to every rank.
    if (localPet == 0) then
      call ESMF_FieldGet(field, farrayPtr=farrayPtr1, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

      ! Fail loudly on a size mismatch rather than read out of bounds. Must
      ! NOT return early (PET-0-only check): other PETs are about to enter
      ! the collective ESMF_VMBroadcast below and would hang waiting on it.
      if (size(farrayPtr1,1) /= nfltype) then
        write(message,'(A,I0,A,I0)') 'SCHISM_ImportFlux: river_flux_segment import '// &
          'size mismatch -- got ', size(farrayPtr1,1), ', expected nfltype=', nfltype
        call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=trim(message), ESMF_CONTEXT, rcToReturn=localrc)
        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
      end if

      ! v1 scope is inflow only. dnwm clamps negatives upstream, but defensively
      ! clamp here too so a negative value can never flip this into an outflow
      ! boundary via qthcon's sign convention.
      q_flux(1:nfltype) = max(0.0_ESMF_KIND_R8, farrayPtr1(1:nfltype))
    end if

    call ESMF_VMBroadcast(vm, bcstData=q_flux, count=nfltype, rootPet=0, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  else if (keep_flux_field) then
    ! nwm_flux_coupling=true but the field is absent/unconnected: abort
    ! rather than silently substitute a fabricated discharge.
    write(message,'(A)') 'SCHISM_ImportFlux: nwm_flux_coupling=true but river_flux_segment '// &
      'is absent/unconnected -- NWM provider not wired. Refusing to fabricate flow. '// &
      'Wire the NWM component+connector, or set flux_stub=true for the constant stub.'
    call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=trim(message), ESMF_CONTEXT, rcToReturn=rc)
    deallocate(q_flux)
    return

  else if (flux_stub_mode) then
    ! Explicit stub mode (flux_stub=true, no NWM provider): uniform constant
    ! discharge on every ifltype=1 segment, computed identically (no communication
    ! needed) on every rank.
    q_flux(1:nfltype) = max(0.0_ESMF_KIND_R8, flux_stub_q)

  else
    ! Unreachable: nuopc_flux_active (checked above) is defined as
    ! keep_flux_field .or. flux_stub_mode. Kept defensive in case that invariant
    ! changes.
    deallocate(q_flux)
    return
  end if

  ! Negate into SCHISM's outward-normal qthcon convention; identical on every rank.
  ! Entry j -> the j-th global boundary k with ifltype(k)==1, in bctides.in scan
  ! order -- EXACTLY the flux.th column order.
  j = 0
  do k = 1, nope_global
    if (ifltype(k) == 1) then
      j = j + 1
      qthcon(k) = -q_flux(j)
    end if
  end do

  ! Per-window diagnostic, positive-reported (the analog of the existing
  ! injected-discharge log in SCHISM_ImportRiver).
  write(message,'(A,I0,A,F0.3)') 'SCHISM_ImportFlux: injected inflow over ', nfltype, &
    ' ifltype=1 segment(s), sum(q) [m^3/s] = ', sum(q_flux(1:nfltype))
  call ESMF_LogWrite(trim(message), ESMF_LOGMSG_INFO)

  deallocate(q_flux)

end subroutine SCHISM_ImportFlux

#undef ESMF_METHOD
#define ESMF_METHOD "CheckImportRiver"
!> @description Re-implements NUOPC's at-currTime staleness check, exempting
!> only river_volume_flux/river_flux_segment (non-mediator-brokered dnwm).
subroutine CheckImportRiver(comp, rc)

  use NUOPC,       only: NUOPC_IsAtTime, NUOPC_IsConnected
  use NUOPC_Model, only: NUOPC_ModelGet

  type(ESMF_GridComp)  :: comp
  integer, intent(out) :: rc

  integer(ESMF_KIND_I4)                   :: localrc
  type(ESMF_Clock)                        :: clock
  type(ESMF_Time)                         :: currTime
  type(ESMF_State)                        :: importState
  type(ESMF_Field)                        :: field
  type(ESMF_StateItem_Flag), allocatable  :: itemTypeList(:)
  character(len=ESMF_MAXSTR), allocatable :: itemNameList(:)
  integer(ESMF_KIND_I4)                   :: i, itemCount
  logical                                 :: connected, atTime
  character(len=ESMF_MAXSTR)              :: message

  rc = ESMF_SUCCESS

  call NUOPC_ModelGet(comp, modelClock=clock, importState=importState, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  call ESMF_ClockGet(clock, currTime=currTime, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  call ESMF_StateGet(importState, itemCount=itemCount, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
  if (itemCount <= 0) return

  allocate(itemTypeList(itemCount), itemNameList(itemCount))
  call ESMF_StateGet(importState, itemTypeList=itemTypeList, &
    itemNameList=itemNameList, rc=localrc)
  _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)

  do i = 1, itemCount
    if (itemTypeList(i) /= ESMF_STATEITEM_FIELD) cycle
    ! Exempted from the at-currTime abort, but still tracked for staleness via
    ! CheckRiverProviderAdvancing (separate state per field).
    if (trim(itemNameList(i)) == 'river_volume_flux') then
      call ESMF_StateGet(importState, trim(itemNameList(i)), field=field, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
      connected = NUOPC_IsConnected(field, rc=localrc)
      if (localrc /= ESMF_SUCCESS) then
        localrc = ESMF_SUCCESS   ! indeterminate connection -> stub path, nothing to track
      else if (connected) then
        call CheckRiverProviderAdvancing(field, trim(itemNameList(i)), currTime, &
          river_stamp_seen, river_prev_stamp, river_prev_time, river_stall_count, rc=localrc)
        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
      end if
      cycle
    end if
    if (trim(itemNameList(i)) == 'river_flux_segment') then
      call ESMF_StateGet(importState, trim(itemNameList(i)), field=field, rc=localrc)
      _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
      connected = NUOPC_IsConnected(field, rc=localrc)
      if (localrc /= ESMF_SUCCESS) then
        localrc = ESMF_SUCCESS   ! indeterminate connection -> stub path, nothing to track
      else if (connected) then
        call CheckRiverProviderAdvancing(field, trim(itemNameList(i)), currTime, &
          flux_stamp_seen, flux_prev_stamp, flux_prev_time, flux_stall_count, rc=localrc)
        _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
      end if
      cycle
    end if
    call ESMF_StateGet(importState, trim(itemNameList(i)), field=field, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    connected = NUOPC_IsConnected(field, rc=localrc)
    if (localrc /= ESMF_SUCCESS) then
      localrc = ESMF_SUCCESS   ! treat an indeterminate connection state as unconnected
      cycle
    end if
    if (.not. connected) cycle
    atTime = NUOPC_IsAtTime(field, currTime, rc=localrc)
    _SCHISM_LOG_AND_FINALIZE_ON_ERROR_(rc)
    if (.not. atTime) then
      write(message,'(A)') 'CheckImportRiver: import field "'//trim(itemNameList(i))// &
        '" is NOT at the current time (stale coupling input)'
      call ESMF_LogSetError(ESMF_RC_VAL_WRONG, msg=trim(message), ESMF_CONTEXT, rcToReturn=rc)
      deallocate(itemTypeList, itemNameList)
      return
    end if
  end do

  deallocate(itemTypeList, itemNameList)

end subroutine CheckImportRiver

#undef ESMF_METHOD
#define ESMF_METHOD "CheckRiverProviderAdvancing"
!> @description Frozen-provider tripwire: warns (never aborts) if a field's
!> NUOPC TimeStamp stalls for RIVER_STALL_LIMIT windows while the clock advances.
subroutine CheckRiverProviderAdvancing(field, fieldName, currTime, &
  stampSeen, prevStamp, prevTime, stallCount, rc)

  type(ESMF_Field), intent(in)  :: field
  character(len=*), intent(in)  :: fieldName
  type(ESMF_Time),  intent(in)  :: currTime
  logical,          intent(inout) :: stampSeen
  integer(ESMF_KIND_I4), allocatable, intent(inout) :: prevStamp(:)
  type(ESMF_Time),  intent(inout) :: prevTime
  integer,          intent(inout) :: stallCount
  integer,          intent(out) :: rc

  integer(ESMF_KIND_I4)              :: localrc
  integer(ESMF_KIND_I4)              :: stampCount
  integer(ESMF_KIND_I4), allocatable :: stamp(:)
  logical                            :: timeAdvanced, stampAdvanced
  character(len=ESMF_MAXSTR)         :: message

  rc = ESMF_SUCCESS

  ! Number of integers in the NUOPC timestamp attribute (yy,mm,dd,h,m,s,...).
  call ESMF_AttributeGet(field, name="TimeStamp", convention="NUOPC", &
    purpose="Instance", itemCount=stampCount, rc=localrc)
  if (localrc /= ESMF_SUCCESS .or. stampCount <= 0) return   ! not stampable -> skip
  allocate(stamp(stampCount))
  call ESMF_AttributeGet(field, name="TimeStamp", convention="NUOPC", &
    purpose="Instance", valueList=stamp, rc=localrc)
  if (localrc /= ESMF_SUCCESS) then
    deallocate(stamp)
    return
  end if

  if (stampSeen) then
    timeAdvanced = (currTime > prevTime)
    if (allocated(prevStamp)) then
      if (size(prevStamp) == size(stamp)) then
        stampAdvanced = any(stamp /= prevStamp)
      else
        stampAdvanced = .true.   ! shape changed -- cannot compare, assume advanced
      end if
    else
      stampAdvanced = .true.
    end if
    if (timeAdvanced .and. .not. stampAdvanced) then
      stallCount = stallCount + 1
      if (stallCount >= RIVER_STALL_LIMIT) then
        write(message,'(A,I0,A)') &
          'CheckImportRiver: '//trim(fieldName)//' timestamp has not advanced for ', &
          stallCount, &
          ' coupling window(s) while the clock did -- SCHISM is holding a frozen '// &
          'import; verify the provider/stream.'
        call ESMF_LogWrite(trim(message), ESMF_LOGMSG_WARNING)
      end if
    else
      stallCount = 0
    end if
  end if

  ! record this window's observation for the next comparison
  if (allocated(prevStamp)) deallocate(prevStamp)
  allocate(prevStamp(size(stamp)))
  prevStamp = stamp
  prevTime  = currTime
  stampSeen = .true.

  deallocate(stamp)

end subroutine CheckRiverProviderAdvancing
#endif /*USE_NUOPC_RIVER*/

#undef ESMF_METHOD
#define ESMF_METHOD "SCHISM_StateWriteVTK"
subroutine SCHISM_StateWriteVTK(state, prefix, rc)

  implicit none

  !> Input/output variables
  type(ESMF_State), intent(in) :: state
  character(len=*), intent(in) :: prefix
  integer, intent(out), optional :: rc

  !> local variables
  integer :: i, itemCount
  type(ESMF_Field) :: field
  character(ESMF_MAXSTR), allocatable :: itemNameList(:)
  character(len=*),parameter :: subname='(SCHISM_StateWriteVTK)'
  !--------------------------------

  rc = ESMF_SUCCESS
  call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

  !> Get number of fields in the state
  call ESMF_StateGet(state, itemCount=itemCount, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  !> Get item names
  if (.not. allocated(itemNameList)) allocate(itemNameList(itemCount))

  call ESMF_StateGet(state, itemNameList=itemNameList, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  !> Loop over fields and write them
  do i = 1, itemCount
     !> Get field
     call ESMF_StateGet(state, itemName=trim(itemNameList(i)), field=field, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

     !> Write field
     call ESMF_FieldWriteVTK(field, trim(prefix)//'_'//trim(itemNameList(i)), rc=rc)
     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out
  end do

  !> Clean temporary variables
  if (allocated(itemNameList)) deallocate(itemNameList)

  call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

end subroutine SCHISM_StateWriteVTK

end module schism_nuopc_cap
