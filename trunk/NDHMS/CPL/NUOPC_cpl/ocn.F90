!==============================================================================
! Earth System Modeling Framework
! Copyright 2002-2019, University Corporation for Atmospheric Research, 
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics 
! Laboratory, University of Michigan, National Centers for Environmental 
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory, 
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

module OCN

  !-----------------------------------------------------------------------------
  ! OCN Component.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS      => SetServices, &
    model_label_SetClock  => label_SetClock, &
    model_label_Advance   => label_Advance
  
  use NWM_ESMF_Utility

  implicit none
  
  private
  
  public SetServices
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine SetServices(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(model, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeP2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! attach specializing method(s)
    call NUOPC_CompSpecialize(model, specLabel=model_label_SetClock, &
      specRoutine=SetClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP1(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! Beheen
    logical                    :: isPresent


    rc = ESMF_SUCCESS


    ! importable field: air_pressure_at_sea_level
    call NUOPC_Advertise(importState, &
      StandardName="air_pressure_at_sea_level", name="pmsl", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! importable field: surface_net_downward_shortwave_flux
    call NUOPC_Advertise(importState, &
      StandardName="surface_net_downward_shortwave_flux", name="rsns", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! importable field: 
    call NUOPC_Advertise(importState, &
      StandardName="flow_rate", name="streamflow", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable field: sea_surface_temperature
    call NUOPC_Advertise(exportState, &
      StandardName="sea_surface_temperature", name="sst", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Beheen waterlevel mockup
    ! exportable field: waterlevel
    isPresent = NUOPC_FieldDictionaryHasEntry( "water_level", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (.not.isPresent) then
        call NUOPC_FieldDictionaryAddEntry("water_level", "m", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
    endif

    call NUOPC_Advertise(exportState, &
      StandardName="water_level", name="wl", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP2(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables    
    type(ESMF_TimeInterval) :: stabilityTimeStep
    type(ESMF_VM)        :: vm
    integer              :: localPet
    integer              :: petCount
    integer,parameter    :: gblElementCnt = 185
    integer              :: gblElementExt
    integer              :: gblElementDiv
    integer              :: locElementBeg
    integer              :: locElementCnt
    integer,allocatable  :: arbSeqIndexList(:)
    integer              :: i,j
    type(ESMF_DistGrid)  :: distgrid
    type(ESMF_Field)     :: field
    type(ESMF_LocStream) :: locStreamIn
    type(ESMF_LocStream) :: locStreamOut
    
    ! Beheen
    type(ESMF_Field)    :: sstField
    type(ESMF_Field)    :: rsnsField
    type(ESMF_Field)    :: pmslField
    type(ESMF_Field)    :: streamflowField
    real, allocatable   :: rsnsarray(:)
    real, allocatable   :: pmslarray(:)
    real, allocatable   :: sstarray(:)
    real(ESMF_KIND_R8), pointer    :: streamflowarray(:)

    ! test for waterlevel
    type(ESMF_Grid)     :: grid
    type(ESMF_Mesh)     :: mesh
    type(ESMF_Field)    :: waterlevelField
    real(ESMF_KIND_R8),pointer :: dataWL(:)

    real(ESMF_KIND_R8), allocatable   :: waterlevelarray(:)
    real(ESMF_KIND_R8), dimension(:), pointer  :: wlPtr => null()
    character(160)      :: msgString
    integer             :: dimCount, numOwnedElements, numOwnedNodes
    integer             :: verbosity
    integer :: parametricDim
    integer :: spatialDim
    integer :: nodeCnt
    integer, allocatable :: nodeIds(:)
    real(ESMF_KIND_R8), allocatable :: nodeCoords(:)
    integer, allocatable :: nodeOwners(:)
    logical :: nodeMaskIsPresent
    integer, allocatable :: nodeMask(:)
    integer :: elementCnt
    integer, allocatable :: elementIds(:)
    integer, allocatable :: elementTypes(:)
    integer :: elementConnCount
    integer, allocatable :: elementConn(:)
    logical :: elementMaskIsPresent
    integer, allocatable :: elementMask(:)
    logical :: elementAreaIsPresent
    real(ESMF_KIND_R8), allocatable :: elementArea(:)
    logical :: elementCoordsIsPresent
    real(ESMF_KIND_R8), allocatable :: elementCoords(:)
    logical :: nodalDistgridIsPresent
    type(ESMF_DistGrid) :: nodalDistgrid
    logical :: elementDistgridIsPresent
    type(ESMF_DistGrid) :: elementDistgrid
    real(ESMF_KIND_R8), allocatable :: ownedNodeCoords(:)
    real(ESMF_KIND_R8), allocatable :: ownedElemCoords(:)
    logical :: isMemFreed
    type(ESMF_Array) :: elemMaskArray
    type(ESMF_Array) :: elemAreaArray
    type(ESMF_CoordSys_Flag):: coordSys
    type(ESMF_MeshStatus_Flag) :: status
    ! end test

    rc = ESMF_SUCCESS


    call ESMF_VMGetCurrent(vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(vm=vm, localPet=localPet, petCount=petCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out


    ! calculate local element count and first local element
    gblElementDiv = gblElementCnt/petCount
    gblElementExt = MOD(gblElementCnt,petCount)
    if (localPet.eq.(petCount-1)) then
      locElementCnt = gblElementDiv + gblElementExt
    else
      locElementCnt = gblElementDiv
    endif
    locElementBeg = 1 + (gblElementDiv*localPet)

    ! create local element list
    allocate(arbSeqIndexList(locElementCnt))
    allocate(streamflowarray(locElementCnt))
    allocate(rsnsarray(locElementCnt))
    allocate(pmslarray(locElementCnt))
    allocate(sstarray(locElementCnt))

    do i=1, locElementCnt
      arbSeqIndexList(i) = locElementBeg + (i - 1)
      streamflowarray(i) = -9.0
      sstarray(i) = 3.0 * arbSeqIndexList(i)
    enddo

    !print *,"OCN: ",localPet,"arbSeqIndices=", &
    !  arbSeqIndexList(1),arbSeqIndexList(locElementCnt)

    ! create DistGrid
    distgrid = ESMF_DistGridCreate(arbSeqIndexList=arbSeqIndexList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out


    ! create a LocationObject object for Fields
    locStreamIn = ESMF_LocStreamCreate(distgrid=distgrid, &
      coordSys=ESMF_COORDSYS_CART, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    locStreamOut = locStreamIn ! for now out same as in

    ! importable field: streamflow
    print*, "Beheen OCN size of streamflowarray", locElementCnt, localPet
    streamflowField = ESMF_FieldCreate(locStreamIn, &
                                   streamflowarray, &
                                ESMF_INDEX_DELOCAL, &
              datacopyflag=ESMF_DATACOPY_REFERENCE, &
                                 name="streamflow", &
                                               rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_Realize(importState, field=streamflowField, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out



    ! importable field: air_pressure_at_sea_level
    pmslField = ESMF_FieldCreate(locStreamIn, &
                                   pmslarray, &
                          ESMF_INDEX_DELOCAL, &
                                 name="pmsl", &
                                         rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(importState, field=pmslField, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! importable field: surface_net_downward_shortwave_flux
    !field = ESMF_FieldCreate(name="rsns", locStream=locStreamIn, &
    !  typekind=ESMF_TYPEKIND_R8, rc=rc)
    rsnsField = ESMF_FieldCreate(locStreamIn, &
                                   rsnsarray, &
                          ESMF_INDEX_DELOCAL, &
                                 name="rsns", &
                                         rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(importState, field=rsnsField, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    

    ! exportable field: sea_surface_temperature
    !field = ESMF_FieldCreate(name="sst", locStream=locStreamOut, &
    !  typekind=ESMF_TYPEKIND_R8, rc=rc)
    sstField = ESMF_FieldCreate(locStreamOut, &
                                    sstarray, &
                          ESMF_INDEX_DELOCAL, &
                                  name="sst", &
                                         rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(exportState, field=sstField, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out


    ! create a Grid object for export - waterlevel
    ! The order of the coordinates is (longitude, latitude) in symmetry with
    ! (x,y) for coordSys=ESMF_COORDSYS_CART.
    ! Also, the longitude is measured in degrees in the eastward direction, so
    ! the longitudes in your map should be -120 to -70.
  !        -------------------------------------------------
  ! 50.0 |   |    |    |    |    |    |    |    |    |       
  !      -------------------------------------------------
  ! 48.4 |   |    |    |    |    |    |    |    |    |       
  !      -------------------------------------------------
  ! 46.8 |   |    |    |    |    |    |    |    |    |       
  !      -------------------------------------------------
  ! 45.2 |   |    |    |    |    |    |    |    |    |       
  !      -------------------------------------------------
  ! 44.2 |   |    |    |    |    |    |    |    |    |       
  !      -------------------------------------------------
  ! 42.6 |   |    |    |    |    |    |    |    |    |       
  !      -------------------------------------------------
  ! 41.0 |   |    |    |    |    |    |    |    |    |       
  !      -------------------------------------------------
  ! 39.4 |   |    |    |    |    |    |    |    |    |       
  !      -------------------------------------------------
  ! 37.8 |   |    |    |    |    |    |    |    |    |       
  !      -------------------------------------------------
  ! 36.2 |   |    |    |    |    |    |    |    |    |       
  !      -------------------------------------------------
  ! 34.6 |   |    |    |    |    |    |    |    |    |       
  !      -------------------------------------------------
  ! 33.0 |   |    |    |    |    |    |    |    |    |       
  !      -------------------------------------------------
  ! 31.4 |   |    |    |    |    |    |    |    |    |       
  !      -------------------------------------------------
  ! 29.8 |   |    |    |    |    |    |    |    |    |       
  !      -------------------------------------------------
  ! 28.2 |   |    |    |    |    |    |    |    |    |       
  !      -------------------------------------------------
  ! 26.6 |   |    |    |    |    |    |    |    |    |       
  !      -------------------------------------------------
  ! 25.0 |   |    |    |    |    |    |    |    |    |       
  !      -------------------------------------------------
  ! -130.0 -123 -116 -109 -102  -95  -88  -81  -74  -67   -60

      ! create a Grid object for Fields
      !gridIn = ESMF_GridCreateNoPeriDimUfrm(maxIndex=(/5, 5/), &
      !minCornerCoord=(/60._ESMF_KIND_R8, 20._ESMF_KIND_R8/), &
      !maxCornerCoord=(/150._ESMF_KIND_R8, 200._ESMF_KIND_R8/), &
      !coordSys=ESMF_COORDSYS_CART,
      !staggerLocList=(/ESMF_STAGGERLOC_CENTER,ESMF_STAGGERLOC_CORNER/), &
      !rc=rc)
      !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      !    line=__LINE__, &
      !    file=__FILE__)) &
      !    return  ! bail out
      !gridOut = gridIn ! for now out same as in


    ! shortcut - Create a uniform Grid with no periodic dim and a regular
    ! distribution
    ! function ESMF_GridCreateNoPeriDimUfrmR(minIndex, maxIndex, &
    !     minCornerCoord, maxCornerCoord, &
    !     regDecomp, decompFlag, &
    !     coordSys, staggerLocList, petMap, name, rc)
    grid = ESMF_GridCreateNoPeriDimUfrm(maxIndex=(/10, 15/), &
      minCornerCoord=(/-96.0_ESMF_KIND_R8, 6.2_ESMF_KIND_R8/), &
      maxCornerCoord=(/-53.0_ESMF_KIND_R8, 47.2_ESMF_KIND_R8/), &
      staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
      name="OCN-GridIn", rc=rc)

    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! call NWM_GridGet(gridIn,ESMF_STAGGERLOC_CENTER,0,2)
    
    ! convert the Grid into a Mesh
    mesh = ESMF_MeshCreate(grid=grid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! get mesh dimenssion
    call ESMF_MeshGet(mesh, spatialDim=dimCount, &
              numOwnedElements=numOwnedElements, &
               numOwnedNodes=numOwnedNodes, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    write(msgString,*) "OCN mesh:   numOwnedElements=", numOwnedElements, &
      "numOwnedNodes=", numOwnedNodes, "dimCount=", dimCount, "pet", localPet
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !allocate(nodeOwners(nodeCount))
    !allocate(nodeCoords(2*nodeOwners))
    !allocate(elementCount())
    !call ESMF_MeshGet(meshOut, nodeOwners=nodeOwners, &
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    !print*,"Beheen in ocn:", "spatialDim=",dimCount, &
    !  "parametricDim=", parametricDim, &
      ! "nodeCount=", nodeCount, "elementCount=", elementCount, &
    !   "numOwnedElements=", numOwnedElements, "numOwnedNodes=", numOwnedNodes

    ! Create a Field from Mesh and Fortran array pointer
    ! function ESMF_FieldCreateMeshDataPtr<rank><type><kind>(mesh, & 
    !            farrayPtr, datacopyflag, meshloc, gridToFieldMap, & 
    ! exportable field: water level
    allocate(dataWL(numOwnedNodes))
    dataWL = 30.0
    field = ESMF_FieldCreate(name="wl", mesh=mesh, &
                                 farrayPtr=dataWL, &
             datacopyflag=ESMF_DATACOPY_REFERENCE, &
                        meshloc=ESMF_MESHLOC_NODE, &
                                              rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

 
    deallocate(arbSeqIndexList)

  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine SetClock(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_TimeInterval)       :: stabilityTimeStep

    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! initialize internal clock
    ! here: parent Clock and stability timeStep determine actual model timeStep
    !TODO: stabilityTimeStep should be read in from configuation
    !TODO: or computed from internal Grid information
    call ESMF_TimeIntervalSet(stabilityTimeStep, h=1, rc=rc) ! 5 minute steps
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetClock(model, clock, stabilityTimeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine

  subroutine printStreamflow(importState)

    type(ESMF_State)               :: importState
    type(ESMF_LocStream)           :: locstream
    character(len=ESMF_MAXSTR)     :: locName

    real(ESMF_KIND_R8), pointer    :: latArrayPtr(:)
    real(ESMF_KIND_R8), pointer    :: lonArrayPtr(:)
    real(ESMF_KIND_R8), pointer    :: flowRatePtr(:)
    integer(ESMF_KIND_I4), pointer :: linkArrayPtr(:)

    character(len=ESMF_MAXSTR), allocatable      :: itemNames(:)
    character (len=ESMF_MAXSTR)                  :: itemName
    type(ESMF_Field)                             :: itemField
    integer :: i, localPet, petCount, j, esmf_comm, k, itemCnt, rc
    type(ESMF_VM) :: vm


    call ESMF_StateGet(importState, itemCount=itemCnt, rc=rc)
    if (rc/=ESMF_SUCCESS) return
    allocate(itemNames(itemCnt))

    call ESMF_StateGet(importState, itemNameList=itemNames, rc=rc)

    do i=1, itemCnt
        itemName = trim(itemNames(i))
        call ESMF_StateGet(importState, itemName, itemField, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                               line=__LINE__, &
                               file=__FILE__)) &
                               return  ! bail out

        SELECT CASE (trim(itemName))
            CASE ('streamflow')
              ! Get a DE-local Fortran array pointer from a Field
              call ESMF_FieldGet(itemField, farrayPtr=flowRatePtr, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                                     line=__LINE__, &
                                     file=__FILE__)) &
                                     return  ! bail out

              call ESMF_FieldGet(itemField, locstream=locstream, vm=vm, rc=rc)
              ! get the vm of this field
              call ESMF_VMGet(vm=vm, localPet=localPet, petCount=petCount, &
                              mpiCommunicator=esmf_comm, rc=rc)

              ! field values 
              !call ESMF_LocStreamGetKey(locstream, "Lat", farray=latArrayPtr, rc=rc)
              !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              !                       line=__LINE__, &
              !                       file=__FILE__)) &
              !                       return  ! bail out

              !call ESMF_LocStreamGetKey(locstream, "Lon", farray=lonArrayPtr, rc=rc)
              !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              !                       line=__LINE__, &
              !                       file=__FILE__)) &
              !                       return  ! bail out

              !call ESMF_LocStreamGetKey(locstream, "link", farray=linkArrayPtr, rc=rc)
              !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              !                       line=__LINE__, &
              !                       file=__FILE__)) &
              !                       return  ! bail out

              do j=0, petCount
                  if (localPet == j) then
                    print*, "OCN getting flowrate for DE:", localPet
                    print*, "flowrate: ",flowRatePtr
                    print*, "element count:", size(flowRatePtr)  
                  endif
                  call MPI_Barrier(esmf_comm, rc)
                  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                      line=__LINE__, &
                      file=__FILE__)) &
                  return  ! bail out
              enddo

        END SELECT
    enddo
    deallocate(itemNames)

  end subroutine


  subroutine printWaterlevel(exportState, dt)
    type(ESMF_State)     :: exportState
    real(ESMF_KIND_R8)   :: dt

    character(len=ESMF_MAXSTR), allocatable      :: itemNames(:)
    character (len=ESMF_MAXSTR)                  :: itemName
    type(ESMF_Field)                             :: itemField
    integer                                      :: itemCnt
    integer :: i, localPet, petCount, j, esmf_comm, k, rc
    integer :: dimCount, numOwnedElements, numOwnedNodes
    type(ESMF_VM)                                :: vm
    type(ESMF_Mesh)                              :: mesh
    real(ESMF_KIND_R8), dimension(:), pointer    :: wlPtr => null()

    rc = ESMF_SUCCESS

    call ESMF_StateGet(exportState, itemCount=itemCnt, rc=rc)
    if (rc/=ESMF_SUCCESS) return

    allocate(itemNames(itemCnt))
    call ESMF_StateGet(exportState, itemNameList=itemNames, rc=rc)

    do i=1, itemCnt
      itemName = trim(itemNames(i))
      call ESMF_StateGet(exportState, itemName, itemField, rc=rc)

      SELECT CASE (itemName)

        CASE ('wl')

          call ESMF_FieldGet(itemField, mesh=mesh, vm=vm, rc=rc)
          ! get mesh dimenssion
          call ESMF_MeshGet(mesh, spatialDim=dimCount, &
                    numOwnedElements=numOwnedElements, &
                     numOwnedNodes=numOwnedNodes, rc=rc)
         
          call ESMF_VMGet(vm=vm, localPet=localPet, petCount=petCount, &
                                     mpiCommunicator=esmf_comm, rc=rc)

          nullify(wlPtr)
          call ESMF_FieldGet(itemField, farrayPtr=wlPtr, rc=rc)
          do j=0, petCount
            if (localPet == j) then
              print*, "OCN setting watelevel,nodes, pet:", numOwnedNodes,localPet
              print*, wlPtr
            endif
            call MPI_Barrier(esmf_comm, rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, &
                file=__FILE__)) &
                return  ! bail out
          enddo

      END SELECT
    enddo
    deallocate(wlPtr)
    deallocate(itemNames)

    ! end test
  end subroutine
    

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Time)             :: currTime
    type(ESMF_TimeInterval)     :: timeStep
    character(len=160)          :: msgString

    real(ESMF_KIND_R8)          :: dt

#define NUOPC_TRACE__OFF
#ifdef NUOPC_TRACE
    call ESMF_TraceRegionEnter("OCN:ModelAdvance")
#endif
    
    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    
    ! Because of the way that the internal Clock was set in SetClock(),
    ! its timeStep is likely smaller than the parent timeStep. As a consequence
    ! the time interval covered by a single parent timeStep will result in 
    ! multiple calls to the ModelAdvance() routine. Every time the currTime
    ! will come in by one internal timeStep advanced. This goes until the
    ! stopTime of the internal Clock has been reached.
    
    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing OCN from: ", unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
   
    ! print waterlevel
    dt = TimeIntervalGetReal(timeStep)
    call printWaterlevel(exportState,dt)

    ! print streamflow 
    call printStreamflow(importState)

      
    call ESMF_TimePrint(currTime + timeStep, &
      preString="---------------------> to: ", unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

#ifdef NUOPC_TRACE
    call ESMF_TraceRegionExit("OCN:ModelAdvance")
#endif
  end subroutine

end module
