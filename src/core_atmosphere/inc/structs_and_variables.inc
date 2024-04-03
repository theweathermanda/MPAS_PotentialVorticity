   subroutine atm_generate_pool_mesh(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_io, only : MPAS_REAL_FILLVAL, MPAS_INT_FILLVAL, MPAS_CHAR_FILLVAL
      implicit none
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), dimension(:), pointer :: r0Ptr
      type (field1DReal), dimension(:), pointer :: r1Ptr
      type (field2DReal), dimension(:), pointer :: r2Ptr
      type (field3DReal), dimension(:), pointer :: r3Ptr
      type (field4DReal), dimension(:), pointer :: r4Ptr
      type (field5DReal), dimension(:), pointer :: r5Ptr
      type (field0DInteger), dimension(:), pointer :: i0Ptr
      type (field1DInteger), dimension(:), pointer :: i1Ptr
      type (field2DInteger), dimension(:), pointer :: i2Ptr
      type (field3DInteger), dimension(:), pointer :: i3Ptr
      type (field0DChar), dimension(:), pointer :: c0Ptr
      type (field1DChar), dimension(:), pointer :: c1Ptr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: mp_kessler_inActive
      logical, pointer :: mp_thompson_inActive
      logical, pointer :: mp_wsm6_inActive
      logical, pointer :: cu_grell_freitas_inActive
      logical, pointer :: cu_kain_fritsch_inActive
      logical, pointer :: cu_tiedtke_inActive
      logical, pointer :: bl_ysu_inActive
      logical, pointer :: bl_mynn_inActive
      logical, pointer :: iauActive
      logical, pointer :: limited_areaActive
      logical, pointer :: jedi_daActive
      logical, pointer :: tendenciesActive


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'mp_kessler_inActive', mp_kessler_inActive)
      call mpas_pool_get_package(packagePool, 'mp_thompson_inActive', mp_thompson_inActive)
      call mpas_pool_get_package(packagePool, 'mp_wsm6_inActive', mp_wsm6_inActive)
      call mpas_pool_get_package(packagePool, 'cu_grell_freitas_inActive', cu_grell_freitas_inActive)
      call mpas_pool_get_package(packagePool, 'cu_kain_fritsch_inActive', cu_kain_fritsch_inActive)
      call mpas_pool_get_package(packagePool, 'cu_tiedtke_inActive', cu_tiedtke_inActive)
      call mpas_pool_get_package(packagePool, 'bl_ysu_inActive', bl_ysu_inActive)
      call mpas_pool_get_package(packagePool, 'bl_mynn_inActive', bl_mynn_inActive)
      call mpas_pool_get_package(packagePool, 'iauActive', iauActive)
      call mpas_pool_get_package(packagePool, 'limited_areaActive', limited_areaActive)
      call mpas_pool_get_package(packagePool, 'jedi_daActive', jedi_daActive)
      call mpas_pool_get_package(packagePool, 'tendenciesActive', tendenciesActive)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'mesh', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'mesh', newSubPool)

! Define variable latCell
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'latCell'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'rad')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Latitude of cells')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'latCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'latCell', r1Ptr)

! Define variable lonCell
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'lonCell'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'rad')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Longitude of cells')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lonCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lonCell', r1Ptr)

! Define variable xCell
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'xCell'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Cartesian x-coordinate of cells')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'xCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'xCell', r1Ptr)

! Define variable yCell
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'yCell'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Cartesian y-coordinate of cells')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'yCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'yCell', r1Ptr)

! Define variable zCell
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'zCell'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Cartesian z-coordinate of cells')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'zCell', r1Ptr)

! Define variable indexToCellID
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'indexToCellID'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .false.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'Mapping from local array index to global cell ID')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'indexToCellID', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'indexToCellID', i1Ptr)

! Define variable latEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'latEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'rad')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Latitude of edges')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'latEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'latEdge', r1Ptr)

! Define variable lonEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'lonEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'rad')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Longitude of edges')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lonEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lonEdge', r1Ptr)

! Define variable xEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'xEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Cartesian x-coordinate of edges')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'xEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'xEdge', r1Ptr)

! Define variable yEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'yEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Cartesian y-coordinate of edges')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'yEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'yEdge', r1Ptr)

! Define variable zEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'zEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Cartesian z-coordinate of edges')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'zEdge', r1Ptr)

! Define variable indexToEdgeID
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'indexToEdgeID'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .false.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nEdges'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'Mapping from local array index to global edge ID')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'indexToEdgeID', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'indexToEdgeID', i1Ptr)

! Define variable latVertex
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'latVertex'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'rad')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Latitude of vertices')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'latVertex', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'latVertex', r1Ptr)

! Define variable lonVertex
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'lonVertex'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'rad')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Longitude of vertices')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lonVertex', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lonVertex', r1Ptr)

! Define variable xVertex
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'xVertex'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Cartesian x-coordinate of vertices')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'xVertex', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'xVertex', r1Ptr)

! Define variable yVertex
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'yVertex'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Cartesian y-coordinate of vertices')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'yVertex', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'yVertex', r1Ptr)

! Define variable zVertex
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'zVertex'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Cartesian z-coordinate of vertices')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zVertex', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'zVertex', r1Ptr)

! Define variable indexToVertexID
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'indexToVertexID'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .false.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nVertices'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'Mapping from local array index to global vertex ID')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'indexToVertexID', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'indexToVertexID', i1Ptr)

! Define variable cellsOnEdge
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'cellsOnEdge'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .false.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'TWO'
      i2Ptr(1) % dimNames(2) = 'nEdges'
      i2Ptr(1) % defaultValue = 0
      allocate(i2Ptr(1) % attLists(1))
      allocate(i2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'long_name', 'IDs of cells divided by an edge')
      i2Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cellsOnEdge', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'cellsOnEdge', i2Ptr)

! Define variable nEdgesOnCell
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'nEdgesOnCell'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .false.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'Number of edges forming the boundary of a cell')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'nEdgesOnCell', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'nEdgesOnCell', i1Ptr)

! Define variable nEdgesOnEdge
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'nEdgesOnEdge'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .false.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nEdges'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'Number of edges involved in reconstruction of tangential '&
//'velocity for an edge')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'nEdgesOnEdge', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'nEdgesOnEdge', i1Ptr)

! Define variable edgesOnCell
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'edgesOnCell'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .false.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'maxEdges'
      i2Ptr(1) % dimNames(2) = 'nCells'
      i2Ptr(1) % defaultValue = 0
      allocate(i2Ptr(1) % attLists(1))
      allocate(i2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'long_name', 'IDs of edges forming the boundary of a cell')
      i2Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'edgesOnCell', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'edgesOnCell', i2Ptr)

! Define variable edgesOnEdge
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'edgesOnEdge'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .false.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'maxEdges2'
      i2Ptr(1) % dimNames(2) = 'nEdges'
      i2Ptr(1) % defaultValue = 0
      allocate(i2Ptr(1) % attLists(1))
      allocate(i2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'long_name', 'IDs of edges involved in reconstruction of tangential '&
//'velocity for an edge')
      i2Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'edgesOnEdge', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'edgesOnEdge', i2Ptr)

! Define variable weightsOnEdge
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'weightsOnEdge'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'maxEdges2'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Weights used in reconstruction of tangential velocity for '&
//'an edge')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'weightsOnEdge', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'weightsOnEdge', r2Ptr)

! Define variable dvEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'dvEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Spherical distance between vertex endpoints of an edge')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dvEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'dvEdge', r1Ptr)

! Define variable invDvEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'invDvEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Inverse distance between vertex endpoints of an edge')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'invDvEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'invDvEdge', r1Ptr)

! Define variable dcEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'dcEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Spherical distance between cells separated by an edge')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dcEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'dcEdge', r1Ptr)

! Define variable invDcEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'invDcEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Inverse distance between cells separated by an edge')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'invDcEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'invDcEdge', r1Ptr)

! Define variable angleEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'angleEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'rad')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Angle between local north and the positive tangential '&
//'direction of an edge')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'angleEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'angleEdge', r1Ptr)

! Define variable areaCell
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'areaCell'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm^2')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Spherical area of a Voronoi cell')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'areaCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'areaCell', r1Ptr)

! Define variable invAreaCell
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'invAreaCell'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Inverse of Voronoi cell area')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'invAreaCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'invAreaCell', r1Ptr)

! Define variable areaTriangle
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'areaTriangle'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm^2')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Spherical area of a Delaunay triangle')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'areaTriangle', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'areaTriangle', r1Ptr)

! Define variable invAreaTriangle
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'invAreaTriangle'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Inverse area of a Delaunay triangle')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'invAreaTriangle', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'invAreaTriangle', r1Ptr)

! Define variable edgeNormalVectors
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'edgeNormalVectors'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'R3'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Cartesian components of the vector normal to an edge and '&
//'tangential to the surface of the sphere')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'edgeNormalVectors', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'edgeNormalVectors', r2Ptr)

! Define variable localVerticalUnitVectors
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'localVerticalUnitVectors'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'R3'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Cartesian components of the vector pointing in the local '&
//'vertical direction for a cell')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'localVerticalUnitVectors', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'localVerticalUnitVectors', r2Ptr)

! Define variable cellTangentPlane
      allocate(r3Ptr(1))

! Setting up time level 1
      r3Ptr(1) % fieldName = 'cellTangentPlane'
      r3Ptr(1) % isVarArray = .false.
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .false.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.
! Setting up dimensions
      r3Ptr(1) % dimNames(1) = 'R3'
      r3Ptr(1) % dimNames(2) = 'TWO'
      r3Ptr(1) % dimNames(3) = 'nCells'
      r3Ptr(1) % defaultValue = 0.0
      allocate(r3Ptr(1) % attLists(1))
      allocate(r3Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'long_name', 'Components of a pair of vectors defining the tangent plane '&
//'at a cell')
      r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r3Ptr(1) % block => block

      r3Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cellTangentPlane', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'cellTangentPlane', r3Ptr)

! Define variable cellsOnCell
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'cellsOnCell'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .false.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'maxEdges'
      i2Ptr(1) % dimNames(2) = 'nCells'
      i2Ptr(1) % defaultValue = 0
      allocate(i2Ptr(1) % attLists(1))
      allocate(i2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'long_name', 'IDs of cells neighboring a cell')
      i2Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cellsOnCell', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'cellsOnCell', i2Ptr)

! Define variable verticesOnCell
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'verticesOnCell'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .false.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'maxEdges'
      i2Ptr(1) % dimNames(2) = 'nCells'
      i2Ptr(1) % defaultValue = 0
      allocate(i2Ptr(1) % attLists(1))
      allocate(i2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'long_name', 'IDs of vertices (corner points) of a cell')
      i2Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'verticesOnCell', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'verticesOnCell', i2Ptr)

! Define variable verticesOnEdge
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'verticesOnEdge'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .false.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'TWO'
      i2Ptr(1) % dimNames(2) = 'nEdges'
      i2Ptr(1) % defaultValue = 0
      allocate(i2Ptr(1) % attLists(1))
      allocate(i2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'long_name', 'IDs of the two vertex endpoints of an edge')
      i2Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'verticesOnEdge', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'verticesOnEdge', i2Ptr)

! Define variable edgesOnVertex
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'edgesOnVertex'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .false.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'vertexDegree'
      i2Ptr(1) % dimNames(2) = 'nVertices'
      i2Ptr(1) % defaultValue = 0
      allocate(i2Ptr(1) % attLists(1))
      allocate(i2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'long_name', 'IDs of the edges that meet at a vertex')
      i2Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'edgesOnVertex', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'edgesOnVertex', i2Ptr)

! Define variable cellsOnVertex
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'cellsOnVertex'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .false.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'vertexDegree'
      i2Ptr(1) % dimNames(2) = 'nVertices'
      i2Ptr(1) % defaultValue = 0
      allocate(i2Ptr(1) % attLists(1))
      allocate(i2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'long_name', 'IDs of the cells that meet at a vertex')
      i2Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cellsOnVertex', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'cellsOnVertex', i2Ptr)

! Define variable kiteAreasOnVertex
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'kiteAreasOnVertex'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'vertexDegree'
      r2Ptr(1) % dimNames(2) = 'nVertices'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm^2')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Intersection areas between primal (Voronoi) and dual '&
//'(triangular) mesh cells')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'kiteAreasOnVertex', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'kiteAreasOnVertex', r2Ptr)

! Define variable fEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'fEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Coriolis parameter at an edge')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'fEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'fEdge', r1Ptr)

! Define variable fVertex
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'fVertex'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Coriolis parameter at a vertex')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'fVertex', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'fVertex', r1Ptr)

! Define variable meshDensity
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'meshDensity'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Mesh density function (used when generating the mesh) '&
//'evaluated at a cell')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'meshDensity', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'meshDensity', r1Ptr)

! Define variable meshScalingDel2
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'meshScalingDel2'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Scaling coefficient for $\nabla^2$ eddy viscosity')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'meshScalingDel2', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'meshScalingDel2', r1Ptr)

! Define variable meshScalingDel4
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'meshScalingDel4'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Scaling coefficient for $\nabla^4$ eddy hyper-viscosity')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'meshScalingDel4', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'meshScalingDel4', r1Ptr)

! Define variable meshScalingRegionalCell
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'meshScalingRegionalCell'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Cell-centered Scaling coefficient for relaxation zone')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'meshScalingRegionalCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'meshScalingRegionalCell', r1Ptr)

! Define variable meshScalingRegionalEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'meshScalingRegionalEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Edge-centered Scaling coefficient for relaxation zone')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'meshScalingRegionalEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'meshScalingRegionalEdge', r1Ptr)

! Define variable edgesOnVertex_sign
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'edgesOnVertex_sign'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'vertexDegree'
      r2Ptr(1) % dimNames(2) = 'nVertices'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Sign for edges incident with a vertex: positive for '&
//'positive inward tengential velocity')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'edgesOnVertex_sign', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'edgesOnVertex_sign', r2Ptr)

! Define variable edgesOnCell_sign
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'edgesOnCell_sign'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'maxEdges'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Sign for edges surrounding a cell: positive for positive '&
//'outward normal velocity')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'edgesOnCell_sign', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'edgesOnCell_sign', r2Ptr)

! Define variable kiteForCell
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'kiteForCell'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .false.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'maxEdges'
      i2Ptr(1) % dimNames(2) = 'nCells'
      i2Ptr(1) % defaultValue = 0
      allocate(i2Ptr(1) % attLists(1))
      allocate(i2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'long_name', 'Index of kite in kiteAreasOnVertex that lies within a cell '&
//'for each of verticesOnCell')
      i2Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'kiteForCell', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'kiteForCell', i2Ptr)

! Define variable cf1
      allocate(r0Ptr(1))

! Setting up time level 1
      r0Ptr(1) % fieldName = 'cf1'
      r0Ptr(1) % isVarArray = .false.
      r0Ptr(1) % isDecomposed = .false.
      r0Ptr(1) % hasTimeDimension = .false.
      r0Ptr(1) % defaultValue = 0.0
      r0Ptr(1) % scalar = 0.0
      allocate(r0Ptr(1) % attLists(1))
      allocate(r0Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r0Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r0Ptr(1) % attLists(1) % attList, 'long_name', 'Surface interpolation weight for level k=1 value')
      r0Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r0Ptr(1) % block => block

      r0Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cf1', r0Ptr)
      call mpas_pool_add_field(block % allFields, 'cf1', r0Ptr)

! Define variable cf2
      allocate(r0Ptr(1))

! Setting up time level 1
      r0Ptr(1) % fieldName = 'cf2'
      r0Ptr(1) % isVarArray = .false.
      r0Ptr(1) % isDecomposed = .false.
      r0Ptr(1) % hasTimeDimension = .false.
      r0Ptr(1) % defaultValue = 0.0
      r0Ptr(1) % scalar = 0.0
      allocate(r0Ptr(1) % attLists(1))
      allocate(r0Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r0Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r0Ptr(1) % attLists(1) % attList, 'long_name', 'Surface interpolation weight for level k=2 value')
      r0Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r0Ptr(1) % block => block

      r0Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cf2', r0Ptr)
      call mpas_pool_add_field(block % allFields, 'cf2', r0Ptr)

! Define variable cf3
      allocate(r0Ptr(1))

! Setting up time level 1
      r0Ptr(1) % fieldName = 'cf3'
      r0Ptr(1) % isVarArray = .false.
      r0Ptr(1) % isDecomposed = .false.
      r0Ptr(1) % hasTimeDimension = .false.
      r0Ptr(1) % defaultValue = 0.0
      r0Ptr(1) % scalar = 0.0
      allocate(r0Ptr(1) % attLists(1))
      allocate(r0Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r0Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r0Ptr(1) % attLists(1) % attList, 'long_name', 'Surface interpolation weight for level k=3 value')
      r0Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r0Ptr(1) % block => block

      r0Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cf3', r0Ptr)
      call mpas_pool_add_field(block % allFields, 'cf3', r0Ptr)

! Define variable zgrid
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'zgrid'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm MSL')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Geometric height of layer interfaces')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zgrid', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'zgrid', r2Ptr)

! Define variable rdzw
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'rdzw'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .false.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertLevels'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Reciprocal dzw')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rdzw', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'rdzw', r1Ptr)

! Define variable dzu
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'dzu'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .false.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertLevels'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'd(zeta) at w levels')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dzu', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'dzu', r1Ptr)

! Define variable rdzu
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'rdzu'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .false.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertLevels'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Reciprocal dzu')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rdzu', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'rdzu', r1Ptr)

! Define variable fzm
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'fzm'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .false.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertLevels'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Weight for linear interpolation to w(k) point for u(k) '&
//'level variable')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'fzm', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'fzm', r1Ptr)

! Define variable fzp
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'fzp'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .false.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertLevels'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Weight for linear interpolation to w(k) point for u(k-1) '&
//'level variable')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'fzp', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'fzp', r1Ptr)

! Define variable zxu
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'zxu'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'dz/dx on horizontal coordinate surfaces at u levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zxu', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'zxu', r2Ptr)

! Define variable zz
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'zz'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'd(zeta)/dz, vertical metric term')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zz', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'zz', r2Ptr)

! Define variable zb
      allocate(r3Ptr(1))

! Setting up time level 1
      r3Ptr(1) % fieldName = 'zb'
      r3Ptr(1) % isVarArray = .false.
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .false.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.
! Setting up dimensions
      r3Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r3Ptr(1) % dimNames(2) = 'TWO'
      r3Ptr(1) % dimNames(3) = 'nEdges'
      r3Ptr(1) % defaultValue = 0.0
      allocate(r3Ptr(1) % attLists(1))
      allocate(r3Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'long_name', 'Coefficients for contribution from u to omega diagnosis, '&
//'edge-oriented')
      r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r3Ptr(1) % block => block

      r3Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zb', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'zb', r3Ptr)

! Define variable zb3
      allocate(r3Ptr(1))

! Setting up time level 1
      r3Ptr(1) % fieldName = 'zb3'
      r3Ptr(1) % isVarArray = .false.
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .false.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.
! Setting up dimensions
      r3Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r3Ptr(1) % dimNames(2) = 'TWO'
      r3Ptr(1) % dimNames(3) = 'nEdges'
      r3Ptr(1) % defaultValue = 0.0
      allocate(r3Ptr(1) % attLists(1))
      allocate(r3Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'long_name', 'Coefficients for 3rd-order correction to contribution from '&
//'u to omega diagnosis, edge-oriented')
      r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r3Ptr(1) % block => block

      r3Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zb3', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'zb3', r3Ptr)

! Define variable zb_cell
      allocate(r3Ptr(1))

! Setting up time level 1
      r3Ptr(1) % fieldName = 'zb_cell'
      r3Ptr(1) % isVarArray = .false.
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .false.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.
! Setting up dimensions
      r3Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r3Ptr(1) % dimNames(2) = 'maxEdges'
      r3Ptr(1) % dimNames(3) = 'nCells'
      r3Ptr(1) % defaultValue = 0.0
      allocate(r3Ptr(1) % attLists(1))
      allocate(r3Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'long_name', 'Coefficients for contribution from u to omega diagnosis, '&
//'cell-oriented')
      r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r3Ptr(1) % block => block

      r3Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zb_cell', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'zb_cell', r3Ptr)

! Define variable zb3_cell
      allocate(r3Ptr(1))

! Setting up time level 1
      r3Ptr(1) % fieldName = 'zb3_cell'
      r3Ptr(1) % isVarArray = .false.
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .false.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.
! Setting up dimensions
      r3Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r3Ptr(1) % dimNames(2) = 'maxEdges'
      r3Ptr(1) % dimNames(3) = 'nCells'
      r3Ptr(1) % defaultValue = 0.0
      allocate(r3Ptr(1) % attLists(1))
      allocate(r3Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'long_name', 'Coefficients for 3rd-order correction to contribution from '&
//'u to omega diagnosis, cell-oriented')
      r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r3Ptr(1) % block => block

      r3Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zb3_cell', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'zb3_cell', r3Ptr)

! Define variable dss
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dss'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'w-damping coefficient')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dss', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dss', r2Ptr)

! Define variable u_init
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'u_init'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .false.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertLevels'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'u reference profile')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'u_init', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'u_init', r1Ptr)

! Define variable v_init
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'v_init'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .false.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertLevels'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'v reference profile')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'v_init', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'v_init', r1Ptr)

! Define variable t_init
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 't_init'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'theta reference profile')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 't_init', r2Ptr)
      call mpas_pool_add_field(block % allFields, 't_init', r2Ptr)

! Define variable qv_init
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'qv_init'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .false.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertLevels'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'qv reference profile')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'qv_init', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'qv_init', r1Ptr)

! Define variable deriv_two
      allocate(r3Ptr(1))

! Setting up time level 1
      r3Ptr(1) % fieldName = 'deriv_two'
      r3Ptr(1) % isVarArray = .false.
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .false.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.
! Setting up dimensions
      r3Ptr(1) % dimNames(1) = 'FIFTEEN'
      r3Ptr(1) % dimNames(2) = 'TWO'
      r3Ptr(1) % dimNames(3) = 'nEdges'
      r3Ptr(1) % defaultValue = 0.0
      allocate(r3Ptr(1) % attLists(1))
      allocate(r3Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'long_name', 'weights for cell-centered second derivative, normal to '&
//'edge, for transport scheme')
      r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r3Ptr(1) % block => block

      r3Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'deriv_two', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'deriv_two', r3Ptr)

! Define variable adv_coefs
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'adv_coefs'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'FIFTEEN'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Weighting coefficients used for reconstructing cell-based '&
//'fields at edges')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'adv_coefs', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'adv_coefs', r2Ptr)

! Define variable adv_coefs_3rd
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'adv_coefs_3rd'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'FIFTEEN'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Weighting coefficients used for reconstructing cell-based '&
//'fields at edges')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'adv_coefs_3rd', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'adv_coefs_3rd', r2Ptr)

! Define variable advCellsForEdge
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'advCellsForEdge'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .false.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'FIFTEEN'
      i2Ptr(1) % dimNames(2) = 'nEdges'
      i2Ptr(1) % defaultValue = 0
      allocate(i2Ptr(1) % attLists(1))
      allocate(i2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'long_name', 'Cells used to reconstruct a cell-based field at an edge')
      i2Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'advCellsForEdge', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'advCellsForEdge', i2Ptr)

! Define variable nAdvCellsForEdge
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'nAdvCellsForEdge'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .false.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nEdges'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'Number of cells used to reconstruct a cell-based field at '&
//'an edge')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'nAdvCellsForEdge', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'nAdvCellsForEdge', i1Ptr)

! Define variable defc_a
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'defc_a'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'maxEdges'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coefficients for computing the off-diagonal components of '&
//'the horizontal deformation')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'defc_a', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'defc_a', r2Ptr)

! Define variable defc_b
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'defc_b'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'maxEdges'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coefficients for computing the diagonal components of the '&
//'horizontal deformation')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'defc_b', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'defc_b', r2Ptr)

! Define variable coeffs_reconstruct
      allocate(r3Ptr(1))

! Setting up time level 1
      r3Ptr(1) % fieldName = 'coeffs_reconstruct'
      r3Ptr(1) % isVarArray = .false.
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .false.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.
! Setting up dimensions
      r3Ptr(1) % dimNames(1) = 'R3'
      r3Ptr(1) % dimNames(2) = 'maxEdges'
      r3Ptr(1) % dimNames(3) = 'nCells'
      r3Ptr(1) % defaultValue = 0.0
      allocate(r3Ptr(1) % attLists(1))
      allocate(r3Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'long_name', 'Coefficients to reconstruct velocity vectors at cell '&
//'centers')
      r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r3Ptr(1) % block => block

      r3Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'coeffs_reconstruct', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'coeffs_reconstruct', r3Ptr)

! Define variable east
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'east'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'R3'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Cartesian components of the local unit vector pointing '&
//'east')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'east', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'east', r2Ptr)

! Define variable north
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'north'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'R3'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Cartesian components of the unit vector pointing north')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'north', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'north', r2Ptr)

! Define variable bdyMaskCell
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'bdyMaskCell'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .false.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'Limited-area specified/relaxation zone index for cells')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'bdyMaskCell', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'bdyMaskCell', i1Ptr)

! Define variable bdyMaskEdge
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'bdyMaskEdge'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .false.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nEdges'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'Limited-area specified/relaxation zone index for edges')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'bdyMaskEdge', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'bdyMaskEdge', i1Ptr)

! Define variable bdyMaskVertex
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'bdyMaskVertex'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .false.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nVertices'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'Limited-area specified/relaxation zone index for vertices')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'bdyMaskVertex', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'bdyMaskVertex', i1Ptr)

! Define variable specZoneMaskCell
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'specZoneMaskCell'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', '0/1 mask on cells, defined as 1 for cells in the '&
//'limited-area specified zone')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'specZoneMaskCell', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'specZoneMaskCell', r1Ptr)

! Define variable specZoneMaskEdge
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'specZoneMaskEdge'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nEdges'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', '0/1 mask on edges, defined as 1 for edges in the '&
//'limited-area specified zone')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'specZoneMaskEdge', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'specZoneMaskEdge', r1Ptr)

! Define variable specZoneMaskVertex
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'specZoneMaskVertex'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', '0/1 mask on vertices, defined as 1 for vertices in the '&
//'limited-area specified zone')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'specZoneMaskVertex', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'specZoneMaskVertex', r1Ptr)

! Define variable nearestRelaxationCell
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'nearestRelaxationCell'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .false.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'For cells in the specified zone, gives the index of the '&
//'nearest cell in the relaxation zone')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'nearestRelaxationCell', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'nearestRelaxationCell', i1Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
         call mpas_pool_add_config(newSubPool, 'is_periodic', block % domain % is_periodic)
         call mpas_pool_add_config(newSubPool, 'x_period', block % domain % x_period)
         call mpas_pool_add_config(newSubPool, 'y_period', block % domain % y_period)
      end if

   end subroutine atm_generate_pool_mesh


   subroutine atm_generate_pool_state(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_io, only : MPAS_REAL_FILLVAL, MPAS_INT_FILLVAL, MPAS_CHAR_FILLVAL
      implicit none
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), dimension(:), pointer :: r0Ptr
      type (field1DReal), dimension(:), pointer :: r1Ptr
      type (field2DReal), dimension(:), pointer :: r2Ptr
      type (field3DReal), dimension(:), pointer :: r3Ptr
      type (field4DReal), dimension(:), pointer :: r4Ptr
      type (field5DReal), dimension(:), pointer :: r5Ptr
      type (field0DInteger), dimension(:), pointer :: i0Ptr
      type (field1DInteger), dimension(:), pointer :: i1Ptr
      type (field2DInteger), dimension(:), pointer :: i2Ptr
      type (field3DInteger), dimension(:), pointer :: i3Ptr
      type (field0DChar), dimension(:), pointer :: c0Ptr
      type (field1DChar), dimension(:), pointer :: c1Ptr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: mp_kessler_inActive
      logical, pointer :: mp_thompson_inActive
      logical, pointer :: mp_wsm6_inActive
      logical, pointer :: cu_grell_freitas_inActive
      logical, pointer :: cu_kain_fritsch_inActive
      logical, pointer :: cu_tiedtke_inActive
      logical, pointer :: bl_ysu_inActive
      logical, pointer :: bl_mynn_inActive
      logical, pointer :: iauActive
      logical, pointer :: limited_areaActive
      logical, pointer :: jedi_daActive
      logical, pointer :: tendenciesActive


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'mp_kessler_inActive', mp_kessler_inActive)
      call mpas_pool_get_package(packagePool, 'mp_thompson_inActive', mp_thompson_inActive)
      call mpas_pool_get_package(packagePool, 'mp_wsm6_inActive', mp_wsm6_inActive)
      call mpas_pool_get_package(packagePool, 'cu_grell_freitas_inActive', cu_grell_freitas_inActive)
      call mpas_pool_get_package(packagePool, 'cu_kain_fritsch_inActive', cu_kain_fritsch_inActive)
      call mpas_pool_get_package(packagePool, 'cu_tiedtke_inActive', cu_tiedtke_inActive)
      call mpas_pool_get_package(packagePool, 'bl_ysu_inActive', bl_ysu_inActive)
      call mpas_pool_get_package(packagePool, 'bl_mynn_inActive', bl_mynn_inActive)
      call mpas_pool_get_package(packagePool, 'iauActive', iauActive)
      call mpas_pool_get_package(packagePool, 'limited_areaActive', limited_areaActive)
      call mpas_pool_get_package(packagePool, 'jedi_daActive', jedi_daActive)
      call mpas_pool_get_package(packagePool, 'tendenciesActive', tendenciesActive)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'state', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'state', newSubPool)

! Define var array scalars
      allocate(r3Ptr(2))
      index_counter = 0
      group_counter = -1
      group_start = -1
      group_started = .false.

! Starting group moist
! Define constituent var qv
! My Packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_qv', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
         end if
         group_started = .true.
      end if
! Define constituent var qc
! My packages are bl_mynn_in;bl_ysu_in;cu_tiedtke_in;mp_kessler_in;mp_thompson_in;mp_wsm6_in
      if (bl_mynn_inActive .or. bl_ysu_inActive .or. cu_tiedtke_inActive .or. mp_kessler_inActive .or. mp_thompson_inActive .or. &
mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qc', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qc', -1)
         end if
      end if
! Define constituent var qr
! My packages are mp_kessler_in;mp_thompson_in;mp_wsm6_in
      if (mp_kessler_inActive .or. mp_thompson_inActive .or. mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qr', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qr', -1)
         end if
      end if
! Define constituent var qi
! My packages are bl_mynn_in;bl_ysu_in;cu_tiedtke_in;mp_thompson_in;mp_wsm6_in
      if (bl_mynn_inActive .or. bl_ysu_inActive .or. cu_tiedtke_inActive .or. mp_thompson_inActive .or. mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qi', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qi', -1)
         end if
      end if
! Define constituent var qs
! My packages are mp_thompson_in;mp_wsm6_in
      if (mp_thompson_inActive .or. mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qs', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qs', -1)
         end if
      end if
! Define constituent var qg
! My packages are mp_thompson_in;mp_wsm6_in
      if (mp_thompson_inActive .or. mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qg', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qg', -1)
         end if
      end if
         if (.not. group_started) then
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', -1)
               call mpas_pool_add_dimension(newSubPool, 'moist_end', -1)
            end if
         else
            group_started = .false.
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_end', index_counter)
            end if
         end if
! End of group       
! Starting group number
! Define constituent var ni
! My Packages are bl_mynn_in;mp_thompson_in
      if (bl_mynn_inActive .or. mp_thompson_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_ni', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'number_start', group_start)
            end if
            group_started = .true.
         end if
      else
           if (associated(newSubPool)) then
              call mpas_pool_add_dimension(newSubPool, 'index_ni', -1)
           end if
      end if
! Define constituent var nr
! My packages are mp_thompson_in
      if (mp_thompson_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_nr', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'number_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_nr', -1)
         end if
      end if
         if (.not. group_started) then
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'number_start', -1)
               call mpas_pool_add_dimension(newSubPool, 'number_end', -1)
            end if
         else
            group_started = .false.
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'number_end', index_counter)
            end if
         end if
! End of group       

      numConstituents = index_counter
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'num_scalars', numConstituents)
      end if
! Defining time level 1
      allocate( r3Ptr(1) % constituentNames(numConstituents) )
      r3Ptr(1) % fieldName = 'scalars'
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .true.
      r3Ptr(1) % isVarArray = .true.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qv', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qv'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qc', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qc'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qr', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qi', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qi'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qs', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qs'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qg', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qg'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_ni', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'ni'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_nr', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'nr'
      end if

! Setup dimensions for       
      r3Ptr(1) % dimNames(1) = 'num_scalars'
      r3Ptr(1) % dimNames(2) = 'nVertLevels'
      r3Ptr(1) % dimNames(3) = 'nCells'

      r3Ptr(1) % defaultValue = 0.0
      allocate(r3Ptr(1) % attLists(size(r3Ptr(1) % constituentNames, dim=1)))
      do index_counter = 1, size(r3Ptr(1) % constituentNames, dim=1)
         allocate(r3Ptr(1) % attLists(index_counter) % attList)
      end do
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qv', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Water vapor mixing ratio')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'qv'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qc', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Cloud water mixing ratio')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'qc'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qr', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Rain water mixing ratio')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'qr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qi', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Ice mixing ratio')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'qi'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qs', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Snow mixing ratio')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'qs'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qg', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Graupel mixing ratio')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'qg'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_ni', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Cloud ice number concentration')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'nb kg^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'ni'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_nr', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Rain number concentration')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'nb kg^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'nr'
      end if
      r3Ptr(1) % block => block
! Defining time level 2
      allocate( r3Ptr(2) % constituentNames(numConstituents) )
      r3Ptr(2) % fieldName = 'scalars'
      r3Ptr(2) % isDecomposed = .true.
      r3Ptr(2) % hasTimeDimension = .true.
      r3Ptr(2) % isVarArray = .true.
      r3Ptr(2) % isPersistent = .true.
      r3Ptr(2) % isActive = .false.

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qv', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'qv'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qc', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'qc'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qr', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'qr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qi', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'qi'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qs', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'qs'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qg', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'qg'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_ni', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'ni'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_nr', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'nr'
      end if

! Setup dimensions for       
      r3Ptr(2) % dimNames(1) = 'num_scalars'
      r3Ptr(2) % dimNames(2) = 'nVertLevels'
      r3Ptr(2) % dimNames(3) = 'nCells'

      r3Ptr(2) % defaultValue = 0.0
      allocate(r3Ptr(2) % attLists(size(r3Ptr(2) % constituentNames, dim=1)))
      do index_counter = 1, size(r3Ptr(2) % constituentNames, dim=1)
         allocate(r3Ptr(2) % attLists(index_counter) % attList)
      end do
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qv', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Water vapor mixing ratio')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'qv'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qc', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Cloud water mixing ratio')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'qc'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qr', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Rain water mixing ratio')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'qr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qi', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Ice mixing ratio')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'qi'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qs', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Snow mixing ratio')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'qs'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qg', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Graupel mixing ratio')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'qg'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_ni', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Cloud ice number concentration')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'nb kg^{-1}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'ni'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_nr', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Rain number concentration')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'nb kg^{-1}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'nr'
      end if
      r3Ptr(2) % block => block

            r3Ptr(1) % isActive = .true.
            r3Ptr(2) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'scalars', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'scalars', r3Ptr)

! Define var array pv_scalars
      allocate(r3Ptr(2))
      index_counter = 0
      group_counter = -1
      group_start = -1
      group_started = .false.

! Starting group passive
! Define constituent var pv_scalar
! My Packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_pv_scalar', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'passive_start', group_start)
         end if
         group_started = .true.
      end if
      if (.not. group_started) then
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'passive_start', -1)
            call mpas_pool_add_dimension(newSubPool, 'passive_end', -1)
         end if
      else
         group_started = .false.
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'passive_end', index_counter)
         end if
      end if
! End of group       

      numConstituents = index_counter
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'num_pv_scalars', numConstituents)
      end if
! Defining time level 1
      allocate( r3Ptr(1) % constituentNames(numConstituents) )
      r3Ptr(1) % fieldName = 'pv_scalars'
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .true.
      r3Ptr(1) % isVarArray = .true.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_pv_scalar', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'pv_scalar'
      end if

! Setup dimensions for       
      r3Ptr(1) % dimNames(1) = 'num_pv_scalars'
      r3Ptr(1) % dimNames(2) = 'nVertLevels'
      r3Ptr(1) % dimNames(3) = 'nCells'

      r3Ptr(1) % defaultValue = 0.0
      allocate(r3Ptr(1) % attLists(size(r3Ptr(1) % constituentNames, dim=1)))
      do index_counter = 1, size(r3Ptr(1) % constituentNames, dim=1)
         allocate(r3Ptr(1) % attLists(index_counter) % attList)
      end do
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_pv_scalar', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Ertel''s potential vorticity to be advected '&
//'by the model (1 PVU = 10^{-6} m^{2} s^{-1} K kg^{-1})')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'PVU')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'pv_scalar'
      end if
      r3Ptr(1) % block => block
! Defining time level 2
      allocate( r3Ptr(2) % constituentNames(numConstituents) )
      r3Ptr(2) % fieldName = 'pv_scalars'
      r3Ptr(2) % isDecomposed = .true.
      r3Ptr(2) % hasTimeDimension = .true.
      r3Ptr(2) % isVarArray = .true.
      r3Ptr(2) % isPersistent = .true.
      r3Ptr(2) % isActive = .false.

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_pv_scalar', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'pv_scalar'
      end if

! Setup dimensions for       
      r3Ptr(2) % dimNames(1) = 'num_pv_scalars'
      r3Ptr(2) % dimNames(2) = 'nVertLevels'
      r3Ptr(2) % dimNames(3) = 'nCells'

      r3Ptr(2) % defaultValue = 0.0
      allocate(r3Ptr(2) % attLists(size(r3Ptr(2) % constituentNames, dim=1)))
      do index_counter = 1, size(r3Ptr(2) % constituentNames, dim=1)
         allocate(r3Ptr(2) % attLists(index_counter) % attList)
      end do
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_pv_scalar', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Ertel''s potential vorticity to be advected '&
//'by the model (1 PVU = 10^{-6} m^{2} s^{-1} K kg^{-1})')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'PVU')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'pv_scalar'
      end if
      r3Ptr(2) % block => block

            r3Ptr(1) % isActive = .true.
            r3Ptr(2) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'pv_scalars', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'pv_scalars', r3Ptr)

! Define var array aerosols
      allocate(r3Ptr(2))
      index_counter = 0
      group_counter = -1
      group_start = -1
      group_started = .false.

! Starting group aer_cam
! Define constituent var sul
! My Packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_sul', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'aer_cam_start', group_start)
         end if
         group_started = .true.
      end if
! Define constituent var sslt
! My packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_sslt', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'aer_cam_start', group_start)
         end if
         group_started = .true.
      end if
! Define constituent var dust1
! My packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_dust1', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'aer_cam_start', group_start)
         end if
         group_started = .true.
      end if
! Define constituent var dust2
! My packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_dust2', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'aer_cam_start', group_start)
         end if
         group_started = .true.
      end if
! Define constituent var dust3
! My packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_dust3', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'aer_cam_start', group_start)
         end if
         group_started = .true.
      end if
! Define constituent var dust4
! My packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_dust4', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'aer_cam_start', group_start)
         end if
         group_started = .true.
      end if
! Define constituent var ocpho
! My packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_ocpho', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'aer_cam_start', group_start)
         end if
         group_started = .true.
      end if
! Define constituent var bcpho
! My packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_bcpho', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'aer_cam_start', group_start)
         end if
         group_started = .true.
      end if
! Define constituent var ocphi
! My packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_ocphi', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'aer_cam_start', group_start)
         end if
         group_started = .true.
      end if
! Define constituent var bcphi
! My packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_bcphi', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'aer_cam_start', group_start)
         end if
         group_started = .true.
      end if
! Define constituent var bg
! My packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_bg', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'aer_cam_start', group_start)
         end if
         group_started = .true.
      end if
! Define constituent var volc
! My packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_volc', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'aer_cam_start', group_start)
         end if
         group_started = .true.
      end if
      if (.not. group_started) then
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'aer_cam_start', -1)
            call mpas_pool_add_dimension(newSubPool, 'aer_cam_end', -1)
         end if
      else
         group_started = .false.
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'aer_cam_end', index_counter)
         end if
      end if
! End of group       

      numConstituents = index_counter
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'num_aerosols', numConstituents)
      end if
! Defining time level 1
      allocate( r3Ptr(1) % constituentNames(numConstituents) )
      r3Ptr(1) % fieldName = 'aerosols'
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .true.
      r3Ptr(1) % isVarArray = .true.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_sul', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'sul'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_sslt', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'sslt'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_dust1', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'dust1'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_dust2', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'dust2'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_dust3', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'dust3'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_dust4', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'dust4'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_ocpho', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'ocpho'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_bcpho', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'bcpho'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_ocphi', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'ocphi'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_bcphi', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'bcphi'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_bg', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'bg'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_volc', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'volc'
      end if

! Setup dimensions for       
      r3Ptr(1) % dimNames(1) = 'num_aerosols'
      r3Ptr(1) % dimNames(2) = 'nAerLevels'
      r3Ptr(1) % dimNames(3) = 'nCells'

      r3Ptr(1) % defaultValue = 0.0
      allocate(r3Ptr(1) % attLists(size(r3Ptr(1) % constituentNames, dim=1)))
      do index_counter = 1, size(r3Ptr(1) % constituentNames, dim=1)
         allocate(r3Ptr(1) % attLists(index_counter) % attList)
      end do
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_sul', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Sulfate soluble (SUL) aerosol concentration')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'sul'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_sslt', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Sea-salt (SSLT) aerosol concentration')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'sslt'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_dust1', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Dust type 1 (DUST1) aerosol concentration')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'dust1'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_dust2', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Dust type 2 (DUST2) aerosol concentration')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'dust2'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_dust3', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Dust type 3 (DUST3) aerosol concentration')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'dust3'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_dust4', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Dust type 4 (DUST4) aerosol concentration')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'dust4'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_ocpho', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Hydrophobic organic carbon (OCPHO) aerosol '&
//'concentration')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'ocpho'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_bcpho', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Hydrophobic black carbon (BCPHO) aerosol '&
//'concentration')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'bcpho'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_ocphi', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Hydrophilic organic carbon (OCPHI) aerosol '&
//'concentration')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'ocphi'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_bcphi', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Hydrophilic black carbon (BCPHI) aerosol '&
//'concentration')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'bcphi'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_bg', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Background (BG) aerosol concentration')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'bg'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_volc', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Volcanic (VOLC) aerosol concentration')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'volc'
      end if
      r3Ptr(1) % block => block
! Defining time level 2
      allocate( r3Ptr(2) % constituentNames(numConstituents) )
      r3Ptr(2) % fieldName = 'aerosols'
      r3Ptr(2) % isDecomposed = .true.
      r3Ptr(2) % hasTimeDimension = .true.
      r3Ptr(2) % isVarArray = .true.
      r3Ptr(2) % isPersistent = .true.
      r3Ptr(2) % isActive = .false.

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_sul', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'sul'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_sslt', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'sslt'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_dust1', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'dust1'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_dust2', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'dust2'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_dust3', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'dust3'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_dust4', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'dust4'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_ocpho', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'ocpho'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_bcpho', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'bcpho'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_ocphi', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'ocphi'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_bcphi', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'bcphi'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_bg', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'bg'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_volc', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'volc'
      end if

! Setup dimensions for       
      r3Ptr(2) % dimNames(1) = 'num_aerosols'
      r3Ptr(2) % dimNames(2) = 'nAerLevels'
      r3Ptr(2) % dimNames(3) = 'nCells'

      r3Ptr(2) % defaultValue = 0.0
      allocate(r3Ptr(2) % attLists(size(r3Ptr(2) % constituentNames, dim=1)))
      do index_counter = 1, size(r3Ptr(2) % constituentNames, dim=1)
         allocate(r3Ptr(2) % attLists(index_counter) % attList)
      end do
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_sul', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Sulfate soluble (SUL) aerosol concentration')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'sul'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_sslt', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Sea-salt (SSLT) aerosol concentration')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'sslt'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_dust1', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Dust type 1 (DUST1) aerosol concentration')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'dust1'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_dust2', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Dust type 2 (DUST2) aerosol concentration')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'dust2'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_dust3', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Dust type 3 (DUST3) aerosol concentration')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'dust3'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_dust4', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Dust type 4 (DUST4) aerosol concentration')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'dust4'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_ocpho', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Hydrophobic organic carbon (OCPHO) aerosol '&
//'concentration')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'ocpho'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_bcpho', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Hydrophobic black carbon (BCPHO) aerosol '&
//'concentration')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'bcpho'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_ocphi', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Hydrophilic organic carbon (OCPHI) aerosol '&
//'concentration')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'ocphi'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_bcphi', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Hydrophilic black carbon (BCPHI) aerosol '&
//'concentration')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'bcphi'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_bg', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Background (BG) aerosol concentration')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'bg'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_volc', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Volcanic (VOLC) aerosol concentration')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg m^{-3}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'volc'
      end if
      r3Ptr(2) % block => block

            r3Ptr(1) % isActive = .true.
            r3Ptr(2) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'aerosols', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'aerosols', r3Ptr)

! Define variable initial_time
      allocate(c0Ptr(2))

! Setting up time level 1
      c0Ptr(1) % fieldName = 'initial_time'
      c0Ptr(1) % isVarArray = .false.
      c0Ptr(1) % isDecomposed = .false.
      c0Ptr(1) % hasTimeDimension = .false.
      c0Ptr(1) % defaultValue = ''
      c0Ptr(1) % scalar = ''
      allocate(c0Ptr(1) % attLists(1))
      allocate(c0Ptr(1) % attLists(1) % attList)
      call mpas_add_att(c0Ptr(1) % attLists(1) % attList, 'units', 'YYYY-MM-DD_hh:mm:ss')
      call mpas_add_att(c0Ptr(1) % attLists(1) % attList, 'long_name', 'Model initialization time')
      c0Ptr(1) % missingValue = MPAS_CHAR_FILLVAL
      c0Ptr(1) % block => block

! Setting up time level 2
      c0Ptr(2) % fieldName = 'initial_time'
      c0Ptr(2) % isVarArray = .false.
      c0Ptr(2) % isDecomposed = .false.
      c0Ptr(2) % hasTimeDimension = .false.
      c0Ptr(2) % defaultValue = ''
      c0Ptr(2) % scalar = ''
      allocate(c0Ptr(2) % attLists(1))
      allocate(c0Ptr(2) % attLists(1) % attList)
      call mpas_add_att(c0Ptr(2) % attLists(1) % attList, 'units', 'YYYY-MM-DD_hh:mm:ss')
      call mpas_add_att(c0Ptr(2) % attLists(1) % attList, 'long_name', 'Model initialization time')
      c0Ptr(2) % missingValue = MPAS_CHAR_FILLVAL
      c0Ptr(2) % block => block

      c0Ptr(1) % isActive = .true.
      c0Ptr(2) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'initial_time', c0Ptr)
      call mpas_pool_add_field(block % allFields, 'initial_time', c0Ptr)

! Define variable xtime
      allocate(c0Ptr(2))

! Setting up time level 1
      c0Ptr(1) % fieldName = 'xtime'
      c0Ptr(1) % isVarArray = .false.
      c0Ptr(1) % isDecomposed = .false.
      c0Ptr(1) % hasTimeDimension = .true.
      c0Ptr(1) % defaultValue = ''
      c0Ptr(1) % scalar = ''
      allocate(c0Ptr(1) % attLists(1))
      allocate(c0Ptr(1) % attLists(1) % attList)
      call mpas_add_att(c0Ptr(1) % attLists(1) % attList, 'units', 'YYYY-MM-DD_hh:mm:ss')
      call mpas_add_att(c0Ptr(1) % attLists(1) % attList, 'long_name', 'Model valid time')
      c0Ptr(1) % missingValue = MPAS_CHAR_FILLVAL
      c0Ptr(1) % block => block

! Setting up time level 2
      c0Ptr(2) % fieldName = 'xtime'
      c0Ptr(2) % isVarArray = .false.
      c0Ptr(2) % isDecomposed = .false.
      c0Ptr(2) % hasTimeDimension = .true.
      c0Ptr(2) % defaultValue = ''
      c0Ptr(2) % scalar = ''
      allocate(c0Ptr(2) % attLists(1))
      allocate(c0Ptr(2) % attLists(1) % attList)
      call mpas_add_att(c0Ptr(2) % attLists(1) % attList, 'units', 'YYYY-MM-DD_hh:mm:ss')
      call mpas_add_att(c0Ptr(2) % attLists(1) % attList, 'long_name', 'Model valid time')
      c0Ptr(2) % missingValue = MPAS_CHAR_FILLVAL
      c0Ptr(2) % block => block

      c0Ptr(1) % isActive = .true.
      c0Ptr(2) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'xtime', c0Ptr)
      call mpas_pool_add_field(block % allFields, 'xtime', c0Ptr)

! Define variable u
      allocate(r2Ptr(2))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'u'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Horizontal normal velocity at edges')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

! Setting up time level 2
      r2Ptr(2) % fieldName = 'u'
      r2Ptr(2) % isVarArray = .false.
      r2Ptr(2) % isDecomposed = .true.
      r2Ptr(2) % hasTimeDimension = .true.
      r2Ptr(2) % isPersistent = .true.
      r2Ptr(2) % isActive = .false.
! Setting up dimensions
      r2Ptr(2) % dimNames(1) = 'nVertLevels'
      r2Ptr(2) % dimNames(2) = 'nEdges'
      r2Ptr(2) % defaultValue = 0.0
      allocate(r2Ptr(2) % attLists(1))
      allocate(r2Ptr(2) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'long_name', 'Horizontal normal velocity at edges')
      r2Ptr(2) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(2) % block => block

      r2Ptr(1) % isActive = .true.
      r2Ptr(2) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'u', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'u', r2Ptr)

! Define variable w
      allocate(r2Ptr(2))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'w'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Vertical velocity at vertical cell faces')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

! Setting up time level 2
      r2Ptr(2) % fieldName = 'w'
      r2Ptr(2) % isVarArray = .false.
      r2Ptr(2) % isDecomposed = .true.
      r2Ptr(2) % hasTimeDimension = .true.
      r2Ptr(2) % isPersistent = .true.
      r2Ptr(2) % isActive = .false.
! Setting up dimensions
      r2Ptr(2) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(2) % dimNames(2) = 'nCells'
      r2Ptr(2) % defaultValue = 0.0
      allocate(r2Ptr(2) % attLists(1))
      allocate(r2Ptr(2) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'long_name', 'Vertical velocity at vertical cell faces')
      r2Ptr(2) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(2) % block => block

      r2Ptr(1) % isActive = .true.
      r2Ptr(2) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'w', r2Ptr)

! Define variable rho_zz
      allocate(r2Ptr(2))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rho_zz'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Dry air density divided by d(zeta)/dz')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

! Setting up time level 2
      r2Ptr(2) % fieldName = 'rho_zz'
      r2Ptr(2) % isVarArray = .false.
      r2Ptr(2) % isDecomposed = .true.
      r2Ptr(2) % hasTimeDimension = .true.
      r2Ptr(2) % isPersistent = .true.
      r2Ptr(2) % isActive = .false.
! Setting up dimensions
      r2Ptr(2) % dimNames(1) = 'nVertLevels'
      r2Ptr(2) % dimNames(2) = 'nCells'
      r2Ptr(2) % defaultValue = 0.0
      allocate(r2Ptr(2) % attLists(1))
      allocate(r2Ptr(2) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'units', 'kg m^{-3}')
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'long_name', 'Dry air density divided by d(zeta)/dz')
      r2Ptr(2) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(2) % block => block

      r2Ptr(1) % isActive = .true.
      r2Ptr(2) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rho_zz', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_zz', r2Ptr)

! Define variable theta_m
      allocate(r2Ptr(2))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'theta_m'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Moist potential temperature: theta*(1+q_v*R_v/R_d)')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

! Setting up time level 2
      r2Ptr(2) % fieldName = 'theta_m'
      r2Ptr(2) % isVarArray = .false.
      r2Ptr(2) % isDecomposed = .true.
      r2Ptr(2) % hasTimeDimension = .true.
      r2Ptr(2) % isPersistent = .true.
      r2Ptr(2) % isActive = .false.
! Setting up dimensions
      r2Ptr(2) % dimNames(1) = 'nVertLevels'
      r2Ptr(2) % dimNames(2) = 'nCells'
      r2Ptr(2) % defaultValue = 0.0
      allocate(r2Ptr(2) % attLists(1))
      allocate(r2Ptr(2) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'long_name', 'Moist potential temperature: theta*(1+q_v*R_v/R_d)')
      r2Ptr(2) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(2) % block => block

      r2Ptr(1) % isActive = .true.
      r2Ptr(2) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'theta_m', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'theta_m', r2Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
         call mpas_pool_add_config(newSubPool, 'is_periodic', block % domain % is_periodic)
         call mpas_pool_add_config(newSubPool, 'x_period', block % domain % x_period)
         call mpas_pool_add_config(newSubPool, 'y_period', block % domain % y_period)
      end if

   end subroutine atm_generate_pool_state


   subroutine atm_generate_pool_diag(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_io, only : MPAS_REAL_FILLVAL, MPAS_INT_FILLVAL, MPAS_CHAR_FILLVAL
      implicit none
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), dimension(:), pointer :: r0Ptr
      type (field1DReal), dimension(:), pointer :: r1Ptr
      type (field2DReal), dimension(:), pointer :: r2Ptr
      type (field3DReal), dimension(:), pointer :: r3Ptr
      type (field4DReal), dimension(:), pointer :: r4Ptr
      type (field5DReal), dimension(:), pointer :: r5Ptr
      type (field0DInteger), dimension(:), pointer :: i0Ptr
      type (field1DInteger), dimension(:), pointer :: i1Ptr
      type (field2DInteger), dimension(:), pointer :: i2Ptr
      type (field3DInteger), dimension(:), pointer :: i3Ptr
      type (field0DChar), dimension(:), pointer :: c0Ptr
      type (field1DChar), dimension(:), pointer :: c1Ptr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: mp_kessler_inActive
      logical, pointer :: mp_thompson_inActive
      logical, pointer :: mp_wsm6_inActive
      logical, pointer :: cu_grell_freitas_inActive
      logical, pointer :: cu_kain_fritsch_inActive
      logical, pointer :: cu_tiedtke_inActive
      logical, pointer :: bl_ysu_inActive
      logical, pointer :: bl_mynn_inActive
      logical, pointer :: iauActive
      logical, pointer :: limited_areaActive
      logical, pointer :: jedi_daActive
      logical, pointer :: tendenciesActive


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'mp_kessler_inActive', mp_kessler_inActive)
      call mpas_pool_get_package(packagePool, 'mp_thompson_inActive', mp_thompson_inActive)
      call mpas_pool_get_package(packagePool, 'mp_wsm6_inActive', mp_wsm6_inActive)
      call mpas_pool_get_package(packagePool, 'cu_grell_freitas_inActive', cu_grell_freitas_inActive)
      call mpas_pool_get_package(packagePool, 'cu_kain_fritsch_inActive', cu_kain_fritsch_inActive)
      call mpas_pool_get_package(packagePool, 'cu_tiedtke_inActive', cu_tiedtke_inActive)
      call mpas_pool_get_package(packagePool, 'bl_ysu_inActive', bl_ysu_inActive)
      call mpas_pool_get_package(packagePool, 'bl_mynn_inActive', bl_mynn_inActive)
      call mpas_pool_get_package(packagePool, 'iauActive', iauActive)
      call mpas_pool_get_package(packagePool, 'limited_areaActive', limited_areaActive)
      call mpas_pool_get_package(packagePool, 'jedi_daActive', jedi_daActive)
      call mpas_pool_get_package(packagePool, 'tendenciesActive', tendenciesActive)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'diag', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'diag', newSubPool)

! Define variable qvapor
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'qvapor'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Water vapor mixing ratio')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'qvapor', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'qvapor', r2Ptr)

! Define variable cofrz
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'cofrz'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .false.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertLevels'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 's m^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'coefficient for implicit contribution of Omega to density '&
//'update')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cofrz', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'cofrz', r1Ptr)

! Define variable cofwr
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'cofwr'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'coefficient for implicit contribution of density to the '&
//'vertical velocity update')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cofwr', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'cofwr', r2Ptr)

! Define variable cofwz
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'cofwz'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1} K^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'coefficient for implicit contribution of density to the '&
//'vertical velocity update')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cofwz', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'cofwz', r2Ptr)

! Define variable coftz
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'coftz'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 's K')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'coefficient for implicit contribution of omega vertical '&
//'derivative to the theta_m update')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'coftz', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'coftz', r2Ptr)

! Define variable cofwt
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'cofwt'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1} K^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'coefficient for implicit contribution of density to the '&
//'vertical velocity update')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cofwt', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'cofwt', r2Ptr)

! Define variable a_tri
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'a_tri'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'implicit tridiagonal solve coefficients')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'a_tri', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'a_tri', r2Ptr)

! Define variable alpha_tri
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'alpha_tri'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'implicit tridiagonal solve coefficients')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'alpha_tri', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'alpha_tri', r2Ptr)

! Define variable gamma_tri
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'gamma_tri'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'implicit tridiagonal solve coefficients')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'gamma_tri', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'gamma_tri', r2Ptr)

! Define variable pressure_p
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'pressure_p'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'Pa')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Perturbation pressure')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'pressure_p', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'pressure_p', r2Ptr)

! Define variable rho
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rho'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Dry air density')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rho', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho', r2Ptr)

! Define variable theta
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'theta'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Potential temperature')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'theta', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'theta', r2Ptr)

! Define variable temperature
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'temperature'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'temperature')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (jedi_daActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'temperature', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'temperature', r2Ptr)

! Define variable relhum
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'relhum'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'percent')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Relative humidity')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'relhum', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'relhum', r2Ptr)

! Define variable spechum
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'spechum'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Specific humidity')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (jedi_daActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'spechum', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'spechum', r2Ptr)

! Define variable v
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'v'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Horizontal tangential velocity at edges')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'v', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'v', r2Ptr)

! Define variable divergence
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'divergence'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 's^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Horizontal velocity divergence at cell center')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'divergence', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'divergence', r2Ptr)

! Define variable vorticity
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'vorticity'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nVertices'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 's^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Relative vorticity at vertices')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'vorticity', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'vorticity', r2Ptr)

! Define variable pv_edge
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'pv_edge'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 's^{-1} kg^{-1} m^3')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'absolute vorticity/rho_zz averaged to the cell edge from '&
//'the vertices')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'pv_edge', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'pv_edge', r2Ptr)

! Define variable rho_edge
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rho_edge'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'rho_zz averaged from cell centers to the cell edge')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rho_edge', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_edge', r2Ptr)

! Define variable ke
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ke'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm^2 s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Kinetic energy at a cell center')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ke', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ke', r2Ptr)

! Define variable pv_vertex
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'pv_vertex'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nVertices'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 's^{-1} kg^{-1} m^3')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'absolute vorticity/rho_zz at a vertex')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'pv_vertex', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'pv_vertex', r2Ptr)

! Define variable pv_cell
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'pv_cell'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 's^{-1} kg^{-1} m^3')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'absolute vorticity/rho_zz averaged to the cell center from '&
//'the vertices')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'pv_cell', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'pv_cell', r2Ptr)

! Define variable uReconstructX
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uReconstructX'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Cartesian x-component of reconstructed horizontal velocity '&
//'at cell centers')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uReconstructX', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructX', r2Ptr)

! Define variable uReconstructY
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uReconstructY'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Cartesian y-component of reconstructed horizontal velocity '&
//'at cell centers')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uReconstructY', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructY', r2Ptr)

! Define variable uReconstructZ
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uReconstructZ'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Cartesian z-component of reconstructed horizontal velocity '&
//'at cell centers')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uReconstructZ', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructZ', r2Ptr)

! Define variable uReconstructZonal
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uReconstructZonal'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Zonal component of reconstructed horizontal velocity at '&
//'cell centers')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uReconstructZonal', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructZonal', r2Ptr)

! Define variable uReconstructMeridional
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uReconstructMeridional'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Meridional component of reconstructed horizontal velocity '&
//'at cell centers')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uReconstructMeridional', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructMeridional', r2Ptr)

! Define variable rv
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rv'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', '???')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', '???')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rv', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rv', r2Ptr)

! Define variable circulation
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'circulation'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nVertices'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm^2 s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Horizontal circulation at vertices')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'circulation', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'circulation', r2Ptr)

! Define variable gradPVt
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'gradPVt'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', '???')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Gradient in the tangential direction of PV at edge '&
//'locations')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'gradPVt', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'gradPVt', r2Ptr)

! Define variable gradPVn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'gradPVn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', '???')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Gradient in the normal direction of PV at edge locations')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'gradPVn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'gradPVn', r2Ptr)

! Define variable h_divergence
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'h_divergence'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', '???')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', '???')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'h_divergence', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'h_divergence', r2Ptr)

! Define variable exner
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'exner'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Exner function')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'exner', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'exner', r2Ptr)

! Define variable exner_base
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'exner_base'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Base-state Exner function')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'exner_base', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'exner_base', r2Ptr)

! Define variable rtheta_base
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rtheta_base'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg K m^{-3}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'reference state rho*theta/zz')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rtheta_base', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rtheta_base', r2Ptr)

! Define variable pressure
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'pressure'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'Pa')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Pressure')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'pressure', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'pressure', r2Ptr)

! Define variable pressure_base
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'pressure_base'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'Pa')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Base state pressure')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'pressure_base', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'pressure_base', r2Ptr)

! Define variable rho_base
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rho_base'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Base state dry air density')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rho_base', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_base', r2Ptr)

! Define variable theta_base
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'theta_base'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Base state potential temperature')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'theta_base', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'theta_base', r2Ptr)

! Define variable rho_zz_old_split
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rho_zz_old_split'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'rho/zz')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rho_zz_old_split', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_zz_old_split', r2Ptr)

! Define variable ruAvg
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ruAvg'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'time-averaged rho*u/zz used in scalar transport')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ruAvg', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ruAvg', r2Ptr)

! Define variable wwAvg
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'wwAvg'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'time-averaged rho*omega/zz used in scalar transport')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'wwAvg', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'wwAvg', r2Ptr)

! Define variable ruAvg_split
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ruAvg_split'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'time-averaged rho*u/zz used in scalar transport')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ruAvg_split', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ruAvg_split', r2Ptr)

! Define variable wwAvg_split
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'wwAvg_split'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'time-averaged rho*omega/zz used in scalar transport')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'wwAvg_split', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'wwAvg_split', r2Ptr)

! Define variable cqu
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'cqu'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'rho_d/rho_m at cell edge (u points)')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cqu', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'cqu', r2Ptr)

! Define variable cqw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'cqw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'rho_d/rho_m at w points')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cqw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'cqw', r2Ptr)

! Define variable ru
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ru'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'horizontal momentum at cell edge (rho*u/zz)')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ru', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ru', r2Ptr)

! Define variable ru_p
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ru_p'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'acoustic perturbation horizontal momentum at cell edge  '&
//'(rho*u/zz)')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ru_p', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ru_p', r2Ptr)

! Define variable ru_save
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ru_save'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'predicted value of horizontal momentum, saved before '&
//'acoustic steps')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ru_save', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ru_save', r2Ptr)

! Define variable rw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'rho*omega/zz carried at w points')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rw', r2Ptr)

! Define variable rw_p
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rw_p'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'acoustic perturbation rho*omega/zz carried at w points')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rw_p', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rw_p', r2Ptr)

! Define variable rw_save
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rw_save'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'predicted value of rho*omega/zz, saved before acoustic '&
//'steps')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rw_save', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rw_save', r2Ptr)

! Define variable rtheta_p
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rtheta_p'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg K m^{-3}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'rho*theta_m/zz perturbation from the reference state '&
//'value')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rtheta_p', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rtheta_p', r2Ptr)

! Define variable rtheta_pp
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rtheta_pp'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg K m^{-3}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'rho*theta_m/zz perturbation from rtheta_p')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rtheta_pp', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rtheta_pp', r2Ptr)

! Define variable rtheta_p_save
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rtheta_p_save'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg K m^{-3}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'predicted value rtheta_p, saved before acoustic steps')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rtheta_p_save', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rtheta_p_save', r2Ptr)

! Define variable rtheta_pp_old
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rtheta_pp_old'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg K m^{-3}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'old time level values of rho*theta_m/zz perturbation from '&
//'rtheta_p, used in 3D divergence damping')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rtheta_pp_old', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rtheta_pp_old', r2Ptr)

! Define variable rho_p
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rho_p'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'rho/zz perturbation from the reference state value, '&
//'advanced over acoustic steps')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rho_p', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_p', r2Ptr)

! Define variable rho_pp
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rho_pp'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'rho/zz perturbation from rho_pp, advanced over acoustic '&
//'steps')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rho_pp', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_pp', r2Ptr)

! Define variable rho_p_save
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rho_p_save'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'predicted value rho_p, saved before acoustic steps')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rho_p_save', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_p_save', r2Ptr)

! Define variable kdiff
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'kdiff'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm^2 s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Smagorinsky horizontal eddy viscosity')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'kdiff', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'kdiff', r2Ptr)

! Define variable tend_rtheta_adv
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_rtheta_adv'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg K m^{-3} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'flux divergence for rho*theta_m/zz, used in the Tiedtke '&
//'convective parameterization')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'tend_rtheta_adv', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_rtheta_adv', r2Ptr)

! Define variable surface_pressure
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'surface_pressure'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'Pa')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Diagnosed surface pressure')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'surface_pressure', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'surface_pressure', r1Ptr)

! Define variable iso_levels
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'iso_levels'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .false.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nIsoLevels'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'Pa')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Levels for vertical interpolation to isobaric surfaces')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'iso_levels', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'iso_levels', r1Ptr)

! Define variable mslp
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'mslp'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'Pa')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Mean sea-level pressure')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'mslp', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'mslp', r1Ptr)

! Define variable temperature_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'temperature_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Temperature vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'temperature_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'temperature_isobaric', r2Ptr)

! Define variable dewpoint_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dewpoint_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Dewpoint temperature vertically interpolated to isobaric '&
//'levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dewpoint_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dewpoint_isobaric', r2Ptr)

! Define variable relhum_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'relhum_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'percent')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Relative humidity vertically interpolated to isobaric '&
//'levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'relhum_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'relhum_isobaric', r2Ptr)

! Define variable qvapor_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'qvapor_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Water vapor mixing ratio vertically interpolated to '&
//'isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'qvapor_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'qvapor_isobaric', r2Ptr)

! Define variable uzonal_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uzonal_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Reconstructed zonal wind at cell centers vertically '&
//'interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uzonal_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uzonal_isobaric', r2Ptr)

! Define variable umeridional_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'umeridional_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Reconstructed meridional wind at cell centers vertically '&
//'interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'umeridional_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'umeridional_isobaric', r2Ptr)

! Define variable height_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'height_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Geometric height interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'height_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'height_isobaric', r2Ptr)

! Define variable geoheight_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'geoheight_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Geopotential height interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'geoheight_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'geoheight_isobaric', r2Ptr)

! Define variable w_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'w_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Vertical velocity vertically interpolated to isobaric '&
//'levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'w_isobaric', r2Ptr)

! Define variable vorticity_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'vorticity_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 's^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Vertical vorticity at cell centers vertically interpolated '&
//'to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'vorticity_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'vorticity_isobaric', r2Ptr)

! Define variable ertel_pv_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ertel_pv_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Ertel''s potential vorticity vertically interpolated to '&
//'isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ertel_pv_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ertel_pv_isobaric', r2Ptr)

! Define variable depv_dt_lw_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_lw_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from longwave radiation vertically '&
//'interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_lw_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_lw_isobaric', r2Ptr)

! Define variable depv_dt_sw_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_sw_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from shortwave radiation vertically '&
//'interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_sw_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_sw_isobaric', r2Ptr)

! Define variable depv_dt_bl_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_bl_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from PBL vertically interpolated to '&
//'isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_bl_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_bl_isobaric', r2Ptr)

! Define variable depv_dt_cu_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_cu_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from cumulus scheme vertically '&
//'interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_cu_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_cu_isobaric', r2Ptr)

! Define variable depv_dt_mp_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_mp_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from microphysics vertically '&
//'interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_mp_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_mp_isobaric', r2Ptr)

! Define variable depv_dt_mix_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_mix_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from explicit numerical mixing '&
//'vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_mix_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_mix_isobaric', r2Ptr)

! Define variable depv_dt_fric_mix_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_fric_mix_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Frictional PV tendency from explicit horizontal mixing '&
//'vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_fric_mix_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_fric_mix_isobaric', r2Ptr)

! Define variable depv_dt_fric_bl_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_fric_bl_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Frictional PV tendency from PBL scheme vertically '&
//'interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_fric_bl_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_fric_bl_isobaric', r2Ptr)

! Define variable depv_dt_fric_cu_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_fric_cu_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Frictional PV tendency from cumulus scheme vertically '&
//'interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_fric_cu_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_fric_cu_isobaric', r2Ptr)

! Define variable depv_dt_dyn_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_dyn_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'PV tendency from dynamics vertically interpolated to '&
//'isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_dyn_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_dyn_isobaric', r2Ptr)

! Define variable depv_dt_diab_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_diab_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Sum of calculated PV tendencies from diabatic processes '&
//'vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_diab_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_diab_isobaric', r2Ptr)

! Define variable depv_dt_fric_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_fric_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Sum of calculated PV tendencies from frictional processes '&
//'vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_fric_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_fric_isobaric', r2Ptr)

! Define variable depv_dt_mp_evap_cw_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_mp_evap_cw_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from microphysics - net '&
//'condensation/evaporation of cloud water - vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_mp_evap_cw_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_mp_evap_cw_isobaric', r2Ptr)

! Define variable depv_dt_mp_evap_rw_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_mp_evap_rw_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from microphysics - evaporation of '&
//'rain water - vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_mp_evap_rw_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_mp_evap_rw_isobaric', r2Ptr)

! Define variable depv_dt_mp_depo_ice_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_mp_depo_ice_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from microphysics - net '&
//'deposition/sublimation of all ice hydrometeors - vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_mp_depo_ice_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_mp_depo_ice_isobaric', r2Ptr)

! Define variable depv_dt_mp_melt_ice_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_mp_melt_ice_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from microphysics - melting of all '&
//'ice hydrometeors - vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_mp_melt_ice_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_mp_melt_ice_isobaric', r2Ptr)

! Define variable depv_dt_mp_frez_ice_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_mp_frez_ice_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from microphysics - freezing of all '&
//'ice hydrometeors - vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_mp_frez_ice_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_mp_frez_ice_isobaric', r2Ptr)

! Define variable depv_dt_mp_allproc_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_mp_allproc_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from summed microphysics processes '&
//'vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_mp_allproc_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_mp_allproc_isobaric', r2Ptr)

! Define variable dtheta_dt_mp_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dtheta_dt_mp_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Potential temperature heating rate from microphysics '&
//'vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dtheta_dt_mp_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dtheta_dt_mp_isobaric', r2Ptr)

! Define variable dtheta_dt_dyn_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dtheta_dt_dyn_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Potential temperature tendency from dynamics vertically '&
//'interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dtheta_dt_dyn_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dtheta_dt_dyn_isobaric', r2Ptr)

! Define variable dtheta_dt_cu_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dtheta_dt_cu_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Potential temperature tendency from cumulus scheme '&
//'vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dtheta_dt_cu_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dtheta_dt_cu_isobaric', r2Ptr)

! Define variable acc_depv_dt_lw_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_lw_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from longwave radiation '&
//'vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_lw_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_lw_isobaric', r2Ptr)

! Define variable acc_depv_dt_sw_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_sw_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from shortwave radiation '&
//'vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_sw_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_sw_isobaric', r2Ptr)

! Define variable acc_depv_dt_bl_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_bl_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from PBL vertically '&
//'interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_bl_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_bl_isobaric', r2Ptr)

! Define variable acc_depv_dt_cu_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_cu_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from cumulus scheme '&
//'vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_cu_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_cu_isobaric', r2Ptr)

! Define variable acc_depv_dt_mp_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_mp_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from microphysics '&
//'vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_mp_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_mp_isobaric', r2Ptr)

! Define variable acc_depv_dt_mix_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_mix_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from explicit horizontal '&
//'mixing vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_mix_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_mix_isobaric', r2Ptr)

! Define variable acc_depv_dt_fric_bl_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_fric_bl_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated frictional PV tendency from PBL scheme '&
//'vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_fric_bl_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_fric_bl_isobaric', r2Ptr)

! Define variable acc_depv_dt_fric_cu_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_fric_cu_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated frictional PV tendency from cumulus scheme '&
//'vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_fric_cu_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_fric_cu_isobaric', r2Ptr)

! Define variable acc_depv_dt_fric_mix_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_fric_mix_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated frictional PV tendency from explicit mixing '&
//'vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_fric_mix_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_fric_mix_isobaric', r2Ptr)

! Define variable acc_depv_dt_diab_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_diab_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated sum of calculated PV tendencies from diabatic '&
//'processes vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_diab_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_diab_isobaric', r2Ptr)

! Define variable acc_depv_dt_fric_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_fric_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated sum of calculated PV tendencies from '&
//'frictional processes vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_fric_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_fric_isobaric', r2Ptr)

! Define variable acc_depv_dt_dyn_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_dyn_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated PV tendency from dynamics vertically '&
//'interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_dyn_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_dyn_isobaric', r2Ptr)

! Define variable acc_depv_dt_mp_evap_cw_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_mp_evap_cw_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from microphysics - net '&
//'condensation/evaporation of cloud water - vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_mp_evap_cw_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_mp_evap_cw_isobaric', r2Ptr)

! Define variable acc_depv_dt_mp_evap_rw_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_mp_evap_rw_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from microphysics - '&
//'evaporation of rain water - vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_mp_evap_rw_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_mp_evap_rw_isobaric', r2Ptr)

! Define variable acc_depv_dt_mp_depo_ice_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_mp_depo_ice_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from microphysics - net '&
//'deposition/sublimation of all ice hydrometeors - vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_mp_depo_ice_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_mp_depo_ice_isobaric', r2Ptr)

! Define variable acc_depv_dt_mp_melt_ice_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_mp_melt_ice_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from microphysics - '&
//'melting of all ice hydrometeors - vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_mp_melt_ice_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_mp_melt_ice_isobaric', r2Ptr)

! Define variable acc_depv_dt_mp_frez_ice_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_mp_frez_ice_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from microphysics - '&
//'freezing of all ice hydrometeors - vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_mp_frez_ice_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_mp_frez_ice_isobaric', r2Ptr)

! Define variable acc_depv_dt_mp_allproc_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_mp_allproc_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from summed microphysics '&
//'processes vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_mp_allproc_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_mp_allproc_isobaric', r2Ptr)

! Define variable acc_dtheta_dt_mp_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_dtheta_dt_mp_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated potential temperature heating rate from '&
//'microphysics vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_dtheta_dt_mp_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_dtheta_dt_mp_isobaric', r2Ptr)

! Define variable acc_dtheta_dt_dyn_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_dtheta_dt_dyn_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated potential temperature tendency from dynamics '&
//'vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_dtheta_dt_dyn_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_dtheta_dt_dyn_isobaric', r2Ptr)

! Define variable acc_dtheta_dt_cu_isobaric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_dtheta_dt_cu_isobaric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nIsoLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated potential temperature tendency from cumulus '&
//'scheme vertically interpolated to isobaric levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_dtheta_dt_cu_isobaric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_dtheta_dt_cu_isobaric', r2Ptr)

! Define variable temperature_200hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'temperature_200hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Temperature vertically interpolated to 200 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'temperature_200hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'temperature_200hPa', r1Ptr)

! Define variable temperature_250hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'temperature_250hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Temperature vertically interpolated to 250 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'temperature_250hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'temperature_250hPa', r1Ptr)

! Define variable temperature_300hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'temperature_300hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Temperature vertically interpolated to 300 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'temperature_300hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'temperature_300hPa', r1Ptr)

! Define variable temperature_500hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'temperature_500hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Temperature vertically interpolated to 500 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'temperature_500hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'temperature_500hPa', r1Ptr)

! Define variable temperature_700hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'temperature_700hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Temperature vertically interpolated to 700 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'temperature_700hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'temperature_700hPa', r1Ptr)

! Define variable temperature_850hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'temperature_850hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Temperature vertically interpolated to 850 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'temperature_850hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'temperature_850hPa', r1Ptr)

! Define variable temperature_925hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'temperature_925hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Temperature vertically interpolated to 925 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'temperature_925hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'temperature_925hPa', r1Ptr)

! Define variable relhum_200hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'relhum_200hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'percent')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Relative humidity vertically interpolated to 200 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'relhum_200hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'relhum_200hPa', r1Ptr)

! Define variable relhum_250hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'relhum_250hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'percent')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Relative humidity vertically interpolated to 250 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'relhum_250hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'relhum_250hPa', r1Ptr)

! Define variable relhum_300hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'relhum_300hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'percent')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Relative humidity vertically interpolated to 300 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'relhum_300hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'relhum_300hPa', r1Ptr)

! Define variable relhum_500hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'relhum_500hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'percent')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Relative humidity vertically interpolated to 500 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'relhum_500hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'relhum_500hPa', r1Ptr)

! Define variable relhum_700hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'relhum_700hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'percent')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Relative humidity vertically interpolated to 700 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'relhum_700hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'relhum_700hPa', r1Ptr)

! Define variable relhum_850hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'relhum_850hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'percent')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Relative humidity vertically interpolated to 850 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'relhum_850hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'relhum_850hPa', r1Ptr)

! Define variable relhum_925hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'relhum_925hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'percent')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Relative humidity vertically interpolated to 925 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'relhum_925hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'relhum_925hPa', r1Ptr)

! Define variable dewpoint_200hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'dewpoint_200hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Dewpoint temperature vertically interpolated to 200 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dewpoint_200hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'dewpoint_200hPa', r1Ptr)

! Define variable dewpoint_250hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'dewpoint_250hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Dewpoint temperature vertically interpolated to 250 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dewpoint_250hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'dewpoint_250hPa', r1Ptr)

! Define variable dewpoint_300hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'dewpoint_300hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Dewpoint temperature vertically interpolated to 300 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dewpoint_300hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'dewpoint_300hPa', r1Ptr)

! Define variable dewpoint_500hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'dewpoint_500hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Dewpoint temperature vertically interpolated to 500 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dewpoint_500hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'dewpoint_500hPa', r1Ptr)

! Define variable dewpoint_700hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'dewpoint_700hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Dewpoint temperature vertically interpolated to 700 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dewpoint_700hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'dewpoint_700hPa', r1Ptr)

! Define variable dewpoint_850hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'dewpoint_850hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Dewpoint temperature vertically interpolated to 850 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dewpoint_850hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'dewpoint_850hPa', r1Ptr)

! Define variable dewpoint_925hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'dewpoint_925hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Dewpoint temperature vertically interpolated to 925 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dewpoint_925hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'dewpoint_925hPa', r1Ptr)

! Define variable height_200hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'height_200hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Geometric height interpolated to 200 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'height_200hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'height_200hPa', r1Ptr)

! Define variable height_250hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'height_250hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Geometric height interpolated to 250 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'height_250hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'height_250hPa', r1Ptr)

! Define variable height_300hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'height_300hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Geometric height interpolated to 300 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'height_300hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'height_300hPa', r1Ptr)

! Define variable height_500hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'height_500hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Geometric height interpolated to 500 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'height_500hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'height_500hPa', r1Ptr)

! Define variable height_700hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'height_700hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Geometric height interpolated to 700 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'height_700hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'height_700hPa', r1Ptr)

! Define variable height_850hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'height_850hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Geometric height interpolated to 850 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'height_850hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'height_850hPa', r1Ptr)

! Define variable height_925hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'height_925hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Geometric height interpolated to 925 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'height_925hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'height_925hPa', r1Ptr)

! Define variable uzonal_200hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'uzonal_200hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Reconstructed zonal wind at cell centers, vertically '&
//'interpolated to 200 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uzonal_200hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'uzonal_200hPa', r1Ptr)

! Define variable uzonal_250hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'uzonal_250hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Reconstructed zonal wind at cell centers, vertically '&
//'interpolated to 250 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uzonal_250hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'uzonal_250hPa', r1Ptr)

! Define variable uzonal_300hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'uzonal_300hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Reconstructed zonal wind at cell centers, vertically '&
//'interpolated to 300 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uzonal_300hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'uzonal_300hPa', r1Ptr)

! Define variable uzonal_500hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'uzonal_500hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Reconstructed zonal wind at cell centers, vertically '&
//'interpolated to 500 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uzonal_500hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'uzonal_500hPa', r1Ptr)

! Define variable uzonal_700hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'uzonal_700hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Reconstructed zonal wind at cell centers, vertically '&
//'interpolated to 700 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uzonal_700hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'uzonal_700hPa', r1Ptr)

! Define variable uzonal_850hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'uzonal_850hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Reconstructed zonal wind at cell centers, vertically '&
//'interpolated to 850 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uzonal_850hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'uzonal_850hPa', r1Ptr)

! Define variable uzonal_925hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'uzonal_925hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Reconstructed zonal wind at cell centers, vertically '&
//'interpolated to 925 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uzonal_925hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'uzonal_925hPa', r1Ptr)

! Define variable umeridional_200hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'umeridional_200hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Reconstructed meridional wind at cell centers, vertically '&
//'interpolated to 200 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'umeridional_200hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'umeridional_200hPa', r1Ptr)

! Define variable umeridional_250hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'umeridional_250hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Reconstructed meridional wind at cell centers, vertically '&
//'interpolated to 250 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'umeridional_250hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'umeridional_250hPa', r1Ptr)

! Define variable umeridional_300hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'umeridional_300hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Reconstructed meridional wind at cell centers, vertically '&
//'interpolated to 300 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'umeridional_300hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'umeridional_300hPa', r1Ptr)

! Define variable umeridional_500hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'umeridional_500hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Reconstructed meridional wind at cell centers, vertically '&
//'interpolated to 500 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'umeridional_500hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'umeridional_500hPa', r1Ptr)

! Define variable umeridional_700hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'umeridional_700hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Reconstructed meridional wind at cell centers, vertically '&
//'interpolated to 700 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'umeridional_700hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'umeridional_700hPa', r1Ptr)

! Define variable umeridional_850hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'umeridional_850hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Reconstructed meridional wind at cell centers, vertically '&
//'interpolated to 850 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'umeridional_850hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'umeridional_850hPa', r1Ptr)

! Define variable umeridional_925hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'umeridional_925hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Reconstructed meridional wind at cell centers, vertically '&
//'interpolated to 925 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'umeridional_925hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'umeridional_925hPa', r1Ptr)

! Define variable w_200hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'w_200hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Vertical velocity vertically interpolated to 200 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w_200hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'w_200hPa', r1Ptr)

! Define variable w_250hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'w_250hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Vertical velocity vertically interpolated to 250 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w_250hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'w_250hPa', r1Ptr)

! Define variable w_300hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'w_300hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Vertical velocity vertically interpolated to 300 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w_300hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'w_300hPa', r1Ptr)

! Define variable w_500hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'w_500hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Vertical velocity vertically interpolated to 500 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w_500hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'w_500hPa', r1Ptr)

! Define variable w_700hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'w_700hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Vertical velocity vertically interpolated to 700 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w_700hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'w_700hPa', r1Ptr)

! Define variable w_850hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'w_850hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Vertical velocity vertically interpolated to 850 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w_850hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'w_850hPa', r1Ptr)

! Define variable w_925hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'w_925hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Vertical velocity vertically interpolated to 925 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w_925hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'w_925hPa', r1Ptr)

! Define variable vorticity_200hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'vorticity_200hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 's^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Relative vorticity vertically interpolated to 200 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'vorticity_200hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'vorticity_200hPa', r1Ptr)

! Define variable vorticity_250hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'vorticity_250hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 's^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Relative vorticity vertically interpolated to 250 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'vorticity_250hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'vorticity_250hPa', r1Ptr)

! Define variable vorticity_300hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'vorticity_300hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 's^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Relative vorticity vertically interpolated to 300 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'vorticity_300hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'vorticity_300hPa', r1Ptr)

! Define variable vorticity_500hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'vorticity_500hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 's^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Relative vorticity vertically interpolated to 500 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'vorticity_500hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'vorticity_500hPa', r1Ptr)

! Define variable vorticity_700hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'vorticity_700hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 's^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Relative vorticity vertically interpolated to 700 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'vorticity_700hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'vorticity_700hPa', r1Ptr)

! Define variable vorticity_850hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'vorticity_850hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 's^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Relative vorticity vertically interpolated to 850 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'vorticity_850hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'vorticity_850hPa', r1Ptr)

! Define variable vorticity_925hPa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'vorticity_925hPa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nVertices'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 's^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Relative vorticity vertically interpolated to 925 hPa')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'vorticity_925hPa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'vorticity_925hPa', r1Ptr)

! Define variable cape
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'cape'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'J kg^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Convective available potential energy')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cape', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'cape', r1Ptr)

! Define variable cin
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'cin'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'J kg^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Convective inhibition')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cin', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'cin', r1Ptr)

! Define variable lcl
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'lcl'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Lifted condensation level')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lcl', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lcl', r1Ptr)

! Define variable lfc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'lfc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Level of free convection')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lfc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lfc', r1Ptr)

! Define variable srh_0_1km
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'srh_0_1km'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm^2 s^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Storm relative helicity, 0-1 km AGL')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'srh_0_1km', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'srh_0_1km', r1Ptr)

! Define variable srh_0_3km
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'srh_0_3km'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm^2 s^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Storm relative helicity, 0-3 km AGL')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'srh_0_3km', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'srh_0_3km', r1Ptr)

! Define variable uzonal_surface
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'uzonal_surface'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Zonal wind component at midpoint of lowest model layer')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uzonal_surface', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'uzonal_surface', r1Ptr)

! Define variable uzonal_1km
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'uzonal_1km'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Zonal wind component at 1 km AGL')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uzonal_1km', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'uzonal_1km', r1Ptr)

! Define variable uzonal_3km
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'uzonal_3km'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Zonal wind component at 3 km AGL')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uzonal_3km', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'uzonal_3km', r1Ptr)

! Define variable uzonal_6km
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'uzonal_6km'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Zonal wind component at 6 km AGL')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uzonal_6km', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'uzonal_6km', r1Ptr)

! Define variable umeridional_surface
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'umeridional_surface'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Meridional wind component at midpoint of lowest model '&
//'layer')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'umeridional_surface', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'umeridional_surface', r1Ptr)

! Define variable umeridional_1km
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'umeridional_1km'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Meridional wind component at 1 km AGL')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'umeridional_1km', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'umeridional_1km', r1Ptr)

! Define variable umeridional_3km
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'umeridional_3km'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Meridional wind component at 3 km AGL')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'umeridional_3km', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'umeridional_3km', r1Ptr)

! Define variable umeridional_6km
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'umeridional_6km'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Meridional wind component at 6 km AGL')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'umeridional_6km', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'umeridional_6km', r1Ptr)

! Define variable temperature_surface
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'temperature_surface'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Temperature at midpoint of lowest model layer')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'temperature_surface', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'temperature_surface', r1Ptr)

! Define variable dewpoint_surface
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'dewpoint_surface'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Dewpoint temperature at midpoint of lowest model layer')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dewpoint_surface', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'dewpoint_surface', r1Ptr)

! Define variable updraft_helicity_max
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'updraft_helicity_max'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm^2 s^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Maximum updraft helicity since last output')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'updraft_helicity_max', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'updraft_helicity_max', r1Ptr)

! Define variable w_velocity_max
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'w_velocity_max'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Maximum column w velocity since last output')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w_velocity_max', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'w_velocity_max', r1Ptr)

! Define variable wind_speed_level1_max
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'wind_speed_level1_max'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Maximum wind speed in lowest model level since last '&
//'output')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'wind_speed_level1_max', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'wind_speed_level1_max', r1Ptr)

! Define variable ertel_pv
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ertel_pv'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Ertel''s potential vorticity (1 PVU = 10^{-6} m^{2} s^{-1} '&
//'K kg^{-1})')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ertel_pv', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ertel_pv', r2Ptr)

! Define variable ertel_pv_prev
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ertel_pv_prev'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Ertel''s PV at beginning of previous time step')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ertel_pv_prev', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ertel_pv_prev', r2Ptr)

! Define variable depv_dt_lw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_lw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from longwave radiation '&
//'parameterization scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_lw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_lw', r2Ptr)

! Define variable depv_dt_sw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_sw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from shortwave radiation '&
//'parameterization scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_sw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_sw', r2Ptr)

! Define variable depv_dt_bl
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_bl'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from PBL parameterization scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_bl', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_bl', r2Ptr)

! Define variable depv_dt_cu
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_cu'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from cumulus parameterization scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_cu', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_cu', r2Ptr)

! Define variable depv_dt_mp
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_mp'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from microphysics parameterization '&
//'scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_mp', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_mp', r2Ptr)

! Define variable depv_dt_mix
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_mix'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from explicit horizontal mixing')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_mix', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_mix', r2Ptr)

! Define variable depv_dt_diab
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_diab'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Sum of calculated PV tendencies from diabatic processes')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_diab', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_diab', r2Ptr)

! Define variable depv_dt_fric_bl
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_fric_bl'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Frictional PV tendency from PBL and GWD parameterization '&
//'schemes')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_fric_bl', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_fric_bl', r2Ptr)

! Define variable depv_dt_fric_cu
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_fric_cu'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Frictional PV tendency from cumulus parameterization '&
//'scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_fric_cu', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_fric_cu', r2Ptr)

! Define variable depv_dt_fric_mix
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_fric_mix'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Frictional PV tendency from explicit horizontal mixing')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_fric_mix', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_fric_mix', r2Ptr)

! Define variable depv_dt_fric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_fric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Sum of calculated PV tendencies from frictional processes')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_fric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_fric', r2Ptr)

! Define variable depv_dt_dyn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_dyn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'PV tendency from dynamics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_dyn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_dyn', r2Ptr)

! Define variable depv_dt_mp_evap_cw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_mp_evap_cw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from microphysics - net '&
//'condensation/evaporation of cloud water')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_mp_evap_cw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_mp_evap_cw', r2Ptr)

! Define variable depv_dt_mp_evap_rw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_mp_evap_rw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from microphysics - evaporation of '&
//'rain water')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_mp_evap_rw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_mp_evap_rw', r2Ptr)

! Define variable depv_dt_mp_depo_ice
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_mp_depo_ice'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from microphysics - net '&
//'deposition/sublimation of all ice hydrometeors')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_mp_depo_ice', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_mp_depo_ice', r2Ptr)

! Define variable depv_dt_mp_melt_ice
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_mp_melt_ice'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from microphysics - melting of all '&
//'ice hydrometeors')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_mp_melt_ice', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_mp_melt_ice', r2Ptr)

! Define variable depv_dt_mp_frez_ice
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_mp_frez_ice'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from microphysics - freezing of all '&
//'ice hydrometeors')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_mp_frez_ice', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_mp_frez_ice', r2Ptr)

! Define variable depv_dt_mp_allproc
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'depv_dt_mp_allproc'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency from summed microphysical processes')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_mp_allproc', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_mp_allproc', r2Ptr)

! Define variable acc_depv_dt_lw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_lw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from longwave radiation '&
//'scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_lw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_lw', r2Ptr)

! Define variable acc_depv_dt_sw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_sw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from shortwave radiation '&
//'scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_sw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_sw', r2Ptr)

! Define variable acc_depv_dt_bl
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_bl'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from PBL scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_bl', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_bl', r2Ptr)

! Define variable acc_depv_dt_cu
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_cu'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from cumulus scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_cu', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_cu', r2Ptr)

! Define variable acc_depv_dt_mp
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_mp'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from microphysics scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_mp', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_mp', r2Ptr)

! Define variable acc_depv_dt_mix
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_mix'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from explicit horizontal '&
//'mixing')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_mix', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_mix', r2Ptr)

! Define variable acc_depv_dt_diab
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_diab'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated sum of calculated PV tendencies from diabatic '&
//'processes')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_diab', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_diab', r2Ptr)

! Define variable acc_depv_dt_fric_bl
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_fric_bl'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated frictional PV tendency from PBL and GWD '&
//'schemes')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_fric_bl', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_fric_bl', r2Ptr)

! Define variable acc_depv_dt_fric_cu
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_fric_cu'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated frictional PV tendency from cumulus scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_fric_cu', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_fric_cu', r2Ptr)

! Define variable acc_depv_dt_fric_mix
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_fric_mix'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated frictional PV tendency from explicit '&
//'horizontal mixing')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_fric_mix', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_fric_mix', r2Ptr)

! Define variable acc_depv_dt_fric
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_fric'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated sum of calculated PV tendencies from '&
//'frictional processes')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_fric', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_fric', r2Ptr)

! Define variable acc_depv_dt_dyn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_dyn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated PV tendency from dynamics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_dyn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_dyn', r2Ptr)

! Define variable acc_depv_dt_mp_evap_cw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_mp_evap_cw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from microphysics - net '&
//'condensation/evaporation of cloud water')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_mp_evap_cw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_mp_evap_cw', r2Ptr)

! Define variable acc_depv_dt_mp_evap_rw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_mp_evap_rw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from microphysics - '&
//'evaporation of rain water')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_mp_evap_rw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_mp_evap_rw', r2Ptr)

! Define variable acc_depv_dt_mp_depo_ice
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_mp_depo_ice'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from microphysics - net '&
//'deposition/sublimation of all ice hydrometeors')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_mp_depo_ice', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_mp_depo_ice', r2Ptr)

! Define variable acc_depv_dt_mp_melt_ice
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_mp_melt_ice'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from microphysics - '&
//'melting of all ice hydrometeors')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_mp_melt_ice', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_mp_melt_ice', r2Ptr)

! Define variable acc_depv_dt_mp_frez_ice
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_mp_frez_ice'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from microphysics - '&
//'freezing of all ice hydrometeors')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_mp_frez_ice', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_mp_frez_ice', r2Ptr)

! Define variable acc_depv_dt_mp_allproc
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_depv_dt_mp_allproc'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency from summed microphysical '&
//'processes')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_mp_allproc', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_mp_allproc', r2Ptr)

! Define variable acc_dtheta_dt_mp
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_dtheta_dt_mp'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated potential temperature tendency due to '&
//'microphysics scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_dtheta_dt_mp', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_dtheta_dt_mp', r2Ptr)

! Define variable acc_dtheta_dt_cu
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_dtheta_dt_cu'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated potential temperature tendency due to cumulus '&
//'scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_dtheta_dt_cu', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_dtheta_dt_cu', r2Ptr)

! Define variable iLev_DT
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'iLev_DT'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'Lowest vertical model level at or above dynamic tropopause '&
//'(< 1 if 2 PVU below column; > nVertLevels if 2 PVU above column)')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'iLev_DT', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'iLev_DT', i1Ptr)

! Define variable iLev_DT_prev
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'iLev_DT_prev'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'Lowest vertical level at or above dynamic tropopause at '&
//'beginning of previous time step')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'iLev_DT_prev', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'iLev_DT_prev', i1Ptr)

! Define variable u_pv
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'u_pv'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Zonal wind on dynamic tropopause diagnosed at current time '&
//'step (i.e., iLev_DT)')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'u_pv', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'u_pv', r1Ptr)

! Define variable v_pv
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'v_pv'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Meridional wind on dynamic tropopause diagnosed at current '&
//'time step (i.e., iLev_DT)')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'v_pv', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'v_pv', r1Ptr)

! Define variable theta_pv
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'theta_pv'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Potential temperature on dynamic tropopause diagnosed at '&
//'current time step (i.e., iLev_DT)')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'theta_pv', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'theta_pv', r1Ptr)

! Define variable vort_pv
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'vort_pv'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 's^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Relative vertical vorticity on dynamic tropopause '&
//'diagnosed at current time step (i.e., iLev_DT)')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'vort_pv', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'vort_pv', r1Ptr)

! Define variable pres_pv
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'pres_pv'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'Pa')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Pressure of dynamic tropopause diagnosed at current time '&
//'step (i.e., iLev_DT)')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'pres_pv', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'pres_pv', r1Ptr)

! Define variable height_pv
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'height_pv'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm MSL')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Geopotential height of dynamic tropopause diagnosed at '&
//'current time step (i.e., iLev_DT)')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'height_pv', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'height_pv', r1Ptr)

! Define variable depv_dt_diab_pv
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'depv_dt_diab_pv'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Diabatic PV tendency on dynamic tropopause diagnosed at '&
//'beginning of previous time step (i.e., iLev_DT_prev)')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_diab_pv', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_diab_pv', r1Ptr)

! Define variable depv_dt_fric_pv
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'depv_dt_fric_pv'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Frictional PV tendency on dynamic tropopause diagnosed at '&
//'beginning of previous time step (i.e., iLev_DT_prev)')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_fric_pv', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_fric_pv', r1Ptr)

! Define variable depv_dt_dyn_pv
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'depv_dt_dyn_pv'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'PV tendency from dynamics on dynamic tropopause diagnosed '&
//'at beginning of previous time step (i.e., iLev_DT_prev)')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'depv_dt_dyn_pv', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'depv_dt_dyn_pv', r1Ptr)

! Define variable acc_depv_dt_diab_pv
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'acc_depv_dt_diab_pv'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated diabatic PV tendency on dynamic tropopause')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_diab_pv', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_diab_pv', r1Ptr)

! Define variable acc_depv_dt_fric_pv
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'acc_depv_dt_fric_pv'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated frictional PV tendency on dynamic tropopause')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_fric_pv', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_fric_pv', r1Ptr)

! Define variable acc_depv_dt_dyn_pv
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'acc_depv_dt_dyn_pv'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'PVU s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated advective PV tendency on dynamic tropopause')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acc_depv_dt_dyn_pv', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_depv_dt_dyn_pv', r1Ptr)

! Define variable wCell
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'wCell'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Vertical velocity interpolated to cell centers')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'wCell', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'wCell', r2Ptr)

! Define variable dtheta_dz
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dtheta_dz'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K m^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Vertical gradient of potential temperature')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dtheta_dz', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dtheta_dz', r2Ptr)

! Define variable uReconstructZonal_prev
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uReconstructZonal_prev'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Reconstructed zonal wind at beginning of previous time '&
//'step')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uReconstructZonal_prev', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructZonal_prev', r2Ptr)

! Define variable uReconstructMeridional_prev
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uReconstructMeridional_prev'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Reconstructed meridional wind at beginning of previous '&
//'time step')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uReconstructMeridional_prev', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructMeridional_prev', r2Ptr)

! Define variable wCell_prev
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'wCell_prev'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Vertical velocity interpolated to cell centers at '&
//'beginning of previous time step')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'wCell_prev', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'wCell_prev', r2Ptr)

! Define variable theta_prev
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'theta_prev'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Potential temperature at beginning of previous time step')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'theta_prev', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'theta_prev', r2Ptr)

! Define variable rho_prev
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rho_prev'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Density at cell center at beginning of previous time step')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rho_prev', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_prev', r2Ptr)

! Define variable qv_prev
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'qv_prev'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Water vapor mixing ratio at beginning of the time step')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'qv_prev', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'qv_prev', r2Ptr)

! Define variable pv_vertex_prev
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'pv_vertex_prev'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nVertices'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 's^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Absolute vertical voricity at cell vertices at beginning '&
//'of previous time step')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'pv_vertex_prev', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'pv_vertex_prev', r2Ptr)

! Define variable dtheta_dt_mix
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dtheta_dt_mix'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Potential temperature heating rate from explicit '&
//'horizontal mixing')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dtheta_dt_mix', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dtheta_dt_mix', r2Ptr)

! Define variable dtheta_dt_dyn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dtheta_dt_dyn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Potential temperature tendency from dynamics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dtheta_dt_dyn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dtheta_dt_dyn', r2Ptr)

! Define variable thmblten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'thmblten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Theta_m tendency from PBL scheme calculated over RK '&
//'timestep integration')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'thmblten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'thmblten', r2Ptr)

! Define variable qvblten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'qvblten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Mixing ratio tendency from PBL scheme calculated over RK '&
//'timestep integration')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'qvblten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'qvblten', r2Ptr)

! Define variable dtheta_dt_pbl
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dtheta_dt_pbl'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Potential temperature tendency from PBL scheme -- derived '&
//'from theta_m and qv tendencies')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dtheta_dt_pbl', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dtheta_dt_pbl', r2Ptr)

! Define variable thmcuten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'thmcuten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Theta_m tendency from cumulus scheme calculated over RK '&
//'timestep integration')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'thmcuten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'thmcuten', r2Ptr)

! Define variable qvcuten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'qvcuten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Mixing ratio tendency from cumulus scheme calculated over '&
//'RK timestep integration')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'qvcuten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'qvcuten', r2Ptr)

! Define variable dtheta_dt_cu
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dtheta_dt_cu'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Potential temperature tendency from cumulus scheme -- '&
//'derived from theta_m and qv tendencies')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dtheta_dt_cu', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dtheta_dt_cu', r2Ptr)

! Define variable thmlwten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'thmlwten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Theta_m tendency from longwave radiation calculated over '&
//'RK timestep integration')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'thmlwten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'thmlwten', r2Ptr)

! Define variable dtheta_dt_lw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dtheta_dt_lw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Potential temperature tendency from longwave radiation '&
//'scheme -- derived from theta_m tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dtheta_dt_lw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dtheta_dt_lw', r2Ptr)

! Define variable thmswten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'thmswten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Theta_m tendency from shortwave radiation calculated over '&
//'RK timestep integration')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'thmswten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'thmswten', r2Ptr)

! Define variable dtheta_dt_sw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dtheta_dt_sw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Potential temperature tendency from shortwave radiation '&
//'scheme -- derived from theta_m tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dtheta_dt_sw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dtheta_dt_sw', r2Ptr)

! Define variable thmmpten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'thmmpten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Theta_m tendency from microphysics scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'thmmpten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'thmmpten', r2Ptr)

! Define variable qvmpten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'qvmpten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Mixing ratio tendency from microphysics scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'qvmpten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'qvmpten', r2Ptr)

! Define variable dtheta_dt_mp
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dtheta_dt_mp'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Potential temperature tendency from microphysics scheme -- '&
//'derived from theta_m and qv tendencies')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dtheta_dt_mp', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dtheta_dt_mp', r2Ptr)

! Define variable tend_u_pbl
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_u_pbl'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Normal wind tendencies from PBL and GWD parameterization '&
//'schemes')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'tend_u_pbl', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_u_pbl', r2Ptr)

! Define variable tend_u_cu
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_u_cu'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Normal wind tendencies from cumulus parameterization '&
//'scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'tend_u_cu', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_u_cu', r2Ptr)

! Define variable uTend_curl_diff
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uTend_curl_diff'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nVertices'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 's^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Curl of the normal wind tendencies from diffusion at the '&
//'cell vertices')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uTend_curl_diff', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uTend_curl_diff', r2Ptr)

! Define variable uTend_curl_pbl
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uTend_curl_pbl'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nVertices'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 's^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Curl of the normal wind tendencies from PBL+GWD at the '&
//'cell vertices')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uTend_curl_pbl', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uTend_curl_pbl', r2Ptr)

! Define variable uTend_curl_cu
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uTend_curl_cu'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nVertices'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 's^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Curl of the normal wind tendencies from cumulus at the '&
//'cell vertices')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uTend_curl_cu', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uTend_curl_cu', r2Ptr)

! Define variable uTend_curl_dyn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uTend_curl_dyn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nVertices'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 's^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Curl of the normal wind tendencies from dynamics at the '&
//'cell vertices')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'uTend_curl_dyn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uTend_curl_dyn', r2Ptr)

! Define variable tend_wCell_diff
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_wCell_diff'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Uncoupled vertical velocity tendency from mixing '&
//'interpolated to mass levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'tend_wCell_diff', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_wCell_diff', r2Ptr)

! Define variable tenddyn_wCell
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tenddyn_wCell'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Uncoupled vertical velocity tendency from dynamics '&
//'interpolated to mass levels')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'tenddyn_wCell', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tenddyn_wCell', r2Ptr)

! Define variable inTropo
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'inTropo'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .true.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'nVertLevels'
      i2Ptr(1) % dimNames(2) = 'nCells'
      i2Ptr(1) % defaultValue = 0
      allocate(i2Ptr(1) % attLists(1))
      allocate(i2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'units', '0/1')
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'long_name', '1 if within troposphere based on PV flood fill')
      i2Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'inTropo', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'inTropo', i2Ptr)

! Define variable candInTropo
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'candInTropo'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .true.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'nVertLevels'
      i2Ptr(1) % dimNames(2) = 'nCells'
      i2Ptr(1) % defaultValue = 0
      allocate(i2Ptr(1) % attLists(1))
      allocate(i2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'units', '0/1')
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'long_name', '1 if candidate of troposphere')
      i2Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'candInTropo', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'candInTropo', i2Ptr)

! Define variable candInStrato
      allocate(i2Ptr(1))

! Setting up time level 1
      i2Ptr(1) % fieldName = 'candInStrato'
      i2Ptr(1) % isVarArray = .false.
      i2Ptr(1) % isDecomposed = .true.
      i2Ptr(1) % hasTimeDimension = .true.
      i2Ptr(1) % isPersistent = .true.
      i2Ptr(1) % isActive = .false.
! Setting up dimensions
      i2Ptr(1) % dimNames(1) = 'nVertLevels'
      i2Ptr(1) % dimNames(2) = 'nCells'
      i2Ptr(1) % defaultValue = 0
      allocate(i2Ptr(1) % attLists(1))
      allocate(i2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'units', '0/1')
      call mpas_add_att(i2Ptr(1) % attLists(1) % attList, 'long_name', '1 if candidate of stratosphere')
      i2Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i2Ptr(1) % block => block

      i2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'candInStrato', i2Ptr)
      call mpas_pool_add_field(block % allFields, 'candInStrato', i2Ptr)

! Define variable rw_tend_dyn_large
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rw_tend_dyn_large'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled vertical momentum tendency due to large-time-step '&
//'dynamics tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rw_tend_dyn_large', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rw_tend_dyn_large', r2Ptr)

! Define variable rw_tend_dyn_small
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rw_tend_dyn_small'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled vertical momentum tendency due to small-time-step '&
//'dynamics tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rw_tend_dyn_small', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rw_tend_dyn_small', r2Ptr)

! Define variable rw_tend_diff
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rw_tend_diff'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled vertical momentum tendency from (horiz) diffusion')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rw_tend_diff', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rw_tend_diff', r2Ptr)

! Define variable w_tend_dcpl
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'w_tend_dcpl'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Vertical velocity tendency component from decoupling from '&
//'mass')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'w_tend_dcpl', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'w_tend_dcpl', r2Ptr)

! Define variable dthetam_dt_dyn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dthetam_dt_dyn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Moist potential temperature dynamic tendency '&
//'(large+small+dcpl)')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'dthetam_dt_dyn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dthetam_dt_dyn', r2Ptr)

! Define variable dthetam_dt_mix
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dthetam_dt_mix'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Moist potential temperature explicit diffusion tendency '&
//'(diff)')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'dthetam_dt_mix', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dthetam_dt_mix', r2Ptr)

! Define variable dqv_dt_dyn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dqv_dt_dyn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Water vapor advective tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'dqv_dt_dyn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dqv_dt_dyn', r2Ptr)

! Define variable du_dt_dyn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'du_dt_dyn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Horizontal winds dynamics tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'du_dt_dyn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'du_dt_dyn', r2Ptr)

! Define variable du_dt_dyn_zonal
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'du_dt_dyn_zonal'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Zonal winds dynamics tendency (reconstructed in PV code)')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'du_dt_dyn_zonal', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'du_dt_dyn_zonal', r2Ptr)

! Define variable du_dt_dyn_meridional
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'du_dt_dyn_meridional'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Meridional winds dynamics tendency (reconstructed in PV '&
//'code)')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'du_dt_dyn_meridional', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'du_dt_dyn_meridional', r2Ptr)

! Define variable dw_dt_dyn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dw_dt_dyn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Vertical winds dynamics tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'dw_dt_dyn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dw_dt_dyn', r2Ptr)

! Define variable ru_tend_dyn_large
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ru_tend_dyn_large'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled horizontal momentum large-time-step dynamic '&
//'tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'ru_tend_dyn_large', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ru_tend_dyn_large', r2Ptr)

! Define variable ru_tend_dyn_small
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ru_tend_dyn_small'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled horizontal momentum acoustic time-step dynamic '&
//'tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'ru_tend_dyn_small', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ru_tend_dyn_small', r2Ptr)

! Define variable ru_tend_diff
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ru_tend_diff'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled horizontal momentum tendency from diffusion')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'ru_tend_diff', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ru_tend_diff', r2Ptr)

! Define variable pv_vertex_dyn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'pv_vertex_dyn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nVertices'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 's^{-1} kg^{-1} m^3')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'absolute vorticity/rho_zz at a vertex after dynamics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'pv_vertex_dyn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'pv_vertex_dyn', r2Ptr)

! Define variable u_tend_diff
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'u_tend_diff'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Uncoupled horizontal momentum tendency from diffusion')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'u_tend_diff', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'u_tend_diff', r2Ptr)

! Define variable u_tend_diff_reconstructZonal
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'u_tend_diff_reconstructZonal'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Uncoupled horizontal momentum tendency from diffusion')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'u_tend_diff_reconstructZonal', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'u_tend_diff_reconstructZonal', r2Ptr)

! Define variable u_tend_diff_reconstructMeridional
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'u_tend_diff_reconstructMeridional'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Uncoupled horizontal momentum tendency from diffusion')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'u_tend_diff_reconstructMeridional', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'u_tend_diff_reconstructMeridional', r2Ptr)

! Define variable w_tend_diff
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'w_tend_diff'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Uncoupled vertical momentum tendency from diffusion')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'w_tend_diff', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'w_tend_diff', r2Ptr)

! Define variable u_dyn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'u_dyn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Uncoupled horizontal velocity field after dynamics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'u_dyn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'u_dyn', r2Ptr)

! Define variable uReconstructX_dyn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uReconstructX_dyn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Uncoupled horizontal velocity field after dynamics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'uReconstructX_dyn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructX_dyn', r2Ptr)

! Define variable uReconstructY_dyn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uReconstructY_dyn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Uncoupled horizontal velocity field after dynamics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'uReconstructY_dyn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructY_dyn', r2Ptr)

! Define variable uReconstructZ_dyn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uReconstructZ_dyn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Uncoupled horizontal velocity field after dynamics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'uReconstructZ_dyn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructZ_dyn', r2Ptr)

! Define variable uReconstructZonal_dyn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uReconstructZonal_dyn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Uncoupled horizontal velocity field after dynamics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'uReconstructZonal_dyn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructZonal_dyn', r2Ptr)

! Define variable uReconstructMeridional_dyn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'uReconstructMeridional_dyn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Uncoupled horizontal velocity field after dynamics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'uReconstructMeridional_dyn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'uReconstructMeridional_dyn', r2Ptr)

! Define variable w_dyn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'w_dyn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Uncoupled vertical velocity field after dynamics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'w_dyn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'w_dyn', r2Ptr)

! Define variable theta_dyn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'theta_dyn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Uncoupled potential temperature field after dynamics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'theta_dyn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'theta_dyn', r2Ptr)

! Define variable rho_dyn
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rho_dyn'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Dry density field after dynamics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rho_dyn', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_dyn', r2Ptr)

! Define variable ru_tend_physics
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ru_tend_physics'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled horizontal momentum tendency from physics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'ru_tend_physics', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ru_tend_physics', r2Ptr)

! Define variable ru_tend_smdiv
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'ru_tend_smdiv'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled horizontal momentum tendency due to divergence '&
//'damping')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'ru_tend_smdiv', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'ru_tend_smdiv', r2Ptr)

! Define variable rublten_tend
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rublten_tend'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled horizontal momentum tendency due to PBL (includes '&
//'GWDO)')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rublten_tend', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rublten_tend', r2Ptr)

! Define variable rugwdo_tend
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rugwdo_tend'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled horizontal momentum tendency due to GWDO scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rugwdo_tend', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rugwdo_tend', r2Ptr)

! Define variable rucuten_tend
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rucuten_tend'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled horizontal momentum tendency due to cumulus '&
//'scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rucuten_tend', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rucuten_tend', r2Ptr)

! Define variable u_tend_dcpl
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'u_tend_dcpl'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Horizontal momentum tendency due to dcoupling from mass')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'u_tend_dcpl', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'u_tend_dcpl', r2Ptr)

! Define variable rth_tend_dyn_large
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rth_tend_dyn_large'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3} K')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled potential temperature large-time-step dynamic '&
//'tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rth_tend_dyn_large', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rth_tend_dyn_large', r2Ptr)

! Define variable rth_tend_dyn_small
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rth_tend_dyn_small'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3} K')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled potential temperature acoustic time-step dynamic '&
//'tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rth_tend_dyn_small', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rth_tend_dyn_small', r2Ptr)

! Define variable rth_tend_diff
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rth_tend_diff'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3} K')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled potential temperature tendency from diffusion')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rth_tend_diff', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rth_tend_diff', r2Ptr)

! Define variable rth_tend_physics
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rth_tend_physics'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3} K')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled potential temperature tendency from physics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rth_tend_physics', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rth_tend_physics', r2Ptr)

! Define variable rthblten_tend
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rthblten_tend'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3} K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled potential temperature tendency due to PBL')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rthblten_tend', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rthblten_tend', r2Ptr)

! Define variable rthcuten_tend
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rthcuten_tend'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3} K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled potential temperature tendency due to cumulus '&
//'scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rthcuten_tend', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rthcuten_tend', r2Ptr)

! Define variable rthratensw_tend
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rthratensw_tend'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3} K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled potential temperature tendency due to cumulus '&
//'scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rthratensw_tend', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rthratensw_tend', r2Ptr)

! Define variable rthratenlw_tend
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rthratenlw_tend'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3} K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled potential temperature tendency due to cumulus '&
//'scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rthratenlw_tend', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rthratenlw_tend', r2Ptr)

! Define variable th_tend_dcpl
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'th_tend_dcpl'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Potential temperature tendency due to dcoupling from mass')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'th_tend_dcpl', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'th_tend_dcpl', r2Ptr)

! Define variable rqv_tend_diff
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rqv_tend_diff'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3} kg kg^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled water vapor tendency from diffusion')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rqv_tend_diff', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rqv_tend_diff', r2Ptr)

! Define variable qvblten_tend
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'qvblten_tend'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled water vapor tendency due to PBL')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'qvblten_tend', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'qvblten_tend', r2Ptr)

! Define variable qvcuten_tend
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'qvcuten_tend'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Coupled water vapor tendency due to cumulus scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'qvcuten_tend', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'qvcuten_tend', r2Ptr)

! Define variable qv_mp_tend
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'qv_mp_tend'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Water vapor tendency due to microphysics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'qv_mp_tend', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'qv_mp_tend', r2Ptr)

! Define variable acc_u_tend_dyn_large
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_dyn_large'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum large-time-step dynamic '&
//'tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_dyn_large', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_dyn_large', r2Ptr)

! Define variable acc_u_tend_dyn_small
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_dyn_small'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated hrizontal momentum acoustic time-step dynamic '&
//'tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_dyn_small', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_dyn_small', r2Ptr)

! Define variable acc_u_tend_diff
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_diff'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency from diffusion')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_diff', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_diff', r2Ptr)

! Define variable acc_u_tend_physics
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_physics'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency from physics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_physics', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_physics', r2Ptr)

! Define variable acc_u_tend_smdiv
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_smdiv'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to divergence '&
//'damping')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_smdiv', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_smdiv', r2Ptr)

! Define variable acc_ublten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_ublten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to PBL scheme '&
//'(includes GWDO) ')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_ublten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_ublten', r2Ptr)

! Define variable acc_ugwdoten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_ugwdoten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to GWDO '&
//'scheme ')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_ugwdoten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_ugwdoten', r2Ptr)

! Define variable acc_ucuten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_ucuten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to cumulus '&
//'scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_ucuten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_ucuten', r2Ptr)

! Define variable acc_u_tend_dcpl
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_dcpl'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to decoupling '&
//'from mass')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_dcpl', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_dcpl', r2Ptr)

! Define variable acc_th_tend_dyn_large
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_th_tend_dyn_large'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated moist potential temperature large-time-step '&
//'dynamic tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_th_tend_dyn_large', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_th_tend_dyn_large', r2Ptr)

! Define variable acc_th_tend_dyn_small
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_th_tend_dyn_small'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated moist potential temperature acoustic time-step '&
//'dynamic tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_th_tend_dyn_small', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_th_tend_dyn_small', r2Ptr)

! Define variable acc_th_tend_diff
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_th_tend_diff'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated potential temperature tendency from diffusion')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_th_tend_diff', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_th_tend_diff', r2Ptr)

! Define variable acc_th_tend_physics
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_th_tend_physics'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated potential temperature tendency from physics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_th_tend_physics', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_th_tend_physics', r2Ptr)

! Define variable acc_thblten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_thblten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated potential temperature tendency due to PBL')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_thblten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_thblten', r2Ptr)

! Define variable acc_thcuten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_thcuten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated potential temperature tendency due to cumulus '&
//'scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_thcuten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_thcuten', r2Ptr)

! Define variable acc_thratensw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_thratensw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated potential temperature tendency due to cumulus '&
//'scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_thratensw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_thratensw', r2Ptr)

! Define variable acc_thratenlw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_thratenlw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated potential temperature tendency due to cumulus '&
//'scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_thratenlw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_thratenlw', r2Ptr)

! Define variable acc_th_tend_diabatic
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_th_tend_diabatic'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated potential temperature tendency due to latent '&
//'heat')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_th_tend_diabatic', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_th_tend_diabatic', r2Ptr)

! Define variable acc_th_tend_dcpl
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_th_tend_dcpl'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated potential temperature tendency due to '&
//'decoupling from mass')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_th_tend_dcpl', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_th_tend_dcpl', r2Ptr)

! Define variable acc_qv_tend_dyn_large
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_qv_tend_dyn_large'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated water vapor large-time-step dynamic tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_qv_tend_dyn_large', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_qv_tend_dyn_large', r2Ptr)

! Define variable acc_qv_tend_diff
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_qv_tend_diff'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated water vapor tendency from diffusion')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_qv_tend_diff', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_qv_tend_diff', r2Ptr)

! Define variable acc_qvblten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_qvblten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated water vapor tendency due to PBL')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_qvblten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_qvblten', r2Ptr)

! Define variable acc_qvcuten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_qvcuten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated water vapor tendency due to cumulus scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_qvcuten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_qvcuten', r2Ptr)

! Define variable acc_qv_mp_tend
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_qv_mp_tend'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated water vapor tendency due to latent heat')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_qv_mp_tend', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_qv_mp_tend', r2Ptr)

! Define variable acc_u_tend_dyn_large_ReconstructZonal
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_dyn_large_ReconstructZonal'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated zonal momentum large-time-step dynamic '&
//'tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_dyn_large_ReconstructZonal', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_dyn_large_ReconstructZonal', r2Ptr)

! Define variable acc_u_tend_dyn_small_ReconstructZonal
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_dyn_small_ReconstructZonal'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated zonal momentum acoustic time-step dynamic '&
//'tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_dyn_small_ReconstructZonal', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_dyn_small_ReconstructZonal', r2Ptr)

! Define variable acc_u_tend_diff_ReconstructZonal
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_diff_ReconstructZonal'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated zonal momentum tendency from diffusion')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_diff_ReconstructZonal', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_diff_ReconstructZonal', r2Ptr)

! Define variable acc_u_tend_physics_ReconstructZonal
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_physics_ReconstructZonal'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated zonal momentum tendency from physics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_physics_ReconstructZonal', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_physics_ReconstructZonal', r2Ptr)

! Define variable acc_u_tend_smdiv_ReconstructZonal
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_smdiv_ReconstructZonal'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated zonal momentum tendency due to divergence '&
//'damping')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_smdiv_ReconstructZonal', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_smdiv_ReconstructZonal', r2Ptr)

! Define variable acc_ublten_ReconstructZonal
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_ublten_ReconstructZonal'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated zonal momentum tendency due to PBL scheme '&
//'(includes GWDO) ')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_ublten_ReconstructZonal', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_ublten_ReconstructZonal', r2Ptr)

! Define variable acc_ugwdoten_ReconstructZonal
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_ugwdoten_ReconstructZonal'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated zonal momentum tendency due to GWDO scheme ')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_ugwdoten_ReconstructZonal', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_ugwdoten_ReconstructZonal', r2Ptr)

! Define variable acc_ucuten_ReconstructZonal
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_ucuten_ReconstructZonal'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated zonal momentum tendency due to cumulus scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_ucuten_ReconstructZonal', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_ucuten_ReconstructZonal', r2Ptr)

! Define variable acc_u_tend_dcpl_ReconstructZonal
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_dcpl_ReconstructZonal'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated zonal momentum tendency due to decoupling from '&
//'mass')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_dcpl_ReconstructZonal', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_dcpl_ReconstructZonal', r2Ptr)

! Define variable acc_u_tend_dyn_large_ReconstructMeridional
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_dyn_large_ReconstructMeridional'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated meridional momentum large-time-step dynamic '&
//'tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_dyn_large_ReconstructMeridional', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_dyn_large_ReconstructMeridional', r2Ptr)

! Define variable acc_u_tend_dyn_small_ReconstructMeridional
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_dyn_small_ReconstructMeridional'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated meridional momentum acoustic time-step dynamic '&
//'tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_dyn_small_ReconstructMeridional', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_dyn_small_ReconstructMeridional', r2Ptr)

! Define variable acc_u_tend_diff_ReconstructMeridional
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_diff_ReconstructMeridional'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated meridional momentum tendency from diffusion')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_diff_ReconstructMeridional', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_diff_ReconstructMeridional', r2Ptr)

! Define variable acc_u_tend_physics_ReconstructMeridional
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_physics_ReconstructMeridional'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated meridional momentum tendency from physics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_physics_ReconstructMeridional', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_physics_ReconstructMeridional', r2Ptr)

! Define variable acc_u_tend_smdiv_ReconstructMeridional
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_smdiv_ReconstructMeridional'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated meridional momentum tendency due to divergence '&
//'damping')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_smdiv_ReconstructMeridional', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_smdiv_ReconstructMeridional', r2Ptr)

! Define variable acc_ublten_ReconstructMeridional
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_ublten_ReconstructMeridional'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated meridional momentum tendency due to PBL scheme '&
//'(includes GWDO)')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_ublten_ReconstructMeridional', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_ublten_ReconstructMeridional', r2Ptr)

! Define variable acc_ugwdoten_ReconstructMeridional
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_ugwdoten_ReconstructMeridional'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated meridional momentum tendency due to GWDO '&
//'scheme ')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_ugwdoten_ReconstructMeridional', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_ugwdoten_ReconstructMeridional', r2Ptr)

! Define variable acc_ucuten_ReconstructMeridional
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_ucuten_ReconstructMeridional'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated meridional momentum tendency due to cumulus '&
//'scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_ucuten_ReconstructMeridional', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_ucuten_ReconstructMeridional', r2Ptr)

! Define variable acc_u_tend_dcpl_ReconstructMeridional
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_dcpl_ReconstructMeridional'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated meridional momentum tendency due to decoupling '&
//'from mass')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_dcpl_ReconstructMeridional', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_dcpl_ReconstructMeridional', r2Ptr)

! Define variable acc_u_tend_dyn_large_ReconstructX
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_dyn_large_ReconstructX'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum large-time-step dynamic '&
//'tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_dyn_large_ReconstructX', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_dyn_large_ReconstructX', r2Ptr)

! Define variable acc_u_tend_dyn_small_ReconstructX
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_dyn_small_ReconstructX'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated hrizontal momentum acoustic time-step dynamic '&
//'tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_dyn_small_ReconstructX', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_dyn_small_ReconstructX', r2Ptr)

! Define variable acc_u_tend_diff_ReconstructX
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_diff_ReconstructX'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency from diffusion')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_diff_ReconstructX', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_diff_ReconstructX', r2Ptr)

! Define variable acc_u_tend_physics_ReconstructX
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_physics_ReconstructX'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency from physics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_physics_ReconstructX', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_physics_ReconstructX', r2Ptr)

! Define variable acc_u_tend_smdiv_ReconstructX
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_smdiv_ReconstructX'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to divergence '&
//'damping')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_smdiv_ReconstructX', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_smdiv_ReconstructX', r2Ptr)

! Define variable acc_ublten_ReconstructX
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_ublten_ReconstructX'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to PBL scheme '&
//'(includes GWDO)')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_ublten_ReconstructX', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_ublten_ReconstructX', r2Ptr)

! Define variable acc_ugwdoten_ReconstructX
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_ugwdoten_ReconstructX'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to GWDO '&
//'scheme ')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_ugwdoten_ReconstructX', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_ugwdoten_ReconstructX', r2Ptr)

! Define variable acc_ucuten_ReconstructX
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_ucuten_ReconstructX'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to cumulus '&
//'scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_ucuten_ReconstructX', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_ucuten_ReconstructX', r2Ptr)

! Define variable acc_u_tend_dcpl_ReconstructX
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_dcpl_ReconstructX'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to decoupling '&
//'from mass')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_dcpl_ReconstructX', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_dcpl_ReconstructX', r2Ptr)

! Define variable acc_u_tend_dyn_large_ReconstructY
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_dyn_large_ReconstructY'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum large-time-step dynamic '&
//'tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_dyn_large_ReconstructY', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_dyn_large_ReconstructY', r2Ptr)

! Define variable acc_u_tend_dyn_small_ReconstructY
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_dyn_small_ReconstructY'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated hrizontal momentum acoustic time-step dynamic '&
//'tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_dyn_small_ReconstructY', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_dyn_small_ReconstructY', r2Ptr)

! Define variable acc_u_tend_diff_ReconstructY
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_diff_ReconstructY'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency from diffusion')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_diff_ReconstructY', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_diff_ReconstructY', r2Ptr)

! Define variable acc_u_tend_physics_ReconstructY
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_physics_ReconstructY'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency from physics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_physics_ReconstructY', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_physics_ReconstructY', r2Ptr)

! Define variable acc_u_tend_smdiv_ReconstructY
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_smdiv_ReconstructY'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to divergence '&
//'damping')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_smdiv_ReconstructY', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_smdiv_ReconstructY', r2Ptr)

! Define variable acc_ublten_ReconstructY
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_ublten_ReconstructY'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to PBL scheme '&
//'(includes GWDO)')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_ublten_ReconstructY', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_ublten_ReconstructY', r2Ptr)

! Define variable acc_ugwdoten_ReconstructY
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_ugwdoten_ReconstructY'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to GWDO '&
//'scheme ')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_ugwdoten_ReconstructY', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_ugwdoten_ReconstructY', r2Ptr)

! Define variable acc_ucuten_ReconstructY
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_ucuten_ReconstructY'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to cumulus '&
//'scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_ucuten_ReconstructY', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_ucuten_ReconstructY', r2Ptr)

! Define variable acc_u_tend_dcpl_ReconstructY
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_dcpl_ReconstructY'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to decoupling '&
//'from mass')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_dcpl_ReconstructY', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_dcpl_ReconstructY', r2Ptr)

! Define variable acc_u_tend_dyn_large_ReconstructZ
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_dyn_large_ReconstructZ'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum large-time-step dynamic '&
//'tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_dyn_large_ReconstructZ', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_dyn_large_ReconstructZ', r2Ptr)

! Define variable acc_u_tend_dyn_small_ReconstructZ
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_dyn_small_ReconstructZ'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated hrizontal momentum acoustic time-step dynamic '&
//'tendency')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_dyn_small_ReconstructZ', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_dyn_small_ReconstructZ', r2Ptr)

! Define variable acc_u_tend_diff_ReconstructZ
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_diff_ReconstructZ'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency from diffusion')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_diff_ReconstructZ', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_diff_ReconstructZ', r2Ptr)

! Define variable acc_u_tend_physics_ReconstructZ
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_physics_ReconstructZ'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency from physics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_physics_ReconstructZ', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_physics_ReconstructZ', r2Ptr)

! Define variable acc_u_tend_smdiv_ReconstructZ
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_smdiv_ReconstructZ'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to divergence '&
//'damping')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_smdiv_ReconstructZ', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_smdiv_ReconstructZ', r2Ptr)

! Define variable acc_ublten_ReconstructZ
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_ublten_ReconstructZ'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to PBL scheme '&
//'(includes GWDO)')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_ublten_ReconstructZ', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_ublten_ReconstructZ', r2Ptr)

! Define variable acc_ugwdoten_ReconstructZ
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_ugwdoten_ReconstructZ'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to GWDO '&
//'scheme ')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_ugwdoten_ReconstructZ', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_ugwdoten_ReconstructZ', r2Ptr)

! Define variable acc_ucuten_ReconstructZ
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_ucuten_ReconstructZ'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to cumulus '&
//'scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_ucuten_ReconstructZ', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_ucuten_ReconstructZ', r2Ptr)

! Define variable acc_u_tend_dcpl_ReconstructZ
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'acc_u_tend_dcpl_ReconstructZ'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Accumulated horizontal momentum tendency due to decoupling '&
//'from mass')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (tendenciesActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'acc_u_tend_dcpl_ReconstructZ', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'acc_u_tend_dcpl_ReconstructZ', r2Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
         call mpas_pool_add_config(newSubPool, 'is_periodic', block % domain % is_periodic)
         call mpas_pool_add_config(newSubPool, 'x_period', block % domain % x_period)
         call mpas_pool_add_config(newSubPool, 'y_period', block % domain % y_period)
      end if

   end subroutine atm_generate_pool_diag


   subroutine atm_generate_pool_tend(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_io, only : MPAS_REAL_FILLVAL, MPAS_INT_FILLVAL, MPAS_CHAR_FILLVAL
      implicit none
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), dimension(:), pointer :: r0Ptr
      type (field1DReal), dimension(:), pointer :: r1Ptr
      type (field2DReal), dimension(:), pointer :: r2Ptr
      type (field3DReal), dimension(:), pointer :: r3Ptr
      type (field4DReal), dimension(:), pointer :: r4Ptr
      type (field5DReal), dimension(:), pointer :: r5Ptr
      type (field0DInteger), dimension(:), pointer :: i0Ptr
      type (field1DInteger), dimension(:), pointer :: i1Ptr
      type (field2DInteger), dimension(:), pointer :: i2Ptr
      type (field3DInteger), dimension(:), pointer :: i3Ptr
      type (field0DChar), dimension(:), pointer :: c0Ptr
      type (field1DChar), dimension(:), pointer :: c1Ptr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: mp_kessler_inActive
      logical, pointer :: mp_thompson_inActive
      logical, pointer :: mp_wsm6_inActive
      logical, pointer :: cu_grell_freitas_inActive
      logical, pointer :: cu_kain_fritsch_inActive
      logical, pointer :: cu_tiedtke_inActive
      logical, pointer :: bl_ysu_inActive
      logical, pointer :: bl_mynn_inActive
      logical, pointer :: iauActive
      logical, pointer :: limited_areaActive
      logical, pointer :: jedi_daActive
      logical, pointer :: tendenciesActive


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'mp_kessler_inActive', mp_kessler_inActive)
      call mpas_pool_get_package(packagePool, 'mp_thompson_inActive', mp_thompson_inActive)
      call mpas_pool_get_package(packagePool, 'mp_wsm6_inActive', mp_wsm6_inActive)
      call mpas_pool_get_package(packagePool, 'cu_grell_freitas_inActive', cu_grell_freitas_inActive)
      call mpas_pool_get_package(packagePool, 'cu_kain_fritsch_inActive', cu_kain_fritsch_inActive)
      call mpas_pool_get_package(packagePool, 'cu_tiedtke_inActive', cu_tiedtke_inActive)
      call mpas_pool_get_package(packagePool, 'bl_ysu_inActive', bl_ysu_inActive)
      call mpas_pool_get_package(packagePool, 'bl_mynn_inActive', bl_mynn_inActive)
      call mpas_pool_get_package(packagePool, 'iauActive', iauActive)
      call mpas_pool_get_package(packagePool, 'limited_areaActive', limited_areaActive)
      call mpas_pool_get_package(packagePool, 'jedi_daActive', jedi_daActive)
      call mpas_pool_get_package(packagePool, 'tendenciesActive', tendenciesActive)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'tend', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'tend', newSubPool)

! Define var array scalars_tend
      allocate(r3Ptr(1))
      index_counter = 0
      group_counter = -1
      group_start = -1
      group_started = .false.

! Starting group moist
! Define constituent var tend_qv
! My Packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_qv', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
         end if
         group_started = .true.
      end if
! Define constituent var tend_qc
! My packages are bl_mynn_in;bl_ysu_in;cu_tiedtke_in;mp_kessler_in;mp_thompson_in;mp_wsm6_in
      if (bl_mynn_inActive .or. bl_ysu_inActive .or. cu_tiedtke_inActive .or. mp_kessler_inActive .or. mp_thompson_inActive .or. &
mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qc', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qc', -1)
         end if
      end if
! Define constituent var tend_qr
! My packages are mp_kessler_in;mp_thompson_in;mp_wsm6_in
      if (mp_kessler_inActive .or. mp_thompson_inActive .or. mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qr', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qr', -1)
         end if
      end if
! Define constituent var tend_qi
! My packages are bl_mynn_in;bl_ysu_in;cu_tiedtke_in;mp_thompson_in;mp_wsm6_in
      if (bl_mynn_inActive .or. bl_ysu_inActive .or. cu_tiedtke_inActive .or. mp_thompson_inActive .or. mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qi', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qi', -1)
         end if
      end if
! Define constituent var tend_qs
! My packages are mp_thompson_in;mp_wsm6_in
      if (mp_thompson_inActive .or. mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qs', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qs', -1)
         end if
      end if
! Define constituent var tend_qg
! My packages are mp_thompson_in;mp_wsm6_in
      if (mp_thompson_inActive .or. mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qg', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qg', -1)
         end if
      end if
         if (.not. group_started) then
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', -1)
               call mpas_pool_add_dimension(newSubPool, 'moist_end', -1)
            end if
         else
            group_started = .false.
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_end', index_counter)
            end if
         end if
! End of group       
! Starting group number
! Define constituent var tend_ni
! My Packages are bl_mynn_in;mp_thompson_in
      if (bl_mynn_inActive .or. mp_thompson_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_ni', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'number_start', group_start)
            end if
            group_started = .true.
         end if
      else
           if (associated(newSubPool)) then
              call mpas_pool_add_dimension(newSubPool, 'index_ni', -1)
           end if
      end if
! Define constituent var tend_nr
! My packages are mp_thompson_in
      if (mp_thompson_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_nr', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'number_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_nr', -1)
         end if
      end if
         if (.not. group_started) then
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'number_start', -1)
               call mpas_pool_add_dimension(newSubPool, 'number_end', -1)
            end if
         else
            group_started = .false.
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'number_end', index_counter)
            end if
         end if
! End of group       

      numConstituents = index_counter
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'num_scalars_tend', numConstituents)
      end if
! Defining time level 1
      allocate( r3Ptr(1) % constituentNames(numConstituents) )
      r3Ptr(1) % fieldName = 'scalars_tend'
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .true.
      r3Ptr(1) % isVarArray = .true.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qv', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_qv'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qc', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_qc'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qr', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_qr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qi', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_qi'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qs', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_qs'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qg', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_qg'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_ni', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_ni'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_nr', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_nr'
      end if

! Setup dimensions for       
      r3Ptr(1) % dimNames(1) = 'num_scalars_tend'
      r3Ptr(1) % dimNames(2) = 'nVertLevels'
      r3Ptr(1) % dimNames(3) = 'nCells'

      r3Ptr(1) % defaultValue = 0.0
      allocate(r3Ptr(1) % attLists(size(r3Ptr(1) % constituentNames, dim=1)))
      do index_counter = 1, size(r3Ptr(1) % constituentNames, dim=1)
         allocate(r3Ptr(1) % attLists(index_counter) % attList)
      end do
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qv', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Tendency of water vapor mass per unit volume '&
//'divided by d(zeta)/dz')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg m^{-3} s^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'tend_qv'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qc', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Tendency of cloud water mass per unit volume '&
//'divided by d(zeta)/dz')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg m^{-3} s^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'tend_qc'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qr', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Tendency of rain water mass per unit volume '&
//'divided by d(zeta)/dz')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg m^{-3} s^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'tend_qr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qi', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Tendency of ice mass per unit volume divided '&
//'by d(zeta)/dz')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg m^{-3} s^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'tend_qi'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qs', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Tendency of snow mass per unit volume divided '&
//'by d(zeta)/dz')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg m^{-3} s^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'tend_qs'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qg', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Tendency of graupel mass per unit volume '&
//'divided by d(zeta)/dz')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg m^{-3} s^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'tend_qg'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_ni', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Tendency of cloud ice number concentration '&
//'multiplied by dry air density divided by d(zeta)/dz')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'nb m^{-3} s^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'tend_ni'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_nr', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Tendency of rain number concentration '&
//'multiplied by dry air density divided by d(zeta)/dz')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'nb m^{-3} s^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'tend_nr'
      end if
      r3Ptr(1) % block => block

            r3Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'scalars_tend', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'scalars_tend', r3Ptr)

! Define var array pv_scalars_tend
      allocate(r3Ptr(1))
      index_counter = 0
      group_counter = -1
      group_start = -1
      group_started = .false.

! Starting group passive
! Define constituent var tend_pv_scalar
! My Packages are (null)
      index_counter = index_counter + 1
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'index_tend_pv_scalar', index_counter)
      end if
      group_counter = group_counter + 1
      if (.not. group_started) then
         group_start = index_counter
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'passive_start', group_start)
         end if
         group_started = .true.
      end if
      if (.not. group_started) then
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'passive_start', -1)
            call mpas_pool_add_dimension(newSubPool, 'passive_end', -1)
         end if
      else
         group_started = .false.
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'passive_end', index_counter)
         end if
      end if
! End of group       

      numConstituents = index_counter
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'num_pv_scalars_tend', numConstituents)
      end if
! Defining time level 1
      allocate( r3Ptr(1) % constituentNames(numConstituents) )
      r3Ptr(1) % fieldName = 'pv_scalars_tend'
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .true.
      r3Ptr(1) % isVarArray = .true.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_tend_pv_scalar', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'tend_pv_scalar'
      end if

! Setup dimensions for       
      r3Ptr(1) % dimNames(1) = 'num_pv_scalars_tend'
      r3Ptr(1) % dimNames(2) = 'nVertLevels'
      r3Ptr(1) % dimNames(3) = 'nCells'

      r3Ptr(1) % defaultValue = 0.0
      allocate(r3Ptr(1) % attLists(size(r3Ptr(1) % constituentNames, dim=1)))
      do index_counter = 1, size(r3Ptr(1) % constituentNames, dim=1)
         allocate(r3Ptr(1) % attLists(index_counter) % attList)
      end do
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_tend_pv_scalar', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Tendency of Ertel''s potential vorticity from '&
//'scalar advection multiplied by dry air density divided by d(zeta)/dz')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'PVU kg m^{-3} s^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'tend_pv_scalar'
      end if
      r3Ptr(1) % block => block

            r3Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'pv_scalars_tend', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'pv_scalars_tend', r3Ptr)

! Define variable tend_u
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_u'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Tendency of u from dynamics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'u', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_u', r2Ptr)

! Define variable tend_w
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_w'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Tendency of w from dynamics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_w', r2Ptr)

! Define variable tend_rho
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_rho'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Tendency of dry density from dynamics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rho_zz', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_rho', r2Ptr)

! Define variable tend_theta
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_theta'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg K m^{-3} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of coupled potential temperature rho*theta_m/zz '&
//'from dynamics and physics, updated each RK step')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'theta_m', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_theta', r2Ptr)

! Define variable rt_diabatic_tend
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rt_diabatic_tend'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Tendency of coupled potential temperature from physics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rt_diabatic_tend', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rt_diabatic_tend', r2Ptr)

! Define variable euler_tend_u
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'euler_tend_u'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Tendency of u from dynamics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'u_euler', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'euler_tend_u', r2Ptr)

! Define variable euler_tend_w
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'euler_tend_w'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Tendency of w from dynamics')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w_euler', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'euler_tend_w', r2Ptr)

! Define variable euler_tend_theta
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'euler_tend_theta'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg K m^{-3} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of coupled potential temperature rho*theta_m/zz '&
//'from dynamics and physics that does not change over a Runge-Kutta timestep (it excludes the flux divergence)')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'theta_euler', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'euler_tend_theta', r2Ptr)

! Define variable tend_w_pgf
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_w_pgf'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Tendency of w due to pressure gradient force')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w_pgf', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_w_pgf', r2Ptr)

! Define variable tend_w_buoy
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_w_buoy'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Tendency of w due to buoyancy force')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'w_buoy', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_w_buoy', r2Ptr)

! Define variable tend_sfc_pressure
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'tend_sfc_pressure'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'Pa s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Tendency of surface pressure')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'tend_sfc_pressure', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_sfc_pressure', r1Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
         call mpas_pool_add_config(newSubPool, 'is_periodic', block % domain % is_periodic)
         call mpas_pool_add_config(newSubPool, 'x_period', block % domain % x_period)
         call mpas_pool_add_config(newSubPool, 'y_period', block % domain % y_period)
      end if

   end subroutine atm_generate_pool_tend


   subroutine atm_generate_pool_lbc(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_io, only : MPAS_REAL_FILLVAL, MPAS_INT_FILLVAL, MPAS_CHAR_FILLVAL
      implicit none
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), dimension(:), pointer :: r0Ptr
      type (field1DReal), dimension(:), pointer :: r1Ptr
      type (field2DReal), dimension(:), pointer :: r2Ptr
      type (field3DReal), dimension(:), pointer :: r3Ptr
      type (field4DReal), dimension(:), pointer :: r4Ptr
      type (field5DReal), dimension(:), pointer :: r5Ptr
      type (field0DInteger), dimension(:), pointer :: i0Ptr
      type (field1DInteger), dimension(:), pointer :: i1Ptr
      type (field2DInteger), dimension(:), pointer :: i2Ptr
      type (field3DInteger), dimension(:), pointer :: i3Ptr
      type (field0DChar), dimension(:), pointer :: c0Ptr
      type (field1DChar), dimension(:), pointer :: c1Ptr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: mp_kessler_inActive
      logical, pointer :: mp_thompson_inActive
      logical, pointer :: mp_wsm6_inActive
      logical, pointer :: cu_grell_freitas_inActive
      logical, pointer :: cu_kain_fritsch_inActive
      logical, pointer :: cu_tiedtke_inActive
      logical, pointer :: bl_ysu_inActive
      logical, pointer :: bl_mynn_inActive
      logical, pointer :: iauActive
      logical, pointer :: limited_areaActive
      logical, pointer :: jedi_daActive
      logical, pointer :: tendenciesActive


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'mp_kessler_inActive', mp_kessler_inActive)
      call mpas_pool_get_package(packagePool, 'mp_thompson_inActive', mp_thompson_inActive)
      call mpas_pool_get_package(packagePool, 'mp_wsm6_inActive', mp_wsm6_inActive)
      call mpas_pool_get_package(packagePool, 'cu_grell_freitas_inActive', cu_grell_freitas_inActive)
      call mpas_pool_get_package(packagePool, 'cu_kain_fritsch_inActive', cu_kain_fritsch_inActive)
      call mpas_pool_get_package(packagePool, 'cu_tiedtke_inActive', cu_tiedtke_inActive)
      call mpas_pool_get_package(packagePool, 'bl_ysu_inActive', bl_ysu_inActive)
      call mpas_pool_get_package(packagePool, 'bl_mynn_inActive', bl_mynn_inActive)
      call mpas_pool_get_package(packagePool, 'iauActive', iauActive)
      call mpas_pool_get_package(packagePool, 'limited_areaActive', limited_areaActive)
      call mpas_pool_get_package(packagePool, 'jedi_daActive', jedi_daActive)
      call mpas_pool_get_package(packagePool, 'tendenciesActive', tendenciesActive)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'lbc', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'lbc', newSubPool)

! Define var array lbc_scalars
      allocate(r3Ptr(2))
      index_counter = 0
      group_counter = -1
      group_start = -1
      group_started = .false.

! Starting group moist
! Define constituent var lbc_qv
! My Packages are limited_area
      if (limited_areaActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qv', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
           if (associated(newSubPool)) then
              call mpas_pool_add_dimension(newSubPool, 'index_qv', -1)
           end if
      end if
! Define constituent var lbc_qc
! My packages are bl_mynn_in;cu_tiedtke_in;mp_kessler_in;mp_thompson_in;mp_wsm6_in
      if (bl_mynn_inActive .or. cu_tiedtke_inActive .or. mp_kessler_inActive .or. mp_thompson_inActive .or. mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qc', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qc', -1)
         end if
      end if
! Define constituent var lbc_qr
! My packages are mp_kessler_in;mp_thompson_in;mp_wsm6_in
      if (mp_kessler_inActive .or. mp_thompson_inActive .or. mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qr', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qr', -1)
         end if
      end if
! Define constituent var lbc_qi
! My packages are bl_mynn_in;cu_tiedtke_in;mp_thompson_in;mp_wsm6_in
      if (bl_mynn_inActive .or. cu_tiedtke_inActive .or. mp_thompson_inActive .or. mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qi', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qi', -1)
         end if
      end if
! Define constituent var lbc_qs
! My packages are mp_thompson_in;mp_wsm6_in
      if (mp_thompson_inActive .or. mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qs', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qs', -1)
         end if
      end if
! Define constituent var lbc_qg
! My packages are mp_thompson_in;mp_wsm6_in
      if (mp_thompson_inActive .or. mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qg', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qg', -1)
         end if
      end if
         if (.not. group_started) then
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', -1)
               call mpas_pool_add_dimension(newSubPool, 'moist_end', -1)
            end if
         else
            group_started = .false.
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_end', index_counter)
            end if
         end if
! End of group       
! Starting group number
! Define constituent var lbc_nr
! My Packages are mp_thompson_in
      if (mp_thompson_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_nr', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'number_start', group_start)
            end if
            group_started = .true.
         end if
      else
           if (associated(newSubPool)) then
              call mpas_pool_add_dimension(newSubPool, 'index_nr', -1)
           end if
      end if
! Define constituent var lbc_ni
! My packages are bl_mynn_in;mp_thompson_in
      if (bl_mynn_inActive .or. mp_thompson_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_ni', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'number_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_ni', -1)
         end if
      end if
         if (.not. group_started) then
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'number_start', -1)
               call mpas_pool_add_dimension(newSubPool, 'number_end', -1)
            end if
         else
            group_started = .false.
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'number_end', index_counter)
            end if
         end if
! End of group       

      numConstituents = index_counter
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'num_lbc_scalars', numConstituents)
      end if
! Defining time level 1
      allocate( r3Ptr(1) % constituentNames(numConstituents) )
      r3Ptr(1) % fieldName = 'lbc_scalars'
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .true.
      r3Ptr(1) % isVarArray = .true.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qv', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'lbc_qv'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qc', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'lbc_qc'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qr', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'lbc_qr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qi', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'lbc_qi'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qs', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'lbc_qs'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qg', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'lbc_qg'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_nr', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'lbc_nr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_ni', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'lbc_ni'
      end if

! Setup dimensions for       
      r3Ptr(1) % dimNames(1) = 'num_lbc_scalars'
      r3Ptr(1) % dimNames(2) = 'nVertLevels'
      r3Ptr(1) % dimNames(3) = 'nCells'

      r3Ptr(1) % defaultValue = 0.0
      allocate(r3Ptr(1) % attLists(size(r3Ptr(1) % constituentNames, dim=1)))
      do index_counter = 1, size(r3Ptr(1) % constituentNames, dim=1)
         allocate(r3Ptr(1) % attLists(index_counter) % attList)
      end do
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qv', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Lateral boundary tendency of water vapor '&
//'mixing ratio')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg kg^{-1} s^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'lbc_qv'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qc', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Lateral boundary tendency of cloud water '&
//'mixing ratio')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg kg^{-1} s^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'lbc_qc'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qr', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Lateral boundary tendency of rain water '&
//'mixing ratio')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg kg^{-1} s^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'lbc_qr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qi', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Lateral boundary tendency of ice mixing '&
//'ratio')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg kg^{-1} s^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'lbc_qi'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qs', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Lateral boundary tendency of snow mixing '&
//'ratio')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg kg^{-1} s^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'lbc_qs'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qg', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Lateral boundary tendency of graupel mixing '&
//'ratio')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg kg^{-1} s^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'lbc_qg'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_nr', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Lateral boundary tendency of rain number '&
//'concentration')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'm^{-3} s^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'lbc_nr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_ni', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Lateral boundary tendency of ice number '&
//'concentration')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'm^{-3} s^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'lbc_ni'
      end if
      r3Ptr(1) % block => block
! Defining time level 2
      allocate( r3Ptr(2) % constituentNames(numConstituents) )
      r3Ptr(2) % fieldName = 'lbc_scalars'
      r3Ptr(2) % isDecomposed = .true.
      r3Ptr(2) % hasTimeDimension = .true.
      r3Ptr(2) % isVarArray = .true.
      r3Ptr(2) % isPersistent = .true.
      r3Ptr(2) % isActive = .false.

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qv', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'lbc_qv'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qc', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'lbc_qc'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qr', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'lbc_qr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qi', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'lbc_qi'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qs', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'lbc_qs'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qg', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'lbc_qg'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_nr', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'lbc_nr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_ni', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(2) % constituentNames(const_index) = 'lbc_ni'
      end if

! Setup dimensions for       
      r3Ptr(2) % dimNames(1) = 'num_lbc_scalars'
      r3Ptr(2) % dimNames(2) = 'nVertLevels'
      r3Ptr(2) % dimNames(3) = 'nCells'

      r3Ptr(2) % defaultValue = 0.0
      allocate(r3Ptr(2) % attLists(size(r3Ptr(2) % constituentNames, dim=1)))
      do index_counter = 1, size(r3Ptr(2) % constituentNames, dim=1)
         allocate(r3Ptr(2) % attLists(index_counter) % attList)
      end do
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qv', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Lateral boundary tendency of water vapor '&
//'mixing ratio')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg kg^{-1} s^{-1}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'lbc_qv'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qc', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Lateral boundary tendency of cloud water '&
//'mixing ratio')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg kg^{-1} s^{-1}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'lbc_qc'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qr', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Lateral boundary tendency of rain water '&
//'mixing ratio')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg kg^{-1} s^{-1}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'lbc_qr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qi', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Lateral boundary tendency of ice mixing '&
//'ratio')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg kg^{-1} s^{-1}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'lbc_qi'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qs', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Lateral boundary tendency of snow mixing '&
//'ratio')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg kg^{-1} s^{-1}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'lbc_qs'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qg', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Lateral boundary tendency of graupel mixing '&
//'ratio')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'kg kg^{-1} s^{-1}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'lbc_qg'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_nr', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Lateral boundary tendency of rain number '&
//'concentration')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'm^{-3} s^{-1}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'lbc_nr'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_ni', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'long_name', 'Lateral boundary tendency of ice number '&
//'concentration')
         call mpas_add_att(r3Ptr(2) % attLists(const_index) % attList, 'units', 'm^{-3} s^{-1}')
         r3Ptr(2) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(2) % constituentNames(const_index) = 'lbc_ni'
      end if
      r3Ptr(2) % block => block

      if (limited_areaActive .or. bl_mynn_inActive .or. cu_tiedtke_inActive .or. mp_kessler_inActive .or. mp_thompson_inActive &
.or. mp_wsm6_inActive) then
         r3Ptr(1) % isActive = .true.
         r3Ptr(2) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'lbc_scalars', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'lbc_scalars', r3Ptr)

! Define variable lbc_u
      allocate(r2Ptr(2))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'lbc_u'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Lateral boundary tendency of u')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

! Setting up time level 2
      r2Ptr(2) % fieldName = 'lbc_u'
      r2Ptr(2) % isVarArray = .false.
      r2Ptr(2) % isDecomposed = .true.
      r2Ptr(2) % hasTimeDimension = .true.
      r2Ptr(2) % isPersistent = .true.
      r2Ptr(2) % isActive = .false.
! Setting up dimensions
      r2Ptr(2) % dimNames(1) = 'nVertLevels'
      r2Ptr(2) % dimNames(2) = 'nEdges'
      r2Ptr(2) % defaultValue = 0.0
      allocate(r2Ptr(2) % attLists(1))
      allocate(r2Ptr(2) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'long_name', 'Lateral boundary tendency of u')
      r2Ptr(2) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(2) % block => block

      if (limited_areaActive) then
         r2Ptr(1) % isActive = .true.
         r2Ptr(2) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'lbc_u', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'lbc_u', r2Ptr)

! Define variable lbc_ru
      allocate(r2Ptr(2))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'lbc_ru'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Lateral boundary tendency of rho_zz-coupled u')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

! Setting up time level 2
      r2Ptr(2) % fieldName = 'lbc_ru'
      r2Ptr(2) % isVarArray = .false.
      r2Ptr(2) % isDecomposed = .true.
      r2Ptr(2) % hasTimeDimension = .true.
      r2Ptr(2) % isPersistent = .true.
      r2Ptr(2) % isActive = .false.
! Setting up dimensions
      r2Ptr(2) % dimNames(1) = 'nVertLevels'
      r2Ptr(2) % dimNames(2) = 'nEdges'
      r2Ptr(2) % defaultValue = 0.0
      allocate(r2Ptr(2) % attLists(1))
      allocate(r2Ptr(2) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-2}')
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'long_name', 'Lateral boundary tendency of rho_zz-coupled u')
      r2Ptr(2) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(2) % block => block

      if (limited_areaActive) then
         r2Ptr(1) % isActive = .true.
         r2Ptr(2) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'lbc_ru', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'lbc_ru', r2Ptr)

! Define variable lbc_w
      allocate(r2Ptr(2))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'lbc_w'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Lateral boundary tendency of w')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

! Setting up time level 2
      r2Ptr(2) % fieldName = 'lbc_w'
      r2Ptr(2) % isVarArray = .false.
      r2Ptr(2) % isDecomposed = .true.
      r2Ptr(2) % hasTimeDimension = .true.
      r2Ptr(2) % isPersistent = .true.
      r2Ptr(2) % isActive = .false.
! Setting up dimensions
      r2Ptr(2) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(2) % dimNames(2) = 'nCells'
      r2Ptr(2) % defaultValue = 0.0
      allocate(r2Ptr(2) % attLists(1))
      allocate(r2Ptr(2) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'long_name', 'Lateral boundary tendency of w')
      r2Ptr(2) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(2) % block => block

      if (limited_areaActive) then
         r2Ptr(1) % isActive = .true.
         r2Ptr(2) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'lbc_w', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'lbc_w', r2Ptr)

! Define variable lbc_rho
      allocate(r2Ptr(2))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'lbc_rho'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Lateral boundary tendency of rho')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

! Setting up time level 2
      r2Ptr(2) % fieldName = 'lbc_rho'
      r2Ptr(2) % isVarArray = .false.
      r2Ptr(2) % isDecomposed = .true.
      r2Ptr(2) % hasTimeDimension = .true.
      r2Ptr(2) % isPersistent = .true.
      r2Ptr(2) % isActive = .false.
! Setting up dimensions
      r2Ptr(2) % dimNames(1) = 'nVertLevels'
      r2Ptr(2) % dimNames(2) = 'nCells'
      r2Ptr(2) % defaultValue = 0.0
      allocate(r2Ptr(2) % attLists(1))
      allocate(r2Ptr(2) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'units', 'kg m^{-3} s^{-1}')
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'long_name', 'Lateral boundary tendency of rho')
      r2Ptr(2) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(2) % block => block

      if (limited_areaActive) then
         r2Ptr(1) % isActive = .true.
         r2Ptr(2) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'lbc_rho', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'lbc_rho', r2Ptr)

! Define variable lbc_rho_zz
      allocate(r2Ptr(2))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'lbc_rho_zz'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Lateral boundary tendency of rho_zz')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

! Setting up time level 2
      r2Ptr(2) % fieldName = 'lbc_rho_zz'
      r2Ptr(2) % isVarArray = .false.
      r2Ptr(2) % isDecomposed = .true.
      r2Ptr(2) % hasTimeDimension = .true.
      r2Ptr(2) % isPersistent = .true.
      r2Ptr(2) % isActive = .false.
! Setting up dimensions
      r2Ptr(2) % dimNames(1) = 'nVertLevels'
      r2Ptr(2) % dimNames(2) = 'nCells'
      r2Ptr(2) % defaultValue = 0.0
      allocate(r2Ptr(2) % attLists(1))
      allocate(r2Ptr(2) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'units', 'kg m^{-3} s^{-1}')
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'long_name', 'Lateral boundary tendency of rho_zz')
      r2Ptr(2) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(2) % block => block

      if (limited_areaActive) then
         r2Ptr(1) % isActive = .true.
         r2Ptr(2) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'lbc_rho_zz', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'lbc_rho_zz', r2Ptr)

! Define variable lbc_rho_edge
      allocate(r2Ptr(2))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'lbc_rho_edge'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Lateral boundary tendency of rho_zz averaged from cell '&
//'centers to cell faces (edges)')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

! Setting up time level 2
      r2Ptr(2) % fieldName = 'lbc_rho_edge'
      r2Ptr(2) % isVarArray = .false.
      r2Ptr(2) % isDecomposed = .true.
      r2Ptr(2) % hasTimeDimension = .true.
      r2Ptr(2) % isPersistent = .true.
      r2Ptr(2) % isActive = .false.
! Setting up dimensions
      r2Ptr(2) % dimNames(1) = 'nVertLevels'
      r2Ptr(2) % dimNames(2) = 'nEdges'
      r2Ptr(2) % defaultValue = 0.0
      allocate(r2Ptr(2) % attLists(1))
      allocate(r2Ptr(2) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'units', 'kg m^{-3} s^{-1}')
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'long_name', 'Lateral boundary tendency of rho_zz averaged from cell '&
//'centers to cell faces (edges)')
      r2Ptr(2) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(2) % block => block

      if (limited_areaActive) then
         r2Ptr(1) % isActive = .true.
         r2Ptr(2) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'lbc_rho_edge', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'lbc_rho_edge', r2Ptr)

! Define variable lbc_theta
      allocate(r2Ptr(2))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'lbc_theta'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Lateral boundary tendency of theta')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

! Setting up time level 2
      r2Ptr(2) % fieldName = 'lbc_theta'
      r2Ptr(2) % isVarArray = .false.
      r2Ptr(2) % isDecomposed = .true.
      r2Ptr(2) % hasTimeDimension = .true.
      r2Ptr(2) % isPersistent = .true.
      r2Ptr(2) % isActive = .false.
! Setting up dimensions
      r2Ptr(2) % dimNames(1) = 'nVertLevels'
      r2Ptr(2) % dimNames(2) = 'nCells'
      r2Ptr(2) % defaultValue = 0.0
      allocate(r2Ptr(2) % attLists(1))
      allocate(r2Ptr(2) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'units', 'K s^{-2}')
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'long_name', 'Lateral boundary tendency of theta')
      r2Ptr(2) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(2) % block => block

      if (limited_areaActive) then
         r2Ptr(1) % isActive = .true.
         r2Ptr(2) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'lbc_theta', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'lbc_theta', r2Ptr)

! Define variable lbc_rtheta_m
      allocate(r2Ptr(2))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'lbc_rtheta_m'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg K m^{-3} s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Lateral boundary tendency of rho_zz-coupled theta_m')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

! Setting up time level 2
      r2Ptr(2) % fieldName = 'lbc_rtheta_m'
      r2Ptr(2) % isVarArray = .false.
      r2Ptr(2) % isDecomposed = .true.
      r2Ptr(2) % hasTimeDimension = .true.
      r2Ptr(2) % isPersistent = .true.
      r2Ptr(2) % isActive = .false.
! Setting up dimensions
      r2Ptr(2) % dimNames(1) = 'nVertLevels'
      r2Ptr(2) % dimNames(2) = 'nCells'
      r2Ptr(2) % defaultValue = 0.0
      allocate(r2Ptr(2) % attLists(1))
      allocate(r2Ptr(2) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'units', 'kg K m^{-3} s^{-2}')
      call mpas_add_att(r2Ptr(2) % attLists(1) % attList, 'long_name', 'Lateral boundary tendency of rho_zz-coupled theta_m')
      r2Ptr(2) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(2) % block => block

      if (limited_areaActive) then
         r2Ptr(1) % isActive = .true.
         r2Ptr(2) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'lbc_rtheta_m', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'lbc_rtheta_m', r2Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
         call mpas_pool_add_config(newSubPool, 'is_periodic', block % domain % is_periodic)
         call mpas_pool_add_config(newSubPool, 'x_period', block % domain % x_period)
         call mpas_pool_add_config(newSubPool, 'y_period', block % domain % y_period)
      end if

   end subroutine atm_generate_pool_lbc


   subroutine atm_generate_pool_diag_physics(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_io, only : MPAS_REAL_FILLVAL, MPAS_INT_FILLVAL, MPAS_CHAR_FILLVAL
      implicit none
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), dimension(:), pointer :: r0Ptr
      type (field1DReal), dimension(:), pointer :: r1Ptr
      type (field2DReal), dimension(:), pointer :: r2Ptr
      type (field3DReal), dimension(:), pointer :: r3Ptr
      type (field4DReal), dimension(:), pointer :: r4Ptr
      type (field5DReal), dimension(:), pointer :: r5Ptr
      type (field0DInteger), dimension(:), pointer :: i0Ptr
      type (field1DInteger), dimension(:), pointer :: i1Ptr
      type (field2DInteger), dimension(:), pointer :: i2Ptr
      type (field3DInteger), dimension(:), pointer :: i3Ptr
      type (field0DChar), dimension(:), pointer :: c0Ptr
      type (field1DChar), dimension(:), pointer :: c1Ptr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: mp_kessler_inActive
      logical, pointer :: mp_thompson_inActive
      logical, pointer :: mp_wsm6_inActive
      logical, pointer :: cu_grell_freitas_inActive
      logical, pointer :: cu_kain_fritsch_inActive
      logical, pointer :: cu_tiedtke_inActive
      logical, pointer :: bl_ysu_inActive
      logical, pointer :: bl_mynn_inActive
      logical, pointer :: iauActive
      logical, pointer :: limited_areaActive
      logical, pointer :: jedi_daActive
      logical, pointer :: tendenciesActive


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'mp_kessler_inActive', mp_kessler_inActive)
      call mpas_pool_get_package(packagePool, 'mp_thompson_inActive', mp_thompson_inActive)
      call mpas_pool_get_package(packagePool, 'mp_wsm6_inActive', mp_wsm6_inActive)
      call mpas_pool_get_package(packagePool, 'cu_grell_freitas_inActive', cu_grell_freitas_inActive)
      call mpas_pool_get_package(packagePool, 'cu_kain_fritsch_inActive', cu_kain_fritsch_inActive)
      call mpas_pool_get_package(packagePool, 'cu_tiedtke_inActive', cu_tiedtke_inActive)
      call mpas_pool_get_package(packagePool, 'bl_ysu_inActive', bl_ysu_inActive)
      call mpas_pool_get_package(packagePool, 'bl_mynn_inActive', bl_mynn_inActive)
      call mpas_pool_get_package(packagePool, 'iauActive', iauActive)
      call mpas_pool_get_package(packagePool, 'limited_areaActive', limited_areaActive)
      call mpas_pool_get_package(packagePool, 'jedi_daActive', jedi_daActive)
      call mpas_pool_get_package(packagePool, 'tendenciesActive', tendenciesActive)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'diag_physics', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'diag_physics', newSubPool)

! Define variable nsteps_accum
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'nsteps_accum'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'number of accumulated time-steps in a day')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'nsteps_accum', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'nsteps_accum', r1Ptr)

! Define variable ndays_accum
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'ndays_accum'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'number of accumulated days in a year')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ndays_accum', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'ndays_accum', r1Ptr)

! Define variable tlag
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tlag'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nLags'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'daily mean surface temperature of prior days')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'tlag', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tlag', r2Ptr)

! Define variable tday_accum
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'tday_accum'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated daily surface temperature for current day')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'tday_accum', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'tday_accum', r1Ptr)

! Define variable tyear_mean
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'tyear_mean'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'annual mean surface temperature')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'tyear_mean', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'tyear_mean', r1Ptr)

! Define variable tyear_accum
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'tyear_accum'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated yearly surface temperature for current year')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'tyear_accum', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'tyear_accum', r1Ptr)

! Define variable refl10cm_3D
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'refl10cm_3D'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'dBZ')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', '10 cm radar reflectivity - 3D')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (mp_thompson_inActive .or. mp_wsm6_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'refl10cm_3D', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'refl10cm_3D', r2Ptr)

! Define variable refl10cm_max
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'refl10cm_max'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'dBZ')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', '10 cm maximum radar reflectivity')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (mp_thompson_inActive .or. mp_wsm6_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'refl10cm_max', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'refl10cm_max', r1Ptr)

! Define variable refl10cm_1km
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'refl10cm_1km'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'dBZ')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'diagnosed 10 cm radar reflectivity at 1 km AGL')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (mp_thompson_inActive .or. mp_wsm6_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'refl10cm_1km', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'refl10cm_1km', r1Ptr)

! Define variable refl10cm_1km_max
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'refl10cm_1km_max'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'dBZ')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'maximum diagnosed 10 cm radar reflectivity at 1 km AGL '&
//'since last output time')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (mp_thompson_inActive .or. mp_wsm6_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'refl10cm_1km_max', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'refl10cm_1km_max', r1Ptr)

! Define variable i_rainnc
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'i_rainnc'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'incidence of accumulated grid-scale precipitation greater '&
//'than config_bucket_rainnc')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      if (mp_kessler_inActive .or. mp_thompson_inActive .or. mp_wsm6_inActive) then
         i1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'i_rainnc', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'i_rainnc', i1Ptr)

! Define variable rainncv
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'rainncv'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'mm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'time-step total grid-scale precipitation')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (mp_kessler_inActive .or. mp_thompson_inActive .or. mp_wsm6_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rainncv', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'rainncv', r1Ptr)

! Define variable snowncv
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'snowncv'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'mm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'time-step grid-scale precipitation of snow')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (mp_thompson_inActive .or. mp_wsm6_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'snowncv', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'snowncv', r1Ptr)

! Define variable graupelncv
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'graupelncv'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'mm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'time-step grid-scale precipitation of graupel')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (mp_thompson_inActive .or. mp_wsm6_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'graupelncv', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'graupelncv', r1Ptr)

! Define variable rainnc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'rainnc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'mm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated total grid-scale precipitation')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (mp_kessler_inActive .or. mp_thompson_inActive .or. mp_wsm6_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rainnc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'rainnc', r1Ptr)

! Define variable snownc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'snownc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'mm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated grid-scale precipitation of snow')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (mp_thompson_inActive .or. mp_wsm6_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'snownc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'snownc', r1Ptr)

! Define variable graupelnc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'graupelnc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'mm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated grid-scale precipitation of graupel')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (mp_thompson_inActive .or. mp_wsm6_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'graupelnc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'graupelnc', r1Ptr)

! Define variable rainprod
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rainprod'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 's^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'rain production rate')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (mp_thompson_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rainprod', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rainprod', r2Ptr)

! Define variable evapprod
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'evapprod'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 's^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'rain evaporation rate')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (mp_thompson_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'evapprod', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'evapprod', r2Ptr)

! Define variable re_cloud
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 're_cloud'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'effective radius of cloud water droplets')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (mp_thompson_inActive .or. mp_wsm6_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 're_cloud', r2Ptr)
      call mpas_pool_add_field(block % allFields, 're_cloud', r2Ptr)

! Define variable re_ice
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 're_ice'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'effective radius of cloud ice crystals')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (mp_thompson_inActive .or. mp_wsm6_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 're_ice', r2Ptr)
      call mpas_pool_add_field(block % allFields, 're_ice', r2Ptr)

! Define variable re_snow
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 're_snow'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'effective radius of snow crystals')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (mp_thompson_inActive .or. mp_wsm6_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 're_snow', r2Ptr)
      call mpas_pool_add_field(block % allFields, 're_snow', r2Ptr)

! Define variable nt_c
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'nt_c'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'nb kg^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'one-moment constant cloud water droplet concentration')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (mp_thompson_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'nt_c', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'nt_c', r1Ptr)

! Define variable mu_c
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'mu_c'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'gamma shape parameter')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (mp_thompson_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'mu_c', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'mu_c', r1Ptr)

! Define variable sr
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'sr'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'time-step ratio of frozen versus total grid-scale '&
//'precipitation')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'sr', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'sr', r1Ptr)

! Define variable precipw
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'precipw'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'precipitable water')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'precipw', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'precipw', r1Ptr)

! Define variable tend_theta_mp_evap_cw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_theta_mp_evap_cw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Net potential temperature heating rate from cloud water '&
//'condensation and evaporation')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (mp_thompson_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'tend_theta_mp_evap_cw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_theta_mp_evap_cw', r2Ptr)

! Define variable tend_theta_mp_evap_rw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_theta_mp_evap_rw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Potential temperature heating rate from rain water '&
//'evaporation')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (mp_thompson_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'tend_theta_mp_evap_rw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_theta_mp_evap_rw', r2Ptr)

! Define variable tend_theta_mp_depo_ice
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_theta_mp_depo_ice'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Net potential temperature heating rate from deposition and '&
//'sublimation of all ice hydrometeors')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (mp_thompson_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'tend_theta_mp_depo_ice', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_theta_mp_depo_ice', r2Ptr)

! Define variable tend_theta_mp_melt_ice
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_theta_mp_melt_ice'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Potential temperature heating rate from melting of all ice '&
//'hydrometeors')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (mp_thompson_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'tend_theta_mp_melt_ice', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_theta_mp_melt_ice', r2Ptr)

! Define variable tend_theta_mp_frez_ice
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tend_theta_mp_frez_ice'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Potential temperature heating rate from freezing/riming of '&
//'all ice hydrometeors')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (mp_thompson_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'tend_theta_mp_frez_ice', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tend_theta_mp_frez_ice', r2Ptr)

! Define variable i_rainc
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'i_rainc'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'incidence of accumulated convective precipitation greater '&
//'than config_bucket_rainc')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      if (cu_grell_freitas_inActive .or. cu_kain_fritsch_inActive .or. cu_tiedtke_inActive) then
         i1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'i_rainc', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'i_rainc', i1Ptr)

! Define variable cuprec
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'cuprec'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'mm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'convective precipitation rate')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (cu_grell_freitas_inActive .or. cu_kain_fritsch_inActive .or. cu_tiedtke_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'cuprec', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'cuprec', r1Ptr)

! Define variable rainc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'rainc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'mm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated convective precipitation')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (cu_grell_freitas_inActive .or. cu_kain_fritsch_inActive .or. cu_tiedtke_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rainc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'rainc', r1Ptr)

! Define variable raincv
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'raincv'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'mm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'time-step convective precipitation')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (cu_grell_freitas_inActive .or. cu_kain_fritsch_inActive .or. cu_tiedtke_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'raincv', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'raincv', r1Ptr)

! Define variable nca
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'nca'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 's')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'relaxation time for KF parameterization of convection')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (cu_kain_fritsch_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'nca', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'nca', r1Ptr)

! Define variable w0avg
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'w0avg'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'time running-averaged vertical velocity')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (cu_kain_fritsch_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'w0avg', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'w0avg', r2Ptr)

! Define variable k22_shallow
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'k22_shallow'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'index of convection at k22 level in Grell-Freitas '&
//'convection scheme')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      if (cu_grell_freitas_inActive) then
         i1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'k22_shallow', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'k22_shallow', i1Ptr)

! Define variable kbcon_shallow
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'kbcon_shallow'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'index of lowest level of shallow convection in '&
//'Grell-Freitas convection scheme')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      if (cu_grell_freitas_inActive) then
         i1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'kbcon_shallow', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'kbcon_shallow', i1Ptr)

! Define variable ktop_shallow
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'ktop_shallow'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'index of highest level of shallow convection in '&
//'Grell-Freitas convection scheme')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      if (cu_grell_freitas_inActive) then
         i1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'ktop_shallow', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'ktop_shallow', i1Ptr)

! Define variable ktop_deep
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'ktop_deep'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'index of highest level of deep convection in Grell-Freitas '&
//'convection scheme')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      if (cu_grell_freitas_inActive) then
         i1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'ktop_deep', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'ktop_deep', i1Ptr)

! Define variable xmb_total
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'xmb_total'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'cloud-base mass flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (cu_grell_freitas_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'xmb_total', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'xmb_total', r1Ptr)

! Define variable xmb_shallow
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'xmb_shallow'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'cloud base mass flux for shallow convection')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (cu_grell_freitas_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'xmb_shallow', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'xmb_shallow', r1Ptr)

! Define variable qc_cu
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'qc_cu'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'in-cloud cloud water mixing ratio in Grell-Freitas cloud '&
//'model')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (cu_grell_freitas_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'qc_cu', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'qc_cu', r2Ptr)

! Define variable qi_cu
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'qi_cu'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'in-cloud cloud ice mixing ratio in Grell-Freitas cloud '&
//'model')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (cu_grell_freitas_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'qi_cu', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'qi_cu', r2Ptr)

! Define variable cubot
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'cubot'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'index of highest level of convection')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (cu_grell_freitas_inActive .or. cu_kain_fritsch_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'cubot', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'cubot', r1Ptr)

! Define variable cutop
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'cutop'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'index of lowest level of convection')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (cu_grell_freitas_inActive .or. cu_kain_fritsch_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'cutop', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'cutop', r1Ptr)

! Define variable kpbl
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'kpbl'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'index level of PBL top')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         i1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'kpbl', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'kpbl', i1Ptr)

! Define variable hpbl
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'hpbl'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Planetary Boundary Layer (PBL) height')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'hpbl', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'hpbl', r1Ptr)

! Define variable delta
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'delta'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'entrainment layer depth from PBL scheme')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'delta', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'delta', r1Ptr)

! Define variable wstar
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'wstar'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'mixed velocity scale from PBL scheme')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'wstar', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'wstar', r1Ptr)

! Define variable kzh
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'kzh'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm^{2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'vertical diffusion coefficient of potential temperature')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'kzh', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'kzh', r2Ptr)

! Define variable kzm
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'kzm'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm^{2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'vertical diffusion coefficient of mommentum')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'kzm', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'kzm', r2Ptr)

! Define variable kzq
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'kzq'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm^{2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'vertical diffusion coefficient of water vapor and cloud '&
//'condensates')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'kzq', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'kzq', r2Ptr)

! Define variable exch_h
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'exch_h'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm^{2} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'exchange coefficient')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_ysu_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'exch_h', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'exch_h', r2Ptr)

! Define variable cov
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'cov'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K kg kg^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'liquid water - liquid water potential temperature '&
//'covariance')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'cov', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'cov', r2Ptr)

! Define variable el_pbl
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'el_pbl'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'mixing length from PBL scheme')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'el_pbl', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'el_pbl', r2Ptr)

! Define variable qke
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'qke'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm^{2} s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'twice turbulent kinetic energy')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'qke', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'qke', r2Ptr)

! Define variable qke_adv
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'qke_adv'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm^{2} s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'twice turbulent kinetic energy')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'qke_adv', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'qke_adv', r2Ptr)

! Define variable qsq
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'qsq'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', '{kg kg^{-1}}^2')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'liquid water variance')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'qsq', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'qsq', r2Ptr)

! Define variable tsq
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tsq'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K^{2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'liquid water potential temperature variance')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'tsq', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tsq', r2Ptr)

! Define variable tke_pbl
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tke_pbl'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm^{2} s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'turbulent kinetic energy from PBL')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'tke_pbl', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tke_pbl', r2Ptr)

! Define variable dqke
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dqke'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm^{2} s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'TKE change')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'dqke', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dqke', r2Ptr)

! Define variable qbuoy
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'qbuoy'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm^{2} s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'TKE production - buoyancy')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'qbuoy', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'qbuoy', r2Ptr)

! Define variable qdiss
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'qdiss'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm^{2} s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'TKE dissipation')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'qdiss', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'qdiss', r2Ptr)

! Define variable qshear
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'qshear'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm^{2} s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'TKE production - shear')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'qshear', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'qshear', r2Ptr)

! Define variable qwt
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'qwt'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm^{2} s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'TKE vertical distribution')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'qwt', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'qwt', r2Ptr)

! Define variable hfx
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'hfx'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'upward heat flux at the surface')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'hfx', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'hfx', r1Ptr)

! Define variable mavail
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'mavail'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'surface moisture availability (between 0 and 1)')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'mavail', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'mavail', r1Ptr)

! Define variable mol
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'mol'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'T* in similarity theory')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'mol', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'mol', r1Ptr)

! Define variable qfx
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'qfx'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'upward moisture flux at the surface')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'qfx', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'qfx', r1Ptr)

! Define variable qsfc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'qsfc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'specific humidity at lower boundary')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'qsfc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'qsfc', r1Ptr)

! Define variable ust
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'ust'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'U* in similarity theory')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'ust', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'ust', r1Ptr)

! Define variable ustm
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'ustm'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'U* in similarity theory without vconv')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'ustm', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'ustm', r1Ptr)

! Define variable zol
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'zol'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'z/L height over Monin-Obukhov length')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'zol', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'zol', r1Ptr)

! Define variable br
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'br'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Richardson number')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'br', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'br', r1Ptr)

! Define variable cd
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'cd'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'drag coefficient at 10-meter')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'cd', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'cd', r1Ptr)

! Define variable cda
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'cda'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'drag coefficient at lowest model level')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'cda', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'cda', r1Ptr)

! Define variable chs
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'chs'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'surface exchange coefficient for heat and moisture')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'chs', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'chs', r1Ptr)

! Define variable chs2
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'chs2'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'surface exchange coefficient for heat at 2-meter')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'chs2', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'chs2', r1Ptr)

! Define variable cqs2
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'cqs2'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'surface exchange coefficient for moisture at 2-meter')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'cqs2', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'cqs2', r1Ptr)

! Define variable ck
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'ck'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'enthalpy exchange coeff at 10-meter')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'ck', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'ck', r1Ptr)

! Define variable cka
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'cka'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'enthalpy exchange coefficient at lowest model level')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'cka', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'cka', r1Ptr)

! Define variable cpm
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'cpm'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'J K^{-1} kg^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'specific heat of dry air at constant pressure at lowest '&
//'model level')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'cpm', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'cpm', r1Ptr)

! Define variable flhc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'flhc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2} K^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'exchange coefficient for heat')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'flhc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'flhc', r1Ptr)

! Define variable flqc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'flqc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2} s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'exchange coefficient for moisture')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'flqc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'flqc', r1Ptr)

! Define variable gz1oz0
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'gz1oz0'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'log(z/z0) where z0 is roughness length')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'gz1oz0', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'gz1oz0', r1Ptr)

! Define variable lh
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'lh'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'latent heat flux at the surface')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'lh', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lh', r1Ptr)

! Define variable psim
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'psim'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'similarity stability function for momentum')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'psim', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'psim', r1Ptr)

! Define variable psih
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'psih'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'similarity stability function for heat')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'psih', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'psih', r1Ptr)

! Define variable qgh
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'qgh'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'lowest level saturation mixing ratio')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'qgh', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'qgh', r1Ptr)

! Define variable regime
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'regime'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'flag indicating the PBL regime (stable,unstable,...)')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'regime', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'regime', r1Ptr)

! Define variable rmol
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'rmol'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', '1./L Monin Obukhov length')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rmol', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'rmol', r1Ptr)

! Define variable wspd
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'wspd'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'wind speed at lowest model level')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'wspd', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'wspd', r1Ptr)

! Define variable fh
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'fh'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'integrated stability function for heat')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'fh', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'fh', r1Ptr)

! Define variable fm
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'fm'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'integrated stability function for moisture')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'fm', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'fm', r1Ptr)

! Define variable ch
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'ch'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'surface exchange coefficient for heat')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'ch', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'ch', r1Ptr)

! Define variable qcg
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'qcg'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'cloud water mixing ratio at the ground surface')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'qcg', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'qcg', r1Ptr)

! Define variable sh3d
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'sh3d'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'stability function for heat')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'sh3d', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'sh3d', r2Ptr)

! Define variable u10
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'u10'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', '10-meter zonal wind')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'u10', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'u10', r1Ptr)

! Define variable v10
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'v10'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', '10-meter meridional wind')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'v10', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'v10', r1Ptr)

! Define variable q2
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'q2'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', '2-meter specific humidity')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'q2', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'q2', r1Ptr)

! Define variable t2m
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 't2m'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', '2-meter temperature')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 't2m', r1Ptr)
      call mpas_pool_add_field(block % allFields, 't2m', r1Ptr)

! Define variable th2m
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'th2m'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', '2-meter potential temperature')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r1Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'th2m', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'th2m', r1Ptr)

! Define variable dusfcg
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'dusfcg'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'Pa m s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'vertically-integrated gravity wave drag over orography '&
//'u-stress')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dusfcg', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'dusfcg', r1Ptr)

! Define variable dvsfcg
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'dvsfcg'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'Pa m s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'vertically-integrated gravity wave drag over orography '&
//'v-stress')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dvsfcg', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'dvsfcg', r1Ptr)

! Define variable dtaux3d
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dtaux3d'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'gravity wave drag over orography u-stress')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dtaux3d', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dtaux3d', r2Ptr)

! Define variable dtauy3d
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dtauy3d'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'gravity wave drag over orography v-stress')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dtauy3d', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dtauy3d', r2Ptr)

! Define variable rubldiff
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rubldiff'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'change in PBL zonal wind tendency due to gravity wave drag '&
//'over orography')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rubldiff', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rubldiff', r2Ptr)

! Define variable rvbldiff
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rvbldiff'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-2}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'change in PBL meridional wind tendency due to gravity wave '&
//'drag over orography')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rvbldiff', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rvbldiff', r2Ptr)

! Define variable i_acswdnb
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'i_acswdnb'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'incidence of accumulated all-sky downward surface '&
//'shortwave radiation greater than config_bucket_radt')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'i_acswdnb', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'i_acswdnb', i1Ptr)

! Define variable i_acswdnbc
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'i_acswdnbc'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'incidence of accumulated clear-sky downward surface '&
//'shortwave radiation greater than config_bucket_radt')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'i_acswdnbc', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'i_acswdnbc', i1Ptr)

! Define variable i_acswdnt
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'i_acswdnt'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'incidence of accumulated all-sky downward '&
//'top-of-atmosphere shortwave radiation greater than config_bucket_radt')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'i_acswdnt', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'i_acswdnt', i1Ptr)

! Define variable i_acswdntc
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'i_acswdntc'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'incidence of accumulated clear-sky downward '&
//'top-of-atmosphere shortwave radiation greater than config_bucket_radt')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'i_acswdntc', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'i_acswdntc', i1Ptr)

! Define variable i_acswupb
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'i_acswupb'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'incidence of accumulated all-sky upward surface shortwave '&
//'radiation greater than config_bucket_radt')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'i_acswupb', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'i_acswupb', i1Ptr)

! Define variable i_acswupbc
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'i_acswupbc'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'incidence of accumulated clear-sky upward surface '&
//'shortwave radiation greater than config_bucket_radt')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'i_acswupbc', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'i_acswupbc', i1Ptr)

! Define variable i_acswupt
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'i_acswupt'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'incidence of accumulated all-sky upward top-of-atmosphere '&
//'shortwave radiation greater than config_bucket_radt')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'i_acswupt', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'i_acswupt', i1Ptr)

! Define variable i_acswuptc
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'i_acswuptc'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'incidence of accumulated clear-sky upward '&
//'top-of-atmosphere shortwave radiation greater than config_bucket_radt')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'i_acswuptc', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'i_acswuptc', i1Ptr)

! Define variable coszr
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'coszr'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'cosine of zenith solar angle')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'coszr', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'coszr', r1Ptr)

! Define variable swcf
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'swcf'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'top-of-atmosphere cloud shortwave radiative forcing')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'swcf', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'swcf', r1Ptr)

! Define variable swdnb
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'swdnb'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'all-sky downward surface shortwave radiation flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'swdnb', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'swdnb', r1Ptr)

! Define variable swdnbc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'swdnbc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'clear-sky downward surface shortwave radiation flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'swdnbc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'swdnbc', r1Ptr)

! Define variable swdnt
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'swdnt'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'all-sky downward top-of-atmosphere shortwave radiation '&
//'flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'swdnt', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'swdnt', r1Ptr)

! Define variable swdntc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'swdntc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'clear-sky downward top-of-atmosphere shortwave radiation '&
//'flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'swdntc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'swdntc', r1Ptr)

! Define variable swupb
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'swupb'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'all-sky upward surface shortwave radiation flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'swupb', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'swupb', r1Ptr)

! Define variable swupbc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'swupbc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'clear-sky upward surface shortwave radiation flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'swupbc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'swupbc', r1Ptr)

! Define variable swupt
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'swupt'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'all-sky upward top-of-atmosphere shortwave radiation flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'swupt', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'swupt', r1Ptr)

! Define variable swuptc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'swuptc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'clear-sky upward top-of-atmosphere shortwave radiation '&
//'flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'swuptc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'swuptc', r1Ptr)

! Define variable acswdnb
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'acswdnb'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated all-sky downward surface shortwave radiation '&
//'flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acswdnb', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'acswdnb', r1Ptr)

! Define variable acswdnbc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'acswdnbc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated clear-sky downward surface shortwave radiation '&
//'flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acswdnbc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'acswdnbc', r1Ptr)

! Define variable acswdnt
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'acswdnt'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated all-sky downward top-of-atmosphere shortwave '&
//'radiation flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acswdnt', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'acswdnt', r1Ptr)

! Define variable acswdntc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'acswdntc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated clear-sky downward top-of-atmosphere shortwave '&
//'radiation flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acswdntc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'acswdntc', r1Ptr)

! Define variable acswupb
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'acswupb'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated all-sky upward surface shortwave radiation '&
//'flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acswupb', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'acswupb', r1Ptr)

! Define variable acswupbc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'acswupbc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated clear-sky upward surface shortwave radiation '&
//'flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acswupbc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'acswupbc', r1Ptr)

! Define variable acswupt
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'acswupt'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated all-sky upward top-of-atmosphere shortwave '&
//'radiation flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acswupt', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'acswupt', r1Ptr)

! Define variable acswuptc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'acswuptc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated clear-sky upward top-of-atmosphere shortwave '&
//'radiation flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acswuptc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'acswuptc', r1Ptr)

! Define variable gsw
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'gsw'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'net surface shortwave radiation flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'gsw', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'gsw', r1Ptr)

! Define variable swdnflx
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'swdnflx'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP2'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', '-')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'swdnflx', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'swdnflx', r2Ptr)

! Define variable swdnflxc
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'swdnflxc'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP2'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', '-')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'swdnflxc', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'swdnflxc', r2Ptr)

! Define variable swupflx
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'swupflx'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP2'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', '-')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'swupflx', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'swupflx', r2Ptr)

! Define variable swupflxc
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'swupflxc'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP2'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', '-')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'swupflxc', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'swupflxc', r2Ptr)

! Define variable nlrad
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'nlrad'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'number of layers added above the model-top in the RRTMG lw '&
//'radiation code')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'nlrad', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'nlrad', i1Ptr)

! Define variable plrad
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'plrad'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'Pa')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'pressure at model-top')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'plrad', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'plrad', r1Ptr)

! Define variable i_aclwdnb
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'i_aclwdnb'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'incidence of accumulated all-sky downward surface longwave '&
//'radiation greater than config_bucket_radt')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'i_aclwdnb', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'i_aclwdnb', i1Ptr)

! Define variable i_aclwdnbc
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'i_aclwdnbc'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'incidence of accumulated clear-sky downward surface '&
//'longwave radiation greater than config_bucket_radt')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'i_aclwdnbc', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'i_aclwdnbc', i1Ptr)

! Define variable i_aclwdnt
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'i_aclwdnt'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'incidence of accumulated all-sky downward '&
//'top-of-atmosphere longwave radiation greater than config_bucket_radt')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'i_aclwdnt', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'i_aclwdnt', i1Ptr)

! Define variable i_aclwdntc
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'i_aclwdntc'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'incidence of accumulated clear-sky downward '&
//'top-of-atmosphere longwave radiation greater than config_bucket_radt')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'i_aclwdntc', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'i_aclwdntc', i1Ptr)

! Define variable i_aclwupb
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'i_aclwupb'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'incidence of accumulated all-sky upward surface longwave '&
//'radiation greater than config_bucket_radt')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'i_aclwupb', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'i_aclwupb', i1Ptr)

! Define variable i_aclwupbc
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'i_aclwupbc'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'incidence of accumulated clear-sky upward surface longwave '&
//'radiation greater than config_bucket_radt')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'i_aclwupbc', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'i_aclwupbc', i1Ptr)

! Define variable i_aclwupt
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'i_aclwupt'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'incidence of accumulated all-sky upward top-of-atmosphere '&
//'longwave radiation greater than config_bucket_radt')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'i_aclwupt', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'i_aclwupt', i1Ptr)

! Define variable i_aclwuptc
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'i_aclwuptc'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .true.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'incidence of accumulated clear-sky upward '&
//'top-of-atmosphere longwave radiation greater than config_bucket_radt')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'i_aclwuptc', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'i_aclwuptc', i1Ptr)

! Define variable lwcf
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'lwcf'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'top-of-atmosphere cloud longwave radiative forcing')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lwcf', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lwcf', r1Ptr)

! Define variable lwdnb
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'lwdnb'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'all-sky downward surface longwave radiation flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lwdnb', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lwdnb', r1Ptr)

! Define variable lwdnbc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'lwdnbc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'clear-sky downward surface longwave radiation flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lwdnbc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lwdnbc', r1Ptr)

! Define variable lwdnt
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'lwdnt'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'all-sky downward top-of-the-atmosphere longwave radiation '&
//'flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lwdnt', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lwdnt', r1Ptr)

! Define variable lwdntc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'lwdntc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'clear-sky downward top-of-the-atmosphere longwave '&
//'radiation flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lwdntc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lwdntc', r1Ptr)

! Define variable lwupb
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'lwupb'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'all-sky upward surface longwave radiation flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lwupb', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lwupb', r1Ptr)

! Define variable lwupbc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'lwupbc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'clear-sky upward surface longwave radiation flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lwupbc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lwupbc', r1Ptr)

! Define variable lwupt
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'lwupt'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'all-sky upward top-of-the-atmosphere longwave radiation '&
//'flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lwupt', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lwupt', r1Ptr)

! Define variable lwuptc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'lwuptc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'clear-sky upward top-of-the-atmosphere longwave radiation '&
//'flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lwuptc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lwuptc', r1Ptr)

! Define variable aclwdnb
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'aclwdnb'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated all-sky downward surface longwave radiation '&
//'flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'aclwdnb', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'aclwdnb', r1Ptr)

! Define variable aclwdnbc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'aclwdnbc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated clear-sky downward surface longwave radiation '&
//'flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'aclwdnbc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'aclwdnbc', r1Ptr)

! Define variable aclwdnt
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'aclwdnt'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated all-sky downward top-of-the-atmosphere '&
//'longwave radiation flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'aclwdnt', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'aclwdnt', r1Ptr)

! Define variable aclwdntc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'aclwdntc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated clear-sky downward top-of-the-atmosphere '&
//'longwave radiation flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'aclwdntc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'aclwdntc', r1Ptr)

! Define variable aclwupb
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'aclwupb'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated all-sky upward surface longwave radiation '&
//'flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'aclwupb', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'aclwupb', r1Ptr)

! Define variable aclwupbc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'aclwupbc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated clear-sky upward surface longwave radiation '&
//'flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'aclwupbc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'aclwupbc', r1Ptr)

! Define variable aclwupt
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'aclwupt'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated all-sky upward top-of-the-atmosphere longwave '&
//'radiation flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'aclwupt', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'aclwupt', r1Ptr)

! Define variable aclwuptc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'aclwuptc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated clear-sky upward top-of-the-atmosphere '&
//'longwave radiation flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'aclwuptc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'aclwuptc', r1Ptr)

! Define variable olrtoa
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'olrtoa'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'all-sky top-of-atmosphere outgoing longwave radiation '&
//'flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'olrtoa', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'olrtoa', r1Ptr)

! Define variable glw
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'glw'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'all-sky downward surface longwave radiation')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'glw', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'glw', r1Ptr)

! Define variable rre_cloud
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rre_cloud'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'microns')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'effective radius of cloud water droplets calculated in '&
//'RRTMG radiation')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rre_cloud', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rre_cloud', r2Ptr)

! Define variable rre_ice
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rre_ice'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'microns')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'effective radius of cloud ice crystals calculated in RRTMG '&
//'radiation')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rre_ice', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rre_ice', r2Ptr)

! Define variable rre_snow
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rre_snow'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'microns')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'effective radius of snow crystals calculated in RRTMG '&
//'radiation')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rre_snow', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rre_snow', r2Ptr)

! Define variable m_hybi
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'm_hybi'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nAerLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Matched hybi (needs to be re-checked)')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'm_hybi', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'm_hybi', r2Ptr)

! Define variable m_ps
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'm_ps'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'Pa')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'Surface pressure from match on MPAS grid')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'm_ps', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'm_ps', r1Ptr)

! Define variable absnxt
      allocate(r3Ptr(1))

! Setting up time level 1
      r3Ptr(1) % fieldName = 'absnxt'
      r3Ptr(1) % isVarArray = .false.
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .true.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.
! Setting up dimensions
      r3Ptr(1) % dimNames(1) = 'nVertLevels'
      r3Ptr(1) % dimNames(2) = 'cam_dim1'
      r3Ptr(1) % dimNames(3) = 'nCells'
      r3Ptr(1) % defaultValue = 0.0
      allocate(r3Ptr(1) % attLists(1))
      allocate(r3Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'long_name', 'Total nearest layer absorptivity')
      r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r3Ptr(1) % block => block

      r3Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'absnxt', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'absnxt', r3Ptr)

! Define variable abstot
      allocate(r3Ptr(1))

! Setting up time level 1
      r3Ptr(1) % fieldName = 'abstot'
      r3Ptr(1) % isVarArray = .false.
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .true.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.
! Setting up dimensions
      r3Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r3Ptr(1) % dimNames(2) = 'nVertLevelsP1'
      r3Ptr(1) % dimNames(3) = 'nCells'
      r3Ptr(1) % defaultValue = 0.0
      allocate(r3Ptr(1) % attLists(1))
      allocate(r3Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'long_name', 'Total absorptivity')
      r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r3Ptr(1) % block => block

      r3Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'abstot', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'abstot', r3Ptr)

! Define variable emstot
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'emstot'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevelsP1'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', '-')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Total emissivity')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'emstot', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'emstot', r2Ptr)

! Define variable o3clim
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'o3clim'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nOznLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'mol mol^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'climatological ozone on prescribed pressure levels at '&
//'current time')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'o3clim', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'o3clim', r2Ptr)

! Define variable o3vmr
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'o3vmr'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'mol mol^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'ozone volume mixing ratio')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'o3vmr', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'o3vmr', r2Ptr)

! Define variable cldfrac
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'cldfrac'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'horizontal cloud fraction')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'cldfrac', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'cldfrac', r2Ptr)

! Define variable acsnom
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'acsnom'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated melted snow')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acsnom', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'acsnom', r1Ptr)

! Define variable acsnow
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'acsnow'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated snow')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'acsnow', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'acsnow', r1Ptr)

! Define variable canwat
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'canwat'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'water in canopy')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'canwat', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'canwat', r1Ptr)

! Define variable chklowq
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'chklowq'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'surface saturation flag')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'chklowq', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'chklowq', r1Ptr)

! Define variable grdflx
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'grdflx'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'ground heat flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'grdflx', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'grdflx', r1Ptr)

! Define variable lai
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'lai'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm^{-2} m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'leaf area index')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'lai', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'lai', r1Ptr)

! Define variable noahres
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'noahres'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'residual of the Noah surface energy budget')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'noahres', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'noahres', r1Ptr)

! Define variable potevp
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'potevp'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'accumulated potential evaporation')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'potevp', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'potevp', r1Ptr)

! Define variable sfc_albedo
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'sfc_albedo'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'surface albedo')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'sfc_albedo', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'sfc_albedo', r1Ptr)

! Define variable sfc_albedo_seaice
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'sfc_albedo_seaice'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'surface albedo over seaice')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'sfc_albedo_seaice', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'sfc_albedo_seaice', r1Ptr)

! Define variable sfc_emiss
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'sfc_emiss'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'surface emissivity')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'sfc_emiss', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'sfc_emiss', r1Ptr)

! Define variable sfc_emibck
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'sfc_emibck'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'background surface emissivity')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'sfc_emibck', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'sfc_emibck', r1Ptr)

! Define variable sfcrunoff
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'sfcrunoff'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'mm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'surface runoff')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'sfcrunoff', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'sfcrunoff', r1Ptr)

! Define variable smstav
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'smstav'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'surface moisture availability')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'smstav', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'smstav', r1Ptr)

! Define variable smstot
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'smstot'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm^{3} m^{-3}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'total soil mositure')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'smstot', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'smstot', r1Ptr)

! Define variable snopcx
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'snopcx'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'W m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'snow phase change heat flux')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'snopcx', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'snopcx', r1Ptr)

! Define variable snotime
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'snotime'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'initial number of time-steps since last snow fall')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'snotime', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'snotime', r1Ptr)

! Define variable sstsk
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'sstsk'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'skin sea-surface temperature')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'sstsk', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'sstsk', r1Ptr)

! Define variable sstsk_dtc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'sstsk_dtc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'skin sea-surface temperature cooling')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'sstsk_dtc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'sstsk_dtc', r1Ptr)

! Define variable sstsk_dtw
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'sstsk_dtw'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'skin sea-surface temperature warming')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'sstsk_dtw', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'sstsk_dtw', r1Ptr)

! Define variable thc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'thc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'Cal cm^{-2} K^{-1} s^{-0.5}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'thermal inertia')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'thc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'thc', r1Ptr)

! Define variable udrunoff
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'udrunoff'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'mm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'underground runoff')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'udrunoff', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'udrunoff', r1Ptr)

! Define variable xicem
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'xicem'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'sea-ice flag from previsous time-step')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'xicem', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'xicem', r1Ptr)

! Define variable z0
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'z0'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'roughness height')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'z0', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'z0', r1Ptr)

! Define variable znt
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'znt'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'roughness length')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'znt', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'znt', r1Ptr)

! Define variable zs
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'zs'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'depth of centers of soil layers')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'zs', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'zs', r1Ptr)

! Define variable t_oml
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 't_oml'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'ocean mixed layer temperature')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 't_oml', r1Ptr)
      call mpas_pool_add_field(block % allFields, 't_oml', r1Ptr)

! Define variable t_oml_initial
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 't_oml_initial'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'ocean mixed layer temperature at initial time')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 't_oml_initial', r1Ptr)
      call mpas_pool_add_field(block % allFields, 't_oml_initial', r1Ptr)

! Define variable t_oml_200m_initial
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 't_oml_200m_initial'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'ocean mixed layer 200 m mean temperature at initial time')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 't_oml_200m_initial', r1Ptr)
      call mpas_pool_add_field(block % allFields, 't_oml_200m_initial', r1Ptr)

! Define variable h_oml
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'h_oml'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'ocean mixed layer depth')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'h_oml', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'h_oml', r1Ptr)

! Define variable h_oml_initial
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'h_oml_initial'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'ocean mixed layer depth at initial time')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'h_oml_initial', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'h_oml_initial', r1Ptr)

! Define variable hu_oml
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'hu_oml'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm^2 s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'ocean mixed layer integrated u (zonal velocity)')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'hu_oml', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'hu_oml', r1Ptr)

! Define variable hv_oml
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'hv_oml'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm^2 s^{-1}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'ocean mixed layer integrated v (meridional velocity)')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'hv_oml', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'hv_oml', r1Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
         call mpas_pool_add_config(newSubPool, 'is_periodic', block % domain % is_periodic)
         call mpas_pool_add_config(newSubPool, 'x_period', block % domain % x_period)
         call mpas_pool_add_config(newSubPool, 'y_period', block % domain % y_period)
      end if

   end subroutine atm_generate_pool_diag_physics


   subroutine atm_generate_pool_tend_physics(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_io, only : MPAS_REAL_FILLVAL, MPAS_INT_FILLVAL, MPAS_CHAR_FILLVAL
      implicit none
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), dimension(:), pointer :: r0Ptr
      type (field1DReal), dimension(:), pointer :: r1Ptr
      type (field2DReal), dimension(:), pointer :: r2Ptr
      type (field3DReal), dimension(:), pointer :: r3Ptr
      type (field4DReal), dimension(:), pointer :: r4Ptr
      type (field5DReal), dimension(:), pointer :: r5Ptr
      type (field0DInteger), dimension(:), pointer :: i0Ptr
      type (field1DInteger), dimension(:), pointer :: i1Ptr
      type (field2DInteger), dimension(:), pointer :: i2Ptr
      type (field3DInteger), dimension(:), pointer :: i3Ptr
      type (field0DChar), dimension(:), pointer :: c0Ptr
      type (field1DChar), dimension(:), pointer :: c1Ptr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: mp_kessler_inActive
      logical, pointer :: mp_thompson_inActive
      logical, pointer :: mp_wsm6_inActive
      logical, pointer :: cu_grell_freitas_inActive
      logical, pointer :: cu_kain_fritsch_inActive
      logical, pointer :: cu_tiedtke_inActive
      logical, pointer :: bl_ysu_inActive
      logical, pointer :: bl_mynn_inActive
      logical, pointer :: iauActive
      logical, pointer :: limited_areaActive
      logical, pointer :: jedi_daActive
      logical, pointer :: tendenciesActive


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'mp_kessler_inActive', mp_kessler_inActive)
      call mpas_pool_get_package(packagePool, 'mp_thompson_inActive', mp_thompson_inActive)
      call mpas_pool_get_package(packagePool, 'mp_wsm6_inActive', mp_wsm6_inActive)
      call mpas_pool_get_package(packagePool, 'cu_grell_freitas_inActive', cu_grell_freitas_inActive)
      call mpas_pool_get_package(packagePool, 'cu_kain_fritsch_inActive', cu_kain_fritsch_inActive)
      call mpas_pool_get_package(packagePool, 'cu_tiedtke_inActive', cu_tiedtke_inActive)
      call mpas_pool_get_package(packagePool, 'bl_ysu_inActive', bl_ysu_inActive)
      call mpas_pool_get_package(packagePool, 'bl_mynn_inActive', bl_mynn_inActive)
      call mpas_pool_get_package(packagePool, 'iauActive', iauActive)
      call mpas_pool_get_package(packagePool, 'limited_areaActive', limited_areaActive)
      call mpas_pool_get_package(packagePool, 'jedi_daActive', jedi_daActive)
      call mpas_pool_get_package(packagePool, 'tendenciesActive', tendenciesActive)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'tend_physics', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'tend_physics', newSubPool)

! Define variable rthcuten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rthcuten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of potential temperature due to cumulus '&
//'convection')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (cu_grell_freitas_inActive .or. cu_kain_fritsch_inActive .or. cu_tiedtke_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rthcuten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rthcuten', r2Ptr)

! Define variable rqvcuten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rqvcuten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of water vapor mixing ratio due to cumulus '&
//'convection')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (cu_grell_freitas_inActive .or. cu_kain_fritsch_inActive .or. cu_tiedtke_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rqvcuten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rqvcuten', r2Ptr)

! Define variable rqccuten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rqccuten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of cloud water mixing ratio due to cumulus '&
//'convection')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (cu_grell_freitas_inActive .or. cu_kain_fritsch_inActive .or. cu_tiedtke_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rqccuten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rqccuten', r2Ptr)

! Define variable rqicuten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rqicuten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of cloud ice mixing ratio due to cumulus '&
//'convection')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (cu_grell_freitas_inActive .or. cu_kain_fritsch_inActive .or. cu_tiedtke_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rqicuten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rqicuten', r2Ptr)

! Define variable rqrcuten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rqrcuten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of rain mixing ratio due to cumulus convection')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (cu_kain_fritsch_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rqrcuten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rqrcuten', r2Ptr)

! Define variable rqscuten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rqscuten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of snow mixing ratio due to cumulus convection')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (cu_kain_fritsch_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rqscuten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rqscuten', r2Ptr)

! Define variable rqvdynten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rqvdynten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of water vapor due to horizontal and vertical '&
//'advections')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rqvdynten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rqvdynten', r2Ptr)

! Define variable rthdynten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rthdynten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of temperature due to horizontal and vertical '&
//'advections')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rthdynten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rthdynten', r2Ptr)

! Define variable rucuten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rucuten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of zonal wind due to cumulus convection')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (cu_grell_freitas_inActive .or. cu_tiedtke_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rucuten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rucuten', r2Ptr)

! Define variable rvcuten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rvcuten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of meridional wind due to cumulus convection')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (cu_grell_freitas_inActive .or. cu_tiedtke_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rvcuten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rvcuten', r2Ptr)

! Define variable rucuten_Edge
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rucuten_Edge'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rucuten_Edge', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rucuten_Edge', r2Ptr)

! Define variable rublten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rublten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of zonal wind due to pbl processes')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rublten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rublten', r2Ptr)

! Define variable rvblten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rvblten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of meridional wind due to pbl processes')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rvblten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rvblten', r2Ptr)

! Define variable rthblten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rthblten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of potential temperature due to pbl processes')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rthblten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rthblten', r2Ptr)

! Define variable rqvblten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rqvblten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of water vapor mixing ratio due to pbl processes')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rqvblten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rqvblten', r2Ptr)

! Define variable rqcblten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rqcblten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of cloud water mixing ratio due to pbl processes')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rqcblten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rqcblten', r2Ptr)

! Define variable rqiblten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rqiblten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of cloud ice mixing ratio due to pbl processes')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive .or. bl_ysu_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rqiblten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rqiblten', r2Ptr)

! Define variable rniblten
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rniblten'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'nb kg^{-1} s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of cloud ice number concentration due to pbl '&
//'processes')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (bl_mynn_inActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rniblten', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rniblten', r2Ptr)

! Define variable rublten_Edge
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rublten_Edge'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rublten_Edge', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rublten_Edge', r2Ptr)

! Define variable rthratensw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rthratensw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of potential temperature due to short wave '&
//'radiation')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rthratensw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rthratensw', r2Ptr)

! Define variable rthratenlw
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rthratenlw'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'tendency of potential temperature due to short wave '&
//'radiation')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'rthratenlw', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rthratenlw', r2Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
         call mpas_pool_add_config(newSubPool, 'is_periodic', block % domain % is_periodic)
         call mpas_pool_add_config(newSubPool, 'x_period', block % domain % x_period)
         call mpas_pool_add_config(newSubPool, 'y_period', block % domain % y_period)
      end if

   end subroutine atm_generate_pool_tend_physics


   subroutine atm_generate_pool_atm_input(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_io, only : MPAS_REAL_FILLVAL, MPAS_INT_FILLVAL, MPAS_CHAR_FILLVAL
      implicit none
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), dimension(:), pointer :: r0Ptr
      type (field1DReal), dimension(:), pointer :: r1Ptr
      type (field2DReal), dimension(:), pointer :: r2Ptr
      type (field3DReal), dimension(:), pointer :: r3Ptr
      type (field4DReal), dimension(:), pointer :: r4Ptr
      type (field5DReal), dimension(:), pointer :: r5Ptr
      type (field0DInteger), dimension(:), pointer :: i0Ptr
      type (field1DInteger), dimension(:), pointer :: i1Ptr
      type (field2DInteger), dimension(:), pointer :: i2Ptr
      type (field3DInteger), dimension(:), pointer :: i3Ptr
      type (field0DChar), dimension(:), pointer :: c0Ptr
      type (field1DChar), dimension(:), pointer :: c1Ptr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: mp_kessler_inActive
      logical, pointer :: mp_thompson_inActive
      logical, pointer :: mp_wsm6_inActive
      logical, pointer :: cu_grell_freitas_inActive
      logical, pointer :: cu_kain_fritsch_inActive
      logical, pointer :: cu_tiedtke_inActive
      logical, pointer :: bl_ysu_inActive
      logical, pointer :: bl_mynn_inActive
      logical, pointer :: iauActive
      logical, pointer :: limited_areaActive
      logical, pointer :: jedi_daActive
      logical, pointer :: tendenciesActive


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'mp_kessler_inActive', mp_kessler_inActive)
      call mpas_pool_get_package(packagePool, 'mp_thompson_inActive', mp_thompson_inActive)
      call mpas_pool_get_package(packagePool, 'mp_wsm6_inActive', mp_wsm6_inActive)
      call mpas_pool_get_package(packagePool, 'cu_grell_freitas_inActive', cu_grell_freitas_inActive)
      call mpas_pool_get_package(packagePool, 'cu_kain_fritsch_inActive', cu_kain_fritsch_inActive)
      call mpas_pool_get_package(packagePool, 'cu_tiedtke_inActive', cu_tiedtke_inActive)
      call mpas_pool_get_package(packagePool, 'bl_ysu_inActive', bl_ysu_inActive)
      call mpas_pool_get_package(packagePool, 'bl_mynn_inActive', bl_mynn_inActive)
      call mpas_pool_get_package(packagePool, 'iauActive', iauActive)
      call mpas_pool_get_package(packagePool, 'limited_areaActive', limited_areaActive)
      call mpas_pool_get_package(packagePool, 'jedi_daActive', jedi_daActive)
      call mpas_pool_get_package(packagePool, 'tendenciesActive', tendenciesActive)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'atm_input', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'atm_input', newSubPool)

! Define variable pin
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'pin'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .false.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nOznLevels'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'Pa')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'fixed pressure levels at which climatological ozone is '&
//'defined')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'pin', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'pin', r1Ptr)

! Define variable ozmixm
      allocate(r3Ptr(1))

! Setting up time level 1
      r3Ptr(1) % fieldName = 'ozmixm'
      r3Ptr(1) % isVarArray = .false.
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .false.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.
! Setting up dimensions
      r3Ptr(1) % dimNames(1) = 'nMonths'
      r3Ptr(1) % dimNames(2) = 'nOznLevels'
      r3Ptr(1) % dimNames(3) = 'nCells'
      r3Ptr(1) % defaultValue = 0.0
      allocate(r3Ptr(1) % attLists(1))
      allocate(r3Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'units', 'mol mol^{-1}')
      call mpas_add_att(r3Ptr(1) % attLists(1) % attList, 'long_name', 'monthly-mean climatological ozone defined at fixed '&
//'pressure levels')
      r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r3Ptr(1) % block => block

      r3Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ozmixm', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'ozmixm', r3Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
         call mpas_pool_add_config(newSubPool, 'is_periodic', block % domain % is_periodic)
         call mpas_pool_add_config(newSubPool, 'x_period', block % domain % x_period)
         call mpas_pool_add_config(newSubPool, 'y_period', block % domain % y_period)
      end if

   end subroutine atm_generate_pool_atm_input


   subroutine atm_generate_pool_sfc_input(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_io, only : MPAS_REAL_FILLVAL, MPAS_INT_FILLVAL, MPAS_CHAR_FILLVAL
      implicit none
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), dimension(:), pointer :: r0Ptr
      type (field1DReal), dimension(:), pointer :: r1Ptr
      type (field2DReal), dimension(:), pointer :: r2Ptr
      type (field3DReal), dimension(:), pointer :: r3Ptr
      type (field4DReal), dimension(:), pointer :: r4Ptr
      type (field5DReal), dimension(:), pointer :: r5Ptr
      type (field0DInteger), dimension(:), pointer :: i0Ptr
      type (field1DInteger), dimension(:), pointer :: i1Ptr
      type (field2DInteger), dimension(:), pointer :: i2Ptr
      type (field3DInteger), dimension(:), pointer :: i3Ptr
      type (field0DChar), dimension(:), pointer :: c0Ptr
      type (field1DChar), dimension(:), pointer :: c1Ptr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: mp_kessler_inActive
      logical, pointer :: mp_thompson_inActive
      logical, pointer :: mp_wsm6_inActive
      logical, pointer :: cu_grell_freitas_inActive
      logical, pointer :: cu_kain_fritsch_inActive
      logical, pointer :: cu_tiedtke_inActive
      logical, pointer :: bl_ysu_inActive
      logical, pointer :: bl_mynn_inActive
      logical, pointer :: iauActive
      logical, pointer :: limited_areaActive
      logical, pointer :: jedi_daActive
      logical, pointer :: tendenciesActive


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'mp_kessler_inActive', mp_kessler_inActive)
      call mpas_pool_get_package(packagePool, 'mp_thompson_inActive', mp_thompson_inActive)
      call mpas_pool_get_package(packagePool, 'mp_wsm6_inActive', mp_wsm6_inActive)
      call mpas_pool_get_package(packagePool, 'cu_grell_freitas_inActive', cu_grell_freitas_inActive)
      call mpas_pool_get_package(packagePool, 'cu_kain_fritsch_inActive', cu_kain_fritsch_inActive)
      call mpas_pool_get_package(packagePool, 'cu_tiedtke_inActive', cu_tiedtke_inActive)
      call mpas_pool_get_package(packagePool, 'bl_ysu_inActive', bl_ysu_inActive)
      call mpas_pool_get_package(packagePool, 'bl_mynn_inActive', bl_mynn_inActive)
      call mpas_pool_get_package(packagePool, 'iauActive', iauActive)
      call mpas_pool_get_package(packagePool, 'limited_areaActive', limited_areaActive)
      call mpas_pool_get_package(packagePool, 'jedi_daActive', jedi_daActive)
      call mpas_pool_get_package(packagePool, 'tendenciesActive', tendenciesActive)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'sfc_input', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'sfc_input', newSubPool)

! Define variable isltyp
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'isltyp'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .false.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'dominant soil category')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'isltyp', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'isltyp', i1Ptr)

! Define variable ivgtyp
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'ivgtyp'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .false.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'dominant vegetation category')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ivgtyp', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'ivgtyp', i1Ptr)

! Define variable mminlu
      allocate(c0Ptr(1))

! Setting up time level 1
      c0Ptr(1) % fieldName = 'mminlu'
      c0Ptr(1) % isVarArray = .false.
      c0Ptr(1) % isDecomposed = .false.
      c0Ptr(1) % hasTimeDimension = .false.
      c0Ptr(1) % defaultValue = ''
      c0Ptr(1) % scalar = ''
      allocate(c0Ptr(1) % attLists(1))
      allocate(c0Ptr(1) % attLists(1) % attList)
      call mpas_add_att(c0Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(c0Ptr(1) % attLists(1) % attList, 'long_name', 'land use classification')
      c0Ptr(1) % missingValue = MPAS_CHAR_FILLVAL
      c0Ptr(1) % block => block

      c0Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'mminlu', c0Ptr)
      call mpas_pool_add_field(block % allFields, 'mminlu', c0Ptr)

! Define variable isice_lu
      allocate(i0Ptr(1))

! Setting up time level 1
      i0Ptr(1) % fieldName = 'isice_lu'
      i0Ptr(1) % isVarArray = .false.
      i0Ptr(1) % isDecomposed = .false.
      i0Ptr(1) % hasTimeDimension = .false.
      i0Ptr(1) % defaultValue = 24
      i0Ptr(1) % scalar = 24
      allocate(i0Ptr(1) % attLists(1))
      allocate(i0Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i0Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i0Ptr(1) % attLists(1) % attList, 'long_name', 'Index category for snow/ice')
      i0Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i0Ptr(1) % block => block

      i0Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'isice', i0Ptr)
      call mpas_pool_add_field(block % allFields, 'isice_lu', i0Ptr)

! Define variable iswater_lu
      allocate(i0Ptr(1))

! Setting up time level 1
      i0Ptr(1) % fieldName = 'iswater_lu'
      i0Ptr(1) % isVarArray = .false.
      i0Ptr(1) % isDecomposed = .false.
      i0Ptr(1) % hasTimeDimension = .false.
      i0Ptr(1) % defaultValue = 16
      i0Ptr(1) % scalar = 16
      allocate(i0Ptr(1) % attLists(1))
      allocate(i0Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i0Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i0Ptr(1) % attLists(1) % attList, 'long_name', 'Index category for water')
      i0Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i0Ptr(1) % block => block

      i0Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'iswater', i0Ptr)
      call mpas_pool_add_field(block % allFields, 'iswater_lu', i0Ptr)

! Define variable landmask
      allocate(i1Ptr(1))

! Setting up time level 1
      i1Ptr(1) % fieldName = 'landmask'
      i1Ptr(1) % isVarArray = .false.
      i1Ptr(1) % isDecomposed = .true.
      i1Ptr(1) % hasTimeDimension = .false.
      i1Ptr(1) % isPersistent = .true.
      i1Ptr(1) % isActive = .false.
! Setting up dimensions
      i1Ptr(1) % dimNames(1) = 'nCells'
      i1Ptr(1) % defaultValue = 0
      allocate(i1Ptr(1) % attLists(1))
      allocate(i1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(i1Ptr(1) % attLists(1) % attList, 'long_name', 'land-ocean mask (1=land ; 0=ocean)')
      i1Ptr(1) % missingValue = MPAS_INT_FILLVAL
      i1Ptr(1) % block => block

      i1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'landmask', i1Ptr)
      call mpas_pool_add_field(block % allFields, 'landmask', i1Ptr)

! Define variable shdmin
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'shdmin'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'minimum fractional coverage of annual green vegetation '&
//'fraction')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'shdmin', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'shdmin', r1Ptr)

! Define variable shdmax
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'shdmax'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'maximum fractional coverage of annual green vegetation '&
//'fraction')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'shdmax', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'shdmax', r1Ptr)

! Define variable snoalb
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'snoalb'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'annual maximum snow albedo')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'snoalb', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'snoalb', r1Ptr)

! Define variable ter
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'ter'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'terrain height')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ter', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'ter', r1Ptr)

! Define variable albedo12m
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'albedo12m'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nMonths'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'monthly-mean climatological surface albedo')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'albedo12m', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'albedo12m', r2Ptr)

! Define variable greenfrac
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'greenfrac'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .false.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nMonths'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'monthly-mean climatological greeness fraction')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'greenfrac', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'greenfrac', r2Ptr)

! Define variable dzs
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'dzs'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nSoilLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'soil layer thickness')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'dzs', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'dzs', r2Ptr)

! Define variable sfc_albbck
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'sfc_albbck'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'background surface albedo')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'sfc_albbck', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'sfc_albbck', r1Ptr)

! Define variable skintemp
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'skintemp'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'ground or water surface temperature')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'skintemp', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'skintemp', r1Ptr)

! Define variable snow
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'snow'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-2}')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'snow water equivalent')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'snow', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'snow', r1Ptr)

! Define variable snowc
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'snowc'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'flag for snow on ground (=0 no snow; =1,otherwise')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'snowc', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'snowc', r1Ptr)

! Define variable snowh
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'snowh'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'physical snow depth')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'snowh', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'snowh', r1Ptr)

! Define variable sst
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'sst'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'sea-surface temperature')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'sst', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'sst', r1Ptr)

! Define variable tmn
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'tmn'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'deep soil temperature')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'tmn', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'tmn', r1Ptr)

! Define variable vegfra
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'vegfra'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'vegetation fraction')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'vegfra', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'vegfra', r1Ptr)

! Define variable seaice
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'seaice'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'sea-ice flag (0=no seaice; =1 otherwise)')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'seaice', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'seaice', r1Ptr)

! Define variable xice
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'xice'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'fractional area coverage of sea-ice')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'xice', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'xice', r1Ptr)

! Define variable xland
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'xland'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .true.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'land-ocean mask (1=land including sea-ice ; 2=ocean)')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'xland', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'xland', r1Ptr)

! Define variable smcrel
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'smcrel'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nSoilLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm3 m^{-3}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'soil moisture threshold below which transpiration begins '&
//'to stress')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'smcrel', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'smcrel', r2Ptr)

! Define variable sh2o
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'sh2o'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nSoilLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm3 m^{-3}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'soil equivalent liquid water ')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'sh2o', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'sh2o', r2Ptr)

! Define variable smois
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'smois'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nSoilLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm3 m^{-3}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'soil moisture')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'smois', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'smois', r2Ptr)

! Define variable tslb
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'tslb'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nSoilLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'soil layer temperature')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      r2Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'tslb', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'tslb', r2Ptr)

! Define variable var2d
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'var2d'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'm')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'standard deviation of subgrid-scale orography')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'var2d', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'var2d', r1Ptr)

! Define variable con
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'con'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'orographic convexity')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'con', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'con', r1Ptr)

! Define variable oa1
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'oa1'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'asymmetry of subgrid-scale orography for westerly flow')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'oa1', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'oa1', r1Ptr)

! Define variable oa2
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'oa2'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'asymmetry of subgrid-scale orography for southerly flow')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'oa2', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'oa2', r1Ptr)

! Define variable oa3
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'oa3'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'asymmetry of subgrid-scale orography for south-westerly '&
//'flow')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'oa3', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'oa3', r1Ptr)

! Define variable oa4
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'oa4'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'asymmetry of subgrid-scale orography for north-westerly '&
//'flow')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'oa4', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'oa4', r1Ptr)

! Define variable ol1
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'ol1'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'effective orographic length for westerly flow')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ol1', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'ol1', r1Ptr)

! Define variable ol2
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'ol2'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'effective orographic length for southerly flow')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ol2', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'ol2', r1Ptr)

! Define variable ol3
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'ol3'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitles')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'effective orographic length for south-westerly flow')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ol3', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'ol3', r1Ptr)

! Define variable ol4
      allocate(r1Ptr(1))

! Setting up time level 1
      r1Ptr(1) % fieldName = 'ol4'
      r1Ptr(1) % isVarArray = .false.
      r1Ptr(1) % isDecomposed = .true.
      r1Ptr(1) % hasTimeDimension = .false.
      r1Ptr(1) % isPersistent = .true.
      r1Ptr(1) % isActive = .false.
! Setting up dimensions
      r1Ptr(1) % dimNames(1) = 'nCells'
      r1Ptr(1) % defaultValue = 0.0
      allocate(r1Ptr(1) % attLists(1))
      allocate(r1Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'units', 'unitless')
      call mpas_add_att(r1Ptr(1) % attLists(1) % attList, 'long_name', 'effective orographic length for north-westerly flow')
      r1Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r1Ptr(1) % block => block

      r1Ptr(1) % isActive = .true.
      call mpas_pool_add_field(newSubPool, 'ol4', r1Ptr)
      call mpas_pool_add_field(block % allFields, 'ol4', r1Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
         call mpas_pool_add_config(newSubPool, 'is_periodic', block % domain % is_periodic)
         call mpas_pool_add_config(newSubPool, 'x_period', block % domain % x_period)
         call mpas_pool_add_config(newSubPool, 'y_period', block % domain % y_period)
      end if

   end subroutine atm_generate_pool_sfc_input


   subroutine atm_generate_pool_tend_iau(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_io, only : MPAS_REAL_FILLVAL, MPAS_INT_FILLVAL, MPAS_CHAR_FILLVAL
      implicit none
      type (block_type), intent(inout), pointer :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      type (field0DReal), dimension(:), pointer :: r0Ptr
      type (field1DReal), dimension(:), pointer :: r1Ptr
      type (field2DReal), dimension(:), pointer :: r2Ptr
      type (field3DReal), dimension(:), pointer :: r3Ptr
      type (field4DReal), dimension(:), pointer :: r4Ptr
      type (field5DReal), dimension(:), pointer :: r5Ptr
      type (field0DInteger), dimension(:), pointer :: i0Ptr
      type (field1DInteger), dimension(:), pointer :: i1Ptr
      type (field2DInteger), dimension(:), pointer :: i2Ptr
      type (field3DInteger), dimension(:), pointer :: i3Ptr
      type (field0DChar), dimension(:), pointer :: c0Ptr
      type (field1DChar), dimension(:), pointer :: c1Ptr

      type (mpas_pool_type), pointer :: newSubPool
      integer :: group_counter
      logical :: group_started
      integer :: group_start
      integer :: index_counter
      integer, pointer :: const_index

      logical, pointer :: mp_kessler_inActive
      logical, pointer :: mp_thompson_inActive
      logical, pointer :: mp_wsm6_inActive
      logical, pointer :: cu_grell_freitas_inActive
      logical, pointer :: cu_kain_fritsch_inActive
      logical, pointer :: cu_tiedtke_inActive
      logical, pointer :: bl_ysu_inActive
      logical, pointer :: bl_mynn_inActive
      logical, pointer :: iauActive
      logical, pointer :: limited_areaActive
      logical, pointer :: jedi_daActive
      logical, pointer :: tendenciesActive


      integer :: numConstituents

      nullify(newSubPool)
      group_counter = -1
      group_started = .false.
      group_start = -1
      call mpas_pool_get_package(packagePool, 'mp_kessler_inActive', mp_kessler_inActive)
      call mpas_pool_get_package(packagePool, 'mp_thompson_inActive', mp_thompson_inActive)
      call mpas_pool_get_package(packagePool, 'mp_wsm6_inActive', mp_wsm6_inActive)
      call mpas_pool_get_package(packagePool, 'cu_grell_freitas_inActive', cu_grell_freitas_inActive)
      call mpas_pool_get_package(packagePool, 'cu_kain_fritsch_inActive', cu_kain_fritsch_inActive)
      call mpas_pool_get_package(packagePool, 'cu_tiedtke_inActive', cu_tiedtke_inActive)
      call mpas_pool_get_package(packagePool, 'bl_ysu_inActive', bl_ysu_inActive)
      call mpas_pool_get_package(packagePool, 'bl_mynn_inActive', bl_mynn_inActive)
      call mpas_pool_get_package(packagePool, 'iauActive', iauActive)
      call mpas_pool_get_package(packagePool, 'limited_areaActive', limited_areaActive)
      call mpas_pool_get_package(packagePool, 'jedi_daActive', jedi_daActive)
      call mpas_pool_get_package(packagePool, 'tendenciesActive', tendenciesActive)

      allocate(newSubPool)
      call mpas_pool_create_pool(newSubPool)
      call mpas_pool_add_subpool(structPool, 'tend_iau', newSubPool)
      call mpas_pool_add_subpool(block % allStructs, 'tend_iau', newSubPool)

! Define var array scalars_amb
      allocate(r3Ptr(1))
      index_counter = 0
      group_counter = -1
      group_start = -1
      group_started = .false.

! Starting group moist
! Define constituent var qv_amb
! My Packages are iau
      if (iauActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qv', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
           if (associated(newSubPool)) then
              call mpas_pool_add_dimension(newSubPool, 'index_qv', -1)
           end if
      end if
! Define constituent var qc_amb
! My packages are bl_mynn_in;bl_ysu_in;cu_tiedtke_in;mp_kessler_in;mp_thompson_in;mp_wsm6_in
      if (bl_mynn_inActive .or. bl_ysu_inActive .or. cu_tiedtke_inActive .or. mp_kessler_inActive .or. mp_thompson_inActive .or. &
mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qc', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qc', -1)
         end if
      end if
! Define constituent var qr_amb
! My packages are mp_kessler_in;mp_thompson_in;mp_wsm6_in
      if (mp_kessler_inActive .or. mp_thompson_inActive .or. mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qr', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qr', -1)
         end if
      end if
! Define constituent var qi_amb
! My packages are bl_mynn_in;bl_ysu_in;cu_tiedtke_in;mp_thompson_in;mp_wsm6_in
      if (bl_mynn_inActive .or. bl_ysu_inActive .or. cu_tiedtke_inActive .or. mp_thompson_inActive .or. mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qi', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qi', -1)
         end if
      end if
! Define constituent var qs_amb
! My packages are mp_thompson_in;mp_wsm6_in
      if (mp_thompson_inActive .or. mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qs', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qs', -1)
         end if
      end if
! Define constituent var qg_amb
! My packages are mp_thompson_in;mp_wsm6_in
      if (mp_thompson_inActive .or. mp_wsm6_inActive) then
         index_counter = index_counter + 1
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qg', index_counter)
         end if
         group_counter = group_counter + 1
         if (.not. group_started) then
            group_start = index_counter
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', group_start)
            end if
            group_started = .true.
         end if
      else
         if (associated(newSubPool)) then
            call mpas_pool_add_dimension(newSubPool, 'index_qg', -1)
         end if
      end if
         if (.not. group_started) then
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_start', -1)
               call mpas_pool_add_dimension(newSubPool, 'moist_end', -1)
            end if
         else
            group_started = .false.
            if (associated(newSubPool)) then
               call mpas_pool_add_dimension(newSubPool, 'moist_end', index_counter)
            end if
         end if
! End of group       

      numConstituents = index_counter
      if (associated(newSubPool)) then
         call mpas_pool_add_dimension(newSubPool, 'num_scalars_amb', numConstituents)
      end if
! Defining time level 1
      allocate( r3Ptr(1) % constituentNames(numConstituents) )
      r3Ptr(1) % fieldName = 'scalars_amb'
      r3Ptr(1) % isDecomposed = .true.
      r3Ptr(1) % hasTimeDimension = .true.
      r3Ptr(1) % isVarArray = .true.
      r3Ptr(1) % isPersistent = .true.
      r3Ptr(1) % isActive = .false.

      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qv', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qv_amb'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qc', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qc_amb'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qr', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qr_amb'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qi', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qi_amb'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qs', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qs_amb'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qg', const_index)
      end if
      if (const_index > 0) then
         r3Ptr(1) % constituentNames(const_index) = 'qg_amb'
      end if

! Setup dimensions for       
      r3Ptr(1) % dimNames(1) = 'num_scalars_amb'
      r3Ptr(1) % dimNames(2) = 'nVertLevels'
      r3Ptr(1) % dimNames(3) = 'nCells'

      r3Ptr(1) % defaultValue = 0.0
      allocate(r3Ptr(1) % attLists(size(r3Ptr(1) % constituentNames, dim=1)))
      do index_counter = 1, size(r3Ptr(1) % constituentNames, dim=1)
         allocate(r3Ptr(1) % attLists(index_counter) % attList)
      end do
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qv', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Water vapor mixing ratio increment')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'qv_amb'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qc', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Cloud water mixing ratio increment')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'qc_amb'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qr', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Rain water mixing ratio increment')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'qr_amb'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qi', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Ice mixing ratio increment')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'qi_amb'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qs', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Snow mixing ratio increment')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'qs_amb'
      end if
      if (associated(newSubPool)) then
         call mpas_pool_get_dimension(newSubPool, 'index_qg', const_index)
      end if
      if (const_index > 0) then
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'long_name', 'Graupel mixing ratio increment')
         call mpas_add_att(r3Ptr(1) % attLists(const_index) % attList, 'units', 'kg kg^{-1}')
         r3Ptr(1) % missingValue = MPAS_REAL_FILLVAL
         r3Ptr(1) % constituentNames(const_index) = 'qg_amb'
      end if
      r3Ptr(1) % block => block

      if (iauActive .or. bl_mynn_inActive .or. bl_ysu_inActive .or. cu_tiedtke_inActive .or. mp_kessler_inActive .or. &
mp_thompson_inActive .or. mp_wsm6_inActive) then
         r3Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'scalars', r3Ptr)
      call mpas_pool_add_field(block % allFields, 'scalars_amb', r3Ptr)

! Define variable u_amb
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'u_amb'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nEdges'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'm s^{-1}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Horizontal normal velocity increment')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (iauActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'u', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'u_amb', r2Ptr)

! Define variable rho_amb
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'rho_amb'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'kg m^{-3}')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Dry air density increment')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (iauActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'rho', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'rho_amb', r2Ptr)

! Define variable theta_amb
      allocate(r2Ptr(1))

! Setting up time level 1
      r2Ptr(1) % fieldName = 'theta_amb'
      r2Ptr(1) % isVarArray = .false.
      r2Ptr(1) % isDecomposed = .true.
      r2Ptr(1) % hasTimeDimension = .true.
      r2Ptr(1) % isPersistent = .true.
      r2Ptr(1) % isActive = .false.
! Setting up dimensions
      r2Ptr(1) % dimNames(1) = 'nVertLevels'
      r2Ptr(1) % dimNames(2) = 'nCells'
      r2Ptr(1) % defaultValue = 0.0
      allocate(r2Ptr(1) % attLists(1))
      allocate(r2Ptr(1) % attLists(1) % attList)
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'units', 'K')
      call mpas_add_att(r2Ptr(1) % attLists(1) % attList, 'long_name', 'Potential temperature increment')
      r2Ptr(1) % missingValue = MPAS_REAL_FILLVAL
      r2Ptr(1) % block => block

      if (iauActive) then
         r2Ptr(1) % isActive = .true.
      end if
      call mpas_pool_add_field(newSubPool, 'theta', r2Ptr)
      call mpas_pool_add_field(block % allFields, 'theta_amb', r2Ptr)



      if (associated(newSubPool)) then
         call mpas_pool_add_config(newSubPool, 'on_a_sphere', block % domain % on_a_sphere)
         call mpas_pool_add_config(newSubPool, 'sphere_radius', block % domain % sphere_radius)
         call mpas_pool_add_config(newSubPool, 'is_periodic', block % domain % is_periodic)
         call mpas_pool_add_config(newSubPool, 'x_period', block % domain % x_period)
         call mpas_pool_add_config(newSubPool, 'y_period', block % domain % y_period)
      end if

   end subroutine atm_generate_pool_tend_iau


   subroutine atm_generate_structs(block, structPool, dimensionPool, packagePool)
      use mpas_derived_types
      use mpas_io_units
      implicit none
      type (block_type), pointer, intent(inout) :: block
      type (mpas_pool_type), intent(inout) :: structPool
      type (mpas_pool_type), intent(inout) :: dimensionPool
      type (mpas_pool_type), intent(in) :: packagePool

      call atm_generate_pool_mesh(block, structPool, dimensionPool, packagePool)

      call atm_generate_pool_state(block, structPool, dimensionPool, packagePool)

      call atm_generate_pool_diag(block, structPool, dimensionPool, packagePool)

      call atm_generate_pool_tend(block, structPool, dimensionPool, packagePool)

      call atm_generate_pool_lbc(block, structPool, dimensionPool, packagePool)

      call atm_generate_pool_diag_physics(block, structPool, dimensionPool, packagePool)

      call atm_generate_pool_tend_physics(block, structPool, dimensionPool, packagePool)

      call atm_generate_pool_atm_input(block, structPool, dimensionPool, packagePool)

      call atm_generate_pool_sfc_input(block, structPool, dimensionPool, packagePool)

      call atm_generate_pool_tend_iau(block, structPool, dimensionPool, packagePool)

   end subroutine atm_generate_structs
