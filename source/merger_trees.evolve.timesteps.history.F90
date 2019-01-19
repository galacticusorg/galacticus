!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!% Implements a merger tree evolution timestepping class which limits the step the next epoch at which to store global history.

  use FGSL               , only : fgsl_interp_accel
  use Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass
  
  !# <mergerTreeEvolveTimestep name="mergerTreeEvolveTimestepHistory" defaultThreadPrivate="yes" autoHook="yes">
  !#  <description>A merger tree evolution timestepping class which limits the step the next epoch at which to store global history.</description>
  !# </mergerTreeEvolveTimestep>
  type, extends(mergerTreeEvolveTimestepClass) :: mergerTreeEvolveTimestepHistory
     !% Implementation of a merger tree evolution timestepping class which limits the step the next epoch at which to store global history.
     private
     class           (cosmologyFunctionsClass), pointer                   :: cosmologyFunctions_
     logical                                                              :: diskActive               , spheroidActive
     integer                                                              :: historyCount
     double precision                                                     :: timeBegin                , timeEnd
     double precision                         , allocatable, dimension(:) :: rateStarFormationDisk    , densityStellarDisk    , &
          &                                                                  expansionFactor          , densityColdGas        , &
          &                                                                  densityHotHaloGas        , densityNode           , &
          &                                                                  rateStarFormationSpheroid, densityStellarSpheroid, &
          &                                                                  rateStarFormation        , densityStellar        , &
          &                                                                  time
     type            (fgsl_interp_accel      )                            :: interpolationAccelerator
   contains
     final     ::                 historyDestructor
     procedure :: timeEvolveTo => historyTimeEvolveTo
     procedure :: autoHook     => historyAutoHook
  end type mergerTreeEvolveTimestepHistory

  interface mergerTreeEvolveTimestepHistory
     !% Constructors for the {\normalfont \ttfamily history} merger tree evolution timestep class.
     module procedure historyConstructorParameters
     module procedure historyConstructorInternal
  end interface mergerTreeEvolveTimestepHistory

contains

  function historyConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily history} merger tree evolution timestep class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (mergerTreeEvolveTimestepHistory)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass        ), pointer       :: cosmologyFunctions_
    integer                                                          :: historyCount
    double precision                                                 :: timeBegin          , timeEnd, &
         &                                                              ageUniverse
    
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    ageUniverse=cosmologyFunctions_%cosmicTime(1.0d0)
    !# <inputParameter>
    !#   <name>timeBegin</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.05d0*ageUniverse</defaultValue>
    !#   <description>The earliest time at which to tabulate the volume averaged history of galaxies (in Gyr).</description>
    !#   <group>timeStepping</group>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>timeEnd</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>ageUniverse</defaultValue>
    !#   <description>The latest time at which to tabulate the volume averaged history of galaxies (in Gyr).</description>
    !#   <group>timeStepping</group>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>historyCount</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>30</defaultValue>
    !#   <description>The number of steps (spaced logarithmically in cosmic time) at which to tabulate the volume averaged history of galaxies.</description>
    !#   <group>timeStepping</group>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>
    self=mergerTreeEvolveTimestepHistory(historyCount,timeBegin,timeEnd,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function historyConstructorParameters

  function historyConstructorInternal(historyCount,timeBegin,timeEnd,cosmologyFunctions_) result(self)
    !% Constructor for the {\normalfont \ttfamily history} merger tree evolution timestep class which takes a parameter set as input.
    use, intrinsic :: ISO_C_Binding
    use            :: Galacticus_Nodes , only : defaultDiskComponent, defaultSpheroidComponent
    use            :: Memory_Management
    use            :: Numerical_Ranges
    implicit none
    type            (mergerTreeEvolveTimestepHistory)                        :: self
    integer                                          , intent(in   )         :: historyCount
    double precision                                 , intent(in   )         :: timeBegin          , timeEnd
    class           (cosmologyFunctionsClass        ), intent(in   ), target :: cosmologyFunctions_    
    integer         (c_size_t                       )                        :: timeIndex
    !# <constructorAssign variables="historyCount, timeBegin, timeEnd, *cosmologyFunctions_"/>

    ! Determine if we have active components that can provide star formation rates.
    self%diskActive    =    defaultDiskComponent%starFormationRateIsGettable()
    self%spheroidActive=defaultSpheroidComponent%starFormationRateIsGettable()
    ! Allocate storage arrays.
    call allocateArray(self%time                     ,[self%historyCount])
    call allocateArray(self%expansionFactor          ,[self%historyCount])
    call allocateArray(self%rateStarFormation        ,[self%historyCount])
    call allocateArray(self%rateStarFormationDisk    ,[self%historyCount])
    call allocateArray(self%rateStarFormationSpheroid,[self%historyCount])
    call allocateArray(self%densityStellar           ,[self%historyCount])
    call allocateArray(self%densityStellarDisk       ,[self%historyCount])
    call allocateArray(self%densityStellarSpheroid   ,[self%historyCount])
    call allocateArray(self%densityColdGas           ,[self%historyCount])
    call allocateArray(self%densityHotHaloGas        ,[self%historyCount])
    call allocateArray(self%densityNode              ,[self%historyCount])
    ! Initialize arrays.
    self%time=Make_Range(self%timeBegin,self%timeEnd,self%historyCount,rangeTypeLogarithmic)
    do timeIndex=1,self%historyCount
       self%expansionFactor(timeIndex)=self%cosmologyFunctions_%expansionFactor(self%time(timeIndex))
    end do
    self%rateStarFormation        =0.0d0
    self%rateStarFormationDisk    =0.0d0
    self%rateStarFormationSpheroid=0.0d0
    self%densityStellar           =0.0d0
    self%densityStellarDisk       =0.0d0
    self%densityStellarSpheroid   =0.0d0
    self%densityColdGas           =0.0d0
    self%densityHotHaloGas        =0.0d0
    self%densityNode              =0.0d0
    return
  end function historyConstructorInternal

  subroutine historyAutoHook(self)
    !% Create a hook to the HDF5 pre-close event to allow us to finalize and write out our data.
    use Events_Hooks
    implicit none
    class(mergerTreeEvolveTimestepHistory), intent(inout) :: self
    
    call hdf5PreCloseEvent%attach(self,historyWrite)
    return
  end subroutine historyAutoHook

  subroutine historyDestructor(self)
    !% Destructor for the {\normalfont \ttfamily history} merger tree evolution timestep class.
    implicit none
    type(mergerTreeEvolveTimestepHistory), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"/>
    return
  end subroutine historyDestructor

  double precision function historyTimeEvolveTo(self,node,task,taskSelf,report,lockNode,lockType)
    !% Determine a suitable timestep for {\normalfont \ttfamily node} using the history method.
    use            :: Galacticus_Nodes       , only : nodeComponentBasic
    use, intrinsic :: ISO_C_Binding
    use            :: Numerical_Interpolation
    use            :: Evolve_To_Time_Reports
    use            :: ISO_Varying_String
    implicit none
    class           (mergerTreeEvolveTimestepHistory), intent(inout), target  :: self
    type            (treeNode                       ), intent(inout), target  :: node
    procedure       (timestepTask                   ), intent(  out), pointer :: task
    class           (*                              ), intent(  out), pointer :: taskSelf
    logical                                          , intent(in   )          :: report
    type            (treeNode                       ), intent(  out), pointer :: lockNode
    type            (varying_string                 ), intent(  out)          :: lockType
    class           (nodeComponentBasic             )               , pointer :: basic
    integer         (c_size_t                       )                         :: timeIndex
    double precision                                                          :: time

    ! Determine how long until next available timestep.
    basic     => node %basic()
    time      =  basic%time ()
    timeIndex =  Interpolate_Locate(self%time,self%interpolationAccelerator,time)
    if (time < self%time(timeIndex+1)) then
       historyTimeEvolveTo =  self%time(timeIndex+1)
       task                => historyStore
       taskSelf            => self
    else
       historyTimeEvolveTo =  huge(0.0d0)
       task                => null(     )
       taskSelf            => null(     )
    end if
    lockNode => node
    lockType =  "history"
    if (report) call Evolve_To_Time_Report("history: ",historyTimeEvolveTo)
    return
  end function historyTimeEvolveTo

  subroutine historyStore(self,tree,node,deadlockStatus)
    !% Store various properties in global arrays.
    use            :: Galacticus_Nodes                  , only : nodeComponentBasic, nodeComponentDisk, nodeComponentSpheroid
    use, intrinsic :: ISO_C_Binding
    use            :: Numerical_Interpolation
    use            :: Galactic_Structure_Options
    use            :: Galactic_Structure_Enclosed_Masses
    use            :: Galacticus_Error
    implicit none
    class           (*                    ), intent(inout)          :: self
    type            (mergerTree           ), intent(in   )          :: tree
    type            (treeNode             ), intent(inout), pointer :: node
    integer                                , intent(inout)          :: deadlockStatus
    class           (nodeComponentBasic   )               , pointer :: basic
    class           (nodeComponentDisk    )               , pointer :: disk
    class           (nodeComponentSpheroid)               , pointer :: spheroid
    integer         (c_size_t             )                         :: timeIndex
    double precision                                                :: rateStarFormationDisk    , massHotGas, &
         &                                                             rateStarFormationSpheroid, time
    !GCC$ attributes unused :: deadlockStatus

    select type (self)
    class is (mergerTreeEvolveTimestepHistory)
       ! Get components.
       basic    => node%basic   ()
       disk     => node%disk    ()
       spheroid => node%spheroid()
       ! Get current cosmic time.
       time=basic%time()
       ! Determine how long until next available timestep.
       if (time == self%time(self%historyCount)) then
          timeIndex=self%historyCount
       else
          timeIndex=Interpolate_Locate(self%time,self%interpolationAccelerator,time)
       end if
       ! Extract disk and spheroid star formation rates.
       rateStarFormationDisk                        =0.0d0
       rateStarFormationSpheroid                    =0.0d0
       if (self%diskActive)     rateStarFormationDisk    =disk    %starFormationRate()
       if (self%spheroidActive) rateStarFormationSpheroid=spheroid%starFormationRate()
       ! Accumulate the properties.
       ! Star formation rate:
       self%rateStarFormation                       (timeIndex)=+self%rateStarFormation          (timeIndex)                                                           &
            &                                                   +(rateStarFormationDisk+rateStarFormationSpheroid)                                                     &
            &                                                   *tree%volumeWeight
       self%rateStarFormationDisk                   (timeIndex)=+self%rateStarFormationDisk      (timeIndex)                                                           &
            &                                                   + rateStarFormationDisk                                                                                &
            &                                                   *tree%volumeWeight
       self%rateStarFormationSpheroid               (timeIndex)=+self%rateStarFormationSpheroid  (timeIndex)                                                           &
            &                                                   +                       rateStarFormationSpheroid                                                      &
            &                                                   *tree%volumeWeight
       ! Stellar densities.
       self%densityStellar                          (timeIndex)=+self%densityStellar             (timeIndex)                                                           &
            &                                                   +  Galactic_Structure_Enclosed_Mass(node,                                    massType=massTypeStellar) &
            &                                                   *tree%volumeWeight
       self%densityStellarDisk                      (timeIndex)=+self%densityStellarDisk         (timeIndex)                                                           &
            &                                                   +  Galactic_Structure_Enclosed_Mass(node,componentType=componentTypeDisk    ,massType=massTypeStellar) &
            &                                                   *tree%volumeWeight
       self%densityStellarSpheroid                  (timeIndex)=+self%densityStellarSpheroid     (timeIndex)                                                           &
            &                                                   +  Galactic_Structure_Enclosed_Mass(node,componentType=componentTypeSpheroid,massType=massTypeStellar) &
            &                                                   *tree%volumeWeight

       ! Hot gas density.
       massHotGas                                              =+  Galactic_Structure_Enclosed_Mass(node,componentType=componentTypeHotHalo                          )
       self%densityHotHaloGas                       (timeIndex)=+self%densityHotHaloGas          (timeIndex)                                                           &
            &                                                   +massHotGas                                                                                            &
            &                                                   *tree%volumeWeight
       ! Galactic gas density.
       self%densityColdGas                          (timeIndex)=+self%densityColdGas             (timeIndex)                                                           &
            &                                                   +(                                                                                                     &
            &                                                     +Galactic_Structure_Enclosed_Mass(node,massType=massTypeGaseous                                    ) &
            &                                                     -massHotGas                                                                                          &
            &                                                    )                                                                                                     &
            &                                                   *tree%volumeWeight
       ! Node density
       if (.not.node%isSatellite()) self%densityNode(timeIndex)=+self%densityNode                (timeIndex)                                                           &
            &                                                   +basic%mass()                                                                                          &
            &                                                   *tree%volumeWeight
    class default
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine historyStore

  subroutine historyWrite(self)
    !% Store the global history data to the \glc\ output file.
    use Galacticus_HDF5
    use Galacticus_Error
    use Numerical_Constants_Astronomical
    implicit none
    class           (*         ), intent(inout)               :: self
    double precision            , allocatable  , dimension(:) :: rateStarFormationDisk , rateStarFormationSpheroid, &
         &                                                       densityStellar        , densityStellarDisk       , &
         &                                                       densityStellarSpheroid, densityColdGas           , &
         &                                                       densityHotHaloGas     , densityNode
    type            (hdf5Object)                              :: historyDataset        , historyGroup
    
    select type (self)
    class is (mergerTreeEvolveTimestepHistory)
      !$ call hdf5Access%set()
       historyGroup=galacticusOutputFile%openGroup('globalHistory','Global (volume averaged) history for this model.')
       if (.not.historyGroup%hasDataset('time')) then
          ! Write non-cumulative datasets on the first write.
          call historyGroup  %writeDataset  (self%time                       ,"time"             ,"Time [Gyr]"                       ,datasetReturned=historyDataset)
          call historyDataset%writeAttribute(gigaYear                        ,"unitsInSI"                                                                           )
          call historyDataset%close         (                                                                                                                       )
          call historyGroup  %writeDataset  (self%expansionFactor            ,"expansionFactor"  ,"Expansion factor []"                                             )
          call historyGroup  %writeDataset  (self%rateStarFormation          ,"rateStarFormation","Star formation rate [M⊙/Gyr/Mpc³]",datasetReturned=historyDataset)
          call historyDataset%writeAttribute(massSolar/gigaYear/megaParsec**3,"unitsInSI"                                                                           )
          call historyDataset%close         (                                                                                                                       )
       else
          ! Read existing datasets, and accumulate.
          call historyGroup%readDataset('rateStarFormationDisk'    ,rateStarFormationDisk    )
          call historyGroup%readDataset('rateStarFormationSpheroid',rateStarFormationSpheroid)
          call historyGroup%readDataset('densityStellar'           ,densityStellar           )
          call historyGroup%readDataset('densityStellarDisk'       ,densityStellarDisk       )
          call historyGroup%readDataset('densityStellarSpheroid'   ,densityStellarSpheroid   )
          call historyGroup%readDataset('densityColdGas'           ,densityColdGas           )
          call historyGroup%readDataset('densityHotHaloGas'        ,densityHotHaloGas        )
          call historyGroup%readDataset('densityNode'              ,densityNode              )
          self%rateStarFormationDisk    =self%rateStarFormationDisk    +rateStarFormationDisk
          self%rateStarFormationSpheroid=self%rateStarFormationSpheroid+rateStarFormationSpheroid
          self%densityStellar           =self%densityStellar           +densityStellar
          self%densityStellarDisk       =self%densityStellarDisk       +densityStellarDisk
          self%densityStellarSpheroid   =self%densityStellarSpheroid   +densityStellarSpheroid
          self%densityColdGas           =self%densityColdGas           +densityColdGas
          self%densityHotHaloGas        =self%densityHotHaloGas        +densityHotHaloGas
          self%densityNode              =self%densityNode              +densityNode
       end if
       ! Accumulate other datasets. (We are doing this inside an hdf5Access lock, so OpenMP threads will not conflict here.)
       call historyGroup  %writeDataset  (self%rateStarFormationDisk      ,"rateStarFormationDisk"    ,"Star formation rate in disks [M⊙/Gyr/Mpc³]"    ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/gigaYear/megaParsec**3,"unitsInSI"                                                                                                )
       call historyDataset%close         (                                                                                                                                            )
       call historyGroup  %writeDataset  (self%rateStarFormationSpheroid  ,"rateStarFormationSpheroid","Star formation rate in spheroids [M⊙/Gyr/Mpc³]",datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/gigaYear/megaParsec**3,"unitsInSI"                                                                                                )
       call historyDataset%close         (                                                                                                                                            )
       call historyGroup  %writeDataset  (self%densityStellar             ,"densityStellar"           ,"Stellar mass density [M⊙/Mpc³]"                ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/megaParsec**3         ,"unitsInSI"                                                                                                )
       call historyDataset%close         (                                                                                                                                            )
       call historyGroup  %writeDataset  (self%densityStellarDisk         ,"densityStellarDisk"       ,"Stellar mass density in disks [M⊙/Mpc³]"       ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/megaParsec**3         ,"unitsInSI"                                                                                                )
       call historyDataset%close         (                                                                                                                                            )
       call historyGroup  %writeDataset  (self%densityStellarSpheroid     ,"densityStellarSpheroid"   ,"Stellar mass density in spheroids [M⊙/Mpc³]"   ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/megaParsec**3         ,"unitsInSI"                                                                                                )
       call historyDataset%close         (                                                                                                                                            )
       call historyGroup  %writeDataset  (self%densityColdGas             ,"densityColdGas"           ,"Gas mass density [M⊙/Mpc³]"                    ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/megaParsec**3         ,"unitsInSI"                                                                                                )
       call historyDataset%close         (                                                                                                                                            )
       call historyGroup  %writeDataset  (self%densityHotHaloGas          ,"densityHotHaloGas"        ,"Hot gas mass density [M⊙/Mpc³]"                ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/megaParsec**3         ,"unitsInSI"                                                                                                )
       call historyDataset%close         (                                                                                                                                            )
       call historyGroup  %writeDataset  (self%densityNode                ,"densityNode"              ,"Node mass density [M⊙/Mpc³]"                   ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/megaParsec**3         ,"unitsInSI"                                                                                                )
       call historyDataset%close         (                                                                                                                                            )
       call historyGroup  %close         (                                                                                                                                            )
       !$ call hdf5Access %unset         (                                                                                                                                            )
    class default
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine historyWrite
