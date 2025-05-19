!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

  !+    Contributions to this file made by: Xiaolong Du, Andrew Benson.

  !!{
  Implements a merger tree operator which creates particle representations of \glc\ halos.
  !!}

  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters    , only : cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Galacticus_Nodes        , only : treeNode
  use :: HDF5                    , only : hsize_t
  use :: ISO_Varying_String      , only : varying_string
  use :: Input_Parameters        , only : inputParameters
  use :: Kind_Numbers            , only : kind_int8
  use :: Tables                  , only : table1D                  , table1DLogarithmicCSpline

  !![
  <enumeration>
   <name>selection</name>
   <description>Options for selection of nodes to particulate.</description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <visibility>private</visibility>
   <entry label="all"        />
   <entry label="hosts"      />
   <entry label="satellites" />
  </enumeration>
  !!]

  !![
  <enumeration>
   <name>particulateKernel</name>
   <description>Options for softening kernel in particulate.</description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <visibility>private</visibility>
   <entry label="delta"  />
   <entry label="plummer"/>
   <entry label="gadget" />
  </enumeration>
  !!]

  !![
  <mergerTreeOperator name="mergerTreeOperatorParticulate">
   <description>Provides a merger tree operator which create particle representations of \glc\ halos.</description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorParticulate
     !!{
     A merger tree operator which create particle representations of \glc\ halos.
     !!}
     private
     class           (cosmologyParametersClass        ), pointer :: cosmologyParameters_  => null()
     class           (cosmologyFunctionsClass         ), pointer :: cosmologyFunctions_   => null()
     class           (darkMatterHaloScaleClass        ), pointer :: darkMatterHaloScale_  => null()
     class           (darkMatterProfileDMOClass       ), pointer :: darkMatterProfileDMO_ => null()
     type            (varying_string                  )          :: outputFileName
     double precision                                            :: massParticle                   , radiusTruncateOverRadiusVirial   , &
          &                                                         timeSnapshot                   , energyDistributionPointsPerDecade, &
          &                                                         lengthSoftening                , toleranceRelativeSmoothing       , &
          &                                                         toleranceMass                  , tolerancePotential
     logical                                                     :: satelliteOffset                , nonCosmological                  , &
          &                                                         positionOffset                 , addHubbleFlow                    , &
          &                                                         haloIdToParticleType           , sampleParticleNumber             , &
          &                                                         subtractRandomOffset
     type            (enumerationParticulateKernelType)          :: kernelSoftening
     type            (enumerationSelectionType        )          :: selection
     integer         (hsize_t                         )          :: chunkSize
     integer         (kind_int8                       )          :: idMultiplier
     ! Pointer to the parameters for this task.
     type            (inputParameters                 )          :: parameters
   contains
     final     ::                        particulateDestructor
     procedure :: operatePreEvolution => particulateOperatePreEvolution
  end type mergerTreeOperatorParticulate

  interface mergerTreeOperatorParticulate
     !!{
     Constructors for the particulate merger tree operator class.
     !!}
     module procedure particulateConstructorParameters
     module procedure particulateConstructorInternal
  end interface mergerTreeOperatorParticulate

  ! Entries in the energy distribution table.
  integer                                           , parameter   :: energyDistributionTableDensity        =1
  integer                                           , parameter   :: energyDistributionTablePotential      =2
  integer                                           , parameter   :: energyDistributionTableDistribution   =3
  integer                                           , parameter   :: energyDistributionTableDensitySmoothed=4
  integer                                           , parameter   :: energyDistributionTableMass           =5

  ! Tables used in construction of distribution functions.
  class           (mergerTreeOperatorParticulate   ), pointer     :: self_
  type            (treeNode                        ), pointer     :: node_
  logical                                                         :: energyDistributionInitialized
  class           (table1D                         ), allocatable :: radiusDistribution
  type            (table1DLogarithmicCSpline       )              :: energyDistribution
  type            (enumerationParticulateKernelType)              :: softeningKernel
  double precision                                                :: energy_                                 , radiusTruncate_, &
       &                                                             height_                                 , radius_        , &
       &                                                             lengthSoftening_
  !$omp threadprivate(self_,node_,radiusDistribution,energyDistribution,energyDistributionInitialized,softeningKernel,energy_,radiusTruncate_,height_,radius_,lengthSoftening_)

contains

  function particulateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the particulate merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeOperatorParticulate   )                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    type            (varying_string                  )                :: outputFileName
    integer         (kind_int8                       )                :: idMultiplier
    double precision                                                  :: massParticle                   , radiusTruncateOverRadiusVirial   , &
         &                                                               timeSnapshot                   , energyDistributionPointsPerDecade, &
         &                                                               lengthSoftening                , toleranceRelativeSmoothing       , &
         &                                                               toleranceMass                  , tolerancePotential
    logical                                                           :: satelliteOffset                , nonCosmological                  , &
         &                                                               positionOffset                 , addHubbleFlow                    , &
         &                                                               haloIdToParticleType           , sampleParticleNumber             , &
         &                                                               subtractRandomOffset
    integer                                                           :: chunkSize
    type            (enumerationParticulateKernelType)                :: kernelSoftening_
    type            (enumerationSelectionType        )                :: selection_
    class           (cosmologyParametersClass        ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass        ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass       ), pointer       :: darkMatterProfileDMO_
    type            (inputParameters                 ), pointer       :: parametersRoot        => null()
    type            (varying_string                  )                :: selection                      , kernelSoftening

    !![
    <inputParameter>
      <name>outputFileName</name>
      <source>parameters</source>
      <description>Name of the file to which particle data should be written.</description>
    </inputParameter>
    <inputParameter>
      <name>idMultiplier</name>
      <source>parameters</source>
      <defaultValue>0_kind_int8</defaultValue>
      <description>If this parameter is greater than zero, particle IDs begin at {\normalfont \ttfamily nodeIndex\*[idMultiplier]} for each node. The multiplier should be chosen to be large enough that duplicate IDs can not occur.</description>
    </inputParameter>
    <inputParameter>
      <name>massParticle</name>
      <source>parameters</source>
      <description>Mass of particles to be used to represent halos.</description>
    </inputParameter>
    <inputParameter>
      <name>timeSnapshot</name>
      <source>parameters</source>
      <description>The time at which to snapshot the tree.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusTruncateOverRadiusVirial</name>
      <source>parameters</source>
      <description>Radius (in units of the virial radius) at which to truncate halo profiles.</description>
    </inputParameter>
    <inputParameter>
      <name>satelliteOffset</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, offset particle representations to the positions/velocities of satellites.</description>
    </inputParameter>
    <inputParameter>
      <name>positionOffset</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, offset particle representations to the positions/velocities of nodes.</description>
    </inputParameter>
    <inputParameter>
      <name>subtractRandomOffset</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, the center-of-mass positions and velocities of the host and satellite are enforced to be matched with the values specified by the offset parameters. If false, due to limited number of particles in the representations, the center-of-mass positions and velocities may deviate slightly from the specified values, i.e. have small random offsets.</description>
    </inputParameter>
    <inputParameter>
      <name>energyDistributionPointsPerDecade</name>
      <source>parameters</source>
      <defaultValue>30.0d0</defaultValue>
      <description>The number of points per decade of radius to use when building the table of energy distribution function.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelativeSmoothing</name>
      <source>parameters</source>
      <defaultValue>1.0d-7</defaultValue>
      <description>The relative tolerance to use in the integrals used in finding the smoothed density profile defined by \cite{barnes_gravitational_2012} to account for gravitational softening.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceMass</name>
      <source>parameters</source>
      <defaultValue>1.0d-8</defaultValue>
      <description>The relative tolerance to use in the integrals over the mass distribution used in finding the smoothed density profile defined by \cite{barnes_gravitational_2012} to account for gravitational softening.</description>
    </inputParameter>
    <inputParameter>
      <name>tolerancePotential</name>
      <source>parameters</source>
      <defaultValue>1.0d-9</defaultValue>
      <description>The relative tolerance to use in the integrals over the potential used in finding the smoothed density profile defined by \cite{barnes_gravitational_2012} to account for gravitational softening.</description>
    </inputParameter>
    <inputParameter>
      <name>selection</name>
      <source>parameters</source>
      <defaultValue>var_str('all')</defaultValue>
      <description>Selects the type of halo to output. Allowed options are ``{\normalfont \ttfamily all}'', ``{\normalfont \ttfamily hosts}'', and ``{\normalfont \ttfamily satellites}''.</description>
    </inputParameter>
    !!]
    selection_=enumerationSelectionEncode(char(selection),includesPrefix=.false.)
    !![
    <inputParameter>
      <name>kernelSoftening</name>
      <source>parameters</source>
      <defaultValue>var_str('plummer')</defaultValue>
      <description>Selects the softening kernel to use. Allowed options are ``{\normalfont \ttfamily plummer}'', and ``{\normalfont \ttfamily gadget}''.</description>
    </inputParameter>
    !!]
    kernelSoftening_=enumerationParticulateKernelEncode(char(kernelSoftening),includesPrefix=.false.)
    !![
    <inputParameter>
      <name>nonCosmological</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, a non-cosmological snapshot file will be created.</description>
    </inputParameter>
    <inputParameter>
      <name>addHubbleFlow</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, Hubble flow will be added to velocity offsets of halos (if applied).</description>
    </inputParameter>
    <inputParameter>
      <name>haloIdToParticleType</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, the halo ID will be used to assign particles from each halo to a different Gadget particle type group.</description>
    </inputParameter>
    <inputParameter>
      <name>sampleParticleNumber</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, the number of particles in each halo will be sampled from a Poisson distribution with the expected mean. Otherwise, the number is set equal to that expectation.</description>
    </inputParameter>
    <inputParameter>
      <name>lengthSoftening</name>
      <source>parameters</source>
      <description>The Plummer-equivalent softening length. That is, the parameter $\epsilon$ in the softening gravitational potential $\phi(r) = -\mathrm{G}m/\sqrt{r^2+\epsilon^2}$. If set to zero, softening is ignored when constructing the particle representation of the halo. For non-zero values softening is accounted for when constructing the velocity distribution following the procedure of \cite{barnes_gravitational_2012}.</description>
    </inputParameter>
    <inputParameter>
      <name>chunkSize</name>
      <source>parameters</source>
      <defaultValue>-1</defaultValue>
      <description>HDF5 dataset chunk size.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"  name="cosmologyParameters_"  source="parameters"/>
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    if (associated(parameters%parent)) then
       parametersRoot => parameters%parent
       do while (associated(parametersRoot%parent))
          parametersRoot => parametersRoot%parent
       end do
       self=mergerTreeOperatorParticulate(outputFileName,idMultiplier,massParticle,radiusTruncateOverRadiusVirial,timeSnapshot,satelliteOffset,positionOffset,subtractRandomOffset,energyDistributionPointsPerDecade,selection_,nonCosmological,addHubbleFlow,haloIdToParticleType,sampleParticleNumber,kernelSoftening_,lengthSoftening,toleranceRelativeSmoothing,toleranceMass,tolerancePotential,chunkSize,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,darkMatterProfileDMO_,parametersRoot)
    else
       self=mergerTreeOperatorParticulate(outputFileName,idMultiplier,massParticle,radiusTruncateOverRadiusVirial,timeSnapshot,satelliteOffset,positionOffset,subtractRandomOffset,energyDistributionPointsPerDecade,selection_,nonCosmological,addHubbleFlow,haloIdToParticleType,sampleParticleNumber,kernelSoftening_,lengthSoftening,toleranceRelativeSmoothing,toleranceMass,tolerancePotential,chunkSize,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,darkMatterProfileDMO_,parameters    )
    end if
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_" />
    <objectDestructor name="cosmologyFunctions_"  />
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function particulateConstructorParameters

  function particulateConstructorInternal(outputFileName,idMultiplier,massParticle,radiusTruncateOverRadiusVirial,timeSnapshot,satelliteOffset,positionOffset,subtractRandomOffset,energyDistributionPointsPerDecade,selection,nonCosmological,addHubbleFlow,haloIdToParticleType,sampleParticleNumber,kernelSoftening,lengthSoftening,toleranceRelativeSmoothing,toleranceMass,tolerancePotential,chunkSize,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,darkMatterProfileDMO_,parameters) result(self)
    !!{
    Internal constructor for the particulate merger tree operator class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (mergerTreeOperatorParticulate   )                        :: self
    type            (varying_string                  ), intent(in   )         :: outputFileName
    integer         (kind_int8                       ), intent(in   )         :: idMultiplier
    double precision                                  , intent(in   )         :: massParticle         , radiusTruncateOverRadiusVirial   , &
         &                                                                       timeSnapshot         , energyDistributionPointsPerDecade, &
         &                                                                       lengthSoftening      , toleranceRelativeSmoothing       , &
         &                                                                       toleranceMass        , tolerancePotential
    logical                                           , intent(in   )         :: satelliteOffset      , nonCosmological                  , &
         &                                                                       positionOffset       , addHubbleFlow                    , &
         &                                                                       haloIdToParticleType , sampleParticleNumber             , &
         &                                                                       subtractRandomOffset
    integer                                           , intent(in   )         :: chunkSize
    type            (enumerationParticulateKernelType), intent(in   )         :: kernelSoftening
    type            (enumerationSelectionType        ), intent(in   )         :: selection
    class           (cosmologyParametersClass        ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass         ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass        ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass       ), intent(in   ), target :: darkMatterProfileDMO_
    type            (inputParameters                 ), intent(in   ), target :: parameters
    !![
    <constructorAssign variables="outputFileName,idMultiplier,massParticle,radiusTruncateOverRadiusVirial,timeSnapshot,satelliteOffset,positionOffset,subtractRandomOffset,energyDistributionPointsPerDecade,selection,nonCosmological,addHubbleFlow,haloIdToParticleType,sampleParticleNumber,kernelSoftening,lengthSoftening,toleranceRelativeSmoothing,toleranceMass,tolerancePotential,chunkSize,*cosmologyParameters_,*cosmologyFunctions_,*darkMatterHaloScale_,*darkMatterProfileDMO_"/>
    !!]

    self%parameters=inputParameters(parameters)
    ! Validate input.
    if (.not.enumerationSelectionIsValid(selection)) call Error_Report('invalid selection type'//{introspection:location})
   return
  end function particulateConstructorInternal

  subroutine particulateDestructor(self)
    !!{
    Destructor for the particulate merger tree operator class.
    !!}
    implicit none
    type(mergerTreeOperatorParticulate), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_" />
    <objectDestructor name="self%cosmologyFunctions_"  />
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine particulateDestructor

  subroutine particulateOperatePreEvolution(self,tree)
    !!{
    Perform a particulation operation on a merger tree (i.e. create a particle representation of the tree).
    !!}
    use    :: Coordinates               , only : assignment(=)                    , coordinateCartesian                , coordinateSpherical
    use    :: Cosmology_Parameters      , only : hubbleUnitsLittleH
    use    :: Display                   , only : displayCounter                   , displayCounterClear                , verbosityLevelStandard, verbosityLevelWorking , &
         &                                       displayGreen                     , displayReset
    use    :: Galactic_Structure_Options, only : massTypeDark
    use    :: Calculations_Resets       , only : Calculations_Reset
    use    :: Error                     , only : Error_Report
    use    :: Galacticus_Nodes          , only : mergerTree                       , nodeComponentBasic                 , nodeComponentPosition , nodeComponentSatellite, &
          &                                      treeNode
    use    :: HDF5_Access               , only : hdf5Access
    use    :: IO_HDF5                   , only : hdf5Object
    use    :: ISO_Varying_String        , only : varying_string                   , var_str
    use    :: Mass_Distributions        , only : massDistributionClass
    use    :: Merger_Tree_Walkers       , only : mergerTreeWalkerAllNodes
    use    :: Node_Components           , only : Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize
    use    :: Numerical_Comparison      , only : Values_Agree
    use    :: Numerical_Constants_Math  , only : Pi
    !$ use :: OMP_Lib                   , only : OMP_Get_Thread_Num
    implicit none
    class           (mergerTreeOperatorParticulate), intent(inout) , target      :: self
    type            (mergerTree                   ), intent(inout) , target      :: tree
    type            (treeNode                     ), pointer                     :: node
    class           (nodeComponentBasic           ), pointer                     :: basic
    class           (nodeComponentPosition        ), pointer                     :: position
    class           (nodeComponentSatellite       ), pointer                     :: satellite
    class           (massDistributionClass        ), pointer                     :: massDistribution_
    double precision                               , parameter                   :: tolerance                 =1.0d-06
    double precision                               , parameter                   :: unitGadgetMass            =1.0d+10
    double precision                               , parameter                   :: unitGadgetLength          =1.0d-03
    double precision                               , parameter                   :: unitGadgetVelocity        =1.0d+00
    double precision                               , parameter                   :: distributionFunctionBuffer=0.3d0
    double precision                               , dimension(3  )              :: positionVector                       , velocityVector             , &
         &                                                                          randomDeviates                       , positionRandomOffset       , &
         &                                                                          velocityRandomOffset
    integer                                        , dimension(6  )              :: particleCounts
    double precision                               , dimension(:,:), allocatable :: particlePosition                     , particleVelocity
    integer         (kind_int8                    ), dimension(  :), allocatable :: particleIDs
    type            (mergerTreeWalkerAllNodes     )                              :: treeWalker
    double precision                                                             :: particleCountMean                    , distributionFunction       , &
         &                                                                          radiusVirial                         , radiusTruncate             , &
         &                                                                          speedPrevious                        , massTruncate               , &
         &                                                                          energy                               , radiusEnergy               , &
         &                                                                          energyPotential                      , speed                      , &
         &                                                                          speedEscape                          , distributionFunctionMaximum
    integer                                                                      :: particleCountActual                  , i                          , &
         &                                                                          j                                    , typeIndex                  , &
         &                                                                          counter
    logical                                                                      :: isNew                                , keepSample                 , &
         &                                                                          firstNode
    type            (coordinateCartesian          )                              :: positionCartesian                    , velocityCartesian
    type            (coordinateSpherical          )                              :: positionSpherical                    , velocitySpherical
    type            (hdf5Object                   )                              :: outputFile                           , header                     , &
         &                                                                          particleGroup
    type            (varying_string               )                              :: message
    character       (len=13                       )                              :: label
    character       (len= 9                       )                              :: groupName

    ! Open the HDF5 file for output.
    !$ call hdf5Access%set  ()
    outputFile=hdf5Object(char(self%outputFileName),overWrite=.true.,readOnly=.false.)
    ! Create the header.
    header=outputFile%openGroup('Header','Group containing Gadget metadata.')
    ! Particle properties.
    call header%writeAttribute(       1                                    ,'NumFilesPerSnapshot'   )
    call header%writeAttribute(spread(0                               ,1,6),'NumPart_ThisFile'      )
    call header%writeAttribute(spread(0                               ,1,6),'NumPart_Total_HighWord')
    call header%writeAttribute(spread(0                               ,1,6),'NumPart_Total'         )
    call header%writeAttribute(spread(self%massParticle/unitGadgetMass,1,6),'MassTable'             )
    ! Time.
    call header%writeAttribute(                                                                                              self%timeSnapshot  ,'Time'    )
    call header%writeAttribute(self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(self%timeSnapshot)),'Redshift')
    ! Cosmology.
    if (self%nonCosmological) then
       call header%writeAttribute(1.0d0                                                        ,'HubbleParam')
       call header%writeAttribute(0.0d0                                                        ,'Omega0'     )
       call header%writeAttribute(0.0d0                                                        ,'OmegaLambda')
    else
       call header%writeAttribute(self%cosmologyParameters_%HubbleConstant (hubbleUnitsLittleH),'HubbleParam')
       call header%writeAttribute(self%cosmologyParameters_%OmegaMatter    (                  ),'Omega0'     )
       call header%writeAttribute(self%cosmologyParameters_%OmegaDarkEnergy(                  ),'OmegaLambda')
    end if
    call header%writeAttribute(0.0d0,'BoxSize'         )
    ! Flags.
    call header%writeAttribute(0    ,'Flag_Cooling'    )
    call header%writeAttribute(0    ,'Flag_Sfr'        )
    call header%writeAttribute(0    ,'Flag_Feedback'   )
    call header%writeAttribute(0    ,'Flag_StellarAge' )
    call header%writeAttribute(0    ,'FlagMetals'      )
    call header%writeAttribute(0    ,'Flag_Entropy_ICs')
    !$ call hdf5Access%unset()
    ! Iterate over nodes.
    firstNode =.true.
    treeWalker=mergerTreeWalkerAllNodes(tree,spanForest=.true.)
    do while (treeWalker%next(node))
       call Calculations_Reset(node)
       ! Check if the node exists at the snapshot time.
       basic => node%basic()
       if     (                                                               &
            &   Values_Agree(basic%time(),self%timeSnapshot,relTol=tolerance) &
            &  .and.                                                          &
            &   (                                                             &
            &            self%selection == selectionAll                       &
            &    .or.                                                         &
            &     (                                                           &
            &            self%selection == selectionHosts                     &
            &      .and.                                                      &
            &       .not.node%isSatellite()                                   &
            &     )                                                           &
            &    .or.                                                         &
            &     (                                                           &
            &            self%selection == selectionSatellites                &
            &      .and.                                                      &
            &            node%isSatellite()                                   &
            &     )                                                           &
            &   )                                                             &
            & ) then
          ! Force the energy distribution tables to be rebuilt.
          energyDistributionInitialized=.false.
          ! Determine the virial radius.
          radiusVirial  =+self%darkMatterHaloScale_%radiusVirial(node)
          ! Determine the truncation radius.
          radiusTruncate=+self%radiusTruncateOverRadiusVirial &
               &         *radiusVirial
          ! Determine the mass within the truncation radius.
          massDistribution_ => node             %massDistribution    (massType=massTypeDark  )
          massTruncate      =  massDistribution_%massEnclosedBySphere(         radiusTruncate)
          !![
	  <objectDestructor name="massDistribution_"/>
	  !!]
          ! Determine the mean number of particles required to represent this node.
          particleCountMean   =+massTruncate      &
               &               /self%massParticle
          ! Determine the actual number of particles to use to represent the node.
          if (self%sampleParticleNumber) then
             particleCountActual=tree%randomNumberGenerator_%poissonSample(particleCountMean)
          else
             particleCountActual=nint(particleCountMean)
          end if
          ! Allocate space for particle data.
          allocate(particlePosition(3,particleCountActual))
          allocate(particleVelocity(3,particleCountActual))
          allocate(particleIDs     (particleCountActual))
          ! Get required components.
          position  => node%position ()
          satellite => node%satellite()
          ! Set pointers to module-scope variables.
          node_            => node
          radiusTruncate_  =  radiusTruncate
          lengthSoftening_ =  self%lengthSoftening
          softeningKernel  =  self%kernelSoftening
          ! Iterate over particles.
          isNew               =.true.
          counter             =0
          positionRandomOffset=0.0d0
          velocityRandomOffset=0.0d0
          !$omp parallel private(i,j,positionSpherical,positionCartesian,velocitySpherical,velocityCartesian,energy,energyPotential,speed,speedEscape,speedPrevious,distributionFunction,distributionFunctionMaximum,keepSample,radiusEnergy,positionVector,velocityVector,randomDeviates) copyin(node_,radiusTruncate_,lengthSoftening_,softeningKernel)
          call Node_Components_Thread_Initialize(self%parameters)
          allocate(self_,mold=self)
          !$omp critical(mergerTreeOperatorsParticulateDeepCopy)
          !![
          <deepCopyReset variables="self"/>
          <deepCopy source="self" destination="self_"/>
          <deepCopyFinalize variables="self_"/>
          !!]
          !$omp end critical(mergerTreeOperatorsParticulateDeepCopy)
          !$omp do reduction(+: positionRandomOffset, velocityRandomOffset)
          do i=1,particleCountActual
             !$ if (OMP_Get_Thread_Num() == 0) then
             call displayCounter(max(1,int(100.0d0*dble(counter)/dble(particleCountActual))),isNew=isNew,verbosity=verbosityLevelStandard)
             isNew=.false.
             !$ end if
             !$omp atomic
             counter=counter+1
             ! Sample particle positions from the halo density distribution. Currently, we assume that halos are spherically
             ! symmetric.
             randomDeviates=-1.0d0
             !$omp critical (mergerTreeOperatorParticulateSample)
             do j=1,3
                ! Ensure that the third random deviate (the enclosed mass fraction) is never precisely zero.
                do while (                                           &
                     &     (j == 3 .and. randomDeviates(j) <= 0.0d0) &
                     &    .or.                                       &
                     &     (j /= 3 .and. randomDeviates(j) <  0.0d0) &
                     &   )
                   randomDeviates(j)=tree%randomNumberGenerator_%uniformSample()
                end do
             end do
             !$omp end critical (mergerTreeOperatorParticulateSample)
             call positionSpherical%  phiSet(     2.0d0*Pi       *randomDeviates(1)       )
             call positionSpherical%thetaSet(acos(2.0d0          *randomDeviates(2)-1.0d0))
             call positionSpherical%    rSet(     radiusTruncate_*randomDeviates(3)       )
             ! Get the corresponding cartesian coordinates.
             positionCartesian=positionSpherical
             ! Construct the energy distribution function encompassing this radius.
             call particulateTabulateEnergyDistribution(positionSpherical%r(),self%energyDistributionPointsPerDecade)
             ! Find the potential energy and escape speed (actually the speed to reach the truncation radius) at this radius.
             energyPotential=+energyDistribution%interpolate(                                        &
                  &                                                positionSpherical%r()           , &
                  &                                          table=energyDistributionTablePotential  &
                  &                                         )
             speedEscape    =+sqrt(                      &
                  &                +2.0d0                &
                  &                *max(                 &
                  &                     energyPotential, &
                  &                     0.0d0            &
                  &                    )                 &
                  &               )
             ! Estimate the maximum of the speed distribution function.
             distributionFunctionMaximum=+0.0d0
             speedPrevious              =-1.0d0
             do j=1,energyDistribution%size()
                ! Find the energy at this radius in the tabulated profile.
                energy=energyDistribution%y(j,table=energyDistributionTablePotential)
                ! Only consider points with positive kinetic energy and speed below the escape speed (i.e. the speed required
                ! to reach the truncation radius).
                if (-energy+energyPotential > 0.0d0 .and. speedPrevious < speedEscape) then
                   ! Compute the speed corresponding to this total energy
                   speed               =+sqrt(                   &
                        &                     +2.0d0             &
                        &                     *(                 &
                        &                       -energy          &
                        &                       +energyPotential &
                        &                      )                 &
                        &                     )
                   speedPrevious       =+speed
                   ! Find the distribution function. This is v²f(E). The tabulated function is the anti-derivative of f(E), so
                   ! we need to take the derivative with respect to energy. We do this by taking the derivative with respect to
                   ! radius and dividing by dE/dr.
                   distributionFunction       =+speed**2                                                                                                  &
                        &                      *energyDistribution%interpolateGradient(energyDistribution%x(j),table=energyDistributionTableDistribution) &
                        &                      /energyDistribution%interpolateGradient(energyDistribution%x(j),table=energyDistributionTablePotential   )
                   ! Track the maximum of this distribution function.
                   distributionFunctionMaximum=+max(                             &
                        &                           distributionFunctionMaximum, &
                        &                           distributionFunction         &
                        &                          )
                end if
             end do
             ! Add some buffer to the distribution function maximum.
             distributionFunctionMaximum=+distributionFunctionMaximum  &
                  &                      *(                            &
                  &                        +1.0d0                      &
                  &                        +distributionFunctionBuffer &
                  &                       )
             ! Draw a speed from the distribution function between zero and the escape speed (to the truncation radius) using
             ! rejection sampling.
             keepSample=.false.
             do while (.not.keepSample)
                ! Draw a speed uniformly at random between zero and the escape velocity.
                !$omp critical (mergerTreeOperatorParticulateSample)
                speed               =+tree%randomNumberGenerator_%uniformSample() &
                     &               *speedEscape
                !$omp end critical (mergerTreeOperatorParticulateSample)
                energy              =+energyPotential    &
                     &               -0.5d0              &
                     &               *speed          **2
                !$omp critical (mergerTreeOperatorParticulateRadius)
                radiusEnergy        =+radiusDistribution%interpolate        (                                           &
                     &                                                             energy                             , &
                     &                                                       table=energyDistributionTablePotential     &
                     &                                                      )
                !$omp end critical (mergerTreeOperatorParticulateRadius)
                distributionFunction=+speed**2                                                                          &
                     &               *energyDistribution%interpolateGradient(                                           &
                     &                                                             radiusEnergy                       , &
                     &                                                       table=energyDistributionTableDistribution  &
                     &                                                      )                                           &
                     &               /energyDistribution%interpolateGradient(                                           &
                     &                                                             radiusEnergy                       , &
                     &                                                       table=energyDistributionTablePotential     &
                     &                                                      )
                if (distributionFunction > distributionFunctionMaximum) then
                   write (label,'(e12.6)') distributionFunction
                   message='distribution function ['//trim(label)//'] exceeds estimated maximum ['
                   write (label,'(e12.6)') distributionFunctionMaximum
                   message=message//trim(label)//']'//char(10)
                   message=message//displayGreen()//'HELP:'//displayReset()//' the issue is probably caused by an inaccurate estimation of the maximum of the distribution function from tabulated values. To resolve this issue, increase the parameter [energyDistributionPointsPerDecade].'//char(10)
                   call Error_Report(message//{introspection:location})
                end if
                !$omp critical (mergerTreeOperatorParticulateSample)
                keepSample=  +tree%randomNumberGenerator_%uniformSample() &
                     &      <                                             &
                     &       +distributionFunction                        &
                     &       /distributionFunctionMaximum
                !$omp end critical (mergerTreeOperatorParticulateSample)
             end do
             ! Choose a velocity vector in spherical coordinates with velocity chosen to give the required kinetic energy.
             !$omp critical (mergerTreeOperatorParticulateSample)
             call velocitySpherical%  phiSet(     2.0d0*Pi*tree%randomNumberGenerator_%uniformSample()       )
             call velocitySpherical%thetaSet(acos(2.0d0   *tree%randomNumberGenerator_%uniformSample()-1.0d0))
             !$omp end critical (mergerTreeOperatorParticulateSample)
             call velocitySpherical%    rSet(speed                                                           )
             ! Get the corresponding cartesian coordinates.
             velocityCartesian=velocitySpherical
             ! Accumulate the particle position and velocity.
             if (self%subtractRandomOffset) then
                positionRandomOffset(1)=+positionRandomOffset(1)+positionCartesian%x()
                positionRandomOffset(2)=+positionRandomOffset(2)+positionCartesian%y()
                positionRandomOffset(3)=+positionRandomOffset(3)+positionCartesian%z()
                velocityRandomOffset(1)=+velocityRandomOffset(1)+velocityCartesian%x()
                velocityRandomOffset(2)=+velocityRandomOffset(2)+velocityCartesian%y()
                velocityRandomOffset(3)=+velocityRandomOffset(3)+velocityCartesian%z()
             end if
             ! Offset position and velocity to the position and velocity of the node.
             if (self%positionOffset) then
                positionVector=position%position()
                velocityVector=position%velocity()
                if (self%addHubbleFlow) velocityVector=+velocityVector                                                          &
                     &                                 +positionVector                                                          &
                     &                                 *self%cosmologyFunctions_%hubbleParameterEpochal(time=self%timeSnapshot)
                call positionCartesian%xSet(positionCartesian%x()+positionVector(1))
                call positionCartesian%ySet(positionCartesian%y()+positionVector(2))
                call positionCartesian%zSet(positionCartesian%z()+positionVector(3))
                call velocityCartesian%xSet(velocityCartesian%x()+velocityVector(1))
                call velocityCartesian%ySet(velocityCartesian%y()+velocityVector(2))
                call velocityCartesian%zSet(velocityCartesian%z()+velocityVector(3))
             end if
             ! Offset position and velocity to position and velocity of satellite.
             if (self%satelliteOffset) then
                positionVector=satellite%position()
                velocityVector=satellite%velocity()
                call positionCartesian%xSet(positionCartesian%x()+positionVector(1))
                call positionCartesian%ySet(positionCartesian%y()+positionVector(2))
                call positionCartesian%zSet(positionCartesian%z()+positionVector(3))
                call velocityCartesian%xSet(velocityCartesian%x()+velocityVector(1))
                call velocityCartesian%ySet(velocityCartesian%y()+velocityVector(2))
                call velocityCartesian%zSet(velocityCartesian%z()+velocityVector(3))
             end if
             ! Store to particle data.
             particlePosition(:,i)=[positionCartesian%x(),positionCartesian%y(),positionCartesian%z()]
             particleVelocity(:,i)=[velocityCartesian%x(),velocityCartesian%y(),velocityCartesian%z()]
             particleIDs     (  i)=i-1
          end do
          !$omp end do
          call Node_Components_Thread_Uninitialize()
          !![
          <objectDestructor name="self_"/>
          !!]
          !$omp end parallel
          call displayCounterClear(verbosity=verbosityLevelWorking)
          ! Subtract random offsets of the center-of-mass position and velocity.
          positionRandomOffset = positionRandomOffset/particleCountActual
          velocityRandomOffset = velocityRandomOffset/particleCountActual
          !$omp parallel do
          do i=1,particleCountActual
             particlePosition(:,i) = particlePosition(:,i)-positionRandomOffset
             particleVelocity(:,i) = particleVelocity(:,i)-velocityRandomOffset
          end do
          !$omp end parallel do
          ! Perform unit conversion.
          particlePosition=particlePosition/unitGadgetLength
          particleVelocity=particleVelocity/unitGadgetVelocity
          ! Accumulate the particle data to file.
          !$ call hdf5Access%set  ()
          outputFile=hdf5Object(char(self%outputFileName),overWrite=.false.,readOnly=.false.,objectsOverwritable=.true.)
          ! Get current count of particles in file.
          header=outputFile%openGroup('Header','Group containing Gadget metadata.')
          call header%readAttributeStatic('NumPart_Total',particleCounts)
          ! Write particle data.
          if (self%haloIdToParticleType) then
             write (groupName,'(a,i1)') 'PartType',node%index()
             typeIndex=int(node%index())+1
          else
             groupName='PartType1'
             typeIndex=2
          end if
          ! Offset particle IDs.
          if (self%idMultiplier > 0) then
             particleIDs=particleIDs+node%index()*self%idMultiplier
          else
             particleIDs=particleIDs+particleCounts(typeIndex)
          end if
          if (.not.firstNode.and.self%chunkSize == -1)                                                                                                        &
               & call Error_Report(                                                                                                                           &
               &                   var_str('can not write multiple halos to output with chunksize=-1')//char(10)//                                            &
               &                   displayGreen()//' HELP: '//displayReset()//                                                                                &
               &                   ' set <chunkSize value="N"/> where N is a non-zero value in the particulate merger tree operator in your parameter file'// &
               &                   {introspection:location}                                                                                                   &
               &                  )
          particleGroup=outputFile%openGroup(groupName,'Group containing particle data for halos',chunkSize=self%chunkSize)
          call particleGroup%writeDataset(particlePosition,'Coordinates','Particle coordinates',appendTo=self%chunkSize /= -1,appendDimension=2)
          call particleGroup%writeDataset(particleVelocity,'Velocities' ,'Particle velocities' ,appendTo=self%chunkSize /= -1,appendDimension=2)
          call particleGroup%writeDataset(particleIDs     ,'ParticleIDs','Particle IDs'        ,appendTo=self%chunkSize /= -1                  )
          firstNode=.false.
          deallocate(particlePosition)
          deallocate(particleVelocity)
          deallocate(particleIDs     )
          ! Update particle counts.
          particleCounts(typeIndex)=particleCounts(typeIndex)+particleCountActual
          call header%writeAttribute(particleCounts,'NumPart_ThisFile')
          call header%writeAttribute(particleCounts,'NumPart_Total'   )
          !$ call hdf5Access%unset()
       end if
    end do
    return
  end subroutine particulateOperatePreEvolution

  subroutine particulateTabulateEnergyDistribution(radius,energyDistributionPointsPerDecade)
    !!{
    Construct the energy distribution function assuming a spherical dark matter halo with
    isotropic velocity dispersion. We solve Eddington's formula
    \citep[][eqn. 4.43a]{binney_galactic_2008}.
    \begin{equation}
     f(E) = {1 \over \sqrt{8} \pi^2} {\mathrm{d}\over \mathrm{d}E} \int_E^0 {\mathrm{d}\Phi \over \sqrt{\Phi-E}} {\mathrm{d}\rho \over \mathrm{d} \Phi}.
    \end{equation}
    In practice, we tabulate:
    \begin{equation}
     F(E) = \int_E^0 {\mathrm{d}\Phi \over \sqrt{\Phi-E}} {\mathrm{d}\rho \over \mathrm{d} \Phi},
    \end{equation}
    which we can then take the derivative of numerically to obtain the distribution function.
    !!}
    use :: Coordinates          , only : coordinateSpherical  , assignment(=)
    use :: Error                , only : Error_Report
    use :: Galacticus_Nodes     , only : nodeComponentBasic
    use :: Mass_Distributions   , only : massDistributionClass
    use :: Numerical_Integration, only : integrator
    use :: Table_Labels         , only : extrapolationTypeFix
    implicit none
    double precision                       , intent(in   ) :: radius
    double precision                       , intent(in   ) :: energyDistributionPointsPerDecade
    class           (nodeComponentBasic   ), pointer       :: basic
    class           (massDistributionClass), pointer       :: massDistribution_
    double precision                       , parameter     :: toleranceTabulation                      =1.0d-6
    double precision                       , parameter     :: toleranceGradient                        =1.0d-6
    ! The largest (absolute) logarithmic gradient dlog(Φ)/dlog(r) at which it is acceptable to have a non-monotonic distribution
    ! function. This allows for numerical inaccuracies that arise in cored density profiles where the central potential has the
    ! form Φ(r) = Φ₀ + k r², such that the potential is very weakly dependent on r at small radii.
    double precision                       , parameter     :: derivativeLogarithmicPotentialTolerance  =1.0d-5
    double precision                                       :: radiusMinimum                                   , derivativeLogarithmicPotential           , &
         &                                                    particulateSmoothingIntegrationRangeLower       , particulateSmoothingIntegrationRangeUpper, &
         &                                                    radiusFactorAsymptote                           , integralAsymptotic                       , &
         &                                                    gradientDensityPotentialLower                   , gradientDensityPotentialUpper            , &
         &                                                    densitySmoothedIntegralLower                    , densitySmoothedIntegralUpper
    logical                                                :: tableRebuild
    integer                                                :: i                                               , j
    integer                                                :: radiusCount
    type            (coordinateSpherical  )                :: coordinates                                     , coordinatesTruncate
    type            (integrator           )                :: intergatorSmoothingZ                            , integratorMass                           , &
         &                                                    integratorPotential                             , integratorEddington

    ! Determine the minimum of the given radius and some small fraction of the virial radius.
    basic             =>                                                      node_%basic       ()
    massDistribution_ => self_     %darkMatterProfileDMO_%get                (node_               )
    radiusMinimum     =  min(                                                                        &
         &                   +0.5d0*                      radius                                   , &
         &                   +      massDistribution_    %radiusEnclosingMass(self_%massParticle  )  &
         &                  )
    ! Rebuild the density vs. potential table to have sufficient range if necessary.
    if (energyDistributionInitialized) then
       tableRebuild=(radiusMinimum < (1.0d0-toleranceTabulation)*energyDistribution%x(1))
       if (tableRebuild) call energyDistribution%destroy()
    else
       tableRebuild                 =.true.
       energyDistributionInitialized=.true.
    end if
    if (tableRebuild) then
       ! Build tables of potential and density.
       radiusCount=int(energyDistributionPointsPerDecade*log10(radiusTruncate_/radiusMinimum))+1
       call energyDistribution%create(                                                                &
            &                                           +radiusMinimum                              , &
            &                                           +radiusTruncate_                            , &
            &                                            radiusCount                                , &
            &                         tableCount       = 5                                          , &
            &                         extrapolationType= [extrapolationTypeFix,extrapolationTypeFix]  &
            &                        )
       do i=1,radiusCount
          radius_    =energyDistribution%x(i)
          coordinates=[radius_,0.0d0,0.0d0]
          call energyDistribution%populate(                                                                                         &
               &                                         +massDistribution_%density  (coordinates)                                , &
               &                                         i                                                                        , &
               &                           table        =energyDistributionTableDensity                                           , &
               &                           computeSpline=i==radiusCount                                                             &
               &                          )
          select case (softeningKernel%ID)
          case (particulateKernelDelta%ID)
             ! No softening is applied, so use the actual density and potential.
             coordinatesTruncate=[radiusTruncate_,0.0d0,0.0d0]
             call energyDistribution%populate(                                                                                      &
                  &                                         energyDistribution%y(i,table=energyDistributionTableDensity)          , &
                  &                                                              i                                                , &
                  &                           table        =energyDistributionTableDensitySmoothed                                , &
                  &                           computeSpline=i==radiusCount                                                          &
                  &                          )
             call energyDistribution%populate(                                                                                      &
                  &                                         massDistribution_%potentialDifference(coordinatesTruncate,coordinates), &
                  &                                         i                                                                     , &
                  &                           table        =energyDistributionTablePotential                                      , &
                  &                           computeSpline=i==radiusCount                                                          &
                  &                          )
          case default
             ! Compute potential from a density field smoothed by the density distribution corresponding to the softened
             ! potential.
             particulateSmoothingIntegrationRangeUpper=+radiusTruncate_-radius_
             particulateSmoothingIntegrationRangeLower=-radiusTruncate_-radius_
             select case (softeningKernel%ID)
             case (particulateKernelGadget %ID)
                particulateSmoothingIntegrationRangeUpper=min(particulateSmoothingIntegrationRangeUpper,+2.8d0*lengthSoftening_)
                particulateSmoothingIntegrationRangeLower=max(particulateSmoothingIntegrationRangeLower,-2.8d0*lengthSoftening_)
             end select
             ! The integral here is performed in two parts - below and above the current radius where the integrand has a
             ! discontinuous gradient - for improved speed and stability.
             intergatorSmoothingZ=integrator(particulateSmoothingIntegrandZ,toleranceRelative=self_%toleranceRelativeSmoothing)
             if (particulateSmoothingIntegrationRangeLower < radius_) then
                densitySmoothedIntegralLower=intergatorSmoothingZ%integrate(                                               &
                     &                                                          particulateSmoothingIntegrationRangeLower, &
                     &                                                      min(                                           &
                     &                                                          radius_                                  , &
                     &                                                          particulateSmoothingIntegrationRangeUpper  &
                     &                                                         )                                           &
                     &                                                     )
             else
                densitySmoothedIntegralLower=0.0d0
             end if
             if (particulateSmoothingIntegrationRangeUpper > radius_) then
                densitySmoothedIntegralUpper=intergatorSmoothingZ%integrate(                                               &
                     &                                                      max(                                           &
                     &                                                          radius_                                  , &
                     &                                                          particulateSmoothingIntegrationRangeLower  &
                     &                                                         )                                         , &
                     &                                                          particulateSmoothingIntegrationRangeUpper  &
                     &                                                      )
             else
                densitySmoothedIntegralUpper=0.0d0
             end if
             call energyDistribution%populate(                                                       &
                  &                                         +densitySmoothedIntegralLower            &
                  &                                         +densitySmoothedIntegralUpper          , &
                  &                                          i                                     , &
                  &                           table        = energyDistributionTableDensitySmoothed, &
                  &                           computeSpline= i==radiusCount                          &
                  &                          )
          end select
       end do
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
       ! If necessary, compute the potential from the smoothed density profile.
       if (softeningKernel /= particulateKernelDelta) then
          integratorMass=integrator(particulateMassIntegrand,toleranceRelative=self_%toleranceMass)
          do i=1,radiusCount
             radius_=energyDistribution%x(i)
             if     (                                                                         &
                  &    i                                                                      &
                  &   >                                                                       &
                  &    1                                                                      &
                  &  .and.                                                                    &
                  &    energyDistribution%y(i  ,table=energyDistributionTableDensitySmoothed) &
                  &   >                                                                       &
                  &    energyDistribution%y(i-1,table=energyDistributionTableDensitySmoothed) &
                  & ) call energyDistribution%populate(energyDistribution%y(i-1,table=energyDistributionTableDensitySmoothed),i,table=energyDistributionTableDensitySmoothed)
             call energyDistribution%populate(                                                  &
                  &                                         +integratorMass%integrate(          &
                  &                                                                   +0.0d0  , &
                  &                                                                   +radius_  &
                  &                                                                  )        , &
                  &                                         i                                 , &
                  &                           table        =energyDistributionTableMass       , &
                  &                           computeSpline=i==radiusCount                      &
                  &                          )
          end do
          integratorPotential=integrator(particulatePotentialIntegrand,toleranceRelative=self_%tolerancePotential)
          do i=1,radiusCount
             radius_=energyDistribution%x(i)
             call energyDistribution%populate(                                                             &
                  &                                         integratorPotential%integrate(                 &
                  &                                                                       radius_        , &
                  &                                                                       radiusTruncate_  &
                  &                                                                      )               , &
                  &                                         i                                            , &
                  &                           table        =energyDistributionTablePotential             , &
                  &                           computeSpline=i==radiusCount                                 &
                  &                          )
          end do
       end if
       ! Evaluate the integral in Eddington's formula. Ignore the normalization as we will simply use rejection sampling to draw
       ! from this distribution.
       call energyDistribution%populate(0.0d0,radiusCount,table=energyDistributionTableDistribution,computeSpline=.false.)
       integratorEddington=integrator(particulateEddingtonIntegrand,toleranceRelative=1.0d-5)
       do i=1,radiusCount-1
          radius_=energyDistribution%x(i)
          energy_=energyDistribution%y(i,table=energyDistributionTablePotential)
          ! As the integrand diverges at radius_ we evaluate the integral analytically in a small region close to that radius. The
          ! size of this region is chosen such that the gradient of density with respect to potential changes very little over that range,
          ! ensuring that the approximation made in the analytic integral is valid.
          gradientDensityPotentialLower=+energyDistribution%interpolateGradient(radius_,table=energyDistributionTableDensity  ) &
               &                        /energyDistribution%interpolateGradient(radius_,table=energyDistributionTablePotential)
          radiusFactorAsymptote        =+energyDistribution%x(i+1) &
               &                        /energyDistribution%x(i  )
          ! Iteratively reduce the size of the region to be integrated analytically until the density vs. potential gradient is
          ! sufficiently constant across that range.
          do while (.true.)
             gradientDensityPotentialUpper=+energyDistribution%interpolateGradient(radius_*radiusFactorAsymptote,table=energyDistributionTableDensity  ) &
                  &                        /energyDistribution%interpolateGradient(radius_*radiusFactorAsymptote,table=energyDistributionTablePotential)
             if     (                                                                                           &
                  &                            abs(gradientDensityPotentialUpper-gradientDensityPotentialLower) &
                  &  <                                                                                          &
                  &   +0.5d0*toleranceGradient*abs(gradientDensityPotentialUpper+gradientDensityPotentialLower) &
                  & ) exit
             radiusFactorAsymptote=sqrt(radiusFactorAsymptote)
          end do
          ! Evaluate the asymptotic part of the integral analytically.
          integralAsymptotic=+2.0d0                                                                                                      &
               &             *gradientDensityPotentialLower                                                                              &
               &             *sqrt(                                                                                                      &
               &                   +energyDistribution%y          (i                            ,table=energyDistributionTablePotential) &
               &                   -energyDistribution%interpolate(radius_*radiusFactorAsymptote,table=energyDistributionTablePotential) &
               &                  )
          ! Evaluate the remainder of the integral numerically and add on the asymptotic part.
          call energyDistribution%populate(                                                                                                     &
               &                                                    +integratorEddington%integrate(                                             &
               &                                                                                   +log(radiusTruncate_                      ), &
               &                                                                                   +log(radius_        *radiusFactorAsymptote)  &
               &                                                                                  )                                             &
               &                                                    +integralAsymptotic                                                       , &
               &                                                    i                                                                         , &
               &                                      table        =energyDistributionTableDistribution                                       , &
               &                                      computeSpline=i==radiusCount-1                                                            &
               &                                     )
          ! Check that the distribution function is monotonically increasing.
          if     (                                                                      &
               &    i                                                                   &
               &   >                                                                    &
               &    1                                                                   &
               &  .and.                                                                 &
               &    energyDistribution%y(i  ,table=energyDistributionTableDistribution) &
               &   >                                                                    &
               &    energyDistribution%y(i-1,table=energyDistributionTableDistribution) &
               & ) then
             ! Find the logarithmic derivative of potential with radius at this point.
             derivativeLogarithmicPotential=+abs(                                                                                                                                    &
                  &                              +(+energyDistribution%y(i,table=energyDistributionTablePotential)-energyDistribution%y(i-1,table=energyDistributionTablePotential)) &
                  &                              /(+energyDistribution%x(i                                       )-energyDistribution%x(i-1                                       )) &
                  &                             )                                                                                                                                    &
                  &                         /       energyDistribution%y(i,table=energyDistributionTablePotential)                                                                   &
                  &                         *       energyDistribution%x(i                                       )
             ! For profiles where the potential asymptotes to a finite value at zero radius we don't need to tabulate to
             ! arbitrarily small radii, just to radii small enough that we have approximately reached this asymptotic value. (This
             ! is useful to avoid numerical inaccuracies in the regime of very small radii where the potential is almost
             ! independent of radius.)
             if (derivativeLogarithmicPotential < derivativeLogarithmicPotentialTolerance) then
                do j=i-1,1,-1
                   call energyDistribution%populate(energyDistribution%y(j+1,table=energyDistributionTablePotential)*(1.0d0+epsilon(0.0d0)),j,table=energyDistributionTableDistribution,computeSpline=i==radiusCount-1)
                end do
             else
                call Error_Report('unphysical distribution function'//{introspection:location})
             end if
          end if
       end do
       ! Construct a reversed (radius vs. potential function) table.
       call energyDistribution%reverse(radiusDistribution,table=energyDistributionTablePotential)
       ! Record that the distribution is initialized.
       energyDistributionInitialized=.true.
    end if
    return
  end subroutine particulateTabulateEnergyDistribution

  double precision function particulateEddingtonIntegrand(radiusLogarithmic)
    !!{
    The integrand appearing in Eddington's formula for the distribution function.
    !!}
    implicit none
    double precision, intent(in   ) :: radiusLogarithmic
    double precision                :: potential        , radius

    radius=exp(radiusLogarithmic)
    potential=energyDistribution%interpolate(radius,table=energyDistributionTablePotential)
    if (potential < energy_) then
       particulateEddingtonIntegrand=radius*energyDistribution%interpolateGradient(radius,table=energyDistributionTableDensity)/sqrt(energy_-potential)
    else
       particulateEddingtonIntegrand=0.0d0
    end if
    return
  end function particulateEddingtonIntegrand

  double precision function particulateSmoothingIntegrandZ(height)
    !!{
    The integrand over cylindrical coordinate $z$ used in finding the smoothed density profile defined by
    \cite{barnes_gravitational_2012} to account for gravitational softening.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    double precision            , intent(in   ) :: height
    type            (integrator)                :: integrator_
    double precision                            :: radiusMaximum, heightOffset, &
         &                                         argumentSqrt

    heightOffset=height-radius_ ! Height in the profile (not the kernel).
    argumentSqrt=+radiusTruncate_**2-heightOffset**2
    if (argumentSqrt <= 0.0d0) then
       particulateSmoothingIntegrandZ=0.0d0
    else
       radiusMaximum    =sqrt(argumentSqrt)
       height_=     height
       select case (softeningKernel%ID)
       case (particulateKernelGadget %ID)
          radiusMaximum=min(radiusMaximum,sqrt(+(2.8d0*lengthSoftening_)**2-height_**2))
       end select
       integrator_                   =integrator           (particulateSmoothingIntegrandR,toleranceRelative=self_%toleranceRelativeSmoothing)
       particulateSmoothingIntegrandZ=integrator_%integrate(0.0d0                         ,                        radiusMaximum             )
    end if
    return
  end function particulateSmoothingIntegrandZ

  double precision function particulateSmoothingIntegrandR(radiusCylindrical)
    !!{
    The integrand over cylindrical coordinate $R$ used in finding the smoothed density profile defined by
    \cite{barnes_gravitational_2012} to account for gravitational softening.
    !!}
    use :: Coordinates       , only : coordinateSpherical  , assignment(=)
    use :: Mass_Distributions, only : massDistributionClass
    implicit none
    double precision                       , intent(in   ) :: radiusCylindrical
    class           (massDistributionClass), pointer       :: massDistribution_
    double precision                                       :: radiusSplineKernel, lengthSplineKernel
    type            (coordinateSpherical  )                :: coordinates
    
    massDistribution_ => self_%darkMatterProfileDMO_%get(node_)
    coordinates       =  [                             &
         &                +sqrt(                       &
         &                      +radiusCylindrical**2  &
         &                      +(                     &
         &                        +height_             &
         &                        -radius_             &
         &                       )                **2  &
         &                     )                     , &
         &                +0.0d0                     , &
         &                +0.0d0                       &
         &               ]
    particulateSmoothingIntegrandR=+                  radiusCylindrical              &
         &                         *massDistribution_%density          (coordinates)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Apply the softening kernel density distribution.
    select case (softeningKernel%ID)
    case (particulateKernelPlummer%ID)
       particulateSmoothingIntegrandR=+particulateSmoothingIntegrandR &
            &                         *1.5d0                          &
            &                         *lengthSoftening_  **2          &
            &                         /(                              &
            &                          +radiusCylindrical**2          &
            &                          +height_          **2          &
            &                          +lengthSoftening_ **2          &
            &                         )**2.5d0
    case (particulateKernelGadget %ID)
       ! Compute the radius in units of the spline kernel length scale.
       lengthSplineKernel=+2.8d0                      &
            &             *lengthSoftening_
       radiusSplineKernel=+sqrt(                      &
            &                   +radiusCylindrical**2 &
            &                   +height_          **2 &
            &                  )                      &
            &             /lengthSplineKernel
       if (radiusSplineKernel <= 0.5d0) then
          particulateSmoothingIntegrandR=+particulateSmoothingIntegrandR &
               &                         *16.0d0                         &
               &                         /lengthSplineKernel**3          &
               &                         *(                              &
               &                           +1.0d0                        &
               &                           -6.0d0*radiusSplineKernel**2  &
               &                           +6.0d0*radiusSplineKernel**3  &
               &                          )
       else if (radiusSplineKernel <= 1.0d0) then
          particulateSmoothingIntegrandR=+particulateSmoothingIntegrandR &
               &                         *32.0d0                         &
               &                         /lengthSplineKernel**3          &
               &                         *(                              &
               &                           +1.0d0                        &
               &                           -radiusSplineKernel           &
               &                          )**3
       else
          particulateSmoothingIntegrandR=0.0d0
       end if
    end select
    return
  end function particulateSmoothingIntegrandR

  double precision function particulateMassIntegrand(radius)
    !!{
    The integrand used to find the enclosed mass in the smoothed density profile defined by \cite{barnes_gravitational_2012} to
    account for gravitational softening.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radius

    if (radius > 0.0d0) then
       particulateMassIntegrand=+4.0d0                                                                               &
            &                   *Pi                                                                                  &
            &                   *radius**2                                                                           &
            &                   *energyDistribution%interpolate(radius,table=energyDistributionTableDensitySmoothed)
    else
       particulateMassIntegrand=+0.0d0
    end if
    return
  end function particulateMassIntegrand

  double precision function particulatePotentialIntegrand(radius)
    !!{
    The integrand used to find the gravitational potential in the smoothed density profile defined by
    \cite{barnes_gravitational_2012} to account for gravitational softening.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    double precision, intent(in   ) :: radius

    ! Evaluate the integrand for gravitational potential. No minus sign here as we actually want the relative potential which will
    ! be positive.
    particulatePotentialIntegrand=+gravitationalConstant_internal                                           &
         &                        *energyDistribution%interpolate(radius,table=energyDistributionTableMass) &
         &                        /radius**2
    return
  end function particulatePotentialIntegrand
