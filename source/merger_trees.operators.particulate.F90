!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which implements a merger tree operator which creates particle representations of \glc\ halos.

  use Cosmology_Parameters
  use Cosmology_Functions
  use ISO_Varying_String
  use Kind_Numbers
  use Tables

  !# <mergerTreeOperator name="mergerTreeOperatorParticulate" defaultThreadPrivate="yes">
  !#  <description>Provides a merger tree operator which create particle representations of \glc\ halos.</description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorParticulate
     !% A merger tree operator which create particle representations of \glc\ halos.
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_
     type            (varying_string          )          :: outputFileName
     double precision                                    :: massParticle        , radiusTruncateOverRadiusVirial, &
          &                                                 timeSnapshot        , lengthSoftening
     logical                                             :: satelliteOffset     , nonCosmological               , &
          &                                                 positionOffset      , addHubbleFlow                 , &
          &                                                 haloIdToParticleType, sampleParticleNumber
     integer                                             :: selection           , kernelSoftening               , &
          &                                                 chunkSize
     integer         (kind_int8               )          :: idMultiplier
   contains
     final     ::            particulateDestructor
     procedure :: operate => particulateOperate
  end type mergerTreeOperatorParticulate

  interface mergerTreeOperatorParticulate
     !% Constructors for the particulate merger tree operator class.
     module procedure particulateConstructorParameters
     module procedure particulateConstructorInternal
  end interface mergerTreeOperatorParticulate
  
  !# <enumeration>
  !#  <name>selection</name>
  !#  <description>Options for selection of nodes to particulate.</description>
  !#  <encodeFunction>yes</encodeFunction>
  !#  <validator>yes</validator>
  !#  <visibility>private</visibility>
  !#  <entry label="all"        />
  !#  <entry label="hosts"      />
  !#  <entry label="satellites" />
  !# </enumeration>
  
  !# <enumeration>
  !#  <name>particulateKernel</name>
  !#  <description>Options for softening kernel in particulate.</description>
  !#  <encodeFunction>yes</encodeFunction>
  !#  <validator>yes</validator>
  !#  <visibility>private</visibility>
  !#  <entry label="delta"  />
  !#  <entry label="plummer"/>
  !#  <entry label="gadget" />
  !# </enumeration>
  
  ! Entries in the energy distribution table.
  integer                                    , parameter   :: energyDistributionTableDensity        =1
  integer                                    , parameter   :: energyDistributionTablePotential      =2
  integer                                    , parameter   :: energyDistributionTableDistribution   =3
  integer                                    , parameter   :: energyDistributionTableDensitySmoothed=4
  integer                                    , parameter   :: energyDistributionTableMass           =5

  ! Tables used in construction of distribution functions.
  type            (treeNode                 ), pointer     :: particulateNode
  logical                                                  :: particularEnergyDistributionInitialized
  class           (table1D                  ), allocatable :: particulateRadiusDistribution
  type            (table1DLogarithmicCSpline)              :: particulateEnergyDistribution
  integer                                                  :: particulateSofteningKernel
  double precision                                         :: particulateEnergy                , particulateRadiusTruncate, &
       &                                                      particulateHeight                , particulateRadius        , &
       &                                                      particulateLengthSoftening
  !$omp threadprivate(particulateRadiusDistribution,particulateEnergyDistribution,particularEnergyDistributionInitialized,particulateEnergy,particulateHeight,particulateRadius)

contains

  function particulateConstructorParameters(parameters)
    !% Constructor for the particulate merger tree operator class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type(mergerTreeOperatorParticulate)                :: particulateConstructorParameters
    type(inputParameters              ), intent(inout) :: parameters
    type(varying_string               )                :: selection                       , kernelSoftening
    
    !# <inputParameter>
    !#   <name>outputFileName</name>
    !#   <source>parameters</source>
    !#   <variable>particulateConstructorParameters%outputFileName</variable>
    !#   <description>Name of the file to which particle data should be written.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>idMultiplier</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0_kind_int8</defaultValue>
    !#   <variable>particulateConstructorParameters%idMultiplier</variable>
    !#   <description>If this parameter is greater than zero, particle IDs begin at {\normalfont \ttfamily nodeIndex\*[idMultiplier]} for each node. The multiplier should be chosen to be large enough that duplicate IDs can not occur.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massParticle</name>
    !#   <source>parameters</source>
    !#   <variable>particulateConstructorParameters%massParticle</variable>
    !#   <description>Mass of particles to be used to represent halos.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>timeSnapshot</name>
    !#   <source>parameters</source>
    !#   <variable>particulateConstructorParameters%timeSnapshot</variable>
    !#   <description>The time at which to snapshot the tree.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>radiusTruncateOverRadiusVirial</name>
    !#   <source>parameters</source>
    !#   <variable>particulateConstructorParameters%radiusTruncateOverRadiusVirial</variable>
    !#   <description>Radius (in units of the virial radius) at which to truncate halo profiles.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>satelliteOffset</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.true.</defaultValue>
    !#   <variable>particulateConstructorParameters%satelliteOffset</variable>
    !#   <description>If true, offset particle representations to the positions/velocities of satellites.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>positionOffset</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.false.</defaultValue>
    !#   <variable>particulateConstructorParameters%positionOffset</variable>
    !#   <description>If true, offset particle representations to the positions/velocities of nodes.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>selection</name>
    !#   <source>parameters</source>
    !#   <defaultValue>var_str('all')</defaultValue>
    !#   <variable>selection</variable>
    !#   <description>Selects the type of halo to ouput. Allowed options are ``{\normalfont \ttfamily all}'', ``{\normalfont \ttfamily hosts}'', and ``{\normalfont \ttfamily satellites}''.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    particulateConstructorParameters%selection=enumerationSelectionEncode(char(selection),includesPrefix=.false.)
    !# <inputParameter>
    !#   <name>kernelSoftening</name>
    !#   <source>parameters</source>
    !#   <defaultValue>var_str('plummer')</defaultValue>
    !#   <variable>kernelSoftening</variable>
    !#   <description>Selects the softening kernel to use. Allowed options are ``{\normalfont \ttfamily plummer}'', and ``{\normalfont \ttfamily gadget}''.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    particulateConstructorParameters%kernelSoftening=enumerationParticulateKernelEncode(char(kernelSoftening),includesPrefix=.false.)
    !# <inputParameter>
    !#   <name>nonCosmological</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.false.</defaultValue>
    !#   <variable>particulateConstructorParameters%nonCosmological</variable>
    !#   <description>If true, a non-cosmological snapshot file will be created.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>addHubbleFlow</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.false.</defaultValue>
    !#   <variable>particulateConstructorParameters%addHubbleFlow</variable>
    !#   <description>If true, Hubble flow will be added to velocity offsets of halos (if applied).</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>haloIdToParticleType</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.false.</defaultValue>
    !#   <variable>particulateConstructorParameters%haloIdToParticleType</variable>
    !#   <description>If true, the halo ID will be used to assign particles from each halo to a different Gadget particle type group.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>sampleParticleNumber</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.false.</defaultValue>
    !#   <variable>particulateConstructorParameters%sampleParticleNumber</variable>
    !#   <description>If true, the number of particles in each halo will be sampled from a Poisson distribution with the expected mean. Otherwise, the number is set equal to that expectation.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>lengthSoftening</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <variable>particulateConstructorParameters%lengthSoftening</variable>
    !#   <description>The Plummer-equivalent softening length. That is, the parameter $\epsilon$ in the softening gravitational potential $\phi(r) = -{\rm G}m/\sqrt{r^2+\epsilon^2}$. If set to zero, softening is ignored when constructing the particle representation of the halo. For non-zero values softening is accounted for when constructing the velocity distribution following the procedure of \cite{barnes_gravitational_2012}.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>chunkSize</name>
    !#   <source>parameters</source>
    !#   <defaultValue>-1</defaultValue>
    !#   <variable>particulateConstructorParameters%chunkSize</variable>
    !#   <description>HDF5 dataset chunk size.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="particulateConstructorParameters%cosmologyParameters_" source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"  name="particulateConstructorParameters%cosmologyFunctions_"  source="parameters"/>
    !# <inputParametersValidate source="parameters"/>
    return
  end function particulateConstructorParameters

  function particulateConstructorInternal(outputFileName,idMultiplier,massParticle,radiusTruncateOverRadiusVirial,timeSnapshot,satelliteOffset,positionOffset,selection,nonCosmological,addHubbleFlow,haloIdToParticleType,sampleParticleNumber,lengthSoftening,chunkSize,cosmologyParameters_,cosmologyFunctions_)
    !% Internal constructor for the particulate merger tree operator class.
    use Galacticus_Error
    implicit none
    type            (mergerTreeOperatorParticulate)                        :: particulateConstructorInternal
    type            (varying_string               ), intent(in   )         :: outputFileName
    integer         (kind_int8                    ), intent(in   )         :: idMultiplier
    double precision                               , intent(in   )         :: massParticle                  , radiusTruncateOverRadiusVirial, &
         &                                                                    timeSnapshot                  , lengthSoftening
    logical                                        , intent(in   )         :: satelliteOffset               , nonCosmological               , &
         &                                                                    positionOffset                , addHubbleFlow                 , &
         &                                                                    haloIdToParticleType          , sampleParticleNumber
    integer                                        , intent(in   )         :: selection                     , chunkSize
    class           (cosmologyParametersClass     ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass      ), intent(in   ), target :: cosmologyFunctions_

    ! Validate input.
    if (.not.enumerationSelectionIsValid(selection)) call Galacticus_Error_Report('invalid selection type'//{introspection:location})
    ! Store properties.
    particulateConstructorInternal%outputFileName                 =  outputFileName
    particulateConstructorInternal%idMultiplier                   =  idMultiplier
    particulateConstructorInternal%massParticle                   =  massParticle
    particulateConstructorInternal%radiusTruncateOverRadiusVirial =  radiusTruncateOverRadiusVirial
    particulateConstructorInternal%timeSnapshot                   =  timeSnapshot
    particulateConstructorInternal%satelliteOffset                =  satelliteOffset
    particulateConstructorInternal%positionOffset                 =  positionOffset
    particulateConstructorInternal%selection                      =  selection
    particulateConstructorInternal%nonCosmological                =  nonCosmological
    particulateConstructorInternal%addHubbleFlow                  =  addHubbleFlow
    particulateConstructorInternal%lengthSoftening                =  lengthSoftening
    particulateConstructorInternal%haloIdToParticleType           =  haloIdToParticleType
    particulateConstructorInternal%sampleParticleNumber           =  sampleParticleNumber
    particulateConstructorInternal%chunkSize                      =  chunkSize
    particulateConstructorInternal%cosmologyParameters_           => cosmologyParameters_  
    particulateConstructorInternal%cosmologyFunctions_            => cosmologyFunctions_
   return
  end function particulateConstructorInternal

  subroutine particulateDestructor(self)
    !% Destructor for the particulate merger tree operator class.
    implicit none
    type(mergerTreeOperatorParticulate), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    !# <objectDestructor name="self%cosmologyFunctions_" />
    return
 end subroutine particulateDestructor

  subroutine particulateOperate(self,tree)
    !% Perform a particulation operation on a merger tree (i.e. create a particle representation of the tree).
    use ISO_Varying_String
    use IO_HDF5
    use Numerical_Constants_Math
    use Pseudo_Random
    use Dark_Matter_Halo_Scales
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Galacticus_Calculations_Resets
    use Coordinates
    use Numerical_Comparison
    use Memory_Management
    use Galacticus_Error
    use Galacticus_Display
    implicit none
    class           (mergerTreeOperatorParticulate), intent(inout)               :: self
    type            (mergerTree                   ), intent(inout) , target      :: tree
    type            (treeNode                     ), pointer                     :: node
    class           (nodeComponentBasic           ), pointer                     :: basic
    class           (nodeComponentPosition        ), pointer                     :: position
    class           (nodeComponentSatellite       ), pointer                     :: satellite
    type            (mergerTree                   ), pointer                     :: treeCurrent
    class           (darkMatterHaloScaleClass     ), pointer                     :: darkMatterHaloScale_
    double precision                               , parameter                   :: tolerance                 =1.0d-06
    double precision                               , parameter                   :: unitGadgetMass            =1.0d+10
    double precision                               , parameter                   :: unitGadgetLength          =1.0d-03
    double precision                               , parameter                   :: unitGadgetVelocity        =1.0d+00
    double precision                               , parameter                   :: distributionFunctionBuffer=0.3d0
    double precision                               , dimension(3  )              :: positionVector                       , velocityVector             , &
         &                                                                          randomDeviates
    integer                                        , dimension(6  )              :: particleCounts
    double precision                               , dimension(:,:), allocatable :: particlePosition                     , particleVelocity
    integer         (kind_int8                    ), dimension(  :), allocatable :: particleIDs
    double precision                                                             :: particleCountMean                    , distributionFunction       , &
         &                                                                          radiusVirial                         , radiusTruncate             , &
         &                                                                          speedPrevious                        , massTruncate               , &
         &                                                                          energy                               , radiusEnergy               , &
         &                                                                          energyPotential                      , speed                      , &
         &                                                                          speedEscape                          , distributionFunctionMaximum
    integer                                                                      :: particleCountActual                  , i                          , &
         &                                                                          j                                    , typeIndex
    logical                                                                      :: isNew                                , keepSample
    type            (coordinateCartesian          )                              :: positionCartesian                    , velocityCartesian
    type            (coordinateSpherical          )                              :: positionSpherical                    , velocitySpherical
    type            (hdf5Object                   )                              :: outputFile                           , header                     , &
         &                                                                          particleGroup
    type            (varying_string               )                              :: message
    character       (len=13                       )                              :: label
    character       (len= 9                       )                              :: groupName
    
    ! Get required objects.
    darkMatterHaloScale_ => darkMatterHaloScale()
    ! Open the HDF5 file for output.
    !$omp critical (HDF5_Access)
    call outputFile%openFile(char(self%outputFileName),overWrite=.false.,readOnly=.false.)
    ! Check if the Gadget header has already been created.
    if (.not.outputFile%hasGroup('Header')) then
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
       call header%close()
    end if
    call outputFile%close()
    !$omp end critical (HDF5_Access)
    ! Iterate over trees.
    treeCurrent             => tree
    do while (associated(treeCurrent))
       ! Get root node of the tree.
       node  => treeCurrent%baseNode
       ! Walk the tree.
       do while (associated(node))
          call Galacticus_Calculations_Reset(node)
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
             particularEnergyDistributionInitialized=.false.
             ! Determine the virial radius.
             radiusVirial=darkMatterHaloScale_%virialRadius(node)
             ! Determine the truncation radius.
             radiusTruncate=+self%radiusTruncateOverRadiusVirial &
                  &         *radiusVirial
             ! Determine the mass within the truncation radius.
             massTruncate=Galactic_Structure_Enclosed_Mass(                             &
                  &                                        node                       , &
                  &                                        radiusTruncate             , &
                  &                                        massType      =massTypeDark, &
                  &                                        haloLoaded    =.true.        &
                  &                                       ) 
             ! Determine the mean number of particles required to represent this node.          
             particleCountMean   =+massTruncate      &
                  &               /self%massParticle
             ! Determine the actual number of particles to use to represent the node.
             if (self%sampleParticleNumber) then
                particleCountActual=tree%randomNumberGenerator%poissonSample(particleCountMean)
             else
                particleCountActual=nint(particleCountMean)
             end if
             ! Allocate space for particle data.
             call allocateArray(particlePosition,[3,particleCountActual])
             call allocateArray(particleVelocity,[3,particleCountActual])
             call allocateArray(particleIDs     ,[  particleCountActual])
             ! Get required components.
             position  => node%position ()
             satellite => node%satellite()
             ! Set pointers to module-scope variables.
             particulateNode            => node
             particulateRadiusTruncate  =  radiusTruncate
             particulateLengthSoftening =  self%lengthSoftening
             particulateSofteningKernel =  self%kernelSoftening
             ! Iterate over particles.
             isNew=.true.
             !$omp parallel do private(i,j,positionSpherical,positionCartesian,velocitySpherical,velocityCartesian,energy,energyPotential,speed,speedEscape,speedPrevious,distributionFunction,distributionFunctionMaximum,keepSample,radiusEnergy,positionVector,velocityVector,randomDeviates)
             do i=1,particleCountActual
                if (OMP_Get_Thread_Num() == 0) then
                   call Galacticus_Display_Counter(max(1,int(100.0d0*dble(i-1)/dble(particleCountActual))),isNew=isNew,verbosity=verbosityStandard)
                   isNew=.false.
                end if
                ! Sample particle positions from the halo density distribution. Currently, we assume that halos are spherically
                ! symmetric.
                !$omp critical (mergerTreeOperatorParticulateSample)
                do j=1,3
                   randomDeviates(j)=tree%randomNumberGenerator%uniformSample()
                end do
                !$omp end critical (mergerTreeOperatorParticulateSample)
                call positionSpherical%  phiSet(     2.0d0*Pi*randomDeviates(1)       )
                call positionSpherical%thetaSet(acos(2.0d0   *randomDeviates(2)-1.0d0))                
                call positionSpherical%    rSet(                                                                        &
                     &                          Galactic_Structure_Radius_Enclosing_Mass(                               &
                     &                                                                   node                         , &
                     &                                                                   mass      =+massTruncate       &
                     &                                                                              *randomDeviates(3), &
                     &                                                                   massType  =massTypeDark      , &
                     &                                                                   haloLoaded=.true.              &
                     &                                                                  )                               &
                     &                         )
                ! Get the corresponding cartesian coordinates.
                positionCartesian=positionSpherical
                ! Construct the energy distribution function encompassing this radius.
                call particulateTabulateEnergyDistribution(positionSpherical%r())
                ! Find the potential energy and escape speed (actually the speed to reach the truncation radius) at this radius.
                energyPotential=+particulateEnergyDistribution%interpolate(                                        &
                     &                                                           positionSpherical%r()           , &
                     &                                                     table=energyDistributionTablePotential  &
                     &                                                     )
                speedEscape    =+sqrt(                 &
                     &                +2.0d0           &
                     &                *energyPotential &
                     &               )
                ! Estimate the maximum of the speed distribution function.
                distributionFunctionMaximum=+0.0d0
                speedPrevious              =-1.0d0
                do j=1,particulateEnergyDistribution%size()
                   ! Find the energy at this radius in the tabulated profile.
                   energy=particulateEnergyDistribution%y(j,table=energyDistributionTablePotential)
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
                      distributionFunction       =+speed**2                                                                                                                        &
                           &                      *particulateEnergyDistribution%interpolateGradient(particulateEnergyDistribution%x(j),table=energyDistributionTableDistribution) &
                           &                      /particulateEnergyDistribution%interpolateGradient(particulateEnergyDistribution%x(j),table=energyDistributionTablePotential   )
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
                   speed               =+tree%randomNumberGenerator%uniformSample() &
                        &               *speedEscape
                   !$omp end critical (mergerTreeOperatorParticulateSample)
                   energy              =+energyPotential    &
                        &               -0.5d0              &
                        &               *speed          **2
                   !$omp critical (mergerTreeOperatorParticulateRadius)
                   radiusEnergy        =+particulateRadiusDistribution%interpolate        (                                           &
                        &                                                                        energy                             , &
                        &                                                                  table=energyDistributionTablePotential     &
                        &                                                                 )
                   !$omp end critical (mergerTreeOperatorParticulateRadius)
                   distributionFunction=+speed**2                                                                                     &
                        &               *particulateEnergyDistribution%interpolateGradient(                                           &
                        &                                                                        radiusEnergy                       , &
                        &                                                                  table=energyDistributionTableDistribution  &
                        &                                                                 )                                           &
                        &               /particulateEnergyDistribution%interpolateGradient(                                           &
                        &                                                                        radiusEnergy                       , &
                        &                                                                  table=energyDistributionTablePotential     &
                        &                                                                 )
                   if (distributionFunction > distributionFunctionMaximum) then
                      write (label,'(e12.6)') distributionFunction
                      message='distribution function ['//trim(label)//'] exceeds estimated maximum ['
                      write (label,'(e12.6)') distributionFunctionMaximum
                      message=message//trim(label)//']'
                      call Galacticus_Error_Report(message//{introspection:location})
                   end if
                   !$omp critical (mergerTreeOperatorParticulateSample)
                   keepSample= +tree%randomNumberGenerator%uniformSample() &
                        &     <                                     &
                        &      +distributionFunction                &
                        &      /distributionFunctionMaximum
                   !$omp end critical (mergerTreeOperatorParticulateSample)
                end do
                ! Choose a velocity vector in spherical coordinates with velocity chosen to give the required kinetic energy.
                !$omp critical (mergerTreeOperatorParticulateSample)
                call velocitySpherical%  phiSet(     2.0d0*Pi*tree%randomNumberGenerator%uniformSample()       )
                call velocitySpherical%thetaSet(acos(2.0d0   *tree%randomNumberGenerator%uniformSample()-1.0d0))
                !$omp end critical (mergerTreeOperatorParticulateSample)
                call velocitySpherical%    rSet(speed                                                                           )
                ! Get the corresponding cartesian coordinates.
                velocityCartesian=velocitySpherical
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
             !$omp end parallel do
             call Galacticus_Display_Counter_Clear(verbosity=verbosityWorking)
             ! Perform unit conversion.
             particlePosition=particlePosition/unitGadgetLength
             particleVelocity=particleVelocity/unitGadgetVelocity
             ! Accumulate the particle data to file.
             !$omp critical (HDF5_Access)
             call outputFile%openFile(char(self%outputFileName),overWrite=.false.,readOnly=.false.,objectsOverwritable=.true.)
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
             particleGroup=outputFile%openGroup(groupName,'Group containing particle data for halos',chunkSize=self%chunkSize)
             call particleGroup%writeDataset(particlePosition,'Coordinates','Particle coordinates',appendTo=self%chunkSize /= -1,appendDimension=2)
             call particleGroup%writeDataset(particleVelocity,'Velocities' ,'Particle velocities' ,appendTo=self%chunkSize /= -1,appendDimension=2)
             call particleGroup%writeDataset(particleIDs     ,'ParticleIDs','Particle IDs'        ,appendTo=self%chunkSize /= -1                  )
             call particleGroup%close()
             call deallocateArray(particlePosition)
             call deallocateArray(particleVelocity)
             call deallocateArray(particleIDs     )
             ! Update particle counts.
             particleCounts(typeIndex)=particleCounts(typeIndex)+particleCountActual
             call header%writeAttribute(particleCounts,'NumPart_ThisFile')
             call header%writeAttribute(particleCounts,'NumPart_Total'   )
             call header%close()
             call outputFile%close()
             !$omp end critical (HDF5_Access)
          end if
          ! Walk to the next node.
          node => node%walkTreeWithSatellites()
       end do
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end subroutine particulateOperate
  
  subroutine particulateTabulateEnergyDistribution(radius)
    !% Construct the energy distribution function assuming a spherical dark matter halo with
    !% isotropic velocity dispersion. We solve Eddington's formula
    !% \citep[][eqn. 4.43a]{binney_galactic_2008}.    
    !% \begin{equation}
    !%  f(E) = {1 \over \sqrt{8} \pi^2} {{\rm d}\over {\rm d}E} \int_E^0 {{\rm d}\Phi \over \sqrt{\Phi-E}} {{\rm d}\rho \over {\rm d} \Phi}.
    !% \end{equation}
    !% In practice, we tabulate:
    !% \begin{equation}
    !%  F(E) = \int_E^0 {{\rm d}\Phi \over \sqrt{\Phi-E}} {{\rm d}\rho \over {\rm d} \Phi},
    !% \end{equation}
    !% which we can then take the derivative of numerically to obtain the distribution function.
    use FGSL
    use Numerical_Integration
    use Dark_Matter_Profiles
    use Dark_Matter_Halo_Scales
    use Table_Labels
    use Galacticus_Error
    implicit none
    double precision                            , intent(in   ) :: radius
    class           (darkMatterProfileClass    ), pointer       :: darkMatterProfile_
    class           (darkMatterHaloScaleClass  ), pointer       :: darkMatterHaloScale_
    double precision                            , parameter     :: radiusVirialFraction                     =1.0d-5
    double precision                            , parameter     :: toleranceTabulation                      =1.0d-6
    double precision                            , parameter     :: toleranceGradient                        =1.0d-6
    double precision                            , parameter     :: energyDistributionPointsPerDecade        =3.0d+1
    double precision                                            :: radiusMinimum                                   , energyPotentialTruncate                  , &
         &                                                         particulateSmoothingIntegrationRangeLower       , particulateSmoothingIntegrationRangeUpper, &
         &                                                         radiusFactorAsymptote                           , integralAsymptotic                       , &
         &                                                         gradientDensityPotentialLower                   , gradientDensityPotentialUpper
    logical                                                     :: tableRebuild
    integer                                                     :: i
    integer                                                     :: radiusCount
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace

    ! Get required objects.
    darkMatterProfile_   => darkMatterProfile  ()
    darkMatterHaloScale_ => darkMatterHaloScale()
    ! Determine the minimum of the given radius and some small fraction of the virial radius.
    radiusMinimum=min(radius/2.0d0,radiusVirialFraction*darkMatterHaloScale_%virialRadius(particulateNode))
    ! Rebuild the density vs. potential table to have sufficient range if necessary.
    if (particularEnergyDistributionInitialized) then
       tableRebuild=(radiusMinimum < (1.0d0-toleranceTabulation)*particulateEnergyDistribution%x(1))
       if (tableRebuild) call particulateEnergyDistribution%destroy()
    else
       tableRebuild                           =.true.
       particularEnergyDistributionInitialized=.true.
    end if
    if (tableRebuild) then
       ! Build tables of potential and density.
       radiusCount=int(energyDistributionPointsPerDecade*log10(particulateRadiusTruncate/radiusMinimum))+1
       call particulateEnergyDistribution%create(                                                          &
            &                                                +radiusMinimum                              , &
            &                                                +particulateRadiusTruncate                  , &
            &                                                 radiusCount                                , &
            &                              tableCount       = 5                                          , &
            &                              extrapolationType= [extrapolationTypeFix,extrapolationTypeFix]  &
            &                             )
       select case (particulateSofteningKernel)
       case (particulateKernelDelta)
          energyPotentialTruncate=+darkMatterProfile_%potential(                            &
               &                                                 particulateNode          , &
               &                                                +particulateRadiusTruncate  &
               &                                               )
       case default
          ! Potential will be computed directly from the smoothed density profile in these cases.
          energyPotentialTruncate=0.0d0
       end select
       do i=1,radiusCount
          particulateRadius=particulateEnergyDistribution%x(i)
          call particulateEnergyDistribution%populate(                                                                          & 
               &                                              +darkMatterProfile_%density  (                                    &
               &                                                                            particulateNode                   , &
               &                                                                            particulateEnergyDistribution%x(i)  &
               &                                                                           )                                  , &
               &                                                                                                            i , &
               &                                table        =energyDistributionTableDensity                                  , &
               &                                computeSpline=i==radiusCount                                                    &
               &                               )            
          select case (particulateSofteningKernel)
          case (particulateKernelDelta)
             ! No softening is applied, so use the actual density and potential.
             call particulateEnergyDistribution%populate(                                                                                 &
                  &                                              particulateEnergyDistribution%y(i,table=energyDistributionTableDensity), &
                  &                                                                                                                   i , &
                  &                                table        =energyDistributionTableDensitySmoothed                                 , &
                  &                                computeSpline=i==radiusCount                                                           &
                  &                               )            
             call particulateEnergyDistribution%populate(                                                                          & 
                  &                                              -darkMatterProfile_%potential(                                    &
                  &                                                                            particulateNode                   , &
                  &                                                                            particulateEnergyDistribution%x(i)  &
                  &                                                                           )                                    &
                  &                                              +energyPotentialTruncate                                        , &
                  &                                                                                                            i , &
                  &                                table        =energyDistributionTablePotential                                , &
                  &                                computeSpline=i==radiusCount                                                    &
                  &                               )
         case default
             ! Compute potential from a density field smoothed by the density distribution corresponding to the softened
             ! potential.
             particulateSmoothingIntegrationRangeUpper=+particulateRadiusTruncate-particulateRadius
             particulateSmoothingIntegrationRangeLower=-particulateRadiusTruncate-particulateRadius
             select case (particulateSofteningKernel)
             case (particulateKernelGadget )
                particulateSmoothingIntegrationRangeUpper=min(particulateSmoothingIntegrationRangeUpper,+2.8d0*particulateLengthSoftening)
                particulateSmoothingIntegrationRangeLower=max(particulateSmoothingIntegrationRangeLower,-2.8d0*particulateLengthSoftening)
             end select
             call particulateEnergyDistribution%populate(                                                                                 &
                  &                                              +Integrate(                                                              &
                  &                                                                           +particulateSmoothingIntegrationRangeLower, &
                  &                                                                           +particulateSmoothingIntegrationRangeUpper, &
                  &                                                                            particulateSmoothingIntegrandZ           , &
                  &                                                                            integrandFunction                        , &
                  &                                                                            integrationWorkspace                     , &
                  &                                                         toleranceAbsolute=+0.0d0                                    , &
                  &                                                         toleranceRelative=+1.0d-9                                     &
                  &                                                       )                                                             , &
                  &                                                                                                                   i , &
                  &                                table        =energyDistributionTableDensitySmoothed                                 , &
                  &                                computeSpline=i==radiusCount                                                           &
                  &                               )            
             call Integrate_Done(integrandFunction,integrationWorkspace)
          end select
       end do
       ! If necessary, compute the potential from the smoothed density profile.
       if (particulateSofteningKernel /= particulateKernelDelta) then
          do i=1,radiusCount
             particulateRadius=particulateEnergyDistribution%x(i)
             if     (                                                                                    &
                  &    i                                                                                 &
                  &   >                                                                                  &
                  &    1                                                                                 &           
                  &  .and.                                                                               &
                  &    particulateEnergyDistribution%y(i  ,table=energyDistributionTableDensitySmoothed) &
                  &   >                                                                                  &
                  &    particulateEnergyDistribution%y(i-1,table=energyDistributionTableDensitySmoothed) &
                  & ) call particulateEnergyDistribution%populate(particulateEnergyDistribution%y(i-1,table=energyDistributionTableDensitySmoothed),i,table=energyDistributionTableDensitySmoothed)
             call particulateEnergyDistribution%populate(                                                                      &
                  &                                              +Integrate(                                                   &
                  &                                                                           +0.0d0                         , &
                  &                                                                           +particulateRadius             , &
                  &                                                                            particulateMassIntegrand      , &
                  &                                                                            integrandFunction             , &
                  &                                                                            integrationWorkspace          , &
                  &                                                         toleranceAbsolute=+0.0d0                         , &
                  &                                                         toleranceRelative=+1.0d-8                          &
                  &                                                       )                                                  , &
                  &                                                                                                        i , &
                  &                                table        =energyDistributionTableMass                                 , &
                  &                                computeSpline=i==radiusCount                                                &
                  &                               )            
             call Integrate_Done(integrandFunction,integrationWorkspace)
          end do
          do i=1,radiusCount
             particulateRadius=particulateEnergyDistribution%x(i)
             call particulateEnergyDistribution%populate(                                                                      &
                  &                                              +Integrate(                                                   &
                  &                                                                           +particulateRadius             , &
                  &                                                                           +particulateRadiusTruncate     , &
                  &                                                                            particulatePotentialIntegrand , &
                  &                                                                            integrandFunction             , &
                  &                                                                            integrationWorkspace          , &
                  &                                                         toleranceAbsolute=+0.0d0                         , &
                  &                                                         toleranceRelative=+1.0d-9                          &
                  &                                                       )                                                  , &
                  &                                                                                                        i , &
                  &                                table        =energyDistributionTablePotential                            , &
                  &                                computeSpline=i==radiusCount                                                &
                  &                               )            
             call Integrate_Done(integrandFunction,integrationWorkspace)
          end do
       end if
       ! Evaluate the integral in Eddington's formula. Ignore the normalization as we will simply use rejection sampling to draw
       ! from this distribution.
       call particulateEnergyDistribution%populate(0.0d0,radiusCount,table=energyDistributionTableDistribution,computeSpline=.false.)
       do i=1,radiusCount-1
          particulateRadius=particulateEnergyDistribution%x(i)
          particulateEnergy=particulateEnergyDistribution%y(i,table=energyDistributionTablePotential)
          ! As the integrand diverges at particulateRadius we evaluate the integral analytically in a small region close to that radius. The
          ! size of this region is chosen such that the gradient of density with respect to potential changes very little over that range,
          ! ensuring that the approximation made in the analytic integral is valid.
          gradientDensityPotentialLower=+particulateEnergyDistribution%interpolateGradient(particulateRadius,table=energyDistributionTableDensity  ) &
               &                        /particulateEnergyDistribution%interpolateGradient(particulateRadius,table=energyDistributionTablePotential)
          radiusFactorAsymptote        =+particulateEnergyDistribution%x(i+1) &
               &                        /particulateEnergyDistribution%x(i  )
          ! Iteratively reduce the size of the region to be integrated analytically until the density vs. potential gradient is
          ! sufficiently constant across that range.
          do while (.true.)
             gradientDensityPotentialUpper=+particulateEnergyDistribution%interpolateGradient(particulateRadius*radiusFactorAsymptote,table=energyDistributionTableDensity  ) &
                  &                        /particulateEnergyDistribution%interpolateGradient(particulateRadius*radiusFactorAsymptote,table=energyDistributionTablePotential)
             if     (                                                                                           &
                  &                            abs(gradientDensityPotentialUpper-gradientDensityPotentialLower) &
                  &  <                                                                                          &
                  &   +0.5d0*toleranceGradient*abs(gradientDensityPotentialUpper+gradientDensityPotentialLower) &
                  & ) exit
             radiusFactorAsymptote=sqrt(radiusFactorAsymptote)
          end do
          ! Evaluate the asymptotic part of the integral analytically.
          integralAsymptotic=+2.0d0                                                                                                                           &
               &             *gradientDensityPotentialLower                                                                                                   &
               &             *sqrt(                                                                                                                           &
               &                   +particulateEnergyDistribution%y          (i                                      ,table=energyDistributionTablePotential) &
               &                   -particulateEnergyDistribution%interpolate(particulateRadius*radiusFactorAsymptote,table=energyDistributionTablePotential) &
               &                  )
          ! Evaluate the remainder of the integral numerically and add on the asympototic part.
          call particulateEnergyDistribution%populate(                                                                                            & 
               &                                              +Integrate(                                                                         &
               &                                                                           +log(particulateRadiusTruncate                      ), &
               &                                                                           +log(particulateRadius        *radiusFactorAsymptote), &
               &                                                                            particulateEddingtonIntegrand                       , &
               &                                                                            integrandFunction                                   , &
               &                                                                            integrationWorkspace                                , &
               &                                                         toleranceAbsolute=+0.0d0                                               , &
               &                                                         toleranceRelative=+1.0d-5                                                &
               &                                                        )                                                                         &
               &                                              +integralAsymptotic                                                               , &
               &                                                                                                                              i , &
               &                                table        =energyDistributionTableDistribution                                               , &
               &                                computeSpline=i==radiusCount-1                                                                    &
               &                               )
          call Integrate_Done(integrandFunction,integrationWorkspace)
          ! Check that the distribution function is monotonically increasing.
          if     (                                                                                 &
               &    i                                                                              &
               &   >                                                                               &
               &    1                                                                              &           
               &  .and.                                                                            &
               &    particulateEnergyDistribution%y(i  ,table=energyDistributionTableDistribution) &
               &   >                                                                               &
               &    particulateEnergyDistribution%y(i-1,table=energyDistributionTableDistribution) &
               & ) call Galacticus_Error_Report('unphysical distribution function'//{introspection:location})
       end do
       ! Construct a reversed (radius vs. potential function) table.
       call particulateEnergyDistribution%reverse(particulateRadiusDistribution,table=energyDistributionTablePotential)
       ! Record that the distribution is initialized.
       particularEnergyDistributionInitialized=.true.
    end if
    return
  end subroutine particulateTabulateEnergyDistribution

  double precision function particulateEddingtonIntegrand(lradius)
    !% The integrand appearing in Eddington's formula for the distribution function.
    implicit none
    double precision, intent(in   ) :: lradius
    double precision                :: potential
    double precision :: radius

    radius=exp(lradius)
    potential=particulateEnergyDistribution%interpolate(radius,table=energyDistributionTablePotential)
    if (potential < particulateEnergy) then
       particulateEddingtonIntegrand=radius*particulateEnergyDistribution%interpolateGradient(radius,table=energyDistributionTableDensity)/sqrt(particulateEnergy-potential)
    else
       particulateEddingtonIntegrand=0.0d0
    end if
    return
  end function particulateEddingtonIntegrand

  double precision function particulateSmoothingIntegrandZ(height)
    !% The integrand over cylindrical coordinate $z$ used in finding the smoothed density profile defined by
    !% \cite{barnes_gravitational_2012} to account for gravitational softening.
    use FGSL
    use Numerical_Integration
    implicit none
    double precision                            , intent(in   ) :: height
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    double precision                                            :: radiusMaximum       , heightOffset, &
         &                                                         argumentSqrt

    heightOffset=height-particulateRadius ! Height in the profile (not the kernel).
    argumentSqrt=+particulateRadiusTruncate**2-heightOffset**2
    if (argumentSqrt <= 0.0d0) then
       particulateSmoothingIntegrandZ=0.0d0
    else
       radiusMaximum    =sqrt(argumentSqrt)
       particulateHeight=     height
       select case (particulateSofteningKernel)
       case (particulateKernelGadget )
          radiusMaximum=min(radiusMaximum,sqrt(+(2.8d0*particulateLengthSoftening)**2-particulateHeight**2))
       end select
       particulateSmoothingIntegrandZ=Integrate(                                                   &
            &                                                     +0.0d0                         , &
            &                                                     +radiusMaximum                 , &
            &                                                      particulateSmoothingIntegrandR, &
            &                                                      integrandFunction             , &
            &                                                      integrationWorkspace          , &
            &                                   toleranceAbsolute=+0.0d0                         , &
            &                                   toleranceRelative=+1.0d-8                          &
            &                                  )
       call Integrate_Done(integrandFunction,integrationWorkspace)
    end if
    return
  end function particulateSmoothingIntegrandZ

  double precision function particulateSmoothingIntegrandR(radiusCylindrical)
    !% The integrand over cylindrical coordinate $z$ used in finding the smoothed density profile defined by
    !% \cite{barnes_gravitational_2012} to account for gravitational softening.
    use Dark_Matter_Profiles
    implicit none
    double precision                        , intent(in   ) :: radiusCylindrical
    class           (darkMatterProfileClass), pointer       :: darkMatterProfile_
    double precision                                        :: radiusSplineKernel, lengthSplineKernel

    darkMatterProfile_ => darkMatterProfile  ()
    particulateSmoothingIntegrandR=+radiusCylindrical                                        &
         &                         *darkMatterProfile_%density(                              &
         &                                                     particulateNode             , &
         &                                                     sqrt(                         &
         &                                                          +radiusCylindrical  **2  &
         &                                                          +(                       &
         &                                                            +particulateHeight     &
         &                                                            -particulateRadius     &
         &                                                           )                  **2  &
         &                                                         )                         &
         &                                                    )
    ! Apply the softening kernel density distribution.
    select case (particulateSofteningKernel)
    case (particulateKernelPlummer)
       particulateSmoothingIntegrandR=+particulateSmoothingIntegrandR &
            &                         *1.5d0                          &
            &                         *particulateLengthSoftening **2 &
            &                         /(                              &
            &                          +radiusCylindrical         **2 &
            &                          +particulateHeight         **2 &
            &                          +particulateLengthSoftening**2 &
            &                         )**2.5d0
    case (particulateKernelGadget )
       ! Compute the radius in units of the spline kernel length scale.
       lengthSplineKernel=+2.8d0                      &
            &             *particulateLengthSoftening
       radiusSplineKernel=+sqrt(                      &
            &                   +radiusCylindrical**2 &
            &                   +particulateHeight**2 &
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
    !% The integrand used to find the enclosed mass in the smoothed density profile defined by \cite{barnes_gravitational_2012} to
    !% account for gravitational softening.
    use Numerical_Constants_Math
    implicit none
    double precision, intent(in   ) :: radius

    if (radius > 0.0d0) then
       particulateMassIntegrand=+4.0d0                                                                                          &
            &                   *Pi                                                                                             &
            &                   *radius**2                                                                                      &
            &                   *particulateEnergyDistribution%interpolate(radius,table=energyDistributionTableDensitySmoothed)
    else
       particulateMassIntegrand=+0.0d0
    end if
    return
  end function particulateMassIntegrand

  double precision function particulatePotentialIntegrand(radius)
    !% The integrand used to find the gravitational potential in the smoothed density profile defined by
    !% \cite{barnes_gravitational_2012} to account for gravitational softening.
    use Numerical_Constants_Physical
    implicit none
    double precision, intent(in   ) :: radius

    ! Evaluate the integrand for gravitational potential. No minus sign here as we actually want the relative potential which will
    ! be positive.
    particulatePotentialIntegrand=+gravitationalConstantGalacticus                                                     &
         &                        *particulateEnergyDistribution%interpolate(radius,table=energyDistributionTableMass) &
         &                        /radius**2
    return
  end function particulatePotentialIntegrand
