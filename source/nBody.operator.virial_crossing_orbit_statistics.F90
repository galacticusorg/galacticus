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

!!{
Implements an N-body data operator which computes virial crossing orbit statistics in bins of separation.
!!}
  
  use            :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use            :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Input_Parameters        , only : inputParameters
  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  !![
  <nbodyOperator name="nbodyOperatorVirialCrossingOrbitStatistics">
   <description>An N-body data operator which computes virial crossing orbit statistics in bins of separation.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorVirialCrossingOrbitStatistics
     !!{
     An N-body data operator which computes virial crossing orbit statistics in bins of separation.
     !!}
     private
     double precision                                      :: velocityMinimum                 , velocityMaximum     , &
          &                                                   separationMinimum               , separationMaximum   , &
          &                                                   bootstrapSampleRate             , time                , &
          &                                                   velocityCut                     , redshift
     integer         (c_size_t                  )          :: velocityCount                   , bootstrapSampleCount
     logical                                               :: includeUnbootstrapped           , crossCount          , &
          &                                                   addHubbleFlow
     class           (randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()
     class           (darkMatterHaloScaleClass  ), pointer :: darkMatterHaloScale_   => null()
     class           (cosmologyFunctionsClass   ), pointer :: cosmologyFunctions_    => null()
     type            (inputParameters           )          :: parameters
   contains
     final     ::            virialCrossingOrbitStatisticsDestructor
     procedure :: operate => virialCrossingOrbitStatisticsOperate
  end type nbodyOperatorVirialCrossingOrbitStatistics

  interface nbodyOperatorVirialCrossingOrbitStatistics
     !!{
     Constructors for the \refClass{nbodyOperatorVirialCrossingOrbitStatistics} N-body operator class.
     !!}
     module procedure virialCrossingOrbitStatisticsConstructorParameters
     module procedure virialCrossingOrbitStatisticsConstructorInternal
  end interface nbodyOperatorVirialCrossingOrbitStatistics

contains

  function virialCrossingOrbitStatisticsConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorVirialCrossingOrbitStatistics} N-body operator class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions, only : cosmologyFunctionsClass
    use :: Input_Parameters   , only : inputParameter
    implicit none
    type            (nbodyOperatorVirialCrossingOrbitStatistics)                :: self
    type            (inputParameters                           ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass                ), pointer       :: randomNumberGenerator_
    class           (darkMatterHaloScaleClass                  ), pointer       :: darkMatterHaloScale_
    class           (cosmologyFunctionsClass                   ), pointer       :: cosmologyFunctions_
    type            (inputParameters                           ), pointer       :: parametersRoot
    double precision                                                            :: velocityMinimum       , velocityMaximum     , &
         &                                                                         separationMinimum     , separationMaximum   , &
         &                                                                         bootstrapSampleRate   , redshift            , &
         &                                                                         time                  , velocityCut
    integer         (c_size_t                                  )                :: velocityCount         , bootstrapSampleCount
    logical                                                                     :: includeUnbootstrapped , crossCount          , &
         &                                                                         addHubbleFlow

    !![
    <inputParameter>
      <name>crossCount</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, compute cross-simulation virial crossing orbit statistics between the first and all simulations. Otherwise, compute virial crossing orbit statistics within each simulation.</description>
    </inputParameter>
    <inputParameter>
      <name>addHubbleFlow</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, add Hubble flow to velocities.</description>
    </inputParameter>
    <inputParameter>
      <name>redshift</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The redshift.</description>
    </inputParameter>
    <inputParameter>
      <name>bootstrapSampleCount</name>
      <source>parameters</source>
      <defaultValue>30_c_size_t</defaultValue>
      <description>The number of bootstrap resamples of the particles that should be used.</description>
    </inputParameter>
    <inputParameter>
      <name>bootstrapSampleRate</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>The sampling rate for particles.</description>
    </inputParameter>
    <inputParameter>
      <name>velocityMinimum</name>
      <source>parameters</source>
      <description>The minimum velocity to consider for virial crossing orbit statistics.</description>
    </inputParameter>
    <inputParameter>
      <name>velocityMaximum</name>
      <source>parameters</source>
      <description>The maximum velocity to consider for virial crossing orbit statistics.</description>
    </inputParameter>
    <inputParameter>
      <name>velocityCount</name>
      <source>parameters</source>
      <description>The number of bins in separation for virial crossing orbit statistics.</description>
    </inputParameter>
    <inputParameter>
      <name>separationMinimum</name>
      <source>parameters</source>
      <description>The minimum separation (in units of virial radii) to consider for virial crossing orbit statistics.</description>
    </inputParameter>
    <inputParameter>
      <name>separationMaximum</name>
      <source>parameters</source>
      <description>The maximum separation (in units of virial radii) to consider for virial crossing orbit statistics.</description>
    </inputParameter>
    <inputParameter>
      <name>includeUnbootstrapped</name>
      <source>parameters</source>
      <description>If true, include results for the unbootstrapped (i.e. original) sample.</description>
      <defaultValue>.true.</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>velocityCut</name>
      <source>parameters</source>
      <description>Pairs with (radial or tangential) velocities in excess of this value (in virial units) are cut from the sample.</description>
      <defaultValue>huge(0.0d0)</defaultValue>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"   name="darkMatterHaloScale_"   source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    !!]
    time=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift))
    if (associated(parameters%parent)) then
       parametersRoot => parameters%parent
       do while (associated(parametersRoot%parent))
          parametersRoot => parametersRoot%parent
       end do
       self=nbodyOperatorVirialCrossingOrbitStatistics(velocityMinimum,velocityMaximum,velocityCount,separationMinimum,separationMaximum,time,crossCount,addHubbleFlow,velocityCut,includeUnbootstrapped,bootstrapSampleCount,bootstrapSampleRate,randomNumberGenerator_,darkMatterHaloScale_,cosmologyFunctions_,parametersRoot)
    else
       self=nbodyOperatorVirialCrossingOrbitStatistics(velocityMinimum,velocityMaximum,velocityCount,separationMinimum,separationMaximum,time,crossCount,addHubbleFlow,velocityCut,includeUnbootstrapped,bootstrapSampleCount,bootstrapSampleRate,randomNumberGenerator_,darkMatterHaloScale_,cosmologyFunctions_,parameters    )
    end if
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    <objectDestructor name="darkMatterHaloScale_"  />
    <objectDestructor name="cosmologyFunctions_"   />
    !!]
    return
  end function virialCrossingOrbitStatisticsConstructorParameters

  function virialCrossingOrbitStatisticsConstructorInternal(velocityMinimum,velocityMaximum,velocityCount,separationMinimum,separationMaximum,time,crossCount,addHubbleFlow,velocityCut,includeUnbootstrapped,bootstrapSampleCount,bootstrapSampleRate,randomNumberGenerator_,darkMatterHaloScale_,cosmologyFunctions_,parameters) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorVirialCrossingOrbitStatistics} N-body operator class.
    !!}
    implicit none
    type            (nbodyOperatorVirialCrossingOrbitStatistics)                        :: self
    double precision                                            , intent(in   )         :: velocityMinimum       , velocityMaximum     , &
         &                                                                                 separationMinimum     , separationMaximum   , &
         &                                                                                 bootstrapSampleRate   , time                , &
         &                                                                                 velocityCut
    integer         (c_size_t                                  ), intent(in   )         :: velocityCount         , bootstrapSampleCount
    logical                                                     , intent(in   )         :: includeUnbootstrapped , crossCount          , &
         &                                                                                 addHubbleFlow
    class           (randomNumberGeneratorClass                ), intent(in   ), target :: randomNumberGenerator_
    class           (darkMatterHaloScaleClass                  ), intent(in   ), target :: darkMatterHaloScale_
    class           (cosmologyFunctionsClass                   ), intent(in   ), target :: cosmologyFunctions_
    type            (inputParameters                           ), intent(in   ), target :: parameters
    !![
    <constructorAssign variables="velocityMinimum, velocityMaximum, velocityCount, separationMinimum, separationMaximum, time, crossCount, addHubbleFlow, velocityCut, includeUnbootstrapped, bootstrapSampleCount, bootstrapSampleRate, *randomNumberGenerator_, *darkMatterHaloScale_, *cosmologyFunctions_"/>
    !!]

    self%parameters=inputParameters(parameters)
    self%redshift  =self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
    return
  end function virialCrossingOrbitStatisticsConstructorInternal

  subroutine virialCrossingOrbitStatisticsDestructor(self)
    !!{
    Destructor for the \refClass{nbodyOperatorVirialCrossingOrbitStatistics} N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorVirialCrossingOrbitStatistics), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    <objectDestructor name="self%darkMatterHaloScale_"  />
    !!]
    return
  end subroutine virialCrossingOrbitStatisticsDestructor

  subroutine virialCrossingOrbitStatisticsOperate(self,simulations)
    !!{
    Compute virial crossing orbit statistics in bins of separation.
    !!}
    use    :: Arrays_Search      , only : searchArray
    use    :: Display            , only : displayCounter                   , displayCounterClear                , displayIndent, displayMessage, &
         &                                displayUnindent                  , verbosityLevelStandard
    use    :: Error              , only : Error_Report
    use    :: Calculations_Resets, only : Calculations_Reset
    use    :: Galacticus_Nodes   , only : nodeComponentBasic               , treeNode
    use    :: HDF5_Access        , only : hdf5Access
    use    :: ISO_Varying_String , only : var_str
#ifdef USEMPI
    use    :: MPI_Utilities      , only : mpiSelf
#endif
    use    :: Nearest_Neighbors  , only : nearestNeighbors
    use    :: Node_Components    , only : Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize
    use    :: Numerical_Ranges   , only : Make_Range                       , rangeTypeLinear
    !$ use :: OMP_Lib            , only : OMP_Get_Thread_Num
    implicit none
    class           (nbodyOperatorVirialCrossingOrbitStatistics), intent(inout)                 :: self
    type            (nBodyData                                 ), intent(inout), dimension(:  ) :: simulations
    double precision                                            , parameter                     :: toleranceZero                =0.0d0
    integer                                                     , allocatable  , dimension(:  ) :: neighborIndex
    double precision                                            , allocatable  , dimension(:  ) :: neighborDistanceSquared            , velocityCentralBin               , &
         &                                                                                         velocityMinimumBin                 , velocityMaximumBin               , &
         &                                                                                         velocitiesRadial                   , velocitiesTangential             , &
         &                                                                                         separationsVirial                  , angularMomentumVirial            , &
         &                                                                                         energyVirial                       , velocityRadialMean               , &
         &                                                                                         velocityTangentialMean             , velocityRadialDispersion         , &
         &                                                                                         velocityTangentialDispersion       , separations
    double precision                                            , pointer      , dimension(:  ) :: massVirial
    double precision                                            , pointer      , dimension(:,:) :: position1                          , velocity1                        , &
         &                                                                                         position2                          , velocity2
    integer         (c_size_t                                  ), allocatable  , dimension(:  ) :: pairCountRadial                    , pairCountTangential
    integer         (c_size_t                                  ), allocatable  , dimension(:,:) :: weight1                            , weight2                          , &
         &                                                                                         weightPair
    logical                                                     , allocatable  , dimension(:  ) :: maskRadial                         , maskTangential
    double precision                                            , allocatable  , dimension(:,:) :: positions                          , velocities                       , &
         &                                                                                         distributionVelocityRadialBin      , distributionVelocityTangentialBin
    type            (treeNode                                  ), pointer                       :: node
    class           (nodeComponentBasic                        ), pointer                       :: basic
    integer         (c_size_t                                  )                                :: i                                  , j                                , &
         &                                                                                         iSample                            , bootstrapSampleCount             , &
         &                                                                                         iSimulation                        , jSimulation                      , &
         &                                                                                         jStart                             , jEnd
    type            (nearestNeighbors                          )                                :: neighborFinder
    integer                                                                                     :: neighborCount
    type            (varying_string                            )                                :: label
    double precision                                                                            :: separationMinimum                  , separationMaximum                , &
         &                                                                                         radiusVirial                       , velocityVirial                   , &
         &                                                                                         velocityWidthBin

#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayIndent('compute virial crossing orbit statistics',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    ! Construct bins of separation.
    bootstrapSampleCount=self%bootstrapSampleCount
    if (self%includeUnbootstrapped) bootstrapSampleCount=bootstrapSampleCount+1
    allocate(velocityCentralBin               (self%velocityCount                     ))
    allocate(velocityMinimumBin               (self%velocityCount                     ))
    allocate(velocityMaximumBin               (self%velocityCount                     ))
    allocate(distributionVelocityRadialBin    (self%velocityCount,bootstrapSampleCount))
    allocate(distributionVelocityTangentialBin(self%velocityCount,bootstrapSampleCount))
    allocate(pairCountRadial                  (                   bootstrapSampleCount))
    allocate(pairCountTangential              (                   bootstrapSampleCount))
    allocate(velocityRadialMean               (                   bootstrapSampleCount))
    allocate(velocityTangentialMean           (                   bootstrapSampleCount))
    allocate(velocityRadialDispersion         (                   bootstrapSampleCount))
    allocate(velocityTangentialDispersion     (                   bootstrapSampleCount))
    velocityCentralBin=Make_Range(self%velocityMinimum,self%velocityMaximum,int(self%velocityCount),rangeTypeLinear,rangeBinned=.true.)
    velocityWidthBin  =velocityCentralBin(2)-velocityCentralBin(1)
    velocityMinimumBin=velocityCentralBin-0.5d0*velocityWidthBin
    velocityMaximumBin=velocityCentralBin+0.5d0*velocityWidthBin
    ! Iterate over simulations.
    do iSimulation=1_c_size_t,size(simulations)
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call displayMessage(var_str('simulation "')//simulations(iSimulation)%label//'"',verbosityLevelStandard)
#ifdef USEMPI
       end if
#endif
       ! Determine reference simulation.
       if (self%crossCount) then
          jSimulation=1
       else
          jSimulation=iSimulation
       end if
       position1 => simulations(jSimulation)%propertiesRealRank1%value('position')
       velocity1 => simulations(jSimulation)%propertiesRealRank1%value('velocity')
       position2 => simulations(iSimulation)%propertiesRealRank1%value('position')
       velocity2 => simulations(iSimulation)%propertiesRealRank1%value('velocity')
       ! Get halo virial masses.
       if (simulations(jSimulation)%propertiesReal%exists('massVirial')) then
          massVirial => simulations(jSimulation)%propertiesReal%value('massVirial')
       else
          call Error_Report('halo virial masses are required, but are not available in the simulation'//{introspection:location})
       end if
       ! Generate bootstrap weights.
       if (self%crossCount) then
          allocate(weight1(bootstrapSampleCount,size(position1,dim=2)))
       else
          allocate(weight1(bootstrapSampleCount,size(position2,dim=2)))
       end if
       allocate   (weight2(bootstrapSampleCount,size(position2,dim=2)))
       do iSample=1,bootstrapSampleCount
          ! Determine weights for particles.
          if (iSample == 1 .and. self%includeUnbootstrapped) then
             weight1(iSample,:)=1_c_size_t
             weight2(iSample,:)=1_c_size_t
          else
             do i=1,size(weight1,dim=2)
                weight1(iSample,i)=int(self%randomNumberGenerator_%poissonSample(self%bootstrapSampleRate),c_size_t)
             end do
             if (self%crossCount) then
                do i=1,size(weight2,dim=2)
                   weight2(iSample,i)=int(self%randomNumberGenerator_%poissonSample(self%bootstrapSampleRate),c_size_t)
                end do
             else
                weight2(iSample,:)=weight1(iSample,:)
             end if
          end if
       end do
       ! Accumulate pairs.
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call displayCounter(0,.true.)
#ifdef USEMPI
       end if
#endif
       distributionVelocityRadialBin    =0.0d0
       distributionVelocityTangentialBin=0.0d0
       pairCountRadial                  =0_c_size_t
       pairCountTangential              =0_c_size_t
       velocityRadialMean               =0.0d0
       velocityTangentialMean           =0.0d0
       velocityRadialDispersion         =0.0d0
       velocityTangentialDispersion     =0.0d0
       !$omp parallel private(iSample,j,jStart,jEnd,node,basic,weightPair,maskRadial,maskTangential,separations,positions,velocities,velocitiesRadial,velocitiesTangential,separationsVirial,angularMomentumVirial,energyVirial,neighborCount,neighborIndex,neighborDistanceSquared,neighborFinder) reduction(+:distributionVelocityRadialBin,distributionVelocityTangentialBin,pairCountRadial,pairCountTangential,velocityRadialMean,velocityTangentialMean,velocityRadialDispersion,velocityTangentialDispersion)
       call Node_Components_Thread_Initialize(self%parameters)
       ! Allocate workspace.
       allocate(weightPair             (bootstrapSampleCount,size(position2,dim=2)))
       allocate(positions              (3                   ,size(position2,dim=2)))
       allocate(velocities             (3                   ,size(position2,dim=2)))
       allocate(separations            (                     size(position2,dim=2)))
       allocate(velocitiesRadial       (                     size(position2,dim=2)))
       allocate(velocitiesTangential   (                     size(position2,dim=2)))
       allocate(neighborIndex          (                     size(position2,dim=2)))
       allocate(neighborDistanceSquared(                     size(position2,dim=2)))
       allocate(maskRadial             (                     size(position2,dim=2)))
       allocate(maskTangential         (                     size(position2,dim=2)))
       allocate(separationsVirial      (                     size(position2,dim=2)))
       allocate(angularMomentumVirial  (                     size(position2,dim=2)))
       allocate(energyVirial           (                     size(position2,dim=2)))
       ! Build a work node.
       node  => treeNode      (                 )
       basic => node    %basic(autoCreate=.true.)
       call basic%timeSet            (self%time)
       call basic%timeLastIsolatedSet(self%time)
       ! Construct the nearest neighbor finder object.
       neighborFinder=nearestNeighbors(transpose(position2))
       ! Iterate over particles.
       !$omp do schedule(dynamic)
       do i=1_c_size_t,size(position1,dim=2,kind=c_size_t)
#ifdef USEMPI
          ! If running under MPI with N processes, process only every Nth particle.
          if (mod(i,mpiSelf%count()) /= mpiSelf%rank()) cycle
#endif
          ! Get virial properties of the target halo.
          call basic%massSet(massVirial(i))
          call Calculations_Reset(node)
          radiusVirial  =self%darkMatterHaloScale_%radiusVirial  (node)
          velocityVirial=self%darkMatterHaloScale_%velocityVirial(node)
          ! Determine radial search range.
          separationMinimum=radiusVirial*self%separationMinimum
          separationMaximum=radiusVirial*self%separationMaximum
          ! Locate particles nearby.
          call neighborFinder%searchFixedRadius(position1(:,i),separationMaximum,toleranceZero,neighborCount,neighborIndex,neighborDistanceSquared)
          ! Find weights of paired particles.
          do j=1_c_size_t,neighborCount
             weightPair(:,j)=weight2  (:,neighborIndex(j))
             positions (:,j)=position2(:,neighborIndex(j))
             velocities(:,j)=velocity2(:,neighborIndex(j))
          end do
          ! Compute velocities.
          do j=1,3
             positions (j,1:neighborCount)=positions (j,1:neighborCount)-position1(j,i)
             velocities(j,1:neighborCount)=velocities(j,1:neighborCount)-velocity1(j,i)
          end do
          separations         (1:neighborCount)=+sqrt(                                      &
               &                                      sum(                                  &
               &                                           positions(:,1:neighborCount)**2, &
               &                                           dim=1                            &
               &                                          )                                 &
               &                                     )
          where (separations(1:neighborCount) > 0.0d0)
             velocitiesRadial    (1:neighborCount)=+    (                                                                                                                                 &
                  &                                        +velocities(1,1:neighborCount)                                  *positions(1,1:neighborCount)                                  &
                  &                                        +velocities(2,1:neighborCount)                                  *positions(2,1:neighborCount)                                  &
                  &                                        +velocities(3,1:neighborCount)                                  *positions(3,1:neighborCount)                                  &
                  &                                     )                                                                                                                                 &
                  &                                /                                                                                                     separations(1:neighborCount)
             velocitiesTangential(1:neighborCount)=+sqrt(                                                                                                                                 &
                  &                                      +(+velocities(1,1:neighborCount)-velocitiesRadial(1:neighborCount)*positions(1,1:neighborCount)/separations(1:neighborCount))**2 &
                  &                                      +(+velocities(2,1:neighborCount)-velocitiesRadial(1:neighborCount)*positions(2,1:neighborCount)/separations(1:neighborCount))**2 &
                  &                                      +(+velocities(3,1:neighborCount)-velocitiesRadial(1:neighborCount)*positions(3,1:neighborCount)/separations(1:neighborCount))**2 &
                  &                                     )
          elsewhere
             velocitiesRadial    (1:neighborCount)=+0.0d0
             velocitiesTangential(1:neighborCount)=+sqrt(sum(velocities(:,1:neighborCount)**2,dim=1))
          end where
          ! Add Hubble flow if requested.
          if (self%addHubbleFlow) velocitiesRadial(1:neighborCount)=velocitiesRadial(1:neighborCount)+self%cosmologyFunctions_%hubbleParameterEpochal(self%time)*separations(1:neighborCount)
          ! Scale velocities to the virial velocity for the halo. Change the sign of radial velocities here so that infalling
          ! halos have positive radial velocity.
          velocitiesRadial    (1:neighborCount)=-velocitiesRadial    (1:neighborCount)/velocityVirial
          velocitiesTangential(1:neighborCount)=+velocitiesTangential(1:neighborCount)/velocityVirial
          ! Propagate velocities to the virial radius.
          !! Separation in virial units.
          separationsVirial(1:neighborCount)=+separations (1:neighborCount) &
               &                             /radiusVirial
          where (separationsVirial(1:neighborCount) > 0.0d0 .and. velocitiesRadial(1:neighborCount) >= 0.0d0)
             !! Angular momentum in virial units.
             angularMomentumVirial(1:neighborCount)=+velocitiesTangential   (1:neighborCount)    &
                  &                                 *separationsVirial      (1:neighborCount)
             !! Energy in virial units.
             energyVirial         (1:neighborCount)=-1.0d0                                       &
                  &                                 /  separationsVirial    (1:neighborCount)    &
                  &                                 +0.5d0                                       &
                  &                                 *(                                           &
                  &                                   +velocitiesRadial     (1:neighborCount)**2 &
                  &                                   +velocitiesTangential (1:neighborCount)**2 &
                  &                                  )
             !! Propagate tangential velocity to virial radius.
             velocitiesTangential (1:neighborCount)=+  angularMomentumVirial(1:neighborCount)
             !! Propagate radial velocity to virial radius. This is actually the square of the radial velocity. For orbits which
             !! have a pericenter outside the virial radius this will be negative.
             velocitiesRadial     (1:neighborCount)=+2.0d0                                       &
                  &                                 *(                                           &
                  &                                   +energyVirial         (1:neighborCount)    &
                  &                                   +1.0d0                                     &
                  &                                   -0.5d0                                     &
                  &                                   *velocitiesTangential (1:neighborCount)**2 &
                  &                                  )
          elsewhere
             velocitiesTangential (1:neighborCount)=+0.0d0
             velocitiesRadial     (1:neighborCount)=+0.0d0
          end where
          ! Convert from square of radial velocity to radial velocity. In cases where the square of the radial velocity is
          ! negative (i.e. the orbit does not cross the virial radius) set both radial and tangential velocities to negative
          ! infinity so that they do not contribute to the statistics.
          where (velocitiesRadial(1:neighborCount) >= 0.0d0)
             velocitiesRadial    (1:neighborCount)=+sqrt(velocitiesRadial(1:neighborCount))
          elsewhere
             velocitiesRadial    (1:neighborCount)=-huge(0.0d0)
             velocitiesTangential(1:neighborCount)=-huge(0.0d0)
          end where
          ! Accumulate particles into bins.
          if (neighborDistanceSquared(1) <= separationMaximum**2 .and. neighborDistanceSquared(neighborCount) >= separationMinimum**2) then
             jStart=searchArray(neighborDistanceSquared(1:neighborCount),separationMinimum**2)
             jEnd  =neighborCount
             ! Accumulate statistics.
             maskRadial    (1:jEnd-jStart+1)= velocitiesRadial       (jStart:jEnd) <  self%velocityCut      &
                  &                          .and.                                                          &
                  &                           velocitiesRadial       (jStart:jEnd) >  0.0d0                 &
                  &                          .and.                                                          &
                  &                           neighborDistanceSquared(jStart:jEnd) >  separationMinimum**2
             maskTangential(1:jEnd-jStart+1)= velocitiesTangential   (jStart:jEnd) <  self%velocityCut      &
                  &                          .and.                                                          &
                  &                           velocitiesRadial       (jStart:jEnd) >  0.0d0                 &
                  &                          .and.                                                          &
                  &                           neighborDistanceSquared(jStart:jEnd) >  separationMinimum**2
             do iSample=1,bootstrapSampleCount
                if (weight1(iSample,i) == 0_c_size_t) cycle
                pairCountRadial             (iSample)=pairCountRadial             (iSample)+     weight1(iSample,i) *sum(     weightPair(iSample,jStart:jEnd)                                      ,maskRadial    (1:jEnd-jStart+1))
                pairCountTangential         (iSample)=pairCountTangential         (iSample)+     weight1(iSample,i) *sum(     weightPair(iSample,jStart:jEnd)                                      ,maskTangential(1:jEnd-jStart+1))
                velocityRadialMean          (iSample)=velocityRadialMean          (iSample)+dble(weight1(iSample,i))*sum(dble(weightPair(iSample,jStart:jEnd))*velocitiesRadial    (jStart:jEnd)   ,maskRadial    (1:jEnd-jStart+1))
                velocityRadialDispersion    (iSample)=velocityRadialDispersion    (iSample)+dble(weight1(iSample,i))*sum(dble(weightPair(iSample,jStart:jEnd))*velocitiesRadial    (jStart:jEnd)**2,maskRadial    (1:jEnd-jStart+1))
                velocityTangentialMean      (iSample)=velocityTangentialMean      (iSample)+dble(weight1(iSample,i))*sum(dble(weightPair(iSample,jStart:jEnd))*velocitiesTangential(jStart:jEnd)   ,maskTangential(1:jEnd-jStart+1))
                velocityTangentialDispersion(iSample)=velocityTangentialDispersion(iSample)+dble(weight1(iSample,i))*sum(dble(weightPair(iSample,jStart:jEnd))*velocitiesTangential(jStart:jEnd)**2,maskTangential(1:jEnd-jStart+1))
             end do
             ! Accumulate distributions.
             do j=1,self%velocityCount
                maskRadial    (1:jEnd-jStart+1)= velocitiesRadial       (jStart:jEnd) >= velocityMinimumBin(j) &
                     &                          .and.                                                          &
                     &                           velocitiesRadial       (jStart:jEnd) <  velocityMaximumBin(j) &
                     &                          .and.                                                          &
                     &                           velocitiesRadial       (jStart:jEnd) <  self%velocityCut      &
                     &                          .and.                                                          &
                     &                           velocitiesRadial       (jStart:jEnd) >  0.0d0                 &
                     &                          .and.                                                          &
                     &                           neighborDistanceSquared(jStart:jEnd) >  separationMinimum**2
                maskTangential(1:jEnd-jStart+1)= velocitiesTangential   (jStart:jEnd) >= velocityMinimumBin(j) &
                     &                          .and.                                                          &
                     &                           velocitiesTangential   (jStart:jEnd) <  velocityMaximumBin(j) &
                     &                          .and.                                                          &
                     &                           velocitiesTangential   (jStart:jEnd) <  self%velocityCut      &
                     &                          .and.                                                          &
                     &                           velocitiesRadial       (jStart:jEnd) >  0.0d0                 &
                     &                          .and.                                                          &
                     &                           neighborDistanceSquared(jStart:jEnd) >  separationMinimum**2
                do iSample=1,bootstrapSampleCount
                   if (weight1(iSample,i) == 0_c_size_t) cycle
                   distributionVelocityRadialBin    (j,iSample)=+distributionVelocityRadialBin    (j,iSample)+dble(weight1(iSample,i))*sum(dble(weightPair(iSample,jStart:jEnd))*velocitiesRadial(jStart:jEnd),maskRadial    (1:jEnd-jStart+1))
                   distributionVelocityTangentialBin(j,iSample)=+distributionVelocityTangentialBin(j,iSample)+dble(weight1(iSample,i))*sum(dble(weightPair(iSample,jStart:jEnd))*velocitiesRadial(jStart:jEnd),maskTangential(1:jEnd-jStart+1))
                end do
             end do
          end if
          ! Update progress.
          !$ if (OMP_Get_Thread_Num() == 0) then
#ifdef USEMPI
          if (mpiSelf%isMaster()) then
#endif
             call displayCounter(                                                 &
                  &                          int(                                             &
                  &                              +100.0d0                                     &
                  &                              *float(i                                  )  &
                  &                              /float(size(position1,dim=2,kind=c_size_t))  &
                  &                             )                                           , &
                  &                          .false.                                          &
                  &                         )
#ifdef USEMPI
          end if
#endif
          !$ end if
       end do
       !$omp end do
       call node%destroy()
       call Node_Components_Thread_Uninitialize()
       !$omp end parallel
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call displayCounterClear()
#ifdef USEMPI
       end if
#endif
       deallocate(weight1   )
       deallocate(weight2   )
#ifdef USEMPI
       ! Reduce across MPI processes.
       distributionVelocityRadialBin    =mpiSelf%sum(distributionVelocityRadialBin    )
       distributionVelocityTangentialBin=mpiSelf%sum(distributionVelocityTangentialBin)
       pairCountRadial                  =mpiSelf%sum(pairCountRadial                  )
       pairCountTangential              =mpiSelf%sum(pairCountTangential              )
       velocityRadialMean               =mpiSelf%sum(velocityRadialMean               )
       velocityRadialDispersion         =mpiSelf%sum(velocityRadialDispersion         )
       velocityTangentialMean           =mpiSelf%sum(velocityTangentialMean           )
       velocityTangentialDispersion     =mpiSelf%sum(velocityTangentialDispersion     )
#endif
       ! Finalize statistics.
       where (pairCountRadial     > 0_c_size_t)
          velocityRadialMean          =               velocityRadialMean          /dble(pairCountRadial    )
          velocityRadialDispersion    =sqrt(max(0.0d0,velocityRadialDispersion    /dble(pairCountRadial    )-velocityRadialMean    **2))
       end where
       where (pairCountTangential > 0_c_size_t)
          velocityTangentialMean      =               velocityTangentialMean      /dble(pairCountTangential)
          velocityTangentialDispersion=sqrt(max(0.0d0,velocityTangentialDispersion/dble(pairCountTangential)-velocityTangentialMean**2))
       end where
       ! Normalize distributions.
       do i=1,bootstrapSampleCount
          if (sum(distributionVelocityRadialBin    (:,i)) > 0.0d0) distributionVelocityRadialBin    (:,i)=distributionVelocityRadialBin    (:,i)/sum(distributionVelocityRadialBin    (:,i))/velocityWidthBin
          if (sum(distributionVelocityTangentialBin(:,i)) > 0.0d0) distributionVelocityTangentialBin(:,i)=distributionVelocityTangentialBin(:,i)/sum(distributionVelocityTangentialBin(:,i))/velocityWidthBin
       end do
       ! Output results.
       if (self%crossCount) then
          jSimulation   = 1
          label         = ":"//simulations(          1)%label// &
               &          ":"//simulations(iSimulation)%label
       else
          jSimulation   = iSimulation
          label         = ":"//simulations(iSimulation)%label
       end if
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          !$ call hdf5Access%  set()
          if (.not.self%crossCount .or. iSimulation == 1 ) &
               & call simulations(jSimulation)%analysis%writeDataset(velocityCentralBin               ,'virialCrossingOrbitVelocity'                           )
          call        simulations(jSimulation)%analysis%writeDataset(distributionVelocityRadialBin    ,'virialCrossingOrbitDistributionRadial'    //char(label))
          call        simulations(jSimulation)%analysis%writeDataset(distributionVelocityTangentialBin,'virialCrossingOrbitDistributionTangential'//char(label))
          call        simulations(jSimulation)%analysis%writeDataset(velocityRadialMean               ,'velocityRadialMean'                       //char(label))
          call        simulations(jSimulation)%analysis%writeDataset(velocityRadialDispersion         ,'velocityRadialDispersion'                 //char(label))
          call        simulations(jSimulation)%analysis%writeDataset(velocityTangentialMean           ,'velocityTangentialMean'                   //char(label))
          call        simulations(jSimulation)%analysis%writeDataset(velocityTangentialDispersion     ,'velocityTangentialDispersion'             //char(label))
          !$ call hdf5Access%unset()
#ifdef USEMPI
       end if
#endif
    end do
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayUnindent('done',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    return
  end subroutine virialCrossingOrbitStatisticsOperate

