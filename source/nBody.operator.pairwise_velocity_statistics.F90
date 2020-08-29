!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which implements an N-body data operator which computes pairwise velocity statistics in bins of separation.

  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
  use            :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use            :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass

  !# <nbodyOperator name="nbodyOperatorPairwiseVelocityStatistics">
  !#  <description>An N-body data operator which computes pairwise velocity statistics in bins of separation.</description>
  !# </nbodyOperator>
  type, extends(nbodyOperatorClass) :: nbodyOperatorPairwiseVelocityStatistics
     !% An N-body data operator which computes pairwise velocity statistics in bins of separation.
     private
     double precision                                      :: separationMinimum               , separationMaximum   , &
          &                                                   bootstrapSampleRate             , time
     integer         (c_size_t                  )          :: separationCount                 , bootstrapSampleCount
     logical                                               :: includeUnbootstrapped           , crossCount          , &
          &                                                   addHubbleFlow
     class           (randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()
     class           (darkMatterHaloScaleClass  ), pointer :: darkMatterHaloScale_   => null()
     class           (cosmologyFunctionsClass   ), pointer :: cosmologyFunctions_    => null()
     type            (inputParameters           )          :: parameters
   contains
     final     ::            pairwiseVelocityStatisticsDestructor
     procedure :: operate => pairwiseVelocityStatisticsOperate
  end type nbodyOperatorPairwiseVelocityStatistics

  interface nbodyOperatorPairwiseVelocityStatistics
     !% Constructors for the ``pairwiseVelocityStatistics'' N-body operator class.
     module procedure pairwiseVelocityStatisticsConstructorParameters
     module procedure pairwiseVelocityStatisticsConstructorInternal
  end interface nbodyOperatorPairwiseVelocityStatistics

contains

  function pairwiseVelocityStatisticsConstructorParameters(parameters) result (self)
    !% Constructor for the ``pairwiseVelocityStatistics'' N-body operator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorPairwiseVelocityStatistics)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass             ), pointer       :: randomNumberGenerator_
    class           (cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass               ), pointer       :: darkMatterHaloScale_
    type            (inputParameters                        ), pointer       :: parametersRoot
    double precision                                                         :: separationMinimum     , separationMaximum   , &
         &                                                                      bootstrapSampleRate   , redshift            , &
         &                                                                      time
    integer         (c_size_t                               )                :: separationCount       , bootstrapSampleCount
    logical                                                                  :: includeUnbootstrapped , crossCount          , &
         &                                                                      addHubbleFlow

    !# <inputParameter>
    !#   <name>crossCount</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>If true, compute cross-simulation pairwise velocity statistics between the first and all simulations. Otherwise, compute pairwise velocity statistics within each simulation.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>addHubbleFlow</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>If true, add Hubble flow to velocities.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>redshift</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The redshift.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>bootstrapSampleCount</name>
    !#   <source>parameters</source>
    !#   <defaultValue>30_c_size_t</defaultValue>
    !#   <description>The number of bootstrap resamples of the particles that should be used.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>bootstrapSampleRate</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The sampling rate for particles.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>separationMinimum</name>
    !#   <source>parameters</source>
    !#   <description>The minimum separation to consider for pairwise velocity statistics.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>separationMaximum</name>
    !#   <source>parameters</source>
    !#   <description>The maximum separation to consider for pairwise velocity statistics.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>separationCount</name>
    !#   <source>parameters</source>
    !#   <description>The number of bins in separation for pairwise velocity statistics.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>includeUnbootstrapped</name>
    !#   <source>parameters</source>
    !#   <description>If true, include results for the unbootstrapped (i.e. original) sample.</description>
    !#   <defaultValue>.true.</defaultValue>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale"   name="darkMatterHaloScale_"   source="parameters"/>
    time=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift))
    if (associated(parameters%parent)) then
       parametersRoot => parameters%parent
       do while (associated(parametersRoot%parent))
          parametersRoot => parametersRoot%parent
       end do
       self=nbodyOperatorPairwiseVelocityStatistics(separationMinimum,separationMaximum,separationCount,time,crossCount,addHubbleFlow,includeUnbootstrapped,bootstrapSampleCount,bootstrapSampleRate,randomNumberGenerator_,cosmologyFunctions_,darkMatterHaloScale_,parametersRoot)
    else
       self=nbodyOperatorPairwiseVelocityStatistics(separationMinimum,separationMaximum,separationCount,time,crossCount,addHubbleFlow,includeUnbootstrapped,bootstrapSampleCount,bootstrapSampleRate,randomNumberGenerator_,cosmologyFunctions_,darkMatterHaloScale_,parameters    )
    end if
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="randomNumberGenerator_"/>
    !# <objectDestructor name="cosmologyFunctions_"   />
    !# <objectDestructor name="darkMatterHaloScale_"  />
    return
  end function pairwiseVelocityStatisticsConstructorParameters

  function pairwiseVelocityStatisticsConstructorInternal(separationMinimum,separationMaximum,separationCount,time,crossCount,addHubbleFlow,includeUnbootstrapped,bootstrapSampleCount,bootstrapSampleRate,randomNumberGenerator_,cosmologyFunctions_,darkMatterHaloScale_,parameters) result (self)
    !% Internal constructor for the ``pairwiseVelocityStatistics'' N-body operator class.
    implicit none
    type            (nbodyOperatorPairwiseVelocityStatistics)                        :: self
    double precision                                         , intent(in   )         :: separationMinimum     , separationMaximum   , &
         &                                                                              bootstrapSampleRate   , time
    integer         (c_size_t                               ), intent(in   )         :: separationCount       , bootstrapSampleCount
    logical                                                  , intent(in   )         :: includeUnbootstrapped , crossCount          , &
         &                                                                              addHubbleFlow
    class           (randomNumberGeneratorClass             ), intent(in   ), target :: randomNumberGenerator_
    class           (darkMatterHaloScaleClass               ), intent(in   ), target :: darkMatterHaloScale_
    class           (cosmologyFunctionsClass                ), intent(in   ), target :: cosmologyFunctions_
    type            (inputParameters                        ), intent(in   ), target :: parameters
    !# <constructorAssign variables="separationMinimum, separationMaximum, separationCount, time, crossCount, addHubbleFlow, includeUnbootstrapped, bootstrapSampleCount, bootstrapSampleRate, *randomNumberGenerator_, *cosmologyFunctions_, *darkMatterHaloScale_"/>

    self%parameters=inputParameters(parameters)
    return
  end function pairwiseVelocityStatisticsConstructorInternal

  subroutine pairwiseVelocityStatisticsDestructor(self)
    !% Destructor for the ``pairwiseVelocityStatistics'' N-body operator class.
    implicit none
    type(nbodyOperatorPairwiseVelocityStatistics), intent(inout) :: self

    !# <objectDestructor name="self%randomNumberGenerator_"/>
    !# <objectDestructor name="self%cosmologyFunctions_"   />
    !# <objectDestructor name="self%darkMatterHaloScale_"  />
    return
  end subroutine pairwiseVelocityStatisticsDestructor
  
  subroutine pairwiseVelocityStatisticsOperate(self,simulations)
    !% Compute pairwise velocity statistics in bins of separation.
    !$ use :: OMP_Lib           , only : OMP_Get_Thread_Num
    use    :: Arrays_Search                 , only : searchArray
    use    :: Galacticus_Calculations_Resets, only : Galacticus_Calculations_Reset
    use    :: Galacticus_Display            , only : Galacticus_Display_Indent        , Galacticus_Display_Unindent, Galacticus_Display_Counter, Galacticus_Display_Counter_Clear, &
         &                                           Galacticus_Display_Message       , verbosityStandard
    use    :: Galacticus_Nodes              , only : treeNode                         , nodeComponentBasic
    use    :: IO_HDF5                       , only : hdf5Access
    use    :: ISO_Varying_String            , only : var_str
    use    :: Memory_Management             , only : deallocateArray
    use    :: Nearest_Neighbors             , only : nearestNeighbors
    use    :: Node_Components               , only : Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize
    use    :: Numerical_Ranges              , only : Make_Range                       , rangeTypeLogarithmic
#ifdef USEMPI
    use    :: MPI_Utilities                 , only : mpiSelf
#endif
    implicit none
    class           (nbodyOperatorPairwiseVelocityStatistics), intent(inout)                 :: self
    type            (nBodyData                              ), intent(inout), dimension(:  ) :: simulations
    double precision                                         , parameter                     :: toleranceZero                          =0.0d0
    integer                                                  , allocatable  , dimension(:  ) :: neighborIndex
    double precision                                         , allocatable  , dimension(:  ) :: neighborDistanceSquared                      , separationCentralBin               , &
         &                                                                                      separationSquaredMinimumBin                  , separationSquaredMaximumBin        , &
         &                                                                                      velocitiesRadial                             , separations                        , &
         &                                                                                      velocitiesTangential                         , massVirial                         , &
         &                                                                                      velocityRadialVirialSquared
    integer         (c_size_t                               ), allocatable  , dimension(:,:) :: weight1                                      , weight2                            , &
         &                                                                                      weightPair
    logical                                                  , allocatable  , dimension(:  ) :: mask
    integer         (c_size_t                               ), allocatable  , dimension(:,:) :: pairCountBin
    double precision                                         , allocatable  , dimension(:,:) :: velocityRadialBin                            , positions                          , &
         &                                                                                      velocities                                   , velocityDispersionRadialBin        , &
         &                                                                                      velocityDispersionTangentialBin              , velocityDispersionRadialWeightedBin, &
         &                                                                                      velocityDispersionTangentialWeightedBin      , velocityRadialWeightedBin          , &
         &                                                                                      pairCountWeightedBin
    type            (treeNode                               ), pointer                       :: node
    class           (nodeComponentBasic                     ), pointer                       :: basic
    double precision                                                                         :: radiusVirial                                 , velocityVirial
    integer         (c_size_t                               )                                :: i                                            , j                                  , &
         &                                                                                      iSample                                      , bootstrapSampleCount               , &
         &                                                                                      iSimulation                                  , jSimulation                        , &
         &                                                                                      jStart                                       , jEnd
    type            (nearestNeighbors                       )                                :: neighborFinder
    integer                                                                                  :: neighborCount
    type            (varying_string                         )                                :: label
    
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call Galacticus_Display_Indent('compute pairwise velocity statistics',verbosityStandard)
#ifdef USEMPI
    end if
#endif
    ! Construct bins of separation.
    bootstrapSampleCount=self%bootstrapSampleCount
    if (self%includeUnbootstrapped) bootstrapSampleCount=bootstrapSampleCount+1
    allocate(separationCentralBin                   (self%separationCount                     ))
    allocate(separationSquaredMinimumBin            (self%separationCount                     ))
    allocate(separationSquaredMaximumBin            (self%separationCount                     ))
    allocate(pairCountBin                           (self%separationCount,bootstrapSampleCount))
    allocate(velocityRadialBin                      (self%separationCount,bootstrapSampleCount))
    allocate(velocityDispersionRadialBin            (self%separationCount,bootstrapSampleCount))
    allocate(velocityDispersionTangentialBin        (self%separationCount,bootstrapSampleCount))
    allocate(pairCountWeightedBin                   (self%separationCount,bootstrapSampleCount))
    allocate(velocityRadialWeightedBin              (self%separationCount,bootstrapSampleCount))
    allocate(velocityDispersionRadialWeightedBin    (self%separationCount,bootstrapSampleCount))
    allocate(velocityDispersionTangentialWeightedBin(self%separationCount,bootstrapSampleCount))    
    separationCentralBin       =Make_Range(self%separationMinimum,self%separationMaximum,int(self%separationCount),rangeTypeLogarithmic,rangeBinned=.true.)
    separationSquaredMinimumBin=(separationCentralBin/sqrt(separationCentralBin(2)/separationCentralBin(1)))**2
    separationSquaredMaximumBin=(separationCentralBin*sqrt(separationCentralBin(2)/separationCentralBin(1)))**2
    ! Iterate over simulations.
    do iSimulation=1_c_size_t,size(simulations)
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call Galacticus_Display_Message(var_str('simulation "')//simulations(iSimulation)%label//'"',verbosityStandard)
#ifdef USEMPI
       end if
#endif
       ! Determine reference simulation.
       if (self%crossCount) then
          jSimulation=1
       else
          jSimulation=iSimulation
       end if
       ! Get halo virial masses.
       if (simulations(jSimulation)%propertiesReal%exists('massVirial')) then
          allocate(massVirial(size(simulations(jSimulation)%position,dim=2)))
          massVirial=simulations(jSimulation)%propertiesReal%value('massVirial')
       else
          call Galacticus_Error_Report('halo virial masses are required, but are not available in the simulation'//{introspection:location})
       end if
       ! Generate bootstrap weights.
       if (self%crossCount) then
          allocate(weight1(bootstrapSampleCount,size(simulations(          1)%position,dim=2)))
       else
          allocate(weight1(bootstrapSampleCount,size(simulations(iSimulation)%position,dim=2)))
       end if
       allocate   (weight2(bootstrapSampleCount,size(simulations(iSimulation)%position,dim=2)))
       ! Generate bootstrap weights.
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
          call Galacticus_Display_Counter(0,.true.)
#ifdef USEMPI
       end if
#endif
       pairCountBin                           =0_c_size_t
       velocityRadialBin                      =0.0d0
       velocityDispersionRadialBin            =0.0d0
       velocityDispersionTangentialBin        =0.0d0
       pairCountWeightedBin                   =0.0d0
       velocityRadialWeightedBin              =0.0d0
       velocityDispersionRadialWeightedBin    =0.0d0
       velocityDispersionTangentialWeightedBin=0.0d0
       !$omp parallel private(iSample,j,jStart,jEnd,weightPair,mask,separations,positions,velocities,velocitiesRadial,velocitiesTangential,velocityRadialVirialSquared,neighborCount,neighborIndex,neighborDistanceSquared,neighborFinder) reduction(+:pairCountBin,velocityRadialBin,velocityDispersionRadialBin,velocityDispersionTangentialBin,pairCountWeightedBin,velocityRadialWeightedBin,velocityDispersionRadialWeightedBin,velocityDispersionTangentialWeightedBin)
       call Node_Components_Thread_Initialize(self%parameters)
       ! Allocate workspace.
       allocate(weightPair                 (bootstrapSampleCount,size(simulations(iSimulation)%position,dim=2)))
       allocate(positions                  (3                   ,size(simulations(iSimulation)%position,dim=2)))
       allocate(velocities                 (3                   ,size(simulations(iSimulation)%position,dim=2)))
       allocate(separations                (                     size(simulations(iSimulation)%position,dim=2)))
       allocate(velocitiesRadial           (                     size(simulations(iSimulation)%position,dim=2)))
       allocate(velocitiesTangential       (                     size(simulations(iSimulation)%position,dim=2)))
       allocate(velocityRadialVirialSquared(                     size(simulations(iSimulation)%position,dim=2)))
       allocate(neighborIndex              (                     size(simulations(iSimulation)%position,dim=2)))
       allocate(neighborDistanceSquared    (                     size(simulations(iSimulation)%position,dim=2)))
       allocate(mask                       (                     size(simulations(iSimulation)%position,dim=2)))
       ! Build a work node.
       node  => treeNode      (                 )
       basic => node    %basic(autoCreate=.true.)
       call basic%timeSet            (self%time)
       call basic%timeLastIsolatedSet(self%time)
       ! Construct the nearest neighbor finder object.
       neighborFinder=nearestNeighbors(transpose(simulations(iSimulation)%position))
       ! Iterate over particles.
       !$omp do schedule(dynamic)
       do i=1_c_size_t,size(simulations(jSimulation)%position,dim=2,kind=c_size_t)
#ifdef USEMPI
          ! If running under MPI with N processes, process only every Nth particle.
          if (mod(i,mpiSelf%count()) /= mpiSelf%rank()) cycle
#endif
          ! Get virial properties of the target halo.
          call basic%massSet(massVirial(i))
          call Galacticus_Calculations_Reset(node)
          radiusVirial  =self%darkMatterHaloScale_%virialRadius  (node)
          velocityVirial=self%darkMatterHaloScale_%virialVelocity(node)
          ! Locate particles nearby.
          call neighborFinder%searchFixedRadius(simulations(jSimulation)%position(:,i),self%separationMaximum,toleranceZero,neighborCount,neighborIndex,neighborDistanceSquared)
          ! Find weights of paired particles.
          do j=1_c_size_t,neighborCount
             weightPair(:,j)=                         weight2 (:,neighborIndex(j))
             positions (:,j)=simulations(iSimulation)%position(:,neighborIndex(j))
             velocities(:,j)=simulations(iSimulation)%velocity(:,neighborIndex(j))
          end do
          ! Compute velocities.
          do j=1,3
             positions (j,1:neighborCount)=positions (j,1:neighborCount)-simulations(jSimulation)%position(j,i)
             velocities(j,1:neighborCount)=velocities(j,1:neighborCount)-simulations(jSimulation)%velocity(j,i)
          end do
          separations(1:neighborCount)=+sqrt(                                      &
               &                             sum(                                  &
               &                                  positions(:,1:neighborCount)**2, &
               &                                  dim=1                            &
               &                                 )                                 &
               &                            )
          where (separations(1:neighborCount) > 0.0d0)
             velocitiesRadial           (1:neighborCount)=+    (                                                                                                                                 &
                  &                                               +velocities(1,1:neighborCount)                                  *positions(1,1:neighborCount)                                  &
                  &                                               +velocities(2,1:neighborCount)                                  *positions(2,1:neighborCount)                                  &
                  &                                               +velocities(3,1:neighborCount)                                  *positions(3,1:neighborCount)                                  &
                  &                                            )                                                                                                                                 &
                  &                                       /                                                                                                     separations(1:neighborCount)
             velocitiesTangential       (1:neighborCount)=+sqrt(                                                                                                                                 &
                  &                                             +(+velocities(1,1:neighborCount)-velocitiesRadial(1:neighborCount)*positions(1,1:neighborCount)/separations(1:neighborCount))**2 &
                  &                                             +(+velocities(2,1:neighborCount)-velocitiesRadial(1:neighborCount)*positions(2,1:neighborCount)/separations(1:neighborCount))**2 &
                  &                                             +(+velocities(3,1:neighborCount)-velocitiesRadial(1:neighborCount)*positions(3,1:neighborCount)/separations(1:neighborCount))**2 &
                  &                                            )
             ! Compute the square of the velocity at the virial radius (in virial units) assuming Keplerian orbits.
             velocityRadialVirialSquared(1:neighborCount)=+2.0d0                                      &
                  &                                       *(                                          &
                  &                                         +1.0d0                                    &
                  &                                         -radiusVirial                             &
                  &                                         /separations         (1:neighborCount)    &
                  &                                        )                                          &
                  &                                       +(                                          &
                  &                                         +velocitiesRadial    (1:neighborCount)**2 &
                  &                                         +velocitiesTangential(1:neighborCount)**2 &
                  &                                         *(                                        &
                  &                                           +1.0d0                                  &
                  &                                           -(                                      &
                  &                                             +separations     (1:neighborCount)    &
                  &                                             /radiusVirial                         &
                  &                                            )**2                                   &
                  &                                          )                                        &
                  &                                        )                                          &
                  &                                       /velocityVirial**2
          elsewhere
             velocitiesRadial           (1:neighborCount)=+0.0d0
             velocitiesTangential       (1:neighborCount)=+sqrt(sum(velocities(:,1:neighborCount)**2,dim=1))
             velocityRadialVirialSquared(1:neighborCount)=-huge(0.0d0)
          end where
          if (self%addHubbleFlow)                                                                                                      &
               & velocityRadialVirialSquared(1:neighborCount)=+                           velocityRadialVirialSquared(1:neighborCount) &
               &                                              +(                                                                       &
               &                                                +self%cosmologyFunctions_%hubbleParameterEpochal     (  self%time    ) &
               &                                                *                         separations                (1:neighborCount) &
               &                                                /                         velocityVirial                               &
               &                                               )**2
          ! Accumulate particles into bins.
          jEnd=0
          do j=1,self%separationCount
             if (neighborDistanceSquared(1) > separationSquaredMaximumBin(j) .or. neighborDistanceSquared(neighborCount) < separationSquaredMinimumBin(j)) cycle
             jStart=jEnd
             jEnd  =searchArray(neighborDistanceSquared(1:neighborCount),separationSquaredMaximumBin(j))
             if (jStart == jEnd) cycle
             do iSample=1,bootstrapSampleCount
                if (weight1(iSample,i) == 0_c_size_t) cycle
                mask(1:jEnd-jStart)= velocitiesRadial           (jStart+1:jEnd) <  self%cosmologyFunctions_%hubbleParameterEpochal(self%time)*separations(jStart+1:jEnd) &
                     &              .and.                                                                                                                                &
                     &               velocityRadialVirialSquared(jStart+1:jEnd) >= 0.0d0
                pairCountBin                           (j,iSample)=+pairCountBin                           (j,iSample)+     weight1(iSample,i) *sum(     weightPair(iSample,jStart+1:jEnd)                                                                                                 )
                velocityRadialBin                      (j,iSample)=+velocityRadialBin                      (j,iSample)+dble(weight1(iSample,i))*sum(dble(weightPair(iSample,jStart+1:jEnd))*velocitiesRadial    (jStart+1:jEnd)                                                            )
                velocityDispersionRadialBin            (j,iSample)=+velocityDispersionRadialBin            (j,iSample)+dble(weight1(iSample,i))*sum(dble(weightPair(iSample,jStart+1:jEnd))*velocitiesRadial    (jStart+1:jEnd)**2                                                         )
                velocityDispersionTangentialBin        (j,iSample)=+velocityDispersionTangentialBin        (j,iSample)+dble(weight1(iSample,i))*sum(dble(weightPair(iSample,jStart+1:jEnd))*velocitiesTangential(jStart+1:jEnd)**2                                                         )
                pairCountWeightedBin                   (j,iSample)=+pairCountWeightedBin                   (j,iSample)+dble(weight1(iSample,i))*sum(dble(weightPair(iSample,jStart+1:jEnd))                                       *abs(velocitiesRadial(jStart+1:jEnd)),mask(1:jEnd-jStart))
                velocityRadialWeightedBin              (j,iSample)=+velocityRadialWeightedBin              (j,iSample)+dble(weight1(iSample,i))*sum(dble(weightPair(iSample,jStart+1:jEnd))*velocitiesRadial    (jStart+1:jEnd)   *abs(velocitiesRadial(jStart+1:jEnd)),mask(1:jEnd-jStart))
                velocityDispersionRadialWeightedBin    (j,iSample)=+velocityDispersionRadialWeightedBin    (j,iSample)+dble(weight1(iSample,i))*sum(dble(weightPair(iSample,jStart+1:jEnd))*velocitiesRadial    (jStart+1:jEnd)**2*abs(velocitiesRadial(jStart+1:jEnd)),mask(1:jEnd-jStart))
                velocityDispersionTangentialWeightedBin(j,iSample)=+velocityDispersionTangentialWeightedBin(j,iSample)+dble(weight1(iSample,i))*sum(dble(weightPair(iSample,jStart+1:jEnd))*velocitiesTangential(jStart+1:jEnd)**2*abs(velocitiesRadial(jStart+1:jEnd)),mask(1:jEnd-jStart))
             end do
          end do
          ! Update progress.
          !$ if (OMP_Get_Thread_Num() == 0) then
#ifdef USEMPI
          if (mpiSelf%isMaster()) then
#endif
             call Galacticus_Display_Counter(                                                               &
                  &                          int(                                                           &
                  &                              +100.0d0                                                   &
                  &                              *float(i                                                )  &
                  &                              /float(size(simulations(1)%position,dim=2,kind=c_size_t))  &
                  &                             )                                                         , &
                  &                          .false.                                                        &
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
          call Galacticus_Display_Counter_Clear()
#ifdef USEMPI
       end if
#endif
       deallocate(weight1   )
       deallocate(weight2   )
       deallocate(massVirial)
#ifdef USEMPI
       ! Reduce across MPI processes.
       pairCountBin                           =mpiSelf%sum(pairCountBin                           )
       velocityRadialBin                      =mpiSelf%sum(velocityRadialBin                      )
       velocityDispersionRadialBin            =mpiSelf%sum(velocityDispersionRadialBin            )
       velocityDispersionTangentialBin        =mpiSelf%sum(velocityDispersionTangentialBin        )
       pairCountWeightedBin                   =mpiSelf%sum(pairCountWeightedBin                   )
       velocityRadialWeightedBin              =mpiSelf%sum(velocityRadialWeightedBin              )
       velocityDispersionRadialWeightedBin    =mpiSelf%sum(velocityDispersionRadialWeightedBin    )
       velocityDispersionTangentialWeightedBin=mpiSelf%sum(velocityDispersionTangentialWeightedBin)
#endif
       ! Compute averages and dispersions.
       where (pairCountBin > 0_c_size_t)
          velocityRadialBin                      =+     velocityRadialBin                        &
               &                                  /dble(pairCountBin             )
          velocityDispersionRadialBin            =+     velocityDispersionRadialBin              &
               &                                  /dble(pairCountBin             )
          velocityDispersionTangentialBin        =+     velocityDispersionTangentialBin          &
               &                                  /dble(pairCountBin             )
       end where
       where (pairCountWeightedBin > 0.0d0)
          velocityRadialWeightedBin              =+     velocityRadialWeightedBin                &
               &                                  /     pairCountWeightedBin      
          velocityDispersionRadialWeightedBin    =+     velocityDispersionRadialWeightedBin      &
               &                                  /     pairCountWeightedBin      
          velocityDispersionTangentialWeightedBin=+     velocityDispersionTangentialWeightedBin  &
               &                                  /     pairCountWeightedBin      
       end where
       velocityDispersionRadialBin        =velocityDispersionRadialBin        -velocityRadialBin        **2
       velocityDispersionRadialWeightedBin=velocityDispersionRadialWeightedBin-velocityRadialWeightedBin**2
       where (velocityDispersionRadialBin         > 0.0d0)
          velocityDispersionRadialBin        =sqrt(velocityDispersionRadialBin        )
       elsewhere
          velocityDispersionRadialBin        =0.0d0
       end where
       where (velocityDispersionRadialWeightedBin > 0.0d0)
          velocityDispersionRadialWeightedBin=sqrt(velocityDispersionRadialWeightedBin)
       elsewhere
          velocityDispersionRadialWeightedBin=0.0d0
       end where
       velocityDispersionTangentialBin        =sqrt(velocityDispersionTangentialBin        )
       velocityDispersionTangentialWeightedBin=sqrt(velocityDispersionTangentialWeightedBin)
       ! Add Hubble flow if requested.
       if (self%addHubbleFlow) then
          do iSample=1,bootstrapSampleCount
             velocityRadialBin        (:,iSample)=+velocityRadialBin       (:,iSample)                        &
                  &                               +separationCentralBin                                       &
                  &                               *self%cosmologyFunctions_%hubbleParameterEpochal(self%time)
             velocityRadialWeightedBin(:,iSample)=+velocityRadialWeightedBin(:,iSample)                       &
                  &                               +separationCentralBin                                       &
                  &                               *self%cosmologyFunctions_%hubbleParameterEpochal(self%time)
          end do
       end if
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
               & call simulations(jSimulation)%analysis%writeDataset(separationCentralBin                   ,'pairwiseVelocitySeparation'                               )
          call        simulations(jSimulation)%analysis%writeDataset(velocityRadialBin                      ,'pairwiseVelocityRadialMean'                  //char(label))
          call        simulations(jSimulation)%analysis%writeDataset(velocityDispersionRadialBin            ,'pairwiseVelocityRadialDispersion'            //char(label))
          call        simulations(jSimulation)%analysis%writeDataset(velocityDispersionTangentialBin        ,'pairwiseVelocityTangentialDispersion'        //char(label))
          call        simulations(jSimulation)%analysis%writeDataset(velocityRadialWeightedBin              ,'pairwiseVelocityRadialWeightedMean'          //char(label))
          call        simulations(jSimulation)%analysis%writeDataset(velocityDispersionRadialWeightedBin    ,'pairwiseVelocityRadialWeightedDispersion'    //char(label))
          call        simulations(jSimulation)%analysis%writeDataset(velocityDispersionTangentialWeightedBin,'pairwiseVelocityTangentialWeightedDispersion'//char(label))
          !$ call hdf5Access%unset()
#ifdef USEMPI
       end if
#endif
    end do
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call Galacticus_Display_Unindent('done',verbosityStandard)
#ifdef USEMPI
    end if
#endif
    return
  end subroutine pairwiseVelocityStatisticsOperate

