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
Implements an N-body data operator which computes pair counts in bins of separation.
!!}

  use, intrinsic :: ISO_C_Binding           , only : c_size_t
  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  !![
  <nbodyOperator name="nbodyOperatorPairCounts">
   <description>An N-body data operator which computes pair counts in bins of separation.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorPairCounts
     !!{
     An N-body data operator which computes pair counts in bins of separation.
     !!}
     private
     double precision                                      :: separationMinimum               , separationMaximum   , &
          &                                                   bootstrapSampleRate
     integer         (c_size_t                  )          :: separationCount                 , bootstrapSampleCount
     logical                                               :: includeUnbootstrapped           , crossCount
     class           (randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()
   contains
     final     ::            pairCountsDestructor
     procedure :: operate => pairCountsOperate
  end type nbodyOperatorPairCounts

  interface nbodyOperatorPairCounts
     !!{
     Constructors for the \refClass{nbodyOperatorPairCounts} N-body operator class.
     !!}
     module procedure pairCountsConstructorParameters
     module procedure pairCountsConstructorInternal
  end interface nbodyOperatorPairCounts

contains

  function pairCountsConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorPairCounts} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorPairCounts   )                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass), pointer       :: randomNumberGenerator_
    double precision                                            :: separationMinimum     , separationMaximum   , &
         &                                                         bootstrapSampleRate
    integer         (c_size_t                  )                :: separationCount       , bootstrapSampleCount
    logical                                                     :: includeUnbootstrapped , crossCount

    !![
    <inputParameter>
      <name>crossCount</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, compute cross-simulation pair counts between the first and all simulations. Otherwise, compute pair counts within each simulation.</description>
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
      <name>separationMinimum</name>
      <source>parameters</source>
      <description>The minimum separation to consider for pair counts.</description>
    </inputParameter>
    <inputParameter>
      <name>separationMaximum</name>
      <source>parameters</source>
      <description>The maximum separation to consider for pair counts.</description>
    </inputParameter>
    <inputParameter>
      <name>separationCount</name>
      <source>parameters</source>
      <description>The number of bins in separation for pair counts.</description>
    </inputParameter>
    <inputParameter>
      <name>includeUnbootstrapped</name>
      <source>parameters</source>
      <description>If true, include results for the unbootstrapped (i.e. original) sample.</description>
      <defaultValue>.true.</defaultValue>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=nbodyOperatorPairCounts(separationMinimum,separationMaximum,separationCount,crossCount,includeUnbootstrapped,bootstrapSampleCount,bootstrapSampleRate,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function pairCountsConstructorParameters

  function pairCountsConstructorInternal(separationMinimum,separationMaximum,separationCount,crossCount,includeUnbootstrapped,bootstrapSampleCount,bootstrapSampleRate,randomNumberGenerator_) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorPairCounts} N-body operator class.
    !!}
    implicit none
    type            (nbodyOperatorPairCounts)                           :: self
    double precision                            , intent(in   )         :: separationMinimum     , separationMaximum   , &
         &                                                                 bootstrapSampleRate
    integer         (c_size_t                  ), intent(in   )         :: separationCount       , bootstrapSampleCount
    logical                                     , intent(in   )         :: includeUnbootstrapped , crossCount
    class           (randomNumberGeneratorClass), intent(in   ), target :: randomNumberGenerator_
    !![
    <constructorAssign variables="separationMinimum, separationMaximum, separationCount, crossCount, includeUnbootstrapped, bootstrapSampleCount, bootstrapSampleRate, *randomNumberGenerator_"/>
    !!]

    return
  end function pairCountsConstructorInternal

  subroutine pairCountsDestructor(self)
    !!{
    Destructor for the \refClass{nbodyOperatorPairCounts} N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorPairCounts), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine pairCountsDestructor
  
  subroutine pairCountsOperate(self,simulations)
    !!{
    Compute pair counts of the particles in bins of separation.
    !!}
    use    :: Arrays_Search     , only : searchArray
    use    :: Display           , only : displayCounter    , displayCounterClear   , displayIndent, displayMessage, &
          &                              displayUnindent   , verbosityLevelStandard
    use    :: HDF5_Access       , only : hdf5Access
    use    :: ISO_Varying_String, only : var_str
#ifdef USEMPI
    use    :: MPI_Utilities     , only : mpiSelf
#endif
    use    :: Nearest_Neighbors , only : nearestNeighbors
    use    :: Numerical_Ranges  , only : Make_Range        , rangeTypeLogarithmic
    !$ use :: OMP_Lib           , only : OMP_Get_Thread_Num
    implicit none
    class           (nbodyOperatorPairCounts   ), intent(inout)                 :: self
    type            (nBodyData                 ), intent(inout), dimension(:  ) :: simulations
    double precision                            , parameter                     :: toleranceZero              =0.0d0
    integer                                     , allocatable  , dimension(:  ) :: neighborIndex
    double precision                            , allocatable  , dimension(:  ) :: neighborDistanceSquared          , separationCentralBin       , &
         &                                                                         separationSquaredMinimumBin      , separationSquaredMaximumBin
    integer         (c_size_t                  ), allocatable  , dimension(:,:) :: weight1                          , weight2                    , &
         &                                                                         weightPair
    integer         (c_size_t                  ), allocatable  , dimension(:,:) :: pairCountBin
    double precision                            , pointer      , dimension(:,:) :: position1                        , position2
    integer         (c_size_t                  )                                :: i                                , j                          , &
         &                                                                         iSample                          , bootstrapSampleCount       , &
         &                                                                         iSimulation                      , jSimulation                , &
         &                                                                         jStart                           , jEnd                       , &
         &                                                                         pairCountTotal
    type            (nearestNeighbors          )                                :: neighborFinder
    integer                                                                     :: neighborCount
    type            (varying_string            )                                :: label
    
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayIndent('compute pair counts',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    ! Construct bins of separation.
    bootstrapSampleCount=self%bootstrapSampleCount
    if (self%includeUnbootstrapped) bootstrapSampleCount=bootstrapSampleCount+1
    allocate(separationCentralBin       (self%separationCount                     ))
    allocate(separationSquaredMinimumBin(self%separationCount                     ))
    allocate(separationSquaredMaximumBin(self%separationCount                     ))
    allocate(pairCountBin               (self%separationCount,bootstrapSampleCount))
    separationCentralBin       =Make_Range(self%separationMinimum,self%separationMaximum,int(self%separationCount),rangeTypeLogarithmic,rangeBinned=.true.)
    separationSquaredMinimumBin=(separationCentralBin/sqrt(separationCentralBin(2)/separationCentralBin(1)))**2
    separationSquaredMaximumBin=(separationCentralBin*sqrt(separationCentralBin(2)/separationCentralBin(1)))**2
    ! Iterate over simulations.
    do iSimulation=1_c_size_t,size(simulations)
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
       call displayMessage(var_str('simulation "')//simulations(iSimulation)%label//'"',verbosityLevelStandard)
#ifdef USEMPI
       end if
#endif
       ! Allocate workspace.
       if (self%crossCount) then
          position1 => simulations(1          )%propertiesRealRank1%value('position')
       else
          position1 => simulations(iSimulation)%propertiesRealRank1%value('position')
       end if
       position2    => simulations(iSimulation)%propertiesRealRank1%value('position')
       allocate(weight1(bootstrapSampleCount,size(position1,dim=2)))
       allocate(weight2(bootstrapSampleCount,size(position2,dim=2)))
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
       ! Accumulate pairs
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call displayCounter(0,.true.)
#ifdef USEMPI
       end if
#endif
       pairCountBin=0_c_size_t
       !$omp parallel private(j,jStart,jEnd,weightPair,neighborCount,neighborIndex,neighborDistanceSquared,neighborFinder) reduction(+:pairCountBin)
       ! Allocate workspace.
       allocate(weightPair             (bootstrapSampleCount,size(position2,dim=2)))
       allocate(neighborIndex          (                     size(position2,dim=2)))
       allocate(neighborDistanceSquared(                     size(position2,dim=2)))
       ! Construct the nearest neighbor finder object.
       neighborFinder=nearestNeighbors(transpose(position2))
       ! Iterate over particles.
       !$omp do schedule(dynamic)
       do i=1_c_size_t,size(position1,dim=2,kind=c_size_t)
#ifdef USEMPI
          ! If running under MPI with N processes, process only every Nth particle.
          if (mod(i,mpiSelf%count()) /= mpiSelf%rank()) cycle
#endif
          ! Locate particles nearby.
          call neighborFinder%searchFixedRadius(position1(:,i),self%separationMaximum,toleranceZero,neighborCount,neighborIndex,neighborDistanceSquared)
          ! Find weights of paired particles.
          do j=1_c_size_t,neighborCount
             weightPair(:,j)=weight2(:,neighborIndex(j))
          end do
          ! Accumulate particles into bins.
          jEnd=0
          do j=1,self%separationCount
             if (neighborDistanceSquared(1) > separationSquaredMaximumBin(j) .or. neighborDistanceSquared(neighborCount) < separationSquaredMinimumBin(j)) cycle
             jStart=jEnd
             jEnd  =searchArray(neighborDistanceSquared(1:neighborCount),separationSquaredMaximumBin(j))
             if (jStart == jEnd) cycle
             pairCountBin(j,:)=+    pairCountBin(j        ,:    )        &
                  &            +    weight1     (:,i            )        &
                  &            *sum(weightPair  (:,jStart+1:jEnd),dim=2)
          end do
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
       !$omp end parallel
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call displayCounterClear()
#ifdef USEMPI
       end if
#endif
       deallocate(weight1)
       deallocate(weight2)
#ifdef USEMPI
       ! Reduce across MPI processes.
       pairCountBin=mpiSelf%sum(pairCountBin)
#endif
       if (self%crossCount) then
          jSimulation   = 1
          pairCountTotal=+size(                         position1,dim=2,kind=c_size_t) &
               &         *size(                         position2,dim=2,kind=c_size_t)
          label         = ":"//simulations(          1)%label// &
               &          ":"//simulations(iSimulation)%label
       else
          jSimulation   = iSimulation
          pairCountTotal=+size(                         position2,dim=2,kind=c_size_t) &
               &         *size(                         position2,dim=2,kind=c_size_t)
          label         = ":"//simulations(iSimulation)%label
       end if
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          !$ call hdf5Access%  set()
          if (.not.self%crossCount .or. iSimulation == 1 ) &
               & call simulations(jSimulation)%analysis%writeDataset  (separationCentralBin,'pairCountSeparation'             )
          call        simulations(jSimulation)%analysis%writeDataset  (pairCountBin        ,'pairCountCount'     //char(label))
          call        simulations(jSimulation)%analysis%writeAttribute(pairCountTotal      ,'pairCountTotal'     //char(label))
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
  end subroutine pairCountsOperate

