!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
Implements an N-body data operator which computes correlations between mass functions at different redshifts.
!!}

  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
  use, intrinsic :: ISO_C_Binding           , only : c_size_t

  !![
  <nbodyOperator name="nbodyOperatorMassFunctionCorrelation">
   <description>An N-body data operator which computes correlations between mass functions at different redshifts.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorMassFunctionCorrelation
     !!{
     An N-body data operator which computes correlations between mass functions at different redshifts.
     !!}
     private
     class           (randomNumberGeneratorClass), pointer :: randomNumberGenerator_   => null()
     double precision                                      :: massMinimum                       , massMaximum
     integer         (c_size_t                  )          :: massCountPerDecade                , countBootstraps
     type            (varying_string            )          :: simulationReference               , simulationURL  , &
          &                                                   description
     logical                                               :: missingRootHalosAreFatal
   contains
     final     ::            massFunctionCorrelationDestructor
     procedure :: operate => massFunctionCorrelationOperate
  end type nbodyOperatorMassFunctionCorrelation

  interface nbodyOperatorMassFunctionCorrelation
     !!{
     Constructors for the \refClass{nbodyOperatorMassFunctionCorrelation} N-body operator class.
     !!}
     module procedure massFunctionCorrelationConstructorParameters
     module procedure massFunctionCorrelationConstructorInternal
  end interface nbodyOperatorMassFunctionCorrelation

  type :: simulationData
     !!{
     Type used to store simulation data for correlation analysis.
     !!}
     double precision                                  :: boxSize
     double precision          , pointer, dimension(:) :: mass
     integer         (c_size_t), pointer, dimension(:) :: particleIDs, hostedRootIDs
  end type simulationData
  
contains

  function massFunctionCorrelationConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorMassFunctionCorrelation} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorMassFunctionCorrelation)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass          ), pointer       :: randomNumberGenerator_
    double precision                                                      :: massMinimum             , massMaximum
    integer         (c_size_t                            )                :: massCountPerDecade      , countBootstraps
    logical                                                               :: missingRootHalosAreFatal
    type            (varying_string                      )                :: simulationReference     , simulationURL, &
         &                                                                   description

    !![
    <inputParameter>
      <name>missingRootHalosAreFatal</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, if a hosted root halo is not found then a fatal error occurs. Otherwise, such cases are ignored and will not contribute to the halo mass function.</description>
    </inputParameter>
    <inputParameter>
      <name>massMinimum</name>
      <source>parameters</source>
      <description>The minimum mass to consider counts.</description>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <source>parameters</source>
      <description>The maximum mass to consider.</description>
    </inputParameter>
    <inputParameter>
      <name>massCountPerDecade</name>
      <source>parameters</source>
      <description>The number of bins per decade of mass.</description>
    </inputParameter>
    <inputParameter>
      <name>countBootstraps</name>
      <source>parameters</source>
      <description>The number of bootstrap samples to perform.</description>
    </inputParameter>
    <inputParameter>
      <name>description</name>
      <source>parameters</source>
      <description>A description of this mass function.</description>
    </inputParameter>
    <inputParameter>
      <name>simulationReference</name>
      <source>parameters</source>
      <description>A reference for the simulation.</description>
    </inputParameter>
    <inputParameter>
      <name>simulationURL</name>
      <source>parameters</source>
      <description>A URL for the simulation.</description>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=nbodyOperatorMassFunctionCorrelation(massMinimum,massMaximum,massCountPerDecade,countBootstraps,missingRootHalosAreFatal,description,simulationReference,simulationURL,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function massFunctionCorrelationConstructorParameters

  function massFunctionCorrelationConstructorInternal(massMinimum,massMaximum,massCountPerDecade,countBootstraps,missingRootHalosAreFatal,description,simulationReference,simulationURL,randomNumberGenerator_) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorMassFunctionCorrelation} N-body operator class.
    !!}
    implicit none
    type            (nbodyOperatorMassFunctionCorrelation)                        :: self
    double precision                                      , intent(in   )         :: massMinimum             , massMaximum
    integer         (c_size_t                            ), intent(in   )         :: massCountPerDecade      , countBootstraps
    logical                                               , intent(in   )         :: missingRootHalosAreFatal
    type            (varying_string                      ), intent(in   )         :: simulationReference     , simulationURL, &
&                                                                                    description
    class           (randomNumberGeneratorClass          ), intent(in   ), target :: randomNumberGenerator_
    !![
    <constructorAssign variables="massMinimum, massMaximum, massCountPerDecade, countBootstraps, missingRootHalosAreFatal, description, simulationReference, simulationURL, *randomNumberGenerator_"/>
    !!]

    return
  end function massFunctionCorrelationConstructorInternal
  
  subroutine massFunctionCorrelationDestructor(self)
    !!{
    Destructor for the \refClass{nbodyOperatorMassFunctionCorrelation} N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorMassFunctionCorrelation), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine massFunctionCorrelationDestructor

  subroutine massFunctionCorrelationOperate(self,simulations)
    !!{
    Compute mass functions of particles.
    !!}
    use    :: Arrays_Search     , only : searchIndexed
    use    :: Dates_and_Times   , only : Formatted_Date_and_Time
    use    :: Display           , only : displayCounter         , displayCounterClear   , displayIndent, displayMessage, &
         &                               displayUnindent        , verbosityLevelStandard
    use    :: Error             , only : Error_Report
    use    :: HDF5_Access       , only : hdf5Access
    use    :: IO_HDF5           , only : hdf5Object
    use    :: ISO_Varying_String, only : var_str
#ifdef USEMPI
    use    :: MPI_Utilities     , only : mpiSelf
#endif
    use    :: Numerical_Ranges  , only : Make_Range             , rangeTypeLogarithmic
    use    :: Output_HDF5       , only : outputFile
    use    :: Sorting           , only : sortIndex
    !$ use :: OMP_Lib           , only : OMP_Get_Thread_Num
    implicit none
    class           (nbodyOperatorMassFunctionCorrelation), intent(inout)                 :: self
    type            (nBodyData                           ), intent(inout), dimension(:  ) :: simulations
    double precision                                      , allocatable  , dimension(:  ) :: massFunction        , massBin        , &
         &                                                                                   massFunctionsAverage
    double precision                                      , allocatable  , dimension(:,:) :: massFunctions       , correlation    , &
         &                                                                                   covariance
    integer         (c_size_t                            ), allocatable  , dimension(:  ) :: countBin            , indexID
    integer                                               , allocatable  , dimension(:  ) :: weight
    type            (simulationData                      ), allocatable  , dimension(:  ) :: simulationData_
    integer         (c_size_t                            )                                :: iSimulation     , jSimulation    , massCount     , &
         &                                                                                   i                   , j              , &
         &                                                                                   k       , ii, jj   , iii, jjj         , iBootstrap
    double precision                                                                      :: binWidthInverse, massMinimum, massMaximum
    type            (hdf5Object                          )                                :: simulationGroup

#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayIndent('compute mass function correlations',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    ! Construct bins in mass.
    massCount=int(log10(self%massMaximum/self%massMinimum)*dble(self%massCountPerDecade),kind=c_size_t)
    allocate(massBin             (massCount                                              ))
    allocate(massFunction        (massCount                                              ))
    allocate(countBin            (massCount                                              ))
    allocate(massFunctions       (massCount*size(simulations),self%countBootstraps       ))
    allocate(massFunctionsAverage(massCount*size(simulations)                            ))
    allocate(covariance          (massCount*size(simulations),massCount*size(simulations)))
    allocate(correlation         (massCount*size(simulations),massCount*size(simulations)))
    massBin        =Make_Range(self%massMinimum,self%massMaximum,int(massCount),rangeTypeLogarithmic,rangeBinned=.true.)
    binWidthInverse=1.0d0/log10(massBin(2)/massBin(1))
    ! Read all simulations.
    allocate(simulationData_(size(simulations)))
    do iSimulation=1_c_size_t,size(simulations)
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
       call displayMessage(var_str('reading simulation "')//simulations(iSimulation)%label//'"',verbosityLevelStandard)
#ifdef USEMPI
       end if
#endif
       ! Get box size.
       if (simulations(iSimulation)%attributesReal%exists('boxSize')) then
          simulationData_(iSimulation)%boxSize=simulations(iSimulation)%attributesReal%value('boxSize')
       else
          simulationData_(iSimulation)%boxSize=0.0d0
          call Error_Report('box size is required, but is not available in the simulation'//{introspection:location})
       end if
       ! Get the mass data.
       if (simulations(iSimulation)%propertiesReal%exists('massVirial')) then
          simulationData_(iSimulation)%mass => simulations(iSimulation)%propertiesReal%value('massVirial')
       else
          call Error_Report('halo virial masses are required, but are not available in the simulation'//{introspection:location})
       end if
       ! Get the ID data.
       if (simulations(iSimulation)%propertiesInteger%exists('particleID')) then
          simulationData_(iSimulation)%particleIDs => simulations(iSimulation)%propertiesInteger%value('particleID')
       else
          call Error_Report('particle IDs are required, but are not available in the simulation'//{introspection:location})
       end if
       if (simulations(iSimulation)%propertiesInteger%exists('hostedRootID')) then
          simulationData_(iSimulation)%hostedRootIDs => simulations(iSimulation)%propertiesInteger%value('hostedRootID')
       else
          call Error_Report('hosted root IDs are required, but are not available in the simulation'//{introspection:location})
       end if
    end do
    ! Build a sort index into hosted IDs of the first simulation.
    indexID=sortIndex(simulationData_(1)%hostedRootIDs)
    ! Iterate over bootstrap resamplings.
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayCounter(0,.true.)
#ifdef USEMPI
    end if
#endif
    allocate(weight(size(simulationData_(1)%mass)))
    do iBootstrap=1_c_size_t,self%countBootstraps
       ! Generate bootstrap weights.
       do i=1_c_size_t,size(simulationData_(1)%mass,kind=c_size_t)
          j        =searchIndexed(simulationData_(1)%hostedRootIDs,indexID,simulationData_(1)%hostedRootIDs(i))
          weight(j)=self%randomNumberGenerator_%poissonSample(1.0d0)
       end do
       ! Iterate over simulations, constructing bootstrap-weighted mass functions.
       do iSimulation=1_c_size_t,size(simulations)
          ! Accumulate counts.
          countBin=0_c_size_t
          !$omp parallel do private(j,k) reduction(+:countBin) schedule(dynamic)
          do i=1_c_size_t,size(simulationData_(iSimulation)%mass,kind=c_size_t)
#ifdef USEMPI
             ! If running under MPI with N processes, process only every Nth particle.
             if (mod(i,mpiSelf%count()) /= mpiSelf%rank()) cycle
#endif
             ! Accumulate particles into bins.
             j=floor(log10(simulationData_(iSimulation)%mass(i)/self%massMinimum)*binWidthInverse)+1
             if (j >= 1 .and. j <= massCount) then
                k=searchIndexed(simulationData_(1)%hostedRootIDs,indexID,simulationData_(iSimulation)%hostedRootIDs(i))
                if     (                    &
                     &   k >= 1_c_size_t    &
                     &  .and.               &
                     &   k <= size(indexID) &
                     & ) then                
                   countBin(j)=+countBin(j) &
                        &      +weight  (k)
                else if (self%missingRootHalosAreFatal) then
                   call Error_Report('can not find hosted root'//{introspection:location})
                end if
             end if
          end do
          !$omp end parallel do
#ifdef USEMPI
          ! Reduce across MPI processes.
          countBin=mpiSelf%sum(countBin)
#endif
          ! Compute mass function.
          massFunction=dble(countBin)*binWidthInverse/log(10.0d0)/simulationData_(iSimulation)%boxSize**3
          ! Accumulate into a single array across all simulations.
          massFunctions((iSimulation-1)*massCount+1:iSimulation*massCount,iBootstrap)=massFunction
       end do
       ! Update progress.
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call displayCounter(                                  &
               &              int(                              &
               &                  +100.0d0                      &
               &                  *float(     iBootstrap     )  &
               &                  /float(self%countBootstraps)  &
               &                 )                            , &
               &              .false.                           &
               &             )
#ifdef USEMPI
       end if
#endif
    end do
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayCounterClear()
#ifdef USEMPI
    end if
#endif
    ! Compute the average mass function across bootstraps.
    massFunctionsAverage=sum(massFunctions,dim=2)/dble(self%countBootstraps)
    ! Subtract the average from all mass functions.
    do iBootstrap=1_c_size_t,self%countBootstraps
       massFunctions(:,iBootstrap)=+massFunctions       (:,iBootstrap) &
            &                      -massFunctionsAverage
    end do
    ! Compute the covariances.
    do i=1_c_size_t,size(massFunctionsAverage)
       do j=i,size(massFunctionsAverage)
          covariance(i,j)=+sum(                       &
               &               +massFunctions(i,:)    &
               &               *massFunctions(j,:)    &
               &              )                       &
               &          /dble(self%countBootstraps)
          covariance(j,i)=+covariance(i,j)
       end do
    end do
    ! Convert to correlation.
    do i=1_c_size_t,size(massFunctionsAverage)
       do j=1_c_size_t,size(massFunctionsAverage)
          if (i == j) then
             correlation(i,j)=1.0d0
          else if (covariance(i,j) /= 0.0d0) then
             correlation(i,j)=+      covariance(i,j) &
                  &           /sqrt(                 &
                  &                 +covariance(i,i) &
                  &                 *covariance(j,j) &
                  &                )
          else
             correlation(i,j)=+0.0d0
          end if
       end do
    end do
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       !$ call hdf5Access%set()
       simulationGroup=outputFile%openGroup('correlation')
       call simulationGroup%writeDataset  (massBin                  ,'mass'       )
       call simulationGroup%writeDataset  (correlation              ,'correlation')
       call simulationGroup%writeAttribute(self%description         ,"description")
       call simulationGroup%writeAttribute(Formatted_Date_and_Time(),"timestamp"  )
       call simulationGroup%close()
       !$ call hdf5Access%unset()
#ifdef USEMPI
    end if
#endif
    ! Clean up
    do iSimulation=1_c_size_t,size(simulations)
       nullify(simulationData_(iSimulation)%mass         )
       nullify(simulationData_(iSimulation)%particleIDs  )
       nullify(simulationData_(iSimulation)%hostedRootIDs)
    end do
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayUnindent('done',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    return
  end subroutine massFunctionCorrelationOperate

