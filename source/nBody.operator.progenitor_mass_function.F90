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
Implements an N-body data operator which computes progenitor mass functions.
!!}

  use            :: Cosmology_Parameters, only : cosmologyParametersClass
  use, intrinsic :: ISO_C_Binding       , only : c_size_t

  !![
  <nbodyOperator name="nbodyOperatorProgenitorMassFunction">
   <description>An N-body data operator which computes progenitor mass functions.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorProgenitorMassFunction
     !!{
     An N-body data operator which computes progenitor mass functions.
     !!}
     private
     class           (cosmologyParametersClass), pointer                   :: cosmologyParameters_       => null()
     integer         (c_size_t                ), allocatable, dimension(:) :: snapshotsProgenitors
     double precision                                                      :: massParentMinimum                   , massParentMaximum                , &
          &                                                                   massRatioProgenitorMinimum          , massRatioProgenitorMaximum
     integer         (c_size_t                )                            :: massParentCountPerDecade            , massRatioProgenitorCountPerDecade, &
          &                                                                   snapshotParents
     type            (varying_string          )                            :: simulationReference                 , simulationURL                    , &
          &                                                                   description
   contains
     final     ::            progenitorMassFunctionDestructor
     procedure :: operate => progenitorMassFunctionOperate
  end type nbodyOperatorProgenitorMassFunction

  interface nbodyOperatorProgenitorMassFunction
     !!{
     Constructors for the \refClass{nbodyOperatorProgenitorMassFunction} N-body operator class.
     !!}
     module procedure progenitorMassFunctionConstructorParameters
     module procedure progenitorMassFunctionConstructorInternal
  end interface nbodyOperatorProgenitorMassFunction

contains

  function progenitorMassFunctionConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorProgenitorMassFunction} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorProgenitorMassFunction)                              :: self
    type            (inputParameters                    ), intent(inout)               :: parameters
    class           (cosmologyParametersClass           ), pointer                     :: cosmologyParameters_
    integer         (c_size_t                           ), allocatable  , dimension(:) :: snapshotsProgenitors
    double precision                                                                   :: massParentMinimum         , massParentMaximum                , &
         &                                                                                massRatioProgenitorMinimum, massRatioProgenitorMaximum
    integer         (c_size_t                           )                              :: massParentCountPerDecade  , massRatioProgenitorCountPerDecade, &
         &                                                                                snapshotParents
    type            (varying_string                     )                              :: simulationReference       , simulationURL                    , &
         &                                                                                description

    !![
    <inputParameter>
      <name>massParentMinimum</name>
      <source>parameters</source>
      <description>The minimum parent mass to consider.</description>
    </inputParameter>
    <inputParameter>
      <name>massParentMaximum</name>
      <source>parameters</source>
      <description>The maximum parent mass to consider.</description>
    </inputParameter>
    <inputParameter>
      <name>massParentCountPerDecade</name>
      <source>parameters</source>
      <description>The number of bins per decade of parent mass.</description>
    </inputParameter>
    <inputParameter>
      <name>massRatioProgenitorMinimum</name>
      <source>parameters</source>
      <description>The minimum mass ratio to consider.</description>
    </inputParameter>
    <inputParameter>
      <name>massRatioProgenitorMaximum</name>
      <source>parameters</source>
      <description>The maximum mass ratio to consider.</description>
    </inputParameter>
    <inputParameter>
      <name>massRatioProgenitorCountPerDecade</name>
      <source>parameters</source>
      <description>The number of bins per decade of mass ratio.</description>
    </inputParameter>
    <inputParameter>
      <name>snapshotParents</name>
      <source>parameters</source>
      <description>The snapshot at which to select parent halos.</description>
    </inputParameter>
    !!]
    allocate(snapshotsProgenitors(parameters%count('snapshotsProgenitors')))
    !![
    <inputParameter>
      <name>snapshotsProgenitors</name>
      <source>parameters</source>
      <description>The snapshots at which to select progenitor halos.</description>
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
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=nbodyOperatorProgenitorMassFunction(massParentMinimum,massParentMaximum,massParentCountPerDecade,massRatioProgenitorMinimum,massRatioProgenitorMaximum,massRatioProgenitorCountPerDecade,snapshotParents,snapshotsProgenitors,description,simulationReference,simulationURL,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function progenitorMassFunctionConstructorParameters

  function progenitorMassFunctionConstructorInternal(massParentMinimum,massParentMaximum,massParentCountPerDecade,massRatioProgenitorMinimum,massRatioProgenitorMaximum,massRatioProgenitorCountPerDecade,snapshotParents,snapshotsProgenitors,description,simulationReference,simulationURL,cosmologyParameters_) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorProgenitorMassFunction} N-body operator class.
    !!}
    implicit none
    type            (nbodyOperatorProgenitorMassFunction)                              :: self
    double precision                                     , intent(in   )               :: massParentMinimum         , massParentMaximum                , &
         &                                                                                massRatioProgenitorMinimum, massRatioProgenitorMaximum
    integer         (c_size_t                           ), intent(in   )               :: massParentCountPerDecade  , massRatioProgenitorCountPerDecade, &
         &                                                                                snapshotParents
    type            (varying_string                     ), intent(in   )               :: simulationReference       , simulationURL                    , &
         &                                                                                description
    integer         (c_size_t                           ), intent(in   ), dimension(:) :: snapshotsProgenitors
    class           (cosmologyParametersClass           ), intent(in   ), target       :: cosmologyParameters_
    !![
    <constructorAssign variables="massParentMinimum, massParentMaximum, massParentCountPerDecade, massRatioProgenitorMinimum, massRatioProgenitorMaximum, massRatioProgenitorCountPerDecade, snapshotParents, snapshotsProgenitors, description, simulationReference, simulationURL, *cosmologyParameters_"/>
    !!]

    return
  end function progenitorMassFunctionConstructorInternal
  
  subroutine progenitorMassFunctionDestructor(self)
    !!{
    Destructor for the \refClass{nbodyOperatorProgenitorMassFunction} N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorProgenitorMassFunction), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine progenitorMassFunctionDestructor

  subroutine progenitorMassFunctionOperate(self,simulations)
    !!{
    Compute mass functions of particles.
    !!}
    use    :: Arrays_Search     , only : searchArray
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
    !$ use :: OMP_Lib           , only : OMP_Get_Thread_Num
    use    :: Sorting           , only : sort
    implicit none
    class           (nbodyOperatorProgenitorMassFunction), intent(inout)                   :: self
    type            (nBodyData                          ), intent(inout), dimension(:    ) :: simulations
    double precision                                     , allocatable  , dimension(:    ) :: massRatioProgenitorBin , massParentBin           , &
         &                                                                                    redshiftProgenitor
    double precision                                     , pointer      , dimension(:    ) :: mass                   , expansionFactor
    integer         (c_size_t                           ), pointer      , dimension(:    ) :: treeID                 , snapshotID
    double precision                                     , allocatable  , dimension(:,:,:) :: progenitorMassFunction
    integer         (c_size_t                           ), allocatable  , dimension(:    ) :: treeIDParents          , indexParents            , &
         &                                                                                    countParentBin
    integer         (c_size_t                           ), allocatable  , dimension(:,:,:) :: countBin
    integer         (c_size_t                           )                                  :: massParentCount        , massRatioProgenitorCount, &
         &                                                                                    iSimulation            , countParents            , &
         &                                                                                    i                      , j                       , &
         &                                                                                    k                      , l                       , &
         &                                                                                    m                      , redshiftProgenitorCount
    double precision                                                                       :: binParentWidthInverse  , binRatioWidthInverse    , &
         &                                                                                    redshiftParent
    logical                                                                                :: snapshotMatched
    type            (hdf5Object                         )                                  :: cosmologyGroup         , simulationGroup

#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayIndent('compute progenitor mass function',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    ! Construct bins in parent mass and in mass ratio.
    massParentCount         =int(log10(self%massParentMaximum         /self%massParentMinimum         )*dble(self%massParentCountPerDecade         ),kind=c_size_t)
    massRatioProgenitorCount=int(log10(self%massRatioProgenitorMaximum/self%massRatioProgenitorMinimum)*dble(self%massRatioProgenitorCountPerDecade),kind=c_size_t)
    redshiftProgenitorCount =size(self%snapshotsProgenitors)
    allocate(massParentBin         (massParentCount                                                 ))
    allocate(massRatioProgenitorBin(                massRatioProgenitorCount                        ))
    allocate(redshiftProgenitor    (                                         redshiftProgenitorCount))
    allocate(countParentBin        (massParentCount                                                 ))
    allocate(countBin              (massParentCount,massRatioProgenitorCount,redshiftProgenitorCount))
    allocate(progenitorMassFunction(massParentCount,massRatioProgenitorCount,redshiftProgenitorCount))
    massParentBin         =Make_Range(self%massParentMinimum         ,self%massParentMaximum         ,int(massParentCount         ),rangeTypeLogarithmic,rangeBinned=.true.)
    massRatioProgenitorBin=Make_Range(self%massRatioProgenitorMinimum,self%massRatioProgenitorMaximum,int(massRatioProgenitorCount),rangeTypeLogarithmic,rangeBinned=.true.)
    binParentWidthInverse =1.0d0/log10(massParentBin         (2)/massParentBin         (1))
    binRatioWidthInverse  =1.0d0/log10(massRatioProgenitorBin(2)/massRatioProgenitorBin(1))
    ! Iterate over simulations.
    do iSimulation=1_c_size_t,size(simulations)
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
       call displayMessage(var_str('simulation "')//simulations(iSimulation)%label//'"',verbosityLevelStandard)
#ifdef USEMPI
       end if
#endif
       ! Get the mass data.
       if (simulations(iSimulation)%propertiesReal   %exists('massVirial'     )) then
          mass            => simulations(iSimulation)%propertiesReal   %value('massVirial'     )
       else
          allocate(mass           (0))
         call Error_Report('halo virial masses are required, but are not available in the simulation'//{introspection:location})
       end if
       ! Get the expansion factor data.
       if (simulations(iSimulation)%propertiesReal   %exists('expansionFactor')) then
          expansionFactor => simulations(iSimulation)%propertiesReal   %value('expansionFactor')
       else
          allocate(expansionFactor(0))
          call Error_Report('expansion factors are required, but are not available in the simulation'//{introspection:location})
       end if
       ! Get the tree ID data.
       if (simulations(iSimulation)%propertiesInteger%exists('hostedRootID'   )) then
          treeID          => simulations(iSimulation)%propertiesInteger%value('hostedRootID'   )
       else
          allocate(treeID         (0))
          call Error_Report('hosted root IDs are required, but are not available in the simulation'  //{introspection:location})
       end if
       ! Get the snapshot ID data.
       if (simulations(iSimulation)%propertiesInteger%exists('snapshotID'     )) then
          snapshotID      => simulations(iSimulation)%propertiesInteger%value('snapshotID'     )
       else
          allocate(snapshotID     (0))
          call Error_Report('snapshot IDs are required, but are not available in the simulation'     //{introspection:location})
       end if
       ! Build a look-up index for parent halo tree IDs.
       countParents      = 0_c_size_t
       redshiftParent    =-1.0d0
       redshiftProgenitor=-1.0d0
       do i=1_c_size_t,size(snapshotID,kind=c_size_t)
          if     (                                         &
               &   snapshotID(i) == self%snapshotParents   &
               &  .and.                                    &
               &   mass      (i) >= self%massParentMinimum &
               &  .and.                                    &
               &   mass      (i) <= self%massParentMaximum &
               & ) countParents=+countParents              &
               &                +1_c_size_t
          if    (snapshotID(i) == self%snapshotParents         .and. redshiftParent        < 0.0d0) redshiftParent       =1.0d0/expansionFactor(i)-1.0d0
          do j=1,redshiftProgenitorCount
             if (snapshotID(i) == self%snapshotsProgenitors(j) .and. redshiftProgenitor(j) < 0.0d0) redshiftProgenitor(j)=1.0d0/expansionFactor(i)-1.0d0
          end do
       end do
       allocate(treeIDParents(countParents))
       allocate(indexParents (countParents))
       countParents  =0_c_size_t
       countParentBin=0_c_size_t
       do i=1_c_size_t,size(snapshotID,kind=c_size_t)
          if     (                                         &
               &   snapshotID(i) == self%snapshotParents   &
               &  .and.                                    &
               &   mass      (i) >= self%massParentMinimum &
               &  .and.                                    &
               &   mass      (i) <= self%massParentMaximum &
               & ) then
             countParents=+countParents &
                  &       +1_c_size_t
             treeIDParents(countParents)=treeID(i)
             indexParents (countParents)=       i
             j=int(log10(mass(i)/self%massParentMinimum)*binParentWidthInverse)+1
             if (j >= 1 .and. j <= massParentCount) countParentBin(j)=countParentBin(j)+1_c_size_t
          end if
       end do
       call sort(treeIDParents,indexParents)
       ! Accumulate counts.
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call displayCounter(0,.true.)
#ifdef USEMPI
       end if
#endif
       countBin=0_c_size_t
       !$omp parallel do private(j,k,l,m,snapshotMatched) reduction(+:countBin) schedule(dynamic)
       do i=1_c_size_t,size(mass,kind=c_size_t)
#ifdef USEMPI
          ! If running under MPI with N processes, process only every Nth particle.
          if (mod(i,mpiSelf%count()) /= mpiSelf%rank()) cycle
#endif
          ! Skip halos not at a requested snapshot.
          snapshotMatched=.false.
          do l=1,redshiftProgenitorCount
             if (self%snapshotsProgenitors(l) == snapshotID(i)) then
                snapshotMatched=.true.
                exit
             end if
          end do
          if (.not.snapshotMatched) cycle
          ! Identify the parent halo, skipping this particle if parent is unacceptable.
          m=searchArray(treeIDParents,treeID(i))
          if (m                <  1_c_size_t .or. m > countParents) cycle
          if (treeIDParents(m) /= treeID(i)                       ) cycle
          ! Accumulate particles into bins.
          j=int(log10(        mass(indexParents(m))/self%massParentMinimum         )*binParentWidthInverse)+1
          if (j < 1 .or. j > massParentCount         ) cycle
          k=int(log10(mass(i)/mass(indexParents(m))/self%massRatioProgenitorMinimum)*binRatioWidthInverse )+1
          if (k < 1 .or. k > massRatioProgenitorCount) cycle
          countBin(j,k,l)=+countBin  (j,k,l) &
               &          +1_c_size_t
          ! Update progress.
          !$ if (OMP_Get_Thread_Num() == 0) then
#ifdef USEMPI
          if (mpiSelf%isMaster()) then
#endif
             call displayCounter(                                      &
                  &                          int(                                  &
                  &                              +100.0d0                          &
                  &                              *float(i                       )  &
                  &                              /float(size(mass,kind=c_size_t))  &
                  &                             )                                , &
                  &                          .false.                               &
                  &                         )
#ifdef USEMPI
          end if
#endif
          !$ end if
       end do
       !$omp end parallel do
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call displayCounterClear()
#ifdef USEMPI
       end if
#endif
#ifdef USEMPI
       ! Reduce across MPI processes.
       countBin=mpiSelf%sum(countBin)
#endif
       ! Compute progenitor mass function.
       progenitorMassFunction=dble(countBin)*binRatioWidthInverse/log(10.0d0)
       do k=1,massRatioProgenitorCount
          do l=1,redshiftProgenitorCount
             where (countParentBin > 0_c_size_t)
                progenitorMassFunction(:,k,l)=progenitorMassFunction(:,k,l)/dble(countParentBin)
             elsewhere
                progenitorMassFunction(:,k,l)=0.0d0
             end where
          end do
       end do
       do j=1,massParentCount
          do l=1,redshiftProgenitorCount
             progenitorMassFunction(j,:,l)=progenitorMassFunction(j,:,l)*massRatioProgenitorBin
          end do
       end do
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call simulations(iSimulation)%analysis%writeDataset  (massParentBin            ,'massParent'            )
          call simulations(iSimulation)%analysis%writeDataset  (massRatioProgenitorBin   ,'massRatioProgenitor'   )
          call simulations(iSimulation)%analysis%writeDataset  (redshiftProgenitor       ,'redshiftProgenitor'    )
          call simulations(iSimulation)%analysis%writeAttribute(redshiftParent           ,'redshiftParent'        )
          call simulations(iSimulation)%analysis%writeDataset  (countBin                 ,'count'                 )
          call simulations(iSimulation)%analysis%writeDataset  (progenitorMassFunction   ,'progenitorMassFunction')
          call simulations(iSimulation)%analysis%writeAttribute(self%description         ,"description"           )
          call simulations(iSimulation)%analysis%writeAttribute(Formatted_Date_and_Time(),"timestamp"             )
          cosmologyGroup=simulations(iSimulation)%analysis%openGroup('cosmology')
          call cosmologyGroup%writeAttribute(self%cosmologyParameters_%OmegaMatter    (),'OmegaMatter'    )
          call cosmologyGroup%writeAttribute(self%cosmologyParameters_%OmegaDarkEnergy(),'OmegaDarkEnergy')
          call cosmologyGroup%writeAttribute(self%cosmologyParameters_%HubbleConstant (),'HubbleConstant' )
          simulationGroup=simulations(iSimulation)%analysis%openGroup('simulation')
          call simulationGroup%writeAttribute(self%simulationReference,'reference')
          call simulationGroup%writeAttribute(self%simulationURL      ,'URL'      )
          if (simulations(iSimulation)%attributesReal%exists('massParticle')) &
               & call simulationGroup%writeAttribute(simulations(iSimulation)%attributesReal%value('massParticle'),'massParticle')
          if (simulations(iSimulation)%attributesReal%exists('boxSize'     )) &
               & call simulationGroup%writeAttribute(simulations(iSimulation)%attributesReal%value('boxSize'     ),'boxSize'     )
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
  end subroutine progenitorMassFunctionOperate

