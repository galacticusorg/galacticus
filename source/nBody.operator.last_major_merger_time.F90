!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Contains a module which implements an N-body data operator which computes progenitor halo order statistics.
!!}

  use            :: Cosmology_Parameters, only : cosmologyParametersClass
  use, intrinsic :: ISO_C_Binding       , only : c_size_t

  !![
  <nbodyOperator name="nbodyOperatorTimeLastMajorMerger">
   <description>An N-body data operator which computes progenitor halo order statistics.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorTimeLastMajorMerger
     !!{
     An N-body data operator which computes progenitor halo order statistics.
     !!}
     private
     class           (cosmologyParametersClass), pointer                   :: cosmologyParameters_     => null()
     integer         (c_size_t                ), allocatable, dimension(:) :: snapshotsProgenitors
     double precision                                                      :: massParentMinimum                 , massParentMaximum         , &
          &                                                                   redshiftMergerMinimum             , redshiftMergerMaximum     , &
          &                                                                   ratioMajor
     integer         (c_size_t                )                            :: massParentCountPerDecade          , redshiftMergerCountPerUnit, &
          &                                                                   snapshotParents
     type            (varying_string          )                            :: simulationReference               , simulationURL             , &
          &                                                                   description
   contains
     final     ::            timeLastMajorMergerDestructor
     procedure :: operate => timeLastMajorMergerOperate
  end type nbodyOperatorTimeLastMajorMerger

  interface nbodyOperatorTimeLastMajorMerger
     !!{
     Constructors for the ``timeLastMajorMerger'' N-body operator class.
     !!}
     module procedure timeLastMajorMergerConstructorParameters
     module procedure timeLastMajorMergerConstructorInternal
  end interface nbodyOperatorTimeLastMajorMerger

contains

  function timeLastMajorMergerConstructorParameters(parameters) result (self)
    !!{
    Constructor for the ``timeLastMajorMerger'' N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorTimeLastMajorMerger)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (cosmologyParametersClass        ), pointer       :: cosmologyParameters_
    double precision                                                  :: massParentMinimum       , massParentMaximum            , &
         &                                                               redshiftMergerMinimum, redshiftMergerMaximum     , &
         &                                                               ratioMajor
    integer         (c_size_t                        )                :: massParentCountPerDecade, redshiftMergerCountPerUnit, &
         &                                                               snapshotParents
    type            (varying_string                  )                :: simulationReference     , simulationURL                , &
         &                                                               description

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
      <name>redshiftMergerMinimum</name>
      <source>parameters</source>
      <description>The minimum merger redshift to consider.</description>
    </inputParameter>
    <inputParameter>
      <name>redshiftMergerMaximum</name>
      <source>parameters</source>
      <description>The maximum merger redshift to consider.</description>
    </inputParameter>
    <inputParameter>
      <name>redshiftMergerCountPerUnit</name>
      <source>parameters</source>
      <description>The number of bins per unit of merger redshift.</description>
    </inputParameter>
    <inputParameter>
      <name>snapshotParents</name>
      <source>parameters</source>
      <description>The snapshot at which to select parent halos.</description>
    </inputParameter>
    <inputParameter>
      <name>ratioMajor</name>
      <source>parameters</source>
      <description>The mass ratio used to define a major merger.</description>
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
    self=nbodyOperatorTimeLastMajorMerger(massParentMinimum,massParentMaximum,massParentCountPerDecade,redshiftMergerMinimum,redshiftMergerMaximum,redshiftMergerCountPerUnit,snapshotParents,ratioMajor,description,simulationReference,simulationURL,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function timeLastMajorMergerConstructorParameters

  function timeLastMajorMergerConstructorInternal(massParentMinimum,massParentMaximum,massParentCountPerDecade,redshiftMergerMinimum,redshiftMergerMaximum,redshiftMergerCountPerUnit,snapshotParents,ratioMajor,description,simulationReference,simulationURL,cosmologyParameters_) result (self)
    !!{
    Internal constructor for the ``timeLastMajorMerger'' N-body operator class.
    !!}
    implicit none
    type            (nbodyOperatorTimeLastMajorMerger)                        :: self
    double precision                                  , intent(in   )         :: massParentMinimum       , massParentMaximum         , &
         &                                                                       redshiftMergerMinimum   , redshiftMergerMaximum     , &
         &                                                                       ratioMajor
    integer         (c_size_t                        ), intent(in   )         :: massParentCountPerDecade, redshiftMergerCountPerUnit, &
         &                                                                       snapshotParents
    type            (varying_string                  ), intent(in   )         :: simulationReference     , simulationURL             , &
         &                                                                       description
    class           (cosmologyParametersClass        ), intent(in   ), target :: cosmologyParameters_
    !![
    <constructorAssign variables="massParentMinimum, massParentMaximum, massParentCountPerDecade, redshiftMergerMinimum, redshiftMergerMaximum, redshiftMergerCountPerUnit, snapshotParents, ratioMajor, description, simulationReference, simulationURL, *cosmologyParameters_"/>
    !!]

    return
  end function timeLastMajorMergerConstructorInternal
  
  subroutine timeLastMajorMergerDestructor(self)
    !!{
    Destructor for the ``timeLastMajorMerger'' N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorTimeLastMajorMerger), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine timeLastMajorMergerDestructor

  subroutine timeLastMajorMergerOperate(self,simulations)
    !!{
    Compute mass functions of particles.
    !!}
    use    :: Arrays_Search     , only : searchArray            , searchIndexed
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
    use    :: Numerical_Ranges  , only : Make_Range             , rangeTypeLogarithmic  , rangeTypeLinear
    !$ use :: OMP_Lib           , only : OMP_Get_Thread_Num
    use    :: Sorting           , only : sort                   , sortIndex
    implicit none
    class           (nbodyOperatorTimeLastMajorMerger), intent(inout)                 :: self
    type            (nBodyData                       ), intent(inout), dimension(:  ) :: simulations
    double precision                                  , allocatable  , dimension(:  ) :: massParentBin        , redshiftMergerBin
    double precision                                  , pointer      , dimension(:  ) :: mass                 , expansionFactor
    integer         (c_size_t                        ), pointer      , dimension(:  ) :: treeID               , snapshotID             , &
         &                                                                               descendantID         , particleID
    double precision                                  , allocatable  , dimension(:,:) :: redshiftMerger
    double precision                                                 , dimension(2  ) :: massesMerger
    integer         (c_size_t                        ), allocatable  , dimension(:  ) :: treeIDParents        , indexParents           , &
         &                                                                               countParentBin       , indexID
    integer         (c_size_t                        ), allocatable  , dimension(:,:) :: countBin
    integer         (c_size_t                        )                                :: massParentCount      , redshiftMergerCount    , &
         &                                                                               iSimulation          , countParents           , &
         &                                                                               i                    , j                      , &
         &                                                                               k                    , l                      , &
         &                                                                               m                    , n
    double precision                                                                  :: binParentWidthInverse, binRedshiftWidthInverse, &
         &                                                                               redshiftParent       , expansionFactorMerger  , &
         &                                                                               redshiftMerger_
    type            (hdf5Object                      )                                :: cosmologyGroup       , simulationGroup        , &
         &                                                                               resultsGroup

#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayIndent('compute last major merger time statistics',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    ! Construct bins in parent mass and in merger redshift.
    massParentCount    =int(log10(self%massParentMaximum    /self%massParentMinimum    )*dble(self%massParentCountPerDecade  ),kind=c_size_t)
    redshiftMergerCount=int(     (self%redshiftMergerMaximum-self%redshiftMergerMinimum)*dble(self%redshiftMergerCountPerUnit),kind=c_size_t)
    allocate(massParentBin    (massParentCount                    ))
    allocate(redshiftMergerBin(                redshiftMergerCount))
    allocate(countParentBin   (massParentCount                    ))
    allocate(countBin         (massParentCount,redshiftMergerCount))
    allocate(redshiftMerger   (massParentCount,redshiftMergerCount))
    massParentBin          =Make_Range(self%massParentMinimum    ,self%massParentMaximum    ,int(massParentCount    ),rangeTypeLogarithmic,rangeBinned=.true.)
    redshiftMergerBin      =Make_Range(self%redshiftMergerMinimum,self%redshiftMergerMaximum,int(redshiftMergerCount),rangeTypeLinear     ,rangeBinned=.true.)
    binParentWidthInverse  =1.0d0/log10(massParentBin    (2)/massParentBin    (1))
    binRedshiftWidthInverse=1.0d0/     (redshiftMergerBin(2)-redshiftMergerBin(1))
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
       ! Get the particle ID data.
       if (simulations(iSimulation)%propertiesInteger%exists('particleID'     )) then
          particleID      => simulations(iSimulation)%propertiesInteger%value('particleID'     )
       else
          allocate(particleID     (0))
          call Error_Report('particle IDs are required, but are not available in the simulation'     //{introspection:location})
       end if
       ! Get the descendant ID data.
       if (simulations(iSimulation)%propertiesInteger%exists('descendantID'   )) then
          descendantID    => simulations(iSimulation)%propertiesInteger%value('descendantID'   )
       else
          allocate(descendantID   (0))
          call Error_Report('descendant IDs are required, but are not available in the simulation'   //{introspection:location})
       end if
       ! Build a look-up index for parent halo tree IDs.
       countParents  = 0_c_size_t
       redshiftParent=-1.0d0
       do i=1_c_size_t,size(snapshotID,kind=c_size_t)
          if     (                                         &
               &   snapshotID(i) == self%snapshotParents   &
               &  .and.                                    &
               &   mass      (i) >= self%massParentMinimum &
               &  .and.                                    &
               &   mass      (i) <= self%massParentMaximum &
               & ) countParents=+countParents              &
               &                +1_c_size_t
          if    (snapshotID(i) == self%snapshotParents .and. redshiftParent < 0.0d0) redshiftParent=1.0d0/expansionFactor(i)-1.0d0
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
       indexID=sortIndex(descendantID)
       ! Accumulate statistics.
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call displayCounter(0,.true.)
#ifdef USEMPI
       end if
#endif
       ! Iterate over parent halos.
       countBin=0_c_size_t
       !$omp parallel do private(j,k,l,m,n,massesMerger,expansionFactorMerger,redshiftMerger_) reduction(+:countBin) schedule(dynamic)
       do i=1_c_size_t,countParents
#ifdef USEMPI
          ! If running under MPI with N processes, process only every Nth particle.
          if (mod(i,mpiSelf%count()) /= mpiSelf%rank()) cycle
#endif
          ! Skip parent halos that are out of range.
          k=int(log10(mass(indexParents(i))/self%massParentMinimum)*binParentWidthInverse)+1
          if (k < 1 .or. k > massParentCount) cycle
          ! Initialize the merger time to an impossible value.
          expansionFactorMerger=-huge(0.0d0)
          ! Initialize to the final halo.
          j=indexParents(i)
          do while (j > 0)
             ! Find halos that descend into this halo.
             l=searchIndexed(descendantID,indexID,particleID(j))
             ! If no progenitors exist, exit the search.
             if (l < 1_c_size_t .or. l > size(descendantID)) exit
             ! Determine the two most massive progenitors and the index of the most massive progenitor.
             massesMerger=-huge(0.0d0     )
             m           =-huge(0_c_size_t)
             do while (l > 0_c_size_t .and. descendantID(indexID(l)) == particleID(j))
                if (mass(indexID(l)) > massesMerger(1)) then
                   massesMerger(2)=massesMerger(        1 )
                   massesMerger(1)=mass        (indexID(l))
                   m              =             indexID(l)
                else if (mass(indexID(l)) > massesMerger(2)) then
                   massesMerger(2)=mass        (indexID(l))
                end if
                l=l-1
             end do
             ! If the merger mass ratio exceeds threshold, record the merger time and exit.
             if (all(massesMerger > 0.0d0) .and. massesMerger(2) >= self%ratioMajor*massesMerger(1)) then
                expansionFactorMerger=expansionFactor(j)
                exit
             end if
             ! Move to the most massive progenitor.
             j=m             
          end do
          ! Accumulate into bin.
          redshiftMerger_=1.0d0/expansionFactorMerger-1.0d0
          j=int((redshiftMerger_-self%redshiftMergerMinimum)*binRedshiftWidthInverse)+1
          if (j < 1 .or. j > redshiftMergerCount) cycle
          countBin(k,j)=+countBin(k,j) &
                     &  +1_c_size_t
          ! Update progress.
          !$ if (OMP_Get_Thread_Num() == 0) then
#ifdef USEMPI
          if (mpiSelf%isMaster()) then
#endif
             call displayCounter(                          &
                  &              int(                      &
                  &                  +100.0d0              &
                  &                  *float(i           )  &
                  &                  /float(countParents)  &
                  &                 )                    , &
                  &              .false.                   &
                  &             )
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
       redshiftMerger=dble(countBin)*binRedshiftWidthInverse
       do k=1,redshiftMergerCount
          where (countParentBin > 0_c_size_t)
             redshiftMerger(:,k)=redshiftMerger(:,k)/dble(countParentBin)
          elsewhere
             redshiftMerger(:,k)=0.0d0
          end where
       end do     
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          resultsGroup=simulations(iSimulation)%analysis%openGroup("timeLastMajorMerger","Last major merger times of halos.")
          call resultsGroup%writeDataset  (massParentBin            ,'massParent'    )
          call resultsGroup%writeDataset  (redshiftMergerBin        ,'redshiftMerger')
          call resultsGroup%writeAttribute(redshiftParent           ,'redshiftParent')
          call resultsGroup%writeDataset  (countBin                 ,'count'         )
          call resultsGroup%writeDataset  (redshiftMerger           ,'distribution'  )
          call resultsGroup%writeAttribute(self%description         ,"description"   )
          call resultsGroup%writeAttribute(Formatted_Date_and_Time(),"timestamp"     )
          call resultsGroup%close         (                                          )
          cosmologyGroup=simulations(iSimulation)%analysis%openGroup('cosmology')
          call cosmologyGroup%writeAttribute(self%cosmologyParameters_%OmegaMatter    (),'OmegaMatter'    )
          call cosmologyGroup%writeAttribute(self%cosmologyParameters_%OmegaDarkEnergy(),'OmegaDarkEnergy')
          call cosmologyGroup%writeAttribute(self%cosmologyParameters_%HubbleConstant (),'HubbleConstant' )
          call cosmologyGroup%close         (                                                             )
          simulationGroup=simulations(iSimulation)%analysis%openGroup('simulation')
          call simulationGroup%writeAttribute(self%simulationReference,'reference')
          call simulationGroup%writeAttribute(self%simulationURL      ,'URL'      )
          if (simulations(iSimulation)%attributesReal%exists('massParticle')) &
               & call simulationGroup%writeAttribute(simulations(iSimulation)%attributesReal%value('massParticle'),'massParticle')
          if (simulations(iSimulation)%attributesReal%exists('boxSize'     )) &
               & call simulationGroup%writeAttribute(simulations(iSimulation)%attributesReal%value('boxSize'     ),'boxSize'     )
          call simulationGroup%close         (                                    )
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
  end subroutine timeLastMajorMergerOperate

