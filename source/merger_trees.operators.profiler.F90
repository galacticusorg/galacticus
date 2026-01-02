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
  Implements a merger tree operator which profiles tree structure.
  !!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass
  use :: Kind_Numbers       , only : kind_int8

  !![
  <mergerTreeOperator name="mergerTreeOperatorProfiler">
   <description>
    A merger tree operator which profiles merger tree structure.
  </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorProfiler
     !!{
     A merger tree operator class which profiles merger tree structure.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer                     :: cosmologyFunctions_ => null()
     integer         (kind_int8              )                              :: nodeCount                  , singleProgenitorCount
     integer         (kind_int8              ), allocatable, dimension(:,:) :: nonPrimaryProgenitorCount
     integer                                                                :: massBinsCount              , timeBinsCount
     double precision                         , allocatable, dimension(:  ) :: mass                       , time
     double precision                                                       :: massMinimumLogarithmic     , timeMinimumLogarithmic     , &
          &                                                                    massLogarithmicDeltaInverse, timeLogarithmicDeltaInverse, &
          &                                                                    massMinimum                , massMaximum                , &
          &                                                                    redshiftMinimum            , redshiftMaximum
     integer                                                                :: massBinsPerDecade          , timeBinsPerDecade
   contains
     final     ::                        profilerDestructor
     procedure :: operatePreEvolution => profilerOperatePreEvolution
     procedure :: finalize            => profilerFinalize
  end type mergerTreeOperatorProfiler

  interface mergerTreeOperatorProfiler
     !!{
     Constructors for the prune-hierarchy merger tree operator class.
     !!}
     module procedure profilerConstructorParameters
     module procedure profilerConstructorInternal
  end interface mergerTreeOperatorProfiler

contains

  function profilerConstructorParameters(parameters) result(self)
    !!{
    Constructor for the information content merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeOperatorProfiler)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass   ), pointer       :: cosmologyFunctions_
    double precision                                            :: massMinimum        , massMaximum      , &
         &                                                         redshiftMinimum    , redshiftMaximum
    integer                                                     :: massBinsPerDecade  , timeBinsPerDecade

    !![
    <inputParameter>
      <name>massMinimum</name>
      <source>parameters</source>
      <defaultValue>1.0d10</defaultValue>
      <description>The minimum mass of non-primary progenitor to count.</description>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <source>parameters</source>
      <defaultValue>1.0d15</defaultValue>
      <description>The maximum mass of non-primary progenitor to count.</description>
    </inputParameter>
    <inputParameter>
      <name>massBinsPerDecade</name>
      <source>parameters</source>
      <defaultValue>10</defaultValue>
      <description>The number of bins per decade of non-primary progenitor mass.</description>
    </inputParameter>
    <inputParameter>
      <name>redshiftMinimum</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The minimum redshift at which to count non-primary progenitors.</description>
    </inputParameter>
    <inputParameter>
      <name>redshiftMaximum</name>
      <source>parameters</source>
      <defaultValue>10.0d0</defaultValue>
      <description>The maximum redshift at which to count non-primary progenitors.</description>
    </inputParameter>
    <inputParameter>
      <name>timeBinsPerDecade</name>
      <source>parameters</source>
      <defaultValue>10</defaultValue>
      <description>The number of bins per decade of time at which to count non-primary progenitors.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=mergerTreeOperatorProfiler(massMinimum,massMaximum,massBinsPerDecade,redshiftMinimum,redshiftMaximum,timeBinsPerDecade,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function profilerConstructorParameters

  function profilerConstructorInternal(massMinimum,massMaximum,massBinsPerDecade,redshiftMinimum,redshiftMaximum,timeBinsPerDecade,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the information content merger tree operator class.
    !!}
    use :: Numerical_Ranges , only : Make_Range   , rangeTypeLogarithmic
    implicit none
    type            (mergerTreeOperatorProfiler)                        :: self
    double precision                            , intent(in   )         :: massMinimum        , massMaximum      , &
         &                                                                 redshiftMinimum    , redshiftMaximum
    integer                                     , intent(in   )         :: massBinsPerDecade  , timeBinsPerDecade
    class           (cosmologyFunctionsClass   ), intent(in   ), target :: cosmologyFunctions_
    double precision                                                    :: timeMinimum        , timeMaximum
    !![
    <constructorAssign variables="*cosmologyFunctions_, massMinimum, massMaximum, redshiftMinimum, redshiftMaximum, massBinsPerDecade, timeBinsPerDecade"/>
    !!]

    ! Construct bins in mass and time.
    timeMinimum                     =self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum))
    timeMaximum                     =self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshiftMinimum))
    self%massBinsCount              =int(dble(massBinsPerDecade)*log10(massMaximum/massMinimum))+1
    self%timeBinsCount              =int(dble(timeBinsPerDecade)*log10(timeMaximum/timeMinimum))+1
    self%mass                       =Make_Range(massMinimum,massMaximum,self%massBinsCount,rangeType=rangeTypeLogarithmic,rangeBinned=.true.)
    self%time                       =Make_Range(timeMinimum,timeMaximum,self%timeBinsCount,rangeType=rangeTypeLogarithmic,rangeBinned=.true.)
    self%massMinimumLogarithmic     =log10(massMinimum)
    self%timeMinimumLogarithmic     =log10(timeMinimum)
    self%massLogarithmicDeltaInverse=dble(self%massBinsCount)/log10(massMaximum/massMinimum)
    self%timeLogarithmicDeltaInverse=dble(self%timeBinsCount)/log10(timeMaximum/timeMinimum)
    ! Allocate bins.
    allocate(self%nonPrimaryProgenitorCount(self%massBinsCount,self%timeBinsCount))
    ! Initialize counts.
    self%nodeCount                =0_kind_int8
    self%singleProgenitorCount    =0_kind_int8
    self%nonPrimaryProgenitorCount=0_kind_int8
    return
  end function profilerConstructorInternal

  subroutine profilerDestructor(self)
    !!{
    Destructor for the profiler merger tree operator function class.
    !!}
    implicit none
    type(mergerTreeOperatorProfiler), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine profilerDestructor

  subroutine profilerOperatePreEvolution(self,tree)
    !!{
    Perform a information content operation on a merger tree.
    !!}
    use :: Galacticus_Nodes   , only : mergerTree                   , nodeComponentBasic, treeNode
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerIsolatedNodes
    implicit none
    class  (mergerTreeOperatorProfiler   ), intent(inout), target   :: self
    type   (mergerTree                   ), intent(inout), target  :: tree
    type   (treeNode                     )               , pointer :: node
    class  (nodeComponentBasic           )               , pointer :: basic
    type   (mergerTreeWalkerIsolatedNodes)                         :: treeWalker
    integer                                                        :: massIndex  , timeIndex

    ! Iterate over nodes.
    treeWalker=mergerTreeWalkerIsolatedNodes(tree,spanForest=.true.)
    do while (treeWalker%next(node))
       ! Count nodes.
       self%nodeCount=self%nodeCount+1_kind_int8
       ! Check for single progenitor nodes.
       if (associated(node%firstChild)) then
          if (associated(node%firstChild%sibling)) then
             ! Count non-primary progenitors as a function of mass and epoch.
             basic     => node%basic()
             massIndex =  int((log10(basic%mass())-self%massMinimumLogarithmic)*self%massLogarithmicDeltaInverse)+1
             timeIndex =  int((log10(basic%time())-self%timeMinimumLogarithmic)*self%timeLogarithmicDeltaInverse)+1
             if     (                                                     &
                  &   massIndex > 0 .and. massIndex <= self%massBinsCount &
                  &  .and.                                                &
                  &   timeIndex > 0 .and. timeIndex <= self%timeBinsCount &
                  & ) self%nonPrimaryProgenitorCount(massIndex,timeIndex)=self%nonPrimaryProgenitorCount(massIndex,timeIndex)+1_kind_int8
          else
             ! Check for single progenitor nodes.
             self%singleProgenitorCount=self%singleProgenitorCount+1_kind_int8
          end if
       end if
     end do
    return
  end subroutine profilerOperatePreEvolution

  subroutine profilerFinalize(self)
    !!{
    Outputs tree information content function.
    !!}
    use :: Output_HDF5, only : outputFile
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Object
    implicit none
    class  (mergerTreeOperatorProfiler), intent(inout)               :: self
    type   (hdf5Object                )                              :: profilerGroup
    integer(kind_int8                 )                              :: nodeCountCurrent                , singleProgenitorCountCurrent
    integer(kind_int8                 ), dimension(:,:), allocatable :: nonPrimaryProgenitorCountCurrent

    !$ call hdf5Access%set()
    ! Output information content information.
    profilerGroup=outputFile%openGroup('treeProfiler','Profiling information on merger trees.',objectsOverwritable=.true.,overwriteOverride=.true.)
    if (profilerGroup%hasAttribute('nodeCount'                )) then
       call profilerGroup%readAttribute('nodeCount'              ,nodeCountCurrent                  )
       self%nodeCount                =self%nodeCount                +nodeCountCurrent
    end if
    if (profilerGroup%hasAttribute('singleProgenitorCount'    )) then
       call profilerGroup%readAttribute('singleProgenitorCount'   ,singleProgenitorCountCurrent     )
       self%singleProgenitorCount    =self%singleProgenitorCount    +singleProgenitorCountCurrent
    end if
   if (profilerGroup%hasDataset   ('nonPrimaryProgenitorCount')) then
       call profilerGroup%readDataset  ('nonPrimaryProgenitorCount',nonPrimaryProgenitorCountCurrent)
       self%nonPrimaryProgenitorCount=self%nonPrimaryProgenitorCount+nonPrimaryProgenitorCountCurrent
    end if
    call profilerGroup%writeAttribute(self%nodeCount                ,'nodeCount'                )
    call profilerGroup%writeAttribute(self%singleProgenitorCount    ,'singleProgenitorCount'    )
    call profilerGroup%writeDataset  (self%nonPrimaryProgenitorCount,'nonPrimaryProgenitorCount')
    call profilerGroup%writeDataset  (self%mass                     ,'nonPrimaryProgenitorMass' )
    call profilerGroup%writeDataset  (self%time                     ,'nonPrimaryProgenitorTime' )
    call profilerGroup%close         (                                                          )
    !$ call hdf5Access%unset()
    return
  end subroutine profilerFinalize
