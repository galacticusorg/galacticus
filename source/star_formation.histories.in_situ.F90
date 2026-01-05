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
Implements a star formation histories class which records \emph{in situ} star formation alongside the star formation history computed by some other class.
!!}

  !![
  <starFormationHistory name="starFormationHistoryInSitu">
   <description>
     A star formation histories class which records \emph{in situ} star formation. Another {\normalfont \ttfamily
     starFormationHistory} object is used to provide the base star formation history. This class tracks a second copy which is
     identical but excludes any star formation from merging galaxies.
   </description>
  </starFormationHistory>
  !!]
  type, extends(starFormationHistoryClass) :: starFormationHistoryInSitu
     !!{
     A star formation histories class which records \emph{in situ} star formation.
     !!}
     private
     class(starFormationHistoryClass), pointer :: starFormationHistory_ => null()
   contains
     final     ::                          inSituDestructor
     procedure :: create                => inSituCreate
     procedure :: rate                  => inSituRate
     procedure :: update                => inSituUpdate
     procedure :: scales                => inSituScales
     procedure :: autoHook              => inSituAutoHook
     procedure :: times                 => inSituTimes
     procedure :: metallicityBoundaries => inSituMetallicityBoundaries
     procedure :: ageDistribution       => inSituAgeDistribution
     procedure :: rangeIsSufficient     => inSituRangeIsSufficient
     procedure :: extend                => inSituExtend
  end type starFormationHistoryInSitu

  interface starFormationHistoryInSitu
     !!{
     Constructors for the \refClass{starFormationHistoryInSitu} star formation history class.
     !!}
     module procedure inSituConstructorParameters
     module procedure inSituConstructorInternal
  end interface starFormationHistoryInSitu

contains

  function inSituConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{starFormationHistoryInSitu} star formation history class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (starFormationHistoryInSitu)                :: self
    type (inputParameters           ), intent(inout) :: parameters
    class(starFormationHistoryClass ), pointer       :: starFormationHistory_

    !![
    <objectBuilder class="starFormationHistory" name="starFormationHistory_" source="parameters"/>
    !!]
    self=starFormationHistoryInSitu(starFormationHistory_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationHistory_"/>
    !!]
    return
  end function inSituConstructorParameters

  function inSituConstructorInternal(starFormationHistory_) result(self)
    !!{
    Internal constructor for the \refClass{starFormationHistoryInSitu} star formation history class.
    !!}
    implicit none
    type (starFormationHistoryInSitu)                        :: self
    class(starFormationHistoryClass ), intent(in   ), target :: starFormationHistory_
    !![
    <constructorAssign variables="*starFormationHistory_"/>
    !!]

    return
  end function inSituConstructorInternal

  subroutine inSituAutoHook(self)
    !!{
    Attach to the satellite merging event hook.
    !!}
    use :: Events_Hooks, only : openMPThreadBindingAllLevels, satelliteMergerEvent
    implicit none
    class(starFormationHistoryInSitu), intent(inout) :: self

    call satelliteMergerEvent%attach(self,inSituSatelliteMerger,openMPThreadBindingAllLevels,label='starFormationHistoryInSitu')
    return
  end subroutine inSituAutoHook

  subroutine inSituDestructor(self)
    !!{
    Destructor for the \refClass{starFormationHistoryInSitu} star formation histories class.
    !!}
    implicit none
    type(starFormationHistoryInSitu), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationHistory_"/>
    !!]
    return
  end subroutine inSituDestructor

  subroutine inSituCreate(self,node,historyStarFormation,timeBegin,timeEnd)
    !!{
    Create the history required for storing star formation history.
    !!}
    implicit none
    class           (starFormationHistoryInSitu), intent(inout)           :: self
    type            (treeNode                  ), intent(inout), target   :: node
    type            (history                   ), intent(inout)           :: historyStarFormation
    double precision                            , intent(in   )           :: timeBegin
    double precision                            , intent(in   ), optional :: timeEnd
    type            (history                   )                          :: history_
    
    call self%starFormationHistory_%create(node,history_,timeBegin,timeEnd)
    historyStarFormation%rangeType=history_%rangeType
    historyStarFormation%time     =history_%time
    allocate(historyStarFormation%data(size(history_%data,dim=1),2*size(history_%data,dim=2)))
    historyStarFormation%data(:,1                          :  size(history_%data,dim=2))=history_%data
    historyStarFormation%data(:,1+size(history_%data,dim=2):2*size(history_%data,dim=2))=history_%data
    return
  end subroutine inSituCreate

  subroutine inSituRate(self,node,historyStarFormation,abundancesFuel,rateStarFormation)
    !!{
    Set the rate the star formation history for {\normalfont \ttfamily node}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (starFormationHistoryInSitu), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    type            (history                   ), intent(inout) :: historyStarFormation
    type            (abundances                ), intent(in   ) :: abundancesFuel
    double precision                            , intent(in   ) :: rateStarFormation
    type            (history                   )                :: history_

    ! Check if history exists.
    if (historyStarFormation%exists()) then
       allocate(history_%time(size(historyStarFormation%data,dim=1)                                        ))
       allocate(history_%data(size(historyStarFormation%data,dim=1),size(historyStarFormation%data,dim=2)/2))
       history_%rangeType=historyStarFormation%rangeType
       history_%time     =historyStarFormation%time
       history_%data=historyStarFormation%data(:,1:size(historyStarFormation%data,dim=2)/2)
       call self%starFormationHistory_%rate(node,history_,abundancesFuel,rateStarFormation)
       historyStarFormation%data(:,1                                        :size(historyStarFormation%data,dim=2)/2)=history_%data
       historyStarFormation%data(:,1+size(historyStarFormation%data,dim=2)/2:size(historyStarFormation%data,dim=2)  )=history_%data
    else
       ! No history exists - this is acceptable only if the star formation rate is zero.
       if (rateStarFormation > 0.0d0) call Error_Report('non-zero star formation rate, but star formation history is uninitialized'//{introspection:location})
    end if
    return
  end subroutine inSituRate
  
  subroutine inSituUpdate(self,node,indexOutput,historyStarFormation)
    !!{
    Update the star formation history after outputting.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class  (starFormationHistoryInSitu), intent(inout)         :: self
    type   (treeNode                  ), intent(inout), target :: node
    type   (history                   ), intent(inout)         :: historyStarFormation
    integer(c_size_t                  ), intent(in   )         :: indexOutput
    type   (history                   )                        :: history1            , history2

    if (.not.historyStarFormation%exists()) return
    allocate(history1%time(size(historyStarFormation%data,dim=1)                                        ))
    allocate(history2%time(size(historyStarFormation%data,dim=1)                                        ))
    allocate(history1%data(size(historyStarFormation%data,dim=1),size(historyStarFormation%data,dim=2)/2))
    allocate(history2%data(size(historyStarFormation%data,dim=1),size(historyStarFormation%data,dim=2)/2))
    history1%rangeType=historyStarFormation%rangeType
    history2%rangeType=historyStarFormation%rangeType
    history1%time     =historyStarFormation%time
    history2%time     =historyStarFormation%time
    history1%data=historyStarFormation%data(:,1                                        :size(historyStarFormation%data,dim=2)/2)
    history2%data=historyStarFormation%data(:,1+size(historyStarFormation%data,dim=2)/2:size(historyStarFormation%data,dim=2)  )
    call self%starFormationHistory_%update(node,indexOutput,history1)
    call self%starFormationHistory_%update(node,indexOutput,history2)
    call historyStarFormation%destroy()
    allocate(historyStarFormation%time(size(history1%data,dim=1)                            ))
    allocate(historyStarFormation%data(size(history1%data,dim=1),size(history1%data,dim=2)*2))
    historyStarFormation%rangeType=history1%rangeType
    historyStarFormation%time     =history1%time
    historyStarFormation%data(:,1                          :  size(history1%data,dim=2))=history1%data
    historyStarFormation%data(:,1+size(history1%data,dim=2):2*size(history1%data,dim=2))=history2%data
    return
  end subroutine inSituUpdate

  subroutine inSituScales(self,historyStarFormation,node,massStellar,massGas,abundancesStellar)
    !!{
    Set the scalings for error control on the absolute values of star formation histories.
    !!}
    implicit none
    class           (starFormationHistoryInSitu), intent(inout) :: self
    double precision                            , intent(in   ) :: massStellar         , massGas
    type            (abundances                ), intent(in   ) :: abundancesStellar
    type            (history                   ), intent(inout) :: historyStarFormation
    type            (treeNode                  ), intent(inout) :: node
    type            (history                   )                :: history_

    allocate(history_%time(size(historyStarFormation%data,dim=1)                                        ))
    allocate(history_%data(size(historyStarFormation%data,dim=1),size(historyStarFormation%data,dim=2)/2))
    history_%rangeType=historyStarFormation%rangeType
    history_%time     =historyStarFormation%time
    history_%data=historyStarFormation%data(:,1:size(historyStarFormation%data,dim=2)/2)
    call self%starFormationHistory_%scales(history_,node,massStellar,massGas,abundancesStellar)
    historyStarFormation%data(:,1                                        :size(historyStarFormation%data,dim=2)/2)=history_%data
    historyStarFormation%data(:,1+size(historyStarFormation%data,dim=2)/2:size(historyStarFormation%data,dim=2)  )=history_%data
    return
  end subroutine inSituScales

  subroutine inSituSatelliteMerger(self,node)
    !!{
    Zero any in-situ star formation history for the galaxy about to merge.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid, treeNode
    implicit none
    class(starFormationHistoryInSitu), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node
    class(nodeComponentDisk         ), pointer       :: disk
    class(nodeComponentSpheroid     ), pointer       :: spheroid
    type (history                   )                :: historyStarFormationDisk, historyStarFormationSpheroid

    select type (self)
    class is (starFormationHistoryInSitu)
       disk                                   => node    %disk                ()
       spheroid                               => node    %spheroid            ()
       historyStarFormationDisk               =  disk    %starFormationHistory()
       historyStarFormationSpheroid           =  spheroid%starFormationHistory()
       if (historyStarFormationDisk    %exists()) then
          historyStarFormationDisk    %data(:,1+size(historyStarFormationDisk    %data,dim=2)/2:size(historyStarFormationDisk    %data,dim=2))=0.0d0
          call disk    %starFormationHistorySet(    historyStarFormationDisk)
       end if
       if (historyStarFormationSpheroid%exists()) then
          historyStarFormationSpheroid%data(:,1+size(historyStarFormationSpheroid%data,dim=2)/2:size(historyStarFormationSpheroid%data,dim=2))=0.0d0
          call spheroid%starFormationHistorySet(historyStarFormationSpheroid)
       end if
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine inSituSatelliteMerger

  function inSituMetallicityBoundaries(self)
    !!{
    Return the boundaries of the metallicities used in this tabulation.
    !!}
    implicit none
    double precision                            , allocatable  , dimension(:) :: inSituMetallicityBoundaries
    class           (starFormationHistoryInSitu), intent(inout)               :: self

    inSituMetallicityBoundaries=self%starFormationHistory_%metallicityBoundaries()
    return
  end function inSituMetallicityBoundaries

  function inSituAgeDistribution(self) result(ageDistribution)
    !!{
    Return true since the tabulation (in time and metallicity) is static (independent of node) per output.
    !!}
    implicit none
    type (enumerationStarFormationHistoryAgesType)                :: ageDistribution
    class(starFormationHistoryInSitu             ), intent(inout) :: self

    ageDistribution=self%starFormationHistory_%ageDistribution()
    return
  end function inSituAgeDistribution

  function inSituTimes(self,node,indexOutput,starFormationHistory,allowTruncation,timeStart)
    !!{
    Return the times used in this tabulation.
    !!}
    implicit none
    double precision                            , allocatable  , dimension(:) :: inSituTimes
    class           (starFormationHistoryInSitu), intent(inout)               :: self
    type            (treeNode                  ), intent(inout), optional     :: node
    integer         (c_size_t                  ), intent(in   ), optional     :: indexOutput
    type            (history                   ), intent(in   ), optional     :: starFormationHistory
    logical                                     , intent(in   ), optional     :: allowTruncation
    double precision                            , intent(  out), optional     :: timeStart

    inSituTimes=self%starFormationHistory_%times(node,indexOutput,starFormationHistory,allowTruncation,timeStart)
    return
  end function inSituTimes

  logical function inSituRangeIsSufficient(self,starFormationHistory,rangeHistory) result(rangeIsSufficient)
    !!{
    Return true if the range of this history is sufficient.
    !!}
    implicit none
    class(starFormationHistoryInSitu), intent(inout) :: self
    type (history                   ), intent(in   ) :: starFormationHistory, rangeHistory

    rangeIsSufficient=self%starFormationHistory_%rangeIsSufficient(starFormationHistory,rangeHistory)
    return
  end function inSituRangeIsSufficient

  subroutine inSituExtend(self,starFormationHistory,times)
    !!{
    Extend this history to span a sufficient range.
    !!}
    implicit none
    class           (starFormationHistoryInSitu), intent(inout)               :: self
    type            (history                   ), intent(inout)               :: starFormationHistory
    double precision                            , intent(in   ), dimension(:) :: times
    type            (history                   )                              :: historyOriginal     , historyInSitu

    allocate(historyOriginal%time(size(starFormationHistory%data,dim=1)                                        ))
    allocate(historyInsitu  %time(size(starFormationHistory%data,dim=1)                                        ))
    allocate(historyOriginal%data(size(starFormationHistory%data,dim=1),size(starFormationHistory%data,dim=2)/2))
    allocate(historyInsitu  %data(size(starFormationHistory%data,dim=1),size(starFormationHistory%data,dim=2)/2))
    historyOriginal%rangeType=starFormationHistory%rangeType
    historyInsitu  %rangeType=starFormationHistory%rangeType
    historyOriginal%time     =starFormationHistory%time
    historyInsitu  %time     =starFormationHistory%time
    historyOriginal%data     =starFormationHistory%data     (:,1:size(starFormationHistory%data,dim=2)/2)
    historyInsitu  %data     =starFormationHistory%data     (:,1:size(starFormationHistory%data,dim=2)/2)
    call self%starFormationHistory_%extend(historyOriginal,times)
    call self%starFormationHistory_%extend(historyInsitu  ,times)
    call starFormationHistory%destroy()
    allocate(starFormationHistory%time(size(historyOriginal%data,dim=1)                                   ))
    allocate(starFormationHistory%data(size(historyOriginal%data,dim=1),size(historyOriginal%data,dim=2)*2))
    starFormationHistory%time                                                                                     =historyOriginal%time
    starFormationHistory%data(:,1                                        :size(starFormationHistory%data,dim=2)/2)=historyOriginal%data
    starFormationHistory%data(:,1+size(starFormationHistory%data,dim=2)/2:size(starFormationHistory%data,dim=2)  )=historyInsitu  %data
    return
  end subroutine inSituExtend
