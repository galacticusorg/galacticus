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
Implements an interval pass filter on halo merger ratio.
!!}

  !![
  <galacticFilter name="galacticFilterMergerRatio">
   <description>An interval pass filter on halo merger ratio.</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterMergerRatio
     !!{
     an interval pass galactic filter class on merger ratio.
     !!}
     private
     double precision :: ratioLow, ratioHigh
   contains
     procedure :: passes => mergerRatioPasses
  end type galacticFilterMergerRatio

  interface galacticFilterMergerRatio
     !!{
     Constructors for the \refClass{galacticFilterMergerRatio} galactic filter class.
     !!}
     module procedure mergerRatioConstructorParameters
     module procedure mergerRatioConstructorInternal
  end interface galacticFilterMergerRatio

contains
  
  function mergerRatioConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterMergerRatio} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterMergerRatio )                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    double precision                                            :: ratioLow  , ratioHigh
    
    !![
    <inputParameter>
      <name>ratioLow</name>
      <source>parameters</source>
      <description>The low ratio value above which to pass.</description>
    </inputParameter>
    <inputParameter>
      <name>ratioHigh</name>
      <source>parameters</source>
      <description>The high ratio value below which to pass.</description>
    </inputParameter>
    !!]
    self=galacticFilterMergerRatio(ratioLow,ratioHigh)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function mergerRatioConstructorParameters

  function mergerRatioConstructorInternal(ratioLow,ratioHigh) result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterMergerRatio} galactic filter class.
    !!}
    implicit none
    type            (galacticFilterMergerRatio)                :: self
    double precision                           , intent(in   ) :: ratioLow, ratioHigh
    !![
    <constructorAssign variables="ratioLow, ratioHigh"/>
    !!]
    
    return
  end function mergerRatioConstructorInternal

  logical function mergerRatioPasses(self,node)
    !!{
    Implement an interval pass galactic filter.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (galacticFilterMergerRatio), intent(inout)          :: self
    type            (treeNode                 ), intent(inout), target  :: node
    class           (nodeComponentBasic       )               , pointer :: basicPrimary, basicSecondary
    double precision                                                    :: mergerRatio

    mergerRatioPasses=.false.
    if (.not.associated(node%firstChild        )) return
    if (.not.associated(node%firstChild%sibling)) return
    basicPrimary      =>  node%firstChild        %basic()
    basicSecondary    =>  node%firstChild%sibling%basic()
    mergerRatio       =  +basicSecondary%mass() &
         &               /basicPrimary  %mass()
    mergerRatioPasses =   mergerRatio >= self%ratioLow  &
         &               .and.                          &
         &                mergerRatio <  self%ratioHigh
    return
  end function mergerRatioPasses
