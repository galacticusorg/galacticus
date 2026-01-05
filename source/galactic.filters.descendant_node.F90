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
Implements a galactic filter which applies another filter to a descendant node of the given node.
!!}
  
  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <galacticFilter name="galacticFilterDescendantNode">
   <description>Applies a filter to a descendant node of the given node.</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterDescendantNode
     !!{
     A galactic filter which applies another filter to a descendant node of the given node.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     class           (galacticFilterClass    ), pointer :: galacticFilter_     => null()
     double precision                                   :: timeDescendant               , redshiftDescendant
     logical                                            :: allowSelf
   contains
     final     ::           descendantNodeDestructor
     procedure :: passes => descendantNodePasses
  end type galacticFilterDescendantNode

  interface galacticFilterDescendantNode
     !!{
     Constructors for the \refClass{galacticFilterDescendantNode} galactic filter class.
     !!}
     module procedure descendantNodeConstructorParameters
     module procedure descendantNodeConstructorInternal
  end interface galacticFilterDescendantNode

contains

  function descendantNodeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterDescendantNode} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters   , only : inputParameter         , inputParameters
    implicit none
    type            (galacticFilterDescendantNode)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass     ), pointer       :: cosmologyFunctions_
    class           (galacticFilterClass         ), pointer       :: galacticFilter_
    double precision                                              :: redshiftDescendant
    logical                                                       :: allowSelf
         
    !![
    <inputParameter>
      <name>redshiftDescendant</name>
      <source>parameters</source>
      <description>The redshift of the descendant node to which to apply the filter.</description>
    </inputParameter>
    <inputParameter>
      <name>allowSelf</name>
      <source>parameters</source>
      <description>If true, the node itself is considered as a possible descendant, otherwise the node itself is excluded from the descendant node search.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="galacticFilter"      name="galacticFilter_"      source="parameters"/>
    !!]
    self=galacticFilterDescendantNode(cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftDescendant)),allowSelf,cosmologyFunctions_,galacticFilter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="galacticFilter_"    />
    !!]
    return
  end function descendantNodeConstructorParameters
  
  function descendantNodeConstructorInternal(timeDescendant,allowSelf,cosmologyFunctions_,galacticFilter_) result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterDescendantNode} galactic filter class.
    !!}
    implicit none
    type            (galacticFilterDescendantNode)                        :: self
    double precision                              , intent(in   )         :: timeDescendant
    logical                                       , intent(in   )         :: allowSelf
    class           (galacticFilterClass         ), intent(in   ), target :: galacticFilter_
    class           (cosmologyFunctionsClass     ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="timeDescendant, allowSelf, *cosmologyFunctions_, *galacticFilter_"/>
    !!]

    self%redshiftDescendant=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeDescendant))
    return
  end function descendantNodeConstructorInternal
  
  subroutine descendantNodeDestructor(self)
    !!{
    Destructor for the \refClass{galacticFilterDescendantNode} galactic filter class.
    !!}
    implicit none
    type(galacticFilterDescendantNode), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    <objectDestructor name="self%galacticFilter_"    />
    !!]
    return
  end subroutine descendantNodeDestructor
  
  logical function descendantNodePasses(self,node)
    !!{
    Implement a filter on descendant node properties.
    !!}
    use :: Error               , only : Error_Report
    use :: Galacticus_Nodes    , only : nodeComponentBasic
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    class           (galacticFilterDescendantNode), intent(inout)         :: self
    type            (treeNode                    ), intent(inout), target :: node
    type            (treeNode                    ), pointer               :: nodeDescendant
    class           (nodeComponentBasic          ), pointer               :: basicDescendant
    double precision                              , parameter             :: timeTolerance  =1.0d-4

    descendantNodePasses=.false.
    if (self%allowSelf) then
       nodeDescendant => node
    else
       nodeDescendant => node%parent
    end if
    do while (associated(nodeDescendant))
       basicDescendant => nodeDescendant%basic()
       if (Values_Agree(basicDescendant%time(),self%timeDescendant,absTol=timeTolerance)) then
          descendantNodePasses=self%galacticFilter_%passes(nodeDescendant)
          return
       end if
       nodeDescendant => nodeDescendant%parent
    end do
    call Error_Report('failed to find descendant node'//{introspection:location})
    return
  end function descendantNodePasses
