!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which implements a galactic filter which applies another filter to a descendent node of the given node.
!!}
  
  !![
  <galacticFilter name="galacticFilterDescendentNode">
   <description>Applies a filter to a descendent node of the given node.</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterDescendentNode
     !!{
     A galactic filter which applies another filter to a descendent node of the given node.
     !!}
     private
     class           (galacticFilterClass), pointer :: galacticFilter_
     double precision                               :: timeDescendent
     logical                                        :: allowSelf
   contains
     final     ::           descendentNodeDestructor
     procedure :: passes => descendentNodePasses
  end type galacticFilterDescendentNode

  interface galacticFilterDescendentNode
     !!{
     Constructors for the ``descendentNode'' galactic filter class.
     !!}
     module procedure descendentNodeConstructorParameters
     module procedure descendentNodeConstructorInternal
  end interface galacticFilterDescendentNode

contains

  function descendentNodeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``descendentNode'' galactic filter class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions, only : cosmologyFunctionsClass
    use :: Input_Parameters   , only : inputParameter         , inputParameters
    implicit none
    type            (galacticFilterDescendentNode)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass     ), pointer       :: cosmologyFunctions_
    class           (galacticFilterClass         ), pointer       :: galacticFilter_
    double precision                                              :: redshiftDescendent
    logical                                                       :: allowSelf
         
    !![
    <inputParameter>
      <name>redshiftDescendent</name>
      <source>parameters</source>
      <description>The redshift of the descendent node to which to apply the filter.</description>
    </inputParameter>
    <inputParameter>
      <name>allowSelf</name>
      <source>parameters</source>
      <description>If true, the node itself is considered as a possible descendent, otherwise the node itself is excldued from the descendent node search.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="galacticFilter"      name="galacticFilter_"      source="parameters"/>
    !!]
    self=galacticFilterDescendentNode(cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftDescendent)),allowSelf,galacticFilter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="galacticFilter_"    />
    !!]
    return
  end function descendentNodeConstructorParameters
  
  function descendentNodeConstructorInternal(timeDescendent,allowSelf,galacticFilter_) result(self)
    !!{
    Internal constructor for the ``descendentNode'' galactic filter class.
    !!}
    implicit none
    type            (galacticFilterDescendentNode)                        :: self
    double precision                              , intent(in   )         :: timeDescendent
    logical                                       , intent(in   )         :: allowSelf
    class           (galacticFilterClass         ), intent(in   ), target :: galacticFilter_
    !![
    <constructorAssign variables="timeDescendent, allowSelf, *galacticFilter_"/>
    !!]

    return
  end function descendentNodeConstructorInternal
  
  subroutine descendentNodeDestructor(self)
    !!{
    Destructor for  the ``descendentNode'' galactic filter class.
    !!}
    implicit none
    type(galacticFilterDescendentNode), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticFilter_"/>
    !!]
    return
  end subroutine descendentNodeDestructor
  
  logical function descendentNodePasses(self,node)
    !!{
    Implement a filter on descendent node properties.
    !!}
    use :: Galacticus_Error    , only : Galacticus_Error_Report
    use :: Galacticus_Nodes    , only : nodeComponentBasic
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    class           (galacticFilterDescendentNode), intent(inout)         :: self
    type            (treeNode                    ), intent(inout), target :: node
    type            (treeNode                    ), pointer               :: nodeDescendent
    class           (nodeComponentBasic          ), pointer               :: basicDescendent
    double precision                              , parameter             :: timeTolerance  =1.0d-4

    descendentNodePasses=.false.
    if (self%allowSelf) then
       nodeDescendent => node
    else
       nodeDescendent => node%parent
    end if
    do while (associated(nodeDescendent))
       basicDescendent => nodeDescendent%basic()
       if (Values_Agree(basicDescendent%time(),self%timeDescendent,absTol=timeTolerance)) then
          descendentNodePasses=self%galacticFilter_%passes(nodeDescendent)
          return
       end if
       nodeDescendent => nodeDescendent%parent
    end do
    call Galacticus_Error_Report('failed to find descendent node'//{introspection:location})
    return
  end function descendentNodePasses
