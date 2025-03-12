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
Implements an output analysis property extractor class that extracts a property from a descendant node of the given node.
!!}
  
  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorDescendantNode">
   <description>An output analysis property extractor class that extracts a property from a descendant node of the given node.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorDescendantNode
     !!{
     A property extractor output analysis class that extracts a property from a descendant node of the given node.
     !!}
     private
     class           (cosmologyFunctionsClass    ), pointer :: cosmologyFunctions_    => null()
     class           (nodePropertyExtractorScalar), pointer :: nodePropertyExtractor_ => null()
     double precision                                       :: timeDescendant                  , redshiftDescendant
   contains
     final     ::                descendantNodeDestructor
     procedure :: extract     => descendantNodeExtract
     procedure :: name        => descendantNodeName
     procedure :: description => descendantNodeDescription
     procedure :: unitsInSI   => descendantNodeUnitsInSI
  end type nodePropertyExtractorDescendantNode

  interface nodePropertyExtractorDescendantNode
     !!{
     Constructors for the ``descendantNode'' node property extractor class.
     !!}
     module procedure descendantNodeConstructorParameters
     module procedure descendantNodeConstructorInternal
  end interface nodePropertyExtractorDescendantNode

contains

  function descendantNodeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``descendantNode'' node property extractor class which takes a parameter set as input.
    !!}
    use :: Error              , only : Error_Report
    use :: Input_Parameters   , only : inputParameters
    implicit none
    type            (nodePropertyExtractorDescendantNode)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass            ), pointer       :: cosmologyFunctions_
    class           (nodePropertyExtractorClass         ), pointer       :: nodePropertyExtractor_
    double precision                                                     :: redshiftDescendant
    
    !![
    <inputParameter>
      <name>redshiftDescendant</name>
      <source>parameters</source>
      <description>The redshift of the descendant node to which to apply the filter.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="nodePropertyExtractor" name="nodePropertyExtractor_" source="parameters"/>
    !!]
    select type (nodePropertyExtractor_)
    class is (nodePropertyExtractorScalar)
       self=nodePropertyExtractorDescendantNode(cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftDescendant)),cosmologyFunctions_,nodePropertyExtractor_)
    class default
       call Error_Report('extracted property must be a real scalar'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="nodePropertyExtractor_"/>
    !!]
    return
  end function descendantNodeConstructorParameters

  function descendantNodeConstructorInternal(timeDescendant,cosmologyFunctions_,nodePropertyExtractor_) result(self)
    !!{
    Internal constructor for the ``descendantNode'' node property extractor class.
    !!}
    implicit none
    type            (nodePropertyExtractorDescendantNode)                        :: self
    double precision                                     , intent(in   )         :: timeDescendant
    class           (cosmologyFunctionsClass            ), intent(in   ), target :: cosmologyFunctions_
    class           (nodePropertyExtractorScalar        ), intent(in   ), target :: nodePropertyExtractor_
    !![
    <constructorAssign variables="timeDescendant, *cosmologyFunctions_, *nodePropertyExtractor_"/>
    !!]

    self%redshiftDescendant=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeDescendant))
    return
  end function descendantNodeConstructorInternal
  
  subroutine descendantNodeDestructor(self)
    !!{
    Destructor for  the ``descendantNode'' node property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorDescendantNode), intent(inout) :: self

    !![
    <objectDestructor name="self%nodePropertyExtractor_"/>
    <objectDestructor name="self%cosmologyFunctions_"   />
    !!]
    return
  end subroutine descendantNodeDestructor
  
  double precision function descendantNodeExtract(self,node,instance)
    !!{
    Implement a descendantNode output analysis.
    !!}
    use :: Error               , only : Error_Report
    use :: Galacticus_Nodes    , only : nodeComponentBasic
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    class           (nodePropertyExtractorDescendantNode), intent(inout), target   :: self
    type            (treeNode                           ), intent(inout), target   :: node
    type            (multiCounter                       ), intent(inout), optional :: instance
    type            (treeNode                           ), pointer                 :: nodeDescendant
    class           (nodeComponentBasic                 ), pointer                 :: basicDescendant
    double precision                                     , parameter               :: timeTolerance  =1.0d-4

    descendantNodeExtract=-huge(0.0d0)
    nodeDescendant => node%parent
    do while (associated(nodeDescendant))
       basicDescendant => nodeDescendant%basic()
       if (Values_Agree(basicDescendant%time(),self%timeDescendant,absTol=timeTolerance)) then
          descendantNodeExtract=self%nodePropertyExtractor_%extract(nodeDescendant,instance)
          return
       end if
       nodeDescendant => nodeDescendant%parent
    end do
    call Error_Report('failed to find descendant node'//{introspection:location})
    return
  end function descendantNodeExtract


  function descendantNodeName(self)
    !!{
    Return the name of the descendantNode property.
    !!}
    use :: String_Handling, only : String_Upper_Case_First
    implicit none
    type (varying_string                     )                :: descendantNodeName
    class(nodePropertyExtractorDescendantNode), intent(inout) :: self

    descendantNodeName=var_str('descendant')//String_Upper_Case_First(char(self%nodePropertyExtractor_%name()))
    return
  end function descendantNodeName

  function descendantNodeDescription(self)
    !!{
    Return a description of the descendantNode property.
    !!}
    implicit none
    type (varying_string                     )                :: descendantNodeDescription
    class(nodePropertyExtractorDescendantNode), intent(inout) :: self

    descendantNodeDescription=self%nodePropertyExtractor_%description()//' (of a descendant node)'
    return
  end function descendantNodeDescription

  double precision function descendantNodeUnitsInSI(self)
    !!{
    Return the units of the descendantNode property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorDescendantNode), intent(inout) :: self

    descendantNodeUnitsInSI=self%nodePropertyExtractor_%unitsInSI()
    return
  end function descendantNodeUnitsInSI
