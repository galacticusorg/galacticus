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
Contains a module which implements an output analysis property extractor class that extracts a property from a descendent node of the given node.
!!}
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorDescendentNode">
   <description>An output analysis property extractor class that extracts a property from a descendent node of the given node.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorDescendentNode
     !!{
     A property extractor output analysis class that extracts a property from a descendent node of the given node.
     !!}
     private
     class           (nodePropertyExtractorScalar), pointer :: nodePropertyExtractor_ => null()
     double precision                                       :: timeDescendent
   contains
     final     ::                descendentNodeDestructor
     procedure :: extract     => descendentNodeExtract
     procedure :: type        => descendentNodeType
     procedure :: name        => descendentNodeName
     procedure :: description => descendentNodeDescription
     procedure :: unitsInSI   => descendentNodeUnitsInSI
  end type nodePropertyExtractorDescendentNode

  interface nodePropertyExtractorDescendentNode
     !!{
     Constructors for the ``descendentNode'' node property extractor class.
     !!}
     module procedure descendentNodeConstructorParameters
     module procedure descendentNodeConstructorInternal
  end interface nodePropertyExtractorDescendentNode

contains

  function descendentNodeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``descendentNode'' node property extractor class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions, only : cosmologyFunctionsClass
    use :: Galacticus_Error   , only : Galacticus_Error_Report
    use :: Input_Parameters   , only : inputParameters
    implicit none
    type            (nodePropertyExtractorDescendentNode)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass            ), pointer       :: cosmologyFunctions_
    class           (nodePropertyExtractorClass         ), pointer       :: nodePropertyExtractor_
    double precision                                                     :: redshiftDescendent
    
    !![
    <inputParameter>
      <name>redshiftDescendent</name>
      <source>parameters</source>
      <description>The redshift of the descendent node to which to apply the filter.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="nodePropertyExtractor" name="nodePropertyExtractor_" source="parameters"/>
    !!]
    select type (nodePropertyExtractor_)
    class is (nodePropertyExtractorScalar)
       self=nodePropertyExtractorDescendentNode(cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftDescendent)),nodePropertyExtractor_)
    class default
       call Galacticus_Error_Report('extracted property must be a real scalar'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="nodePropertyExtractor_"/>
    !!]
    return
  end function descendentNodeConstructorParameters

  function descendentNodeConstructorInternal(timeDescendent,nodePropertyExtractor_) result(self)
    !!{
    Internal constructor for the ``descendentNode'' node property extractor class.
    !!}
    implicit none
    type            (nodePropertyExtractorDescendentNode)                        :: self
    double precision                                     , intent(in   )         :: timeDescendent
    class           (nodePropertyExtractorScalar        ), intent(in   ), target :: nodePropertyExtractor_
    !![
    <constructorAssign variables="timeDescendent, *nodePropertyExtractor_"/>
    !!]

    return
  end function descendentNodeConstructorInternal
  
  subroutine descendentNodeDestructor(self)
    !!{
    Destructor for  the ``descendentNode'' node property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorDescendentNode), intent(inout) :: self

    !![
    <objectDestructor name="self%nodePropertyExtractor_"/>
    !!]
    return
  end subroutine descendentNodeDestructor
  
  double precision function descendentNodeExtract(self,node,instance)
    !!{
    Implement a descendentNode output analysis.
    !!}
    use :: Galacticus_Error    , only : Galacticus_Error_Report
    use :: Galacticus_Nodes    , only : nodeComponentBasic
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    class           (nodePropertyExtractorDescendentNode), intent(inout)           :: self
    type            (treeNode                           ), intent(inout), target   :: node
    type            (multiCounter                       ), intent(inout), optional :: instance
    type            (treeNode                           ), pointer                 :: nodeDescendent
    class           (nodeComponentBasic                 ), pointer                 :: basicDescendent
    double precision                                     , parameter               :: timeTolerance  =1.0d-4

    descendentNodeExtract=-huge(0.0d0)
    nodeDescendent => node%parent
    do while (associated(nodeDescendent))
       basicDescendent => nodeDescendent%basic()
       if (Values_Agree(basicDescendent%time(),self%timeDescendent,absTol=timeTolerance)) then
          descendentNodeExtract=self%nodePropertyExtractor_%extract(nodeDescendent,instance)
          return
       end if
       nodeDescendent => nodeDescendent%parent
    end do
    call Galacticus_Error_Report('failed to find descendent node'//{introspection:location})
    return
  end function descendentNodeExtract

  integer function descendentNodeType(self)
    !!{
    Return the type of the halo mass property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorDescendentNode), intent(inout) :: self

    descendentNodeType=self%nodePropertyExtractor_%type()
    return
  end function descendentNodeType

  function descendentNodeName(self)
    !!{
    Return the name of the descendentNode property.
    !!}
    use :: String_Handling, only : String_Upper_Case_First
    implicit none
    type (varying_string                     )                :: descendentNodeName
    class(nodePropertyExtractorDescendentNode), intent(inout) :: self

    descendentNodeName=var_str('descendent')//String_Upper_Case_First(char(self%nodePropertyExtractor_%name()))
    return
  end function descendentNodeName

  function descendentNodeDescription(self)
    !!{
    Return a description of the descendentNode property.
    !!}
    implicit none
    type (varying_string                     )                :: descendentNodeDescription
    class(nodePropertyExtractorDescendentNode), intent(inout) :: self

    descendentNodeDescription=self%nodePropertyExtractor_%description()//' (of a descendent node)'
    return
  end function descendentNodeDescription

  double precision function descendentNodeUnitsInSI(self)
    !!{
    Return the units of the descendentNode property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorDescendentNode), intent(inout) :: self

    descendentNodeUnitsInSI=self%nodePropertyExtractor_%unitsInSI()
    return
  end function descendentNodeUnitsInSI
