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
Implements an output analysis property extractor class that scalarizes one element from an array node property extractor.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorScalarizer">
   <description>An output analysis property extractor class that scalarizes one element from an array node property extractor.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorScalarizer
     !!{
     A property extractor output analysis class that scalarizes one element from an array node property extractor.
     !!}
     private
     integer                                      :: item                            , element
     class  (nodePropertyExtractorClass), pointer :: nodePropertyExtractor_ => null()
   contains
     final     ::                scalarizerDestructor
     procedure :: extract     => scalarizerExtract
     procedure :: name        => scalarizerName
     procedure :: description => scalarizerDescription
     procedure :: unitsInSI   => scalarizerUnitsInSI
  end type nodePropertyExtractorScalarizer

  interface nodePropertyExtractorScalarizer
     !!{
     Constructors for the \refClass{nodePropertyExtractorScalarizer} output analysis class.
     !!}
     module procedure scalarizerConstructorParameters
     module procedure scalarizerConstructorInternal
  end interface nodePropertyExtractorScalarizer

contains

  function scalarizerConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorScalarizer} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodePropertyExtractorScalarizer)                :: self
    type   (inputParameters                ), intent(inout) :: parameters
    class  (nodePropertyExtractorClass     ), pointer       :: nodePropertyExtractor_
    integer                                                 :: item                  , element

    !![
    <objectBuilder class="nodePropertyExtractor" name="nodePropertyExtractor_" source="parameters"/>
    <inputParameter>
      <name>element</name>
      <description>The element to scalarize from the array.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    select type (nodePropertyExtractor_)
    class is (nodePropertyExtractorArray)
       !![
       <inputParameter>
	 <name>item</name>
	 <description>The item to scalarize from the array.</description>
	 <source>parameters</source>
       </inputParameter>
       !!]
    class default
       ! "item" is not relevant for non-array extractors.
       item=-1
    end select
    self=nodePropertyExtractorScalarizer(item,element,nodePropertyExtractor_)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function scalarizerConstructorParameters

  function scalarizerConstructorInternal(item,element,nodePropertyExtractor_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorScalarizer} property extractor class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type   (nodePropertyExtractorScalarizer)                        :: self
    integer                                 , intent(in   )         :: item                  , element
    class  (nodePropertyExtractorClass     ), intent(in   ), target :: nodePropertyExtractor_
    !![
    <constructorAssign variables="item, element, *nodePropertyExtractor_"/>
    !!]

    select type (nodePropertyExtractor__ => self%nodePropertyExtractor_)
    class is (nodePropertyExtractorArray)
       ! This is as expected.
    class is (nodePropertyExtractorTuple)
       ! This is as expected.
    class default
       call Error_Report('class must be nodePropertyExtractorArray'//{introspection:location})
    end select
    return
  end function scalarizerConstructorInternal

  subroutine scalarizerDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorScalarizer} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorScalarizer), intent(inout) :: self

    !![
    <objectDestructor name="self%nodePropertyExtractor_"/>
    !!]
    return
  end subroutine scalarizerDestructor

  double precision function scalarizerExtract(self,node,instance)
    !!{
    Implement a scalarizer output analysis.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (nodePropertyExtractorScalarizer), intent(inout), target         :: self
    type            (treeNode                       ), intent(inout), target         :: node
    type            (multiCounter                   ), intent(inout), optional       :: instance
    class           (nodeComponentBasic             ), pointer                       :: basic
    double precision                                 , allocatable  , dimension(:,:) :: array
    double precision                                 , allocatable  , dimension(  :) :: tuple

    select type (nodePropertyExtractor__ => self%nodePropertyExtractor_)
    class is (nodePropertyExtractorArray)
       basic => node%basic()
       if (self%item    > nodePropertyExtractor__%size        (basic%time())) call Error_Report('item exceeds size of array'    //{introspection:location})
       if (self%element > nodePropertyExtractor__%elementCount(basic%time())) call Error_Report('element exceeds count of array'//{introspection:location})
       array            =nodePropertyExtractor__%extract(node     ,basic%time   (),instance)
       scalarizerExtract=array                          (self%item,self %element           )
    class is (nodePropertyExtractorTuple)
       basic => node%basic()
       if (self%element > nodePropertyExtractor__%elementCount(basic%time())) call Error_Report('element exceeds count of tuple'//{introspection:location})
       tuple            =nodePropertyExtractor__%extract(node     ,basic%time   (),instance)
       scalarizerExtract=tuple                          (          self %element           )
     class default
       scalarizerExtract=0.0d0
       call Error_Report('class must be nodePropertyExtractorArray'//{introspection:location})
    end select
    return
  end function scalarizerExtract

  function scalarizerName(self)
    !!{
    Return the name of the scalarizer property.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (varying_string                 )                              :: scalarizerName
    class(nodePropertyExtractorScalarizer), intent(inout)               :: self
    type (varying_string                 ), allocatable  , dimension(:) :: names
    
    select type (nodePropertyExtractor__ => self%nodePropertyExtractor_)
    class is (nodePropertyExtractorArray)
       call nodePropertyExtractor__%names(             names)
       scalarizerName=names(self%element)
    class is (nodePropertyExtractorTuple)
       call nodePropertyExtractor__%names(-huge(0.0d0),names)
       scalarizerName=names(self%element)
    class default
       call Error_Report('class must be nodePropertyExtractorArray'//{introspection:location})
    end select
    return
   end function scalarizerName

  function scalarizerDescription(self)
    !!{
    Return a description of the scalarizer property.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (varying_string                 )                              :: scalarizerDescription
    class(nodePropertyExtractorScalarizer), intent(inout)               :: self
    type (varying_string                 ), allocatable  , dimension(:) :: descriptions

    select type (nodePropertyExtractor__ => self%nodePropertyExtractor_)
    class is (nodePropertyExtractorArray)
       call nodePropertyExtractor__%descriptions(             descriptions)
       scalarizerDescription=descriptions(self%element)
    class is (nodePropertyExtractorTuple)
       call nodePropertyExtractor__%descriptions(-huge(0.0d0),descriptions)
       scalarizerDescription=descriptions(self%element)
    class default
       call Error_Report('class must be nodePropertyExtractorArray'//{introspection:location})
    end select
    return
  end function scalarizerDescription

  double precision function scalarizerUnitsInSI(self)
    !!{
    Return the units of the scalarizer property in the SI system.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (nodePropertyExtractorScalarizer), intent(inout)               :: self
    double precision                                 , allocatable  , dimension(:) :: unitsInSI

    select type (nodePropertyExtractor__ => self%nodePropertyExtractor_)
    class is (nodePropertyExtractorArray)
       unitsInSI          =nodePropertyExtractor__%unitsInSI(            )
       scalarizerUnitsInSI=unitsInSI                        (self%element)
    class is (nodePropertyExtractorTuple)
       unitsInSI          =nodePropertyExtractor__%unitsInSI(-huge(0.0d0))
       scalarizerUnitsInSI=unitsInSI                        (self%element)
    class default
       scalarizerUnitsInSI=0.0d0
       call Error_Report('class must be nodePropertyExtractorArray'//{introspection:location})
    end select
    return
  end function scalarizerUnitsInSI
