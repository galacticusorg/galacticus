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
Implements a ratio output analysis property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRatio">
   <description>A ratio output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRatio
     !!{
     A ratio extractor output analysis class. This extractor extracts two other properties and takes their ratio.
     !!}
     private
     class           (nodePropertyExtractorClass), pointer :: propertyNumerator_ => null(), propertyDenominator_ => null()
     type            (varying_string            )          :: name_                       , description_
     double precision                                      :: unitsInSI_
   contains
     final     ::                ratioDestructor
     procedure :: extract     => ratioExtract
     procedure :: name        => ratioName
     procedure :: description => ratioDescription
     procedure :: unitsInSI   => ratioUnitsInSI
  end type nodePropertyExtractorRatio

  interface nodePropertyExtractorRatio
     !!{
     Constructors for the ``ratio'' output analysis class.
     !!}
     module procedure ratioConstructorParameters
     module procedure ratioConstructorInternal
  end interface nodePropertyExtractorRatio

contains

  function ratioConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``ratio'' output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorRatio)                :: self
    type (inputParameters           ), intent(inout) :: parameters
    class(nodePropertyExtractorClass), pointer       :: propertyNumerator_, propertyDenominator_
    type (varying_string            )                :: name              , description

    !![
    <inputParameter>
      <name>name</name>
      <source>parameters</source>
      <description>The name of this property.</description>
    </inputParameter>
    <inputParameter>
      <name>description</name>
      <source>parameters</source>
      <description>A description of this property.</description>
    </inputParameter>
    <objectBuilder class="nodePropertyExtractor" name="propertyNumerator_"   parameterName="nodePropertyExtractorNumerator"   source="parameters"/>
    <objectBuilder class="nodePropertyExtractor" name="propertyDenominator_" parameterName="nodePropertyExtractorDenominator" source="parameters"/>
    !!]
    self=nodePropertyExtractorRatio(char(name),char(description),propertyNumerator_,propertyDenominator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="propertyNumerator_"  />
    <objectDestructor name="propertyDenominator_"/>
    !!]
    return
  end function ratioConstructorParameters

  function ratioConstructorInternal(name,description,propertyNumerator_,propertyDenominator_) result(self)
    !!{
    Internal constructor for the ``ratio'' output analysis property extractor class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type     (nodePropertyExtractorRatio)                        :: self
    class    (nodePropertyExtractorClass), intent(inout), target :: propertyNumerator_, propertyDenominator_
    character(len=*                     ), intent(in   )         :: name              , description
    !![
    <constructorAssign variables="*propertyNumerator_, *propertyDenominator_"/>
    !!]

    self%name_       =name
    self%description_=description
    select type (propertyNumerator_  )
    class is (nodePropertyExtractorScalar)
       self%unitsInSI_=+propertyNumerator_  %unitsInSI()
    class default
       call Error_Report('numerator property must be a scalar'  //{introspection:location})
    end select
    select type (propertyDenominator_)
    class is (nodePropertyExtractorScalar)
       self%unitsInSI_=+self%unitsInSI_                  &
            &          /propertyDenominator_%unitsInSI()
    class default
       call Error_Report('denominator property must be a scalar'//{introspection:location})
    end select
    return
  end function ratioConstructorInternal

  subroutine ratioDestructor(self)
    !!{
    Destructor for the ``ratio'' output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorRatio), intent(inout) :: self

    !![
    <objectDestructor name="self%propertyNumerator_"  />
    <objectDestructor name="self%propertyDenominator_"/>
    !!]
    return
  end subroutine ratioDestructor

  double precision function ratioExtract(self,node,instance)
    !!{
    Implement a ratio output analysis.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (nodePropertyExtractorRatio), intent(inout), target   :: self
    type            (treeNode                  ), intent(inout), target   :: node
    type            (multiCounter              ), intent(inout), optional :: instance
    double precision                                                      :: numerator, denominator
    !$GLC attributes unused :: instance

    select type (propertyNumerator_   => self%propertyNumerator_  )
    class is (nodePropertyExtractorScalar)
       numerator   =propertyNumerator_  %extract(node,instance)
    class default
       numerator   =huge(0.0d0)
    end select
    select type (propertyDenominator_ => self%propertyDenominator_)
    class is (nodePropertyExtractorScalar)
       denominator = propertyDenominator_%extract(node,instance)
    class default
       denominator =     0.0d0
    end select
    if (denominator == 0.0d0) then
       ratioExtract=+0.0d0
       if (numerator /= 0.0d0) call Error_Report('denominator is zero'//{introspection:location})
    else
       ratioExtract=+numerator   &
            &       /denominator
    end if
    return
  end function ratioExtract

  function ratioName(self)
    !!{
    Return the name of this property.
    !!}
    implicit none
    type (varying_string            )                :: ratioName
    class(nodePropertyExtractorRatio), intent(inout) :: self

    ratioName=self%name_
    return
  end function ratioName

  function ratioDescription(self)
    !!{
    Return the description of this property.
    !!}
    implicit none
    type (varying_string            )                :: ratioDescription
    class(nodePropertyExtractorRatio), intent(inout) :: self

    ratioDescription=self%description_
    return
  end function ratioDescription

  double precision function ratioUnitsInSI(self)
    !!{
    Return the description of this property.
    !!}
    implicit none
    class(nodePropertyExtractorRatio), intent(inout) :: self

    ratioUnitsInSI=self%unitsInSI_
    return
  end function ratioUnitsInSI

