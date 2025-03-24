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
Implements a high-pass filter analysis property operator class.
!!}

  !![
  <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorFilterHighPass">
   <description>A high-pass filter analysis property operator class.</description>
  </outputAnalysisPropertyOperator>
  !!]
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorFilterHighPass
     !!{
     A high-pass filter property operator class.
     !!}
     private
     double precision :: filterThreshold, filterWidth
     logical          :: normalized
   contains
     procedure :: operate => filterHighPassOperate
  end type outputAnalysisPropertyOperatorFilterHighPass

  interface outputAnalysisPropertyOperatorFilterHighPass
     !!{
     Constructors for the {\normalfont \ttfamily filterHighPass} output analysis class.
     !!}
     module procedure filterHighPassConstructorParameters
     module procedure filterHighPassConstructorInternal
  end interface outputAnalysisPropertyOperatorFilterHighPass

contains

  function filterHighPassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily filterHighPass} output analysis property operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisPropertyOperatorFilterHighPass)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    double precision                                                              :: filterThreshold, filterWidth
    logical                                                                       :: normalized

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>filterThreshold</name>
      <source>parameters</source>
      <description>Threshold for the high-pass filter distribution operator.</description>
    </inputParameter>
    <inputParameter>
      <name>filterWidth</name>
      <defaultValue>0.0d0</defaultValue>
      <source>parameters</source>
      <description>The width of the filter (0 for a sharp transition; $&gt;0$ for a smoothed transition.</description>
    </inputParameter>
    <inputParameter>
      <name>normalized</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If true, the property value is set to the filter value, otherwise it is multiplied by it.</description>
    </inputParameter>
    !!]
    self=outputAnalysisPropertyOperatorFilterHighPass(filterThreshold,filterWidth,normalized)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function filterHighPassConstructorParameters

  function filterHighPassConstructorInternal(filterThreshold,filterWidth,normalized) result (self)
    !!{
    Internal constructor for the {\normalfont \ttfamily filterHighPass} output analysis distribution operator class.
    !!}
    implicit none
    type            (outputAnalysisPropertyOperatorFilterHighPass)                          :: self
    double precision                                              , intent(in   )           :: filterThreshold
    double precision                                              , intent(in   ), optional :: filterWidth
    logical                                                       , intent(in   ), optional :: normalized
    !![
    <constructorAssign variables="filterThreshold, filterWidth, normalized"/>
    !!]

    if (.not.present(filterWidth)) self%filterWidth=0.0d0
    if (.not.present(normalized )) self%normalized =.false.
    return
  end function filterHighPassConstructorInternal

  double precision function filterHighPassOperate(self,propertyValue,node,propertyType,outputIndex)
    !!{
    Implement an filterHighPass output analysis property operator.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (outputAnalysisPropertyOperatorFilterHighPass), intent(inout)           :: self
    double precision                                              , intent(in   )           :: propertyValue
    type            (treeNode                                    ), intent(inout), optional :: node
    type            (enumerationOutputAnalysisPropertyTypeType   ), intent(inout), optional :: propertyType
    integer         (c_size_t                                    ), intent(in   ), optional :: outputIndex
    double precision                                                                        :: filterValue
    !$GLC attributes unused :: propertyType, outputIndex, node

    if (self%filterWidth <= 0.0d0) then
       ! Sharp filter.
       if (propertyValue > self%filterThreshold) then
          filterValue=+1.0d0
       else
          filterValue=+0.0d0
       end if
    else
       ! Smooth filter using a sigmoid function
       filterValue   =+1.0d0                         &
            &         /(                             &
            &           +1.0d0                       &
            &           +exp(                        &
            &                -(                      &
            &                  +     propertyValue   &
            &                  -self%filterThreshold &
            &                 )                      &
            &                /  self%filterWidth     &
            &               )                        &
            &          )
    end if
    if (self%normalized) then
       filterHighPassOperate=+filterValue
    else
       filterHighPassOperate=+filterValue   &
            &                *propertyValue
    end if
    return
  end function filterHighPassOperate
