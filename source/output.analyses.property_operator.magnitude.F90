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
Implements an output analysis property operator class which converts luminosity to absolute magnitude.
!!}

  !![
  <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorMagnitude">
   <description>An output analysis property operator class.</description>
  </outputAnalysisPropertyOperator>
  !!]
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorMagnitude
     !!{
     An magnitude output property operator class.
     !!}
     private
   contains
     procedure :: operate => magnitudeOperate
  end type outputAnalysisPropertyOperatorMagnitude

  interface outputAnalysisPropertyOperatorMagnitude
     !!{
     Constructors for the ``magnitude'' output analysis class.
     !!}
     module procedure magnitudeConstructorParameters
  end interface outputAnalysisPropertyOperatorMagnitude

contains

  function magnitudeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``magnitude'' output analysis property operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(outputAnalysisPropertyOperatorMagnitude)                :: self
    type(inputParameters                        ), intent(inout) :: parameters

    self=outputAnalysisPropertyOperatorMagnitude()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function magnitudeConstructorParameters

  double precision function magnitudeOperate(self,propertyValue,node,propertyType,outputIndex)
    !!{
    Implement an magnitude output analysis property operator.
    !!}
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    use            :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear, outputAnalysisPropertyTypeMagnitude, outputAnalysisPropertyTypeUnknown
    implicit none
    class           (outputAnalysisPropertyOperatorMagnitude  ), intent(inout)           :: self
    double precision                                           , intent(in   )           :: propertyValue
    type            (treeNode                                 ), intent(inout), optional :: node
    type            (enumerationOutputAnalysisPropertyTypeType), intent(inout), optional :: propertyType
    integer         (c_size_t                                 ), intent(in   ), optional :: outputIndex
    !$GLC attributes unused :: self, outputIndex, node

    if (propertyValue > 0.0d0) then
       magnitudeOperate=-2.5d0*log10(propertyValue)
    else
       magnitudeOperate=+      huge (        0.0d0)
    end if
    ! Change the property type.
    if (present(propertyType)) then
       if (propertyType == outputAnalysisPropertyTypeLinear) then
          propertyType=outputAnalysisPropertyTypeMagnitude
       else
          propertyType=outputAnalysisPropertyTypeUnknown
       end if
    end if
    return
  end function magnitudeOperate
