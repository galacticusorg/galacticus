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
Implements an log10 output analysis property operator class.
!!}

  !![
  <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorLog10">
   <description>An log10 output analysis property operator class.</description>
  </outputAnalysisPropertyOperator>
  !!]
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorLog10
     !!{
     An log10 output property operator class.
     !!}
     private
   contains
     procedure :: operate => log10Operate
  end type outputAnalysisPropertyOperatorLog10

  interface outputAnalysisPropertyOperatorLog10
     !!{
     Constructors for the {\normalfont \ttfamily log10} output analysis class.
     !!}
     module procedure log10ConstructorParameters
  end interface outputAnalysisPropertyOperatorLog10

contains

  function log10ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily log10} output analysis property operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(outputAnalysisPropertyOperatorLog10)                :: self
    type(inputParameters                    ), intent(inout) :: parameters

    self=outputAnalysisPropertyOperatorLog10()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function log10ConstructorParameters

  double precision function log10Operate(self,propertyValue,node,propertyType,outputIndex)
    !!{
    Implement an log10 output analysis property operator.
    !!}
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    use            :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear, outputAnalysisPropertyTypeLog10, outputAnalysisPropertyTypeUnknown
    implicit none
    class           (outputAnalysisPropertyOperatorLog10      ), intent(inout)           :: self
    double precision                                           , intent(in   )           :: propertyValue
    type            (treeNode                                 ), intent(inout), optional :: node
    type            (enumerationOutputAnalysisPropertyTypeType), intent(inout), optional :: propertyType
    integer         (c_size_t                                 ), intent(in   ), optional :: outputIndex
    !$GLC attributes unused :: self, outputIndex, node

    if (propertyValue > 0.0d0) then
       log10Operate=log10(propertyValue)
    else
       log10Operate=-huge(1.0d0)
    end if
    ! Change the property type.
    if (present(propertyType)) then
       if (propertyType == outputAnalysisPropertyTypeLinear) then
          propertyType=outputAnalysisPropertyTypeLog10
       else
          propertyType=outputAnalysisPropertyTypeUnknown
       end if
    end if
    return
  end function log10Operate
