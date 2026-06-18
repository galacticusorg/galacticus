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

!!{RST
Implements a min-max analysis property operator class.
!!}

  !![
  <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorMinMax" docformat="rst">
   <description>
   An output analysis property operator that clamps a galaxy property value to lie within a specified range, replacing values below ``thresholdMinimum`` with the minimum threshold and values above ``thresholdMaximum`` with the maximum threshold.
   </description>
  </outputAnalysisPropertyOperator>
  !!]
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorMinMax
     !!{RST
     A high-pass filter property operator class.
     !!}
     private
     double precision :: thresholdMinimum, thresholdMaximum
   contains
     procedure :: operate => minMaxOperate
  end type outputAnalysisPropertyOperatorMinMax

  interface outputAnalysisPropertyOperatorMinMax
     !!{RST
     Constructors for the :galacticus-class:`outputAnalysisPropertyOperatorMinMax` output analysis property operator class.
     !!}
     module procedure minMaxConstructorParameters
     module procedure minMaxConstructorInternal
  end interface outputAnalysisPropertyOperatorMinMax

contains

  function minMaxConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`outputAnalysisPropertyOperatorMinMax` output analysis property operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(outputAnalysisPropertyOperatorMinMax)                :: self
    type(inputParameters                     ), intent(inout) :: parameters
    double precision                                          :: thresholdMinimum, thresholdMaximum

    ! Check and read parameters.
    !![
    <inputParameter docformat="rst">
      <name>thresholdMinimum</name>
      <source>parameters</source>
      <description>
      Minimum threshold for the min-max property operator.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>thresholdMaximum</name>
      <source>parameters</source>
      <description>
      Maximum threshold for the min-max property operator.
      </description>
    </inputParameter>
    !!]
    self=outputAnalysisPropertyOperatorMinMax(thresholdMinimum,thresholdMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function minMaxConstructorParameters

  function minMaxConstructorInternal(thresholdMinimum,thresholdMaximum) result (self)
    !!{RST
    Internal constructor for the :galacticus-class:`outputAnalysisPropertyOperatorMinMax` output analysis property operator class.
    !!}
    implicit none
    type            (outputAnalysisPropertyOperatorMinMax)                :: self
    double precision                                      , intent(in   ) :: thresholdMinimum, thresholdMaximum
    !![
    <constructorAssign variables="thresholdMinimum, thresholdMaximum"/>
    !!]

    return
  end function minMaxConstructorInternal

  double precision function minMaxOperate(self,propertyValue,node,propertyType,outputIndex)
    !!{RST
    Implement an minMax output analysis property operator.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (outputAnalysisPropertyOperatorMinMax     ), intent(inout)           :: self
    double precision                                           , intent(in   )           :: propertyValue
    type            (treeNode                                 ), intent(inout), optional :: node
    type            (enumerationOutputAnalysisPropertyTypeType), intent(inout), optional :: propertyType
    integer         (c_size_t                                 ), intent(in   ), optional :: outputIndex
    !$GLC attributes unused :: propertyType, outputIndex, node

    minMaxOperate=min(max(propertyValue,self%thresholdMinimum),self%thresholdMaximum)
    return
  end function minMaxOperate
