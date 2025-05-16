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
  Implements a property operator class which converts a metallicity, assumed to be a mass ratio of a
  given element to hydrogen, to $12+\log_{10}(\mathrm{N}/\mathrm{H})$ form.
  !!}

  !![
  <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorMetallicity12LogNH">
   <description>A property operator class which converts a metallicity, assumed to be a mass ratio of a given element to hydrogen, to $12+\log_{10}(\mathrm{N}/\mathrm{H})$ form.</description>
  </outputAnalysisPropertyOperator>
  !!]
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorMetallicity12LogNH
     !!{
     A metallicity property operator class which converts a metallicity, assumed to be a mass ratio of a given element to
     hydrogen, to $12+\log_{10}(\mathrm{N}/\mathrm{H})$ form.
     !!}
     private
     double precision :: massElement
   contains
     procedure :: operate => metallicity12LogNHOperate
  end type outputAnalysisPropertyOperatorMetallicity12LogNH

  interface outputAnalysisPropertyOperatorMetallicity12LogNH
     !!{
     Constructors for the \refClass{outputAnalysisPropertyOperatorMetallicity12LogNH} output analysis class.
     !!}
     module procedure metallicity12LogNHConstructorParameters
     module procedure metallicity12LogNHConstructorInternal
  end interface outputAnalysisPropertyOperatorMetallicity12LogNH

contains

  function metallicity12LogNHConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisPropertyOperatorMetallicity12LogNH} output analysis property operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisPropertyOperatorMetallicity12LogNH)                :: self
    type            (inputParameters                                 ), intent(inout) :: parameters
    double precision                                                                  :: massElement

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>massElement</name>
      <source>parameters</source>
      <description>The atomic mass of the element used to define metallicity.</description>
    </inputParameter>
    !!]
    self=outputAnalysisPropertyOperatorMetallicity12LogNH(massElement)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function metallicity12LogNHConstructorParameters

  function metallicity12LogNHConstructorInternal(massElement) result (self)
    !!{
    Internal constructor for the \refClass{outputAnalysisPropertyOperatorMetallicity12LogNH} output analysis distribution operator class.
    !!}
    implicit none
    type            (outputAnalysisPropertyOperatorMetallicity12LogNH)                :: self
    double precision                                                  , intent(in   ) :: massElement
    !![
    <constructorAssign variables="massElement"/>
    !!]

    return
  end function metallicity12LogNHConstructorInternal

  double precision function metallicity12LogNHOperate(self,propertyValue,node,propertyType,outputIndex)
    !!{
    Implement an metallicity output analysis property operator.
    !!}
    use, intrinsic :: ISO_C_Binding             , only : c_size_t
    use            :: Numerical_Constants_Atomic, only : atomicMassHydrogen
    implicit none
    class           (outputAnalysisPropertyOperatorMetallicity12LogNH), intent(inout)           :: self
    double precision                                                  , intent(in   )           :: propertyValue
    type            (treeNode                                        ), intent(inout), optional :: node
    type            (enumerationOutputAnalysisPropertyTypeType       ), intent(inout), optional :: propertyType
    integer         (c_size_t                                        ), intent(in   ), optional :: outputIndex
    double precision                                                                            :: ratioByNumber
    !$GLC attributes unused :: propertyType, outputIndex, node

    ratioByNumber=+propertyValue      &
            &     *atomicMassHydrogen &
            &     /self%massElement
    if (ratioByNumber > 0.0d0) then
       metallicity12LogNHOperate=+log10(ratioByNumber) &
            &                    +12.0d0
    else
       metallicity12LogNHOperate=-huge (+0.0d0       )
    end if
    return
  end function metallicity12LogNHOperate
