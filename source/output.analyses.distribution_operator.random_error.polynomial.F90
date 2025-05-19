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
  Implements a random error output analysis distribution operator class with an error magnitude that is
  a polynomial function of the property value.
  !!}

  !![
  <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorRandomErrorPlynml">
   <description>A random error output analysis distribution operator class.</description>
  </outputAnalysisDistributionOperator>
  !!]
  type, extends(outputAnalysisDistributionOperatorRandomError) :: outputAnalysisDistributionOperatorRandomErrorPlynml
     !!{
     A random error output distribution operator class which has an error magnitude that is a polynomial function of the
     property value.
     !!}
     private
     double precision                            :: zeroPoint, errorMinimum, &
          &                                         errorMaximum
     double precision, allocatable, dimension(:) :: coefficient
   contains
     procedure :: rootVariance => randomErrorPolynomialRootVariance
  end type outputAnalysisDistributionOperatorRandomErrorPlynml

  interface outputAnalysisDistributionOperatorRandomErrorPlynml
     !!{
     Constructors for the {\normalfont \ttfamily randomErrorPolynomial} output analysis distribution operator class.
     !!}
     module procedure randomErrorPolynomialConstructorParameters
     module procedure randomErrorPolynomialConstructorInternal
  end interface outputAnalysisDistributionOperatorRandomErrorPlynml

contains

  function randomErrorPolynomialConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily randomErrorPolynomial} output analysis distribution operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml)                              :: self
    type            (inputParameters                                    ), intent(inout)               :: parameters
    double precision                                                                                   :: zeroPoint   , errorMinimum, &
         &                                                                                                errorMaximum
    double precision                                                     , allocatable  , dimension(:) :: coefficient

    ! Check and read parameters.
    allocate(coefficient(parameters%count('coefficient')))
    !![
    <inputParameter>
      <name>zeroPoint</name>
      <source>parameters</source>
      <variable>zeroPoint</variable>
      <description>The zero-point of the property value used in the polynomial random error distribution class.</description>
    </inputParameter>
    <inputParameter>
      <name>coefficient</name>
      <source>parameters</source>
      <variable>coefficient</variable>
      <description>The coefficients in the polynomial random error distribution class.</description>
    </inputParameter>
    <inputParameter>
      <name>errorMinimum</name>
      <source>parameters</source>
      <variable>errorMinimum</variable>
      <description>The minimum error in the polynomial random error distribution class.</description>
    </inputParameter>
    <inputParameter>
      <name>errorMaximum</name>
      <source>parameters</source>
      <variable>errorMaximum</variable>
      <description>The maximum error in the polynomial random error distribution class.</description>
    </inputParameter>
    !!]
    ! Construct the object.
    self=outputAnalysisDistributionOperatorRandomErrorPlynml(errorMinimum,errorMaximum,zeroPoint,coefficient)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function randomErrorPolynomialConstructorParameters

  function randomErrorPolynomialConstructorInternal(errorMinimum,errorMaximum,zeroPoint,coefficient) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily randomErrorPolynomial} output analysis distribution operator class.
    !!}
    implicit none
    type            (outputAnalysisDistributionOperatorRandomErrorPlynml)                              :: self
    double precision                                                     , intent(in   )               :: errorMinimum, errorMaximum, &
         &                                                                                                zeroPoint
    double precision                                                     , intent(in   ), dimension(:) :: coefficient
    !![
    <constructorAssign variables="errorMinimum, errorMaximum, zeroPoint, coefficient"/>
    !!]

    return
  end function randomErrorPolynomialConstructorInternal

  double precision function randomErrorPolynomialRootVariance(self,propertyValue,node)
    !!{
    Return the root-variance in the polynomial random error distribution operator.
    !!}
    implicit none
    class           (outputAnalysisDistributionOperatorRandomErrorPlynml), intent(inout) :: self
    double precision                                                     , intent(in   ) :: propertyValue
    type            (treeNode                                           ), intent(inout) :: node
    integer                                                                              :: i
    !$GLC attributes unused :: node

    randomErrorPolynomialRootVariance=0.0d0
    do i=1,size(self%coefficient)
       randomErrorPolynomialRootVariance=+randomErrorPolynomialRootVariance &
            &                            +self%coefficient(i)               &
            &                            *(                                 &
            &                              +     propertyValue              &
            &                              -self%zeroPoint                  &
            &                             )**(i-1)
    end do
    randomErrorPolynomialRootVariance=min(self%errorMaximum,max(self%errorMinimum,randomErrorPolynomialRootVariance))
    return
  end function randomErrorPolynomialRootVariance
