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

!!{
Implements a polynomial systematic shift output analysis property operator class.
!!}

  !![
  <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorSystmtcPolynomial">
   <description>
    A polynomial systematic shift output analysis property operator class. This operator allows for a systematic shift in
    properties (to account for systematic uncertainties in the observational analysis) using a simple model. Specifically,
    properties are mapped by this model as follows \begin{equation} \log_\mathrm{10} x \rightarrow \log_{10} x + \sum_{i=0}^N
    \alpha_i \log^i_{10}(x/x_0), \end{equation} where $x_0=${\normalfont \ttfamily [zeroPoint]} is a zero-point, and the
    coefficients $\alpha_{i=1\ldots N}=${\normalfont \ttfamily [coefficient]} are specified by input parameters.
   </description>
  </outputAnalysisPropertyOperator>
  !!]
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorSystmtcPolynomial
     !!{
     A polynomial systematic shift output property operator class.
     !!}
     private
     double precision                            :: zeroPoint
     double precision, allocatable, dimension(:) :: coefficient
   contains
     procedure :: operate => systmtcPolynomialOperate
  end type outputAnalysisPropertyOperatorSystmtcPolynomial

  interface outputAnalysisPropertyOperatorSystmtcPolynomial
     !!{
     Constructors for the \refClass{outputAnalysisPropertyOperatorSystmtcPolynomial} output analysis class.
     !!}
     module procedure systmtcPolynomialConstructorParameters
     module procedure systmtcPolynomialConstructorInternal
  end interface outputAnalysisPropertyOperatorSystmtcPolynomial

contains

  function systmtcPolynomialConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisPropertyOperatorSystmtcPolynomial} output analysis property operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial)                              :: self
    type            (inputParameters                                ), intent(inout)               :: parameters
    double precision                                                                               :: zeroPoint
    double precision                                                 , allocatable  , dimension(:) :: coefficient

    ! Check and read parameters.
    allocate(coefficient(parameters%count('coefficient')))
    !![
    <inputParameter>
      <name>zeroPoint</name>
      <source>parameters</source>
      <variable>zeroPoint</variable>
      <description>The zero-point of the property value used in the polynomial systematic offset property operator class.</description>
    </inputParameter>
    <inputParameter>
      <name>coefficient</name>
      <source>parameters</source>
      <variable>coefficient</variable>
      <description>The coefficients in the polynomial systematic offset property operator class.</description>
    </inputParameter>
    !!]
    ! Construct the object.
    self=outputAnalysisPropertyOperatorSystmtcPolynomial(zeroPoint,coefficient)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function systmtcPolynomialConstructorParameters

  function systmtcPolynomialConstructorInternal(zeroPoint,coefficient) result(self)
    !!{
    Internal constructor for the \refClass{outputAnalysisPropertyOperatorSystmtcPolynomial} output analysis distribution operator class.
    !!}
    implicit none
    type            (outputAnalysisPropertyOperatorSystmtcPolynomial)                              :: self
    double precision                                                 , intent(in   )               :: zeroPoint
    double precision                                                 , intent(in   ), dimension(:) :: coefficient
    !![
    <constructorAssign variables="zeroPoint, coefficient"/>
    !!]

    return
  end function systmtcPolynomialConstructorInternal

  double precision function systmtcPolynomialOperate(self,propertyValue,node,propertyType,outputIndex)
    !!{
    Implement an systmtcPolynomial output analysis property operator.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (outputAnalysisPropertyOperatorSystmtcPolynomial), intent(inout)           :: self
    double precision                                                 , intent(in   )           :: propertyValue
    type            (treeNode                                       ), intent(inout), optional :: node
    type            (enumerationOutputAnalysisPropertyTypeType      ), intent(inout), optional :: propertyType
    integer         (c_size_t                                       ), intent(in   ), optional :: outputIndex
    integer                                                                                    :: i
    !$GLC attributes unused :: outputIndex, propertyType, node

    systmtcPolynomialOperate=propertyValue
    ! Do not attempt to modify out-of-range values.
    if     (                              &
         &   propertyValue > -huge(0.0d0) &
         &  .and.                         &
         &   propertyValue < +huge(0.0d0) &
         & ) then
       do i=1,size(self%coefficient)
          systmtcPolynomialOperate=+systmtcPolynomialOperate &
               &                   +self%coefficient(i)      &
               &                   *(                        &
               &                     +     propertyValue     & 
               &                     -self%zeroPoint         &
               &                    )**(i-1)
       end do
    end if
    return
  end function systmtcPolynomialOperate
