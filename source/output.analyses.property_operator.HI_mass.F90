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
Implements a conversion of ISM mass to HI mass analysis property operator class.
!!}

  use :: Output_Analysis_Molecular_Ratios, only : outputAnalysisMolecularRatioClass

  !![
  <outputAnalysisPropertyOperator name="outputAnalysisPropertyOperatorHIMass">
   <description>A conversion of ISM mass to HI mass analysis property operator class.</description>
  </outputAnalysisPropertyOperator>
  !!]
  type, extends(outputAnalysisPropertyOperatorClass) :: outputAnalysisPropertyOperatorHIMass
     !!{
     A conversion of ISM mass to HI mass property operator class.
     !!}
     private
     class(outputAnalysisMolecularRatioClass), pointer :: outputAnalysisMolecularRatio_ => null()
   contains
     final     ::            hiMassDestructor
     procedure :: operate => hiMassOperate
  end type outputAnalysisPropertyOperatorHIMass

  interface outputAnalysisPropertyOperatorHIMass
     !!{
     Constructors for the \refClass{outputAnalysisPropertyOperatorHIMass} output analysis class.
     !!}
     module procedure hiMassConstructorParameters
     module procedure hiMassConstructorInternal
  end interface outputAnalysisPropertyOperatorHIMass

contains

  function hiMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisPropertyOperatorHIMass} output analysis property operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (outputAnalysisPropertyOperatorHIMass)                :: self
    type (inputParameters                     ), intent(inout) :: parameters
    class(outputAnalysisMolecularRatioClass   ), pointer       :: outputAnalysisMolecularRatio_

    ! Check and read parameters.
    !![
    <objectBuilder class="outputAnalysisMolecularRatio" name="outputAnalysisMolecularRatio_" source="parameters" />
    !!]
    self=outputAnalysisPropertyOperatorHIMass(outputAnalysisMolecularRatio_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputAnalysisMolecularRatio_"/>
    !!]
    return
  end function hiMassConstructorParameters

  function hiMassConstructorInternal(outputAnalysisMolecularRatio_) result (self)
    !!{
    Internal constructor for the \refClass{outputAnalysisPropertyOperatorHIMass} output analysis distribution operator class.
    !!}
    implicit none
    type (outputAnalysisPropertyOperatorHIMass)                        :: self
    class(outputAnalysisMolecularRatioClass   ), intent(in   ), target :: outputAnalysisMolecularRatio_
    !![
    <constructorAssign variables="*outputAnalysisMolecularRatio_"/>
    !!]

    return
  end function hiMassConstructorInternal

  subroutine hiMassDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisPropertyOperatorHIMass} output analysis distribution operator class.
    !!}
    implicit none
    type (outputAnalysisPropertyOperatorHIMass), intent(inout) :: self

    !![
    <objectDestructor name="self%outputAnalysisMolecularRatio_"/>
    !!]
  return
  end subroutine hiMassDestructor

  double precision function hiMassOperate(self,propertyValue,node,propertyType,outputIndex)
    !!{
    Implement an hiMass output analysis property operator.
    !!}
    use            :: Error        , only : Error_Report
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (outputAnalysisPropertyOperatorHIMass     ), intent(inout)           :: self
    double precision                                           , intent(in   )           :: propertyValue
    type            (treeNode                                 ), intent(inout), optional :: node
    type            (enumerationOutputAnalysisPropertyTypeType), intent(inout), optional :: propertyType
    integer         (c_size_t                                 ), intent(in   ), optional :: outputIndex
    double precision                                                                     :: ratioHydrogenMolecularHydrogenNeutral
    !$GLC attributes unused :: propertyType, outputIndex
    
    if (.not.present(node)) call Error_Report('node must be provided'//{introspection:location})
    ratioHydrogenMolecularHydrogenNeutral=+self%outputAnalysisMolecularRatio_                                 %ratio               (propertyValue,node) &
         &                                *10.0d0**(                                                                                                    &
         &                                          +node%hostTree                     %randomNumberGenerator_%standardNormalSample(                  ) &
         &                                          *self%outputAnalysisMolecularRatio_                       %ratioScatter        (propertyValue,node) &
         &                                         )
    hiMassOperate                        =+propertyValue                                                                                                &
         &                                /(                                                                                                            &
         &                                  +1.0d0                                                                                                      &
         &                                  +ratioHydrogenMolecularHydrogenNeutral                                                                      &
         &                                 )
    return
  end function hiMassOperate

