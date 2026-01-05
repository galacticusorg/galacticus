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
  Implements an identity output analysis distribution normalizer class.
  !!}

  !![
  <outputAnalysisDistributionNormalizer name="outputAnalysisDistributionNormalizerIdentity">
   <description>An identity output analysis distribution normalizer class.</description>
  </outputAnalysisDistributionNormalizer>
  !!]
  type, extends(outputAnalysisDistributionNormalizerClass) :: outputAnalysisDistributionNormalizerIdentity
     !!{
     An identity output distribution normalizer class.
     !!}
     private
   contains
     procedure :: normalize => identityNormalize
  end type outputAnalysisDistributionNormalizerIdentity

  interface outputAnalysisDistributionNormalizerIdentity
     !!{
     Constructors for the \refClass{outputAnalysisDistributionNormalizerIdentity} output analysis distribution normalizer class.
     !!}
     module procedure identityConstructorParameters
  end interface outputAnalysisDistributionNormalizerIdentity

contains

  function identityConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisDistributionNormalizerIdentity} output analysis distribution normalizer class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(outputAnalysisDistributionNormalizerIdentity)                :: self
    type(inputParameters                             ), intent(inout) :: parameters

    self=outputAnalysisDistributionNormalizerIdentity()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function identityConstructorParameters

  subroutine identityNormalize(self,distribution,covariance,propertyValueMinimum,propertyValueMaximum)
    !!{
    Implement a bin width output analysis distribution normalizer.
    !!}
    implicit none
    class           (outputAnalysisDistributionNormalizerIdentity), intent(inout)                           :: self
    double precision                                              , intent(inout), dimension(:  ), optional :: distribution
    double precision                                              , intent(inout), dimension(:,:), optional :: covariance
    double precision                                              , intent(in   ), dimension(:  )           :: propertyValueMinimum, propertyValueMaximum
    !$GLC attributes unused :: self, propertyValueMinimum, propertyValueMaximum, distribution, covariance

    ! Leave everything unchanged.
    return
  end subroutine identityNormalize
