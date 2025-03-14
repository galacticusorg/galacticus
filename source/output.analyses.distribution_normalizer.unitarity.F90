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
  Implements a unitarity output analysis distribution normalizer class.
  !!}

  !![
  <outputAnalysisDistributionNormalizer name="outputAnalysisDistributionNormalizerUnitarity">
   <description>A unitarity output analysis distribution normalizer class.</description>
  </outputAnalysisDistributionNormalizer>
  !!]
  type, extends(outputAnalysisDistributionNormalizerClass) :: outputAnalysisDistributionNormalizerUnitarity
     !!{
     A unitarity output distribution normalizer class.
     !!}
     private
   contains
     procedure :: normalize => unitarityNormalize
  end type outputAnalysisDistributionNormalizerUnitarity

  interface outputAnalysisDistributionNormalizerUnitarity
     !!{
     Constructors for the ``unitarity'' output analysis distribution normalizer class.
     !!}
     module procedure unitarityConstructorParameters
  end interface outputAnalysisDistributionNormalizerUnitarity

contains

  function unitarityConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``unitarity'' output analysis distribution normalizer class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(outputAnalysisDistributionNormalizerUnitarity)                :: self
    type(inputParameters                              ), intent(inout) :: parameters

    self=outputAnalysisDistributionNormalizerUnitarity()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function unitarityConstructorParameters

  subroutine unitarityNormalize(self,distribution,covariance,propertyValueMinimum,propertyValueMaximum)
    !!{
    Implement a unitarity output analysis distribution normalizer.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (outputAnalysisDistributionNormalizerUnitarity), intent(inout)                           :: self
    double precision                                               , intent(inout), dimension(:  ), optional :: distribution
    double precision                                               , intent(inout), dimension(:,:), optional :: covariance
    double precision                                               , intent(in   ), dimension(:  )           :: propertyValueMinimum, propertyValueMaximum
    double precision                                                                                         :: distributionSum
    !$GLC attributes unused :: self, propertyValueMinimum, propertyValueMaximum

    if (.not.present(distribution)) call Error_Report('"distribution" required for this class'//{introspection:location})
    distributionSum=+sum(distribution   )
    if (distributionSum /= 0.0d0) then
       distribution     =+distribution       &
            &            /distributionSum
       if (present(covariance))              &
            & covariance=+covariance         &
            &            /distributionSum**2
    end if
    return
  end subroutine unitarityNormalize
