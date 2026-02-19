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
  Implements a $\log_{10}\rightarrow \log$ output analysis distribution normalizer class.
  !!}

  !![
  <outputAnalysisDistributionNormalizer name="outputAnalysisDistributionNormalizerLog10ToLog">
   <description>A $\log_{10}\rightarrow \log$ output analysis distribution normalizer class.</description>
  </outputAnalysisDistributionNormalizer>
  !!]
  type, extends(outputAnalysisDistributionNormalizerClass) :: outputAnalysisDistributionNormalizerLog10ToLog
     !!{
     A $\log_{10}\rightarrow \log$ output distribution normalizer class.
     !!}
     private
   contains
     procedure :: normalize => log10ToLogNormalize
  end type outputAnalysisDistributionNormalizerLog10ToLog

  interface outputAnalysisDistributionNormalizerLog10ToLog
     !!{
     Constructors for the \refClass{outputAnalysisDistributionNormalizerLog10ToLog} output analysis distribution normalizer class.
     !!}
     module procedure log10ToLogConstructorParameters
  end interface outputAnalysisDistributionNormalizerLog10ToLog

contains

  function log10ToLogConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisDistributionNormalizerLog10ToLog} output analysis distribution normalizer class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(outputAnalysisDistributionNormalizerLog10ToLog)                :: self
    type(inputParameters                               ), intent(inout) :: parameters

    self=outputAnalysisDistributionNormalizerLog10ToLog()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function log10ToLogConstructorParameters

  subroutine log10ToLogNormalize(self,distribution,covariance,propertyValueMinimum,propertyValueMaximum)
    !!{
    Implement a bin width output analysis distribution normalizer.
    !!}
    implicit none
    class           (outputAnalysisDistributionNormalizerLog10ToLog), intent(inout)                           :: self
    double precision                                                , intent(inout), dimension(:  ), optional :: distribution
    double precision                                                , intent(inout), dimension(:,:), optional :: covariance
    double precision                                                , intent(in   ), dimension(:  )           :: propertyValueMinimum, propertyValueMaximum
    !$GLC attributes unused :: self, propertyValueMinimum, propertyValueMaximum

    if (present(distribution)) distribution=distribution/log(10.0d0)
    if (present(covariance  )) covariance  =covariance  /log(10.0d0)**2
    return
  end subroutine log10ToLogNormalize
