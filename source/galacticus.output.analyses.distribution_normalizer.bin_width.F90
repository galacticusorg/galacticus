!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
  Contains a module which implements a bin width output analysis distribution normalizer class.
  !!}

  !![
  <outputAnalysisDistributionNormalizer name="outputAnalysisDistributionNormalizerBinWidth">
   <description>A bin width output analysis distribution normalizer class.</description>
  </outputAnalysisDistributionNormalizer>
  !!]
  type, extends(outputAnalysisDistributionNormalizerClass) :: outputAnalysisDistributionNormalizerBinWidth
     !!{
     A bin width output distribution normalizer class.
     !!}
     private
   contains
     procedure :: normalize => binWidthNormalize
  end type outputAnalysisDistributionNormalizerBinWidth

  interface outputAnalysisDistributionNormalizerBinWidth
     !!{
     Constructors for the ``binWidth'' output analysis distribution normalizer class.
     !!}
     module procedure binWidthConstructorParameters
  end interface outputAnalysisDistributionNormalizerBinWidth

contains

  function binWidthConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``binWidth'' output analysis distribution normalizer class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(outputAnalysisDistributionNormalizerBinWidth)                :: self
    type(inputParameters                             ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=outputAnalysisDistributionNormalizerBinWidth()
    return
  end function binWidthConstructorParameters

  subroutine binWidthNormalize(self,distribution,covariance,propertyValueMinimum,propertyValueMaximum)
    !!{
    Implement a bin width output analysis distribution normalizer.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (outputAnalysisDistributionNormalizerBinWidth), intent(inout)                 :: self
    double precision                                              , intent(inout), dimension(:  ) :: distribution
    double precision                                              , intent(inout), dimension(:,:) :: covariance
    double precision                                              , intent(in   ), dimension(:  ) :: propertyValueMinimum, propertyValueMaximum
    integer         (c_size_t                                    )                                :: i
    !$GLC attributes unused :: self

    distribution=+distribution            &
         &        /(                      &
         &          +propertyValueMaximum &
         &          -propertyValueMinimum &
         &        )
    forall(i=1:size(propertyValueMinimum))
       covariance(:,i)=+covariance(:,i)           &
         &             /(                         &
         &               +propertyValueMaximum(i) &
         &               -propertyValueMinimum(i) &
         &              )                         &
         &             /(                         &
         &               +propertyValueMaximum    &
         &               -propertyValueMinimum    &
         &              )
    end forall
    return
  end subroutine binWidthNormalize
