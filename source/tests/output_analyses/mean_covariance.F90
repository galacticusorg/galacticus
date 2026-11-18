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

!+    Contributions to this file made by: Claude.

!!{RST
Contains a program which tests the covariance of binned means computed by ``Output_Analysis_Mean_Covariance``.
!!}

program Test_Output_Analyses_Mean_Covariance
  !!{RST
  Tests ``Output_Analysis_Mean_Covariance``, which computes the covariance of the mean of a property in each bin of a
  ``meanFunction1D`` output analysis. Weights and covariances are accumulated from a small set of galaxies, as for the Poisson
  covariance model. The covariance is checked against a direct evaluation of the first-order Taylor series expression for
  well-scaled weights, and then shown to be unchanged when the weights in each bin are rescaled by factors down to
  :math:`10^{-150}` (the covariance of a ratio is invariant under such rescaling). Finally, a bin with a weight of
  :math:`10^{-165}`, whose square underflows to zero (as found in issue #1758), and an empty bin are checked to have zero
  covariance, leaving the covariance of all other bins unchanged.
  !!}
  use :: Display                  , only : displayVerbositySet            , verbosityLevelStandard
  use :: Error                    , only : Error_Handler_Register
  use :: Output_Analysis_Utilities, only : Output_Analysis_Mean_Covariance
  use :: Unit_Tests               , only : Assert                         , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  integer         , parameter                           :: countBins           =4, countGalaxies=5
  double precision, parameter                           :: toleranceRelative   =1.0d-10
  ! Distribution of each galaxy across the bins, its weight, and its property value.
  double precision, parameter, dimension(countBins,countGalaxies) :: distribution=reshape(                                     &
       &                                                                                  [                                    &
       &                                                                                   0.70d0, 0.20d0, 0.08d0, 0.02d0,     &
       &                                                                                   0.10d0, 0.60d0, 0.25d0, 0.05d0,     &
       &                                                                                   0.05d0, 0.15d0, 0.50d0, 0.30d0,     &
       &                                                                                   0.30d0, 0.30d0, 0.20d0, 0.20d0,     &
       &                                                                                   0.02d0, 0.08d0, 0.30d0, 0.60d0      &
       &                                                                                  ]                                  , &
       &                                                                                  [countBins,countGalaxies]            &
       &                                                                                 )
  double precision, parameter, dimension(          countGalaxies) :: weight      =[0.8d0, 1.3d0, 0.5d0, 1.1d0, 0.9d0]
  double precision, parameter, dimension(          countGalaxies) :: property    =[1.0d0, 2.5d0, 0.7d0, 3.1d0, 1.8d0]
  double precision           , dimension(countBins            ) :: valueUnweighted     , valueWeighted         , &
       &                                                           valueMean           , scale
  double precision           , dimension(countBins,countBins  ) :: covarianceUnweighted, covarianceWeighted    , &
       &                                                           covarianceCross     , covarianceMean        , &
       &                                                           covarianceReference , covarianceDirect
  double precision                                              :: toleranceAbsolute
  integer                                                       :: i                   , j

  call displayVerbositySet(verbosityLevelStandard)
  call Error_Handler_Register()
  call Unit_Tests_Begin_Group("Output analyses: covariance of binned means")
  ! Well-scaled weights: compare with direct evaluation of the Taylor series expression in the form
  ! [C^W U_i U_j + C^U W_i W_j - C^× (U_i W_j + U_j W_i)] / (U_i² U_j²).
  scale=1.0d0
  call accumulate(scale)
  covarianceReference=Output_Analysis_Mean_Covariance(valueUnweighted,valueWeighted,covarianceUnweighted,covarianceWeighted,covarianceCross)
  valueMean          =valueWeighted/valueUnweighted
  do i=1,countBins
     do j=1,countBins
        covarianceDirect(i,j)=+(                                                                                             &
             &                  +  covarianceWeighted(i,j)* valueUnweighted(i)              *valueUnweighted(j)              &
             &                  +covarianceUnweighted(i,j)*                   valueWeighted(i)                 *valueWeighted(j) &
             &                  -     covarianceCross(i,j)*(valueUnweighted(i)*valueWeighted(j)+valueUnweighted(j)*valueWeighted(i)) &
             &                 )                                                                                             &
             &                /valueUnweighted(i)**2                                                                         &
             &                /valueUnweighted(j)**2
     end do
     covarianceDirect(i,i)=max(covarianceDirect(i,i),0.0d0)
  end do
  toleranceAbsolute=1.0d-12*maxval(abs(covarianceDirect))
  call Assert('well-scaled weights: agrees with direct evaluation'                 ,covarianceReference                       ,covarianceDirect                   ,absTol=toleranceAbsolute,relTol=toleranceRelative)
  call Assert('well-scaled weights: variances are positive'                        ,all([(covarianceReference(i,i) > 0.0d0,i=1,countBins)]),.true.                                                                 )
  ! Rescale the weights in each bin. The covariance of the mean is invariant under this rescaling. Weights down to ~10⁻⁷⁷
  ! cause U_i² U_j² to underflow, so this checks that this product is never formed.
  scale=[1.0d0,1.0d-50,1.0d-100,1.0d-150]
  call accumulate(scale)
  covarianceMean=Output_Analysis_Mean_Covariance(valueUnweighted,valueWeighted,covarianceUnweighted,covarianceWeighted,covarianceCross)
  call Assert('weights rescaled to 1e-150: covariance unchanged'                  ,covarianceMean                            ,covarianceReference                ,absTol=toleranceAbsolute,relTol=toleranceRelative)
  ! A bin with weight ~10⁻¹⁶⁵ (whose square underflows to zero) - modeled on the weights found in issue #1758.
  scale=[1.0d0,1.0d-28,1.0d-78,1.0d-165]
  call accumulate(scale)
  call Assert('weight 1e-165 test case: squared weight underflows'                ,valueUnweighted(countBins)**2             ,0.0d0                                                                                     )
  covarianceMean=Output_Analysis_Mean_Covariance(valueUnweighted,valueWeighted,covarianceUnweighted,covarianceWeighted,covarianceCross)
  call Assert('weight 1e-165: other bins unchanged'                               ,covarianceMean(1:countBins-1,1:countBins-1),covarianceReference(1:countBins-1,1:countBins-1),absTol=toleranceAbsolute,relTol=toleranceRelative)
  call Assert('weight 1e-165: zero row'                                           ,covarianceMean(countBins,:)               ,[(0.0d0,i=1,countBins)]                                                                   )
  call Assert('weight 1e-165: zero column'                                        ,covarianceMean(:,countBins)               ,[(0.0d0,i=1,countBins)]                                                                   )
  ! An empty bin.
  scale=[1.0d0,1.0d0,1.0d0,0.0d0]
  call accumulate(scale)
  covarianceMean=Output_Analysis_Mean_Covariance(valueUnweighted,valueWeighted,covarianceUnweighted,covarianceWeighted,covarianceCross)
  call Assert('empty bin: other bins unchanged'                                    ,covarianceMean(1:countBins-1,1:countBins-1),covarianceReference(1:countBins-1,1:countBins-1),absTol=toleranceAbsolute,relTol=toleranceRelative)
  call Assert('empty bin: zero row'                                                ,covarianceMean(countBins,:)               ,[(0.0d0,i=1,countBins)]                                                                   )
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

contains

  subroutine accumulate(scale)
    !!{RST
    Accumulate weights and covariances (under the Poisson model) from our galaxies, with the weights in each bin multiplied by
    the given scale factors. Products are formed in the same order as in the :galacticus-class:`outputAnalysisCrossCorrelator1D`
    class, such that they underflow as they would there.
    !!}
    implicit none
    double precision, intent(in   ), dimension(countBins) :: scale
    double precision               , dimension(countBins) :: weightUnweighted, weightWeighted
    integer                                               :: k               , i              , j

    valueUnweighted     =0.0d0
    valueWeighted       =0.0d0
    covarianceUnweighted=0.0d0
    covarianceWeighted  =0.0d0
    covarianceCross     =0.0d0
    do k=1,countGalaxies
       weightUnweighted=scale*distribution(:,k)*weight(k)
       weightWeighted  =scale*distribution(:,k)*weight(k)*property(k)
       valueUnweighted =valueUnweighted+weightUnweighted
       valueWeighted   =valueWeighted  +weightWeighted
       do i=1,countBins
          do j=1,countBins
             covarianceUnweighted(i,j)=covarianceUnweighted(i,j)+weightUnweighted(i)*weightUnweighted(j)
             covarianceWeighted  (i,j)=covarianceWeighted  (i,j)+weightWeighted  (i)*weightWeighted  (j)
             covarianceCross     (i,j)=covarianceCross     (i,j)+weightUnweighted(i)*weightWeighted  (j)
          end do
       end do
    end do
    return
  end subroutine accumulate

end program Test_Output_Analyses_Mean_Covariance
