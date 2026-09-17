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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a program which tests the velocity dispersion of the :cite:t:`zhao_analytical_1996` mass distribution against
numerical integration of the Jeans equation.
!!}

program Test_Mass_Distributions_Zhao1996_Dispersion
  !!{RST
  Tests the 1D velocity dispersion of :galacticus-class:`massDistributionZhao1996` against an independent reference
  (``zhao1996DispersionCheck.py`` in `galacticusDevTools <https://github.com/galacticusorg/galacticusDevTools>`_), which
  integrates the isotropic Jeans equation numerically.

  :galacticus-class:`kinematicsDistributionZhao1996` evaluates the dispersion of a self-gravitating profile from closed-form
  solutions for four special cases of :math:`(\alpha,\beta,\gamma)=(1,3,\gamma)`, with :math:`\gamma \in \{0,1/2,1,3/2\}`. Each
  is implemented as three branches - a series for small radii, a full solution, and a series for large radii - switching at
  :math:`r/r_\mathrm{s} = 10^{-3}` and :math:`10^2`. The radii below span all three, so that all twelve branches are exercised.
  ``tests.mass_distributions.Zhao1996_limits`` covers density, enclosed mass and potential, but not the dispersion.

  The reference must be built from the Jeans integral rather than by comparing the series against the full solutions: those lose
  accuracy catastrophically to cancellation at large radius - the :math:`\gamma=1` solution is wrong by a factor of 240 at
  :math:`r/r_\mathrm{s}=10^4` - which is precisely why the large-radius series exist. Comparing against them would condemn the
  series wrongly.

  The dispersion is also checked against the independently-implemented :galacticus-class:`massDistributionNFW` in the
  :math:`\gamma=1` limit, where the two profiles are identical.
  !!}
  use :: Coordinates                     , only : coordinateSpherical           , assignment(=)
  use :: Display                         , only : displayVerbositySet           , verbosityLevelStandard
  use :: Events_Hooks                    , only : eventsHooksInitialize
  use :: Mass_Distributions              , only : massDistributionZhao1996      , massDistributionNFW      , kinematicsDistributionClass, &
       &                                          kinematicsDistributionZhao1996, kinematicsDistributionNFW
  use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
  use :: Unit_Tests                      , only : Assert                        , Unit_Tests_Begin_Group   , Unit_Tests_End_Group       , &
       &                                          Unit_Tests_Finish
  implicit none
  ! These are allocatable, and reallocated for each shape parameter below: assigning a constructor result onto an object which
  ! has already been used double-finalizes its reference-counted components.
  type            (massDistributionZhao1996  ), target , allocatable    :: massDistributionZhao_
  type            (massDistributionNFW       ), target , allocatable    :: massDistributionNFW_
  class           (kinematicsDistributionClass), pointer                :: kinematicsDistribution_ => null(), kinematicsDistributionNFW_ => null()
  type            (coordinateSpherical        )                         :: coordinates
  double precision                             , parameter              :: radiusScale            =2.0d-2, densityNormalization=1.0d15 ! [Mpc], [M☉/Mpc³]
  integer                                      , parameter              :: countShapes            =4     , countRadii          =10
  double precision                             , dimension(countShapes) :: gammas                 =[0.0d0,0.5d0,1.0d0,1.5d0]
  double precision                             , dimension(countRadii ) :: radii                  =[1.00000000000000d-05,1.00000000000000d-04,1.00000000000000d-03,1.00000000000000d-02,1.00000000000000d-01,1.00000000000000d+00,1.00000000000000d+01,1.00000000000000d+02,1.00000000000000d+03,1.00000000000000d+04]
  ! σ²/(Gρ₀r_s²) from the Jeans equation, one row per γ.
  double precision                             , dimension(countRadii,countShapes) :: dispersionSquaredReference=reshape([                                                                                                                                                                                                                  &
       &                                                                                                                 2.95709806505266d-01,2.95789633807025d-01,2.96586832453910d-01,3.04452563341716d-01,3.73553029987541d-01,6.24428405300340d-01,4.26770200751873d-01,1.06825383224977d-01,1.77929378699372d-02,2.50104389285257d-03, &
       &                                                                                                                 4.78736670332368d-03,1.51378652858839d-02,4.78348244157399d-02,1.50159642036469d-01,4.42508551921485d-01,8.33865255833240d-01,4.82986604119177d-01,1.13522589422858d-01,1.84801669022655d-02,2.57000721519792d-03, &
       &                                                                                                                 6.20952177751750d-04,4.76378179319818d-03,3.32459425452428d-02,1.92781294060178d-01,7.40078984189940d-01,1.17419414595177d+00,5.57687568152596d-01,1.22122095449207d-01,1.93581693328800d-02,2.65805340367372d-03, &
       &                                                                                                                 2.64857183186662d-02,8.36159875560162d-02,2.61326490718637d-01,7.69054687449379d-01,1.72050279960631d+00,1.77299687367832d+00,6.64935353629876d-01,1.34041748233601d-01,2.05688596415887d-02,2.77937522435531d-03  &
       &                                                                                                                ],[countRadii,countShapes])
  double precision                             , dimension(countRadii ) :: dispersionSquared             , dispersionSquaredNFW
  ! The reference is converged to better than 10⁻¹¹; the closed forms are exact, so the comparison is limited by the quadrature.
  double precision                             , parameter              :: tolerance              =1.0d-7
  integer                                                               :: i                             , j

  call displayVerbositySet  (verbosityLevelStandard)
  call eventsHooksInitialize(                      )
  call Unit_Tests_Begin_Group("Zhao (1996) velocity dispersion")
  do j=1,countShapes
     if (allocated(massDistributionZhao_)) deallocate(massDistributionZhao_)
     allocate(massDistributionZhao_)
     massDistributionZhao_  =massDistributionZhao1996(alpha=1.0d0,beta=3.0d0,gamma=gammas(j),scaleLength=radiusScale,densityNormalization=densityNormalization)
     if (associated(kinematicsDistribution_)) deallocate(kinematicsDistribution_)
     allocate(kinematicsDistributionZhao1996 :: kinematicsDistribution_)
     select type (kinematicsDistribution_)
     type is (kinematicsDistributionZhao1996)
        kinematicsDistribution_=kinematicsDistributionZhao1996()
     end select
     do i=1,countRadii
        coordinates          =[radii(i)*radiusScale,0.0d0,0.0d0]
        dispersionSquared(i)=kinematicsDistribution_%velocityDispersion1D(coordinates,massDistributionZhao_,massDistributionZhao_)**2
     end do
     call Assert('σ²(r) against the Jeans equation',dispersionSquared,dispersionSquaredReference(:,j)*gravitationalConstant_internal*densityNormalization*radiusScale**2,relTol=tolerance)
  end do
  ! The γ=1 case is NFW, which is implemented independently.
  if (allocated(massDistributionZhao_)) deallocate(massDistributionZhao_)
  allocate(massDistributionZhao_)
  allocate(massDistributionNFW_ )
  massDistributionZhao_     =massDistributionZhao1996(alpha=1.0d0,beta=3.0d0,gamma=1.0d0,scaleLength=radiusScale,densityNormalization=densityNormalization)
  massDistributionNFW_      =massDistributionNFW     (                                   scaleLength=radiusScale,densityNormalization=densityNormalization)
  if (associated(kinematicsDistribution_)) deallocate(kinematicsDistribution_)
  allocate(kinematicsDistributionZhao1996 :: kinematicsDistribution_   )
  allocate(kinematicsDistributionNFW      :: kinematicsDistributionNFW_)
  select type (kinematicsDistribution_)
  type is (kinematicsDistributionZhao1996)
     kinematicsDistribution_   =kinematicsDistributionZhao1996()
  end select
  select type (kinematicsDistributionNFW_)
  type is (kinematicsDistributionNFW     )
     kinematicsDistributionNFW_=kinematicsDistributionNFW     (useSeriesApproximation=.false.)
  end select
  do i=1,countRadii
     coordinates             =[radii(i)*radiusScale,0.0d0,0.0d0]
     dispersionSquared   (i)=kinematicsDistribution_   %velocityDispersion1D(coordinates,massDistributionZhao_,massDistributionZhao_)**2
     dispersionSquaredNFW(i)=kinematicsDistributionNFW_%velocityDispersion1D(coordinates,massDistributionNFW_ ,massDistributionNFW_ )**2
  end do
  call Assert('(α,β,γ)=(1,3,1) → the NFW dispersion',dispersionSquared,dispersionSquaredNFW,relTol=tolerance)
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Mass_Distributions_Zhao1996_Dispersion
