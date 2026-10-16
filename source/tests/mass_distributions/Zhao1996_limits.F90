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
Contains a program which tests that the :cite:t:`zhao_analytical_1996` mass distribution reduces to the NFW and Hernquist
mass distributions in the appropriate limits.
!!}

program Test_Mass_Distributions_Zhao1996_Limits
  !!{RST
  Tests that the :galacticus-class:`massDistributionZhao1996` mass distribution, with density :math:`\rho(r) = \rho_0
  x^{-\gamma} (1+x^\alpha)^{-(\beta-\gamma)/\alpha}` for :math:`x=r/r_\mathrm{s}`, reduces to the independently-implemented
  :galacticus-class:`massDistributionNFW` for :math:`(\alpha,\beta,\gamma)=(1,3,1)` and to
  :galacticus-class:`massDistributionHernquist` for :math:`(\alpha,\beta,\gamma)=(1,4,1)`, in density, enclosed mass, and
  potential. All three classes share the same definition of :math:`\rho_0`.

  The Zhao profile evaluates :math:`(1,3,1)` using a dedicated NFW solution, but :math:`(1,4,1)` using its general solution
  in terms of hypergeometric functions. The general solution is also tested against NFW by choosing :math:`\gamma=1+\delta`
  with :math:`\delta` just large enough to avoid the NFW special case. The density then differs from NFW by a factor
  :math:`[x/(1+x)]^{-\delta}`, and so by at most :math:`\delta |\ln[x/(1+x)]|` over the radii tested, with comparable
  differences in enclosed mass and potential.

  For :math:`\beta>3` the total mass is finite,

  .. math::

     M_\mathrm{total} = \frac{4\pi}{\alpha} \frac{\Gamma([3-\gamma]/\alpha) \Gamma([\beta-3]/\alpha)}{\Gamma([\beta-\gamma]/\alpha)} \rho_0 r_\mathrm{s}^3,

  which is :math:`2\pi\rho_0 r_\mathrm{s}^3` for :math:`(1,4,1)` (the Hernquist profile), and :math:`(4\pi/3)\rho_0
  r_\mathrm{s}^3` for :math:`(1,4,0)`. The latter case has :math:`\gamma \ne \alpha`, and so tests all of the :math:`\Gamma`
  function arguments.

  The profiles adopt different zero points for the potential, so potentials are compared as differences relative to the
  potential at the scale radius.
  !!}
  use :: Coordinates             , only : coordinateSpherical   , assignment(=)
  use :: Display                 , only : displayVerbositySet   , verbosityLevelStandard
  use :: Events_Hooks            , only : eventsHooksInitialize
  use :: Mass_Distributions      , only : massDistributionNFW   , massDistributionHernquist, massDistributionZhao1996
  use :: Numerical_Constants_Math, only : Pi
  use :: Unit_Tests              , only : Assert                , Unit_Tests_Begin_Group   , Unit_Tests_End_Group    , Unit_Tests_Finish
  implicit none
  type            (massDistributionNFW      )                         :: massDistributionNFW_
  type            (massDistributionHernquist)                         :: massDistributionHernquist_
  type            (massDistributionZhao1996 )                         :: massDistributionZhaoNFW_          , massDistributionZhaoNFWGeneral_, &
       &                                                                 massDistributionZhaoHernquist_    , massDistributionZhaoCored_
  double precision                           , parameter              :: radiusScale             =2.0d-2   , densityNormalization  =1.0d15 ! [Mpc], [M☉/Mpc³]
  ! Offset in γ used to force the general solution. This must exceed the absolute tolerance (10⁻⁶) used by the Zhao profile to
  ! detect the NFW special case.
  double precision                           , parameter              :: gammaOffset             =1.5d-6
  integer                                    , parameter              :: countRadii              =5
  double precision                           , dimension(countRadii)  :: radiusScaleFree         =[1.0d-2,1.0d-1,1.0d0,1.0d1,1.0d2]
  double precision                           , dimension(countRadii)  :: densityZhao                       , densityReference              , &
       &                                                                 massZhao                          , massReference                 , &
       &                                                                 potentialZhao                     , potentialReference
  ! Tolerances. The NFW and Hernquist comparisons are between closed-form expressions, except that the general Zhao solution
  ! involves hypergeometric and Γ functions, evaluated by GSL to near machine precision. The comparison of the general solution
  ! with NFW is limited by the offset in γ, which gives differences no larger than δ ln(101) ≈ 7×10⁻⁶ at the radii tested.
  double precision                           , parameter              :: toleranceClosedForm     =1.0d-10  , toleranceHypergeometric=1.0d-8 , &
       &                                                                 toleranceOffset         =1.0d-4

  call displayVerbositySet  (verbosityLevelStandard)
  call eventsHooksInitialize(                      )
  call Unit_Tests_Begin_Group("Zhao (1996) mass distribution limits")
  massDistributionNFW_           =massDistributionNFW      (                                      scaleLength=radiusScale,densityNormalization=densityNormalization)
  massDistributionHernquist_     =massDistributionHernquist(                                      scaleLength=radiusScale,densityNormalization=densityNormalization)
  massDistributionZhaoNFW_       =massDistributionZhao1996 (alpha=1.0d0,beta=3.0d0,gamma=1.0d0            ,scaleLength=radiusScale,densityNormalization=densityNormalization)
  massDistributionZhaoNFWGeneral_=massDistributionZhao1996 (alpha=1.0d0,beta=3.0d0,gamma=1.0d0+gammaOffset,scaleLength=radiusScale,densityNormalization=densityNormalization)
  massDistributionZhaoHernquist_ =massDistributionZhao1996 (alpha=1.0d0,beta=4.0d0,gamma=1.0d0            ,scaleLength=radiusScale,densityNormalization=densityNormalization)
  massDistributionZhaoCored_     =massDistributionZhao1996 (alpha=1.0d0,beta=4.0d0,gamma=0.0d0            ,scaleLength=radiusScale,densityNormalization=densityNormalization)
  ! NFW special case.
  call Unit_Tests_Begin_Group("(α,β,γ)=(1,3,1) → NFW")
  call evaluate(massDistributionZhaoNFW_       ,massDistributionNFW_      )
  call Assert('density'      ,densityZhao  ,densityReference  ,relTol=toleranceClosedForm    )
  call Assert('enclosed mass',massZhao     ,massReference     ,relTol=toleranceClosedForm    )
  call Assert('potential'    ,potentialZhao,potentialReference,relTol=toleranceClosedForm    )
  call Unit_Tests_End_Group()
  ! General solution close to NFW.
  call Unit_Tests_Begin_Group("(α,β,γ)=(1,3,1+δ) → NFW (general solution)")
  call evaluate(massDistributionZhaoNFWGeneral_,massDistributionNFW_      )
  call Assert('density'      ,densityZhao  ,densityReference  ,relTol=toleranceOffset        )
  call Assert('enclosed mass',massZhao     ,massReference     ,relTol=toleranceOffset        )
  call Assert('potential'    ,potentialZhao,potentialReference,relTol=toleranceOffset        )
  call Unit_Tests_End_Group()
  ! Hernquist, via the general solution.
  call Unit_Tests_Begin_Group("(α,β,γ)=(1,4,1) → Hernquist (general solution)")
  call evaluate(massDistributionZhaoHernquist_ ,massDistributionHernquist_)
  call Assert('density'      ,densityZhao  ,densityReference  ,relTol=toleranceHypergeometric)
  call Assert('enclosed mass',massZhao     ,massReference     ,relTol=toleranceHypergeometric)
  call Assert('potential'    ,potentialZhao,potentialReference,relTol=toleranceHypergeometric)
  call Assert('total mass'   ,massDistributionZhaoHernquist_%massTotal(),massDistributionHernquist_%massTotal(),relTol=toleranceHypergeometric)
  call Unit_Tests_End_Group()
  ! Total mass for a case with γ≠α.
  call Unit_Tests_Begin_Group("(α,β,γ)=(1,4,0) (general solution)")
  call Assert('total mass'   ,massDistributionZhaoCored_%massTotal(),4.0d0*Pi/3.0d0*densityNormalization*radiusScale**3,relTol=toleranceHypergeometric)
  call Unit_Tests_End_Group()
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

contains

  subroutine evaluate(massDistributionZhao,massDistributionReference)
    !!{RST
    Evaluate density, enclosed mass, and potential (relative to that at the scale radius) in a Zhao profile and in the
    reference profile which it should match.
    !!}
    use :: Mass_Distributions, only : massDistributionClass
    implicit none
    class           (massDistributionClass), intent(inout) :: massDistributionZhao, massDistributionReference
    type            (coordinateSpherical  )                :: coordinates         , coordinatesScale
    integer                                                :: i

    coordinatesScale=[radiusScale,0.0d0,0.0d0]
    do i=1,countRadii
       coordinates          =[radiusScaleFree(i)*radiusScale,0.0d0,0.0d0]
       densityZhao       (i)=+massDistributionZhao     %density             (coordinates                  )
       densityReference  (i)=+massDistributionReference%density             (coordinates                  )
       massZhao          (i)=+massDistributionZhao     %massEnclosedBySphere(radiusScaleFree(i)*radiusScale)
       massReference     (i)=+massDistributionReference%massEnclosedBySphere(radiusScaleFree(i)*radiusScale)
       potentialZhao     (i)=+massDistributionZhao     %potential           (coordinates                  ) &
            &                -massDistributionZhao     %potential           (coordinatesScale             )
       potentialReference(i)=+massDistributionReference%potential           (coordinates                  ) &
            &                -massDistributionReference%potential           (coordinatesScale             )
    end do
    ! At the scale radius itself the potential difference is identically zero, which can not be compared with a relative
    ! tolerance, so replace it by the potential at the scale radius relative to that at ten scale radii.
    coordinates          =[1.0d1*radiusScale,0.0d0,0.0d0]
    potentialZhao     (3)=+massDistributionZhao     %potential(coordinatesScale)-massDistributionZhao     %potential(coordinates)
    potentialReference(3)=+massDistributionReference%potential(coordinatesScale)-massDistributionReference%potential(coordinates)
    return
  end subroutine evaluate

end program Test_Mass_Distributions_Zhao1996_Limits
