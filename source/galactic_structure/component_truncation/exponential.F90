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
  Implementation of a galactic component truncation class which applies an exponential truncation.
  !!}

  use :: Numerical_Interpolation, only : interpolator

  !![
  <galacticComponentTruncation name="galacticComponentTruncationExponential" docformat="rst">
    <description>
    A galactic component truncation class which truncates the mass distributions of galactic components exponentially (see
    :galacticus-class:`massDistributionSphericalTruncatedExponential`). The truncation begins at radius
    :math:`r_\mathrm{t}`, beyond which the density decays exponentially with scale length
    :math:`r_\mathrm{decay}=`\ ``[radiusDecayOverRadiusTruncate]``\ :math:`\times r_\mathrm{t}`. The truncation radius is chosen
    such that the truncated distribution retains the required fraction, :math:`f`, of the mass of the untruncated distribution.

    The relation between :math:`f` and :math:`r_\mathrm{t}` depends only on the shape of the dimensionless mass distribution, so
    it is tabulated once (as a function of :math:`x=-\log f`) and interpolated thereafter. Truncations with
    :math:`x &lt; 10^{-4}` are negligible and are ignored, while those beyond the tabulated range (in which the component
    retains a fraction of order :math:`10^{-11}` of its original mass) are limited to the smallest tabulated truncation radius.
    </description>
  </galacticComponentTruncation>
  !!]
  type, extends(galacticComponentTruncationClass) :: galacticComponentTruncationExponential
     !!{RST
     A galactic component truncation class which applies an exponential truncation.
     !!}
     private
     double precision               :: radiusDecayOverRadiusTruncate
     logical                        :: tabulated                    =.false.
     double precision               :: truncationLogarithmicMaximum
     type            (interpolator) :: interpolatorRadius
   contains
     !![
     <methods>
       <method method="tabulate" description="Tabulate the truncation radius as a function of the degree of truncation."/>
     </methods>
     !!]
     procedure :: isActive         => exponentialIsActive
     procedure :: massDistribution => exponentialMassDistribution
     procedure :: tabulate         => exponentialTabulate
  end type galacticComponentTruncationExponential

  interface galacticComponentTruncationExponential
     !!{RST
     Constructors for the :galacticus-class:`galacticComponentTruncationExponential` galactic component truncation class.
     !!}
     module procedure exponentialConstructorParameters
     module procedure exponentialConstructorInternal
  end interface galacticComponentTruncationExponential

  ! Limits on the tabulated degree of truncation, x=-log(f). Truncations smaller than the lower limit are negligible. The upper
  ! limit corresponds to a component retaining a fraction of order 10⁻¹¹ of its original mass.
  double precision, parameter :: exponentialTruncationNegligible=1.0d-4, exponentialTruncationLarge=25.0d0

contains

  function exponentialConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`galacticComponentTruncationExponential` galactic component truncation class which
    takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticComponentTruncationExponential)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    double precision                                                        :: radiusDecayOverRadiusTruncate

    !![
    <inputParameter docformat="rst">
      <name>radiusDecayOverRadiusTruncate</name>
      <defaultValue>0.1d0</defaultValue>
      <description>
      The scale length of the exponential decay of the density beyond the truncation radius, in units of the truncation radius.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=galacticComponentTruncationExponential(radiusDecayOverRadiusTruncate)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function exponentialConstructorParameters

  function exponentialConstructorInternal(radiusDecayOverRadiusTruncate) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`galacticComponentTruncationExponential` galactic component truncation class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (galacticComponentTruncationExponential)                :: self
    double precision                                        , intent(in   ) :: radiusDecayOverRadiusTruncate
    !![
    <constructorAssign variables="radiusDecayOverRadiusTruncate"/>
    !!]

    if (radiusDecayOverRadiusTruncate <= 0.0d0) call Error_Report('`radiusDecayOverRadiusTruncate` must be positive'//{introspection:location})
    return
  end function exponentialConstructorInternal

  logical function exponentialIsActive(self) result(isActive)
    !!{RST
    Return true, as this class truncates mass distributions.
    !!}
    implicit none
    class(galacticComponentTruncationExponential), intent(inout) :: self
    !$GLC attributes unused :: self

    isActive=.true.
    return
  end function exponentialIsActive

  function exponentialMassDistribution(self,massDistribution_,fractionMassRetained) result(massDistributionTruncated)
    !!{RST
    Return the exponentially-truncated mass distribution.
    !!}
    use :: Mass_Distributions, only : massDistributionSphericalTruncatedExponential, nonAnalyticSolversNumerical
    implicit none
    class           (massDistributionSpherical                    ), pointer               :: massDistributionTruncated
    class           (galacticComponentTruncationExponential       ), intent(inout)         :: self
    class           (massDistributionSpherical                    ), intent(inout), target :: massDistribution_
    double precision                                               , intent(in   )         :: fractionMassRetained
    type            (massDistributionSphericalTruncatedExponential), pointer               :: massDistributionTruncatedExponential
    double precision                                                                       :: truncation                          , radiusTruncate

    truncation=-log(max(fractionMassRetained,tiny(0.0d0)))
    if (truncation < exponentialTruncationNegligible) then
       ! Truncation is negligible - return the untruncated distribution.
       massDistributionTruncated => massDistribution_
       !![
       <referenceCountIncrement object="massDistributionTruncated"/>
       !!]
    else
       if (.not.self%tabulated) call self%tabulate(massDistribution_)
       radiusTruncate=exp(self%interpolatorRadius%interpolate(min(log(truncation),self%truncationLogarithmicMaximum)))
       allocate(massDistributionTruncatedExponential)
       !![
       <referenceConstruct object="massDistributionTruncatedExponential" constructor="massDistributionSphericalTruncatedExponential(radiusTruncateMinimum=radiusTruncate,radiusTruncateDecay=radiusTruncate*self%radiusDecayOverRadiusTruncate,nonAnalyticSolver=nonAnalyticSolversNumerical,massDistribution_=massDistribution_)"/>
       !!]
       massDistributionTruncated => massDistributionTruncatedExponential
    end if
    return
  end function exponentialMassDistribution

  subroutine exponentialTabulate(self,massDistribution_)
    !!{RST
    Tabulate the truncation radius as a function of the degree of truncation, :math:`x=-\log f`. Truncation radii are stepped
    logarithmically outward and inward from the length scale of the dimensionless distribution until the degree of truncation
    spans the range of interest.
    !!}
    use :: Error             , only : Error_Report
    use :: Mass_Distributions, only : massDistributionSphericalTruncatedExponential, nonAnalyticSolversNumerical
    use :: Table_Labels      , only : extrapolationTypeFix
    implicit none
    class           (galacticComponentTruncationExponential), intent(inout)                        :: self
    class           (massDistributionSpherical             ), intent(inout), target                :: massDistribution_
    integer                                                 , parameter                            :: countPerDecade =50, countMaximum   =1000
    double precision                                        , dimension(countMaximum)              :: truncationOuter   , truncationInner     , &
         &                                                                                            radiusOuter       , radiusInner
    double precision                                        , dimension(:           ), allocatable :: truncations       , radii
    double precision                                                                               :: radiusStep        , radius              , &
         &                                                                                            massUntruncated
    integer                                                                                        :: countOuter        , countInner          , &
         &                                                                                            i

    massUntruncated=massDistribution_%massTotal()
    ! Step outward from the length scale of the distribution until truncation becomes negligible.
    radiusStep=10.0d0**(1.0d0/dble(countPerDecade))
    radius    = 1.0d0
    countOuter= 0
    do while (countOuter < countMaximum)
       countOuter=countOuter+1
       radiusOuter    (countOuter)=radius
       truncationOuter(countOuter)=truncationEvaluate(radius)
       if (truncationOuter(countOuter) < exponentialTruncationNegligible) exit
       radius=radius*radiusStep
    end do
    ! Step inward from the length scale of the distribution until truncation becomes large.
    radius    =1.0d0/radiusStep
    countInner=0
    do while (countInner < countMaximum)
       countInner=countInner+1
       radiusInner    (countInner)=radius
       truncationInner(countInner)=truncationEvaluate(radius)
       if (truncationInner(countInner) > exponentialTruncationLarge     ) exit
       radius=radius/radiusStep
    end do
    if     (                                                                       &
         &   truncationOuter(countOuter) >= exponentialTruncationNegligible        &
         &  .or.                                                                   &
         &   truncationInner(countInner) <= exponentialTruncationLarge             &
         & ) call Error_Report('failed to span the required range of truncation'//{introspection:location})
    ! Assemble the table in order of increasing truncation (i.e. decreasing truncation radius).
    truncations=[truncationOuter(countOuter:1:-1),truncationInner(1:countInner)]
    radii      =[radiusOuter    (countOuter:1:-1),radiusInner    (1:countInner)]
    do i=2,size(truncations)
       if (truncations(i) <= truncations(i-1)) call Error_Report('truncation is not a monotonic function of truncation radius'//{introspection:location})
    end do
    self%truncationLogarithmicMaximum=log(truncations(size(truncations)))
    self%interpolatorRadius          =interpolator(log(truncations),log(radii),extrapolationType=extrapolationTypeFix)
    self%tabulated                   =.true.
    return

  contains

    double precision function truncationEvaluate(radiusTruncate) result(truncation)
      !!{RST
      Evaluate the degree of truncation for the given truncation radius.
      !!}
      implicit none
      double precision                                               , intent(in   ) :: radiusTruncate
      type            (massDistributionSphericalTruncatedExponential), pointer       :: massDistributionTruncatedExponential

      allocate(massDistributionTruncatedExponential)
      !![
      <referenceConstruct object="massDistributionTruncatedExponential" constructor="massDistributionSphericalTruncatedExponential(radiusTruncateMinimum=radiusTruncate,radiusTruncateDecay=radiusTruncate*self%radiusDecayOverRadiusTruncate,nonAnalyticSolver=nonAnalyticSolversNumerical,massDistribution_=massDistribution_)"/>
      !!]
      truncation=-log(massDistributionTruncatedExponential%massTotal()/massUntruncated)
      !![
      <objectDestructor name="massDistributionTruncatedExponential"/>
      !!]
      return
    end function truncationEvaluate

  end subroutine exponentialTabulate
