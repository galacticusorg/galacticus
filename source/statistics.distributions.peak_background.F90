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
  Implementation of a peak-background split density 1D distribution function.
  !!}

  use :: Tables, only : table1D, table1DLinearLinear

  !![
  <distributionFunction1D name="distributionFunction1DPeakBackground">
   <description>A peakBackground 1D distribution function class.</description>
  </distributionFunction1D>
  !!]
  type, extends(distributionFunction1DClass) :: distributionFunction1DPeakBackground
     !!{
     Implementation of a peakBackground 1D distribution function.
     !!}
     private
     double precision                                   :: varianceBackground, thresholdCollapse, &
          &                                                normalization
     type            (table1DLinearLinear)              :: cdf
     class           (table1D            ), allocatable :: cdfInverse
   contains
     final     ::               peakBackgroundDestructor
     procedure :: density    => peakBackgroundDensity
     procedure :: cumulative => peakBackgroundCumulative
     procedure :: inverse    => peakBackgroundInverse
     procedure :: minimum    => peakBackgroundMinimum
     procedure :: maximum    => peakBackgroundMaximum
  end type distributionFunction1DPeakBackground

  interface distributionFunction1DPeakBackground
     !!{
     Constructors for the \refClass{distributionFunction1DPeakBackground} 1D distribution function class.
     !!}
     module procedure peakBackgroundConstructorParameters
     module procedure peakBackgroundConstructorInternal
  end interface distributionFunction1DPeakBackground

contains

  function peakBackgroundConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{distributionFunction1DPeakBackground} 1D distribution function class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (distributionFunction1DPeakBackground)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass          ), pointer       :: randomNumberGenerator_
    double precision                                                      :: varianceBackground    , thresholdCollapse

    !![
    <inputParameter>
      <name>varianceBackground</name>
      <description>The variance in the background density field.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>thresholdCollapse</name>
      <description>The threshold for collapse of density perturbations.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=distributionFunction1DPeakBackground(varianceBackground,thresholdCollapse,randomNumberGenerator_)
    !![
    <objectDestructor name="randomNumberGenerator_"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function peakBackgroundConstructorParameters

  function peakBackgroundConstructorInternal(varianceBackground,thresholdCollapse,randomNumberGenerator_) result(self)
    !!{
    Constructor for the \refClass{distributionFunction1DPeakBackground} 1D distribution function class.
    !!}
    use :: Error_Functions, only : Error_Function
    implicit none
    type            (distributionFunction1DPeakBackground)                                  :: self
    double precision                                      , intent(in   )                   :: varianceBackground          , thresholdCollapse
    class           (randomNumberGeneratorClass          ), intent(in   ), target, optional :: randomNumberGenerator_
    double precision                                      , parameter                       :: rangeTable            =7.0d0
    integer                                               , parameter                       :: cdfCount              =1000
    integer                                                                                 :: i
    !![
    <constructorAssign variables="varianceBackground, thresholdCollapse, *randomNumberGenerator_"/>
    !!]

    ! Compute normalization of the distribution.
    self%normalization     =+1.0d0                                                                  &
         &                  /Error_Function(thresholdCollapse/sqrt(2.0d0)/sqrt(varianceBackground))
    ! Tabulate the cumulative distribution function.
    call self%cdf%create(-rangeTable*sqrt(varianceBackground),thresholdCollapse,cdfCount)
    do i=1,cdfCount
       call self%cdf%populate(                                                                                            &
            &                 +self%normalization                                                                         &
            &                 *0.50d0                                                                                     &
            &                 *(                                                                                          &
            &                   +Error_Function((+self%cdf%x(i)                        )/sqrt(2.0d0*varianceBackground))  &
            &                   -Error_Function((+self%cdf%x(i)-2.0d0*thresholdCollapse)/sqrt(2.0d0*varianceBackground))  &
            &                  )                                                                                        , &
            &                 i                                                                                           &
            &                )
       if (i > 1) then
          ! Test for monotonicity. If monotonicity fails enforce it by adding a small increase in the CDF. This will make
          ! negligible difference to our results.
          if (self%cdf%y(i) <= self%cdf%y(i-1)) call self%cdf%populate(self%cdf%y(i-1)*(1.0d0+epsilon(0.0d0)),i)
       end if
    end do
    call self%cdf%reverse(self%cdfInverse)
    return

  contains

    double precision function cdfIntegrand(x)
      !!{
      The integrand for the cumulative distribution function.
      !!}
      implicit none
      double precision, intent(in   ) :: x

      cdfIntegrand=self%density(x)
      return
    end function cdfIntegrand

  end function peakBackgroundConstructorInternal

  subroutine peakBackgroundDestructor(self)
    !!{
    Destructor for the \refClass{distributionFunction1DPeakBackground} 1D distribution function class.
    !!}
    implicit none
    type(distributionFunction1DPeakBackground), intent(inout) :: self

    call    self%cdf       %destroy()
    if (allocated(self%cdfInverse)) then
       call self%cdfInverse%destroy()
       deallocate(self%cdfInverse)
    end if
    return
  end subroutine peakBackgroundDestructor

  double precision function peakBackgroundMinimum(self)
    !!{
    Return the minimum possible value of a peak-background split distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(distributionFunction1DPeakBackground), intent(inout) :: self
    !$GLC attributes unused :: self

    peakBackgroundMinimum=0.0d0
    call Error_Report('no minimum exists'//{introspection:location})
    return
  end function peakBackgroundMinimum

  double precision function peakBackgroundMaximum(self)
    !!{
    Return the maximum possible value of a peak-background split distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(distributionFunction1DPeakBackground), intent(inout) :: self
    !$GLC attributes unused :: self

    peakBackgroundMaximum=0.0d0
    call Error_Report('no maximum exists'//{introspection:location})
    return
  end function peakBackgroundMaximum

  double precision function peakBackgroundDensity(self,x)
    !!{
    Return the density of a normal distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (distributionFunction1DPeakBackground), intent(inout) :: self
    double precision                                      , intent(in   ) :: x

    if (x < self%thresholdCollapse) then
       peakBackgroundDensity=+self%normalization                  &
            &                /sqrt(                               &
            &                      +2.0d0                         &
            &                      *Pi                            &
            &                     )                               &
            &                *(                                   &
            &                  +exp(                              &
            &                       -0.5d0                        &
            &                       *x                        **2 &
            &                       /  self%varianceBackground    &
            &                      )                              &
            &                  -exp(                              &
            &                       -0.5d0                        &
            &                       *(                            &
            &                         +x                          &
            &                         -2.0d0                      &
            &                         *self%thresholdCollapse     &
            &                        )                        **2 &
            &                       /  self%varianceBackground    &
            &                      )                              &
            &                 )                                   &
            &                /sqrt(                               &
            &                      +   self%varianceBackground    &
            &                     )
    else
       peakBackgroundDensity=+0.0d0
    end if
    return
  end function peakBackgroundDensity

  double precision function peakBackgroundCumulative(self,x)
    !!{
    Return the cumulative probability of a normal distribution.
    !!}
    implicit none
    class           (distributionFunction1DPeakBackground), intent(inout) :: self
    double precision                                      , intent(in   ) :: x

    if      (x < self%cdf%x(+1)) then
       peakBackgroundCumulative=0.0d0
    else if (x > self%cdf%x(-1)) then
       peakBackgroundCumulative=1.0d0
    else
       peakBackgroundCumulative=self%cdf%interpolate(x)
    end if
    return
  end function peakBackgroundCumulative

  double precision function peakBackgroundInverse(self,p)
    !!{
    Return the inverse of a normal distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (distributionFunction1DPeakBackground), intent(inout), target :: self
    double precision                                      , intent(in   )         :: p

    if (p < 0.0d0 .or. p > 1.0d0)                                    &
         & call Error_Report(                             &
         &                              'probability out of range'// &
         &                              {introspection:location}     &
         &                             )
    peakBackgroundInverse=self%cdfInverse%interpolate(p)
    return
  end function peakBackgroundInverse
