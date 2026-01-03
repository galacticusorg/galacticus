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
  An implementation of cosmological density field mass variance which scales that of another class.
  !!}

  !![
  <cosmologicalMassVariance name="cosmologicalMassVarianceScaled">
   <description>
    The mass variance of cosmological density fields is computed by scaling that from another class by a factor of {\normalfont \ttfamily [scale]}.
   </description>
  </cosmologicalMassVariance>
  !!]
  type, extends(cosmologicalMassVarianceClass) :: cosmologicalMassVarianceScaled
     !!{
     A cosmological mass variance class that scales the variance from another class.
     !!}
     private
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     double precision                                         :: scale                              , rootScale
   contains
     final     ::                                        scaledDestructor
     procedure :: sigma8                              => scaledSigma8
     procedure :: powerNormalization                  => scaledPowerNormalization
     procedure :: rootVariance                        => scaledRootVariance
     procedure :: rootVarianceLogarithmicGradient     => scaledRootVarianceLogarithmicGradient
     procedure :: rootVarianceLogarithmicGradientTime => scaledRootVarianceLogarithmicGradientTime
     procedure :: rootVarianceAndLogarithmicGradient  => scaledRootVarianceAndLogarithmicGradient
     procedure :: mass                                => scaledMass   
     procedure :: growthIsMassDependent               => scaledGrowthIsMassDependent
  end type cosmologicalMassVarianceScaled

  interface cosmologicalMassVarianceScaled
     !!{
     Constructors for the \refClass{cosmologicalMassVarianceScaled} cosmological mass variance class.
     !!}
     module procedure scaledConstructorParameters
     module procedure scaledConstructorInternal
  end interface cosmologicalMassVarianceScaled

contains

  function scaledConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{cosmologicalMassVarianceScaled} cosmological mass variance class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    use :: Error           , only : Error_Report
    implicit none
    type            (cosmologicalMassVarianceScaled)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    class           (cosmologicalMassVarianceClass ), pointer       :: cosmologicalMassVariance_
    double precision                                                :: scale


    !![
    <inputParameter>
      <name>scale</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>The factor by which to scale the variance.</description>
    </inputParameter>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=cosmologicalMassVarianceScaled(scale,cosmologicalMassVariance_)
    !![
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]    
    return
  end function scaledConstructorParameters

  function scaledConstructorInternal(scale,cosmologicalMassVariance_) result(self)
    !!{
    Internal constructor for the \refClass{cosmologicalMassVarianceScaled} linear growth class.
    !!}
    implicit none
    type            (cosmologicalMassVarianceScaled)                        :: self
    class           (cosmologicalMassVarianceClass ), intent(in   ), target :: cosmologicalMassVariance_
    double precision                                , intent(in   )         :: scale
    !![
    <constructorAssign variables="scale, *cosmologicalMassVariance_"/>
    !!]

    self%rootScale=sqrt(scale)
    return
  end function scaledConstructorInternal

  subroutine scaledDestructor(self)
    !!{
    Destructor for the \refClass{cosmologicalMassVarianceScaled} linear growth class.
    !!}
    implicit none
    type   (cosmologicalMassVarianceScaled), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine scaledDestructor

  double precision function scaledPowerNormalization(self)
    !!{
    Return the normalization of the power spectrum.
    !!}
    implicit none
    class(cosmologicalMassVarianceScaled), intent(inout) :: self

    scaledPowerNormalization=self%cosmologicalMassVariance_%powerNormalization()
    return
  end function scaledPowerNormalization

  double precision function scaledSigma8(self)
    !!{
    Return the value of $\sigma_8$.
    !!}
    implicit none
    class(cosmologicalMassVarianceScaled), intent(inout) :: self

    scaledSigma8=self%cosmologicalMassVariance_%sigma8()
    return
  end function scaledSigma8

  double precision function scaledRootVariance(self,mass,time)
    !!{
    Return the root-variance of the cosmological density field in a spherical region containing the given {\normalfont
    \ttfamily mass} on average.
    !!}
    implicit none
    class           (cosmologicalMassVarianceScaled), intent(inout) :: self
    double precision                                , intent(in   ) :: mass, time

    scaledRootVariance=+self%cosmologicalMassVariance_%rootVariance(mass,time) &
         &             *self                          %rootScale
    return
  end function scaledRootVariance

  double precision function scaledRootVarianceLogarithmicGradient(self,mass,time)
    !!{
    Return the logarithmic gradient with respect to mass of the root-variance of the cosmological density field in a spherical
    region containing the given {\normalfont \ttfamily mass} on average.
    !!}
    implicit none
    class           (cosmologicalMassVarianceScaled), intent(inout) :: self
    double precision                                , intent(in   ) :: mass, time

    scaledRootVarianceLogarithmicGradient=self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient(mass,time)
    return
  end function scaledRootVarianceLogarithmicGradient

  double precision function scaledRootVarianceLogarithmicGradientTime(self,mass,time)
    !!{
    Return the logarithmic gradient with respect to time of the root-variance of the cosmological density field in a spherical
    region containing the given {\normalfont \ttfamily mass} on average.
    !!}
    implicit none
    class           (cosmologicalMassVarianceScaled), intent(inout) :: self
    double precision                                , intent(in   ) :: mass, time

    scaledRootVarianceLogarithmicGradientTime=self%cosmologicalMassVariance_%rootVarianceLogarithmicGradientTime(mass,time)
   return
  end function scaledRootVarianceLogarithmicGradientTime

  subroutine scaledRootVarianceAndLogarithmicGradient(self,mass,time,rootVariance,rootVarianceLogarithmicGradient)
    !!{
    Return the value and logarithmic gradient with respect to mass of the root-variance of the cosmological density field in a
    spherical region containing the given {\normalfont \ttfamily mass} on average.
    !!}
    implicit none
    class           (cosmologicalMassVarianceScaled), intent(inout) :: self
    double precision                                , intent(in   ) :: mass        , time
    double precision                                , intent(  out) :: rootVariance, rootVarianceLogarithmicGradient

    call self%cosmologicalMassVariance_%rootVarianceAndLogarithmicGradient(mass,time,rootVariance,rootVarianceLogarithmicGradient)
    rootVariance=+     rootVariance &
         &       *self%rootScale
    return
  end subroutine scaledRootVarianceAndLogarithmicGradient

  double precision function scaledMass(self,rootVariance,time)
    !!{
    Return the mass corresponding to the given {\normalfont \ttfamily } root-variance of the cosmological density field.
    !!}
    implicit none
    class           (cosmologicalMassVarianceScaled), intent(inout) :: self
    double precision                                , intent(in   ) :: rootVariance, time

    scaledMass=self%cosmologicalMassVariance_%mass(rootVariance*self%rootScale,time)
    return
  end function scaledMass

  logical function scaledGrowthIsMassDependent(self)
    !!{
    Return true if the growth rate of the variance is mass-dependent.
    !!}
    implicit none
    class(cosmologicalMassVarianceScaled), intent(inout) :: self

    scaledGrowthIsMassDependent=self%cosmologicalMassVariance_%growthIsMassDependent()
    return
  end function scaledGrowthIsMassDependent
