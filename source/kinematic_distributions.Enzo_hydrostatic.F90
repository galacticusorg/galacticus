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
  Implementation of a kinematic distribution class that mimics the ``hydrostatic'' solution from the Enzo code.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionEnzoHydrostatic">
    <description>
      A kinematic class that implements the ``hydrostatic'' temperature profile available in the \gls{enzo}
      code. Specifically,
      \begin{equation}
        T(r) = \hbox{max}\left( {\mathrm{G} M(&lt;r) \mu m_\mathrm{H} \over 3 \mathrm{k_B} r} , T_\mathrm{min} \right),
      \end{equation}
      where $M(&lt;r)$ is the total mass enclosed within radius $r$, $\mu$ is the primordial mean atomic mass, and
      $T_\mathrm{min}=100$~K is a temperature floor introduced so as to avoid the temperature reaching arbitrarily low values.
    </description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionClass) :: kinematicsDistributionEnzoHydrostatic
     !!{
     A enzoHydrostatic kinematic distribution.
     !!}
     class(massDistributionClass), pointer :: massDistribution_ => null()
   contains
     final     ::                                   enzoHydrostaticDestructor
     procedure :: isCollisional                  => enzoHydrostaticIsCollisional
     procedure :: temperature                    => enzoHydrostaticTemperature
     procedure :: temperatureGradientLogarithmic => enzoHydrostaticTemperature
  end type kinematicsDistributionEnzoHydrostatic

  interface kinematicsDistributionEnzoHydrostatic
     !!{
     Constructors for the \refClass{kinematicsDistributionEnzoHydrostatic} kinematic distribution class.
     !!}
     module procedure enzoHydrostaticKinematicsConstructorParameters
     module procedure enzoHydrostaticKinematicsConstructorInternal
  end interface kinematicsDistributionEnzoHydrostatic

  ! Minimum temperature allowed in this distribution.
  double precision, parameter :: temperatureMinimum=1.0d2

contains

  function enzoHydrostaticKinematicsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{kinematicsDistributionEnzoHydrostatic} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (kinematicsDistributionEnzoHydrostatic)                :: self
    type (inputParameters                      ), intent(inout) :: parameters
    class(massDistributionClass                ), pointer       :: massDistribution_

    !![
    <objectBuilder class="massDistribution" name="massDistribution_" source="parameters"/>
    !!]
    self=kinematicsDistributionEnzoHydrostatic(massDistribution_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function enzoHydrostaticKinematicsConstructorParameters
  
  function enzoHydrostaticKinematicsConstructorInternal(massDistribution_) result(self)
    !!{
    Constructor for the \refClass{kinematicsDistributionEnzoHydrostatic} kinematics distribution class.
    !!}
    implicit none
    type (kinematicsDistributionEnzoHydrostatic)                        :: self
    class(massDistributionClass                ), intent(in   ), target :: massDistribution_
    !![
    <constructorAssign variables="*massDistribution_"/>
    !!]

    return
  end function enzoHydrostaticKinematicsConstructorInternal

  subroutine enzoHydrostaticDestructor(self)
    !!{
    Destructor for the \refClass{kinematicsDistributionEnzoHydrostatic} kinematic distribution class.
    !!}
    type(kinematicsDistributionEnzoHydrostatic), intent(inout) :: self
    implicit none

    !![
    <objectDestructor name="self%massDistribution_"/>
    !!]
    return
  end subroutine enzoHydrostaticDestructor
  
  logical function enzoHydrostaticIsCollisional(self)
    !!{
    Return false indicating that the enzoHydrostatic kinematic distribution represents collisionless particles.
    !!}
    implicit none
    class(kinematicsDistributionEnzoHydrostatic), intent(inout) :: self
    
    enzoHydrostaticIsCollisional=.true.
    return
  end function enzoHydrostaticIsCollisional

  double precision function enzoHydrostaticTemperature(self,coordinates) result(temperature)
    !!{
    Return the temperature at the specified {\normalfont \ttfamily coordinates} in an Enzo hydrostatic kinematic distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : meanAtomicMassPrimordial, gravitationalConstant_internal
    use :: Numerical_Constants_Atomic      , only : massHydrogenAtom
    use :: Numerical_Constants_Physical    , only : boltzmannsConstant
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    class           (kinematicsDistributionEnzoHydrostatic), intent(inout) :: self
    class           (coordinate                           ), intent(in   ) :: coordinates
    double precision                                                       :: massEnclosed
    
    if (coordinates%rSpherical() == 0.0d0) then
       temperature =temperatureMinimum
    else
       massEnclosed=self%massDistribution_%massEnclosedBySphere(coordinates%rSpherical())
       temperature =max(                                    &
            &           +kilo                          **2  &
            &           *gravitationalConstant_internal     &
            &           *massEnclosed                       &
            &           *meanAtomicMassPrimordial           &
            &           *massHydrogenAtom                   &
            &           /3.0d0                              &
            &           /boltzmannsConstant                 &
            &           /coordinates%rSpherical()         , &
            &           temperatureMinimum                  &
            &          )
    end if
    return
  end function enzoHydrostaticTemperature

  double precision function enzoHydrostaticTemperatureGradientLogarithmic(self,coordinates) result(temperatureGradientLogarithmic)
    !!{
    Return the logarithmic gradient of the temperature at the specified {\normalfont \ttfamily coordinates} in an Enzo hydrostatic kinematic distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (kinematicsDistributionEnzoHydrostatic), intent(inout) :: self
    class           (coordinate                           ), intent(in   ) :: coordinates
    double precision                                                       :: massEnclosed, density
    
    if (self%temperature(coordinates) <= temperatureMinimum) then
       temperatureGradientLogarithmic=temperatureMinimum
    else
       density                   =self%massDistribution_%density             (coordinates             )
       massEnclosed              =self%massDistribution_%massEnclosedBySphere(coordinates%rSpherical())
       temperatureGradientLogarithmic=+            4.0d0             &
            &                         *            Pi                &
            &                         *coordinates%rSpherical  ()**3 &
            &                         *            density           &
            &                         /            massEnclosed      &
            &                         -            1.0d0
    end if
    return
  end function enzoHydrostaticTemperatureGradientLogarithmic
