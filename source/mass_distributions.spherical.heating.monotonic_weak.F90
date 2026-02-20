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

  !+    Contributions to this file made by: Xiaolong Du.

  !!{
  Implements a mass distribution heating class which takes another heating source and enforces monotonic heating energy perturbation.
  !!}

  !![
  <massDistributionHeating name="massDistributionHeatingMonotonicWeak">
    <description>
      A mass distribution heating class which takes another heating source and enforces monotonic heating energy
      perturbation. This class enforces a weaker condition (compared to \refClass{massDistributionHeatingMonotonic}):      
      \begin{equation}
      \frac{\mathrm{d}r_\mathrm{f}}{\mathrm{d}r_\mathrm{i}} > 0,
      \end{equation}
      where $r_\mathrm{i}$ and $r_\mathrm{f}$ are the initial and final radii of the shell respectively.
      
      Note that this condition does not ensure that the gradient of the specific heating energy is continuous through the
      shell-crossing radius. As such, the heated density profile may be discontinuous at this radius also.      

      Using the fact that
      \begin{equation}
      -\frac{\mathrm{G}M}{2 r_\mathrm{f}} = -\frac{\mathrm{G}M}{2 r_\mathrm{i}} + \epsilon(r_\mathrm{i}),
      \end{equation}

      where $\epsilon(r)$ is the specific heating energy as a function of radius, and $M$ is the mass enclosed by the shell, we
      can re-write the above condition as
      \begin{equation}
      \frac{\mathrm{d}\epsilon}{\mathrm{d}r_\mathrm{i}} r_\mathrm{i} - \epsilon(r_\mathrm{i}) \frac{4 \pi r_\mathrm{i}^3 \rho_\mathrm{i}(r_\mathrm{i})}{M} + \frac{\mathrm{G}M}{2r_\mathrm{i}} > \xi \frac{\mathrm{G}M}{r_\mathrm{i}},
      \end{equation}
      where $\rho_\mathrm{i}(r_\mathrm{i})$ is the density in the unheated profile. Here, $\xi$ should equal zero to precisely
      match the criterion for no shell-crossing. However, it is often useful to allow $\xi$ to be a small positive number---this
      avoids getting too close to the boundary of the shell crossing region (where the density can diverge as there is, by
      definition, a caustic in density at this point).
    </description>
  </massDistributionHeating>
  !!]
  type, extends(massDistributionHeatingMonotonic) :: massDistributionHeatingMonotonicWeak
     !!{
     Implementation of a mass distribution heating class which takes another heating source and enforces monotonic heating energy perturbation.
     !!}
     private
     double precision :: toleranceShellCrossing
   contains
     procedure :: noShellCrossingIsValid  => monotonicWeakNoShellCrossingIsValid
     procedure :: radiusShellCrossingRoot => monotonicWeakRadiusShellCrossingRoot
  end type massDistributionHeatingMonotonicWeak

  interface massDistributionHeatingMonotonicWeak
     !!{
     Constructors for the \refClass{massDistributionHeatingMonotonicWeak} mass distribution class.
     !!}
     module procedure monotonicWeakConstructorParameters
     module procedure monotonicWeakConstructorInternal
  end interface massDistributionHeatingMonotonicWeak

contains

  function monotonicWeakConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionHeatingMonotonicWeak} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (massDistributionHeatingMonotonicWeak)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (massDistributionHeatingClass        ), pointer       :: massDistributionHeating_
    double precision                                                      :: toleranceShellCrossing

    !![
    <inputParameter>
      <name>toleranceShellCrossing</name>
      <defaultValue>1.0d-3</defaultValue>
      <source>parameters</source>
      <description>The tolerance adopted in determining if the no-shell-crossing assumption is valid.</description>
    </inputParameter>
    <objectBuilder class="massDistributionHeating" name="massDistributionHeating_" source="parameters"/>
    !!]
    self=massDistributionHeatingMonotonicWeak(toleranceShellCrossing,massDistributionHeating_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistributionHeating_"/>
    !!]
    return
  end function monotonicWeakConstructorParameters
  
  function monotonicWeakConstructorInternal(toleranceShellCrossing,massDistributionHeating_) result(self)
    !!{
    Constructor for the \refClass{massDistributionHeatingMonotonicWeak} mass distribution heating class.
    !!}
    implicit none
    type            (massDistributionHeatingMonotonicWeak)                        :: self
    class           (massDistributionHeatingClass        ), target, intent(in   ) :: massDistributionHeating_
    double precision                                              , intent(in   ) :: toleranceShellCrossing

    self%massDistributionHeatingMonotonic=massDistributionHeatingMonotonic(massDistributionHeating_)
    self%toleranceShellCrossing          =toleranceShellCrossing
    return
  end function monotonicWeakConstructorInternal

  logical function monotonicWeakNoShellCrossingIsValid(self,radius,massDistribution_) result(isValid)
    !!{
    Determines if the no shell crossing assumption is valid.
    !!}
    use :: Coordinates                     , only : coordinateSpherical           , assignment(=)
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionHeatingMonotonicWeak), intent(inout) :: self
    class           (massDistributionClass               ), intent(inout) :: massDistribution_
    double precision                                      , intent(in   ) :: radius
    double precision                                                      :: energySpecific     , energySpecificGradient, &
         &                                                                   energySpecificScale, massEnclosed
    type            (coordinateSpherical                 )                :: coordinates

    massEnclosed             = +massDistribution_%massEnclosedBySphere(radius)
    if (massEnclosed > 0.0d0) then
       coordinates           = [radius,0.0d0,0.0d0]
       energySpecific        = +self%massDistributionHeating_%specificEnergy        (                   &
            &                                                                        radius           , &
            &                                                                        massDistribution_  &
            &                                                                       )
       energySpecificGradient= +self%massDistributionHeating_%specificEnergyGradient(                   &
            &                                                                        radius           , &
            &                                                                        massDistribution_  &     
            &                                                                       )
       energySpecificScale   = +gravitationalConstant_internal &
            &                  *massEnclosed                   &
            &                  /radius
       isValid               = +energySpecificGradient                    &
            &                  *                          radius          &
            &                  -energySpecific                            &
            &                  *4.0d0                                     &
            &                  *Pi                                        &
            &                  *                          radius      **3 &
            &                  *massDistribution_%density(coordinates)    &
            &                  /massEnclosed                              &
            &                  +0.5d0                                     &
            &                  *energySpecificScale                       &
            &                 >                                           &
            &                  +self%toleranceShellCrossing               &
            &                  *energySpecificScale
    else
       isValid=.true.
    end if
    return
  end function monotonicWeakNoShellCrossingIsValid

  double precision function monotonicWeakRadiusShellCrossingRoot(self,radius,massDistribution_) result(root)
    !!{
    Root function used in finding the radius where shell crossing happens.
    !!}
    use :: Coordinates                     , only : coordinateSpherical           , assignment(=)
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionHeatingMonotonicWeak), intent(inout) :: self
    double precision                                      , intent(in   ) :: radius
    class           (massDistributionClass               ), intent(inout) :: massDistribution_
    double precision                                                      :: massEnclosed
    type            (coordinateSpherical)                                 :: coordinates

    massEnclosed  =+massDistribution_%massEnclosedBySphere(radius)
    if (massEnclosed > 0.0d0) then
       coordinates=[radius,0.0d0,0.0d0]
       root       =+self_%massDistributionHeating_%specificEnergyGradient(                   &
            &                                                             radius           , &
            &                                                             massDistribution_  &
            &                                                            )                   &
            &      *                                                      radius             &
            &      -self_%massDistributionHeating_%specificEnergy        (                   &
            &                                                             radius           , &
            &                                                             massDistribution_  &
            &                                                            )                   &
            &      *4.0d0                                                                    &
            &      *Pi                                                                       &
            &      *                                                      radius**3          &
            &      *massDistribution__          %density                 (                   &
            &                                                             coordinates        &
            &                                                            )                   &
            &      /massEnclosed                                                             &
            &      +(                                                                        &
            &        +0.5d0                                                                  &
            &        -self%toleranceShellCrossing                                            &
            &       )                                                                        &
            &      *gravitationalConstant_internal                                           &
            &      *massEnclosed                                                             &
            &      /radius
    else
       root       =+0.0d0
    end if
    return
  end function monotonicWeakRadiusShellCrossingRoot
