!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
  A dark matter halo profile heating class which takes another heating source and enforces monotonic heating energy perturbation.
  !!}

  !![
  <darkMatterProfileHeating name="darkMatterProfileHeatingMonotonicWeak">
    <description>
      A dark matter profile heating model which takes another heating source and enforces monotonic heating energy
      perturbation. This classes enforces a weaker condition (compared to \refClass{darkMatterProfileHeatingMonotonic}):      
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
  </darkMatterProfileHeating>
  !!]
  type, extends(darkMatterProfileHeatingMonotonic) :: darkMatterProfileHeatingMonotonicWeak
     !!{
     A dark matter profile heating class which takes another heating source and enforces monotonic heating energy perturbation.
     !!}
     private
     double precision :: toleranceShellCrossing
   contains
     procedure :: noShellCrossingIsValid  => monotonicWeakNoShellCrossingIsValid
     procedure :: radiusShellCrossingRoot => monotonicWeakRadiusShellCrossingRoot
  end type darkMatterProfileHeatingMonotonicWeak

  interface darkMatterProfileHeatingMonotonicWeak
     !!{
     Constructors for the {\normalfont \ttfamily monotonicWeak} dark matter profile heating class.
     !!}
     module procedure monotonicWeakConstructorParameters
     module procedure monotonicWeakConstructorInternal
  end interface darkMatterProfileHeatingMonotonicWeak

contains

  function monotonicWeakConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily monotonicWeak} dark matter profile heating class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileHeatingMonotonicWeak), target        :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (darkMatterProfileHeatingClass        ), pointer       :: darkMatterProfileHeating_
    double precision                                                       :: toleranceShellCrossing
    
    !![
    <inputParameter>
      <name>toleranceShellCrossing</name>
      <defaultValue>1.0d-3</defaultValue>
      <source>parameters</source>
      <description>The tolerance adopted in determining if the no-shell-crossing assumption is valid.</description>
    </inputParameter>
    <objectBuilder class="darkMatterProfileHeating" name="darkMatterProfileHeating_" source="parameters"/>
    !!]
    self=darkMatterProfileHeatingMonotonicWeak(toleranceShellCrossing,darkMatterProfileHeating_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileHeating_"/>
    !!]
    return
  end function monotonicWeakConstructorParameters

  function monotonicWeakConstructorInternal(toleranceShellCrossing,darkMatterProfileHeating_) result(self)
    !!{
    Internal constructor for the ``monotonicWeak'' dark matter profile heating class.
    !!}
    implicit none
    type            (darkMatterProfileHeatingMonotonicWeak)                        :: self
    class           (darkMatterProfileHeatingClass        ), target, intent(in   ) :: darkMatterProfileHeating_
    double precision                                               , intent(in   ) :: toleranceShellCrossing

    self%darkMatterProfileHeatingMonotonic=darkMatterProfileHeatingMonotonic(darkMatterProfileHeating_)
    self%toleranceShellCrossing           =toleranceShellCrossing
    return
  end function monotonicWeakConstructorInternal

  logical function monotonicWeakNoShellCrossingIsValid(self,node,radius,darkMatterProfileDMO_) result(shellCrossingIsValid)
    !!{
    Determines if the no shell crossing assumption is valid.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileHeatingMonotonicWeak), intent(inout) :: self
    type            (treeNode                             ), intent(inout) :: node
    class           (darkMatterProfileDMOClass            ), intent(inout) :: darkMatterProfileDMO_
    double precision                                       , intent(in   ) :: radius
    double precision                                                       :: energySpecific       , energySpecificGradient, &
         &                                                                    energySpecificScale  , massEnclosed

    massEnclosed             = darkMatterProfileDMO_%enclosedMass(node,radius)
    if (massEnclosed > 0.0d0) then
       energySpecific        = +self%darkMatterProfileHeating_%specificEnergy        (                       &
            &                                                                         node                 , &
            &                                                                         radius               , &
            &                                                                         darkMatterProfileDMO_  &     
            &                                                                        )
       energySpecificGradient= +self%darkMatterProfileHeating_%specificEnergyGradient(                       &
            &                                                                         node                 , &
            &                                                                         radius               , &
            &                                                                         darkMatterProfileDMO_  &
            &                                                                        )
       energySpecificScale   = +gravitationalConstantGalacticus &
            &                  *massEnclosed                    &
            &                  /radius
       shellCrossingIsValid  = +energySpecificGradient                   &
            &                  *                                radius   &
            &                  -energySpecific                           &
            &                  *4.0d0                                    &
            &                  *Pi                                       &
            &                  *                              radius**3  &
            &                  *darkMatterProfileDMO_%density(           &
            &                                                 node     , &
            &                                                 radius     &
            &                                                )           &
            &                  /massEnclosed                             &
            &                  +0.5d0                                    &
            &                  *energySpecificScale                      &
            &                 >                                          &
            &                  +self%toleranceShellCrossing              &
            &                  *energySpecificScale
    else
       shellCrossingIsValid=.true.
    end if
    return
  end function monotonicWeakNoShellCrossingIsValid

  double precision function monotonicWeakRadiusShellCrossingRoot(self,node,radius,darkMatterProfileDMO_) result(root)
    !!{
    Root function used in finding the radius where shell crossing happens.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileHeatingMonotonicWeak), intent(inout) :: self
    type            (treeNode                             ), intent(inout) :: node
    double precision                                       , intent(in   ) :: radius
    class           (darkMatterProfileDMOClass            ), intent(inout) :: darkMatterProfileDMO_
    double precision                                                       :: massEnclosed

    massEnclosed= darkMatterProfileDMO_%enclosedMass(node,radius)
    if (massEnclosed > 0.0d0) then
       root   =+self%darkMatterProfileHeating_%specificEnergyGradient(                        &
            &                                                          node                 , &
            &                                                          radius               , &
            &                                                          darkMatterProfileDMO_  &
            &                                                         )                       &
            &  *                                                       radius                 &
            &  -self%darkMatterProfileHeating_%specificEnergy         (                       &
            &                                                          node                 , &
            &                                                          radius               , &
            &                                                          darkMatterProfileDMO_  &
            &                                                         )                       &
            &  *4.0d0                                                                         &
            &  *Pi                                                                            &
            &  *                                                       radius**3              &
            &  *      darkMatterProfileDMO_       %density            (                       &
            &                                                          node                 , &
            &                                                          radius                 &
            &                                                         )                       &
            &  /massEnclosed                                                                  &
            &  +(                                                                             &
            &    +0.5d0                                                                       &
            &    -self%toleranceShellCrossing                                                 &
            &   )                                                                             &
            &  *gravitationalConstantGalacticus                                               &
            &  *massEnclosed                                                                  &
            &  /radius
    else
       root   =+0.0d0
    end if
    return
  end function monotonicWeakRadiusShellCrossingRoot
