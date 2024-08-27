!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  <massDistributionHeating name="massDistributionHeatingMonotonic">
    <description>
      A mass distribution heating class which takes another heating source and enforces monotonic heating energy perturbation.
    </description>
  </massDistributionHeating>
  !!]
  type, extends(massDistributionHeatingClass) :: massDistributionHeatingMonotonic
     !!{
     Implementation of a mass distribution heating class which takes another heating source and enforces monotonic heating energy perturbation.
     !!}
     private
     class           (massDistributionHeatingClass), pointer :: massDistributionHeating_ => null()
     type            (rootFinder                  )          :: finder
     double precision                                        :: radiusShellCrossing               , energyPerturbationShellCrossing
   contains
     !![
     <methods>
       <method description="Return true if the no shell crossing assumption is valid locally." method="noShellCrossingIsValid"    />
       <method description="Compute the radius where shell crossing happens."                  method="computeRadiusShellCrossing"/>
     </methods>
     !!]
     procedure :: specificEnergy                 => monotonicSpecificEnergy
     procedure :: specificEnergyGradient         => monotonicSpecificEnergyGradient
     procedure :: specificEnergyIsEveryWhereZero => monotonicSpecificEnergyIsEverywhereZero
     procedure :: noShellCrossingIsValid         => monotonicNoShellCrossingIsValid
     procedure :: computeRadiusShellCrossing     => monotonicComputeRadiusShellCrossing
  end type massDistributionHeatingMonotonic

  interface massDistributionHeatingMonotonic
     !!{
     Constructors for the {\normalfont \ttfamily monotonic} mass distribution class.
     !!}
     module procedure monotonicConstructorParameters
     module procedure monotonicConstructorInternal
  end interface massDistributionHeatingMonotonic

  ! Global variables used in root solving.
  type (massDistributionHeatingMonotonic), pointer :: self_
  class(massDistributionClass           ), pointer :: massDistribution__
  !$omp threadprivate(self_,massDistribution__)

contains

  function monotonicConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily monotonic} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (massDistributionHeatingMonotonic)                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(massDistributionHeatingClass    ), pointer       :: massDistributionHeating_
    
    !![
    <objectBuilder class="massDistributionHeating" name="massDistributionHeating_" source="parameters"/>
    !!]
    self=massDistributionHeatingMonotonic(massDistributionHeating_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistributionHeating_"/>
    !!]
    return
  end function monotonicConstructorParameters
  
  function monotonicConstructorInternal(massDistributionHeating_) result(self)
    !!{
    Constructor for ``monotonic'' mass distribution heating class.
    !!}
    implicit none
    type            (massDistributionHeatingMonotonic)                        :: self
    class           (massDistributionHeatingClass    ), target, intent(in   ) :: massDistributionHeating_
    double precision                                  , parameter             :: toleranceAbsolute       =0.0d0, toleranceRelative=1.0d-6
    !![
    <constructorAssign variables="*massDistributionHeating_"/>
    !!]

    self%radiusShellCrossing            =-1.0d0
    self%energyPerturbationShellCrossing=-1.0d0
    self%finder                         =rootFinder(                                                    &
         &                                          rootFunction     =monotonicRadiusShellCrossingRoot, &
         &                                          toleranceAbsolute=toleranceAbsolute               , &
         &                                          toleranceRelative=toleranceRelative                 &
         &                                         ) 
    return
  end function monotonicConstructorInternal

  subroutine monotonicDestructor(self)
    !!{
    Destructor for the ``monotonic'' mass distribution heating class.
    !!}
    implicit none
    type(massDistributionHeatingMonotonic), intent(inout) :: self

    !![
    <objectDestructor name="self%massDistributionHeating_"/>
    !!]
    return
  end subroutine monotonicDestructor

  double precision function monotonicSpecificEnergy(self,radius,massDistribution_) result(energySpecific)
    !!{
    Compute the specific energy in a monotonicly-heated mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (massDistributionHeatingMonotonic), intent(inout) :: self
    double precision                                  , intent(in   ) :: radius
    class           (massDistributionClass           ), intent(inout) :: massDistribution_

    if (self%noShellCrossingIsValid(radius,massDistribution_)) then
       energySpecific=self%massDistributionHeating_%specificEnergy(                   &
            &                                                      radius           , &
            &                                                      massDistribution_  &
            &                                                     )
    else
       if (self%energyPerturbationShellCrossing < 0.0d0)            &
            call self%computeRadiusShellCrossing(                   &
            &                                    radius           , &
            &                                    massDistribution_  &
            &                                   )
       energySpecific=+self%energyPerturbationShellCrossing                         &
            &         *0.5d0                                                        &
            &         *gravitationalConstantGalacticus                              &
            &         *massDistribution_              %massEnclosedBySphere(radius) &
            &         /radius
    end if
    return
  end function monotonicSpecificEnergy

  double precision function monotonicSpecificEnergyGradient(self,radius,massDistribution_) result(energySpecificGradient)
    !!{
    Returns the gradient of the specific energy of heating.
    !!}
    use :: Coordinates                     , only : coordinateSpherical            , assignment(=)
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (massDistributionHeatingMonotonic), intent(inout) :: self
    double precision                                  , intent(in   ) :: radius
    class           (massDistributionClass           ), intent(inout) :: massDistribution_
    type            (coordinateSpherical             )                :: coordinates

    if (self%noShellCrossingIsValid(radius,massDistribution_)) then
       energySpecificGradient=self%massDistributionHeating_%specificEnergyGradient(                   &
            &                                                                      radius           , &
            &                                                                      massDistribution_  &
            &                                                                     )
    else
       if (self%energyPerturbationShellCrossing < 0.0d0)              &
            & call self%computeRadiusShellCrossing(                   &
            &                                      radius           , &
            &                                      massDistribution_  &
            &                                     )
       coordinates                    =[radius,0.0d0,0.0d0]
       energySpecificGradient=+self%energyPerturbationShellCrossing                  &
            &                 *0.5d0                                                 &
            &                 *gravitationalConstantGalacticus                       &
            &                 *(                                                     &
            &                   +4.0d0                                               &
            &                   *Pi                                                  &
            &                   *                                        radius      &
            &                   *massDistribution_%density             (coordinates) &
            &                   -massDistribution_%massEnclosedBySphere(radius     ) &
            &                   /radius**2                                           &
            &                  )
    end if
    return
  end function monotonicSpecificEnergyGradient

  logical function monotonicSpecificEnergyIsEverywhereZero(self) result(energySpecificIsEverywhereZero)
    !!{
    Returns true if the specific energy is everywhere zero.
    !!}
    implicit none
    class(massDistributionHeatingMonotonic), intent(inout) :: self

    energySpecificIsEverywhereZero=self%massDistributionHeating_%specificEnergyIsEverywhereZero()
    return
  end function monotonicSpecificEnergyIsEverywhereZero

  logical function monotonicNoShellCrossingIsValid(self,radius,massDistribution_)
    !!{
    Determines if the no shell crossing assumption is valid.
    !!}
    use :: Coordinates             , only : coordinateSpherical, assignment(=)
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionHeatingMonotonic), intent(inout) :: self
    class           (massDistributionClass           ), intent(inout) :: massDistribution_
    double precision                                  , intent(in   ) :: radius
    double precision                                                  :: massEnclosed
    type            (coordinateSpherical             )                :: coordinates

    massEnclosed                      = massDistribution_%massEnclosedBySphere(radius)
    if (massEnclosed > 0.0d0) then
       coordinates                    =[radius,0.0d0,0.0d0]
       monotonicNoShellCrossingIsValid=+self%massDistributionHeating_%specificEnergyGradient(                   &
            &                                                                                radius           , &
            &                                                                                massDistribution_  &
            &                                                                               )                   &
            &                          *                                                     radius             &
            &                          +self%massDistributionHeating_%specificEnergy        (                   &
            &                                                                                radius           , &
            &                                                                                massDistribution_  &     
            &                                                                               )                   &
            &                          *(                                                                       &
            &                            +1.0d0                                                                 &
            &                            -4.0d0                                                                 &
            &                            *Pi                                                                    &
            &                            *                                                   radius**3          &
            &                            *massDistribution_          %density               (                   &
            &                                                                                coordinates        &
            &                                                                               )                   &
            &                            /massEnclosed                                                          &
            &                           )                                                                       &
            &                          >=0.0d0
    else
       monotonicNoShellCrossingIsValid=.true.
    end if
    return
  end function monotonicNoShellCrossingIsValid

  subroutine monotonicComputeRadiusShellCrossing(self,radius,massDistribution_)
    !!{
    Determines if the no shell crossing assumption is valid.
    !!}
    use :: Root_Finder                     , only : rangeExpandMultiplicative      , rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (massDistributionHeatingMonotonic), intent(inout), target :: self
    class           (massDistributionClass           ), intent(inout), target :: massDistribution_
    double precision                                  , intent(in   )         :: radius

    if (self%energyPerturbationShellCrossing < 0.0d0) then
       self_              => self
       massDistribution__ => massDistribution_
       call self%finder%rangeExpand(                                                             &
            &                       rangeExpandUpward            =2.0d0                        , &
            &                       rangeExpandDownward          =1.0d0                        , &
            &                       rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
            &                       rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &                       rangeExpandType              =rangeExpandMultiplicative      &
            &                      )
       self%radiusShellCrossing             =+self%finder%find(rootGuess=radius)
       self%energyPerturbationShellCrossing =+self%massDistributionHeating_%specificEnergy      (self%radiusShellCrossing,massDistribution_) &
            &                                /(                                                                                              &
            &                                  +0.5d0                                                                                        &
            &                                  *gravitationalConstantGalacticus                                                              &
            &                                  *massDistribution_          %massEnclosedBySphere(self%radiusShellCrossing                  ) &
            &                                  /                                                 self%radiusShellCrossing                    &
            &                                 )
    end if
    return
  end subroutine monotonicComputeRadiusShellCrossing

  double precision function monotonicRadiusShellCrossingRoot(radius)
    !!{
    Root function used in finding the radius where shell crossing happens.
    !!}
    use :: Coordinates             , only : coordinateSpherical, assignment(=)
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision                     , intent(in   ) :: radius
    double precision                                     :: massEnclosed
    type            (coordinateSpherical)                :: coordinates

    massEnclosed                       = massDistribution__%massEnclosedBySphere(radius)
    if (massEnclosed > 0.0d0) then
       coordinates                     =[radius,0.0d0,0.0d0]
       monotonicRadiusShellCrossingRoot=+self_%massDistributionHeating_%specificEnergyGradient(                    &
            &                                                                                  radius            , &
            &                                                                                  massDistribution__  &
            &                                                                                 )                    &
            &                           *                                                      radius              &
            &                           +self_%massDistributionHeating_%specificEnergy        (                    &
            &                                                                                  radius            , &
            &                                                                                  massDistribution__  &
            &                                                                                 )                    &
            &                           *(                                                                         &
            &                             +1.0d0                                                                   &
            &                             -4.0d0                                                                   &
            &                             *Pi                                                                      &
            &                             *                                                    radius**3           &
            &                             *massDistribution__          %density               (                    &
            &                                                                                  coordinates         &
            &                                                                                 )                    &
            &                             /massEnclosed                                                            &
            &                            )
    else
       monotonicRadiusShellCrossingRoot=0.0d0
    end if
    return
  end function monotonicRadiusShellCrossingRoot
