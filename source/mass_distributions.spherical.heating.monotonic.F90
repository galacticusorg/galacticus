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
  <massDistributionHeating name="massDistributionHeatingMonotonic">
    <description>
      A mass distribution heating class which takes another heating source and enforces monotonic heating energy
      perturbation. This is achieved by enforcing the constraint that      
      \begin{equation}
      \mathrm{d}/\mathrm{d}r \left[ \frac{\epsilon}{\mathrm{G} M(r) / r} \right] > 0,
      \end{equation}
      where $\epsilon$ is the specific heating energy \citep{du_tidal_2024}. At radii smaller than the shell-crossing radius
      defined by the above condition the specific energy is assumed to be proportional to $\mathrm{G} M(r) / r$, with a smooth
      transition through this radius.
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
     double precision                                        :: radiusShellCrossing               , energyPerturbationShellCrossing, &
          &                                                      radiusCheckedMinimum
     logical                                                 :: shellCrossingAllRadii
   contains
     !![
     <methods>
       <method description="Return true if the no shell crossing assumption is valid locally." method="noShellCrossingIsValid"    />
       <method description="Compute the radius where shell crossing happens."                  method="computeRadiusShellCrossing"/>
       <method description="Root function used in finding the radius of shell crossing."       method="radiusShellCrossingRoot"   />
     </methods>
     !!]
     final     ::                                   monotonicDestructor
     procedure :: specificEnergy                 => monotonicSpecificEnergy
     procedure :: specificEnergyGradient         => monotonicSpecificEnergyGradient
     procedure :: specificEnergyIsEveryWhereZero => monotonicSpecificEnergyIsEverywhereZero
     procedure :: noShellCrossingIsValid         => monotonicNoShellCrossingIsValid
     procedure :: computeRadiusShellCrossing     => monotonicComputeRadiusShellCrossing
     procedure :: radiusShellCrossingRoot        => monotonicRadiusShellCrossingRoot
  end type massDistributionHeatingMonotonic

  interface massDistributionHeatingMonotonic
     !!{
     Constructors for the \refClass{massDistributionHeatingMonotonic} mass distribution class.
     !!}
     module procedure monotonicConstructorParameters
     module procedure monotonicConstructorInternal
  end interface massDistributionHeatingMonotonic

  ! Global variables used in root solving.
  class(massDistributionHeatingMonotonic), pointer :: self_
  class(massDistributionClass           ), pointer :: massDistribution__
  !$omp threadprivate(self_,massDistribution__)

contains

  function monotonicConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionHeatingMonotonic} mass distribution class which builds the object from a parameter
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
    Constructor for the \refClass{massDistributionHeatingMonotonic} mass distribution heating class.
    !!}
    implicit none
    type            (massDistributionHeatingMonotonic)                        :: self
    class           (massDistributionHeatingClass    ), target, intent(in   ) :: massDistributionHeating_
    double precision                                  , parameter             :: toleranceAbsolute       =0.0d0, toleranceRelative=1.0d-6
    !![
    <constructorAssign variables="*massDistributionHeating_"/>
    !!]

    self%radiusShellCrossing            =-huge(0.0d0)
    self%radiusCheckedMinimum           =+huge(0.0d0)
    self%energyPerturbationShellCrossing=-huge(0.0d0)
    self%shellCrossingAllRadii          =.false.
    self%finder                         =rootFinder(                                                     &
         &                                          rootFunction     =monotonicRadiusShellCrossingRoot_, &
         &                                          toleranceAbsolute=toleranceAbsolute                , &
         &                                          toleranceRelative=toleranceRelative                  &
         &                                         ) 
    return
  end function monotonicConstructorInternal

  subroutine monotonicDestructor(self)
    !!{
    Destructor for the \refClass{massDistributionHeatingMonotonic} mass distribution heating class.
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
    Compute the specific energy in a monotonically-heated mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionHeatingMonotonic), intent(inout) :: self
    double precision                                  , intent(in   ) :: radius
    class           (massDistributionClass           ), intent(inout) :: massDistribution_

    call self%computeRadiusShellCrossing(                   &
         &                               radius           , &
         &                               massDistribution_  &
         &                              )
    if (radius > self%radiusShellCrossing .and. .not.self%shellCrossingAllRadii) then
       energySpecific=self%massDistributionHeating_%specificEnergy(                   &
            &                                                      radius           , &
            &                                                      massDistribution_  &
            &                                                     )
    else
       energySpecific=+self%energyPerturbationShellCrossing                        &
            &         *0.5d0                                                       &
            &         *gravitationalConstant_internal                              &
            &         *massDistribution_             %massEnclosedBySphere(radius) &
            &         /radius
    end if
    return
  end function monotonicSpecificEnergy

  double precision function monotonicSpecificEnergyGradient(self,radius,massDistribution_) result(energySpecificGradient)
    !!{
    Returns the gradient of the specific energy of heating.
    !!}
    use :: Coordinates                     , only : coordinateSpherical           , assignment(=)
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (massDistributionHeatingMonotonic), intent(inout) :: self
    double precision                                  , intent(in   ) :: radius
    class           (massDistributionClass           ), intent(inout) :: massDistribution_
    type            (coordinateSpherical             )                :: coordinates

    call self%computeRadiusShellCrossing(                   &
         &                               radius           , &
         &                               massDistribution_  &
         &                              )
    if (radius > self%radiusShellCrossing .and. .not.self%shellCrossingAllRadii) then
       energySpecificGradient=self%massDistributionHeating_%specificEnergyGradient(                   &
            &                                                                      radius           , &
            &                                                                      massDistribution_  &
            &                                                                     )
    else
       coordinates           =[radius,0.0d0,0.0d0]
       energySpecificGradient=+self%energyPerturbationShellCrossing                  &
            &                 *0.5d0                                                 &
            &                 *gravitationalConstant_internal                        &
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
    use :: Coordinates                     , only : coordinateSpherical           , assignment(=)
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionHeatingMonotonic), intent(inout) :: self
    class           (massDistributionClass           ), intent(inout) :: massDistribution_
    double precision                                  , intent(in   ) :: radius
    double precision                                  , parameter     :: toleranceRelative=1.0d-12
    double precision                                                  :: massEnclosed             , energySpecificScale   , &
         &                                                               energySpecific           , energySpecificGradient
    type            (coordinateSpherical             )                :: coordinates

    massEnclosed                      = massDistribution_%massEnclosedBySphere(radius)
    if (massEnclosed > 0.0d0) then
       coordinates                    =[radius,0.0d0,0.0d0]
       energySpecific                 =+self%massDistributionHeating_%specificEnergyGradient(                   &
            &                                                                                radius           , &
            &                                                                                massDistribution_  &
            &                                                                               )
       energySpecificGradient         =+self%massDistributionHeating_%specificEnergy        (                   &
            &                                                                                radius           , &
            &                                                                                massDistribution_  &     
            &                                                                               )
       energySpecificScale            =+gravitationalConstant_internal &
            &                          *massEnclosed                   &
            &                          /radius
       monotonicNoShellCrossingIsValid=+energySpecific                              &
            &                          *                            radius          &
            &                          +energySpecificGradient                      &
            &                          *(                                           &
            &                            +1.0d0                                     &
            &                            -4.0d0                                     &
            &                            *Pi                                        &
            &                            *                          radius      **3 &
            &                            *massDistribution_%density(coordinates)    &
            &                            /massEnclosed                              &
            &                           )                                           &
            &                          >=                                           &
            &                           -toleranceRelative                          &
            &                           *energySpecificScale
    else
       monotonicNoShellCrossingIsValid=.true.
    end if
    return
  end function monotonicNoShellCrossingIsValid

  subroutine monotonicComputeRadiusShellCrossing(self,radius,massDistribution_)
    !!{
    Determines if the no shell crossing assumption is valid.
    !!}
    use :: Root_Finder                     , only : rangeExpandMultiplicative     , rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionHeatingMonotonic), intent(inout), target :: self
    class           (massDistributionClass           ), intent(inout), target :: massDistribution_
    double precision                                  , intent(in   )         :: radius
    double precision                                  , parameter             :: radiusSearchMaximum  =10.0d0, factorRadius              =1.1d0
    double precision                                                          :: radiusSearch                , radiusShellCrossingMinimum
    logical                                                                   :: foundShellCrossing          , isValidPrevious                 , &
         &                                                                       isValid

    if (self%energyPerturbationShellCrossing < 0.0d0) then
       ! Search for the largest radius at which shell-crossing occurs.
       radiusSearch              =radius
       radiusShellCrossingMinimum=radius
       foundShellCrossing        =.false.
       isValidPrevious           =.true.
       do while (radiusSearch <= radiusSearchMaximum)
          isValid=self%noShellCrossingIsValid(radiusSearch,massDistribution_)
          if (isValid .and. .not.isValidPrevious) foundShellCrossing=.true.
          if (.not.isValid) radiusShellCrossingMinimum=radiusSearch
          isValidPrevious= isValid
          radiusSearch   =+factorRadius &
               &          *radiusSearch
       end do
       ! Determine if a the shell crossing radius was found.
       if (foundShellCrossing .and. isValid) then
          ! Seek the exact radius at which shell-crossing first occurs. Use an expansion step matched to that in our prior search
          ! since we know that the root should be within this range.
          self_              => self
          massDistribution__ => massDistribution_
          call self%finder%rangeExpand(                                                             &
               &                       rangeExpandUpward            =1.0d0*factorRadius           , &
               &                       rangeExpandDownward          =1.0d0/factorRadius           , &
               &                       rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
               &                       rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
               &                       rangeExpandType              =rangeExpandMultiplicative      &
               &                      )
          self%radiusShellCrossing  =+self%finder%find(rootGuess=radiusShellCrossingMinimum)
          self%shellCrossingAllRadii=.false.
          ! If we exceeded the radius for which the specific energy is non-zero, back up to the prior radius.
          if (self%massDistributionHeating_%specificEnergy(self%radiusShellCrossing,massDistribution_) <= 0.0d0) &
               & self%radiusShellCrossing=radiusShellCrossingMinimum
       else if (foundShellCrossing) then
          self%radiusShellCrossing  =+radiusSearchMaximum
          self%shellCrossingAllRadii=.true.
       else
          self%radiusShellCrossing  =-huge(0.0d0)
          self%shellCrossingAllRadii=.false.
       end if
       if (self%radiusShellCrossing > 0.0d0) then
          self%energyPerturbationShellCrossing =+self%massDistributionHeating_%specificEnergy      (self%radiusShellCrossing,massDistribution_) &
               &                                /(                                                                                              &
               &                                  +0.5d0                                                                                        &
               &                                  *gravitationalConstant_internal                                                               &
               &                                  *massDistribution_          %massEnclosedBySphere(self%radiusShellCrossing                  ) &
               &                                  /                                                 self%radiusShellCrossing                    &
               &                                 )
       else
          self%energyPerturbationShellCrossing=-huge(0.0d0)
       end if
    end if
    return
  end subroutine monotonicComputeRadiusShellCrossing

  double precision function monotonicRadiusShellCrossingRoot_(radius) result(root)
    !!{
    Root function used in finding the radius where shell crossing happens.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    root=self_%radiusShellCrossingRoot(radius,massDistribution__)
    return
  end function monotonicRadiusShellCrossingRoot_
  
  double precision function monotonicRadiusShellCrossingRoot(self,radius,massDistribution_)
    !!{
    Root function used in finding the radius where shell crossing happens.
    !!}
    use :: Coordinates             , only : coordinateSpherical, assignment(=)
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionHeatingMonotonic), intent(inout) :: self
    double precision                                  , intent(in   ) :: radius
    class           (massDistributionClass           ), intent(inout) :: massDistribution_
    double precision                                                  :: massEnclosed
    type            (coordinateSpherical)                             :: coordinates

    massEnclosed                       = massDistribution_%massEnclosedBySphere(radius)
    if (massEnclosed > 0.0d0) then
       coordinates                     =[radius,0.0d0,0.0d0]
       monotonicRadiusShellCrossingRoot=+self_%massDistributionHeating_%specificEnergyGradient(                   &
            &                                                                                  radius           , &
            &                                                                                  massDistribution_  &
            &                                                                                 )                   &
            &                           *                                                      radius             &
            &                           +self_%massDistributionHeating_%specificEnergy        (                   &
            &                                                                                  radius           , &
            &                                                                                  massDistribution_  &
            &                                                                                 )                   &
            &                           *(                                                                        &
            &                             +1.0d0                                                                  &
            &                             -4.0d0                                                                  &
            &                             *Pi                                                                     &
            &                             *                                                    radius**3          &
            &                             *massDistribution__          %density               (                   &
            &                                                                                  coordinates        &
            &                                                                                 )                   &
            &                             /massEnclosed                                                           &
            &                            )
    else
       monotonicRadiusShellCrossingRoot=0.0d0
    end if
    return
  end function monotonicRadiusShellCrossingRoot
