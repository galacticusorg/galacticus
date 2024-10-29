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

  !+    Contributions to this file made by: Xiaolong Du.

  !!{
  A dark matter halo profile heating class which takes another heating source and enforces monotonic heating energy perturbation.
  !!}

  use :: Root_Finder, only : rootFinder

  !![
  <darkMatterProfileHeating name="darkMatterProfileHeatingMonotonic">
    <description>
      A dark matter profile heating model which takes another heating source and enforces monotonic heating energy
      perturbation. This is achieved by enforcing the constraint that
      \begin{equation}
      \mathrm{d}/\mathrm{d}r \left[ \frac{\epsilon}{\mathrm{G} M(r) / r} \right] > 0,
      \end{equation}
      where $\epsilon$ is the specific heating energy \citep{du_tidal_2024}. At radii smaller than the shell-crossing radius
      defined by the above condition the specific energy is assumed to be proportional to $\mathrm{G} M(r) / r$, with a smooth
      transition through this radius.
    </description>
  </darkMatterProfileHeating>
  !!]

  type, extends(darkMatterProfileHeatingClass) :: darkMatterProfileHeatingMonotonic
     !!{
     A dark matter profile heating class which takes another heating source and enforces monotonic heating energy perturbation.
     !!}
     private
     class           (darkMatterProfileHeatingClass), pointer :: darkMatterProfileHeating_ => null()
     type            (rootFinder                   )          :: finder
     double precision                                         :: radiusShellCrossing                , energyPerturbationShellCrossing
     integer         (kind_int8                    )          :: lastUniqueID
   contains
     !![
     <methods>
       <method description="Reset memoized calculations."                                      method="calculationReset"          />
       <method description="Return true if the no shell crossing assumption is valid locally." method="noShellCrossingIsValid"    />
       <method description="Compute the radius where shell crossing happens."                  method="computeRadiusShellCrossing"/>
       <method description="Root function used in finding the radius of shell crossing."       method="radiusShellCrossingRoot"   />
     </methods>
     !!]
     final     ::                                   monotonicDestructor
     procedure :: autoHook                       => monotonicAutoHook
     procedure :: calculationReset               => monotonicCalculationReset
     procedure :: specificEnergy                 => monotonicSpecificEnergy
     procedure :: specificEnergyGradient         => monotonicSpecificEnergyGradient
     procedure :: specificEnergyIsEverywhereZero => monotonicSpecificEnergyIsEverywhereZero
     procedure :: noShellCrossingIsValid         => monotonicNoShellCrossingIsValid
     procedure :: computeRadiusShellCrossing     => monotonicComputeRadiusShellCrossing
     procedure :: radiusShellCrossingRoot        => monotonicRadiusShellCrossingRoot
  end type darkMatterProfileHeatingMonotonic

  interface darkMatterProfileHeatingMonotonic
     !!{
     Constructors for the {\normalfont \ttfamily monotonic} dark matter profile heating class.
     !!}
     module procedure monotonicConstructorParameters
     module procedure monotonicConstructorInternal
  end interface darkMatterProfileHeatingMonotonic

  ! Global variables used in root solving.
  type (treeNode                         ), pointer :: node_
  class(darkMatterProfileHeatingMonotonic), pointer :: self_
  class(darkMatterProfileDMOClass        ), pointer :: darkMatterProfileDMO__
  !$omp threadprivate(node_,self_,darkMatterProfileDMO__)

contains

  function monotonicConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily monotonic} dark matter profile heating class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (darkMatterProfileHeatingMonotonic), target        :: self
    type   (inputParameters                  ), intent(inout) :: parameters
    class  (darkMatterProfileHeatingClass    ), pointer       :: darkMatterProfileHeating_

    !![
    <objectBuilder class="darkMatterProfileHeating" name="darkMatterProfileHeating_" source="parameters"/>
    !!]
    self=darkMatterProfileHeatingMonotonic(darkMatterProfileHeating_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileHeating_"/>
    !!]
    return
  end function monotonicConstructorParameters

  function monotonicConstructorInternal(darkMatterProfileHeating_) result(self)
    !!{
    Internal constructor for the ``monotonic'' dark matter profile heating class.
    !!}
    implicit none
    type            (darkMatterProfileHeatingMonotonic)                        :: self
    class           (darkMatterProfileHeatingClass    ), target, intent(in   ) :: darkMatterProfileHeating_
    double precision                                   , parameter             :: toleranceAbsolute        =0.0d0, toleranceRelative=1.0d-6
    !![
    <constructorAssign variables="*darkMatterProfileHeating_"/>
    !!]

    self%radiusShellCrossing            =-1.0d0
    self%energyPerturbationShellCrossing=-1.0d0
    self%lastUniqueID                   =-1_kind_int8
    self%finder                         =rootFinder(                                                     &
         &                                          rootFunction     =monotonicRadiusShellCrossingRoot_, &
         &                                          toleranceAbsolute=toleranceAbsolute                , &
         &                                          toleranceRelative=toleranceRelative                  &
         &                                         )
    return
  end function monotonicConstructorInternal

  subroutine monotonicAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileHeatingMonotonic), intent(inout) :: self

    call calculationResetEvent%attach(self,monotonicCalculationReset,openMPThreadBindingAllLevels,label='darkMatterProfileHeatingMonotonic')
    return 
  end subroutine monotonicAutoHook
    
  subroutine monotonicCalculationReset(self,node,uniqueID)
    !!{ 
    Reset the stored shell crossing radius.
    !!}
    implicit none   
    class  (darkMatterProfileHeatingMonotonic), intent(inout) :: self
    type   (treeNode                         ), intent(inout) :: node
    integer(kind_int8                        ), intent(in   ) :: uniqueID

    self%radiusShellCrossing            =-1.0d0
    self%energyPerturbationShellCrossing=-1.0d0
    self%lastUniqueID                   =uniqueID
    return
  end subroutine monotonicCalculationReset

  subroutine monotonicDestructor(self)
    !!{
    Destructor for the ``monotonic'' dark matter profile heating class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(darkMatterProfileHeatingMonotonic), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileHeating_"/>
    !!]
    if (calculationResetEvent%isAttached(self,monotonicCalculationReset)) call calculationResetEvent%detach(self,monotonicCalculationReset)
    return
  end subroutine monotonicDestructor

  double precision function monotonicSpecificEnergy(self,node,radius,darkMatterProfileDMO_)
    !!{
    Returns the specific energy of heating in the given {\normalfont \ttfamily node}.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileHeatingMonotonic), intent(inout) :: self
    type            (treeNode                         ), intent(inout) :: node
    class           (darkMatterProfileDMOClass        ), intent(inout) :: darkMatterProfileDMO_
    double precision                                   , intent(in   ) :: radius
    integer         (kind_int8                        )                :: uniqueID

    if (self%noShellCrossingIsValid(node,radius,darkMatterProfileDMO_)) then
       monotonicSpecificEnergy=self%darkMatterProfileHeating_%specificEnergy(                       &
            &                                                                node                 , &
            &                                                                radius               , &
            &                                                                darkMatterProfileDMO_  &
            &                                                               )
    else
       uniqueID=node%uniqueID()
       if (uniqueID /= self%lastUniqueID) call self%calculationReset(node,uniqueID)
       if (self%energyPerturbationShellCrossing < 0.0d0) then
          call self%computeRadiusShellCrossing                              (                       &
               &                                                             node                 , &
               &                                                             radius               , &
               &                                                             darkMatterProfileDMO_  &
               &                                                            )
       end if
       monotonicSpecificEnergy=+self%energyPerturbationShellCrossing                      &
            &                  *0.5d0                                                     &
            &                  *gravitationalConstantGalacticus                           &
            &                  *darkMatterProfileDMO_          %enclosedMass(node,radius) &
            &                  /radius
    end if
    return
  end function monotonicSpecificEnergy

  double precision function monotonicSpecificEnergyGradient(self,node,radius,darkMatterProfileDMO_)
    !!{
    Returns the gradient of the specific energy of heating in the given {\normalfont \ttfamily node}.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (darkMatterProfileHeatingMonotonic), intent(inout) :: self
    type            (treeNode                         ), intent(inout) :: node
    class           (darkMatterProfileDMOClass        ), intent(inout) :: darkMatterProfileDMO_
    double precision                                   , intent(in   ) :: radius
    integer         (kind_int8                        )                :: uniqueID

    if (self%noShellCrossingIsValid(node,radius,darkMatterProfileDMO_)) then
       monotonicSpecificEnergyGradient=self%darkMatterProfileHeating_%specificEnergyGradient(                       &
            &                                                                                node                 , &
            &                                                                                radius               , &
            &                                                                                darkMatterProfileDMO_  &
            &                                                                               )
    else
       uniqueID=node%uniqueID()
       if (uniqueID /= self%lastUniqueID) call self%calculationReset(node,uniqueID)
       if (self%energyPerturbationShellCrossing < 0.0d0) then
          call self%computeRadiusShellCrossing                                              (                       &
               &                                                                             node                 , &
               &                                                                             radius               , &
               &                                                                             darkMatterProfileDMO_  &
               &                                                                            )
       end if
       monotonicSpecificEnergyGradient=+self%energyPerturbationShellCrossing              &
            &                          *0.5d0                                             &
            &                          *gravitationalConstantGalacticus                   &
            &                          *(                                                 &
            &                            +4.0d0                                           &
            &                            *Pi                                              &
            &                            *                                        radius  &
            &                            *darkMatterProfileDMO_%density     (node,radius) &
            &                            -darkMatterProfileDMO_%enclosedMass(node,radius) &
            &                            /radius**2                                       &
            &                           )
    end if
    return
  end function monotonicSpecificEnergyGradient

  logical function monotonicSpecificEnergyIsEverywhereZero(self,node,darkMatterProfileDMO_)
    !!{
    Returns true if the specific energy is everywhere zero in the given {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileHeatingMonotonic), intent(inout) :: self
    type (treeNode                         ), intent(inout) :: node
    class(darkMatterProfileDMOClass        ), intent(inout) :: darkMatterProfileDMO_

    monotonicSpecificEnergyIsEverywhereZero=self%darkMatterProfileHeating_%specificEnergyIsEverywhereZero(node,darkMatterProfileDMO_)
    return
  end function monotonicSpecificEnergyIsEverywhereZero

  logical function monotonicNoShellCrossingIsValid(self,node,radius,darkMatterProfileDMO_)
    !!{
    Determines if the no shell crossing assumption is valid.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileHeatingMonotonic), intent(inout) :: self
    type            (treeNode                         ), intent(inout) :: node
    class           (darkMatterProfileDMOClass        ), intent(inout) :: darkMatterProfileDMO_
    double precision                                   , intent(in   ) :: radius
    double precision                                   , parameter     :: toleranceRelative    =1.0d-12
    double precision                                                   :: massEnclosed                 , energySpecificScale   , &
         &                                                                energySpecific               , energySpecificGradient

    massEnclosed                      = darkMatterProfileDMO_%enclosedMass(node,radius)
    if (massEnclosed > 0.0d0) then
       energySpecific                 =+self%darkMatterProfileHeating_%specificEnergy        (                       &
            &                                                                                 node                 , &
            &                                                                                 radius               , &
            &                                                                                 darkMatterProfileDMO_  &     
            &                                                                                )
       energySpecificGradient         =+self%darkMatterProfileHeating_%specificEnergyGradient(                       &
            &                                                                                 node                 , &
            &                                                                                 radius               , &
            &                                                                                 darkMatterProfileDMO_  &
            &                                                                                )
       energySpecificScale            =+gravitationalConstantGalacticus &
            &                          *massEnclosed                    &
            &                          /radius
       monotonicNoShellCrossingIsValid=+energySpecificGradient                     &
            &                          *                                radius     &
            &                          +energySpecific                             &
            &                          *(                                          &
            &                            +1.0d0                                    &
            &                            -4.0d0                                    &
            &                            *Pi                                       &
            &                            *                              radius**3  &
            &                            *darkMatterProfileDMO_%density(           &
            &                                                           node     , &
            &                                                           radius     &
            &                                                          )           &
            &                            /massEnclosed                             &
            &                           )                                          &
            &                          >=                                          &
            &                           -toleranceRelative                         &
            &                           *energySpecificScale
    else
       monotonicNoShellCrossingIsValid=.true.
    end if
    return
  end function monotonicNoShellCrossingIsValid

  subroutine monotonicComputeRadiusShellCrossing(self,node,radius,darkMatterProfileDMO_)
    !!{
    Determines if the no shell crossing assumption is valid.
    !!}
    use :: Root_Finder                     , only : rangeExpandMultiplicative      , rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileHeatingMonotonic), intent(inout), target :: self
    type            (treeNode                         ), intent(inout), target :: node
    class           (darkMatterProfileDMOClass        ), intent(inout), target :: darkMatterProfileDMO_
    double precision                                   , intent(in   )         :: radius
    double precision                                   , parameter             :: radiusSearchMaximum  =10.0d0, factorRadius              =1.1d0
    integer         (kind_int8                        )                        :: uniqueID
    double precision                                                           :: radiusSearch                , radiusShellCrossingMinimum

    uniqueID=node%uniqueID()
    if (uniqueID /= self%lastUniqueID) call self%calculationReset(node,uniqueID)
    if (self%energyPerturbationShellCrossing < 0.0d0) then
       ! Search for the largest radius at which shell-crossing occurs.
       radiusSearch              =radius
       radiusShellCrossingMinimum=radius
       do while (radiusSearch <= radiusSearchMaximum)
          if (.not.self%noShellCrossingIsValid(node,radiusSearch,darkMatterProfileDMO_)) radiusShellCrossingMinimum=radiusSearch
          radiusSearch=factorRadius*radiusSearch
       end do
       ! Seek the exact radius at which shell-crossing first occurs. Use an expansion step matched to that in our prior search
       ! since we know that the root should be within this range.
       self_                  => self
       node_                  => node
       darkMatterProfileDMO__ => darkMatterProfileDMO_
       call self%finder%rangeExpand(                                                             &
            &                       rangeExpandUpward            =1.0d0*factorRadius           , &
            &                       rangeExpandDownward          =1.0d0/factorRadius           , &
            &                       rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
            &                       rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &                       rangeExpandType              =rangeExpandMultiplicative      &
            &                      )
       self%radiusShellCrossing             =+self%finder%find(rootGuess=radiusShellCrossingMinimum)
       ! If we exceeded the radius for which the specific energy is non-zero, back up to the prior radius.
       if (self%darkMatterProfileHeating_%specificEnergy(node,self%radiusShellCrossing,darkMatterProfileDMO_) <= 0.0d0) &
            & self%radiusShellCrossing=radiusShellCrossingMinimum
       self%energyPerturbationShellCrossing =+self%darkMatterProfileHeating_%specificEnergy(node,self%radiusShellCrossing,darkMatterProfileDMO_) &
            &                                /(                                                                                                  &
            &                                  +0.5d0                                                                                            &
            &                                  *gravitationalConstantGalacticus                                                                  &
            &                                  *darkMatterProfileDMO_       %enclosedMass  (node,self%radiusShellCrossing                      ) &
            &                                  /                                                 self%radiusShellCrossing                        &
            &                                 )
    end if
    return
  end subroutine monotonicComputeRadiusShellCrossing

  double precision function monotonicRadiusShellCrossingRoot_(radius) result(root)
    !!{
    Root function used in finding the radius where shell crossing happens.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    root=self_%radiusShellCrossingRoot(node_,radius,darkMatterProfileDMO__)
    return
  end function monotonicRadiusShellCrossingRoot_
  
  double precision function monotonicRadiusShellCrossingRoot(self,node,radius,darkMatterProfileDMO_)
    !!{
    Root function used in finding the radius where shell crossing happens.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileHeatingMonotonic), intent(inout) :: self
    type            (treeNode                         ), intent(inout) :: node
    double precision                                   , intent(in   ) :: radius
    class           (darkMatterProfileDMOClass        ), intent(inout) :: darkMatterProfileDMO_
    double precision                                                   :: massEnclosed

    massEnclosed                       = darkMatterProfileDMO_%enclosedMass(node,radius)
    if (massEnclosed > 0.0d0) then
       monotonicRadiusShellCrossingRoot=+self%darkMatterProfileHeating_%specificEnergyGradient(                       &
            &                                                                                  node                 , &
            &                                                                                  radius               , &
            &                                                                                  darkMatterProfileDMO_  &
            &                                                                                 )                       &
            &                           *                                                      radius                 &
            &                           +self%darkMatterProfileHeating_%specificEnergy        (                       &
            &                                                                                  node                 , &
            &                                                                                  radius               , &
            &                                                                                  darkMatterProfileDMO_  &
            &                                                                                 )                       &
            &                           *(                                                                            &
            &                             +1.0d0                                                                      &
            &                             -4.0d0                                                                      &
            &                             *Pi                                                                         &
            &                             *                                                    radius**3              &
            &                             *darkMatterProfileDMO_        %density              (                       &
            &                                                                                  node                 , &
            &                                                                                  radius                 &
            &                                                                                 )                       &
            &                             /massEnclosed                                                               &
            &                            )
    else
       monotonicRadiusShellCrossingRoot=+0.0d0
    end if
    return
  end function monotonicRadiusShellCrossingRoot
