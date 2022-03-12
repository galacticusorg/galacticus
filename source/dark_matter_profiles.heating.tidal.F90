!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  A dark matter halo profile heating class which accounts for heating from tidal shocking.
  !!}

  use :: Kind_Numbers, only : kind_int8

  !![
  <darkMatterProfileHeating name="darkMatterProfileHeatingTidal">
   <description>
    A dark matter profile heating model which accounts for heating due to tidal shocking. The model follows the general
    approach of \cite{gnedin_tidal_1999}. The change in the specific energy of particles at radius $r$ in a halo is gievn by
    $\Delta \epsilon = \Delta \epsilon_1 + \Delta \epsilon_2$, where $\Delta \epsilon_1$, and $\Delta \epsilon_2$ are the first
    and second order perturbations respectively. The first order term is given by $\Delta \epsilon_1 = Q r^2$ where $Q$ is the
    tidal tensor integrated along the orbital path (see, for example, \citealt{taylor_dynamics_2001}), while the second order
    term is given by $\Delta \epsilon_2 = (2/3) f \sigma_\mathrm{rms} (1+\chi_\mathrm{r,v}) \sqrt{\Delta \epsilon_1}$
    \citep[][eqn.~20, see also \protect\citealt{gnedin_self-consistent_1999}; eqn.~18a,b]{gnedin_tidal_1999}. For the particle
    velocity dispersion, $v_\mathrm{rms}$, we use $\sqrt{3} \sigma_\mathrm{r}(r)$, the radial velocity dispersion in the dark
    matter profile scaled to the total velocity dispersion assuming an isotropic velocity distribution. The position-velocity
    correlation function, $\chi_\mathrm{r,v}$, is taken to be a constant given by the parameter {\normalfont \ttfamily
    [correlationVelocityRadius]}. The coefficient, $f=${\normalfont \ttfamily [coefficientSecondOrder]} is introduced to allow
    some freedom to adjust the contribution of the second order term. It is degenerate with the value of $\chi_\mathrm{r,v}$
    but is introduced to allow for possible future promotion of $\chi_\mathrm{r,v}$ from a constant to a function of the dark
    matter profile potential \citep[see, for example,][appendix~B]{gnedin_self-consistent_1999}.
   </description>
  </darkMatterProfileHeating>
  !!]
  type, extends(darkMatterProfileHeatingClass) :: darkMatterProfileHeatingTidal
     !!{
     A dark matter profile heating class which accounts for heating due to tidal shocking.
     !!}
     private
     double precision            :: specificEnergyOverRadiusSquared_, specificEnergyOverRadiusSquaredParent_, &
          &                         correlationVelocityRadius       , coefficientSecondOrder
     integer         (kind_int8) :: lastUniqueID                    , parentUniqueID
  contains
     !![
     <methods>
       <method description="Reset memoized calculations."                             method="calculationReset"               />
       <method description="Compute $Q = E / r^2$."                                   method="specificEnergyOverRadiusSquared"/>
       <method description="Compute the first and second order energy perturbations." method="specificEnergyTerms"            />
     </methods>
     !!]
     final     ::                                    tidalDestructor
     procedure :: autoHook                        => tidalAutoHook
     procedure :: calculationReset                => tidalCalculationReset
     procedure :: specificEnergy                  => tidalSpecificEnergy
     procedure :: specificEnergyGradient          => tidalSpecificEnergyGradient
     procedure :: specificEnergyIsEverywhereZero  => tidalSpecificEnergyIsEverywhereZero
     procedure :: specificEnergyOverRadiusSquared => tidalSpecificEnergyOverRadiusSquared
     procedure :: specificEnergyTerms             => tidalSpecificEnergyTerms
  end type darkMatterProfileHeatingTidal

  interface darkMatterProfileHeatingTidal
     !!{
     Constructors for the {\normalfont \ttfamily tidal} dark matter profile heating class.
     !!}
     module procedure tidalConstructorParameters
     module procedure tidalConstructorInternal
  end interface darkMatterProfileHeatingTidal

contains

  function tidalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily tidal} dark matter profile heating scales class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (darkMatterProfileHeatingTidal), target        :: self
    type            (inputParameters              ), intent(inout) :: parameters
    double precision                                               :: coefficientSecondOrder, correlationVelocityRadius
    
    !![
    <inputParameter>
      <name>coefficientSecondOrder</name>
      <defaultValue>0.0d0</defaultValue>
      <source>parameters</source>
      <description>A multiplicative coefficient for the second-order heating term.</description>
    </inputParameter>
    <inputParameter>
      <name>correlationVelocityRadius</name>
      <defaultValue>-1.0d0</defaultValue>
      <source>parameters</source>
      <description>The velocity-position correlation function, $\chi_\mathrm{r,v}$, as defined by \cite[][eqn.~B1]{gnedin_self-consistent_1999} which controls the strength of the second order heating term.</description>
    </inputParameter>
    !!]
    self=darkMatterProfileHeatingTidal(coefficientSecondOrder,correlationVelocityRadius)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function tidalConstructorParameters

  function tidalConstructorInternal(coefficientSecondOrder,correlationVelocityRadius) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily tidal} dark matter profile heating scales class.
    !!}
    implicit none
    type            (darkMatterProfileHeatingTidal)                :: self
    double precision                               , intent(in   ) :: coefficientSecondOrder, correlationVelocityRadius
    !![
    <constructorAssign variables="coefficientSecondOrder, correlationVelocityRadius"/>
    !!]
    
    self%specificEnergyOverRadiusSquared_      =-1.0d0
    self%specificEnergyOverRadiusSquaredParent_=-1.0d0
    self%lastUniqueID                          =-1_kind_int8
    self%parentUniqueID                        =-1_kind_int8
    return
  end function tidalConstructorInternal

  subroutine tidalAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileHeatingTidal), intent(inout) :: self

    call calculationResetEvent%attach(self,tidalCalculationReset,openMPThreadBindingAllLevels)
    return
  end subroutine tidalAutoHook

  subroutine tidalCalculationReset(self,node)
    !!{
    Reset the stored tidal radii.
    !!}
    implicit none
    class(darkMatterProfileHeatingTidal), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node

    self   %specificEnergyOverRadiusSquared_      =-1.0d0
    self   %specificEnergyOverRadiusSquaredParent_=-1.0d0
    self   %lastUniqueID                          =node       %uniqueID()
    if (associated(node%parent)) then
       self%parentUniqueID                        =node%parent%uniqueID()
    else
       self%parentUniqueID                        =-1_kind_int8
    end if
    return
  end subroutine tidalCalculationReset

  subroutine tidalDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily tidal} dark matter profile heating class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(darkMatterProfileHeatingTidal), intent(inout) :: self

    if (calculationResetEvent%isAttached(self,tidalCalculationReset)) call calculationResetEvent%detach(self,tidalCalculationReset)
    return
  end subroutine tidalDestructor

  double precision function tidalSpecificEnergy(self,node,darkMatterProfileDMO_,radius)
    !!{
    Returns the specific energy of heating in the given {\normalfont \ttfamily node}.
    !!}
    implicit none
    class           (darkMatterProfileHeatingTidal), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    class           (darkMatterProfileDMOClass    ), intent(inout) :: darkMatterProfileDMO_
    double precision                               , intent(in   ) :: radius
    double precision                                               :: energyPerturbationFirstOrder, energyPerturbationSecondOrder

    call self%specificEnergyTerms(node,darkMatterProfileDMO_,radius,energyPerturbationFirstOrder,energyPerturbationSecondOrder)
    tidalSpecificEnergy=+energyPerturbationFirstOrder  &
         &              +energyPerturbationSecondOrder
    return
  end function tidalSpecificEnergy

  double precision function tidalSpecificEnergyGradient(self,node,darkMatterProfileDMO_,radius)
    !!{
    Returns the gradient of the specific energy of heating in the given {\normalfont \ttfamily node}.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileHeatingTidal), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    class           (darkMatterProfileDMOClass    ), intent(inout) :: darkMatterProfileDMO_
    double precision                               , intent(in   ) :: radius
    double precision                                               :: energyPerturbationFirstOrder, energyPerturbationSecondOrder

    
    if (radius > 0.0d0) then
       call self%specificEnergyTerms(node,darkMatterProfileDMO_,radius,energyPerturbationFirstOrder,energyPerturbationSecondOrder)
       if (energyPerturbationSecondOrder > 0.0d0) then
          tidalSpecificEnergyGradient=+(                                                                                                &
               &                        +energyPerturbationFirstOrder *  2.0d0                                                          & !   dlog[r²    ]/dlog(r) term
               &                        +energyPerturbationSecondOrder*(                                                                &
               &                                                        -0.5d0                                                          & ! ⎧ dlog[σ_r(r)]/dlog[r] term
               &                                                        *darkMatterProfileDMO_%densityLogSlope         (node,radius)    & ! ⎥
               &                                                        -0.5d0                                                          & ! ⎥ Assumes the Jeans equation in
               &                                                        *gravitationalConstantGalacticus                                & ! ⎥ spherical symmetry with anisotropy
               &                                                        *darkMatterProfileDMO_%enclosedMass            (node,radius)    & ! ⎥ parameter β=0. Would be better to
               &                                                        /                                                    radius     & ! ⎥ have this provided by the
               &                                                        /darkMatterProfileDMO_%radialVelocityDispersion(node,radius)**2 & ! ⎩ darkMatterProfileDMO class.
               &                                                        +1.0d0                                                          & !   dlog[r     ]/dlog(r) term
               &                                                       )                                                                &
               &                       )                                                                                                &
               &                      /radius
       else
          tidalSpecificEnergyGradient=+  energyPerturbationFirstOrder *  2.0d0                                                          & !   dlog[r²    ]/dlog(r) term
               &                      /radius
       end if
    else
       tidalSpecificEnergyGradient=+0.0d0
    end if
    return
  end function tidalSpecificEnergyGradient

  subroutine tidalSpecificEnergyTerms(self,node,darkMatterProfileDMO_,radius,energyPerturbationFirstOrder,energyPerturbationSecondOrder)
    !!{
    Compute the first and second order perturbations to the energy.
    !!}
    implicit none
    class           (darkMatterProfileHeatingTidal), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    class           (darkMatterProfileDMOClass    ), intent(inout) :: darkMatterProfileDMO_
    double precision                               , intent(in   ) :: radius
    double precision                               , intent(  out) :: energyPerturbationFirstOrder, energyPerturbationSecondOrder

    energyPerturbationFirstOrder    =+self%specificEnergyOverRadiusSquared(node)                  &
         &                           *radius                                    **2
    if (self%coefficientSecondOrder > 0.0d0) then
       energyPerturbationSecondOrder=+sqrt(2.0d0)                                                 &
            &                        *self%coefficientSecondOrder                                 &
            &                        *(                                                           &
            &                          +1.0d0                                                     &
            &                          +self%correlationVelocityRadius                            &
            &                         )                                                           &
            &                        *sqrt(energyPerturbationFirstOrder)                          &
            &                        *darkMatterProfileDMO_%radialVelocityDispersion(node,radius)
    else
       energyPerturbationSecondOrder=+0.0d0
    end if
    return
  end subroutine tidalSpecificEnergyTerms
  
  logical function tidalSpecificEnergyIsEverywhereZero(self,node,darkMatterProfileDMO_)
    !!{
    Returns true if the specific energy is everywhere zero in the given {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileHeatingTidal), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node
    class(darkMatterProfileDMOClass    ), intent(inout) :: darkMatterProfileDMO_
    !$GLC attributes unused :: darkMatterProfileDMO_

    tidalSpecificEnergyIsEverywhereZero=self%specificEnergyOverRadiusSquared(node) <= 0.0d0
    return
  end function tidalSpecificEnergyIsEverywhereZero

  double precision function tidalSpecificEnergyOverRadiusSquared(self,node)
    !!{
    Compute $Q = E / r^2$.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite, treeNode
    implicit none
    class(darkMatterProfileHeatingTidal), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node
    class(nodeComponentSatellite       ), pointer       :: satellite

    if     (                                        &
         &   node%uniqueID() /= self%parentUniqueID &
         &  .and.                                   &
         &   node%uniqueID() /= self%lastUniqueID   &
         & ) call self%calculationReset(node)
    if (node%uniqueID() == self%parentUniqueID) then
       if (self%specificEnergyOverRadiusSquaredParent_ < 0.0d0) then
          satellite                                   =>      node     %satellite             ()
          self%specificEnergyOverRadiusSquaredParent_ =  max(                                     &
               &                                             +0.0d0                             , &
               &                                             +satellite%tidalHeatingNormalized()  &
               &                                            )
       end if
       tidalSpecificEnergyOverRadiusSquared=self%specificEnergyOverRadiusSquaredParent_
    else
       if (self%specificEnergyOverRadiusSquared_       < 0.0d0) then
          satellite                                   =>      node     %satellite             ()
          self%specificEnergyOverRadiusSquared_       =  max(                                     &
               &                                             +0.0d0                             , &
               &                                             +satellite%tidalHeatingNormalized()  &
               &                                            )
       end if
       tidalSpecificEnergyOverRadiusSquared=self%specificEnergyOverRadiusSquared_
    end if
    return
  end function tidalSpecificEnergyOverRadiusSquared
