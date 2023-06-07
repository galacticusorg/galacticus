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

  !!{
  Implements a tidal mass distribution heating class.
  !!}

  !![
  <massDistributionHeating name="massDistributionHeatingTidal">
    <description>A tidal mass distribution heating class.</description>
  </massDistributionHeating>
  !!]
  type, extends(massDistributionHeatingClass) :: massDistributionHeatingTidal
     !!{
     Implementation of a tidal mass distribution heating class.
     !!}
     private
     double precision :: correlationVelocityRadius, coefficientSecondOrder0, &
          &              coefficientSecondOrder1  , coefficientSecondOrder2, &
          &              heatSpecificNormalized
   contains
     !![
     <methods>
       <method description="Compute the first and second order energy perturbations." method="specificEnergyTerms"/>
     </methods>
     !!]
     procedure :: specificEnergy                 => tidalSpecificEnergy
     procedure :: specificEnergyGradient         => tidalSpecificEnergyGradient
     procedure :: specificEnergyIsEveryWhereZero => tidalSpecificEnergyIsEverywhereZero
     procedure :: specificEnergyTerms            => tidalSpecificEnergyTerms
  end type massDistributionHeatingTidal

  interface massDistributionHeatingTidal
     !!{
     Constructors for the {\normalfont \ttfamily tidal} mass distribution class.
     !!}
     module procedure tidalConstructorParameters
     module procedure tidalConstructorInternal
  end interface massDistributionHeatingTidal

contains

  function tidalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily tidal} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (massDistributionHeatingTidal)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    double precision                                              :: heatSpecificNormalized , correlationVelocityRadius, &
         &                                                           coefficientSecondOrder0, coefficientSecondOrder1  , &
         &                                                           coefficientSecondOrder2

    !![
    <inputParameter>
      <name>heatSpecificNormalized</name>
      <description>The normalized specific tidal heating, $Q = \epsilon / r^2$.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>coefficientSecondOrder0</name>
      <defaultValue>0.0d0</defaultValue>
      <source>parameters</source>
      <description>The coefficient, $a_0$, appearing in the second-order heating term, $f_2 = a_0 + a_1 \mathrm{d}\log \rho/\mathrm{d} \log r + a_2 (\mathrm{d}\log \rho/\mathrm{d} \log r)^2$.</description>
    </inputParameter>
    <inputParameter>
      <name>coefficientSecondOrder1</name>
      <defaultValue>0.0d0</defaultValue>
      <source>parameters</source>
      <description>The coefficient, $a_1$, appearing in the second-order heating term, $f_2 = a_0 + a_1 \mathrm{d}\log \rho/\mathrm{d} \log r + a_2 (\mathrm{d}\log \rho/\mathrm{d} \log r)^2$.</description>
    </inputParameter>
    <inputParameter>
      <name>coefficientSecondOrder2</name>
      <defaultValue>0.0d0</defaultValue>
      <source>parameters</source>
      <description>The coefficient, $a_2$, appearing in the second-order heating term, $f_2 = a_0 + a_1 \mathrm{d}\log \rho/\mathrm{d} \log r + a_2 (\mathrm{d}\log \rho/\mathrm{d} \log r)^2$.</description>
    </inputParameter>
    <inputParameter>
      <name>correlationVelocityRadius</name>
      <defaultValue>-1.0d0</defaultValue>
      <source>parameters</source>
      <description>The velocity-position correlation function, $\chi_\mathrm{r,v}$, as defined by \cite[][eqn.~B1]{gnedin_self-consistent_1999} which controls the strength of the second order heating term.</description>
    </inputParameter>
    !!]
    self=massDistributionHeatingTidal(heatSpecificNormalized,coefficientSecondOrder0,coefficientSecondOrder1,coefficientSecondOrder2,correlationVelocityRadius)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function tidalConstructorParameters
  
  function tidalConstructorInternal(heatSpecificNormalized,coefficientSecondOrder0,coefficientSecondOrder1,coefficientSecondOrder2,correlationVelocityRadius) result(self)
    !!{
    Constructor for ``tidal'' convergence class.
    !!}
    implicit none
    type             (massDistributionHeatingTidal)                :: self
    double precision                               , intent(in   ) :: heatSpecificNormalized   , coefficientSecondOrder0, &
         &                                                            coefficientSecondOrder1  , coefficientSecondOrder2, &
         &                                                            correlationVelocityRadius
    !![
    <constructorAssign variables="heatSpecificNormalized, coefficientSecondOrder0, coefficientSecondOrder1, coefficientSecondOrder2, correlationVelocityRadius"/>
    !!]
 
    return
  end function tidalConstructorInternal

  double precision function tidalSpecificEnergy(self,radius,massDistribution_) result(energySpecific)
    !!{
    Compute the specific energy in a tidally-heated mass distribution.
    !!}
    implicit none
    class           (massDistributionHeatingTidal), intent(inout) :: self
    double precision                              , intent(in   ) :: radius
    class           (massDistributionClass       ), intent(inout) :: massDistribution_
    double precision                                              :: energyPerturbationFirstOrder, energyPerturbationSecondOrder

    call self%specificEnergyTerms(radius,massDistribution_,energyPerturbationFirstOrder,energyPerturbationSecondOrder)
    energySpecific=+energyPerturbationFirstOrder  &
         &         +energyPerturbationSecondOrder
    return
  end function tidalSpecificEnergy

  double precision function tidalSpecificEnergyGradient(self,radius,massDistribution_) result(energySpecificGradient)
    !!{
    Returns the gradient of the specific energy of heating.
    !!}
    use :: Coordinates                     , only : coordinateSpherical            , assignment(=)
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (massDistributionHeatingTidal), intent(inout) :: self
    double precision                              , intent(in   ) :: radius
    class           (massDistributionClass       ), intent(inout) :: massDistribution_
    double precision                                              :: energyPerturbationFirstOrder, energyPerturbationSecondOrder
    type            (coordinateSpherical         )                :: coordinates

    if (radius > 0.0d0) then
       call self%specificEnergyTerms(radius,massDistribution_,energyPerturbationFirstOrder,energyPerturbationSecondOrder)
       if (energyPerturbationSecondOrder > 0.0d0) then
       coordinates=[radius,0.0d0,0.0d0]
          energySpecificGradient=+(                                                                                                                                              &
               &                   +energyPerturbationFirstOrder *  2.0d0                                                                                                        & !   dlog[r²    ]/dlog(r) term
               &                   +energyPerturbationSecondOrder*(                                                                                                              &
               &                                                   -0.5d0                                                                                                        & ! ⎧ dlog[σ_r(r)]/dlog[r] term
               &                                                   *massDistribution_%densityGradientRadial                       (coordinates,logarithmic=.true.           )    & ! ⎥
               &                                                   -0.5d0                                                                                                        & ! ⎥ Assumes the Jeans equation in
               &                                                   *gravitationalConstantGalacticus                                                                              & ! ⎥ spherical symmetry with anisotropy
               &                                                   *massDistribution_%massEnclosedBySphere                        (radius                                   )    & ! ⎥ parameter β=0. Would be better to
               &                                                   /                                                               radius                                        & ! ⎥ have this provided by the
               &                                                   /massDistribution_%kinematicsDistribution_%velocityDispersion1D(coordinates,            massDistribution_)**2 & ! ⎩ darkMatterProfileDMO class.
               &                                                   +1.0d0                                                                                                        & !   dlog[r     ]/dlog(r) term
               &                                                  )                                                                                                              &
               &                  )                                                                                                                                              &
               &                 /radius
       else
          energySpecificGradient=+  energyPerturbationFirstOrder *  2.0d0                                                                                                        & !   dlog[r²    ]/dlog(r) term
               &                 /radius
       end if
    else
       energySpecificGradient   =+0.0d0
    end if
    return
  end function tidalSpecificEnergyGradient

  subroutine tidalSpecificEnergyTerms(self,radius,massDistribution_,energyPerturbationFirstOrder,energyPerturbationSecondOrder)
    !!{
    Compute the first and second order perturbations to the energy.
    !!}
    use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    class           (massDistributionHeatingTidal), intent(inout) :: self
    double precision                              , intent(in   ) :: radius
    class           (massDistributionClass       ), intent(inout) :: massDistribution_
    double precision                              , intent(  out) :: energyPerturbationFirstOrder, energyPerturbationSecondOrder
    double precision                                              :: coefficientSecondOrder      , densityLogSlope
    type            (coordinateSpherical         )                :: coordinates

    energyPerturbationFirstOrder=+self%heatSpecificNormalized    &
         &                       *radius                     **2
    if     (                                       &
         &   self%coefficientSecondOrder0 /= 0.0d0 &
         &  .or.                                   &
         &   self%coefficientSecondOrder1 /= 0.0d0 &
         &  .or.                                   &
         &   self%coefficientSecondOrder2 /= 0.0d0 &
         & ) then
       ! Compute the coefficient for the second order term.
       coordinates=[radius,0.0d0,0.0d0]
       densityLogSlope       =+massDistribution_%densityGradientRadial(coordinates,logarithmic=.true.)
       coefficientSecondOrder=+self             %coefficientSecondOrder0                    &
            &                 +self             %coefficientSecondOrder1*densityLogSlope    &
            &                 +self             %coefficientSecondOrder2*densityLogSlope**2
       ! Compute the second order energy perturbation.
       energyPerturbationSecondOrder=+sqrt(2.0d0)                                                                                   &
            &                        *coefficientSecondOrder                                                                        &
            &                        *(                                                                                             &
            &                          +1.0d0                                                                                       &
            &                          +self%correlationVelocityRadius                                                              &
            &                         )                                                                                             &
            &                        *sqrt(energyPerturbationFirstOrder)                                                            &
            &                        *massDistribution_%kinematicsDistribution_%velocityDispersion1D(coordinates,massDistribution_)
    else
       energyPerturbationSecondOrder=+0.0d0
    end if
    return
  end subroutine tidalSpecificEnergyTerms

  logical function tidalSpecificEnergyIsEverywhereZero(self) result(energySpecificIsEverywhereZero)
    !!{
    Returns true if the specific energy is everywhere zero.
    !!}
    implicit none
    class(massDistributionHeatingTidal), intent(inout) :: self

    energySpecificIsEverywhereZero=self%heatSpecificNormalized <= 0.0d0
    return
  end function tidalSpecificEnergyIsEverywhereZero
