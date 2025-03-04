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
  An implementation of atomic recombination recombination rates which computes rates from the corresponding photoionization
  cross section using detailed balance (i.e. the Milne relation---see for example ``Astrophysics of the Diffuse Universe'' by
  Dopita \& Sutherland, 2004, Springer Science \& Business Media, section 5.3.3).
  !!}

  use :: Atomic_Cross_Sections_Ionization_Photo , only : atomicCrossSectionIonizationPhotoClass
  use :: Atomic_Ionization_Potentials           , only : atomicIonizationPotentialClass

  !![
  <atomicRecombinationRateRadiative name="atomicRecombinationRateRadiativeComputed">
   <description>Atomic radiative recombination rates computed from the corresponding photoionization cross section using detailed balance (i.e. the Milne relation---see for example ``Astrophysics of the Diffuse Universe'' by Dopita \&amp; Sutherland, 2004, Springer Science \&amp; Business Media, section 5.3.3).</description>
  </atomicRecombinationRateRadiative>
  !!]
  type, extends(atomicRecombinationRateRadiativeClass) :: atomicRecombinationRateRadiativeComputed
     !!{
     A recombination recombination rate class assuming a thermal electron distribution.
     !!}
     private
     class(atomicCrossSectionIonizationPhotoClass), pointer :: atomicCrossSectionIonizationPhoto_ => null()
     class(atomicIonizationPotentialClass        ), pointer :: atomicIonizationPotential_         => null()
   contains
     final     ::         computedDestructor
     procedure :: rate => computedRate
  end type atomicRecombinationRateRadiativeComputed

  interface atomicRecombinationRateRadiativeComputed
     !!{
     Constructors for the {\normalfont \ttfamily computed} atomic radiative recombination class.
     !!}
     module procedure computedConstructorParameters
     module procedure computedConstructorInternal
  end interface atomicRecombinationRateRadiativeComputed

contains

  function computedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily computed} atomic radiative recombination class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (atomicRecombinationRateRadiativeComputed)                :: self
    type (inputParameters                         ), intent(inout) :: parameters
    class(atomicCrossSectionIonizationPhotoClass  ), pointer       :: atomicCrossSectionIonizationPhoto_
    class(atomicIonizationPotentialClass          ), pointer       :: atomicIonizationPotential_


    !![
    <objectBuilder class="atomicCrossSectionIonizationPhoto" name="atomicCrossSectionIonizationPhoto_" source="parameters"/>
    <objectBuilder class="atomicIonizationPotential"         name="atomicIonizationPotential_"         source="parameters"/>
    !!]
    self=atomicRecombinationRateRadiativeComputed(atomicCrossSectionIonizationPhoto_,atomicIonizationPotential_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="atomicCrossSectionIonizationPhoto_"/>
    <objectDestructor name="atomicIonizationPotential_"        />
    !!]
    return
  end function computedConstructorParameters

  function computedConstructorInternal(atomicCrossSectionIonizationPhoto_,atomicIonizationPotential_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily computed} atomic radiative recombination class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (atomicRecombinationRateRadiativeComputed)                        :: self
    class(atomicCrossSectionIonizationPhotoClass  ), intent(in   ), target :: atomicCrossSectionIonizationPhoto_
    class(atomicIonizationPotentialClass          ), intent(in   ), target :: atomicIonizationPotential_
    !![
    <constructorAssign variables="*atomicCrossSectionIonizationPhoto_, *atomicIonizationPotential_"/>
    !!]

    return
  end function computedConstructorInternal

  subroutine computedDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily computed} recombination recombination class.
    !!}
    implicit none
    type(atomicRecombinationRateRadiativeComputed), intent(inout) :: self

    !![
    <objectDestructor name="self%atomicCrossSectionIonizationPhoto_"/>
    <objectDestructor name="self%atomicIonizationPotential_"        />
    !!]
    return
  end subroutine computedDestructor

  double precision function computedRate(self,atomicNumber,ionizationState,temperature,level)
    !!{
    Returns the recombination rate coefficient.
    !!}
    use :: Error                       , only : Error_Report
    use :: Numerical_Integration       , only : integrator
    use :: Numerical_Constants_Physical, only : boltzmannsConstant, electronMass
    use :: Numerical_Constants_Atomic  , only : atomicMassUnit
    implicit none
    class           (atomicRecombinationRateRadiativeComputed), intent(inout)           :: self
    integer                                                   , intent(in   )           :: atomicNumber                 , ionizationState
    double precision                                          , intent(in   )           :: temperature
    type            (enumerationRecombinationCaseType        ), intent(in   ), optional :: level
    double precision                                                                    :: velocityMaximumFactor=100.0d0
    type            (integrator                              )                          :: integrator_
    double precision                                                                    :: velocityMaximum

    ! A level is required.
    if (present(level)) then
       ! Find a maximum velocity to integrate to which is some large multiple of the typical thermal velocity.
       velocityMaximum=+velocityMaximumFactor    &
            &          *sqrt(                    &
            &                +boltzmannsConstant &
            &                *temperature        &
            &                /electronMass       &
            &               )
       ! Compute the recombination rate by integration of the recombination cross section over the thermal distribution of electron
       ! velocities.
       integrator_ =integrator           (rateIntegrand,toleranceRelative=1.0d-4         )
       computedRate=integrator_%integrate(0.0d0        ,                  velocityMaximum)
    else
       computedRate=0.0d0
       call Error_Report('calculation is only supported to specific levels'//{introspection:location})
    end if
    return

  contains

    double precision function rateIntegrand(velocity)
      !!{
      Integrand for the recombination coefficient due to recombination, see equation (2.5) of
      \cite{osterbrock_astrophysics_2006}, but note that their expression for the Maxwell-Boltzmann distribution function in
      equation (2.6) is missing a factor (1/2) in the exponential term.
      !!}
      use :: Numerical_Constants_Math    , only : Pi
      use :: Numerical_Constants_Units   , only : metersToAngstroms, electronVolt
      use :: Numerical_Constants_Physical, only : plancksConstant  , speedLight  , electronMass, fineStructure
      use :: Numerical_Constants_Prefixes, only : centi
      implicit none
      double precision, intent(in   ) :: velocity
      double precision                :: energyPhoton    , wavelengthPhoton         , &
           &                             energyIonization, crossSectionRecombination
      
      energyIonization=self%atomicIonizationPotential_%potential(atomicNumber,ionizationState)*electronVolt
      energyPhoton=0.5d0*electronMass*velocity**2+energyIonization
      wavelengthPhoton=metersToAngstroms*speedLight*plancksConstant/energyPhoton
      if (velocity > 0.0d0) then
         ! Recombination cross section is computed from the corresponding photoionization cross section using the result from
         ! detailed balance (Arfken, 1961, Ionization of the Interplanetary Gas, Rept. LAMS-2596, Los Alamos Scientific Laboratory
         ! of the University of California). The cross section used here is summed over all states corresponding to the given
         ! principal quantum number (the argument "level").
         crossSectionRecombination=+(+energyPhoton/energyIonization      )**2                                                                                          &
              &                    /(+energyPhoton/energyIonization-1.0d0)                                                                                             &
              &                    *fineStructure                         **2                                                                                          &
              &                    /2.0d0                                                                                                                              &
              &                    *self%atomicCrossSectionIonizationPhoto_%crossSection(atomicNumber,ionizationState,shellNumber=level%ID,wavelength=wavelengthPhoton)
         rateIntegrand=+4.0d0                     &
              &        /sqrt(Pi)                  &
              &        *(                         &
              &          +electronMass            &
              &          /2.0d0                   &
              &          /boltzmannsConstant      &
              &          /temperature             &
              &        )**1.5d0                   &
              &        *velocity**3               &
              &        *exp(                      &
              &             -0.5d0                &
              &             *electronMass         &
              &             *velocity**2          &
              &             /boltzmannsConstant   &
              &             /temperature          &
              &            )                      &
              &        *crossSectionRecombination &
              &        /centi
      else
         rateIntegrand=0.0d0
      end if
      return
    end function rateIntegrand

  end function computedRate
