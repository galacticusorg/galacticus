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
  Implements a supernovae type Ia class based on \cite{nagashima_metal_2005}.
  !!}

  use :: Stellar_Astrophysics , only : stellarAstrophysicsClass
  use :: Numerical_Integration, only : integrator
  
  !![
  <supernovaeTypeIa name="supernovaeTypeIaNagashima2005">
   <description>
    A supernovae type Ia class which uses the prescriptions from \cite{nagashima_metal_2005} to compute the numbers and yields
    of Type Ia supernovae.
   </description>
  </supernovaeTypeIa>
  !!]
  type, extends(supernovaeTypeIaFixedYield) :: supernovaeTypeIaNagashima2005
     !!{
     A supernovae type Ia class based on \cite{nagashima_metal_2005}.
     !!}
     private
     class(stellarAstrophysicsClass), pointer     :: stellarAstrophysics_ => null()
     type (integrator              ), allocatable :: integrator_
   contains    
     final     ::                     nagashima2005Destructor
     procedure :: massInitialRange => nagashima2005MassInitialRange
     procedure :: number           => nagashima2005Number
  end type supernovaeTypeIaNagashima2005

  interface supernovaeTypeIaNagashima2005
     !!{
     Constructors for the \refClass{supernovaeTypeIaNagashima2005} supernovae type Ia class.
     !!}
     module procedure nagashima2005ConstructorParameters
     module procedure nagashima2005ConstructorInternal
  end interface supernovaeTypeIaNagashima2005

  ! Sub-module-scope variables used in integration.
  class           (initialMassFunctionClass), pointer :: initialMassFunction__
  double precision                                    :: massSecondary
  !$omp threadprivate(initialMassFunction__,massSecondary)

  ! Parameters of the distribution of binaries from Nagashima et al. (2005; MNRAS; 358; 1427; eqn. 17).
  double precision, parameter :: primaryMassMaximum=6.0d0
  double precision, parameter :: binaryMassMaximum =2.0d0*primaryMassMaximum, binaryMassMinimum  =3.00d0
  double precision, parameter :: gamma             =2.0d0                   , typeIaNormalization=0.07d0
  
contains

  function nagashima2005ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{supernovaeTypeIaNagashima2005} supernovae type Ia class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (supernovaeTypeIaNagashima2005)                :: self
    type (inputParameters              ), intent(inout) :: parameters
    class(stellarAstrophysicsClass     ), pointer       :: stellarAstrophysics_

    !![
    <objectBuilder class="stellarAstrophysics" name="stellarAstrophysics_" source="parameters"/>
    !!]
    self=supernovaeTypeIaNagashima2005(stellarAstrophysics_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="stellarAstrophysics_"/>
    !!]
    return
  end function nagashima2005ConstructorParameters

  function nagashima2005ConstructorInternal(stellarAstrophysics_) result(self)
    !!{
    Internal constructor for the \refClass{supernovaeTypeIaNagashima2005} supernovae type Ia class.
    !!}
    implicit none
    type (supernovaeTypeIaNagashima2005)                        :: self
    class(stellarAstrophysicsClass     ), intent(in   ), target :: stellarAstrophysics_
    !![
    <constructorAssign variables="*stellarAstrophysics_"/>
    !!]

    self%initialized=.false.
    allocate(self%integrator_)
    self%integrator_=integrator(nagashima2005NumberIntegrand,toleranceRelative=1.0d-3)
    return
  end function nagashima2005ConstructorInternal

  subroutine nagashima2005Destructor(self)
    !!{
    Destructor for the \refClass{supernovaeTypeIaNagashima2005} supernovae type Ia class.
    !!}
    implicit none
    type(supernovaeTypeIaNagashima2005), intent(inout) :: self

    !![
    <objectDestructor name="self%stellarAstrophysics_"/>
    !!]
    return
  end subroutine nagashima2005Destructor

  subroutine nagashima2005MassInitialRange(self,initialMassFunction_,age,metallicity,massInitialMinimum,massInitialMaximum)
    !!{
    Return the range of initial stellar masses contributing to the Type Ia population.
    !!}
    implicit none
     class           (supernovaeTypeIaNagashima2005), intent(inout) :: self
     class           (initialMassFunctionClass     ), intent(inout) :: initialMassFunction_
     double precision                               , intent(in   ) :: age                 , metallicity
     double precision                               , intent(  out) :: massInitialMinimum  , massInitialMaximum
     !$GLC attributes unused :: initialMassFunction_

     ! The minimum initial mass is set by the requirement that the secondary has evolved off of the main sequence at this age. The
     ! maximum mass is the largest single-star mass for which the endpoint is a C-O white dwarf (Nagashima et al. 2005).
     massInitialMinimum=self%stellarAstrophysics_%massInitial       (age,metallicity)
     massInitialMaximum=                          primaryMassMaximum
    return
  end subroutine nagashima2005MassInitialRange
  
  double precision function nagashima2005Number(self,initialMassFunction_,initialMass,age,metallicity)
    !!{
    Compute the cumulative number of Type Ia supernovae originating per unit interval of secondary star mass with given
    {\normalfont \ttfamily initialMass} and {\normalfont \ttfamily metallicity} after a time {\normalfont \ttfamily age}. The
    calculation is based on that of \cite{nagashima_metal_2005}. This function is expected to be integrated over the initial mass
    function of secondary stars.
    !!}
    implicit none
    class           (supernovaeTypeIaNagashima2005), intent(inout), target :: self
    class           (initialMassFunctionClass     ), intent(inout), target :: initialMassFunction_
    double precision                               , intent(in   )         :: age                         , initialMass, &
         &                                                                    metallicity
    double precision                                                       :: massFractionSecondaryMinimum
    
    ! Check if the secondary has evolved off of the main sequence at this age.
    if (initialMass > self%stellarAstrophysics_%massInitial(age,metallicity)) then
       ! The secondary has evolved off of the main sequence. Integrate over all possible secondary mass fractions, μ=m₂/m_b. The
       ! minimum secondary mass fraction is determined by the maximum possible primary mass (i.e. the maximum mass in the initial
       ! mass function).
       initialMassFunction__        => initialMassFunction_
       massSecondary                =  initialMass
       massFractionSecondaryMinimum =  +1.0d0                                  &
            &                          /(                                      &
            &                            +1.0d0                                &
            &                            +initialMassFunction_%massMaximum  () &
            &                            /                     massSecondary   &
            &                           )
       nagashima2005Number          =  self%integrator_%integrate(massFractionSecondaryMinimum,0.5d0)       
    else
       ! Secondary has not yet evolved off of the main sequence - no SNIa occurs as yet.
       nagashima2005Number=0.0d0
    end if
    return
  end function nagashima2005Number
  
  double precision function nagashima2005NumberIntegrand(massFractionSecondary) result(integrand)
    !!{
    Integrand used in computing the number of Type Ia supernovae.
    !!}
    implicit none
    double precision, intent(in   ) :: massFractionSecondary
    double precision                :: massBinary           , massPrimary
    
    ! Check if the binary initial mass is within the range of binary masses that lead to Type Ia supernovae, and that the primary
    ! mass is below the maximum single-star mass to produce a C-O white dwarf.
    massBinary=+massSecondary         &
         &     /massFractionSecondary
    massPrimary=+massBinary           &
         &      -massSecondary
    if     (                                  &
         &   massBinary  > binaryMassMinimum  &
         &  .and.                             &
         &   massBinary  < binaryMassMaximum  &
         &  .and.                             &
         &   massPrimary < primaryMassMaximum &
         & ) then
       ! Evaluate the integrand. Nagashima et al. (2005) give this integrand in the 2D space of (m_b,μ). Here, the quantity we
       ! compute will be integrated over m₂, weighted by the initial mass function of for single-stars with masses corresponding
       ! to our secondaries. As such, we must:
       !   1. divide the integrand given by Nagashima et al. (2005) by φ(m₂) (as it will later be multiplied by it), and;
       !   2. multiply the integrand by dm_b/dm₂ = 1/υ to convert from integration over m_b to integration over m₂.
       integrand=+                      typeIaNormalization                             &
            &    *initialMassFunction__%phi                     (massBinary           ) &
            &    /initialMassFunction__%phi                     (massSecondary        ) &
            &    *                      massFractionDistribution(massFractionSecondary) &
            &    /                      massFractionSecondary
    else
       integrand=0.0d0
    end if
    return
    
  contains

    double precision function massFractionDistribution(massFractionSecondary)
      !!{
      The distribution function for secondary mass fractions in binary star systems from \cite{nagashima_metal_2005}.
      !!}
      implicit none
      double precision, intent(in   ) :: massFractionSecondary

      massFractionDistribution=+2.0d0**(1.0d0+gamma)         &
           &                   *       (1.0d0+gamma)         &
           &                   *massFractionSecondary**gamma
      return
    end function massFractionDistribution

  end function nagashima2005NumberIntegrand
