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
  Implements a stellar winds class based on \cite{leitherer_deposition_1992}.
  !!}

  use :: Numerical_Constants_Astronomical, only : metallicitySolar
  use :: Stellar_Astrophysics_Tracks     , only : stellarTracksClass

  !![
  <stellarWinds name="stellarWindsLeitherer1992">
   <description>
    A stellar winds class using the fitting formulae of \cite{leitherer_deposition_1992} to compute stellar wind energy input
    from the luminosity and effective temperature of a star.
   </description>
  </stellarWinds>
  !!]
  type, extends(stellarWindsClass) :: stellarWindsLeitherer1992
     !!{
     A stellar winds class based on \cite{leitherer_deposition_1992}.
     !!}
     private
     class(stellarTracksClass), pointer :: stellarTracks_ => null()
   contains
     final     ::                     leitherer1992Destructor
     procedure :: rateMassLoss     => leitherer1992RateMassLoss
     procedure :: velocityTerminal => leitherer1992VelocityTerminal
  end type stellarWindsLeitherer1992

  interface stellarWindsLeitherer1992
     !!{
     Constructors for the \refClass{stellarWindsLeitherer1992} stellar winds class.
     !!}
     module procedure leitherer1992ConstructorParameters
     module procedure leitherer1992ConstructorInternal
  end interface stellarWindsLeitherer1992

  ! Minimum metallicity to which we trust Leitherer et al.'s metallicity scaling.
  double precision, parameter :: metallicityMinimum=1.0d-4*metallicitySolar

contains

  function leitherer1992ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{stellarWindsLeitherer1992} stellar winds class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (stellarWindsLeitherer1992)                :: self
    type (inputParameters          ), intent(inout) :: parameters
    class(stellarTracksClass       ), pointer       :: stellarTracks_

    !![
    <objectBuilder class="stellarTracks" name="stellarTracks_" source="parameters"/>
    !!]
    self=stellarWindsLeitherer1992(stellarTracks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="stellarTracks_"/>
    !!]
    return
  end function leitherer1992ConstructorParameters

  function leitherer1992ConstructorInternal(stellarTracks_) result(self)
    !!{
    Internal constructor for the \refClass{stellarWindsLeitherer1992} stellar winds class.
    !!}
    implicit none
    type (stellarWindsLeitherer1992)                        :: self
    class(stellarTracksClass       ), intent(in   ), target :: stellarTracks_
    !![
    <constructorAssign variables="*stellarTracks_"/>
    !!]

    return
  end function leitherer1992ConstructorInternal

  subroutine leitherer1992Destructor(self)
    !!{
    Destructor for the \refClass{stellarWindsLeitherer1992} stellar winds class.
    !!}
    implicit none
    type(stellarWindsLeitherer1992), intent(inout) :: self

    !![
    <objectDestructor name="self%stellarTracks_"/>
    !!]
    return
  end subroutine leitherer1992Destructor

  double precision function leitherer1992RateMassLoss(self,initialMass,age,metallicity)
    !!{
    Compute the mass loss rate (in $M_\odot$/Gyr) from a star of given {\normalfont \ttfamily initialMass}, {\normalfont
    \ttfamily age} and {\normalfont \ttfamily metallicity} using the fitting formula of \cite{leitherer_deposition_1992}.
    !!}
    implicit none
    class           (stellarWindsLeitherer1992), intent(inout) :: self
    double precision                           , intent(in   ) :: age                        , initialMass       , &
         &                                                        metallicity
    double precision                                           :: stellarEffectiveTemperature, stellarLuminosity
    !$GLC attributes unused :: self

    ! Get luminosity and effective temperature of the star.
    stellarLuminosity          =self%stellarTracks_%luminosity          (initialMass,metallicity,age)
    stellarEffectiveTemperature=self%stellarTracks_%temperatureEffective(initialMass,metallicity,age)
    ! Compute mass loss rate using fitting formula. (Initial constant 9 converts Leitherer's mass loss rate from per year to per Gyr.)
    if     (                                     &
         &   stellarLuminosity           > 0.0d0 &
         &  .and.                                &
         &   stellarEffectiveTemperature > 0.0d0 &
         & ) then
       leitherer1992RateMassLoss=+10.0d0**(                                                                                                   &
            &                              + 9.00d0                                                                                           &
            &                              -24.06d0                                                                                           &
            &                              + 2.45d0*log10(    stellarLuminosity                                                             ) &
            &                              - 1.10d0*log10(    initialMass                                                                   ) &
            &                              + 1.31d0*log10(    stellarEffectiveTemperature                                                   ) &
            &                              + 0.80d0*log10(max(metallicity                 ,metallicityMinimum)/metallicitySolar) &
            &                             )
    else
       leitherer1992RateMassLoss=+ 0.0d0
    end if
    return
  end function leitherer1992RateMassLoss

  double precision function leitherer1992VelocityTerminal(self,initialMass,age,metallicity)
    !!{
    Compute the terminal velocity (in km/s) from a star of given {\normalfont \ttfamily initialMass}, {\normalfont \ttfamily age} and {\normalfont \ttfamily metallicity} using
    the fitting formula of \cite{leitherer_deposition_1992}.
    !!}
    implicit none
     class           (stellarWindsLeitherer1992), intent(inout) :: self
     double precision                           , intent(in   ) :: age                        , initialMass       , &
         &                                                         metallicity
    double precision                                            :: stellarEffectiveTemperature, stellarLuminosity
    !$GLC attributes unused :: self

    ! Get luminosity and effective temperature of the star.
    stellarLuminosity          =self%stellarTracks_%luminosity          (initialMass,metallicity,age)
    stellarEffectiveTemperature=self%stellarTracks_%temperatureEffective(initialMass,metallicity,age)
    ! Compute mass loss rate using fitting formula.
    if (stellarLuminosity > 0.0d0 .and. stellarEffectiveTemperature > 0.0d0) then
       leitherer1992VelocityTerminal=+10.0d0**(                                                                                    &
            &                                  +1.23d0                                                                             &
            &                                  -0.30d0*log10(    stellarLuminosity                                               ) &
            &                                  +0.55d0*log10(    initialMass                                                     ) &
            &                                  +0.64d0*log10(    stellarEffectiveTemperature                                     ) &
            &                                  +0.12d0*log10(max(metallicity                ,metallicityMinimum)/metallicitySolar) &
            &                                 )
    else
       leitherer1992VelocityTerminal=+ 0.0d0
    end if
    return
  end function leitherer1992VelocityTerminal
