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
Implements the geometry of the SDSS survey used by \cite{montero-dorta_sdss_2009}.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <surveyGeometry name="surveyGeometryMonteroDorta2009SDSS">
   <description>Implements the geometry of the SDSS survey of \cite{montero-dorta_sdss_2009}.</description>
  </surveyGeometry>
  !!]
  type, extends(surveyGeometryBernardi2013SDSS) :: surveyGeometryMonteroDorta2009SDSS
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     character       (len=1                  )          :: band
     double precision                                   :: redshiftMinimum         , redshiftMaximum         , &
          &                                                magnitudeApparentMinimum, magnitudeApparentMaximum
   contains
     final     ::                    monteroDorta2009SDSSDestructor
     procedure :: distanceMinimum => monteroDorta2009SDSSDistanceMinimum
     procedure :: distanceMaximum => monteroDorta2009SDSSDistanceMaximum
  end type surveyGeometryMonteroDorta2009SDSS

  interface surveyGeometryMonteroDorta2009SDSS
     !!{
     Constructors for the \cite{montero-dorta_sdss_2009} survey geometry class.
     !!}
     module procedure monteroDorta2009SDSSConstructorParameters
     module procedure monteroDorta2009SDSSConstructorInternal
  end interface surveyGeometryMonteroDorta2009SDSS

contains

  function monteroDorta2009SDSSConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \cite{montero-dorta_sdss_2009} conditional mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type     (surveyGeometryMonteroDorta2009SDSS)                :: self
    type     (inputParameters                   ), intent(inout) :: parameters
    class    (cosmologyFunctionsClass           ), pointer       :: cosmologyFunctions_
    character(len=1                             )                :: band

    !![
    <inputParameter>
      <name>band</name>
      <source>parameters</source>
      <description>The band for which the survey geometry should be computed.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=surveyGeometryMonteroDorta2009SDSS(band,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function monteroDorta2009SDSSConstructorParameters

  function monteroDorta2009SDSSConstructorInternal(band,cosmologyFunctions_,redshiftMinimum,redshiftMaximum) result (self)
    !!{
    Default constructor for the \cite{montero-dorta_sdss_2009} survey geometry class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (surveyGeometryMonteroDorta2009SDSS)                          :: self
    character       (len=1                             ), intent(in   )           :: band
    class           (cosmologyFunctionsClass           ), intent(in   ), target   :: cosmologyFunctions_
    double precision                                    , intent(in   ), optional :: redshiftMinimum    , redshiftMaximum
    !![
    <constructorAssign variables="band, *cosmologyFunctions_"/>
    !!]

    ! Set survey limits.
    select case (band)
    case ('u')
       self%redshiftMinimum         = 0.02d0
       self%redshiftMaximum         = 0.19d0
       self%magnitudeApparentMinimum=16.45d0
       self%magnitudeApparentMaximum=19.00d0
    case ('g')
       self%redshiftMinimum         = 0.02d0
       self%redshiftMaximum         = 0.16d0
       self%magnitudeApparentMinimum=14.55d0
       self%magnitudeApparentMaximum=17.91d0
    case ('r')
       self%redshiftMinimum         = 0.02d0
       self%redshiftMaximum         = 0.22d0
       self%magnitudeApparentMinimum=13.93d0
       self%magnitudeApparentMaximum=17.77d0
    case ('i')
       self%redshiftMinimum         = 0.02d0
       self%redshiftMaximum         = 0.22d0
       self%magnitudeApparentMinimum=13.55d0
       self%magnitudeApparentMaximum=17.24d0
    case ('z')
       self%redshiftMinimum         = 0.02d0
       self%redshiftMaximum         = 0.23d0
       self%magnitudeApparentMinimum=13.40d0
       self%magnitudeApparentMaximum=16.97d0
    case default
       call Error_Report('band ∈ {u,g,r,i,z} is required'//{introspection:location})
    end select
    ! If redshift ranges are provided, override the defaults.
    if (present(redshiftMinimum)) self%redshiftMinimum=redshiftMinimum
    if (present(redshiftMaximum)) self%redshiftMaximum=redshiftMaximum
    ! Initialize geometry.
    call self%initialize()
    return
  end function monteroDorta2009SDSSConstructorInternal

  subroutine monteroDorta2009SDSSDestructor(self)
    !!{
    Destructor for the \refClass{surveyGeometryMonteroDorta2009SDSS} survey geometry class.
    !!}
    implicit none
    type(surveyGeometryMonteroDorta2009SDSS), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine monteroDorta2009SDSSDestructor

  double precision function monteroDorta2009SDSSDistanceMinimum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (surveyGeometryMonteroDorta2009SDSS), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                                           luminosity, starFormationRate
    integer                                             , intent(in   ), optional :: field
    !$GLC attributes unused :: field, mass, luminosity, starFormationRate

    ! Validate input.
    if (.not.present(magnitudeAbsolute)) call Error_Report('absolute magnitude must be supplied '//{introspection:location})
    ! Compute limiting distances. This is due only to the redshift limit.
    monteroDorta2009SDSSDistanceMinimum=self   %cosmologyFunctions_%distanceComoving           (                        &
         &                               self  %cosmologyFunctions_%cosmicTime                  (                       &
         &                                self %cosmologyFunctions_%expansionFactorFromRedshift  (                      &
         &                                                                                        +self%redshiftMinimum &
         &                                                                                       )                      &
         &                                                                                      )                       &
         &                                                                                     )
    return
  end function monteroDorta2009SDSSDistanceMinimum

  double precision function monteroDorta2009SDSSDistanceMaximum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    use :: Error                      , only : Error_Report
    implicit none
    class           (surveyGeometryMonteroDorta2009SDSS), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: mass                   , magnitudeAbsolute       , &
         &                                                                           luminosity             , starFormationRate
    integer                                             , intent(in   ), optional :: field
    double precision                                                              :: distanceMaximumRedshift, distanceMaximumMagnitude
    !$GLC attributes unused :: field

    ! Validate arguments.
    if (present(mass             )) call Error_Report(             '`mass` is not supported'//{introspection:location})
    if (present(luminosity       )) call Error_Report(       '`luminosity` is not supported'//{introspection:location})
    if (present(starFormationRate)) call Error_Report('`starFormationRate` is not supported'//{introspection:location})
    ! Compute limiting distances. Note that for the magnitudes we have:
    !  m = M₀.₁ + D(z) - 2.5log₁₀(1+z) - K
    ! where D(z)=25+5log(Dₗ) is the regular distance modulus, Dₗ(z)=(1+z)Dᵪ(z) is the luminosity distance, Dᵪ(z) is the comoving
    ! distance, the -2.5log₁₀(1+z) terms accounts for compression of photon frequencies due to redshifting, and K is the
    ! k-correction. Since we do knot know K for each galaxy we neglect it. As such:
    !  D(z) - 2.5log₁₀(1+z) = 25 + 5 log₁₀[Dᵪ(z)] + 2.5log₁₀(1+z) = 25 + 5 log₁₀[Dᵪ(z) √(1+z)] = m - M₀.₁
    ! This is converted to a comoving distance by calling the relevant cosmological conversion function, with the input difference
    ! of magnitudes explicitly noted to contain the -2.5log₁₀(1+z) K-correction factor.
    distanceMaximumRedshift    =self   %cosmologyFunctions_%distanceComoving           (                                                          &
         &                       self  %cosmologyFunctions_%cosmicTime                  (                                                         &
         &                        self %cosmologyFunctions_%expansionFactorFromRedshift  (                                                        &
         &                                                                                                        +self%redshiftMaximum           &
         &                                                                               )                                                        &
         &                                                                              )                                                         &
         &                                                                             )
    if (present(magnitudeAbsolute)) then
       distanceMaximumMagnitude=self   %cosmologyFunctions_%distanceComovingConvert    (                                                          &
            &                                                                           output                   =      distanceTypeComoving    , &
            &                                                                           distanceModulusKCorrected=+self%magnitudeApparentMaximum  &
            &                                                                                                     -     magnitudeAbsolute         &
            &                                                                          )
    else
       distanceMaximumMagnitude=huge(0.0d0)
    end if
    ! Take the smaller of the two distances.
    monteroDorta2009SDSSDistanceMaximum=min(                          &
         &                                  distanceMaximumRedshift , &
         &                                  distanceMaximumMagnitude  &
         &                                 )
    return
  end function monteroDorta2009SDSSDistanceMaximum
