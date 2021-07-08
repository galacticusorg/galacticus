!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Implements the geometry of the SDSS survey used by \cite{gunawardhana_galaxy_2013}.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <surveyGeometry name="surveyGeometryGunawardhana2013SDSS">
   <description>Implements the geometry of the SDSS survey of \cite{gunawardhana_galaxy_2013}.</description>
  </surveyGeometry>
  !!]
  type, extends(surveyGeometryBernardi2013SDSS) :: surveyGeometryGunawardhana2013SDSS
     private
     class(cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
   contains
     final     ::                    gunawardhana2013SDSSDestructor
     procedure :: distanceMinimum => gunawardhana2013SDSSDistanceMinimum
     procedure :: distanceMaximum => gunawardhana2013SDSSDistanceMaximum
  end type surveyGeometryGunawardhana2013SDSS

  interface surveyGeometryGunawardhana2013SDSS
     !!{
     Constructors for the \cite{gunawardhana_galaxy_2013} survey geometry class.
     !!}
     module procedure gunawardhana2013SDSSConstructorParameters
     module procedure gunawardhana2013SDSSConstructorInternal
  end interface surveyGeometryGunawardhana2013SDSS

contains

  function gunawardhana2013SDSSConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \cite{gunawardhana_galaxy_2013} conditional mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (surveyGeometryGunawardhana2013SDSS)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(cosmologyFunctionsClass           ), pointer       :: cosmologyFunctions_

    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=surveyGeometryGunawardhana2013SDSS(cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function gunawardhana2013SDSSConstructorParameters

  function gunawardhana2013SDSSConstructorInternal(cosmologyFunctions_) result (self)
    !!{
    Default constructor for the \cite{gunawardhana_galaxy_2013} survey geometry class.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type (surveyGeometryGunawardhana2013SDSS)                        :: self
    class(cosmologyFunctionsClass           ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="*cosmologyFunctions_"/>
    !!]

    ! Initialize geometry.
    call self%initialize()
    return
  end function gunawardhana2013SDSSConstructorInternal

  subroutine gunawardhana2013SDSSDestructor(self)
    !!{
    Destructor for the ``gunawardhana2013SDSS'' survey geometry class.
    !!}
    implicit none
    type(surveyGeometryGunawardhana2013SDSS), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine gunawardhana2013SDSSDestructor

  double precision function gunawardhana2013SDSSDistanceMinimum(self,mass,magnitudeAbsolute,luminosity,field)
    !!{
    Compute the minimum distance at which a galaxy is visible.
    !!}
    implicit none
    class           (surveyGeometryGunawardhana2013SDSS), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: mass , magnitudeAbsolute, luminosity
    integer                                             , intent(in   ), optional :: field
    !$GLC attributes unused :: field, mass, luminosity, magnitudeAbsolute

    ! Compute limiting distances. This is due only to the redshift limit.
    gunawardhana2013SDSSDistanceMinimum=self   %cosmologyFunctions_%distanceComoving           (          &
         &                               self  %cosmologyFunctions_%cosmicTime                  (         &
         &                                self %cosmologyFunctions_%expansionFactorFromRedshift  (        &
         &                                                                                        +1.0d-3 &
         &                                                                                       )        &
         &                                                                                      )         &
         &                                                                                     )
    return
  end function gunawardhana2013SDSSDistanceMinimum

  double precision function gunawardhana2013SDSSDistanceMaximum(self,mass,magnitudeAbsolute,luminosity,field)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    use :: Cosmology_Functions_Options     , only : distanceTypeComoving
    use :: Galacticus_Error                , only : Galacticus_Error_Report
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Units       , only : ergs
    implicit none
    class           (surveyGeometryGunawardhana2013SDSS), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: mass                           , magnitudeAbsolute        , &
         &                                                                           luminosity
    integer                                             , intent(in   ), optional :: field
    double precision                                    , parameter               :: fluxLimiting           =1.0d-18 ! W m⁻².
    double precision                                                              :: distanceMaximumRedshift        , distanceMaximumLuminosity, &
         &                                                                           distanceLuminosity
    !$GLC attributes unused :: field, mass, magnitudeAbsolute

    ! Validate input.
    if (.not.present(luminosity)) call Galacticus_Error_Report('luminosity must be supplied '//{introspection:location})
    ! Compute limiting distances. We find the luminosity distance from the supplied luminosity and the limiting flux of the survey.
    distanceLuminosity      =sqrt(                 &
         &                        +luminosity      &
         &                        *ergs            &
         &                        /4.0d0           &
         &                        /Pi              &
         &                        /fluxLimiting    &
         &                        /megaParsec  **2 &
         &                       )
    distanceMaximumRedshift =self   %cosmologyFunctions_%distanceComoving           (                                         &
         &                    self  %cosmologyFunctions_%cosmicTime                  (                                        &
         &                     self %cosmologyFunctions_%expansionFactorFromRedshift  (                                       &
         &                                                                                              1.0d-1                &
         &                                                                            )                                       &
         &                                                                           )                                        &
         &                                                                          )
    distanceMaximumLuminosity=self   %cosmologyFunctions_%distanceComovingConvert    (                                        &
         &                                                                           output            =distanceTypeComoving, &
         &                                                                           distanceLuminosity=distanceLuminosity    &
         &                                                                          )
    ! Take the smaller of the two distances.
    gunawardhana2013SDSSDistanceMaximum=min(                           &
         &                                  distanceMaximumRedshift  , &
         &                                  distanceMaximumLuminosity  &
         &                                 )
    return
  end function gunawardhana2013SDSSDistanceMaximum
