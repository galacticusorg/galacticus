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

!!{
Implements the geometry of the GAMAnear survey used by \cite{kelvin_galaxy_2014-1}.
!!}

  !![
  <surveyGeometry name="surveyGeometryKelvin2014GAMAnear">
   <description>Implements the geometry of the GAMAnear survey of \cite{kelvin_galaxy_2014-1}.</description>
  </surveyGeometry>
  !!]

  type, extends(surveyGeometryBaldry2012GAMA) :: surveyGeometryKelvin2014GAMAnear
     private
     double precision :: distanceMinimumSurvey
   contains
     procedure :: distanceMinimum => kelvin2014GAMAnearDistanceMinimum
     procedure :: distanceMaximum => kelvin2014GAMAnearDistanceMaximum
  end type surveyGeometryKelvin2014GAMAnear

  interface surveyGeometryKelvin2014GAMAnear
     !!{
     Constructors for the \cite{kelvin_galaxy_2014-1} survey geometry class.
     !!}
     module procedure kelvin2014GAMAnearConstructorParameters
     module procedure kelvin2014GAMAnearConstructorInternal
  end interface surveyGeometryKelvin2014GAMAnear

contains

  function kelvin2014GAMAnearConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the \cite{kelvin_galaxy_2014-1} conditional mass function class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (surveyGeometryKelvin2014GAMAnear)                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_

    ! Check and read parameters.
    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    ! Build the object.
    self=surveyGeometryKelvin2014GAMAnear(cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function kelvin2014GAMAnearConstructorParameters

  function kelvin2014GAMAnearConstructorInternal(cosmologyFunctions_) result (self)
    !!{
    Internal constructor for the \cite{kelvin_galaxy_2014-1} conditional mass function class.
    !!}
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    implicit none
    type            (surveyGeometryKelvin2014GAMAnear)                        :: self
    class           (cosmologyFunctionsClass         ), intent(in   ), target :: cosmologyFunctions_
    double precision                                  , parameter             :: redshiftMinimum    =0.025d0
    double precision                                  , parameter             :: redshiftMaximum    =0.060d0
    !![
    <constructorAssign variables="*cosmologyFunctions_"/>
    !!]

    call self%initialize()
    self%distanceMinimumSurvey=self%cosmologyFunctions_%distanceComovingConvert(distanceTypeComoving,redshift=redshiftMinimum)
    self%distanceMaximumSurvey=self%cosmologyFunctions_%distanceComovingConvert(distanceTypeComoving,redshift=redshiftMaximum)
    return
  end function kelvin2014GAMAnearConstructorInternal

  double precision function kelvin2014GAMAnearDistanceMinimum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the minimum distance at which a galaxy is included in the survey.
    !!}
    implicit none
    class           (surveyGeometryKelvin2014GAMAnear), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                                         luminosity, starFormationRate
    integer                                           , intent(in   ), optional :: field
    !$GLC attributes unused :: mass, field, magnitudeAbsolute, luminosity, starFormationRate

    kelvin2014GAMAnearDistanceMinimum=self%distanceMinimumSurvey
    return
  end function kelvin2014GAMAnearDistanceMinimum

  double precision function kelvin2014GAMAnearDistanceMaximum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (surveyGeometryKelvin2014GAMAnear), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: mass           , magnitudeAbsolute, &
         &                                                                         luminosity     , starFormationRate
    integer                                           , intent(in   ), optional :: field
    double precision                                                            :: logarithmicMass

    ! Validate field.
    if (present(field).and.(field < 1 .or. field > 3)) call Error_Report('1 ≤ field ≤ 3 required'//{introspection:location})
    ! Validate arguments.
    if (present(magnitudeAbsolute)) call Error_Report('`magnitudeAbsolute` is not supported'//{introspection:location})
    if (present(luminosity       )) call Error_Report(       '`luminosity` is not supported'//{introspection:location})
    if (present(starFormationRate)) call Error_Report('`starFormationRate` is not supported'//{introspection:location})
    ! Compute the limiting distance. For the GAMAnear sample, all fields are limited to r=19.4
    if (present(mass)) then
       logarithmicMass=log10(mass)
       kelvin2014GAMAnearDistanceMaximum                      &
            & =10.0d0**(                                      &
            &           -0.521147071716417d0                  &
            &           +0.318557607893107d0*logarithmicMass  &
            &          )
       kelvin2014GAMAnearDistanceMaximum=min(kelvin2014GAMAnearDistanceMaximum,self%distanceMaximumSurvey)
    else
       kelvin2014GAMAnearDistanceMaximum=                                      self%distanceMaximumSurvey
    end if
    return
  end function kelvin2014GAMAnearDistanceMaximum
