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
Implements the survey geometry of the SDSS sample used by \cite{hearin_dark_2013}.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <surveyGeometry name="surveyGeometryHearin2014SDSS">
   <description>Implements the survey geometry of the SDSS sample used by \cite{hearin_dark_2013}.</description>
  </surveyGeometry>
  !!]
  type, extends(surveyGeometryBernardi2013SDSS) :: surveyGeometryHearin2014SDSS
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_  => null()
     double precision                                   :: distanceMinimumLimit          , distanceMaximumLimit   , &
          &                                                massPrevious                  , distanceMaximumPrevious
   contains
     final     ::                      hearin2014SDSSDestructor
     procedure :: distanceMinimum   => hearin2014SDSSDistanceMinimum
     procedure :: distanceMaximum   => hearin2014SDSSDistanceMaximum
  end type surveyGeometryHearin2014SDSS

  interface surveyGeometryHearin2014SDSS
     !!{
     Constructors for the \cite{hearin_dark_2013} survey geometry class.
     !!}
     module procedure hearin2014SDSSConstructorParameters
     module procedure hearin2014SDSSConstructorInternal
  end interface surveyGeometryHearin2014SDSS

  ! Redshift limits.
  double precision, parameter :: hearing2014RedshiftMinimum=0.020d0
  double precision, parameter :: hearing2014RedshiftMaximum=0.068d0

contains

  function hearin2014SDSSConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the \cite{hearin_dark_2013} conditional mass function class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (surveyGeometryHearin2014SDSS)                :: self
    type (inputParameters             ), intent(inout) :: parameters
    class(cosmologyFunctionsClass     ), pointer       :: cosmologyFunctions_

    ! Check and read parameters.
    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    ! Build the object.
    self=surveyGeometryHearin2014SDSS(cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function hearin2014SDSSConstructorParameters

  function hearin2014SDSSConstructorInternal(cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \cite{hearin_dark_2013} conditional mass function class.
    !!}
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    implicit none
    type (surveyGeometryHearin2014SDSS)                        :: self
    class(cosmologyFunctionsClass     ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="*cosmologyFunctions_"/>
    !!]

    self%solidAnglesInitialized =.false.
    self%angularPowerInitialized=.false.
    self%windowInitialized      =.false.
    self%massPrevious           =-1.0d0
    ! Get the default cosmology functions object.
    self%distanceMinimumLimit                                                                     &
         & =self%cosmologyFunctions_%distanceComovingConvert(                                     &
         &                                                   output  =distanceTypeComoving      , &
         &                                                   redshift=hearing2014RedshiftMinimum  &
         &                                                  )
    self%distanceMaximumLimit                                                                     &
         & =self%cosmologyFunctions_%distanceComovingConvert(                                     &
         &                                                   output  =distanceTypeComoving      , &
         &                                                   redshift=hearing2014RedshiftMaximum  &
         &                                                  )
    return
  end function hearin2014SDSSConstructorInternal

  subroutine hearin2014SDSSDestructor(self)
    !!{
    Destructor for the \refClass{surveyGeometryHearin2014SDSS} survey geometry class.
    !!}
    implicit none
    type(surveyGeometryHearin2014SDSS), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine hearin2014SDSSDestructor

  double precision function hearin2014SDSSDistanceMinimum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the minimum distance at which a galaxy is visible.
    !!}
    implicit none
    class           (surveyGeometryHearin2014SDSS), intent(inout)           :: self
    double precision                              , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                                     luminosity, starFormationRate
    integer                                       , intent(in   ), optional :: field
    !$GLC attributes unused :: mass, field, magnitudeAbsolute, luminosity, starFormationRate

    hearin2014SDSSDistanceMinimum=self%distanceMinimumLimit
    return
  end function hearin2014SDSSDistanceMinimum

  double precision function hearin2014SDSSDistanceMaximum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    implicit none
    class           (surveyGeometryHearin2014SDSS), intent(inout)           :: self
    double precision                              , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                                     luminosity, starFormationRate
    integer                                       , intent(in   ), optional :: field
    !$GLC attributes unused :: field
        ! Validate arguments.
    if (present(magnitudeAbsolute)) call Error_Report('`magnitudeAbsolute` is not supported'//{introspection:location})
    if (present(luminosity       )) call Error_Report(       '`luminosity` is not supported'//{introspection:location})
    if (present(starFormationRate)) call Error_Report('`starFormationRate` is not supported'//{introspection:location})
    if (present(mass)) then
       if (mass /= self%massPrevious)                                                                                 &
            & self%distanceMaximumPrevious=min(                                                                       &
            &                                  self%surveyGeometryBernardi2013SDSS%distanceMaximum(mass,field=field), &
            &                                  self%distanceMaximumLimit                                              &
            &                                 )
       hearin2014SDSSDistanceMaximum=self%distanceMaximumPrevious
    else
       hearin2014SDSSDistanceMaximum=self%distanceMaximumLimit
    end if
    return
  end function hearin2014SDSSDistanceMaximum
