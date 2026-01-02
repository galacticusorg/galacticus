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
Implements the geometry of the SDSS survey with a depth for Local Group dwarf detection.
!!}

  !![
  <surveyGeometry name="surveyGeometryLocalGroupSDSS">
   <description>Implements the geometry of the SDSS survey with a depth for Local Group dwarf detection.</description>
  </surveyGeometry>
  !!]
  type, extends(surveyGeometryBernardi2013SDSS) :: surveyGeometryLocalGroupSDSS
     private
     double precision :: distanceMaximumSurvey
   contains
     procedure :: distanceMaximum => localGroupSDSSDistanceMaximum
  end type surveyGeometryLocalGroupSDSS

  interface surveyGeometryLocalGroupSDSS
     !!{
     Constructors for the \cite{bernardi_massive_2013} survey geometry class.
     !!}
     module procedure localGroupSDSSConstructorParameters
     module procedure localGroupSDSSConstructorInternal
  end interface surveyGeometryLocalGroupSDSS

contains

  function localGroupSDSSConstructorParameters(parameters) result (self)
    !!{
    Constructor for the Local Group SDSS survey geometry class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (surveyGeometryLocalGroupSDSS)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    double precision                                              :: distanceMaximumSurvey

    !![
    <inputParameter>
      <name>distanceMaximumSurvey</name>
      <source>parameters</source>
      <defaultValue>300.0d-3</defaultValue>
      <description>The maximum distance at which galaxies are to be included in the survey.</description>
    </inputParameter>
    !!]
    self=surveyGeometryLocalGroupSDSS(distanceMaximumSurvey)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function localGroupSDSSConstructorParameters

  function localGroupSDSSConstructorInternal(distanceMaximumSurvey) result (self)
    !!{
    Internal constructor for the Local Group SDSS survey geometry class
    !!}
    implicit none
    type            (surveyGeometryLocalGroupSDSS)                :: self
    double precision                              , intent(in   ) :: distanceMaximumSurvey
    !![
    <constructorAssign variables="distanceMaximumSurvey"/>
    !!]

    call self%initialize()
    return
  end function localGroupSDSSConstructorInternal

  double precision function localGroupSDSSDistanceMaximum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    implicit none
    class           (surveyGeometryLocalGroupSDSS), intent(inout)           :: self
    double precision                              , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                                     luminosity, starFormationRate
    integer                                       , intent(in   ), optional :: field
    !$GLC attributes unused :: field

    ! Validate arguments.
    if (present(magnitudeAbsolute)) call Error_Report('`magnitudeAbsolute` is not supported'//{introspection:location})
    if (present(luminosity       )) call Error_Report(       '`luminosity` is not supported'//{introspection:location})
    if (present(starFormationRate)) call Error_Report('`starFormationRate` is not supported'//{introspection:location})
    ! Find the limiting distance for this mass completeness limits. We adopt the model of Kim, Peter & Hargis (2018,
    ! Phys. Rev. Lett. 121, 211302; https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.211302), assuming a
    ! luminosity-to-mass ratio of 1 in Solar units. (Their model is in turn based on the results of Walsh, Willman & Jerjen, 2008,
    ! AJ, 137, 1; https://iopscience.iop.org/article/10.1088/0004-6256/137/1/450/meta.)
    if (present(mass)) then
       localGroupSDSSDistanceMaximum=min(15.7d-3*(mass/100.0d0)**0.51d0,self%distanceMaximumSurvey)
    else
       localGroupSDSSDistanceMaximum=                                   self%distanceMaximumSurvey
    end if
    return
  end function localGroupSDSSDistanceMaximum
