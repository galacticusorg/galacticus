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
Implements a geometry corresponding to the detectability of classical Local Group galaxies.
!!}


  !![
  <surveyGeometry name="surveyGeometryLocalGroupClassical">
   <description>Implements a geometry corresponding to the detectability of classical Local Group galaxies.</description>
  </surveyGeometry>
  !!]
  type, extends(surveyGeometryMangle) :: surveyGeometryLocalGroupClassical
     private
     double precision :: distanceMaximumSurvey, massThreshold
   contains
     procedure :: fieldCount                => localGroupClassicalFieldCount
     procedure :: distanceMaximum           => localGroupClassicalDistanceMaximum
     procedure :: angularPowerMaximumDegree => localGroupClassicalAngularPowerMaximumDegree
     procedure :: mangleDirectory           => localGroupClassicalMangleDirectory
     procedure :: mangleFiles               => localGroupClassicalMangleFiles
  end type surveyGeometryLocalGroupClassical

  interface surveyGeometryLocalGroupClassical
     !!{
     Constructors for the \cite{baldry_galaxy_2012} survey geometry class.
     !!}
     module procedure localGroupClassicalConstructorParameters
     module procedure localGroupClassicalConstructorInternal
  end interface surveyGeometryLocalGroupClassical

  ! Number of fields.
  integer, parameter :: countFields         =  1

  ! Maximum degree for angular power spectrum
  integer, parameter :: angularPowerMaximumL=360

contains

  function localGroupClassicalConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{surveyGeometryLocalGroupClassical} survey geometry class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (surveyGeometryLocalGroupClassical)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    double precision                                                   :: distanceMaximumSurvey, massThreshold

    !![
    <inputParameter>
      <name>distanceMaximumSurvey</name>
      <source>parameters</source>
      <defaultValue>300.0d-3</defaultValue>
      <description>The maximum distance for the sample of classical Local Group galaxies.</description>
    </inputParameter>
    <inputParameter>
      <name>massThreshold</name>
      <source>parameters</source>
      <defaultValue>1.0d5</defaultValue>
      <description>The minimum stellar mass for a classical Local Group dwarf galaxy.</description>
    </inputParameter>
    !!]
    self=surveyGeometryLocalGroupClassical(distanceMaximumSurvey,massThreshold)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function localGroupClassicalConstructorParameters

  function localGroupClassicalConstructorInternal(distanceMaximumSurvey,massThreshold) result (self)
    !!{
    Internal constructor for the \cite{baldry_galaxy_2012} conditional mass function class.
    !!}
    implicit none
    type            (surveyGeometryLocalGroupClassical)                :: self
    double precision                                   , intent(in   ) :: distanceMaximumSurvey, massThreshold
    !![
    <constructorAssign variables="distanceMaximumSurvey, massThreshold"/>
    !!]

    call self%initialize()
   return
  end function localGroupClassicalConstructorInternal

  integer function localGroupClassicalFieldCount(self)
    !!{
    Return the number of fields in this sample.
    !!}
    implicit none
    class(surveyGeometryLocalGroupClassical), intent(inout) :: self
    !$GLC attributes unused :: self

    localGroupClassicalFieldCount=countFields
    return
  end function localGroupClassicalFieldCount

  double precision function localGroupClassicalDistanceMaximum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    implicit none
    class           (surveyGeometryLocalGroupClassical), intent(inout)           :: self
    double precision                                   , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                                          luminosity, starFormationRate
    integer                                            , intent(in   ), optional :: field
    !$GLC attributes unused :: field

    ! Validate arguments.
    if (present(magnitudeAbsolute)) call Error_Report('`magnitudeAbsolute` is not supported'//{introspection:location})
    if (present(luminosity       )) call Error_Report(       '`luminosity` is not supported'//{introspection:location})
    if (present(starFormationRate)) call Error_Report('`starFormationRate` is not supported'//{introspection:location})
    ! For galaxies above the mass threshold, assume they can be detected out to the maximum specified distance. Galaxies below the
    ! threshold are never detected.
    if (present(mass)) then
       if (mass > self%massThreshold) then
          localGroupClassicalDistanceMaximum=self%distanceMaximumSurvey
       else
          localGroupClassicalDistanceMaximum=0.0d0
       end if
    else
       localGroupClassicalDistanceMaximum   =self%distanceMaximumSurvey
    end if
    return
  end function localGroupClassicalDistanceMaximum

  function localGroupClassicalMangleDirectory(self)
    !!{
    Return the path to the directory containing \gls{mangle} files.
    !!}
    use :: Input_Paths, only : inputPath, pathTypeDataStatic
    implicit none
    class(surveyGeometryLocalGroupClassical), intent(inout) :: self
    type (varying_string                   )                :: localGroupClassicalMangleDirectory
    !$GLC attributes unused :: self

    localGroupClassicalMangleDirectory=inputPath(pathTypeDataStatic)//"surveyGeometry/localGroupClassical/"
    return
  end function localGroupClassicalMangleDirectory

  subroutine localGroupClassicalMangleFiles(self,mangleFiles)
    !!{
    Return a list of \gls{mangle} files.
    !!}
    implicit none
    class(surveyGeometryLocalGroupClassical)                           , intent(inout) :: self
    type (varying_string                   ), allocatable, dimension(:), intent(inout) :: mangleFiles

    allocate(mangleFiles(1))
    mangleFiles(1)=self%mangleDirectory()//"zoneOfAvoidance.ply"
    return
  end subroutine localGroupClassicalMangleFiles

  integer function localGroupClassicalAngularPowerMaximumDegree(self)
    !!{
    Return the maximum degree for which angular power is computed for the Local Group Classical galaxies survey.
    !!}
    implicit none
    class(surveyGeometryLocalGroupClassical), intent(inout) :: self
    !$GLC attributes unused :: self

    localGroupClassicalAngularPowerMaximumDegree=angularPowerMaximumL
    return
  end function localGroupClassicalAngularPowerMaximumDegree

