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
Implements the geometry of the DES survey for Local Group dwarfs.
!!}


  !![
  <surveyGeometry name="surveyGeometryLocalGroupDES">
   <description>Implements the geometry of the DES survey for Local Group dwarfs.</description>
  </surveyGeometry>
  !!]
  type, extends(surveyGeometryMangle) :: surveyGeometryLocalGroupDES
     private
     double precision :: distanceMaximumSurvey
   contains
     procedure :: fieldCount                => localGroupDESFieldCount
     procedure :: distanceMaximum           => localGroupDESDistanceMaximum
     procedure :: angularPowerMaximumDegree => localGroupDESAngularPowerMaximumDegree
     procedure :: mangleDirectory           => localGroupDESMangleDirectory
     procedure :: mangleFiles               => localGroupDESMangleFiles
  end type surveyGeometryLocalGroupDES

  interface surveyGeometryLocalGroupDES
     !!{
     Constructors for the \refClass{surveyGeometryLocalGroupDES} survey geometry class.
     !!}
     module procedure localGroupDESConstructorParameters
     module procedure localGroupDESConstructorInternal
  end interface surveyGeometryLocalGroupDES

  ! Angular power spectra.
  integer, parameter :: angularPowerMaximumL=360

contains

  function localGroupDESConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{surveyGeometryLocalGroupDES} conditional mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (surveyGeometryLocalGroupDES)                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    double precision                                             :: distanceMaximumSurvey

    !![
    <inputParameter>
      <name>distanceMaximumSurvey</name>
      <source>parameters</source>
      <defaultValue>300.0d-3</defaultValue>
      <description>The maximum distance at which galaxies are to be included in the survey.</description>
    </inputParameter>
    !!]
    self=surveyGeometryLocalGroupDES(distanceMaximumSurvey)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function localGroupDESConstructorParameters

  function localGroupDESConstructorInternal(distanceMaximumSurvey) result (self)
    !!{
    Internal constructor for the \refClass{surveyGeometryLocalGroupDES} survey geometry class.
    !!}
    implicit none
    type            (surveyGeometryLocalGroupDES)                :: self
    double precision                             , intent(in   ) :: distanceMaximumSurvey
    !![
    <constructorAssign variables="distanceMaximumSurvey"/>
    !!]

    call self%initialize()
    return
  end function localGroupDESConstructorInternal

  integer function localGroupDESFieldCount(self)
    !!{
    Return the number of fields in this sample.
    !!}
    implicit none
    class(surveyGeometryLocalGroupDES), intent(inout) :: self
    !$GLC attributes unused :: self

    localGroupDESFieldCount=1
    return
  end function localGroupDESFieldCount

  double precision function localGroupDESDistanceMaximum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    implicit none
    class           (surveyGeometryLocalGroupDES), intent(inout)           :: self
    double precision                             , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                                    luminosity, starFormationRate
    integer                                      , intent(in   ), optional :: field
    !$GLC attributes unused :: field

    ! Validate arguments.
    if (present(magnitudeAbsolute)) call Error_Report('`magnitudeAbsolute` is not supported'//{introspection:location})
    if (present(luminosity       )) call Error_Report(       '`luminosity` is not supported'//{introspection:location})
    if (present(starFormationRate)) call Error_Report('`starFormationRate` is not supported'//{introspection:location})
    ! Find the limiting distance for this mass completeness limits. The following functional form should be considered to be
    ! approximate at best. It was derived by fitting a polynomial in stellar mass (assuming a mass-to-light ratio of 1 in Solar
    ! units in the V band) to the detection efficiencies reported by Drlica-Wagner et al. (2015, ApJ, 813, 109;
    ! https://ui.adsabs.harvard.edu/abs/2015ApJ...813..109D) - with galaxies where the efficiency was 1 excluded (since by
    ! necessity the efficiency is truncated at 1), and then solving for distance as a function of mass for the point at which the
    ! detection efficiency is 90%.
    if (present(mass)) then
       localGroupDESDistanceMaximum=min(11.3d-3*(mass/100.0d0)**0.48d0,self%distanceMaximumSurvey)
    else
       localGroupDESDistanceMaximum=                                   self%distanceMaximumSurvey
    end if
    return
  end function localGroupDESDistanceMaximum

  integer function localGroupDESAngularPowerMaximumDegree(self)
    !!{
    Return the maximum degree for which angular power is computed for the {\normalfont \ttfamily localGroupDES} survey.
    !!}
    implicit none
    class(surveyGeometryLocalGroupDES), intent(inout) :: self
    !$GLC attributes unused :: self

    localGroupDESAngularPowerMaximumDegree=angularPowerMaximumL
    return
  end function localGroupDESAngularPowerMaximumDegree

  function localGroupDESMangleDirectory(self)
    !!{
    Return the path to the directory containing \gls{mangle} files.
    !!}
    use :: Input_Paths, only : inputPath, pathTypeDataStatic
    implicit none
    class(surveyGeometryLocalGroupDES), intent(inout) :: self
    type (varying_string             )                :: localGroupDESMangleDirectory
    !$GLC attributes unused :: self

    localGroupDESMangleDirectory=inputPath(pathTypeDataStatic)//"surveyGeometry/darkEnergySurvey/"
    return
  end function localGroupDESMangleDirectory

  subroutine localGroupDESMangleFiles(self,mangleFiles)
    !!{
    Return a list of \gls{mangle} files.
    !!}
    implicit none
    class(surveyGeometryLocalGroupDES)                           , intent(inout) :: self
    type (varying_string             ), allocatable, dimension(:), intent(inout) :: mangleFiles

    allocate(mangleFiles(1))
    mangleFiles(1)=self%mangleDirectory()//"darkEnergySurvey.ply"
    return
  end subroutine localGroupDESMangleFiles
