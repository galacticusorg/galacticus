!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Implements the geometry of the GAMAnear survey used by \cite{kelvin_galaxy_2014-1}.
  
  !# <surveyGeometry name="surveyGeometryKelvin2014GAMAnear">
  !#  <description>Implements the geometry of the GAMAnear survey of \cite{kelvin_galaxy_2014-1}.</description>
  !# </surveyGeometry>

  use Galacticus_Input_Paths

  type, extends(surveyGeometryBaldry2012GAMA) :: surveyGeometryKelvin2014GAMAnear
     double precision :: distanceMinimumSurvey
   contains
     procedure :: distanceMinimum => kelvin2014GAMAnearDistanceMinimum
     procedure :: distanceMaximum => kelvin2014GAMAnearDistanceMaximum
  end type surveyGeometryKelvin2014GAMAnear

  interface surveyGeometryKelvin2014GAMAnear
     !% Constructors for the \cite{kelvin_galaxy_2014-1} survey geometry class.
     module procedure kelvin2014GAMAnearDefaultConstructor
  end interface surveyGeometryKelvin2014GAMAnear

contains

  function kelvin2014GAMAnearDefaultConstructor()
    !% Default constructor for the \cite{kelvin_galaxy_2014-1} conditional mass function class.
    use Cosmology_Functions
    use Cosmology_Functions_Options
    implicit none
    type            (surveyGeometryKelvin2014GAMAnear)            :: kelvin2014GAMAnearDefaultConstructor
    class           (cosmologyFunctionsClass         ), pointer   :: cosmologyFunctions_
    double precision                                  , parameter :: redshiftMinimum                     =0.025d0
    double precision                                  , parameter :: redshiftMaximum                     =0.060d0

    cosmologyFunctions_=> cosmologyFunctions()
    kelvin2014GAMAnearDefaultConstructor%distanceMinimumSurvey   =  &
         & cosmologyFunctions_%distanceComovingConvert(distanceTypeComoving,redshift=redshiftMinimum)
    kelvin2014GAMAnearDefaultConstructor%distanceMaximumSurvey   =  &
         & cosmologyFunctions_%distanceComovingConvert(distanceTypeComoving,redshift=redshiftMaximum)
    ! Initialize state.
    kelvin2014GAMAnearDefaultConstructor%solidAnglesInitialized  =.false.
    kelvin2014GAMAnearDefaultConstructor%angularPowerInitialized =.false.
    kelvin2014GAMAnearDefaultConstructor%windowInitialized       =.false.
   return
  end function kelvin2014GAMAnearDefaultConstructor

  double precision function kelvin2014GAMAnearDistanceMinimum(self,mass,field)
    !% Compute the minimum distance at which a galaxy is included in the survey.
    use Galacticus_Error
    implicit none
    class           (surveyGeometryKelvin2014GAMAnear), intent(inout)           :: self
    double precision                                  , intent(in   )           :: mass
    integer                                           , intent(in   ), optional :: field
    !GCC$ attributes unused :: mass, field
    
    kelvin2014GAMAnearDistanceMinimum=self%distanceMinimumSurvey
    return
  end function kelvin2014GAMAnearDistanceMinimum
  
  double precision function kelvin2014GAMAnearDistanceMaximum(self,mass,field)
    !% Compute the maximum distance at which a galaxy is visible.
    use Galacticus_Error
    implicit none
    class           (surveyGeometryKelvin2014GAMAnear), intent(inout)           :: self
    double precision                                  , intent(in   )           :: mass
    integer                                           , intent(in   ), optional :: field
    double precision                                                            :: logarithmicMass
    
    ! Validate field.
    if (present(field).and.(field < 1 .or. field > 3)) call Galacticus_Error_Report('kelvin2014GAMAnearDistanceMaximum','1 ≤ field ≤ 3 required')
    ! Compute the limiting distance. For the GAMAnear sample, all fields are limited to r=19.4
    logarithmicMass=log10(mass)
    kelvin2014GAMAnearDistanceMaximum                      &
         & =10.0d0**(                                      &
         &           -0.521147071716417d0                  &
         &           +0.318557607893107d0*logarithmicMass  &
         &          )    
    kelvin2014GAMAnearDistanceMaximum=min(kelvin2014GAMAnearDistanceMaximum,self%distanceMaximumSurvey)
    return
  end function kelvin2014GAMAnearDistanceMaximum
