!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Implements the survey geometry of the SDSS sample used by \cite{hearin_dark_2013}.

  !# <surveyGeometry name="surveyGeometryHearin2014SDSS">
  !#  <description>Implements the survey geometry of the SDSS sample used by \cite{hearin_dark_2013}.</description>
  !# </surveyGeometry>

  type, extends(surveyGeometryBernardi2013SDSS) :: surveyGeometryHearin2014SDSS
     double precision :: distanceMinimumLimit, distanceMaximumLimit, massPrevious, distanceMaximumPrevious
   contains
     procedure :: distanceMinimum   => hearin2014SDSSDistanceMinimum
     procedure :: distanceMaximum   => hearin2014SDSSDistanceMaximum
  end type surveyGeometryHearin2014SDSS

  interface surveyGeometryHearin2014SDSS
     !% Constructors for the \cite{hearin_dark_2013} survey geometry class.
     module procedure hearin2014SDSSDefaultConstructor
  end interface surveyGeometryHearin2014SDSS

  ! Redshift limits.
  double precision, parameter :: hearing2014RedshiftMinimum=0.020d0
  double precision, parameter :: hearing2014RedshiftMaximum=0.068d0

contains

  function hearin2014SDSSDefaultConstructor()
    !% Default constructor for the \cite{hearin_dark_2013} conditional mass function class.
    use Cosmology_Functions
    use Cosmology_Functions_Options
    implicit none
    type (surveyGeometryHearin2014SDSS)          :: hearin2014SDSSDefaultConstructor
    class(cosmologyFunctionsClass     ), pointer :: cosmologyFunctions_
    
    hearin2014SDSSDefaultConstructor%solidAnglesInitialized =.false.
    hearin2014SDSSDefaultConstructor%angularPowerInitialized=.false.
    hearin2014SDSSDefaultConstructor%windowInitialized      =.false.
    hearin2014SDSSDefaultConstructor%massPrevious           =-1.0d0
    ! Get the default cosmology functions object.
    cosmologyFunctions_ => cosmologyFunctions()
    hearin2014SDSSDefaultConstructor%distanceMinimumLimit                                    &
         & =cosmologyFunctions_%distanceComovingConvert(                                     &
         &                                              output  =distanceTypeComoving      , &
         &                                              redshift=hearing2014RedshiftMinimum  &
         &                                             )
    hearin2014SDSSDefaultConstructor%distanceMaximumLimit                                    &
         & =cosmologyFunctions_%distanceComovingConvert(                                     &
         &                                              output  =distanceTypeComoving      , &
         &                                              redshift=hearing2014RedshiftMaximum  &
         &                                             )
    return
  end function hearin2014SDSSDefaultConstructor

  double precision function hearin2014SDSSDistanceMinimum(self,mass,field)
    !% Compute the minimum distance at which a galaxy is visible.
    implicit none
    class           (surveyGeometryHearin2014SDSS), intent(inout)           :: self
    double precision                              , intent(in   )           :: mass
    integer                                       , intent(in   ), optional :: field

    hearin2014SDSSDistanceMinimum=self%distanceMinimumLimit
    return
  end function hearin2014SDSSDistanceMinimum

  double precision function hearin2014SDSSDistanceMaximum(self,mass,field)
    !% Compute the maximum distance at which a galaxy is visible.
    implicit none
    class           (surveyGeometryHearin2014SDSS), intent(inout)           :: self
    double precision                              , intent(in   )           :: mass
    integer                                       , intent(in   ), optional :: field

    if (mass /= self%massPrevious)                                                                           &
         & self%distanceMaximumPrevious=min(                                                                 &
         &                                  self%surveyGeometryBernardi2013SDSS%distanceMaximum(mass,field), &
         &                                  self%distanceMaximumLimit                                        &
         &                                 )
    hearin2014SDSSDistanceMaximum=self%distanceMaximumPrevious
    return
  end function hearin2014SDSSDistanceMaximum
