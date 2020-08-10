!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Implements the geometry of the GAMA survey used by \cite{baldry_galaxy_2012}.

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !# <surveyGeometry name="surveyGeometryBaldry2012GAMA">
  !#  <description>Implements the geometry of the GAMA survey of \cite{baldry_galaxy_2012}.</description>
  !# </surveyGeometry>
  type, extends(surveyGeometryMangle) :: surveyGeometryBaldry2012GAMA
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_   => null()
     double precision                                   :: distanceMaximumSurvey
   contains
     final     ::                              baldry2012GAMADestructor
     procedure :: fieldCount                => baldry2012GAMAFieldCount
     procedure :: distanceMaximum           => baldry2012GAMADistanceMaximum
     procedure :: angularPowerMaximumDegree => baldry2012GAMAAngularPowerMaximumDegree
     procedure :: mangleDirectory           => baldry2012GAMAMangleDirectory
     procedure :: mangleFiles               => baldry2012GAMAMangleFiles
  end type surveyGeometryBaldry2012GAMA

  interface surveyGeometryBaldry2012GAMA
     !% Constructors for the \cite{baldry_galaxy_2012} survey geometry class.
     module procedure baldry2012GAMAConstructorParameters
     module procedure baldry2012GAMAConstructorInternal
  end interface surveyGeometryBaldry2012GAMA

  ! Number of fields.
  integer, parameter :: baldry2012GAMAFields              =  3

  ! Maximum degree for angular power spectrum
  integer, parameter :: baldry2012GAMAAngularPowerMaximumL=360

contains

  function baldry2012GAMAConstructorParameters(parameters) result (self)
    !% Constructor for the \cite{baldry_galaxy_2012} conditional mass function class which takes a parameter set as input.
    use :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass
    use :: Input_Parameters   , only : inputParameter    , inputParameters
    implicit none
    type (surveyGeometryBaldry2012GAMA)                :: self
    type (inputParameters             ), intent(inout) :: parameters
    class(cosmologyFunctionsClass     ), pointer       :: cosmologyFunctions_

    ! Check and read parameters.
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    ! Build the object.
    self=surveyGeometryBaldry2012GAMA(cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"/>
    return
  end function baldry2012GAMAConstructorParameters

  function baldry2012GAMAConstructorInternal(cosmologyFunctions_) result (self)
    !% Internal constructor for the \cite{baldry_galaxy_2012} conditional mass function class.
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    implicit none
    type            (surveyGeometryBaldry2012GAMA)                        :: self
    class           (cosmologyFunctionsClass     ), intent(in   ), target :: cosmologyFunctions_
    double precision                              , parameter             :: redshiftMaximum    =0.06d0
    !# <constructorAssign variables="*cosmologyFunctions_"/>

    call self%initialize()
    self%distanceMaximumSurvey=self%cosmologyFunctions_%distanceComovingConvert(distanceTypeComoving,redshift=redshiftMaximum)
   return
  end function baldry2012GAMAConstructorInternal

  subroutine baldry2012GAMADestructor(self)
    !% Destructor for the ``baldry2012GAMA'' survey geometry class.
    implicit none
    type(surveyGeometryBaldry2012GAMA), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"/>
    return
  end subroutine baldry2012GAMADestructor

  integer function baldry2012GAMAFieldCount(self)
    !% Return the number of fields in this sample.
    implicit none
    class(surveyGeometryBaldry2012GAMA), intent(inout) :: self
    !$GLC attributes unused :: self

    baldry2012GAMAFieldCount=baldry2012GAMAFields
    return
  end function baldry2012GAMAFieldCount

  double precision function baldry2012GAMADistanceMaximum(self,mass,magnitudeAbsolute,luminosity,field)
    !% Compute the maximum distance at which a galaxy is visible.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (surveyGeometryBaldry2012GAMA), intent(inout)           :: self
    double precision                              , intent(in   ), optional :: mass           , magnitudeAbsolute, luminosity
    integer                                       , intent(in   ), optional :: field
    double precision                                                        :: logarithmicMass
    !$GLC attributes unused :: magnitudeAbsolute, luminosity

    ! Validate field.
    if (.not.present(field)) call Galacticus_Error_Report('field must be specified'//{introspection:location})
    ! Compute the limiting distance.
    logarithmicMass=log10(mass)
    select case (field)
    case (1,3) ! Fields G09 and G15.
       baldry2012GAMADistanceMaximum                          &
            & =10.0d0**(                                      &
            &           -0.521147071716417d0                  &
            &           +0.318557607893107d0*logarithmicMass  &
            &          )
    case (2)
       baldry2012GAMADistanceMaximum                          &
            & =10.0d0**(                                      &
            &           -0.361147071716369d0                  &
            &           +0.318557607893101d0*logarithmicMass  &
            &          )
    case default
       baldry2012GAMADistanceMaximum=0.0d0
       call Galacticus_Error_Report('1 ≤ field ≤ 3 required'//{introspection:location})
    end select
    baldry2012GAMADistanceMaximum=min(baldry2012GAMADistanceMaximum,self%distanceMaximumSurvey)
    return
  end function baldry2012GAMADistanceMaximum

  function baldry2012GAMAMangleDirectory(self)
    !% Return the path to the directory containing \gls{mangle} files.
    use :: Galacticus_Paths, only : galacticusPath, pathTypeDataStatic
    implicit none
    class(surveyGeometryBaldry2012GAMA), intent(inout) :: self
    type (varying_string              )                :: baldry2012GAMAMangleDirectory
    !$GLC attributes unused :: self

    baldry2012GAMAMangleDirectory=galacticusPath(pathTypeDataStatic)//"surveyGeometry/GAMA/"
    return
  end function baldry2012GAMAMangleDirectory

  subroutine baldry2012GAMAMangleFiles(self,mangleFiles)
    !% Return a list of \gls{mangle} files.
    implicit none
    class(surveyGeometryBaldry2012GAMA)                           , intent(inout) :: self
    type (varying_string              ), allocatable, dimension(:), intent(inout) :: mangleFiles

    allocate(mangleFiles(3))
    mangleFiles=                                                   &
         &      [                                                  &
         &       self%mangleDirectory()//"angularGeometryG09.ply", &
         &       self%mangleDirectory()//"angularGeometryG12.ply", &
         &       self%mangleDirectory()//"angularGeometryG15.ply"  &
         &      ]
    return
  end subroutine baldry2012GAMAMangleFiles

  integer function baldry2012GAMAAngularPowerMaximumDegree(self)
    !% Return the maximum degree for which angular power is computed for the \cite{bernardi_massive_2013} survey.
    implicit none
    class(surveyGeometryBaldry2012GAMA), intent(inout) :: self
    !$GLC attributes unused :: self

    baldry2012GAMAAngularPowerMaximumDegree=baldry2012GAMAAngularPowerMaximumL
    return
  end function baldry2012GAMAAngularPowerMaximumDegree

