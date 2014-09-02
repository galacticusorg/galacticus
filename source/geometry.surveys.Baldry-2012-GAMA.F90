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

!% Implements the geometry of the GAMA survey used by \cite{baldry_galaxy_2012}.
  
  !# <surveyGeometry name="surveyGeometryBaldry2012GAMA">
  !#  <description>Implements the geometry of the GAMA survey of \cite{baldry_galaxy_2012}.</description>
  !# </surveyGeometry>

  use Galacticus_Input_Paths

  type, extends(surveyGeometryMangle) :: surveyGeometryBaldry2012GAMA
     double precision :: distanceMaximumSurvey
   contains
     procedure :: fieldCount                => baldry2012GAMAFieldCount
     procedure :: distanceMaximum           => baldry2012GAMADistanceMaximum
     procedure :: angularPowerMaximumDegree => baldry2012GAMAAngularPowerMaximumDegree
     procedure :: mangleDirectory           => baldry2012GAMAMangleDirectory
     procedure :: mangleFiles               => baldry2012GAMAMangleFiles
  end type surveyGeometryBaldry2012GAMA

  interface surveyGeometryBaldry2012GAMA
     !% Constructors for the \cite{baldry_galaxy_2012} survey geometry class.
     module procedure baldry2012GAMADefaultConstructor
  end interface surveyGeometryBaldry2012GAMA

  ! Number of fields.
  integer, parameter :: baldry2012GAMAFields              =  3
 
  ! Maximum degree for angular power spectrum
  integer, parameter :: baldry2012GAMAAngularPowerMaximumL=360

contains

  function baldry2012GAMADefaultConstructor()
    !% Default constructor for the \cite{baldry_galaxy_2012} conditional mass function class.
    use Cosmology_Functions
    use Cosmology_Functions_Options
    implicit none
    type            (surveyGeometryBaldry2012GAMA)            :: baldry2012GAMADefaultConstructor
    class           (cosmologyFunctionsClass     ), pointer   :: cosmologyFunctions_
    double precision                              , parameter :: redshiftMaximum                 =0.06d0

    cosmologyFunctions_                                      => cosmologyFunctions                         (                                             )
    baldry2012GAMADefaultConstructor%distanceMaximumSurvey   =  cosmologyFunctions_%distanceComovingConvert(distanceTypeComoving,redshift=redshiftMaximum)
    baldry2012GAMADefaultConstructor%solidAnglesInitialized  =.false.
    baldry2012GAMADefaultConstructor%angularPowerInitialized =.false.
    return
  end function baldry2012GAMADefaultConstructor

  integer function baldry2012GAMAFieldCount(self)
    !% Return the number of fields in this sample.
    implicit none
    class(surveyGeometryBaldry2012GAMA), intent(inout) :: self

    baldry2012GAMAFieldCount=baldry2012GAMAFields
    return
  end function baldry2012GAMAFieldCount

  double precision function baldry2012GAMADistanceMaximum(self,mass,field)
    !% Compute the maximum distance at which a galaxy is visible.
    use Galacticus_Error
    implicit none
    class           (surveyGeometryBaldry2012GAMA), intent(inout)           :: self
    double precision                              , intent(in   )           :: mass
    integer                                       , intent(in   ), optional :: field
    double precision                                                        :: logarithmicMass
    
   ! Validate field.
    if (.not.present(field)) call Galacticus_Error_Report('baldry2012GAMADistanceMaximum','field must be specified')
    if (field < 1 .or. field > 3) call Galacticus_Error_Report('baldry2012GAMADistanceMaximum','1 ≤ field ≤ 3 required')
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
    end select
    baldry2012GAMADistanceMaximum=min(baldry2012GAMADistanceMaximum,self%distanceMaximumSurvey)
    return
  end function baldry2012GAMADistanceMaximum

  function baldry2012GAMAMangleDirectory(self)
    !% Return the path to the directory containing \gls{mangle} files.
    implicit none
    class(surveyGeometryBaldry2012GAMA), intent(inout) :: self
    type (varying_string              )                :: baldry2012GAMAMangleDirectory

    baldry2012GAMAMangleDirectory=Galacticus_Input_Path()//"constraints/dataAnalysis/stellarMassFunction_GAMA_z0.03/"
    return
  end function baldry2012GAMAMangleDirectory
  
  subroutine baldry2012GAMAMangleFiles(self,mangleFiles)
    !% Return a list of \gls{mangle} files.
    implicit none
    class(surveyGeometryBaldry2012GAMA)                           , intent(inout) :: self
    type (varying_string              ), allocatable, dimension(:), intent(  out) :: mangleFiles

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

    baldry2012GAMAAngularPowerMaximumDegree=baldry2012GAMAAngularPowerMaximumL
    return
  end function baldry2012GAMAAngularPowerMaximumDegree
  
