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

!% Implements the geometry of the VIPERS survey used by \cite{davidzon_vimos_2013}.

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !# <surveyGeometry name="surveyGeometryDavidzon2013VIPERS">
  !#  <description>Implements the geometry of the VIPERS survey of \cite{davidzon_vimos_2013}.</description>
  !# </surveyGeometry>
  type, extends(surveyGeometryMangle) :: surveyGeometryDavidzon2013VIPERS
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     integer                                            :: redshiftBin
     double precision                                   :: binDistanceMinimum , binDistanceMaximum
   contains
     final     ::                              davidzon2013VIPERSDestructor
     procedure :: fieldCount                => davidzon2013VIPERSFieldCount
     procedure :: distanceMinimum           => davidzon2013VIPERSDistanceMinimum
     procedure :: distanceMaximum           => davidzon2013VIPERSDistanceMaximum
     procedure :: volumeMaximum             => davidzon2013VIPERSVolumeMaximum
     procedure :: angularPowerMaximumDegree => davidzon2013VIPERSAngularPowerMaximumDegree
     procedure :: mangleDirectory           => davidzon2013VIPERSMangleDirectory
     procedure :: mangleFiles               => davidzon2013VIPERSMangleFiles
   end type surveyGeometryDavidzon2013VIPERS

  interface surveyGeometryDavidzon2013VIPERS
     !% Constructors for the \cite{davidzon_vimos_2013} survey geometry class.
     module procedure davidzon2013VIPERSConstructorParameters
     module procedure davidzon2013VIPERSConstructorInternal
  end interface surveyGeometryDavidzon2013VIPERS

  ! Number of fields.
  integer, parameter :: davidzon2013VIPERSFields        =  1

  ! Angular power spectra.
  integer, parameter :: davidzon2013AngularPowerMaximumL=720

contains

  function davidzon2013VIPERSConstructorParameters(parameters) result(self)
    !% Constructor for the \cite{davidzon_vimos_2013} conditional mass function class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (surveyGeometryDavidzon2013VIPERS)                :: self
    type   (inputParameters                 ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    integer                                                  :: redshiftBin

    ! Check and read parameters.
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !# <inputParameter>
    !#   <name>redshiftBin</name>
    !#   <source>parameters</source>
    !#   <description>The redshift bin (0, 1, 2) of the \cite{davidzon_vimos_2013} mass function to use.</description>
    !# </inputParameter>
    self=surveyGeometryDavidzon2013VIPERS(redshiftBin,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"/>
    return
  end function davidzon2013VIPERSConstructorParameters

  function davidzon2013VIPERSConstructorInternal(redshiftBin,cosmologyFunctions_) result(self)
    !% Generic constructor for the \cite{davidzon_vimos_2013} mass function class.
    use :: Cosmology_Functions        , only : cosmologyFunctionsClass
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    use :: Galacticus_Error           , only : Galacticus_Error_Report
    implicit none
    type            (surveyGeometryDavidzon2013VIPERS)                        :: self
    integer                                           , intent(in   )         :: redshiftBin
    class           (cosmologyFunctionsClass         ), intent(in   ), target :: cosmologyFunctions_
    double precision                                                          :: redshiftMinimum    , redshiftMaximum
    !# <constructorAssign variables="redshiftBin, *cosmologyFunctions_"/>

    ! Find distance limits for this redshift bin.
    select case (redshiftBin)
    case(0)
       redshiftMinimum=0.50d0
       redshiftMaximum=0.60d0
    case(1)
       redshiftMinimum=0.60d0
       redshiftMaximum=0.80d0
    case(2)
       redshiftMinimum=0.80d0
       redshiftMaximum=1.00d0
    case default
       call Galacticus_Error_Report('0≤redshiftBin≤3 is required'//{introspection:location})
    end select
    self%binDistanceMinimum                                                                 &
         & =self%cosmologyFunctions_%distanceComovingConvert(                               &
         &                                                   output  =distanceTypeComoving, &
         &                                                   redshift=redshiftMinimum       &
         &                                                  )
    self%binDistanceMaximum                                                                 &
         & =self%cosmologyFunctions_%distanceComovingConvert(                               &
         &                                                   output  =distanceTypeComoving, &
         &                                                   redshift=redshiftMaximum       &
         &                                                  )
    call self%initialize()
    return
  end function davidzon2013VIPERSConstructorInternal

  subroutine davidzon2013VIPERSDestructor(self)
    !% Destructor for the ``baldry2012GAMA'' survey geometry class.
    implicit none
    type(surveyGeometryDavidzon2013VIPERS), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"/>
    return
  end subroutine davidzon2013VIPERSDestructor

  integer function davidzon2013VIPERSFieldCount(self)
    !% Return the number of fields in this sample.
    implicit none
    class(surveyGeometryDavidzon2013VIPERS), intent(inout) :: self
    !$GLC attributes unused :: self

    davidzon2013VIPERSFieldCount=davidzon2013VIPERSFields
    return
  end function davidzon2013VIPERSFieldCount

  double precision function davidzon2013VIPERSDistanceMinimum(self,mass,magnitudeAbsolute,luminosity,field)
    !% Compute the minimum distance at which a galaxy is included.
    implicit none
    class           (surveyGeometryDavidzon2013VIPERS), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: mass , magnitudeAbsolute, luminosity
    integer                                           , intent(in   ), optional :: field
    !$GLC attributes unused :: mass, field, magnitudeAbsolute, luminosity

    davidzon2013VIPERSDistanceMinimum=self%binDistanceMinimum
    return
  end function davidzon2013VIPERSDistanceMinimum

  double precision function davidzon2013VIPERSDistanceMaximum(self,mass,magnitudeAbsolute,luminosity,field)
    !% Compute the maximum distance at which a galaxy is visible.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (surveyGeometryDavidzon2013VIPERS), intent(inout)           :: self
    double precision                                  , intent(in   ), optional :: mass , magnitudeAbsolute, luminosity
    integer                                           , intent(in   ), optional :: field
    double precision                                                            :: logarithmicMass
    !$GLC attributes unused :: field, magnitudeAbsolute, luminosity

    ! Find the limiting distance for this mass. (See
    ! constraints/dataAnalysis/stellarMassFunctions_VIPERS_z0_1/massDistanceRelation.pl for details.)
    logarithmicMass=log10(mass)
    select case (self%redshiftBin)
    case (0)
       davidzon2013VIPERSDistanceMaximum=3.20663737335189d0+logarithmicMass*(0.0124101908903665d0)
    case (1)
       davidzon2013VIPERSDistanceMaximum=3.14840402683405d0+logarithmicMass*(0.0268494389098537d0)
    case (2)
       davidzon2013VIPERSDistanceMaximum=3.20688538211015d0+logarithmicMass*(0.0273132827274515d0)
    case default
       davidzon2013VIPERSDistanceMaximum=0.0d0
       call Galacticus_Error_Report('invalid redshift bin'//{introspection:location})
    end select
    ! Limit the maximum distance.
    davidzon2013VIPERSDistanceMaximum=min(10.0d0**davidzon2013VIPERSDistanceMaximum,self%binDistanceMaximum)
    return
  end function davidzon2013VIPERSDistanceMaximum

  double precision function davidzon2013VIPERSVolumeMaximum(self,mass,field)
    !% Compute the maximum volume within which a galaxy is visible.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (surveyGeometryDavidzon2013VIPERS), intent(inout)           :: self
    double precision                                  , intent(in   )           :: mass
    integer                                           , intent(in   ), optional :: field

     ! Compute the volume.
    davidzon2013VIPERSVolumeMaximum                              &
         & =max(                                                 &
         &       0.0d0                                         , &
         &       self%solidAngle()                               &
         &      *(                                               &
         &        +self%distanceMaximum   (mass,field=field)**3  &
         &        -self%binDistanceMinimum                  **3  &
         &       )                                               &
         &      /3.0d0                                           &
         &     )
    return
  end function davidzon2013VIPERSVolumeMaximum

  function davidzon2013VIPERSMangleDirectory(self)
    !% Return the path to the directory containing \gls{mangle} files.
    use :: Galacticus_Paths, only : galacticusPath, pathTypeDataStatic
    implicit none
    class(surveyGeometryDavidzon2013VIPERS), intent(inout) :: self
    type (varying_string                  )                :: davidzon2013VIPERSMangleDirectory
    !$GLC attributes unused :: self

    davidzon2013VIPERSMangleDirectory=galacticusPath(pathTypeDataStatic)//"surveyGeometry/VIPERS/"
    return
  end function davidzon2013VIPERSMangleDirectory

  subroutine davidzon2013VIPERSMangleFiles(self,mangleFiles)
    !% Return a list of \gls{mangle} files.
    implicit none
    class(surveyGeometryDavidzon2013VIPERS)                           , intent(inout) :: self
    type (varying_string                  ), allocatable, dimension(:), intent(inout) :: mangleFiles

    allocate(mangleFiles(3))
    mangleFiles=                                                       &
         &      [                                                      &
         &       "+"//self%mangleDirectory()//"/maskCombinedBU.ply:"// &
         &       "-"//self%mangleDirectory()//"/photoW1.BU2.ply:"   // &
         &       "-"//self%mangleDirectory()//"/photoW4.BU2.ply"       &
         &      ]
    return
  end subroutine davidzon2013VIPERSMangleFiles

  integer function davidzon2013VIPERSAngularPowerMaximumDegree(self)
    !% Return the maximum degree for which angular power is computed for the \cite{davidzon_vimos_2013} survey.
    implicit none
    class(surveyGeometryDavidzon2013VIPERS), intent(inout) :: self
    !$GLC attributes unused :: self

    davidzon2013VIPERSAngularPowerMaximumDegree=davidzon2013AngularPowerMaximumL
    return
  end function davidzon2013VIPERSAngularPowerMaximumDegree

