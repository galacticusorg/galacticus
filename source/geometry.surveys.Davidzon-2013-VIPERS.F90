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

!% Implements the geometry of the VIPERS survey used by \cite{davidzon_vimos_2013}.
  
  !# <surveyGeometry name="surveyGeometryDavidzon2013VIPERS">
  !#  <description>Implements the geometry of the VIPERS survey of \cite{davidzon_vimos_2013}.</description>
  !# </surveyGeometry>

  use Galacticus_Input_Paths

  type, extends(surveyGeometryMangle) :: surveyGeometryDavidzon2013VIPERS
     integer          :: redshiftBin
     double precision :: binDistanceMinimum, binDistanceMaximum
   contains
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
     module procedure davidzon2013VIPERSConstructor
     module procedure davidzon2013VIPERSDefaultConstructor
  end interface surveyGeometryDavidzon2013VIPERS

  ! Paths and file names for mangle polygon files.
  integer, parameter :: davidzon2013VIPERSFields        =1
  logical            :: davidzon2013VIPERSBinInitialized=.false.
  integer            :: davidzon2013VIPERSRedshiftBin

  ! Angular power spectra.
  integer, parameter :: davidzon2013AngularPowerMaximumL=720

contains

  function davidzon2013VIPERSDefaultConstructor()
    !% Default constructor for the \cite{davidzon_vimos_2013} conditional mass function class.
    use Input_Parameters
    implicit none
    type(surveyGeometryDavidzon2013VIPERS) :: davidzon2013VIPERSDefaultConstructor

    if (.not.davidzon2013VIPERSBinInitialized) then
       !$omp critical(davidzon2013VIPERSBinInitialize)
       if (.not.davidzon2013VIPERSBinInitialized) then
          ! Get the redshift bin to use.
          !@ <inputParameter>
          !@   <name>davidzon2013VIPERSRedshiftBin</name>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The redshift bin (0, 1, 2) of the \cite{davidzon_vimos_2013} mass function to use.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>starFormation</group>
          !@ </inputParameter>
          call Get_Input_Parameter('davidzon2013VIPERSRedshiftBin',davidzon2013VIPERSRedshiftBin)
          davidzon2013VIPERSBinInitialized=.true.
       end if
       !$omp end critical(davidzon2013VIPERSBinInitialize)
    end if
    davidzon2013VIPERSDefaultConstructor=davidzon2013VIPERSConstructor(davidzon2013VIPERSRedshiftBin)
   return
  end function davidzon2013VIPERSDefaultConstructor

  function davidzon2013VIPERSConstructor(redshiftBin)
    !% Generic constructor for the \cite{davidzon_vimos_2013} mass function class.
    use Galacticus_Error
    use Input_Parameters
    use Cosmology_Functions
    use Cosmology_Functions_Options
    implicit none
    type            (surveyGeometryDavidzon2013VIPERS)                :: davidzon2013VIPERSConstructor
    integer                                           , intent(in   ) :: redshiftBin
    class           (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    double precision                                                  :: redshiftMinimum               , redshiftMaximum

    ! Find distance limits for this redshift bin.
    davidzon2013VIPERSConstructor%redshiftBin=redshiftBin
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
       call Galacticus_Error_Report('davidzon2013VIPERSConstructor','0≤redshiftBin≤3 is required')
    end select
    cosmologyFunctions_ => cosmologyFunctions()
    davidzon2013VIPERSConstructor%binDistanceMinimum                                     &
         & =cosmologyFunctions_%distanceComovingConvert(                                  &
         &                                                 output  =distanceTypeComoving, &
         &                                                 redshift=redshiftMinimum       &
         &                                                )
    davidzon2013VIPERSConstructor%binDistanceMaximum                                     &
         & =cosmologyFunctions_%distanceComovingConvert(                                  &
         &                                                 output  =distanceTypeComoving, &
         &                                                 redshift=redshiftMaximum       &
         &                                                )
    davidzon2013VIPERSConstructor%solidAnglesInitialized =.false.
    davidzon2013VIPERSConstructor%angularPowerInitialized=.false.
    davidzon2013VIPERSConstructor%windowInitialized      =.false.
   return
  end function davidzon2013VIPERSConstructor
  
  integer function davidzon2013VIPERSFieldCount(self)
    !% Return the number of fields in this sample.
    implicit none
    class(surveyGeometryDavidzon2013VIPERS), intent(inout) :: self

    davidzon2013VIPERSFieldCount=davidzon2013VIPERSFields
    return
  end function davidzon2013VIPERSFieldCount

  double precision function davidzon2013VIPERSDistanceMinimum(self,mass,field)
    !% Compute the minimum distance at which a galaxy is included.
    implicit none
    class           (surveyGeometryDavidzon2013VIPERS), intent(inout)           :: self
    double precision                                   , intent(in   )           :: mass
    integer                                            , intent(in   ), optional :: field

    davidzon2013VIPERSDistanceMinimum=self%binDistanceMinimum
    return
  end function davidzon2013VIPERSDistanceMinimum

  double precision function davidzon2013VIPERSDistanceMaximum(self,mass,field)
    !% Compute the maximum distance at which a galaxy is visible.
    use Cosmology_Functions
    use Cosmology_Functions_Options
    use Galacticus_Error
    implicit none
    class           (surveyGeometryDavidzon2013VIPERS), intent(inout)           :: self
    double precision                                  , intent(in   )           :: mass
    integer                                           , intent(in   ), optional :: field
    class           (cosmologyFunctionsClass         ), pointer                 :: cosmologyFunctions_
    double precision                                                            :: redshift           , logarithmicMass
    
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
    end select
    ! Limit the maximum distance.
    davidzon2013VIPERSDistanceMaximum=min(10.0d0**davidzon2013VIPERSDistanceMaximum,self%binDistanceMaximum)
    return
  end function davidzon2013VIPERSDistanceMaximum

  double precision function davidzon2013VIPERSVolumeMaximum(self,mass,field)
    !% Compute the maximum volume within which a galaxy is visible.
    use Galacticus_Error
    implicit none
    class           (surveyGeometryDavidzon2013VIPERS), intent(inout)           :: self
    double precision                                  , intent(in   )           :: mass
    integer                                           , intent(in   ), optional :: field

     ! Compute the volume.
    davidzon2013VIPERSVolumeMaximum                        &
         & =max(                                           &
         &       0.0d0                                   , &
         &       self%solidAngle()                         &
         &      *(                                         &
         &        +self%distanceMaximum   (mass,field)**3  &
         &        -self%binDistanceMinimum            **3  &
         &       )                                         &
         &      /3.0d0                                     &
         &     )
    return
  end function davidzon2013VIPERSVolumeMaximum

  function davidzon2013VIPERSMangleDirectory(self)
    !% Return the path to the directory containing \gls{mangle} files.
    implicit none
    class(surveyGeometryDavidzon2013VIPERS), intent(inout) :: self
    type (varying_string                  )                :: davidzon2013VIPERSMangleDirectory

    davidzon2013VIPERSMangleDirectory=Galacticus_Input_Path()//"constraints/dataAnalysis/stellarMassFunctions_VIPERS_z0_1/"
    return
  end function davidzon2013VIPERSMangleDirectory
  
  subroutine davidzon2013VIPERSMangleFiles(self,mangleFiles)
    !% Return a list of \gls{mangle} files.
    implicit none
    class(surveyGeometryDavidzon2013VIPERS)                           , intent(inout) :: self
    type (varying_string                  ), allocatable, dimension(:), intent(  out) :: mangleFiles
    
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

    davidzon2013VIPERSAngularPowerMaximumDegree=davidzon2013AngularPowerMaximumL
    return
  end function davidzon2013VIPERSAngularPowerMaximumDegree
  
