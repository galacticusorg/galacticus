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

!% Implements the geometry of the ZFOURGE survey used by \cite{tomczak_galaxy_2014}.
  
  use Galacticus_Input_Paths
  use Cosmology_Functions

  !# <surveyGeometry name="surveyGeometryTomczak2014ZFOURGE">
  !#  <description>Implements the geometry of the ZFOURGE survey of \cite{tomczak_galaxy_2014}.</description>
  !# </surveyGeometry>
  type, extends(surveyGeometryMangle) :: surveyGeometryTomczak2014ZFOURGE
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_
     integer                                            :: redshiftBin
     double precision                                   :: binDistanceMinimum, binDistanceMaximum
   contains
     final     ::                              tomczak2014ZFOURGEDestructor
     procedure :: fieldCount                => tomczak2014ZFOURGEFieldCount
     procedure :: distanceMinimum           => tomczak2014ZFOURGEDistanceMinimum
     procedure :: distanceMaximum           => tomczak2014ZFOURGEDistanceMaximum
     procedure :: volumeMaximum             => tomczak2014ZFOURGEVolumeMaximum
     procedure :: angularPowerMaximumDegree => tomczak2014ZFOURGEAngularPowerMaximumDegree
     procedure :: mangleDirectory           => tomczak2014ZFOURGEMangleDirectory
     procedure :: mangleFiles               => tomczak2014ZFOURGEMangleFiles
   end type surveyGeometryTomczak2014ZFOURGE

  interface surveyGeometryTomczak2014ZFOURGE
     !% Constructors for the \cite{tomczak_galaxy_2014} survey geometry class.
     module procedure tomczak2014ZFOURGEConstructorParameters
     module procedure tomczak2014ZFOURGEConstructorInternal
  end interface surveyGeometryTomczak2014ZFOURGE

  ! Paths and file names for mangle polygon files.
  integer, parameter :: tomczak2014ZFOURGEFields       =   2

  ! Angular power spectra.
  integer, parameter :: tomczak2014AngularPowerMaximumL=5000

contains

  function tomczak2014ZFOURGEConstructorParameters(parameters) result(self)
    !% Default constructor for the \cite{tomczak_galaxy_2014} conditional mass function class.
    use Input_Parameters2
    implicit none
    type   (surveyGeometryTomczak2014ZFOURGE)                :: self
    type   (inputParameters                 ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    integer                                                  :: redshiftBin
    
    ! Check and read parameters.
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !# <inputParameter>
    !#   <name>redshiftBin</name>
    !#   <source>parameters</source>
    !#   <description>The redshift bin (0, 1, 2, 3, 4, 5, 6, or 7) of the \cite{tomczak_galaxy_2014} mass function to use.</description>
    !#   <type>integer</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    self=surveyGeometryTomczak2014ZFOURGE(redshiftBin,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function tomczak2014ZFOURGEConstructorParameters

  function tomczak2014ZFOURGEConstructorInternal(redshiftBin,cosmologyFunctions_) result(self)
    !% Generic constructor for the \cite{tomczak_galaxy_2014} mass function class.
    use Galacticus_Error
    use Cosmology_Functions_Options
    implicit none
    type            (surveyGeometryTomczak2014ZFOURGE)                        :: self
    integer                                           , intent(in   )         :: redshiftBin
    class           (cosmologyFunctionsClass         ), intent(in   ), target :: cosmologyFunctions_
    double precision                                                          :: redshiftMinimum    , redshiftMaximum
    !# <constructorAssign variables="redshiftBin, *cosmologyFunctions_"/>

    ! Find distance limits for this redshift bin.
    select case (redshiftBin)
    case(0)
       redshiftMinimum=0.20d0
       redshiftMaximum=0.50d0
    case(1)
       redshiftMinimum=0.50d0
       redshiftMaximum=0.75d0
    case(2)
       redshiftMinimum=0.75d0
       redshiftMaximum=1.00d0
    case(3)
       redshiftMinimum=1.00d0
       redshiftMaximum=1.25d0
    case(4)
       redshiftMinimum=1.25d0
       redshiftMaximum=1.50d0
    case(5)
       redshiftMinimum=1.50d0
       redshiftMaximum=2.00d0
    case(6)
       redshiftMinimum=2.00d0
       redshiftMaximum=2.50d0
    case(7)
       redshiftMinimum=2.50d0
       redshiftMaximum=3.00d0
    case default
       call Galacticus_Error_Report('tomczak2014ZFOURGEConstructorInternal','0≤redshiftBin≤7 is required')
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
  end function tomczak2014ZFOURGEConstructorInternal
  
  subroutine tomczak2014ZFOURGEDestructor(self)
    !% Destructor for the ``tomczak2014ZFOURGE'' survey geometry class.
    implicit none
    type(surveyGeometryTomczak2014ZFOURGE), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"/>
    return
  end subroutine tomczak2014ZFOURGEDestructor
  
  integer function tomczak2014ZFOURGEFieldCount(self)
    !% Return the number of fields in this sample.
    implicit none
    class(surveyGeometryTomczak2014ZFOURGE), intent(inout) :: self
    !GCC$ attributes unused :: self
    
    tomczak2014ZFOURGEFieldCount=tomczak2014ZFOURGEFields
    return
  end function tomczak2014ZFOURGEFieldCount

  double precision function tomczak2014ZFOURGEDistanceMinimum(self,mass,field)
    !% Compute the minimum distance at which a galaxy is included.
    implicit none
    class           (surveyGeometryTomczak2014ZFOURGE), intent(inout)           :: self
    double precision                                  , intent(in   )           :: mass
    integer                                           , intent(in   ), optional :: field
    !GCC$ attributes unused :: mass, field
    
    tomczak2014ZFOURGEDistanceMinimum=self%binDistanceMinimum
    return
  end function tomczak2014ZFOURGEDistanceMinimum

  double precision function tomczak2014ZFOURGEDistanceMaximum(self,mass,field)
    !% Compute the maximum distance at which a galaxy is visible.
    use Cosmology_Functions
    use Cosmology_Functions_Options
    use Galacticus_Error
    implicit none
    class           (surveyGeometryTomczak2014ZFOURGE), intent(inout)           :: self
    double precision                                  , intent(in   )           :: mass
    integer                                           , intent(in   ), optional :: field
    class           (cosmologyFunctionsClass         ), pointer                 :: cosmologyFunctions_
    double precision                                                            :: redshift           , logarithmicMass
    
    ! Validate field.
    if (.not.present(field)) call Galacticus_Error_Report('tomczak2014ZFOURGEDistanceMaximum','field must be specified')
    if (field < 1 .or. field > 2) call Galacticus_Error_Report('tomczak2014ZFOURGEDistanceMaximum','1 ≤ field ≤ 2 required')
    ! Find the limiting redshift for this mass. (See
    ! constraints/dataAnalysis/stellarMassFunctions_ZFOURGE_z0.2_2.5/massRedshiftRelation.pl for details.)
    logarithmicMass=log10(mass)
    select case (field)
    case (1)
       redshift=-114.659703302477d0+logarithmicMass*(45.9008293873578d0+logarithmicMass*(-6.16172903321726d0+logarithmicMass*(0.278223082977791d0)))
    case (2)
       redshift=-58.4827675323933d0+logarithmicMass*(20.2501133330335d0+logarithmicMass*(-2.35628286307041d0+logarithmicMass*(0.0927047000361006d0)))
    end select
    ! Get the default cosmology functions object.
    cosmologyFunctions_ => cosmologyFunctions()    
    ! Convert from redshift to comoving distance.
    tomczak2014ZFOURGEDistanceMaximum                                                 &
         &=cosmologyFunctions_%distanceComovingConvert(                               &
         &                                             output  =distanceTypeComoving, &
         &                                             redshift=redshift              &
         &                                            )
    ! Limit the maximum distance.
    tomczak2014ZFOURGEDistanceMaximum=min(tomczak2014ZFOURGEDistanceMaximum,self%binDistanceMaximum)
    return
  end function tomczak2014ZFOURGEDistanceMaximum

  double precision function tomczak2014ZFOURGEVolumeMaximum(self,mass,field)
    !% Compute the maximum volume within which a galaxy is visible.
    use Galacticus_Error
    implicit none
    class           (surveyGeometryTomczak2014ZFOURGE), intent(inout)           :: self
    double precision                                  , intent(in   )           :: mass
    integer                                           , intent(in   ), optional :: field

     ! Compute the volume.
    tomczak2014ZFOURGEVolumeMaximum                        &
         & =max(                                           &
         &       0.0d0                                   , &
         &       self%solidAngle(field)                    &
         &      *(                                         &
         &        +self%distanceMaximum   (mass,field)**3  &
         &        -self%binDistanceMinimum            **3  &
         &       )                                         &
         &      /3.0d0                                     &
         &     )
    return
  end function tomczak2014ZFOURGEVolumeMaximum

  function tomczak2014ZFOURGEMangleDirectory(self)
    !% Return the path to the directory containing \gls{mangle} files.
    implicit none
    class(surveyGeometryTomczak2014ZFOURGE), intent(inout) :: self
    type (varying_string                  )                :: tomczak2014ZFOURGEMangleDirectory
    !GCC$ attributes unused :: self
    
    tomczak2014ZFOURGEMangleDirectory=Galacticus_Input_Path()//"constraints/dataAnalysis/stellarMassFunctions_ZFOURGE_z0.2_2.5/"
    return
  end function tomczak2014ZFOURGEMangleDirectory
  
  subroutine tomczak2014ZFOURGEMangleFiles(self,mangleFiles)
    !% Return a list of \gls{mangle} files.
    implicit none
    class(surveyGeometryTomczak2014ZFOURGE)                           , intent(inout) :: self
    type (varying_string                  ), allocatable, dimension(:), intent(  out) :: mangleFiles
    
    mangleFiles=                                                       &
         &      [                                                      &
         &       "+"//self%mangleDirectory()//"/ZFOURGE-CDFS.ply:"  // &
         &       "+"//self%mangleDirectory()//"/ZFOURGE-COSMOS.ply:"// &
         &       "+"//self%mangleDirectory()//"/ZFOURGE-UDS.ply"    ,  &
         &       "+"//self%mangleDirectory()//"/NMBS-COSMOS.ply:"   // &
         &       "+"//self%mangleDirectory()//"/NMBS-AEGIS.ply"        &
         &      ]
    return
  end subroutine tomczak2014ZFOURGEMangleFiles

  integer function tomczak2014ZFOURGEAngularPowerMaximumDegree(self)
    !% Return the maximum degree for which angular power is computed for the \cite{tomczak_galaxy_2014} survey.
    implicit none
    class(surveyGeometryTomczak2014ZFOURGE), intent(inout) :: self
    !GCC$ attributes unused :: self
    
    tomczak2014ZFOURGEAngularPowerMaximumDegree=tomczak2014AngularPowerMaximumL
    return
  end function tomczak2014ZFOURGEAngularPowerMaximumDegree
  
