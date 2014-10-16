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

!% Implements the geometry of the ULTRAVISTA survey used by \cite{muzzin_evolution_2013}.
  
  !# <surveyGeometry name="surveyGeometryMuzzin2013ULTRAVISTA">
  !#  <description>Implements the geometry of the ULTRAVISTA survey of \cite{muzzin_evolution_2013}.</description>
  !# </surveyGeometry>

  use Galacticus_Input_Paths

  type, extends(surveyGeometryMangle) :: surveyGeometryMuzzin2013ULTRAVISTA
     integer          :: redshiftBin
     double precision :: binDistanceMinimum, binDistanceMaximum
   contains
     procedure :: fieldCount                => muzzin2013ULTRAVISTAFieldCount
     procedure :: distanceMinimum           => muzzin2013ULTRAVISTADistanceMinimum
     procedure :: distanceMaximum           => muzzin2013ULTRAVISTADistanceMaximum
     procedure :: volumeMaximum             => muzzin2013ULTRAVISTAVolumeMaximum
     procedure :: angularPowerMaximumDegree => muzzin2013ULTRAVISTAAngularPowerMaximumDegree
     procedure :: mangleDirectory           => muzzin2013ULTRAVISTAMangleDirectory
     procedure :: mangleFiles               => muzzin2013ULTRAVISTAMangleFiles
   end type surveyGeometryMuzzin2013ULTRAVISTA

  interface surveyGeometryMuzzin2013ULTRAVISTA
     !% Constructors for the \cite{muzzin_evolution_2013} survey geometry class.
     module procedure muzzin2013ULTRAVISTAConstructor
     module procedure muzzin2013ULTRAVISTADefaultConstructor
  end interface surveyGeometryMuzzin2013ULTRAVISTA

  ! Paths and file names for mangle polygon files.
  integer, parameter :: muzzin2013ULTRAVISTAFields        =1
  logical            :: muzzin2013ULTRAVISTABinInitialized=.false.
  integer            :: muzzin2013ULTRAVISTARedshiftBin

  ! Angular power spectra.
  integer, parameter :: muzzin2013AngularPowerMaximumL=3600

contains

  function muzzin2013ULTRAVISTADefaultConstructor()
    !% Default constructor for the \cite{muzzin_evolution_2013} conditional mass function class.
    use Input_Parameters
    implicit none
    type(surveyGeometryMuzzin2013ULTRAVISTA) :: muzzin2013ULTRAVISTADefaultConstructor

    if (.not.muzzin2013ULTRAVISTABinInitialized) then
       !$omp critical(muzzin2013ULTRAVISTABinInitialize)
       if (.not.muzzin2013ULTRAVISTABinInitialized) then
          ! Get the redshift bin to use.
          !@ <inputParameter>
          !@   <name>muzzin2013ULTRAVISTARedshiftBin</name>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The redshift bin (0, 1, 2, 3, 4, 5, 6) of the \cite{muzzin_evolution_2013} mass function to use.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('muzzin2013ULTRAVISTARedshiftBin',muzzin2013ULTRAVISTARedshiftBin)
          muzzin2013ULTRAVISTABinInitialized=.true.
       end if
       !$omp end critical(muzzin2013ULTRAVISTABinInitialize)
    end if
    muzzin2013ULTRAVISTADefaultConstructor=muzzin2013ULTRAVISTAConstructor(muzzin2013ULTRAVISTARedshiftBin)
   return
  end function muzzin2013ULTRAVISTADefaultConstructor

  function muzzin2013ULTRAVISTAConstructor(redshiftBin)
    !% Generic constructor for the \cite{muzzin_evolution_2013} mass function class.
    use Galacticus_Error
    use Input_Parameters
    use Cosmology_Functions
    use Cosmology_Functions_Options
    implicit none
    type            (surveyGeometryMuzzin2013ULTRAVISTA)                :: muzzin2013ULTRAVISTAConstructor
    integer                                             , intent(in   ) :: redshiftBin
    class           (cosmologyFunctionsClass           ), pointer       :: cosmologyFunctions_
    double precision                                                    :: redshiftMinimum               , redshiftMaximum

    ! Find distance limits for this redshift bin.
    muzzin2013ULTRAVISTAConstructor%redshiftBin=redshiftBin
    select case (redshiftBin)
    case(0)
       redshiftMinimum=0.20d0
       redshiftMaximum=0.50d0
    case(1)
       redshiftMinimum=0.50d0
       redshiftMaximum=1.00d0
    case(2)
       redshiftMinimum=1.00d0
       redshiftMaximum=1.50d0
    case(3)
       redshiftMinimum=1.50d0
       redshiftMaximum=2.00d0
    case(4)
       redshiftMinimum=2.00d0
       redshiftMaximum=2.50d0
    case(5)
       redshiftMinimum=2.50d0
       redshiftMaximum=3.00d0
    case(6)
       redshiftMinimum=3.00d0
       redshiftMaximum=4.00d0
 
    case default
       call Galacticus_Error_Report('muzzin2013ULTRAVISTAConstructor','0≤redshiftBin≤6 is required')
    end select
    cosmologyFunctions_ => cosmologyFunctions()
    muzzin2013ULTRAVISTAConstructor%binDistanceMinimum                                 &
         & =cosmologyFunctions_%distanceComovingConvert(                               &
         &                                              output  =distanceTypeComoving, &
         &                                              redshift=redshiftMinimum       &
         &                                             )
    muzzin2013ULTRAVISTAConstructor%binDistanceMaximum                                 &
         & =cosmologyFunctions_%distanceComovingConvert(                               &
         &                                              output  =distanceTypeComoving, &
         &                                              redshift=redshiftMaximum       &
         &                                             )
    muzzin2013ULTRAVISTAConstructor%solidAnglesInitialized =.false.
    muzzin2013ULTRAVISTAConstructor%angularPowerInitialized=.false.
    muzzin2013ULTRAVISTAConstructor%windowInitialized      =.false.
    return
  end function muzzin2013ULTRAVISTAConstructor
  
  integer function muzzin2013ULTRAVISTAFieldCount(self)
    !% Return the number of fields in this sample.
    implicit none
    class(surveyGeometryMuzzin2013ULTRAVISTA), intent(inout) :: self

    muzzin2013ULTRAVISTAFieldCount=muzzin2013ULTRAVISTAFields
    return
  end function muzzin2013ULTRAVISTAFieldCount

  double precision function muzzin2013ULTRAVISTADistanceMinimum(self,mass,field)
    !% Compute the minimum distance at which a galaxy is included.
    implicit none
    class           (surveyGeometryMuzzin2013ULTRAVISTA), intent(inout)           :: self
    double precision                                    , intent(in   )           :: mass
    integer                                             , intent(in   ), optional :: field

    muzzin2013ULTRAVISTADistanceMinimum=self%binDistanceMinimum
    return
  end function muzzin2013ULTRAVISTADistanceMinimum

  double precision function muzzin2013ULTRAVISTADistanceMaximum(self,mass,field)
    !% Compute the maximum distance at which a galaxy is visible.
    use Cosmology_Functions
    use Cosmology_Functions_Options
    use Galacticus_Error
    implicit none
    class           (surveyGeometryMuzzin2013ULTRAVISTA), intent(inout)           :: self
    double precision                                    , intent(in   )           :: mass
    integer                                             , intent(in   ), optional :: field
    class           (cosmologyFunctionsClass           ), pointer                 :: cosmologyFunctions_
    double precision                                                              :: redshift           , logarithmicMass
    
    ! Find the limiting redshift for this mass. (See
    ! constraints/dataAnalysis/stellarMassFunctions_ULTRAVISTA_z0.2_4.0/massRedshiftRelation.pl for details.)
    logarithmicMass=log10(mass)
    redshift=-6076.22869161192d0+logarithmicMass*(3231.43806947672d0+logarithmicMass*(-686.81594922437d0+logarithmicMass*(72.9147764759627d0+logarithmicMass*(-3.86638121388152d0+logarithmicMass*(0.0819398410572916d0)))))
    if (logarithmicMass < 11.24d0) then
       redshift=min(4.0d0,redshift/(1.0d0-exp((logarithmicMass-11.24d0)/0.02d0)))
    else
       redshift=    4.0d0
    end if
    ! Get the default cosmology functions object.
    cosmologyFunctions_ => cosmologyFunctions()    
    ! Convert from redshift to comoving distance.
    muzzin2013ULTRAVISTADistanceMaximum                                               &
         &=cosmologyFunctions_%distanceComovingConvert(                               &
         &                                             output  =distanceTypeComoving, &
         &                                             redshift=redshift              &
         &                                            )
    ! Limit the maximum distance.
    muzzin2013ULTRAVISTADistanceMaximum=min(muzzin2013ULTRAVISTADistanceMaximum,self%binDistanceMaximum)
    return
  end function muzzin2013ULTRAVISTADistanceMaximum

  double precision function muzzin2013ULTRAVISTAVolumeMaximum(self,mass,field)
    !% Compute the maximum volume within which a galaxy is visible.
    use Galacticus_Error
    implicit none
    class           (surveyGeometryMuzzin2013ULTRAVISTA), intent(inout)           :: self
    double precision                                    , intent(in   )           :: mass
    integer                                             , intent(in   ), optional :: field

     ! Compute the volume.
    muzzin2013ULTRAVISTAVolumeMaximum                      &
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
  end function muzzin2013ULTRAVISTAVolumeMaximum

  function muzzin2013ULTRAVISTAMangleDirectory(self)
    !% Return the path to the directory containing \gls{mangle} files.
    implicit none
    class(surveyGeometryMuzzin2013ULTRAVISTA), intent(inout) :: self
    type (varying_string                    )                :: muzzin2013ULTRAVISTAMangleDirectory

    muzzin2013ULTRAVISTAMangleDirectory=Galacticus_Input_Path()//"constraints/dataAnalysis/stellarMassFunctions_ULTRAVISTA_z0.2_4.0/"
    return
  end function muzzin2013ULTRAVISTAMangleDirectory
  
  subroutine muzzin2013ULTRAVISTAMangleFiles(self,mangleFiles)
    !% Return a list of \gls{mangle} files.
    implicit none
    class(surveyGeometryMuzzin2013ULTRAVISTA)                           , intent(inout) :: self
    type (varying_string                    ), allocatable, dimension(:), intent(  out) :: mangleFiles
    
   mangleFiles=                  &
        &      [                 &
        &       'surveyMask.ply' &
        &      ]
    return
  end subroutine muzzin2013ULTRAVISTAMangleFiles

  integer function muzzin2013ULTRAVISTAAngularPowerMaximumDegree(self)
    !% Return the maximum degree for which angular power is computed for the \cite{muzzin_evolution_2013} survey.
    implicit none
    class(surveyGeometryMuzzin2013ULTRAVISTA), intent(inout) :: self

    muzzin2013ULTRAVISTAAngularPowerMaximumDegree=muzzin2013AngularPowerMaximumL
    return
  end function muzzin2013ULTRAVISTAAngularPowerMaximumDegree
  
