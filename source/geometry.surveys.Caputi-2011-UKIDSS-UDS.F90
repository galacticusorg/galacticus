!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Implements the survey geometry used by \cite{caputi_stellar_2011}.
  
  use Cosmology_Functions

  !# <surveyGeometry name="surveyGeometryCaputi2011UKIDSSUDS">
  !#  <description>Implements the survey geometry of the SDSS sample used by \cite{caputi_stellar_2011}.</description>
  !# </surveyGeometry>
  type, extends(surveyGeometryRandomPoints) :: surveyGeometryCaputi2011UKIDSSUDS
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     double precision                                   :: binDistanceMinimum           , binDistanceMaximum, &
          &                                                redshiftMinimum              , redshiftMaximum
   contains
     final     ::                      caputi2011UKIDSSUDSDestructor
     procedure :: distanceMinimum   => caputi2011UKIDSSUDSDistanceMinimum
     procedure :: distanceMaximum   => caputi2011UKIDSSUDSDistanceMaximum
     procedure :: volumeMaximum     => caputi2011UKIDSSUDSVolumeMaximum
     procedure :: solidAngle        => caputi2011UKIDSSUDSSolidAngle
     procedure :: randomsInitialize => caputi2011UKIDSSUDSRandomsInitialize
  end type surveyGeometryCaputi2011UKIDSSUDS

  interface surveyGeometryCaputi2011UKIDSSUDS
     !% Constructors for the \cite{caputi_stellar_2011} survey geometry class.
     module procedure caputi2011UKIDSSUDSConstructorParameters
     module procedure caputi2011UKIDSSUDSConstructorInternal
  end interface surveyGeometryCaputi2011UKIDSSUDS

contains

  function caputi2011UKIDSSUDSConstructorParameters(parameters) result(self)
    !% Default constructor for the \cite{caputi_stellar_2011} conditional mass function class.
    use Input_Parameters
    implicit none
    type   (surveyGeometryCaputi2011UKIDSSUDS)                :: self
    type   (inputParameters                  ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    integer                                                   :: redshiftBin
    
    ! Check and read parameters.
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !# <inputParameter>
    !#   <name>redshiftBin</name>
    !#   <source>parameters</source>
    !#   <description>The redshift bin (0, 1, or 2) of the \cite{caputi_stellar_2011} to use.</description>
    !#   <type>integer</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    self=surveyGeometryCaputi2011UKIDSSUDS(redshiftBin,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"/>
    return
  end function caputi2011UKIDSSUDSConstructorParameters

  function caputi2011UKIDSSUDSConstructorInternal(redshiftBin,cosmologyFunctions_) result(self)
    !% Internal constructor for the \cite{caputi_stellar_2011} conditional mass function class.
    use Galacticus_Error
    use Cosmology_Functions_Options
    implicit none
    type            (surveyGeometryCaputi2011UKIDSSUDS)                        :: self
    integer                                            , intent(in   )         :: redshiftBin
    class           (cosmologyFunctionsClass          ), intent(in   ), target :: cosmologyFunctions_
    !# <constructorAssign variables="*cosmologyFunctions_"/>

    ! Find distance limits for this redshift bin.
    select case (redshiftBin)
    case(0)
       self%redshiftMinimum=3.00d0
       self%redshiftMaximum=3.50d0
    case(1)
       self%redshiftMinimum=3.50d0
       self%redshiftMaximum=4.25d0
    case(2)
       self%redshiftMinimum=4.25d0
       self%redshiftMaximum=5.00d0
    case default
       call Galacticus_Error_Report('0≤redshiftBin≤2 is required'//{introspection:location})
    end select
    self%binDistanceMinimum                                                                 &
         & =self%cosmologyFunctions_%distanceComovingConvert(                               &
         &                                                   output  =distanceTypeComoving, &
         &                                                   redshift=self%redshiftMinimum  &
         &                                                  )
    self%binDistanceMaximum                                                                 &
         & =self%cosmologyFunctions_%distanceComovingConvert(                               &
         &                                                   output  =distanceTypeComoving, &
         &                                                   redshift=self%redshiftMaximum  &
         &                                                  )
    self%geometryInitialized=.false.
    return
  end function caputi2011UKIDSSUDSConstructorInternal
  
  subroutine caputi2011UKIDSSUDSDestructor(self)
    !% Destructor for the ``caputi2011UKIDSSUDS'' survey geometry class.
    implicit none
    type(surveyGeometryCaputi2011UKIDSSUDS), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"/>
    return
  end subroutine caputi2011UKIDSSUDSDestructor
  
  double precision function caputi2011UKIDSSUDSDistanceMinimum(self,mass,magnitudeAbsolute,luminosity,field)
    !% Compute the minimum distance at which a galaxy is included.
    implicit none
    class           (surveyGeometryCaputi2011UKIDSSUDS), intent(inout)           :: self
    double precision                                   , intent(in   ), optional :: mass , magnitudeAbsolute, luminosity
    integer                                            , intent(in   ), optional :: field
    !GCC$ attributes unused :: mass, field, magnitudeAbsolute, luminosity
    
    caputi2011UKIDSSUDSDistanceMinimum=self%binDistanceMinimum
    return
  end function caputi2011UKIDSSUDSDistanceMinimum

  double precision function caputi2011UKIDSSUDSDistanceMaximum(self,mass,magnitudeAbsolute,luminosity,field)
    !% Compute the maximum distance at which a galaxy is visible.
    use Cosmology_Functions_Options
    use Galacticus_Error
    implicit none
    class           (surveyGeometryCaputi2011UKIDSSUDS), intent(inout)           :: self
    double precision                                   , intent(in   ), optional :: mass    , magnitudeAbsolute, luminosity
    integer                                            , intent(in   ), optional :: field
    double precision                                                             :: redshift, logarithmicMass
    !GCC$ attributes unused :: magnitudeAbsolute, luminosity

    ! Validate field.
    if (present(field).and.field /= 1) call Galacticus_Error_Report('field = 1 required'//{introspection:location})
    ! Find the limiting redshift for this mass using a fit derived from Millennium Simulation SAMs. (See
    ! constraints/dataAnalysis/stellarMassFunctions_UKIDSS_UDS_z3_5/massLuminosityRelation.pl for details.)
    logarithmicMass=log10(mass)
    redshift=-56.247426278132d0+logarithmicMass*(5.88091022342758d0)
    ! Convert from redshift to comoving distance.
    caputi2011UKIDSSUDSDistanceMaximum                                                                                             &
         &=self%cosmologyFunctions_%distanceComovingConvert(                                                                       &
         &                                                  output  =distanceTypeComoving                                        , &
         &                                                  redshift=min(max(redshift,self%redshiftMinimum),self%redshiftMaximum)  &
         &                                                 )
    ! Limit the maximum distance.
    caputi2011UKIDSSUDSDistanceMaximum=min(caputi2011UKIDSSUDSDistanceMaximum,self%binDistanceMaximum)
    return
  end function caputi2011UKIDSSUDSDistanceMaximum

  double precision function caputi2011UKIDSSUDSVolumeMaximum(self,mass,field)
    !% Compute the maximum volume within which a galaxy is visible.
    use Galacticus_Error
    implicit none
    class           (surveyGeometryCaputi2011UKIDSSUDS), intent(inout)           :: self
    double precision                                   , intent(in   )           :: mass
    integer                                            , intent(in   ), optional :: field
    
    ! Validate field.
    if (present(field).and.field /= 1) call Galacticus_Error_Report('field = 1 required'//{introspection:location})
    ! Compute the volume.
    caputi2011UKIDSSUDSVolumeMaximum                             &
         & =max(                                                 &
         &       0.0d0                                         , &
         &       self%solidAngle(field)                          &
         &      *(                                               &
         &        +self%distanceMaximum   (mass,field=field)**3  &
         &        -self%binDistanceMinimum                  **3  &
         &       )                                               &
         &      /3.0d0                                           &
         &     )
    return
  end function caputi2011UKIDSSUDSVolumeMaximum

  double precision function caputi2011UKIDSSUDSSolidAngle(self,field)
    !% Return the solid angle of the \cite{caputi_stellar_2011} sample. Computed from survey mask (see {\normalfont \ttfamily
    !% constraints/dataAnalysis/stellarMassFunctions\_UKIDSS\_UDS\_z3\_5/surveyGeometryRandoms.pl}).
    use Galacticus_Error
    implicit none
    class           (surveyGeometryCaputi2011UKIDSSUDS), intent(inout)           :: self
    integer                                            , intent(in   ), optional :: field
    double precision                                   , parameter               :: solidAngleSurvey=1.59233703487973d-4
    !GCC$ attributes unused :: self
    
    ! Validate field.
    if (present(field).and.field /= 1) call Galacticus_Error_Report('field = 1 required'//{introspection:location})
    caputi2011UKIDSSUDSSolidAngle=solidAngleSurvey
    return
  end function caputi2011UKIDSSUDSSolidAngle

  subroutine caputi2011UKIDSSUDSRandomsInitialize(self)
    !% Load random points for the survey.
    use Numerical_Constants_Math
    use Memory_Management
    use Galacticus_Display
    use ISO_Varying_String
    use String_Handling
    use Galacticus_Paths
    use System_Command
    use Galacticus_Error
    use IO_HDF5
    use File_Utilities
    implicit none
    class(surveyGeometryCaputi2011UKIDSSUDS), intent(inout) :: self
    type (hdf5Object                       )                :: surveyGeometryRandomsFile

    ! Generate the randoms file if necessary.
    if (.not.File_Exists(galacticusPath(pathTypeExec)//&
         &"constraints/dataAnalysis/stellarMassFunctions_UKIDSS_UDS_z3_5/data/surveyGeometryRandoms.hdf5")) then
       call System_Command_Do(galacticusPath(pathTypeExec)//"constraints/dataAnalysis/stellarMassFunctions_UKIDSS_UDS_z3_5/surveyGeometryRandoms.pl")
       if (.not.File_Exists(galacticusPath(pathTypeExec)//"constraints/dataAnalysis/stellarMassFunctions_UKIDSS_UDS_z3_5/data/surveyGeometryRandoms.hdf5")) call Galacticus_Error_Report('unable to create survey geometry randoms file'//{introspection:location})
    end if
    ! Read the distribution of random points from file.
    !$ call hdf5Access%set()
    call surveyGeometryRandomsFile%openFile(char(galacticusPath(pathTypeExec)//&
         &'constraints/dataAnalysis/stellarMassFunctions_UKIDSS_UDS_z3_5/data/surveyGeometryRandoms.hdf5')&
         &,readOnly=.true.)
    call surveyGeometryRandomsFile%readDataset('theta',self%randomTheta)
    call surveyGeometryRandomsFile%readDataset('phi'  ,self%randomPhi  )
    call surveyGeometryRandomsFile%close()
    !$ call hdf5Access%unset()    
    return
  end subroutine caputi2011UKIDSSUDSRandomsInitialize
