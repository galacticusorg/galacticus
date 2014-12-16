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

!% Implements the survey geometry used by \cite{caputi_stellar_2011}.
  
  !# <surveyGeometry name="surveyGeometryCaputi2011UKIDSSUDS">
  !#  <description>Implements the survey geometry of the SDSS sample used by \cite{caputi_stellar_2011}.</description>
  !# </surveyGeometry>

  type, extends(surveyGeometryRandomPoints) :: surveyGeometryCaputi2011UKIDSSUDS
     double precision :: binDistanceMinimum, binDistanceMaximum
   contains
     procedure :: distanceMinimum   => caputi2011UKIDSSUDSDistanceMinimum
     procedure :: distanceMaximum   => caputi2011UKIDSSUDSDistanceMaximum
     procedure :: volumeMaximum     => caputi2011UKIDSSUDSVolumeMaximum
     procedure :: solidAngle        => caputi2011UKIDSSUDSSolidAngle
     procedure :: randomsInitialize => caputi2011UKIDSSUDSRandomsInitialize
  end type surveyGeometryCaputi2011UKIDSSUDS

  interface surveyGeometryCaputi2011UKIDSSUDS
     !% Constructors for the \cite{caputi_stellar_2011} survey geometry class.
     module procedure caputi2011UKIDSSUDSDefaultConstructor
     module procedure caputi2011UKIDSSUDSConstructor
  end interface surveyGeometryCaputi2011UKIDSSUDS

  ! Default redshift bin to use.
  logical :: caputi2011UKIDSSUDSInitialized=.false.
  integer :: caputi2011UKIDSSUDSRedshiftBin

contains

  function caputi2011UKIDSSUDSDefaultConstructor()
    !% Default constructor for the \cite{caputi_stellar_2011} conditional mass function class.
    use Input_Parameters
    implicit none
    type(surveyGeometryCaputi2011UKIDSSUDS) :: caputi2011UKIDSSUDSDefaultConstructor

    if (.not.caputi2011UKIDSSUDSInitialized) then
       !$omp critical(caputi2011UKIDSSUDSInitialize)
       if (.not.caputi2011UKIDSSUDSInitialized) then
          ! Get the redshift bin to use.
          !@ <inputParameter>
          !@   <name>caputi2011UKIDSSUDSRedshiftBin</name>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The redshift bin (0, 1, or 2) of the \cite{caputi_stellar_2011} to use.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>starFormation</group>
          !@ </inputParameter>
          call Get_Input_Parameter('caputi2011UKIDSSUDSRedshiftBin',caputi2011UKIDSSUDSRedshiftBin)
          ! Record that the survey is initialized.
          caputi2011UKIDSSUDSInitialized=.true.
       end if
       !$omp end critical(caputi2011UKIDSSUDSInitialize)
    end if
    caputi2011UKIDSSUDSDefaultConstructor=caputi2011UKIDSSUDSConstructor(caputi2011UKIDSSUDSRedshiftBin)
    return
  end function caputi2011UKIDSSUDSDefaultConstructor

  function caputi2011UKIDSSUDSConstructor(redshiftBin)
    !% Default constructor for the \cite{caputi_stellar_2011} conditional mass function class.
    use Galacticus_Error
    use Input_Parameters
    use Cosmology_Functions
    use Cosmology_Functions_Options
    implicit none
    type            (surveyGeometryCaputi2011UKIDSSUDS)                :: caputi2011UKIDSSUDSConstructor
    integer                                            , intent(in   ) :: redshiftBin
    class           (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    double precision                                                   :: redshiftMinimum               , redshiftMaximum

    ! Find distance limits for this redshift bin.
    select case (redshiftBin)
    case(0)
       redshiftMinimum=3.00d0
       redshiftMaximum=3.50d0
    case(1)
       redshiftMinimum=3.50d0
       redshiftMaximum=4.25d0
    case(2)
       redshiftMinimum=4.25d0
       redshiftMaximum=5.00d0
    case default
       call Galacticus_Error_Report('caputi2011UKIDSSUDSConstructor','0≤redshiftBin≤2 is required')
    end select
    cosmologyFunctions_ => cosmologyFunctions()
    caputi2011UKIDSSUDSConstructor%binDistanceMinimum                                     &
         & =cosmologyFunctions_%distanceComovingConvert(                                  &
         &                                                 output  =distanceTypeComoving, &
         &                                                 redshift=redshiftMinimum       &
         &                                                )
    caputi2011UKIDSSUDSConstructor%binDistanceMaximum                                     &
         & =cosmologyFunctions_%distanceComovingConvert(                                  &
         &                                                 output  =distanceTypeComoving, &
         &                                                 redshift=redshiftMaximum       &
         &                                                )
    caputi2011UKIDSSUDSConstructor%geometryInitialized=.false.
    return
  end function caputi2011UKIDSSUDSConstructor
  
  double precision function caputi2011UKIDSSUDSDistanceMinimum(self,mass,field)
    !% Compute the minimum distance at which a galaxy is included.
    implicit none
    class           (surveyGeometryCaputi2011UKIDSSUDS), intent(inout)           :: self
    double precision                                   , intent(in   )           :: mass
    integer                                            , intent(in   ), optional :: field
    
    caputi2011UKIDSSUDSDistanceMinimum=self%binDistanceMinimum
    return
  end function caputi2011UKIDSSUDSDistanceMinimum

  double precision function caputi2011UKIDSSUDSDistanceMaximum(self,mass,field)
    !% Compute the maximum distance at which a galaxy is visible.
    use Cosmology_Functions
    use Cosmology_Functions_Options
    use Galacticus_Error
    implicit none
    class           (surveyGeometryCaputi2011UKIDSSUDS), intent(inout)           :: self
    double precision                                   , intent(in   )           :: mass
    integer                                            , intent(in   ), optional :: field
    class           (cosmologyFunctionsClass          ), pointer                 :: cosmologyFunctions_
    double precision                                                             :: redshift           , logarithmicMass
    
    ! Validate field.
    if (present(field).and.field /= 1) call Galacticus_Error_Report('caputi2011UKIDSSUDSDistanceMaximum','field = 1 required')
    ! Find the limiting redshift for this mass using a fit derived from Millennium Simulation SAMs. (See
    ! constraints/dataAnalysis/stellarMassFunctions_UKIDSS_UDS_z3_5/massLuminosityRelation.pl for details.)
    logarithmicMass=log10(mass)
    redshift=-56.247426278132d0+logarithmicMass*(5.88091022342758d0)
    ! Get the default cosmology functions object.
    cosmologyFunctions_ => cosmologyFunctions()    
    ! Convert from redshift to comoving distance.
    caputi2011UKIDSSUDSDistanceMaximum                                                &
         &=cosmologyFunctions_%distanceComovingConvert(                               &
         &                                             output  =distanceTypeComoving, &
         &                                             redshift=redshift              &
         &                                            )
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
    if (present(field).and.field /= 1) call Galacticus_Error_Report('caputi2011UKIDSSUDSVolumeMaximum','field = 1 required')
    ! Compute the volume.
    caputi2011UKIDSSUDSVolumeMaximum                       &
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
  end function caputi2011UKIDSSUDSVolumeMaximum

  double precision function caputi2011UKIDSSUDSSolidAngle(self,field)
    !% Return the solid angle of the \cite{caputi_stellar_2011} sample. Computed from survey mask (see {\tt
    !% constraints/dataAnalysis/stellarMassFunctions\_UKIDSS\_UDS\_z3\_5/surveyGeometryRandoms.pl}).
    use Galacticus_Error
    implicit none
    class           (surveyGeometryCaputi2011UKIDSSUDS), intent(inout)           :: self
    integer                                            , intent(in   ), optional :: field
    double precision                                   , parameter               :: solidAngleSurvey=1.59233703487973d-4
    
    ! Validate field.
    if (present(field).and.field /= 1) call Galacticus_Error_Report('caputi2011UKIDSSUDSSolidAngle','field = 1 required')
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
    use Galacticus_Input_Paths
    use System_Command
    use Galacticus_Error
    use IO_HDF5
    use File_Utilities
    implicit none
    class(surveyGeometryCaputi2011UKIDSSUDS), intent(inout) :: self
    type (hdf5Object                       )                :: surveyGeometryRandomsFile

    ! Generate the randoms file if necessary.
    if (.not.File_Exists(Galacticus_Input_Path()//&
         &"constraints/dataAnalysis/stellarMassFunctions_UKIDSS_UDS_z3_5/data/surveyGeometryRandoms.hdf5")) then
       call System_Command_Do(Galacticus_Input_Path()//"constraints/dataAnalysis/stellarMassFunctions_UKIDSS_UDS_z3_5/surveyGeometryRandoms.pl")
       if (.not.File_Exists(Galacticus_Input_Path()//"constraints/dataAnalysis/stellarMassFunctions_UKIDSS_UDS_z3_5/data/surveyGeometryRandoms.hdf5")) call Galacticus_Error_Report('caputi2011UKIDSSUDSWindowFunctions','unable to create survey geometry randoms file')
    end if
    ! Read the distribution of random points from file.
    !$omp critical(HDF5_Access)
    call surveyGeometryRandomsFile%openFile(char(Galacticus_Input_Path()//&
         &'constraints/dataAnalysis/stellarMassFunctions_UKIDSS_UDS_z3_5/data/surveyGeometryRandoms.hdf5')&
         &,readOnly=.true.)
    call surveyGeometryRandomsFile%readDataset('theta',self%randomTheta)
    call surveyGeometryRandomsFile%readDataset('phi'  ,self%randomPhi  )
    call surveyGeometryRandomsFile%close()
    !$omp end critical(HDF5_Access)    
    return
  end subroutine caputi2011UKIDSSUDSRandomsInitialize
