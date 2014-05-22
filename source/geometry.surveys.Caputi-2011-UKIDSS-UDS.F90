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

  type, extends(surveyGeometryClass) :: surveyGeometryCaputi2011UKIDSSUDS
     double precision :: binDistanceMinimum, binDistanceMaximum
   contains
     procedure :: windowFunctionAvailable => caputi2011UKIDSSUDSWindowFunctionAvailable
     procedure :: angularPowerAvailable   => caputi2011UKIDSSUDSAngularPowerAvailable
     procedure :: distanceMaximum         => caputi2011UKIDSSUDSDistanceMaximum
     procedure :: volumeMaximum           => caputi2011UKIDSSUDSVolumeMaximum
     procedure :: solidAngle              => caputi2011UKIDSSUDSSolidAngle
     procedure :: windowFunctions         => caputi2011UKIDSSUDSWindowFunctions
     procedure :: angularPower            => caputi2011UKIDSSUDSAngularPower
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
    return
  end function caputi2011UKIDSSUDSConstructor
  
  logical function caputi2011UKIDSSUDSWindowFunctionAvailable(self)
    !% Return true to indicate that survey window function is available.
    implicit none
    class(surveyGeometryCaputi2011UKIDSSUDS), intent(inout) :: self

    caputi2011UKIDSSUDSWindowFunctionAvailable=.true.
    return
  end function caputi2011UKIDSSUDSWindowFunctionAvailable

  logical function caputi2011UKIDSSUDSAngularPowerAvailable(self)
    !% Return false to indicate that survey angular power is not available.
    implicit none
    class(surveyGeometryCaputi2011UKIDSSUDS), intent(inout) :: self

    caputi2011UKIDSSUDSAngularPowerAvailable=.false.
    return
  end function caputi2011UKIDSSUDSAngularPowerAvailable

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

  subroutine caputi2011UKIDSSUDSWindowFunctions(self,mass1,mass2,gridCount,boxLength,windowFunction1,windowFunction2)
    !% Compute the window function for the survey.
    use FFTW3
    use Vectors
    use Pseudo_Random
    use FGSL
    use Meshes
    use, intrinsic :: ISO_C_Binding
    use Numerical_Constants_Math
    use Memory_Management
    use Galacticus_Display
    use ISO_Varying_String
    use Cosmology_Functions
    use String_Handling
    use Galacticus_Input_Paths
    use System_Command
    use Galacticus_Error
    use IO_HDF5
    use File_Utilities
    implicit none
    class           (surveyGeometryCaputi2011UKIDSSUDS), intent(inout) :: self
    double precision                                   , intent(in   )                                               :: mass1,mass2
    integer                                            , intent(in   )                                               :: gridCount
    double precision                                   , intent(  out)                                               :: boxLength
    complex         (c_double_complex                 ), intent(  out),     dimension(gridCount,gridCount,gridCount) :: windowFunction1,windowFunction2
    double precision                                   ,                    dimension(3                            ) :: origin,position1,position2
    integer                                                                                                          :: i,j
    double precision                                                                                                 :: comovingDistanceMaximum1, comovingDistanceMaximum2, &
         &                                                                                                              comovingDistanceMinimum1, comovingDistanceMinimum2
    type            (c_ptr                            )                                                              :: plan
    complex         (c_double_complex                 ),                    dimension(gridCount,gridCount,gridCount) :: selectionFunction1,selectionFunction2
    complex         (c_double_complex                 )                                                              :: normalization
    logical                                            , save                                                        :: geometryInitialized=.false.
    double precision                                                                                                 :: rightAscension,declination,distance1,distance2
    double precision                                   , save, allocatable, dimension(:                            ) :: randomTheta,randomPhi
    integer                                            , save                                                        :: randomsCount
    type            (fgsl_rng                         ), save                                                        :: pseudoSequenceObject
    logical                                            , save                                                        :: reset=.true.
    type            (varying_string                   )                                                              :: message
    double precision                                   , save                                                        :: surveyDistanceMinimum, surveyDistanceMaximum
    type          (hdf5Object                         )                                                              :: surveyGeometryRandomsFile
    class         (cosmologyFunctionsClass            ), pointer                                                     :: cosmologyFunctions_

    ! Initialize geometry if necessary.
    if (.not.geometryInitialized) then
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
       call surveyGeometryRandomsFile%readDataset('theta',randomTheta)
       call surveyGeometryRandomsFile%readDataset('phi'  ,randomPhi  )
       call surveyGeometryRandomsFile%close()
       randomsCount=size(randomTheta)
       !$omp end critical(HDF5_Access)
       
       ! Get the default cosmology functions object.
       cosmologyFunctions_ => cosmologyFunctions()
       ! Compute the distances corresponding to the minimum and maximum redshifts.
       surveyDistanceMinimum=cosmologyFunctions_%distanceComoving(cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(0.001d0)))
       surveyDistanceMaximum=cosmologyFunctions_%distanceComoving(cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(0.500d0)))
       geometryInitialized=.true.
    end if

    ! Find the comoving distance corresponding to this distance module.
    comovingDistanceMaximum1=min(self%distanceMaximum(mass1),surveyDistanceMaximum)
    comovingDistanceMaximum2=min(self%distanceMaximum(mass2),surveyDistanceMaximum)
    comovingDistanceMinimum1=surveyDistanceMinimum
    comovingDistanceMinimum2=surveyDistanceMinimum

    ! Determine a suitable box length for the window function calculation.
    boxLength=3.0d0*max(comovingDistanceMaximum1,comovingDistanceMaximum2)

    ! Set up origin.
    origin=([0.5d0,0.5d0,0.0d0]/dble(gridCount)+[0.5d0,0.5d0,0.5d0])*dble(boxLength)

    ! Populate the cube with the survey selection function.
    selectionFunction1=cmplx(0.0d0,0.0d0)
    selectionFunction2=cmplx(0.0d0,0.0d0)

    ! Loop over randoms.
    do i=1,randomsCount
       ! Choose random distances.
       distance1=(Pseudo_Random_Get(pseudoSequenceObject,reset)*(comovingDistanceMaximum1**3-comovingDistanceMinimum1**3)+comovingDistanceMinimum1**3)**(1.0d0/3.0d0)
       distance2=(Pseudo_Random_Get(pseudoSequenceObject,reset)*(comovingDistanceMaximum2**3-comovingDistanceMinimum2**3)+comovingDistanceMinimum2**3)**(1.0d0/3.0d0)
       ! Convert to Cartesian coordinates.
       position1=distance1*[sin(randomTheta(i))*cos(randomPhi(i)),sin(randomTheta(i))*sin(randomPhi(i)),cos(randomTheta(i))]+origin
       position2=distance2*[sin(randomTheta(i))*cos(randomPhi(i)),sin(randomTheta(i))*sin(randomPhi(i)),cos(randomTheta(i))]+origin
       ! Apply to the mesh.
       call Meshes_Apply_Point(selectionFunction1,boxLength,position1,pointWeight=cmplx(1.0d0,0.0d0,kind=c_double_complex),cloudType=cloudTypeTriangular)
       call Meshes_Apply_Point(selectionFunction2,boxLength,position2,pointWeight=cmplx(1.0d0,0.0d0,kind=c_double_complex),cloudType=cloudTypeTriangular)
    end do
    ! Take the Fourier transform of the selection function.
    plan=fftw_plan_dft_3d(gridCount,gridCount,gridCount,selectionFunction1,windowFunction1,FFTW_FORWARD,FFTW_ESTIMATE)
    call fftw_execute_dft(plan,selectionFunction1,windowFunction1)
    call fftw_execute_dft(plan,selectionFunction2,windowFunction2)
    call fftw_destroy_plan(plan)

    ! Normalize the window function.
    normalization=windowFunction1(1,1,1)
    if (real(normalization) > 0.0d0) windowFunction1=windowFunction1/normalization
    normalization=windowFunction2(1,1,1)
    if (real(normalization) > 0.0d0) windowFunction2=windowFunction2/normalization

    return
  end subroutine caputi2011UKIDSSUDSWindowFunctions

  double precision function caputi2011UKIDSSUDSAngularPower(self,i,j,l)
    !% Angular power is not available, so simply aborts.
    use Galacticus_Error
    implicit none
    class  (surveyGeometryCaputi2011UKIDSSUDS), intent(inout) :: self
    integer                                   , intent(in   ) :: i   , j, l

    call Galacticus_Error_Report('caputi2011UKIDSSUDSAngularPower','angular power is not available')
    return
  end function caputi2011UKIDSSUDSAngularPower
