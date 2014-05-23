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

!% Implements the geometry of the PRIMUS survey used by \cite{moustakas_primus:_2013}.
  
  !# <surveyGeometry name="surveyGeometryMoustakas2013PRIMUS">
  !#  <description>Implements the geometry of the PRIMUS survey of \cite{moustakas_primus:_2013}.</description>
  !# </surveyGeometry>

  use Galacticus_Input_Paths

  type, extends(surveyGeometryClass) :: surveyGeometryMoustakas2013PRIMUS
     double precision :: binDistanceMinimum, binDistanceMaximum
   contains
     procedure :: windowFunctionAvailable   => moustakas2013PRIMUSWindowFunctionAvailable
     procedure :: angularPowerAvailable     => moustakas2013PRIMUSAngularPowerAvailable
     procedure :: fieldCount                => moustakas2013PRIMUSFieldCount
     procedure :: distanceMinimum           => moustakas2013PRIMUSDistanceMinimum
     procedure :: distanceMaximum           => moustakas2013PRIMUSDistanceMaximum
     procedure :: volumeMaximum             => moustakas2013PRIMUSVolumeMaximum
     procedure :: solidAngle                => moustakas2013PRIMUSSolidAngle
     procedure :: windowFunctions           => moustakas2013PRIMUSWindowFunctions
     procedure :: angularPower              => moustakas2013PRIMUSAngularPower
     procedure :: angularPowerMaximumDegree => moustakas2013PRIMUSAngularPowerMaximumDegree
  end type surveyGeometryMoustakas2013PRIMUS

  interface surveyGeometryMoustakas2013PRIMUS
     !% Constructors for the \cite{moustakas_primus:_2013} survey geometry class.
     module procedure moustakas2013PRIMUSConstructor
     module procedure moustakas2013PRIMUSDefaultConstructor
  end interface surveyGeometryMoustakas2013PRIMUS

  ! Paths and file names for mangle polygon files.
  integer                         , parameter                                                                     :: moustakas2013PRIMUSFields                =5
  integer                         , parameter                                                                     :: moustakas2013PRIMUSFieldPairs            =  moustakas2013PRIMUSFields     &
       &                                                                                                                                                       *(moustakas2013PRIMUSFields+1)  &
       &                                                                                                                                                       /2
  logical                                                                                                         :: moustakas2013PRIMUSInitialized           =.false.                       , &
       &                                                                                                             moustakas2013PRIMUSBinInitialized        =.false.
  integer                                                                                                         :: moustakas2013PRIMUSRedshiftBin
  type            (varying_string)                                                                                :: moustakas2013PRIMUSPath
  type            (varying_string), dimension(                                     moustakas2013PRIMUSFields    ) :: moustakas2013PRIMUSMangleFiles

  ! Solid angles.
  logical                                                                                                         :: moustakas2013PRIMUSSolidAnglesInitialized=.false.
  double precision                , dimension(                                     moustakas2013PRIMUSFields    ) :: moustakas2013PRIMUSSolidAngles

  ! Angular power spectra.
  integer                         , parameter                                                                     :: moustaskas2013AngularPowerMaximumL        =3000
  logical                                                                                                         :: moustakas2013PRIMUSAngularPowerInitialized=.false.
  double precision                , dimension(moustaskas2013AngularPowerMaximumL+1,moustakas2013PRIMUSFieldPairs) :: moustakas2013PRIMUSAngularPowerSpectra

contains

  function moustakas2013PRIMUSDefaultConstructor()
    !% Default constructor for the \cite{moustakas_primus:_2013} conditional mass function class.
    use Input_Parameters
    implicit none
    type(surveyGeometryMoustakas2013PRIMUS) :: moustakas2013PRIMUSDefaultConstructor

    call moustakas2013PRIMUSInitialize()
    if (.not.moustakas2013PRIMUSBinInitialized) then
       !$omp critical(moustakas2013PRIMUSBinInitialize)
       if (.not.moustakas2013PRIMUSBinInitialized) then
          ! Get the redshift bin to use.
          !@ <inputParameter>
          !@   <name>moustakas2013PRIMUSRedshiftBin</name>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The redshift bin (0, 1, or 2) of the \cite{caputi_stellar_2011} to use.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>starFormation</group>
          !@ </inputParameter>
          call Get_Input_Parameter('moustakas2013PRIMUSRedshiftBin',moustakas2013PRIMUSRedshiftBin)
          moustakas2013PRIMUSBinInitialized=.true.
       end if
       !$omp end critical(moustakas2013PRIMUSBinInitialize)
    end if
    moustakas2013PRIMUSDefaultConstructor=moustakas2013PRIMUSConstructor(moustakas2013PRIMUSRedshiftBin)
   return
  end function moustakas2013PRIMUSDefaultConstructor

  function moustakas2013PRIMUSConstructor(redshiftBin)
    !% Generic constructor for the \cite{moustakas_primus:_2013} mass function class.
    use Galacticus_Error
    use Input_Parameters
    use Cosmology_Functions
    use Cosmology_Functions_Options
    implicit none
    type            (surveyGeometryMoustakas2013PRIMUS)                :: moustakas2013PRIMUSConstructor
    integer                                            , intent(in   ) :: redshiftBin
    class           (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    double precision                                                   :: redshiftMinimum               , redshiftMaximum

    call moustakas2013PRIMUSInitialize()
    ! Find distance limits for this redshift bin.
    select case (redshiftBin)
    case(0)
       redshiftMinimum=0.20d0
       redshiftMaximum=0.30d0
    case(1)
       redshiftMinimum=0.30d0
       redshiftMaximum=0.40d0
    case(2)
       redshiftMinimum=0.40d0
       redshiftMaximum=0.50d0
    case(3)
       redshiftMinimum=0.50d0
       redshiftMaximum=0.65d0
    case(4)
       redshiftMinimum=0.65d0
       redshiftMaximum=0.80d0
    case(5)
       redshiftMinimum=0.80d0
       redshiftMaximum=1.00d0
    case default
       call Galacticus_Error_Report('moustakas2013PRIMUSConstructor','0≤redshiftBin≤5 is required')
    end select
    cosmologyFunctions_ => cosmologyFunctions()
    moustakas2013PRIMUSConstructor%binDistanceMinimum                                     &
         & =cosmologyFunctions_%distanceComovingConvert(                                  &
         &                                                 output  =distanceTypeComoving, &
         &                                                 redshift=redshiftMinimum       &
         &                                                )
    moustakas2013PRIMUSConstructor%binDistanceMaximum                                     &
         & =cosmologyFunctions_%distanceComovingConvert(                                  &
         &                                                 output  =distanceTypeComoving, &
         &                                                 redshift=redshiftMaximum       &
         &                                                )
    return
  end function moustakas2013PRIMUSConstructor
  
  integer function moustakas2013PRIMUSFieldCount(self)
    !% Return the number of fields in this sample.
    implicit none
    class(surveyGeometryMoustakas2013PRIMUS), intent(inout) :: self

    moustakas2013PRIMUSFieldCount=moustakas2013PRIMUSFields
    return
  end function moustakas2013PRIMUSFieldCount
    
  logical function moustakas2013PRIMUSWindowFunctionAvailable(self)
    !% Return false to indicate that survey window function is available.
    implicit none
    class(surveyGeometryMoustakas2013PRIMUS), intent(inout) :: self

    moustakas2013PRIMUSWindowFunctionAvailable=.false.
    return
  end function moustakas2013PRIMUSWindowFunctionAvailable

  logical function moustakas2013PRIMUSAngularPowerAvailable(self)
    !% Return true to indicate that survey angular power is not available.
    implicit none
    class(surveyGeometryMoustakas2013PRIMUS), intent(inout) :: self

    moustakas2013PRIMUSAngularPowerAvailable=.true.
    return
  end function moustakas2013PRIMUSAngularPowerAvailable

  double precision function moustakas2013PRIMUSDistanceMinimum(self,mass,field)
    !% Compute the minimum distance at which a galaxy is included.
    implicit none
    class           (surveyGeometryMoustakas2013PRIMUS), intent(inout)           :: self
    double precision                                   , intent(in   )           :: mass
    integer                                            , intent(in   ), optional :: field

    moustakas2013PRIMUSDistanceMinimum=self%binDistanceMinimum
    return
  end function moustakas2013PRIMUSDistanceMinimum

  double precision function moustakas2013PRIMUSDistanceMaximum(self,mass,field)
    !% Compute the maximum distance at which a galaxy is visible.
    use Cosmology_Functions
    use Cosmology_Functions_Options
    use Galacticus_Error
    implicit none
    class           (surveyGeometryMoustakas2013PRIMUS), intent(inout)           :: self
    double precision                                   , intent(in   )           :: mass
    integer                                            , intent(in   ), optional :: field
    class           (cosmologyFunctionsClass          ), pointer                 :: cosmologyFunctions_
    double precision                                                             :: redshift           , logarithmicMass
    
    ! Validate field.
    if (.not.present(field)) call Galacticus_Error_Report('moustakas2013PRIMUSDistanceMaximum','field must be specified')
    if (field < 1 .or. field > 5) call Galacticus_Error_Report('moustakas2013PRIMUSDistanceMaximum','1 ≤ field ≤ 5 required')
    ! Find the limiting redshift for this mass completeness limits from Moustakas et al. (2013; Table 2). (See
    ! constraints/dataAnalysis/stellarMassFunctions_PRIMUS_z0_1/massRedshiftRelation.pl for details.)
    logarithmicMass=log10(mass)
    select case (field)
    case (1) ! COSMOS
       redshift=                  +3.51240871481968000d0  &
            &   +logarithmicMass*(-0.94131511297034200d0  &
            &   +logarithmicMass*(+0.06507208866075860d0) &
            &                                           )
    case (2) ! XMM-SXDS
       redshift=                  +2.46068289817352000d0  &
            &   +logarithmicMass*(-0.72960045705258400d0  &
            &   +logarithmicMass*(+0.05422457500058130d0) &
            &                                           )
    case (3) ! XMM-CFHTLS
       redshift=                  -3.60001396783385000d0  &
            &   +logarithmicMass*(+0.50007933123305700d0  &
            &   +logarithmicMass*(-0.00781044013508246d0) &
            &                                           )
    case (4) ! CDFS
       redshift=                  +5.86929723910240000d0  &
            &   +logarithmicMass*(-1.52816338306828000d0  &
            &   +logarithmicMass*(+0.09815249638377170d0) &
            &                                           )
    case (5) ! ELAIS-S1
       redshift=                  +6.87489619768651000d0  &
            &   +logarithmicMass*(-1.65556365363183000d0  &
            &   +logarithmicMass*(+0.10030052053225000d0) &
            &                                           )
    end select
    ! Get the default cosmology functions object.
    cosmologyFunctions_ => cosmologyFunctions()    
    ! Convert from redshift to comoving distance.
    moustakas2013PRIMUSDistanceMaximum                                                &
         &=cosmologyFunctions_%distanceComovingConvert(                               &
         &                                             output  =distanceTypeComoving, &
         &                                             redshift=redshift              &
         &                                            )
    ! Limit the maximum distance.
    moustakas2013PRIMUSDistanceMaximum=min(moustakas2013PRIMUSDistanceMaximum,self%binDistanceMaximum)
    return
  end function moustakas2013PRIMUSDistanceMaximum

  double precision function moustakas2013PRIMUSVolumeMaximum(self,mass,field)
    !% Compute the maximum volume within which a galaxy is visible.
    use Galacticus_Error
    implicit none
    class           (surveyGeometryMoustakas2013PRIMUS), intent(inout)           :: self
    double precision                                   , intent(in   )           :: mass
    integer                                            , intent(in   ), optional :: field
    
    ! Validate field.
    if (.not.present(field)) call Galacticus_Error_Report('moustakas2013PRIMUSDistanceMaximum','field must be specified')
    if (field < 1 .or. field > 5) call Galacticus_Error_Report('moustakas2013PRIMUSDistanceMaximum','1 ≤ field ≤ 5 required')
    ! Compute the volume.
    moustakas2013PRIMUSVolumeMaximum                       &
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
  end function moustakas2013PRIMUSVolumeMaximum

  subroutine moustakas2013PRIMUSInitialize()
    !% Establish data paths and file names for the {\tt moustakas2013PRIMUS} survey geometry class.
    implicit none
    
    if (.not.moustakas2013PRIMUSInitialized) then
       !$omp critical(moustakas2013PRIMUSInitialize)
       if (.not.moustakas2013PRIMUSInitialized) then
          ! Specify input paths and files.
          moustakas2013PRIMUSPath=Galacticus_Input_Path()                                     // &
               &                  "constraints/dataAnalysis/stellarMassFunctions_PRIMUS_z0_1/"
          moustakas2013PRIMUSMangleFiles(1)=moustakas2013PRIMUSPath//"cosmos_field_galex_window_2mask.ply"
          moustakas2013PRIMUSMangleFiles(2)=moustakas2013PRIMUSPath//"xmm_swire_field_galex_window_2mask.ply"
          moustakas2013PRIMUSMangleFiles(3)=moustakas2013PRIMUSPath//"cfhtls_xmm_field_galex_window_2mask.ply"
          moustakas2013PRIMUSMangleFiles(4)=moustakas2013PRIMUSPath//"cdfs_field_galex_window_2mask.ply"
          moustakas2013PRIMUSMangleFiles(5)=moustakas2013PRIMUSPath//"es1_field_galex_window_2mask.ply"
          moustakas2013PRIMUSInitialized   =.true.
       end if
       !$omp end critical(moustakas2013PRIMUSInitialize)
    end if
    return
  end subroutine moustakas2013PRIMUSInitialize

  double precision function moustakas2013PRIMUSSolidAngle(self,field)
    !% Return the solid angle of the \cite{moustakas_primus:_2013} sample.
    use Galacticus_Error
    use File_Utilities
    use String_Handling
    use System_Command
    use IO_HDF5
    implicit none
    class  (surveyGeometryMoustakas2013PRIMUS), intent(inout)           :: self
    integer                                   , intent(in   ), optional :: field
    type   (hdf5Object                       )                          :: solidAngleFile
    
    ! Validate field.
    if (.not.present(field)) call Galacticus_Error_Report('moustakas2013PRIMUSDistanceMaximum','field must be specified')
    if (field < 1 .or. field > 5) &
         & call Galacticus_Error_Report('moustakas2013PRIMUSDistanceMaximum','1 ≤ field ≤ 5 required')
    ! Read solid angles for the fields.
    if (.not.moustakas2013PRIMUSSolidAnglesInitialized) then
       !$omp critical(moustakas2013PRIMUSSolidAnglesInitialize)
       if (.not.moustakas2013PRIMUSSolidAnglesInitialized) then
          ! Construct a solid angle file if one does not already exist.
          if (.not.File_Exists(moustakas2013PRIMUSPath//"solidAngles.hdf5")) then
             call System_Command_Do(Galacticus_Input_Path()//"scripts/aux/mangleRansack.pl "//String_Join(moustakas2013PRIMUSMangleFiles," ")//" "//moustakas2013PRIMUSPath//"solidAngles.hdf5 0")
             if (.not.File_Exists(moustakas2013PRIMUSPath//"solidAngles.hdf5")) &
                  & call Galacticus_Error_Report('moustakas2013PRIMUSSolidAngle','unable to generate solid angles from mangle files')
          end if
          ! Read the solid angles.
          !$omp critical(HDF5_Access)
          call solidAngleFile%openFile         (char(moustakas2013PRIMUSPath//"solidAngles.hdf5"))
          call solidAngleFile%readDatasetStatic('solidAngle',moustakas2013PRIMUSSolidAngles      )
          call solidAngleFile%close            (                                                 )
          !$omp end critical(HDF5_Access)
          ! Record that solid angles are now initialized.
          moustakas2013PRIMUSSolidAnglesInitialized=.true.
       end if
       !$omp end critical(moustakas2013PRIMUSSolidAnglesInitialize)
    end if
    moustakas2013PRIMUSSolidAngle=moustakas2013PRIMUSSolidAngles(field)
    return
  end function moustakas2013PRIMUSSolidAngle

  subroutine moustakas2013PRIMUSWindowFunctions(self,mass1,mass2,gridCount,boxLength,windowFunction1,windowFunction2)
    !% Compute the window function for the survey.
    use, intrinsic :: ISO_C_Binding
    use Galacticus_Error
    implicit none
    class           (surveyGeometryMoustakas2013PRIMUS), intent(inout) :: self
    double precision                                   , intent(in   )                                               :: mass1,mass2
    integer                                            , intent(in   )                                               :: gridCount
    double precision                                   , intent(  out)                                               :: boxLength
    complex         (c_double_complex                 ), intent(  out),     dimension(gridCount,gridCount,gridCount) :: windowFunction1,windowFunction2

    call Galacticus_Error_Report('moustakas2013PRIMUSWindowFunctions','window function construction is not supported')
    return
  end subroutine moustakas2013PRIMUSWindowFunctions

  double precision function moustakas2013PRIMUSAngularPower(self,i,j,l)
    !% Return the angular power $C^{ij}_\ell$ for the \cite{moustakas_primus:_2013} survey.
    use Galacticus_Error
    use File_Utilities
    use String_Handling
    use System_Command
    use IO_HDF5
    implicit none
    class           (surveyGeometryMoustakas2013PRIMUS), intent(inout)             :: self
    integer                                            , intent(in   )             :: i               , j, &
         &                                                                            l
    integer                                                                        :: m               , n
    double precision                                   , allocatable, dimension(:) :: degree
    type            (varying_string                   )                            :: datasetName
    type            (hdf5Object                       )                            :: angularPowerFile

    ! Validate fields.
    if     (                                                                                        &
         &   i < 1 .or. i > 5                                                                       &
         &  .or.                                                                                    &
         &   j < 1 .or. j > 5                                                                       &
         & )                                                                                        &
         & call Galacticus_Error_Report('moustakas2013PRIMUSAngularPower','1 ≤ field ≤ 5 required')
    ! Read angular power spectra.
    if (.not.moustakas2013PRIMUSAngularPowerInitialized) then
       !$omp critical(moustakas2013PRIMUSAngularPowerInitialize)
       if (.not.moustakas2013PRIMUSAngularPowerInitialized) then
          ! Construct an angular power file if one does not already exist.
          if (.not.File_Exists(moustakas2013PRIMUSPath//"angularPower.hdf5")) then
             call System_Command_Do(Galacticus_Input_Path()//"scripts/aux/mangleHarmonize.pl "//String_Join(moustakas2013PRIMUSMangleFiles," ")//" "//moustakas2013PRIMUSPath//"angularPower.hdf5 "//moustaskas2013AngularPowerMaximumL)
             if (.not.File_Exists(moustakas2013PRIMUSPath//"angularPower.hdf5")) &
                  & call Galacticus_Error_Report('moustakas2013PRIMUSAngularPower','unable to generate angular power spectra from mangle files')
          end if
          ! Read the angular power spectra.
          !$omp critical(HDF5_Access)
          call angularPowerFile%openFile(char(moustakas2013PRIMUSPath//"angularPower.hdf5"))
          call angularPowerFile%readDataset('l',degree)
          if (degree(1) /= 0 .or. degree(size(degree)) /= moustaskas2013AngularPowerMaximumL)        &
               & call Galacticus_Error_Report(                                                       &
               &                              'moustakas2013PRIMUSAngularPower'                    , &
               &                              'power spectra do not span expected range of degrees'  &
               &                             )
          do m=1,moustakas2013PRIMUSFields
             do n=m,moustakas2013PRIMUSFields
                datasetName="Cl_"
                datasetName=datasetName//(m-1)//"_"//(n-1)
                call angularPowerFile%readDatasetStatic(char(datasetName),moustakas2013PRIMUSAngularPowerSpectra(:,moustakas2013PRIMUSFieldPairIndex(m,n)))
             end do
          end do
          call angularPowerFile%close()
          !$omp end critical(HDF5_Access)
          ! Record that solid angles are now initialized.
          moustakas2013PRIMUSAngularPowerInitialized=.true.
       end if
       !$omp end critical(moustakas2013PRIMUSAngularPowerInitialize)
    end if
    ! Return the appropriate angular power.
    if (l > moustaskas2013AngularPowerMaximumL) then
       moustakas2013PRIMUSAngularPower=0.0d0
    else
       moustakas2013PRIMUSAngularPower=moustakas2013PRIMUSAngularPowerSpectra(l+1,moustakas2013PRIMUSFieldPairIndex(i,j))
    end if
    return
  end function moustakas2013PRIMUSAngularPower
  
  integer function moustakas2013PRIMUSFieldPairIndex(i,j)
    !% Compute the index of a pair of fields in the \cite{moustakas_primus:_2013} survey.
    implicit none
    integer, intent(in   ) ::  i,  j
    integer                :: ii, jj

    ii=min(i,j)
    jj=max(i,j)
    moustakas2013PRIMUSFieldPairIndex=(ii-1)*(2*moustakas2013PRIMUSFields-ii+2)/2+(jj-ii+1)
    return
  end function moustakas2013PRIMUSFieldPairIndex
  
  integer function moustakas2013PRIMUSAngularPowerMaximumDegree(self)
    !% Return the maximum degree for which angular power is computed for the \cite{moustakas_primus:_2013} survey.
    implicit none
    class(surveyGeometryMoustakas2013PRIMUS), intent(inout) :: self

    moustakas2013PRIMUSAngularPowerMaximumDegree=moustaskas2013AngularPowerMaximumL
    return
  end function moustakas2013PRIMUSAngularPowerMaximumDegree
  
