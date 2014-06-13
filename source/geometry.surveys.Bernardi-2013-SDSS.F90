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

!% Implements the geometry of the PRIMUS survey used by \cite{bernardi_massive_2013}.
  
  !# <surveyGeometry name="surveyGeometryBernardi2013SDSS">
  !#  <description>Implements the geometry of the PRIMUS survey of \cite{bernardi_massive_2013}.</description>
  !# </surveyGeometry>

  use Galacticus_Input_Paths

  type, extends(surveyGeometryClass) :: surveyGeometryBernardi2013SDSS
   contains
     procedure :: windowFunctionAvailable   => bernardi2013SDSSWindowFunctionAvailable
     procedure :: angularPowerAvailable     => bernardi2013SDSSAngularPowerAvailable
     procedure :: fieldCount                => bernardi2013SDSSFieldCount
     procedure :: distanceMaximum           => bernardi2013SDSSDistanceMaximum
     procedure :: solidAngle                => bernardi2013SDSSSolidAngle
     procedure :: windowFunctions           => bernardi2013SDSSWindowFunctions
     procedure :: angularPower              => bernardi2013SDSSAngularPower
     procedure :: angularPowerMaximumDegree => bernardi2013SDSSAngularPowerMaximumDegree
  end type surveyGeometryBernardi2013SDSS

  interface surveyGeometryBernardi2013SDSS
     !% Constructors for the \cite{bernardi_massive_2013} survey geometry class.
     module procedure bernardi2013SDSSDefaultConstructor
  end interface surveyGeometryBernardi2013SDSS

  ! Paths and file names for mangle polygon files.
  logical                                                                             :: bernardi2013SDSSInitialized            =.false. 
  type            (varying_string)                                                    :: bernardi2013SDSSPath
  type            (varying_string)                                                    :: bernardi2013SDSSMangleFile

  ! Solid angles.
  logical                                                                             :: bernardi2013SDSSSolidAngleInitialized  =.false.
  double precision                , dimension(                                     1) :: bernardi2013SDSSSolidAngles

  ! Angular power spectra.
  integer                         , parameter                                         :: bernardi2013SDSSAngularPowerMaximumL   =360
  logical                                                                             :: bernardi2013SDSSAngularPowerInitialized=.false.
  double precision                , dimension(bernardi2013SDSSAngularPowerMaximumL+1) :: bernardi2013SDSSAngularPowerSpectra

contains

  function bernardi2013SDSSDefaultConstructor()
    !% Default constructor for the \cite{bernardi_massive_2013} conditional mass function class.
    use Input_Parameters
    implicit none
    type(surveyGeometryBernardi2013SDSS) :: bernardi2013SDSSDefaultConstructor

    call bernardi2013SDSSInitialize()
    return
  end function bernardi2013SDSSDefaultConstructor

  integer function bernardi2013SDSSFieldCount(self)
    !% Return the number of fields in this sample.
    implicit none
    class(surveyGeometryBernardi2013SDSS), intent(inout) :: self

    bernardi2013SDSSFieldCount=1
    return
  end function bernardi2013SDSSFieldCount
    
  logical function bernardi2013SDSSWindowFunctionAvailable(self)
    !% Return false to indicate that survey window function is available.
    implicit none
    class(surveyGeometryBernardi2013SDSS), intent(inout) :: self

    bernardi2013SDSSWindowFunctionAvailable=.false.
    return
  end function bernardi2013SDSSWindowFunctionAvailable

  logical function bernardi2013SDSSAngularPowerAvailable(self)
    !% Return true to indicate that survey angular power is not available.
    implicit none
    class(surveyGeometryBernardi2013SDSS), intent(inout) :: self

    bernardi2013SDSSAngularPowerAvailable=.true.
    return
  end function bernardi2013SDSSAngularPowerAvailable

  double precision function bernardi2013SDSSDistanceMaximum(self,mass,field)
    !% Compute the maximum distance at which a galaxy is visible.
    use Cosmology_Functions
    use Cosmology_Functions_Options
    use Galacticus_Error
    implicit none
    class           (surveyGeometryBernardi2013SDSS), intent(inout)           :: self
    double precision                                , intent(in   )           :: mass
    integer                                         , intent(in   ), optional :: field
    class           (cosmologyFunctionsClass       ), pointer                 :: cosmologyFunctions_
    double precision                                                          :: redshift           , logarithmicMass
    

    ! Find the limiting redshift for this mass completeness limits from Moustakas et al. (2013; Table 2). (See
    ! constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07_Bernardi/massDistanceRelation.pl for details.)
    logarithmicMass=log10(mass)
    bernardi2013SDSSDistanceMaximum                                                         &
         &                         =10.0d0**(                                               &
         &                                                         +1282.1065495948200000d0 &
         &                                   +logarithmicMass*(    - 626.6442739444630000d0 &
         &                                   +logarithmicMass* (   + 122.0914916099620000d0 &
         &                                   +logarithmicMass*  (  -  11.8431000301984000d0 &
         &                                   +logarithmicMass*   ( +   0.5723990783953920d0 &
         &                                   +logarithmicMass*    (-   0.0110301089727899d0 &
         &                                                        )                         &
         &                                                       )                          &
         &                                                      )                           &
         &                                                     )                            &
         &                                                    )                             &
         &                                  )
    return
  end function bernardi2013SDSSDistanceMaximum

  subroutine bernardi2013SDSSInitialize()
    !% Establish data paths and file names for the {\tt bernardi2013SDSS} survey geometry class.
    implicit none
    
    if (.not.bernardi2013SDSSInitialized) then
       !$omp critical(bernardi2013SDSSInitialize)
       if (.not.bernardi2013SDSSInitialized) then
          ! Specify input paths and files.
          bernardi2013SDSSPath=Galacticus_Input_Path()                                               // &
               &                  "constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07_Bernardi/"
          bernardi2013SDSSMangleFile =bernardi2013SDSSPath//"sdss_dr72safe0_res6d.pol"
          bernardi2013SDSSInitialized=.true.
       end if
       !$omp end critical(bernardi2013SDSSInitialize)
    end if
    return
  end subroutine bernardi2013SDSSInitialize

  double precision function bernardi2013SDSSSolidAngle(self,field)
    !% Return the solid angle of the \cite{bernardi_massive_2013} sample.
    use Galacticus_Error
    use File_Utilities
    use String_Handling
    use System_Command
    use IO_HDF5
    implicit none
    class  (surveyGeometryBernardi2013SDSS), intent(inout)           :: self
    integer                                , intent(in   ), optional :: field
    type   (hdf5Object                    )                          :: solidAngleFile
    
    ! Read solid angles for the field.
    if (.not.bernardi2013SDSSSolidAngleInitialized) then
       !$omp critical(bernardi2013SDSSSolidAngleInitialize)
       if (.not.bernardi2013SDSSSolidAngleInitialized) then
          ! Construct a solid angle file if one does not already exist.
          if (.not.File_Exists(bernardi2013SDSSPath//"solidAngle.hdf5")) then
             call System_Command_Do(Galacticus_Input_Path()//"scripts/aux/mangleRansack.pl "//bernardi2013SDSSMangleFile//" "//bernardi2013SDSSPath//"solidAngle.hdf5 0")
             if (.not.File_Exists(bernardi2013SDSSPath//"solidAngle.hdf5")) &
                  & call Galacticus_Error_Report('bernardi2013SDSSSolidAngle','unable to generate solid angle from mangle file')
          end if
          ! Read the solid angle.
          !$omp critical(HDF5_Access)
          call solidAngleFile%openFile         (char(bernardi2013SDSSPath//"solidAngle.hdf5"))
          call solidAngleFile%readDatasetStatic('solidAngle',bernardi2013SDSSSolidAngles     )
          call solidAngleFile%close            (                                             )
          !$omp end critical(HDF5_Access)
          ! Record that solid angle are now initialized.
          bernardi2013SDSSSolidAngleInitialized=.true.
       end if
       !$omp end critical(bernardi2013SDSSSolidAngleInitialize)
    end if
    bernardi2013SDSSSolidAngle=bernardi2013SDSSSolidAngles(1)
    return
  end function bernardi2013SDSSSolidAngle

  subroutine bernardi2013SDSSWindowFunctions(self,mass1,mass2,gridCount,boxLength,windowFunction1,windowFunction2)
    !% Compute the window function for the survey.
    use, intrinsic :: ISO_C_Binding
    use Galacticus_Error
    implicit none
    class           (surveyGeometryBernardi2013SDSS), intent(inout)                                               :: self
    double precision                                , intent(in   )                                               :: mass1,mass2
    integer                                         , intent(in   )                                               :: gridCount
    double precision                                , intent(  out)                                               :: boxLength
    complex         (c_double_complex              ), intent(  out),     dimension(gridCount,gridCount,gridCount) :: windowFunction1,windowFunction2

    call Galacticus_Error_Report('bernardi2013SDSSWindowFunctions','window function construction is not supported')
    return
  end subroutine bernardi2013SDSSWindowFunctions

  double precision function bernardi2013SDSSAngularPower(self,i,j,l)
    !% Return the angular power $C^{ij}_\ell$ for the \cite{bernardi_massive_2013} survey.
    use Galacticus_Error
    use File_Utilities
    use String_Handling
    use System_Command
    use IO_HDF5
    implicit none
    class           (surveyGeometryBernardi2013SDSS), intent(inout)             :: self
    integer                                         , intent(in   )             :: i               , j, &
         &                                                                         l
    integer                                                                     :: m               , n
    double precision                                , allocatable, dimension(:) :: degree
    type            (varying_string                )                            :: datasetName
    type            (hdf5Object                    )                            :: angularPowerFile

    ! Read angular power spectra.
    if (.not.bernardi2013SDSSAngularPowerInitialized) then
       !$omp critical(bernardi2013SDSSAngularPowerInitialize)
       if (.not.bernardi2013SDSSAngularPowerInitialized) then
          ! Construct an angular power file if one does not already exist.
          if (.not.File_Exists(bernardi2013SDSSPath//"angularPower.hdf5")) then
             call System_Command_Do(Galacticus_Input_Path()//"scripts/aux/mangleHarmonize.pl "//bernardi2013SDSSMangleFile//" "//bernardi2013SDSSPath//"angularPower.hdf5 "//bernardi2013SDSSAngularPowerMaximumL)
             if (.not.File_Exists(bernardi2013SDSSPath//"angularPower.hdf5")) &
                  & call Galacticus_Error_Report('bernardi2013SDSSAngularPower','unable to generate angular power spectra from mangle files')
          end if
          ! Read the angular power spectra.
          !$omp critical(HDF5_Access)
          call angularPowerFile%openFile(char(bernardi2013SDSSPath//"angularPower.hdf5"))
          call angularPowerFile%readDataset('l',degree)
          if (degree(1) /= 0 .or. degree(size(degree)) /= bernardi2013SDSSAngularPowerMaximumL)      &
               & call Galacticus_Error_Report(                                                       &
               &                              'bernardi2013SDSSAngularPower'                       , &
               &                              'power spectra do not span expected range of degrees'  &
               &                             )
          call angularPowerFile%readDatasetStatic("Cl_0_0",bernardi2013SDSSAngularPowerSpectra(:))
          call angularPowerFile%close()
          !$omp end critical(HDF5_Access)
          ! Record that solid angles are now initialized.
          bernardi2013SDSSAngularPowerInitialized=.true.
       end if
       !$omp end critical(bernardi2013SDSSAngularPowerInitialize)
    end if
    ! Return the appropriate angular power.
    if (l > bernardi2013SDSSAngularPowerMaximumL) then
       bernardi2013SDSSAngularPower=0.0d0
    else
       bernardi2013SDSSAngularPower=bernardi2013SDSSAngularPowerSpectra(l+1)
    end if
    return
  end function bernardi2013SDSSAngularPower
  
  integer function bernardi2013SDSSAngularPowerMaximumDegree(self)
    !% Return the maximum degree for which angular power is computed for the \cite{bernardi_massive_2013} survey.
    implicit none
    class(surveyGeometryBernardi2013SDSS), intent(inout) :: self

    bernardi2013SDSSAngularPowerMaximumDegree=bernardi2013SDSSAngularPowerMaximumL
    return
  end function bernardi2013SDSSAngularPowerMaximumDegree
  
