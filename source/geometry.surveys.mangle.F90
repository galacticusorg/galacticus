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

!% Implements an abstract survey geometry using \gls{mangle} polygons.
  
  !# <surveyGeometry name="surveyGeometryMangle">
  !#  <description>Implements an abstract survey geometry using \gls{mangle} polygons.</description>
  !#  <abstract>yes</abstract>
  !# </surveyGeometry>

  use Galacticus_Input_Paths

  type, abstract, extends(surveyGeometryClass) :: surveyGeometryMangle
     logical                                       :: solidAnglesInitialized, angularPowerInitialized
     double precision, allocatable, dimension(:  ) :: solidAngles
     double precision, allocatable, dimension(:,:) :: angularPowerSpectra
   contains
     !@ <objectMethods>
     !@   <object>surveyGeometryMangle</object>
     !@   <objectMethod>
     !@     <method>mangleDirectory</method>
     !@     <type>\textcolor{red}{\textless type(varying\_string) \textgreater}</type>
     !@     <arguments></arguments>
     !@     <description>Return the directory containing \gls{mangle} files for this survey geometry.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>mangleFiles</method>
     !@     <type>\textcolor{red}{\textless type(varying\_string)(:) \textgreater}</type>
     !@     <arguments></arguments>
     !@     <description>Return array of \gls{mangle} filenames for this survey geometry.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                                  :: windowFunctionAvailable => mangleWindowFunctionAvailable
     procedure                                  :: angularPowerAvailable   => mangleAngularPowerAvailable
     procedure                                  :: solidAngle              => mangleSolidAngle
     procedure                                  :: windowFunctions         => mangleWindowFunctions
     procedure                                  :: angularPower            => mangleAngularPower
     procedure(mangleMangleDirectory), deferred :: mangleDirectory
     procedure(mangleMangleFiles    ), deferred :: mangleFiles
  end type surveyGeometryMangle

  abstract interface
     function mangleMangleDirectory(self)
       import varying_string, surveyGeometryMangle
       type (varying_string      )                :: mangleMangleDirectory
       class(surveyGeometryMangle), intent(inout) :: self
     end function mangleMangleDirectory
  end interface

  abstract interface
     subroutine mangleMangleFiles(self,mangleFiles)
       import varying_string, surveyGeometryMangle
       type (varying_string      ), allocatable, dimension(:), intent(  out) :: mangleFiles
       class(surveyGeometryMangle)                           , intent(inout) :: self
     end subroutine mangleMangleFiles
  end interface

  abstract interface
     integer function mangleAngularPowerMaximumDegree(self)
       import surveyGeometryMangle
       class(surveyGeometryMangle), intent(inout) :: self
     end function mangleAngularPowerMaximumDegree
  end interface

contains

  logical function mangleWindowFunctionAvailable(self)
    !% Return false to indicate that survey window function is not available.
    implicit none
    class(surveyGeometryMangle), intent(inout) :: self

    mangleWindowFunctionAvailable=.false.
    return
  end function mangleWindowFunctionAvailable

  logical function mangleAngularPowerAvailable(self)
    !% Return true to indicate that survey angular power is available.
    implicit none
    class(surveyGeometryMangle), intent(inout) :: self

    mangleAngularPowerAvailable=.true.
    return
  end function mangleAngularPowerAvailable

  double precision function mangleSolidAngle(self,field)
    !% Return the survey solid angle computed from \gls{mangle} polygons.
    use Galacticus_Error
    use File_Utilities
    use String_Handling
    use System_Command
    use IO_HDF5
    implicit none
    class  (surveyGeometryMangle), intent(inout)               :: self
    integer                      , intent(in   ), optional     :: field
    type   (varying_string      ), allocatable  , dimension(:) :: mangleFiles
    type   (hdf5Object          )                              :: solidAngleFile
    integer                                                    :: fieldActual
    type   (varying_string      )                              :: message
    
    ! Validate field.
    if (.not.present(field)) then
       if (self%fieldCount() > 1) call Galacticus_Error_Report('mangleDistanceMaximum','field must be specified')
       fieldActual=1
    else
       fieldActual=field
    end if
    if (fieldActual < 1 .or. fieldActual > self%fieldCount()) then
       message='1 ≤ field ≤ '
       message=message//self%fieldCount()//' required'
       call Galacticus_Error_Report('mangleSolidAngle',message)
    end if
    ! Read solid angles for the fields.
    if (.not.self%solidAnglesInitialized) then
       !$omp critical(mangleSolidAnglesInitialize)
       if (.not.self%solidAnglesInitialized) then
          ! Construct a solid angle file if one does not already exist.
          if (.not.File_Exists(self%mangleDirectory()//"solidAngles.hdf5")) then
             call self%mangleFiles(mangleFiles)
             call System_Command_Do(Galacticus_Input_Path()//"scripts/aux/mangleRansack.pl "//String_Join(mangleFiles," ")//" "//self%mangleDirectory()//"solidAngles.hdf5 0")
             if (.not.File_Exists(self%mangleDirectory()//"solidAngles.hdf5")) &
                  & call Galacticus_Error_Report('mangleSolidAngle','unable to generate solid angles from mangle files')
          end if
          ! Read the solid angles.
          !$omp critical(HDF5_Access)
          call solidAngleFile%openFile   (char(self%mangleDirectory()//"solidAngles.hdf5"))
          call solidAngleFile%readDataset('solidAngle',self%solidAngles                   )
          call solidAngleFile%close      (                                                )
          !$omp end critical(HDF5_Access)
          ! Record that solid angles are now initialized.
          self%solidAnglesInitialized=.true.
       end if
       !$omp end critical(mangleSolidAnglesInitialize)
    end if
    mangleSolidAngle=self%solidAngles(fieldActual)
    return
  end function mangleSolidAngle

  subroutine mangleWindowFunctions(self,mass1,mass2,gridCount,boxLength,windowFunction1,windowFunction2)
    !% Provides window functions for \gls{mangle}-based survey geometries.
    use, intrinsic :: ISO_C_Binding
    use Galacticus_Error
    implicit none
    class           (surveyGeometryMangle), intent(inout)                                           :: self
    double precision                      , intent(in   )                                           :: mass1,mass2
    integer                               , intent(in   )                                           :: gridCount
    double precision                      , intent(  out)                                           :: boxLength
    complex         (c_double_complex    ), intent(  out), dimension(gridCount,gridCount,gridCount) :: windowFunction1,windowFunction2

    call Galacticus_Error_Report('mangleWindowFunctions','window function construction is not supported')
    return
  end subroutine mangleWindowFunctions

  double precision function mangleAngularPower(self,i,j,l)
    !% Return the survey angular power $C^{ij}_\ell$ from \gls{mangle} polygons.
    use Galacticus_Error
    use File_Utilities
    use String_Handling
    use System_Command
    use IO_HDF5
    use Memory_Management
    implicit none
    class           (surveyGeometryMangle), intent(inout)               :: self
    integer                               , intent(in   )               :: i               , j, &
         &                                                                 l
    integer                                                             :: m               , n
    double precision                      , allocatable  , dimension(:) :: degree
    type            (varying_string      ), allocatable  , dimension(:) :: mangleFiles
    type            (varying_string      )                              :: datasetName     , message
    type            (hdf5Object          )                              :: angularPowerFile

    ! Validate fields.
    if     (                                  &
         &   i < 1 .or. i > self%fieldCount() &
         &  .or.                              &
         &   j < 1 .or. j > self%fieldCount() &
         & ) then
       message='1 ≤ field ≤ '
       message=message//self%fieldCount()//' required'
       call Galacticus_Error_Report('mangleAngularPower',message)
    end if
    ! Read angular power spectra.
    if (.not.self%angularPowerInitialized) then
       !$omp critical(mangleAngularPowerInitialize)
       if (.not.self%angularPowerInitialized) then
          ! Construct an angular power file if one does not already exist.
          if (.not.File_Exists(self%mangleDirectory()//"angularPower.hdf5")) then
             call self%mangleFiles(mangleFiles)
             call System_Command_Do(Galacticus_Input_Path()//"scripts/aux/mangleHarmonize.pl "//String_Join(mangleFiles," ")//" "//self%mangleDirectory()//"angularPower.hdf5 "//self%angularPowerMaximumDegree())
             if (.not.File_Exists(self%mangleDirectory()//"angularPower.hdf5")) &
                  & call Galacticus_Error_Report('mangleAngularPower','unable to generate angular power spectra from mangle files')
          end if
          ! Read the angular power spectra.
          !$omp critical(HDF5_Access)
          call angularPowerFile%openFile(char(self%mangleDirectory()//"angularPower.hdf5"))
          call angularPowerFile%readDataset('l',degree)
          if (degree(1) /= 0 .or. degree(size(degree)) /= self%angularPowerMaximumDegree())               &
               & call Galacticus_Error_Report(                                                       &
               &                              'mangleAngularPower'                                 , &
               &                              'power spectra do not span expected range of degrees'  &
               &                             )
          call Alloc_Array(self%angularPowerSpectra,[self%angularPowerMaximumDegree()+1,self%fieldCount()*(self%fieldCount()+1)])
          do m=1,self%fieldCount()
             do n=m,self%fieldCount()
                datasetName="Cl_"
                datasetName=datasetName//(m-1)//"_"//(n-1)
                call angularPowerFile%readDatasetStatic(char(datasetName),self%angularPowerSpectra(:,mangleFieldPairIndex(self,m,n)))
             end do
          end do
          call angularPowerFile%close()
          !$omp end critical(HDF5_Access)
          ! Record that solid angles are now initialized.
          self%angularPowerInitialized=.true.
       end if
       !$omp end critical(mangleAngularPowerInitialize)
    end if
    ! Return the appropriate angular power.
    if (l > self%angularPowerMaximumDegree()) then
       mangleAngularPower=0.0d0
    else
       mangleAngularPower=self%angularPowerSpectra(l+1,mangleFieldPairIndex(self,i,j))
    end if
    return
  end function mangleAngularPower
  
  integer function mangleFieldPairIndex(self,i,j)
    !% Compute the index of a pair of fields in \gls{mangle}-based survey geometries.
    implicit none
    class  (surveyGeometryMangle), intent(inout) :: self
    integer                      , intent(in   ) ::  i,  j
    integer                                      :: ii, jj

    ii=min(i,j)
    jj=max(i,j)
    mangleFieldPairIndex=(ii-1)*(2*self%fieldCount()-ii+2)/2+(jj-ii+1)
    return
  end function mangleFieldPairIndex

