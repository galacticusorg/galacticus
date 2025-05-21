!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

!!{
Implements an abstract survey geometry using \gls{mangle} polygons.
!!}

  use :: Geometry_Mangle, only : window

  !![
  <surveyGeometry name="surveyGeometryMangle" abstract="yes">
   <description>Implements an abstract survey geometry using \gls{mangle} polygons.</description>
  </surveyGeometry>
  !!]
  type, abstract, extends(surveyGeometryClass) :: surveyGeometryMangle
     private
     logical                                               :: solidAnglesInitialized, angularPowerInitialized, windowInitialized
     double precision        , allocatable, dimension(:  ) :: solidAngles
     double precision        , allocatable, dimension(:,:) :: angularPowerSpectra
     type            (window)                              :: mangleWindow
   contains
     !![
     <methods>
       <method description="Return the directory containing \gls{mangle} files for this survey geometry." method="mangleDirectory" />
       <method description="Return array of \gls{mangle} filenames for this survey geometry." method="mangleFiles" />
       <method description="Initialize an instance of the \gls{mangle} survey geometry class." method="initialize" />
     </methods>
     !!]
     procedure                                  :: windowFunctionAvailable => mangleWindowFunctionAvailable
     procedure                                  :: angularPowerAvailable   => mangleAngularPowerAvailable
     procedure                                  :: solidAngle              => mangleSolidAngle
     procedure                                  :: windowFunctions         => mangleWindowFunctions
     procedure                                  :: angularPower            => mangleAngularPower
     procedure                                  :: pointIncluded           => manglePointIncluded
     procedure                                  :: initialize              => mangleInitialize
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
       type (varying_string      ), allocatable, dimension(:), intent(inout) :: mangleFiles
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

  subroutine mangleInitialize(self)
    !!{
    Internal constructor for the \refClass{surveyGeometryMangle} survey geometry class.
    !!}
    implicit none
    class(surveyGeometryMangle), intent(inout) :: self

    self%solidAnglesInitialized  =.false.
    self%angularPowerInitialized =.false.
    self%windowInitialized       =.false.
    return
  end subroutine mangleInitialize

  logical function mangleWindowFunctionAvailable(self)
    !!{
    Return false to indicate that survey window function is not available.
    !!}
    implicit none
    class(surveyGeometryMangle), intent(inout) :: self
    !$GLC attributes unused :: self

    mangleWindowFunctionAvailable=.false.
    return
  end function mangleWindowFunctionAvailable

  logical function mangleAngularPowerAvailable(self)
    !!{
    Return true to indicate that survey angular power is available.
    !!}
    implicit none
    class(surveyGeometryMangle), intent(inout) :: self
    !$GLC attributes unused :: self

    mangleAngularPowerAvailable=.true.
    return
  end function mangleAngularPowerAvailable

  double precision function mangleSolidAngle(self,field)
    !!{
    Return the survey solid angle computed from \gls{mangle} polygons.
    !!}
    use :: Error          , only : Error_Report
    use :: Geometry_Mangle, only : geometryMangleSolidAngle
    use :: String_Handling, only : operator(//)
    implicit none
    class  (surveyGeometryMangle), intent(inout)               :: self
    integer                      , intent(in   ), optional     :: field
    type   (varying_string      ), allocatable  , dimension(:) :: mangleFiles
    integer                                                    :: fieldActual
    type   (varying_string      )                              :: message

    ! Validate field.
    if (.not.present(field)) then
       if (self%fieldCount() > 1) call Error_Report('field must be specified'//{introspection:location})
       fieldActual=1
    else
       fieldActual=field
    end if
    if (fieldActual < 1 .or. fieldActual > self%fieldCount()) then
       message='1 ≤ field ≤ '
       message=message//self%fieldCount()//' required'
       call Error_Report(message//{introspection:location})
    end if
    ! Read solid angles for the fields.
    if (.not.self%solidAnglesInitialized) then
       if (.not.self%solidAnglesInitialized) then
          call self%mangleFiles(mangleFiles)
          self%solidAngles           =geometryMangleSolidAngle(mangleFiles,char(self%mangleDirectory()//"solidAngles.hdf5"))
          self%solidAnglesInitialized=.true.
       end if
    end if
    mangleSolidAngle=self%solidAngles(fieldActual)
    return
  end function mangleSolidAngle

  subroutine mangleWindowFunctions(self,mass1,mass2,gridCount,boxLength,windowFunction1,windowFunction2)
    !!{
    Provides window functions for \gls{mangle}-based survey geometries.
    !!}
    use            :: Error        , only : Error_Report
    use, intrinsic :: ISO_C_Binding, only : c_double_complex
    implicit none
    class           (surveyGeometryMangle), intent(inout)                                           :: self
    double precision                      , intent(in   )                                           :: mass1,mass2
    integer                               , intent(in   )                                           :: gridCount
    double precision                      , intent(  out)                                           :: boxLength
    complex         (c_double_complex    ), intent(  out), dimension(gridCount,gridCount,gridCount) :: windowFunction1,windowFunction2
    !$GLC attributes unused :: self, mass1, mass2, gridCount, boxLength, windowFunction1, windowFunction2

    call Error_Report('window function construction is not supported'//{introspection:location})
    return
  end subroutine mangleWindowFunctions

  double precision function mangleAngularPower(self,i,j,l)
    !!{
    Return the survey angular power $C^{ij}_\ell$ from \gls{mangle} polygons.
    !!}
    use :: Error          , only : Error_Report
    use :: Geometry_Mangle, only : geometryMangleAngularPower
    use :: String_Handling, only : operator(//)
    implicit none
    class           (surveyGeometryMangle), intent(inout)               :: self
    integer                               , intent(in   )               :: i          , j, &
         &                                                                 l
    type            (varying_string      ), allocatable  , dimension(:) :: mangleFiles
    type            (varying_string      )                              :: message

    ! Validate fields.
    if     (                                  &
         &   i < 1 .or. i > self%fieldCount() &
         &  .or.                              &
         &   j < 1 .or. j > self%fieldCount() &
         & ) then
       message='1 ≤ field ≤ '
       message=message//self%fieldCount()//' required'
       call Error_Report(message//{introspection:location})
    end if
    ! Read angular power spectra.
    if (.not.self%angularPowerInitialized) then
       if (.not.self%angularPowerInitialized) then
          call self%mangleFiles(mangleFiles)
          self%angularPowerSpectra    =geometryMangleAngularPower(mangleFiles,self%angularPowerMaximumDegree(),char(self%mangleDirectory()//"angularPower.hdf5"))
          self%angularPowerInitialized=.true.
       end if
    end if
    ! Return the appropriate angular power.
    if (l > self%angularPowerMaximumDegree()) then
       mangleAngularPower=0.0d0
    else
       mangleAngularPower=self%angularPowerSpectra(mangleFieldPairIndex(self,i,j),l+1)
    end if
    return
  end function mangleAngularPower

  integer function mangleFieldPairIndex(self,i,j)
    !!{
    Compute the index of a pair of fields in \gls{mangle}-based survey geometries.
    !!}
    implicit none
    class  (surveyGeometryMangle), intent(inout) :: self
    integer                      , intent(in   ) ::  i,  j
    integer                                      :: ii, jj

    ii=min(i,j)
    jj=max(i,j)
    mangleFieldPairIndex=(ii-1)*(2*self%fieldCount()-ii+2)/2+(jj-ii+1)
    return
  end function mangleFieldPairIndex

  logical function manglePointIncluded(self,point,mass)
    !!{
    Return true if a point is included in the survey geometry.
    !!}
    use :: Error  , only : Error_Report
    use :: Vectors, only : Vector_Magnitude
    implicit none
    class           (surveyGeometryMangle), intent(inout)               :: self
    double precision                      , intent(in   ), dimension(3) :: point
    double precision                      , intent(in   )               :: mass
    type            (varying_string      ), allocatable  , dimension(:) :: mangleFiles
    double precision                                                    :: pointDistance

    ! Initialize the mangle window if necessary.
    if (.not.self%windowInitialized) then
       if (.not.self%windowInitialized) then
          if (self%fieldCount() > 1) call Error_Report('only single field surveys are supported'//{introspection:location})
          call self%mangleFiles(mangleFiles)
          call self%mangleWindow%read(char(mangleFiles(1)))
          self%windowInitialized=.true.
       end if
    end if
    ! Get the distance to the point.
    pointDistance=Vector_Magnitude(point)
    ! Compute if point lies within survey bounds.
    manglePointIncluded=                              &
         & pointDistance > self%distanceMinimum(mass) &
         & .and.                                      &
         & pointDistance < self%distanceMaximum(mass)
    if (manglePointIncluded) manglePointIncluded=self%mangleWindow%pointIncluded(point)
    return
  end function manglePointIncluded
