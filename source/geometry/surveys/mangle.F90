!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!!{RST
Implements an abstract survey geometry using :term:`mangle` polygons.
!!}

  use :: Geometry_Mangle         , only : window
  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  !![
  <surveyGeometry name="surveyGeometryMangle" abstract="yes" docformat="rst">
   <description>
   Implements an abstract survey geometry using :term:`mangle` polygons.
   </description>
  </surveyGeometry>
  !!]
  type, abstract, extends(surveyGeometryClass) :: surveyGeometryMangle
     private
     logical                                                                   :: solidAnglesInitialized          , angularPowerInitialized, &
          &                                                                       windowInitialized
     double precision                            , allocatable, dimension(:  ) :: solidAngles
     double precision                            , allocatable, dimension(:,:) :: angularPowerSpectra
     type            (window                    )                              :: mangleWindow
     class           (randomNumberGeneratorClass), pointer                     :: randomNumberGenerator_ => null()
   contains
     !![
     <methods docformat="rst">
       <method description="Return the directory containing :term:`mangle` files for this survey geometry." method="mangleDirectory" />
       <method description="Return array of :term:`mangle` filenames for this survey geometry." method="mangleFiles" />
       <method description="Initialize an instance of the :term:`mangle` survey geometry class." method="initialize" />
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
    !!{RST
    Internal constructor for the :galacticus-class:`surveyGeometryMangle` survey geometry class.
    !!}
    implicit none
    class(surveyGeometryMangle), intent(inout) :: self

    self%solidAnglesInitialized  =.false.
    self%angularPowerInitialized =.false.
    self%windowInitialized       =.false.
    return
  end subroutine mangleInitialize

  logical function mangleWindowFunctionAvailable(self)
    !!{RST
    Window functions are available whenever a random number generator was supplied (they are constructed by Monte Carlo sampling of the :term:`mangle` mask).
    !!}
    implicit none
    class(surveyGeometryMangle), intent(inout) :: self

    mangleWindowFunctionAvailable=associated(self%randomNumberGenerator_)
    return
  end function mangleWindowFunctionAvailable

  logical function mangleAngularPowerAvailable(self)
    !!{RST
    Return true to indicate that survey angular power is available.
    !!}
    implicit none
    class(surveyGeometryMangle), intent(inout) :: self
    !$GLC attributes unused :: self

    mangleAngularPowerAvailable=.true.
    return
  end function mangleAngularPowerAvailable

  double precision function mangleSolidAngle(self,field)
    !!{RST
    Return the survey solid angle computed from :term:`mangle` polygons.
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
    !!{RST
    Provides window functions for :term:`mangle`-based survey geometries. The survey selection function is Monte Carlo
    sampled: directions drawn uniformly on the sphere are rejection-sampled against the :term:`mangle` angular mask, and
    each accepted direction is assigned depths drawn uniformly in volume within the survey's distance limits (following
    the random-points implementation). A random number generator must have been supplied to the survey geometry's
    constructor.
    !!}
#ifdef FFTW3AVAIL
    use            :: FFTW3                   , only : fftw_plan_dft_3d  , FFTW_FORWARD       , FFTW_ESTIMATE, fftw_execute_dft, &
         &                                             fftw_destroy_plan
#endif
    use            :: Error                   , only : Error_Report
    use, intrinsic :: ISO_C_Binding           , only : c_double_complex  , c_ptr
    use            :: Meshes                  , only : Meshes_Apply_Point, cloudTypeTriangular
    implicit none
    class           (surveyGeometryMangle), intent(inout)                                           :: self
    double precision                      , intent(in   )                                           :: mass1                   , mass2
    integer                               , intent(in   )                                           :: gridCount
    double precision                      , intent(  out)                                           :: boxLength
    complex         (c_double_complex    ), intent(  out), dimension(gridCount,gridCount,gridCount) :: windowFunction1         , windowFunction2
    ! Monte Carlo sample count: at 250,000 samples the per-cell noise of the window on typical grid sizes is at the
    ! percent level, while keeping the cost of the polygon-sampling loop modest (tens of seconds for the SDSS mask).
    integer                               , parameter                                               :: sampleCount=250000
    type            (varying_string      ), allocatable  , dimension(:)                             :: mangleFiles
    complex         (c_double_complex    )               , dimension(gridCount,gridCount,gridCount) :: selectionFunction1      , selectionFunction2
    double precision                                     , dimension(3)                             :: origin                  , direction               , &
         &                                                                                             position1               , position2
    integer                                                                                         :: i
    double precision                                                                                :: comovingDistanceMaximum1, comovingDistanceMaximum2, &
         &                                                                                             comovingDistanceMinimum1, comovingDistanceMinimum2, &
         &                                                                                             distance1               , distance2
#ifdef FFTW3AVAIL
    type            (c_ptr               )                                                          :: plan
#endif
    complex         (c_double_complex    )                                                          :: normalization

#ifdef FFTW3UNAVAIL
    call Error_Report('FFTW3 library is required but was not found'//{introspection:location})
#endif
    if (.not.associated(self%randomNumberGenerator_)) &
         & call Error_Report('window function construction requires a random number generator - supply one to the survey geometry constructor'//{introspection:location})
    ! Initialize the mangle window if necessary.
    if (.not.self%windowInitialized) then
       if (self%fieldCount() > 1) call Error_Report('only single field surveys are supported'//{introspection:location})
       call self%mangleFiles(mangleFiles)
       call self%mangleWindow%read(char(mangleFiles(1)))
       self%windowInitialized=.true.
    end if
    ! Find the comoving distances corresponding to the two masses.
    comovingDistanceMaximum1=self%distanceMaximum(mass1)
    comovingDistanceMaximum2=self%distanceMaximum(mass2)
    comovingDistanceMinimum1=self%distanceMinimum(mass1)
    comovingDistanceMinimum2=self%distanceMinimum(mass2)
    ! Determine a suitable box length for the window function calculation.
    boxLength=3.0d0*max(comovingDistanceMaximum1-comovingDistanceMinimum1,comovingDistanceMaximum2-comovingDistanceMinimum2)
    ! Set up origin.
    origin=([0.5d0,0.5d0,0.0d0]/dble(gridCount)+[0.5d0,0.5d0,0.5d0])*dble(boxLength)
    ! Populate the cube with the survey selection function.
    selectionFunction1=dcmplx(0.0d0,0.0d0)
    selectionFunction2=dcmplx(0.0d0,0.0d0)
    do i=1,sampleCount
       ! Draw a direction uniformly from within the mask (polygons sampled proportional to solid angle; balkanized
       ! mangle windows have disjoint polygons so this is exactly uniform over the mask, at O(1) cost per sample
       ! rather than the O(polygonCount) cost of rejection sampling the full sphere against the mask).
       direction=self%mangleWindow%sampleDirection(self%randomNumberGenerator_)
       ! Choose random distances, uniform in volume within the survey limits.
       distance1=+(                                             &
            &      +self%randomNumberGenerator_%uniformSample() &
            &      *(                                           &
            &        +comovingDistanceMaximum1**3               &
            &        -comovingDistanceMinimum1**3               &
            &       )                                           &
            &        +comovingDistanceMinimum1**3               &
            &     )**(1.0d0/3.0d0)                              &
            &    -min(                                          &
            &         comovingDistanceMinimum1                , &
            &         comovingDistanceMinimum2                  &
            &        )
       distance2=+(                                             &
            &      +self%randomNumberGenerator_%uniformSample() &
            &      *(                                           &
            &        +comovingDistanceMaximum2**3               &
            &        -comovingDistanceMinimum2**3               &
            &       )                                           &
            &        +comovingDistanceMinimum2**3               &
            &     )**(1.0d0/3.0d0)                              &
            &    -min(                                          &
            &         comovingDistanceMinimum1                , &
            &         comovingDistanceMinimum2                  &
            &        )
       position1=distance1*direction+origin
       position2=distance2*direction+origin
       ! Apply to the mesh.
       call Meshes_Apply_Point(selectionFunction1,boxLength,position1,pointWeight=cmplx(1.0d0,0.0d0,kind=c_double_complex),cloudType=cloudTypeTriangular)
       call Meshes_Apply_Point(selectionFunction2,boxLength,position2,pointWeight=cmplx(1.0d0,0.0d0,kind=c_double_complex),cloudType=cloudTypeTriangular)
    end do
    ! Take the Fourier transform of the selection function.
#ifdef FFTW3AVAIL
    plan=fftw_plan_dft_3d(gridCount,gridCount,gridCount,selectionFunction1,windowFunction1,FFTW_FORWARD,FFTW_ESTIMATE)
    call fftw_execute_dft(plan,selectionFunction1,windowFunction1)
    call fftw_execute_dft(plan,selectionFunction2,windowFunction2)
    call fftw_destroy_plan(plan)
#endif
    ! Normalize the window function.
    normalization=windowFunction1(1,1,1)
    if (real(normalization) > 0.0d0) windowFunction1=windowFunction1/normalization
    normalization=windowFunction2(1,1,1)
    if (real(normalization) > 0.0d0) windowFunction2=windowFunction2/normalization
    return
  end subroutine mangleWindowFunctions

  double precision function mangleAngularPower(self,i,j,l)
    !!{RST
    Return the survey angular power :math:`C^{ij}_\ell` from :term:`mangle` polygons.
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
    !!{RST
    Compute the index of a pair of fields in :term:`mangle`-based survey geometries.
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
    !!{RST
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
