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

!!{
Implements survey geometries defined by random points.
!!}

  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  !![
  <surveyGeometry name="surveyGeometryRandomPoints" abstract="yes">
   <description>Implements survey geometries defined by random points.</description>
  </surveyGeometry>
  !!]
  type, abstract, extends(surveyGeometryClass) :: surveyGeometryRandomPoints
     private
     logical                                                                 :: geometryInitialized
     double precision                            , allocatable, dimension(:) :: randomTheta                     , randomPhi
     class           (randomNumberGeneratorClass), pointer                   :: randomNumberGenerator_ => null()
   contains
     !![
     <methods>
       <method description="Initialize arrays of random points to define the survey angular geometry." method="randomsInitialize" />
     </methods>
     !!]
     procedure                                          :: windowFunctionAvailable => randomPointsWindowFunctionAvailable
     procedure                                          :: angularPowerAvailable   => randomPointsAngularPowerAvailable
     procedure                                          :: windowFunctions         => randomPointsWindowFunctions
     procedure                                          :: angularPower            => randomPointsAngularPower
     procedure                                          :: pointIncluded           => randomPointIncluded
     procedure(randomPointsRandomsInitialize), deferred :: randomsInitialize
  end type surveyGeometryRandomPoints

  abstract interface
     subroutine randomPointsRandomsInitialize(self)
       import surveyGeometryRandomPoints
       class(surveyGeometryRandomPoints), intent(inout) :: self
     end subroutine randomPointsRandomsInitialize
  end interface

contains

  logical function randomPointsWindowFunctionAvailable(self)
    !!{
    Return true to indicate that survey window function is available.
    !!}
    implicit none
    class(surveyGeometryRandomPoints), intent(inout) :: self
    !$GLC attributes unused :: self

    randomPointsWindowFunctionAvailable=.true.
    return
  end function randomPointsWindowFunctionAvailable

  logical function randomPointsAngularPowerAvailable(self)
    !!{
    Return false to indicate that survey angular power is not available.
    !!}
    implicit none
    class(surveyGeometryRandomPoints), intent(inout) :: self
    !$GLC attributes unused :: self

    randomPointsAngularPowerAvailable=.false.
    return
  end function randomPointsAngularPowerAvailable

  subroutine randomPointsWindowFunctions(self,mass1,mass2,gridCount,boxLength,windowFunction1,windowFunction2)
    !!{
    Compute the window function for the survey.
    !!}
#ifdef FFTW3AVAIL
    use            :: FFTW3        , only : fftw_plan_dft_3d       , FFTW_FORWARD       , FFTW_ESTIMATE, fftw_execute_dft, &
         &                                  fftw_destroy_plan
#endif
    use            :: Error        , only : Error_Report
    use, intrinsic :: ISO_C_Binding, only : c_double_complex       , c_ptr
    use            :: Meshes       , only : Meshes_Apply_Point     , cloudTypeTriangular
    implicit none
    class           (surveyGeometryRandomPoints), intent(inout)                                           :: self
    double precision                            , intent(in   )                                           :: mass1                   , mass2
    integer                                     , intent(in   )                                           :: gridCount
    double precision                            , intent(  out)                                           :: boxLength
    complex         (c_double_complex          ), intent(  out), dimension(gridCount,gridCount,gridCount) :: windowFunction1         , windowFunction2
    double precision                            ,                dimension(3                            ) :: origin                  , position1              , &
         &                                                                                                   position2
    integer                                                                                               :: i
    double precision                                                                                      :: comovingDistanceMaximum1, comovingDistanceMaximum2, &
         &                                                                                                   comovingDistanceMinimum1, comovingDistanceMinimum2, &
         &                                                                                                   distance1               , distance2
#ifdef FFTW3AVAIL
    type            (c_ptr                     )                                                          :: plan
#endif
    complex         (c_double_complex          ),                dimension(gridCount,gridCount,gridCount) :: selectionFunction1      , selectionFunction2
    complex         (c_double_complex          )                                                          :: normalization

#ifdef FFTW3UNAVAIL
    call Error_Report('FFTW3 library is required but was not found'//{introspection:location})
#endif
    if (.not.associated(self%randomNumberGenerator_)) call Error_Report('no random number generator supplied'//{introspection:location})
    ! Initialize geometry if necessary.
    if (.not.self%geometryInitialized) then
       call self%randomsInitialize()
       self%geometryInitialized=.true.
    end if
    ! Find the comoving distance corresponding to this distance module.
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
    ! Loop over randoms.
    do i=1,size(self%randomPhi)
       ! Choose random distances.
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
       ! Convert to Cartesian coordinates.
       position1=distance1*[sin(self%randomTheta(i))*cos(self%randomPhi(i)),sin(self%randomTheta(i))*sin(self%randomPhi(i)),cos(self%randomTheta(i))]+origin
       position2=distance2*[sin(self%randomTheta(i))*cos(self%randomPhi(i)),sin(self%randomTheta(i))*sin(self%randomPhi(i)),cos(self%randomTheta(i))]+origin
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
  end subroutine randomPointsWindowFunctions

  double precision function randomPointsAngularPower(self,i,j,l)
    !!{
    Angular power is not available, so simply aborts.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (surveyGeometryRandomPoints), intent(inout) :: self
    integer                            , intent(in   ) :: i   , j, l
    !$GLC attributes unused :: self, i, j, l

    randomPointsAngularPower=0.0d0
    call Error_Report('angular power is not available'//{introspection:location})
    return
  end function randomPointsAngularPower

  logical function randomPointIncluded(self,point,mass)
    !!{
    Return true if a point is included in the survey geometry.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (surveyGeometryRandomPoints), intent(inout)               :: self
    double precision                            , intent(in   ), dimension(3) :: point
    double precision                            , intent(in   )               :: mass
    !$GLC attributes unused :: self, point, mass

    randomPointIncluded=.false.
    call Error_Report('point inclusion is not supported for window functions defined by random points'//{introspection:location})
    return
  end function randomPointIncluded

