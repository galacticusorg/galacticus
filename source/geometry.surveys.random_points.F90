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

!% Implements survey geometries defined by random points.
  
  !# <surveyGeometry name="surveyGeometryRandomPoints">
  !#  <description>Implements survey geometries defined by random points.</description>
  !#  <abstract>yes</abstract>
  !# </surveyGeometry>

  type, abstract, extends(surveyGeometryClass) :: surveyGeometryRandomPoints
     logical                                     :: geometryInitialized
     double precision, allocatable, dimension(:) :: randomTheta        , randomPhi
   contains
     !@ <objectMethods>
     !@   <object>surveyGeometryRandomPoints</object>
     !@   <objectMethod>
     !@     <method>randomsInitialize</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Initialize arrays of random points to define the survey angular geometry.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                                          :: windowFunctionAvailable => randomPointsWindowFunctionAvailable
     procedure                                          :: angularPowerAvailable   => randomPointsAngularPowerAvailable
     procedure                                          :: windowFunctions         => randomPointsWindowFunctions
     procedure                                          :: angularPower            => randomPointsAngularPower
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
    !% Return true to indicate that survey window function is available.
    implicit none
    class(surveyGeometryRandomPoints), intent(inout) :: self

    randomPointsWindowFunctionAvailable=.true.
    return
  end function randomPointsWindowFunctionAvailable

  logical function randomPointsAngularPowerAvailable(self)
    !% Return false to indicate that survey angular power is not available.
    implicit none
    class(surveyGeometryRandomPoints), intent(inout) :: self

    randomPointsAngularPowerAvailable=.false.
    return
  end function randomPointsAngularPowerAvailable

  subroutine randomPointsWindowFunctions(self,mass1,mass2,gridCount,boxLength,windowFunction1,windowFunction2)
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
    use String_Handling
    use Galacticus_Error
    implicit none
    class           (surveyGeometryRandomPoints), intent(inout)                                           :: self
    double precision                            , intent(in   )                                           :: mass1                   , mass2
    integer                                     , intent(in   )                                           :: gridCount
    double precision                            , intent(  out)                                           :: boxLength
    complex         (c_double_complex          ), intent(  out), dimension(gridCount,gridCount,gridCount) :: windowFunction1         , windowFunction2
    double precision                            ,                dimension(3                            ) :: origin                  , position1              , &
         &                                                                                                   position2
    type            (fgsl_rng                  ), save                                                    :: pseudoSequenceObject
    type            (varying_string)                                                                      :: message
    logical                                     , save                                                    :: reset=.true.
    integer                                                                                               :: i                       , j
    double precision                                                                                      :: comovingDistanceMaximum1, comovingDistanceMaximum2, &
         &                                                                                                   comovingDistanceMinimum1, comovingDistanceMinimum2, &
         &                                                                                                   distance1               , distance2
    type            (c_ptr                     )                                                          :: plan
    complex         (c_double_complex          ),                dimension(gridCount,gridCount,gridCount) :: selectionFunction1      , selectionFunction2
    complex         (c_double_complex          )                                                          :: normalization
    
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
    selectionFunction1=cmplx(0.0d0,0.0d0)
    selectionFunction2=cmplx(0.0d0,0.0d0)
    ! Loop over randoms.
    do i=1,size(self%randomPhi)
       ! Choose random distances.
       distance1=+(                                               &
            &       Pseudo_Random_Get(pseudoSequenceObject,reset) &
            &      *(                                             &
            &        +comovingDistanceMaximum1**3                 &
            &        -comovingDistanceMinimum1**3                 &
            &       )                                             &
            &        +comovingDistanceMinimum1**3                 &
            &     )**(1.0d0/3.0d0)                                &
            &    -min(                                            &
            &         comovingDistanceMinimum1   ,                &
            &         comovingDistanceMinimum2                    &
            &        )
       distance2=+(                                               &
            &       Pseudo_Random_Get(pseudoSequenceObject,reset) &
            &      *(                                             &
            &        +comovingDistanceMaximum2**3                 &
            &        -comovingDistanceMinimum2**3                 &
            &       )                                             &
            &        +comovingDistanceMinimum2**3                 &
            &     )**(1.0d0/3.0d0)                                &
            &    -min(                                            &
            &         comovingDistanceMinimum1   ,                &
            &         comovingDistanceMinimum2                    &
            &        )
       ! Convert to Cartesian coordinates.
       position1=distance1*[sin(self%randomTheta(i))*cos(self%randomPhi(i)),sin(self%randomTheta(i))*sin(self%randomPhi(i)),cos(self%randomTheta(i))]+origin
       position2=distance2*[sin(self%randomTheta(i))*cos(self%randomPhi(i)),sin(self%randomTheta(i))*sin(self%randomPhi(i)),cos(self%randomTheta(i))]+origin
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
  end subroutine randomPointsWindowFunctions

  double precision function randomPointsAngularPower(self,i,j,l)
    !% Angular power is not available, so simply aborts.
    use Galacticus_Error
    implicit none
    class  (surveyGeometryRandomPoints), intent(inout) :: self
    integer                            , intent(in   ) :: i   , j, l

    call Galacticus_Error_Report('randomPointsAngularPower','angular power is not available')
    return
  end function randomPointsAngularPower
