!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Implements survey geometries over the full sky.
  
  !# <surveyGeometry name="surveyGeometryFullSky">
  !#  <description>Implements survey geometries over the full sky.</description>
  !# </surveyGeometry>

  type, extends(surveyGeometryClass) :: surveyGeometryFullSky
     double precision :: limitDistanceMinimum, limitDistanceMaximum
   contains
     procedure :: distanceMinimum         => fullSkyDistanceMinimum
     procedure :: distanceMaximum         => fullSkyDistanceMaximum
     procedure :: solidAngle              => fullSkySolidAngle
     procedure :: windowFunctionAvailable => fullSkyWindowFunctionAvailable
     procedure :: angularPowerAvailable   => fullSkyAngularPowerAvailable
     procedure :: windowFunctions         => fullSkyWindowFunctions
     procedure :: angularPower            => fullSkyAngularPower
     procedure :: pointIncluded           => fullSkyPointIncluded
  end type surveyGeometryFullSky

  interface surveyGeometryFullSky
     !% Constructors for the full sky survey geometry class.
     module procedure fullSkyDefaultConstructor
     module procedure fullSkyRedshiftConstructor
  end interface surveyGeometryFullSky

contains

  function fullSkyDefaultConstructor()
    !% Default constructor for the full sky survey geometry.
    use Galacticus_Error
    implicit none
    type(surveyGeometryFullSky) :: fullSkyDefaultConstructor

    fullSkyDefaultConstructor%limitDistanceMinimum=-1.0d0
    fullSkyDefaultConstructor%limitDistanceMaximum=-1.0d0
    call Galacticus_Error_Report('fullSkyDefaultConstructor','default constructor not available')
    return
  end function fullSkyDefaultConstructor

  function fullSkyRedshiftConstructor(redshiftMinimum,redshiftMaximum)
    !% Constructor for the full sky survey class which allows specification of minimum and maximum redshifts.
    use Cosmology_Functions
    use Cosmology_Functions_Options
    implicit none
    type            (surveyGeometryFullSky  )                :: fullSkyRedshiftConstructor
    double precision                         , intent(in   ) :: redshiftMinimum           , redshiftMaximum
    class           (cosmologyFunctionsClass), pointer       :: cosmologyFunctions_
    
    cosmologyFunctions_                             => cosmologyFunctions()     
    fullSkyRedshiftConstructor%limitDistanceMinimum =  cosmologyFunctions_%distanceComovingConvert(                               &
         &                                                                                         output  =distanceTypeComoving, &
         &                                                                                         redshift=redshiftMinimum       &
         &                                                                                        )
    fullSkyRedshiftConstructor%limitDistanceMaximum =  cosmologyFunctions_%distanceComovingConvert(                               &
         &                                                                                         output  =distanceTypeComoving, &
         &                                                                                         redshift=redshiftMaximum       &
         &                                                                                        )
    return
  end function fullSkyRedshiftConstructor

  double precision function fullSkyDistanceMinimum(self,mass,field)
    !% Compute the minimum distance at which a galaxy is visible.
    implicit none
    class           (surveyGeometryFullSky), intent(inout)           :: self
    double precision                       , intent(in   )           :: mass
    integer                                , intent(in   ), optional :: field
    !GCC$ attributes unused :: mass, field
    
    fullSkyDistanceMinimum=self%limitDistanceMinimum
    return
  end function fullSkyDistanceMinimum

  double precision function fullSkyDistanceMaximum(self,mass,field)
    !% Compute the maximum distance at which a galaxy is visible.
    implicit none
    class           (surveyGeometryFullSky), intent(inout)           :: self
    double precision                       , intent(in   )           :: mass
    integer                                , intent(in   ), optional :: field
    !GCC$ attributes unused :: mass, field

    fullSkyDistanceMaximum=self%limitDistanceMaximum
    return
  end function fullSkyDistanceMaximum

  double precision function fullSkySolidAngle(self,field)
    !% Return the solid angle of the \cite{li_distribution_2009} sample.
    use Numerical_Constants_Math
    implicit none
    class           (surveyGeometryFullSky), intent(inout)           :: self
    integer                                , intent(in   ), optional :: field
    !GCC$ attributes unused :: self, field

    fullSkySolidAngle=4.0d0*Pi
    return
  end function fullSkySolidAngle

  logical function fullSkyWindowFunctionAvailable(self)
    !% Return true to indicate that survey window function is available.
    implicit none
    class(surveyGeometryFullSky), intent(inout) :: self
    !GCC$ attributes unused :: self

    fullSkyWindowFunctionAvailable=.true.
    return
  end function fullSkyWindowFunctionAvailable

  logical function fullSkyAngularPowerAvailable(self)
    !% Return true to indicate that survey angular power is available.
    implicit none
    class(surveyGeometryFullSky), intent(inout) :: self
    !GCC$ attributes unused :: self

    fullSkyAngularPowerAvailable=.true.
    return
  end function fullSkyAngularPowerAvailable

  subroutine fullSkyWindowFunctions(self,mass1,mass2,gridCount,boxLength,windowFunction1,windowFunction2)
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
    class           (surveyGeometryFullSky), intent(inout)                                           :: self
    double precision                       , intent(in   )                                           :: mass1                   , mass2
    integer                                , intent(in   )                                           :: gridCount
    double precision                       , intent(  out)                                           :: boxLength
    complex         (c_double_complex     ), intent(  out), dimension(gridCount,gridCount,gridCount) :: windowFunction1         , windowFunction2
    complex         (c_double_complex     ),                dimension(gridCount,gridCount,gridCount) :: selectionFunction1      , selectionFunction2
    double precision                       ,                dimension(3                            ) :: origin                  , position
    integer                                                                                          :: i                       , j                       , &
         &                                                                                              k
    double precision                                                                                 :: comovingDistanceMaximum1, comovingDistanceMaximum2, &
         &                                                                                              comovingDistanceMinimum1, comovingDistanceMinimum2, &
         &                                                                                              distance
    type            (c_ptr                )                                                          :: plan
    complex         (c_double_complex     )                                                          :: normalization
    
    ! Find the comoving distance corresponding to this distance module.
    comovingDistanceMaximum1=self%distanceMaximum(mass1)
    comovingDistanceMaximum2=self%distanceMaximum(mass2)
    comovingDistanceMinimum1=self%distanceMinimum(mass1)
    comovingDistanceMinimum2=self%distanceMinimum(mass2)
    ! Determine a suitable box length for the window function calculation.
    boxLength=3.0d0*max(comovingDistanceMaximum1-comovingDistanceMinimum1,comovingDistanceMaximum2-comovingDistanceMinimum2)
    ! Set up origin.
    origin=[0.5d0,0.5d0,0.5d0]*boxLength
    ! Populate the cube with the survey selection function.
    selectionFunction1=dcmplx(0.0d0,0.0d0)
    selectionFunction2=dcmplx(0.0d0,0.0d0)
    ! Iterate over grid cells.
    do i=1,gridCount
       do j=1,gridCount
          do k=1,gridCount
             ! Compute distance to this point.
             position=([dble(i),dble(j),dble(k)]-[0.5d0,0.5d0,0.5d0])*boxLength-origin
             distance=Vector_Magnitude(position)
             ! Apply to the mesh.
             if (distance > comovingDistanceMinimum1 .and. distance < comovingDistanceMaximum1) &
                  & call Meshes_Apply_Point(selectionFunction1,boxLength,position,pointWeight=cmplx(1.0d0,0.0d0,kind=c_double_complex),cloudType=cloudTypeTriangular)
             if (distance > comovingDistanceMinimum2 .and. distance < comovingDistanceMaximum2) &
                  & call Meshes_Apply_Point(selectionFunction2,boxLength,position,pointWeight=cmplx(1.0d0,0.0d0,kind=c_double_complex),cloudType=cloudTypeTriangular)
          end do
       end do
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
  end subroutine fullSkyWindowFunctions

  double precision function fullSkyAngularPower(self,i,j,l)
    !% Return the angular power for the full sky.
    use NUmerical_Constants_Math
    implicit none
    class  (surveyGeometryFullSky), intent(inout) :: self
    integer                       , intent(in   ) :: i   , j, l
    !GCC$ attributes unused :: self, i, j

    if (l == 0) then
       fullSkyAngularPower=4.0d0*Pi
    else
       fullSkyAngularPower=0.0d0
    end if
    return
  end function fullSkyAngularPower

  logical function fullSkyPointIncluded(self,point,mass)
    !% Return true if a point is included in the survey geometry.
    use Vectors
    implicit none
    class           (surveyGeometryFullSky), intent(inout)               :: self
    double precision                       , intent(in   ), dimension(3) :: point
    double precision                       , intent(in   )               :: mass
    double precision                                                     :: distance
    !GCC$ attributes unused :: mass

    distance            =Vector_Magnitude(point)
    fullSkyPointIncluded=(distance >= self%limitDistanceMinimum .and. distance <= self%limitDistanceMaximum)
    return
  end function fullSkyPointIncluded

