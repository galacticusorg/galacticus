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
Implements survey geometries over the full sky.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <surveyGeometry name="surveyGeometryFullSky">
   <description>Implements survey geometries over the full sky.</description>
  </surveyGeometry>
  !!]
  type, extends(surveyGeometryClass) :: surveyGeometryFullSky
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_  => null()
     double precision                                   :: limitDistanceMinimum          , limitDistanceMaximum, &
          &                                                redshiftMinimum               , redshiftMaximum
   contains
     final     ::                            fullSkyDestructor
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
     !!{
     Constructors for the full sky survey geometry class.
     !!}
     module procedure fullSkyConstructorParameters
     module procedure fullSkyConstructorInternal
  end interface surveyGeometryFullSky

contains

  function fullSkyConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the full sky survey geometry.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (surveyGeometryFullSky  )                :: self
    type            (inputParameters        ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass), pointer       :: cosmologyFunctions_
    double precision                                         :: redshiftMinimum    , redshiftMaximum

    ! Check and read parameters.
    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <inputParameter>
      <name>redshiftMinimum</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The minimum redshift for the survey.</description>
    </inputParameter>
    <inputParameter>
      <name>redshiftMaximum</name>
      <defaultValue>huge(1.0d0)</defaultValue>
      <source>parameters</source>
      <description>The maximum redshift for the survey.</description>
    </inputParameter>
    !!]
    ! Build the object.
    self=surveyGeometryFullSky(                                         &
         &                     redshiftMinimum    =redshiftMinimum    , &
         &                     redshiftMaximum    =redshiftMaximum    , &
         &                     cosmologyFunctions_=cosmologyFunctions_  &
         &                    )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function fullSkyConstructorParameters

  function fullSkyConstructorInternal(redshiftMinimum,redshiftMaximum,distanceMinimum,distanceMaximum,cosmologyFunctions_) result(self)
    !!{
    Constructor for the full sky survey class which allows specification of minimum and maximum redshifts.
    !!}
    use :: Error                      , only : Error_Report
    use :: Cosmology_Functions_Options, only : distanceTypeComoving
    implicit none
    type            (surveyGeometryFullSky  )                                  :: self
    double precision                         , intent(in   ), optional         :: redshiftMinimum    , redshiftMaximum, &
         &                                                                        distanceMinimum    , distanceMaximum
    class           (cosmologyFunctionsClass), intent(in   ), optional, target :: cosmologyFunctions_
    !![
    <constructorAssign variables="*cosmologyFunctions_"/>
    !!]

    if (present(redshiftMinimum)) then
       if (     present(distanceMinimum    )) call Error_Report('ambiguous minimum distance'         //{introspection:location})
       if (.not.present(cosmologyFunctions_)) call Error_Report('cosmology function must be supplied'//{introspection:location})
       self   %limitDistanceMinimum=self%cosmologyFunctions_%distanceComovingConvert(                               &
            &                                                                        output  =distanceTypeComoving, &
            &                                                                        redshift=redshiftMinimum       &
            &                                                                       )
    else if (present(distanceMinimum)) then
       self   %limitDistanceMinimum=distanceMinimum
    else
       call Error_Report('no minimum distance was specified'//{introspection:location})
    end if
    if (present(redshiftMinimum)) then
       if (present(distanceMaximum)) call Error_Report('ambiguous maximum distance'//{introspection:location})
       if (redshiftMaximum < huge(1.0d0)) then
       if (.not.present(cosmologyFunctions_)) call Error_Report('cosmology function must be supplied'//{introspection:location})
          self%limitDistanceMaximum=self%cosmologyFunctions_%distanceComovingConvert(                               &
               &                                                                     output  =distanceTypeComoving, &
               &                                                                     redshift=redshiftMaximum       &
               &                                                                    )
       else
          self%limitDistanceMaximum=huge(1.0d0)
       end if
    else if (present(distanceMaximum)) then
       self   %limitDistanceMaximum=distanceMaximum
    else
       call Error_Report('no maximum distance was specified'//{introspection:location})
    end if
    return
  end function fullSkyConstructorInternal

  subroutine fullSkyDestructor(self)
    !!{
    Destructor for the \refClass{surveyGeometryFullSky} survey geometry class.
    !!}
    implicit none
    type(surveyGeometryFullSky), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine fullSkyDestructor

  double precision function fullSkyDistanceMinimum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the minimum distance at which a galaxy is visible.
    !!}
    implicit none
    class           (surveyGeometryFullSky), intent(inout)           :: self
    double precision                       , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                              luminosity, starFormationRate
    integer                                , intent(in   ), optional :: field
    !$GLC attributes unused :: mass, field, magnitudeAbsolute, luminosity, starFormationRate

    fullSkyDistanceMinimum=self%limitDistanceMinimum
    return
  end function fullSkyDistanceMinimum

  double precision function fullSkyDistanceMaximum(self,mass,magnitudeAbsolute,luminosity,starFormationRate,field)
    !!{
    Compute the maximum distance at which a galaxy is visible.
    !!}
    implicit none
    class           (surveyGeometryFullSky), intent(inout)           :: self
    double precision                       , intent(in   ), optional :: mass      , magnitudeAbsolute, &
         &                                                              luminosity, starFormationRate
    integer                                , intent(in   ), optional :: field
    !$GLC attributes unused :: mass, magnitudeAbsolute, field, luminosity, starFormationRate

    fullSkyDistanceMaximum=self%limitDistanceMaximum
    return
  end function fullSkyDistanceMaximum

  double precision function fullSkySolidAngle(self,field)
    !!{
    Return the solid angle of the \cite{li_distribution_2009} sample.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (surveyGeometryFullSky), intent(inout)           :: self
    integer                                , intent(in   ), optional :: field
    !$GLC attributes unused :: self, field

    fullSkySolidAngle=4.0d0*Pi
    return
  end function fullSkySolidAngle

  logical function fullSkyWindowFunctionAvailable(self)
    !!{
    Return true to indicate that survey window function is available.
    !!}
    implicit none
    class(surveyGeometryFullSky), intent(inout) :: self
    !$GLC attributes unused :: self

    fullSkyWindowFunctionAvailable=.true.
    return
  end function fullSkyWindowFunctionAvailable

  logical function fullSkyAngularPowerAvailable(self)
    !!{
    Return true to indicate that survey angular power is available.
    !!}
    implicit none
    class(surveyGeometryFullSky), intent(inout) :: self
    !$GLC attributes unused :: self

    fullSkyAngularPowerAvailable=.true.
    return
  end function fullSkyAngularPowerAvailable

  subroutine fullSkyWindowFunctions(self,mass1,mass2,gridCount,boxLength,windowFunction1,windowFunction2)
    !!{
    Compute the window function for the survey.
    !!}
#ifdef FFTW3AVAIL
    use            :: FFTW3        , only : fftw_plan_dft_3d  , FFTW_FORWARD       , FFTW_ESTIMATE, fftw_execute_dft, &
         &                                  fftw_destroy_plan
#else
    use            :: Error        , only : Error_Report
#endif
    use, intrinsic :: ISO_C_Binding, only : c_ptr             , c_double_complex
    use            :: Meshes       , only : Meshes_Apply_Point, cloudTypeTriangular
    use            :: Vectors      , only : Vector_Magnitude
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
#ifdef FFTW3AVAIL
    type            (c_ptr                )                                                          :: plan
#endif
    complex         (c_double_complex     )                                                          :: normalization

#ifdef FFTW3UNAVAIL
    call Error_Report('FFTW3 library is required but was not found'//{introspection:location})
#endif
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
  end subroutine fullSkyWindowFunctions

  double precision function fullSkyAngularPower(self,i,j,l)
    !!{
    Return the angular power for the full sky.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class  (surveyGeometryFullSky), intent(inout) :: self
    integer                       , intent(in   ) :: i   , j, l
    !$GLC attributes unused :: self, i, j

    if (l == 0) then
       fullSkyAngularPower=4.0d0*Pi
    else
       fullSkyAngularPower=0.0d0
    end if
    return
  end function fullSkyAngularPower

  logical function fullSkyPointIncluded(self,point,mass)
    !!{
    Return true if a point is included in the survey geometry.
    !!}
    use :: Vectors, only : Vector_Magnitude
    implicit none
    class           (surveyGeometryFullSky), intent(inout)               :: self
    double precision                       , intent(in   ), dimension(3) :: point
    double precision                       , intent(in   )               :: mass
    double precision                                                     :: distance
    !$GLC attributes unused :: mass

    distance            =Vector_Magnitude(point)
    fullSkyPointIncluded=(distance >= self%limitDistanceMinimum .and. distance <= self%limitDistanceMaximum)
    return
  end function fullSkyPointIncluded

