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
  Implementation of a Gaussian ellipsoid mass distribution class.
  !!}

  use :: Linear_Algebra, only : matrix
  
  !![
  <massDistribution name="massDistributionGaussianEllipsoid">
   <description>A mass distribution class for Gaussian ellipsoids.</description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionClass) :: massDistributionGaussianEllipsoid
     !!{
     A Gaussian ellipsoid mass distribution. The formulation and calculations follow the parameterizations and approach of
     \cite[][chapter 3]{chandrasekhar_ellipsoidal_1987}. The ellipsoidal has scale lengths $a_\mathrm{i}$ along each of the
     three perpendicular axes (which are assumed to be aligned with the Cartesian $(x,y,z)$ axes).
     !!}
     double precision                     , dimension(3          ) :: scaleLength
     double precision                                              :: mass                             , density_
     logical                                                       :: accelerationInitialized
     double precision        , allocatable, dimension(:          ) :: accelerationX                    , accelerationScaleLength
     double precision        , allocatable, dimension(:,:,:,:,:,:) :: accelerationVector
     double precision                                              :: accelerationXMinimumLog          , accelerationXMaximumLog               , &
          &                                                           accelerationScaleLengthMinimumLog, accelerationScaleLengthMaximumLog     , &
          &                                                           accelerationXInverseInterval     , accelerationScaleLengthInverseInterval, &
          &                                                           scaleLengthMaximum
     integer                              , dimension(3          ) :: axesMapIn                        , axesMapOut
     double precision                     , dimension(2          ) :: axisRatio
     type            (matrix), allocatable                         :: rotationIn                       , rotationOut
     double precision                     , dimension(3)           :: axis1                            , axis2                                 , &
          &                                                           axis3
   contains
     !![
     <methods>
       <method description="Compute the density on the isodensity surface defined by the parameter $m^2$2." method="densityEllipsoidal"     />
       <method description="Tabulate the gravitational acceleration due to the ellipsoid."                  method="accelerationTabulate"   />
       <method description="Interpolate in the tabulated gravitational acceleration due to the ellipsoid."  method="accelerationInterpolate"/>
       <method description="(Re)initialize the structural properties of the Gaussian ellispoid."            method="initialize"             />
     </methods>
     !!]
     procedure :: density                 => gaussianEllipsoidDensity
     procedure :: densityEllipsoidal      => gaussianEllipsoidDensityEllipsoidal
     procedure :: acceleration            => gaussianEllipsoidAcceleration
     procedure :: accelerationTabulate    => gaussianEllipsoidAccelerationTabulate
     procedure :: accelerationInterpolate => gaussianEllipsoidAccelerationInterpolate
     procedure :: initialize              => gaussianEllipsoidInitialize
  end type massDistributionGaussianEllipsoid

  interface massDistributionGaussianEllipsoid
     !!{
     Constructors for the {\normalfont \ttfamily gaussianEllipsoid} mass distribution class.
     !!}
     module procedure gaussianEllipsoidConstructorParameters
     module procedure gaussianEllipsoidConstructorInternal
  end interface massDistributionGaussianEllipsoid

contains

  function gaussianEllipsoidConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily gaussianEllipsoid} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Linear_Algebra            , only : assignment(=)                 , vector
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionGaussianEllipsoid)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    double precision                                   , dimension(3)  :: scaleLength  , axis1, &
         &                                                                axis2        , axis3
    double precision                                                   :: mass
    logical                                                            :: dimensionless
    type            (vector                           ), dimension(3)  :: axes
    type            (varying_string                   )                :: componentType
    type            (varying_string                   )                :: massType

    !![
    <inputParameter>
      <name>mass</name>
      <description>The mass of the ellipsoid.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>scaleLength</name>
      <description>The scale lengths of the ellipsoid along each axis.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>axis1</name>
      <defaultValue>[1.0d0,0.0d0,0.0d0]</defaultValue>
      <description>The unit vector defining the first axis of the ellipsoid.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>axis2</name>
      <defaultValue>[0.0d0,1.0d0,0.0d0]</defaultValue>
      <description>The unit vector defining the second axis of the ellipsoid.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>axis3</name>
      <defaultValue>[0.0d0,0.0d0,1.0d0]</defaultValue>
      <description>The unit vector defining the third axis of the ellipsoid.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>If true the Gaussian ellipsoid profile is considered to be dimensionless.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>componentType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The component type that this mass distribution represents.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The mass type that this mass distribution represents.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    axes(1)=axis1
    axes(2)=axis2
    axes(3)=axis3    
    !![
    <conditionalCall>
     <call>self=massDistributionGaussianEllipsoid(scaleLength=scaleLength,axes=axes,componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})</call>
     <argument name="mass"          value="mass"          parameterPresent="parameters"/>
     <argument name="dimensionless" value="dimensionless" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function gaussianEllipsoidConstructorParameters
  
  function gaussianEllipsoidConstructorInternal(scaleLength,axes,rotation,mass,dimensionless,componentType,massType) result(self)
    !!{
    Constructor for ``gaussianEllipsoid'' convergence class.
    !!}
    use :: Error               , only : Error_Report
    use :: Linear_Algebra      , only : vector       , assignment(=)
    use :: Numerical_Comparison, only : Values_Differ
    implicit none
    type            (massDistributionGaussianEllipsoid)                                        :: self
    double precision                                   , intent(in   ), dimension(3)           :: scaleLength
    type            (vector                           ), intent(in   ), dimension(3), optional :: axes
    type            (matrix                           ), intent(in   )              , optional :: rotation
    double precision                                   , intent(in   )              , optional :: mass
    logical                                            , intent(in   )              , optional :: dimensionless
    type            (enumerationComponentTypeType     ), intent(in   )              , optional :: componentType
    type            (enumerationMassTypeType          ), intent(in   )              , optional :: massType
    !![
    <constructorAssign variables="mass, dimensionless, componentType, massType"/>
    !!]

    ! Determine if profile is dimensionless.
    if (present(dimensionless)) then
       self%dimensionless=dimensionless
    else
       self%dimensionless=.false.
    end if
    ! If dimensionless, then maximum mass should be 1.0.
    if (self%dimensionless) then
       if (present(mass)) then
          if (Values_Differ(mass,1.0d0,absTol=1.0d-6)) call Error_Report('mass should be unity for a dimensionless profile (or simply do not specify a mass)'//{introspection:location})
       end if
       self%mass=1.0d0
    end if
    ! Initialize structural properties.
    call self%initialize(scaleLength,axes,rotation)
    ! Set acceleration as uninitialized.
    self%accelerationInitialized=.false.
    return
  end function gaussianEllipsoidConstructorInternal

  subroutine gaussianEllipsoidInitialize(self,scaleLength,axes,rotation)
    !!{
    (Re)initialize properties of a Gaussian ellipsoid mass distribution.
    !!}
    use :: Error                   , only : Error_Report
    use :: Linear_Algebra          , only : assignment(=), matrixRotation, vector
    use :: Numerical_Comparison    , only : Values_Differ
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionGaussianEllipsoid), intent(inout)                         :: self
    double precision                                   , intent(in   ), dimension(3)           :: scaleLength
    type            (vector                           ), intent(in   ), dimension(3), optional :: axes
    type            (matrix                           ), intent(in   )              , optional :: rotation
    type            (vector                           )               , dimension(3)           :: axesPrinciple
    integer                                                                                    :: i

    ! If dimensionless, then maximum scale length should be 1.0.
    if (self%dimensionless) then
       if (Values_Differ(maxval(scaleLength),1.0d0,absTol=1.0d-6)) call Error_Report('maximum scaleRadius should be unity for a dimensionless profile'                   //{introspection:location})
    end if
    ! Assign scalelengths.
    self%scaleLength=scaleLength
    ! Assert that axis vectors are unit vectors and are orthogonal.
    if (present(axes)) then
       do i=1,3
          if (Values_Differ(axes(i)%magnitude(),1.0d0,absTol=1.0d-6)) call Error_Report('axis vectors must be unit vectors'//{introspection:location})
       end do
       if (Values_Differ(axes(1).dot.axes(2),0.0d0,absTol=1.0d-6)) call Error_Report('axis vectors must be orthogonal'//{introspection:location})
       if (Values_Differ(axes(1).dot.axes(3),0.0d0,absTol=1.0d-6)) call Error_Report('axis vectors must be orthogonal'//{introspection:location})
       if (Values_Differ(axes(2).dot.axes(3),0.0d0,absTol=1.0d-6)) call Error_Report('axis vectors must be orthogonal'//{introspection:location})
    end if
    ! Compute the central density.
    self%density_=+        self%mass         &
         &        /product(self%scaleLength) &
         &        /(                         &
         &          +2.0d0                   &
         &          *Pi                      &
         &         )**1.5d0
    ! Compute the scaling and axes maps that convert to an ellipsoid where the largest semi-axis is a₃=1.
    self%scaleLengthMaximum=maxval(                 scaleLength         )
    self%axesMapIn         =mod   ([0,1,2]+3-maxloc(scaleLength,dim=1),3)+1
    do i=1,3
       self%axesMapOut(self%axesMapIn(i))=i
    end do
    do i=1,2
       self%axisRatio(i)=+self%scaleLength(self%axesMapOut(i)) &
            &            /self%scaleLengthMaximum
    end do
    ! Compute rotation matrices required to rotate the ellipsoid to be aligned with the principle Cartesian axes, and back again.
    if (allocated(self%rotationIn )) deallocate(self%rotationIn )
    if (allocated(self%rotationOut)) deallocate(self%rotationOut)
    allocate(self%rotationIn )
    allocate(self%rotationOut)
    if (present(axes)) then
       self%axis1      =axes(1)
       self%axis2      =axes(2)
       self%axis3      =axes(3)
       axesPrinciple(1)=vector([1.0d0,0.0d0,0.0d0])
       axesPrinciple(2)=vector([0.0d0,1.0d0,0.0d0])
       axesPrinciple(3)=vector([0.0d0,0.0d0,1.0d0])
       self%rotationIn =matrixRotation(axes,axesPrinciple)
    else if (present(rotation)) then
       self%rotationIn=matrix(rotation)
    else
       call Error_Report('either principle axes or a rotation matrix must be supplied'//{introspection:location})
    end if
    self%rotationOut=self%rotationIn%inverse()
    return
  end subroutine gaussianEllipsoidInitialize
  
  double precision function gaussianEllipsoidDensity(self,coordinates)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a Gaussian ellipsoid mass distribution.
    !!}
    use :: Coordinates   , only : assignment(=), coordinateCartesian
    use :: Linear_Algebra, only : assignment(=), operator(*)        , vector
    implicit none
    class           (massDistributionGaussianEllipsoid), intent(inout) :: self
    class           (coordinate                       ), intent(in   ) :: coordinates
    type            (coordinateCartesian              )                :: position
    double precision                                   , dimension(3)  :: positionComponents
    double precision                                                   :: mSquared
    type           (vector                            )                :: positionVectorUnrotated, positionVector

    ! Rotate the position to the frame where the ellipsoid is aligned with the principle Cartesian axes.
    position                = coordinates
    positionComponents      = position
    positionVectorUnrotated = positionComponents
    positionVector          = self%rotationIn         &
         &                   *positionVectorUnrotated
    positionComponents      = positionVector
    position                = positionComponents
    ! Evaluate the density.
    mSquared                =+(position%x()/self%scaleLength(1))**2 &
         &                   +(position%y()/self%scaleLength(2))**2 &
         &                   +(position%z()/self%scaleLength(3))**2
    gaussianEllipsoidDensity=+self%densityEllipsoidal(mSquared)
    return
  end function gaussianEllipsoidDensity

  double precision function gaussianEllipsoidDensityEllipsoidal(self,mSquared)
    !!{
    Return the density on the isodensity surface defined by $m^2 = \sum_{\mathrm{i}=1}^3 (x_\mathrm{i}/a_\mathrm{i})^2$.
    !!}
    use :: Coordinates, only : assignment(=), coordinateCartesian
    implicit none
    class           (massDistributionGaussianEllipsoid), intent(inout) :: self
    double precision                                   , intent(in   ) :: mSquared

    gaussianEllipsoidDensityEllipsoidal=+self%density_ &
         &                              *exp(          &
         &                                   -0.5d0    &
         &                                   *mSquared &
         &                                  )
    return
  end function gaussianEllipsoidDensityEllipsoidal

  function gaussianEllipsoidAcceleration(self,coordinates)
    !!{
    Computes the gravitational acceleration at {\normalfont \ttfamily coordinates} for Gaussian ellipsoid mass distributions.
    !!}
    use :: Coordinates                     , only : assignment(=)                 , coordinateCartesian
    use :: Linear_Algebra                  , only : assignment(=)                 , operator(*)        , vector
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    double precision                                   , dimension(3)  :: gaussianEllipsoidAcceleration
    class           (massDistributionGaussianEllipsoid), intent(inout) :: self
    class           (coordinate                       ), intent(in   ) :: coordinates
    type            (coordinateCartesian              )                :: coordinatesCartesian
    double precision                                   , dimension(3)  :: positionCartesian            , positionCartesianScaleFree , &
         &                                                                accelerationScaleFree
    integer                                                            :: i
    type            (vector                            )               :: positionVector               , positionVectorUnrotated    , &
         &                                                                accelerationVector           , accelerationVectorUnrotated
    
    ! Ensure that acceleration is tabulated.
    call self%accelerationTabulate()
    ! Construct the scale-free (and rotated) position.
    coordinatesCartesian    = coordinates
    positionCartesian       = coordinatesCartesian
    positionVectorUnrotated = positionCartesian
    positionVector          = self%rotationIn         &
         &                   *positionVectorUnrotated
    positionCartesian       = positionVector
    do i=1,3
       positionCartesianScaleFree(self%axesMapIn(i))=positionCartesian(i)
    end do
    positionCartesianScaleFree=+positionCartesianScaleFree &
         &                     /self%scaleLengthMaximum
    ! Interpolate to get the scale-free acceleration.
    accelerationScaleFree=self%accelerationInterpolate(positionCartesianScaleFree)
    ! De-rotate and scale the acceleration.
    do i=1,3
       gaussianEllipsoidAcceleration(self%axesMapOut(i))=accelerationScaleFree(i)
    end do
    accelerationVector           = gaussianEllipsoidAcceleration
    accelerationVectorUnrotated  = self%rotationOut              &
         &                        *accelerationVector
    gaussianEllipsoidAcceleration= accelerationVectorUnrotated
    if (.not.self%isDimensionless())                                     &
         & gaussianEllipsoidAcceleration=+gaussianEllipsoidAcceleration  &
         &                               *gravitationalConstant_internal &
         &                               *self%mass                      &
         &                               /self%scaleLengthMaximum**2
    return
  end function gaussianEllipsoidAcceleration
  
  function gaussianEllipsoidAccelerationInterpolate(self,positionCartesian)
    !!{
    Interpolate gravitational acceleration in the tabulated solutions for a Gaussian ellipsoid.
    !!}
    use :: Trigonometric_Functions, only : hypotenuse
    implicit none
    double precision                                   , dimension(    3)                :: gaussianEllipsoidAccelerationInterpolate
    class           (massDistributionGaussianEllipsoid)                  , intent(inout) :: self
    double precision                                   , dimension(    3), intent(in   ) :: positionCartesian
    double precision                                   , dimension(0:1,5)                :: h
    integer                                            , dimension(    5)                :: i
    double precision                                                                     :: positionCartesianLog                    , axisRatioLog, &
         &                                                                                  radiusSpherical
    integer                                                                              :: j1                                      , j2          , &
         &                                                                                  j3                                      , j4          , &
         &                                                                                  j5                                      , k

    ! Find interpolating factors.
    if (any(abs(positionCartesian) > self%accelerationX(size(self%accelerationX)))) then
       ! Use spherical approximation.
       radiusSpherical=hypotenuse(positionCartesian)
       if (3*exponent(radiusSpherical) >= maxexponent(0.0d0)) then
          gaussianEllipsoidAccelerationInterpolate=+0.0d0
       else
          gaussianEllipsoidAccelerationInterpolate=-positionCartesian    &
               &                                   /radiusSpherical  **3
       end if
    else
       ! Interpolate in tabulated solution.
       !! Find interpolating factors in position.
       do k=1,3
          if (abs(positionCartesian(k)) < self%accelerationX(1)) then
             ! The acceleration scales as x for sufficiently small x.
             i(  k)= 1
             h(:,k)=[0.0d0,1.0d0]
          else
             positionCartesianLog=log(abs(positionCartesian(k)))
             h(0,k)=+(                              &
                  &   +     positionCartesianLog    &
                  &   -self%accelerationXMinimumLog &
                  &  )                              &
                  & *  self%accelerationXInverseInterval
             i(  k)=+int(h(0,k))+          1
             h(0,k)=+    h(0,k) -dble(i(k)-1    )
             h(1,k)=-    h(0,k) +          1.0d0
          end if
       end do
       !! Find interpolating factors in axis ratios.
       do k=1,2
          if (self%axisRatio(k) < self%accelerationScaleLength(1)) then
             i(  k+3)= 1
             h(:,k+3)=[0.0d0,1.0d0]
          else
             axisRatioLog=log(self%axisRatio(k))
             h(0,k+3)=+(                                             &
                  &     +     axisRatioLog                           &
                  &     -self%accelerationScaleLengthMinimumLog      &
                  &    )                                             &
                  &   *  self%accelerationScaleLengthInverseInterval
             i(  k+3)=+int(h(0,k+3))+            1
             h(0,k+3)=+    h(0,k+3) -dble(i(k+3)-1    )
             h(1,k+3)=-    h(0,k+3) +            1.0d0
          end if
       end do
       !! Do the interpolation.
       gaussianEllipsoidAccelerationInterpolate=0.0d0
       do j1=0,1
          do j2=0,1
             do j3=0,1
                do j4=0,1
                   do j5=0,1
                      gaussianEllipsoidAccelerationInterpolate=+gaussianEllipsoidAccelerationInterpolate &
                           &                                   +self%accelerationVector(                 &
                           &                                                            :      ,         &
                           &                                                            i(1)+j1,         &
                           &                                                            i(2)+j2,         &
                           &                                                            i(3)+j3,         &
                           &                                                            i(4)+j4,         &
                           &                                                            i(5)+j5          &
                           &                                                           )                 &
                           &                                   *(1.0d0-h(j1,1))                          &
                           &                                   *(1.0d0-h(j2,2))                          &
                           &                                   *(1.0d0-h(j3,3))                          &
                           &                                   *(1.0d0-h(j4,4))                          &
                           &                                   *(1.0d0-h(j5,5))
                   end do
                end do
             end do
          end do
       end do
       ! Extrapolate to small radii.
       do k=1,3
          if (abs(positionCartesian(k)) < self%accelerationX(1)) gaussianEllipsoidAccelerationInterpolate(k)=gaussianEllipsoidAccelerationInterpolate(k)*abs(positionCartesian(k))/self%accelerationX(1)
       end do
       ! Reverse direction of acceleration for negative coordinates.
       gaussianEllipsoidAccelerationInterpolate=gaussianEllipsoidAccelerationInterpolate*sign(1.0d0,positionCartesian)
    end if
    return
  end function gaussianEllipsoidAccelerationInterpolate

  subroutine gaussianEllipsoidAccelerationTabulate(self)
    !!{
    Tabulate the gravitational acceleration due to a Gaussian ellipsoid. Follows the approach of \cite[][chapter
    3]{chandrasekhar_ellipsoidal_1987}. Specifically, using his eqn.~(93) for the gravitational potential
    of a heterogeneous ellipsoid:
    \begin{equation}
     \mathfrak{B} = \pi \mathrm{G} a_1 a_2 a_3 \int_0^\infty {\mathrm{d} u \over \Delta(u)} \left[ \psi(1)-\psi(m^2(u)) \right],
    \end{equation}
    where $a_1 \le a_2 \le a_3$ are the semi-axes of the ellipsoid, $\Delta^2(u) = \prod_{i=1}^3 (a_i^2 + u)$,
    $m^2(u)=\sum_{i=1}^3 x_i^2/(a_i^2+u)$, and
    \begin{equation}
     \psi(m^2) = \int_1^{m^2} \mathrm{d}m^2 \rho(m^2),
    \end{equation}
    such that
    \begin{equation}
     {\mathrm{d} \psi(m^2) \over \mathrm{d} m^2} = \rho(m^2).
    \end{equation}
    Then the acceleration is given by
    \begin{equation}
     a_i = - 2 \pi \mathrm{G} a_1 a_2 a_3 x_i \int_0^\infty {\mathrm{d} u \over \Delta(u)} {\rho(m^2) \over (a_i^2+u)}.
    \end{equation}
    We take $a_3=1$ in all cases in order to construct the scale-free solution.
    !!}
    use :: Display                 , only : displayCounter       , displayCounterClear , displayIndent , displayUnindent, &
          &                                 verbosityLevelWorking
    use :: File_Utilities          , only : Directory_Make       , File_Exists         , File_Lock     , File_Path      , &
          &                                 File_Unlock          , lockDescriptor
    use :: Input_Paths             , only : inputPath            , pathTypeDataDynamic
    use :: HDF5_Access             , only : hdf5Access
    use :: IO_HDF5                 , only : hdf5Object
    use :: ISO_Varying_String      , only : char                 , operator(//)        , varying_string
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator
    use :: Numerical_Ranges        , only : Make_Range           , rangeTypeLogarithmic
    implicit none
    class           (massDistributionGaussianEllipsoid), intent(inout)       :: self
    double precision                                   , dimension(3) , save :: positionCartesian          , scaleLength
    double precision                                   , parameter           :: uHighFactor         =1.0d+3
    double precision                                   , parameter           :: xMinimum            =1.0d-3, xMaximum            =5.0d1
    double precision                                   , parameter           :: scaleLengthMinimum  =1.0d-2, scaleLengthMaximum  =1.0d0
    double precision                                   , parameter           :: xPerDecade          =1.0d+1, scaleLengthPerDecade=1.0d1
    type            (integrator                       )               , save :: integrator_
    double precision                                                  , save :: uLow                       , uHigh
    integer                                                           , save :: iAxis                      , l                       , &
         &                                                                      m
    integer                                                                  :: countX                     , countScaleLength        , &
         &                                                                      i                          , j                       , &
         &                                                                      k                          , countWork
    !$omp threadprivate(positionCartesian,scaleLength,integrator_,uLow,uHigh,iAxis,l,m)
    
    ! Return if acceleration is initialized.
    if (self%accelerationInitialized) return
    block
      type(varying_string) :: fileName
      type(hdf5Object    ) :: file
      type(lockDescriptor) :: fileLock
      
      ! Construct a file name for the table.
      fileName=inputPath(pathTypeDataDynamic)// &
           &   'galacticStructure/'          // &
           &   self%objectType()             // &
           &   '.hdf5'
      call Directory_Make(char(File_Path(char(fileName))))
      ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
      call File_Lock(char(fileName),fileLock,lockIsShared=.true.)
      if (File_Exists(fileName)) then
         !$ call hdf5Access%set()
         call file%openFile    (char(fileName      ),readOnly=.true.             )
         call file%readDataset(      'x'            ,self%accelerationX          )
         call file%readDataset(      'scaleLength'  ,self%accelerationScaleLength)
         call file%readDataset(      'acceleration' ,self%accelerationVector     )
         call file%close      (                                                  )
         !$ call hdf5Access%unset()
      else
         ! Generate grid in position and scale length.
         countX          =int(log10(          xMaximum/          xMinimum)*          xPerDecade)+1
         countScaleLength=int(log10(scaleLengthMaximum/scaleLengthMinimum)*scaleLengthPerDecade)+1
         allocate(self%accelerationX          (  countX                                                ))
         allocate(self%accelerationScaleLength(                       countScaleLength                 ))
         allocate(self%accelerationVector     (3,countX,countX,countX,countScaleLength,countScaleLength))
         self%accelerationX          =Make_Range(          xMinimum,          xMaximum,countX          ,rangeTypeLogarithmic)
         self%accelerationScaleLength=Make_Range(scaleLengthMinimum,scaleLengthMaximum,countScaleLength,rangeTypeLogarithmic)
         ! Iterate over all positions and scale lengths, computing the accelerations.
         call displayIndent("tabulating gravitational accelerations for Gaussian ellipsoids",verbosityLevelWorking)
         countWork=0
         do i=1,countX
            do j=1,countX
               !$omp parallel
               integrator_=integrator(accelerationIntegrand,toleranceRelative=1.0d-3)
               !$omp do schedule(dynamic)
               do k=1,countX
                  !$omp atomic
                  countWork=countWork+1
                  call displayCounter(int(100.0d0*dble(countWork)/dble(countX**3)),i == 1 .and. j == 1 .and. k == 1,verbosityLevelWorking)
                  positionCartesian=[self%accelerationX          (i),self%accelerationX          (j),self%accelerationX(k)]
                  do l=1,countScaleLength
                     do m=1,countScaleLength
                        scaleLength=[self%accelerationScaleLength(l),self%accelerationScaleLength(m),1.0d0                ]
                        uLow =0.0d0
                        uHigh=max(maxval(self%scaleLength),maxval(abs(positionCartesian)))*uHighFactor
                        do iAxis=1,3
                           ! Note that the factor of ∏ᵢ₌₁³ aᵢ has been canceled with that appearing in the density normalization.
                           self%accelerationVector(iAxis,i,j,k,l,m)=-2.0d0                         &
                                &                                   *Pi                            &
                                &                                   *positionCartesian(iAxis)      &
                                &                                   *integrator_%integrate(        &
                                &                                                          uLow  , &
                                &                                                          uHigh   &
                                &                                                         )
                        end do
                     end do
                  end do
               end do
               !$omp end do
               !$omp end parallel
            end do
         end do
         call displayCounterClear(       verbosityLevelWorking)
         call displayUnindent     ("done",verbosityLevelWorking)
         !$ call hdf5Access%set()
         call file%openFile    (char   (fileName                    )               ,overWrite=.true.,readOnly=.false.)
         call file%writeDataset(        self%accelerationX           ,'x'                                             )
         call file%writeDataset(        self%accelerationScaleLength ,'scaleLength'                                   )
         call file%writeDataset(        self%accelerationVector      ,'acceleration'                                  )
         call file%close       (                                                                                      )
         !$ call hdf5Access%unset()
      end if
      call File_Unlock(fileLock)
      ! Compute factors needed for interpolation.
      self%accelerationXMinimumLog               =      log(self%accelerationX          (                                1 ))
      self%accelerationXMaximumLog               =      log(self%accelerationX          (size(self%accelerationX          )))
      self%accelerationScaleLengthMinimumLog     =      log(self%accelerationScaleLength(                                1 ))
      self%accelerationScaleLengthMaximumLog     =      log(self%accelerationScaleLength(size(self%accelerationScaleLength)))
      self%accelerationXInverseInterval          =1.0d0/log(self%accelerationX          (2)/self%accelerationX          (1))
      self%accelerationScaleLengthInverseInterval=1.0d0/log(self%accelerationScaleLength(2)/self%accelerationScaleLength(1))
      ! Record that the acceleration table is initialized.
      self%accelerationInitialized=.true.
    end block
    return

  contains

    double precision function accelerationIntegrand(u)
      !!{
      Integrand for the gravitational acceleration.
      !!}
      implicit none
      double precision, intent(in   ) :: u

      accelerationIntegrand=+  densityMSquared(mSquared(u,positionCartesian))    &
           &                /                  Delta   (u                  )     &
           &                /(                                                   &
           &                  +scaleLength    (           iAxis             )**2 &
           &                  +                         u                        &
           &                 )
      return
    end function accelerationIntegrand

    double precision function densityMSquared(mSquared)
      !!{
      Density as a function of $m^2$.
      !!}
      implicit none
      double precision, intent(in   ) :: mSquared
      
      ! Note that the factor of ∏ᵢ₌₁³ aᵢ has been canceled with that appearing in the expression for the gravitational potential.
      densityMSquared=exp(-0.5d0*mSquared)/(2.0d0*Pi)**1.5d0
      return
    end function densityMSquared

    double precision function Delta(u)
      !!{
      The function $\Delta^2(u) = \prod_{i=1}^3 (a_i^2 + u)$.
      !!}
      implicit none
      double precision, intent(in   ) :: u

      Delta=sqrt(product(scaleLength**2+u))
      return
    end function Delta

    double precision function mSquared(u,x)
      !!{
      The function $m^2(u)=\sum_{i=1}^3 x_i^2/(a_i^2+u)$.
      !!}
      implicit none
      double precision, intent(in   )               :: u
      double precision, intent(in   ), dimension(3) :: x

      mSquared=sum(x**2/(scaleLength**2+u))
      return
    end function mSquared

  end subroutine gaussianEllipsoidAccelerationTabulate
  
