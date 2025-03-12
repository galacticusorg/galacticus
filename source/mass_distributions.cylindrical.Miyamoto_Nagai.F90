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
  Implementation of an Miyamoto-Nagai model \citep{miyamoto_three-dimensional_1975} mass distribution class.
  !!}

  use :: Tables, only : table1DLogarithmicLinear

  !![
  <massDistribution name="massDistributionMiyamotoNagai">
   <description>An Miyamoto-Nagai model \citep{miyamoto_three-dimensional_1975} mass distribution class.</description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionCylindrical) :: massDistributionMiyamotoNagai
     !!{
     The Miyamoto-Nagai model \citep{miyamoto_three-dimensional_1975} mass distribution.
     !!}
     private
     double precision                           :: a                        , b                      , &
          &                                        densityNormalization     , mass                   , &
          &                                        shape                    , radiusHalfMassValue
     logical                                    :: surfaceDensityInitialized, massEnclosedInitialized
     type            (table1DLogarithmicLinear) :: surfaceDensityTable      , massEnclosedTable
   contains
     !![
     <methods>
       <method description="Initialize the surface density tabulation." method="surfaceDensityTabulate"/>
       <method description="Initialize the enclosed mass tabulation."   method="massEnclosedTabulate"  />
     </methods>
     !!]
     procedure :: density                    => miyamotoNagaiDensity
     procedure :: densitySphericalAverage    => miyamotoNagaiDensitySphericalAverage
     procedure :: surfaceDensity             => miyamotoNagaiSurfaceDensity
     procedure :: surfaceDensityTabulate     => miyamotoNagaiSurfaceDensityTabulate
     procedure :: surfaceDensityRadialMoment => miyamotoNagaiSurfaceDensityRadialMoment
     procedure :: massEnclosedBySphere       => miyamotoNagaiMassEnclosedBySphere
     procedure :: massEnclosedTabulate       => miyamotoNagaiMassEnclosedTabulate
     procedure :: potentialIsAnalytic        => miyamotoNagaiPotentialIsAnalytic
     procedure :: potential                  => miyamotoNagaiPotential
     procedure :: rotationCurve              => miyamotoNagaiRotationCurve
     procedure :: rotationCurveGradient      => miyamotoNagaiRotationCurveGradient
     procedure :: radiusHalfMass             => miyamotoNagaiRadiusHalfMass
  end type massDistributionMiyamotoNagai

  interface massDistributionMiyamotoNagai
     !!{
     Constructors for the {\normalfont \ttfamily miyamotoNagai} mass distribution class.
     !!}
     module procedure miyamotoNagaiConstructorParameters
     module procedure miyamotoNagaiConstructorInternal
  end interface massDistributionMiyamotoNagai

contains

  function miyamotoNagaiConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily miyamotoNagai} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionMiyamotoNagai)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    double precision                                               :: mass         , a, &
         &                                                            b
    logical                                                        :: dimensionless
    type            (varying_string               )                :: componentType
    type            (varying_string               )                :: massType

    !![
    <inputParameter>
      <name>a</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The $a$ parameter of the MiyamotoNagai profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>b</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The $b$ parameter of the MiyamotoNagai profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>mass</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The mass of the MiyamotoNagai profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>If true the MiyamotoNagai profile is considered to be dimensionless.</description>
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
    <conditionalCall>
     <call>self=massDistributionMiyamotoNagai(componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})</call>
     <argument name="mass"          value="mass"          parameterPresent="parameters"/>
     <argument name="a"             value="a"             parameterPresent="parameters"/>
     <argument name="b"             value="b"             parameterPresent="parameters"/>
     <argument name="dimensionless" value="dimensionless" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function miyamotoNagaiConstructorParameters

  function miyamotoNagaiConstructorInternal(a,b,mass,dimensionless,componentType,massType) result(self)
    !!{
    Internal constructor for ``miyamotoNagai'' mass distribution class.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Comparison    , only : Values_Differ
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionMiyamotoNagai)                          :: self
    double precision                               , intent(in   ), optional :: a            , b, &
         &                                                                      mass
    logical                                        , intent(in   ), optional :: dimensionless
    type            (enumerationComponentTypeType ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType      ), intent(in   ), optional :: massType

    ! Determine if profile is dimensionless.
    if (present(dimensionless)) then
       self%dimensionless=dimensionless
    else
       self%dimensionless=.false.
    end if
    ! If dimensionless, then set scale length and mass to unity.
    if (self%dimensionless) then
       if (.not.present(b))                            call Error_Report('"b" must be specified for dimensionless profiles'                                  //{introspection:location})
       if (present(a   )) then
          if (Values_Differ(a   ,1.0d0,absTol=1.0d-6)) call Error_Report('"a" should be unity for a dimensionless profile (or simply do not specify it)'     //{introspection:location})
       end if
       if (present(mass)) then
          if (Values_Differ(mass,1.0d0,absTol=1.0d-6)) call Error_Report('mass should be unity for a dimensionless profile (or simply do not specify a mass)'//{introspection:location})
       end if
       self%a                   =1.0d0
       self%b                   =b
       self%mass                =1.0d0
    else
       if (.not.present(a   )) call Error_Report('"a" must be specified for dimensionful profiles' //{introspection:location})
       if (.not.present(b   )) call Error_Report('"b" must be specified for dimensionful profiles' //{introspection:location})
       if (.not.present(mass)) call Error_Report('mass must be specified for dimensionful profiles'//{introspection:location})
       self%a                   =a
       self%b                   =b
       self%mass                =mass
    end if
    ! Compute shape parameter and density normalization.
    self%shape               =+self%b       &
         &                    /self%a
    self%densityNormalization=              &
         &                    +self%mass    &
         &                    /self%a   **2 &
         &                    /self%b       &
         &                    /4.0d0        &
         &                    /Pi
    ! Set table initialization states.
    self%surfaceDensityInitialized=.false.
    self%  massEnclosedInitialized=.false.
    return
  end function miyamotoNagaiConstructorInternal

  double precision function miyamotoNagaiDensity(self,coordinates)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in an \citep{miyamoto_three-dimensional_1975} disk mass distribution.
    !!}
    use :: Coordinates, only : assignment(=), coordinateCylindrical
    implicit none
    class           (massDistributionMiyamotoNagai), intent(inout) :: self
    class           (coordinate                   ), intent(in   ) :: coordinates
    type            (coordinateCylindrical        )                :: position
    double precision                                               :: r          , z

    ! Get position in cylindrical coordinate system.
    position=coordinates
    ! Compute density.
    r=    position%r() /self%a
    z=abs(position%z())/self%b
    ! Evaluate the density.
    miyamotoNagaiDensity=+self%densityNormalization &
         &               *(                         &
         &                 +r**2                    &
         &                 +(                       &
         &                   +1.0d0                 &
         &                   +self%shape            &
         &                   *sqrt(                 &
         &                         +1.0d0           &
         &                         +z**2            &
         &                        )                 &
         &                  )**2                    &
         &                 *(                       &
         &                   +1.0d0                 &
         &                   +3.0d0                 &
         &                   *self%shape            &
         &                   *sqrt(                 &
         &                         +1.0d0           &
         &                         +z**2            &
         &                        )                 &
         &                  )                       &
         &                )                         &
         &               /(                         &
         &                 +r**2                    &
         &                 +(                       &
         &                   +1.0d0                 &
         &                   +self%shape            &
         &                   *sqrt(                 &
         &                         +1.0d0           &
         &                         +z**2            &
         &                        )                 &
         &                  )**2                    &
         &                )**2.5d0                  &
         &               /(                         &
         &                 +1.0d0                   &
         &                 +z**2                    &
         &                )**1.5d0
    return
  end function miyamotoNagaiDensity

  double precision function miyamotoNagaiDensitySphericalAverage(self,radius)
    !!{
    Return the spherically-averaged density at the specified {\normalfont \ttfamily radius} in an \citep{miyamoto_three-dimensional_1975} disk mass distribution.
    !!}
    implicit none
    class           (massDistributionMiyamotoNagai), intent(inout) :: self
    double precision                               , intent(in   ) :: radius
 
    miyamotoNagaiDensitySphericalAverage=0.0d0
    call Error_Report('spherically-averaged density profile is not implemented'//{introspection:location})
    return
  end function miyamotoNagaiDensitySphericalAverage

  double precision function miyamotoNagaiMassEnclosedBySphere(self,radius)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for \citep{miyamoto_three-dimensional_1975} disk mass distributions.
    !!}
    implicit none
    class           (massDistributionMiyamotoNagai), intent(inout), target :: self
    double precision                               , intent(in   )         :: radius
    
    ! Ensure mass enclosed profile is tabulated.
    call self%massEnclosedTabulate()
    ! Evaluate the mass enclosed.
    miyamotoNagaiMassEnclosedBySphere=self%massEnclosedTable%interpolate(radius)
    return
  end function miyamotoNagaiMassEnclosedBySphere

  double precision function miyamotoNagaiRadiusHalfMass(self)
    !!{
    Return the half-mass radius in a Miyamoto-Nagai mass distribution.
    !!}
    implicit none
    class(massDistributionMiyamotoNagai), intent(inout) :: self

    ! Ensure mass enclosed profile is tabulated.
    call self%massEnclosedTabulate()
    ! Return the half-mass radius.
    miyamotoNagaiRadiusHalfMass=self%radiusHalfMassValue
    return
  end function miyamotoNagaiRadiusHalfMass

  subroutine miyamotoNagaiMassEnclosedTabulate(self)
    !!{
    Construct a tabulation of the mass enclosed by a sphere in a Miyamoto-Nagai mass distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator
    use :: Table_Labels            , only : extrapolationTypeExtrapolate
    implicit none
    class           (massDistributionMiyamotoNagai), intent(inout) :: self
    double precision                               , parameter     :: radiusMinimum       =1.0d-6, radiusMaximum =1.0d3
    double precision                               , parameter     :: heightStart         =0.0d+0, heightInfinity=1.0d3
    integer                                        , parameter     :: radiusCount         =300
    type            (integrator                   )                :: integrator_
    integer                                                        :: i
    double precision                                               :: radius

    ! Mass enclosed is not analytically tractable. Tabulate the results of numerical integration.
    if (.not.self%massEnclosedInitialized) then
       ! Initialize a table of radii.
       call self%massEnclosedTable%create(radiusMinimum*self%a,radiusMaximum*self%a,radiusCount,1,extrapolationType=spread(extrapolationTypeExtrapolate,1,2))
       ! Compute enclosed mass at each radius.
       integrator_=integrator(integrandR,toleranceRelative=1.0d-3)
       do i=1,radiusCount
          call self                                                             &
               & %massEnclosedTable                                             &
               &  %populate(                                                    &
               &            +4.0d0                                              &
               &            *Pi                                                 &
               &            *integrator_%integrate(                             &
               &                                   0.0d0                      , &
               &                                   self%massEnclosedTable%x(i)  &
               &                                  )                           , &
               &                                                            i   &
               &           )
       end do
       ! Find the half-mass radius.
       do i=2,radiusCount
          if     (                                                                     &
               &   self%massEnclosedTable%y(i  ) >= 0.5d0*self%massEnclosedTable%y(-1) &
               &  .and.                                                                &
               &   self%massEnclosedTable%y(i-1) <  0.5d0*self%massEnclosedTable%y(-1) &
               & ) self%radiusHalfMassValue=+  self%massEnclosedTable%x(i-1)           &
               &                            +(                                         &
               &                              +self%massEnclosedTable%x(i  )           &
               &                              -self%massEnclosedTable%x(i-1)           &
               &                             )                                         &
               &                            *(                                         &
               &                              +0.5d0                                   &
               &                              *self%massEnclosedTable%y( -1)           &
               &                              -self%massEnclosedTable%y(i-1)           &
               &                             )                                         &
               &                            /(                                         &
               &                              +self%massEnclosedTable%y(i  )           &
               &                              -self%massEnclosedTable%y(i-1)           &
               &                            )
       end do
       ! Record that surface density is now initialized.
       self%massEnclosedInitialized=.true.
    end if
    return

  contains

    double precision function integrandR(r)
      !!{
      Integrand function used for finding the mass enclosed by a sphere in Miyamoto-Nagai disks.
      !!}
      implicit none
      double precision            , intent(in   ) :: r
      type            (integrator)                :: integrator1_

      radius    =+r
      integrator1_=integrator(integrandZ,toleranceRelative=1.0d-3)
      integrandR=+                             radius                          &
           &     *integrator1_%integrate(                                      &
           &                             0.0d0                               , &
           &                             sqrt(                                 &
           &                                  +self%massEnclosedTable%x(i)**2  &
           &                                  -radius                     **2  &
           &                                 )                                 &
           &                            )
      return
    end function integrandR

    double precision function integrandZ(z)
      !!{
      Integrand function used for finding the mass enclosed by a sphere in Miyamoto-Nagai disks.
      !!}
      use :: Coordinates, only : assignment(=), coordinateCylindrical
      implicit none
      double precision                       , intent(in   ) :: z
      type            (coordinateCylindrical)                :: position

      position  =[radius,0.0d0,z]
      integrandZ=self%density(position)
      return
    end function integrandZ

  end subroutine miyamotoNagaiMassEnclosedTabulate

  subroutine miyamotoNagaiSurfaceDensityTabulate(self)
    !!{
    Construct a tabulation of the surface density profile in a Miyamoto-Nagai mass distribution.
    !!}
    use :: Numerical_Integration, only : integrator
    use :: Table_Labels         , only : extrapolationTypeExtrapolate
    implicit none
    class           (massDistributionMiyamotoNagai), intent(inout) :: self
    double precision                               , parameter     :: radiusMinimum=1.0d-6, radiusMaximum =1.0d3
    double precision                               , parameter     :: heightStart  =0.0d+0, heightInfinity=1.0d3
    integer                                        , parameter     :: radiusCount  =300
    type            (integrator                   )                :: integrator_
    integer                                                        :: i

    ! Surface density is not analytically tractable. Tabulate the results of numerical integration.
    if (.not.self%surfaceDensityInitialized) then
       ! Initialize a table of radii.
       call self%surfaceDensityTable%create(radiusMinimum*self%a,radiusMaximum*self%a,radiusCount,1,extrapolationType=spread(extrapolationTypeExtrapolate,1,2))
       ! Compute surface density at each radius.
       integrator_=integrator(integrandSurfaceDensity,toleranceRelative=1.0d-3)
       do i=1,radiusCount
          call self                                                       &
               & %surfaceDensityTable                                     &
               &  %populate(                                              &
               &            +2.0d0                                        &
               &            *integrator_%integrate(                       &
               &                                   self%b*heightStart   , &
               &                                   self%b*heightInfinity  &
               &                                  )                     , &
               &            i                                             &
               &           )
       end do
       ! Record that surface density is now initialized.
       self%surfaceDensityInitialized=.true.
    end if
    return

  contains

    double precision function integrandSurfaceDensity(height)
      !!{
      Integrand function used for finding the surface density of Miyamoto-Nagai disks.
      !!}
      use :: Coordinates, only : assignment(=), coordinateCylindrical
      implicit none
      double precision                       , intent(in   ) :: height
      type            (coordinateCylindrical)                :: position

      position                =[self%surfaceDensityTable%x(i),0.0d0,height]
      integrandSurfaceDensity=self%density(position)
      return
    end function integrandSurfaceDensity

  end subroutine miyamotoNagaiSurfaceDensityTabulate

  double precision function miyamotoNagaiSurfaceDensity(self,coordinates)
    !!{
    Return the surface density at the specified {\normalfont \ttfamily coordinates} in a Miyamoto-Nagai mass distribution.
    !!}
    use :: Coordinates, only : assignment(=), coordinateCylindrical
    implicit none
    class(massDistributionMiyamotoNagai), intent(inout) :: self
    class(coordinate                   ), intent(in   ) :: coordinates
    type (coordinateCylindrical        )                :: position

    ! Ensure surface density profile is tabulated.
    call self%surfaceDensityTabulate()
    ! Evaluate the surface density.
    position=coordinates
    miyamotoNagaiSurfaceDensity=self%surfaceDensityTable%interpolate(position%r())
    return
  end function miyamotoNagaiSurfaceDensity

  double precision function miyamotoNagaiSurfaceDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite)
    !!{
    Compute radial moments of the Miyamoto-Nagai mass distribution surface density profile.
    !!}
    use :: Error , only : Error_Report
    use :: Tables, only : tablesIntegrationWeightFunction
    implicit none
    class           (massDistributionMiyamotoNagai  ), intent(inout)           :: self
    double precision                                 , intent(in   )           :: moment
    double precision                                 , intent(in   ), optional :: radiusMinimum          , radiusMaximum
    logical                                          , intent(  out), optional :: isInfinite
    procedure       (tablesIntegrationWeightFunction), pointer                 :: integrandWeightFunction
    double precision                                                           :: radiusMinimumActual    , radiusMaximumActual

    ! Set infinity status.
    if (present(isInfinite)) isInfinite=.false.
    ! Ensure surface density profile is tabulated.
    call self%surfaceDensityTabulate()
    ! Determine the radii to use.
    if (present(radiusMinimum)) then
       radiusMinimumActual=radiusMinimum
    else
       radiusMinimumActual=self%surfaceDensityTable%x(+1)
    end if
    if (present(radiusMaximum)) then
       radiusMaximumActual=radiusMaximum
    else
       radiusMaximumActual=self%surfaceDensityTable%x(-1)
    end if
    ! Check for finite moments.
    if     (                                                                     &
         &   (moment <= -1.0d0 .and.              radiusMinimumActual  <= 0.0d0) &
         &  .or.                                                                 &
         &   (moment >=  2.0d0 .and. .not.present(radiusMaximum      )         ) &
         & ) then
       miyamotoNagaiSurfaceDensityRadialMoment=0.0d0
       if (present(isInfinite)) then
          isInfinite=.true.
          return
       else
          call Error_Report('-1 < momemnt < 2 is required for finite moment'//{introspection:location})
       end if
    end if
    ! Perform the integration.
    integrandWeightFunction => momentWeight
    miyamotoNagaiSurfaceDensityRadialMoment=                                                   &
         & +sum(                                                                               &
         &      +         self%surfaceDensityTable%integrationWeights(                         &
         &                                                            radiusMinimumActual    , &
         &                                                            radiusMaximumActual    , &
         &                                                            integrandWeightFunction  &
         &                                                           )                         &
         &      *reshape(                                                                      &
         &                self%surfaceDensityTable%ys                (                         &
         &                                                           )                       , &
         &               [                                                                     &
         &                self%surfaceDensityTable%size              (                         &
         &                                                           )                         &
         &               ]                                                                     &
         &              )                                                                      &
         &     )
    return

  contains

    double precision function momentWeight(radius)
      !!{
      The weight function used in computing radial moments of the Miyamoto-Nagai mass distribution surface density.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      momentWeight=radius**moment
      return
    end function momentWeight

  end function miyamotoNagaiSurfaceDensityRadialMoment

  double precision function miyamotoNagaiRotationCurve(self,radius)
    !!{
    Return the mid-plane rotation curve for a Miyamoto-Nagai mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionMiyamotoNagai), intent(inout) :: self
    double precision                               , intent(in   ) :: radius
    double precision                                               :: r

    ! Get dimensionless radius.
    r=radius/self%a
    ! Evaluate the rotation curve.
    miyamotoNagaiRotationCurve=+sqrt(           &
         &                           +self%mass &
         &                           /self%a    &
         &                          )           &
         &                     *r               &
         &                     /(               &
         &                       +r**2          &
         &                       +(             &
         &                         +1.0d0       &
         &                         +self%shape  &
         &                        )**2          &
         &                      )**0.75d0
    ! Make dimensionfull if necessary.
    if (.not.self%dimensionless)                                      &
         & miyamotoNagaiRotationCurve=+gravitationalConstant_internal &
         &                            *miyamotoNagaiRotationCurve
    return
  end function miyamotoNagaiRotationCurve

  double precision function miyamotoNagaiRotationCurveGradient(self,radius)
    !!{
    Return the mid-plane rotation curve gradient for an \citep{miyamoto_three-dimensional_1975} disk.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionMiyamotoNagai), intent(inout) :: self
    double precision                               , intent(in   ) :: radius
    double precision                                               :: r

    ! Get dimensionless radius.
    r=radius/self%a
    ! Evaluate the rotation curve.
    miyamotoNagaiRotationCurveGradient=+self%mass        &
         &                             /self%a   **2     &
         &                             *(                &
         &                               +2.0d0          &
         &                               *r              &
         &                               -3.0d0          &
         &                               *r      **3     &
         &                               /(              &
         &                                 +r    **2     &
         &                                 +(            &
         &                                   +1.0d0      &
         &                                   +self%shape &
         &                                  )**2         &
         &                                )              &
         &                              )                &
         &                             /(                &
         &                               +r**2           &
         &                               +(              &
         &                                 +1.0d0        &
         &                                 +self%shape   &
         &                                )**2           &
         &                              )**1.5d0
    ! Make dimensionfull if necessary.
    if (.not.self%dimensionless)                                                 &
         & miyamotoNagaiRotationCurveGradient=+gravitationalConstant_internal    &
         &                                    *miyamotoNagaiRotationCurveGradient
    return
  end function miyamotoNagaiRotationCurveGradient

  logical function miyamotoNagaiPotentialIsAnalytic(self) result(isAnalytic)
    !!{
    Return that the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionMiyamotoNagai), intent(inout) :: self

    isAnalytic=.true.
    return
  end function miyamotoNagaiPotentialIsAnalytic

  double precision function miyamotoNagaiPotential(self,coordinates,status)
    !!{
    Return the gravitational potential for an \citep{miyamoto_three-dimensional_1975} disk.
    !!}
    use :: Coordinates                     , only : assignment(=)                 , coordinateCylindrical
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionMiyamotoNagai    ), intent(inout), target   :: self
    class           (coordinate                       ), intent(in   )           :: coordinates
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status
    type            (coordinateCylindrical            )                          :: position
    double precision                                                             :: r          , z

    if (present(status)) status=structureErrorCodeSuccess
    ! Get position in cylindrical coordinate system.
    position=coordinates
    ! Compute density.
    r=    position%r() /self%a
    z=abs(position%z())/self%b
    ! Evaluate the potential.
    miyamotoNagaiPotential=-self%mass                       &
         &                 /self%a                          &
         &                 /sqrt(                           &
         &                       +r**2                      &
         &                       +(                         &
         &                         +1.0d0                   &
         &                         +self%shape              &
         &                         *sqrt(                   &
         &                               +1.0d0             &
         &                               +z**2              &
         &                              )                   &
         &                        )**2                      &
         &                      )
    ! Make dimensionfull if necessary.
    if (.not.self%dimensionless)                                  &
         & miyamotoNagaiPotential=+gravitationalConstant_internal &
         &                        *miyamotoNagaiPotential
    return
  end function miyamotoNagaiPotential
