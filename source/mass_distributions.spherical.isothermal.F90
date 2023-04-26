!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  Implementation of an isothermal mass distribution class.
  !!}

  !![
  <massDistribution name="massDistributionIsothermal">
    <description>
      An isothermal mass distribution class in which the density profile is given by:
      \begin{equation}
      \rho_\mathrm{dark matter}(r) \propto r^{-2}.
      \end{equation}
   </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionIsothermal
     !!{
     The isothermal mass distribution.
     !!}
     private
     double precision :: densityNormalization, lengthReference
   contains
     procedure :: massTotal             => isothermalMassTotal
     procedure :: density               => isothermalDensity
     procedure :: densityGradientRadial => isothermalDensityGradientRadial
     procedure :: densityRadialMoment   => isothermalDensityRadialMoment
     procedure :: massEnclosedBySphere  => isothermalMassEnclosedBySphere
     procedure :: potential             => isothermalPotential
     procedure :: descriptor            => isothermalDescriptor
  end type massDistributionIsothermal

  interface massDistributionIsothermal
     !!{
     Constructors for the {\normalfont \ttfamily isothermal} mass distribution class.
     !!}
     module procedure massDistributionIsothermalConstructorParameters
     module procedure massDistributionIsothermalConstructorInternal
  end interface massDistributionIsothermal

contains

  function massDistributionIsothermalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily isothermal} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    use :: Numerical_Constants_Math  , only : Pi
    implicit none
    type            (massDistributionIsothermal)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    double precision                                            :: mass                , lengthReference, &
         &                                                         densityNormalization
    logical                                                     :: dimensionless
    type            (varying_string            )                :: componentType
    type            (varying_string            )                :: massType

    !![
    <inputParameter>
      <name>densityNormalization</name>
      <defaultValue>0.25d0/Pi</defaultValue>
      <description>The density normalization of the isothermal profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>lengthReference</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The scale radius of the isothermal profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>mass</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The mass of the isothermal profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>If true the isothermal profile is considered to be dimensionless.</description>
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
     <call>self=massDistributionIsothermal(componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})</call>
     <argument name="densityNormalization" value="densityNormalization" parameterPresent="parameters"/>
     <argument name="mass"                 value="mass"                 parameterPresent="parameters"/>
     <argument name="lengthReference"      value="lengthReference"      parameterPresent="parameters"/>
     <argument name="dimensionless"        value="dimensionless"        parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massDistributionIsothermalConstructorParameters

  function massDistributionIsothermalConstructorInternal(densityNormalization,mass,lengthReference,dimensionless,componentType,massType) result(self)
    !!{
    Internal constructor for ``isothermal'' mass distribution class.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Comparison    , only : Values_Differ
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionIsothermal   )                          :: self
    double precision                               , intent(in   ), optional :: densityNormalization, mass, &
         &                                                                      lengthReference
    logical                                        , intent(in   ), optional :: dimensionless
    type            (enumerationComponentTypeType ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType      ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="componentType, massType"/>
    !!]

    ! Determine if profile is dimensionless.
    self%dimensionless=.false.
    if (present(dimensionless)) self%dimensionless=dimensionless
    ! If dimensionless, then set scale length and mass to unity.
    if (self%dimensionless) then
       if (present(lengthReference     )) then
          if (Values_Differ(lengthReference     ,1.0d0    ,absTol=1.0d-6)) call Error_Report('scaleLength should be unity for a dimensionless profile (or simply do not specify a scale length)'               //{introspection:location})
       end if
       if (present(mass                )) then
          if (Values_Differ(mass                ,1.0d0    ,absTol=1.0d-6)) call Error_Report('mass should be unity for a dimensionless profile (or simply do not specify a mass)'                              //{introspection:location})
       end if
       if (present(densityNormalization)) then
          if (Values_Differ(densityNormalization,0.25d0/Pi,absTol=1.0d-6)) call Error_Report('densityNormalization should be π/4 for a dimensionless profile (or simply do not specify a densityNormalization)'//{introspection:location})
       end if
       self%lengthReference     =1.00d0
       self%densityNormalization=0.25d0/Pi
    else
       if      (present(lengthReference     )) then
          self%lengthReference=lengthReference
       else
          call Error_Report('"lengthReference" must be specified'//{introspection:location})
       end if
       if      (present(densityNormalization)) then
          self%densityNormalization=densityNormalization
       else if (present(mass                )) then
          self%densityNormalization=mass                /4.0d0/Pi/lengthReference**3
       else
          call Error_Report('one of "densityNormalization" or "mass" must be specified'//{introspection:location})
       end if
    end if
    return
  end function massDistributionIsothermalConstructorInternal

  double precision function isothermalMassTotal(self,componentType,massType)
    !!{
    Return the total mass in an Isothermal mass distribution.
    !!}
    implicit none
    class(massDistributionIsothermal  ), intent(inout)           :: self
    type (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type (enumerationMassTypeType     ), intent(in   ), optional :: massType

    if (self%matches(componentType,massType)) then
       isothermalMassTotal=huge(0.0d0)
    else
       isothermalMassTotal=0.0d0
    end if
    return
  end function isothermalMassTotal

  double precision function isothermalDensity(self,coordinates,componentType,massType)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in an isothermal mass distribution.
    !!}
    use :: Coordinates, only : assignment(=), coordinateSpherical
    implicit none
    class(massDistributionIsothermal  ), intent(inout)           :: self
    class(coordinate                  ), intent(in   )           :: coordinates
    type (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type (enumerationMassTypeType     ), intent(in   ), optional :: massType

    if (.not.self%matches(componentType,massType)) then
       isothermalDensity=0.0d0
       return
    end if
    isothermalDensity=+ self       %densityNormalization   &
         &            /(                                   &
         &             +coordinates%rSpherical          () &
         &             /self       %lengthReference        &
         &            )**2
    return
  end function isothermalDensity

  double precision function isothermalDensityGradientRadial(self,coordinates,logarithmic,componentType,massType)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in an isothermal mass distribution.
    !!}
    implicit none
    class           (massDistributionIsothermal  ), intent(inout)           :: self
    class           (coordinate                  ), intent(in   )           :: coordinates
    logical                                       , intent(in   ), optional :: logarithmic
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: radius
    logical                                                                 :: logarithmicActual

    if (.not.self%matches(componentType,massType)) then
       isothermalDensityGradientRadial=0.0d0
       return
    end if
    ! Set default options.
    logarithmicActual=.false.
    if (present(logarithmic)) logarithmicActual=logarithmic
    ! Get position in spherical coordinate system.
    radius=coordinates%rSpherical()
    ! Compute density gradient.
    if (logarithmicActual) then
       isothermalDensityGradientRadial=-2.0d0
    else
       isothermalDensityGradientRadial=-2.0d0                        &
            &                          *self%densityNormalization    &
            &                          *self%lengthReference     **2 &
            &                          /radius                   **3
    end if
    return
  end function isothermalDensityGradientRadial

  double precision function isothermalMassEnclosedBySphere(self,radius,componentType,massType)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for isothermal mass distributions.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionIsothermal  ), intent(inout), target   :: self
    double precision                              , intent(in   )           :: radius
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType

    if (.not.self%matches(componentType,massType)) then
       isothermalMassEnclosedBySphere=0.0d0
       return
    end if
    isothermalMassEnclosedBySphere=+4.0d0                        &
         &                         *Pi                           &
         &                         *self%densityNormalization    &
         &                         *self%lengthReference     **2 &
         &                         *radius
    return
  end function isothermalMassEnclosedBySphere

  double precision function isothermalDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite,componentType,massType)
    !!{
    Returns a radial density moment for the Isothermal mass distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (massDistributionIsothermal  ), intent(inout)           :: self
    double precision                              , intent(in   )           :: moment
    double precision                              , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                                       , intent(  out), optional :: isInfinite
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: momentMinimum, momentMaximum

    if (.not.self%matches(componentType,massType)) then
       isothermalDensityRadialMoment=0.0d0
       return
    end if
    isothermalDensityRadialMoment=0.0d0
    if (present(isInfinite)) isInfinite=.false.
    if (present(radiusMinimum)) then
       if (moment == 1.0d0) then
          momentMinimum=+log(+radiusMinimum                )
       else
          momentMinimum=     +radiusMinimum**(moment-1.0d0)
       end if
    else if (moment <= 1.0d0) then
       momentMinimum=+0.0d0
       if (present(isInfinite)) then
          isInfinite=.true.
          return
       else
          call Error_Report('radial moment is infinite'//{introspection:location})
       end if
    else
       momentMinimum=0.0d0
    end if
    if (present(radiusMaximum)) then
       if (moment == 1.0d0) then
          momentMaximum=+log(+radiusMaximum                )
       else
          momentMaximum=     +radiusMaximum**(moment-1.0d0)
       end if
    else if (moment >= 1.0d0) then
       momentMaximum=+0.0d0
       if (present(isInfinite)) then
          isInfinite=.true.
          return
       else
          call Error_Report('radial moment is infinite'//{introspection:location})
       end if
    else
       momentMaximum=0.0d0
    end if
    isothermalDensityRadialMoment=+momentMaximum &
         &                        -momentMinimum
    return
  end function isothermalDensityRadialMoment

  double precision function isothermalPotential(self,coordinates,componentType,massType,status)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in an isothermal mass distribution.
    !!}
    use :: Coordinates                     , only : assignment(=)
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess      , structureErrorCodeInfinite
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Error                           , only : Error_Report
    implicit none
    class(massDistributionIsothermal       ), intent(inout)           :: self
    class(coordinate                       ), intent(in   )           :: coordinates
    type (enumerationComponentTypeType     ), intent(in   ), optional :: componentType
    type (enumerationMassTypeType          ), intent(in   ), optional :: massType
    type (enumerationStructureErrorCodeType), intent(  out), optional :: status

    if (present(status)) status=structureErrorCodeSuccess
    if (.not.self%matches(componentType,massType)) then
       isothermalPotential=0.0d0
       return
    end if
    ! Compute the potential at this position.
    if (coordinates%rSpherical() <= 0.0d0) then
       isothermalPotential=0.0d0
       if (present(status)) then
          status=structureErrorCodeInfinite
          return
       else
          call Error_Report('potential is divergent at zero radius'//{introspection:location})
       end if
    end if
    isothermalPotential=-     self       %densityNormalization      &
         &              *     self       %lengthReference       **2 &
         &              *log(                                       &
         &                   +coordinates%rSpherical          ()    &
         &                   /self       %lengthReference           &
         &                  )
    if (.not.self%isDimensionless()) isothermalPotential=+gravitationalConstantGalacticus &
         &                                               *isothermalPotential
    return
  end function isothermalPotential

  subroutine isothermalDescriptor(self,descriptor,includeClass)
    !!{
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class    (massDistributionIsothermal), intent(inout)           :: self
    type     (inputParameters           ), intent(inout)           :: descriptor
    logical                              , intent(in   ), optional :: includeClass
    character(len=18                    )                          :: parameterLabel
    type     (inputParameters           )                          :: parameters

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('massDistribution','isothermal')
    parameters=descriptor%subparameters('massDistribution')
    write (parameterLabel,'(e17.10)') self%densityNormalization
    call parameters%addParameter('densityNormalization',trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%lengthReference
    call parameters%addParameter('lengthReference'     ,trim(adjustl(parameterLabel)))
    return
  end subroutine isothermalDescriptor
