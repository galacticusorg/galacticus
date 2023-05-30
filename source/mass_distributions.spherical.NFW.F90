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
  Implementation of an NFW \citep{navarro_structure_1996} mass distribution class.
  !!}

  !![
  <massDistribution name="massDistributionNFW">
   <description>An NFW \citep{navarro_structure_1996} mass distribution class.</description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionNFW
     !!{
     The NFW \citep{navarro_structure_1996} mass distribution.
     !!}
     private
     double precision :: densityNormalization, scaleLength
   contains
     procedure :: density               => nfwDensity
     procedure :: densityGradientRadial => nfwDensityGradientRadial
     procedure :: densityRadialMoment   => nfwDensityRadialMoment
     procedure :: descriptor            => nfwDescriptor
  end type massDistributionNFW

  interface massDistributionNFW
     !!{
     Constructors for the {\normalfont \ttfamily nfw} mass distribution class.
     !!}
     module procedure nfwConstructorParameters
     module procedure nfwConstructorInternal
  end interface massDistributionNFW

contains

  function nfwConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily nfw} mass distribution class which builds the object from a parameter
    set.
    !!}
     use :: Input_Parameters          , only : inputParameter                , inputParameters
     use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
     use :: Numerical_Constants_Math  , only : Pi
    implicit none
    type            (massDistributionNFW)                :: self
    type            (inputParameters    ), intent(inout) :: parameters
    double precision                                     :: mass                , scaleLength  , &
         &                                                  densityNormalization, concentration, &
         &                                                  virialRadius
    logical                                              :: dimensionless
    type            (varying_string     )                :: componentType
    type            (varying_string     )                :: massType

    !![
    <inputParameter>
      <name>densityNormalization</name>
      <defaultValue>0.5d0/Pi</defaultValue>
      <description>The density normalization of the NFW profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>scaleLength</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The scale radius of the NFW profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>mass</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The mass of the NFW profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>concentration</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The concentration of the NFW profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>virialRadius</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The virial radius of the NFW profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>If true the NFW profile is considered to be dimensionless.</description>
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
     <call>self=massDistributionNFW(componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})</call>
     <argument name="densityNormalization" value="densityNormalization" parameterPresent="parameters"/>
     <argument name="mass"                 value="mass"                 parameterPresent="parameters"/>
     <argument name="scaleLength"          value="scaleLength"          parameterPresent="parameters"/>
     <argument name="virialRadius"         value="virialRadius"         parameterPresent="parameters"/>
     <argument name="concentration"        value="concentration"        parameterPresent="parameters"/>
     <argument name="dimensionless"        value="dimensionless"        parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nfwConstructorParameters

  function nfwConstructorInternal(scaleLength,concentration,densityNormalization,mass,virialRadius,dimensionless,componentType,massType) result(self)
    !!{
    Internal constructor for ``nfw'' mass distribution class.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionNFW         )                          :: self
    double precision                              , intent(in   ), optional :: scaleLength         , concentration, &
         &                                                                     densityNormalization, mass         , &
         &                                                                     virialRadius
    logical                                       , intent(in   ), optional :: dimensionless
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: radiusScaleFree
    !![
    <constructorAssign variables="componentType, massType"/>
    !!]

    ! Determine scale length
    if      (                            &
         &   present(scaleLength  )      &
         &  ) then
       self%scaleLength         =scaleLength
    else if (                            &
         &   present(concentration).and. &
         &   present(virialRadius )      &
         &  ) then
       self%scaleLength=virialRadius/concentration
    else
       call Error_Report('no means to determine scale length'//{introspection:location})
    end if
    ! Determine density normalization.
    if      (                                   &
         &   present(densityNormalization)      &
         &  ) then
       self%densityNormalization=densityNormalization
    else if (                                   &
         &   present(mass                ).and. &
         &   present(virialRadius        )      &
         &  ) then
       radiusScaleFree          =+virialRadius/self%scaleLength
       self%densityNormalization=+mass/4.0d0/Pi/self%scaleLength**3/(log(1.0d0+radiusScaleFree)-radiusScaleFree/(1.0d0+radiusScaleFree))
    else
       call Error_Report('either "densityNormalization", or "mass" and "virialRadius" must be specified'//{introspection:location})
    end if
    ! Determine if profile is dimensionless.
    if      (present(dimensionless     )) then
       self%dimensionless=dimensionless
    else
       self%dimensionless=.false.
    end if
    return
  end function nfwConstructorInternal

  double precision function nfwDensity(self,coordinates,componentType,massType)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in an NFW mass distribution.
    !!}
    implicit none
    class           (massDistributionNFW         ), intent(inout)           :: self
    class           (coordinate                  ), intent(in   )           :: coordinates
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: radiusScaleFree

    if (.not.self%matches(componentType,massType)) then
       nfwDensity=0.0d0
       return
    end if    
    ! Compute the density at this position.
    radiusScaleFree=+coordinates%rSpherical          () &
         &          /self       %scaleLength
    nfwDensity     =+self       %densityNormalization   &
         &          /       radiusScaleFree             &
         &          /(1.0d0+radiusScaleFree)**2
    return
  end function nfwDensity

  double precision function nfwDensityGradientRadial(self,coordinates,logarithmic,componentType,massType) result(densityGradientRadial)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in an NFW \citep{navarro_structure_1996} mass distribution.
    !!}
    implicit none
    class           (massDistributionNFW         ), intent(inout), target   :: self
    class           (coordinate                  ), intent(in   )           :: coordinates
    logical                                       , intent(in   ), optional :: logarithmic
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: radiusScaleFree
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    densityGradientRadial=0.0d0
    if (.not.self%matches(componentType,massType)) return
    radiusScaleFree      =+coordinates%rSpherical()         &
         &                /self       %scaleLength
    densityGradientRadial=-self       %densityNormalization &
         &                /self       %scaleLength          &
         &                /             radiusScaleFree **2 &
         &                *(1.0d0+3.0d0*radiusScaleFree)    &
         &                /(1.0d0+      radiusScaleFree)**3
    if (logarithmic_) densityGradientRadial=+            densityGradientRadial              &
         &                                  /self       %density              (coordinates) &
         &                                  *coordinates%rSpherical           (           )
    return
  end function nfwDensityGradientRadial

  double precision function nfwDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite,componentType,massType) result(densityRadialMoment)
    !!{
    Computes radial moments of the density in an NFW \citep{navarro_structure_1996} mass distribution.
    !!}
    implicit none
    class           (massDistributionNFW         ), intent(inout)           :: self
    double precision                              , intent(in   )           :: moment
    double precision                              , intent(in   ), optional :: radiusMinimum      , radiusMaximum
    logical                                       , intent(  out), optional :: isInfinite
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: radialMomentMinimum, radialMomentMaximum

    densityRadialMoment=0.0d0
    if (present(isInfinite)) isInfinite=.false.
    if (.not.self%matches(componentType,massType)) return
    if (present(radiusMinimum)) then
       radialMomentMinimum=radialMomentScaleFree(radiusMinimum/self%scaleLength)
    else
       radialMomentMinimum=radialMomentScaleFree(                         0.0d0)
    end if
    if (present(radiusMaximum)) then
       radialMomentMaximum=radialMomentScaleFree(radiusMaximum/self%scaleLength)
    else
       radialMomentMaximum=0.0d0
       if (moment >= 3.0d0) then
          if (present(isInfinite)) then
             isInfinite=.true.
             return
          else
             call Error_Report('moment is infinite'//{introspection:location})
          end if
       end if
    end if
    densityRadialMoment=+self%densityNormalization                 &
         &              *self%scaleLength         **(moment+1.0d0) &
         &              *(                                         &
         &                +radialMomentMaximum                     &
         &                -radialMomentMinimum                     &
         &               )    
    return

  contains

    double precision function radialMomentScaleFree(radius)
      !!{
      Provides the scale-free part of the radial moment of the NFW density profile.
      !!}
      use :: Hypergeometric_Functions, only : Hypergeometric_2F1
      use :: Numerical_Comparison    , only : Values_Agree
      implicit none
      double precision, intent(in   ) :: radius

      if (Values_Agree(moment,0.0d0,absTol=1.0d-6)) then
         ! Take the real part of this improper integral. The imaginary parts must cancel when taking differences to compute a
         ! proper integral.
         radialMomentScaleFree=+1.0d0/                 (1.0d0+      radius        ) &
              &                -2.0d0*real(atanh(dcmplx(1.0d0+2.0d0*radius,0.0d0)))
      else if (Values_Agree(moment,1.0d0,absTol=1.0d-6)) then
         radialMomentScaleFree=-1.0d0/                 (1.0d0      +radius        )
      else if (Values_Agree(moment,2.0d0,absTol=1.0d-6)) then
         radialMomentScaleFree=+1.0d0/                 (1.0d0      +radius        ) &
              &                +      log              (1.0d0      +radius        )
      else if (Values_Agree(moment,3.0d0,absTol=1.0d-6)) then
         radialMomentScaleFree=+                                    radius          &
              &                -1.0d0/                 (1.0d0      +radius        ) &
              &                -2.0d0*log              (1.0d0      +radius        )
      else
         radialMomentScaleFree=+(1.0d0+radius)**(moment-1.0d0)                                                     &
              &                /moment                                                                             &
              &                /                (moment-1.0d0)                                                     &
              &                *(                                                                                  &
              &                  - moment                                                                          &
              &                  *  Hypergeometric_2F1([1.0d0-moment,-moment],[2.0d0-moment],1.0d0/(1.0d0+radius)) &
              &                  +(1.0d0+radius)                                                                   &
              &                  *(moment-1.0d0)                                                                   &
              &                  *(                                                                                &
              &                    +(radius/(1.0d0+radius))**moment                                                &
              &                    -Hypergeometric_2F1([     -moment,-moment],[1.0d0-moment],1.0d0/(1.0d0+radius)) &
              &                  )                                                                                 &
              &                 )
      end if
      return
    end function radialMomentScaleFree

  end function nfwDensityRadialMoment

  subroutine nfwDescriptor(self,descriptor,includeClass)
    !!{
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class    (massDistributionNFW), intent(inout)           :: self
    type     (inputParameters    ), intent(inout)           :: descriptor
    logical                       , intent(in   ), optional :: includeClass
    character(len=18             )                          :: parameterLabel
    type     (inputParameters    )                          :: parameters

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('massDistribution','NFW')
    parameters=descriptor%subparameters('massDistribution')
    write (parameterLabel,'(e17.10)') self%densityNormalization
    call parameters%addParameter('densityNormalization',trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)') self%scaleLength
    call parameters%addParameter('scaleLength'         ,trim(adjustl(parameterLabel)))
    return
  end subroutine nfwDescriptor

