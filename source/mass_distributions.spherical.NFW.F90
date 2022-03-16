!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
     procedure :: density => nfwDensity
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
    use :: Input_Parameters        , only : inputParameter, inputParameters
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionNFW)                :: self
    type            (inputParameters    ), intent(inout) :: parameters
    double precision                                     :: mass                , scaleLength, &
         &                                                  densityNormalization, concentration, &
         &                                                  virialRadius
    logical                                              :: dimensionless

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
    <conditionalCall>
     <call>self=massDistributionNFW({conditions})</call>
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

  function nfwConstructorInternal(scaleLength,concentration,densityNormalization,mass,virialRadius,dimensionless) result(self)
    !!{
    Internal constructor for ``nfw'' mass distribution class.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionNFW)                          :: self
    double precision                     , intent(in   ), optional :: scaleLength         , concentration, &
         &                                                            densityNormalization, mass         , &
         &                                                            virialRadius
    logical                              , intent(in   ), optional :: dimensionless
    double precision                                               :: r

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
       r=virialRadius/self%scaleLength
       self%densityNormalization=mass/4.0d0/Pi/self%scaleLength**3/(log(1.0d0+r)-r/(1.0d0+r))
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

  double precision function nfwDensity(self,coordinates)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in an NFW mass distribution.
    !!}
    use :: Coordinates, only : assignment(=), coordinate, coordinateSpherical
    implicit none
    class           (massDistributionNFW), intent(inout) :: self
    class           (coordinate         ), intent(in   ) :: coordinates
    type            (coordinateSpherical)                :: position
    double precision                                     :: r

    ! Get position in spherical coordinate system.
    position  = coordinates
    ! Compute the density at this position.
    r         =+position%r                   () &
         &     /self    %scaleLength
    nfwDensity=+self    %densityNormalization   &
         &     /       r                        &
         &     /(1.0d0+r)**2
    return
  end function nfwDensity
