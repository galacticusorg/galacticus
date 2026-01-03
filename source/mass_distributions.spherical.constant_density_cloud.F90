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
  Implementation of a constant density spherical cloud mass distribution class.
  !!}

  !![
  <massDistribution name="massDistributionConstantDensityCloud">
   <description>A mass distribution class for constant density spherical clouds.</description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionConstantDensityCloud
     !!{
     A constant density spherical cloud mass distribution
     !!}
     double precision :: mass    , radius       , &
          &              density_, radiusSquared
   contains
     procedure :: density               => constantDensityCloudDensity
     procedure :: densityGradientRadial => constantDensityCloudDensityGradientRadial
     procedure :: densityRadialMoment   => constantDensityCloudDensityRadialMoment
     procedure :: massEnclosedBySphere  => constantDensityCloudMassEnclosedBySphere
     procedure :: potentialIsAnalytic   => constantDensityCloudPotentialIsAnalytic
    procedure :: potential              => constantDensityCloudPotential
  end type massDistributionConstantDensityCloud

  interface massDistributionConstantDensityCloud
     !!{
     Constructors for the \refClass{massDistributionConstantDensityCloud} mass distribution class.
     !!}
     module procedure constantDensityCloudConstructorParameters
     module procedure constantDensityCloudConstructorInternal
  end interface massDistributionConstantDensityCloud

contains

  function constantDensityCloudConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionConstantDensityCloud} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionConstantDensityCloud)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    double precision                                                      :: mass         , radius
    type            (varying_string                      )                :: componentType
    type            (varying_string                      )                :: massType

    !![
    <inputParameter>
      <name>mass</name>
      <description>The mass of the cloud.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radius</name>
      <description>The radius of the cloud.</description>
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
    self=massDistributionConstantDensityCloud(mass,radius,enumerationComponentTypeEncode(componentType,includesPrefix=.false.),enumerationMassTypeEncode(massType,includesPrefix=.false.))
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function constantDensityCloudConstructorParameters
  
  function constantDensityCloudConstructorInternal(mass,radius,componentType,massType) result(self)
    !!{
    Constructor for the \refClass{massDistributionConstantDensityCloud} convergence class.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionConstantDensityCloud)                          :: self
    double precision                                      , intent(in   )           :: mass         , radius
    type            (enumerationComponentTypeType        ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType             ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="mass, radius, componentType, massType"/>
    !!]

    self%radiusSquared=+radius**2
    self%density_     =+3.0d0     &
         &             *mass      &
         &             /4.0d0     &
         &             /Pi        &
         &             /radius**3
    return
  end function constantDensityCloudConstructorInternal

  double precision function constantDensityCloudDensity(self,coordinates)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a $\beta$-profile mass distribution.
    !!}
    implicit none
    class(massDistributionConstantDensityCloud), intent(inout) :: self
    class(coordinate                          ), intent(in   ) :: coordinates
    
    if (coordinates%rSphericalSquared() < self%radiusSquared) then
       constantDensityCloudDensity=self%density_
    else
       constantDensityCloudDensity=0.0d0
    end if
    return
  end function constantDensityCloudDensity

  double precision function constantDensityCloudDensityGradientRadial(self,coordinates,logarithmic)
    !!{
    Return the density gradient in the radial direction in a constant density cloud.
    !!}
    implicit none
    class  (massDistributionConstantDensityCloud), intent(inout), target   :: self
    class  (coordinate                          ), intent(in   )           :: coordinates
    logical                                      , intent(in   ), optional :: logarithmic
    !$GLC attributes unused :: self, coordinates, logarithmic
    
    constantDensityCloudDensityGradientRadial=0.0d0
    return
  end function constantDensityCloudDensityGradientRadial

  double precision function constantDensityCloudMassEnclosedBySphere(self,radius)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for a constant density cloud.
    !!}
    implicit none
    class           (massDistributionConstantDensityCloud), intent(inout), target :: self
    double precision                                      , intent(in   )         :: radius

    if (radius > self%radius) then
       constantDensityCloudMassEnclosedBySphere=+self%mass
    else
       constantDensityCloudMassEnclosedBySphere=+  self%mass   &
            &                                   *(             &
            &                                     +     radius &
            &                                     /self%radius &
            &                                    )**3
    end if
    return
  end function constantDensityCloudMassEnclosedBySphere

  logical function constantDensityCloudPotentialIsAnalytic(self) result(isAnalytic)
    !!{
    Return that the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionConstantDensityCloud), intent(inout) :: self

    isAnalytic=.true.
    return
  end function constantDensityCloudPotentialIsAnalytic

  double precision function constantDensityCloudPotential(self,coordinates,status)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} in a constant density cloud.
    !!}
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionConstantDensityCloud), intent(inout), target   :: self
    class           (coordinate                          ), intent(in   )           :: coordinates
    type            (enumerationStructureErrorCodeType   ), intent(  out), optional :: status
    double precision                                                                :: radius

    if (present(status)) status=structureErrorCodeSuccess
    radius=coordinates%rSpherical()
    if (radius > self%radius) then
       constantDensityCloudPotential=-gravitationalConstant_internal &
            &                        *self%mass                      &
            &                        /radius
    else
       constantDensityCloudPotential=-gravitationalConstant_internal &
            &                        *self%mass                      &
            &                        /self%radius                    &
            &                        /2.0d0                          &
            &                        *(                              &
            &                          +1.0d0                        &
            &                          +(                            &
            &                            +     radius                &
            &                            /self%radius                &
            &                           )**2                         &
            &                         )
    end if
    return
  end function constantDensityCloudPotential

  double precision function constantDensityCloudDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite)
    !!{
    Computes radial moments of the density in a constant density cloud.
    !!}
    use :: Numerical_Comparison, only : Values_Agree
    use :: Error               , only : Error_Report
    implicit none
    class           (massDistributionConstantDensityCloud), intent(inout)           :: self
    double precision                                      , intent(in   )           :: moment
    double precision                                      , intent(in   ), optional :: radiusMinimum , radiusMaximum
    logical                                               , intent(  out), optional :: isInfinite
    double precision                                                                :: radiusMaximum_
    
    constantDensityCloudDensityRadialMoment=+0.0d0
    if (present(isInfinite)) isInfinite=.false.
    radiusMaximum_=min(radiusMaximum,self%radius)
    if (radiusMinimum < radiusMaximum_) then
       if (radiusMinimum <= 0.0d0 .and. moment <= -1.0d0) then
          if (present(isInfinite)) then
             isInfinite=.true.
             return
          else
             call Error_Report('radial moment is infinite'//{introspection:location})
          end if
       end if
       if (Values_Agree(moment,-1.0d0,absTol=1.0d-6)) then
          constantDensityCloudDensityRadialMoment=+self%density_                    &
               &                                  *log(                             &
               &                                       +radiusMaximum_              &
               &                                       /radiusMinimum               &
               &                                      )
       else
          constantDensityCloudDensityRadialMoment=+self%density_                    &
               &                                  *(                                &
               &                                    +radiusMaximum_**(1.0d0+moment) &
               &                                    -radiusMinimum **(1.0d0+moment) &
               &                                   )                                &
               &                                  /                  (1.0d0+moment)
       end if
    end if
    return    
  end function constantDensityCloudDensityRadialMoment
