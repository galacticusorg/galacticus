!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Implementation of a constant density spherical cloud mass distribution class.

  !# <massDistribution name="massDistributionConstantDensityCloud">
  !#  <description>A mass distribution class for constant density spherical clouds.</description>
  !# </massDistribution>
  type, public, extends(massDistributionSpherical) :: massDistributionConstantDensityCloud
     !% A constant density spherical cloud mass distribution
     double precision :: mass    , radius       , &
          &              density_, radiusSquared
   contains
     procedure :: density               => constantDensityCloudDensity
     procedure :: densityGradientRadial => constantDensityCloudDensityGradientRadial
     procedure :: densityRadialMoment   => constantDensityCloudDensityRadialMoment
     procedure :: massEnclosedBySphere  => constantDensityCloudMassEnclosedBySphere
     procedure :: potential             => constantDensityCloudPotential
  end type massDistributionConstantDensityCloud

  interface massDistributionConstantDensityCloud
     !% Constructors for the {\normalfont \ttfamily constantDensityCloud} mass distribution class.
     module procedure constantDensityCloudConstructorParameters
     module procedure constantDensityCloudConstructorInternal
  end interface massDistributionConstantDensityCloud

contains

  function constantDensityCloudConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily constantDensityCloud} mass distribution class which builds the object from a parameter
    !% set.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (massDistributionConstantDensityCloud)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    double precision                                                      :: mass      , radius

    !# <inputParameter>
    !#   <name>mass</name>
    !#   <description>The mass of the cloud.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>radius</name>
    !#   <description>The radius of the cloud.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    self=massDistributionConstantDensityCloud(mass,radius)
    !# <inputParametersValidate source="parameters"/>
    return
  end function constantDensityCloudConstructorParameters
  
  function constantDensityCloudConstructorInternal(mass,radius) result(self)
    !% Constructor for ``constantDensityCloud'' convergence class.
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionConstantDensityCloud)                :: self
    double precision                                      , intent(in   ) :: mass, radius
    !# <constructorAssign variables="mass, radius"/>

    self%radiusSquared=+radius**2
    self%density_     =+3.0d0     &
         &             *mass      &
         &             /4.0d0     &
         &             /Pi        &
         &             /radius**3
    return
  end function constantDensityCloudConstructorInternal

  double precision function constantDensityCloudDensity(self,coordinates)
    !% Return the density at the specified {\normalfont \ttfamily coordinates} in a $\beta$-profile mass distribution.
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
    !% Return the density gradient in the radial direction in a constant density cloud.
    implicit none
    class  (massDistributionConstantDensityCloud), intent(inout)           :: self
    class  (coordinate                          ), intent(in   )           :: coordinates
    logical                                      , intent(in   ), optional :: logarithmic
    !$GLC attributes unused :: self, coordinates, logarithmic
    
    constantDensityCloudDensityGradientRadial=0.0d0
    return
  end function constantDensityCloudDensityGradientRadial

  double precision function constantDensityCloudMassEnclosedBySphere(self,radius)
    !% Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for a constant density cloud.
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

  double precision function constantDensityCloudPotential(self,coordinates)
    !% Return the potential at the specified {\normalfont \ttfamily coordinates} in a $\beta$-profile mass distribution. Calculated using
    !% \href{http://www.wolframalpha.com/input/?i=integrate+4\%2F3+\%CF\%80+r+\%CF\%81+2F1\%283\%2F2\%2C+\%283+\%CE\%B2\%29\%2F2\%2C+5\%2F2\%2C+-r^2\%29}{Wolfram
    !% Alpha}.
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (massDistributionConstantDensityCloud), intent(inout) :: self
    class           (coordinate                          ), intent(in   ) :: coordinates
    double precision                                                      :: radius

    radius=coordinates%rSpherical()
    if (radius > self%radius) then
       constantDensityCloudPotential=-gravitationalConstantGalacticus &
            &                        *self%mass                       &
            &                        /radius
    else
       constantDensityCloudPotential=-gravitationalConstantGalacticus &
            &                        *self%mass                       &
            &                        /self%radius                     &
            &                        /2.0d0                           &
            &                        *(                               &
            &                          +1.0d0                         &
            &                          +(                             &
            &                            +     radius                 &
            &                            /self%radius                 &
            &                           )**2                          &
            &                         )
    end if
    return
  end function constantDensityCloudPotential

  double precision function constantDensityCloudDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite)
    !% Computes radial moments of the density in a constant density cloud.
    use :: Numerical_Comparison, only : Values_Agree
    use :: Galacticus_Error    , only : Galacticus_Error_Report
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
             call Galacticus_Error_Report('radial moment is infinite'//{introspection:location})
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
