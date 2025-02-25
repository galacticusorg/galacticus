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
  Implementation of a local kinematic distribution class.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionLocal">
    <description>
      A local kinematic distribution class in which the 1D velocity dispersion is given by
      \begin{equation}
      \sigma_\mathrm{1D}(\mathbf{r}) = \alpha V_\mathrm{c}(r),
      \end{equation}      
      where $r = |\mathbf{r}|$ and $V_\mathrm{c}(r)$ is the rotation curve. Here, $\alpha=${\normalfont \ttfamily []} is a
      parameter.
    </description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionClass) :: kinematicsDistributionLocal
     !!{
     A local kinematic distribution.
     !!}
     double precision :: alpha
   contains
     procedure :: isCollisional        => localIsCollisional
     procedure :: velocityDispersion1D => localVelocityDispersion1D
  end type kinematicsDistributionLocal

  interface kinematicsDistributionLocal
     !!{
     Constructors for the {\normalfont \ttfamily local} kinematic distribution class.
     !!}
     module procedure localConstructorParameters
     module procedure localConstructorInternal
  end interface kinematicsDistributionLocal

contains

  function localConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily isothermal} kinematic distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (kinematicsDistributionLocal)                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    double precision                                             :: alpha

    !![
    <inputParameter>
      <name>alpha</name>
      <description>The parameter $\alpha$ in the relation $\sigma_\mathrm{1D}(\mathbf{r}) = \alpha V_\mathrm{c}(r)$.</description>
      <defaultValue>1.0d0/sqrt(2.0d0)</defaultValue>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=kinematicsDistributionLocal(alpha)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function localConstructorParameters
  
  function localConstructorInternal(alpha) result(self)
    !!{
    Constructor for {\normalfont \ttfamily local} kinematics distribution class.
    !!}
    implicit none
    type            (kinematicsDistributionLocal)                :: self
    double precision                             , intent(in   ) :: alpha
    !![
    <constructorAssign variables="alpha"/>
    !!]

    return
  end function localConstructorInternal

  logical function localIsCollisional(self)
    !!{
    Return false indicating that the local kinematic distribution represents collisionless particles.
    !!}
    implicit none
    class(kinematicsDistributionLocal), intent(inout) :: self
    
    localIsCollisional=.false.
    return
  end function localIsCollisional

  double precision function localVelocityDispersion1D(self,coordinates,massDistribution_,massDistributionEmbedding)
    !!{
    Return the 1D velocity dispersion at the specified {\normalfont \ttfamily coordinates} in an local kinematic distribution.
    !!}
    implicit none
    class(kinematicsDistributionLocal), intent(inout)         :: self
    class(coordinate                 ), intent(in   )         :: coordinates
    class(massDistributionClass      ), intent(inout), target :: massDistribution_, massDistributionEmbedding
    !$GLC attributes unused :: massDistribution_
    
    localVelocityDispersion1D=+self                     %alpha                                   &
         &                    *massDistributionEmbedding%rotationCurve(coordinates%rSpherical())
    return
  end function localVelocityDispersion1D
