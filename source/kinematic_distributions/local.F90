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

  !!{RST
  Implementation of a local kinematic distribution class.
  !!}

  !![
  <kinematicsDistribution name="kinematicsDistributionLocal" docformat="rst">
    <description>
    A local kinematic distribution class in which the 1D velocity dispersion is given by

    .. math::

       \sigma_\mathrm{1D}(\mathbf{r}) = \alpha V_\mathrm{c}(r),

    where :math:`r = |\mathbf{r}|` and :math:`V_\mathrm{c}(r)` is the rotation curve. Here, :math:`\alpha=`\ ``[]`` is a parameter.
    </description>
  </kinematicsDistribution>
  !!]
  type, public, extends(kinematicsDistributionClass) :: kinematicsDistributionLocal
     !!{RST
     A local kinematic distribution.
     !!}
     double precision :: alpha
   contains
     procedure :: isCollisional        => localIsCollisional
     procedure :: velocityDispersion1D => localVelocityDispersion1D
  end type kinematicsDistributionLocal

  interface kinematicsDistributionLocal
     !!{RST
     Constructors for the :galacticus-class:`kinematicsDistributionLocal` kinematic distribution class.
     !!}
     module procedure localConstructorParameters
     module procedure localConstructorInternal
  end interface kinematicsDistributionLocal

contains

  function localConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`kinematicsDistributionLocal` kinematic distribution class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (kinematicsDistributionLocal)                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    double precision                                             :: alpha

    !![
    <inputParameter docformat="rst">
      <name>alpha</name>
      <description>
      The parameter :math:`\alpha` in the relation :math:`\sigma_\mathrm{1D}(\mathbf{r}) = \alpha V_\mathrm{c}(r)`.
      </description>
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
    !!{RST
    Constructor for the :galacticus-class:`kinematicsDistributionLocal` kinematic distribution class.
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
    !!{RST
    Return false indicating that the local kinematic distribution represents collisionless particles.
    !!}
    implicit none
    class(kinematicsDistributionLocal), intent(inout) :: self
    
    localIsCollisional=.false.
    return
  end function localIsCollisional

  double precision function localVelocityDispersion1D(self,coordinates,massDistribution_,massDistributionEmbedding)
    !!{RST
    Return the 1D velocity dispersion at the specified ``coordinates`` in an local kinematic distribution.
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
