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
  Implementation of a satellite dynamical friction class which scales the acceleration returned by another satellite dynamical friction class.
  !!}

  !![
  <satelliteDynamicalFriction name="satelliteDynamicalFrictionScale" docformat="rst">
   <description>
    A satellite dynamical friction class which multiplies the acceleration returned by another satellite dynamical friction class by a constant scale factor.
   </description>
  </satelliteDynamicalFriction>
  !!]
  type, extends(satelliteDynamicalFrictionClass) :: satelliteDynamicalFrictionScale
     !!{RST
     Implementation of a satellite dynamical friction class which multiplies the acceleration returned by another satellite dynamical friction class by a constant scale factor.
     !!}
     private
     class           (satelliteDynamicalFrictionClass), pointer :: satelliteDynamicalFriction_ => null()
     double precision                                           :: scaleFactor
   contains
     final     ::                 scaleDestructor
     procedure :: acceleration => scaleAcceleration
  end type satelliteDynamicalFrictionScale

  interface satelliteDynamicalFrictionScale
     !!{RST
     Constructors for the scale satellite dynamical friction class.
     !!}
     module procedure scaleConstructorParameters
     module procedure scaleConstructorInternal
  end interface satelliteDynamicalFrictionScale

  ! Submodule-scope variables used in root finding.
  class(satelliteDynamicalFrictionScale), pointer :: self_
  !$omp threadprivate(self_)
  
contains

  function scaleConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`satelliteDynamicalFrictionScale` satellite dynamical friction class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (satelliteDynamicalFrictionScale)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(satelliteDynamicalFrictionClass), pointer       :: satelliteDynamicalFriction_
    double precision                                      :: scaleFactor

    !![
    <inputParameter docformat="rst">
      <name>scaleFactor</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The multiplicative scale factor applied to the dynamical friction acceleration. Values less than unity reduce the strength of dynamical friction, while values greater than unity increase it.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="satelliteDynamicalFriction" name="satelliteDynamicalFriction_" source="parameters"/>
    !!]
    self=satelliteDynamicalFrictionScale(satelliteDynamicalFriction_,scaleFactor)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteDynamicalFriction_"/>
    !!]
    return
  end function scaleConstructorParameters

  function scaleConstructorInternal(satelliteDynamicalFriction_,scaleFactor) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`satelliteDynamicalFrictionScale` satellite dynamical friction class.
    !!}
    use :: Table_Labels, only : extrapolationTypeExtrapolate
    implicit none
    type (satelliteDynamicalFrictionScale)                     :: self
    class(satelliteDynamicalFrictionClass), intent(in), target :: satelliteDynamicalFriction_
    double precision                      , intent(in)         :: scaleFactor
    !![
    <constructorAssign variables="*satelliteDynamicalFriction_, scaleFactor"/>
    !!]
    return
  end function scaleConstructorInternal

  subroutine scaleDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`satelliteDynamicalFrictionScale` satellite dynamical friction class.
    !!}
    implicit none
    type(satelliteDynamicalFrictionScale), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteDynamicalFriction_"/>
    !!]
    return
  end subroutine scaleDestructor

  function scaleAcceleration(self,node)
    !!{RST
    Return the acceleration due to dynamical friction for a satellite by scaling the acceleration returned by the child dynamical friction class.
    !!}
    implicit none
    double precision                                 , dimension(3)          :: scaleAcceleration
    class           (satelliteDynamicalFrictionScale), intent(inout), target :: self
    type            (treeNode                       ), intent(inout)         :: node
    
    ! Compute the base acceleration.    
    scaleAcceleration=+self%satelliteDynamicalFriction_%acceleration(node)
    ! Evaluate the dynamical friction acceleration.
    scaleAcceleration=+     scaleAcceleration &
         &            *self%scaleFactor
    return
  end function scaleAcceleration
