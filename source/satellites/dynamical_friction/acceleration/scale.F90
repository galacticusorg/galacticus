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

  !+ Contributions to this file made by:  Anthony Pullen, Andrew Benson.

  !!{
  Implementation of a satellite dynamical friction class which applies the core-stalling model of \cite{kaur_stalling_2018} to another dynamical
  friction class.
  !!}

  !![
  <satelliteDynamicalFriction name="satelliteDynamicalFrictionScale">
   <description>
    A satellite dynamical friction class which applies the core-stalling model of \cite{kaur_stalling_2018} to another dynamical
    friction class. Specifically, the acceleration due to dynamical friction is given by:
    \begin{equation}
     \mathbf{a}_\mathrm{df} = \mathbf{a}^\prime_\mathrm{df} [S_\mathrm{trail}(r/r_*)+S_\mathrm{lead}(r/r_*)],
    \end{equation}
    where $\mathbf{a}^\prime_\mathrm{df}$ is the acceleration provided by the child dynamical friction class, and
    $S_\mathrm{trail}(x)$ and $S_\mathrm{lead}(x)$ are the torque suppression factors defined by
    \cite[][eqn.~79]{kaur_stalling_2018}. The functional forms for $S(x)$ are taken from Figure~8 of
    \cite{kaur_stalling_2018}---extrapolation in $\log(S)$ vs. $x$ is used to extend $S(x)$ to lower $x$ than shown in that figure,
    while for values of $x$ higher than that shown in that figure $S(x)$ is held constant at the maximum value shown.
  
    The stalling radius, $r_*$, is computed using equation~(55) of \cite{kaur_stalling_2018}. For dark matter halo profiles with a
    central cusp (for which that equation has no solution), the acceleration provided by the child dynamical friction class is
    returned unmodified.
   </description>
  </satelliteDynamicalFriction>
  !!]
  type, extends(satelliteDynamicalFrictionClass) :: satelliteDynamicalFrictionScale
     !!{
     Implementation of a satellite dynamical friction class which applies the core-stalling model of \cite{kaur_stalling_2018}
     to another dynamical friction class.
     !!}
     private
     class           (satelliteDynamicalFrictionClass), pointer :: satelliteDynamicalFriction_         => null()
     double precision                                           :: scaleFactor
   contains
     final     ::                 scaleDestructor
     procedure :: acceleration => scaleAcceleration
  end type satelliteDynamicalFrictionScale

  interface satelliteDynamicalFrictionScale
     !!{
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
    !!{
    Constructor for the \refClass{satelliteDynamicalFrictionScale} satellite dynamical friction class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (satelliteDynamicalFrictionScale)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(satelliteDynamicalFrictionClass), pointer       :: satelliteDynamicalFriction_
    double precision                                      :: scaleFactor

    !![
    <inputParameter>
      <name>scaleFactor</name>
      <defaultValue>0.3d0</defaultValue>
      <source>parameters</source>
      <description>The fractional scatter in the solitonic core-halo mass relation (default corresponds to a 50\% fractional scatter).</description>
    </inputParameter>
    <objectBuilder class="satelliteDynamicalFriction" name="satelliteDynamicalFriction_" source="parameters"/>
    !!]
    self=satelliteDynamicalFrictionScale(satelliteDynamicalFriction_, scaleFactor)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteDynamicalFriction_"/>
    !!]
    return
  end function scaleConstructorParameters

  function scaleConstructorInternal(satelliteDynamicalFriction_,scaleFactor) result(self)
    !!{
    Internal constructor for the \refClass{satelliteDynamicalFrictionScale} satellite dynamical friction class.
    !!}
    use :: Table_Labels, only : extrapolationTypeExtrapolate
    implicit none
    type (satelliteDynamicalFrictionScale)                        :: self
    class(satelliteDynamicalFrictionClass), intent(in   ), target :: satelliteDynamicalFriction_
    double precision                      , intent(in)          :: scaleFactor
    !![
    <constructorAssign variables="*satelliteDynamicalFriction_,scaleFactor"/>
    !!]
    return
  end function scaleConstructorInternal

  subroutine scaleDestructor(self)
    !!{
    Destructor for the simple cooling radius class.
    !!}
    implicit none
    type(satelliteDynamicalFrictionScale), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteDynamicalFriction_"/>
    !!]
    return
  end subroutine scaleDestructor

  function scaleAcceleration(self,node)
    !!{
    Return an acceleration for satellites due to dynamical friction using the core-stalling model of \cite{kaur_stalling_2018}.
    !!}
    implicit none
    double precision                                    , dimension(3)          :: scaleAcceleration
    class           (satelliteDynamicalFrictionScale), intent(inout), target :: self
    type            (treeNode                          ), intent(inout)         :: node
    
    ! Compute the base acceleration.    
    scaleAcceleration=+self%satelliteDynamicalFriction_%acceleration(node)
    ! Evaluate the dynamical friction acceleration.
    scaleAcceleration=+scaleAcceleration &
         &               *self%scaleFactor
    return
  end function scaleAcceleration
