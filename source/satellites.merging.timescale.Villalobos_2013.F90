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
  Implements a satellite merging timescale class which applies the \cite{villalobos_improved_2013} modifier to another selected
  satellite merging time class.
  !!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <satelliteMergingTimescales name="satelliteMergingTimescalesVillalobos2013">
   <description>
    A satellite merging timescale class which computes merging timescales using the modifier of \cite{villalobos_improved_2013}
    as
    \begin{equation}
    \tau_\mathrm{merge} = (1+z)^\alpha \tau^\prime_\mathrm{merge},
    \end{equation}
    where $\alpha=${\normalfont \ttfamily [exponent]} and $\tau^\prime_\mathrm{merge}$ is the merging timescale computed by
    another satellite merging timescale.
   </description>
  </satelliteMergingTimescales>
  !!]
  type, extends(satelliteMergingTimescalesClass) :: satelliteMergingTimescalesVillalobos2013
     !!{
     A class implementing calculations of satellite merging times by applying the \cite{villalobos_improved_2013} modifier to
     another selected satellite merging time method.
     !!}
     private
     class           (cosmologyFunctionsClass        ), pointer :: cosmologyFunctions_         => null()
     class           (satelliteMergingTimescalesClass), pointer :: satelliteMergingTimescales_ => null()
     double precision                                           :: exponent
   contains
     final     ::                     villalobos2013Destructor
     procedure :: timeUntilMerging => villalobos2013TimeUntilMerging
  end type satelliteMergingTimescalesVillalobos2013

  interface satelliteMergingTimescalesVillalobos2013
     !!{
     Constructors for the \cite{lacey_merger_1993} merging timescale class.
     !!}
     module procedure villalobos2013ConstructorParameters
     module procedure villalobos2013ConstructorInternal
  end interface satelliteMergingTimescalesVillalobos2013

contains

  function villalobos2013ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \cite{lacey_merger_1993} merging timescale class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (satelliteMergingTimescalesVillalobos2013)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (satelliteMergingTimescalesClass         ), pointer       :: satelliteMergingTimescales_
    class           (cosmologyFunctionsClass                 ), pointer       :: cosmologyFunctions_
    double precision                                                          :: exponent

    !![
    <inputParameter>
      <name>exponent</name>
      <defaultSource>\citep{villalobos_improved_2013}</defaultSource>
      <defaultValue>0.44d0</defaultValue>
      <description>The exponent of $1+z$ appearing in the \cite{villalobos_improved_2013} modifier for satellite merging timescales.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="satelliteMergingTimescales" name="satelliteMergingTimescales_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"         name="cosmologyFunctions_"         source="parameters"/>
    !!]
    self=satelliteMergingTimescalesVillalobos2013(exponent,satelliteMergingTimescales_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteMergingTimescales_"/>
    <objectDestructor name="cosmologyFunctions_"        />
    !!]
    return
  end function villalobos2013ConstructorParameters

  function villalobos2013ConstructorInternal(exponent,satelliteMergingTimescales_,cosmologyFunctions_) result(self)
    !!{
    Constructor for the \cite{lacey_merger_1993} merging timescale class.
    !!}
    implicit none
    type            (satelliteMergingTimescalesVillalobos2013)                        :: self
    double precision                                          , intent(in   )         :: exponent
    class           (satelliteMergingTimescalesClass         ), intent(in   ), target :: satelliteMergingTimescales_
    class           (cosmologyFunctionsClass                 ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="exponent, *satelliteMergingTimescales_, *cosmologyFunctions_"/>
    !!]

    return
  end function villalobos2013ConstructorInternal

  subroutine villalobos2013Destructor(self)
    !!{
    Destructor for the \refClass{satelliteMergingTimescalesVillalobos2013} satellite merging timescale class.
    !!}
    implicit none
    type(satelliteMergingTimescalesVillalobos2013), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteMergingTimescales_"/>
    <objectDestructor name="self%cosmologyFunctions_"        />
    !!]
    return
  end subroutine villalobos2013Destructor

  double precision function villalobos2013TimeUntilMerging(self,node,orbit)
    !!{
    Return the timescale for merging satellites using the \cite{villalobos_improved_2013} method.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    use :: Kepler_Orbits   , only : keplerOrbit
    implicit none
    class           (satelliteMergingTimescalesVillalobos2013), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    type            (keplerOrbit                             ), intent(inout) :: orbit
    class           (nodeComponentBasic                      ), pointer       :: basic
    double precision                                                          :: expansionFactor

    ! Compute expansion factor.
    basic           => node                    %basic          (            )
    expansionFactor =  self%cosmologyFunctions_%expansionFactor(basic%time())
    ! Compute dynamical friction timescale.
    villalobos2013TimeUntilMerging=+self%satelliteMergingTimescales_%timeUntilMerging(node,orbit) &
         &                         /expansionFactor**self%exponent
    return
  end function villalobos2013TimeUntilMerging
