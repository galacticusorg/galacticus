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

!+    Contributions to this file made by:  Stéphane Mangeon, Andrew Benson.

  !!{
  Implements a class for black hole binary initial separation based on the model of \cite{volonteri_assembly_2003}.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <blackHoleBinaryInitialSeparation name="blackHoleBinaryInitialSeparationVolonteri2003">
   <description>
    A black hole binary initial separation class that assumes that the initial separation follows the relationship described in
    \cite{volonteri_assembly_2003}
    \begin{equation}
     r_\mathrm{initial} = { \mathrm{G} (M_{\bullet,1} + M_{\bullet, 2}) \over 2
    \sigma_\mathrm{DM}^2 }
    \end{equation}
    where $M_{\bullet, 1}$ and $M_{\bullet, 2}$ are the masses of the black holes
    and $\sigma_\mathrm{DM}$ is the velocity dispersion of the dark matter, which we
    assume to equal the virial velocity of the dark matter halo.
   </description>
  </blackHoleBinaryInitialSeparation>
  !!]
  type, extends(blackHoleBinaryInitialSeparationClass) :: blackHoleBinaryInitialSeparationVolonteri2003
     !!{
     A black hole binary initial separation class based on the model of \cite{volonteri_assembly_2003}.
     !!}
     private
     class(darkMatterHaloScaleClass), pointer:: darkMatterHaloScale_ => null()
   contains
     final     ::                      volonteri2003Destructor
     procedure :: separationInitial => volonteri2003SeparationInitial
  end type blackHoleBinaryInitialSeparationVolonteri2003

  interface blackHoleBinaryInitialSeparationVolonteri2003
     !!{
     Constructors for the {\normalfont \ttfamily volonteri2003} black hole binary initial radius class.
     !!}
     module procedure volonteri2003ConstructorParameters
     module procedure volonteri2003ConstructorInternal
  end interface blackHoleBinaryInitialSeparationVolonteri2003

contains

  function volonteri2003ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily volonteri2003} black hole binary initial radius class which takes a parameter
    list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (blackHoleBinaryInitialSeparationVolonteri2003)                :: self
    type (inputParameters                              ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass                     ), pointer       :: darkMatterHaloScale_

    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=blackHoleBinaryInitialSeparationVolonteri2003(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function volonteri2003ConstructorParameters

  function volonteri2003ConstructorInternal(darkMatterHaloScale_) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily volonteri2003} black hole binary initial radius class which takes a parameter
    list as input.
    !!}
    implicit none
    type (blackHoleBinaryInitialSeparationVolonteri2003)                        :: self
    class(darkMatterHaloScaleClass                     ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    !!]

    return
  end function volonteri2003ConstructorInternal

  subroutine volonteri2003Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily volonteri2003} black hole binary initial separation class.
    !!}
    implicit none
    type(blackHoleBinaryInitialSeparationVolonteri2003), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine volonteri2003Destructor

  double precision function volonteri2003SeparationInitial(self,node,nodeHost)
    !!{
    Returns an initial separation for binary black holes using the method of \cite{volonteri_assembly_2003}, with the
    assumption that the local velocity dispersion is approximately the dark matter halo virial velocity.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBlackHole        , treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class(blackHoleBinaryInitialSeparationVolonteri2003), intent(inout), target :: self
    type (treeNode                                     ), intent(inout), target :: nodeHost     , node
    class(nodeComponentBlackHole                       ), pointer               :: blackHoleHost, blackHole

    blackHole                      =>  node    %blackHole()
    blackHoleHost                  =>  nodeHost%blackHole()
    volonteri2003SeparationInitial =  +gravitationalConstant_internal                        &
         &                            *(                                                     &
         &                              +blackHole    %mass()                                &
         &                              +blackHoleHost%mass()                                &
         &                             )                                                     &
         &                            /2.0d0                                                 &
         &                            /self%darkMatterHaloScale_%velocityVirial(nodeHost)**2
    return
  end function volonteri2003SeparationInitial
