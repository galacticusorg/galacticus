!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which implements a galactic high-pass filter for halo mass under a given definition.
!!}

  use :: Virial_Density_Contrast, only : virialDensityContrast, virialDensityContrastClass

  !![
  <galacticFilter name="galacticFilterHaloMass">
   <description>
   A high-pass filter for basic mass. Halos with a halo mass greater than or equal to a fixed threshold,
   $M_0=${\normalfont \ttfamily [massThreshold]}.
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterHaloMass
     !!{
     A galactic high-pass filter class for halo mass.
     !!}
     private
     double precision                                      :: massThreshold
     class           (virialDensityContrastClass), pointer :: virialDensityContrast_ => null()
   contains
     final     ::           haloMassDestructor
     procedure :: passes => haloMassPasses
  end type galacticFilterHaloMass

  interface galacticFilterHaloMass
     !!{
     Constructors for the ``haloMass'' galactic filter class.
     !!}
     module procedure haloMassConstructorParameters
     module procedure haloMassConstructorInternal
  end interface galacticFilterHaloMass

contains

  function haloMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``haloMass'' galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterHaloMass    )                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (virialDensityContrastClass), pointer       :: virialDensityContrast_
    double precision                                            :: massThreshold

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>massThreshold</name>
      <source>parameters</source>
      <description>The parameter $M_0$ (in units of $M_\odot$) appearing in the mass threshold for the halo mass galactic filter class.</description>
    </inputParameter>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    !!]
    self=galacticFilterHaloMass(massThreshold,virialDensityContrast_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="virialDensityContrast_"/>
    !!]
    return
  end function haloMassConstructorParameters

  function haloMassConstructorInternal(massThreshold,virialDensityContrast_) result(self)
    !!{
    Internal constructor for the ``haloMass'' galactic filter class.
    !!}
    implicit none
    type            (galacticFilterHaloMass    )                        :: self
    double precision                            , intent(in   )         :: massThreshold
    class           (virialDensityContrastClass), intent(in   ), target :: virialDensityContrast_
    !![
    <constructorAssign variables="massThreshold, *virialDensityContrast_"/>
    !!]

    return
  end function haloMassConstructorInternal

  subroutine haloMassDestructor(self)
    !!{
    Destructor for the ``haloMass'' galactic filter class.
    !!}
    implicit none
    type(galacticFilterHaloMass), intent(inout) :: self

    !![
    <objectDestructor name="self%virialDensityContrast_"/>
    !!]
    return
  end subroutine haloMassDestructor

  logical function haloMassPasses(self,node)
    !!{
    Implement a halo mass high-pass galactic filter.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    class(galacticFilterHaloMass), intent(inout)         :: self
    type (treeNode              ), intent(inout), target :: node
    class(nodeComponentBasic    ), pointer               :: basic

    basic          =>                                                                                  node %basic()
    haloMassPasses =   Dark_Matter_Profile_Mass_Definition(                                                           &
         &                                                 node                                                     , &
         &                                                 self%virialDensityContrast_%densityContrast(               &
         &                                                                                             basic%mass (), &
         &                                                                                             basic%time ()  &
         &                                                                                            )               &
         &                                                )                                                           &
         &            >=                                                                                              &
         &             self %massThreshold
    return
  end function haloMassPasses
