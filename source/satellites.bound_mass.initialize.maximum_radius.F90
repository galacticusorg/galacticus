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
  Implementation of a satellite bound mass initializor class that sets the initial bound mass by integrating the density profile up to a maximum radius.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <satelliteMassBoundInitializor name="satelliteMassBoundInitializorMaximumRadius" docformat="rst">
   <description>
   A satellite bound mass initializor class that sets the initial bound mass of the satellite halo by integrating its density profile up to a maximum radius, :math:`r_\mathrm{max} =`\ ``[radiusMaximumOverRadiusVirial]`` :math:`\times r_\mathrm{virial}`. The density profile is assumed to be zero beyond this maximum radius.
   </description>
  </satelliteMassBoundInitializor>
  !!]
  type, extends(satelliteMassBoundInitializorClass) :: satelliteMassBoundInitializorMaximumRadius
     !!{RST
     Implementation of a satellite bound mass initializor class that sets the initial bound mass by integrating the density profile up to a maximum radius.
     !!}
     private
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_          => null()
     double precision                                    :: radiusMaximumOverRadiusVirial
   contains
     final     ::              maximumRadiusDestructor
     procedure :: massBound => maximumRadiusMassBound
  end type satelliteMassBoundInitializorMaximumRadius

  interface satelliteMassBoundInitializorMaximumRadius
     !!{RST
     Constructors for the :galacticus-class:`satelliteMassBoundInitializorMaximumRadius` satellite bound mass initializor class.
     !!}
     module procedure maximumRadiusConstructorParameters
     module procedure maximumRadiusConstructorInternal
  end interface satelliteMassBoundInitializorMaximumRadius

contains

  function maximumRadiusConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`satelliteMassBoundInitializorMaximumRadius` satellite bound mass initializor class which builds the object from a parameter set.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (satelliteMassBoundInitializorMaximumRadius)                :: self
    type            (inputParameters                           ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass                  ), pointer       :: darkMatterHaloScale_
    double precision                                                            :: radiusMaximumOverRadiusVirial

    !![
    <inputParameter docformat="rst">
      <name>radiusMaximumOverRadiusVirial</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The maximum radius of the satellite halo in units of its virial radius. This value will be used to compute the initial bound mass of the satellite halo by integrating the density profile up to this maximum radius, assuming that the density profile is zero beyond this radius.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    if (radiusMaximumOverRadiusVirial <= 0.0d0) call Error_Report('specify a positive maximum radius for the satellite'//{introspection:location})
    self=satelliteMassBoundInitializorMaximumRadius(radiusMaximumOverRadiusVirial,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function maximumRadiusConstructorParameters

  function maximumRadiusConstructorInternal(radiusMaximumOverRadiusVirial,darkMatterHaloScale_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`satelliteMassBoundInitializorMaximumRadius` satellite bound mass initializor class.
    !!}
    implicit none
    type            (satelliteMassBoundInitializorMaximumRadius)                        :: self
    class           (darkMatterHaloScaleClass                  ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                            , intent(in   )         :: radiusMaximumOverRadiusVirial
    !![
    <constructorAssign variables="radiusMaximumOverRadiusVirial, *darkMatterHaloScale_"/>
    !!]

    return
  end function maximumRadiusConstructorInternal

  subroutine maximumRadiusDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`satelliteMassBoundInitializorMaximumRadius` satellite bound mass initializor class.
    !!}
    implicit none
    type(satelliteMassBoundInitializorMaximumRadius), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine maximumRadiusDestructor

  double precision function maximumRadiusMassBound(self,node) result(massBound)
    !!{RST
    Returns the initial bound mass of a satellite halo by integrating the density profile up to a maximum radius.
    !!}
    use :: Mass_Distributions, only : massDistributionClass
    implicit none
    class           (satelliteMassBoundInitializorMaximumRadius), intent(inout) :: self
    type            (treeNode                                  ), intent(inout) :: node
    class           (massDistributionClass                     ), pointer       :: massDistribution_
    double precision                                                            :: radiusMaximum

    radiusMaximum     =  +self                                  %radiusMaximumOverRadiusVirial                &
         &               *self             %darkMatterHaloScale_%radiusVirial                 (node         )
    massDistribution_ =>  node                                  %massDistribution             (             )
    massBound         =   massDistribution_                     %massEnclosedBySphere         (radiusMaximum)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function maximumRadiusMassBound

