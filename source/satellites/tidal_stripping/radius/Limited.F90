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

!+    Contributions to this file made by:  Yu Zhao.

  !!{RST
  Implementation of a satellite tidal radius class which follows the method of :cite:t:`king_structure_1962`.
  !!}
  
  !![
  <satelliteTidalStrippingRadius name="satelliteTidalStrippingRadiusLimited" docformat="rst">
   <description>
   A satellite tidal radius class which uses the method of :cite:t:`king_structure_1962`:

   .. math::

      r_\mathrm{tidal}=\left(\frac{GM_\mathrm{sat}}{\gamma_\mathrm{c} \omega^2-d^2\Phi/dr^2}\right)^{1/3},

   where :math:`\omega` is the orbital angular velocity of the satellite, :math:`\Phi(r)` is the gravitational potential due to the host, and :math:`\gamma_\mathrm{c}=`\ ``[efficiencyCentrifugal]`` is the a model parameter that controls the efficiency of centrifugal force. The calculation is based on the dark matter only density profile of the satellite---no accounting is made for the baryonic components.
   </description>
  </satelliteTidalStrippingRadius>
  !!]
  type, extends(satelliteTidalStrippingRadiusClass) :: satelliteTidalStrippingRadiusLimited
     !!{RST
     Implementation of a satellite tidal radius class which follows the method of :cite:t:`king_structure_1962`.
     !!}
     private
     class  (satelliteTidalStrippingRadiusClass), pointer :: satelliteTidalStrippingRadius_ => null()
     integer                                              :: radiusSolitonID
   contains
     final     ::           limitedDestructor
     procedure :: radius => limitedRadius
  end type satelliteTidalStrippingRadiusLimited

  interface satelliteTidalStrippingRadiusLimited
     !!{RST
     Constructors for the :galacticus-class:`satelliteTidalStrippingRadiusLimited` satellite tidal stripping class.
     !!}
     module procedure limitedConstructorParameters
     module procedure limitedConstructorInternal
  end interface satelliteTidalStrippingRadiusLimited

contains

  function limitedConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`satelliteTidalStrippingRadiusLimited` satellite tidal stripping class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (satelliteTidalStrippingRadiusLimited)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (satelliteTidalStrippingRadiusClass  ), pointer       :: satelliteTidalStrippingRadius_

    !![
    <objectBuilder class="satelliteTidalStrippingRadius" name="satelliteTidalStrippingRadius_" source="parameters"/>
    !!]
    self=satelliteTidalStrippingRadiusLimited(satelliteTidalStrippingRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteTidalStrippingRadius_"/>
    !!]
    return
  end function limitedConstructorParameters

  function limitedConstructorInternal(satelliteTidalStrippingRadius_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`satelliteTidalStrippingRadiusLimited` satellite tidal stripping class.
    !!}
    implicit none
    type            (satelliteTidalStrippingRadiusLimited)                     :: self
    class           (satelliteTidalStrippingRadiusClass  ), intent(in), target :: satelliteTidalStrippingRadius_
    !![
    <constructorAssign variables="*satelliteTidalStrippingRadius_"/>
    <addMetaProperty component="darkMatterProfile" name="solitonRadiusSoliton"  id="self%radiusSolitonID"  isEvolvable="no"  isCreator="no"/>
    !!]
    return
  end function limitedConstructorInternal

  subroutine limitedDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`satelliteTidalStrippingRadiusLimited` satellite tidal stripping class.
    !!}
    implicit none
    type(satelliteTidalStrippingRadiusLimited), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteTidalStrippingRadius_"/>
    !!]
    return
  end subroutine limitedDestructor

  double precision function limitedRadius(self,node)
    !!{RST
    Return the tidal radius using the formulation of :cite:t:`king_structure_1962`.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentDarkMatterProfile, treeNode
    use :: Mass_Distributions              , only : massDistributionClass         , massDistributionSolitonNFWHeated
    implicit none
    class           (satelliteTidalStrippingRadiusLimited), intent(inout), target :: self
    type            (treeNode                            ), intent(inout), target :: node
    class           (nodeComponentDarkMatterProfile      ), pointer               :: darkMatterProfile
    class           (massDistributionClass               ), pointer               :: massDistribution_
    double precision                                                              :: radiusTidal      , radiusSoliton

    massDistribution_  => node%massDistribution()
    darkMatterProfile  => node%darkMatterProfile()

    radiusTidal        =  self%satelliteTidalStrippingRadius_%radius(node)

    select type (massDistribution_)
        type is (massDistributionSolitonNFWHeated)
            radiusSoliton = darkMatterProfile%floatRank0MetaPropertyGet(self%radiusSolitonID)
            limitedRadius = max(radiusTidal, radiusSoliton)
    class default
       limitedRadius = radiusTidal
    end select

    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    
    return
  end function limitedRadius
