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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
An implementation of the hot halo outflow stripping class using a simple estimate based on the outer radius of the :term:`CGM`.
!!}

  use :: Dark_Matter_Halo_Scales      , only : darkMatterHaloScaleClass
  use :: Hot_Halo_Mass_Distributions  , only : hotHaloMassDistributionClass
  use :: Hot_Halo_Temperature_Profiles, only : hotHaloTemperatureProfileClass

  !![
  <hotHaloOutflowStripping name="hotHaloOutflowStrippingStandard" docformat="rst">
   <description>
   An implementation of the hot halo outflow stripping class using a simple estimate based on the outer radius of the :term:`CGM`,
   following the approach of :cite:t:`font_colours_2008`. For satellite nodes the fraction of outflowing material which is
   stripped from the halo is

   .. math::

      f_\mathrm{outflow, stripped} = \epsilon_\mathrm{strip} \left( 1 - {M_\mathrm{hot}(&lt;r_\mathrm{outer}) \over M_\mathrm{hot}(&lt;r_\mathrm{virial})} \right),

   where :math:`\epsilon_\mathrm{strip}=`\ ``[efficiency]``, :math:`r_\mathrm{outer}` is the current outer radius of the
   :term:`CGM` (limited to be no larger than the virial radius, :math:`r_\mathrm{virial}`), and :math:`M_\mathrm{hot}(&lt;r)` is
   the mass enclosed within radius :math:`r` in the *untruncated* hot halo mass distribution (see
   :galacticus-class:`hotHaloMassDistribution`), i.e. the distribution, with its current shape, that would be present had no
   stripping occurred. Following :cite:t:`font_colours_2008`, the stripped fraction is therefore :math:`\epsilon_\mathrm{strip}`
   times the fraction of that unstripped hot halo lying beyond the stripping radius, growing from zero for an unstripped halo
   toward :math:`\epsilon_\mathrm{strip}` as the halo is stripped down. Only the shape of the distribution matters, so this is
   well-defined even if the hot halo is empty. Note that the shape is computed using the current properties of the node, so for
   hot halo mass distributions whose shape depends on the current hot halo mass (e.g. through the
   :galacticus-class:`hotHaloMassDistributionCoreRadiusGrowing` core radius) it may differ from that of the halo prior to
   stripping. If no hot halo mass distribution is available no stripping occurs. For nodes which are not satellites no stripping
   occurs.
   </description>
  </hotHaloOutflowStripping>
  !!]
  type, extends(hotHaloOutflowStrippingClass) :: hotHaloOutflowStrippingStandard
     !!{RST
     An implementation of the hot halo outflow stripping class using a simple estimate based on the outer radius of the :term:`CGM`.
     !!}
     private
     class           (darkMatterHaloScaleClass      ), pointer :: darkMatterHaloScale_       => null()
     class           (hotHaloMassDistributionClass  ), pointer :: hotHaloMassDistribution_   => null()
     class           (hotHaloTemperatureProfileClass), pointer :: hotHaloTemperatureProfile_ => null()
     double precision                                          :: efficiency
   contains
     final     ::                     standardDestructor
     procedure :: neverStripped    => standardNeverStripped
     procedure :: fractionStripped => standardFractionStripped
  end type hotHaloOutflowStrippingStandard

  interface hotHaloOutflowStrippingStandard
     !!{RST
     Constructors for the :galacticus-class:`hotHaloOutflowStrippingStandard` hot halo outflow stripping class.
     !!}
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface hotHaloOutflowStrippingStandard

contains

  function standardConstructorParameters(parameters) result(self)
    !!{RST
    Default constructor for the ``standard`` hot halo outflow stripping class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (hotHaloOutflowStrippingStandard)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass       ), pointer       :: darkMatterHaloScale_
    class           (hotHaloMassDistributionClass   ), pointer       :: hotHaloMassDistribution_
    class           (hotHaloTemperatureProfileClass ), pointer       :: hotHaloTemperatureProfile_
    double precision                                                 :: efficiency
    
    !![
    <inputParameter docformat="rst">
      <name>efficiency</name>
      <source>parameters</source>
      <defaultValue>0.1d0</defaultValue>
      <description>
      Specifies the efficiency with which outflowing gas is stripped from the hot halo, following the prescription of :cite:author:`font_colours_2008` (:cite:year:`font_colours_2008`; i.e. this is the parameter :math:`\epsilon_\mathrm{strip}` in their eqn. 6).
      </description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"       name="darkMatterHaloScale_"       source="parameters"/>
    <objectBuilder class="hotHaloMassDistribution"   name="hotHaloMassDistribution_"   source="parameters"/>
    <objectBuilder class="hotHaloTemperatureProfile" name="hotHaloTemperatureProfile_" source="parameters"/>
    !!]
    self=hotHaloOutflowStrippingStandard(efficiency,darkMatterHaloScale_,hotHaloMassDistribution_,hotHaloTemperatureProfile_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"      />
    <objectDestructor name="hotHaloMassDistribution_"  />
    <objectDestructor name="hotHaloTemperatureProfile_"/>
    !!]
    return
  end function standardConstructorParameters

  function standardConstructorInternal(efficiency,darkMatterHaloScale_,hotHaloMassDistribution_,hotHaloTemperatureProfile_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`hotHaloOutflowStrippingStandard` hot halo outflow stripping class.
    !!}
    implicit none
    type            (hotHaloOutflowStrippingStandard)                        :: self
    class           (darkMatterHaloScaleClass       ), intent(in   ), target :: darkMatterHaloScale_
    class           (hotHaloMassDistributionClass   ), intent(in   ), target :: hotHaloMassDistribution_
    class           (hotHaloTemperatureProfileClass ), intent(in   ), target :: hotHaloTemperatureProfile_
    double precision                                 , intent(in   )         :: efficiency
    !![
    <constructorAssign variables="efficiency, *darkMatterHaloScale_, *hotHaloMassDistribution_, *hotHaloTemperatureProfile_"/>
    !!]
    
    return
  end function standardConstructorInternal

  subroutine standardDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`hotHaloOutflowStrippingStandard` hot halo outflow stripping class.
    !!}
    implicit none
    type(hotHaloOutflowStrippingStandard), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"      />
    <objectDestructor name="self%hotHaloMassDistribution_"  />
    <objectDestructor name="self%hotHaloTemperatureProfile_"/>
    !!]
    return
  end subroutine standardDestructor

  logical function standardNeverStripped(self,node) result(neverStripped)
    !!{RST
    Indicate if outflowing mass is never stripped in the hot halo.
    !!}
    implicit none
    class(hotHaloOutflowStrippingStandard), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    neverStripped=.not.node%isSatellite()
    return
  end function standardNeverStripped

  double precision function standardFractionStripped(self,node) result(fractionStripped)
    !!{RST
    Return the fraction of outflowing mass stripped in the hot halo.
    !!}
    use :: Mass_Distributions, only : massDistributionClass, kinematicsDistributionClass
    use :: Galacticus_Nodes  , only : nodeComponentHotHalo
    implicit none
    class           (hotHaloOutflowStrippingStandard), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    class           (massDistributionClass          ), pointer       :: massDistribution_
    class           (kinematicsDistributionClass    ), pointer       :: kinematicsDistribution_
    class           (nodeComponentHotHalo           ), pointer       :: hotHalo
    double precision                                                 :: massOuter              , massVirial  , &
         &                                                              radiusOuter            , radiusVirial

    ! For non-satellite galaxies, no stripping occurs.
    fractionStripped=+0.0d0
    if (.not.node%isSatellite()) return
    ! If the outer radius is at or beyond the virial radius, no stripping has occurred.
    hotHalo      => node                     %hotHalo     (    )
    radiusOuter  =  hotHalo                  %outerRadius (    )
    radiusVirial =  self%darkMatterHaloScale_%radiusVirial(node)
    if (radiusOuter >= radiusVirial) return
    ! Get the untruncated hot halo mass distribution - i.e. the distribution as it would be had no stripping occurred. This is
    ! obtained directly from the hot halo mass distribution class (rather than from the node) as it is not the distribution that
    ! is physically present in the node.
    massDistribution_ => self%hotHaloMassDistribution_%get(node,untruncated=.true.)
    if (.not.associated(massDistribution_)) return
    kinematicsDistribution_ => self%hotHaloTemperatureProfile_%get(node)
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    massOuter =massDistribution_%massEnclosedBySphere(radiusOuter )
    massVirial=massDistribution_%massEnclosedBySphere(radiusVirial)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    <objectDestructor name="massDistribution_"      />
    !!]
    if (massVirial > 0.0d0)                  &
         & fractionStripped=+self%efficiency &
         &                  *(               &
         &                    +1.0d0         &
         &                    -massOuter     &
         &                    /massVirial    &
         &                   )
    return
  end function standardFractionStripped
