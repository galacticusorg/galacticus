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
An implementation of the hot halo outflow stripping class using a simple estimate based on the outer radius of the \gls{cgm}.
!!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <hotHaloOutflowStripping name="hotHaloOutflowStrippingStandard">
   <description>An implementation of the hot halo outflow stripping class using a simple estimate based on the outer radius of the \gls{cgm}.</description>
  </hotHaloOutflowStripping>
  !!]
  type, extends(hotHaloOutflowStrippingClass) :: hotHaloOutflowStrippingStandard
     !!{
     An implementation of the hot halo outflow stripping class using a simple estimate based on the outer radius of the \gls{cgm}.
     !!}
     private
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     double precision                                    :: efficiency
   contains
     final     ::                     standardDestructor
     procedure :: neverStripped    => standardNeverStripped
     procedure :: fractionStripped => standardFractionStripped
  end type hotHaloOutflowStrippingStandard

  interface hotHaloOutflowStrippingStandard
     !!{
     Constructors for the \refClass{hotHaloOutflowStrippingStandard} hot halo outflow stripping class.
     !!}
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface hotHaloOutflowStrippingStandard

contains

  function standardConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily standard} hot halo outflow stripping class which takes a parameter set
    as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (hotHaloOutflowStrippingStandard)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass       ), pointer       :: darkMatterHaloScale_
    double precision                                                 :: efficiency
    
    !![
    <inputParameter>
      <name>efficiency</name>
      <source>parameters</source>
      <defaultValue>0.1d0</defaultValue>
      <description>Specifies the efficiency with which outflowing gas is stripped from the hot halo, following the prescription of \citeauthor{font_colours_2008}~(\citeyear{font_colours_2008}; i.e. this is the parameter $\epsilon_\mathrm{strip}$ in their eqn.~6).</description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=hotHaloOutflowStrippingStandard(efficiency,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function standardConstructorParameters

  function standardConstructorInternal(efficiency,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{hotHaloOutflowStrippingStandard} class.
    !!}
    implicit none
    type            (hotHaloOutflowStrippingStandard)                        :: self
    class           (darkMatterHaloScaleClass       ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                 , intent(in   )         :: efficiency
    !![
    <constructorAssign variables="efficiency, *darkMatterHaloScale_"/>
    !!]
    
    return
  end function standardConstructorInternal

  subroutine standardDestructor(self)
    !!{
    Destructor for the \refClass{hotHaloOutflowStrippingStandard} class.
    !!}
    implicit none
    type(hotHaloOutflowStrippingStandard), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine standardDestructor

  logical function standardNeverStripped(self,node) result(neverStripped)
    !!{
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
    !!{
    Return the fraction of outflowing mass stripped in the hot halo.
    !!}
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Galactic_Structure_Options, only : componentTypeHotHalo , massTypeGaseous
    use :: Galacticus_Nodes          , only : nodeComponentHotHalo
    implicit none
    class           (hotHaloOutflowStrippingStandard), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    class           (massDistributionClass          ), pointer       :: massDistribution_
    class           (nodeComponentHotHalo           ), pointer       :: hotHalo
    double precision                                                 :: massOuter        , massVirial  , &
         &                                                              radiusOuter      , radiusVirial
    
    if (node%isSatellite()) then
       hotHalo           => node                                  %hotHalo             (                                    )
       massDistribution_ => node                                  %massDistribution    (componentTypeHotHalo,massTypeGaseous)
       radiusOuter       =  hotHalo                               %outerRadius         (                                    )
       radiusVirial      =  self             %darkMatterHaloScale_%radiusVirial        (node                                )
       massOuter         =  massDistribution_                     %massEnclosedBySphere(radiusOuter                         )
       massVirial        =  massDistribution_                     %massEnclosedBySphere(radiusVirial                        )
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
       if (massVirial > 0.0d0) then
          fractionStripped=+self%efficiency &
               &           *(               &
               &             +1.0d0         &
               &             -massOuter     &
               &             /massVirial    &
               &           )
       else
          fractionStripped=+self%efficiency
       end if
    else
       ! For non-satellite galaxies, no stripping occurs.
       fractionStripped   =+0.0d0
    end if
    return
  end function standardFractionStripped
