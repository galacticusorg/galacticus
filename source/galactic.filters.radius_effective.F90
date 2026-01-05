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

!!{
Implements a galactic high-pass filter for stellar mass effective radius.
!!}

  !![
  <galacticFilter name="galacticFilterRadiusEffective">
   <description>
   A galactic high-pass filter for stellar mass effective radius. Galaxies with a stellar mass effective radius greater than or
   equal to a fixed threshold, $R_{\mathrm{eff},0}=${\normalfont \ttfamily [radiusEffectiveThreshold]}, are passed.
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterRadiusEffective
     !!{
     A galactic high-pass filter for stellar mass effective radius.
     !!}
     private
     double precision :: radiusThreshold
   contains
     procedure :: passes => radiusEffectivePasses
  end type galacticFilterRadiusEffective

  interface galacticFilterRadiusEffective
     !!{
     Constructors for the \refClass{galacticFilterRadiusEffective} galactic filter class.
     !!}
     module procedure radiusEffectiveConstructorParameters
     module procedure radiusEffectiveConstructorInternal
  end interface galacticFilterRadiusEffective

contains

  function radiusEffectiveConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterRadiusEffective} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterRadiusEffective)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    double precision                                               :: radiusThreshold

    !![
    <inputParameter>
      <name>radiusThreshold</name>
      <source>parameters</source>
      <description>The parameter $R_{\mathrm{eff},0}$ (in units of Mpc) appearing in the stellar mass effective radius threshold.</description>
    </inputParameter>
    !!]
    self=galacticFilterRadiusEffective(radiusThreshold)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function radiusEffectiveConstructorParameters

  function radiusEffectiveConstructorInternal(radiusThreshold) result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterRadiusEffective} galactic filter class.
    !!}
    implicit none
    type            (galacticFilterRadiusEffective)                :: self
    double precision                               , intent(in   ) :: radiusThreshold
    !![
    <constructorAssign variables="radiusThreshold"/>
    !!]
    
    return
  end function radiusEffectiveConstructorInternal

  logical function radiusEffectivePasses(self,node) result(passes)
    !!{
    Implement a stellar mass high-pass galactic filter.
    !!}
    use :: Mass_Distributions        , only : massDistributionClass, massDistributionComposite
    use :: Galactic_Structure_Options, only : massTypeStellar
    implicit none
    class(galacticFilterRadiusEffective), intent(inout)          :: self
    type (treeNode                     ), intent(inout), target  :: node
    class(massDistributionClass        )               , pointer :: massDistribution_

    massDistribution_ =>  node             %massDistribution              (massType      =massTypeStellar)
    passes            =   massDistribution_%radiusCylindricalEnclosingMass(massFractional=0.5d0          ) &
         &               >=                                                                                &
         &                self             %radiusThreshold
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function radiusEffectivePasses
