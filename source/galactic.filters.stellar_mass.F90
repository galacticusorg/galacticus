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
Contains a module which implements a galactic high-pass filter for total stellar mass.
!!}

  !![
  <galacticFilter name="galacticFilterStellarMass">
   <description>
   A galactic high-pass filter for stellar mass. Galaxies with a combined disk, spheroid, plus \gls{nsc} stellar mass greater than or equal
   to a fixed threshold, $M_{\star,0}=${\normalfont \ttfamily [massThreshold]}.
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterStellarMass
     !!{
     A galactic high-pass filter class for stellar mass.
     !!}
     private
     double precision :: massThreshold
   contains
     procedure :: passes => stellarMassPasses
  end type galacticFilterStellarMass

  interface galacticFilterStellarMass
     !!{
     Constructors for the ``stellarMass'' galactic filter class.
     !!}
     module procedure stellarMassConstructorParameters
     module procedure stellarMassConstructorInternal
  end interface galacticFilterStellarMass

contains

  function stellarMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``stellarMass'' galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterStellarMass)                :: self
    type            (inputParameters          ), intent(inout) :: parameters
    double precision                                           :: massThreshold

    !![
    <inputParameter>
      <name>massThreshold</name>
      <source>parameters</source>
      <description>The parameter $M_0$ (in units of $M_\odot$) appearing in the stellar mass threshold for the stellar mass galactic filter class.</description>
    </inputParameter>
    !!]
    self=galacticFilterStellarMass(massThreshold)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function stellarMassConstructorParameters

  function stellarMassConstructorInternal(massThreshold) result(self)
    !!{
    Internal constructor for the ``stellarMass'' galactic filter class.
    !!}
    implicit none
    type            (galacticFilterStellarMass)                :: self
    double precision                           , intent(in   ) :: massThreshold
    !![
    <constructorAssign variables="massThreshold"/>
    !!]
    
    return
  end function stellarMassConstructorInternal

  logical function stellarMassPasses(self,node) result(passes)
    !!{
    Implement a stellar mass high-pass galactic filter.
    !!}
    use :: Mass_Distributions        , only : massDistributionClass, massDistributionComposite
    use :: Galactic_Structure_Options, only : massTypeStellar
    implicit none
    class(galacticFilterStellarMass), intent(inout)          :: self
    type (treeNode                 ), intent(inout), target  :: node
    class(massDistributionClass    )               , pointer :: massDistribution_

    massDistribution_ =>  node             %massDistribution(massType=massTypeStellar)
    passes            =   massDistribution_%massTotal       (                        ) &
         &               >=                                                            &
         &                self             %massThreshold
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function stellarMassPasses
