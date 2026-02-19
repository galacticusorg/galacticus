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
  An implementation of accretion from the \gls{igm} onto halos using simple truncation to
  mimic the effects of reionization, and the Bertschinger mass to define available mass.
  !!}

  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !![
  <accretionHalo name="accretionHaloBertschinger">
   <description>Accretion onto halos using simple truncation to mimic the effects of reionization, and the Bertschinger mass to define available mass.</description>
  </accretionHalo>
  !!]
  type, extends(accretionHaloSimple) :: accretionHaloBertschinger
     !!{
     A halo accretion class using simple truncation to mimic the effects of reionization, and the Bertschinger mass to define
     available mass.
     !!}
     private
     class(darkMatterProfileDMOClass), pointer:: darkMatterProfileDMO_ => null()
   contains
     final     ::                  bertschingerDestructor
     procedure :: velocityScale => bertschingerVelocityScale
  end type accretionHaloBertschinger

  interface accretionHaloBertschinger
     !!{
     Constructors for the \refClass{accretionHaloBertschinger} halo accretion class.
     !!}
     module procedure bertschingerConstructorParameters
     module procedure bertschingerConstructorInternal
  end interface accretionHaloBertschinger

contains

  function bertschingerConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily bertschinger} halo accretion class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (accretionHaloBertschinger)                :: self
    type (inputParameters          ), intent(inout) :: parameters

    self%accretionHaloSimple=accretionHaloSimple(parameters)
    !![
    <objectBuilder class="darkMatterProfileDMO" name="self%darkMatterProfileDMO_" source="parameters"/>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function bertschingerConstructorParameters

  function bertschingerConstructorInternal(timeReionization,velocitySuppressionReionization,accretionNegativeAllowed,accretionNewGrowthOnly,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,accretionHaloTotal_,chemicalState_,intergalacticMediumState_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{accretionHaloBertschinger} halo accretion class.
    !!}
    implicit none
    type            (accretionHaloBertschinger    )                        :: self
    double precision                               , intent(in   )         :: timeReionization        , velocitySuppressionReionization
    logical                                        , intent(in   )         :: accretionNegativeAllowed, accretionNewGrowthOnly
    class           (cosmologyParametersClass     ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass      ), intent(in   ), target :: cosmologyFunctions_
    class           (accretionHaloTotalClass      ), intent(in   ), target :: accretionHaloTotal_
    class           (darkMatterHaloScaleClass     ), intent(in   ), target :: darkMatterHaloScale_
    class           (chemicalStateClass           ), intent(in   ), target :: chemicalState_
    class           (darkMatterProfileDMOClass    ), intent(in   ), target :: darkMatterProfileDMO_
    class           (intergalacticMediumStateClass), intent(in   ), target :: intergalacticMediumState_
    !![
    <constructorAssign variables="*darkMatterProfileDMO_"/>
    !!]

    self%accretionHaloSimple=accretionHaloSimple(timeReionization,velocitySuppressionReionization,accretionNegativeAllowed,accretionNewGrowthOnly,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,accretionHaloTotal_,chemicalState_,intergalacticMediumState_)
    return
  end function bertschingerConstructorInternal

  subroutine bertschingerDestructor(self)
    !!{
    Destructor for the \refClass{accretionHaloBertschinger} halo accretion class.
    !!}
    implicit none
    type(accretionHaloBertschinger), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine bertschingerDestructor

  double precision function bertschingerVelocityScale(self,node) result(velocityScale)
    !!{
    Returns the velocity scale to use for {\normalfont \ttfamily node}. Use the maximum circular velocity.
    !!}
    use :: Mass_Distributions, only : massDistributionClass
    implicit none
    class(accretionHaloBertschinger), intent(inout) :: self
    type (treeNode                 ), intent(inout) :: node
    class(massDistributionClass    ), pointer       :: massDistribution_
    
    massDistribution_ => self             %darkMatterProfileDMO_       %get(node)
    velocityScale     =  massDistribution_%velocityRotationCurveMaximum    (    )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function bertschingerVelocityScale
