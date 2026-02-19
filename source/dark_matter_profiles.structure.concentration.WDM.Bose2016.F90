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

!+    Contributions to this file made by:  Xiaolong Du, Andrew Benson.

  !!{
  An implementation of warm dark matter halo profile concentrations using the
  \cite{bose_copernicus_2016} modifier.
  !!}

  use :: Transfer_Functions , only : transferFunctionClass
  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <darkMatterProfileConcentration name="darkMatterProfileConcentrationWDMBose2016">
   <description>
    A dark matter profile concentration class in which the concentration is computed by applying the correction factor of
    \cite{bose_copernicus_2016}:
    \begin{equation}
    c_\mathrm{WDM} = c_\mathrm{CDM} \left( 1 + \gamma_1 {M_\mathrm{1/2} \over M_\mathrm{halo}} \right)^{-\gamma_2} (1+z)^{\beta(z)},
    \end{equation}
    where $\gamma_1=60$, $\gamma_2=0.17$, $M_\mathrm{1/2}$ is the mass corresponding to the wavenumber at which the WDM transfer
    function is suppressed below the CDM transfer function by a factor of 2, $M_\mathrm{halo}$ is the mass of the dark matter
    halo, and $\beta(z)=0.026 z-0.04$, to a CDM concentration algorithm as specified by {\normalfont \ttfamily [cdmConcentration]}.
   </description>
  </darkMatterProfileConcentration>
  !!]
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationWDMBose2016
     !!{
     A dark matter halo profile concentration class implementing the modifier of
     \cite{bose_copernicus_2016}.
     !!}
     private
     class(darkMatterProfileConcentrationClass), pointer :: cdmConcentration    => null()
     class(transferFunctionClass              ), pointer :: transferFunction_   => null()
     class(cosmologyFunctionsClass            ), pointer :: cosmologyFunctions_ => null()
   contains
     final     ::                                   wdmBose2016Destructor
     procedure :: concentration                  => wdmBose2016Concentration
     procedure :: densityContrastDefinition      => wdmBose2016DensityContrastDefinition
     procedure :: darkMatterProfileDMODefinition => wdmBose2016DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationWDMBose2016

  interface darkMatterProfileConcentrationWDMBose2016
     !!{
     Constructors for the \refClass{darkMatterProfileConcentrationWDMBose2016} dark matter halo profile concentration
     class.
     !!}
     module procedure wdmBose2016ConstructorParameters
     module procedure wdmBose2016ConstructorInternal
  end interface darkMatterProfileConcentrationWDMBose2016

contains

  function wdmBose2016ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily wdm} dark matter halo profile concentration class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterProfileConcentrationWDMBose2016)                :: self
    type (inputParameters                          ), intent(inout) :: parameters
    class(darkMatterProfileConcentrationClass      ), pointer       :: cdmConcentration
    class(transferFunctionClass                    ), pointer       :: transferFunction_
    class(cosmologyFunctionsClass                  ), pointer       :: cosmologyFunctions_

    !![
    <objectBuilder class="darkMatterProfileConcentration" name="cdmConcentration"     source="parameters"/>
    <objectBuilder class="transferFunction"               name="transferFunction_"    source="parameters"/>
    <objectBuilder class="cosmologyFunctions"             name="cosmologyFunctions_"  source="parameters"/>
    !!]
    self=darkMatterProfileConcentrationWDMBose2016(cdmConcentration,transferFunction_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cdmConcentration"   />
    <objectDestructor name="transferFunction_"  />
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function wdmBose2016ConstructorParameters

  function wdmBose2016ConstructorInternal(cdmConcentration,transferFunction_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileConcentrationWDMBose2016} dark matter halo concentration class.
    !!}
    implicit none
    type (darkMatterProfileConcentrationWDMBose2016)                        :: self
    class(darkMatterProfileConcentrationClass      ), intent(in   ), target :: cdmConcentration
    class(transferFunctionClass                    ), intent(in   ), target :: transferFunction_
    class(cosmologyFunctionsClass                  ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="*cdmConcentration, *transferFunction_, *cosmologyFunctions_"/>
    !!]

    return
  end function wdmBose2016ConstructorInternal

  subroutine wdmBose2016Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileConcentrationWDMBose2016} dark matter halo profile concentration class.
    !!}
    implicit none
    type(darkMatterProfileConcentrationWDMBose2016), intent(inout) :: self

    !![
    <objectDestructor name="self%cdmConcentration"   />
    <objectDestructor name="self%transferFunction_"  />
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine wdmBose2016Destructor

  double precision function wdmBose2016Concentration(self,node)
    !!{
    Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node}
    using the warm dark matter modifier of \cite{bose_copernicus_2016}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (darkMatterProfileConcentrationWDMBose2016), intent(inout), target  :: self
    type            (treeNode                                 ), intent(inout), target  :: node
    class           (nodeComponentBasic                       )               , pointer :: basic
    ! Parameters of Bose et al. (2016)'s fitting formula.
    double precision                                           , parameter              :: gamma1  =60.0d0, gamma2=0.17d0
    double precision                                                                    :: redshift

    ! Get required objects.
    basic    => node %basic()
    redshift =  1.0d0/self%cosmologyFunctions_%expansionFactor(basic%time())-1.0d0
    ! Get WDM concentration
    wdmBose2016Concentration=+self%cdmConcentration%concentration(node)              &
         &                   /(                                                      &
         &                     +1.0d0                                                &
         &                     +gamma1                                               &
         &                     *(self%transferFunction_%halfModeMass()/basic%mass()) &
         &                    )**gamma2                                              &
         &                   *(                                                      &
         &                     1.0d0+redshift                                        &
         &                    )**(0.026d0*redshift-0.04d0)
    return
  end function wdmBose2016Concentration

  function wdmBose2016DensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of
    concentration in the warm dark matter modifier of \cite{bose_copernicus_2016}.
    !!}
    implicit none
    class(virialDensityContrastClass               ), pointer       :: wdmBose2016DensityContrastDefinition
    class(darkMatterProfileConcentrationWDMBose2016), intent(inout) :: self

    wdmBose2016DensityContrastDefinition => self%cdmConcentration%densityContrastDefinition()
    return
  end function wdmBose2016DensityContrastDefinition

  function wdmBose2016DarkMatterProfileDefinition(self)
    !!{
    Return a dark matter density profile object defining that used in the definition of concentration in the
    warm dark matter modifier of \cite{bose_copernicus_2016}.
    !!}
    implicit none
    class(darkMatterProfileDMOClass                ), pointer       :: wdmBose2016DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationWDMBose2016), intent(inout) :: self

    wdmBose2016DarkMatterProfileDefinition => self%cdmConcentration%darkMatterProfileDMODefinition()
    return
  end function wdmBose2016DarkMatterProfileDefinition
