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
An implementation of the hot halo outflow reincorporation class in which implements the model of
\cite{henriques_simulations_2013}.
!!}

  use :: Cosmology_Functions    , only : cosmologyFunctions , cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScale, darkMatterHaloScaleClass

  !![
  <hotHaloOutflowReincorporation name="hotHaloOutflowReincorporationHenriques2013">
   <description>
    An implementation of the hot halo outflow reincorporation class which implements the model of
    \cite{henriques_simulations_2013}. Specifically, outflowed gas is returned at a rate:
    \begin{equation}
     \dot{M}_\mathrm{return} = \gamma (1+z)^{-\delta_1} \left({V_\mathrm{vir}\over 200\hbox{km/s}}\right)^{\delta_2} {M_\mathrm{outflowed} \over \tau_\mathrm{dyn}}
    \end{equation}
   </description>
  </hotHaloOutflowReincorporation>
  !!]
  type, extends(hotHaloOutflowReincorporationClass) :: hotHaloOutflowReincorporationHenriques2013
     !!{
     An implementation of the hot halo outflow reincorporation class which implements the model of \cite{henriques_simulations_2013}.
     !!}
     private
     double precision                                    :: gamma                         , delta1, &
          &                                                 delta2
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     final     ::         henriques2013Destructor
     procedure :: rate => henriques2013Rate
  end type hotHaloOutflowReincorporationHenriques2013

  interface hotHaloOutflowReincorporationHenriques2013
     !!{
     Constructors for the \refClass{hotHaloOutflowReincorporationHenriques2013} hot halo outflow reincorporation class.
     !!}
     module procedure henriques2013ConstructorParameters
     module procedure henriques2013ConstructorInternal
  end interface hotHaloOutflowReincorporationHenriques2013

contains

  function henriques2013ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{hotHaloOutflowReincorporationHenriques2013} hot halo outflow reincorporation class which takes a parameter set
    as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (hotHaloOutflowReincorporationHenriques2013)                :: self
    type            (inputParameters                           ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                   ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass                  ), pointer       :: darkMatterHaloScale_
    double precision                                                            :: gamma               , delta1, &
          &                                                                        delta2

    !![
    <inputParameter>
      <name>gamma</name>
      <defaultValue>5.0d0</defaultValue>
      <description>The dimensionless parameter $\gamma$ which multiplier the rate at which reheated mass is returned to the hot phase.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>delta1</name>
      <defaultValue>2.40d0</defaultValue>
      <defaultSource>\cite{henriques_simulations_2013}</defaultSource>
      <description>The exponent of the $(1+z)$ term, $\delta_1$.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>delta2</name>
      <defaultValue>3.07d0</defaultValue>
      <defaultSource>\cite{henriques_simulations_2013}</defaultSource>
      <description>The exponent of the $V_\mathrm{vir}$ term, $\delta_2$.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !!]
    self=hotHaloOutflowReincorporationHenriques2013(gamma,delta1,delta2,cosmologyFunctions_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function henriques2013ConstructorParameters

  function henriques2013ConstructorInternal(gamma,delta1,delta2,cosmologyFunctions_,darkMatterHaloScale_) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily henriques2013} hot halo outflow reincorporation class.
    !!}
    implicit none
    type            (hotHaloOutflowReincorporationHenriques2013)                        :: self
    double precision                                            , intent(in   )         :: gamma               , delta1, &
          &                                                                                delta2
    class           (cosmologyFunctionsClass                   ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass                  ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="gamma, delta1, delta2, *cosmologyFunctions_, *darkMatterHaloScale_"/>
    !!]

    return
  end function henriques2013ConstructorInternal

  subroutine henriques2013Destructor(self)
    !!{
    Destructor for the \refClass{hotHaloOutflowReincorporationHenriques2013} hot halo outflow reincorporation class.
    !!}
    implicit none
    type(hotHaloOutflowReincorporationHenriques2013), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine henriques2013Destructor

  double precision function henriques2013Rate(self,node)
    !!{
    Return the rate of mass reincorporation for outflowed gas in the hot halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentHotHalo, treeNode
    implicit none
    class           (hotHaloOutflowReincorporationHenriques2013), intent(inout) :: self
    type            (treeNode                                  ), intent(inout) :: node
    class           (nodeComponentBasic                        ), pointer       :: basic
    class           (nodeComponentHotHalo                      ), pointer       :: hotHalo
    double precision                                            , parameter     :: velocityNormalization=200.0d0 ! km/s.

    basic             =>  node   %basic                                    (            )
    hotHalo           =>  node   %hotHalo                                  (            )
    henriques2013Rate =  +hotHalo%outflowedMass                            (            )              &
         &               *  self%gamma                                                                 &
         &               *  self%cosmologyFunctions_ %expansionFactor      (basic%time())**self%delta1 &
         &               *(                                                                            &
         &                 +self%darkMatterHaloScale_%velocityVirial       (node        )              &
         &                 /                          velocityNormalization                            &
         &                )                                                              **self%delta2 &
         &               /  self%darkMatterHaloScale_%timescaleDynamical   (node        )
    return
  end function henriques2013Rate
