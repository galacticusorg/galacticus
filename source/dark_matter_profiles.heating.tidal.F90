!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
  A dark matter halo profile heating class which accounts for heating from tidal shocking.
  !!}

  !![
  <darkMatterProfileHeating name="darkMatterProfileHeatingTidal">
   <description>
    A dark matter profile heating model which accounts for heating due to tidal shocking. The model follows the general
    approach of \cite{gnedin_tidal_1999}. The change in the specific energy of particles at radius $r$ in a halo is given by
    $\Delta \epsilon = \Delta \epsilon_1 + \Delta \epsilon_2$, where $\Delta \epsilon_1$, and $\Delta \epsilon_2$ are the first
    and second order perturbations respectively. The first order term is given by $\Delta \epsilon_1 = Q r^2$ where $Q$ is the
    tidal tensor integrated along the orbital path (see, for example, \citealt{taylor_dynamics_2001}), while the second order
    term is given by $\Delta \epsilon_2 = (2/3) f \sigma_\mathrm{rms} (1+\chi_\mathrm{r,v}) \sqrt{\Delta \epsilon_1}$
    \citep[][eqn.~20, see also \protect\citealt{gnedin_self-consistent_1999}; eqn.~18a,b]{gnedin_tidal_1999}. For the particle
    velocity dispersion, $v_\mathrm{rms}$, we use $\sqrt{3} \sigma_\mathrm{r}(r)$, the radial velocity dispersion in the dark
    matter profile scaled to the total velocity dispersion assuming an isotropic velocity distribution. The position-velocity
    correlation function, $\chi_\mathrm{r,v}$, is taken to be a constant given by the parameter {\normalfont \ttfamily
    [correlationVelocityRadius]}. The coefficient, $f=${\normalfont \ttfamily [coefficientSecondOrder]} is introduced to allow
    some freedom to adjust the contribution of the second order term. It is degenerate with the value of $\chi_\mathrm{r,v}$
    but is introduced to allow for possible future promotion of $\chi_\mathrm{r,v}$ from a constant to a function of the dark
    matter profile potential \citep[see, for example,][appendix~B]{gnedin_self-consistent_1999}.
   </description>
  </darkMatterProfileHeating>
  !!]
  type, extends(darkMatterProfileHeatingClass) :: darkMatterProfileHeatingTidal
     !!{
     A dark matter profile heating class which accounts for heating due to tidal shocking.
     !!}
     private
     double precision :: correlationVelocityRadius, coefficientSecondOrder0, &
          &              coefficientSecondOrder1  , coefficientSecondOrder2
   contains
     procedure :: get => tidalGet
  end type darkMatterProfileHeatingTidal

  interface darkMatterProfileHeatingTidal
     !!{
     Constructors for the {\normalfont \ttfamily tidal} dark matter profile heating class.
     !!}
     module procedure tidalConstructorParameters
     module procedure tidalConstructorInternal
  end interface darkMatterProfileHeatingTidal

contains

  function tidalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily tidal} dark matter profile heating scales class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (darkMatterProfileHeatingTidal), target        :: self
    type            (inputParameters              ), intent(inout) :: parameters
    double precision                                               :: correlationVelocityRadius, coefficientSecondOrder0, &
         &                                                            coefficientSecondOrder1  , coefficientSecondOrder2
         
    !![
    <inputParameter>
      <name>coefficientSecondOrder0</name>
      <defaultValue>0.0d0</defaultValue>
      <source>parameters</source>
      <description>The coefficient, $a_0$, appearing in the second-order heating term, $f_2 = a_0 + a_1 \mathrm{d}\log \rho/\mathrm{d} \log r + a_2 (\mathrm{d}\log \rho/\mathrm{d} \log r)^2$.</description>
    </inputParameter>
    <inputParameter>
      <name>coefficientSecondOrder1</name>
      <defaultValue>0.0d0</defaultValue>
      <source>parameters</source>
      <description>The coefficient, $a_1$, appearing in the second-order heating term, $f_2 = a_0 + a_1 \mathrm{d}\log \rho/\mathrm{d} \log r + a_2 (\mathrm{d}\log \rho/\mathrm{d} \log r)^2$.</description>
    </inputParameter>
    <inputParameter>
      <name>coefficientSecondOrder2</name>
      <defaultValue>0.0d0</defaultValue>
      <source>parameters</source>
      <description>The coefficient, $a_2$, appearing in the second-order heating term, $f_2 = a_0 + a_1 \mathrm{d}\log \rho/\mathrm{d} \log r + a_2 (\mathrm{d}\log \rho/\mathrm{d} \log r)^2$.</description>
    </inputParameter>
    <inputParameter>
      <name>correlationVelocityRadius</name>
      <defaultValue>-1.0d0</defaultValue>
      <source>parameters</source>
      <description>The velocity-position correlation function, $\chi_\mathrm{r,v}$, as defined by \cite[][eqn.~B1]{gnedin_self-consistent_1999} which controls the strength of the second order heating term.</description>
    </inputParameter>
    !!]
    self=darkMatterProfileHeatingTidal(coefficientSecondOrder0,coefficientSecondOrder1,coefficientSecondOrder2,correlationVelocityRadius)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function tidalConstructorParameters

  function tidalConstructorInternal(coefficientSecondOrder0,coefficientSecondOrder1,coefficientSecondOrder2,correlationVelocityRadius) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily tidal} dark matter profile heating scales class.
    !!}
    implicit none
    type            (darkMatterProfileHeatingTidal)                :: self
    double precision                               , intent(in   ) :: correlationVelocityRadius, coefficientSecondOrder0, &
         &                                                            coefficientSecondOrder1  , coefficientSecondOrder2
    !![
    <constructorAssign variables="coefficientSecondOrder0, coefficientSecondOrder1, coefficientSecondOrder2, correlationVelocityRadius"/>
    !!]
    
    return
  end function tidalConstructorInternal

  function tidalGet(self,node) result(massDistributionHeating_)
    !!{
    Return the dark matter mass distribution heating for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes  , only : nodeComponentSatellite
    use :: Mass_Distributions, only : massDistributionHeatingTidal
    implicit none
    class           (massDistributionHeatingClass ), pointer       :: massDistributionHeating_
    class           (darkMatterProfileHeatingTidal), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    class           (nodeComponentSatellite       ), pointer       :: satellite
    double precision                                               :: heatSpecificNormalized
 
    ! Create the mass distribution.
    allocate(massDistributionHeatingTidal :: massDistributionHeating_)
    select type(massDistributionHeating_)
    type is (massDistributionHeatingTidal)
       satellite              =>      node     %satellite             ()
       heatSpecificNormalized =  max(                                     &
            &                        +0.0d0                             , &
            &                        +satellite%tidalHeatingNormalized()  &
            &                       )
       !![
       <referenceConstruct object="massDistributionHeating_">
	 <constructor>
           massDistributionHeatingTidal(                                                          &amp;
            &amp;                       heatSpecificNormalized   =heatSpecificNormalized        , &amp;
            &amp;                       coefficientSecondOrder0  =self%coefficientSecondOrder0  , &amp;
            &amp;                       coefficientSecondOrder1  =self%coefficientSecondOrder1  , &amp;
            &amp;                       coefficientSecondOrder2  =self%coefficientSecondOrder2  , &amp;
            &amp;                       correlationVelocityRadius=self%correlationVelocityRadius  &amp;
            &amp;                      )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    return
  end function tidalGet

