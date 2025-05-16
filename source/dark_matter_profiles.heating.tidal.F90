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
  A dark matter halo profile heating class which accounts for heating from tidal shocking.
  !!}

  !![
  <darkMatterProfileHeating name="darkMatterProfileHeatingTidal">
   <description>
    A dark matter profile heating class that constructs \refClass{massDistributionHeatingTidal} objects to compute heating due to
    tidal shocks.
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
     Constructors for the \refClass{darkMatterProfileHeatingTidal} dark matter profile heating class.
     !!}
     module procedure tidalConstructorParameters
     module procedure tidalConstructorInternal
  end interface darkMatterProfileHeatingTidal

contains

  function tidalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileHeatingTidal} dark matter profile heating scales class which takes a parameter set as input.
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
    Internal constructor for the \refClass{darkMatterProfileHeatingTidal} dark matter profile heating scales class.
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

