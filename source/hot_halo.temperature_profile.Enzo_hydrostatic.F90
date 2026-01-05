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
An implementation of the hot halo temperature class which uses the ``hydrostatic'' solution from the Enzo code.
!!}

  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !![
  <hotHaloTemperatureProfile name="hotHaloTemperatureProfileEnzoHydrostatic">
   <description>
    A hot halo temperature profile class that implements the ``hydrostatic'' temperature profile available in the \gls{enzo}
    code. Specifically,
    \begin{equation}
      T(r) = \hbox{max}\left( {\mathrm{G} M(&lt;r) \mu m_\mathrm{H} \over 3 \mathrm{k_B} r} , T_\mathrm{min} \right),
    \end{equation}
    where $M(&lt;r)$ is the total mass enclosed within radius $r$, $\mu$ is the primordial mean atomic mass, and
    $T_\mathrm{min}=100$~K is a temperature floor introduced so as to avoid the temperature reaching arbitrarily low values.
   </description>
  </hotHaloTemperatureProfile>
  !!]
  type, extends(hotHaloTemperatureProfileClass) :: hotHaloTemperatureProfileEnzoHydrostatic
     !!{
     An implementation of the hot halo temperature profile class which uses the ``hydrostatic'' solution from the Enzo code.
     !!}
     private
     class(darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
   contains
     final     ::        enzoHydrostaticDestructor
     procedure :: get => enzoHydrostaticGet
  end type hotHaloTemperatureProfileEnzoHydrostatic

  interface hotHaloTemperatureProfileEnzoHydrostatic
     !!{
     Constructors for the enzoHydrostatic hot halo temperature profile class.
     !!}
     module procedure enzoHydrostaticConstructorParameters
     module procedure enzoHydrostaticConstructorInternal
  end interface hotHaloTemperatureProfileEnzoHydrostatic

  double precision, parameter :: temperatureMinimum=1.0d2

contains

  function enzoHydrostaticConstructorParameters(parameters) result(self)
    !!{
    Constructor for the enzoHydrostatic cooling rate class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (hotHaloTemperatureProfileEnzoHydrostatic)                :: self
    type (inputParameters                         ), intent(inout) :: parameters
    class(darkMatterProfileDMOClass               ), pointer       :: darkMatterProfileDMO_

    !![
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=hotHaloTemperatureProfileEnzoHydrostatic(darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function enzoHydrostaticConstructorParameters

  function enzoHydrostaticConstructorInternal(darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the enzoHydrostatic cooling rate class.
    !!}
    implicit none
    type (hotHaloTemperatureProfileEnzoHydrostatic)                        :: self
    class(darkMatterProfileDMOClass               ), intent(in   ), target :: darkMatterProfileDMO_
    !![
    <constructorAssign variables="*darkMatterProfileDMO_"/>
    !!]

    return
  end function enzoHydrostaticConstructorInternal

  subroutine enzoHydrostaticDestructor(self)
    !!{
    Destructor for the \refClass{hotHaloTemperatureProfileEnzoHydrostatic} hot halo temperature profile class.
    !!}
    implicit none
    type(hotHaloTemperatureProfileEnzoHydrostatic), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine enzoHydrostaticDestructor

  function enzoHydrostaticGet(self,node) result(kinematicsDistribution_)
    !!{
    Return the virial hot halo temperature distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Mass_Distributions, only : massDistributionClass, kinematicsDistributionEnzoHydrostatic
    implicit none
    class(kinematicsDistributionClass             ), pointer       :: kinematicsDistribution_
    class(hotHaloTemperatureProfileEnzoHydrostatic), intent(inout) :: self
    type (treeNode                                ), intent(inout) :: node
    class(massDistributionClass                   ), pointer       :: massDistribution_
    
    ! Create an isothermal kinematics distribution.
    allocate(kinematicsDistributionEnzoHydrostatic :: kinematicsDistribution_)
    select type(kinematicsDistribution_)
    type is (kinematicsDistributionEnzoHydrostatic)
       massDistribution_ => self%darkMatterProfileDMO_%get(node)
       !![
       <referenceConstruct object="kinematicsDistribution_">
	 <constructor>
           kinematicsDistributionEnzoHydrostatic(massDistribution_=massDistribution_)
	 </constructor>
       </referenceConstruct>
       <objectDestructor name="massDistribution_"/>
       !!]
    end select
    return
  end function enzoHydrostaticGet
