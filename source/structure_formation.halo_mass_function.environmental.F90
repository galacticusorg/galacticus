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
Implements a dark matter halo mass function class which handles the transition through the environment
mass scale.
!!}


  !![
  <haloMassFunction name="haloMassFunctionEnvironmental">
   <description>
    The halo mass function is computed by handling the transition though the environment mass scale.
   </description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionEnvironmentAveraged) :: haloMassFunctionEnvironmental
     !!{
     A halo mass function class which implements a dark matter halo mass function class which handles the transition through
     the environment mass scale.
     !!}
     private
   contains
     procedure :: differential => environmentalDifferential
  end type haloMassFunctionEnvironmental

  interface haloMassFunctionEnvironmental
     !!{
     Constructors for the {\normalfont \ttfamily environmental} halo mass function class.
     !!}
     module procedure environmentalConstructorParameters
     module procedure environmentalConstructorInternal
  end interface haloMassFunctionEnvironmental

contains

  function environmentalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily environmental} halo mass function class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (haloMassFunctionEnvironmental)                :: self
    type (inputParameters              ), intent(inout) :: parameters
    class(haloMassFunctionClass        ), pointer       :: haloMassFunctionConditioned_, haloMassFunctionUnconditioned_
    class(haloEnvironmentClass         ), pointer       :: haloEnvironment_
    class(cosmologyParametersClass     ), pointer       :: cosmologyParameters_

    !![
    <objectBuilder class="haloMassFunction"    name="haloMassFunctionConditioned_"   source="parameters" parameterName="haloMassFunctionConditioned"  />
    <objectBuilder class="haloMassFunction"    name="haloMassFunctionUnconditioned_" source="parameters" parameterName="haloMassFunctionUnconditioned"/>
    <objectBuilder class="haloEnvironment"     name="haloEnvironment_"               source="parameters"                                              />
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_"           source="parameters"                                              />
    !!]
    self=haloMassFunctionEnvironmental(haloMassFunctionConditioned_,haloMassFunctionUnconditioned_,haloEnvironment_,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="haloMassFunctionConditioned_"  />
    <objectDestructor name="haloMassFunctionUnconditioned_"/>
    <objectDestructor name="haloEnvironment_"              />
    <objectDestructor name="cosmologyParameters_"          />
    !!]
    return
  end function environmentalConstructorParameters

  function environmentalConstructorInternal(haloMassFunctionConditioned_,haloMassFunctionUnconditioned_,haloEnvironment_,cosmologyParameters_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily environmental} halo mass function class.
    !!}
    implicit none
    type (haloMassFunctionEnvironmental)                        :: self
    class(haloMassFunctionClass        ), target, intent(in   ) :: haloMassFunctionConditioned_,haloMassFunctionUnconditioned_
    class(haloEnvironmentClass         ), target, intent(in   ) :: haloEnvironment_
    class(cosmologyParametersClass     ), target, intent(in   ) :: cosmologyParameters_
    !![
    <constructorAssign variables="*haloMassFunctionConditioned_, *haloMassFunctionUnconditioned_, *haloEnvironment_, *cosmologyParameters_"/>
    !!]

    self%timeMatching=-1.0d0
    return
  end function environmentalConstructorInternal

  double precision function environmentalDifferential(self,time,mass,node)
    !!{
    Return the differential halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionEnvironmental), intent(inout), target   :: self
    double precision                               , intent(in   )           :: time          , mass
    type            (treeNode                     ), intent(inout), optional :: node
    double precision                                                         :: massBackground

    massBackground=self%haloEnvironment_%environmentMass()
    if (mass >= massBackground) then
       ! If the halo mass is equal to or greater than the mass of the environment, we simply use the unconditioned mass function.
       environmentalDifferential=+self%haloMassFunctionEnvironmentAveraged%differential(time,mass,node)
    else
       ! Halo mass is less than the mass of the environment. Use the conditioned mass function
       environmentalDifferential=+self%haloMassFunctionConditioned_       %differential(time,mass,node)
    end if
    return
  end function environmentalDifferential
