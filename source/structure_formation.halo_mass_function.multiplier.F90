!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Contains a module which implements a dark matter halo mass function class which modifies another mass function by multiplying it
by a constant factor.
!!}

   use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <haloMassFunction name="haloMassFunctionMultiplier">
   <description>
    The halo mass function is computed by multiplying another halo mass function by a constant factor.
   </description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionClass) :: haloMassFunctionMultiplier
     !!{
     A halo mass function class that modifies another mass function by multiplying it by a constant factor.
     !!}
     private
     double precision                                   :: multiplier                   , exponentRedshift
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     class           (haloMassFunctionClass  ), pointer :: massFunction_       => null()
   contains
     final     ::                 multiplierDestructor
     procedure :: differential => multiplierDifferential
     procedure :: integrated   => multiplierIntegrated
  end type haloMassFunctionMultiplier

  interface haloMassFunctionMultiplier
     !!{
     Constructors for the {\normalfont \ttfamily multiplier} halo mass function class.
     !!}
     module procedure multiplierConstructorParameters
     module procedure multiplierConstructorInternal
  end interface haloMassFunctionMultiplier

contains

  function multiplierConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily multiplier} halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloMassFunctionMultiplier)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (haloMassFunctionClass     ), pointer       :: massFunction_
    class           (cosmologyFunctionsClass   ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass  ), pointer       :: cosmologyParameters_
    double precision                                            :: multiplier          , exponentRedshift

    !![
    <inputParameter>
      <name>multiplier</name>
      <source>parameters</source>
      <description>The factor by which to multiply the halo mass function.</description>
    </inputParameter>
    <inputParameter>
      <name>exponentRedshift</name>
      <source>parameters</source>
      <description>The exponent of $(1+z)$ in the multiplier of the halo mass function.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="haloMassFunction"    name="massFunction_"        source="parameters"/>
    !!]
    self=haloMassFunctionMultiplier(multiplier,exponentRedshift,massFunction_,cosmologyParameters_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massFunction_"       />
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function multiplierConstructorParameters

  function multiplierConstructorInternal(multiplier,exponentRedshift,massFunction_,cosmologyParameters_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily multiplier} halo mass function class.
    !!}
    implicit none
    type            (haloMassFunctionMultiplier)                        :: self
    class           (haloMassFunctionClass     ), target, intent(in   ) :: massFunction_
    class           (cosmologyParametersClass  ), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass   ), target, intent(in   ) :: cosmologyFunctions_
    double precision                                    , intent(in   ) :: multiplier          , exponentRedshift
    !![
    <constructorAssign variables="multiplier, exponentRedshift, *cosmologyParameters_, *cosmologyFunctions_, *massFunction_"/>
    !!]

    return
  end function multiplierConstructorInternal

  subroutine multiplierDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily multiplier} halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionMultiplier), intent(inout) :: self

    !![
    <objectDestructor name="self%massFunction_"       />
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    return
  end subroutine multiplierDestructor

  double precision function multiplierDifferential(self,time,mass,node) result(massFunction)
    !!{
    Return the differential halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionMultiplier), intent(inout), target   :: self
    double precision                            , intent(in   )           :: time, mass
    type            (treeNode                  ), intent(inout), optional :: node

    massFunction=+(                                                                                   &
         &         +(                                                                                 &
         &           +self                    %multiplier                                             &
         &           -1.0d0                                                                           &
         &          )                                                                                 &
         &         /  self%cosmologyFunctions_%expansionFactor(time          )**self%exponentRedshift &
         &         +  1.0d0                                                                           &
         &        )                                                                                   &
         &       *    self%massFunction_      %differential   (time,mass,node)
    return
  end function multiplierDifferential

  double precision function multiplierIntegrated(self,time,massLow,massHigh,node,status) result(massFunction)
    !!{
    Return the integrated halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionMultiplier), intent(inout), target           :: self
    double precision                            , intent(in   )                   :: time    , massLow, &
         &                                                                           massHigh
    type            (treeNode                  ), intent(inout), target, optional :: node 
    integer                                     , intent(  out)        , optional :: status

    massFunction=+(                                                                                   &
         &         +(                                                                                 &
         &           +self                    %multiplier                                             &
         &           -1.0d0                                                                           &
         &          )                                                                                 &
         &         /  self%cosmologyFunctions_%expansionFactor(time          )**self%exponentRedshift &
         &         +  1.0d0                                                                           &
         &        )                                                                                   &
         &       *self%massFunction_%integrated(time,massLow,massHigh,node,status)
    return
  end function multiplierIntegrated
  
