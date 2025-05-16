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
Implements a dark matter halo mass function class which modifies another mass function using a simple model for systematics.
!!}

  !![
  <haloMassFunction name="haloMassFunctionSimpleSystematic">
   <description>
  !!]

  !![
    A halo mass function class in which the mass function is computed by modifying another halo mass function (specified as a
    subparameter) using a simple model for systematic errors as follows:
    \begin{equation}
      {\mathrm{d} n\over \mathrm{d}M}(M) \rightarrow {\mathrm{d} n\over \mathrm{d}M}(M) \left( 1 + \alpha + \beta
      \log_{10}\left[ {M \over 10^{12}M_\odot} \right] \right)
    \end{equation}
    where $\alpha=${\normalfont \ttfamily [alpha]}, and $\beta=${\normalfont \ttfamily [beta]}.
   </description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionClass) :: haloMassFunctionSimpleSystematic
     !!{
     A halo mass function class which modifies another mass function using a simple model for systematics.
     !!}
     private
     double precision                                 :: alpha                          , beta
     class           (haloMassFunctionClass), pointer :: referenceMassFunction => null()
    contains
     final     ::                 simpleSystematicDestructor
     procedure :: differential => simpleSystematicDifferential
  end type haloMassFunctionSimpleSystematic

  interface haloMassFunctionSimpleSystematic
     !!{
     Constructors for the \refClass{haloMassFunctionSimpleSystematic} halo mass function class.
     !!}
     module procedure simpleSystematicConstructorParameters
     module procedure simpleSystematicConstructorInternal
  end interface haloMassFunctionSimpleSystematic

contains

  function simpleSystematicConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{haloMassFunctionSimpleSystematic} halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloMassFunctionSimpleSystematic)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (cosmologyParametersClass        ), pointer       :: cosmologyParameters_
    class           (haloMassFunctionClass           ), pointer       :: referenceMassFunction
    double precision                                                  :: alpha                , beta

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>alpha</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>Parameter $\alpha$ appearing in model for simple systematic shift in the halo mass function.</description>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>Parameter $\beta$ appearing in model for simple systematic shift in the halo mass function.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_"  source="parameters"/>
    <objectBuilder class="haloMassFunction"    name="referenceMassFunction" source="parameters"/>
    !!]
    self=haloMassFunctionSimpleSystematic(alpha,beta,cosmologyParameters_,referenceMassFunction)
    !![
    <inputParametersValidate source="parameters"/>
   <objectDestructor name="cosmologyParameters_" />
   <objectDestructor name="referenceMassFunction"/>
   !!]
   return
  end function simpleSystematicConstructorParameters

  function simpleSystematicConstructorInternal(alpha,beta,cosmologyParameters_,referenceMassFunction) result(self)
    !!{
    Internal constructor for the \refClass{haloMassFunctionSimpleSystematic} halo mass function class.
    !!}
    implicit none
    type            (haloMassFunctionSimpleSystematic)                        :: self
    class           (cosmologyParametersClass        ), target, intent(in   ) :: cosmologyParameters_
    class           (haloMassFunctionClass           ), target, intent(in   ) :: referenceMassFunction
    double precision                                          , intent(in   ) :: alpha                , beta
    !![
    <constructorAssign variables="alpha, beta, *cosmologyParameters_, *referenceMassFunction"/>
    !!]

    return
  end function simpleSystematicConstructorInternal

  subroutine simpleSystematicDestructor(self)
    !!{
    Destructor for the \refClass{haloMassFunctionSimpleSystematic} halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionSimpleSystematic), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"  />
    <objectDestructor name="self%referenceMassFunction" />
    !!]
    return
  end subroutine simpleSystematicDestructor

  double precision function simpleSystematicDifferential(self,time,mass,node)
    !!{
    Return the differential halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionSimpleSystematic), intent(inout), target   :: self
    double precision                                  , intent(in   )           :: time                , mass
    type            (treeNode                        ), intent(inout), optional :: node
    double precision                                  , parameter               :: massZeroPoint=1.0d12

    ! Compute the systematic shift in the halo mass function.
    simpleSystematicDifferential=+self%referenceMassFunction%differential(time,mass,node) &
         &                       *(                                                       &
         &                         +1.0d0                                                 &
         &                         +self%alpha                                            &
         &                         +self%beta                                             &
         &                         *log10(                                                &
         &                                 mass                                           &
         &                                /massZeroPoint                                  &
         &                               )                                                &
         &                        )
    return
  end function simpleSystematicDifferential
