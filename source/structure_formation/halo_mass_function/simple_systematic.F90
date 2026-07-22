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

!!{RST
Implements a dark matter halo mass function class which modifies another mass function using a simple model for systematics.
!!}

  !![
  <haloMassFunction name="haloMassFunctionSimpleSystematic" docformat="rst">
   <description>
   A halo mass function class in which the mass function is computed by modifying another halo mass function (specified as a subparameter) using a simple model for systematic errors as follows:

   .. math::

      {\mathrm{d} n\over \mathrm{d}M}(M) \rightarrow {\mathrm{d} n\over \mathrm{d}M}(M) \left( 1 + \alpha + \beta
      \log_{10}\left[ {M \over 10^{12}\mathrm{M}_\odot} \right] \right)

   where :math:`\alpha=`\ ``[alpha]``, and :math:`\beta=`\ ``[beta]``.
   </description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionClass) :: haloMassFunctionSimpleSystematic
     !!{RST
     A halo mass function class which modifies another mass function using a simple model for systematics.
     !!}
     private
     double precision                                 :: alpha                          , beta
     class           (haloMassFunctionClass), pointer :: referenceMassFunction => null()
    contains
     final     ::                 simpleSystematicDestructor
     procedure :: differential => simpleSystematicDifferential
     procedure :: isCriticalOverdensityDependent => simpleSystematicIsCriticalOverdensityDependent
  end type haloMassFunctionSimpleSystematic

  interface haloMassFunctionSimpleSystematic
     !!{RST
     Constructors for the :galacticus-class:`haloMassFunctionSimpleSystematic` halo mass function class.
     !!}
     module procedure simpleSystematicConstructorParameters
     module procedure simpleSystematicConstructorInternal
  end interface haloMassFunctionSimpleSystematic

contains

  function simpleSystematicConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`haloMassFunctionSimpleSystematic` halo mass function class which takes a parameter set as input.
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
    <inputParameter docformat="rst">
      <name>alpha</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>
      Parameter :math:`\alpha` appearing in model for simple systematic shift in the halo mass function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>beta</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>
      Parameter :math:`\beta` appearing in model for simple systematic shift in the halo mass function.
      </description>
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
    !!{RST
    Internal constructor for the :galacticus-class:`haloMassFunctionSimpleSystematic` halo mass function class.
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
    !!{RST
    Destructor for the :galacticus-class:`haloMassFunctionSimpleSystematic` halo mass function class.
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
    !!{RST
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

  logical function simpleSystematicIsCriticalOverdensityDependent(self)
    !!{RST
    Return whether the differential halo mass function depends on the critical overdensity for
    collapse, by forwarding the query to the wrapped halo mass function.
    !!}
    implicit none
    class(haloMassFunctionSimpleSystematic), intent(inout) :: self

    simpleSystematicIsCriticalOverdensityDependent=self%referenceMassFunction%isCriticalOverdensityDependent()
    return
  end function simpleSystematicIsCriticalOverdensityDependent
