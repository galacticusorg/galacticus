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
  Implements a :cite:t:`tinker_towardhalo_2008` dark matter halo mass function class with user-specified parameters.
  !!}

  !![
  <haloMassFunction name="haloMassFunctionTinker2008Generic" docformat="rst">
   <description>
   The dark matter halo mass function is computed using the empirical fitting function of :cite:t:`tinker_towardhalo_2008`, calibrated against N-body simulations over a wide range of halo masses and redshifts. The normalization :math:`A` and shape parameters :math:`a`, :math:`b`, :math:`c` of the fit can each be specified directly via input parameters.
   </description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionTinker2008Form) :: haloMassFunctionTinker2008Generic
     !!{RST
     A halo mass function class using the fitting function of :cite:t:`tinker_towardhalo_2008` with user-specified parameters.
     !!}
     private
     double precision :: normalization_, a_, &
          &              b_            , c_
   contains
     final     ::                  tinker2008GenericDestructor
     procedure :: normalization => tinker2008GenericNormalization
     procedure :: a             => tinker2008GenericA
     procedure :: b             => tinker2008GenericB
     procedure :: c             => tinker2008GenericC
  end type haloMassFunctionTinker2008Generic

  interface haloMassFunctionTinker2008Generic
     !!{RST
     Constructors for the :galacticus-class:`haloMassFunctionTinker2008Generic` halo mass function class.
     !!}
     module procedure tinker2008GenericConstructorParameters
     module procedure tinker2008GenericConstructorInternal
  end interface haloMassFunctionTinker2008Generic

contains

  function tinker2008GenericConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`haloMassFunctionTinker2008Generic` halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloMassFunctionTinker2008Generic)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    class           (cosmologyParametersClass         ), pointer       :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass    ), pointer       :: cosmologicalMassVariance_
    class           (linearGrowthClass                ), pointer       :: linearGrowth_
    class           (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    double precision                                                   :: normalization            , a, &
         &                                                                b                        , c

    !![
    <inputParameter docformat="rst">
      <name>normalization</name>
      <defaultValue>0.150d0</defaultValue>
      <defaultSource>
      :cite:p:`trac_scorch_2015`
      </defaultSource>
      <description>
      The normalization parameter, :math:`A`, for the :cite:t:`tinker_towardhalo_2008` halo mass function.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>a</name>
      <defaultValue>1.36d0</defaultValue>
      <defaultSource>
      :cite:p:`trac_scorch_2015`
      </defaultSource>
      <description>
      The parameter :math:`a` for the :cite:t:`tinker_towardhalo_2008` halo mass function.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>b</name>
      <defaultValue>2.54d0</defaultValue>
      <defaultSource>
      :cite:p:`trac_scorch_2015`
      </defaultSource>
      <description>
      The parameter :math:`b` for the :cite:t:`tinker_towardhalo_2008` halo mass function.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>c</name>
      <defaultValue>1.14d0</defaultValue>
      <defaultSource>
      :cite:p:`trac_scorch_2015`
      </defaultSource>
      <description>
      The parameter :math:`c` for the :cite:t:`tinker_towardhalo_2008` halo mass function.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    !!]
    self=haloMassFunctionTinker2008Generic(                           &
         &                                 normalization            , &
         &                                 a                        , &
         &                                 b                        , &
         &                                 c                        , &
         &                                 cosmologyParameters_     , &
         &                                 cosmologicalMassVariance_, &
         &                                 linearGrowth_            , &
         &                                 cosmologyFunctions_        &
         &                                )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="linearGrowth_"            />
    <objectDestructor name="cosmologyFunctions_"      />
    !!]
    return
  end function tinker2008GenericConstructorParameters

  function tinker2008GenericConstructorInternal(normalization,a,b,c,cosmologyParameters_,cosmologicalMassVariance_,linearGrowth_,cosmologyFunctions_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`haloMassFunctionTinker2008Generic` halo mass function class.
    !!}
    implicit none
    type            (haloMassFunctionTinker2008Generic)                             :: self
    class           (cosmologyParametersClass         ), target     , intent(in   ) :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass    ), target     , intent(in   ) :: cosmologicalMassVariance_
    class           (linearGrowthClass                ), target     , intent(in   ) :: linearGrowth_
    class           (cosmologyFunctionsClass          ), target     , intent(in   ) :: cosmologyFunctions_
    double precision                                                , intent(in   ) :: normalization            , a, &
         &                                                                             b                        , c
    !![
    <constructorAssign variables="*cosmologyParameters_, *cosmologicalMassVariance_, *linearGrowth_, *cosmologyFunctions_"/>
    !!]

    self%time              =-1.0d0
    self%mass              =-1.0d0
    self%normalization_=normalization
    self%            a_=a
    self%            b_=b
    self%            c_=c
    return
  end function tinker2008GenericConstructorInternal

  subroutine tinker2008GenericDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`haloMassFunctionTinker2008Generic` halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionTinker2008Generic), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%linearGrowth_"            />
    <objectDestructor name="self%cosmologyParameters_"     />
    !!]
    return
  end subroutine tinker2008GenericDestructor

  double precision function tinker2008GenericNormalization(self,time,mass)
    !!{RST
    Return the normalization for the ``tinker2008Generic`` halo mass function class.
    !!}
    implicit none
    class           (haloMassFunctionTinker2008Generic), intent(inout) :: self
    double precision                                   , intent(in   ) :: time, mass
    !$GLC attributes unused :: time, mass

    tinker2008GenericNormalization=self%normalization_
    return
  end function tinker2008GenericNormalization

  double precision function tinker2008GenericA(self,time,mass)
    !!{RST
    Return the :math:`a` parameter for the ``tinker2008Generic`` halo mass function class.
    !!}
    implicit none
    class           (haloMassFunctionTinker2008Generic), intent(inout) :: self
    double precision                                   , intent(in   ) :: time, mass
    !$GLC attributes unused :: time, mass

    tinker2008GenericA=self%a_
    return
  end function tinker2008GenericA

  double precision function tinker2008GenericB(self,time,mass)
    !!{RST
    Return the :math:`b` parameter for the ``tinker2008Generic`` halo mass function class.
    !!}
    implicit none
    class           (haloMassFunctionTinker2008Generic), intent(inout) :: self
    double precision                                   , intent(in   ) :: time, mass
    !$GLC attributes unused :: time, mass

    tinker2008GenericB=self%b_
    return
  end function tinker2008GenericB

  double precision function tinker2008GenericC(self,time,mass)
    !!{RST
    Return the :math:`c` parameter for the ``tinker2008Generic`` halo mass function class.
    !!}
    implicit none
    class           (haloMassFunctionTinker2008Generic), intent(inout) :: self
    double precision                                   , intent(in   ) :: time, mass
    !$GLC attributes unused :: time, mass

    tinker2008GenericC=self%c_
    return
  end function tinker2008GenericC

