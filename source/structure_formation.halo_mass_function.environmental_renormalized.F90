!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which implements a dark matter halo mass function class which renormalizes to account for environmental dependence.
!!}


  !![
  <haloMassFunction name="haloMassFunctionEnvironmentalRenormalized">
   <description>
    The halo mass function is renormalized to account for environmental dependence.
   </description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionClass) :: haloMassFunctionEnvironmentalRenormalized
     !!{
     A halo mass function class which renormalizes to account for environmental dependence.
     !!}
     private
     class           (haloMassFunctionClass), pointer :: haloMassFunctionConditioned_     => null(), haloMassFunctionUnconditioned_ => null(), &
          &                                              haloMassFunction_                => null()
     class           (haloEnvironmentClass ), pointer :: haloEnvironment_                 => null()
     double precision                                 :: normalization                             , timeNormalized                          , &
          &                                              fractionMassEnvironmentNormalize
   contains
     final     ::                 environmentalRenormalizedDestructor
     procedure :: differential => environmentalRenormalizedDifferential
  end type haloMassFunctionEnvironmentalRenormalized

  interface haloMassFunctionEnvironmentalRenormalized
     !!{
     Constructors for the {\normalfont \ttfamily environmentalRenormalized} halo mass function class.
     !!}
     module procedure environmentalRenormalizedConstructorParameters
     module procedure environmentalRenormalizedConstructorInternal
  end interface haloMassFunctionEnvironmentalRenormalized

contains

  function environmentalRenormalizedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily environmentalRenormalized} halo mass function class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloMassFunctionEnvironmentalRenormalized)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (haloMassFunctionClass                    ), pointer       :: haloMassFunctionConditioned_    , haloMassFunctionUnconditioned_, &
         &                                                                        haloMassFunction_
    class           (haloEnvironmentClass                     ), pointer       :: haloEnvironment_
    class           (cosmologyParametersClass                 ), pointer       :: cosmologyParameters_
    double precision                                                           :: fractionMassEnvironmentNormalize

    !![
    <inputParameter>
      <name>fractionMassEnvironmentNormalize</name>
      <source>parameters</source>
      <description>The fraction of the environment mass at which to compute the renormalization.</description>
    </inputParameter>
    <objectBuilder class="haloMassFunction"    name="haloMassFunction_"              source="parameters"                                              />
    <objectBuilder class="haloMassFunction"    name="haloMassFunctionConditioned_"   source="parameters" parameterName="haloMassFunctionConditioned"  />
    <objectBuilder class="haloMassFunction"    name="haloMassFunctionUnconditioned_" source="parameters" parameterName="haloMassFunctionUnconditioned"/>
    <objectBuilder class="haloEnvironment"     name="haloEnvironment_"               source="parameters"                                              />
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_"           source="parameters"                                              />
    !!]
    self=haloMassFunctionEnvironmentalRenormalized(fractionMassEnvironmentNormalize,haloMassFunction_,haloMassFunctionConditioned_,haloMassFunctionUnconditioned_,haloEnvironment_,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="haloMassFunction_"             />
    <objectDestructor name="haloMassFunctionConditioned_"  />
    <objectDestructor name="haloMassFunctionUnconditioned_"/>
    <objectDestructor name="haloEnvironment_"              />
    <objectDestructor name="cosmologyParameters_"          />
    !!]
    return
  end function environmentalRenormalizedConstructorParameters

  function environmentalRenormalizedConstructorInternal(fractionMassEnvironmentNormalize,haloMassFunction_,haloMassFunctionConditioned_,haloMassFunctionUnconditioned_,haloEnvironment_,cosmologyParameters_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily environmentalRenormalized} halo mass function class.
    !!}
    implicit none
    type            (haloMassFunctionEnvironmentalRenormalized)                        :: self
    class           (haloMassFunctionClass                    ), target, intent(in   ) :: haloMassFunctionConditioned_    , haloMassFunctionUnconditioned_, &
         &                                                                                haloMassFunction_
    class           (haloEnvironmentClass                     ), target, intent(in   ) :: haloEnvironment_
    class           (cosmologyParametersClass                 ), target, intent(in   ) :: cosmologyParameters_
    double precision                                                   , intent(in   ) :: fractionMassEnvironmentNormalize
    !![
    <constructorAssign variables="fractionMassEnvironmentNormalize, *haloMassFunction_, *haloMassFunctionConditioned_, *haloMassFunctionUnconditioned_, *haloEnvironment_, *cosmologyParameters_"/>
    !!]

    self%timeNormalized=-1.0d0
    return
  end function environmentalRenormalizedConstructorInternal

  subroutine environmentalRenormalizedDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily shethTormen} halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionEnvironmentalRenormalized), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"          />
    <objectDestructor name="self%haloMassFunction_"             />
    <objectDestructor name="self%haloMassFunctionConditioned_"  />
    <objectDestructor name="self%haloMassFunctionUnconditioned_"/>
    <objectDestructor name="self%haloEnvironment_"              />
    !!]
    return
  end subroutine environmentalRenormalizedDestructor

  double precision function environmentalRenormalizedDifferential(self,time,mass,node)
    !!{
    Return the differential halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionEnvironmentalRenormalized), intent(inout), target   :: self
    double precision                                           , intent(in   )           :: time         , mass
    type            (treeNode                                 ), intent(inout), optional :: node
    double precision                                                                     :: massNormalize

    if (time /= self%timeNormalized) then
       massNormalize      =+self%haloEnvironment_              %environmentMass                 (                       ) &
            &              *self%                               fractionMassEnvironmentNormalize
       self%normalization =+self%haloMassFunctionUnconditioned_%differential                    (time,massNormalize,node) &
            &              /self%haloMassFunctionConditioned_  %differential                    (time,massNormalize,node)
       self%timeNormalized=+                                                                     time
    end if
    environmentalRenormalizedDifferential=+self%haloMassFunction_%differential(time,mass,node) &
         &                                *self%normalization
    return
  end function environmentalRenormalizedDifferential
