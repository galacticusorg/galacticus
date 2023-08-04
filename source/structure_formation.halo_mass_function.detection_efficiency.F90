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
Implements a dark matter halo mass function class which modifies another mass function by account for halo-finder detection efficiency.
!!}

  !![
  <haloMassFunction name="haloMassFunctionDetectionEfficiency">
   <description>
    The halo mass function is computed by modifying another mass function by the halo-finder detection efficiency.
   </description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionClass) :: haloMassFunctionDetectionEfficiency
     !!{
     A halo mass function class that modifies another mass function by the halo-finder detection efficiency.
     !!}
     private
     double precision                                 :: massMinimum            , efficiencyAtMassMinimum, &
          &                                              exponentMass
     class           (haloMassFunctionClass), pointer :: massFunction_ => null()
   contains
     final     ::                 detectionEfficiencyDestructor
     procedure :: differential => detectionEfficiencyDifferential
  end type haloMassFunctionDetectionEfficiency

  interface haloMassFunctionDetectionEfficiency
     !!{
     Constructors for the {\normalfont \ttfamily detectionEfficiency} halo mass function class.
     !!}
     module procedure detectionEfficiencyConstructorParameters
     module procedure detectionEfficiencyConstructorInternal
  end interface haloMassFunctionDetectionEfficiency

contains

  function detectionEfficiencyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily detectionEfficiency} halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloMassFunctionDetectionEfficiency)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (haloMassFunctionClass              ), pointer       :: massFunction_
    class           (cosmologyParametersClass           ), pointer       :: cosmologyParameters_
    double precision                                                     :: massMinimum         , efficiencyAtMassMinimum, &
         &                                                                  exponentMass

    !![
    <inputParameter>
      <name>massMinimum</name>
      <source>parameters</source>
      <description>The minimum mass at which halos are detected.</description>
    </inputParameter>
    <inputParameter>
      <name>efficiencyAtMassMinimum</name>
      <source>parameters</source>
      <description>The efficiency of detection at the minimum mass.</description>
    </inputParameter>
    <inputParameter>
      <name>exponentMass</name>
      <source>parameters</source>
      <description>The exponent $\alpha$ in the detection efficiency, $f(M) = 1 - (1-\epsilon) (M/M_\mathrm{min})^\alpha$.</description>
    </inputParameter>
    <objectBuilder class="haloMassFunction"    name="massFunction_"        source="parameters"/>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=haloMassFunctionDetectionEfficiency(massMinimum,efficiencyAtMassMinimum,exponentMass,massFunction_,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massFunction_"       />
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function detectionEfficiencyConstructorParameters

  function detectionEfficiencyConstructorInternal(massMinimum,efficiencyAtMassMinimum,exponentMass,massFunction_,cosmologyParameters_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily detectionEfficiency} halo mass function class.
    !!}
    implicit none
    type            (haloMassFunctionDetectionEfficiency)                        :: self
    class           (haloMassFunctionClass              ), target, intent(in   ) :: massFunction_
    class           (cosmologyParametersClass           ), target, intent(in   ) :: cosmologyParameters_
    double precision                                             , intent(in   ) :: massMinimum         , efficiencyAtMassMinimum, &
         &                                                                          exponentMass
    !![
    <constructorAssign variables="massMinimum, efficiencyAtMassMinimum, exponentMass, *cosmologyParameters_, *massFunction_"/>
    !!]

    return
  end function detectionEfficiencyConstructorInternal

  subroutine detectionEfficiencyDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily detectionEfficiency} halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionDetectionEfficiency), intent(inout) :: self

    !![
    <objectDestructor name="self%massFunction_"       />
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine detectionEfficiencyDestructor

  double precision function detectionEfficiencyDifferential(self,time,mass,node) result(massFunction)
    !!{
    Return the differential halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionDetectionEfficiency), intent(inout), target   :: self
    double precision                                     , intent(in   )           :: time, mass
    type            (treeNode                           ), intent(inout), optional :: node

    if (mass >= self%massMinimum) then
       massFunction=+(                                               &
            &         +  1.0d0                                       &
            &         -(                                             &
            &           +1.0d0                                       &
            &           -self%efficiencyAtMassMinimum                &
            &          )                                             &
            &         *(                                             &
            &           +     mass                                   &
            &           /self%massMinimum                            &
            &          )**self%exponentMass                          &
            &        )                                               &
            &       *self%massFunction_%differential(time,mass,node)
    else
       massFunction=+0.0d0
    end if
    return
  end function detectionEfficiencyDifferential
