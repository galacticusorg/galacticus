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
Implements a galactic filter for halo mass under a given definition.
!!}

  use :: Cosmology_Parameters   , only : cosmologyParametersClass
  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Virial_Density_Contrast, only : virialDensityContrastClass

  !![
  <galacticFilter name="galacticFilterHaloMassRange">
   <description>
    Halos with a halo mass {\normalfont \ttfamily massLow} $\le M &lt;$ {\normalfont \ttfamily [massHigh]} are passed.
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterHaloMassRange
     !!{
     A galactic filter class on halo mass.
     !!}
     private
     double precision                                      :: massLow                         , massHigh
     class           (virialDensityContrastClass), pointer :: virialDensityContrast_ => null(), virialDensityContrastDefinition_ => null()
     class           (cosmologyParametersClass  ), pointer :: cosmologyParameters_   => null()
     class           (cosmologyFunctionsClass   ), pointer :: cosmologyFunctions_    => null()
   contains
     final     ::           haloMassRangeDestructor
     procedure :: passes => haloMassRangePasses
  end type galacticFilterHaloMassRange

  interface galacticFilterHaloMassRange
     !!{
     Constructors for the ``haloMassRange'' galactic filter class.
     !!}
     module procedure haloMassRangeConstructorParameters
     module procedure haloMassRangeConstructorInternal
  end interface galacticFilterHaloMassRange

contains

  function haloMassRangeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``haloMassRange'' galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterHaloMassRange)                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass    ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass   ), pointer       :: cosmologyParameters_
    class           (virialDensityContrastClass ), pointer       :: virialDensityContrast_, virialDensityContrastDefinition_
    double precision                                             :: massLow               , massHigh

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>massLow</name>
      <source>parameters</source>
      <description>The lowest halo mass to pass.</description>
    </inputParameter>
    <inputParameter>
      <name>massHigh</name>
      <source>parameters</source>
      <description>The highest halo mass to pass.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"              source="parameters"                                                />
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"             source="parameters"                                                />
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_"           source="parameters"                                                />
    <objectBuilder class="virialDensityContrast" name="virialDensityContrastDefinition_" source="parameters" parameterName="virialDensityContrastDefinition"/>
    !!]
    self=galacticFilterHaloMassRange(massLow,massHigh,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,virialDensityContrastDefinition_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"             />
    <objectDestructor name="cosmologyParameters_"            />
    <objectDestructor name="virialDensityContrast_"          />
    <objectDestructor name="virialDensityContrastDefinition_"/>
    !!]
    return
  end function haloMassRangeConstructorParameters

  function haloMassRangeConstructorInternal(massLow,massHigh,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,virialDensityContrastDefinition_) result(self)
    !!{
    Internal constructor for the ``haloMassRange'' galactic filter class.
    !!}
    implicit none
    type            (galacticFilterHaloMassRange)                        :: self
    double precision                             , intent(in   )         :: massLow               , massHigh
    class           (cosmologyParametersClass   ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass    ), intent(in   ), target :: cosmologyFunctions_
    class           (virialDensityContrastClass ), intent(in   ), target :: virialDensityContrast_, virialDensityContrastDefinition_
    !![
    <constructorAssign variables="massLow, massHigh, *cosmologyFunctions_, *cosmologyParameters_, *virialDensityContrast_, *virialDensityContrastDefinition_"/>
    !!]

    return
  end function haloMassRangeConstructorInternal

  subroutine haloMassRangeDestructor(self)
    !!{
    Destructor for the ``haloMassRange'' galactic filter class.
    !!}
    implicit none
    type(galacticFilterHaloMassRange), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%virialDensityContrast_"          />
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    !!]
    return
  end subroutine haloMassRangeDestructor

  logical function haloMassRangePasses(self,node)
    !!{
    Implement a halo mass high-pass galactic filter.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    class           (galacticFilterHaloMassRange), intent(inout)         :: self
    type            (treeNode                   ), intent(inout), target :: node
    class           (nodeComponentBasic         ), pointer               :: basic
    double precision                                                     :: massHalo

    basic               =>                                                                                                                  node %basic()
    massHalo            =  Dark_Matter_Profile_Mass_Definition(                                                                                            &
         &                                                                            node                                                               , &
         &                                                                            self%virialDensityContrastDefinition_%densityContrast(               &
         &                                                                                                                                  basic%mass (), &
         &                                                                                                                                  basic%time ()  &
         &                                                                                                                                 )             , &
         &                                                     cosmologyParameters_  =self%cosmologyParameters_                                          , &
         &                                                     cosmologyFunctions_   =self%cosmologyFunctions_                                           , &
         &                                                     virialDensityContrast_=self%virialDensityContrast_                                          &
         &                                                    ) 
    haloMassRangePasses =  massHalo >= self%massLow  &
         &                .and.                      &
         &                 massHalo <  self%massHigh
    return
  end function haloMassRangePasses
