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
  Implementation of a satellite bound mass initializor class that sets the initial bound mass using a specified density contrast
  definition.
  !!}

  use :: Cosmology_Functions   , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters  , only : cosmologyParametersClass
  use :: Virial_Density_Contrast, only : virialDensityContrastClass

  !![
  <satelliteMassBoundInitializor name="satelliteMassBoundInitializorDensityContrast">
   <description>
    A satellite bound mass initializor class that sets the initial bound mass of the satellite halo to the mass enclosed within
    the radius defined by a given density contrast. The density contrast is evaluated using the
    \refClass{virialDensityContrastClass} object specified by \mono{[virialDensityContrastDefinition]}, and the mass is computed
    using the \refFunction{Dark\_Matter\_Profile\_Mass\_Definition} function, which accounts for the halo density profile.
   </description>
  </satelliteMassBoundInitializor>
  !!]
  type, extends(satelliteMassBoundInitializorClass) :: satelliteMassBoundInitializorDensityContrast
     !!{
     Implementation of a satellite bound mass initializor class that sets the initial bound mass using a specified density
     contrast definition.
     !!}
     private
     class(cosmologyFunctionsClass  ), pointer :: cosmologyFunctions_             => null()
     class(cosmologyParametersClass ), pointer :: cosmologyParameters_            => null()
     class(virialDensityContrastClass), pointer :: virialDensityContrast_          => null()
     class(virialDensityContrastClass), pointer :: virialDensityContrastDefinition_ => null()
   contains
     final     ::             densityContrastDestructor
     procedure :: massBound => densityContrastMassBound
  end type satelliteMassBoundInitializorDensityContrast

  interface satelliteMassBoundInitializorDensityContrast
     !!{
     Constructors for the \refClass{satelliteMassBoundInitializorDensityContrast} satellite bound mass initializor class.
     !!}
     module procedure densityContrastConstructorParameters
     module procedure densityContrastConstructorInternal
  end interface satelliteMassBoundInitializorDensityContrast

contains

  function densityContrastConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{satelliteMassBoundInitializorDensityContrast} satellite bound mass initializor class which
    builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (satelliteMassBoundInitializorDensityContrast)                :: self
    type (inputParameters                             ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                     ), pointer       :: cosmologyFunctions_
    class(cosmologyParametersClass                    ), pointer       :: cosmologyParameters_
    class(virialDensityContrastClass                  ), pointer       :: virialDensityContrast_
    class(virialDensityContrastClass                  ), pointer       :: virialDensityContrastDefinition_

    !![
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"             source="parameters"                                          />
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"            source="parameters"                                          />
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_"          source="parameters"                                          />
    <objectBuilder class="virialDensityContrast" name="virialDensityContrastDefinition_" source="parameters" parameterName="virialDensityContrastDefinition"/>
    !!]
    self=satelliteMassBoundInitializorDensityContrast(cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,virialDensityContrastDefinition_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"            />
    <objectDestructor name="cosmologyParameters_"           />
    <objectDestructor name="virialDensityContrast_"         />
    <objectDestructor name="virialDensityContrastDefinition_"/>
    !!]
    return
  end function densityContrastConstructorParameters

  function densityContrastConstructorInternal(cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,virialDensityContrastDefinition_) result(self)
    !!{
    Internal constructor for the \refClass{satelliteMassBoundInitializorDensityContrast} satellite bound mass initializor class.
    !!}
    implicit none
    type (satelliteMassBoundInitializorDensityContrast)                        :: self
    class(cosmologyFunctionsClass                     ), intent(in   ), target :: cosmologyFunctions_
    class(cosmologyParametersClass                    ), intent(in   ), target :: cosmologyParameters_
    class(virialDensityContrastClass                  ), intent(in   ), target :: virialDensityContrast_
    class(virialDensityContrastClass                  ), intent(in   ), target :: virialDensityContrastDefinition_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *cosmologyParameters_, *virialDensityContrast_, *virialDensityContrastDefinition_"/>
    !!]

    return
  end function densityContrastConstructorInternal

  subroutine densityContrastDestructor(self)
    !!{
    Destructor for the \refClass{satelliteMassBoundInitializorDensityContrast} satellite bound mass initializor class.
    !!}
    implicit none
    type(satelliteMassBoundInitializorDensityContrast), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"            />
    <objectDestructor name="self%cosmologyParameters_"           />
    <objectDestructor name="self%virialDensityContrast_"         />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    !!]
    return
  end subroutine densityContrastDestructor

  double precision function densityContrastMassBound(self,node)
    !!{
    Returns the initial bound mass of a satellite halo based on a specified density contrast definition.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    class(satelliteMassBoundInitializorDensityContrast), intent(inout) :: self
    type (treeNode                                    ), intent(inout) :: node
    class(nodeComponentBasic                          ), pointer       :: basic

    basic                      => node%basic()
    densityContrastMassBound   =  Dark_Matter_Profile_Mass_Definition(                                                                                                              &
         &                                                                         node                                                                                           , &
         &                                                             densityContrast       =self%virialDensityContrastDefinition_%densityContrast(basic%mass(),basic%time())    , &
         &                                                             cosmologyParameters_  =self%cosmologyParameters_                                                           , &
         &                                                             cosmologyFunctions_   =self%cosmologyFunctions_                                                            , &
         &                                                             virialDensityContrast_=self%virialDensityContrast_                                                           &
         &                                                            )
    return
  end function densityContrastMassBound

