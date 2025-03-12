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
Implements a virial radius output analysis property extractor class.
!!}

  use :: Cosmology_Parameters    , only : cosmologyParametersClass
  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Virial_Density_Contrast , only : virialDensityContrastClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRadiusVirial">
   <description>A virial radius output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRadiusVirial
     !!{
     A virial radius property extractor output analysis class. The property extracted is the ''\gls{dmou}'' virial radius,
     defined as the radius enclosing a density contrast as defined by the supplied {\normalfont \ttfamily
     virialDensityContrast} class object. Note that the density contrast is defined here at the time at which the halo
     presently exists, \emph{not} at the time at which is was last isolated (as is used for standard definition of virial
     radius).
     !!}
     private
     class  (virialDensityContrastClass), pointer :: virialDensityContrast_ => null(), virialDensityContrastDefinition_ => null()
     class  (cosmologyParametersClass  ), pointer :: cosmologyParameters_   => null()
     class  (cosmologyFunctionsClass   ), pointer :: cosmologyFunctions_    => null()
     class  (darkMatterProfileDMOClass ), pointer :: darkMatterProfileDMO_  => null()
     logical                                      :: useLastIsolatedTime
   contains
     final     ::                radiusVirialDestructor
     procedure :: extract     => radiusVirialExtract
     procedure :: name        => radiusVirialName
     procedure :: description => radiusVirialDescription
     procedure :: unitsInSI   => radiusVirialUnitsInSI
  end type nodePropertyExtractorRadiusVirial

  interface nodePropertyExtractorRadiusVirial
     !!{
     Constructors for the ``radiusVirial'' output analysis class.
     !!}
     module procedure radiusVirialConstructorParameters
     module procedure radiusVirialConstructorInternal
  end interface nodePropertyExtractorRadiusVirial

contains

  function radiusVirialConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``radiusVirial'' output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorRadiusVirial)                :: self
    type   (inputParameters                  ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    class  (cosmologyParametersClass         ), pointer       :: cosmologyParameters_
    class  (darkMatterProfileDMOClass        ), pointer       :: darkMatterProfileDMO_
    class  (virialDensityContrastClass       ), pointer       :: virialDensityContrast_, virialDensityContrastDefinition_
    logical                                                   :: useLastIsolatedTime

    !![
    <inputParameter>
     <name>useLastIsolatedTime</name>
     <source>parameters</source>
     <defaultValue>.false.</defaultValue>
     <description>If true, evaluate the halo virial radius using a the virial density definition at the last isolated time of the halo.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"              source="parameters"                                                />
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"             source="parameters"                                                />
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"            source="parameters"                                                />
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_"           source="parameters"                                                />
    <objectBuilder class="virialDensityContrast" name="virialDensityContrastDefinition_" source="parameters" parameterName="virialDensityContrastDefinition"/>
    !!]
    self=nodePropertyExtractorRadiusVirial(useLastIsolatedTime,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"             />
    <objectDestructor name="cosmologyParameters_"            />
    <objectDestructor name="darkMatterProfileDMO_"           />
    <objectDestructor name="virialDensityContrast_"          />
    <objectDestructor name="virialDensityContrastDefinition_"/>
    !!]
    return
  end function radiusVirialConstructorParameters

  function radiusVirialConstructorInternal(useLastIsolatedTime,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_) result(self)
    !!{
    Internal constructor for the ``radiusVirial'' output analysis property extractor class.
    !!}
    implicit none
    type   (nodePropertyExtractorRadiusVirial)                        :: self
    class  (cosmologyParametersClass         ), intent(in   ), target :: cosmologyParameters_
    class  (cosmologyFunctionsClass          ), intent(in   ), target :: cosmologyFunctions_
    class  (darkMatterProfileDMOClass        ), intent(in   ), target :: darkMatterProfileDMO_
    class  (virialDensityContrastClass       ), intent(in   ), target :: virialDensityContrast_, virialDensityContrastDefinition_
    logical                                   , intent(in   )         :: useLastIsolatedTime
    !![
    <constructorAssign variables="useLastIsolatedTime, *cosmologyFunctions_, *cosmologyParameters_, *darkMatterProfileDMO_, *virialDensityContrast_, *virialDensityContrastDefinition_"/>
    !!]

    return
  end function radiusVirialConstructorInternal

  subroutine radiusVirialDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily mass} output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorRadiusVirial), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%darkMatterProfileDMO_"           />
    <objectDestructor name="self%virialDensityContrast_"          />
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    !!]
    return
  end subroutine radiusVirialDestructor

  double precision function radiusVirialExtract(self,node,instance)
    !!{
    Implement a radiusVirial output analysis.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    class           (nodePropertyExtractorRadiusVirial), intent(inout), target   :: self
    type            (treeNode                         ), intent(inout), target   :: node
    type            (multiCounter                     ), intent(inout), optional :: instance
    class           (nodeComponentBasic               ), pointer                 :: basic
    double precision                                                             :: massVirial, time
    !$GLC attributes unused :: instance
    
    basic      => node %basic           ()
    if (self%useLastIsolatedTime) then
       time    =  basic%timeLastIsolated()
    else
       time    =  basic%time            ()
    end if
    massVirial =  Dark_Matter_Profile_Mass_Definition(                                                                                                &
         &                                                                  node                                                                    , &
         &                                                                  self%virialDensityContrastDefinition_%densityContrast(basic%mass(),time), &
         &                                            radius                =     radiusVirialExtract                                               , &
         &                                            cosmologyParameters_  =self%cosmologyParameters_                                              , &
         &                                            cosmologyFunctions_   =self%cosmologyFunctions_                                               , &
         &                                            virialDensityContrast_=self%virialDensityContrast_                                            , &
         &                                            darkMatterProfileDMO_ =self%darkMatterProfileDMO_                                             , &
         &                                            useLastIsolatedTime   =self%useLastIsolatedTime                                                 &
         &                                           )
    return
  end function radiusVirialExtract


  function radiusVirialName(self)
    !!{
    Return the name of the radiusVirial property.
    !!}
    implicit none
    type (varying_string                   )                :: radiusVirialName
    class(nodePropertyExtractorRadiusVirial), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusVirialName=var_str('radiusVirialCurrent')
    return
  end function radiusVirialName

  function radiusVirialDescription(self)
    !!{
    Return a description of the radiusVirial property.
    !!}
    implicit none
    type (varying_string                   )                :: radiusVirialDescription
    class(nodePropertyExtractorRadiusVirial), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusVirialDescription=var_str('The current virial radius of the dark-matter-only halo, defined as the radius enclosing a specified density contrast.')
    return
  end function radiusVirialDescription

  double precision function radiusVirialUnitsInSI(self)
    !!{
    Return the units of the radiusVirial property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    implicit none
    class(nodePropertyExtractorRadiusVirial), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusVirialUnitsInSI=megaParsec
    return
  end function radiusVirialUnitsInSI
