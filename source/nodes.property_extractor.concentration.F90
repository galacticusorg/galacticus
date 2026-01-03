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
Implements a concentration output analysis property extractor class.
!!}

  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters    , only : cosmologyParametersClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Virial_Density_Contrast , only : virialDensityContrastClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorConcentration">
   <description>A concentration output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorConcentration
     !!{
     A concentration property extractor output analysis class. The property extracted is the ''\gls{dmou}'' concentration of
     the halo with a radius defined by a density contrast as given by the supplied {\normalfont \ttfamily
     virialDensityContrast} class object. Note that the density contrast is defined here at the time at which the halo
     presently exists, \emph{not} at the time at which is was last isolated (as is used for standard definition of
     concentration).
     !!}
     private
     class  (cosmologyParametersClass  ), pointer :: cosmologyParameters_   => null()
     class  (cosmologyFunctionsClass   ), pointer :: cosmologyFunctions_    => null()
     class  (virialDensityContrastClass), pointer :: virialDensityContrast_ => null(), virialDensityContrastDefinition_ => null()
     class  (darkMatterProfileDMOClass ), pointer :: darkMatterProfileDMO_  => null()
     logical                                      :: useLastIsolatedTime
   contains
     final     ::                concentrationDestructor
     procedure :: extract     => concentrationExtract
     procedure :: name        => concentrationName
     procedure :: description => concentrationDescription
     procedure :: unitsInSI   => concentrationUnitsInSI
  end type nodePropertyExtractorConcentration

  interface nodePropertyExtractorConcentration
     !!{
     Constructors for the \refClass{nodePropertyExtractorConcentration} output analysis class.
     !!}
     module procedure concentrationConstructorParameters
     module procedure concentrationConstructorInternal
  end interface nodePropertyExtractorConcentration

contains

  function concentrationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorConcentration} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorConcentration)                :: self
    type   (inputParameters                   ), intent(inout) :: parameters
    class  (cosmologyParametersClass          ), pointer       :: cosmologyParameters_
    class  (cosmologyFunctionsClass           ), pointer       :: cosmologyFunctions_
    class  (darkMatterProfileDMOClass         ), pointer       :: darkMatterProfileDMO_
    class  (virialDensityContrastClass        ), pointer       :: virialDensityContrast_, virialDensityContrastDefinition_
    logical                                                    :: useLastIsolatedTime

    !![
    <inputParameter>
     <name>useLastIsolatedTime</name>
     <source>parameters</source>
     <defaultValue>.false.</defaultValue>
     <description>If true, evaluate the concentration using a the virial density definition at the last isolated time of the halo.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"             source="parameters"                                                />
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"              source="parameters"                                                />
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"            source="parameters"                                                />
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_"           source="parameters"                                                />
    <objectBuilder class="virialDensityContrast" name="virialDensityContrastDefinition_" source="parameters" parameterName="virialDensityContrastDefinition"/>
    !!]
    self=nodePropertyExtractorConcentration(useLastIsolatedTime,cosmologyParameters_,cosmologyFunctions_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"            />
    <objectDestructor name="cosmologyFunctions_"             />
    <objectDestructor name="darkMatterProfileDMO_"           />
    <objectDestructor name="virialDensityContrast_"          />
    <objectDestructor name="virialDensityContrastDefinition_"/>
    !!]
    return
  end function concentrationConstructorParameters

  function concentrationConstructorInternal(useLastIsolatedTime,cosmologyParameters_,cosmologyFunctions_,darkMatterProfileDMO_,virialDensityContrast_,virialDensityContrastDefinition_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorConcentration} output analysis property extractor class.
    !!}
    implicit none
    type   (nodePropertyExtractorConcentration)                        :: self
    class  (cosmologyParametersClass          ), intent(in   ), target :: cosmologyParameters_
    class  (cosmologyFunctionsClass           ), intent(in   ), target :: cosmologyFunctions_
    class  (darkMatterProfileDMOClass         ), intent(in   ), target :: darkMatterProfileDMO_
    class  (virialDensityContrastClass        ), intent(in   ), target :: virialDensityContrast_, virialDensityContrastDefinition_
    logical                                    , intent(in   )         :: useLastIsolatedTime
    !![
    <constructorAssign variables="useLastIsolatedTime, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterProfileDMO_, *virialDensityContrast_, *virialDensityContrastDefinition_"/>
    !!]

    return
  end function concentrationConstructorInternal

  subroutine concentrationDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorConcentration} output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorConcentration), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%darkMatterProfileDMO_"           />
    <objectDestructor name="self%virialDensityContrast_"          />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    !!]
    return
  end subroutine concentrationDestructor

  double precision function concentrationExtract(self,node,instance)
    !!{
    Implement a concentration output analysis.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (nodePropertyExtractorConcentration), intent(inout), target   :: self
    type            (treeNode                          ), intent(inout), target   :: node
    type            (multiCounter                      ), intent(inout), optional :: instance
    class           (nodeComponentBasic                ), pointer                 :: basic
    class           (nodeComponentDarkMatterProfile    ), pointer                 :: darkMatterProfile
    double precision                                                              :: massHalo         , radiusHalo, &
         &                                                                           time
    !$GLC attributes unused :: instance

    basic                =>  node%basic            ()
    if (self%useLastIsolatedTime) then
       time              =   basic%timeLastIsolated()
    else
       time              =   basic%time            ()
    end if
    darkMatterProfile    =>  node%darkMatterProfile()
    massHalo             =   Dark_Matter_Profile_Mass_Definition(                                                                                                 &
         &                                                                              node                                                                    , &
         &                                                                              self%virialDensityContrastDefinition_%densityContrast(basic%mass(),time), &
         &                                                       radius                =     radiusHalo                                                         , &
         &                                                       cosmologyParameters_  =self%cosmologyParameters_                                               , &
         &                                                       cosmologyFunctions_   =self%cosmologyFunctions_                                                , &
         &                                                       virialDensityContrast_=self%virialDensityContrast_                                             , &
         &                                                       darkMatterProfileDMO_ =self%darkMatterProfileDMO_                                              , &
         &                                                       useLastIsolatedTime   =self%useLastIsolatedTime                                                  &
         &                                                      )
    concentrationExtract =  +                  radiusHalo   &
         &                  /darkMatterProfile%scale     ()
    return
  end function concentrationExtract

  function concentrationName(self)
    !!{
    Return the name of the concentration property.
    !!}
    implicit none
    type (varying_string                    )                :: concentrationName
    class(nodePropertyExtractorConcentration), intent(inout) :: self
    !$GLC attributes unused :: self

    concentrationName=var_str('concentration')
    return
  end function concentrationName

  function concentrationDescription(self)
    !!{
    Return a description of the concentration property.
    !!}
    implicit none
    type (varying_string                    )                :: concentrationDescription
    class(nodePropertyExtractorConcentration), intent(inout) :: self
    !$GLC attributes unused :: self

    concentrationDescription=var_str('The concentration parameter of the dark matter halos.')
    return
  end function concentrationDescription

  double precision function concentrationUnitsInSI(self)
    !!{
    Return the units of the concentration property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorConcentration), intent(inout) :: self
    !$GLC attributes unused :: self

    concentrationUnitsInSI=0.0d0
    return
  end function concentrationUnitsInSI

