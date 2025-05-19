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
  Implements a stellar feedback mass outflow rate property extractor class.
  !!}

  use :: Stellar_Feedback_Outflows     , only : stellarFeedbackOutflowsClass
  use :: Star_Formation_Rates_Disks    , only : starFormationRateDisksClass
  use :: Star_Formation_Rates_Spheroids, only : starFormationRateSpheroidsClass
  use :: Stellar_Population_Properties , only : stellarPopulationPropertiesClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorStellarFeedbackOutflowRate">
   <description>
    A node property extractor which extracts the stellar feedback-driven mass outflow rate from a galaxy. The type of mass outflow rate is controlled
    by the {\normalfont \ttfamily [component]} parameter, which can be either ``{\normalfont \ttfamily disk}'', ``{\normalfont
    \ttfamily spheroid}'', or ``{\normalfont \ttfamily total}''. The corresponding mass outflow rate is extracted as
    {\normalfont \ttfamily \textless\ component\textgreater\ StellarFeedbackOutflowRate} in units of $M_\odot$/Gyr.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorStellarFeedbackOutflowRate
     !!{
     A stellar feedback-driven mass outflow rate property extractor class.
     !!}
     private
     class(starFormationRateDisksClass     ), pointer :: starFormationRateDisks_      => null()
     class(starFormationRateSpheroidsClass ), pointer :: starFormationRateSpheroids_  => null()
     class(stellarPopulationPropertiesClass), pointer :: stellarPopulationProperties_ => null()
     class(stellarFeedbackOutflowsClass    ), pointer :: stellarFeedbackOutflows_     => null()
     type (varying_string                  )          :: name_                                 , description_
     type (enumerationGalacticComponentType)          :: component
   contains
     final     ::                stellarFeedbackOutflowRateDestructor
     procedure :: extract     => stellarFeedbackOutflowRateExtract
     procedure :: name        => stellarFeedbackOutflowRateName
     procedure :: description => stellarFeedbackOutflowRateDescription
     procedure :: unitsInSI   => stellarFeedbackOutflowRateUnitsInSI
  end type nodePropertyExtractorStellarFeedbackOutflowRate

  interface nodePropertyExtractorStellarFeedbackOutflowRate
     !!{
     Constructors for the {\normalfont \ttfamily stellarFeedbackOutflowRate} output analysis class.
     !!}
     module procedure stellarFeedbackOutflowRateConstructorParameters
     module procedure stellarFeedbackOutflowRateConstructorInternal
  end interface nodePropertyExtractorStellarFeedbackOutflowRate

contains

  function stellarFeedbackOutflowRateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily stellarFeedbackOutflowRate} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorStellarFeedbackOutflowRate)                :: self
    type (inputParameters                                ), intent(inout) :: parameters
    class(starFormationRateDisksClass                    ), pointer       :: starFormationRateDisks_
    class(starFormationRateSpheroidsClass                ), pointer       :: starFormationRateSpheroids_
    class(stellarPopulationPropertiesClass               ), pointer       :: stellarPopulationProperties_
    class(stellarFeedbackOutflowsClass                   ), pointer       :: stellarFeedbackOutflows_
    type (varying_string                                 )                :: component

    !![
    <inputParameter>
      <name>component</name>
      <source>parameters</source>
      <description>The component from which to extract star formation rate.</description>
    </inputParameter>
    <objectBuilder class="starFormationRateDisks"      name="starFormationRateDisks_"      source="parameters"/>
    <objectBuilder class="starFormationRateSpheroids" name="starFormationRateSpheroids_"   source="parameters"/>
    <objectBuilder class="stellarPopulationProperties" name="stellarPopulationProperties_" source="parameters"/>
    <objectBuilder class="stellarFeedbackOutflows"     name="stellarFeedbackOutflows_"     source="parameters"/>
    !!]
    self=nodePropertyExtractorStellarFeedbackOutflowRate(enumerationGalacticComponentEncode(char(component),includesPrefix=.false.),starFormationRateDisks_,starFormationRateSpheroids_,stellarPopulationProperties_,stellarFeedbackOutflows_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateDisks_"     />
    <objectDestructor name="starFormationRateSpheroids_" />
    <objectDestructor name="stellarPopulationProperties_"/>
    <objectDestructor name="stellarFeedbackOutflows_"    />
    !!]
    return
  end function stellarFeedbackOutflowRateConstructorParameters

  function stellarFeedbackOutflowRateConstructorInternal(component,starFormationRateDisks_,starFormationRateSpheroids_,stellarPopulationProperties_,stellarFeedbackOutflows_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily stellarFeedbackOutflowRate} property extractor class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (nodePropertyExtractorStellarFeedbackOutflowRate)                        :: self
    type (enumerationGalacticComponentType               ), intent(in   )         :: component
    class(starFormationRateDisksClass                    ), intent(in   ), target :: starFormationRateDisks_
    class(starFormationRateSpheroidsClass                ), intent(in   ), target :: starFormationRateSpheroids_
    class(stellarPopulationPropertiesClass               ), intent(in   ), target :: stellarPopulationProperties_
    class(stellarFeedbackOutflowsClass                   ), intent(in   ), target :: stellarFeedbackOutflows_
    !![
    <constructorAssign variables="component, *starFormationRateDisks_, *starFormationRateSpheroids_, *stellarPopulationProperties_, *stellarFeedbackOutflows_"/>
    !!]

    select case (component%ID)
    case (galacticComponentTotal   %ID)
       self%name_       ="totalStellarFeedbackOutflowRate"
       self%description_="Total (disk + spheroid) stellar feedback-driven outflow rate [M☉ Gyr⁻¹]."
    case (galacticComponentDisk    %ID)
       self%name_       ="diskStellarFeedbackOutflowRate"
       self%description_="Disk stellar feedback-driven outflow rate [M☉ Gyr⁻¹]."
    case (galacticComponentSpheroid%ID)
       self%name_       ="spheroidStellarFeedbackOutflowRate"
       self%description_="Spheroid stellar feedback-driven outflow rate [M☉ Gyr⁻¹]."
    case default
       call Error_Report('Unknown component.'//{introspection:location})
    end select
    return
  end function stellarFeedbackOutflowRateConstructorInternal

  subroutine stellarFeedbackOutflowRateDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily stellarFeedbackOutflowRate} property extractor class.
    !!}
    implicit none
    type   (nodePropertyExtractorStellarFeedbackOutflowRate), intent(inout) :: self
  
    !![
    <objectDestructor name="self%starFormationRateDisks_"     />
    <objectDestructor name="self%starFormationRateSpheroids_"     />
    <objectDestructor name="self%stellarPopulationProperties_"/>
    <objectDestructor name="self%stellarFeedbackOutflows_"    />
    !!]
    return
  end subroutine stellarFeedbackOutflowRateDestructor

  double precision function stellarFeedbackOutflowRateExtract(self,node,instance)
    !!{
    Implement an emission line output analysis property extractor.
    !!}
    use :: Abundances_Structure          , only : abundances
    use :: Galacticus_Nodes              , only : nodeComponentDisk  , nodeComponentSpheroid
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : stellarLuminosities
    implicit none
    class           (nodePropertyExtractorStellarFeedbackOutflowRate), intent(inout), target   :: self
    type            (treeNode                                       ), intent(inout), target   :: node
    type            (multiCounter                                   ), intent(inout), optional :: instance
    class           (nodeComponentDisk                              ), pointer                 :: disk
    class           (nodeComponentSpheroid                          ), pointer                 :: spheroid
    double precision                                                                           :: rateStarFormation          , rateEnergyInput             , &
         &                                                                                        rateOutflowEjectiveDisk    , rateOutflowExpulsiveDisk    , &
         &                                                                                        rateOutflowEjectiveSpheroid, rateOutflowExpulsiveSpheroid, &
         &                                                                                        rateMassStellar            , rateMassFuel                , &
         &                                                                                        massGas
    type            (abundances                                      )                         :: abundancesGas              , rateAbundancesFuels         , &
         &                                                                                        rateAbundancesStellar
    type            (history                                         )                         :: ratePropertiesStellar
    type            (stellarLuminosities                             )                         :: rateLuminositiesStellar    
    !$GLC attributes unused :: instance

    rateOutflowEjectiveSpheroid =0.0d0
    rateOutflowExpulsiveSpheroid=0.0d0
    rateOutflowEjectiveDisk     =0.0d0
    rateOutflowExpulsiveDisk    =0.0d0
    if (self%component == galacticComponentDisk     .or. self%component == galacticComponentTotal) then
       disk     => node%disk    ()
       if     (      disk    %angularMomentum() >= 0.0d0 &
            &  .and. disk    %radius         () >= 0.0d0 &
            &  .and. disk    %massGas        () >= 0.0d0 &
            & ) then
          ! Get the star formation rate.
          rateStarFormation=self%starFormationRateDisks_    %rate(node)   
          ! Compute abundances of star forming gas.
          massGas      =disk    %massGas      ()
          abundancesGas=disk    %abundancesGas()
          call abundancesGas%massToMassFraction(massGas)
          ! Find rates of change of stellar mass, gas mass, abundances and luminosities.
          ratePropertiesStellar=disk    %stellarPropertiesHistory()
          call self%stellarPopulationProperties_%rates(                                                      &
               &                                                                    rateStarFormation      , &
               &                                                                    abundancesGas          , &
               &                                                                    disk                   , &
               &                                                                    node                   , &
               &                                                                    ratePropertiesStellar  , &
               &                                                                    rateMassStellar        , &
               &                                                                    rateMassFuel           , &
               &                                                                    rateEnergyInput        , &
               &                                                                    rateAbundancesFuels    , &
               &                                                                    rateAbundancesStellar  , &
               &                                                                    rateLuminositiesStellar, &
               &                                       computeRateLuminosityStellar=.false.                  &
               &                                      )
          
          call self%stellarFeedbackOutflows_%outflowRate(disk    ,rateStarFormation,rateEnergyInput,rateOutflowEjectiveDisk    ,rateOutflowExpulsiveDisk    )
       end if
    end if
    if (self%component == galacticComponentSpheroid .or. self%component == galacticComponentTotal) then
       spheroid => node%spheroid()
       if     (      spheroid%angularMomentum() >= 0.0d0 &
            &  .and. spheroid%radius         () >= 0.0d0 &
            &  .and. spheroid%massGas        () >= 0.0d0 &
            & ) then
          ! Get the star formation rate.
          rateStarFormation=self%starFormationRateSpheroids_%rate(node)   
          ! Compute abundances of star forming gas.
          massGas      =spheroid%massGas      ()
          abundancesGas=spheroid%abundancesGas()
          call abundancesGas%massToMassFraction(massGas)
          ! Find rates of change of stellar mass, gas mass, abundances and luminosities.
          ratePropertiesStellar=spheroid%stellarPropertiesHistory()
          call self%stellarPopulationProperties_%rates(                                                      &
               &                                                                    rateStarFormation      , &
               &                                                                    abundancesGas          , &
               &                                                                    spheroid               , &
               &                                                                    node                   , &
               &                                                                    ratePropertiesStellar  , &
               &                                                                    rateMassStellar        , &
               &                                                                    rateMassFuel           , &
               &                                                                    rateEnergyInput        , &
               &                                                                    rateAbundancesFuels    , &
               &                                                                    rateAbundancesStellar  , &
               &                                                                    rateLuminositiesStellar, &
               &                                       computeRateLuminosityStellar=.false.                  &
               &                                      )
          
          call self%stellarFeedbackOutflows_%outflowRate(spheroid,rateStarFormation,rateEnergyInput,rateOutflowEjectiveSpheroid,rateOutflowExpulsiveSpheroid)
       end if
    end if
    ! Sum the rates.
    stellarFeedbackOutflowRateExtract=+rateOutflowEjectiveDisk    +rateOutflowExpulsiveDisk     &
         &                            +rateOutflowEjectiveSpheroid+rateOutflowExpulsiveSpheroid
    return
  end function stellarFeedbackOutflowRateExtract


  function stellarFeedbackOutflowRateName(self)
    !!{
    Return the name of the stellarFeedbackOutflowRate property.
    !!}
    implicit none
    type (varying_string                                 )                :: stellarFeedbackOutflowRateName
    class(nodePropertyExtractorStellarFeedbackOutflowRate), intent(inout) :: self

    stellarFeedbackOutflowRateName=self%name_
    return
  end function stellarFeedbackOutflowRateName

  function stellarFeedbackOutflowRateDescription(self)
    !!{
    Return a description of the stellarFeedbackOutflowRate property.
    !!}
    implicit none
    type (varying_string                                 )                  :: stellarFeedbackOutflowRateDescription
    class(nodePropertyExtractorStellarFeedbackOutflowRate), intent(inout) :: self

    stellarFeedbackOutflowRateDescription=self%description_
    return
  end function stellarFeedbackOutflowRateDescription

  double precision function stellarFeedbackOutflowRateUnitsInSI(self)
    !!{
    Return the units of the stellarFeedbackOutflowRate property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, gigaYear
    implicit none
    class(nodePropertyExtractorStellarFeedbackOutflowRate), intent(inout) :: self
    !$GLC attributes unused :: self

    stellarFeedbackOutflowRateUnitsInSI=massSolar/gigaYear
    return
  end function stellarFeedbackOutflowRateUnitsInSI
