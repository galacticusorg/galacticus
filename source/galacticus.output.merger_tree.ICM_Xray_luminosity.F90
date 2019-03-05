!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which handles outputting of rotation curve data to the \glc\ output file.

module Galacticus_Output_Trees_ICM_Xray_Luminosity
  !% Handles outputting of ICM X-ray luminosity data to the \glc\ output file.
  use ISO_Varying_String
  implicit none
  private
  public :: Galacticus_Output_Tree_ICM_Xray_Luminosity      , Galacticus_Output_Tree_ICM_Xray_Luminosity_Property_Count, &
       &    Galacticus_Output_Tree_ICM_Xray_Luminosity_Names

  ! Flag indicating whether or not ICM X-ray luminosity information is to be output.
  logical :: outputICMXrayLuminosities

  ! Flag indicating whether or not this module has been initialized.
  logical :: outputICMXrayLuminositiesInitialized=.false.

contains

  subroutine Galacticus_Output_Tree_ICM_Xray_Luminosity_Initialize()
    !% Initializes the module by determining whether or not ICM X-ray luminosities should be output.
    use Input_Parameters
    implicit none

    if (.not.outputICMXrayLuminositiesInitialized) then
       !$omp critical(Galacticus_Output_Tree_ICM_Xray_Luminosity_Initialize)
       if (.not.outputICMXrayLuminositiesInitialized) then
          !# <inputParameter>
          !#   <name>outputICMXrayLuminosities</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>.false.</defaultValue>
          !#   <description>Specifies whether or not ICM X-ray luminosities should be included in the output file.</description>
          !#   <group>output</group>
          !#   <source>globalParameters</source>
          !#   <type>boolean</type>
          !# </inputParameter>
          ! Flag that module is now initialized.
          outputICMXrayLuminositiesInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_ICM_Xray_Luminosity_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_ICM_Xray_Luminosity_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_ICM_Xray_Luminosity_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_ICM_Xray_Luminosity</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_ICM_Xray_Luminosity_Names(node,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of projected density properties to be written to the \glc\ output file.
    use Galacticus_Nodes            , only : treeNode
    use Numerical_Constants_Units   , only : ergs    , electronVolt
    use Numerical_Constants_Prefixes, only : kilo
    implicit none
    type            (treeNode)              , intent(inout) :: node
    double precision                        , intent(in   ) :: time
    integer                                 , intent(inout) :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout) :: doublePropertyComments , doublePropertyNames   , &
         &                                                     integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout) :: doublePropertyUnitsSI  , integerPropertyUnitsSI
    !GCC$ attributes unused :: node, integerProperty, integerPropertyNames, integerPropertyComments, integerPropertyUnitsSI, time
    
    ! Initialize the module.
    call Galacticus_Output_Tree_ICM_Xray_Luminosity_Initialize()

    ! Return property names if we are outputting ICM X-ray luminosities.
    if (outputICMXrayLuminosities) then
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='icmXrayLuminosity'
       doublePropertyComments(doubleProperty)='X-ray luminosity of the ICM [ergs/s]'
       doublePropertyUnitsSI (doubleProperty)=ergs
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='icmXrayTemperature'
       doublePropertyComments(doubleProperty)='X-ray luminosity-weighted temperature of the ICM [keV]'
       doublePropertyUnitsSI (doubleProperty)=kilo*electronVolt
    end if
    return
  end subroutine Galacticus_Output_Tree_ICM_Xray_Luminosity_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_ICM_Xray_Luminosity_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_ICM_Xray_Luminosity</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_ICM_Xray_Luminosity_Property_Count(node,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of projected density properties to be written to the \glc\ output file.
    use Galacticus_Nodes, only : treeNode
    implicit none
    type            (treeNode), intent(inout) :: node
    double precision          , intent(in   ) :: time
    integer                   , intent(inout) :: doublePropertyCount, integerPropertyCount
    !GCC$ attributes unused :: node, time, integerPropertyCount
    
    ! Initialize the module.
    call Galacticus_Output_Tree_ICM_Xray_Luminosity_Initialize()

    ! Increment property count if we are outputting ICM X-ray luminosities.
    if (outputICMXrayLuminosities) doublePropertyCount=doublePropertyCount+2
    return
  end subroutine Galacticus_Output_Tree_ICM_Xray_Luminosity_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_ICM_Xray_Luminosity</unitName>
  !#  <sortName>Galacticus_Output_Tree_ICM_Xray_Luminosity</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_ICM_Xray_Luminosity(node,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time,instance)
    !% Store projected density properties in the \glc\ output file buffers.
    use FGSL                             , only : fgsl_function            , fgsl_integration_workspace
    use Numerical_Integration
    use Galacticus_Nodes                 , only : treeNode                 , nodeComponentHotHalo
    use Kind_Numbers
    use Multi_Counters
    use Dark_Matter_Halo_Scales          , only : darkMatterHaloScale      , darkMatterHaloScaleClass
    use Hot_Halo_Mass_Distributions      , only : hotHaloMassDistribution  , hotHaloMassDistributionClass
    use Hot_Halo_Temperature_Profiles    , only : hotHaloTemperatureProfile, hotHaloTemperatureProfileClass
    use Cooling_Functions                , only : coolingFunction          , coolingFunctionClass
    use Cosmology_Functions              , only : cosmologyFunctions       , cosmologyFunctionsClass
    use Chemical_Reaction_Rates_Utilities
    use Radiation_Fields
    use Numerical_Constants_Prefixes     , only : kilo
    use Numerical_Constants_Physical     , only : boltzmannsConstant
    use Numerical_Constants_Units        , only : electronVolt
    implicit none
    double precision                                         , intent(in   )                 :: time
    type            (treeNode                               ), intent(inout)                 :: node
    integer                                                  , intent(inout)                 :: doubleBufferCount         , doubleProperty, &
         &                                                                                      integerBufferCount        , integerProperty
    integer         (kind=kind_int8                         ), intent(inout), dimension(:,:) :: integerBuffer
    double precision                                         , intent(inout), dimension(:,:) :: doubleBuffer
    type            (multiCounter                           ), intent(inout)                 :: instance
    class           (darkMatterHaloScaleClass               ), pointer                       :: darkMatterHaloScale_
    class           (hotHaloMassDistributionClass           ), pointer                       :: hotHaloMassDistribution_
    class           (hotHaloTemperatureProfileClass         ), pointer                       :: hotHaloTemperatureProfile_
    class           (coolingFunctionClass                   ), pointer                       :: coolingFunction_
    class           (cosmologyFunctionsClass                ), pointer                       :: cosmologyFunctions_
    type            (radiationFieldCosmicMicrowaveBackground), pointer                       :: radiation_
    type            (fgsl_function                          )                                :: integrandFunction
    type            (fgsl_integration_workspace             )                                :: integrationWorkspace
    logical                                                                                  :: integrationReset
    double precision                                                                         :: luminosity                , temperature
    !GCC$ attributes unused :: time, integerProperty, integerBufferCount, integerBuffer, instance
    
    ! Initialize the module.
    call Galacticus_Output_Tree_ICM_Xray_Luminosity_Initialize()
    ! Store property data if we are outputting ICM X-ray luminosities.
    if (outputICMXrayLuminosities) then
       ! Compute required quantities.
       darkMatterHaloScale_       => darkMatterHaloScale      ()
       hotHaloMassDistribution_   => hotHaloMassDistribution  ()
       hotHaloTemperatureProfile_ => hotHaloTemperatureProfile()
       coolingFunction_           => coolingFunction          ()
       cosmologyFunctions_        => cosmologyFunctions       ()
       ! Initialize radiation field.
       allocate(radiation_)
       !# <referenceConstruct object="radiation_" constructor="radiationFieldCosmicMicrowaveBackground(cosmologyFunctions_)"/>
       ! Compute luminosity and temperature.
       integrationReset=.true.
       luminosity      =Integrate(                                                           &
            &                                       0.0d0                                  , &
            &                                       darkMatterHaloScale_%virialRadius(node), &
            &                                       integrandLuminosityXray                , &
            &                                       integrandFunction                      , &
            &                                       integrationWorkspace                   , &
            &                     reset            =integrationReset                       , &
            &                     toleranceAbsolute=0.0d+0                                 , &
            &                     toleranceRelative=1.0d-3                                   &
            &                    )     
       call Integrate_Done(integrandFunction,integrationWorkspace)
       integrationReset=.true.
       temperature     =Integrate(                                                           &
            &                                       0.0d0                                  , &
            &                                       darkMatterHaloScale_%virialRadius(node), &
            &                                       integrandTemperatureXray               , &
            &                                       integrandFunction                      , &
            &                                       integrationWorkspace                   , &
            &                     reset            =integrationReset                       , &
            &                     toleranceAbsolute=0.0d+0                                 , &
            &                     toleranceRelative=1.0d-3                                   &
            &                    )     
       call Integrate_Done(integrandFunction,integrationWorkspace)
       if (luminosity > 0.0d0) then
          temperature=+temperature        &
               &      /luminosity         &
               &      *boltzmannsConstant &
               &      /kilo               &
               &      /electronVolt
       else
          luminosity =+0.0d0
          temperature=+0.0d0
       end if
       ! Store the results.
       doubleProperty  =doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=luminosity
       doubleProperty  =doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=temperature
       !# <objectDestructor name="radiation_"/>
    end if
    return
    
  contains

    double precision function integrandLuminosityXray(radius)
      !% Integrand function used for computing ICM X-ray luminosities.
      use Numerical_Constants_Math        , only : Pi
      use Numerical_Constants_Astronomical
      use Numerical_Constants_Physical
      use Abundances_Structure
      use Chemical_Abundances_Structure
      implicit none
      double precision                      , intent(in   ) :: radius
      class           (nodeComponentHotHalo), pointer       :: hotHalo
      double precision                                      :: density                , temperature       , &
           &                                                   numberDensityHydrogen  , massICM           , &
           &                                                   massToDensityConversion
      type            (abundances          )                :: abundancesICM
      type            (chemicalAbundances  )                :: massChemicalICM        , densityChemicalICM

      ! Get the density of the ICM.
      density    =hotHaloMassDistribution_  %density    (node,radius)
      ! Get the temperature of the ICM.
      temperature=hotHaloTemperatureProfile_%temperature(node,radius)
      ! Get abundances and chemistry of the ICM.
      hotHalo         => node   %hotHalo   ()
      massICM         =  hotHalo%mass      ()
      abundancesICM   =  hotHalo%abundances()
      massChemicalICM =  hotHalo%chemicals ()
      call abundancesICM  %massToMassFraction(           massICM)
      call massChemicalICM%massToNumber      (densityChemicalICM)
      ! Compute factor converting mass of chemicals in (M☉/Mₐₘᵤ) to number density in cm⁻³.
      if (hotHalo%outerRadius() > 0.0d0) then
         massToDensityConversion=Chemicals_Mass_To_Density_Conversion(hotHalo%outerRadius())
      else
         massToDensityConversion=0.0d0
      end if
      ! Convert to number density.
      densityChemicalICM=densityChemicalICM*massToDensityConversion
      ! Compute number density of hydrogen (in cm⁻³).
      numberDensityHydrogen  =+density                                    &
           &                  *abundancesICM   %hydrogenMassFraction()    &
           &                  *massSolar                                  &
           &                  /massHydrogenAtom                           &
           &                  /hecto                                  **3 &
           &                  /megaParsec                             **3
      ! Evaluate the integrand.    
      integrandLuminosityXray=+4.0d0                                                                                                           &
           &                  *Pi                                                                                                              &
           &                  *radius**2                                                                                                       &
           &                  *coolingFunction_%coolingFunction(numberDensityHydrogen,temperature,abundancesICM,densityChemicalICM,radiation_) &
           &                  *(                                                                                                               &
           &                    +megaParsec                                                                                                    &
           &                    /centi                                                                                                         &
           &                   )**3
      return
    end function integrandLuminosityXray

    double precision function integrandTemperatureXray(radius)
      !% Integrand function used for computing ICM X-ray luminosity-weighted temperatures.
      implicit none
      double precision                      , intent(in   ) :: radius

      integrandTemperatureXray=+integrandLuminosityXray                (     radius) &
           &                    *hotHaloTemperatureProfile_%temperature(node,radius)
      return
    end function integrandTemperatureXray
      
  end subroutine Galacticus_Output_Tree_ICM_Xray_Luminosity

end module Galacticus_Output_Trees_ICM_Xray_Luminosity
