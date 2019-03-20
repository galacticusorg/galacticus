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

!% Contains a module which handles outputting of Sunyaev-Zeldovich data to the \glc\ output file.

module Galacticus_Output_Trees_ICM_SZ
  !% Handles outputting of ICM Sunyaev-Zeldovich data to the \glc\ output file.
  use ISO_Varying_String
  implicit none
  private
  public :: Galacticus_Output_Tree_ICM_SZ      , Galacticus_Output_Tree_ICM_SZ_Property_Count, &
       &    Galacticus_Output_Tree_ICM_SZ_Names

  ! Flag indicating whether or not ICM SZ luminosity information is to be output.
  logical :: outputICMSZ

  ! Flag indicating whether or not this module has been initialized.
  logical :: outputICMSZInitialized=.false.

contains

  subroutine Galacticus_Output_Tree_ICM_SZ_Initialize()
    !% Initializes the module by determining whether or not ICM SZ data should be output.
    use Input_Parameters
    implicit none

    if (.not.outputICMSZInitialized) then
       !$omp critical(Galacticus_Output_Tree_ICM_SZ_Initialize)
       if (.not.outputICMSZInitialized) then
          !# <inputParameter>
          !#   <name>outputICMSZ</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>.false.</defaultValue>
          !#   <description>Specifies whether or not ICM SZ data should be included in the output file.</description>
          !#   <group>output</group>
          !#   <source>globalParameters</source>
          !#   <type>boolean</type>
          !# </inputParameter>
          ! Flag that module is now initialized.
          outputICMSZInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_ICM_SZ_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_ICM_SZ_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_ICM_SZ_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_ICM_SZ</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_ICM_SZ_Names(node,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
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
    call Galacticus_Output_Tree_ICM_SZ_Initialize()

    ! Return property names if we are outputting ICM SZ data.
    if (outputICMSZ) then
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='icmSZComptonYMean'
       doublePropertyComments(doubleProperty)='Mean thermal Sunyaev-Zeldovich Compton y-parameter of the ICM within the virial radius []'
       doublePropertyUnitsSI (doubleProperty)=1.0d0
    end if
    return
  end subroutine Galacticus_Output_Tree_ICM_SZ_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_ICM_SZ_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_ICM_SZ</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_ICM_SZ_Property_Count(node,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of projected density properties to be written to the \glc\ output file.
    use Galacticus_Nodes, only : treeNode
    implicit none
    type            (treeNode), intent(inout) :: node
    double precision          , intent(in   ) :: time
    integer                   , intent(inout) :: doublePropertyCount, integerPropertyCount
    !GCC$ attributes unused :: node, time, integerPropertyCount
    
    ! Initialize the module.
    call Galacticus_Output_Tree_ICM_SZ_Initialize()

    ! Increment property count if we are outputting ICM SZ data.
    if (outputICMSZ) doublePropertyCount=doublePropertyCount+1
    return
  end subroutine Galacticus_Output_Tree_ICM_SZ_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_ICM_SZ</unitName>
  !#  <sortName>Galacticus_Output_Tree_ICM_SZ</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_ICM_SZ(node,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time,instance)
    !% Store projected density properties in the \glc\ output file buffers.
    use FGSL                             , only : fgsl_function                          , fgsl_integration_workspace
    use Numerical_Integration            , only : Integrate                              , Integrate_Done
    use Galacticus_Nodes                 , only : treeNode                               , nodeComponentHotHalo
    use Kind_Numbers                     , only : kind_int8
    use Multi_Counters                   , only : multiCounter
    use Dark_Matter_Halo_Scales          , only : darkMatterHaloScale                    , darkMatterHaloScaleClass
    use Hot_Halo_Mass_Distributions      , only : hotHaloMassDistribution                , hotHaloMassDistributionClass
    use Hot_Halo_Temperature_Profiles    , only : hotHaloTemperatureProfile              , hotHaloTemperatureProfileClass
    use Cosmology_Functions              , only : cosmologyFunctions                     , cosmologyFunctionsClass
    use Chemical_States                  , only : chemicalState                          , chemicalStateClass
    use Radiation_Fields                 , only : radiationFieldCosmicMicrowaveBackground
    use Numerical_Constants_Math         , only : Pi
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
    class           (cosmologyFunctionsClass                ), pointer                       :: cosmologyFunctions_
    class           (chemicalStateClass                     ), pointer                       :: chemicalState_
    type            (radiationFieldCosmicMicrowaveBackground), pointer                       :: radiation_
    type            (fgsl_function                          )                                :: integrandFunction
    type            (fgsl_integration_workspace             )                                :: integrationWorkspace
    logical                                                                                  :: integrationReset
    !GCC$ attributes unused :: time, integerProperty, integerBufferCount, integerBuffer, instance
    
    ! Initialize the module.
    call Galacticus_Output_Tree_ICM_SZ_Initialize()
    ! Store property data if we are outputting ICM SZ data.
    if (outputICMSZ) then
       ! Compute required quantities.
       darkMatterHaloScale_       => darkMatterHaloScale      ()
       hotHaloMassDistribution_   => hotHaloMassDistribution  ()
       hotHaloTemperatureProfile_ => hotHaloTemperatureProfile()
       cosmologyFunctions_        => cosmologyFunctions       ()
       chemicalState_             => chemicalState            ()
       ! Initialize radiation field.
       allocate(radiation_)
       !# <referenceConstruct object="radiation_" constructor="radiationFieldCosmicMicrowaveBackground(cosmologyFunctions_)"/>
       integrationReset=.true.
       doubleProperty                                  =+doubleProperty+1
       doubleBuffer  (doubleBufferCount,doubleProperty)=+Integrate(                                                             &
            &                                                                       0.0d0                                     , &
            &                                                                       darkMatterHaloScale_%virialRadius(node)   , &
            &                                                                       integrandComptionY                        , &
            &                                                                       integrandFunction                         , &
            &                                                                       integrationWorkspace                      , &
            &                                                     reset            =integrationReset                          , &
            &                                                     toleranceAbsolute=0.0d+0                                    , &
            &                                                     toleranceRelative=1.0d-3                                      &
            &                                                    )                                                              &
            &                                           /Pi                                                                     &
            &                                           /                           darkMatterHaloScale_%virialRadius(node)**2
       call Integrate_Done(integrandFunction,integrationWorkspace)
       !# <objectDestructor name="radiation_"/>
    end if
    return
    
  contains

    double precision function integrandComptionY(radius)
      !% Integrand function used for computing ICM SZ properties.
      use Numerical_Constants_Astronomical, only : megaParsec        , massSolar
      use Numerical_Constants_Physical    , only : boltzmannsConstant, thomsonCrossSection, speedLight, electronMass
      use Numerical_Constants_Prefixes    , only : centi             , hecto
      use Numerical_Constants_Atomic      , only : massHydrogenAtom
      use Abundances_Structure            , only : abundances
      implicit none
      double precision                      , intent(in   ) :: radius
      class           (nodeComponentHotHalo), pointer       :: hotHalo
      double precision                                      :: density              , temperature, &
           &                                                   numberDensityHydrogen, massICM
      type            (abundances          )                :: abundancesICM

      ! Get the density of the ICM.
      density    =hotHaloMassDistribution_  %density    (node,radius)
      ! Get the temperature of the ICM.
      temperature=hotHaloTemperatureProfile_%temperature(node,radius)
      ! Get abundances and chemistry of the ICM.
      hotHalo         => node   %hotHalo   ()
      massICM         =  hotHalo%mass      ()
      abundancesICM   =  hotHalo%abundances()
      call abundancesICM%massToMassFraction(massICM)
      ! Compute number density of hydrogen (in cm⁻³).
      numberDensityHydrogen  =+density                                    &
           &                  *abundancesICM   %hydrogenMassFraction()    &
           &                  *massSolar                                  &
           &                  /massHydrogenAtom                           &
           &                  /hecto                                  **3 &
           &                  /megaParsec                             **3
      ! Evaluate the integrand.    
      integrandComptionY=+4.0d0                                                                                         &
           &             *Pi                                                                                            &
           &             *radius                                                                                    **2 &
           &             *boltzmannsConstant                                                                            &
           &             *thomsonCrossSection                                                                           &
           &             /electronMass                                                                                  &
           &             /speedLight                                                                                **2 &
           &             *chemicalState_%electronDensity(numberDensityHydrogen,temperature,abundancesICM,radiation_)    &
           &             *temperature                                                                                   &
           &             *megaParsec                                                                                    &
           &             /centi                                                                                     **3
      return
    end function integrandComptionY

  end subroutine Galacticus_Output_Tree_ICM_SZ  

end module Galacticus_Output_Trees_ICM_SZ
