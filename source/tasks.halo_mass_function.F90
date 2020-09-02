!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  use :: Cosmological_Density_Field      , only : cosmologicalMassVarianceClass    , criticalOverdensityClass         , haloEnvironmentClass
  use :: Cosmology_Functions             , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters            , only : cosmologyParametersClass
  use :: Dark_Matter_Halo_Biases         , only : darkMatterHaloBiasClass
  use :: Dark_Matter_Halo_Scales         , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profile_Scales      , only : darkMatterProfileScaleRadius     , darkMatterProfileScaleRadiusClass
  use :: Dark_Matter_Profiles_DMO        , only : darkMatterProfileDMOClass
  use :: Halo_Mass_Functions             , only : haloMassFunctionClass
  use :: Linear_Growth                   , only : linearGrowthClass
  use :: Output_Times                    , only : outputTimesClass
  use :: Transfer_Functions              , only : transferFunctionClass
  use :: Unevolved_Subhalo_Mass_Functions, only : unevolvedSubhaloMassFunctionClass
  use :: Virial_Density_Contrast         , only : virialDensityContrastClass

  !# <task name="taskHaloMassFunction">
  !#  <description>A task which computes and outputs the halo mass function and related quantities.</description>
  !# </task>
  type, extends(taskClass) :: taskHaloMassFunction
     !% Implementation of a task which computes and outputs the halo mass function and related quantities.
     private
     class           (cosmologyParametersClass         ), pointer :: cosmologyParameters_                => null()
     class           (cosmologyFunctionsClass          ), pointer :: cosmologyFunctions_                 => null()
     class           (virialDensityContrastClass       ), pointer :: virialDensityContrast_              => null()
     class           (darkMatterProfileDMOClass        ), pointer :: darkMatterProfileDMO_               => null()
     class           (criticalOverdensityClass         ), pointer :: criticalOverdensity_                => null()
     class           (linearGrowthClass                ), pointer :: linearGrowth_                       => null()
     class           (haloMassFunctionClass            ), pointer :: haloMassFunction_                   => null()
     class           (haloEnvironmentClass             ), pointer :: haloEnvironment_                    => null()
     class           (unevolvedSubhaloMassFunctionClass), pointer :: unevolvedSubhaloMassFunction_       => null()
     class           (darkMatterHaloScaleClass         ), pointer :: darkMatterHaloScale_                => null()
     class           (cosmologicalMassVarianceClass    ), pointer :: cosmologicalMassVariance_           => null()
     class           (darkMatterHaloBiasClass          ), pointer :: darkMatterHaloBias_                 => null()
     class           (transferFunctionClass            ), pointer :: transferFunction_                   => null()
     class           (outputTimesClass                 ), pointer :: outputTimes_                        => null()
     class           (darkMatterProfileScaleRadiusClass), pointer :: darkMatterProfileScaleRadius_       => null()
     double precision                                             :: haloMassMinimum                              , haloMassMaximum
     integer                                                      :: pointsPerDecade
     type            (varying_string                   )          :: outputGroup
     logical                                                      :: includeUnevolvedSubhaloMassFunction
     ! Pointer to the parameters for this task.
     type            (inputParameters                  )          :: parameters
  contains
     final     ::            haloMassFunctionDestructor
     procedure :: perform => haloMassFunctionPerform
  end type taskHaloMassFunction

  interface taskHaloMassFunction
     !% Constructors for the {\normalfont \ttfamily haloMassFunction} task.
     module procedure haloMassFunctionConstructorParameters
     module procedure haloMassFunctionConstructorInternal
  end interface taskHaloMassFunction

contains

  function haloMassFunctionConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily haloMassFunction} task class which takes a parameter set as input.
    use :: Galacticus_Nodes, only : nodeClassHierarchyInitialize, treeNode
    use :: Input_Parameters, only : inputParameter              , inputParameters
    use :: Node_Components , only : Node_Components_Initialize
    implicit none
    type            (taskHaloMassFunction             )                        :: self
    type            (inputParameters                  ), intent(inout), target :: parameters
    class           (cosmologyParametersClass         ), pointer               :: cosmologyParameters_
    class           (cosmologyFunctionsClass          ), pointer               :: cosmologyFunctions_
    class           (virialDensityContrastClass       ), pointer               :: virialDensityContrast_
    class           (darkMatterProfileDMOClass        ), pointer               :: darkMatterProfileDMO_
    class           (criticalOverdensityClass         ), pointer               :: criticalOverdensity_
    class           (linearGrowthClass                ), pointer               :: linearGrowth_
    class           (haloMassFunctionClass            ), pointer               :: haloMassFunction_
    class           (haloEnvironmentClass             ), pointer               :: haloEnvironment_
    class           (unevolvedSubhaloMassFunctionClass), pointer               :: unevolvedSubhaloMassFunction_
    class           (darkMatterHaloScaleClass         ), pointer               :: darkMatterHaloScale_
    class           (cosmologicalMassVarianceClass    ), pointer               :: cosmologicalMassVariance_
    class           (darkMatterHaloBiasClass          ), pointer               :: darkMatterHaloBias_
    class           (transferFunctionClass            ), pointer               :: transferFunction_
    class           (outputTimesClass                 ), pointer               :: outputTimes_
    class           (darkMatterProfileScaleRadiusClass), pointer               :: darkMatterProfileScaleRadius_
    type            (inputParameters                  ), pointer               :: parametersRoot
    type            (varying_string                   )                        :: outputGroup
    double precision                                                           :: haloMassMinimum                    , haloMassMaximum
    integer                                                                    :: pointsPerDecade
    logical                                                                    :: includeUnevolvedSubhaloMassFunction
    
    ! Ensure the nodes objects are initialized.
    if (associated(parameters%parent)) then
       parametersRoot => parameters%parent
       do while (associated(parametersRoot%parent))
          parametersRoot => parametersRoot%parent
       end do
       call nodeClassHierarchyInitialize(parametersRoot)
       call Node_Components_Initialize  (parametersRoot)
    else
       parametersRoot => parameters
       call nodeClassHierarchyInitialize(parameters    )
       call Node_Components_Initialize  (parameters    )
    end if
    !# <inputParameter>
    !#   <name>haloMassMinimum</name>
    !#   <defaultValue>1.0d10</defaultValue>
    !#   <description>The minimum mass at which to tabulate halo mass functions.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>haloMassMaximum</name>
    !#   <defaultValue>1.0d15</defaultValue>
    !#   <description>The maximum mass at which to tabulate halo mass functions.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>pointsPerDecade</name>
    !#   <defaultValue>10</defaultValue>
    !#   <description>The number of points per decade of halo mass at which to tabulate halo mass functions.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>outputGroup</name>
    !#   <defaultValue>var_str('.')</defaultValue>
    !#   <description>The HDF5 output group within which to write mass function data.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>includeUnevolvedSubhaloMassFunction</name>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>If true then also compute and output the unevolved subhalo mass function.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters"          name="cosmologyParameters_"          source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"           name="cosmologyFunctions_"           source="parameters"/>
    !# <objectBuilder class="virialDensityContrast"        name="virialDensityContrast_"        source="parameters"/>
    !# <objectBuilder class="darkMatterProfileDMO"         name="darkMatterProfileDMO_"         source="parameters"/>
    !# <objectBuilder class="criticalOverdensity"          name="criticalOverdensity_"          source="parameters"/>
    !# <objectBuilder class="linearGrowth"                 name="linearGrowth_"                 source="parameters"/>
    !# <objectBuilder class="haloMassFunction"             name="haloMassFunction_"             source="parameters"/>
    !# <objectBuilder class="haloEnvironment"              name="haloEnvironment_"              source="parameters"/>
    !# <objectBuilder class="unevolvedSubhaloMassFunction" name="unevolvedSubhaloMassFunction_" source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance"     name="cosmologicalMassVariance_"     source="parameters"/>
    !# <objectBuilder class="darkMatterHaloBias"           name="darkMatterHaloBias_"           source="parameters"/>
    !# <objectBuilder class="transferFunction"             name="transferFunction_"             source="parameters"/>
    !# <objectBuilder class="outputTimes"                  name="outputTimes_"                  source="parameters"/>
    !# <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    self=taskHaloMassFunction(                                     &
         &                    haloMassMinimum                    , &
         &                    haloMassMaximum                    , &
         &                    pointsPerDecade                    , &
         &                    outputGroup                        , &
         &                    includeUnevolvedSubhaloMassFunction, &
         &                    cosmologyParameters_               , &
         &                    cosmologyFunctions_                , &
         &                    virialDensityContrast_             , &
         &                    darkMatterProfileDMO_              , &
         &                    criticalOverdensity_               , &
         &                    linearGrowth_                      , &
         &                    haloMassFunction_                  , &
         &                    haloEnvironment_                   , &
         &                    unevolvedSubhaloMassFunction_      , &
         &                    darkMatterHaloScale_               , &
         &                    darkMatterProfileScaleRadius_      , &
         &                    cosmologicalMassVariance_          , &
         &                    darkMatterHaloBias_                , &
         &                    transferFunction_                  , &
         &                    outputTimes_                       , &
         &                    parametersRoot                       &
         &                   )
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"         />
    !# <objectDestructor name="cosmologyFunctions_"          />
    !# <objectDestructor name="virialDensityContrast_"       />
    !# <objectDestructor name="darkMatterProfileDMO_"        />
    !# <objectDestructor name="criticalOverdensity_"         />
    !# <objectDestructor name="linearGrowth_"                />
    !# <objectDestructor name="haloMassFunction_"            />
    !# <objectDestructor name="haloEnvironment_"             />
    !# <objectDestructor name="unevolvedSubhaloMassFunction_"/>
    !# <objectDestructor name="darkMatterHaloScale_"         />
    !# <objectDestructor name="cosmologicalMassVariance_"    />
    !# <objectDestructor name="darkMatterHaloBias_"          />
    !# <objectDestructor name="transferFunction_"            />
    !# <objectDestructor name="outputTimes_"                 />
    !# <objectDestructor name="darkMatterProfileScaleRadius_"/>
    return
  end function haloMassFunctionConstructorParameters

  function haloMassFunctionConstructorInternal(                                     &
       &                                       haloMassMinimum                    , &
       &                                       haloMassMaximum                    , &
       &                                       pointsPerDecade                    , &
       &                                       outputGroup                        , &
       &                                       includeUnevolvedSubhaloMassFunction, &
       &                                       cosmologyParameters_               , &
       &                                       cosmologyFunctions_                , &
       &                                       virialDensityContrast_             , &
       &                                       darkMatterProfileDMO_              , &
       &                                       criticalOverdensity_               , &
       &                                       linearGrowth_                      , &
       &                                       haloMassFunction_                  , &
       &                                       haloEnvironment_                   , &
       &                                       unevolvedSubhaloMassFunction_      , &
       &                                       darkMatterHaloScale_               , &
       &                                       darkMatterProfileScaleRadius_      , &
       &                                       cosmologicalMassVariance_          , &
       &                                       darkMatterHaloBias_                , &
       &                                       transferFunction_                  , &
       &                                       outputTimes_                       , &
       &                                       parameters                           &
       &                                      ) result(self)
    !% Constructor for the {\normalfont \ttfamily haloMassFunction} task class which takes a parameter set as input.
    implicit none
    type            (taskHaloMassFunction             )                        :: self
    class           (cosmologyParametersClass         ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass          ), intent(in   ), target :: cosmologyFunctions_
    class           (virialDensityContrastClass       ), intent(in   ), target :: virialDensityContrast_
    class           (darkMatterProfileDMOClass        ), intent(in   ), target :: darkMatterProfileDMO_
    class           (criticalOverdensityClass         ), intent(in   ), target :: criticalOverdensity_
    class           (linearGrowthClass                ), intent(in   ), target :: linearGrowth_
    class           (haloMassFunctionClass            ), intent(in   ), target :: haloMassFunction_
    class           (haloEnvironmentClass             ), intent(in   ), target :: haloEnvironment_
    class           (unevolvedSubhaloMassFunctionClass), intent(in   ), target :: unevolvedSubhaloMassFunction_
    class           (darkMatterHaloScaleClass         ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileScaleRadiusClass), intent(in   ), target :: darkMatterProfileScaleRadius_
    class           (cosmologicalMassVarianceClass    ), intent(in   ), target :: cosmologicalMassVariance_
    class           (darkMatterHaloBiasClass          ), intent(in   ), target :: darkMatterHaloBias_
    class           (transferFunctionClass            ), intent(in   ), target :: transferFunction_
    class           (outputTimesClass                 ), intent(in   ), target :: outputTimes_
    type            (varying_string                   ), intent(in   )         :: outputGroup
    double precision                                   , intent(in   )         :: haloMassMinimum                    , haloMassMaximum
    integer                                            , intent(in   )         :: pointsPerDecade
    logical                                            , intent(in   )         :: includeUnevolvedSubhaloMassFunction
    type            (inputParameters                  ), intent(in   ), target :: parameters
    !# <constructorAssign variables="haloMassMinimum,haloMassMaximum,pointsPerDecade,outputGroup,includeUnevolvedSubhaloMassFunction,*cosmologyParameters_,*cosmologyFunctions_,*virialDensityContrast_,*darkMatterProfileDMO_,*criticalOverdensity_,*linearGrowth_,*haloMassFunction_,*haloEnvironment_,*unevolvedSubhaloMassFunction_,*darkMatterHaloScale_, *darkMatterProfileScaleRadius_, *cosmologicalMassVariance_,*darkMatterHaloBias_,*transferFunction_, *outputTimes_"/>

    self%parameters=inputParameters(parameters)
    call self%parameters%parametersGroupCopy(parameters)
    return
  end function haloMassFunctionConstructorInternal

  subroutine haloMassFunctionDestructor(self)
    !% Destructor for the {\normalfont \ttfamily haloMassFunction} task class.
    use :: Node_Components, only : Node_Components_Uninitialize
    implicit none
    type(taskHaloMassFunction), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"         />
    !# <objectDestructor name="self%cosmologyFunctions_"          />
    !# <objectDestructor name="self%virialDensityContrast_"       />
    !# <objectDestructor name="self%darkMatterProfileDMO_"        />
    !# <objectDestructor name="self%criticalOverdensity_"         />
    !# <objectDestructor name="self%linearGrowth_"                />
    !# <objectDestructor name="self%haloMassFunction_"            />
    !# <objectDestructor name="self%haloEnvironment_"             />
    !# <objectDestructor name="self%unevolvedSubhaloMassFunction_"/>
    !# <objectDestructor name="self%darkMatterHaloScale_"         />
    !# <objectDestructor name="self%darkMatterProfileScaleRadius_"/>
    !# <objectDestructor name="self%cosmologicalMassVariance_"    />
    !# <objectDestructor name="self%darkMatterHaloBias_"          />
    !# <objectDestructor name="self%transferFunction_"            />
    !# <objectDestructor name="self%outputTimes_"                 />
    call Node_Components_Uninitialize()
    return
  end subroutine haloMassFunctionDestructor

  subroutine haloMassFunctionPerform(self,status)
    !% Compute and output the halo mass function.
    use            :: Galacticus_Calculations_Resets  , only : Galacticus_Calculations_Reset
    use            :: Galacticus_Display              , only : Galacticus_Display_Indent        , Galacticus_Display_Unindent
    use            :: Galacticus_Error                , only : errorStatusSuccess
    use            :: Galacticus_HDF5                 , only : galacticusOutputFile
    use            :: Galacticus_Nodes                , only : mergerTree                       , nodeComponentBasic                 , nodeComponentDarkMatterProfile, treeNode
    use            :: IO_HDF5                         , only : hdf5Object
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Memory_Management               , only : allocateArray
    use            :: Node_Components                 , only : Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize
    use            :: Numerical_Constants_Astronomical, only : massSolar                        , megaParsec
    use            :: Numerical_Constants_Math        , only : Pi
    use            :: Numerical_Constants_Prefixes    , only : kilo
    use            :: Numerical_Integration           , only : integrator                       , GSL_Integ_Gauss15
    use            :: Numerical_Ranges                , only : Make_Range                       , rangeTypeLogarithmic
    use            :: String_Handling                 , only : operator(//)
    implicit none
    class           (taskHaloMassFunction          ), intent(inout), target         :: self
    integer                                         , intent(  out), optional       :: status
    double precision                                , allocatable  , dimension(:,:) :: massFunctionDifferentialLogarithmicBinAveraged       , biasHalo                     , &
         &                                                                             massFunctionCumulative                               , massFunctionDifferential     , &
         &                                                                             massFunctionDifferentialLogarithmic                  , massFunctionMassFraction     , &
         &                                                                             peakHeight                                           , densityFieldRootVariance     , &
         &                                                                             radiusVirial                                         , temperatureVirial            , &
         &                                                                             velocityVirial                                       , darkMatterProfileRadiusScale , &
         &                                                                             velocityMaximum                                      , peakHeightMassFunction       , &
         &                                                                             densityFieldRootVarianceGradientLogarithmic          , massFunctionCumulativeSubhalo
    double precision                                , allocatable  , dimension(:  ) :: outputCharacteristicMass                             , outputCriticalOverdensities  , &
         &                                                                             outputExpansionFactors                               , outputGrowthFactors          , &
         &                                                                             outputRedshifts                                      , outputTimes                  , &
         &                                                                             outputVirialDensityContrast                          , outputTurnAroundRadius       , &
         &                                                                             massHalo
    ! The upper limit to halo mass used when computing cumulative mass functions.
    double precision                                , parameter                     :: haloMassEffectiveInfinity                     =1.0d16
    ! Largest subhalo mass (in units of host mass) for which we expect significant unevolved subhalo mass function.
    double precision                                , parameter                     :: subhaloMassMaximum                            =1.0d02
    class           (nodeComponentBasic            ), pointer                       :: basic
    class           (nodeComponentDarkMatterProfile), pointer                       :: darkMatterProfileHalo
    type            (mergerTree                    ), target                        :: tree
    type            (integrator                    )                                :: integrator_
    integer         (c_size_t                      )                                :: iOutput                                              , outputCount                  , &
         &                                                                             iMass                                                , massCount
    double precision                                                                :: massHaloBinMinimum                                   , massHaloBinMaximum           , &
         &                                                                             massHaloLogarithmicInterval                          , massHalfMode                 , &
         &                                                                             wavenumberComoving                                   , growthFactor
    type            (hdf5Object                    )                                :: outputsGroup                                         , outputGroup                  , &
         &                                                                             containerGroup                                       , powerSpectrumGroup           , &
         &                                                                             cosmologyGroup                                       , dataset
    integer                                                                         :: statusHalfModeMass
    type            (varying_string                )                                :: groupName                                            , commentText

    call Galacticus_Display_Indent('Begin task: halo mass function')
    ! Call routines to perform initializations which must occur for all threads if run in parallel.
    call Node_Components_Thread_Initialize(self%parameters)
    ! Get the requested output redshifts.
    outputCount=self%outputTimes_%count()
    call allocateArray(outputTimes                                   ,[          outputCount])
    call allocateArray(outputRedshifts                               ,[          outputCount])
    call allocateArray(outputExpansionFactors                        ,[          outputCount])
    call allocateArray(outputGrowthFactors                           ,[          outputCount])
    call allocateArray(outputCriticalOverdensities                   ,[          outputCount])
    call allocateArray(outputVirialDensityContrast                   ,[          outputCount])
    call allocateArray(outputTurnaroundRadius                        ,[          outputCount])
    call allocateArray(outputCharacteristicMass                      ,[          outputCount])
    ! Compute number of tabulation points.
    massCount=int(log10(self%haloMassMaximum/self%haloMassMinimum)*dble(self%pointsPerDecade))+1
    call allocateArray(massHalo                                      ,[massCount            ])
    call allocateArray(massFunctionDifferential                      ,[massCount,outputCount])
    call allocateArray(massFunctionDifferentialLogarithmic           ,[massCount,outputCount])
    call allocateArray(massFunctionDifferentialLogarithmicBinAveraged,[massCount,outputCount])
    call allocateArray(massFunctionCumulative                        ,[massCount,outputCount])
    call allocateArray(massFunctionCumulativeSubhalo                 ,[massCount,outputCount])
    call allocateArray(massFunctionMassFraction                      ,[massCount,outputCount])
    call allocateArray(biasHalo                                      ,[massCount,outputCount])
    call allocateArray(densityFieldRootVariance                      ,[massCount,outputCount])
    call allocateArray(densityFieldRootVarianceGradientLogarithmic   ,[massCount,outputCount])
    call allocateArray(peakHeight                                    ,[massCount,outputCount])
    call allocateArray(peakHeightMassFunction                        ,[massCount,outputCount])
    call allocateArray(velocityVirial                                ,[massCount,outputCount])
    call allocateArray(temperatureVirial                             ,[massCount,outputCount])
    call allocateArray(radiusVirial                                  ,[massCount,outputCount])
    call allocateArray(darkMatterProfileRadiusScale                  ,[massCount,outputCount])
    call allocateArray(velocityMaximum                               ,[massCount,outputCount])
    ! Compute output time properties.
    do iOutput=1,outputCount
       outputTimes                   (iOutput)=self%outputTimes_%time                                 (                                                                       iOutput )
       outputExpansionFactors        (iOutput)=self%cosmologyFunctions_   %expansionFactor            (                                                outputTimes           (iOutput))
       outputRedshifts               (iOutput)=self%cosmologyFunctions_   %redshiftFromExpansionFactor(                                                outputExpansionFactors(iOutput))
       outputGrowthFactors           (iOutput)=self%linearGrowth_         %value                      (                                                outputTimes           (iOutput))
       outputCharacteristicMass      (iOutput)=self%criticalOverdensity_  %collapsingMass             (                                                outputTimes           (iOutput))
       if (outputCharacteristicMass(iOutput) > 0.0d0) then
          outputCriticalOverdensities(iOutput)=self%criticalOverdensity_  %value                      (mass=outputCharacteristicMass    (iOutput),time=outputTimes           (iOutput))
       else
          outputCriticalOverdensities(iOutput)=0.0d0
       end if
       outputVirialDensityContrast   (iOutput)=self%virialDensityContrast_%densityContrast            (mass=self%haloMassMinimum                 ,time=outputTimes           (iOutput))
       outputTurnaroundRadius        (iOutput)=self%virialDensityContrast_%turnAroundOverVirialRadii  (mass=self%haloMassMinimum                 ,time=outputTimes           (iOutput))
    end do
    ! Create a node object, assume zero environmental overdensity.
    tree%baseNode          => treeNode()
    tree%baseNode%hostTree => tree
    call tree%properties      %initialize          (                               )
    call self%haloEnvironment_%overdensityLinearSet(tree%baseNode,overdensity=0.0d0)
    ! Build a range of halo masses.
    massHalo                   =Make_Range(self%haloMassMinimum,self%haloMassMaximum,int(massCount),rangeTypeLogarithmic)
    massHaloLogarithmicInterval=log(self%haloMassMaximum/self%haloMassMinimum)/dble(massCount-1)
    ! Get the basic and dark matter profile components.
    basic                 => tree%baseNode%basic            (autoCreate=.true.)
    darkMatterProfileHalo => tree%baseNode%darkMatterProfile(autoCreate=.true.)
    ! Build an integrator.
    integrator_=integrator(subhaloMassFunctionIntegrand,toleranceRelative=1.0d-3,integrationRule=GSL_Integ_Gauss15)
    ! Iterate over all output times.    
    do iOutput=outputCount,1,-1
       ! Set the time in the node.
       call basic%timeSet(outputTimes(iOutput))
       ! Loop over all halo masses.
       do iMass=1,massCount
          ! Reset calculations.
          call Galacticus_Calculations_Reset(tree%baseNode)
          ! Set the mass in the node.
          call basic                %massSet (massHalo                                 (iMass        ))
          ! Set the node scale radius.
          call darkMatterProfileHalo%scaleSet(self%darkMatterProfileScaleRadius_%radius(tree%baseNode))
          ! Compute bin interval.
          massHaloBinMinimum=massHalo(iMass)*exp(-0.5*massHaloLogarithmicInterval)
          massHaloBinMaximum=massHalo(iMass)*exp(+0.5*massHaloLogarithmicInterval)
          ! Get the growth factor and this epoch and mass.
          wavenumberComoving=2.0d0*Pi/(3.0d0*massHalo(iMass)/4.0d0/Pi/self%cosmologyFunctions_%matterDensityEpochal(outputTimes(iOutput)))**(1.0d0/3.0d0)*self%cosmologyFunctions_%expansionFactor(outputTimes(iOutput))
          growthFactor      =self%linearGrowth_%value(outputTimes(iOutput),wavenumber=wavenumberComoving)
          ! Compute halo properties.
          densityFieldRootVariance                      (iMass,iOutput)=+self%cosmologicalMassVariance_    %rootVariance                   (mass   =massHalo          (iMass)                                   ,time=outputTimes(iOutput)                   )
          densityFieldRootVarianceGradientLogarithmic   (iMass,iOutput)=+self%cosmologicalMassVariance_    %rootVarianceLogarithmicGradient(mass   =massHalo          (iMass)                                   ,time=outputTimes(iOutput)                   )
          peakHeight                                    (iMass,iOutput)=+self%criticalOverdensity_         %value                          (mass   =massHalo          (iMass)                                   ,time=outputTimes(iOutput)                   )    &
               &                                                        /densityFieldRootVariance                                          (                           iMass                                    ,                 iOutput)
          massFunctionDifferentialLogarithmicBinAveraged(iMass,iOutput)=+self%haloMassFunction_            %integrated                     (massLow=massHaloBinMinimum       ,massHigh=massHaloBinMaximum       ,time=outputTimes(iOutput),node=tree%baseNode)    &
               &                                                        /massHaloLogarithmicInterval
          massFunctionDifferential                      (iMass,iOutput)=+self%haloMassFunction_            %differential                   (mass   =massHalo          (iMass)                                   ,time=outputTimes(iOutput),node=tree%baseNode)
          massFunctionCumulative                        (iMass,iOutput)=+self%haloMassFunction_            %integrated                     (massLow=massHalo          (iMass),massHigh=haloMassEffectiveInfinity,time=outputTimes(iOutput),node=tree%baseNode)
          massFunctionMassFraction                      (iMass,iOutput)=+self%haloMassFunction_            %massFraction                   (massLow=massHalo          (iMass),massHigh=haloMassEffectiveInfinity,time=outputTimes(iOutput),node=tree%baseNode)
          peakHeightMassFunction                        (iMass,iOutput)=+massHalo                                                          (                           iMass                                                                                 )**2 &
               &                                                        *massFunctionDifferential                                          (                           iMass                                ,                     iOutput                    )    &
               &                                                        /self%cosmologyParameters_         %densityCritical                (                                                                                                                 )    &
               &                                                        /self%cosmologyParameters_         %OmegaMatter                    (                                                                                                                 )    &
               &                                                        /abs(                                                                                                                                                                                     &
               &                                                             self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient(mass   =massHalo          (iMass)                                   ,time=outputTimes(iOutput)                   )    &
               &                                                            )
          biasHalo                                      (iMass,iOutput)=self%darkMatterHaloBias_           %bias                           (                                                                                               node=tree%baseNode)
          velocityVirial                                (iMass,iOutput)=self%darkMatterHaloScale_          %virialVelocity                 (                                                                                               node=tree%baseNode)
          temperatureVirial                             (iMass,iOutput)=self%darkMatterHaloScale_          %virialTemperature              (                                                                                               node=tree%baseNode)
          radiusVirial                                  (iMass,iOutput)=self%darkMatterHaloScale_          %virialRadius                   (                                                                                               node=tree%baseNode)
          velocityMaximum                               (iMass,iOutput)=self%darkMatterProfileDMO_         %circularVelocityMaximum        (                                                                                               node=tree%baseNode)
          darkMatterProfileRadiusScale                  (iMass,iOutput)=darkMatterProfileHalo              %scale                          (                                                                                                                 )
          ! Integrate the unevolved subhalo mass function over the halo mass function to get the total subhalo mass function.
          if (self%includeUnevolvedSubhaloMassFunction)                                                                                   &
               & massFunctionCumulativeSubhalo          (iMass,iOutput)=integrator_%integrate(                                            &
               &                                                                              log(massHalo(1)/subhaloMassMaximum       ), &
               &                                                                              log(            haloMassEffectiveInfinity)  &
               &                                                                             )
       end do
       massFunctionDifferentialLogarithmic(:,iOutput)=+massFunctionDifferential(:,iOutput) &
            &                                         *massHalo
    end do
    ! Open the group for output time information.
    if (self%outputGroup == ".") then
       outputsGroup  =galacticusOutputFile%openGroup(     'Outputs'        ,'Group containing datasets relating to output times.')
    else
       containerGroup=galacticusOutputFile%openGroup(char(self%outputGroup),'Group containing halo mass function data.'          )
       outputsGroup  =containerGroup      %openGroup(     'Outputs'        ,'Group containing datasets relating to output times.')
    end if
    ! Store half-mode mass if possible.
    massHalfMode=self%transferFunction_%halfModeMass(statusHalfModeMass)
    if (statusHalfModeMass == errorStatusSuccess) then
       if (self%outputGroup == ".") then
          powerSpectrumGroup=galacticusOutputFile%openGroup('powerSpectrum','Group containing data relating to the power spectrum.')
       else
          powerSpectrumGroup=containerGroup      %openGroup('powerSpectrum','Group containing data relating to the power spectrum.')
       end if
       call powerSpectrumGroup%writeAttribute(massHalfMode,'massHalfMode')
       call powerSpectrumGroup%close         (                           )
    end if
    ! Store other usual information.
    if (self%outputGroup == ".") then
       cosmologyGroup=galacticusOutputFile%openGroup('cosmology','Group containing data relating to cosmology.')
    else
       cosmologyGroup=containerGroup      %openGroup('cosmology','Group containing data relating to cosmology.')
    end if
    call cosmologyGroup%writeAttribute(self%cosmologyParameters_%densityCritical(),'densityCritical')
    call cosmologyGroup%close()
    ! Iterate over output times and output data.
    do iOutput=1,outputCount
       groupName  ='Output'
       commentText='Data for output number '
       groupName  =groupName  //iOutput
       commentText=commentText//iOutput
       outputGroup=outputsGroup%openGroup(char(groupName),char(commentText))
       call    outputGroup%writeAttribute(outputTimes                                   (  iOutput),'outputTime'                                                                                                                          )
       call    outputGroup%writeAttribute(outputRedshifts                               (  iOutput),'outputRedshift'                                                                                                                      )
       call    outputGroup%writeAttribute(outputExpansionFactors                        (  iOutput),'outputExpansionFactor'                                                                                                               )
       call    outputGroup%writeAttribute(outputGrowthFactors                           (  iOutput),'growthFactor'                                                                                                                        )
       call    outputGroup%writeAttribute(outputCriticalOverdensities                   (  iOutput),'criticalOverdensity'                                                                                                                 )
       call    outputGroup%writeAttribute(outputVirialDensityContrast                   (  iOutput),'virialDensityContrast'                                                                                                               )
       call    outputGroup%writeAttribute(outputTurnaroundRadius                        (  iOutput),'turnaroundToVirialRadiusRatio'                                                                                                       )
       call    outputGroup%writeAttribute(outputCharacteristicMass                      (  iOutput),'massHaloCharacteristic'                                                                                                              )
       call    outputGroup%writeDataset  (massHalo                                      (:        ),'haloMass'                      ,'The mass of the halo.'                                                      ,datasetReturned=dataset)
       call    dataset    %writeAttribute(massSolar                                                ,'unitsInSI'                                                                                                                           )
       call    dataset    %close         (                                                                                                                                                                                                )
       call    outputGroup%writeDataset  (massFunctionDifferential                      (:,iOutput),'haloMassFunctionM'             ,'The halo mass function (per unit halo mass).'                               ,datasetReturned=dataset)
       call    dataset    %writeAttribute(1.0d0/megaParsec**3/massSolar                            ,'unitsInSI'                                                                                                                           )
       call    dataset    %close         (                                                                                                                                                                                                )
       call    outputGroup%writeDataset  (massFunctionDifferentialLogarithmic           (:,iOutput),'haloMassFunctionLnM'           ,'The halo mass function (per logarithmic halo mass).'                        ,datasetReturned=dataset)
       call    dataset    %writeAttribute(1.0d0/megaParsec**3                                      ,'unitsInSI'                                                                                                                           )
       call    dataset    %close         (                                                                                                                                                                                                )
       call    outputGroup%writeDataset  (massFunctionDifferentialLogarithmicBinAveraged(:,iOutput),'haloMassFunctionLnMBinAveraged','The halo mass function (per logarithmic halo mass averaged across the bin).',datasetReturned=dataset)
       call    dataset    %writeAttribute(1.0d0/megaParsec**3                                      ,'unitsInSI'                                                                                                                           )
       call    dataset    %close         (                                                                                                                                                                                                )
       call    outputGroup%writeDataset  (massFunctionCumulative                        (:,iOutput),'haloMassFunctionCumulative'    ,'The halo cumulative mass function.'                                         ,datasetReturned=dataset)
       call    dataset    %writeAttribute(1.0d0/megaParsec**3                                      ,'unitsInSI'                                                                                                                           )
       call    dataset    %close         (                                                                                                                                                                                                )
       if (self%includeUnevolvedSubhaloMassFunction) then
          call outputGroup%writeDataset  (massFunctionCumulativeSubhalo                 (:,iOutput),'subhaloMassFunctionCumulative' ,'The subhalo cumulative mass function.'                                      ,datasetReturned=dataset)
          call dataset    %writeAttribute(1.0d0/megaParsec**3                                      ,'unitsInSI'                                                                                                                           )
          call dataset    %close         (                                                                                                                                                                                                )
       end if
       call    outputGroup%writeDataset  (velocityVirial                                (:,iOutput),'haloVirialVelocity'            ,'The virial velocity of halos.'                                              ,datasetReturned=dataset)
       call    dataset    %writeAttribute(kilo                                                     ,'unitsInSI'                                                                                                                           )
       call    dataset    %close         (                                                                                                                                                                                                )
       call    outputGroup%writeDataset  (temperatureVirial                             (:,iOutput),'haloVirialTemperature'         ,'The virial temperature of halos.'                                           ,datasetReturned=dataset)
       call    dataset    %writeAttribute(1.0d0                                                    ,'unitsInSI'                                                                                                                           )
       call    dataset    %close         (                                                                                                                                                                                                )
       call    outputGroup%writeDataset  (radiusVirial                                  (:,iOutput),'haloVirialRadius'              ,'The virial radius of halos.'                                                ,datasetReturned=dataset)
       call    dataset    %writeAttribute(megaParsec                                               ,'unitsInSI'                                                                                                                           )
       call    dataset    %close         (                                                                                                                                                                                                )
       call    outputGroup%writeDataset  (darkMatterProfileRadiusScale                  (:,iOutput),'haloScaleRadius'               ,'The scale radius of halos.'                                                 ,datasetReturned=dataset)
       call    dataset    %writeAttribute(megaParsec                                               ,'unitsInSI'                                                                                                                           )
       call    dataset    %close         (                                                                                                                                                                                                )
       call    outputGroup%writeDataset  (velocityMaximum                               (:,iOutput),'haloVelocityMaximum'           ,'The maximum circular velocity of halos.'                                    ,datasetReturned=dataset)
       call    dataset    %writeAttribute(kilo                                                     ,'unitsInSI'                                                                                                                           )
       call    dataset    %close         (                                                                                                                                                                                                )
       call    outputGroup%writeDataset  (biasHalo                                      (:,iOutput),'haloBias'                      ,'The large scale linear bias of halos.'                                                              )
       call    outputGroup%writeDataset  (densityFieldRootVariance                      (:,iOutput),'haloSigma'                     ,'The mass fluctuation on the scale of the halo.'                                                     )
       call    outputGroup%writeDataset  (densityFieldRootVarianceGradientLogarithmic   (:,iOutput),'haloAlpha'                     ,'dlog(σ)/dlog(m).'                                                                                   )
       call    outputGroup%writeDataset  (peakHeight                                    (:,iOutput),'haloPeakHeightNu'              ,'The peak height, ν, of the halo.'                                                                   )
       call    outputGroup%writeDataset  (massFunctionMassFraction                      (:,iOutput),'haloMassFractionCumulative'    ,'The halo cumulative mass fraction.'                                                                 )
       call    outputGroup%writeDataset  (peakHeightMassFunction                        (:,iOutput),'haloMassFunctionNuFNu'         ,'The halo mass fraction function as a function of the ν parameter, νF(ν).'                           )
       call    outputGroup%close         (                                                                                                                                                                                                )
    end do
    call outputsGroup%close()
    if (containerGroup%isOpen()) call containerGroup%close()
    if (present(status)) status=errorStatusSuccess
    call Node_Components_Thread_Uninitialize()
    call Galacticus_Display_Unindent('Done task: halo mass function' )
    return

  contains

    double precision function subhaloMassFunctionIntegrand(logMass)
      !% Integrand function used to find the cumulative subhalo mass function.
      implicit none
      double precision, intent(in   ) :: logMass
      double precision                :: mass

      ! Extract integrand parameters.
      mass=exp(logMass)
      ! Return the differential halo mass function multiplied by the integrated unevolved subhalo mass function in such hosts.
      subhaloMassFunctionIntegrand=+                                                                                                               mass                     &
           &                       *self%haloMassFunction_            %differential(outputTimes(iOutput)                                          ,mass,node=tree%baseNode) &
           &                       *self%unevolvedSubhaloMassFunction_%integrated  (outputTimes(iOutput),massHalo(iMass),haloMassEffectiveInfinity,mass                   )
      return
    end function subhaloMassFunctionIntegrand

  end subroutine haloMassFunctionPerform
