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

  use Dark_Matter_Profile_Scales      , only : darkMatterProfileScaleRadius, darkMatterProfileScaleRadiusClass
  use Halo_Mass_Functions
  use Unevolved_Subhalo_Mass_Functions
  use Dark_Matter_Halo_Biases
  use Dark_Matter_Profiles
  use Dark_Matter_Halo_Scales
  use Cosmology_Parameters
  use Cosmology_Functions
  use Linear_Growth
  use Virial_Density_Contrast
  use Transfer_Functions
  use Cosmological_Density_Field
  use Output_Times

  !# <task name="taskHaloMassFunction">
  !#  <description>A task which computes and outputs the halo mass function and related quantities.</description>
  !# </task>
  type, extends(taskClass) :: taskHaloMassFunction
     !% Implementation of a task which computes and outputs the halo mass function and related quantities.
     private
     class           (cosmologyParametersClass         ), pointer :: cosmologyParameters_ => null()
     class           (cosmologyFunctionsClass          ), pointer :: cosmologyFunctions_ => null()
     class           (virialDensityContrastClass       ), pointer :: virialDensityContrast_ => null()
     class           (darkMatterProfileClass           ), pointer :: darkMatterProfile_ => null()
     class           (criticalOverdensityClass         ), pointer :: criticalOverdensity_ => null()
     class           (linearGrowthClass                ), pointer :: linearGrowth_ => null()
     class           (haloMassFunctionClass            ), pointer :: haloMassFunction_ => null()
     class           (haloEnvironmentClass             ), pointer :: haloEnvironment_ => null()
     class           (unevolvedSubhaloMassFunctionClass), pointer :: unevolvedSubhaloMassFunction_ => null()
     class           (darkMatterHaloScaleClass         ), pointer :: darkMatterHaloScale_ => null()
     class           (cosmologicalMassVarianceClass    ), pointer :: cosmologicalMassVariance_ => null()
     class           (darkMatterHaloBiasClass          ), pointer :: darkMatterHaloBias_ => null()
     class           (transferFunctionClass            ), pointer :: transferFunction_ => null()
     class           (outputTimesClass                 ), pointer :: outputTimes_ => null()
     class           (darkMatterProfileScaleRadiusClass), pointer :: darkMatterProfileScaleRadius_ => null()
     double precision                                             :: haloMassMinimum              , haloMassMaximum
     integer                                                      :: pointsPerDecade
     type            (varying_string                   )          :: outputGroup
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
    use Galacticus_Nodes
    use Node_Components
    use Input_Parameters
    implicit none
    type            (taskHaloMassFunction             )                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    class           (cosmologyParametersClass         ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    class           (virialDensityContrastClass       ), pointer       :: virialDensityContrast_
    class           (darkMatterProfileClass           ), pointer       :: darkMatterProfile_
    class           (criticalOverdensityClass         ), pointer       :: criticalOverdensity_
    class           (linearGrowthClass                ), pointer       :: linearGrowth_
    class           (haloMassFunctionClass            ), pointer       :: haloMassFunction_
    class           (haloEnvironmentClass             ), pointer       :: haloEnvironment_
    class           (unevolvedSubhaloMassFunctionClass), pointer       :: unevolvedSubhaloMassFunction_
    class           (darkMatterHaloScaleClass         ), pointer       :: darkMatterHaloScale_
    class           (cosmologicalMassVarianceClass    ), pointer       :: cosmologicalMassVariance_
    class           (darkMatterHaloBiasClass          ), pointer       :: darkMatterHaloBias_
    class           (transferFunctionClass            ), pointer       :: transferFunction_
    class           (outputTimesClass                 ), pointer       :: outputTimes_
    class           (darkMatterProfileScaleRadiusClass), pointer       :: darkMatterProfileScaleRadius_
    type            (varying_string                   )                :: outputGroup
    double precision                                                   :: haloMassMinimum              , haloMassMaximum
    integer                                                            :: pointsPerDecade

    call nodeClassHierarchyInitialize     (parameters)
    call Node_Components_Initialize       (parameters)
    call Node_Components_Thread_Initialize(parameters)
    !# <inputParameter>
    !#   <name>haloMassMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d10</defaultValue>
    !#   <description>The minimum mass at which to tabulate halo mass functions.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>haloMassMaximum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d15</defaultValue>
    !#   <description>The maximum mass at which to tabulate halo mass functions.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>pointsPerDecade</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>10</defaultValue>
    !#   <description>The number of points per decade of halo mass at which to tabulate halo mass functions.</description>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>    
    !# <inputParameter>
    !#   <name>outputGroup</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>var_str('.')</defaultValue>
    !#   <description>The HDF5 output group within which to write mass function data.</description>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>    
    !# <objectBuilder class="cosmologyParameters"          name="cosmologyParameters_"          source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"           name="cosmologyFunctions_"           source="parameters"/>
    !# <objectBuilder class="virialDensityContrast"        name="virialDensityContrast_"        source="parameters"/>
    !# <objectBuilder class="darkMatterProfile"            name="darkMatterProfile_"            source="parameters"/>
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
    self=taskHaloMassFunction(                               &
         &                    haloMassMinimum              , &
         &                    haloMassMaximum              , &
         &                    pointsPerDecade              , &
         &                    outputGroup                  , &
         &                    cosmologyParameters_         , &
         &                    cosmologyFunctions_          , &
         &                    virialDensityContrast_       , &
         &                    darkMatterProfile_           , &
         &                    criticalOverdensity_         , &
         &                    linearGrowth_                , &
         &                    haloMassFunction_            , &
         &                    haloEnvironment_             , &
         &                    unevolvedSubhaloMassFunction_, &
         &                    darkMatterHaloScale_         , &
         &                    darkMatterProfileScaleRadius_, &
         &                    cosmologicalMassVariance_    , &
         &                    darkMatterHaloBias_          , &
         &                    transferFunction_            , &
         &                    outputTimes_                   &
         &                   )
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"         />
    !# <objectDestructor name="cosmologyFunctions_"          />
    !# <objectDestructor name="virialDensityContrast_"       />
    !# <objectDestructor name="darkMatterProfile_"           />
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

  function haloMassFunctionConstructorInternal(                               &
       &                                       haloMassMinimum              , &
       &                                       haloMassMaximum              , &
       &                                       pointsPerDecade              , &
       &                                       outputGroup                  , &
       &                                       cosmologyParameters_         , &
       &                                       cosmologyFunctions_          , &
       &                                       virialDensityContrast_       , &
       &                                       darkMatterProfile_           , &
       &                                       criticalOverdensity_         , &
       &                                       linearGrowth_                , &
       &                                       haloMassFunction_            , &
       &                                       haloEnvironment_             , &
       &                                       unevolvedSubhaloMassFunction_, &
       &                                       darkMatterHaloScale_         , &
       &                                       darkMatterProfileScaleRadius_, &
       &                                       cosmologicalMassVariance_    , &
       &                                       darkMatterHaloBias_          , &
       &                                       transferFunction_            , &
       &                                       outputTimes_                   &
       &                                      ) result(self)
    !% Constructor for the {\normalfont \ttfamily haloMassFunction} task class which takes a parameter set as input.
    implicit none
    type            (taskHaloMassFunction             )                        :: self
    class           (cosmologyParametersClass         ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass          ), intent(in   ), target :: cosmologyFunctions_
    class           (virialDensityContrastClass       ), intent(in   ), target :: virialDensityContrast_
    class           (darkMatterProfileClass           ), intent(in   ), target :: darkMatterProfile_
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
    double precision                                   , intent(in   )         :: haloMassMinimum              , haloMassMaximum
    integer                                            , intent(in   )         :: pointsPerDecade
    !# <constructorAssign variables="haloMassMinimum,haloMassMaximum,pointsPerDecade,outputGroup,*cosmologyParameters_,*cosmologyFunctions_,*virialDensityContrast_,*darkMatterProfile_,*criticalOverdensity_,*linearGrowth_,*haloMassFunction_,*haloEnvironment_,*unevolvedSubhaloMassFunction_,*darkMatterHaloScale_, *darkMatterProfileScaleRadius_, *cosmologicalMassVariance_,*darkMatterHaloBias_,*transferFunction_, *outputTimes_"/>
    
    return
  end function haloMassFunctionConstructorInternal
  
  subroutine haloMassFunctionDestructor(self)
    !% Destructor for the {\normalfont \ttfamily haloMassFunction} task class.
    use Node_Components, only : Node_Components_Uninitialize, Node_Components_Thread_Uninitialize
    implicit none
    type(taskHaloMassFunction), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"         />
    !# <objectDestructor name="self%cosmologyFunctions_"          />
    !# <objectDestructor name="self%virialDensityContrast_"       />
    !# <objectDestructor name="self%darkMatterProfile_"           />
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
    call Node_Components_Uninitialize       ()
    call Node_Components_Thread_Uninitialize()
    return
  end subroutine haloMassFunctionDestructor

  subroutine haloMassFunctionPerform(self,status)
    !% Compute and output the halo mass function.
    use, intrinsic :: ISO_C_Binding
    use               Galacticus_Error
    use               Galacticus_Display    
    use               Galacticus_Nodes
    use               Galacticus_Calculations_Resets
    use               Galacticus_HDF5
    use               Numerical_Constants_Astronomical
    use               Memory_Management
    use               Numerical_Ranges
    use               Numerical_Integration
    use               Dark_Matter_Profile_Scales
    use               IO_HDF5
    use               String_Handling
    use               FGSL                            , only : fgsl_function, fgsl_integration_workspace, FGSL_Integ_Gauss15
    implicit none
    class           (taskHaloMassFunction          ), intent(inout)                 :: self
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
         &                                                                             outputVirialDensityContrast                          , massHalo
    ! The upper limit to halo mass used when computing cumulative mass functions.
    double precision                                , parameter                     :: haloMassEffectiveInfinity                     =1.0d16
    ! Largest subhalo mass (in units of host mass) for which we expect significant unevolved subhalo mass function.
    double precision                                , parameter                     :: subhaloMassMaximum                            =1.0d02
    class           (nodeComponentBasic            ), pointer                       :: basic
    class           (nodeComponentDarkMatterProfile), pointer                       :: darkMatterProfileHalo
    type            (mergerTree                    ), target                        :: tree
    type            (fgsl_function                 )                                :: integrandFunction
    type            (fgsl_integration_workspace    )                                :: integrationWorkspace
    integer         (c_size_t                      )                                :: iOutput                                              , outputCount                  , &
         &                                                                             iMass                                                , massCount
    double precision                                                                :: massHaloBinMinimum                                   , massHaloBinMaximum           , &
         &                                                                             massHaloLogarithmicInterval                          , massHalfMode
    type            (hdf5Object                    )                                :: outputsGroup                                         , outputGroup                  , &
         &                                                                             containerGroup                                       , powerSpectrumGroup           , &
         &                                                                             cosmologyGroup                                       , dataset
    integer                                                                         :: statusHalfModeMass
    type            (varying_string                )                                :: groupName                                            , commentText
   
    call Galacticus_Display_Indent('Begin task: halo mass function')
    ! Get the requested output redshifts.
    outputCount=self%outputTimes_%count()
    call allocateArray(outputTimes                                   ,[          outputCount])
    call allocateArray(outputRedshifts                               ,[          outputCount])
    call allocateArray(outputExpansionFactors                        ,[          outputCount])
    call allocateArray(outputGrowthFactors                           ,[          outputCount])
    call allocateArray(outputCriticalOverdensities                   ,[          outputCount])
    call allocateArray(outputVirialDensityContrast                   ,[          outputCount])
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
    ! Iterate over all output times.
    do iOutput=1,outputCount
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
          ! Compute halo properties.
          massFunctionDifferential                      (iMass,iOutput)=+self%haloMassFunction_            %differential                   (outputTimes(iOutput),massHalo          (iMass        )                          ,node=tree%baseNode)
          massFunctionCumulative                        (iMass,iOutput)=+self%haloMassFunction_            %integrated                     (outputTimes(iOutput),massHalo          (iMass        ),haloMassEffectiveInfinity,node=tree%baseNode)
          massFunctionMassFraction                      (iMass,iOutput)=+self%haloMassFunction_            %massFraction                   (outputTimes(iOutput),massHalo          (iMass        ),haloMassEffectiveInfinity,node=tree%baseNode)
          massFunctionDifferentialLogarithmicBinAveraged(iMass,iOutput)=+self%haloMassFunction_            %integrated                     (outputTimes(iOutput),massHaloBinMinimum               ,massHaloBinMaximum       ,node=tree%baseNode)    &
               &                                                        /massHaloLogarithmicInterval
          densityFieldRootVariance                      (iMass,iOutput)=+self%cosmologicalMassVariance_    %rootVariance                   (                     massHalo          (iMass        )                                             )
          densityFieldRootVarianceGradientLogarithmic   (iMass,iOutput)=+self%cosmologicalMassVariance_    %rootVarianceLogarithmicGradient(                     massHalo          (iMass        )                                             )
          peakHeightMassFunction                        (iMass,iOutput)=+massHalo                                                                                                  (iMass        )                                              **2 &
               &                                                        *massFunctionDifferential                                                                                  (iMass,iOutput)                                                  &
               &                                                        /self%cosmologyParameters_         %densityCritical                ()                                                                                                       &
               &                                                        /self%cosmologyParameters_         %OmegaMatter                    ()                                                                                                       &
               &                                                        /abs(                                                                                                                                                                       &
               &                                                             self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient(                                                                                                        &
               &                                                                                                                                                 massHalo          (iMass        )                                                  &
               &                                                                                                                           )                                                                                                        &
               &                                                            )
          biasHalo                                      (iMass,iOutput)=self%darkMatterHaloBias_           %bias                           (                                                                                      tree%baseNode)
          velocityVirial                                (iMass,iOutput)=self%darkMatterHaloScale_          %virialVelocity                 (                                                                                      tree%baseNode)
          temperatureVirial                             (iMass,iOutput)=self%darkMatterHaloScale_          %virialTemperature              (                                                                                      tree%baseNode)
          radiusVirial                                  (iMass,iOutput)=self%darkMatterHaloScale_          %virialRadius                   (                                                                                      tree%baseNode)
          velocityMaximum                               (iMass,iOutput)=self%darkMatterProfile_            %circularVelocityMaximum        (                                                                                      tree%baseNode)
          darkMatterProfileRadiusScale                  (iMass,iOutput)=darkMatterProfileHalo              %scale                          (                                                                                                   )
          ! Integrate the unevolved subhalo mass function over the halo mass function to get the total subhalo mass function.
          massFunctionCumulativeSubhalo                 (iMass,iOutput)=Integrate(                                                              &
               &                                                                                    log(massHalo(1)/subhaloMassMaximum       ), &
               &                                                                                    log(            haloMassEffectiveInfinity), &
               &                                                                                    subhaloMassFunctionIntegrand              , &
               &                                                                                    integrandFunction                         , &
               &                                                                                    integrationWorkspace                      , &
               &                                                                  toleranceAbsolute=0.0d+0                                    , &
               &                                                                  toleranceRelative=1.0d-3                                    , &
               &                                                                  integrationRule  =FGSL_Integ_Gauss15                          &
               &                                                                 )
          call Integrate_Done(integrandFunction,integrationWorkspace)
       end do
       peakHeight                         (:,iOutput)=+outputCriticalOverdensities(  iOutput) &
            &                                         /densityFieldRootVariance   (:,iOutput)
       massFunctionDifferentialLogarithmic(:,iOutput)=+massFunctionDifferential   (:,iOutput) &
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
       call outputGroup%writeAttribute(outputTimes                                   (  iOutput),'outputTime'                                                                                                                          )
       call outputGroup%writeAttribute(outputRedshifts                               (  iOutput),'outputRedshift'                                                                                                                      )
       call outputGroup%writeAttribute(outputExpansionFactors                        (  iOutput),'outputExpansionFactor'                                                                                                               )
       call outputGroup%writeAttribute(outputGrowthFactors                           (  iOutput),'growthFactor'                                                                                                                        )
       call outputGroup%writeAttribute(outputCriticalOverdensities                   (  iOutput),'criticalOverdensity'                                                                                                                 )
       call outputGroup%writeAttribute(outputVirialDensityContrast                   (  iOutput),'virialDensityContrast'                                                                                                               )    
       call outputGroup%writeAttribute(outputCharacteristicMass                      (  iOutput),'massHaloCharacteristic'                                                                                                              )
       call outputGroup%writeDataset  (massHalo                                      (:        ),'haloMass'                      ,'The mass of the halo.'                                                      ,datasetReturned=dataset)
       call dataset    %writeAttribute(massSolar                                                ,'unitsInSI'                                                                                                                           )
       call dataset    %close         (                                                                                                                                                                                                )
       call outputGroup%writeDataset  (massFunctionDifferential                      (:,iOutput),'haloMassFunctionM'             ,'The halo mass function (per unit halo mass).'                               ,datasetReturned=dataset)
       call dataset    %writeAttribute(1.0d0/megaParsec**3/massSolar                            ,'unitsInSI'                                                                                                                           )
       call dataset    %close         (                                                                                                                                                                                                )
       call outputGroup%writeDataset  (massFunctionDifferentialLogarithmic           (:,iOutput),'haloMassFunctionLnM'           ,'The halo mass function (per logarithmic halo mass).'                        ,datasetReturned=dataset)
       call dataset    %writeAttribute(1.0d0/megaParsec**3                                      ,'unitsInSI'                                                                                                                           )
       call dataset    %close         (                                                                                                                                                                                                )
       call outputGroup%writeDataset  (massFunctionDifferentialLogarithmicBinAveraged(:,iOutput),'haloMassFunctionLnMBinAveraged','The halo mass function (per logarithmic halo mass averaged across the bin).',datasetReturned=dataset)
       call dataset    %writeAttribute(1.0d0/megaParsec**3                                      ,'unitsInSI'                                                                                                                           )
       call dataset    %close         (                                                                                                                                                                                                )
       call outputGroup%writeDataset  (massFunctionCumulative                        (:,iOutput),'haloMassFunctionCumulative'    ,'The halo cumulative mass function.'                                         ,datasetReturned=dataset)
       call dataset    %writeAttribute(1.0d0/megaParsec**3                                      ,'unitsInSI'                                                                                                                           )
       call dataset    %close         (                                                                                                                                                                                                )
       call outputGroup%writeDataset  (massFunctionCumulativeSubhalo                 (:,iOutput),'subhaloMassFunctionCumulative' ,'The subhalo cumulative mass function.'                                      ,datasetReturned=dataset)
       call dataset    %writeAttribute(1.0d0/megaParsec**3                                      ,'unitsInSI'                                                                                                                           )
       call dataset    %close         (                                                                                                                                                                                                )
       call outputGroup%writeDataset  (velocityVirial                                (:,iOutput),'haloVirialVelocity'            ,'The virial velocity of halos.'                                              ,datasetReturned=dataset)
       call dataset    %writeAttribute(kilo                                                     ,'unitsInSI'                                                                                                                           )
       call dataset    %close         (                                                                                                                                                                                                )
       call outputGroup%writeDataset  (temperatureVirial                             (:,iOutput),'haloVirialTemperature'         ,'The virial temperature of halos.'                                           ,datasetReturned=dataset)
       call dataset    %writeAttribute(1.0d0                                                    ,'unitsInSI'                                                                                                                           )
       call dataset    %close         (                                                                                                                                                                                                )
       call outputGroup%writeDataset  (radiusVirial                                  (:,iOutput),'haloVirialRadius'              ,'The virial radius of halos.'                                                ,datasetReturned=dataset)
       call dataset    %writeAttribute(megaParsec                                               ,'unitsInSI'                                                                                                                           )
       call dataset    %close         (                                                                                                                                                                                                )
       call outputGroup%writeDataset  (darkMatterProfileRadiusScale                  (:,iOutput),'haloScaleRadius'               ,'The scale radius of halos.'                                                 ,datasetReturned=dataset)
       call dataset    %writeAttribute(megaParsec                                               ,'unitsInSI'                                                                                                                           )
       call dataset    %close         (                                                                                                                                                                                                )
       call outputGroup%writeDataset  (velocityMaximum                               (:,iOutput),'haloVelocityMaximum'           ,'The maximum circular velocity of halos.'                                    ,datasetReturned=dataset)
       call dataset    %writeAttribute(kilo                                                     ,'unitsInSI'                                                                                                                           )
       call dataset    %close         (                                                                                                                                                                                                )
       call outputGroup%writeDataset  (biasHalo                                      (:,iOutput),'haloBias'                      ,'The large scale linear bias of halos.'                                                              )
       call outputGroup%writeDataset  (densityFieldRootVariance                      (:,iOutput),'haloSigma'                     ,'The mass fluctuation on the scale of the halo.'                                                     )
       call outputGroup%writeDataset  (densityFieldRootVarianceGradientLogarithmic   (:,iOutput),'haloAlpha'                     ,'dlog(σ)/dlog(m).'                                                                                   )
       call outputGroup%writeDataset  (peakHeight                                    (:,iOutput),'haloPeakHeightNu'              ,'The peak height, ν, of the halo.'                                                                   )
       call outputGroup%writeDataset  (massFunctionMassFraction                      (:,iOutput),'haloMassFractionCumulative'    ,'The halo cumulative mass fraction.'                                                                 )
       call outputGroup%writeDataset  (peakHeightMassFunction                        (:,iOutput),'haloMassFunctionNuFNu'         ,'The halo mass fraction function as a function of the ν parameter, νF(ν).'                           )
       call outputGroup%close         (                                                                                                                                                                                                )
    end do
    call outputsGroup%close()
    if (containerGroup%isOpen()) call containerGroup%close()
    if (present(status)) status=errorStatusSuccess
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
