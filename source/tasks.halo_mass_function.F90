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

  use :: Cosmological_Density_Field               , only : cosmologicalMassVarianceClass          , criticalOverdensityClass         , haloEnvironmentClass
  use :: Cosmology_Functions                      , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters                     , only : cosmologyParametersClass
  use :: Dark_Matter_Halo_Biases                  , only : darkMatterHaloBiasClass
  use :: Dark_Matter_Halo_Scales                  , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Halo_Mass_Accretion_Histories, only : darkMatterHaloMassAccretionHistoryClass
  use :: Dark_Matter_Profile_Scales               , only : darkMatterProfileScaleRadius           , darkMatterProfileScaleRadiusClass
  use :: Dark_Matter_Profiles_Shape               , only : darkMatterProfileShapeClass
  use :: Dark_Matter_Profiles_DMO                 , only : darkMatterProfileDMOClass
  use :: Halo_Mass_Functions                      , only : haloMassFunctionClass
  use :: Linear_Growth                            , only : linearGrowthClass
  use :: Output_Times                             , only : outputTimesClass
  use :: Numerical_Random_Numbers                 , only : randomNumberGeneratorClass
  use :: Transfer_Functions                       , only : transferFunctionClass
  use :: Unevolved_Subhalo_Mass_Functions         , only : unevolvedSubhaloMassFunctionClass
  use :: Virial_Density_Contrast                  , only : virialDensityContrastClass

  type :: virialDensityContrastList
     !!{
     Type used to store a list of virial density contrasts.
     !!}
     type (varying_string            )          :: label
     class(virialDensityContrastClass), pointer :: virialDensityContrast_ => null()
  end type virialDensityContrastList
  
  !![
  <task name="taskHaloMassFunction">
    <description>A task which computes and outputs the halo mass function and related quantities.</description>
    <descriptorSpecial>descriptorSpecial</descriptorSpecial>
  </task>
  !!]
  type, extends(taskClass) :: taskHaloMassFunction
     !!{
     Implementation of a task which computes and outputs the halo mass function and related quantities.
     !!}
     private
     class           (cosmologyParametersClass               ), pointer                   :: cosmologyParameters_                => null()
     class           (cosmologyFunctionsClass                ), pointer                   :: cosmologyFunctions_                 => null()
     class           (virialDensityContrastClass             ), pointer                   :: virialDensityContrast_              => null()
     class           (criticalOverdensityClass               ), pointer                   :: criticalOverdensity_                => null()
     class           (linearGrowthClass                      ), pointer                   :: linearGrowth_                       => null()
     class           (haloMassFunctionClass                  ), pointer                   :: haloMassFunction_                   => null()
     class           (haloEnvironmentClass                   ), pointer                   :: haloEnvironment_                    => null()
     class           (unevolvedSubhaloMassFunctionClass      ), pointer                   :: unevolvedSubhaloMassFunction_       => null()
     class           (darkMatterHaloScaleClass               ), pointer                   :: darkMatterHaloScale_                => null()
     class           (cosmologicalMassVarianceClass          ), pointer                   :: cosmologicalMassVariance_           => null()
     class           (darkMatterHaloBiasClass                ), pointer                   :: darkMatterHaloBias_                 => null()
     class           (transferFunctionClass                  ), pointer                   :: transferFunction_                   => null(), transferFunctionReference => null(), &
          &                                                                                  transferFunctionRelative            => null()
     class           (outputTimesClass                       ), pointer                   :: outputTimes_                        => null()
     class           (darkMatterProfileScaleRadiusClass      ), pointer                   :: darkMatterProfileScaleRadius_       => null()
     class           (darkMatterProfileShapeClass            ), pointer                   :: darkMatterProfileShape_             => null()
     class           (darkMatterHaloMassAccretionHistoryClass), pointer                   :: darkMatterHaloMassAccretionHistory_ => null()
     class           (randomNumberGeneratorClass             ), pointer                   :: randomNumberGenerator_              => null()
     double precision                                                                     :: haloMassMinimum                              , haloMassMaximum                     , &
          &                                                                                  pointsPerDecade
     type            (varying_string                         )                            :: outputGroup
     logical                                                                              :: includeUnevolvedSubhaloMassFunction          , includeMassAccretionRate            , &
          &                                                                                  massesRelativeToHalfModeMass                 , nodeComponentsInitialized =  .false., &
          &                                                                                  errorsAreFatal
     double precision                                         , allocatable, dimension(:) :: fractionModeMasses
     type            (virialDensityContrastList              ), allocatable, dimension(:) :: virialDensityContrasts
     ! Pointer to the parameters for this task.
     type            (inputParameters                        ), pointer                   :: parameters                          => null()
   contains
     !![
     <methods>
       <method method="descriptorSpecial" description="Handle adding special parameters to the descriptor."/>
     </methods>
     !!]
     final     ::                      haloMassFunctionDestructor
     procedure :: perform           => haloMassFunctionPerform
     procedure :: descriptorSpecial => haloMassFunctionDescriptorSpecial
  end type taskHaloMassFunction

  interface taskHaloMassFunction
     !!{
     Constructors for the \refClass{taskHaloMassFunction} task.
     !!}
     module procedure haloMassFunctionConstructorParameters
     module procedure haloMassFunctionConstructorInternal
  end interface taskHaloMassFunction

contains

  function haloMassFunctionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskHaloMassFunction} task class which takes a parameter set as input.
    !!}
    use :: Galacticus_Nodes, only : nodeClassHierarchyInitialize, treeNode
    use :: Input_Parameters, only : inputParameter              , inputParameters
    use :: Node_Components , only : Node_Components_Initialize
    implicit none
    type            (taskHaloMassFunction                   )                              :: self
    type            (inputParameters                        ), intent(inout), target       :: parameters
    class           (cosmologyParametersClass               ), pointer                     :: cosmologyParameters_
    class           (cosmologyFunctionsClass                ), pointer                     :: cosmologyFunctions_
    class           (virialDensityContrastClass             ), pointer                     :: virialDensityContrast_
    class           (criticalOverdensityClass               ), pointer                     :: criticalOverdensity_
    class           (linearGrowthClass                      ), pointer                     :: linearGrowth_
    class           (haloMassFunctionClass                  ), pointer                     :: haloMassFunction_
    class           (haloEnvironmentClass                   ), pointer                     :: haloEnvironment_
    class           (unevolvedSubhaloMassFunctionClass      ), pointer                     :: unevolvedSubhaloMassFunction_
    class           (darkMatterHaloScaleClass               ), pointer                     :: darkMatterHaloScale_
    class           (darkMatterProfileShapeClass            ), pointer                     :: darkMatterProfileShape_
    class           (cosmologicalMassVarianceClass          ), pointer                     :: cosmologicalMassVariance_
    class           (darkMatterHaloBiasClass                ), pointer                     :: darkMatterHaloBias_
    class           (transferFunctionClass                  ), pointer                     :: transferFunction_                  , transferFunctionReference, &
         &                                                                                    transferFunctionRelative
    class           (outputTimesClass                       ), pointer                     :: outputTimes_
    class           (darkMatterProfileScaleRadiusClass      ), pointer                     :: darkMatterProfileScaleRadius_
    class           (darkMatterHaloMassAccretionHistoryClass), pointer                     :: darkMatterHaloMassAccretionHistory_
    class           (randomNumberGeneratorClass             ), pointer                     :: randomNumberGenerator_
    type            (inputParameters                        ), pointer                     :: parametersRoot
    type            (inputParameters                        ),                             :: parametersMassDefinitions
    type            (virialDensityContrastList              ), allocatable  , dimension(:) :: virialDensityContrasts
    type            (varying_string                         ), allocatable  , dimension(:) :: labels
    double precision                                         , allocatable  , dimension(:) :: fractionModeMasses
    type            (varying_string                         )                              :: outputGroup
    double precision                                                                       :: haloMassMinimum                    , haloMassMaximum           , &
         &                                                                                    pointsPerDecade
    logical                                                                                :: includeUnevolvedSubhaloMassFunction, includeMassAccretionRate  , &
          &                                                                                   massesRelativeToHalfModeMass       , errorsAreFatal
    integer                                                                                :: i
    
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
    self%nodeComponentsInitialized=.true.
    !![
    <inputParameter>
      <name>haloMassMinimum</name>
      <defaultValue>1.0d10</defaultValue>
      <description>The minimum mass at which to tabulate halo mass functions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>haloMassMaximum</name>
      <defaultValue>1.0d15</defaultValue>
      <description>The maximum mass at which to tabulate halo mass functions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>pointsPerDecade</name>
      <defaultValue>10.0d0</defaultValue>
      <description>The number of points per decade of halo mass at which to tabulate halo mass functions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>outputGroup</name>
      <defaultValue>var_str('.')</defaultValue>
      <description>The HDF5 output group within which to write mass function data.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>includeUnevolvedSubhaloMassFunction</name>
      <defaultValue>.false.</defaultValue>
      <description>If true then also compute and output the unevolved subhalo mass function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>includeMassAccretionRate</name>
      <defaultValue>.true.</defaultValue>
      <description>If true then also compute and output the mass accretion rate of the halos.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massesRelativeToHalfModeMass</name>
      <defaultValue>.false.</defaultValue>
      <description>If true then masses are interpreted (and output) relative to the half-mode mass. (If the half-mode mass is undefined an error will occur.) If false, masses are absolute.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>errorsAreFatal</name>
      <defaultValue>.true.</defaultValue>
      <description>If true then errors in evaluating the halo mass function are considered to be fatal.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder    class="cosmologyParameters"                name="cosmologyParameters_"                source="parameters"                                          />
    <objectBuilder    class="cosmologyFunctions"                 name="cosmologyFunctions_"                 source="parameters"                                          />
    <objectBuilder    class="virialDensityContrast"              name="virialDensityContrast_"              source="parameters"                                          />
    <objectBuilder    class="criticalOverdensity"                name="criticalOverdensity_"                source="parameters"                                          />
    <objectBuilder    class="linearGrowth"                       name="linearGrowth_"                       source="parameters"                                          />
    <objectBuilder    class="haloMassFunction"                   name="haloMassFunction_"                   source="parameters"                                          />
    <objectBuilder    class="haloEnvironment"                    name="haloEnvironment_"                    source="parameters"                                          />
    <objectBuilder    class="unevolvedSubhaloMassFunction"       name="unevolvedSubhaloMassFunction_"       source="parameters"                                          />
    <objectBuilder    class="darkMatterHaloScale"                name="darkMatterHaloScale_"                source="parameters"                                          />
    <objectBuilder    class="cosmologicalMassVariance"           name="cosmologicalMassVariance_"           source="parameters"                                          />
    <objectBuilder    class="darkMatterHaloBias"                 name="darkMatterHaloBias_"                 source="parameters"                                          />
    <objectBuilder    class="transferFunction"                   name="transferFunction_"                   source="parameters"                                          />
    <objectBuilder    class="outputTimes"                        name="outputTimes_"                        source="parameters"                                          />
    <objectBuilder    class="darkMatterProfileScaleRadius"       name="darkMatterProfileScaleRadius_"       source="parameters"                                          />
    <objectBuilder    class="darkMatterProfileShape"             name="darkMatterProfileShape_"             source="parameters"                                          />
    <objectBuilder    class="darkMatterHaloMassAccretionHistory" name="darkMatterHaloMassAccretionHistory_" source="parameters"                                          />
    <objectBuilder    class="randomNumberGenerator"              name="randomNumberGenerator_"              source="parameters"                                          />
    !!]
    if (parameters%isPresent('transferFunctionReference')) then
       !![
       <objectBuilder class="transferFunction"                   name="transferFunctionReference"           source="parameters" parameterName="transferFunctionReference"/>
       !!]
    end if
    if (parameters%isPresent('transferFunctionRelative' )) then
       !![
       <objectBuilder class="transferFunction"                   name="transferFunctionRelative"            source="parameters" parameterName="transferFunctionRelative" />
       !!]
    end if
    if (parameters%isPresent('massDefinitions',requireValue=.false.)) then
       parametersMassDefinitions=parameters%subParameters('massDefinitions',requireValue=.false.)
       if (parametersMassDefinitions%copiesCount('virialDensityContrast') /= parametersMassDefinitions%count('labels')) &
            & call Error_Report('number of labels must match number of virial density contrasts'//{introspection:location})
       allocate(virialDensityContrasts(parametersMassDefinitions%copiesCount('virialDensityContrast')))
       allocate(labels                (parametersMassDefinitions%copiesCount('virialDensityContrast')))
       !![
       <inputParameter>
         <name>labels</name>
         <description>Labels for virial density contrast mass definitions.</description>
         <source>parametersMassDefinitions</source>
       </inputParameter>
       !!]
       do i=1,parametersMassDefinitions%copiesCount('virialDensityContrast')
          !![
          <objectBuilder class="virialDensityContrast" name="virialDensityContrasts(i)%virialDensityContrast_" source="parametersMassDefinitions" copy="i" />
          !!]
          virialDensityContrasts(i)%label=labels(i)
       end do
    else
       allocate(virialDensityContrasts(0))
    end if
    allocate(fractionModeMasses(parameters%count('fractionModeMasses',zeroIfNotPresent=.true.)))
    if (size(fractionModeMasses) > 0) then
       !![
       <inputParameter>
         <name>fractionModeMasses</name>
         <description>List of suppression fractions at which to compute the fractional mode mass.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
    end if    
    !![
    <conditionalCall>
      <call>
	self=taskHaloMassFunction(                                     &amp;
         &amp;                    haloMassMinimum                    , &amp;
         &amp;                    haloMassMaximum                    , &amp;
         &amp;                    pointsPerDecade                    , &amp;
         &amp;                    outputGroup                        , &amp;
         &amp;                    includeUnevolvedSubhaloMassFunction, &amp;
         &amp;                    includeMassAccretionRate           , &amp;
         &amp;                    massesRelativeToHalfModeMass       , &amp;
	 &amp;                    errorsAreFatal                     , &amp;
	 &amp;                    fractionModeMasses                 , &amp;
         &amp;                    cosmologyParameters_               , &amp;
         &amp;                    cosmologyFunctions_                , &amp;
         &amp;                    virialDensityContrast_             , &amp;
         &amp;                    criticalOverdensity_               , &amp;
         &amp;                    linearGrowth_                      , &amp;
         &amp;                    haloMassFunction_                  , &amp;
         &amp;                    haloEnvironment_                   , &amp;
         &amp;                    unevolvedSubhaloMassFunction_      , &amp;
         &amp;                    darkMatterHaloScale_               , &amp;
         &amp;                    darkMatterProfileScaleRadius_      , &amp;
         &amp;                    darkMatterProfileShape_            , &amp;
         &amp;                    darkMatterHaloMassAccretionHistory_, &amp;
         &amp;                    cosmologicalMassVariance_          , &amp;
         &amp;                    darkMatterHaloBias_                , &amp;
         &amp;                    transferFunction_                  , &amp;
         &amp;                    outputTimes_                       , &amp;
         &amp;                    randomNumberGenerator_             , &amp;
         &amp;                    virialDensityContrasts             , &amp;
         &amp;                    parametersRoot                       &amp;
         &amp;                    {conditions}                         &amp;
         &amp;                   )
     </call>
     <argument name="transferFunctionReference" value="transferFunctionReference" parameterPresent="parameters"/>
     <argument name="transferFunctionRelative"  value="transferFunctionRelative"  parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"               />
    <objectDestructor name="cosmologyFunctions_"                />
    <objectDestructor name="virialDensityContrast_"             />
    <objectDestructor name="criticalOverdensity_"               />
    <objectDestructor name="linearGrowth_"                      />
    <objectDestructor name="haloMassFunction_"                  />
    <objectDestructor name="haloEnvironment_"                   />
    <objectDestructor name="unevolvedSubhaloMassFunction_"      />
    <objectDestructor name="darkMatterHaloScale_"               />
    <objectDestructor name="cosmologicalMassVariance_"          />
    <objectDestructor name="darkMatterHaloBias_"                />
    <objectDestructor name="transferFunction_"                  />
    <objectDestructor name="outputTimes_"                       />
    <objectDestructor name="darkMatterProfileScaleRadius_"      />
    <objectDestructor name="darkMatterProfileShape_"            />
    <objectDestructor name="darkMatterHaloMassAccretionHistory_"/>
    <objectDestructor name="randomNumberGenerator_"             />
    !!]
    if (parameters%isPresent('transferFunctionReference')) then
       !![
       <objectDestructor name="transferFunctionReference"/>
       !!]
    end if
    if (parameters%isPresent('transferFunctionRelative' )) then
       !![
       <objectDestructor name="transferFunctionRelative" />
       !!]
    end if
    if (size(virialDensityContrasts) > 0) then
       do i=1,parametersMassDefinitions%copiesCount('virialDensityContrast')
          !![
          <objectDestructor name="virialDensityContrasts(i)%virialDensityContrast_"/>
          !!]
       end do
    end if
    return
  end function haloMassFunctionConstructorParameters

  function haloMassFunctionConstructorInternal(                                     &
       &                                       haloMassMinimum                    , &
       &                                       haloMassMaximum                    , &
       &                                       pointsPerDecade                    , &
       &                                       outputGroup                        , &
       &                                       includeUnevolvedSubhaloMassFunction, &
       &                                       includeMassAccretionRate           , &
       &                                       massesRelativeToHalfModeMass       , &
       &                                       errorsAreFatal                     , &
       &                                       fractionModeMasses                 , &
       &                                       cosmologyParameters_               , &
       &                                       cosmologyFunctions_                , &
       &                                       virialDensityContrast_             , &
       &                                       criticalOverdensity_               , &
       &                                       linearGrowth_                      , &
       &                                       haloMassFunction_                  , &
       &                                       haloEnvironment_                   , &
       &                                       unevolvedSubhaloMassFunction_      , &
       &                                       darkMatterHaloScale_               , &
       &                                       darkMatterProfileScaleRadius_      , &
       &                                       darkMatterProfileShape_            , &
       &                                       darkMatterHaloMassAccretionHistory_, &
       &                                       cosmologicalMassVariance_          , &
       &                                       darkMatterHaloBias_                , &
       &                                       transferFunction_                  , &
       &                                       outputTimes_                       , &
       &                                       randomNumberGenerator_             , &
       &                                       virialDensityContrasts             , &
       &                                       parameters                         , &
       &                                       transferFunctionReference          , &
       &                                       transferFunctionRelative             &
       &                                      ) result(self)
    !!{
    Constructor for the \refClass{taskHaloMassFunction} task class which takes a parameter set as input.
    !!}
    implicit none
    type            (taskHaloMassFunction                   )                                        :: self
    class           (cosmologyParametersClass               ), intent(in   ), target                 :: cosmologyParameters_
    class           (cosmologyFunctionsClass                ), intent(in   ), target                 :: cosmologyFunctions_
    class           (virialDensityContrastClass             ), intent(in   ), target                 :: virialDensityContrast_
    class           (criticalOverdensityClass               ), intent(in   ), target                 :: criticalOverdensity_
    class           (linearGrowthClass                      ), intent(in   ), target                 :: linearGrowth_
    class           (haloMassFunctionClass                  ), intent(in   ), target                 :: haloMassFunction_
    class           (haloEnvironmentClass                   ), intent(in   ), target                 :: haloEnvironment_
    class           (unevolvedSubhaloMassFunctionClass      ), intent(in   ), target                 :: unevolvedSubhaloMassFunction_
    class           (darkMatterHaloScaleClass               ), intent(in   ), target                 :: darkMatterHaloScale_
    class           (darkMatterProfileScaleRadiusClass      ), intent(in   ), target                 :: darkMatterProfileScaleRadius_
    class           (darkMatterProfileShapeClass            ), intent(in   ), target                 :: darkMatterProfileShape_
    class           (darkMatterHaloMassAccretionHistoryClass), intent(in   ), target                 :: darkMatterHaloMassAccretionHistory_
    class           (cosmologicalMassVarianceClass          ), intent(in   ), target                 :: cosmologicalMassVariance_
    class           (darkMatterHaloBiasClass                ), intent(in   ), target                 :: darkMatterHaloBias_
    class           (transferFunctionClass                  ), intent(in   ), target                 :: transferFunction_
    class           (transferFunctionClass                  ), intent(in   ), target      , optional :: transferFunctionReference          , transferFunctionRelative
    class           (outputTimesClass                       ), intent(in   ), target                 :: outputTimes_
    class           (randomNumberGeneratorClass             ), intent(in   ), target                 :: randomNumberGenerator_
    type            (virialDensityContrastList              ), intent(in   ), dimension(:)           :: virialDensityContrasts
    type            (varying_string                         ), intent(in   )                         :: outputGroup
    double precision                                         , intent(in   )                         :: haloMassMinimum                    , haloMassMaximum           , &
         &                                                                                              pointsPerDecade
    logical                                                  , intent(in   )                         :: includeUnevolvedSubhaloMassFunction, includeMassAccretionRate  , &
         &                                                                                              massesRelativeToHalfModeMass       , errorsAreFatal
    double precision                                         , intent(in   ), dimension(:)           :: fractionModeMasses
    type            (inputParameters                        ), intent(in   ), target                 :: parameters
    integer                                                                                          :: i
    !![
    <constructorAssign variables="haloMassMinimum, haloMassMaximum, pointsPerDecade, outputGroup, includeUnevolvedSubhaloMassFunction, includeMassAccretionRate, massesRelativeToHalfModeMass, errorsAreFatal, fractionModeMasses, *cosmologyParameters_, *cosmologyFunctions_, *virialDensityContrast_, *criticalOverdensity_, *linearGrowth_, *haloMassFunction_, *haloEnvironment_, *unevolvedSubhaloMassFunction_, *darkMatterHaloScale_, *darkMatterProfileScaleRadius_, *darkMatterProfileShape_, *darkMatterHaloMassAccretionHistory_, *cosmologicalMassVariance_, *darkMatterHaloBias_, *transferFunction_, *transferFunctionReference, *transferFunctionRelative, *outputTimes_, *randomNumberGenerator_"/>
    !!]

    self%parameters  => parameters
    allocate(self%virialDensityContrasts(size(virialDensityContrasts)))
    do i=1,size(virialDensityContrasts)
       self%virialDensityContrasts(i)%label=virialDensityContrasts(i)%label
       !![
       <referenceAcquire isResult="yes" owner="self" target="virialDensityContrasts(i)%virialDensityContrast_" source="virialDensityContrasts(i)%virialDensityContrast_"/>
       !!]
    end do
    return
  end function haloMassFunctionConstructorInternal

  subroutine haloMassFunctionDestructor(self)
    !!{
    Destructor for the \refClass{taskHaloMassFunction} task class.
    !!}
    use :: Node_Components, only : Node_Components_Uninitialize
    implicit none
    type   (taskHaloMassFunction), intent(inout) :: self
    integer                                      :: i

    !![
    <objectDestructor name="self%cosmologyParameters_"               />
    <objectDestructor name="self%cosmologyFunctions_"                />
    <objectDestructor name="self%virialDensityContrast_"             />
    <objectDestructor name="self%criticalOverdensity_"               />
    <objectDestructor name="self%linearGrowth_"                      />
    <objectDestructor name="self%haloMassFunction_"                  />
    <objectDestructor name="self%haloEnvironment_"                   />
    <objectDestructor name="self%unevolvedSubhaloMassFunction_"      />
    <objectDestructor name="self%darkMatterHaloScale_"               />
    <objectDestructor name="self%darkMatterProfileScaleRadius_"      />
    <objectDestructor name="self%darkMatterProfileShape_"            />
    <objectDestructor name="self%darkMatterHaloMassAccretionHistory_"/>
    <objectDestructor name="self%cosmologicalMassVariance_"          />
    <objectDestructor name="self%darkMatterHaloBias_"                />
    <objectDestructor name="self%transferFunction_"                  />
    <objectDestructor name="self%outputTimes_"                       />
    <objectDestructor name="self%randomNumberGenerator_"             />
    !!]
    if (associated(self%transferFunctionReference)) then
       !![
       <objectDestructor name="self%transferFunctionReference"/>
       !!]
    end if
    if (associated(self%transferFunctionRelative )) then
       !![
       <objectDestructor name="self%transferFunctionRelative" />
       !!]
    end if
    if (allocated(self%virialDensityContrasts)) then
       do i=1,size(self%virialDensityContrasts)
          !![
	  <objectDestructor name="self%virialDensityContrasts(i)%virialDensityContrast_"/>
          !!]
       end do
    end if
    if (self%nodeComponentsInitialized) call Node_Components_Uninitialize()
    return
  end subroutine haloMassFunctionDestructor

  subroutine haloMassFunctionPerform(self,status)
    !!{
    Compute and output the halo mass function.
    !!}
    use            :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use            :: Display                             , only : displayIndent                      , displayUnindent
    use            :: Calculations_Resets                 , only : Calculations_Reset
    use            :: Error                               , only : errorStatusSuccess                 , Error_Report                       , Warn
    use            :: Output_HDF5                         , only : outputFile
    use            :: Galacticus_Nodes                    , only : mergerTree                         , nodeComponentBasic                 , nodeComponentDarkMatterProfile, treeNode
    use            :: Galactic_Structure_Options          , only : componentTypeDarkMatterOnly        , massTypeDark
    use            :: IO_HDF5                             , only : hdf5Object
    use, intrinsic :: ISO_C_Binding                       , only : c_size_t
    use            :: Mass_Distributions                  , only : massDistributionClass
    use            :: Node_Components                     , only : Node_Components_Thread_Initialize  , Node_Components_Thread_Uninitialize
    use            :: Numerical_Constants_Astronomical    , only : gigaYear                           , massSolar                          , megaParsec
    use            :: Numerical_Constants_Math            , only : Pi
    use            :: Numerical_Constants_Prefixes        , only : kilo
    use            :: Numerical_Integration               , only : GSL_Integ_Gauss15                  , integrator
    use            :: Numerical_Ranges                    , only : Make_Range                         , rangeTypeLogarithmic
    use            :: String_Handling                     , only : String_Upper_Case_First            , operator(//)
    implicit none
    class           (taskHaloMassFunction                   ), intent(inout), target                 :: self
    integer                                                  , intent(  out), optional               :: status
    double precision                                         , allocatable  , dimension(:,:  )       :: massFunctionDifferentialLogarithmicBinAveraged         , biasHalo                     , &
         &                                                                                              massFunctionCumulative                                 , massFunctionDifferential     , &
         &                                                                                              massFunctionDifferentialLogarithmic                    , massFunctionMassFraction     , &
         &                                                                                              peakHeight                                             , densityFieldRootVariance     , &
         &                                                                                              radiusVirial                                           , temperatureVirial            , &
         &                                                                                              velocityVirial                                         , darkMatterProfileRadiusScale , &
         &                                                                                              velocityMaximum                                        , peakHeightMassFunction       , &
         &                                                                                              densityFieldRootVarianceGradientLogarithmic            , massFunctionCumulativeSubhalo, &
         &                                                                                              massAccretionRate
    double precision                                         , allocatable  , dimension(:    )       :: outputCharacteristicMass                               , outputCriticalOverdensities  , &
         &                                                                                              outputExpansionFactors                                 , outputGrowthFactors          , &
         &                                                                                              outputRedshifts                                        , outputTimes                  , &
         &                                                                                              outputVirialDensityContrast                            , outputTurnAroundRadius       , &
         &                                                                                              massHalo                                               , massHaloOutput               , &
         &                                                                                              massFractionMode                                       , slopeFractionMode            , &
         &                                                                                              wavenumberFractionMode
    double precision                                         , allocatable  , dimension(:,:,:)       :: massAlternate                                          , radiusAlternate
    integer                                                  , allocatable  , dimension(:    )       :: statusFractionMode
    ! The upper limit to halo mass used when computing cumulative mass functions.
    double precision                                         , parameter                             :: haloMassEffectiveInfinity                     =1.0d16
    ! Largest subhalo mass (in units of host mass) for which we expect significant unevolved subhalo mass function.
    double precision                                         , parameter                             :: subhaloMassMaximum                            =1.0d02
    class           (nodeComponentBasic                     ), pointer                               :: basic
    class           (nodeComponentDarkMatterProfile         ), pointer                               :: darkMatterProfileHalo
    class           (cosmologyParametersClass               ), pointer                               :: cosmologyParameters_                          => null()
    class           (cosmologyFunctionsClass                ), pointer                               :: cosmologyFunctions_                           => null()
    class           (criticalOverdensityClass               ), pointer                               :: criticalOverdensity_                          => null()
    class           (haloMassFunctionClass                  ), pointer                               :: haloMassFunction_                             => null()
    class           (haloEnvironmentClass                   ), pointer                               :: haloEnvironment_                              => null()
    class           (unevolvedSubhaloMassFunctionClass      ), pointer                               :: unevolvedSubhaloMassFunction_                 => null()
    class           (darkMatterHaloScaleClass               ), pointer                               :: darkMatterHaloScale_                          => null()
    class           (cosmologicalMassVarianceClass          ), pointer                               :: cosmologicalMassVariance_                     => null()
    class           (darkMatterHaloBiasClass                ), pointer                               :: darkMatterHaloBias_                           => null()
    class           (darkMatterProfileScaleRadiusClass      ), pointer                               :: darkMatterProfileScaleRadius_                 => null()
    class           (darkMatterProfileShapeClass            ), pointer                               :: darkMatterProfileShape_                       => null()
    class           (darkMatterHaloMassAccretionHistoryClass), pointer                               :: darkMatterHaloMassAccretionHistory_           => null()
    class           (virialDensityContrastClass             ), pointer                               :: virialDensityContrast_                        => null()
    !$omp threadprivate(haloEnvironment_,cosmologyFunctions_,cosmologyParameters_,cosmologicalMassVariance_,haloMassFunction_,darkMatterHaloScale_,unevolvedSubhaloMassFunction_,darkMatterHaloBias_,darkMatterProfileScaleRadius_,darkMatterProfileShape_,darkMatterHaloMassAccretionHistory_,virialDensityContrast_,criticalOverdensity_)
    class           (massDistributionClass                  ), pointer                        , save :: massDistribution_                             => null()
    !$omp threadprivate(massDistribution_)
    type            (virialDensityContrastList              ), allocatable   , dimension(:   )       :: virialDensityContrasts
    type            (mergerTree                             ), allocatable   , target         , save :: tree
    !$omp threadprivate(tree)
    type            (integrator                             ), allocatable                           :: integrator_
    type            (inputParameters                        ), allocatable                    , save :: parameters
    !$omp threadprivate(parameters)
    integer         (c_size_t                               )                                        :: iOutput                                                , outputCount                  , &
         &                                                                                              iMass                                                  , massCount                    , &
         &                                                                                              iAlternate
    double precision                                                                                 :: massHaloBinMinimum                                     , massHaloBinMaximum           , &
         &                                                                                              massHaloLogarithmicInterval                            , massHalfMode                 , &
         &                                                                                              massCritical                                           , massQuarterMode              , &
         &                                                                                              densityMean                                            , densityCritical              , &
         &                                                                                              massHaloMinimum                                        , massHaloMaximum              , &
         &                                                                                              massHalfModeReference                                  , slopeHalfMode                , &
         &                                                                                              slopeQuarterMode                                       , wavenumberHalfMode           , &
         &                                                                                              wavenumberQuarterMode
    type            (hdf5Object                             )                                        :: outputsGroup                                           , outputGroup                  , &
         &                                                                                              containerGroup                                         , powerSpectrumGroup           , &
         &                                                                                              cosmologyGroup                                         , dataset
    integer                                                                                          :: statusHalfModeMass                                     , statusQuarterModeMass        , &
         &                                                                                              statusHalfModeMassReference                            , statusIntegrated
    type            (varying_string                         )                                        :: groupName                                              , description
    character       (len=32                                 )                                        :: label
    logical                                                                                          :: scaleIsSettable                                        , shapeIsSettable              , &
         &                                                                                              warnedIntegratedFailure
    
    call displayIndent('Begin task: halo mass function')
    ! Call routines to perform initialization which must occur for all threads if run in parallel.
    call Node_Components_Thread_Initialize(self%parameters)
    ! Get the requested output redshifts.
    outputCount=self%outputTimes_%count()
    allocate(outputTimes                                   (outputCount))
    allocate(outputRedshifts                               (outputCount))
    allocate(outputExpansionFactors                        (outputCount))
    allocate(outputGrowthFactors                           (outputCount))
    allocate(outputCriticalOverdensities                   (outputCount))
    allocate(outputVirialDensityContrast                   (outputCount))
    allocate(outputTurnaroundRadius                        (outputCount))
    allocate(outputCharacteristicMass                      (outputCount))
    ! Compute number of tabulation points.
    massCount=int(log10(self%haloMassMaximum/self%haloMassMinimum)*self%pointsPerDecade)+1
    allocate(massHalo                                      (massCount            ))
    allocate(massHaloOutput                                (massCount            ))
    allocate(massFunctionDifferential                      (massCount,outputCount))
    allocate(massFunctionDifferentialLogarithmic           (massCount,outputCount))
    allocate(massFunctionDifferentialLogarithmicBinAveraged(massCount,outputCount))
    allocate(massFunctionCumulative                        (massCount,outputCount))
    allocate(massFunctionCumulativeSubhalo                 (massCount,outputCount))
    allocate(massFunctionMassFraction                      (massCount,outputCount))
    allocate(massAccretionRate                             (massCount,outputCount))
    allocate(biasHalo                                      (massCount,outputCount))
    allocate(densityFieldRootVariance                      (massCount,outputCount))
    allocate(densityFieldRootVarianceGradientLogarithmic   (massCount,outputCount))
    allocate(peakHeight                                    (massCount,outputCount))
    allocate(peakHeightMassFunction                        (massCount,outputCount))
    allocate(velocityVirial                                (massCount,outputCount))
    allocate(temperatureVirial                             (massCount,outputCount))
    allocate(radiusVirial                                  (massCount,outputCount))
    allocate(darkMatterProfileRadiusScale                  (massCount,outputCount))
    allocate(velocityMaximum                               (massCount,outputCount))
    allocate(massAlternate                                 (size(self%virialDensityContrasts),massCount,outputCount))
    allocate(radiusAlternate                               (size(self%virialDensityContrasts),massCount,outputCount))
    ! Compute output time properties.
    do iOutput=1,outputCount
       outputTimes                (iOutput)=self%outputTimes_%time                                 (                                                      iOutput )
       outputExpansionFactors     (iOutput)=self%cosmologyFunctions_   %expansionFactor            (                               outputTimes           (iOutput))
       outputRedshifts            (iOutput)=self%cosmologyFunctions_   %redshiftFromExpansionFactor(                               outputExpansionFactors(iOutput))
       outputGrowthFactors        (iOutput)=self%linearGrowth_         %value                      (                               outputTimes           (iOutput))
       outputCharacteristicMass   (iOutput)=self%criticalOverdensity_  %collapsingMass             (                               outputTimes           (iOutput))
       if (outputCharacteristicMass(iOutput) > 0.0d0) then
          massCritical=outputCharacteristicMass(iOutput)
       else
          massCritical=self%haloMassMinimum
       end if
       outputCriticalOverdensities(iOutput)=self%criticalOverdensity_  %value                      (mass=     massCritical   ,time=outputTimes           (iOutput))
       outputVirialDensityContrast(iOutput)=self%virialDensityContrast_%densityContrast            (mass=self%haloMassMinimum,time=outputTimes           (iOutput))
       outputTurnaroundRadius     (iOutput)=self%virialDensityContrast_%turnAroundOverVirialRadii  (mass=self%haloMassMinimum,time=outputTimes           (iOutput))
    end do
    ! Get half- and quarter-mode masses.
    massHalfMode   =self%transferFunction_%halfModeMass   (statusHalfModeMass   )
    massQuarterMode=self%transferFunction_%quarterModeMass(statusQuarterModeMass)
    if (size(self%fractionModeMasses) > 0) then
       allocate(massFractionMode      (size(self%fractionModeMasses)))
       allocate(statusFractionMode    (size(self%fractionModeMasses)))
       allocate(wavenumberFractionMode(size(self%fractionModeMasses)))
       allocate(slopeFractionMode     (size(self%fractionModeMasses)))
       do iMass=1,size(self%fractionModeMasses)
          massFractionMode(iMass)=self%transferFunction_%fractionModeMass(self%fractionModeMasses(iMass),statusFractionMode(iMass))
       end do
    else
       allocate(massFractionMode      (                            0))
       allocate(statusFractionMode    (                            0))
       allocate(wavenumberFractionMode(                            0))
       allocate(slopeFractionMode     (                            0))
    end if
    ! If a relative transfer function is provided, compute the relative logarithmic slope of the transfer function at the mode masses.
    wavenumberHalfMode   =-huge(0.0d0)
    wavenumberQuarterMode=-huge(0.0d0)
    slopeHalfMode        =-huge(0.0d0)
    slopeQuarterMode     =-huge(0.0d0)
    if (size(self%fractionModeMasses) > 0) then
       wavenumberFractionMode=-huge(0.0d0)
       slopeFractionMode     =-huge(0.0d0)
    end if
    if (associated(self%transferFunctionRelative )) then
       if (statusHalfModeMass    == errorStatusSuccess) then
          wavenumberHalfMode   =+self%transferFunction_       %wavenumberFromMass   (         massHalfMode)
          slopeHalfMode        =+self%transferFunction_       %logarithmicDerivative(   wavenumberHalfMode) &
               &                -self%transferFunctionRelative%logarithmicDerivative(   wavenumberHalfMode)
       end if
       if (statusQuarterModeMass == errorStatusSuccess) then
          wavenumberQuarterMode=+self%transferFunction_       %wavenumberFromMass   (      massQuarterMode)
          slopeQuarterMode     =+self%transferFunction_       %logarithmicDerivative(wavenumberQuarterMode) &
               &                -self%transferFunctionRelative%logarithmicDerivative(wavenumberQuarterMode)
       end if
       do iMass=1,size(self%fractionModeMasses)
          if (statusFractionMode(iMass) == errorStatusSuccess) then
             wavenumberFractionMode(iMass)=+self%transferFunction_       %wavenumberFromMass   (      massFractionMode(iMass))
             slopeFractionMode     (iMass)=+self%transferFunction_       %logarithmicDerivative(wavenumberFractionMode(iMass)) &
                  &                        -self%transferFunctionRelative%logarithmicDerivative(wavenumberFractionMode(iMass))
          end if
       end do
    end if
    ! If a reference transfer function is provided from which to derive a half-mode mass for mass scaling, get that now.
    if (associated(self%transferFunctionReference)) then
       massHalfModeReference      =self%transferFunctionReference%halfModeMass(statusHalfModeMassReference)
    else
       statusHalfModeMassReference=statusHalfModeMass
       massHalfModeReference      =massHalfMode
    end if
    ! Initialize warnings state.
    warnedIntegratedFailure=.false.
    ! Build a range of halo masses.
    massHaloMinimum            =self%haloMassMinimum
    massHaloMaximum            =self%haloMassMaximum
    if (self%massesRelativeToHalfModeMass) then
       if (statusHalfModeMassReference /= errorStatusSuccess) call Error_Report('half-mode mass is not defined'//{introspection:location})
       massHaloMinimum=massHaloMinimum*massHalfModeReference
       massHaloMaximum=massHaloMaximum*massHalfModeReference
    end if
    massHalo                   =Make_Range(massHaloMinimum,massHaloMaximum,int(massCount),rangeTypeLogarithmic)
    massHaloLogarithmicInterval=log(massHaloMaximum/massHaloMinimum)/dble(massCount-1)
    !$omp parallel private(iOutput,iMass,densityMean,densityCritical,basic,darkMatterProfileHalo,scaleIsSettable,shapeIsSettable,massHaloBinMinimum,massHaloBinMaximum,virialDensityContrasts,integrator_)
    allocate(haloEnvironment_                   ,mold=self%haloEnvironment_                   )
    allocate(cosmologyFunctions_                ,mold=self%cosmologyFunctions_                )
    allocate(cosmologyParameters_               ,mold=self%cosmologyParameters_               )
    allocate(virialDensityContrast_             ,mold=self%virialDensityContrast_             )
    allocate(cosmologicalMassVariance_          ,mold=self%cosmologicalMassVariance_          )
    allocate(criticalOverdensity_               ,mold=self%criticalOverdensity_               )
    allocate(haloMassFunction_                  ,mold=self%haloMassFunction_                  )
    allocate(darkMatterHaloScale_               ,mold=self%darkMatterHaloScale_               )
    allocate(unevolvedSubhaloMassFunction_      ,mold=self%unevolvedSubhaloMassFunction_      )
    allocate(darkMatterHaloBias_                ,mold=self%darkMatterHaloBias_                )
    allocate(darkMatterProfileScaleRadius_      ,mold=self%darkMatterProfileScaleRadius_      )
    allocate(darkMatterProfileShape_            ,mold=self%darkMatterProfileShape_            )
    allocate(darkMatterHaloMassAccretionHistory_,mold=self%darkMatterHaloMassAccretionHistory_)
    allocate(virialDensityContrasts(size(self%virialDensityContrasts)))
    do iAlternate=1,size(self%virialDensityContrasts)
       allocate(virialDensityContrasts(iAlternate)%virialDensityContrast_,mold=self%virialDensityContrasts(iAlternate)%virialDensityContrast_)
    end do
    !$omp critical(taskHaloMassFunctionDeepCopy)
    !![
    <deepCopyReset variables="self%haloEnvironment_ self%cosmologyFunctions_ self%cosmologyParameters_ self%virialDensityContrast_ self%cosmologicalMassVariance_ self%criticalOverdensity_ self%haloMassFunction_ self%darkMatterHaloScale_ self%unevolvedSubhaloMassFunction_ self%darkMatterHaloBias_ self%darkMatterProfileScaleRadius_ self%darkMatterProfileShape_ self%darkMatterHaloMassAccretionHistory_"/>
    <deepCopy source="self%haloEnvironment_                   " destination="haloEnvironment_                   "/>
    <deepCopy source="self%virialDensityContrast_             " destination="virialDensityContrast_             "/>
    <deepCopy source="self%cosmologyParameters_               " destination="cosmologyParameters_               "/>
    <deepCopy source="self%cosmologyFunctions_                " destination="cosmologyFunctions_                "/>
    <deepCopy source="self%cosmologicalMassVariance_          " destination="cosmologicalMassVariance_          "/>
    <deepCopy source="self%criticalOverdensity_               " destination="criticalOverdensity_               "/>
    <deepCopy source="self%haloMassFunction_                  " destination="haloMassFunction_                  "/>
    <deepCopy source="self%darkMatterHaloScale_               " destination="darkMatterHaloScale_               "/>
    <deepCopy source="self%unevolvedSubhaloMassFunction_      " destination="unevolvedSubhaloMassFunction_      "/>
    <deepCopy source="self%darkMatterHaloBias_                " destination="darkMatterHaloBias_                "/>
    <deepCopy source="self%darkMatterProfileScaleRadius_      " destination="darkMatterProfileScaleRadius_      "/>
    <deepCopy source="self%darkMatterProfileShape_            " destination="darkMatterProfileShape_            "/>
    <deepCopy source="self%darkMatterHaloMassAccretionHistory_" destination="darkMatterHaloMassAccretionHistory_"/>
    <deepCopyFinalize variables="haloEnvironment_ cosmologyFunctions_ cosmologyParameters_ virialDensityContrast_ cosmologicalMassVariance_ criticalOverdensity_ haloMassFunction_ darkMatterHaloScale_ unevolvedSubhaloMassFunction_ darkMatterHaloBias_ darkMatterProfileScaleRadius_ darkMatterProfileShape_ darkMatterHaloMassAccretionHistory_"/>
    !!]
    do iAlternate=1,size(self%virialDensityContrasts)
       !![
       <deepCopyReset variables="self%virialDensityContrasts(iAlternate)%virialDensityContrast_"/>
       <deepCopy source="self%virialDensityContrasts(iAlternate)%virialDensityContrast_" destination="virialDensityContrasts(iAlternate)%virialDensityContrast_"/>
       <deepCopyFinalize variables="virialDensityContrasts(iAlternate)%virialDensityContrast_"/>
       !!]
    end do
    !$omp end critical(taskHaloMassFunctionDeepCopy)
    ! Call routines to perform initialization which must occur for all threads if run in parallel.
    allocate(parameters)
    parameters=inputParameters(self%parameters)
    call parameters%parametersGroupCopy(self%parameters)
    call Node_Components_Thread_Initialize(parameters)
    ! Build an integrator.
    allocate(integrator_)
    integrator_=integrator(subhaloMassFunctionIntegrand,toleranceRelative=1.0d-3,integrationRule=GSL_Integ_Gauss15)
    ! Create a node object, assume zero environmental overdensity.
    allocate(tree)
    tree                   =  mergerTree()
    tree%nodeBase          => treeNode  ()
    tree%nodeBase%hostTree => tree
    call tree                   %properties%initialize(                               )
    if (haloEnvironment_%overdensityIsSettable())                                       &
         & call haloEnvironment_%overdensityLinearSet (tree%nodeBase,overdensity=0.0d0)
    allocate(tree%randomNumberGenerator_,mold=self%randomNumberGenerator_)
    !$omp critical(taskHaloMassFunctionDeepCopy)
    !![
    <deepCopyReset variables="self%randomNumberGenerator_"/>
    <deepCopy source="self%randomNumberGenerator_" destination="tree%randomNumberGenerator_"/>
    <deepCopyFinalize variables="tree%randomNumberGenerator_"/>
    !!]
    !$omp end critical(taskHaloMassFunctionDeepCopy)
    ! Get the basic and dark matter profile components.
    basic                 => tree%nodeBase%basic            (autoCreate=.true.)
    darkMatterProfileHalo => tree%nodeBase%darkMatterProfile(autoCreate=.true.)
    ! Test dark matter profile property attributes.
    scaleIsSettable       =  darkMatterProfileHalo%scaleIsSettable()
    shapeIsSettable       =  darkMatterProfileHalo%shapeIsSettable()
    ! Iterate over all output times.    
    do iOutput=outputCount,1,-1
       ! Compute characteristic densities.
       densityMean    =+cosmologyFunctions_%matterDensityEpochal(outputTimes(iOutput))
       densityCritical=+densityMean                                                    &
            &          /cosmologyFunctions_%omegaMatterEpochal  (outputTimes(iOutput))
       ! Set the time in the node.
       call basic%timeSet(outputTimes(iOutput))
       ! Loop over all halo masses.
       !$omp do
       do iMass=1,massCount
          ! Reset calculations.
          call Calculations_Reset(tree%nodeBase)
          ! Set the mass in the node.
          call                      basic                %massSet (massHalo                            (iMass        ))
          ! Set the node scale radius.
          if (scaleIsSettable) call darkMatterProfileHalo%scaleSet(darkMatterProfileScaleRadius_%radius(tree%nodeBase))
          ! Set the node shape parameter.
          if (shapeIsSettable) call darkMatterProfileHalo%shapeSet(darkMatterProfileShape_      %shape (tree%nodeBase))
          ! Get the mass distribution.
          massDistribution_ => tree%nodeBase%massDistribution(componentTypeDarkMatterOnly,massTypeDark)
          ! Compute bin interval.
          massHaloBinMinimum=massHalo(iMass)*exp(-0.5*massHaloLogarithmicInterval)
          massHaloBinMaximum=massHalo(iMass)*exp(+0.5*massHaloLogarithmicInterval)
          ! Compute halo properties.
          densityFieldRootVariance                      (iMass,iOutput)=+cosmologicalMassVariance_         %rootVariance                   (mass   =massHalo          (iMass)                                   ,time=outputTimes(iOutput)                                           )
          densityFieldRootVarianceGradientLogarithmic   (iMass,iOutput)=+cosmologicalMassVariance_         %rootVarianceLogarithmicGradient(mass   =massHalo          (iMass)                                   ,time=outputTimes(iOutput)                                           )
          if (densityFieldRootVariance(iMass,iOutput) > 0.0d0) then
             peakHeight                                 (iMass,iOutput)=+criticalOverdensity_              %value                          (mass   =massHalo          (iMass)                                   ,time=outputTimes(iOutput)                                           )    &
                  &                                                     /densityFieldRootVariance                                          (                           iMass                                    ,                 iOutput)
          else
             peakHeight                                 (iMass,iOutput)=+0.0d0
          end if
          massFunctionDifferentialLogarithmicBinAveraged(iMass,iOutput)=+haloMassFunction_                 %integrated                     (massLow=massHaloBinMinimum       ,massHigh=massHaloBinMaximum       ,time=outputTimes(iOutput),node=tree%nodeBase,status=statusIntegrated)    &
               &                                                        /massHaloLogarithmicInterval
          if (statusIntegrated /= errorStatusSuccess) then
             if (self%errorsAreFatal) then
                call    Error_Report('integrated halo mass function failed'//{introspection:location})
             else
                if (.not.warnedIntegratedFailure) then
                   warnedIntegratedFailure=.true.
                   call Warn        ('integrated halo mass function failed'                          )
                end if
             end if
          end if
          massFunctionDifferential                      (iMass,iOutput)=+haloMassFunction_                 %differential                   (mass   =massHalo          (iMass)                                   ,time=outputTimes(iOutput),node=tree%nodeBase                        )
          massFunctionCumulative                        (iMass,iOutput)=+haloMassFunction_                 %integrated                     (massLow=massHalo          (iMass),massHigh=haloMassEffectiveInfinity,time=outputTimes(iOutput),node=tree%nodeBase,status=statusIntegrated)
          if (statusIntegrated /= errorStatusSuccess) then
             if (self%errorsAreFatal) then
                call    Error_Report('integrated halo mass function failed'//{introspection:location})
             else
                if (.not.warnedIntegratedFailure) then
                   warnedIntegratedFailure=.true.
                   call Warn        ('integrated halo mass function failed'                          )
                end if
             end if
          end if
          massFunctionMassFraction                      (iMass,iOutput)=+haloMassFunction_                 %massFraction                   (massLow=massHalo          (iMass),massHigh=haloMassEffectiveInfinity,time=outputTimes(iOutput),node=tree%nodeBase                        )
          if     (                                                 &
               &   massFunctionDifferential(iMass,iOutput) > 0.0d0 &
               &  .and.                                            &
               &   densityFieldRootVariance(iMass,iOutput) > 0.0d0 &
               & ) then
             peakHeightMassFunction                     (iMass,iOutput)=+massHalo                                                          (                           iMass                                                                                                         )**2 &
                  &                                                     *massFunctionDifferential                                          (                           iMass                                    ,                 iOutput                                            )    &
                  &                                                     /cosmologyParameters_              %densityCritical                (                                                                                                                                         )    &
                  &                                                     /cosmologyParameters_              %OmegaMatter                    (                                                                                                                                         )    &
                  &                                                     /abs(densityFieldRootVariance                                      (                           iMass                                    ,                 iOutput                                            ))
          else
             peakHeightMassFunction                     (iMass,iOutput)=+0.0d0
          end if
          biasHalo                                      (iMass,iOutput)=darkMatterHaloBias_                %bias                           (                                                                                               node=tree%nodeBase                        )
          velocityVirial                                (iMass,iOutput)=darkMatterHaloScale_               %velocityVirial                 (                                                                                               node=tree%nodeBase                        )
          temperatureVirial                             (iMass,iOutput)=darkMatterHaloScale_               %temperatureVirial              (                                                                                               node=tree%nodeBase                        )
          radiusVirial                                  (iMass,iOutput)=darkMatterHaloScale_               %radiusVirial                   (                                                                                               node=tree%nodeBase                        )
          velocityMaximum                               (iMass,iOutput)=massDistribution_                  %velocityRotationCurveMaximum   (                                                                                                                 )
          darkMatterProfileRadiusScale                  (iMass,iOutput)=darkMatterProfileHalo              %scale                          (                                                                                                                                         )
          if (self%includeMassAccretionRate) &
               & massAccretionRate                      (iMass,iOutput)=darkMatterHaloMassAccretionHistory_%massAccretionRate              (                                                                     time=outputTimes(iOutput),node=tree%nodeBase                        )
          ! Compute alternate mass definitions for halos.
          do iAlternate=1,size(self%virialDensityContrasts)
             massAlternate(iAlternate,iMass,iOutput)=Dark_Matter_Profile_Mass_Definition(tree%nodeBase,virialDensityContrasts(iAlternate)%virialDensityContrast_%densityContrast(mass=massHalo(iMass),time=outputTimes(iOutput)),radius=radiusAlternate(iAlternate,iMass,iOutput),cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_,virialDensityContrast_=virialDensityContrast_)
          end do
          ! Integrate the unevolved subhalo mass function over the halo mass function to get the total subhalo mass function.
          if (self%includeUnevolvedSubhaloMassFunction)                                                                                   &
               & massFunctionCumulativeSubhalo          (iMass,iOutput)=integrator_%integrate(                                            &
               &                                                                              log(massHalo(1)/subhaloMassMaximum       ), &
               &                                                                              log(            haloMassEffectiveInfinity)  &
               &                                                                             )
          !![
	  <objectDestructor name="massDistribution_"/>
          !!]
       end do
       !$omp end do
       !$omp single
       massFunctionDifferentialLogarithmic(:,iOutput)=+massFunctionDifferential(:,iOutput) &
            &                                         *massHalo
       !$omp end single
    end do    
    call tree%destroy()
    nullify   (basic                )
    nullify   (darkMatterProfileHalo)
    deallocate(tree                 )
    deallocate(integrator_          )
    !![
    <objectDestructor name="haloEnvironment_                   "/>
    <objectDestructor name="virialDensityContrast_             "/>
    <objectDestructor name="cosmologyParameters_               "/>
    <objectDestructor name="cosmologyFunctions_                "/>
    <objectDestructor name="criticalOverdensity_               "/>
    <objectDestructor name="cosmologicalMassVariance_          "/>
    <objectDestructor name="haloMassFunction_                  "/>
    <objectDestructor name="darkMatterHaloScale_               "/>
    <objectDestructor name="unevolvedSubhaloMassFunction_      "/>
    <objectDestructor name="darkMatterHaloBias_                "/>
    <objectDestructor name="darkMatterProfileScaleRadius_      "/>
    <objectDestructor name="darkMatterProfileShape_            "/>
    <objectDestructor name="darkMatterHaloMassAccretionHistory_"/>
    !!]
    do iAlternate=1,size(self%virialDensityContrasts)
       !![
       <objectDestructor name="virialDensityContrasts(iAlternate)%virialDensityContrast_"/>
       !!]
    end do
    deallocate(virialDensityContrasts)
    call Node_Components_Thread_Uninitialize()
    !$omp end parallel
    ! Open the group for output time information.
    if (self%outputGroup == ".") then
       outputsGroup  =outputFile    %openGroup(     'Outputs'        ,'Group containing datasets relating to output times.')
    else
       containerGroup=outputFile    %openGroup(char(self%outputGroup),'Group containing halo mass function data.'          )
       outputsGroup  =containerGroup%openGroup(     'Outputs'        ,'Group containing datasets relating to output times.')
    end if
    ! Store half- and quarter-mode masses if possible.
    if     (                                             &
         &   statusHalfModeMass    == errorStatusSuccess &
         &  .or.                                         &
         &   statusQuarterModeMass == errorStatusSuccess &
         & ) then
       if (self%outputGroup == ".") then
          powerSpectrumGroup=outputFile%openGroup('powerSpectrum','Group containing data relating to the power spectrum.')
       else
          powerSpectrumGroup=containerGroup      %openGroup('powerSpectrum','Group containing data relating to the power spectrum.')
       end if
       if (statusHalfModeMass    == errorStatusSuccess) then
          call powerSpectrumGroup%writeAttribute(         massHalfMode,                  'massHalfMode')
          call powerSpectrumGroup%writeAttribute(   wavenumberHalfMode,            'wavenumberHalfMode')
          call powerSpectrumGroup%writeAttribute(        slopeHalfMode,   'logarithmicGradientHalfMode')
       end if
       if (statusQuarterModeMass == errorStatusSuccess) then
          call powerSpectrumGroup%writeAttribute(      massQuarterMode,               'massQuarterMode')
          call powerSpectrumGroup%writeAttribute(wavenumberQuarterMode,         'wavenumberQuarterMode')
          call powerSpectrumGroup%writeAttribute(     slopeQuarterMode,'logarithmicGradientQuarterMode')
       end if
       if (size(self%fractionModeMasses) > 0) then
          do iMass=1,size(self%fractionModeMasses)
             write (label,'(e8.2)') self%fractionModeMasses(iMass)
             if (statusFractionMode(iMass) == errorStatusSuccess) then
                call powerSpectrumGroup%writeAttribute(      massFractionMode(iMass),               'massFractionMode_'//trim(label))
                call powerSpectrumGroup%writeAttribute(wavenumberFractionMode(iMass),         'wavenumberFractionMode_'//trim(label))
                call powerSpectrumGroup%writeAttribute(     slopeFractionMode(iMass),'logarithmicGradientFractionMode_'//trim(label))
             end if
          end do
       end if
       call                                                            powerSpectrumGroup%close         (                                   )
    end if
    ! Store σ₈.
    if (self%outputGroup == ".") then
       powerSpectrumGroup=outputFile%openGroup('powerSpectrum','Group containing data relating to the power spectrum.')
    else
       powerSpectrumGroup=containerGroup      %openGroup('powerSpectrum','Group containing data relating to the power spectrum.')
    end if
    call powerSpectrumGroup%writeAttribute(self%cosmologicalMassVariance_%sigma8(),'sigma8')
    call powerSpectrumGroup%close         (                                                )
    ! Store other usual information.
    if (self%outputGroup == ".") then
       cosmologyGroup=outputFile    %openGroup('cosmology','Group containing data relating to cosmology.')
    else
       cosmologyGroup=containerGroup%openGroup('cosmology','Group containing data relating to cosmology.')
    end if
    call cosmologyGroup%writeAttribute(self%cosmologyParameters_%densityCritical(),'densityCritical')
    call cosmologyGroup%close()
    ! Iterate over output times and output data.
    do iOutput=1,outputCount
       groupName  ='Output'
       description='Data for output number '
       groupName  =groupName  //iOutput
       description=description//iOutput
       outputGroup=outputsGroup%openGroup(char(groupName),char(description))
       call    outputGroup%writeAttribute(outputTimes                                   (  iOutput),'outputTime'                                                                                                                          )
       call    outputGroup%writeAttribute(outputRedshifts                               (  iOutput),'outputRedshift'                                                                                                                      )
       call    outputGroup%writeAttribute(outputExpansionFactors                        (  iOutput),'outputExpansionFactor'                                                                                                               )
       call    outputGroup%writeAttribute(outputGrowthFactors                           (  iOutput),'growthFactor'                                                                                                                        )
       call    outputGroup%writeAttribute(outputCriticalOverdensities                   (  iOutput),'criticalOverdensity'                                                                                                                 )
       call    outputGroup%writeAttribute(outputVirialDensityContrast                   (  iOutput),'virialDensityContrast'                                                                                                               )
       call    outputGroup%writeAttribute(outputTurnaroundRadius                        (  iOutput),'turnaroundToVirialRadiusRatio'                                                                                                       )
       call    outputGroup%writeAttribute(outputCharacteristicMass                      (  iOutput),'massHaloCharacteristic'                                                                                                              )
       if (self%massesRelativeToHalfModeMass) then
          massHaloOutput=massHalo/massHalfModeReference
          call outputGroup%writeDataset  (massHaloOutput                                (:        ),'haloMass'                      ,'The mass of the halo relative to the half-mode mass.'                                               )
       else
          call outputGroup%writeDataset  (massHalo                                      (:        ),'haloMass'                      ,'The mass of the halo.'                                                      ,datasetReturned=dataset)
          call dataset    %writeAttribute(massSolar                                                ,'unitsInSI'                                                                                                                           )
          call dataset    %close         (                                                                                                                                                                                                )
       end if
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
       if (self%includeMassAccretionRate) then
          call outputGroup%writeDataset  (massAccretionRate                             (:,iOutput),'haloMassAccretionRate'         ,'The mass accretion rate of halos.'                                          ,datasetReturned=dataset)
          call dataset    %writeAttribute(massSolar/gigaYear                                       ,'unitsInSI'                                                                                                                           )
          call dataset    %close         (                                                                                                                                                                                                )
       end if
       call    outputGroup%writeDataset  (biasHalo                                      (:,iOutput),'haloBias'                      ,'The large scale linear bias of halos.'                                                              )
       call    outputGroup%writeDataset  (densityFieldRootVariance                      (:,iOutput),'haloSigma'                     ,'The mass fluctuation on the scale of the halo.'                                                     )
       call    outputGroup%writeDataset  (densityFieldRootVarianceGradientLogarithmic   (:,iOutput),'haloAlpha'                     ,'dlog(σ)/dlog(m).'                                                                                   )
       call    outputGroup%writeDataset  (peakHeight                                    (:,iOutput),'haloPeakHeightNu'              ,'The peak height, ν, of the halo.'                                                                   )
       call    outputGroup%writeDataset  (massFunctionMassFraction                      (:,iOutput),'haloMassFractionCumulative'    ,'The halo cumulative mass fraction.'                                                                 )
       call    outputGroup%writeDataset  (peakHeightMassFunction                        (:,iOutput),'haloMassFunctionNuFNu'         ,'The halo mass fraction function as a function of the ν parameter, νF(ν).'                           )
       if (size(self%virialDensityContrasts) > 0) then
          do iAlternate=1,size(self%virialDensityContrasts)
             if (self%massesRelativeToHalfModeMass) then
                massHaloOutput=massAlternate(iAlternate,:,iOutput)/massHalfModeReference
                call outputGroup%writeDataset  (massHaloOutput                       ,'haloMass'  //String_Upper_Case_First(char(self%virialDensityContrasts(iAlternate)%label)),'The mass of the halo under the "'  //char(self%virialDensityContrasts(iAlternate)%label)//'" definition relative to the half-mode mass.'                        )
           else
                call outputGroup%writeDataset  (massAlternate  (iAlternate,:,iOutput),'haloMass'  //String_Upper_Case_First(char(self%virialDensityContrasts(iAlternate)%label)),'The mass of the halo under the "'  //char(self%virialDensityContrasts(iAlternate)%label)//'" definition.'                               ,datasetReturned=dataset)
                call dataset    %writeAttribute(massSolar                            ,'unitsInSI'                                                                                                                                                                                                                                                 )
                call dataset    %close         (                                                                                                                                                                                                                                                                                                  )
             end if
             call    outputGroup%writeDataset  (radiusAlternate(iAlternate,:,iOutput),'haloRadius'//String_Upper_Case_First(char(self%virialDensityContrasts(iAlternate)%label)),'The radius of the halo under the "'//char(self%virialDensityContrasts(iAlternate)%label)//'" definition.'                               ,datasetReturned=dataset)
             call    dataset    %writeAttribute(megaParsec                           ,'unitsInSI'                                                                                                                                                                                                                                                 )
             call    dataset    %close         (                                                                                                                                                                                                                                                                                                  )
          end do
       end if
       call outputGroup%close()
    end do
    call outputsGroup%close()
    if (containerGroup%isOpen()) call containerGroup%close()
    if (present(status)) status=errorStatusSuccess
    call Node_Components_Thread_Uninitialize()
    call displayUnindent('Done task: halo mass function' )
    return

  contains

    double precision function subhaloMassFunctionIntegrand(logMass)
      !!{
      Integrand function used to find the cumulative subhalo mass function.
      !!}
      implicit none
      double precision, intent(in   ) :: logMass
      double precision                :: mass

      ! Extract integrand parameters.
      mass=exp(logMass)
      ! Return the differential halo mass function multiplied by the integrated unevolved subhalo mass function in such hosts.
      subhaloMassFunctionIntegrand=+                                                                                                          mass                     &
           &                       *haloMassFunction_            %differential(outputTimes(iOutput)                                          ,mass,node=tree%nodeBase) &
           &                       *unevolvedSubhaloMassFunction_%integrated  (outputTimes(iOutput),massHalo(iMass),haloMassEffectiveInfinity,mass                   )
      return
    end function subhaloMassFunctionIntegrand

  end subroutine haloMassFunctionPerform

  subroutine haloMassFunctionDescriptorSpecial(self,descriptor)
    !!{
    Add special parameters to the descriptor.
    !!}
    use :: Input_Parameters  , only : inputParameters
    use :: ISO_Varying_String, only : char
    use :: String_Handling   , only : String_Join
    implicit none
    class  (taskHaloMassFunction), intent(inout) :: self
    type   (inputParameters     ), intent(inout) :: descriptor
    type   (inputParameters     )                :: subParameters
    integer                                      :: i
    
    if (allocated(self%virialDensityContrasts)) then
       subParameters=descriptor%subparameters('massDefinitions')
       call subParameters%addParameter('labels',char(String_Join(self%virialDensityContrasts%label," ")))
       do i=1,size(self%virialDensityContrasts)
          call self%virialDensityContrasts(i)%virialDensityContrast_%descriptor(subParameters)
       end do
    end if
    return
  end subroutine haloMassFunctionDescriptorSpecial
