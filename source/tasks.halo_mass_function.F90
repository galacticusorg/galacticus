!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
  use :: Dark_Matter_Profiles_DMO                 , only : darkMatterProfileDMOClass
  use :: Halo_Mass_Functions                      , only : haloMassFunctionClass
  use :: Linear_Growth                            , only : linearGrowthClass
  use :: Output_Times                             , only : outputTimesClass
  use :: Transfer_Functions                       , only : transferFunctionClass
  use :: Unevolved_Subhalo_Mass_Functions         , only : unevolvedSubhaloMassFunctionClass
  use :: Virial_Density_Contrast                  , only : virialDensityContrastClass

  type :: virialDensityContrastList
     !!{
     Type used to store a list of virial density contrasts.
     !!}
     type (varying_string            )          :: label
     class(virialDensityContrastClass), pointer :: virialDensityContrast_
  end type virialDensityContrastList
  
  !![
  <task name="taskHaloMassFunction">
   <description>A task which computes and outputs the halo mass function and related quantities.</description>
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
     class           (darkMatterProfileDMOClass              ), pointer                   :: darkMatterProfileDMO_               => null()
     class           (criticalOverdensityClass               ), pointer                   :: criticalOverdensity_                => null()
     class           (linearGrowthClass                      ), pointer                   :: linearGrowth_                       => null()
     class           (haloMassFunctionClass                  ), pointer                   :: haloMassFunction_                   => null()
     class           (haloEnvironmentClass                   ), pointer                   :: haloEnvironment_                    => null()
     class           (unevolvedSubhaloMassFunctionClass      ), pointer                   :: unevolvedSubhaloMassFunction_       => null()
     class           (darkMatterHaloScaleClass               ), pointer                   :: darkMatterHaloScale_                => null()
     class           (cosmologicalMassVarianceClass          ), pointer                   :: cosmologicalMassVariance_           => null()
     class           (darkMatterHaloBiasClass                ), pointer                   :: darkMatterHaloBias_                 => null()
     class           (transferFunctionClass                  ), pointer                   :: transferFunction_                   => null(), transferFunctionReference => null()
     class           (outputTimesClass                       ), pointer                   :: outputTimes_                        => null()
     class           (darkMatterProfileScaleRadiusClass      ), pointer                   :: darkMatterProfileScaleRadius_       => null()
     class           (darkMatterHaloMassAccretionHistoryClass), pointer                   :: darkMatterHaloMassAccretionHistory_ => null()
     double precision                                                                     :: haloMassMinimum                              , haloMassMaximum                    , &
          &                                                                                  pointsPerDecade
     type            (varying_string                         )                            :: outputGroup
     logical                                                                              :: includeUnevolvedSubhaloMassFunction          , includeMassAccretionRate           , &
          &                                                                                  massesRelativeToHalfModeMass
     double precision                                         , allocatable, dimension(:) :: fractionModeMasses
     type            (virialDensityContrastList              ), allocatable, dimension(:) :: virialDensityContrasts
     ! Pointer to the parameters for this task.
     type            (inputParameters                        )                            :: parameters
  contains
     final     ::            haloMassFunctionDestructor
     procedure :: perform => haloMassFunctionPerform
  end type taskHaloMassFunction

  interface taskHaloMassFunction
     !!{
     Constructors for the {\normalfont \ttfamily haloMassFunction} task.
     !!}
     module procedure haloMassFunctionConstructorParameters
     module procedure haloMassFunctionConstructorInternal
  end interface taskHaloMassFunction

contains

  function haloMassFunctionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily haloMassFunction} task class which takes a parameter set as input.
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
    class           (darkMatterProfileDMOClass              ), pointer                     :: darkMatterProfileDMO_
    class           (criticalOverdensityClass               ), pointer                     :: criticalOverdensity_
    class           (linearGrowthClass                      ), pointer                     :: linearGrowth_
    class           (haloMassFunctionClass                  ), pointer                     :: haloMassFunction_
    class           (haloEnvironmentClass                   ), pointer                     :: haloEnvironment_
    class           (unevolvedSubhaloMassFunctionClass      ), pointer                     :: unevolvedSubhaloMassFunction_
    class           (darkMatterHaloScaleClass               ), pointer                     :: darkMatterHaloScale_
    class           (cosmologicalMassVarianceClass          ), pointer                     :: cosmologicalMassVariance_
    class           (darkMatterHaloBiasClass                ), pointer                     :: darkMatterHaloBias_
    class           (transferFunctionClass                  ), pointer                     :: transferFunction_                  , transferFunctionReference
    class           (outputTimesClass                       ), pointer                     :: outputTimes_
    class           (darkMatterProfileScaleRadiusClass      ), pointer                     :: darkMatterProfileScaleRadius_
    class           (darkMatterHaloMassAccretionHistoryClass), pointer                     :: darkMatterHaloMassAccretionHistory_
    type            (inputParameters                        ), pointer                     :: parametersRoot
    type            (inputParameters                        ),                             :: parametersMassDefinitions
    type            (virialDensityContrastList              ), allocatable  , dimension(:) :: virialDensityContrasts
    type            (varying_string                         ), allocatable  , dimension(:) :: labels
    double precision                                         , allocatable  , dimension(:) :: fractionModeMasses
    type            (varying_string                         )                              :: outputGroup
    double precision                                                                       :: haloMassMinimum                    , haloMassMaximum           , &
         &                                                                                    pointsPerDecade
    logical                                                                                :: includeUnevolvedSubhaloMassFunction, includeMassAccretionRate  , &
          &                                                                                   massesRelativeToHalfModeMass
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
    <objectBuilder    class="cosmologyParameters"                name="cosmologyParameters_"                source="parameters"                                          />
    <objectBuilder    class="cosmologyFunctions"                 name="cosmologyFunctions_"                 source="parameters"                                          />
    <objectBuilder    class="virialDensityContrast"              name="virialDensityContrast_"              source="parameters"                                          />
    <objectBuilder    class="darkMatterProfileDMO"               name="darkMatterProfileDMO_"               source="parameters"                                          />
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
    <objectBuilder    class="darkMatterHaloMassAccretionHistory" name="darkMatterHaloMassAccretionHistory_" source="parameters"                                          />
    !!]
    if (parameters%isPresent('transferFunctionReference')) then
       !![
       <objectBuilder class="transferFunction"                   name="transferFunction_"                   source="parameters" parameterName="transferFunctionReference"/>
       !!]
    end if
    if (parameters%isPresent('massDefinitions',requireValue=.false.)) then
       parametersMassDefinitions=parameters%subParameters('massDefinitions',requireValue=.false.)
       if (parametersMassDefinitions%copiesCount('virialDensityContrast') /= parametersMassDefinitions%count('labels')) &
            & call Galacticus_Error_Report('number of labels must match number of virial density contrasts'//{introspection:location})
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
	 &amp;                    fractionModeMasses                 , &amp;
         &amp;                    cosmologyParameters_               , &amp;
         &amp;                    cosmologyFunctions_                , &amp;
         &amp;                    virialDensityContrast_             , &amp;
         &amp;                    darkMatterProfileDMO_              , &amp;
         &amp;                    criticalOverdensity_               , &amp;
         &amp;                    linearGrowth_                      , &amp;
         &amp;                    haloMassFunction_                  , &amp;
         &amp;                    haloEnvironment_                   , &amp;
         &amp;                    unevolvedSubhaloMassFunction_      , &amp;
         &amp;                    darkMatterHaloScale_               , &amp;
         &amp;                    darkMatterProfileScaleRadius_      , &amp;
         &amp;                    darkMatterHaloMassAccretionHistory_, &amp;
         &amp;                    cosmologicalMassVariance_          , &amp;
         &amp;                    darkMatterHaloBias_                , &amp;
         &amp;                    transferFunction_                  , &amp;
         &amp;                    outputTimes_                       , &amp;
         &amp;                    virialDensityContrasts             , &amp;
         &amp;                    parametersRoot                       &amp;
         &amp;                    {conditions}                         &amp;
         &amp;                   )
     </call>
     <argument name="transferFunctionReference" value="transferFunctionReference" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"               />
    <objectDestructor name="cosmologyFunctions_"                />
    <objectDestructor name="virialDensityContrast_"             />
    <objectDestructor name="darkMatterProfileDMO_"              />
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
    <objectDestructor name="darkMatterHaloMassAccretionHistory_"/>
    !!]
    if (parameters%isPresent('transferFunctionReference')) then
       !![
       <objectDestructor name="transferFunctionReference"/>
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
       &                                       fractionModeMasses                 , &
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
       &                                       darkMatterHaloMassAccretionHistory_, &
       &                                       cosmologicalMassVariance_          , &
       &                                       darkMatterHaloBias_                , &
       &                                       transferFunction_                  , &
       &                                       outputTimes_                       , &
       &                                       virialDensityContrasts             , &
       &                                       parameters                         , &
       &                                       transferFunctionReference            &
       &                                      ) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily haloMassFunction} task class which takes a parameter set as input.
    !!}
    implicit none
    type            (taskHaloMassFunction                   )                                        :: self
    class           (cosmologyParametersClass               ), intent(in   ), target                 :: cosmologyParameters_
    class           (cosmologyFunctionsClass                ), intent(in   ), target                 :: cosmologyFunctions_
    class           (virialDensityContrastClass             ), intent(in   ), target                 :: virialDensityContrast_
    class           (darkMatterProfileDMOClass              ), intent(in   ), target                 :: darkMatterProfileDMO_
    class           (criticalOverdensityClass               ), intent(in   ), target                 :: criticalOverdensity_
    class           (linearGrowthClass                      ), intent(in   ), target                 :: linearGrowth_
    class           (haloMassFunctionClass                  ), intent(in   ), target                 :: haloMassFunction_
    class           (haloEnvironmentClass                   ), intent(in   ), target                 :: haloEnvironment_
    class           (unevolvedSubhaloMassFunctionClass      ), intent(in   ), target                 :: unevolvedSubhaloMassFunction_
    class           (darkMatterHaloScaleClass               ), intent(in   ), target                 :: darkMatterHaloScale_
    class           (darkMatterProfileScaleRadiusClass      ), intent(in   ), target                 :: darkMatterProfileScaleRadius_
    class           (darkMatterHaloMassAccretionHistoryClass), intent(in   ), target                 :: darkMatterHaloMassAccretionHistory_
    class           (cosmologicalMassVarianceClass          ), intent(in   ), target                 :: cosmologicalMassVariance_
    class           (darkMatterHaloBiasClass                ), intent(in   ), target                 :: darkMatterHaloBias_
    class           (transferFunctionClass                  ), intent(in   ), target                 :: transferFunction_
    class           (transferFunctionClass                  ), intent(in   ), target      , optional :: transferFunctionReference
    class           (outputTimesClass                       ), intent(in   ), target                 :: outputTimes_
    type            (virialDensityContrastList              ), intent(in   ), dimension(:)           :: virialDensityContrasts
    type            (varying_string                         ), intent(in   )                         :: outputGroup
    double precision                                         , intent(in   )                         :: haloMassMinimum                    , haloMassMaximum           , &
         &                                                                                              pointsPerDecade
    logical                                                  , intent(in   )                         :: includeUnevolvedSubhaloMassFunction, includeMassAccretionRate  , &
         &                                                                                              massesRelativeToHalfModeMass
    double precision                                         , intent(in   ), dimension(:)           :: fractionModeMasses
    type            (inputParameters                        ), intent(in   ), target                 :: parameters
    integer                                                                                          :: i
    !![
    <constructorAssign variables="haloMassMinimum, haloMassMaximum, pointsPerDecade, outputGroup, includeUnevolvedSubhaloMassFunction, includeMassAccretionRate, massesRelativeToHalfModeMass, fractionModeMasses, *cosmologyParameters_, *cosmologyFunctions_, *virialDensityContrast_, *darkMatterProfileDMO_, *criticalOverdensity_, *linearGrowth_, *haloMassFunction_, *haloEnvironment_, *unevolvedSubhaloMassFunction_, *darkMatterHaloScale_, *darkMatterProfileScaleRadius_, *darkMatterHaloMassAccretionHistory_, *cosmologicalMassVariance_, *darkMatterHaloBias_, *transferFunction_, *transferFunctionReference, *outputTimes_"/>
    !!]

    self%parameters=inputParameters(parameters)
    call self%parameters%parametersGroupCopy(parameters)
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
    Destructor for the {\normalfont \ttfamily haloMassFunction} task class.
    !!}
    use :: Node_Components, only : Node_Components_Uninitialize
    implicit none
    type   (taskHaloMassFunction), intent(inout) :: self
    integer                                      :: i

    !![
    <objectDestructor name="self%cosmologyParameters_"               />
    <objectDestructor name="self%cosmologyFunctions_"                />
    <objectDestructor name="self%virialDensityContrast_"             />
    <objectDestructor name="self%darkMatterProfileDMO_"              />
    <objectDestructor name="self%criticalOverdensity_"               />
    <objectDestructor name="self%linearGrowth_"                      />
    <objectDestructor name="self%haloMassFunction_"                  />
    <objectDestructor name="self%haloEnvironment_"                   />
    <objectDestructor name="self%unevolvedSubhaloMassFunction_"      />
    <objectDestructor name="self%darkMatterHaloScale_"               />
    <objectDestructor name="self%darkMatterProfileScaleRadius_"      />
    <objectDestructor name="self%darkMatterHaloMassAccretionHistory_"/>
    <objectDestructor name="self%cosmologicalMassVariance_"          />
    <objectDestructor name="self%darkMatterHaloBias_"                />
    <objectDestructor name="self%transferFunction_"                  />
    <objectDestructor name="self%outputTimes_"                       />
    !!]
    if (associated(self%transferFunctionReference)) then
       !![
       <objectDestructor name="self%transferFunctionReference"/>
       !!]
    end if
    if (allocated(elf%virialDensityContrasts)) then
       do i=1,size(self%virialDensityContrasts)
          !![
	  <objectDestructor name="self%virialDensityContrasts(i)%virialDensityContrast_"/>
          !!]
       end do
    end if
    call Node_Components_Uninitialize()
    return
  end subroutine haloMassFunctionDestructor

  subroutine haloMassFunctionPerform(self,status)
    !!{
    Compute and output the halo mass function.
    !!}
    use            :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use            :: Display                             , only : displayIndent                      , displayUnindent
    use            :: Galacticus_Calculations_Resets      , only : Galacticus_Calculations_Reset
    use            :: Galacticus_Error                    , only : errorStatusSuccess
    use            :: Galacticus_HDF5                     , only : galacticusOutputFile
    use            :: Galacticus_Nodes                    , only : mergerTree                         , nodeComponentBasic                 , nodeComponentDarkMatterProfile, treeNode
    use            :: IO_HDF5                             , only : hdf5Object
    use, intrinsic :: ISO_C_Binding                       , only : c_size_t
    use            :: Memory_Management                   , only : allocateArray
    use            :: Node_Components                     , only : Node_Components_Thread_Initialize  , Node_Components_Thread_Uninitialize
    use            :: Numerical_Constants_Astronomical    , only : massSolar                          , megaParsec                         , gigaYear
    use            :: Numerical_Constants_Math            , only : Pi
    use            :: Numerical_Constants_Prefixes        , only : kilo
    use            :: Numerical_Integration               , only : GSL_Integ_Gauss15                  , integrator
    use            :: Numerical_Ranges                    , only : Make_Range                         , rangeTypeLogarithmic
    use            :: String_Handling                     , only : operator(//)                       , String_Upper_Case_First
    implicit none
    class           (taskHaloMassFunction          ), intent(inout), target           :: self
    integer                                         , intent(  out), optional         :: status
    double precision                                , allocatable  , dimension(:,:)   :: massFunctionDifferentialLogarithmicBinAveraged       , biasHalo                     , &
         &                                                                               massFunctionCumulative                               , massFunctionDifferential     , &
         &                                                                               massFunctionDifferentialLogarithmic                  , massFunctionMassFraction     , &
         &                                                                               peakHeight                                           , densityFieldRootVariance     , &
         &                                                                               radiusVirial                                         , temperatureVirial            , &
         &                                                                               velocityVirial                                       , darkMatterProfileRadiusScale , &
         &                                                                               velocityMaximum                                      , peakHeightMassFunction       , &
         &                                                                               densityFieldRootVarianceGradientLogarithmic          , massFunctionCumulativeSubhalo, &
         &                                                                               massAccretionRate
    double precision                                , allocatable  , dimension(:  )   :: outputCharacteristicMass                             , outputCriticalOverdensities  , &
         &                                                                               outputExpansionFactors                               , outputGrowthFactors          , &
         &                                                                               outputRedshifts                                      , outputTimes                  , &
         &                                                                               outputVirialDensityContrast                          , outputTurnAroundRadius       , &
         &                                                                               massHalo                                             , massHaloOutput               , &
         &                                                                               massFractionMode
    double precision                                , allocatable  , dimension(:,:,:) :: massAlternate                                        , radiusAlternate
    integer                                         , allocatable  , dimension(:    ) :: statusFractionMode
    ! The upper limit to halo mass used when computing cumulative mass functions.
    double precision                                , parameter                       :: haloMassEffectiveInfinity                     =1.0d16
    ! Largest subhalo mass (in units of host mass) for which we expect significant unevolved subhalo mass function.
    double precision                                , parameter                       :: subhaloMassMaximum                            =1.0d02
    class           (nodeComponentBasic            ), pointer                         :: basic
    class           (nodeComponentDarkMatterProfile), pointer                         :: darkMatterProfileHalo
    type            (mergerTree                    ), target                          :: tree
    type            (integrator                    )                                  :: integrator_
    integer         (c_size_t                      )                                  :: iOutput                                              , outputCount                  , &
         &                                                                               iMass                                                , massCount                    , &
         &                                                                               iAlternate
    double precision                                                                  :: massHaloBinMinimum                                   , massHaloBinMaximum           , &
         &                                                                               massHaloLogarithmicInterval                          , massHalfMode                 , &
         &                                                                               wavenumberComoving                                   , growthFactor                 , &
         &                                                                               massCritical                                         , massQuarterMode              , &
         &                                                                               densityMean                                          , densityCritical              , &
         &                                                                               massHaloMinimum                                      , massHaloMaximum              , &
         &                                                                               massHalfModeReference
    type            (hdf5Object                    )                                  :: outputsGroup                                         , outputGroup                  , &
         &                                                                               containerGroup                                       , powerSpectrumGroup           , &
         &                                                                               cosmologyGroup                                       , dataset
    integer                                                                           :: statusHalfModeMass                                   , statusQuarterModeMass        , &
         &                                                                               statusHalfModeMassReference
    type            (varying_string                )                                  :: groupName                                            , commentText
    character       (len=32                        )                                  :: label
    
    call displayIndent('Begin task: halo mass function')
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
    massCount=int(log10(self%haloMassMaximum/self%haloMassMinimum)*self%pointsPerDecade)+1
    call allocateArray(massHalo                                      ,[                                                massCount            ])
    call allocateArray(massHaloOutput                                ,[                                                massCount            ])
    call allocateArray(massFunctionDifferential                      ,[                                                massCount,outputCount])
    call allocateArray(massFunctionDifferentialLogarithmic           ,[                                                massCount,outputCount])
    call allocateArray(massFunctionDifferentialLogarithmicBinAveraged,[                                                massCount,outputCount])
    call allocateArray(massFunctionCumulative                        ,[                                                massCount,outputCount])
    call allocateArray(massFunctionCumulativeSubhalo                 ,[                                                massCount,outputCount])
    call allocateArray(massFunctionMassFraction                      ,[                                                massCount,outputCount])
    call allocateArray(massAccretionRate                             ,[                                                massCount,outputCount])
    call allocateArray(biasHalo                                      ,[                                                massCount,outputCount])
    call allocateArray(densityFieldRootVariance                      ,[                                                massCount,outputCount])
    call allocateArray(densityFieldRootVarianceGradientLogarithmic   ,[                                                massCount,outputCount])
    call allocateArray(peakHeight                                    ,[                                                massCount,outputCount])
    call allocateArray(peakHeightMassFunction                        ,[                                                massCount,outputCount])
    call allocateArray(velocityVirial                                ,[                                                massCount,outputCount])
    call allocateArray(temperatureVirial                             ,[                                                massCount,outputCount])
    call allocateArray(radiusVirial                                  ,[                                                massCount,outputCount])
    call allocateArray(darkMatterProfileRadiusScale                  ,[                                                massCount,outputCount])
    call allocateArray(velocityMaximum                               ,[                                                massCount,outputCount])
    call allocateArray(massAlternate                                 ,[size(self%virialDensityContrasts,kind=c_size_t),massCount,outputCount])
    call allocateArray(radiusAlternate                               ,[size(self%virialDensityContrasts,kind=c_size_t),massCount,outputCount])
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
    ! Create a node object, assume zero environmental overdensity.
    tree%baseNode          => treeNode()
    tree%baseNode%hostTree => tree
    call tree%properties      %initialize          (                               )
    call self%haloEnvironment_%overdensityLinearSet(tree%baseNode,overdensity=0.0d0)
    ! Get half- and quarter-mode masses.
    massHalfMode   =self%transferFunction_%halfModeMass   (statusHalfModeMass   )
    massQuarterMode=self%transferFunction_%quarterModeMass(statusQuarterModeMass)
    if (size(self%fractionModeMasses) > 0) then
       allocate(massFractionMode  (size(self%fractionModeMasses)))
       allocate(statusFractionMode(size(self%fractionModeMasses)))
       do iMass=1,size(self%fractionModeMasses)
          massFractionMode(iMass)=self%transferFunction_%fractionModeMass(self%fractionModeMasses(iMass),statusFractionMode(iMass))
       end do
    end if
    ! If a reference transfer function is provided from which to derive a half-mode mass for mass scaling, get that now.
    if (associated(self%transferFunctionReference)) then
       massHalfModeReference      =self%transferFunctionReference%halfModeMass(statusHalfModeMassReference)
    else
       statusHalfModeMassReference=statusHalfModeMass
       massHalfModeReference      =massHalfMode
    end if
    ! Build a range of halo masses.
    massHaloMinimum            =self%haloMassMinimum
    massHaloMaximum            =self%haloMassMaximum
    if (self%massesRelativeToHalfModeMass) then
       if (statusHalfModeMassReference /= errorStatusSuccess) call Galacticus_Error_Report('half-mode mass is not defined'//{introspection:location})
       massHaloMinimum=massHaloMinimum*massHalfModeReference
       massHaloMaximum=massHaloMaximum*massHalfModeReference
    end if
    massHalo                   =Make_Range(massHaloMinimum,massHaloMaximum,int(massCount),rangeTypeLogarithmic)
    massHaloLogarithmicInterval=log(massHaloMaximum/massHaloMinimum)/dble(massCount-1)
    ! Get the basic and dark matter profile components.
    basic                 => tree%baseNode%basic            (autoCreate=.true.)
    darkMatterProfileHalo => tree%baseNode%darkMatterProfile(autoCreate=.true.)
    ! Build an integrator.
    integrator_=integrator(subhaloMassFunctionIntegrand,toleranceRelative=1.0d-3,integrationRule=GSL_Integ_Gauss15)
    ! Iterate over all output times.    
    do iOutput=outputCount,1,-1
       ! Compute characteristic densities.
       densityMean    =+self%cosmologyFunctions_%matterDensityEpochal(outputTimes(iOutput))
       densityCritical=+densityMean                                                         &
            &          /self%cosmologyFunctions_%omegaMatterEpochal  (outputTimes(iOutput))
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
          densityFieldRootVariance                      (iMass,iOutput)=+self%cosmologicalMassVariance_         %rootVariance                   (mass   =massHalo          (iMass)                                   ,time=outputTimes(iOutput)                   )
          densityFieldRootVarianceGradientLogarithmic   (iMass,iOutput)=+self%cosmologicalMassVariance_         %rootVarianceLogarithmicGradient(mass   =massHalo          (iMass)                                   ,time=outputTimes(iOutput)                   )
          peakHeight                                    (iMass,iOutput)=+self%criticalOverdensity_              %value                          (mass   =massHalo          (iMass)                                   ,time=outputTimes(iOutput)                   )    &
               &                                                        /densityFieldRootVariance                                               (                           iMass                                    ,                 iOutput)
          massFunctionDifferentialLogarithmicBinAveraged(iMass,iOutput)=+self%haloMassFunction_                 %integrated                     (massLow=massHaloBinMinimum       ,massHigh=massHaloBinMaximum       ,time=outputTimes(iOutput),node=tree%baseNode)    &
               &                                                        /massHaloLogarithmicInterval
          massFunctionDifferential                      (iMass,iOutput)=+self%haloMassFunction_                 %differential                   (mass   =massHalo          (iMass)                                   ,time=outputTimes(iOutput),node=tree%baseNode)
          massFunctionCumulative                        (iMass,iOutput)=+self%haloMassFunction_                 %integrated                     (massLow=massHalo          (iMass),massHigh=haloMassEffectiveInfinity,time=outputTimes(iOutput),node=tree%baseNode)
          massFunctionMassFraction                      (iMass,iOutput)=+self%haloMassFunction_                 %massFraction                   (massLow=massHalo          (iMass),massHigh=haloMassEffectiveInfinity,time=outputTimes(iOutput),node=tree%baseNode)
          peakHeightMassFunction                        (iMass,iOutput)=+massHalo                                                               (                           iMass                                                                                 )**2 &
               &                                                        *massFunctionDifferential                                               (                           iMass                                    ,                 iOutput                    )    &
               &                                                        /self%cosmologyParameters_              %densityCritical                (                                                                                                                 )    &
               &                                                        /self%cosmologyParameters_              %OmegaMatter                    (                                                                                                                 )    &
               &                                                        /abs(                                                                                                                                                                                          &
               &                                                             self%cosmologicalMassVariance_     %rootVarianceLogarithmicGradient(mass   =massHalo          (iMass)                                   ,time=outputTimes(iOutput)                   )    &
               &                                                            )
          biasHalo                                      (iMass,iOutput)=self%darkMatterHaloBias_                %bias                           (                                                                                               node=tree%baseNode)
          velocityVirial                                (iMass,iOutput)=self%darkMatterHaloScale_               %virialVelocity                 (                                                                                               node=tree%baseNode)
          temperatureVirial                             (iMass,iOutput)=self%darkMatterHaloScale_               %virialTemperature              (                                                                                               node=tree%baseNode)
          radiusVirial                                  (iMass,iOutput)=self%darkMatterHaloScale_               %virialRadius                   (                                                                                               node=tree%baseNode)
          velocityMaximum                               (iMass,iOutput)=self%darkMatterProfileDMO_              %circularVelocityMaximum        (                                                                                               node=tree%baseNode)
          darkMatterProfileRadiusScale                  (iMass,iOutput)=darkMatterProfileHalo                   %scale                          (                                                                                                                 )
          if (self%includeMassAccretionRate) &
               & massAccretionRate                      (iMass,iOutput)=self%darkMatterHaloMassAccretionHistory_%massAccretionRate              (                                                                     time=outputTimes(iOutput),node=tree%baseNode)
          ! Compute alternate mass definitions for halos.
          do iAlternate=1,size(self%virialDensityContrasts)
             massAlternate(iAlternate,iMass,iOutput)=Dark_Matter_Profile_Mass_Definition(tree%baseNode,self%virialDensityContrasts(iAlternate)%virialDensityContrast_%densityContrast(mass=massHalo(iMass),time=outputTimes(iOutput)),radius=radiusAlternate(iAlternate,iMass,iOutput))
          end do
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
    ! Store half- and quarter-mode masses if possible.
    if     (                                             &
         &   statusHalfModeMass    == errorStatusSuccess &
         &  .or.                                         &
         &   statusQuarterModeMass == errorStatusSuccess &
         & ) then
       if (self%outputGroup == ".") then
          powerSpectrumGroup=galacticusOutputFile%openGroup('powerSpectrum','Group containing data relating to the power spectrum.')
       else
          powerSpectrumGroup=containerGroup      %openGroup('powerSpectrum','Group containing data relating to the power spectrum.')
       end if
       if (statusHalfModeMass    == errorStatusSuccess) call powerSpectrumGroup%writeAttribute(massHalfMode   ,'massHalfMode'   )
       if (statusQuarterModeMass == errorStatusSuccess) call powerSpectrumGroup%writeAttribute(massQuarterMode,'massQuarterMode')
       if (size(self%fractionModeMasses) > 0) then
          do iMass=1,size(self%fractionModeMasses)
             write (label,'(a,e8.2)') 'massFractionMode_',self%fractionModeMasses(iMass)
             if (statusFractionMode(iMass) == errorStatusSuccess) call powerSpectrumGroup%writeAttribute(massFractionMode(iMass),trim(label))
          end do          
       end if
       call                                                  powerSpectrumGroup%close         (                                 )
    end if
    ! Store sigma8.
    if (self%outputGroup == ".") then
       powerSpectrumGroup=galacticusOutputFile%openGroup('powerSpectrum','Group containing data relating to the power spectrum.')
    else
       powerSpectrumGroup=containerGroup      %openGroup('powerSpectrum','Group containing data relating to the power spectrum.')
    end if
    call powerSpectrumGroup%writeAttribute(self%cosmologicalMassVariance_%sigma8(),'sigma8')
    call powerSpectrumGroup%close         (                                                )
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
      subhaloMassFunctionIntegrand=+                                                                                                               mass                     &
           &                       *self%haloMassFunction_            %differential(outputTimes(iOutput)                                          ,mass,node=tree%baseNode) &
           &                       *self%unevolvedSubhaloMassFunction_%integrated  (outputTimes(iOutput),massHalo(iMass),haloMassEffectiveInfinity,mass                   )
      return
    end function subhaloMassFunctionIntegrand

  end subroutine haloMassFunctionPerform
