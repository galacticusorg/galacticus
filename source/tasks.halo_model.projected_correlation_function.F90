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

  use :: Conditional_Mass_Functions, only : conditionalMassFunction     , conditionalMassFunctionClass
  use :: Cosmology_Functions       , only : cosmologyFunctions          , cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Biases   , only : darkMatterHaloBias          , darkMatterHaloBiasClass
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScale         , darkMatterHaloScaleClass
  use :: Dark_Matter_Profile_Scales, only : darkMatterProfileScaleRadius, darkMatterProfileScaleRadiusClass
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMO        , darkMatterProfileDMOClass
  use :: Geometry_Surveys          , only : surveyGeometry              , surveyGeometryClass
  use :: Halo_Mass_Functions       , only : haloMassFunction            , haloMassFunctionClass
  use :: Linear_Growth             , only : linearGrowth                , linearGrowthClass
  use :: Power_Spectra             , only : powerSpectrum               , powerSpectrumClass

  !![
  <task name="taskHaloModelProjectedCorrelationFunction">
   <description>A task which generates a mock catalog of galaxies based on a simple halo model approach.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskHaloModelProjectedCorrelationFunction
     !!{
     Implementation of a task which generates a mock catalog of galaxies based on a simple halo model approach.
     !!}
     private
     class           (conditionalMassFunctionClass     ), pointer                   :: conditionalMassFunction_      => null()
     class           (powerSpectrumClass               ), pointer                   :: powerSpectrum_                => null()
     class           (cosmologyFunctionsClass          ), pointer                   :: cosmologyFunctions_           => null()
     class           (surveyGeometryClass              ), pointer                   :: surveyGeometry_               => null()
     class           (darkMatterHaloScaleClass         ), pointer                   :: darkMatterHaloScale_          => null()
     class           (haloMassFunctionClass            ), pointer                   :: haloMassFunction_             => null()
     class           (darkMatterProfileDMOClass        ), pointer                   :: darkMatterProfileDMO_         => null()
     class           (darkMatterHaloBiasClass          ), pointer                   :: darkMatterHaloBias_           => null()
     class           (darkMatterProfileScaleRadiusClass), pointer                   :: darkMatterProfileScaleRadius_ => null()
      double precision                                  , allocatable, dimension(:) :: separationProjectedBinned               , correlationProjectedBinned
     double precision                                                               :: separationMinimum                       , separationMaximum         , &
          &                                                                            massMinimum                             , massMaximum               , &
          &                                                                            massHaloMinimum                         , massHaloMaximum           , &
          &                                                                            depthLineOfSight
     integer                                                                        :: countSeparations
     logical                                                                        :: nodeComponentsInitialized     =  .false., halfIntegral
     type            (varying_string                   )                            :: outputGroup
     ! Pointer to the parameters for this task.
     type            (inputParameters                  )                            :: parameters
   contains
     final     ::                       haloModelProjectedCorrelationFunctionDestructor
     procedure :: perform            => haloModelProjectedCorrelationFunctionPerform
     procedure :: requiresOutputFile => haloModelProjectedCorrelationFunctionRequiresOutputFile
  end type taskHaloModelProjectedCorrelationFunction

  interface taskHaloModelProjectedCorrelationFunction
     !!{
     Constructors for the \refClass{taskHaloModelProjectedCorrelationFunction} task.
     !!}
     module procedure haloModelProjectedCorrelationFunctionConstructorParameters
     module procedure haloModelProjectedCorrelationFunctionConstructorInternal
  end interface taskHaloModelProjectedCorrelationFunction

contains

  function haloModelProjectedCorrelationFunctionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskHaloModelProjectedCorrelationFunction} task class which takes a parameter set as input.
    !!}
    use :: Galacticus_Nodes, only : nodeClassHierarchyInitialize
    use :: Input_Parameters, only : inputParameter              , inputParameters
    use :: Node_Components , only : Node_Components_Initialize
    implicit none
    type            (taskHaloModelProjectedCorrelationFunction)                        :: self
    type            (inputParameters                          ), intent(inout), target :: parameters
    class           (conditionalMassFunctionClass             ), pointer               :: conditionalMassFunction_
    class           (powerSpectrumClass                       ), pointer               :: powerSpectrum_
    class           (cosmologyFunctionsClass                  ), pointer               :: cosmologyFunctions_
    class           (surveyGeometryClass                      ), pointer               :: surveyGeometry_
    class           (darkMatterHaloScaleClass                 ), pointer               :: darkMatterHaloScale_
    class           (haloMassFunctionClass                    ), pointer               :: haloMassFunction_
    class           (darkMatterProfileDMOClass                ), pointer               :: darkMatterProfileDMO_
    class           (darkMatterHaloBiasClass                  ), pointer               :: darkMatterHaloBias_
    class           (darkMatterProfileScaleRadiusClass        ), pointer               :: darkMatterProfileScaleRadius_
    type            (inputParameters                          ), pointer               :: parametersRoot
    double precision                                                                   :: separationMinimum            , separationMaximum, &
         &                                                                                massMinimum                  , massMaximum      , &
         &                                                                                massHaloMinimum              , massHaloMaximum  , &
         &                                                                                depthLineOfSight
    integer                                                                            :: countSeparations
    logical                                                                            :: halfIntegral
    type            (varying_string                           )                        :: outputGroup

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
      <name>separationMinimum</name>
      <description>The minimum separation at which to compute the projected correlation function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>separationMaximum</name>
      <description>The maximum separation at which to compute the projected correlation function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>countSeparations</name>
      <description>The number of separations at which to compute the projected correlation function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>depthLineOfSight</name>
      <description>The maximum line of sight depth to which to integrate when computing the projected correlation function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>halfIntegral</name>
      <defaultValue>.false.</defaultValue>
      <description>Set to {\normalfont \ttfamily true} if the projected correlation function is computed as $w_\mathrm{p}(r_\mathrm{p})=\int_0^{+\pi_\mathrm{max}} \xi(r_\mathrm{p},\pi) \mathrm{d} \pi$, instead of the usual $w_\mathrm{p}(r_\mathrm{p})=\int_{-\pi_\mathrm{max}}^{+\pi_\mathrm{max}} \xi(r_\mathrm{p},\pi) \mathrm{d} \pi$.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMinimum</name>
      <defaultValue>1.0d8</defaultValue>
      <description>The minimum mass of galaxies to include in the projected correlation function calculation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <defaultValue>1.0d12</defaultValue>
      <description>The maximum mass of galaxies to include in the projected correlation function calculation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massHaloMinimum</name>
      <defaultValue>1.0d6</defaultValue>
      <description>The minimum halo mass to use when integrating over the halo mass function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massHaloMaximum</name>
      <defaultValue>1.0d16</defaultValue>
      <description>The maximum halo mass to use when integrating over the halo mass function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>outputGroup</name>
      <defaultValue>var_str('projectedCorrelationFunction')</defaultValue>
      <description>The HDF5 output group within which to write the projected correlation function.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="conditionalMassFunction"      name="conditionalMassFunction_"      source="parameters"/>
    <objectBuilder class="powerSpectrum"                name="powerSpectrum_"                source="parameters"/>
    <objectBuilder class="cosmologyFunctions"           name="cosmologyFunctions_"           source="parameters"/>
    <objectBuilder class="surveyGeometry"               name="surveyGeometry_"               source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters"/>
    <objectBuilder class="haloMassFunction"             name="haloMassFunction_"             source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"         name="darkMatterProfileDMO_"         source="parameters"/>
    <objectBuilder class="darkMatterHaloBias"           name="darkMatterHaloBias_"           source="parameters"/>
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    !!]
    self=taskHaloModelProjectedCorrelationFunction(separationMinimum,separationMaximum,countSeparations,massMinimum,massMaximum,massHaloMinimum,massHaloMaximum,depthLineOfSight,halfIntegral,outputGroup,conditionalMassFunction_,powerSpectrum_,cosmologyFunctions_,surveyGeometry_,darkMatterHaloScale_,haloMassFunction_,darkMatterProfileDMO_,darkMatterHaloBias_,darkMatterProfileScaleRadius_,parametersRoot)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="conditionalMassFunction_"     />
    <objectDestructor name="powerSpectrum_"               />
    <objectDestructor name="cosmologyFunctions_"          />
    <objectDestructor name="surveyGeometry_"              />
    <objectDestructor name="darkMatterHaloScale_"         />
    <objectDestructor name="haloMassFunction_"            />
    <objectDestructor name="darkMatterProfileDMO_"        />
    <objectDestructor name="darkMatterHaloBias_"          />
    <objectDestructor name="darkMatterProfileScaleRadius_"/>
    !!]
    return
  end function haloModelProjectedCorrelationFunctionConstructorParameters

  function haloModelProjectedCorrelationFunctionConstructorInternal(separationMinimum,separationMaximum,countSeparations,massMinimum,massMaximum,massHaloMinimum,massHaloMaximum,depthLineOfSight,halfIntegral,outputGroup,conditionalMassFunction_,powerSpectrum_,cosmologyFunctions_,surveyGeometry_,darkMatterHaloScale_,haloMassFunction_,darkMatterProfileDMO_,darkMatterHaloBias_,darkMatterProfileScaleRadius_,parameters) result(self)
    !!{
    Constructor for the \refClass{taskHaloModelProjectedCorrelationFunction} task class which takes a parameter set as input.
    !!}
    use :: Numerical_Ranges , only : Make_Range   , rangeTypeLogarithmic
    implicit none
    type            (taskHaloModelProjectedCorrelationFunction)                        :: self
    double precision                                           , intent(in   )         :: separationMinimum       , separationMaximum, &
         &                                                                                massMinimum             , massMaximum      , &
         &                                                                                massHaloMinimum         , massHaloMaximum  , &
         &                                                                                depthLineOfSight
    integer                                                    , intent(in   )         :: countSeparations
    logical                                                    , intent(in   )         :: halfIntegral
    type            (varying_string                           ), intent(in   )         :: outputGroup
    class           (conditionalMassFunctionClass             ), intent(in   ), target :: conditionalMassFunction_
    class           (powerSpectrumClass                       ), intent(in   ), target :: powerSpectrum_
    class           (cosmologyFunctionsClass                  ), intent(in   ), target :: cosmologyFunctions_
    class           (surveyGeometryClass                      ), intent(in   ), target :: surveyGeometry_
    class           (darkMatterHaloScaleClass                 ), intent(in   ), target :: darkMatterHaloScale_
    class           (haloMassFunctionClass                    ), intent(in   ), target :: haloMassFunction_
    class           (darkMatterProfileDMOClass                ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloBiasClass                  ), intent(in   ), target :: darkMatterHaloBias_
    class           (darkMatterProfileScaleRadiusClass        ), intent(in   ), target :: darkMatterProfileScaleRadius_
    type            (inputParameters                          ), intent(in   ), target :: parameters
    !![
    <constructorAssign variables="separationMinimum, separationMaximum, massMinimum, massMaximum, massHaloMinimum, massHaloMaximum, depthLineOfSight, countSeparations, halfIntegral, outputGroup, *conditionalMassFunction_, *powerSpectrum_, *cosmologyFunctions_, *surveyGeometry_, *darkMatterHaloScale_, *haloMassFunction_, *darkMatterProfileDMO_, *darkMatterHaloBias_, *darkMatterProfileScaleRadius_"/>
    !!]

    self%parameters=inputParameters(parameters)
    call self%parameters%parametersGroupCopy(parameters)
    allocate(self%separationProjectedBinned (self%countSeparations))
    allocate(self%correlationProjectedBinned(self%countSeparations))
    self%separationProjectedBinned=Make_Range(self%separationMinimum,self%separationMaximum,self%countSeparations,rangeTypeLogarithmic)
    return
  end function haloModelProjectedCorrelationFunctionConstructorInternal

  subroutine haloModelProjectedCorrelationFunctionDestructor(self)
    !!{
    Destructor for the \refClass{taskHaloModelProjectedCorrelationFunction} task class.
    !!}
    use :: Node_Components, only : Node_Components_Uninitialize
    implicit none
    type(taskHaloModelProjectedCorrelationFunction), intent(inout) :: self

    !![
    <objectDestructor name="self%conditionalMassFunction_"     />
    <objectDestructor name="self%powerSpectrum_"               />
    <objectDestructor name="self%cosmologyFunctions_"          />
    <objectDestructor name="self%surveyGeometry_"              />
    <objectDestructor name="self%darkMatterHaloScale_"         />
    <objectDestructor name="self%haloMassFunction_"            />
    <objectDestructor name="self%darkMatterProfileDMO_"        />
    <objectDestructor name="self%darkMatterHaloBias_"          />
    <objectDestructor name="self%darkMatterProfileScaleRadius_"/>
    !!]
    if (self%nodeComponentsInitialized) call Node_Components_Uninitialize()
    return
  end subroutine haloModelProjectedCorrelationFunctionDestructor

  subroutine haloModelProjectedCorrelationFunctionPerform(self,status)
    !!{
    Generate a mock galaxy catalog using a simple halo model approach.
    !!}
    use :: Display                          , only : displayIndent                    , displayUnindent
    use :: Error                            , only : errorStatusSuccess
    use :: Output_HDF5                      , only : outputFile
    use :: Halo_Model_Projected_Correlations, only : Halo_Model_Projected_Correlation
    use :: IO_HDF5                          , only : hdf5Object
    use :: Node_Components                  , only : Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize
    implicit none
    class  (taskHaloModelProjectedCorrelationFunction), intent(inout), target   :: self
    integer                                           , intent(  out), optional :: status
    type   (hdf5Object                               )                          :: outputGroup

    call displayIndent('Begin task: halo model projected correlation function')
    ! Call routines to perform initialization which must occur for all threads if run in parallel.
    call Node_Components_Thread_Initialize(self%parameters)
    call Halo_Model_Projected_Correlation(                                    &
         &                                self%conditionalMassFunction_     , &
         &                                self%powerSpectrum_               , &
         &                                self%cosmologyFunctions_          , &
         &                                self%surveyGeometry_              , &
         &                                self%darkMatterHaloScale_         , &
         &                                self%haloMassFunction_            , &
         &                                self%darkMatterProfileDMO_        , &
         &                                self%darkMatterHaloBias_          , &
         &                                self%darkMatterProfileScaleRadius_, &
         &                                self%separationProjectedBinned    , &
         &                                self%massMinimum                  , &
         &                                self%massMaximum                  , &
         &                                self%massHaloMinimum              , &
         &                                self%massHaloMaximum              , &
         &                                self%depthLineOfSight             , &
         &                                self%halfIntegral                 , &
         &                                self%correlationProjectedBinned     &
         &                               )
    outputGroup=outputFile%openGroup(char(self%outputGroup),'Group containing halo mass function data.')
    call outputGroup%writeDataset(self%separationProjectedBinned ,"separation"          ,comment="Projected separation [Mpc]." )
    call outputGroup%writeDataset(self%correlationProjectedBinned,"projectedCorrelation",comment="Projected correlation [Mpc].")
    call outputGroup%close       (                                                                                                 )
    call Node_Components_Thread_Uninitialize()
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: halo model projected correlation function' )
  end subroutine haloModelProjectedCorrelationFunctionPerform

  logical function haloModelProjectedCorrelationFunctionRequiresOutputFile(self)
    !!{
    Specifies that this task does not requires the main output file.
    !!}
    implicit none
    class(taskHaloModelProjectedCorrelationFunction), intent(inout) :: self
    !$GLC attributes unused :: self

    haloModelProjectedCorrelationFunctionRequiresOutputFile=.true.
    return
  end function haloModelProjectedCorrelationFunctionRequiresOutputFile
