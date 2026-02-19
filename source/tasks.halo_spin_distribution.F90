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

  use :: Cosmology_Functions     , only : cosmologyFunctions      , cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Halo_Spin_Distributions , only : haloSpinDistribution    , haloSpinDistributionClass
  use :: Output_Times            , only : outputTimes             , outputTimesClass

  !![
  <task name="taskHaloSpinDistribution">
   <description>A task which computes and outputs the halo spin distribution.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskHaloSpinDistribution
     !!{
     Implementation of a task which computes and outputs the halo spin distribution.
     !!}
     private
     class           (haloSpinDistributionClass), pointer :: haloSpinDistribution_     => null()
     class           (outputTimesClass         ), pointer :: outputTimes_              => null()
     class           (cosmologyFunctionsClass  ), pointer :: cosmologyFunctions_       => null()
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_      => null()
     double precision                                     :: spinMinimum                         , spinMaximum    , &
          &                                                  spinPointsPerDecade                 , haloMassMinimum
     type            (varying_string           )          :: outputGroup
     ! Pointer to the parameters for this task.
     type            (inputParameters          )          :: parameters
     logical                                              :: nodeComponentsInitialized =  .false.
   contains
     final     ::            haloSpinDistributionDestructor
     procedure :: perform => haloSpinDistributionPerform
  end type taskHaloSpinDistribution

  interface taskHaloSpinDistribution
     !!{
     Constructors for the \refClass{taskHaloSpinDistribution} task.
     !!}
     module procedure haloSpinDistributionConstructorParameters
     module procedure haloSpinDistributionConstructorInternal
  end interface taskHaloSpinDistribution

contains

  function haloSpinDistributionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskHaloSpinDistribution} task class which takes a parameter set as input.
    !!}
    use :: Galacticus_Nodes, only : nodeClassHierarchyInitialize, treeNode
    use :: Input_Parameters, only : inputParameter              , inputParameters
    use :: Node_Components , only : Node_Components_Initialize
    implicit none
    type            (taskHaloSpinDistribution )                        :: self
    type            (inputParameters          ), intent(inout), target :: parameters
    class           (haloSpinDistributionClass), pointer               :: haloSpinDistribution_
    class           (outputTimesClass         ), pointer               :: outputTimes_
    class           (cosmologyFunctionsClass  ), pointer               :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass ), pointer               :: darkMatterHaloScale_
    type            (inputParameters          ), pointer               :: parametersRoot
    type            (varying_string           )                        :: outputGroup
    double precision                                                   :: spinMinimum          , spinMaximum    , &
         &                                                                spinPointsPerDecade  , haloMassMinimum

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
      <name>spinMinimum</name>
      <variable>spinMinimum</variable>
      <defaultValue>3.0d-4</defaultValue>
      <description>Minimum spin for which the distribution function should be calculated.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>spinMaximum</name>
      <variable>spinMaximum</variable>
      <defaultValue>0.5d0</defaultValue>
      <description>Maximum spin for which the distribution function should be calculated.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>spinPointsPerDecade</name>
      <variable>spinPointsPerDecade</variable>
      <defaultValue>10.0d0</defaultValue>
      <description>Number of points per decade of spin at which to calculate the distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>haloMassMinimum</name>
      <variable>haloMassMinimum</variable>
      <defaultValue>0.0d0</defaultValue>
      <description>Minimum halo mass above which spin distribution should be averaged.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>outputGroup</name>
      <defaultValue>var_str('.')</defaultValue>
      <description>The HDF5 output group within which to write spin distribution data.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="haloSpinDistribution" name="haloSpinDistribution_" source="parameters"/>
    <objectBuilder class="outputTimes"          name="outputTimes_"          source="parameters"/>
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !!]
    self=taskHaloSpinDistribution(spinMinimum,spinMaximum,spinPointsPerDecade,haloMassMinimum,outputGroup,darkMatterHaloScale_,haloSpinDistribution_,outputTimes_,cosmologyFunctions_,parametersRoot)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="haloSpinDistribution_"/>
    <objectDestructor name="outputTimes_"         />
    <objectDestructor name="cosmologyFunctions_"  />
    <objectDestructor name="darkMatterHaloScale_" />
    !!]
    return
  end function haloSpinDistributionConstructorParameters

  function haloSpinDistributionConstructorInternal(spinMinimum,spinMaximum,spinPointsPerDecade,haloMassMinimum,outputGroup,darkMatterHaloScale_,haloSpinDistribution_,outputTimes_,cosmologyFunctions_,parameters) result(self)
    !!{
    Constructor for the \refClass{taskHaloSpinDistribution} task class which takes a parameter set as input.
    !!}
    implicit none
    type            (taskHaloSpinDistribution )                        :: self
    class           (haloSpinDistributionClass), intent(in   ), target :: haloSpinDistribution_
    class           (outputTimesClass         ), intent(in   ), target :: outputTimes_
    class           (cosmologyFunctionsClass  ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass ), intent(in   ), target :: darkMatterHaloScale_
    type            (varying_string           ), intent(in   )         :: outputGroup
    double precision                           , intent(in   )         :: spinMinimum          , spinMaximum    , &
         &                                                                spinPointsPerDecade  , haloMassMinimum
    type            (inputParameters          ), intent(in   ), target :: parameters
    !![
    <constructorAssign variables="spinMinimum, spinMaximum, spinPointsPerDecade, haloMassMinimum, outputGroup, *darkMatterHaloScale_, *haloSpinDistribution_, *outputTimes_, *cosmologyFunctions_"/>
    !!]

    self%parameters=inputParameters(parameters)
    call self%parameters%parametersGroupCopy(parameters)
    return
  end function haloSpinDistributionConstructorInternal

  subroutine haloSpinDistributionDestructor(self)
    !!{
    Destructor for the \refClass{taskHaloSpinDistribution} task class.
    !!}
    use :: Node_Components, only : Node_Components_Uninitialize
    implicit none
    type(taskHaloSpinDistribution), intent(inout) :: self

    !![
    <objectDestructor name="self%haloSpinDistribution_"/>
    <objectDestructor name="self%outputTimes_"         />
    <objectDestructor name="self%cosmologyFunctions_"  />
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    if (self%nodeComponentsInitialized) call Node_Components_Uninitialize()
    return
  end subroutine haloSpinDistributionDestructor

  subroutine haloSpinDistributionPerform(self,status)
    !!{
    Compute and output the halo spin distribution.
    !!}
    use            :: Dark_Matter_Halo_Spins , only : Dark_Matter_Halo_Angular_Momentum_Scale
    use            :: Display                , only : displayIndent                          , displayUnindent
    use            :: Error                  , only : Error_Report                           , errorStatusSuccess
    use            :: Output_HDF5            , only : outputFile
    use            :: Galacticus_Nodes       , only : nodeComponentBasic                     , nodeComponentDarkMatterProfile     , nodeComponentSpin, treeNode
    use            :: Halo_Spin_Distributions, only : haloSpinDistributionNbodyErrors
    use            :: IO_HDF5                , only : hdf5Object
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    use            :: Node_Components        , only : Node_Components_Thread_Initialize      , Node_Components_Thread_Uninitialize
    use            :: String_Handling        , only : operator(//)
    implicit none
    class           (taskHaloSpinDistribution      ), intent(inout), target       :: self
    integer                                         , intent(  out), optional     :: status
    type            (treeNode                      ), pointer                     :: node
    class           (nodeComponentBasic            ), pointer                     :: nodeBasic
    class           (nodeComponentSpin             ), pointer                     :: nodeSpin
    class           (nodeComponentDarkMatterProfile), pointer                     :: nodeDarkMatterProfile
    double precision                                , allocatable  , dimension(:) :: spin                 , spinDistribution
    integer         (c_size_t                      )                              :: iOutput
    integer                                                                       :: iSpin                , spinCount
    type            (hdf5Object                    )                              :: outputsGroup         , outputGroup     , &
         &                                                                           containerGroup
    type            (varying_string                )                              :: groupName            , description

    call displayIndent('Begin task: halo spin distribution')
    ! Call routines to perform initialization which must occur for all threads if run in parallel.
    call Node_Components_Thread_Initialize(self%parameters)
    ! Create a tree node.
    node                  => treeNode                  (                 )
    nodeBasic             => node    %basic            (autoCreate=.true.)
    nodeSpin              => node    %spin             (autoCreate=.true.)
    nodeDarkMatterProfile => node    %darkMatterProfile(autoCreate=.true.)
    ! Compute number of spin points to tabulate.
    spinCount=int(log10(self%spinMaximum/self%spinMinimum)*self%spinPointsPerDecade)+1
    allocate(spin            (spinCount))
    allocate(spinDistribution(spinCount))
    ! Open the group for output time information.
    if (self%outputGroup == ".") then
       outputsGroup  =outputFile    %openGroup(     'Outputs'        ,'Group containing datasets relating to output times.')
    else
       containerGroup=outputFile    %openGroup(char(self%outputGroup),'Group containing halo mass function data.'          )
       outputsGroup  =containerGroup%openGroup(     'Outputs'        ,'Group containing datasets relating to output times.')
    end if
    ! Iterate over output redshifts.
    do iOutput=1,self%outputTimes_%count()
       ! Set the epoch for the node.
       call nodeBasic%timeSet(self%outputTimes_%time(iOutput))
       ! Iterate over spins.
       do iSpin=1,spinCount
          spin(iSpin)=exp(log(self%spinMinimum)+log(self%spinMaximum/self%spinMinimum)*dble(iSpin-1)/dble(spinCount-1))
          call nodeSpin%angularMomentumSet(spin(iSpin)*Dark_Matter_Halo_Angular_Momentum_Scale(node,self%darkMatterHaloScale_))
          ! Evaluate the distribution.
          if (self%haloMassMinimum <= 0.0d0) then
             ! No minimum halo mass specified - simply evaluate the spin distribution.
             spinDistribution(iSpin)=self%haloSpinDistribution_%distribution(node)
          else
             ! A minimum halo mass was specified. Evaluate the spin distribution averaged over halo mass, if the distribution class
             ! supports this.
             select type (haloSpinDistribution_ => self%haloSpinDistribution_)
             class is (haloSpinDistributionNbodyErrors)
                spinDistribution(iSpin)=haloSpinDistribution_%distributionAveraged(node,self%haloMassMinimum)
             class default
                call Error_Report('halo spin distribution class does not support averaging over halo mass'//{introspection:location})
             end select
          end if
       end do
       ! Open the output group.
       groupName  ='Output'
       description='Data for output number '
       groupName  =groupName  //iOutput
       description=description//iOutput
       outputGroup=outputsGroup%openGroup(char(groupName),char(description))
       ! Store the distribution, redshifts, and spins.
       call outputGroup%writeAttribute(                                                                                              self%outputTimes_%time(iOutput)  ,'outputTime'                                                                )
       call outputGroup%writeAttribute(                                                     self%cosmologyFunctions_%expansionFactor(self%outputTimes_%time(iOutput)) ,'outputExpansionFactor'                                                     )
       call outputGroup%writeAttribute(self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(self%outputTimes_%time(iOutput))),'outputRedshift'                                                            )
       call outputGroup%writeDataset  (spin                                                                                                                           ,'spin'                 ,'Spins at which the spin distribution is tabulated.')
       call outputGroup%writeDataset  (spinDistribution                                                                                                               ,'spinDistribution'     ,'Spin parameter distribution.'                      )
       call outputGroup%close         (                                                                                                                                                                                                            )
    end do
    call outputsGroup%close()
    if (containerGroup%isOpen()) call containerGroup%close()
    call Node_Components_Thread_Uninitialize()
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: halo spin distribution' )
    return
  end subroutine haloSpinDistributionPerform
