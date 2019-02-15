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
  
  !# <task name="taskHaloModelProjectedCorrelationFunction">
  !#  <description>A task which generates a mock catalog of galaxies based on a simple halo model approach.</description>
  !# </task>
  type, extends(taskClass) :: taskHaloModelProjectedCorrelationFunction
     !% Implementation of a task which generates a mock catalog of galaxies based on a simple halo model approach.
     private 
     class           (conditionalMassFunctionClass), pointer                   :: conditionalMassFunction_ => null()
     double precision                              , allocatable, dimension(:) :: separationProjectedBinned, correlationProjectedBinned
     double precision                                                          :: separationMinimum        , separationMaximum         , &
          &                                                                       massMinimum              , massMaximum               , &
          &                                                                       massHaloMinimum          , massHaloMaximum           , &
          &                                                                       depthLineOfSight
     integer                                                                   :: countSeparations
     logical                                                                   :: halfIntegral
     type            (varying_string              )                            :: outputGroup
   contains
     final     ::                       haloModelProjectedCorrelationFunctionDestructor
     procedure :: perform            => haloModelProjectedCorrelationFunctionPerform
     procedure :: requiresOutputFile => haloModelProjectedCorrelationFunctionRequiresOutputFile
  end type taskHaloModelProjectedCorrelationFunction

  interface taskHaloModelProjectedCorrelationFunction
     !% Constructors for the {\normalfont \ttfamily haloModelProjectedCorrelationFunction} task.
     module procedure haloModelProjectedCorrelationFunctionConstructorParameters
     module procedure haloModelProjectedCorrelationFunctionConstructorInternal
  end interface taskHaloModelProjectedCorrelationFunction

contains

  function haloModelProjectedCorrelationFunctionConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily haloModelProjectedCorrelationFunction} task class which takes a parameter set as input.
    use Input_Parameters, only : inputParameters             , inputParameter
    use Galacticus_Nodes, only : nodeClassHierarchyInitialize
    use Node_Components , only : Node_Components_Initialize  , Node_Components_Thread_Initialize
    implicit none
    type            (taskHaloModelProjectedCorrelationFunction)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (conditionalMassFunctionClass             ), pointer       :: conditionalMassFunction_
    double precision                                                           :: separationMinimum       , separationMaximum, &
         &                                                                        massMinimum             , massMaximum      , &
         &                                                                        massHaloMinimum         , massHaloMaximum  , &
         &                                                                        depthLineOfSight
    integer                                                                    :: countSeparations
    logical                                                                    :: halfIntegral
    type            (varying_string                           )                :: outputGroup

    call nodeClassHierarchyInitialize     (parameters)
    call Node_Components_Initialize       (parameters)
    call Node_Components_Thread_Initialize(parameters)
    !# <inputParameter>
    !#   <name>separationMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The minimum separation at which to compute the projected correlation function.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>separationMaximum</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The maximum separation at which to compute the projected correlation function.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>countSeparations</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The number of separations at which to compute the projected correlation function.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>depthLineOfSight</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The maximum line of sight depth to which to integrate when computing the projected correlation function.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>halfIntegral</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>Set to {\normalfont \ttfamily true} if the projected correlation function is computed as $w_\mathrm{p}(r_\mathrm{p})=\int_0^{+\pi_\mathrm{max}} \xi(r_\mathrm{p},\pi) \mathrm{d} \pi$, instead of the usual $w_\mathrm{p}(r_\mathrm{p})=\int_{-\pi_\mathrm{max}}^{+\pi_\mathrm{max}} \xi(r_\mathrm{p},\pi) \mathrm{d} \pi$.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d8</defaultValue>
    !#   <description>The minimum mass of galaxies to include in the projected correlation function calculation.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massMaximum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d12</defaultValue>
    !#   <description>The maximum mass of galaxies to include in the projected correlation function calculation.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massHaloMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d6</defaultValue>
    !#   <description>The minimum halo mass to use when integrating over the halo mass function.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massHaloMaximum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d16</defaultValue>
    !#   <description>The maximum halo mass to use when integrating over the halo mass function.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>outputGroup</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>var_str('projectedCorrelationFunction')</defaultValue>
    !#   <description>The HDF5 output group within which to write the projected correlation function.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>    
    !# <objectBuilder class="conditionalMassFunction" name="conditionalMassFunction_" source="parameters"/>
    self=taskHaloModelProjectedCorrelationFunction(separationMinimum,separationMaximum,countSeparations,massMinimum,massMaximum,massHaloMinimum,massHaloMaximum,depthLineOfSight,halfIntegral,outputGroup,conditionalMassFunction_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="conditionalMassFunction_"/>
    return
  end function haloModelProjectedCorrelationFunctionConstructorParameters

  function haloModelProjectedCorrelationFunctionConstructorInternal(separationMinimum,separationMaximum,countSeparations,massMinimum,massMaximum,massHaloMinimum,massHaloMaximum,depthLineOfSight,halfIntegral,outputGroup,conditionalMassFunction_) result(self)
    !% Constructor for the {\normalfont \ttfamily haloModelProjectedCorrelationFunction} task class which takes a parameter set as input.
    use Numerical_Ranges , only : Make_Range   , rangeTypeLogarithmic
    use Memory_Management, only : allocateArray
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
    !# <constructorAssign variables="separationMinimum, separationMaximum, massMinimum, massMaximum, massHaloMinimum, massHaloMaximum, depthLineOfSight, countSeparations, halfIntegral, outputGroup, *conditionalMassFunction_"/>

    call allocateArray(self%separationProjectedBinned ,[self%countSeparations])
    call allocateArray(self%correlationProjectedBinned,[self%countSeparations])
    self%separationProjectedBinned=Make_Range(self%separationMinimum,self%separationMaximum,self%countSeparations,rangeTypeLogarithmic)
    return
  end function haloModelProjectedCorrelationFunctionConstructorInternal

  subroutine haloModelProjectedCorrelationFunctionDestructor(self)
    !% Destructor for the {\normalfont \ttfamily haloModelProjectedCorrelationFunction} task class.
    use Node_Components, only : Node_Components_Uninitialize, Node_Components_Thread_Uninitialize
    implicit none
    type(taskHaloModelProjectedCorrelationFunction), intent(inout) :: self

    !# <objectDestructor name="self%conditionalMassFunction_"/>
    call Node_Components_Uninitialize       ()
    call Node_Components_Thread_Uninitialize()
    return
  end subroutine haloModelProjectedCorrelationFunctionDestructor
  
  subroutine haloModelProjectedCorrelationFunctionPerform(self,status)
    !% Generate a mock galaxy catalog using a simple halo model approach.   
    use IO_HDF5                          , only : hdf5Object
    use Galacticus_HDF5                  , only : galacticusOutputFile
    use Galacticus_Display               , only : Galacticus_Display_Indent       , Galacticus_Display_Unindent
    use Galacticus_Error                 , only : errorStatusSuccess
    use Halo_Model_Projected_Correlations, only : Halo_Model_Projected_Correlation
    implicit none
    class  (taskHaloModelProjectedCorrelationFunction), intent(inout)           :: self
    integer                                           , intent(  out), optional :: status
    type   (hdf5Object                               )                          :: outputGroup

    call Galacticus_Display_Indent('Begin task: halo model projected correlation function')
    call Halo_Model_Projected_Correlation(                                 &
         &                                self%conditionalMassFunction_  , &
         &                                self%separationProjectedBinned , &
         &                                self%massMinimum               , &
         &                                self%massMaximum               , &
         &                                self%massHaloMinimum           , &
         &                                self%massHaloMaximum           , &
         &                                self%depthLineOfSight          , &
         &                                self%halfIntegral              , &
         &                                self%correlationProjectedBinned  &
         &                               )
    outputGroup=galacticusOutputFile%openGroup(char(self%outputGroup),'Group containing halo mass function data.'          )
    call outputGroup%writeDataset(self%separationProjectedBinned ,"separation"          ,commentText="Projected separation [Mpc]." )
    call outputGroup%writeDataset(self%correlationProjectedBinned,"projectedCorrelation",commentText="Projected correlation [Mpc].")
    call outputGroup%close       (                                                                                                 )
    if (present(status)) status=errorStatusSuccess
    call Galacticus_Display_Unindent('Done task: halo model projected correlation function' )
  end subroutine haloModelProjectedCorrelationFunctionPerform

  logical function haloModelProjectedCorrelationFunctionRequiresOutputFile(self)
    !% Specifies that this task does not requires the main output file.
    implicit none
    class(taskHaloModelProjectedCorrelationFunction), intent(inout) :: self    
    !GCC$ attributes unused :: self

    haloModelProjectedCorrelationFunctionRequiresOutputFile=.true.
    return
  end function haloModelProjectedCorrelationFunctionRequiresOutputFile
