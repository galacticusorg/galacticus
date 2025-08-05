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

  use :: Cosmology_Functions                  , only : cosmologyFunctions              , cosmologyFunctionsClass
  use :: Intergalactic_Medium_Filtering_Masses, only : intergalacticMediumFilteringMass, intergalacticMediumFilteringMassClass
  use :: Intergalactic_Medium_State           , only : intergalacticMediumState        , intergalacticMediumStateClass
  use :: Output_Times                         , only : outputTimes                     , outputTimesClass

  !![
  <task name="taskIntergalacticMediumState">
   <description>A task which outputs the state of the intergalactic medium.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskIntergalacticMediumState
     !!{
     Implementation of a task which computes and outputs the halo mass function and related quantities.
     !!}
     private
     class(cosmologyFunctionsClass              ), pointer :: cosmologyFunctions_               => null()
     class(outputTimesClass                     ), pointer :: outputTimes_                      => null()
     class(intergalacticMediumStateClass        ), pointer :: intergalacticMediumState_         => null()
     class(intergalacticMediumFilteringMassClass), pointer :: intergalacticMediumFilteringMass_ => null()
     type (varying_string                       )          :: outputGroup
   contains
     final     ::            intergalacticMediumStateDestructor
     procedure :: perform => intergalacticMediumStatePerform
  end type taskIntergalacticMediumState

  interface taskIntergalacticMediumState
     !!{
     Constructors for the \refClass{taskIntergalacticMediumState} task.
     !!}
     module procedure intergalacticMediumStateConstructorParameters
     module procedure intergalacticMediumStateConstructorInternal
  end interface taskIntergalacticMediumState

contains

  function intergalacticMediumStateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskIntergalacticMediumState} task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (taskIntergalacticMediumState         )                :: self
    type (inputParameters                      ), intent(inout) :: parameters
    class(cosmologyFunctionsClass              ), pointer       :: cosmologyFunctions_
    class(outputTimesClass                     ), pointer       :: outputTimes_
    class(intergalacticMediumStateClass        ), pointer       :: intergalacticMediumState_
    class(intergalacticMediumFilteringMassClass), pointer       :: intergalacticMediumFilteringMass_
    type (varying_string                       )                :: outputGroup

    !![
    <inputParameter>
      <name>outputGroup</name>
      <defaultValue>var_str('.')</defaultValue>
      <description>The HDF5 output group within which to write intergalactic medium state data.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"               name="cosmologyFunctions_"               source="parameters"/>
    <objectBuilder class="outputTimes"                      name="outputTimes_"                      source="parameters"/>
    <objectBuilder class="intergalacticMediumState"         name="intergalacticMediumState_"         source="parameters"/>
    <objectBuilder class="intergalacticMediumFilteringMass" name="intergalacticMediumFilteringMass_" source="parameters"/>
    !!]
    self=taskIntergalacticMediumState(                                   &
         &                            outputGroup                      , &
         &                            cosmologyFunctions_              , &
         &                            outputTimes_                     , &
         &                            intergalacticMediumState_        , &
         &                            intergalacticMediumFilteringMass_  &
         &                           )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"              />
    <objectDestructor name="outputTimes_"                     />
    <objectDestructor name="intergalacticMediumState_"        />
    <objectDestructor name="intergalacticMediumFilteringMass_"/>
    !!]
    return
  end function intergalacticMediumStateConstructorParameters

  function intergalacticMediumStateConstructorInternal(                                   &
       &                                               outputGroup                      , &
       &                                               cosmologyFunctions_              , &
       &                                               outputTimes_                     , &
       &                                               intergalacticMediumState_        , &
       &                                               intergalacticMediumFilteringMass_  &
       &                                              ) result(self)
    !!{
    Constructor for the \refClass{taskIntergalacticMediumState} task class which takes a parameter set as input.
    !!}
    implicit none
    type (taskIntergalacticMediumState         )                        :: self
    class(cosmologyFunctionsClass              ), intent(in   ), target :: cosmologyFunctions_
    class(outputTimesClass                     ), intent(in   ), target :: outputTimes_
    class(intergalacticMediumStateClass        ), intent(in   ), target :: intergalacticMediumState_
    class(intergalacticMediumFilteringMassClass), intent(in   ), target :: intergalacticMediumFilteringMass_
    type (varying_string                       ), intent(in   )         :: outputGroup
    !![
    <constructorAssign variables="outputGroup, *cosmologyFunctions_, *outputTimes_, *intergalacticMediumState_, *intergalacticMediumFilteringMass_"/>
    !!]

    return
  end function intergalacticMediumStateConstructorInternal

  subroutine intergalacticMediumStateDestructor(self)
    !!{
    Destructor for the \refClass{taskIntergalacticMediumState} task class.
    !!}
    implicit none
    type(taskIntergalacticMediumState), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"              />
    <objectDestructor name="self%outputTimes_"                     />
    <objectDestructor name="self%intergalacticMediumState_"        />
    <objectDestructor name="self%intergalacticMediumFilteringMass_"/>
    !!]
    return
  end subroutine intergalacticMediumStateDestructor

  subroutine intergalacticMediumStatePerform(self,status)
    !!{
    Output \gls{igm} state to the \glc\ output file.
    !!}
    use            :: Display         , only : displayIndent     , displayUnindent
    use            :: Error           , only : errorStatusSuccess
    use            :: Output_HDF5     , only : outputFile
    use            :: HDF5_Access     , only : hdf5Access
    use            :: IO_HDF5         , only : hdf5Object
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    use            :: String_Handling , only : operator(//)
    implicit none
    class           (taskIntergalacticMediumState), intent(inout), target       :: self
    integer                                       , intent(  out), optional     :: status
    double precision                              , allocatable  , dimension(:) :: time           , redshift     , &
         &                                                                         temperature    , massFiltering, &
         &                                                                         expansionFactor, massJeans
    integer         (c_size_t                    )                              :: outputCount    , i
    type            (hdf5Object                  )                              :: outputsGroup   , outputGroup  , &
         &                                                                         containerGroup
    type            (varying_string              )                              :: groupName      , description

    call displayIndent('Begin task: intergalactic medium state')
    ! Get the requested output redshifts.
    outputCount=self%outputTimes_%count()
    allocate(time           (outputCount))
    allocate(redshift       (outputCount))
    allocate(expansionFactor(outputCount))
    allocate(temperature    (outputCount))
    allocate(massFiltering  (outputCount))
    allocate(massJeans      (outputCount))
    ! Populate arrays with IGM state.
    do i=1,outputCount
       time           (i)=self%outputTimes_                     %time                                (i)
       redshift       (i)=self%outputTimes_                     %redshift                            (i)
       expansionFactor(i)=self%cosmologyFunctions_              %expansionFactorFromRedshift(redshift(i))
       temperature    (i)=self%intergalacticMediumState_        %temperature                (time    (i))
       massJeans      (i)=self%intergalacticMediumState_        %massJeans                  (time    (i))
       massFiltering  (i)=self%intergalacticMediumFilteringMass_%massFiltering              (time    (i))
    end do
    !$ call hdf5Access%set()
    ! Open the group for output time information.
    if (self%outputGroup == ".") then
       outputsGroup  =outputFile    %openGroup(     'Outputs'        ,'Group containing datasets relating to output times.')
    else
       containerGroup=outputFile    %openGroup(char(self%outputGroup),'Group containing intergalactic medium state data.'  )
       outputsGroup  =containerGroup%openGroup(     'Outputs'        ,'Group containing datasets relating to output times.')
    end if
    do i=1,outputCount
       groupName  ='Output'
       description='Data for output number '
       groupName  =groupName  //i
       description=description//i
       outputGroup=outputsGroup%openGroup(char(groupName),char(description))
       call outputGroup%writeAttribute(time           (i),'outputTime'           )
       call outputGroup%writeAttribute(redshift       (i),'outputRedshift'       )
       call outputGroup%writeAttribute(expansionFactor(i),'outputExpansionFactor')
       call outputGroup%writeAttribute(temperature    (i),'temperature'          )
       call outputGroup%writeAttribute(massJeans      (i),'massJeans'            )
       call outputGroup%writeAttribute(massFiltering  (i),'massFiltering'        )
    end do
    !$ call hdf5Access%unset()
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: intergalactic medium state' )
    return
  end subroutine intergalacticMediumStatePerform
