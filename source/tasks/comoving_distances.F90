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

  use :: Cosmology_Functions, only : cosmologyFunctionsClass
  use :: Output_Times       , only : outputTimesClass

  !![
  <task name="taskComovingDistances" docformat="rst">
   <description>
   A task which computes and outputs the comoving distance to each output.
   </description>
  </task>
  !!]
  type, extends(taskClass) :: taskComovingDistances
     !!{RST
     Implementation of a task which computes and outputs the comoving distance to each output.
     !!}
     private
     class(cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     class(outputTimesClass       ), pointer :: outputTimes_        => null()
   contains
     final     ::            comovingDistancesDestructor
     procedure :: perform => comovingDistancesPerform
  end type taskComovingDistances

  interface taskComovingDistances
     !!{RST
     Constructors for the :galacticus-class:`taskComovingDistances` task.
     !!}
     module procedure comovingDistancesConstructorParameters
     module procedure comovingDistancesConstructorInternal
  end interface taskComovingDistances

contains

  function comovingDistancesConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`taskComovingDistances` task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (taskComovingDistances  )                :: self
    type (inputParameters        ), intent(inout) :: parameters
    class(cosmologyFunctionsClass), pointer       :: cosmologyFunctions_
    class(outputTimesClass       ), pointer       :: outputTimes_
    
    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <objectBuilder class="outputTimes"        name="outputTimes_"        source="parameters"/>
    !!]
    self=taskComovingDistances(cosmologyFunctions_,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="outputTimes_"       />
    !!]
    return
  end function comovingDistancesConstructorParameters

  function comovingDistancesConstructorInternal(cosmologyFunctions_,outputTimes_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`taskComovingDistances` task class.
    !!}
    implicit none
    type (taskComovingDistances  )                        :: self
    class(cosmologyFunctionsClass), intent(in   ), target :: cosmologyFunctions_
    class(outputTimesClass       ), intent(in   ), target :: outputTimes_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *outputTimes_"/>
    !!]

    return
  end function comovingDistancesConstructorInternal

  subroutine comovingDistancesDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`taskComovingDistances` task class.
    !!}
    implicit none
    type(taskComovingDistances), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"         />
    <objectDestructor name="self%outputTimes_"                />
    !!]
    return
  end subroutine comovingDistancesDestructor

  subroutine comovingDistancesPerform(self,status)
    !!{RST
    Compute and output the comoving distances to each output.
    !!}
    use            :: Display         , only : displayIndent     , displayUnindent
    use            :: Error           , only : errorStatusSuccess
    use            :: Output_HDF5     , only : outputFile
    use            :: IO_HDF5         , only : hdf5Object, hdf5Group
    use            :: HDF5_Access     , only : hdf5Access
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    use            :: String_Handling , only : operator(//)
    implicit none
    class  (taskComovingDistances), intent(inout), target   :: self
    integer                       , intent(  out), optional :: status
    integer(c_size_t             )                          :: output
    type   (hdf5Group            )                          :: outputsGroup, outputGroup
    type   (varying_string       )                          :: groupName   , description

    call displayIndent('Begin task: comoving distances')
    ! Open the group for output time information.
    !$ call hdf5Access%set()
    outputsGroup  =outputFile%openGroup('Outputs','Group containing datasets relating to output times.')
    !$ call hdf5Access%unset()
    ! Iterate over output times and output data.
    do output=1,self%outputTimes_%count()
       groupName  ='Output'
       description='Data for output number '
       groupName  =groupName  //output
       description=description//output
       outputGroup=outputsGroup%openGroup(char(groupName),char(description))
       !$ call hdf5Access%set()
       call outputGroup%writeAttribute(self%cosmologyFunctions_%distanceComoving(self%outputTimes_%time(output)),'distanceComoving')
       !$ call hdf5Access%unset()
    end do
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: comoving distances')
    return
 end subroutine comovingDistancesPerform
