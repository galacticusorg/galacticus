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

  !!{RST
  Implements a stellar population spectra postprocessor builder which simply looks up postprocessors by name.
  !!}

  !![
  <stellarPopulationSpectraPostprocessorBuilder name="stellarPopulationSpectraPostprocessorBuilderLookup" docformat="rst">
   <description>
   A stellar population spectra postprocessor builder which simply looks up postprocessors by name.
   </description>
   <descriptorSpecial>descriptorSpecial</descriptorSpecial>
   <assignment forceArrayAssign="names"/>
  </stellarPopulationSpectraPostprocessorBuilder>
  !!]
  type, extends(stellarPopulationSpectraPostprocessorBuilderClass) :: stellarPopulationSpectraPostprocessorBuilderLookup
     !!{RST
     A stellar population spectra postprocessor builder which simply looks up postprocessors by name.
     !!}
     private
     !![
     <workaround type="gfortran" docformat="rst">
       <description>
       Having `names` be `allocatable` in the following (which would be preferable) results in a memory leak on assignment. I have not been able to figure out why this happens. So, we make it a `pointer` array, and handle deallocation manually (and set `forceArrayAssign` above so that assignment is done, rather than just pointer assignment).
       </description>
     </workaround>
     !!]
     type(varying_string                           ), pointer    , dimension(:) :: names          => null()
     type(stellarPopulationSpectraPostprocessorList), allocatable, dimension(:) :: postprocessors
   contains
     !![
     <methods>
       <method method="descriptorSpecial" description="Handle adding special parameters to the descriptor."/>
     </methods>
     !!]
     final     ::                      lookupDestructor
     procedure :: build             => lookupBuild
     procedure :: descriptorSpecial => lookupDescriptorSpecial
  end type stellarPopulationSpectraPostprocessorBuilderLookup

  interface stellarPopulationSpectraPostprocessorBuilderLookup
     !!{RST
     Constructors for the ``stellarPopulationSpectraPostprocessorBuilderLookup`` stellar population spectra postprocessor builder class.
     !!}
     module procedure lookupConstructorParameters
     module procedure lookupConstructorInternal
  end interface stellarPopulationSpectraPostprocessorBuilderLookup

contains

  function lookupConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``stellarPopulationSpectraPostprocessorBuilderLookup`` stellar population spectra postprocessor builder class which takes a parameter list as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (stellarPopulationSpectraPostprocessorBuilderLookup)                              :: self
    type   (inputParameters                                   ), intent(inout)               :: parameters
    type   (varying_string                                    ), allocatable  , dimension(:) :: names
    type   (stellarPopulationSpectraPostprocessorList         ), allocatable  , dimension(:) :: postprocessors
    type   (varying_string                                    )               , dimension(1) :: namesDefault
    integer                                                                                  :: countPostprocessors, countPostprocessorNames, &
         &                                                                                      i

    countPostprocessors    =max(parameters%copiesCount('stellarPopulationSpectraPostprocessor',zeroIfNotPresent=.true.),1)
    countPostprocessorNames=max(parameters%      count('names'                                ,zeroIfNotPresent=.true.),1)
    if (countPostprocessors /= countPostprocessorNames) call Error_Report('number of names must match number of postprocessors'//{introspection:location})
    allocate(names         (countPostprocessors))
    allocate(postprocessors(countPostprocessors))
    !![
    <workaround type="gfortran" PR="37336" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=37336" docformat="rst">
      <description>
      Missing finalization of array constructors after their use in gfortran (PR#37336); workaround initializes the default names array using ``var_str`` before passing to ``inputParameter``.
      </description>
    </workaround>
    !!]   
    namesDefault(1)=var_str('default')
    !![
    <inputParameter docformat="rst">
      <name>names</name>
      <defaultValue>namesDefault</defaultValue>
      <description>
      The string names (e.g., ``default``, ``intrinsic``) assigned to each stellar spectra postprocessor, used as lookup keys when the builder is asked to return the postprocessor appropriate for a given filter or descriptor.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="stellarPopulationSpectraPostprocessor" name="postprocessors(i)%stellarPopulationSpectraPostprocessor_" source="parameters" copy="i=1,countPostprocessors"/>
    !!]
    self=stellarPopulationSpectraPostprocessorBuilderLookup(names,postprocessors)
    !![
    <inputParametersValidate source="parameters" multiParameters="stellarPopulationSpectraPostprocessor"/>
    !!]
    do i=1,countPostprocessors
       !![
       <objectDestructor name="postprocessors(i)%stellarPopulationSpectraPostprocessor_"/>
       !!]
    end do
    return
  end function lookupConstructorParameters

  function lookupConstructorInternal(names,postprocessors) result(self)
    !!{RST
    Internal constructor for the ``stellarPopulationSpectraPostprocessorBuilderLookup`` stellar population spectra postprocessor builder.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type   (stellarPopulationSpectraPostprocessorBuilderLookup)                              :: self
    type   (varying_string                                    ), intent(in   ), dimension(:) :: names
    type   (stellarPopulationSpectraPostprocessorList         ), intent(in   ), dimension(:) :: postprocessors
    integer                                                                                  :: i
    !![
    <constructorAssign variables="postprocessors"/>
    !!]

    if (size(names) /= size(postprocessors)) call Error_Report('number of names must match number of postprocessors'//{introspection:location})
    allocate(self%names(size(names)))
    do i=1,size(postprocessors)
       self%names(i)=names(i)
       !![
       <referenceCountIncrement owner="self%postprocessors(i)" object="stellarPopulationSpectraPostprocessor_"/>
       !!]
    end do
    return
  end function lookupConstructorInternal

  subroutine lookupDestructor(self)
    !!{RST
    Destructor for the ``stellarPopulationSpectraPostprocessorBuilderLookup`` stellar population spectra postprocessor builder.
    !!}
    implicit none
    type   (stellarPopulationSpectraPostprocessorBuilderLookup), intent(inout) :: self
    integer                                                                    :: i

    if (associated(self%names)) deallocate(self%names)
    if (allocated(self%postprocessors)) then
       do i=1,size(self%postprocessors)
          !![
	  <objectDestructor name="self%postprocessors(i)%stellarPopulationSpectraPostprocessor_"/>
          !!]
       end do
    end if
    return
  end subroutine lookupDestructor

  function lookupBuild(self,descriptor) result(postprocessor)
    !!{RST
    Return a stellar population spectra postprocessor by lookup via name.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (stellarPopulationSpectraPostprocessorClass        ), pointer       :: postprocessor
    class  (stellarPopulationSpectraPostprocessorBuilderLookup), intent(inout) :: self
    type   (varying_string                                    ), intent(in   ) :: descriptor
    integer                                                                    :: i

    postprocessor => null()
    do i=1,size(self%names)
       if (self%names(i) == descriptor) postprocessor => self%postprocessors(i)%stellarPopulationSpectraPostprocessor_
    end do
    if (.not.associated(postprocessor)) call Error_Report('unable to located postprocessor "'//descriptor//'"'//{introspection:location})
    return
  end function lookupBuild

  subroutine lookupDescriptorSpecial(self,descriptor)
    !!{RST
    Add special parameters to the descriptor.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class  (stellarPopulationSpectraPostprocessorBuilderLookup), intent(inout) :: self
    type   (inputParameters                                   ), intent(inout) :: descriptor
    integer                                                                    :: i
    
    if (allocated(self%postprocessors)) then
       do i=1,size(self%postprocessors)
          call self%postprocessors(i)%stellarPopulationSpectraPostprocessor_%descriptor(descriptor)
       end do
    end if
    return
  end subroutine lookupDescriptorSpecial
