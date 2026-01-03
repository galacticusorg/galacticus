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

  !!{
  Implements a stellar population spectra postprocessor builder which simply looks up postprocessors by name.
  !!}

  !![
  <stellarPopulationSpectraPostprocessorBuilder name="stellarPopulationSpectraPostprocessorBuilderLookup">
   <description>A stellar population spectra postprocessor builder which simply looks up postprocessors by name.</description>
   <descriptorSpecial>descriptorSpecial</descriptorSpecial>
  </stellarPopulationSpectraPostprocessorBuilder>
  !!]
  type, extends(stellarPopulationSpectraPostprocessorBuilderClass) :: stellarPopulationSpectraPostprocessorBuilderLookup
     !!{
     A stellar population spectra postprocessor builder which simply looks up postprocessors by name.
     !!}
     private
     type(varying_string                           ), allocatable, dimension(:) :: names
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
     !!{
     Constructors for the \refClass{stellarPopulationSpectraPostprocessorBuilderLookup} stellar population spectra postprocessor builder class.
     !!}
     module procedure lookupConstructorParameters
     module procedure lookupConstructorInternal
  end interface stellarPopulationSpectraPostprocessorBuilderLookup

contains

  function lookupConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{stellarPopulationSpectraPostprocessorBuilderLookup} stellar population spectra postprocessor builder class which takes a
    parameter list as input.
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
    <workaround type="gfortran" PR="37336" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=37336">
      <description>Missing finalization of array constructors after their use.</description>
    </workaround>
    !!]   
    namesDefault=var_str('default')
    !![
    <inputParameter>
      <name>names</name>
      <defaultValue>namesDefault</defaultValue>
      <description>The names assigned to stellar spectra postprocessors.</description>
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
    !!{
    Internal constructor for the \refClass{stellarPopulationSpectraPostprocessorBuilderLookup} stellar population spectra postprocessor builder.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type   (stellarPopulationSpectraPostprocessorBuilderLookup)                              :: self
    type   (varying_string                                    ), intent(in   ), dimension(:) :: names
    type   (stellarPopulationSpectraPostprocessorList         ), intent(in   ), dimension(:) :: postprocessors
    integer                                                                                  :: i
    !![
    <constructorAssign variables="names, postprocessors"/>
    !!]

    if (size(names) /= size(postprocessors)) call Error_Report('number of names must match number of postprocessors'//{introspection:location})
    do i=1,size(postprocessors)
       !![
       <referenceCountIncrement owner="self%postprocessors(i)" object="stellarPopulationSpectraPostprocessor_"/>
       !!]
    end do
    return
  end function lookupConstructorInternal

  subroutine lookupDestructor(self)
    !!{
    Destructor for the \refClass{stellarPopulationSpectraPostprocessorBuilderLookup} stellar population spectra postprocessor builder.
    !!}
    implicit none
    type   (stellarPopulationSpectraPostprocessorBuilderLookup), intent(inout) :: self
    integer                                                                    :: i

    if (.not.allocated(self%postprocessors)) return
    do i=1,size(self%postprocessors)
       !![
       <objectDestructor name="self%postprocessors(i)%stellarPopulationSpectraPostprocessor_"/>
       !!]
    end do
    return
  end subroutine lookupDestructor

  !![
  <workaround type="gfortran" PR="93422" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=93422">
   <description>
    If the function name is used as the result variable, instead of using "result(postprocessor)", this PR is triggered.
   </description>
  </workaround>
  !!]
  function lookupBuild(self,descriptor) result(postprocessor)
    !!{
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
    !!{
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
