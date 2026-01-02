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

  use, intrinsic :: ISO_C_Binding, only : c_size_t

  !![
  <evolveForestsWorkShare name="evolveForestsWorkShareStride">
   <description>A forest evolution work sharing class in which forests are assigned by another work sharing class, but then strode over in steps of a specified size.</description>
  </evolveForestsWorkShare>
  !!]
  type, extends(evolveForestsWorkShareClass) :: evolveForestsWorkShareStride
     !!{
     Implementation of a forest evolution work sharing class in which forests are assigned by cycling through processes.
     !!}
     private
     integer(c_size_t                   )          :: stride                           , offset
     class  (evolveForestsWorkShareClass), pointer :: evolveForestsWorkShare_ => null()
   contains
     final     ::                 strideDestructor
     procedure :: forestNumber => strideForestNumber
  end type evolveForestsWorkShareStride

  interface evolveForestsWorkShareStride
     !!{
     Constructors for the \refClass{evolveForestsWorkShareStride} forest evolution work sharing class.
     !!}
     module procedure strideConstructorParameters
     module procedure strideConstructorInternal
  end interface evolveForestsWorkShareStride

contains

  function strideConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{evolveForestsWorkShareStride} forest evolution work sharing class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (evolveForestsWorkShareStride)                :: self
    type   (inputParameters             ), intent(inout) :: parameters
    class  (evolveForestsWorkShareClass ), pointer       :: evolveForestsWorkShare_
    integer(c_size_t                    )                :: stride                 , offset

    !![
    <inputParameter>
      <name>stride</name>
      <description>The size of the stride to take over forests.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>offset</name>
      <description>The offset of the stride to take over forests.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="evolveForestsWorkShare" name="evolveForestsWorkShare_" source="parameters"/>
    !!]
    self=evolveForestsWorkShareStride(stride,offset,evolveForestsWorkShare_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="evolveForestsWorkShare_"/>
    !!]
    return
  end function strideConstructorParameters

  function strideConstructorInternal(stride,offset,evolveForestsWorkShare_) result(self)
    !!{
    Internal constructor for the \refClass{evolveForestsWorkShareStride} forest evolution work sharing class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type   (evolveForestsWorkShareStride)                        :: self
    integer(c_size_t                    ), intent(in   )         :: stride                 , offset
    class  (evolveForestsWorkShareClass ), intent(in   ), target :: evolveForestsWorkShare_
    !![
    <constructorAssign variables="stride, offset, *evolveForestsWorkShare_"/>
    !!]

    if (offset < 0_c_size_t .or. offset >= stride) call Error_Report('0 ≤ [offset] < [stride]'//{introspection:location})
    return
  end function strideConstructorInternal

  subroutine strideDestructor(self)
    !!{
    Destructor for the \refClass{evolveForestsWorkShareStride} forest evolution work sharing class.
    !!}
    implicit none
    type(evolveForestsWorkShareStride), intent(inout) :: self

    !![
    <objectDestructor name="self%evolveForestsWorkShare_"/>
    !!]
    return
  end subroutine strideDestructor

  function strideForestNumber(self,utilizeOpenMPThreads)
    !!{
    Return the number of the next forest to process.
    !!}
    implicit none
    integer(c_size_t                    )                :: strideForestNumber
    class  (evolveForestsWorkShareStride), intent(inout) :: self
    logical                              , intent(in   ) :: utilizeOpenMPThreads

    strideForestNumber=+(                                                                 &
         &               +self%evolveForestsWorkShare_%forestNumber(utilizeOpenMPThreads) &
         &               -1_c_size_t                                                      &
         &              )                                                                 &
         &             *self%stride                                                       &
         &             +self%offset                                                       &
         &             +1_c_size_t
    return
  end function strideForestNumber
