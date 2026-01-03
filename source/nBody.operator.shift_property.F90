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
  Implements an N-body data operator which shifts values of a property by an integer amount.
  !!}
  
  !![
  <nbodyOperator name="nbodyOperatorShiftProperty">
   <description>An N-body data operator which shifts values of a property by an integer amount.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorShiftProperty
     !!{
     An N-body data operator which shifts values of a property by an integer amount.
     !!}
     private
     type   (varying_string) :: propertyName
     integer(c_size_t      ) :: shiftBy
   contains
     procedure :: operate => shiftPropertyOperate
  end type nbodyOperatorShiftProperty

  interface nbodyOperatorShiftProperty
     !!{
     Constructors for the \refClass{nbodyOperatorShiftProperty} N-body operator class.
     !!}
     module procedure shiftPropertyConstructorParameters
     module procedure shiftPropertyConstructorInternal
  end interface nbodyOperatorShiftProperty

contains

  function shiftPropertyConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorShiftProperty} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nbodyOperatorShiftProperty)                :: self
    type   (inputParameters           ), intent(inout) :: parameters
    integer(c_size_t                  )                :: shiftBy
    type   (varying_string            )                :: propertyName
    
    !![
    <inputParameter>
      <name>propertyName</name>
      <source>parameters</source>
      <description>A named property on which to select.</description>
    </inputParameter>
    <inputParameter>
      <name>shiftBy</name>
      <source>parameters</source>
      <description>The amount by which to shift the property.</description>
    </inputParameter>
    !!]
    self=nbodyOperatorShiftProperty(propertyName,shiftBy)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function shiftPropertyConstructorParameters

  function shiftPropertyConstructorInternal(propertyName,shiftBy) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorShiftProperty} N-body operator class.
    !!}
    implicit none
    type   (nbodyOperatorShiftProperty)                   :: self
    integer(c_size_t                  ), intent(in   ) :: shiftBy
    type   (varying_string            ), intent(in   ) :: propertyName
    !![
    <constructorAssign variables="propertyName, shiftBy"/>
    !!]

    return
  end function shiftPropertyConstructorInternal

  subroutine shiftPropertyOperate(self,simulations)
    !!{
    Select particles matching a list of integer properties. 
    !!}
    use :: Error  , only : Error_Report
    use :: Display, only : displayIndent, displayMessage, displayUnindent, verbosityLevelStandard
    implicit none
    class  (nbodyOperatorShiftProperty), intent(inout)               :: self
    type   (nBodyData                 ), intent(inout), dimension(:) :: simulations
    integer(c_size_t                  ), pointer      , dimension(:) :: propertyInteger
    integer                                                          :: i
    
    call displayIndent('shift a property values',verbosityLevelStandard)
    do i=1,size(simulations)
       if (simulations(i)%propertiesInteger%exists(self%propertyName)) then
          propertyInteger =>  simulations(i)%propertiesInteger%value(self%propertyName)
          propertyInteger =  +propertyInteger &
               &             +self%shiftBy
          nullify(propertyInteger)
       else
          call Error_Report('property "'//self%propertyName//'"does not exist'//{introspection:location})
       end if
    end do
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine shiftPropertyOperate
