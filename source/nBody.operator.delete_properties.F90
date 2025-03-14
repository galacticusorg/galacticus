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

!!{
Implements an N-body data operator which deletes named properties from the simulation.
!!}

  !![
  <nbodyOperator name="nbodyOperatorDeleteProperties">
   <description>An N-body data operator which deletes named properties from the simulation.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorDeleteProperties
     !!{
     An N-body data operator which deletes named properties from the simulation.
     !!}
     private
    type(varying_string), allocatable  , dimension(:) :: propertyNames
   contains
     procedure :: operate => deletePropertiesOperate
  end type nbodyOperatorDeleteProperties

  interface nbodyOperatorDeleteProperties
     !!{
     Constructors for the {\normalfont \ttfamily deleteProperties} N-body operator class.
     !!}
     module procedure deletePropertiesConstructorParameters
     module procedure deletePropertiesConstructorInternal
  end interface nbodyOperatorDeleteProperties

contains

  function deletePropertiesConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily deleteProperties} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nbodyOperatorDeleteProperties)                              :: self
    type(inputParameters              ), intent(inout)               :: parameters
    type(varying_string               ), allocatable  , dimension(:) :: propertyNames

    allocate(propertyNames(parameters%count('propertyNames',zeroIfNotPresent=.true.)))
    !![
    <inputParameter>
      <name>propertyNames</name>
      <source>parameters</source>
      <description>A list of named properties to be deleted from the simulation.</description>
    </inputParameter>
    !!]
    self=nbodyOperatorDeleteProperties(propertyNames)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function deletePropertiesConstructorParameters

  function deletePropertiesConstructorInternal(propertyNames) result (self)
    !!{
    Internal constructor for the {\normalfont \ttfamily deleteProperties} N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorDeleteProperties)                              :: self
    type(varying_string               ), intent(in   ), dimension(:) :: propertyNames
    !![
    <constructorAssign variables="propertyNames"/>
    !!]
    
    return
  end function deletePropertiesConstructorInternal

  subroutine deletePropertiesOperate(self,simulations)
    !!{
    Identify and flag particles which have been always isolated.
    !!}
    use :: Display, only : displayIndent, displayMessage, displayUnindent, verbosityLevelStandard
    use :: Error  , only : Error_Report
    implicit none
    class           (nbodyOperatorDeleteProperties), intent(inout)                 :: self
    type            (nBodyData                    ), intent(inout), dimension(  :) :: simulations
    integer         (c_size_t                     ), pointer      , dimension(  :) :: propertyInteger
    double precision                               , pointer      , dimension(  :) :: propertyReal
    integer         (c_size_t                     ), pointer      , dimension(:,:) :: propertyIntegerRank1
    double precision                               , pointer      , dimension(:,:) :: propertyRealRank1
    logical                                                                        :: propertyFound
    integer                                                                        :: i            , j

    call displayIndent('delete named properties',verbosityLevelStandard)
    do i=1,size(self%propertyNames)
       do j=1,size(simulations)
          propertyFound=.false.
          if (simulations(j)%propertiesInteger     %exists(self%propertyNames(i))) then
             propertyFound        =  .true.
             propertyInteger      => simulations(j)%propertiesInteger     %value(self%propertyNames(i))
             deallocate(propertyInteger)
             call simulations(j)%propertiesInteger     %delete(self%propertyNames(i))
          end if
          if (simulations(j)%propertiesReal        %exists(self%propertyNames(i))) then
             propertyFound        =  .true.
             propertyReal         => simulations(j)%propertiesReal        %value(self%propertyNames(i))
             deallocate(propertyReal   )
             call simulations(j)%propertiesReal        %delete(self%propertyNames(i))
          end if
          if (simulations(j)%propertiesIntegerRank1%exists(self%propertyNames(i))) then
             propertyFound        =  .true.
             propertyIntegerRank1 => simulations(j)%propertiesIntegerRank1%value(self%propertyNames(i))
             deallocate(propertyIntegerRank1)
             call simulations(j)%propertiesIntegerRank1%delete(self%propertyNames(i))
          end if
          if (simulations(j)%propertiesRealRank1   %exists(self%propertyNames(i))) then
             propertyFound        =  .true.
             propertyRealRank1    => simulations(j)%propertiesRealRank1   %value(self%propertyNames(i))
             deallocate(propertyRealRank1   )
             call simulations(j)%propertiesRealRank1   %delete(self%propertyNames(i))
          end if
          if (propertyFound) then
             call displayMessage('deleted "' //self%propertyNames(i)//'"'                                    )
          else
             call Error_Report  ('property "'//self%propertyNames(i)//'" not found'//{introspection:location})
          end if
       end do
    end do
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine deletePropertiesOperate
