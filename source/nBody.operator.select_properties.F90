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
  Implements an N-body data operator which selects particles matching a list of integer properties.
  !!}
  
  !![
  <nbodyOperator name="nbodyOperatorSelectProperties">
   <description>An N-body data operator which selects particles matching a list of integer properties.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorSelectProperties
     !!{
     An N-body data operator which selects particles matching a list of integer properties.
     !!}
     private
     type   (varying_string)                            :: propertyName
     integer(c_size_t      ), allocatable, dimension(:) :: selectedValues
   contains
     procedure :: operate => selectPropertiesOperate
  end type nbodyOperatorSelectProperties

  interface nbodyOperatorSelectProperties
     !!{
     Constructors for the \refClass{nbodyOperatorSelectProperties} N-body operator class.
     !!}
     module procedure selectPropertiesConstructorParameters
     module procedure selectPropertiesConstructorInternal
  end interface nbodyOperatorSelectProperties

contains

  function selectPropertiesConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorSelectProperties} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nbodyOperatorSelectProperties)                              :: self
    type   (inputParameters              ), intent(inout)               :: parameters
    integer(c_size_t                     ), allocatable  , dimension(:) :: selectedValues
    type   (varying_string               )                              :: propertyName
    
    !![
    <inputParameter>
      <name>propertyName</name>
      <source>parameters</source>
      <description>A named property on which to select.</description>
    </inputParameter>
    !!]
    allocate(selectedValues(parameters%count('selectedValues')))
    !![
    <inputParameter>
      <name>selectedValues</name>
      <source>parameters</source>
      <description>A list of allowed values for the property</description>
    </inputParameter>
    !!]
    self=nbodyOperatorSelectProperties(propertyName,selectedValues)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function selectPropertiesConstructorParameters

  function selectPropertiesConstructorInternal(propertyName,selectedValues) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorSelectProperties} N-body operator class.
    !!}
    implicit none
    type   (nbodyOperatorSelectProperties)                              :: self
    integer(c_size_t                     ), intent(in   ), dimension(:) :: selectedValues
    type   (varying_string               ), intent(in   )               :: propertyName
    !![
    <constructorAssign variables="propertyName, selectedValues"/>
    !!]

    return
  end function selectPropertiesConstructorInternal

  subroutine selectPropertiesOperate(self,simulations)
    !!{
    Select particles matching a list of integer properties. 
    !!}
    use :: Display, only : displayIndent, displayMessage, displayUnindent, verbosityLevelStandard
    use :: Error  , only : Error_Report
    implicit none
    class           (nbodyOperatorSelectProperties), intent(inout)                 :: self
    type            (nBodyData                    ), intent(inout), dimension(  :) :: simulations
    logical                                        , allocatable  , dimension(  :) :: mask
    integer         (c_size_t                     ), pointer      , dimension(  :) :: propertyInteger     , propertyIntegerFiltered
    double precision                               , pointer      , dimension(  :) :: propertyReal        , propertyRealFiltered
    integer         (c_size_t                     ), pointer      , dimension(:,:) :: propertyIntegerRank1, propertyIntegerRank1Filtered
    double precision                               , pointer      , dimension(:,:) :: propertyRealRank1   , propertyRealRank1Filtered
    integer                                                                        :: i                   , j                           , &
         &                                                                            k
    integer         (c_size_t                     )                                :: countFiltered
    
    call displayIndent('select on property values',verbosityLevelStandard)
    do i=1,size(simulations)
       if (simulations(i)%propertiesInteger%exists(self%propertyName)) then
          propertyInteger => simulations(i)%propertiesInteger%value(self%propertyName)
          allocate(mask(size(propertyInteger)))
          do j=1,size(mask)
             mask(j)=any(propertyInteger(j) == self%selectedValues)
          end do
          nullify(propertyInteger)
       else
          call Error_Report('property "'//self%propertyName//'"does not exist'//{introspection:location})
       end if
       countFiltered=count(mask)
       ! Filter properties.
       !! Integer properties.
       do j=1,simulations(i)%propertiesInteger    %size()
          propertyInteger      => simulations(i)%propertiesInteger    %value(j)
          allocate(propertyIntegerFiltered    (                                  countFiltered))
          propertyIntegerFiltered             =pack(propertyInteger          ,mask)
          call simulations(i)%propertiesInteger    %set(simulations(i)%propertiesInteger        %key(j),propertyIntegerFiltered    )
          deallocate(propertyInteger            )
          nullify   (propertyIntegerFiltered    )
       end do
       !! Real properties.
       do j=1,simulations(i)%propertiesReal   %size()
          propertyReal         => simulations(i)%propertiesReal       %value(j)
          allocate(propertyRealFiltered   (                                      countFiltered))
          propertyRealFiltered                =pack(propertyReal             ,mask)
          call simulations(i)%propertiesReal       %set(simulations(i)%propertiesReal           %key(j),propertyRealFiltered       )
          deallocate(propertyReal               )
          nullify   (propertyRealFiltered       )
       end do
       !! Integer rank-1 properties.
       do j=1,simulations(i)%propertiesIntegerRank1%size()
          propertyIntegerRank1 => simulations(i)%propertiesIntegerRank1%value(j)
          allocate(propertyIntegerRank1Filtered(size(propertyIntegerRank1,dim=1),countFiltered))
          do k=1,size(propertyIntegerRank1,dim=1)
             propertyIntegerRank1Filtered(k,:)=pack(propertyIntegerRank1(k,:),mask)
          end do
          call simulations(i)%propertiesIntegerRank1%set(simulations(i)%propertiesIntegerRank1%key(j),propertyIntegerRank1Filtered)
          deallocate(propertyIntegerRank1        )
          nullify   (propertyIntegerRank1Filtered)
       end do
       !! Real rank-1 properties.
       do j=1,simulations(i)%propertiesRealRank1   %size()
          propertyRealRank1    => simulations(i)%propertiesRealRank1   %value(j)
          allocate(propertyRealRank1Filtered   (size(propertyRealRank1   ,dim=1),countFiltered))
          do k=1,size(propertyRealRank1   ,dim=1)
             propertyRealRank1Filtered  (k,:)=  pack(propertyRealRank1   (k,:),mask)
          end do
          call simulations(i)%propertiesRealRank1   %set(simulations(i)%propertiesRealRank1   %key(j),propertyRealRank1Filtered   )
          deallocate(propertyRealRank1           )
          nullify   (propertyRealRank1Filtered   )
       end do
       deallocate(mask)
    end do
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine selectPropertiesOperate
