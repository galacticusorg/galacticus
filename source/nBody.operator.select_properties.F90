!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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
  
  !% Contains a module which implements an N-body data operator which selects particles matching a list of integer properties.
  
  !# <nbodyOperator name="nbodyOperatorSelectProperties">
  !#  <description>An N-body data operator which selects particles matching a list of integer properties.</description>
  !# </nbodyOperator>
  type, extends(nbodyOperatorClass) :: nbodyOperatorSelectProperties
     !% An N-body data operator which selects particles matching a list of integer properties.
     private
     type   (varying_string)                            :: propertyName
     integer(c_size_t      ), allocatable, dimension(:) :: selectedValues
   contains
     procedure :: operate => selectPropertiesOperate
  end type nbodyOperatorSelectProperties

  interface nbodyOperatorSelectProperties
     !% Constructors for the {\normalfont \ttfamily selectProperties} N-body operator class.
     module procedure selectPropertiesConstructorParameters
     module procedure selectPropertiesConstructorInternal
  end interface nbodyOperatorSelectProperties

contains

  function selectPropertiesConstructorParameters(parameters) result (self)
    !% Constructor for the {\normalfont \ttfamily selectProperties} N-body operator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nbodyOperatorSelectProperties)                              :: self
    type   (inputParameters              ), intent(inout)               :: parameters
    integer(c_size_t                     ), allocatable  , dimension(:) :: selectedValues
    type   (varying_string               )                              :: propertyName
    
    !# <inputParameter>
    !#   <name>propertyName</name>
    !#   <source>parameters</source>
    !#   <description>A named property on which to select.</description>
    !# </inputParameter>
    allocate(selectedValues(parameters%count('selectedValues')))
    !# <inputParameter>
    !#   <name>selectedValues</name>
    !#   <source>parameters</source>
    !#   <description>A list of allowed values for the property</description>
    !# </inputParameter>
    self=nbodyOperatorSelectProperties(propertyName,selectedValues)
    !# <inputParametersValidate source="parameters"/>
    return
  end function selectPropertiesConstructorParameters

  function selectPropertiesConstructorInternal(propertyName,selectedValues) result (self)
    !% Internal constructor for the {\normalfont \ttfamily selectProperties} N-body operator class.
    implicit none
    type   (nbodyOperatorSelectProperties)                              :: self
    integer(c_size_t                     ), intent(in   ), dimension(:) :: selectedValues
    type   (varying_string               ), intent(in   )               :: propertyName
    !# <constructorAssign variables="propertyName, selectedValues"/>

    return
  end function selectPropertiesConstructorInternal

  subroutine selectPropertiesOperate(self,simulations)
    !% Select particles matching a list of integer properties. 
    use :: Galacticus_Display, only : Galacticus_Display_Indent, Galacticus_Display_Unindent, Galacticus_Display_Message, verbosityStandard
    implicit none
    class           (nbodyOperatorSelectProperties), intent(inout)                 :: self
    type            (nBodyData                    ), intent(inout), dimension(  :) :: simulations
    logical                                        , allocatable  , dimension(  :) :: mask
    integer         (c_size_t                     ), allocatable  , dimension(  :) :: propertyInteger, propertyIntegerFiltered
    double precision                               , allocatable  , dimension(  :) :: propertyReal   , propertyRealFiltered
    double precision                               , allocatable  , dimension(:,:) :: position
    integer                                                                        :: i              , j
    integer         (c_size_t                     )                                :: countFiltered
    
    call Galacticus_Display_Indent('select on property values',verbosityStandard)
    do i=1,size(simulations)
       allocate(mask(size(simulations(i)%particleIDs)))
       if (simulations(i)%propertiesInteger%exists(self%propertyName)) then
          propertyInteger=simulations(i)%propertiesInteger%value(self%propertyName)
          do j=1,size(mask)
             mask(j)=any(propertyInteger(j) == self%selectedValues)
          end do
          deallocate(propertyInteger)
       else
          call Galacticus_Error_Report('property "'//self%propertyName//'"does not exist'//{introspection:location})
       end if
       ! Filter default properties.
       countFiltered=count(mask)
       !! Position
       allocate(position(3,countFiltered))
       do j=1,3
          position(j,:)=pack(simulations(i)%position(j,:),mask)
       end do
       deallocate(simulations(i)%position)
       call move_alloc(position,simulations(i)%position)
       !! Velocity
       allocate(position(3,countFiltered))
       do j=1,3
          position(j,:)=pack(simulations(i)%velocity(j,:),mask)
       end do
       deallocate(simulations(i)%velocity)
       call move_alloc(position,simulations(i)%velocity)
       !! IDs
       allocate(propertyInteger(countFiltered))
       propertyInteger=pack(simulations(i)%particleIDs,mask)
       deallocate(simulations(i)%particleIDs)
       call move_alloc(propertyInteger,simulations(i)%particleIDs)
       ! Filter other properties.
       !! Integer properties.
       allocate(propertyIntegerFiltered(countFiltered))
       do j=1,simulations(i)%propertiesInteger%size()
          propertyInteger        =simulations(i)%propertiesInteger%value(j)
          propertyIntegerFiltered=pack(propertyInteger,mask)
          call simulations(i)%propertiesInteger%set(simulations(i)%propertiesInteger%key(j),propertyIntegerFiltered)
       end do
       deallocate(propertyIntegerFiltered)
       !! Real properties.
       allocate(propertyRealFiltered   (countFiltered))
       do j=1,simulations(i)%propertiesReal   %size()
          propertyReal           =simulations(i)%propertiesReal   %value(j)
          propertyRealFiltered   =pack(propertyReal   ,mask)
          call simulations(i)%propertiesReal   %set(simulations(i)%propertiesReal   %key(j),propertyRealFiltered   )
       end do
       deallocate(propertyRealFiltered   )
       deallocate(mask                   )
    end do
    call Galacticus_Display_Unindent('done',verbosityStandard)
    return
  end subroutine selectPropertiesOperate
