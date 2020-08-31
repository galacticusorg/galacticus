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

!% Contains a module which implements an N-body data operator which deletes named properties from the simulation.

  !# <nbodyOperator name="nbodyOperatorDeleteProperties">
  !#  <description>An N-body data operator which deletes named properties from the simulation.</description>
  !# </nbodyOperator>
  type, extends(nbodyOperatorClass) :: nbodyOperatorDeleteProperties
     !% An N-body data operator which deletes named properties from the simulation.
     private
    type(varying_string), allocatable  , dimension(:) :: propertyNames
   contains
     procedure :: operate => deletePropertiesOperate
  end type nbodyOperatorDeleteProperties

  interface nbodyOperatorDeleteProperties
     !% Constructors for the {\normalfont \ttfamily deleteProperties} N-body operator class.
     module procedure deletePropertiesConstructorParameters
     module procedure deletePropertiesConstructorInternal
  end interface nbodyOperatorDeleteProperties

contains

  function deletePropertiesConstructorParameters(parameters) result (self)
    !% Constructor for the {\normalfont \ttfamily deleteProperties} N-body operator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nbodyOperatorDeleteProperties)                              :: self
    type(inputParameters              ), intent(inout)               :: parameters
    type(varying_string               ), allocatable  , dimension(:) :: propertyNames

    allocate(propertyNames(parameters%count('propertyNames',zeroIfNotPresent=.true.)))
    !# <inputParameter>
    !#   <name>propertyNames</name>
    !#   <source>parameters</source>
    !#   <description>A list of named properties to be deleted from the simulation.</description>
    !#   <type>string</type>
    !#   <cardinality>0..*</cardinality>
    !# </inputParameter>
    self=nbodyOperatorDeleteProperties(propertyNames)
    !# <inputParametersValidate source="parameters"/>
    return
  end function deletePropertiesConstructorParameters

  function deletePropertiesConstructorInternal(propertyNames) result (self)
    !% Internal constructor for the {\normalfont \ttfamily deleteProperties} N-body operator class.
    implicit none
    type(nbodyOperatorDeleteProperties)                              :: self
    type(varying_string               ), intent(in   ), dimension(:) :: propertyNames
    !# <constructorAssign variables="propertyNames"/>
    
    return
  end function deletePropertiesConstructorInternal

  subroutine deletePropertiesOperate(self,simulations)
    !% Identify and flag particles which have been always isolated.
    use :: Galacticus_Display, only : Galacticus_Display_Indent, Galacticus_Display_Unindent, Galacticus_Display_Message, verbosityStandard
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    implicit none
    class  (nbodyOperatorDeleteProperties), intent(inout)               :: self
    type   (nBodyData                    ), intent(inout), dimension(:) :: simulations
    logical                                                             :: propertyFound
    integer                                                             :: i            , j

    call Galacticus_Display_Indent('delete named properties',verbosityStandard)
    do i=1,size(self%propertyNames)
       do j=1,size(simulations)
          propertyFound=.false.
          if (simulations(j)%propertiesInteger%exists(self%propertyNames(i))) then
             propertyFound=.true.
             call simulations(j)%propertiesInteger%delete(self%propertyNames(i))
          end if
          if (simulations(j)%propertiesReal   %exists(self%propertyNames(i))) then
             propertyFound=.true.
             call simulations(j)%propertiesReal   %delete(self%propertyNames(i))
          end if
          if (propertyFound) then
             call Galacticus_Display_Message('deleted "' //self%propertyNames(i)//'"'                                    )
          else
             call Galacticus_Error_Report   ('property "'//self%propertyNames(i)//'" not found'//{introspection:location})
          end if
       end do
    end do
    call Galacticus_Display_Unindent('done',verbosityStandard)
    return
  end subroutine deletePropertiesOperate
