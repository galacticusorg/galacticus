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

!% Contains a module which implements an N-body data operator which filters out particles based on a property range.
  
  type :: propertyRange
     !% Type used to store filter ranges.
     type            (varying_string) :: name
     double precision                 :: rangeLowReal   , rangeHighReal
     integer         (c_size_t      ) :: rangeLowInteger, rangeHighInteger
     integer                          :: type
  end type propertyRange
  
  !# <nbodyOperator name="nbodyOperatorFilterProperties">
  !#  <description>An N-body data operator which filters out particles based on a property range.</description>
  !# </nbodyOperator>
  type, extends(nbodyOperatorClass) :: nbodyOperatorFilterProperties
     !% An N-body data operator which filters out particles based on a property range.
     private
     type(varying_string), allocatable, dimension(:) :: rangeLow      , rangeHigh, &
          &                                             propertyNames
     type(propertyRange ), allocatable, dimension(:) :: propertyRanges
   contains
     procedure :: operate => filterPropertiesOperate
  end type nbodyOperatorFilterProperties

  interface nbodyOperatorFilterProperties
     !% Constructors for the {\normalfont \ttfamily filterProperties} N-body operator class.
     module procedure filterPropertiesConstructorParameters
     module procedure filterPropertiesConstructorInternal
  end interface nbodyOperatorFilterProperties

contains

  function filterPropertiesConstructorParameters(parameters) result (self)
    !% Constructor for the {\normalfont \ttfamily filterProperties} N-body operator class which takes a parameter set as input.
    use :: Input_Parameters     , only : inputParameters
    use :: NBody_Simulation_Data, only : propertyTypeInteger, propertyTypeReal, propertyTypeUnknown, nBodyDataPropertyType
    implicit none
    type     (nbodyOperatorFilterProperties)                              :: self
    type     (inputParameters              ), intent(inout)               :: parameters
    type     (varying_string               ), allocatable  , dimension(:) :: rangeLow      , rangeHigh, &
         &                                                                   propertyNames
    type     (propertyRange                ), allocatable  , dimension(:) :: propertyRanges
    character(len=64                       )                              :: range
    integer                                                               :: i
    
    allocate(propertyNames(parameters%count('propertyNames',zeroIfNotPresent=.true.)))
    allocate(rangeLow     (parameters%count('rangeLow'     ,zeroIfNotPresent=.true.)))
    allocate(rangeHigh    (parameters%count('rangeHigh'    ,zeroIfNotPresent=.true.)))
    if (size(rangeLow ) /= size(propertyNames)) call Galacticus_Error_Report('[rangeLow] must have same cardinality as [propertyNames]' //{introspection:location})
    if (size(rangeHigh) /= size(propertyNames)) call Galacticus_Error_Report('[rangeHigh] must have same cardinality as [propertyNames]'//{introspection:location})
    !# <inputParameter>
    !#   <name>propertyNames</name>
    !#   <source>parameters</source>
    !#   <description>A list of named properties on which to filter.</description>
    !#   <type>string</type>
    !#   <cardinality>0..*</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>rangeLow</name>
    !#   <source>parameters</source>
    !#   <description>The lowest value of each property to pass (``{\normalfont \ttfamily -infinity}'' is interpreted as the lowest possible value for the property.</description>
    !#   <type>string</type>
    !#   <cardinality>0..*</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>rangeHigh</name>
    !#   <source>parameters</source>
    !#   <description>The highest value of each property to pass (``{\normalfont \ttfamily +infinity}'' is interpreted as the lowest possible value for the property.</description>
    !#   <type>string</type>
    !#   <cardinality>0..*</cardinality>
    !# </inputParameter>
    allocate(propertyRanges(size(propertyNames)))
    do i=1,size(propertyNames)
       propertyRanges(i)%name=                           propertyNames(i)
       propertyRanges(i)%type=nBodyDataPropertyType(char(propertyNames(i)))
       select case (propertyRanges(i)%type)
       case (propertyTypeInteger)
          if (rangeLow( i) == "-infinity") then
             propertyRanges(i)%rangeLowInteger =-huge(0_c_size_t)
          else
             range=rangeLow (i)
             read (range,*) propertyRanges(i)%rangeLowInteger
          end if
          if (rangeHigh(i) == "+infinity") then
             propertyRanges(i)%rangeHighInteger=+huge(0_c_size_t)
          else
             range=rangeHigh(i)
             read (range,*) propertyRanges(i)%rangeHighInteger
          end if
       case (propertyTypeReal   )
          if (rangeLow( i) == "-infinity") then
             propertyRanges(i)%rangeLowReal    =-huge(0.0d0   )
          else
             range=rangeLow (i)
             read (range,*) propertyRanges(i)%rangeLowReal
          end if
          if (rangeHigh(i) == "+infinity") then
             propertyRanges(i)%rangeHighReal   =+huge(0.0d0   )
          else
             range=rangeHigh(i)
             read (range,*) propertyRanges(i)%rangeHighReal
          end if
       case (propertyTypeUnknown)
          call Galacticus_Error_Report('unknown property "'//char(propertyNames(i))//'"'//{introspection:location})
       end select
    end do
    self=nbodyOperatorFilterProperties(propertyRanges)
    !# <inputParametersValidate source="parameters"/>
    return
  end function filterPropertiesConstructorParameters

  function filterPropertiesConstructorInternal(propertyRanges) result (self)
    !% Internal constructor for the {\normalfont \ttfamily filterProperties} N-body operator class.
    use :: Galacticus_Error     , only : Galacticus_Error_Report
    use :: NBody_Simulation_Data, only : propertyTypeInteger    , propertyTypeReal
    implicit none
    type   (nbodyOperatorFilterProperties)                              :: self
    type   (propertyRange                ), intent(in   ), dimension(:) :: propertyRanges
    integer                                                             :: i
    !# <constructorAssign variables="propertyRanges"/>

    ! Validate ranges.
    do i=1,size(propertyRanges)
       select case (propertyRanges(i)%type)
       case (propertyTypeInteger)
          if (propertyRanges(i)%rangeLowInteger > propertyRanges(i)%rangeHighInteger) &
               & call Galacticus_Error_Report('range for property "'//char(propertyRanges(i)%name)//'" will exclude all'//{introspection:location})
       case (propertyTypeReal   )
          if (propertyRanges(i)%rangeLowReal    > propertyRanges(i)%rangeHighReal   ) &
               & call Galacticus_Error_Report('range for property "'//char(propertyRanges(i)%name)//'" will exclude all'//{introspection:location})
       end select
    end do
    return
  end function filterPropertiesConstructorInternal

  subroutine filterPropertiesOperate(self,simulations)
    !% Identify and flag particles which have been always isolated.
    use :: Galacticus_Display   , only : Galacticus_Display_Indent, Galacticus_Display_Unindent, Galacticus_Display_Message, verbosityStandard
    use :: Galacticus_Error     , only : Galacticus_Error_Report
    use :: NBody_Simulation_Data, only : propertyTypeInteger      , propertyTypeReal
    implicit none
    class           (nbodyOperatorFilterProperties), intent(inout)                 :: self
    type            (nBodyData                    ), intent(inout), dimension(  :) :: simulations
    logical                                        , allocatable  , dimension(  :) :: mask
    integer         (c_size_t                     ), allocatable  , dimension(  :) :: propertyInteger, propertyIntegerFiltered
    double precision                               , allocatable  , dimension(  :) :: propertyReal   , propertyRealFiltered
    double precision                               , allocatable  , dimension(:,:) :: position
    integer                                                                        :: i              , j
    integer         (c_size_t                     )                                :: countFiltered
    
    call Galacticus_Display_Indent('filter on property ranges',verbosityStandard)
    do i=1,size(simulations)
       allocate(mask(size(simulations(i)%particleIDs)))
       mask=.true.
       do j=1,size(self%propertyRanges)
          select case (self%propertyRanges(j)%type)
          case (propertyTypeInteger)
             if (simulations(i)%propertiesInteger%exists(self%propertyRanges(j)%name)) then
                propertyInteger=simulations(i)%propertiesInteger%value(self%propertyRanges(j)%name)
                mask           = mask                                                       &
                     &          .and.                                                       &
                     &           propertyInteger >= self%propertyRanges(j)%rangeLowInteger  &
                     &          .and.                                                       &
                     &           propertyInteger <= self%propertyRanges(j)%rangeHighInteger
                deallocate(propertyInteger)
             else
                call Galacticus_Error_Report('property "'//self%propertyRanges(j)%name//'"does not exist'//{introspection:location})
             end if
          case (propertyTypeReal   )
             if (simulations(i)%propertiesReal   %exists(self%propertyRanges(j)%name)) then
                propertyReal   =simulations(i)%propertiesReal   %value(self%propertyRanges(j)%name)
                mask           = mask                                                       &
                     &          .and.                                                       &
                     &           propertyReal    >= self%propertyRanges(j)%rangeLowReal     &
                     &          .and.                                                       &
                     &           propertyReal    <= self%propertyRanges(j)%rangeHighReal   
                deallocate(propertyReal   )
             else
                call Galacticus_Error_Report('property "'//self%propertyRanges(j)%name//'"does not exist'//{introspection:location})
             end if
          case default
             call Galacticus_Error_Report('unsupported property type'//{introspection:location})
          end select
       end do
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
  end subroutine filterPropertiesOperate
