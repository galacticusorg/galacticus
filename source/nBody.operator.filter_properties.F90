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
Implements an N-body data operator which filters out particles based on a property range.
!!}

  use :: NBody_Simulation_Data, only : enumerationPropertyTypeType

  type :: propertyRange
     !!{
     Type used to store filter ranges.
     !!}
     type            (varying_string             ) :: name
     double precision                              :: rangeLowReal   , rangeHighReal
     integer         (c_size_t                   ) :: rangeLowInteger, rangeHighInteger
     type            (enumerationPropertyTypeType) :: type
  end type propertyRange
  
  !![
  <nbodyOperator name="nbodyOperatorFilterProperties">
   <description>An N-body data operator which filters out particles based on a property range.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorFilterProperties
     !!{
     An N-body data operator which filters out particles based on a property range.
     !!}
     private
     type(varying_string), allocatable, dimension(:) :: rangeLow      , rangeHigh, &
          &                                             propertyNames
     type(propertyRange ), allocatable, dimension(:) :: propertyRanges
   contains
     procedure :: operate => filterPropertiesOperate
  end type nbodyOperatorFilterProperties

  interface nbodyOperatorFilterProperties
     !!{
     Constructors for the \refClass{nbodyOperatorFilterProperties} N-body operator class.
     !!}
     module procedure filterPropertiesConstructorParameters
     module procedure filterPropertiesConstructorInternal
  end interface nbodyOperatorFilterProperties

contains

  function filterPropertiesConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorFilterProperties} N-body operator class which takes a parameter set as input.
    !!}
    use :: Error                , only : Error_Report
    use :: Input_Parameters     , only : inputParameters
    use :: NBody_Simulation_Data, only : nBodyDataPropertyType, propertyTypeInteger, propertyTypeReal, propertyTypeUnknown
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
    if (size(rangeLow ) /= size(propertyNames)) call Error_Report('[rangeLow] must have same cardinality as [propertyNames]' //{introspection:location})
    if (size(rangeHigh) /= size(propertyNames)) call Error_Report('[rangeHigh] must have same cardinality as [propertyNames]'//{introspection:location})
    !![
    <inputParameter>
      <name>propertyNames</name>
      <source>parameters</source>
      <description>A list of named properties on which to filter.</description>
    </inputParameter>
    <inputParameter>
      <name>rangeLow</name>
      <source>parameters</source>
      <description>The lowest value of each property to pass (``{\normalfont \ttfamily -infinity}'' is interpreted as the lowest possible value for the property.</description>
    </inputParameter>
    <inputParameter>
      <name>rangeHigh</name>
      <source>parameters</source>
      <description>The highest value of each property to pass (``{\normalfont \ttfamily +infinity}'' is interpreted as the lowest possible value for the property.</description>
    </inputParameter>
    !!]
    allocate(propertyRanges(size(propertyNames)))
    do i=1,size(propertyNames)
       propertyRanges(i)%name=                           propertyNames(i)
       propertyRanges(i)%type=nBodyDataPropertyType(char(propertyNames(i)))
       select case (propertyRanges(i)%type%ID)
       case (propertyTypeInteger%ID)
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
       case (propertyTypeReal   %ID)
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
       case (propertyTypeUnknown%ID)
          call Error_Report('unknown property "'//char(propertyNames(i))//'"'//{introspection:location})
       end select
    end do
    self=nbodyOperatorFilterProperties(propertyRanges)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function filterPropertiesConstructorParameters

  function filterPropertiesConstructorInternal(propertyRanges) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorFilterProperties} N-body operator class.
    !!}
    use :: Error                , only : Error_Report
    use :: NBody_Simulation_Data, only : propertyTypeInteger, propertyTypeReal
    implicit none
    type   (nbodyOperatorFilterProperties)                              :: self
    type   (propertyRange                ), intent(in   ), dimension(:) :: propertyRanges
    integer                                                             :: i
    !![
    <constructorAssign variables="propertyRanges"/>
    !!]

    ! Validate ranges.
    do i=1,size(propertyRanges)
       select case (propertyRanges(i)%type%ID)
       case (propertyTypeInteger%ID)
          if (propertyRanges(i)%rangeLowInteger > propertyRanges(i)%rangeHighInteger) &
               & call Error_Report('range for property "'//char(propertyRanges(i)%name)//'" will exclude all'//{introspection:location})
       case (propertyTypeReal   %ID)
          if (propertyRanges(i)%rangeLowReal    > propertyRanges(i)%rangeHighReal   ) &
               & call Error_Report('range for property "'//char(propertyRanges(i)%name)//'" will exclude all'//{introspection:location})
       end select
    end do
    return
  end function filterPropertiesConstructorInternal

  subroutine filterPropertiesOperate(self,simulations)
    !!{
    Identify and flag particles which have been always isolated.
    !!}
    use :: Display              , only : displayIndent      , displayMessage  , displayUnindent, verbosityLevelStandard
    use :: Error                , only : Error_Report
    use :: NBody_Simulation_Data, only : propertyTypeInteger, propertyTypeReal
    implicit none
    class           (nbodyOperatorFilterProperties), intent(inout)                 :: self
    type            (nBodyData                    ), intent(inout), dimension(  :) :: simulations
    logical                                        , allocatable  , dimension(  :) :: mask
    integer         (c_size_t                     ), pointer      , dimension(  :) :: propertyInteger     , propertyIntegerFiltered
    double precision                               , pointer      , dimension(  :) :: propertyReal        , propertyRealFiltered
    integer         (c_size_t                     ), pointer      , dimension(:,:) :: propertyIntegerRank1, propertyIntegerRank1Filtered
    double precision                               , pointer      , dimension(:,:) :: propertyRealRank1   , propertyRealRank1Filtered
    integer                                                                        :: i                   , j                           , &
         &                                                                            k
    integer         (c_size_t                     )                                :: countFiltered
    
    call displayIndent('filter on property ranges',verbosityLevelStandard)
    do i=1,size(simulations)
       do j=1,size(self%propertyRanges)
          select case (self%propertyRanges(j)%type%ID)
          case (propertyTypeInteger%ID)
             if (simulations(i)%propertiesInteger%exists(self%propertyRanges(j)%name)) then
                propertyInteger => simulations(i)%propertiesInteger%value(self%propertyRanges(j)%name)
                if (.not.allocated(mask)) then
                   allocate(mask(size(propertyInteger)))
                   mask=.true.
                end if
                mask           = mask                                                       &
                     &          .and.                                                       &
                     &           propertyInteger >= self%propertyRanges(j)%rangeLowInteger  &
                     &          .and.                                                       &
                     &           propertyInteger <= self%propertyRanges(j)%rangeHighInteger
                nullify(propertyInteger)
             else
                call Error_Report('property "'//self%propertyRanges(j)%name//'"does not exist'//{introspection:location})
             end if
          case (propertyTypeReal   %ID)
             if (simulations(i)%propertiesReal   %exists(self%propertyRanges(j)%name)) then
                propertyReal   => simulations(i)%propertiesReal   %value(self%propertyRanges(j)%name)
                if (.not.allocated(mask)) then
                   allocate(mask(size(propertyReal   )))
                   mask=.true.
                end if
                mask           = mask                                                       &
                     &          .and.                                                       &
                     &           propertyReal    >= self%propertyRanges(j)%rangeLowReal     &
                     &          .and.                                                       &
                     &           propertyReal    <= self%propertyRanges(j)%rangeHighReal   
               nullify(propertyReal   )
             else
                call Error_Report('property "'//self%propertyRanges(j)%name//'"does not exist'//{introspection:location})
             end if
          case default
             call Error_Report('unsupported property type'//{introspection:location})
          end select
       end do
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
  end subroutine filterPropertiesOperate
