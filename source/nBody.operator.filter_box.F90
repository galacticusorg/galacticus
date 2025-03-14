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
Implements an N-body data operator which filters particles outside a cuboid region.
!!}
  
  !![
  <nbodyOperator name="nbodyOperatorFilterBox">
   <description>An N-body data operator which filters particles outside a cuboid region.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorFilterBox
     !!{
     An N-body data operator which filters particles outside a cuboid region.
     !!}
     private
     double precision, dimension(3) :: boundLow, boundHigh
   contains
     procedure :: operate => filterBoxOperate
  end type nbodyOperatorFilterBox

  interface nbodyOperatorFilterBox
     !!{
     Constructors for the {\normalfont \ttfamily filterBox} N-body operator class.
     !!}
     module procedure filterBoxConstructorParameters
     module procedure filterBoxConstructorInternal
  end interface nbodyOperatorFilterBox

contains

  function filterBoxConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily filterBox} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nbodyOperatorFilterBox)                :: self
    type            (inputParameters       ), intent(inout) :: parameters
    double precision                        , dimension(3)  :: boundLow  , boundHigh
    
    !![
    <inputParameter>
      <name>boundLow</name>
      <source>parameters</source>
      <description>The lower boundary of the cuboid region.</description>
    </inputParameter>
    <inputParameter>
      <name>boundHigh</name>
      <source>parameters</source>
      <description>The upper boundary of the cuboid region.</description>
    </inputParameter>
    !!]
    self=nbodyOperatorFilterBox(boundLow,boundHigh)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function filterBoxConstructorParameters

  function filterBoxConstructorInternal(boundLow,boundHigh) result (self)
    !!{
    Internal constructor for the {\normalfont \ttfamily filterBox} N-body operator class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (nbodyOperatorFilterBox)                              :: self
    double precision                        , intent(in   ), dimension(3) :: boundLow, boundHigh
    !![
    <constructorAssign variables="boundLow, boundHigh"/>
    !!]
    
    if (any(boundLow >= boundHigh))                                                          &
         & call Error_Report('filter will exclude all'//{introspection:location})
    return
  end function filterBoxConstructorInternal

  subroutine filterBoxOperate(self,simulations)
    !!{
    Filter particles outside of a cuboid region.
    !!}
    use :: Display, only : displayIndent, displayMessage  , displayUnindent, verbosityLevelStandard
    implicit none
    class           (nbodyOperatorFilterBox), intent(inout)                 :: self
    type            (nBodyData             ), intent(inout), dimension(  :) :: simulations
    logical                                 , allocatable  , dimension(  :) :: mask
    integer         (c_size_t              ), pointer      , dimension(  :) :: propertyInteger     , propertyIntegerFiltered
    double precision                        , pointer      , dimension(  :) :: propertyReal        , propertyRealFiltered
    integer         (c_size_t              ), pointer      , dimension(:,:) :: propertyIntegerRank1, propertyIntegerRank1Filtered
    double precision                        , pointer      , dimension(:,:) :: propertyRealRank1   , propertyRealRank1Filtered
    integer                                                                 :: i                   , j                           , &
         &                                                                     k
    integer         (c_size_t              )                                :: countFiltered
    
    call displayIndent('filter on cuboid region',verbosityLevelStandard)
    do i=1,size(simulations)
       propertyRealRank1 => simulations(i)%propertiesRealRank1%value('position')
       allocate(mask(size(propertyRealRank1,dim=2)))
       !$omp workshare
       mask   = propertyRealRank1(1,:) >= self%boundLow (1) &
            &  .and.                                        &
            &   propertyRealRank1(1,:) <  self%boundHigh(1) &
            &  .and.                                        &
            &   propertyRealRank1(2,:) >= self%boundLow (2) &
            &  .and.                                        &
            &   propertyRealRank1(2,:) <  self%boundHigh(2) &
            &  .and.                                        &
            &   propertyRealRank1(3,:) >= self%boundLow (3) &
            &  .and.                                        &
            &   propertyRealRank1(3,:) <  self%boundHigh(3)
       !$omp end workshare
       nullify(propertyRealRank1)
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
  end subroutine filterBoxOperate
