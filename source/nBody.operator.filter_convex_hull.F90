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
Implements an N-body data operator which filters particles outside of a convex hull.
!!}
  
  !![
  <nbodyOperator name="nbodyOperatorFilterConvexHull">
   <description>An N-body data operator which filters particles outside of a convex hull.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorFilterConvexHull
     !!{
     An N-body data operator which filters particles outside of a convex hull.
     !!}
     private
     integer :: hullFromSimulation
   contains
     procedure :: operate => filterConvexHullOperate
  end type nbodyOperatorFilterConvexHull

  interface nbodyOperatorFilterConvexHull
     !!{
     Constructors for the \refClass{nbodyOperatorFilterConvexHull} N-body operator class.
     !!}
     module procedure filterConvexHullConstructorParameters
     module procedure filterConvexHullConstructorInternal
  end interface nbodyOperatorFilterConvexHull

contains

  function filterConvexHullConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorFilterConvexHull} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nbodyOperatorFilterConvexHull)                :: self
    type   (inputParameters              ), intent(inout) :: parameters
    integer                                               :: hullFromSimulation
   
    !![
    <inputParameter>
      <name>hullFromSimulation</name>
      <source>parameters</source>
      <description>The index of the simulation from which to take the convex hull.</description>
    </inputParameter>
    !!]
    self=nbodyOperatorFilterConvexHull(hullFromSimulation)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function filterConvexHullConstructorParameters

  function filterConvexHullConstructorInternal(hullFromSimulation) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorFilterConvexHull} N-body operator class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nbodyOperatorFilterConvexHull)               :: self
    integer                               ,intent(in   ) :: hullFromSimulation
    !![
    <constructorAssign variables="hullFromSimulation"/>
    !!]
   
    return
  end function filterConvexHullConstructorInternal

  subroutine filterConvexHullOperate(self,simulations)
    !!{
    Filter particles outside of a convex hull.
    !!}
    use :: Display           , only : displayIndent , displayMessage     , displayUnindent, verbosityLevelStandard, &
         &                            displayCounter, displayCounterClear
    use :: Error             , only : Error_Report
    use :: Points_Convex_Hull, only : convexHull
    implicit none
    class           (nbodyOperatorFilterConvexHull), intent(inout)                 :: self
    type            (nBodyData                    ), intent(inout), dimension(  :) :: simulations
    logical                                        , allocatable  , dimension(  :) :: mask
    integer         (c_size_t                     ), pointer      , dimension(  :) :: propertyInteger     , propertyIntegerFiltered
    double precision                               , pointer      , dimension(  :) :: propertyReal        , propertyRealFiltered
    integer         (c_size_t                     ), pointer      , dimension(:,:) :: propertyIntegerRank1, propertyIntegerRank1Filtered
    double precision                               , pointer      , dimension(:,:) :: propertyRealRank1   , propertyRealRank1Filtered   , &
         &                                                                            position
    class           (*                            ), pointer                       :: hull
    integer                                                                        :: i                   , j                           , &
         &                                                                            k
    integer         (c_size_t                     )                                :: countFiltered
    
    call displayIndent('filter on convex hull',verbosityLevelStandard)
    hull => simulations(self%hullFromSimulation)%attributesGeneric%value('convexHull')
    select type (hull)
    class is (convexHull)
       do i=1,size(simulations)
          position => simulations(i)%propertiesRealRank1%value('position')
          allocate(mask(size(position,dim=2)))
          call displayCounter(0,.true.)
          do j=1,size(position,dim=2)
             if (mod(j,size(position,dim=2)/100) == 0)                                  &
                  & call displayCounter(                                                &
                  &                     int(                                            &
                  &                         +100.0d0                                    &
                  &                         *float(j                                 )  &
                  &                         /float(size(position,dim=2,kind=c_size_t))  &
                  &                        )                                          , &
                  &                     .false.                                         &
                  &                    )
             mask(j)=hull%pointIsInHull(position(:,j))
          end do
          call displayCounterClear()
          nullify(position)
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
          ! Set the convex hull volume for this simulation.
          call simulations(i)%attributesReal%set('convexHullVolume',hull%volume())
       end do
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    nullify(hull)
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine filterConvexHullOperate
  
