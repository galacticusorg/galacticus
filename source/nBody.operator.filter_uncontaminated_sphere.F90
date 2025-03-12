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
Implements an N-body data operator which filters particles to select those within a sphere such that the
contamination by particles of non-preferred type is below a specified level.
!!}
  
  !![
  <nbodyOperator name="nbodyOperatorFilterUncontaminatedSphere">
    <description>
      An N-body data operator which filters particles to select those within a sphere such that the contamination by particles of
      non-preferred type is below a specified level.
    </description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorFilterUncontaminatedSphere
     !!{
     An N-body data operator which filters particles to select those within a sphere such that the contamination by particles of
     non-preferred type is below a specified level.
     !!}
     private
     double precision, dimension(3) :: point
     double precision               :: fractionContamination
     integer                        :: particleType
     logical                        :: massWeighted
   contains
     procedure :: operate => filterUncontaminatedSphereOperate
  end type nbodyOperatorFilterUncontaminatedSphere

  interface nbodyOperatorFilterUncontaminatedSphere
     !!{
     Constructors for the {\normalfont \ttfamily filterUncontaminatedSphere} N-body operator class.
     !!}
     module procedure filterUncontaminatedSphereConstructorParameters
     module procedure filterUncontaminatedSphereConstructorInternal
  end interface nbodyOperatorFilterUncontaminatedSphere

contains

  function filterUncontaminatedSphereConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily filterUncontaminatedSphere} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nbodyOperatorFilterUncontaminatedSphere)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    double precision                                         , dimension(3)  :: point
    double precision                                                         :: fractionContamination
    integer                                                                  :: particleType
    logical                                                                  :: massWeighted

    !![
    <inputParameter>
      <name>point</name>
      <source>parameters</source>
      <description>The point at which to center the filtered sphere.</description>
    </inputParameter>
    <inputParameter>
      <name>fractionContamination</name>
      <source>parameters</source>
      <description>The contamination fraction allowed within the sphere.</description>
    </inputParameter>
    <inputParameter>
      <name>particleType</name>
      <source>parameters</source>
      <description>The preferred particle type for filtering.</description>
    </inputParameter>
    <inputParameter>
      <name>massWeighted</name>
      <source>parameters</source>
      <description>If true, contamination is weighted by particle mass, otherwise it is unweighted (i.e. depends on the number of particles only).</description>
    </inputParameter>
    !!]
    self=nbodyOperatorFilterUncontaminatedSphere(point,fractionContamination,particleType,massWeighted)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function filterUncontaminatedSphereConstructorParameters

  function filterUncontaminatedSphereConstructorInternal(point,fractionContamination,particleType,massWeighted) result (self)
    !!{
    Internal constructor for the {\normalfont \ttfamily filterUncontaminatedSphere} N-body operator class.
    !!}
    implicit none
    type            (nbodyOperatorFilterUncontaminatedSphere)                              :: self
    double precision                                         , intent(in   ), dimension(3) :: point
    double precision                                         , intent(in   )               :: fractionContamination
    integer                                                  , intent(in   )               :: particleType
    logical                                                  , intent(in   )               :: massWeighted
    !![
    <constructorAssign variables="point, fractionContamination, particleType, massWeighted"/>
    !!]

    return
  end function filterUncontaminatedSphereConstructorInternal

  subroutine filterUncontaminatedSphereOperate(self,simulations)
    !!{
    Filter particles outside of an uncontaminated region.
    !!}
    use :: Sorting, only : sortIndex
    use :: Display, only : displayIndent, displayMessage, displayUnindent, verbosityLevelStandard
    implicit none
    class           (nbodyOperatorFilterUncontaminatedSphere), intent(inout)                 :: self
    type            (nBodyData                              ), intent(inout), dimension(  :) :: simulations
    double precision                                         , pointer      , dimension(:,:) :: position
    double precision                                         , pointer      , dimension(:  ) :: distanceFromPoint   , mass
    logical                                                  , allocatable  , dimension(  :) :: mask
    integer         (c_size_t                               ), allocatable  , dimension(  :) :: order
    integer         (c_size_t                               ), pointer      , dimension(  :) :: propertyInteger     , propertyIntegerFiltered
    double precision                                         , pointer      , dimension(  :) :: propertyReal        , propertyRealFiltered
    integer         (c_size_t                               ), pointer      , dimension(:,:) :: propertyIntegerRank1, propertyIntegerRank1Filtered
    double precision                                         , pointer      , dimension(:,:) :: propertyRealRank1   , propertyRealRank1Filtered
    double precision                                                                         :: distanceMaximum     , massAll                     , &
         &                                                                                      massContaminated
    integer                                                                                  :: i                   , j                           , &
         &                                                                                      k
    integer         (c_size_t                               )                                :: countFiltered       , countContaminated           , &
         &                                                                                      countAll
    character       (len=16                                 )                                :: label
    
    call displayIndent('filter on contamination',verbosityLevelStandard)
    do i=1,size(simulations)
       ! Retrieve required properties.
       position => simulations(i)%propertiesRealRank1%value('position')
       ! Allocate workspace.
       allocate(distanceFromPoint(size(position,dim=2)))
       allocate(order            (size(position,dim=2)))
       ! Compute distances.
       !$omp workshare
       distanceFromPoint=sqrt(                                  &
            &                 +(position(1,:)-self%point(1))**2 &
            &                 +(position(2,:)-self%point(2))**2 &
            &                 +(position(3,:)-self%point(3))**2 &
            &                )
       !$omp end workshare
       nullify(position)
       ! Get order index into distances.
       order=sortIndex(distanceFromPoint)
       ! Get particle type data.
       propertyInteger => simulations(i)%propertiesInteger%value('particleType')
       ! Find maximum distance which maintains required contamination fraction.
       distanceMaximum  =distanceFromPoint(order(size(order)))
       if (self%massWeighted) then
          mass             => simulations(i)%propertiesReal%value('massParticle')
          massAll          =  0_c_size_t
          massContaminated =  0_c_size_t
          do j=1,size(order)
             massAll=massAll+mass(order(j))
             if (propertyInteger(order(j)) /= self%particleType) massContaminated=massContaminated+mass(order(j))
             if (massContaminated > self%fractionContamination*massAll) then
                distanceMaximum=distanceFromPoint(order(j))
                exit
             end if
          end do
          nullify(mass)
       else
          countAll         =0_c_size_t
          countContaminated=0_c_size_t
          do j=1,size(order)
             countAll=countAll+1_c_size_t
             if (propertyInteger(order(j)) /= self%particleType) countContaminated=countContaminated+1_c_size_t
             if (countContaminated > int(self%fractionContamination*dble(countAll),kind=c_size_t)) then
                distanceMaximum=distanceFromPoint(order(j))
                exit
             end if
          end do
       end if
       ! Construct the filter mask.
       allocate(mask(size(propertyInteger)))
       !$omp workshare
       mask=distanceFromPoint < distanceMaximum
       !$omp end workshare
       countFiltered=count(mask)
       nullify   (propertyInteger  )
       deallocate(distanceFromPoint)
       deallocate(order            )
       write (label,'(e12.6)') distanceMaximum
       call displayMessage('filtering particles within a radius of '//trim(adjustl(label))//' Mpc',verbosityLevelStandard)
       call simulations(i)%attributesReal%set(keyCH='radiusUncontaminated',value=distanceMaximum)
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
  end subroutine filterUncontaminatedSphereOperate
