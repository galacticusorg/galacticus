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
Contains a module which provieds tools for working with sets of points.
!!}

module Points
  !!{
  Provide tools for working with sets of points.
  !!}
  private
  public :: Points_Prune, Points_Translate, Points_Rotate, Points_Replicate, Points_Survey_Geometry

contains

  subroutine Points_Prune(points,mask)
    !!{
    Prune a set of points.
    !!}
    implicit none
    double precision, allocatable, dimension(:,:), intent(inout) :: points
    logical                      , dimension(:  ), intent(in   ) :: mask
    double precision, allocatable, dimension(:,:)                :: pointsTmp
    integer                                                      :: pointsCount, i

    pointsCount=count(mask)
    call Move_Alloc(points,pointsTmp)
    allocate(points(3,pointsCount))
    pointsCount=0
    do i=1,size(pointsTmp,dim=2)
       if (mask(i)) then
          pointsCount=pointsCount+1
          points(:,pointsCount)=pointsTmp(:,i)
       end if
    end do
    deallocate(pointsTmp)
    return
  end subroutine Points_Prune

  subroutine Points_Translate(points,shift,periodicLength)
    !!{
    Apply a simple translation to a set of points.
    !!}
    implicit none
    double precision, dimension(:,:), intent(inout) :: points
    double precision, dimension(:  ), intent(in   ) :: shift
    double precision, optional      , intent(in   ) :: periodicLength
    integer                                         :: i             , j

    forall(i=1:size(points,dim=2))
       points(:,i)=points(:,i)+shift
    end forall
    if (present(periodicLength)) then
       do i=1,size(points,dim=2)
          do j=1,3
             do while (points(j,i) < 0.0d0)
                points(j,i)=points(j,i)+periodicLength
             end do
             do while (points(j,i) >= periodicLength)
                points(j,i)=points(j,i)-periodicLength
             end do
          end do
       end do
    end if
    return
  end subroutine Points_Translate

  subroutine Points_Replicate(points,periodicLength,replicantStart,replicantEnd)
    !!{
    Apply a simple translation to a set of points.
    !!}
    implicit none
    double precision, dimension(:,:), allocatable, intent(inout) :: points
    integer         , dimension(3  )             , intent(in   ) :: replicantStart  , replicantEnd
    double precision                             , intent(in   ) :: periodicLength
    double precision, dimension(:,:), allocatable                :: pointsTmp
    integer                                                      :: i               , j           , &
         &                                                          k               , l           , &
         &                                                          replicationCount, pointsCount

    pointsCount=size(points,dim=2)
    call move_alloc(points,  pointsTmp                                          )
    allocate       (points(3,pointsCount*product(replicantEnd-replicantStart+1)))
    replicationCount=0
    do i=replicantStart(1),replicantEnd(1)
       do j=replicantStart(2),replicantEnd(2)
          do k=replicantStart(3),replicantEnd(3)
             points(:,1+replicationCount*pointsCount:(replicationCount+1)*pointsCount)=pointsTmp
             forall(l=1:pointsCount)
                points(:,l+replicationCount*pointsCount)=points(:,l+replicationCount*pointsCount)+periodicLength*dble([i,j,k])
             end forall
             replicationCount=replicationCount+1
          end do
       end do
    end do
    deallocate(pointsTmp)
    return
  end subroutine Points_Replicate

  subroutine Points_Rotate(points,axis,angle)
    !!{
    Apply a rotation to a set of points.
    !!}
    use :: Vectors, only : Vector_Magnitude
    implicit none
    double precision, dimension(:,:), intent(inout) :: points
    double precision, dimension(3  ), intent(in   ) :: axis
    double precision                , intent(in   ) :: angle
    double precision, dimension(3  )                :: unitAxis      , point
    double precision, dimension(3,3)                :: rotationMatrix
    integer                                         :: i             , j

    ! Construct the unit axis vector.
    unitAxis=axis/Vector_Magnitude(axis)
    ! Construct the rotation matrix.
    rotationMatrix=                                                                    &
         & reshape(                                                                    &
         &         [                                                                   &
         &          unitAxis(1)*unitAxis(1)*(1.0d0-cos(angle))+            cos(angle), &
         &          unitAxis(1)*unitAxis(2)*(1.0d0-cos(angle))-unitAxis(3)*sin(angle), &
         &          unitAxis(1)*unitAxis(3)*(1.0d0-cos(angle))+unitAxis(2)*sin(angle), &
         &          unitAxis(2)*unitAxis(1)*(1.0d0-cos(angle))+unitAxis(3)*sin(angle), &
         &          unitAxis(2)*unitAxis(2)*(1.0d0-cos(angle))+            cos(angle), &
         &          unitAxis(2)*unitAxis(3)*(1.0d0-cos(angle))-unitAxis(1)*sin(angle), &
         &          unitAxis(3)*unitAxis(1)*(1.0d0-cos(angle))-unitAxis(2)*sin(angle), &
         &          unitAxis(3)*unitAxis(2)*(1.0d0-cos(angle))+unitAxis(1)*sin(angle), &
         &          unitAxis(3)*unitAxis(3)*(1.0d0-cos(angle))+            cos(angle)  &
         &         ]                                                                 , &
         &         [                                                                   &
         &          3                                                                , &
         &          3                                                                  &
         &         ]                                                                   &
         &        )
    ! Rotate each point.
    do i=1,size(points,dim=2)
       forall(j=1:3)
          point(j)=sum(rotationMatrix(j,:)*points(:,i))
       end forall
       points(:,i)=point
    end do
    return
  end subroutine Points_Rotate

  subroutine Points_Survey_Geometry(points,surveyGeometry_,mass)
    !!{
    Select a set of points that lie within a given survey geometry
    !!}
    use :: Display          , only : displayCounter     , displayCounterClear, displayIndent, displayUnindent
    use :: Geometry_Surveys , only : surveyGeometryClass
    implicit none
    double precision                     , allocatable, dimension(:,:), intent(inout) :: points
    class           (surveyGeometryClass)                             , intent(inout) :: surveyGeometry_
    double precision                                                  , intent(in   ) :: mass
    logical                              , allocatable, dimension(:  )                :: pointInclusionMask
    double precision                     , allocatable, dimension(:,:)                :: pointsTmp
    integer                                                                           :: iPoint            , pointInclusionCount, &
         &                                                                               pointTestCount

    allocate(pointInclusionMask(size(points,dim=2)))
    pointTestCount=0
    call displayIndent('Finding points in survey geometry')
    !$omp parallel do
    do iPoint=1,size(points,dim=2)
       !$omp atomic
       pointTestCount=pointTestCount+1
       call displayCounter(int(100.0d0*dble(pointTestCount-1)/dble(size(points,dim=2))),isNew=(pointTestCount==1))
       pointInclusionMask(iPoint)=surveyGeometry_%pointIncluded(points(:,iPoint),mass)
    end do
    !$omp end parallel do
    call displayCounterClear()
    call displayUnindent('done')
    pointInclusionCount=count(pointInclusionMask)
    call Move_Alloc(points,pointsTmp)
    allocate(points(3,pointInclusionCount))
    pointInclusionCount=0
    do iPoint=1,size(pointsTmp,dim=2)
       if (pointInclusionMask(iPoint)) then
          pointInclusionCount=pointInclusionCount+1
          points(:,pointInclusionCount)=pointsTmp(:,iPoint)
       end if
    end do
    deallocate(pointsTmp)
    return
  end subroutine Points_Survey_Geometry

end module Points
