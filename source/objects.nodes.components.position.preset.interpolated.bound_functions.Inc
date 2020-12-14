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

!% Contains custom functions for the preset position component.

function PositionPresetInterpolatedPosition(self) result(position)
  !% Return the position of the node.
  use :: Node_Component_Position_Preset_Interpolated_Data, only : cosmologyFunctions_, positionPresetBoxLength, isPeriodic
  use :: Numerical_Interpolation                         , only : interpolator
  implicit none
  double precision                                         , allocatable  , dimension(:     ) :: position
  class           (nodeComponentPositionPresetInterpolated), intent(inout)                    :: self
  type            (treeNode                               ), pointer                          :: node
  class           (nodeComponentBasic                     ), pointer                          :: basic                , basicParent
  class           (nodeComponentPosition                  ), pointer                          :: positionParent
  double precision                                                        , dimension(4,3   ) :: coefficientsCubic
  double precision                                                        , dimension(    20) :: coefficientsSpiral
  double precision                                                        , dimension(  2, 2) :: coefficientsAngle    , coefficientsLogRadius
  double precision                                                        , dimension(2,2, 3) :: vectorInPlaneNormal
  integer                                                                 , parameter         :: historyBegin       =1, historyEnd           =3
  double precision                                                                            :: time                 , angle                  , &
       &                                                                                         logRadius            , timeParent             , &
       &                                                                                         boxLength
  integer         (c_size_t                               )                                   :: i
  logical                                                                                     :: usingHistory
  type            (interpolator                           )                                   :: interpolator_
  type            (history                                )                                   :: history_
  
  allocate(position(3))
  node  => self %host ()
  basic => node %basic()
  time  =  basic%time ()
  ! Use the size of the interpolation coefficient array to determine the type of interpolation. A more elegant/robust solution
  ! would be preferable.
  if (size(self%interpolationCoefficientsData) == 0) then
     ! Interpolation coefficients have not yet been computed or are not defined. Use the non-interpolated position.
     usingHistory=.false.
     history_    =self%positionHistory()
     if (history_%exists()) then
        if (history_%time(1) <= basic%time()) then
           interpolator_=interpolator        (history_%time                                  )
           i            =interpolator_%locate(basic   %time(),closest=.true.                 )
           position     =history_     %data  (i              ,        historyBegin:historyEnd)
           usingHistory =.true.
        end if
     end if
     if (.not.usingHistory) position=self%positionData
  else if (size(self%interpolationCoefficientsData) == 12) then
     ! Using cubic polynomial interpolation.
     coefficientsCubic=reshape(self%interpolationCoefficients(),[4,3])
     do i=1,3
        position(i)=+coefficientsCubic(1,i)*time**3 &
             &      +coefficientsCubic(2,i)*time**2 &
             &      +coefficientsCubic(3,i)*time    &
             &      +coefficientsCubic(4,i)
     end do
     ! Convert from comoving back to physical position.
     position=position*cosmologyFunctions_%expansionFactor(time)     
  else if (size(self%interpolationCoefficientsData) == 20) then
     ! Use logarithmic spiral interoplation in physical position.
     coefficientsSpiral   = self%interpolationCoefficients()
     vectorInPlaneNormal  = reshape(coefficientsSpiral( 1:12),[2,2,3])
     coefficientsAngle    = reshape(coefficientsSpiral(13:16),[2,2  ])
     coefficientsLogRadius= reshape(coefficientsSpiral(17:20),[2,2  ])
     angle                = coefficientsAngle    (1,1)+coefficientsAngle    (1,2)*time
     logRadius            = coefficientsLogRadius(1,1)+coefficientsLogRadius(1,2)*time
     position             =+(                                       &
          &                  +vectorInPlaneNormal(1,1,:)*cos(angle) &
          &                  +vectorInPlaneNormal(1,2,:)*sin(angle) &
          &                 )                                       &
          &                *exp(logRadius)
     ! Our interpolation is relative to the host halo center - add that on now.
     basicParent    => node       %parent%basic   ()
     positionParent => node       %parent%position()
     timeParent     =  basicParent%time           ()
     call basicParent%timeSet(time      )
     position=+               position   &
          &   +positionParent%position()
     call basicParent%timeSet(timeParent)     
  end if
  ! Handle periodic boundaries.
  if (isPeriodic) then
     boxLength=positionPresetBoxLength*cosmologyFunctions_%expansionFactor(time)
     do i=1,3
        do while (position(i) <  0.0d0    )
           position(i)=position(i)+boxLength
        end do
        do while (position(i) >= boxLength)
           position(i)=position(i)-boxLength
        end do
     end do
  end if
  return
end function PositionPresetInterpolatedPosition

function PositionPresetInterpolatedVelocity(self) result (velocity)
  !% Return the velocity of the node.
  use :: Node_Component_Position_Preset_Interpolated_Data, only : cosmologyFunctions_
  use :: Numerical_Interpolation                         , only : interpolator
  implicit none
  double precision                                         , allocatable  , dimension(:     ) :: velocity
  class           (nodeComponentPositionPresetInterpolated), intent(inout)                    :: self
  type            (treeNode                               ), pointer                          :: node
  class           (nodeComponentBasic                     ), pointer                          :: basic                , basicParent
  class           (nodeComponentPosition                  ), pointer                          :: positionParent
  double precision                                                        , dimension(4,3   ) :: coefficientsCubic
  double precision                                                        , dimension(    20) :: coefficientsSpiral
  double precision                                                        , dimension(  2, 2) :: coefficientsAngle    , coefficientsLogRadius
  double precision                                                        , dimension(2,2, 3) :: vectorInPlaneNormal
  integer                                                                 , parameter         :: historyBegin       =1, historyEnd           =3
  double precision                                                                            :: time                 , angle                  , &
       &                                                                                         logRadius            , timeParent
  integer         (c_size_t                               )                                   :: i
  logical                                                                                     :: usingHistory
  type            (interpolator                           )                                   :: interpolator_
  type            (history                                )                                   :: history_
  
  allocate(velocity(3))
  node  => self %host ()
  basic => node %basic()
  time  =  basic%time ()
  ! Use the size of the interpolation coefficient array to determine the type of interpolation. A more elegant/robust solution
  ! would be preferable.
  if (size(self%interpolationCoefficientsData) == 0) then
     ! Interpolation coefficients have not yet been computed or are not defined. Use the non-interpolated velocity.
     usingHistory=.false.
     history_    =self%positionHistory()
     if (history_%exists()) then
        if (history_%time(1) <= basic%time()) then
           interpolator_=interpolator        (history_%time                                  )
           i            =interpolator_%locate(basic   %time(),closest=.true.                 )
           velocity     =history_     %data  (i              ,        historyBegin:historyEnd)
           usingHistory =.true.
        end if
     end if
     if (.not.usingHistory) velocity=self%velocityData
  else if (size(self%interpolationCoefficientsData) == 12) then
     ! Using cubic polynomial interpolation.
     coefficientsCubic=reshape(self%interpolationCoefficients(),[4,3])
     do i=1,3
        velocity(i)=+3.0d0*coefficientsCubic(1,i)*time**2 &
             &      +2.0d0*coefficientsCubic(2,i)*time    &
             &      +      coefficientsCubic(3,i)
     end do
     ! Convert from comoving back to physical velocity.
     velocity=velocity*cosmologyFunctions_%expansionFactor(time)     
  else if (size(self%interpolationCoefficientsData) == 20) then
     ! Use logarithmic spiral interoplation in physical velocity.
     coefficientsSpiral   = self%interpolationCoefficients()
     vectorInPlaneNormal  = reshape(coefficientsSpiral( 1:12),[2,2,3])
     coefficientsAngle    = reshape(coefficientsSpiral(13:16),[2,2  ])
     coefficientsLogRadius= reshape(coefficientsSpiral(17:20),[2,2  ])
     angle                = coefficientsAngle    (2,1)+coefficientsAngle    (2,2)*time
     logRadius            = coefficientsLogRadius(2,1)+coefficientsLogRadius(2,2)*time
     velocity             =+(                                       &
          &                  +vectorInPlaneNormal(2,1,:)*cos(angle) &
          &                  +vectorInPlaneNormal(2,2,:)*sin(angle) &
          &                 )                                       &
          &                *exp(logRadius)
     ! Our interpolation is relative to the host halo center - add that on now.
     basicParent    => node       %parent%basic   ()
     positionParent => node       %parent%position()
     timeParent     =  basicParent%time           ()
     call basicParent%timeSet(time      )
     velocity=+               velocity   &
          &   +positionParent%velocity()
     call basicParent%timeSet(timeParent)     
  end if
  return
end function PositionPresetInterpolatedVelocity
