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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a program to test multi-dimensional interpolation.
!!}

program Test_Interpolation_MultiD
  !!{RST
  Tests that the multi-dimensional interpolation code works correctly.

  Multilinear interpolation is exact for functions which are multilinear in the coordinates, so the primary tests
  build such a function, tabulate it, and check that interpolation recovers it exactly. Interpolation at the grid
  nodes of a function which is *not* multilinear is also checked, since that catches errors in the flattened
  indexing which a smooth test function could mask.
  !!}
  use            :: Display                       , only : displayVerbositySet, verbosityLevelStandard
  use, intrinsic :: ISO_C_Binding                 , only : c_size_t
  use            :: Numerical_Interpolation       , only : interpolator
  use            :: Numerical_Interpolation_MultiD, only : interpolatorMultiD
  use            :: Unit_Tests                    , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (interpolatorMultiD)                             :: interpolator2_                          , interpolator3_
  type            (interpolator      ), dimension(2)               :: interpolators2
  type            (interpolator      ), dimension(3)               :: interpolators3
  double precision                    , dimension(4)               :: x1            =[0.0d0,1.0d0,2.0d0,3.0d0]
  double precision                    , dimension(3)               :: x2            =[0.0d0,2.0d0,4.0d0]
  double precision                    , dimension(2)               :: y1            =[0.0d0,1.0d0]
  double precision                    , dimension(2)               :: y2            =[0.0d0,1.0d0]
  double precision                    , dimension(2)               :: y3            =[0.0d0,2.0d0]
  double precision                    , dimension(4,3)             :: values2
  double precision                    , dimension(2,2,2)           :: values3
  double precision                    , dimension(4,3)             :: valuesLumpy
  integer         (c_size_t          ), allocatable, dimension(:)  :: shape_
  integer         (c_size_t          ), dimension(4)               :: indices2
  double precision                    , dimension(4)               :: weights2
  integer                                                          :: i                                       , j             , &
       &                                                              k
  double precision                                                 :: y

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Multi-dimensional interpolation")

  ! Build a 2-D interpolator over a (4,3) grid.
  interpolators2(1)=interpolator(x1)
  interpolators2(2)=interpolator(x2)
  interpolator2_   =interpolatorMultiD(interpolators2)

  ! Tabulate a bilinear function, for which bilinear interpolation is exact, plus a "lumpy" function which is not
  ! bilinear and so only agrees at the grid nodes.
  do i=1,size(x1)
     do j=1,size(x2)
        values2    (i,j)=+1.0d0             &
             &           +2.0d0*x1(i)       &
             &           +3.0d0      *x2(j) &
             &           +4.0d0*x1(i)*x2(j)
        valuesLumpy(i,j)=+dble(i  )**3      &
             &           -dble(  j)**2      &
             &           +dble(i*j)
     end do
  end do

  ! Test the reported geometry of the interpolator.
  call Unit_Tests_Begin_Group("geometry")
  call Assert("dimension count",     interpolator2_%countDimensions() , 2   )
  call Assert("corner count"   ,     interpolator2_%countCorners   () , 4   )
  call Assert("table count"    ,int (interpolator2_%countTable     ()),12   )
  shape_=interpolator2_%shapeTable()
  call Assert("table shape"    ,int (shape_                          ),[4,3])
  call Unit_Tests_End_Group()

  ! Test that interpolation of a bilinear function is exact.
  call Unit_Tests_Begin_Group("bilinear function")
  y=interpolator2_%interpolate(values2,[1.5d0,3.0d0])
  call Assert("interior point"      ,y,31.0d0,relTol=1.0d-12)
  y=interpolator2_%interpolate(values2,[0.0d0,0.0d0])
  call Assert("lower corner"        ,y, 1.0d0,relTol=1.0d-12)
  y=interpolator2_%interpolate(values2,[3.0d0,4.0d0])
  call Assert("upper corner"        ,y,+1.0d0+2.0d0*3.0d0+3.0d0*4.0d0+4.0d0*3.0d0*4.0d0,relTol=1.0d-12)
  y=interpolator2_%interpolate(values2,[2.0d0,1.0d0])
  call Assert("on a grid line"      ,y,+1.0d0+2.0d0*2.0d0+3.0d0*1.0d0+4.0d0*2.0d0*1.0d0,relTol=1.0d-12)
  call Unit_Tests_End_Group()

  ! Test that interpolation at grid nodes recovers the tabulated values even for a function which is not bilinear.
  ! This verifies that the mapping from grid indices to positions in the flattened table is correct.
  call Unit_Tests_Begin_Group("grid nodes")
  do i=1,size(x1)
     do j=1,size(x2)
        y=interpolator2_%interpolate(valuesLumpy,[x1(i),x2(j)])
        call Assert("node value",y,valuesLumpy(i,j),relTol=1.0d-12)
     end do
  end do
  call Unit_Tests_End_Group()

  ! Test that the corner weights form a partition of unity.
  call Unit_Tests_Begin_Group("corner weights")
  call interpolator2_%corners([1.5d0,3.0d0],indices2,weights2)
  call Assert("weights sum to unity"       ,sum(weights2),1.0d0,relTol=1.0d-12)
  call Assert("corner indices within table",all(indices2 >= 1_c_size_t .and. indices2 <= interpolator2_%countTable()),.true.)
  call interpolator2_%corners([0.0d0,0.0d0],indices2,weights2)
  call Assert("weights sum to unity at node",sum(weights2),1.0d0,relTol=1.0d-12)
  call Unit_Tests_End_Group()

  ! Build a 3-D interpolator and check that interpolation of a trilinear function is exact.
  call Unit_Tests_Begin_Group("trilinear function")
  interpolators3(1)=interpolator(y1)
  interpolators3(2)=interpolator(y2)
  interpolators3(3)=interpolator(y3)
  interpolator3_   =interpolatorMultiD(interpolators3)
  do i=1,size(y1)
     do j=1,size(y2)
        do k=1,size(y3)
           values3(i,j,k)=+1.0d0                   &
                &         +      y1(i)             &
                &         +2.0d0*y2(j)             &
                &         +3.0d0*y3(k)             &
                &         +      y1(i)*y2(j)       &
                &         +      y1(i)*      y3(k) &
                &         +            y2(j)*y3(k) &
                &         +      y1(i)*y2(j)*y3(k)
        end do
     end do
  end do
  call Assert("corner count"  ,interpolator3_%countCorners(),8)
  y=interpolator3_%interpolate(values3,[0.5d0,0.5d0,1.0d0])
  call Assert("interior point",y,7.0d0,relTol=1.0d-12)
  y=interpolator3_%interpolate(values3,[1.0d0,1.0d0,2.0d0])
  call Assert("upper corner"  ,y,values3(2,2,2),relTol=1.0d-12)
  y=interpolator3_%interpolate(values3,[0.0d0,1.0d0,0.0d0])
  call Assert("mixed corner"  ,y,values3(1,2,1),relTol=1.0d-12)
  call Unit_Tests_End_Group()

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Interpolation_MultiD
