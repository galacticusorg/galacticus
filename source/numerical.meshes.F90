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
Contains a module which provieds tools for working with grids.
!!}

module Meshes
  !!{
  Provide tools for working with grids.
  !!}
  private
  public :: Meshes_Apply_Point

  ! Cloud types.
  integer, parameter, public :: cloudTypePoint     =0
  integer, parameter, public :: cloudTypeCubic     =1
  integer, parameter, public :: cloudTypeTriangular=2

contains

  subroutine Meshes_Apply_Point(mesh,boxLength,pointPosition,pointWeight,cloudType)
    !!{
    Apply a point to a mesh.
    !!}
    use            :: Error        , only : Error_Report
    use, intrinsic :: ISO_C_Binding, only : c_double_complex
    implicit none
    complex(c_double_complex), intent(inout), dimension(:,:,:) :: mesh
    double precision                  , intent(in   ) :: boxLength
    double precision, dimension(3    ), intent(in   ) :: pointPosition
    complex(c_double_complex), intent(in   ), optional         :: pointWeight
    integer         , intent(in   )   , optional :: cloudType
    integer         , dimension(3    )           :: gridIndices
    double precision, dimension(3    )           :: gridPosition
    integer                                      :: cloudTypeActual, i, j, k
    complex(c_double_complex)                                  :: pointWeightActual,meshFraction

    ! Find the particle weight.
    if (present(pointWeight)) then
       pointWeightActual=pointWeight
    else
       pointWeightActual=1.0d0
    end if

     ! Find the cloud type to use.
     if (present(cloudType)) then
        cloudTypeActual=cloudType
     else
        cloudTypeActual=cloudTypePoint
     end if

    ! Compute the position in the grid.
    gridPosition=pointPosition*dble(shape(mesh))/boxLength+1.0d0
    gridIndices =int(gridPosition)

    ! Apply to the grid.
    select case (cloudTypeActual)
    case (cloudTypePoint)
       if (all(gridIndices > 0 .and. gridIndices <= shape(mesh)))  &
            &  mesh(gridIndices(1),gridIndices(2),gridIndices(3))  &
            & =mesh(gridIndices(1),gridIndices(2),gridIndices(3))  &
            & +pointWeightActual
    case (cloudTypeCubic)
       do i      =-1,1
          do j   =-1,1
             do k=-1,1
                meshFraction=product(                                                                    &
                     &               +min(max(gridPosition-dble(gridIndices+[i,j,k])+0.5d0,0.0d0),1.0d0) &
                     &               -min(max(gridPosition-dble(gridIndices+[i,j,k])-0.5d0,0.0d0),1.0d0) &
                     &              )
                if (all(gridIndices+[i,j,k] > 0 .and. gridIndices+[i,j,k] <= shape(mesh)))  &
                     &  mesh(gridIndices(1)+i,gridIndices(2)+j,gridIndices(3)+k)            &
                     & =mesh(gridIndices(1)+i,gridIndices(2)+j,gridIndices(3)+k)            &
                     & +pointWeightActual*meshFraction
             end do
          end do
       end do
    case (cloudTypeTriangular)
       do i      =-1,1
          do j   =-1,1
             do k=-1,1
                meshFraction=product(                                                                     &
                     &               Triangular_Shaped_Cloud_Integral(                                    &
                     &                +min(max(gridPosition-dble(gridIndices+[i,j,k])+0.5d0,0.0d0),1.0d0) &
                     &                -min(max(gridPosition-dble(gridIndices+[i,j,k])-0.5d0,0.0d0),1.0d0) &
                     &               )                                                                    &
                     &              )
                if (all(gridIndices+[i,j,k] > 0 .and. gridIndices+[i,j,k] <= shape(mesh)))  &
                     &  mesh(gridIndices(1)+i,gridIndices(2)+j,gridIndices(3)+k)            &
                     & =mesh(gridIndices(1)+i,gridIndices(2)+j,gridIndices(3)+k)            &
                     & +pointWeightActual*meshFraction
             end do
          end do
       end do
    case default
       call Error_Report('unrecognized cloud type'//{introspection:location})
    end select

    return
  end subroutine Meshes_Apply_Point

  elemental double precision function Triangular_Shaped_Cloud_Integral(cellFraction)
    !!{
    Return the integral over a triangular shaped cloud given the fraction of the cloud length in a cell.
    !!}
    implicit none
    double precision, intent(in   ) :: cellFraction

    if (cellFraction < 0.5d0) then
       Triangular_Shaped_Cloud_Integral=2.0d0*              cellFraction **2
    else
       Triangular_Shaped_Cloud_Integral=2.0d0*(1.0d0-(1.0d0-cellFraction)**2)-1.0d0
    end if
    return
  end function Triangular_Shaped_Cloud_Integral

end module Meshes
