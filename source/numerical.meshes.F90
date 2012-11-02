!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which provides tools for working with grids.

module Meshes
  !% Provide tools for working with grids.
  private
  public :: Meshes_Apply_Point

  ! Cloud types.
  integer, public, parameter :: cloudTypePoint     =0
  integer, public, parameter :: cloudTypeCubic     =1
  integer, public, parameter :: cloudTypeTriangular=2

contains

  subroutine Meshes_Apply_Point(mesh,boxLength,pointPosition,pointWeight,cloudType)
    !% Apply a point to a mesh.
    use, intrinsic :: ISO_C_Binding
    use Galacticus_Error
    implicit none
    complex(c_double_complex), intent(inout), dimension(:,:,:) :: mesh
    double precision         , intent(in   )                   :: boxLength
    double precision         , intent(in   ), dimension(3    ) :: pointPosition
    complex(c_double_complex), intent(in   ), optional         :: pointWeight
    integer                  , intent(in   ), optional         :: cloudType
    integer                  ,                dimension(3    ) :: gridIndices
    double precision         ,                dimension(3    ) :: gridPosition
    integer                                                    :: cloudTypeActual,i,j,k
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
       call Galacticus_Error_Report('Meshes_Apply_Point','unrecognized cloud type')
    end select

    return
  end subroutine Meshes_Apply_Point

  elemental double precision function Triangular_Shaped_Cloud_Integral(cellFraction)
    !% Return the integral over a triangular shaped cloud given the fraction of the cloud length in a cell.
    implicit none
    double precision, intent(in) :: cellFraction

    if (cellFraction < 0.5d0) then
       Triangular_Shaped_Cloud_Integral=2.0d0*              cellFraction **2
    else
       Triangular_Shaped_Cloud_Integral=2.0d0*(1.0d0-(1.0d0-cellFraction)**2)-1.0d0
    end if
    return
  end function Triangular_Shaped_Cloud_Integral

end module Meshes
