!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a black hole binary initial separation which is a fixed fraction of the scale radius of the
!% larger of the host and satellite spheroids.

module Black_Hole_Binary_Initial_Radii_Spheroid_Size
  !% Implements a black hole binary initial separation which is a fixed fraction of the scale radius of the larger of the host and
  !% satellite spheroids.
  implicit none
  private
  public :: Black_Hole_Binary_Initial_Radii_Spheroid_Size_Initialize

  ! The fraction of the spheroid radius at which merging black holes should be placed.
  double precision :: blackHoleInitialRadiusSpheroidRadiusRatio  
                                                              
contains

  !# <blackHoleBinaryInitialRadiiMethod>
  !#  <unitName>Black_Hole_Binary_Initial_Radii_Spheroid_Size_Initialize</unitName>
  !# </blackHoleBinaryInitialRadiiMethod>
  subroutine Black_Hole_Binary_Initial_Radii_Spheroid_Size_Initialize(blackHoleBinaryInitialRadiiMethod,Black_Hole_Binary_Initial_Radius_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string                                ), intent(in   )          :: blackHoleBinaryInitialRadiiMethod     
    procedure(Black_Hole_Binary_Initial_Radius_Spheroid_Size), intent(inout), pointer :: Black_Hole_Binary_Initial_Radius_Get  
                                                                                                                            
    if (blackHoleBinaryInitialRadiiMethod == 'spheroidRadiusFraction') then
       Black_Hole_Binary_Initial_Radius_Get => Black_Hole_Binary_Initial_Radius_Spheroid_Size
       !@ <inputParameter>
       !@   <name>blackHoleInitialRadiusSpheroidRadiusRatio</name>
       !@   <defaultValue>0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The fraction of the spheroid radius at which merging black holes will be initially placed.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('blackHoleInitialRadiusSpheroidRadiusRatio',blackHoleInitialRadiusSpheroidRadiusRatio,defaultValue=0.0d0)
    end if
    return
  end subroutine Black_Hole_Binary_Initial_Radii_Spheroid_Size_Initialize

  double precision function Black_Hole_Binary_Initial_Radius_Spheroid_Size(thisNode,hostNode)
    !% Returns an initial separation for a binary black holes that is a fixed fraction of the scale radius of the larger of the
    !% host and satellite spheroids.
    use Galacticus_Nodes
    implicit none
    type (treeNode             ), intent(inout), pointer :: hostNode             , thisNode               
    class(nodeComponentSpheroid)               , pointer :: hostSpheroidComponent, thisSpheroidComponent  
    
    ! Get the spheroid components.                                                                                                   
    thisSpheroidComponent => thisNode%spheroid()
    hostSpheroidComponent => hostNode%spheroid()
    ! Compute the initial radius.
    Black_Hole_Binary_Initial_Radius_Spheroid_Size=blackHoleInitialRadiusSpheroidRadiusRatio*max(                                &
         &                                                                                       thisSpheroidComponent%radius(), &
         &                                                                                       hostSpheroidComponent%radius()  &
         &                                                                                      )
    return
  end function Black_Hole_Binary_Initial_Radius_Spheroid_Size

end module Black_Hole_Binary_Initial_Radii_Spheroid_Size
