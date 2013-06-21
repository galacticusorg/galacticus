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

!% Contains a module which provides various internal option codes for the galactic structure functions.

module Galactic_Structure_Options
  !% Provides various internal option codes for the galactic structure functions.
  implicit none
  public

  ! Values used to represent different mass types.
  !@ <enumeration>
  !@  <name>massType</name>
  !@  <description>Used to specify the mass type(s) to be queried in galactic structure functions.</description>
  !@  <entry label="massTypeUnknown"   />
  !@  <entry label="massTypeAll"       />
  !@  <entry label="massTypeDark"      />
  !@  <entry label="massTypeBaryonic"  />
  !@  <entry label="massTypeGalactic"  />
  !@  <entry label="massTypeGaseous"   />
  !@  <entry label="massTypeStellar"   />
  !@  <entry label="massTypeBlackHole" />
  !@ </enumeration>
  integer, parameter :: massTypeCount    =7  
  integer, parameter :: massTypeUnknown  =-1 
  integer, parameter :: massTypeAll      =0  
  integer, parameter :: massTypeDark     =1  
  integer, parameter :: massTypeBaryonic =2  
  integer, parameter :: massTypeGalactic =3  
  integer, parameter :: massTypeGaseous  =4  
  integer, parameter :: massTypeStellar  =5  
  integer, parameter :: massTypeBlackHole=6  
  character(len=9), parameter, dimension(0:massTypeCount-1) :: massTypesName    = &
       &                                                        [                 &
       &                                                         'all      ',     &
       &                                                         'dark     ',     &
       &                                                         'baryonic ',     &
       &                                                         'galactic ',     &
       &                                                         'gaseous  ',     &
       &                                                         'stellar  ',     &
       &                                                         'blackHole'      &
       &                                                        ]
  
  ! Values used to represent different component types.
  !@ <enumeration>
  !@  <name>componentType</name>
  !@  <description>Used to specify the component(s) to be queried in galactic structure functions.</description>
  !@  <entry label="componentTypeUnknown"   />
  !@  <entry label="componentTypeAll"       />
  !@  <entry label="componentTypeDisk"      />
  !@  <entry label="componentTypeSpheroid"  />
  !@  <entry label="componentTypeHotHalo"   />
  !@  <entry label="componentTypeDarkHalo"  />
  !@  <entry label="componentTypeBlackHole" />
  !@ </enumeration>
  integer, parameter :: componentTypeCount    =6  
  integer, parameter :: componentTypeUnknown  =-1 
  integer, parameter :: componentTypeAll      =0  
  integer, parameter :: componentTypeDisk     =1  
  integer, parameter :: componentTypeSpheroid =2  
  integer, parameter :: componentTypeHotHalo  =3  
  integer, parameter :: componentTypeDarkHalo =4  
  integer, parameter :: componentTypeBlackHole=5  
  character(len=9), parameter, dimension(0:componentTypeCount-1) :: componentTypesName    = &
       &                                                        [                           &
       &                                                         'all      ',               &
       &                                                         'disk     ',               &
       &                                                         'spheroid ',               &
       &                                                         'hotHalo  ',               &
       &                                                         'darkHalo ',               &
       &                                                         'blackHole'                &
       &                                                         ]

  ! Coordinate system options.
  !@ <enumeration>
  !@  <name>coordinateSystem</name>
  !@  <description>Used to specify the coordinate system of the input coordinates in galactic structure functions.</description>
  !@  <entry label="coordinateSystemSpherical"   />
  !@  <entry label="coordinateSystemCylindrical" />
  !@  <entry label="coordinateSystemCartesian"   />
  !@ </enumeration>
  integer         , parameter :: coordinateSystemSpherical  =1      
  integer         , parameter :: coordinateSystemCylindrical=2      
  integer         , parameter :: coordinateSystemCartesian  =3      
  
  ! Weighting options.
  !@ <enumeration>
  !@  <name>weightBy</name>
  !@  <description>Used to specify by which quantity to weight the results in galactic structure functions.</description>
  !@  <entry label="weightByMass"       />
  !@  <entry label="weightByLuminosity" />
  !@ </enumeration>
  integer         , parameter :: weightByMass               =0      
  integer         , parameter :: weightByLuminosity         =1      
  integer         , parameter :: weightIndexNull            =0      
  
  ! Suitably large value to represent infinite radius.
  double precision, parameter :: radiusLarge                =1.0d10 
  
contains

  integer function Galactic_Structure_Mass_Type_Decode(massTypeName)
    !% Decode a mass type from a string, returning the appropriate identifier.
    use Galacticus_Error
    implicit none
    character(len=*), intent(in   ) :: massTypeName 
    integer                         :: i            
    
    Galactic_Structure_Mass_Type_Decode=massTypeUnknown
    do i=0,massTypeCount-1
       if (trim(massTypeName) == trim(massTypesName(i))) then
          Galactic_Structure_Mass_Type_Decode=i
          exit
       end if
    end do
    if (Galactic_Structure_Mass_Type_Decode == massTypeUnknown) &
         & call Galacticus_Error_Report('Galactic_Structure_Mass_Type_Decode','unrecognized mass specifier ['//trim(massTypeName)//']')
    return
  end function Galactic_Structure_Mass_Type_Decode

  integer function Galactic_Structure_Component_Type_Decode(componentTypeName)
    !% Decode a component type from a string, returning the appropriate identifier.
    use Galacticus_Error
    implicit none
    character(len=*), intent(in   ) :: componentTypeName 
    integer                         :: i                 
    
    Galactic_Structure_Component_Type_Decode=componentTypeUnknown
    do i=0,componentTypeCount-1
       if (trim(componentTypeName) == trim(componentTypesName(i))) then
          Galactic_Structure_Component_Type_Decode=i
          exit
       end if
    end do
    if (Galactic_Structure_Component_Type_Decode == componentTypeUnknown) &
         & call Galacticus_Error_Report('Galactic_Structure_Component_Type_Decode','unrecognized component specifier ['//trim(componentTypeName)//']')
    return
  end function Galactic_Structure_Component_Type_Decode

end module Galactic_Structure_Options
