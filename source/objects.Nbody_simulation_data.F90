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
  Contains a module which provides a class to store N-body simulation data.
  !!}

module NBody_Simulation_Data
  !!{
  Provides a class to store N-body simulation data.
  !!}
  use :: IO_HDF5           , only : hdf5Object
  use :: ISO_Varying_String, only : varying_string
  use :: Hashes            , only : rank1IntegerSizeTPtrHash, rank2IntegerSizeTPtrHash, rank1DoublePtrHash, rank2DoublePtrHash, &
       &                            integerSizeTHash        , doubleHash              , varyingStringHash , genericHash
  implicit none
  private
  public :: nBodyData, nBodyDataPropertyType

  type :: nBodyData
     !!{
     A class to store N-body simulation data.
     !!}
     type(varying_string          ) :: label
     type(hdf5Object              ) :: analysis
     type(integerSizeTHash        ) :: attributesInteger
     type(doubleHash              ) :: attributesReal
     type(varyingStringHash       ) :: attributesText
     type(genericHash             ) :: attributesGeneric
     type(rank1IntegerSizeTPtrHash) :: propertiesInteger
     type(rank1DoublePtrHash      ) :: propertiesReal
     type(rank2IntegerSizeTPtrHash) :: propertiesIntegerRank1
     type(rank2DoublePtrHash      ) :: propertiesRealRank1
   contains
     final :: nBodyDataDestructorScalar, nBodyDataDestructorRank1
  end type nBodyData

  !![
  <enumeration>
   <name>propertyType</name>
   <description>Enumeration of property types for N-body data properties.</description>
   <visibility>public</visibility>
   <entry label="unknown"/>
   <entry label="integer"/>
   <entry label="real"   />
  </enumeration>
  !!]
  
contains

  subroutine nBodyDataDestructorScalar(self)
    !!{
    Destruct for scalar {\normalfont \ttfamily nBodyData} objects.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    use iso_varying_string
    implicit none
    type            (nBodyData), intent(inout)                 :: self
    integer         (c_size_t ), pointer      , dimension(:  ) :: propertyInteger
    integer         (c_size_t ), pointer      , dimension(:,:) :: propertyIntegerRank1
    double precision           , pointer      , dimension(:  ) :: propertyReal
    double precision           , pointer      , dimension(:,:) :: propertyRealRank1
    class           (*        ), pointer                       :: attributeGeneric
    integer                                                    :: i

    do i=1,self%propertiesInteger%size()
       propertyInteger      => self%propertiesInteger     %value(i)
       deallocate(propertyInteger     )
    end do
    do i=1,self%propertiesIntegerRank1%size()
       propertyIntegerRank1 => self%propertiesIntegerRank1%value(i)
       deallocate(propertyIntegerRank1)
    end do
    do i=1,self%propertiesReal%size()
       propertyReal         => self%propertiesReal        %value(i)
       deallocate(propertyReal        )
    end do
    do i=1,self%propertiesRealRank1%size()
       propertyRealRank1    => self%propertiesRealRank1   %value(i)
       deallocate(propertyRealRank1   )
    end do       
    do i=1,self%attributesGeneric%size()
       attributeGeneric     => self%attributesGeneric     %value(i)
       deallocate(attributeGeneric    )
    end do       
    return
  end subroutine nBodyDataDestructorScalar
  
  subroutine nBodyDataDestructorRank1(self)
    !!{
    Destruct for rank-1 {\normalfont \ttfamily nBodyData} objects.
    !!}
    implicit none
    type   (nBodyData), intent(inout), dimension(:) :: self
    integer                                         :: i

    do i=1,size(self)
       call nBodyDataDestructorScalar(self(i))
    end do
    return
  end subroutine nBodyDataDestructorRank1
  
  function nBodyDataPropertyType(propertyName)
    !!{
    Returns the type of the named property.
    !!}
    implicit none
    type     (enumerationPropertyTypeType)                :: nBodyDataPropertyType
    character(len=*                      ), intent(in   ) :: propertyName

    select case (propertyName)
    case('particleID'               )
       nBodyDataPropertyType=propertyTypeInteger
    case('descendantID'             )
       nBodyDataPropertyType=propertyTypeInteger
    case('progenitorCount'          )
       nBodyDataPropertyType=propertyTypeInteger
    case('hostID'                   )
       nBodyDataPropertyType=propertyTypeInteger
    case('isolatedHostID'           )
       nBodyDataPropertyType=propertyTypeInteger
    case('descendantHostID'         )
       nBodyDataPropertyType=propertyTypeInteger
    case('isPhantom'                )
       nBodyDataPropertyType=propertyTypeInteger
    case('alwaysIsolated'           )
       nBodyDataPropertyType=propertyTypeInteger
    case('isFlyby'                  )
       nBodyDataPropertyType=propertyTypeInteger
    case('expansionFactor'          )
       nBodyDataPropertyType=propertyTypeReal
    case('descendantExpansionFactor')
       nBodyDataPropertyType=propertyTypeReal
    case('massVirial'               )
       nBodyDataPropertyType=propertyTypeReal
    case('radiusVirial'             )
       nBodyDataPropertyType=propertyTypeReal
    case('radiusScale'              )
       nBodyDataPropertyType=propertyTypeReal
    case('spin'                     )
       nBodyDataPropertyType=propertyTypeReal
    case('virialRatio'              )
       nBodyDataPropertyType=propertyTypeReal
    case('distanceFromPoint'        )
       nBodyDataPropertyType=propertyTypeReal
    case default
       nBodyDataPropertyType=propertyTypeUnknown
    end select
    return
  end function nBodyDataPropertyType

end module NBody_Simulation_Data
