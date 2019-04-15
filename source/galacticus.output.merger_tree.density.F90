!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which handles outputting of density profile data to the \glc\ output file.

module Galacticus_Output_Trees_Density
  !% Handles outputting of density profile data to the \glc\ output file.
  use ISO_Varying_String
  use Galactic_Structure_Radii_Definitions
  implicit none
  private
  public :: Galacticus_Output_Tree_Density      , Galacticus_Output_Tree_Density_Property_Count, &
       &    Galacticus_Output_Tree_Density_Names

  ! Flag indicating whether or not projected density information is to be output.
  logical :: outputDensityData, outputDensityIncludeRadii

  ! Flag indicating whether or not this module has been initialized.
  logical :: outputDensityDataInitialized=.false.

  ! Radii at which projected density is to be output.
  integer                                             :: radiiCount  
  type   (radiusSpecifier), allocatable, dimension(:) :: radii
  type   (varying_string ), allocatable, dimension(:) :: outputDensityRadii
  logical                                             :: darkMatterScaleRadiusIsNeeded, diskIsNeeded        , &
       &                                                 spheroidIsNeeded             , virialRadiusIsNeeded

contains

  subroutine Galacticus_Output_Tree_Density_Initialize()
    !% Initializes the module by determining whether or not projected density data should be output.
    use Input_Parameters
    use Galacticus_Error
    implicit none

    if (.not.outputDensityDataInitialized) then
       !$omp critical(Galacticus_Output_Tree_Density_Initialize)
       if (.not.outputDensityDataInitialized) then
          !# <inputParameter>
          !#   <name>outputDensityData</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>.false.</defaultValue>
          !#   <description>Specifies whether or not density data should be included in the output file.</description>
          !#   <group>output</group>
          !#   <source>globalParameters</source>
          !#   <type>boolean</type>
          !# </inputParameter>
          !# <inputParameter>
          !#   <name>outputDensityIncludeRadii</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>.false.</defaultValue>
          !#   <description>Specifies whether or not the radii at which density data are output should also be included in the output file.</description>
          !#   <group>output</group>
          !#   <source>globalParameters</source>
          !#   <type>boolean</type>
          !# </inputParameter>
          ! Read and parse parameter controlling radii of output.
          if (outputDensityData) then
             if (.not.globalParameters%isPresent('outputDensityRadii'))                                      &
                  & call Galacticus_Error_Report(                                                            &
                  &                              'outputDensityRadii must be specified for density output'// &
                  &                              {introspection:location}                                    &
                  &                             )
             radiiCount=globalParameters%count('outputDensityRadii')
             allocate(outputDensityRadii(radiiCount))
             !# <inputParameter>
             !#   <name>outputDensityRadii</name>
             !#   <cardinality>1..*</cardinality>
             !#   <description>Specifies the radii at which the density will be output.</description>
             !#   <group>output</group>
             !#   <source>globalParameters</source>
             !#   <type>string</type>
             !# </inputParameter>
             ! Parse the radii definitions.
             call Galactic_Structure_Radii_Definition_Decode(outputDensityRadii,radii,diskIsNeeded,spheroidIsNeeded,virialRadiusIsNeeded,darkMatterScaleRadiusIsNeeded)
             deallocate(outputDensityRadii)
          end if
          ! Flag that module is now initialized.
          outputDensityDataInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Density_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Density_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Density_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Density</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Density_Names(node,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of density properties to be written to the \glc\ output file.
    use Galacticus_Nodes, only : treeNode
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode)              , intent(inout) :: node
    double precision                        , intent(in   ) :: time
    integer                                 , intent(inout) :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout) :: doublePropertyComments , doublePropertyNames   , &
         &                                                     integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout) :: doublePropertyUnitsSI  , integerPropertyUnitsSI
    integer                                                 :: i
    !GCC$ attributes unused :: node, integerProperty, integerPropertyNames, integerPropertyComments, integerPropertyUnitsSI, time
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Density_Initialize

    ! Return property names if we are outputting density data.
    if (outputDensityData) then
       do i=1,radiiCount
          doubleProperty=doubleProperty+1
          doublePropertyNames   (doubleProperty)='densityProfile:'//char(radii(i)%name)
          doublePropertyComments(doubleProperty)='density at a given radius'
          doublePropertyUnitsSI (doubleProperty)=+massSolar     &
               &                                 /megaParsec**3
          if (outputDensityIncludeRadii) then
             doubleProperty=doubleProperty+1
             doublePropertyNames   (doubleProperty)='densityProfileRadius:'//char(radii(i)%name)
             doublePropertyComments(doubleProperty)='radius at which density is output'
             doublePropertyUnitsSI (doubleProperty)=+megaParsec
          end if
       end do
    end if
    return
  end subroutine Galacticus_Output_Tree_Density_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Density_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Density</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Density_Property_Count(node,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of density properties to be written to the \glc\ output file.
    use Galacticus_Nodes, only : treeNode
    implicit none
    type            (treeNode), intent(inout) :: node
    double precision          , intent(in   ) :: time
    integer                   , intent(inout) :: doublePropertyCount, integerPropertyCount
    !GCC$ attributes unused :: node, time, integerPropertyCount
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Density_Initialize

    ! Increment property count if we are outputting density data.
    if (outputDensityData) then
       doublePropertyCount=doublePropertyCount+radiiCount
       if (outputDensityIncludeRadii)doublePropertyCount=doublePropertyCount+radiiCount
    end if
    return
  end subroutine Galacticus_Output_Tree_Density_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Density</unitName>
  !#  <sortName>Galacticus_Output_Tree_Density</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Density(node,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time,instance)
    !% Store density properties in the \glc\ output file buffers.
    use Galacticus_Nodes                  , only : treeNode, nodeComponentDisk, nodeComponentSpheroid, nodeComponentDarkMatterProfile
    use Kind_Numbers
    use Galactic_Structure_Densities
    use Dark_Matter_Halo_Scales
    use Galactic_Structure_Options
    use Galactic_Structure_Enclosed_Masses
    use Multi_Counters
    implicit none
    double precision                                , intent(in   ) :: time
    type            (treeNode                      ), intent(inout) :: node
    integer                                         , intent(inout) :: doubleBufferCount         , doubleProperty, integerBufferCount, &
         &                                                             integerProperty
    integer         (kind=kind_int8                ), intent(inout) :: integerBuffer        (:,:)
    double precision                                , intent(inout) :: doubleBuffer         (:,:)
    type            (multiCounter                  ), intent(inout) :: instance
    class           (nodeComponentDisk             ), pointer       :: disk
    class           (nodeComponentSpheroid         ), pointer       :: spheroid
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    class           (darkMatterHaloScaleClass      ), pointer       :: darkMatterHaloScale_
    integer                                                         :: i
    double precision                                                :: radius                    , radiusVirial
    !GCC$ attributes unused :: node, time, integerProperty, integerBufferCount, integerBuffer, instance
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Density_Initialize()
    ! Store property data if we are outputting density data.
    if (outputDensityData) then
       ! Compute required quantities.
       darkMatterHaloScale_ => darkMatterHaloScale()
       radiusVirial         = 0.0d0
       if (         virialRadiusIsNeeded) radiusVirial      =  darkMatterHaloScale_%virialRadius(node                    )
       if (                 diskIsNeeded) disk              =>                                   node%disk             ()
       if (             spheroidIsNeeded) spheroid          =>                                   node%spheroid         ()
       if (darkMatterScaleRadiusIsNeeded) darkMatterProfile =>                                   node%darkMatterProfile()
       do i=1,radiiCount
          ! Find the projected radius.
          radius=radii(i)%value
          select case (radii(i)%type)
          case   (radiusTypeRadius                )
             ! Nothing to do.
          case   (radiusTypeVirialRadius          )
             radius=+radius*radiusVirial
          case   (radiusTypeDarkMatterScaleRadius )
             radius=+radius*darkMatterProfile%         scale()
          case   (radiusTypeDiskRadius            )
             radius=+radius*disk             %        radius()
          case   (radiusTypeSpheroidRadius        )
             radius=+radius*spheroid         %        radius()
          case   (radiusTypeDiskHalfMassRadius    )
             radius=+radius*disk             %halfMassRadius()
          case   (radiusTypeSpheroidHalfMassRadius)
             radius=+radius*spheroid         %halfMassRadius()
          case   (radiusTypeGalacticMassFraction ,  &
               &  radiusTypeGalacticLightFraction )
             radius=+radius                                   &
                  & *Galactic_Structure_Radius_Enclosing_Mass &
                  &  (                                        &
                  &   node                                 ,  &
                  &   fractionalMass=radii(i)%fraction     ,  &
                  &   massType      =massTypeGalactic      ,  &
                  &   componentType =componentTypeAll      ,  &
                  &   weightBy      =radii(i)%weightBy     ,  &
                  &   weightIndex   =radii(i)%weightByIndex   &
                  &  )
          end select
          ! Store the density.
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=Galactic_Structure_Density(                                  &
              &                                                                     node                            , &
              &                                                                     [                                 &
              &                                                                      radius                         , &
              &                                                                      0.0d0                          , &
              &                                                                      0.0d0                            &
              &                                                                     ]                               , &
              &                                                                     componentType=radii(i)%component, &
              &                                                                     massType     =radii(i)%mass       &
              &                                                                    )
          if (outputDensityIncludeRadii) then
             doubleProperty=doubleProperty+1
             doubleBuffer(doubleBufferCount,doubleProperty)=radius
          end if
       end do
    end if
    return    
  end subroutine Galacticus_Output_Tree_Density

end module Galacticus_Output_Trees_Density
