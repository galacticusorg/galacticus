!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which implements a property extractor class for the mass enclosed by a set of radii.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMassProfile">
   <description>
    A node property extractor class which extracts the mass enclosed within specified radii. A list of radii, $r$ (in Mpc),
    must be specified via the {\normalfont \ttfamily [radii]} parameter. For each specified radius, the total (dark + baryonic)
    mass will be extracted as {\normalfont \ttfamily massProfile}$r$.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorArray) :: nodePropertyExtractorMassProfile
     !!{
     A property extractor class for the mass enclosed within a specified list of radii.
     !!}
     private
     integer         (c_size_t)                            :: radiiCount
     double precision          , allocatable, dimension(:) :: radii
   contains
     procedure :: columnDescriptions => massProfileColumnDescriptions
     procedure :: size               => massProfileSize
     procedure :: elementCount       => massProfileElementCount
     procedure :: extract            => massProfileExtract
     procedure :: names              => massProfileNames
     procedure :: descriptions       => massProfileDescriptions
     procedure :: unitsInSI          => massProfileUnitsInSI
     procedure :: type               => massProfileType
  end type nodePropertyExtractorMassProfile

  interface nodePropertyExtractorMassProfile
     !!{
     Constructors for the ``massProfile'' output analysis class.
     !!}
     module procedure massProfileConstructorParameters
     module procedure massProfileConstructorInternal
  end interface nodePropertyExtractorMassProfile

contains

  function massProfileConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily massProfile} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nodePropertyExtractorMassProfile)                              :: self
    type            (inputParameters                 ), intent(inout)               :: parameters
    double precision                                  , allocatable  , dimension(:) :: radii

    allocate(radii(parameters%count('radii')))
    !![
    <inputParameter>
      <name>radii</name>
      <description>A list of radii at which to output the mass profile.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=nodePropertyExtractorMassProfile(radii)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massProfileConstructorParameters

  function massProfileConstructorInternal(radii) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily massProfile} property extractor class.
    !!}
    use :: Galactic_Structure_Options, only : massTypeAll, massTypeDark
    implicit none
    type            (nodePropertyExtractorMassProfile)                              :: self
    double precision                                  , intent(in   ), dimension(:) :: radii
    !![
    <constructorAssign variables="radii"/>
    !!]

    self%radiiCount=size(radii)
    return
  end function massProfileConstructorInternal

  integer function massProfileElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily massProfile} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorMassProfile), intent(inout) :: self
    double precision                                  , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    massProfileElementCount=1
    return
  end function massProfileElementCount

  function massProfileSize(self,time)
    !!{
    Return the number of array alements in the {\normalfont \ttfamily massProfile} property extractors.
    !!}
    implicit none
    integer         (c_size_t                        )                :: massProfileSize
    class           (nodePropertyExtractorMassProfile), intent(inout) :: self
    double precision                                  , intent(in   ) :: time
    !$GLC attributes unused :: time

    massProfileSize=self%radiiCount
    return
  end function massProfileSize

  function massProfileExtract(self,node,time,instance)
    !!{
    Implement a last isolated redshift output analysis.
    !!}
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Enclosed_Mass
    use :: Galactic_Structure_Options        , only : componentTypeAll                , massTypeAll
    implicit none
    double precision                                  , dimension(:,:), allocatable :: massProfileExtract
    class           (nodePropertyExtractorMassProfile), intent(inout) , target      :: self
    type            (treeNode                        ), intent(inout) , target      :: node
    double precision                                  , intent(in   )               :: time
    type            (multiCounter                    ), intent(inout) , optional    :: instance
    integer         (c_size_t                        )                              :: i
    !$GLC attributes unused :: time, instance

    allocate(massProfileExtract(self%radiiCount,1_c_size_t))
    do i=1,self%radiiCount
       massProfileExtract(i,1)=Galactic_Structure_Enclosed_Mass(node,self%radii(i),componentType=componentTypeAll,massType=massTypeAll)
    end do
    return
  end function massProfileExtract

  function massProfileNames(self,time)
    !!{
    Return the names of the {\normalfont \ttfamily massProfile} properties.
    !!}
    implicit none
    type            (varying_string                  ), dimension(:) , allocatable :: massProfileNames
    class           (nodePropertyExtractorMassProfile), intent(inout)              :: self
    double precision                                  , intent(in   )              :: time
    !$GLC attributes unused :: time

    allocate(massProfileNames(1))
    massProfileNames(1)='massProfile'
    return
  end function massProfileNames

  function massProfileDescriptions(self,time)
    !!{
    Return descriptions of the {\normalfont \ttfamily massProfile} property.
    !!}
    implicit none
    type            (varying_string                  ), dimension(:) , allocatable :: massProfileDescriptions
    class           (nodePropertyExtractorMassProfile), intent(inout)              :: self
    double precision                                  , intent(in   )              :: time
    !$GLC attributes unused :: time

    allocate(massProfileDescriptions(1))
    massProfileDescriptions(1)='Mass enclosed within each radius'
    return
  end function massProfileDescriptions

  function massProfileColumnDescriptions(self,time)
    !!{
    Return column descriptions of the {\normalfont \ttfamily massProfile} property.
    !!}
    implicit none
    type            (varying_string                  ), dimension(:) , allocatable :: massProfileColumnDescriptions
    class           (nodePropertyExtractorMassProfile), intent(inout)              :: self
    double precision                                  , intent(in   )              :: time
    integer         (c_size_t                        )                             :: i
    character       (len=22                          )                             :: name
    !$GLC attributes unused :: time
    
    allocate(massProfileColumnDescriptions(self%radiiCount))
    do i=1,self%radiiCount
       write (name,'(a,e9.3,a)') 'r = ',self%radii(i),' Mpc'       
       massProfileColumnDescriptions(i)=trim(adjustl(name))
    end do
    return
  end function massProfileColumnDescriptions
  
  function massProfileUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily massProfile} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    double precision                                  , allocatable  , dimension(:) :: massProfileUnitsInSI
    class           (nodePropertyExtractorMassProfile), intent(inout)               :: self
    double precision                                  , intent(in   )               :: time
    !$GLC attributes unused :: time

    allocate(massProfileUnitsInSI(1))
    massProfileUnitsInSI(1)=massSolar
    return
  end function massProfileUnitsInSI

  integer function massProfileType(self)
    !!{
    Return the type of the {\normalfont \ttfamily massProfile} properties.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorMassProfile), intent(inout) :: self
    !$GLC attributes unused :: self

    massProfileType=outputAnalysisPropertyTypeLinear
    return
  end function massProfileType
