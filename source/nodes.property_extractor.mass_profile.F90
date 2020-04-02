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

!% Contains a module which implements a property extractor class for the mass enclosed by a set of radii.

  !# <nodePropertyExtractor name="nodePropertyExtractorMassProfile">
  !#  <description>A property extractor class for the mass enclosed by a set of radii.</description>
  !# </nodePropertyExtractor>
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorMassProfile
     !% A property extractor class for the mass and radii of spheres are specified density contrast..
     private
     integer                                     :: elementCount_
     double precision, allocatable, dimension(:) :: radii
   contains
     procedure :: elementCount => massProfileElementCount
     procedure :: extract      => massProfileExtract
     procedure :: names        => massProfileNames
     procedure :: descriptions => massProfileDescriptions
     procedure :: unitsInSI    => massProfileUnitsInSI
     procedure :: type         => massProfileType
  end type nodePropertyExtractorMassProfile

  interface nodePropertyExtractorMassProfile
     !% Constructors for the ``massProfile'' output analysis class.
     module procedure massProfileConstructorParameters
     module procedure massProfileConstructorInternal
  end interface nodePropertyExtractorMassProfile

contains

  function massProfileConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily massProfile} property extractor class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nodePropertyExtractorMassProfile)                              :: self
    type            (inputParameters                 ), intent(inout)               :: parameters
    double precision                                  , allocatable  , dimension(:) :: radii

    allocate(radii(parameters%count('radii')))
    !# <inputParameter>
    !#   <name>radii</name>
    !#   <cardinality>1..*</cardinality>
    !#   <description>A list of radii at which to output the mass profile.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    self=nodePropertyExtractorMassProfile(radii)
    !# <inputParametersValidate source="parameters"/>
    return
  end function massProfileConstructorParameters

  function massProfileConstructorInternal(radii) result(self)
    !% Internal constructor for the {\normalfont \ttfamily massProfile} property extractor class.
    use :: Galactic_Structure_Options, only : massTypeAll, massTypeDark
    implicit none
    type            (nodePropertyExtractorMassProfile)                              :: self
    double precision                                  , intent(in   ), dimension(:) :: radii
    !# <constructorAssign variables="radii"/>

    self%elementCount_=size(radii)
    return
  end function massProfileConstructorInternal

  integer function massProfileElementCount(self,time)
    !% Return the number of elements in the {\normalfont \ttfamily massProfile} property extractors.
    implicit none
    class           (nodePropertyExtractorMassProfile), intent(inout) :: self
    double precision                                  , intent(in   ) :: time
    !$GLC attributes unused :: time

    massProfileElementCount=self%elementCount_
    return
  end function massProfileElementCount

  function massProfileExtract(self,node,time,instance)
    !% Implement a last isolated redshift output analysis.
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Enclosed_Mass
    use :: Galactic_Structure_Options        , only : componentTypeAll                , massTypeAll
    implicit none
    double precision                                  , dimension(:) , allocatable :: massProfileExtract
    class           (nodePropertyExtractorMassProfile), intent(inout), target      :: self
    type            (treeNode                        ), intent(inout), target      :: node
    double precision                                  , intent(in   )              :: time
    type            (multiCounter                    ), intent(inout), optional    :: instance
    integer                                                                        :: i
    !$GLC attributes unused :: time, instance

    allocate(massProfileExtract(self%elementCount_))
    do i=1,self%elementCount_
       massProfileExtract(i)=Galactic_Structure_Enclosed_Mass(node,self%radii(i),componentType=componentTypeAll,massType=massTypeAll)
    end do
    return
  end function massProfileExtract

  function massProfileNames(self,time)
    !% Return the names of the {\normalfont \ttfamily massProfile} properties.
    implicit none
    type            (varying_string                  ), dimension(:) , allocatable :: massProfileNames
    class           (nodePropertyExtractorMassProfile), intent(inout)              :: self
    double precision                                  , intent(in   )              :: time
    integer                                                                        :: i
    character       (len=22                          )                             :: name
    !$GLC attributes unused :: time

    allocate(massProfileNames(self%elementCount_))
    do i=1,self%elementCount_
       write (name,'(a,e9.3)') 'massProfile',self%radii(i)
       massProfileNames(i)=trim(adjustl(name))
    end do
    return
  end function massProfileNames

  function massProfileDescriptions(self,time)
    !% Return descriptions of the {\normalfont \ttfamily massProfile} property.
    implicit none
    type            (varying_string                  ), dimension(:) , allocatable :: massProfileDescriptions
    class           (nodePropertyExtractorMassProfile), intent(inout)              :: self
    double precision                                  , intent(in   )              :: time
    integer                                                                        :: i
    character       (len=64                          )                             :: description
    !$GLC attributes unused :: time

    allocate(massProfileDescriptions(self%elementCount_))
    do i=1,self%elementCount_
       write (description,'(a,e9.3,a)') 'Mass enclosed within a radius of ',self%radii(i),' Mpc [M☉].'
       massProfileDescriptions(i)=trim(adjustl(description))
    end do
    return
  end function massProfileDescriptions

  function massProfileUnitsInSI(self,time)
    !% Return the units of the {\normalfont \ttfamily massProfile} properties in the SI system.
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    double precision                                  , allocatable  , dimension(:) :: massProfileUnitsInSI
    class           (nodePropertyExtractorMassProfile), intent(inout)               :: self
    double precision                                  , intent(in   )               :: time
    integer                                                                         :: i
    !$GLC attributes unused :: time

    allocate(massProfileUnitsInSI(self%elementCount_))
    do i=1,self%elementCount_
       massProfileUnitsInSI(i)=massSolar
    end do
    return
  end function massProfileUnitsInSI

  integer function massProfileType(self)
    !% Return the type of the {\normalfont \ttfamily massProfile} properties.
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorMassProfile), intent(inout) :: self
    !$GLC attributes unused :: self

    massProfileType=outputAnalysisPropertyTypeLinear
    return
  end function massProfileType
