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
  Provides a module storing object references needed for galaxy merger trees.
  !!}
  
  module Node_Property_Extractor_Galaxy_Merger_Trees
    !!{
    A module storing object references needed for galaxy merger trees.
    !!}
    private
    
    ! NOTE: This is an unpleasant solution - essentially functioning as a communicator between the relevant `nodeOperator` and
    ! `nodePropertyExtractor` classes. It is necessitated by the fact that we don't have a good way for objects defined in the
    ! parameter file to communicate their behavior.

    type, public :: nodePropertyExtractorListWrapper
       class(*), pointer :: extractor_ => null()
    end type nodePropertyExtractorListWrapper
     
    ! Public-scope pointer to the extractor used in galaxy merger trees.
    integer                                                             , public :: nodePropertyExtractorGalaxyMergerTreeCount
    type   (nodePropertyExtractorListWrapper), allocatable, dimension(:), public :: nodePropertyExtractorGalaxyMergerTree_
    !$omp threadprivate(nodePropertyExtractorGalaxyMergerTreeCount,nodePropertyExtractorGalaxyMergerTree_)
  
  end module Node_Property_Extractor_Galaxy_Merger_Trees
  
