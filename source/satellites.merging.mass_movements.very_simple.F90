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

  !% Implements a merger mass movements class which uses a simple calculation.
  
  !# <mergerMassMovements name="mergerMassMovementsVerySimple">
  !#  <description>A merger mass movements class which uses a simple calculation.</description>
  !# </mergerMassMovements>
  type, extends(mergerMassMovementsClass) :: mergerMassMovementsVerySimple
     !% A merger mass movements class which uses a simple calculation.
     private
     double precision :: massRatioMajorMerger
   contains
     procedure :: get => verySimpleGet
  end type mergerMassMovementsVerySimple

  interface mergerMassMovementsVerySimple
     !% Constructors for the {\normalfont \ttfamily verySimple} merger mass movements class.
     module procedure verySimpleConstructorParameters
     module procedure verySimpleConstructorInternal
  end interface mergerMassMovementsVerySimple

contains

  function verySimpleConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily verySimple} merger mass movements class which takes a parameter list as input.
    use Input_Parameters
    implicit none
    type            (mergerMassMovementsVerySimple)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    double precision                                               :: massRatioMajorMerger

    !# <inputParameter>
    !#   <name>massRatioMajorMerger</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.25d0</defaultValue>
    !#   <description>The mass ratio above which mergers are considered to be ``major''.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    self=mergerMassMovementsVerySimple(massRatioMajorMerger)
    !# <inputParametersValidate source="parameters"/>
    return
  end function verySimpleConstructorParameters

  function verySimpleConstructorInternal(massRatioMajorMerger) result(self)
    !% Internal constructor for the {\normalfont \ttfamily verySimple} merger mass movements.
    implicit none
    type            (mergerMassMovementsVerySimple)                :: self
    double precision                               , intent(in   ) :: massRatioMajorMerger
    !# <constructorAssign variables="massRatioMajorMerger"/>
    
    return
  end function verySimpleConstructorInternal

  subroutine verySimpleGet(self,node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    !% Determine where stars and gas move as the result of a merger event using a very simple algorithm.
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    class           (mergerMassMovementsVerySimple), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    integer                                        , intent(  out) :: destinationGasSatellite, destinationGasHost       , &
         &                                                            destinationStarsHost   , destinationStarsSatellite
    logical                                        , intent(  out) :: mergerIsMajor
    type            (treeNode                     ), pointer       :: nodeHost
    double precision                                               :: massHost               , massSatellite

    if      (self%massRatioMajorMerger <= 0.0d0) then
       mergerIsMajor=.true.
    else if (self%massRatioMajorMerger >  1.0d0) then
       mergerIsMajor=.false.
    else
       nodeHost      => node%mergesWith()
       massSatellite =  Galactic_Structure_Enclosed_Mass(node    ,massType=massTypeGalactic)
       massHost      =  Galactic_Structure_Enclosed_Mass(nodeHost,massType=massTypeGalactic)    
       mergerIsMajor =  massSatellite >= self%massRatioMajorMerger*massHost
    end if
    destinationGasSatellite  =destinationMergerDisk
    destinationStarsSatellite=destinationMergerDisk
    destinationGasHost       =destinationMergerUnmoved
    destinationStarsHost     =destinationMergerUnmoved
    return
  end subroutine verySimpleGet
