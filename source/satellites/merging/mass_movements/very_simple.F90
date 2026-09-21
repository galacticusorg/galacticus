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

  !!{RST
  Implements a merger mass movements class which uses a simple calculation.
  !!}

  !![
  <mergerMassMovements name="mergerMassMovementsVerySimple" docformat="rst">
   <description>
   A merger mass movements class which assumes that the satellite material is always added to the disk of the host, while the host mass is not moved.
   </description>
  </mergerMassMovements>
  !!]
  type, extends(mergerMassMovementsMemoized) :: mergerMassMovementsVerySimple
     !!{RST
     A merger mass movements class which uses a simple calculation.
     !!}
     private
     double precision :: massRatioMajorMerger
   contains
     final     ::              verySimpleDestructor
     procedure :: calculate => verySimpleCalculate
  end type mergerMassMovementsVerySimple

  interface mergerMassMovementsVerySimple
     !!{RST
     Constructors for the :galacticus-class:`mergerMassMovementsVerySimple` merger mass movements class.
     !!}
     module procedure verySimpleConstructorParameters
     module procedure verySimpleConstructorInternal
  end interface mergerMassMovementsVerySimple

contains

  function verySimpleConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`mergerMassMovementsVerySimple` merger mass movements class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerMassMovementsVerySimple)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    double precision                                               :: massRatioMajorMerger

    !![
    <inputParameter docformat="rst">
      <name>massRatioMajorMerger</name>
      <defaultValue>0.25d0</defaultValue>
      <description>
      The mass ratio above which mergers are considered to be "major".
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=mergerMassMovementsVerySimple(massRatioMajorMerger)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function verySimpleConstructorParameters

  function verySimpleConstructorInternal(massRatioMajorMerger) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`mergerMassMovementsVerySimple` merger mass movements.
    !!}
    implicit none
    type            (mergerMassMovementsVerySimple)                :: self
    double precision                               , intent(in   ) :: massRatioMajorMerger
    !![
    <constructorAssign variables="massRatioMajorMerger"/>
    !!]

    return
  end function verySimpleConstructorInternal

  subroutine verySimpleDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`mergerMassMovementsVerySimple` merger mass movements class.
    !!}
    implicit none
    type(mergerMassMovementsVerySimple), intent(inout) :: self

    call self%detachHooks()
    return
  end subroutine verySimpleDestructor

  subroutine verySimpleCalculate(self,node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    !!{RST
    Determine where stars and gas move as the result of a merger event using a very simple algorithm.
    !!}
    use :: Galactic_Structure_Options, only : massTypeGalactic
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class           (mergerMassMovementsVerySimple   ), intent(inout)         :: self
    type            (treeNode                        ), intent(inout), target :: node
    type            (enumerationDestinationMergerType), intent(  out)         :: destinationGasSatellite  , destinationGasHost       , &
         &                                                                       destinationStarsHost     , destinationStarsSatellite
    logical                                           , intent(  out)         :: mergerIsMajor
    type            (treeNode                        ), pointer               :: nodeHost
    class           (massDistributionClass           ), pointer               :: massDistributionSatellite, massDistributionHost
    double precision                                                          :: massHost                 , massSatellite
    
    if      (self%massRatioMajorMerger <= 0.0d0) then
       mergerIsMajor=.true.
    else if (self%massRatioMajorMerger >  1.0d0) then
       mergerIsMajor=.false.
    else
       nodeHost                  => node                     %mergesWith      (                         )
       massDistributionHost      => nodeHost                 %massDistribution(massType=massTypeGalactic)
       massDistributionSatellite => node                     %massDistribution(massType=massTypeGalactic)
       massSatellite             =  massDistributionSatellite%massTotal       (                         )
       massHost                  =  massDistributionHost     %massTotal       (                         )
       mergerIsMajor             =  massSatellite >= self%massRatioMajorMerger*massHost
       !![
	  <objectDestructor name="massDistributionHost"     />
	  <objectDestructor name="massDistributionSatellite"/>
	  !!]
    end if
    destinationGasSatellite     =  destinationMergerDisk
    destinationStarsSatellite   =  destinationMergerDisk
    destinationGasHost          =  destinationMergerUnmoved
    destinationStarsHost        =  destinationMergerUnmoved
    return
  end subroutine verySimpleCalculate
