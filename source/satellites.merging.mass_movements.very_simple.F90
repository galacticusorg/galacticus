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
  Implements a merger mass movements class which uses a simple calculation.
  !!}

  use :: Kind_Numbers, only : kind_int8

  !![
  <mergerMassMovements name="mergerMassMovementsVerySimple">
   <description>
    A merger mass movements class which assumes that the satellite material is always added to the disk of the host, while the
    host mass is not moved.
   </description>
  </mergerMassMovements>
  !!]
  type, extends(mergerMassMovementsClass) :: mergerMassMovementsVerySimple
     !!{
     A merger mass movements class which uses a simple calculation.
     !!}
     private
     double precision                          :: massRatioMajorMerger
     integer         (kind=kind_int8)          :: lastUniqueID
     logical                                   :: mergerIsMajor       , movementsCalculated
   contains
     final     ::             verySimpleDestructor
     procedure :: autoHook => verySimpleAutoHook
     procedure :: get      => verySimpleGet
  end type mergerMassMovementsVerySimple

  interface mergerMassMovementsVerySimple
     !!{
     Constructors for the \refClass{mergerMassMovementsVerySimple} merger mass movements class.
     !!}
     module procedure verySimpleConstructorParameters
     module procedure verySimpleConstructorInternal
  end interface mergerMassMovementsVerySimple

contains

  function verySimpleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerMassMovementsVerySimple} merger mass movements class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerMassMovementsVerySimple)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    double precision                                               :: massRatioMajorMerger

    !![
    <inputParameter>
      <name>massRatioMajorMerger</name>
      <defaultValue>0.25d0</defaultValue>
      <description>The mass ratio above which mergers are considered to be ``major''.</description>
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
    !!{
    Internal constructor for the \refClass{mergerMassMovementsVerySimple} merger mass movements.
    !!}
    implicit none
    type            (mergerMassMovementsVerySimple)                :: self
    double precision                               , intent(in   ) :: massRatioMajorMerger
    !![
    <constructorAssign variables="massRatioMajorMerger"/>
    !!]

    self%lastUniqueID       =-huge(0_kind_int8)
    self%mergerIsMajor      =.false.
    self%movementsCalculated=.false.
    return
  end function verySimpleConstructorInternal

  subroutine verySimpleAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels, satelliteMergerEvent
    implicit none
    class(mergerMassMovementsVerySimple), intent(inout) :: self

    call calculationResetEvent%attach(self,verySimpleCalculationReset,openMPThreadBindingAllLevels,label='remnantStructure:massMovementsVerySimple')
    call satelliteMergerEvent %attach(self,verySimpleGetHook         ,openMPThreadBindingAllLevels,label='remnantStructure:massMovementsVerySimple')
    return
  end subroutine verySimpleAutoHook

  subroutine verySimpleDestructor(self)
    !!{
    Destructor for the \refClass{mergerMassMovementsVerySimple} dark matter halo profile class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, satelliteMergerEvent
    implicit none
    type(mergerMassMovementsVerySimple), intent(inout) :: self

    if (calculationResetEvent%isAttached(self,verySimpleCalculationReset)) call calculationResetEvent%detach(self,verySimpleCalculationReset)
    if (satelliteMergerEvent %isAttached(self,verySimpleGetHook         )) call satelliteMergerEvent %detach(self,verySimpleGetHook         )
    return
  end subroutine verySimpleDestructor

  subroutine verySimpleCalculationReset(self,node,uniqueID)
    !!{
    Reset the dark matter profile calculation.
    !!}
    use :: Error       , only : Error_Report
    use :: Kind_Numbers, only : kind_int8
    implicit none
    class  (*        ), intent(inout) :: self
    type   (treeNode ), intent(inout) :: node
    integer(kind_int8), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    select type (self)
    class is (mergerMassMovementsVerySimple)
       self%movementsCalculated=.false.
       self%lastUniqueID       =uniqueID
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine verySimpleCalculationReset

  subroutine verySimpleGetHook(self,node)
    !!{
    Hookable wrapper around the get function.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (*                               ), intent(inout)         :: self
    type   (treeNode                        ), intent(inout), target :: node
    type   (enumerationDestinationMergerType)                        :: destinationGasSatellite, destinationGasHost       , &
         &                                                              destinationStarsHost   , destinationStarsSatellite
    logical                                                          :: mergerIsMajor

    select type (self)
    type is (mergerMassMovementsVerySimple)
       call self%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine verySimpleGetHook

  subroutine verySimpleGet(self,node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    !!{
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
    
    ! The calculation of how mass moves as a result of the merger is computed when first needed and then stored. This ensures that
    ! the results are determined by the properties of the merge target prior to any modification that will occur as node
    ! components are modified in response to the merger.
    if (node%uniqueID() /= self%lastUniqueID) call verySimpleCalculationReset(self,node,node%uniqueID())
    if (.not.self%movementsCalculated) then
       self%movementsCalculated=.true.
       if      (self%massRatioMajorMerger <= 0.0d0) then
          self%mergerIsMajor=.true.
       else if (self%massRatioMajorMerger >  1.0d0) then
          self%mergerIsMajor=.false.
       else
          nodeHost                  => node                     %mergesWith      (                         )
          massDistributionHost      => nodeHost                 %massDistribution(massType=massTypeGalactic)
          massDistributionSatellite => node                     %massDistribution(massType=massTypeGalactic)
          massSatellite             =  massDistributionSatellite%massTotal       (                         )
          massHost                  =  massDistributionHost     %massTotal       (                         )
          self%mergerIsMajor        =  massSatellite >= self%massRatioMajorMerger*massHost
          !![
	  <objectDestructor name="massDistributionHost"     />
	  <objectDestructor name="massDistributionSatellite"/>
	  !!]
       end if
    end if
    mergerIsMajor            =self%mergerIsMajor
    destinationGasSatellite  =     destinationMergerDisk
    destinationStarsSatellite=     destinationMergerDisk
    destinationGasHost       =     destinationMergerUnmoved
    destinationStarsHost     =     destinationMergerUnmoved
    return
  end subroutine verySimpleGet
