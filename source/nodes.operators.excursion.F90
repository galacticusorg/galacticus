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
  Implements a node operator class that collects and stores the (infimum of the) excursion corresponding to the mass accretion history for each node.
  !!}
  
  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Linear_Growth             , only : linearGrowthClass 
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  
  !![
  <nodeOperator name="nodeOperatorExcursion">
   <description>A node operator class that collects and stores the (infimum of the) excursion corresponding to the mass accretion history for each node.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorExcursion
     !!{
     A node operator class that collects and stores the mass accretion history of each node.
     !!}
     private
     class           (cosmologyFunctionsClass       ), pointer :: cosmologyFunctions_       => null()
     class           (criticalOverdensityClass      ), pointer :: criticalOverdensity_      => null()
     class           (cosmologicalMassVarianceClass ), pointer :: cosmologicalMassVariance_ => null()
     class           (linearGrowthClass             ), pointer :: linearGrowth_             => null()
     integer                                                   :: excursionOverdensityID             , excursionVarianceID
     double precision                                          :: timePresent
   contains
     final     ::                   excursionDestructor
     procedure :: nodeInitialize => excursionNodeInitialize
  end type nodeOperatorExcursion
  
  interface nodeOperatorExcursion
     !!{
     Constructors for the \refClass{nodeOperatorExcursion} node operator class.
     !!}
     module procedure excursionConstructorParameters
     module procedure excursionConstructorInternal
  end interface nodeOperatorExcursion
  
contains

  function excursionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorExcursion} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorExcursion        )                :: self
    type (inputParameters              ), intent(inout) :: parameters
    class(cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class(criticalOverdensityClass     ), pointer       :: criticalOverdensity_
    class(cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    class(linearGrowthClass            ), pointer       :: linearGrowth_
    
    !![
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=nodeOperatorExcursion(cosmologyFunctions_,cosmologicalMassVariance_,criticalOverdensity_,linearGrowth_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="linearGrowth_"            />
    !!]
    return
  end function excursionConstructorParameters

  function excursionConstructorInternal(cosmologyFunctions_,cosmologicalMassVariance_,criticalOverdensity_,linearGrowth_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorExcursion} node operator class.
    !!}
    use :: Galacticus_Nodes, only : defaultBasicComponent
    implicit none
    type (nodeOperatorExcursion        )                        :: self
    class(criticalOverdensityClass     ), intent(in   ), target :: criticalOverdensity_
    class(cosmologyFunctionsClass      ), intent(in   ), target :: cosmologyFunctions_
    class(cosmologicalMassVarianceClass), intent(in   ), target :: cosmologicalMassVariance_
    class(linearGrowthClass            ), intent(in   ), target :: linearGrowth_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_, *linearGrowth_"/>
    !!]
    
    self%timePresent=self%cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0)
    !![
    <addMetaProperty component="basic" name="excursionTime" id="self%excursionOverdensityID" rank="1" isCreator="yes"/>
    <addMetaProperty component="basic" name="excursionMass" id="self%excursionVarianceID"    rank="1" isCreator="yes"/>
    !!]
    return
  end function excursionConstructorInternal

  subroutine excursionDestructor(self)
    !!{
    Destructor for the critical overdensity excursion set barrier class.
    !!}
    implicit none
    type(nodeOperatorExcursion), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%linearGrowth_"            />
    !!]
    return
  end subroutine excursionDestructor

  subroutine excursionNodeInitialize(self,node)
    !!{
    Record the mass accretion history of the node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (nodeOperatorExcursion), intent(inout), target      :: self
    type            (treeNode             ), intent(inout), target      :: node
    double precision                       , dimension(:) , allocatable :: overdensities, variances
    type            (treeNode             )               , pointer     :: nodeWork
    class           (nodeComponentBasic   )               , pointer     :: basic
    integer                                                             :: countNodes
    
    if (associated(node%firstChild)) return
    nodeWork   => node
    countNodes =  1
    do while (nodeWork%isPrimaryProgenitor())
       countNodes =  countNodes       +1
       nodeWork   => nodeWork  %parent
    end do
    allocate(overdensities(countNodes))
    allocate(variances    (countNodes))
    nodeWork                  => node
    countNodes                =  1
    basic                     => nodeWork%basic()
    overdensities(countNodes) =  +self      %criticalOverdensity_     %value       (time=basic%time       (),mass=basic%mass())    &
            &                    /self      %linearGrowth_            %value       (time=basic%time       ()                  )
    variances    (countNodes) =  +self      %cosmologicalMassVariance_%rootVariance(time=self %timePresent  ,mass=basic%mass())**2
    do while (nodeWork%isPrimaryProgenitor())
       countNodes                =   countNodes+1
       nodeWork                  =>  nodeWork                            %parent
       basic                     =>  nodeWork                            %basic       (                                          )
       overdensities(countNodes) =  +self      %criticalOverdensity_     %value       (time=basic%time       (),mass=basic%mass())    &
            &                       /self      %linearGrowth_            %value       (time=basic%time       ()                  )
       variances    (countNodes) =  +self      %cosmologicalMassVariance_%rootVariance(time=self %timePresent  ,mass=basic%mass())**2
    end do
    basic => node%basic()
    call basic%floatRank1MetaPropertySet(self%excursionOverdensityID,overdensities)
    call basic%floatRank1MetaPropertySet(self%excursionVarianceID   ,variances    )
    return
  end subroutine excursionNodeInitialize
  
