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

  !+    Contributions to this file made by: Matías Liempi

  !!{
  Implements a node operator class that computes the growth of the \gls{nsc} from spheroid gas.
  !!}

  use :: Nuclear_Star_Cluster_Growth_Rates, only : nuclearStarClusterGrowthRatesClass
  
  !![
  <nodeOperator name="nodeOperatorNuclearStarClusterGrowth">
    <description>
      A node operator class that handle the gas inflow rate ontothe nuclear star cluster.
    </description>
  </nodeOperator>
  !!]

  type, extends(nodeOperatorClass) :: nodeOperatorNuclearStarClusterGrowth
     !!{
     A node operator class that computes the growth of the \gls{nsc} from spheroid gas.
     !!}
     private
     class  (nuclearStarClusterGrowthRatesClass), pointer :: nuclearStarClusterGrowthRates_ => null()
   contains
     final     ::                          nuclearStarClusterGrowthDestructor
     procedure :: differentialEvolution => nuclearStarClusterGrowthDifferentialEvolution
  end type nodeOperatorNuclearStarClusterGrowth
  
  interface nodeOperatorNuclearStarClusterGrowth
     !!{
     Constructors for the {\normalfont \ttfamily gasMassRateNSC} node operator class.
     !!}
     module procedure nuclearStarClusterGrowthConstructorParameters
     module procedure nuclearStarClusterGrowthConstructorInternal
  end interface nodeOperatorNuclearStarClusterGrowth
  
contains

  function nuclearStarClusterGrowthConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily gasMassRateNSC} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorNuclearStarClusterGrowth)                :: self
    type (inputParameters                     ), intent(inout) :: parameters
    class(nuclearStarClusterGrowthRatesClass  ), pointer       :: nuclearStarClusterGrowthRates_

    !![
    <objectBuilder class="nuclearStarClusterGrowthRates" name="nuclearStarClusterGrowthRates_" source="parameters"/>
    !!]
    self=nodeOperatorNuclearStarClusterGrowth(nuclearStarClusterGrowthRates_)
    !![
    <inputParametersValidate source="parameters"      />
    <objectDestructor name="nuclearStarClusterGrowthRates_"/>
    !!]
    return
  end function nuclearStarClusterGrowthConstructorParameters
  
  function nuclearStarClusterGrowthConstructorInternal(nuclearStarClusterGrowthRates_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily NuclearStarClusterGrowth} node operator class.
    !!}
    implicit none
    type (nodeOperatorNuclearStarClusterGrowth)                        :: self
    class(nuclearStarClusterGrowthRatesClass  ), intent(in   ), target :: nuclearStarClusterGrowthRates_
    !![
    <constructorAssign variables="*nuclearStarClusterGrowthRates_"/>
    !!]

    return
  end function nuclearStarClusterGrowthConstructorInternal
  
  subroutine nuclearStarClusterGrowthDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily NuclearStarClusterGrowth} node operator class.
    !!}
    implicit none
    type(nodeOperatorNuclearStarClusterGrowth), intent(inout) :: self
    
    !![
    <objectDestructor name="self%nuclearStarClusterGrowthRates_"/>
    !!]
    return
  end subroutine nuclearStarClusterGrowthDestructor

  subroutine nuclearStarClusterGrowthDifferentialEvolution(self,node,interrupt,functioninterrupt,propertyType)
    !!{
    Compute the growth rate of \gls{nsc} due to accretion from the spheroid.
    !!}
    use :: Galacticus_Nodes    , only : interruptTask   , nodeComponentNSC, nodeComponentNSCStandard, nodeComponentSpheroid, &
       &                                propertyInactive, treeNode
    use :: Abundances_Structure, only : operator(*)
    implicit none
    class           (nodeOperatorNuclearStarClusterGrowth), intent(inout), target :: self
    type            (treeNode                            ), intent(inout), target :: node
    logical                                               , intent(inout)         :: interrupt
    procedure       (interruptTask                       ), intent(inout), pointer:: functioninterrupt
    integer                                               , intent(in   )         :: propertyType
    class           (nodeComponentNSC                    ),                pointer:: nuclearStarCluster
    class           (nodeComponentSpheroid               ),                pointer:: spheroid
    double precision                                                              :: rateMassAccretionGas

    ! Return immediately if inactive variables are requested.
    if (propertyInactive(propertyType)) return
    ! Find the rate of gas accretion.
    rateMassAccretionGas=self%nuclearStarClusterGrowthRates_%rate(node)
    ! Finish if there is no gas inflow onto the nuclear star cluster.
    if (rateMassAccretionGas <= 0.0d0) return
    ! Get the spheroid and nuclear star cluster component.
    spheroid           => node%spheroid()
    nuclearStarCluster => node%     NSC()
    ! Detect nuclear star cluster component type.
    select type (nuclearStarCluster)
    type is (nodeComponentNSC)
       ! Generic type - interrupt and create a nuclear star cluster.
       interrupt         =  .true.
       functionInterrupt => nuclearStarClusterCreate
       return
    class default
       ! A nuclear star cluster exists - continue processing.
       ! Remove gas from the spheroid component and add to the nuclear star cluster component.
       call spheroid          %        massGasRate(-rateMassAccretionGas)
       call spheroid          %angularMomentumRate(-rateMassAccretionGas*spheroid%angularMomentum()/(spheroid%massGas()+spheroid%massStellar()))
       call spheroid          %  abundancesGasRate(-rateMassAccretionGas*spheroid%abundancesGas  ()/ spheroid%massGas()                        )
       call nuclearStarCluster%        massGasRate(+rateMassAccretionGas)
       call nuclearStarCluster%angularMomentumRate(+rateMassAccretionGas*spheroid%angularMomentum()/(spheroid%massGas()+spheroid%massStellar()))
       call nuclearStarCluster%  abundancesGasRate(+rateMassAccretionGas*spheroid%abundancesGas  ()/ spheroid%massGas()                        )
    end select
    return
  end subroutine nuclearStarClusterGrowthDifferentialEvolution

  subroutine nuclearStarClusterCreate(node,timeEnd)
    !!{
    Creates the nuclear star cluster via interrupt.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC, treeNode
    implicit none
    type             (treeNode        ), intent(inout), target  :: node
    double precision                   , intent(in   ), optional:: timeEnd
    class            (nodeComponentNSC),                pointer :: nuclearStarCluster
    !$GLC attributes unused :: timeEnd

    nuclearStarCluster => node%NSC(autoCreate=.true.)
    return 
  end subroutine nuclearStarClusterCreate
  
