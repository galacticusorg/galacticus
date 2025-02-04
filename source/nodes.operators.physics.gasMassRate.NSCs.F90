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
  Implements a node operator class that handle the gas mass rate in the nuclear star cluster.
  !!}
  use :: Star_Formation_Rates_Spheroids , only : starFormationRateSpheroidsClass

  !![
  <nodeOperator name="nodeOperatorgasMassRateNSC">
   <description>A node operator class that handle the gass mass rate in the nuclear star cluster.</description>
  </nodeOperator>
  !!]

  type, extends(nodeOperatorClass) :: nodeOperatorgasMassRateNSC
     !!{
     A node operator class that handle the gas mass rate in the nuclear star cluster.
     !!}
     private
     class(starFormationRateSpheroidsClass), pointer :: starFormationRateSpheroids_ => null()
     double precision                                :: Ares
   contains
     final     ::                          gasMassRateNSCDestructor
     procedure :: differentialEvolution => gasMassRateNSCDifferentialEvolution
  end type nodeOperatorgasMassRateNSC
  
  interface nodeOperatorgasMassRateNSC
     !!{
     Constructors for the {\normalfont \ttfamily gasMassRateNSC} node operator class.
     !!}
     module procedure gasMassRateNSCConstructorParameters
     module procedure gasMassRateNSCConstructorInternal
  end interface nodeOperatorgasMassRateNSC
  
contains

  function gasMassRateNSCConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily gasMassRateNSC} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorgasMassRateNSC     )                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(starFormationRateSpheroidsClass), pointer       :: starFormationRateSpheroids_
    double precision                                      :: Ares
    !![
    <inputParameter>
    <name>Ares</name>
    <defaultValue>6.0d-3</defaultValue>
    <description> Free parameter controlling the transference of gas from the spheroid to the NSC gas reservoir</description>
    <source>parameters</source>
    </inputParameter>
    <objectBuilder class="starFormationRateSpheroids" name="starFormationRateSpheroids_" source="parameters"/>
    !!]
    self=nodeOperatorgasMassRateNSC(Ares,starFormationRateSpheroids_)
    !![
    <inputParametersValidate source="parameters"        />
    <objectDestructor name="starFormationRateSpheroids_"/>
    !!]
    return
  end function gasMassRateNSCConstructorParameters
  
  function gasMassRateNSCConstructorInternal(Ares, starFormationRateSpheroids_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily gasMassRateNSC} node operator class.
    !!}
    implicit none
    type (nodeOperatorgasMassRateNSC     )                        :: self
    class(starFormationRateSpheroidsClass), intent(in   ), target :: starFormationRateSpheroids_
    double precision                      , intent(in   ), target :: Ares
    !![
    <constructorAssign variables="Ares"                        />
    <constructorAssign variables="*starFormationRateSpheroids_"/>
    !!]
    return
  end function gasMassRateNSCConstructorInternal
  
  subroutine gasMassRateNSCDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily gasMassRateNSC} node operator class.
    !!}
    implicit none
    type(nodeOperatorgasMassRateNSC), intent(inout) :: self
    !![
    <objectDestructor name="self%starFormationRateSpheroids_"/>
    !!]
    return
  end subroutine gasMassRateNSCDestructor

  subroutine gasMassRateNSCDifferentialEvolution(self,node,interrupt,functioninterrupt,propertyType)
      !!{
        Compute the nuclear star cluster gas mass rate change.
      !!}
    use :: Galacticus_Nodes , only : interruptTask   , nodeComponentNSC, nodeComponentNSCStandard, nodeComponentSpheroid, &
                                 &   propertyInactive, treeNode
    implicit none
    class(nodeOperatorgasMassRateNSC), intent(inout), target :: self
    type (treeNode                  ), intent(inout), target :: node
    logical                          , intent(inout)         :: interrupt
    procedure (interruptTask        ), intent(inout), pointer:: functioninterrupt
    integer                          , intent(in   )         :: propertyType
    class (nodeComponentNSC         ),                pointer:: NSC
    class (nodeComponentSpheroid    ),                pointer:: spheroid
    double precision                                         :: gasMassAccretionRate , rateStarFormationSpheroid

    gasMassAccretionRate = 0.0d0

    ! Return immediately if inactive variables are requested.
    if (propertyInactive(propertyType)) return

    ! Get the spheroid component.
    spheroid => node%spheroid()
  
    ! Get the star formation rate of the spheroid component.
    rateStarFormationSpheroid =  self%starFormationRateSpheroids_%rate(node)   

    ! Find the rate of mass accretion onto the nuclear star cluster.
    if  (rateStarFormationSpheroid <= 0.0d0) then
      gasMassAccretionRate = 0.0d0
    else
      ! Gas mass accretion rate model from F. Antonini, E. Barausse & J. Silk (2015; https://ui.adsabs.harvard.edu/abs/2015ApJ...812...72A/abstract)
      gasMassAccretionRate = self%Ares*rateStarFormationSpheroid
    end if 
    
    ! Finish if there is no accretion.
    if (gasMassAccretionRate <= 0.0d0) return

    ! Get the nuclear stars cluster component.
    NSC   => node%NSC()

    ! Detect nuclear star cluster component type.
    select type (NSC)
    type is (nodeComponentNSC)
      ! Generic type - interrupt and create a standard nuclear star cluster if accretion rate is non-zero.
      if (gasmassAccretionRate /= 0.0d0) then
        interrupt=.true.
        functionInterrupt => NSCCreate
      end if
      return
      ! Standard type - continue processing.
      class is (nodeComponentNSCStandard)
        call NSC%massGasRate(gasMassAccretionRate)
    end select
    return
  end subroutine gasMassRateNSCDifferentialEvolution

  subroutine NSCCreate(node,timeEnd)
    !!{
       Creates the nuclear star cluster via interrupt.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC, treeNode
    implicit none
    type             (treeNode        ), intent(inout), target  :: node
    double precision                   , intent(in   ), optional:: timeEnd
    class            (nodeComponentNSC),                pointer :: NSC
    !$GLC attributes unused :: timeEnd

    ! Creates the nuclear star cluster.
    NSC => node%NSC(autoCreate=.true.)
    return 
  end subroutine NSCCreate
  