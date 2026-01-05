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
  Implements a node operator class that implements an an empirical power law relationship between disk stellar mass and
  stellar radius.
  !!}

  !![
  <nodeOperator name="nodeOperatorDiskRadiusPowerLaw">
   <description>
    A node operator that implements an an empirical power law relationship between disk stellar radius and stellar mass.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorDiskRadiusPowerLaw
     !!{
     Implements a power law prescription for the stellar mass--stellar radius relation of disks. Specifically:
     \begin{equation}
       r_\mathrm{s} = \gamma \left( \frac{M_\star}{M_\odot} \right)^\alpha \left( 1 + \frac{M_\star}{M_0}^{\beta-\alpha} \right), 
     \end{equation}
     where $r_\mathrm{s}$ is the disk scale radius, $M_\star$ is the stellar mass of the disk, and $M_0$, $\alpha$, $\beta$, and $\gamma$
     are free parameters.
     !!}
     private
     double precision :: alpha, beta     , &
          &              gamma, massPivot
   contains
     !![
     <methods>
       <method method="update" description="Update the disk radius to be consistent with its stellar mass."/>
     </methods>
     !!]
     procedure :: update                              => diskRadiusPowerLawUpdate
     procedure :: nodeInitialize                      => diskRadiusPowerLawNodeInitialize
     procedure :: differentialEvolutionSolveAnalytics => diskRadiusPowerLawSolveAnalytics
     procedure :: nodesMerge                          => diskRadiusPowerLawNodesMerge
  end type nodeOperatorDiskRadiusPowerLaw
  
  interface nodeOperatorDiskRadiusPowerLaw
     !!{
     Constructors for the \refClass{nodeOperatorDiskRadiusPowerLaw} node operator class.
     !!}
     module procedure diskRadiusPowerLawConstructorParameters
     module procedure diskRadiusPowerLawConstructorInternal
  end interface nodeOperatorDiskRadiusPowerLaw
  
contains

  function diskRadiusPowerLawConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorDiskRadiusPowerLaw} class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorDiskRadiusPowerLaw)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    double precision                                                :: alpha     , beta     , &
         &                                                             gamma     , massPivot

    !![
    <inputParameter>
      <name>alpha</name>
      <source>parameters</source>
      <description>Exponent $\alpha$ in the power law fit.</description> 
      <defaultValue>0.14d0</defaultValue>
      <defaultSource>\cite[][table 1: Parameter $\alpha$, for late type galaxies]{shen_size_2003}</defaultSource>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <source>parameters</source>
      <description>Exponent $\beta$ in the power law fit.</description>
      <defaultValue>0.39d0</defaultValue>
      <defaultSource>\cite[][table 1: Parameter $\beta$, for late type galaxies]{shen_size_2003}</defaultSource>
    </inputParameter>
    <inputParameter>
      <name>gamma</name>
      <source>parameters</source>
      <description>Coefficient $\gamma$ in the power law fit.</description>
      <defaultValue>5.958d-5</defaultValue>
      <defaultSource>\cite[][table 1: Parameter $\gamma$, for late type galaxies and re-scaled from half light radius to exponential radius]{shen_size_2003}</defaultSource>
    </inputParameter>
    <inputParameter>
      <name>massPivot</name>
      <source>parameters</source>
      <description>Pivot mass $M_0$ in the power law fit.</description>
      <defaultValue>3.98d10</defaultValue>
      <defaultSource>\cite[][table 1: Parameter $M_0$, for late type galaxies]{shen_size_2003}</defaultSource>
    </inputParameter>
    !!]
    self=diskRadiusPowerLawConstructorInternal(alpha,beta,gamma,massPivot)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
  end function diskRadiusPowerLawConstructorParameters

  function diskRadiusPowerLawConstructorInternal(alpha,beta,gamma,massPivot) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorDiskRadiusPowerLaw} node operator class.
    !!}
    implicit none
    type            (nodeOperatorDiskRadiusPowerLaw)             :: self
    double precision                                , intent(in) :: alpha, beta     , &
         &                                                          gamma, massPivot
    !![
    <constructorAssign variables="alpha, beta, gamma, massPivot"/>
    !!]
    return
  end function diskRadiusPowerLawConstructorInternal

  subroutine diskRadiusPowerLawUpdate(self,node)
    !!{
    Update radius of the disk.
    !!} 
    use :: Galacticus_Nodes, only : nodeComponentDisk
    implicit none
    class           (nodeOperatorDiskRadiusPowerLaw), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    class           (nodeComponentDisk             ), pointer       :: disk
    double precision                                                :: radiusStellar

    if (.not.node%isOnMainBranch()) return 
    disk          =>  node    %disk   ()
    radiusStellar =  +self    %gamma                                                &
         &           *  (disk%massStellar()               )**           self%alpha  &
         &           *(                                                             &
         &             +1.0d0                                                       &
         &             +(disk%massStellar()/self%massPivot)**(self%beta-self%alpha) &
         &            )
    call disk%radiusSet(radiusStellar)
    return  
  end subroutine diskRadiusPowerLawUpdate
   
  subroutine diskRadiusPowerLawNodeInitialize(self,node)
    !!{
    Initialize the disk.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk
    implicit none
    class(nodeOperatorDiskRadiusPowerLaw), intent(inout), target  :: self
    type (treeNode                      ), intent(inout), target  :: node
    class(nodeComponentDisk             ),                pointer :: disk

    ! Initialize disk.
    if (.not.node%isOnMainBranch().or.associated(node%firstChild)) return 
    disk => node%disk(autoCreate=.true.)
    call self%update(node)
    return
  end subroutine diskRadiusPowerLawNodeInitialize

  subroutine diskRadiusPowerLawSolveAnalytics(self,node,time)
    !!{
    Set the radius of the disk.
    !!}
    implicit none
    class           (nodeOperatorDiskRadiusPowerLaw), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(in   ) :: time
    !$GLC attributes unused :: time

    call self%update(node)
    return
  end subroutine diskRadiusPowerLawSolveAnalytics

  subroutine diskRadiusPowerLawNodesMerge(self,node)
    !!{
    Update the radius of the disk after a merger.
    !!}
    implicit none
    class(nodeOperatorDiskRadiusPowerLaw), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node

    call self%update(node)
    return
  end subroutine diskRadiusPowerLawNodesMerge
