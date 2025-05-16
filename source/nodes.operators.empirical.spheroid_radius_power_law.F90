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

!+    Contributions to this file made by: Charles Gannon, Andrew Benson.

  !!{
  Implements a node operator class that implements an an empirical power law relationship between spheroid stellar mass and
  stellar radius.
  !!}

  !![
  <nodeOperator name="nodeOperatorSpheroidRadiusPowerLaw">
   <description>
    A node operator that implements an an empirical power law relationship between spheroid stellar radius and stellar mass.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSpheroidRadiusPowerLaw
     !!{
     Implements a power law prescription for the stellar mass--stellar radius relation of spheroids. Specificially:
     \begin{equation}
       r_\mathrm{s} = \beta \left( \frac{M_\star}{M_\odot} \right)^\alpha, 
     \end{equation}
     where $r_\mathrm{s}$ is the spheroid scale radius, $M_\star$ is the stellar mass of the spheroid, and $\alpha$ and $\beta$
     are free parameters.
     !!}
     private
     double precision :: alpha, beta
   contains
     !![
     <methods>
       <method method="update" description="Update the spheroid radius to be consistent with its stellar mass."/>
     </methods>
     !!]
     procedure :: update                              => spheroidRadiusPowerLawUpdate
     procedure :: nodeInitialize                      => spheroidRadiusPowerLawNodeInitialize
     procedure :: differentialEvolutionSolveAnalytics => spheroidRadiusPowerLawSolveAnalytics
     procedure :: nodesMerge                          => spheroidRadiusPowerLawNodesMerge
  end type nodeOperatorSpheroidRadiusPowerLaw
  
  interface nodeOperatorSpheroidRadiusPowerLaw
     !!{
     Constructors for the \refClass{nodeOperatorSpheroidRadiusPowerLaw} node operator class.
     !!}
     module procedure spheroidRadiusPowerLawConstructorParameters
     module procedure spheroidRadiusPowerLawConstructorInternal
  end interface nodeOperatorSpheroidRadiusPowerLaw
  
contains

  function spheroidRadiusPowerLawConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorSpheroidRadiusPowerLaw} which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorSpheroidRadiusPowerLaw)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    double precision                                                    :: alpha     , beta    

    !![
    <inputParameter>
      <name>alpha</name>
      <source>parameters</source>
      <description>Exponent $\alpha$ in the power law fit.</description> 
      <defaultValue>0.56d0</defaultValue>
      <defaultSource>\cite[][table J1: Parameter $a$, for early type galaxies]{shen_erratum_2007}</defaultSource>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <source>parameters</source>
      <description>Coefficient $\beta$ in the power law fit.</description>
      <defaultValue>1.19d-9</defaultValue>
      <defaultSource>\cite[][table J1: Parameter $b$, for early type galaxies and re-scaled from half light radius to Hernquist radius \protect\citep{hernquist_analytical_1990}---Note: there was a typo in the originally provided value, see \protect\cite{shen_erratum_2007} for the corrected value]{shen_size_2003}</defaultSource>
    </inputParameter>
    !!]
    self=spheroidRadiusPowerLawConstructorInternal(alpha, beta)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
  end function spheroidRadiusPowerLawConstructorParameters

  function spheroidRadiusPowerLawConstructorInternal(alpha, beta) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorSpheroidRadiusPowerLaw} node operator class.
    !!}
    implicit none
    type            (nodeOperatorSpheroidRadiusPowerLaw)             :: self
    double precision                                    , intent(in) :: alpha, beta
    !![
    <constructorAssign variables="alpha, beta"/>
    !!]
    return
  end function spheroidRadiusPowerLawConstructorInternal

  subroutine spheroidRadiusPowerLawUpdate(self, node)
    !!{
    Update radius of the spheroid.
    !!} 
    use :: Galacticus_Nodes, only : nodeComponentSpheroid
    implicit none
    class           (nodeOperatorSpheroidRadiusPowerLaw), intent(inout) :: self
    type            (treeNode                          ), intent(inout) :: node
    class           (nodeComponentSpheroid             ), pointer       :: spheroid
    double precision                                                    :: radiusStellar

    if (.not.node%isOnMainBranch()) return 
    spheroid      =>  node    %spheroid   ()
    radiusStellar =  +self    %beta                      &
      &              *spheroid%massStellar()**self%alpha 
    call spheroid%radiusSet(radiusStellar)
    return  
  end subroutine spheroidRadiusPowerLawUpdate
   
  subroutine spheroidRadiusPowerLawNodeInitialize(self,node)
    !!{
    Initialize the spheroid.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpheroid
    implicit none
    class(nodeOperatorSpheroidRadiusPowerLaw), intent(inout), target  :: self
    type (treeNode                          ), intent(inout), target  :: node
    class           (nodeComponentSpheroid  ),                pointer :: spheroid

    ! Initialize spheroid.
    if (.not.node%isOnMainBranch().or.associated(node%firstChild)) return 
    spheroid => node%spheroid(autoCreate=.true.)
    call self%update(node)
    return
  end subroutine spheroidRadiusPowerLawNodeInitialize

  subroutine spheroidRadiusPowerLawSolveAnalytics(self,node,time)
    !!{
    Set the radius of the spheroid.
    !!}
    implicit none
    class           (nodeOperatorSpheroidRadiusPowerLaw), intent(inout) :: self
    type            (treeNode                          ), intent(inout) :: node
    double precision                                    , intent(in   ) :: time
    !$GLC attributes unused :: time

    call self%update(node)
    return
  end subroutine spheroidRadiusPowerLawSolveAnalytics

  subroutine spheroidRadiusPowerLawNodesMerge(self,node)
    !!{
    Update the radius of the spheroid after a merger.
    !!}
    implicit none
    class(nodeOperatorSpheroidRadiusPowerLaw), intent(inout) :: self
    type (treeNode                          ), intent(inout) :: node

    call self%update(node)
    return
  end subroutine spheroidRadiusPowerLawNodesMerge
