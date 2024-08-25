!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
  Implements a node operator class that inserts an empirical model of the formation history of a massive elliptical galaxy.
  !!}

  !![
  <nodeOperator name="nodeOperatorspheroidRadiusPowerLaw">
   <description>
    An empirical power law for the Stellar Mass - Stellar Radius of a spheroid.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorspheroidRadiusPowerLaw
     !!{
      Node operator that sets a power law prescription for the Stellar Mass - Stellar Radius of spheroid.
      The power law is given by:
      \begin{equation}
          r_{hq} = \beta \left( \frac{M_\star}{M_\odot} \right)^\alpha
      \end{equation}
     !!}
     private
     double precision :: alpha, beta
   contains
     procedure :: nodeInitialize                      => spheroidRadiusPowerLawNodeInitialize
     procedure :: differentialEvolutionSolveAnalytics => spheroidRadiusPowerLawSolveAnalytics
     procedure :: nodesMerge                          => spheroidRadiusPowerLawNodesMerge
  end type nodeOperatorspheroidRadiusPowerLaw
  
  interface nodeOperatorspheroidRadiusPowerLaw
     !!{
     Constructors for the {\normalfont \ttfamily spheroidRadiusPowerLaw} node operator class.
     !!}
     module procedure spheroidRadiusPowerLawConstructorParameters
     module procedure spheroidRadiusPowerLawConstructorInternal
  end interface nodeOperatorspheroidRadiusPowerLaw
  
contains

  function spheroidRadiusPowerLawConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily starFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorspheroidRadiusPowerLaw)                    :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    double precision                                                        :: alpha     , beta    

    !![
    <inputParameter>
      <name>alpha</name>
      <source>parameters</source>
      <description>
        Exponent $\alpha$ in the power law fit. 
        Best fit value taken from \cite{2003MNRAS.343..978S} (table J1: Parameter a, for early type galaxies). 
      </description> 
      <defaultValue>0.56d0</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <source>parameters</source>
      <description>
        Coefitient $\beta$ in the power law fit.
        Best fit value taken from \cite{2003MNRAS.343..978S} (table J1: Parameter b, for early type galaxies) and rescaled from half light radius to Hernquist radius \citep(hernquist_analytical_1990).
        Note: there was a typo in the originally provided value see \cite{hernquist_analytical_1990} for the corrected value.        
        Additionaly, 
      </description>
    
      <defaultValue>1.19d-6</defaultValue>
    </inputParameter>
    !!]
  end function spheroidRadiusPowerLawConstructorParameters

  function spheroidRadiusPowerLawConstructorInternal(alpha, beta) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily spheroidRadiusPowerLaw} node operator class.
    !!}
    implicit none
    type            (nodeOperatorspheroidRadiusPowerLaw)             :: self
    double precision                                    , intent(in) :: alpha, beta
    !![
    <constructorAssign variables="alpha, beta"/>
    !!]
    return
  end function spheroidRadiusPowerLawConstructorInternal

  subroutine spheroidRadiusPowerLawNodeUpdate(self, node)
    !!{
    Update radius of galaxy
    !!} 
    use :: Galacticus_Nodes, only : nodeComponentSpheroid
    implicit none
    class (nodeOperatorspheroidRadiusPowerLaw), intent(inout)          :: self
    type  (treeNode                          ), intent(inout)          :: node
    class (nodeComponentSpheroid             )               , pointer :: spheroid
    double precision                                                   :: radiusStellar

    if (.not.node%isOnMainBranch()) return 

    spheroid     =>node%spheroid()
    
    radiusStellar= + self%beta              &
      &            + spheroid%massStellar() &
      &            **self%alpha 

    call spheroid%radiusSet(radiusStellar)

    return  
  end subroutine spheroidRadiusPowerLawNodeUpdate
   
  subroutine spheroidRadiusPowerLawNodeInitialize(self,node)
    !!{
    Initialize radii of galaxy
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpheroid
    implicit none
    class           (nodeOperatorspheroidRadiusPowerLaw), intent(inout), target  :: self
    type            (treeNode                          ), intent(inout), target  :: node
    class           (nodeComponentSpheroid             )               , pointer :: spheroid

    if (.not.node%isOnMainBranch()) return 

    spheroid=>node%spheroid(autoCreate=.true.)

    call spheroidRadiusPowerLawNodeUpdate(self, node)

    return
  end subroutine spheroidRadiusPowerLawNodeInitialize

  subroutine spheroidRadiusPowerLawSolveAnalytics(self,node,time)
    !!{
    Set radius of galaxy
    !!}
    implicit none
    
    class           (nodeOperatorspheroidRadiusPowerLaw), intent(inout) :: self
    type            (treeNode                          ), intent(inout) :: node
    double precision                                    , intent(in   ) :: time

    call spheroidRadiusPowerLawNodeUpdate(self, node)

    return
  end subroutine spheroidRadiusPowerLawSolveAnalytics

  subroutine spheroidRadiusPowerLawNodesMerge(self,node)
    !!{
    Update radius of galaxy after merges
    !!}
    implicit none
    class           (nodeOperatorspheroidRadiusPowerLaw), intent(inout) :: self
    type            (treeNode                          ), intent(inout) :: node

    call spheroidRadiusPowerLawNodeUpdate(self, node)

    return
  end subroutine spheroidRadiusPowerLawNodesMerge
