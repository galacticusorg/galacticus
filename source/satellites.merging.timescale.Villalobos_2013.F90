!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% Implements calculations of satellite merging times by applying the \cite{villalobos_improved_2013} modifier to another
  !% selected satellite merging time method.

  !# <satelliteMergingTimescales name="satelliteMergingTimescalesVillalobos2013" />

  type, extends(satelliteMergingTimescalesClass) :: satelliteMergingTimescalesVillalobos2013
     !% A class implementing calculations of satellite merging times by applying the \cite{villalobos_improved_2013} modifier to
     !% another selected satellite merging time method.
     private
     class           (satelliteMergingTimescalesClass), pointer :: baseMethod
     double precision                                           :: expansionFactorExponent
   contains
     final     ::                     villalobos2013Destructor
     procedure :: timeUntilMerging => villalobos2013TimeUntilMerging
  end type satelliteMergingTimescalesVillalobos2013

  interface satelliteMergingTimescalesVillalobos2013
     !% Constructors for the \cite{villalobos_improved_2013} merging timescale class.
     module procedure villalobos2013DefaultConstructor
     module procedure villalobos2013GenericConstructor
  end interface satelliteMergingTimescalesVillalobos2013

  ! Record of whether method is initialized.
  logical                          :: villalobos2013Initialized=.false.
  ! Default settings for this method.
  type            (varying_string) :: satelliteMergingTimescaleVillalobos2013BaseMethod
  double precision                 :: satelliteMergingTimescaleVillalobos2013Exponent

contains

  function villalobos2013DefaultConstructor()
    !% Default constructor for the \cite{villalobos_improved_2013} merging timescale class.
    use Galacticus_Display
    use Input_Parameters
    implicit none
    type(satelliteMergingTimescalesVillalobos2013) :: villalobos2013DefaultConstructor

    ! Initialize if necessary.
    if (.not.villalobos2013Initialized) then
       !$omp critical(villalobos2013Initialization)
       if (.not.villalobos2013Initialized) then
          ! Read parameter controlling the default base method.
          !@ <inputParameter>
          !@   <name>satelliteMergingTimescaleVillalobos2013BaseMethod</name>
          !@   <defaultValue>boylanKolchin2008</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The base {\tt satelliteMergingTimescales} method to which the \cite{villalobos_improved_2013} modifier for satellite merging timescales should be applied.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>starFormation</group>
          !@ </inputParameter>
          call Get_Input_Parameter('satelliteMergingTimescaleVillalobos2013BaseMethod',satelliteMergingTimescaleVillalobos2013BaseMethod,defaultValue="boylanKolchin2008")
          ! Read parameter controlling the default scaling exponent.
          !@ <inputParameter>
          !@   <name>satelliteMergingTimescaleVillalobos2013Exponent</name>
          !@   <defaultValue>$0.44$ \citep{villalobos_improved_2013}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The exponent of $1+z$ appearing in the \cite{villalobos_improved_2013} modifier for satellite merging timescales.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>starFormation</group>
          !@ </inputParameter>
          call Get_Input_Parameter('satelliteMergingTimescaleVillalobos2013Exponent',satelliteMergingTimescaleVillalobos2013Exponent,defaultValue=0.44d0)
          ! Record that this method is now initialized
          villalobos2013Initialized=.true.
       end if
       !$omp end critical(villalobos2013Initialization)
    end if
    ! Construct the default object.
    villalobos2013DefaultConstructor%baseMethod              => satelliteMergingTimescales(char(satelliteMergingTimescaleVillalobos2013BaseMethod))
    villalobos2013DefaultConstructor%expansionFactorExponent =  satelliteMergingTimescaleVillalobos2013Exponent
    return
  end function villalobos2013DefaultConstructor

  function villalobos2013GenericConstructor(baseMethod,expansionFactorExponent)
    !% Generic constructor for the \cite{villalobos_improved_2013} merging timescale class.
    implicit none
    type            (satelliteMergingTimescalesVillalobos2013)                        :: villalobos2013GenericConstructor
    class           (satelliteMergingTimescalesClass         ), intent(in   ), target :: baseMethod
    double precision                                          , intent(in   )         :: expansionFactorExponent

    ! Construct the object.
    villalobos2013GenericConstructor%baseMethod              => baseMethod
    villalobos2013GenericConstructor%expansionFactorExponent =  expansionFactorExponent
    return
  end function villalobos2013GenericConstructor

  elemental subroutine villalobos2013Destructor(self)
    !% Default constructor for the \cite{villalobos_improved_2013} merging timescale class.
    implicit none
    type(satelliteMergingTimescalesVillalobos2013), intent(inout) :: self

    ! Deallocate the base method.
    deallocate(self%baseMethod)
    return
  end subroutine villalobos2013Destructor

  double precision function villalobos2013TimeUntilMerging(self,thisNode,thisOrbit)
    !% Return the timescale for merging satellites using the \cite{villalobos_improved_2013} method.
    use Galacticus_Nodes
    use Dynamical_Friction_Timescale_Utilities
    use Cosmology_Functions
    use Kepler_Orbits
    implicit none
    class           (satelliteMergingTimescalesVillalobos2013 ), intent(inout)          :: self
    type            (treeNode                                 ), intent(inout), pointer :: thisNode
    type            (keplerOrbit                              ), intent(inout)          :: thisOrbit
    class           (nodeComponentBasic                       )               , pointer :: thisBasic
    class           (cosmologyFunctionsClass                  )               , pointer :: cosmologyFunctionsDefault
    double precision                                                                    :: expansionFactor

    ! Get the default cosmology functions object.
    cosmologyFunctionsDefault => cosmologyFunctions()
    ! Compute expansion factor.
    thisBasic       => thisNode%basic()
    expansionFactor =  cosmologyFunctionsDefault%expansionFactor(thisBasic%time())
    ! Compute dynamical friction timescale.
    villalobos2013TimeUntilMerging=self%baseMethod%timeUntilMerging(thisNode,thisOrbit)/expansionFactor**self%expansionFactorExponent
    return
  end function villalobos2013TimeUntilMerging
