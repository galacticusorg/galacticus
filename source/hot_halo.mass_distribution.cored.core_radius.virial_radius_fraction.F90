!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% An implementation of the hot halo mass distribution core radius class which sets the core radius to a fraction of the virial radius.

  double precision :: virialFractionCoreRadiusOverVirialRadius
  logical          :: virialFractionInitialized               =.false.

  !# <hotHaloMassDistributionCoreRadius name="hotHaloMassDistributionCoreRadiusVirialFraction">
  !#  <description>Provides an implementation of the hot halo mass distribution core radius class which sets the core radius to a fraction of the virial radius.</description>
  !# </hotHaloMassDistributionCoreRadius>
  type, extends(hotHaloMassDistributionCoreRadiusClass) :: hotHaloMassDistributionCoreRadiusVirialFraction
     !% An implementation of the hot halo mass distribution core radius class which sets the core radius to a fraction of the virial radius.
     private
     double precision :: coreRadiusOverVirialRadius
   contains
     final     ::           virialFractionDestructor
     procedure :: radius => virialFractionRadius
  end type hotHaloMassDistributionCoreRadiusVirialFraction

  interface hotHaloMassDistributionCoreRadiusVirialFraction
     !% Constructors for the {\tt virialFraction} hot halo mass distribution core radius class.
     module procedure virialFractionDefaultConstructor
     module procedure virialFractionConstructor
  end interface hotHaloMassDistributionCoreRadiusVirialFraction

contains

  function virialFractionDefaultConstructor()
    !% Default constructor for the {\tt virialFraction} hot halo mass distribution core radius class.
    use Input_Parameters
    implicit none
    type(hotHaloMassDistributionCoreRadiusVirialFraction) :: virialFractionDefaultConstructor
    
    if (.not.virialFractionInitialized) then
       !$omp critical (hotHaloMassDistributionCoreRadiusVirialFractionInitialize)
       if (.not.virialFractionInitialized) then
          !@ <inputParameter>
          !@   <name>hotHaloCoreRadiusOverVirialRadius</name>
          !@   <defaultValue>0.3</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The core radius in the hot halo density profile in units of the virial radius.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('hotHaloCoreRadiusOverVirialRadius',virialFractionCoreRadiusOverVirialRadius,defaultValue=0.3d0)
          ! Record that this implementation is now initialized.
          virialFractionInitialized=.true.
       end if
       !$omp end critical (hotHaloMassDistributionCoreRadiusVirialFractionInitialize)
    end if
    virialFractionDefaultConstructor=virialFractionConstructor(virialFractionCoreRadiusOverVirialRadius)
    return
  end function virialFractionDefaultConstructor

  function virialFractionConstructor(coreRadiusOverVirialRadius)
    !% Default constructor for the {\tt virialFraction} hot halo mass distribution core radius class.
    use Input_Parameters
    implicit none
    type            (hotHaloMassDistributionCoreRadiusVirialFraction)                :: virialFractionConstructor
    double precision                                                 , intent(in   ) :: coreRadiusOverVirialRadius

    virialFractionConstructor%coreRadiusOverVirialRadius=coreRadiusOverVirialRadius
   return
  end function virialFractionConstructor

  elemental subroutine virialFractionDestructor(self)
    !% Destructor for the {\tt virialFraction} hot halo mass distribution class.
    implicit none
    type(hotHaloMassDistributionCoreRadiusVirialFraction), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine virialFractionDestructor

  double precision function virialFractionRadius(self,node)
    !% Return the core radius of the hot halo mass distribution.
    use Dark_Matter_Halo_Scales
    implicit none
    class           (hotHaloMassDistributionCoreRadiusVirialFraction), intent(inout)          :: self
    type            (treeNode                                       ), intent(inout), pointer :: node

    virialFractionRadius=self%coreRadiusOverVirialRadius*Dark_Matter_Halo_Virial_Radius(node)
    return
  end function virialFractionRadius
