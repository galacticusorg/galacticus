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

!% Contains a module which implements calculations of dark matter halo angular momentum.

module Dark_Matter_Halo_Spins
  !% Implements calculations of dark matter halo angular momentum.
  implicit none
  private
  public :: Dark_Matter_Halo_Angular_Momentum, Dark_Matter_Halo_Angular_Momentum_Growth_Rate

  ! Record of whether the module has been initialized.
  logical :: moduleInitialized=.false.

contains

  subroutine Dark_Matter_Halo_Spins_Initialize()
    !% Initialize the halo spins module.
    use Galacticus_Nodes
    use Galacticus_Error
    implicit none

    if (.not.moduleInitialized) then
       !$omp critical(Dark_Matter_Halo_Spins_Initialize)
       if (.not.moduleInitialized) then
          ! Ensure that the spin property is available.
          if (.not. defaultSpinComponent%spinIsGettable()) call Galacticus_Error_Report('Dark_Matter_Halo_Spins_Initialize','spin property of spin component must be gettable')
          if (.not.defaultBasicComponent%massIsGettable()) call Galacticus_Error_Report('Dark_Matter_Halo_Spins_Initialize','mass property of basic component must be gettable')
          ! Record that the module is now initialized.
          moduleInitialized=.true.
       end if
       !$omp end critical(Dark_Matter_Halo_Spins_Initialize)
    end if
    return
  end subroutine Dark_Matter_Halo_Spins_Initialize

  double precision function Dark_Matter_Halo_Angular_Momentum(thisNode)
    !% Returns the total anuglar momentum of {\tt thisNode} based on its mass, energy and spin parameter.
    use Galacticus_Nodes
    use Numerical_Constants_Physical
    use Dark_Matter_Profiles
    implicit none
    type (treeNode          ), pointer, intent(inout) :: thisNode
    class(nodeComponentBasic), pointer                :: thisBasicComponent
    class(nodeComponentSpin ), pointer                :: thisSpinComponent

    ! Ensure that the module is initialized.
    call Dark_Matter_Halo_Spins_Initialize

    thisBasicComponent => thisNode%basic(                 )
    thisSpinComponent  => thisNode%spin (autoCreate=.true.)
    Dark_Matter_Halo_Angular_Momentum=thisSpinComponent%spin()*gravitationalConstantGalacticus*thisBasicComponent%mass()**2.5d0 &
         &/dsqrt(dabs(Dark_Matter_Profile_Energy(thisNode)))
    return
  end function Dark_Matter_Halo_Angular_Momentum

  double precision function Dark_Matter_Halo_Angular_Momentum_Growth_Rate(thisNode)
    !% Returns the rate of change of the total anuglar momentum of {\tt thisNode} based on its mass, energy and spin parameter.
    use Galacticus_Nodes
    use Dark_Matter_Profiles
    implicit none
    type (treeNode          ), pointer, intent(inout) :: thisNode
    class(nodeComponentBasic), pointer                :: thisBasicComponent
    class(nodeComponentSpin ), pointer                :: thisSpinComponent

    ! Ensure that the module is initialized.
    call Dark_Matter_Halo_Spins_Initialize

    thisBasicComponent => thisNode%basic(                 )
    thisSpinComponent  => thisNode%spin (autoCreate=.true.)

    Dark_Matter_Halo_Angular_Momentum_Growth_Rate=Dark_Matter_Halo_Angular_Momentum(thisNode)&
         &*(thisSpinComponent%spinGrowthRate()/thisSpinComponent%spin()+2.5d0*thisBasicComponent%accretionRate()&
         &/thisBasicComponent%mass()-0.5d0*Dark_Matter_Profile_Energy_Growth_Rate(thisNode)/Dark_Matter_Profile_Energy(thisNode))

    return
  end function Dark_Matter_Halo_Angular_Momentum_Growth_Rate

end module Dark_Matter_Halo_Spins
