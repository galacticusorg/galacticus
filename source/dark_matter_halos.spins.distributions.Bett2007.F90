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

!% Contains a module which implements the \cite{bett_spin_2007} halo spin distribution.

module Halo_Spin_Distributions_Bett2007
  !% Implements the \cite{bett_spin_2007} halo spin distribution.
  use FGSL
  use Galacticus_Nodes
  use Tables
  implicit none
  private
  public :: Halo_Spin_Distribution_Bett2007_Initialize, Halo_Spin_Distribution_Bett2007_Snapshot,&
       & Halo_Spin_Distribution_Bett2007_State_Store, Halo_Spin_Distribution_Bett2007_State_Retrieve

  ! Parameters of the spin distribution.
  double precision                                        :: spinDistributionBett2007Alpha           , spinDistributionBett2007Lambda0

  ! Tabulation of the spin distribution.
  integer                                   , parameter   :: spinDistributionTableNumberPoints=1000
  double precision                          , parameter   :: spinDistributionTableSpinMaximum =0.2d0                                   !   Maximum spin to tabulate.
  double precision                          , parameter   :: spinDistributionTableMinimum     =1.0d-6                                  !   Minimum spin in units of lambda_.
  double precision                                        :: spinDistributionTableMaximum
  type            (table1DLogarithmicLinear)              :: spinDistributionTable
  class           (table1D                 ), allocatable :: spinDistributionTableInverse

  ! Random number objects.
  type            (fgsl_rng                )              :: clonedPseudoSequenceObject              , randomSequenceObject
  logical                                                 :: resetRandomSequence              =.true., resetRandomSequenceSnapshot
  !$omp threadprivate(resetRandomSequence,randomSequenceObject,resetRandomSequenceSnapshot,clonedPseudoSequenceObject)

contains

  !# <haloSpinDistributionMethod>
  !#  <unitName>Halo_Spin_Distribution_Bett2007_Initialize</unitName>
  !# </haloSpinDistributionMethod>
  subroutine Halo_Spin_Distribution_Bett2007_Initialize(haloSpinDistributionMethod,Halo_Spin_Sample_Get)
    !% Initializes the ``Bett2007'' halo spin distribution module.
    use ISO_Varying_String
    use Input_Parameters
    use Memory_Management
    use Numerical_Ranges
    use Gamma_Functions
    implicit none
    type     (varying_string                 ), intent(in   )          :: haloSpinDistributionMethod
    procedure(Halo_Spin_Distribution_Bett2007), intent(inout), pointer :: Halo_Spin_Sample_Get
    integer                                                            :: iSpin
    double precision                                                   :: spinDimensionless

    if (haloSpinDistributionMethod == 'Bett2007') then
       Halo_Spin_Sample_Get => Halo_Spin_Distribution_Bett2007
       !@ <inputParameter>
       !@   <name>spinDistributionBett2007Lambda0</name>
       !@   <defaultValue>0.04326 \citep{bett_spin_2007}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The median in a lognormal halo spin distribution.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('spinDistributionBett2007Lambda0',spinDistributionBett2007Lambda0,defaultValue=0.04326d0)
       !@ <inputParameter>
       !@   <name>spinDistributionBett2007Alpha</name>
       !@   <defaultValue>2.509 \citep{bett_spin_2007}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The dispersion in a lognormal halo spin distribution.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('spinDistributionBett2007Alpha',spinDistributionBett2007Alpha,defaultValue=2.509d0)

       ! Maximum value of x=(lambda/lambda_0)^(3/alpha) to tabulate.
       spinDistributionTableMaximum=(spinDistributionTableSpinMaximum/spinDistributionBett2007Lambda0)**(3.0d0/spinDistributionBett2007Alpha)
       ! Tabulate the cumulative distribution.
       call spinDistributionTable%destroy()
       call spinDistributiontable%create(&
            &                            spinDistributionBett2007Lambda0*spinDistributionTableMinimum**(spinDistributionBett2007Alpha/3.0d0), &
            &                            spinDistributionBett2007Lambda0*spinDistributionTableMaximum**(spinDistributionBett2007Alpha/3.0d0), &
            &                            spinDistributionTableNumberPoints                                                                  , &
            &                            extrapolationType=extrapolationTypeFix                                                               &
            &                           )
       ! Compute the cumulative probability distribution.
       do iSpin=1,spinDistributionTableNumberPoints
          spinDimensionless=                                         &
               &            (                                        &
               &              spinDistributionTable%x(iSpin)         &
               &             /spinDistributionBett2007Lambda0        &
               &            )**(3.0d0/spinDistributionBett2007Alpha)
          call spinDistributionTable%populate(                                                                        &
               &                              Gamma_Function_Incomplete_Complementary                                 &
               &                                                                     (                                &
               &                                                                       spinDistributionBett2007Alpha, &
               &                                                                       spinDistributionBett2007Alpha  &
               &                                                                      *spinDimensionless              &
               &                                                                     )                              , &
               &                              iSpin                                                                   &
               &                             )
       end do
       call spinDistributionTable%reverse(spinDistributionTableInverse)
    end if
    return
  end subroutine Halo_Spin_Distribution_Bett2007_Initialize

  double precision function Halo_Spin_Distribution_Bett2007(thisNode)
    !% Return a halo spin from a lognormal distribution.
    use Galacticus_Nodes
    use Pseudo_Random
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision                                   :: randomDeviate

    randomDeviate=Pseudo_Random_Get(randomSequenceObject,resetRandomSequence)
    Halo_Spin_Distribution_Bett2007=spinDistributionTableInverse%interpolate(randomDeviate)
    return
  end function Halo_Spin_Distribution_Bett2007

  !# <galacticusStateSnapshotTask>
  !#  <unitName>Halo_Spin_Distribution_Bett2007_Snapshot</unitName>
  !# </galacticusStateSnapshotTask>
  subroutine Halo_Spin_Distribution_Bett2007_Snapshot
    !% Store a snapshot of the random number generator internal state.
    implicit none

    if (.not.resetRandomSequence) clonedPseudoSequenceObject=FGSL_Rng_Clone(randomSequenceObject)
    resetRandomSequenceSnapshot=resetRandomSequence
    return
  end subroutine Halo_Spin_Distribution_Bett2007_Snapshot

  !# <galacticusStateStoreTask>
  !#  <unitName>Halo_Spin_Distribution_Bett2007_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Halo_Spin_Distribution_Bett2007_State_Store(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    use Pseudo_Random
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    write (stateFile) resetRandomSequenceSnapshot
    if (.not.resetRandomSequenceSnapshot) call Pseudo_Random_Store(clonedPseudoSequenceObject,fgslStateFile)
    return
  end subroutine Halo_Spin_Distribution_Bett2007_State_Store

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Halo_Spin_Distribution_Bett2007_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Halo_Spin_Distribution_Bett2007_State_Retrieve(stateFile,fgslStateFile)
    !% Write the stored snapshot of the random number state to file.
    use FGSL
    use Pseudo_Random
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    read (stateFile) resetRandomSequence
    if (.not.resetRandomSequence) call Pseudo_Random_Retrieve(randomSequenceObject,fgslStateFile)
    return
  end subroutine Halo_Spin_Distribution_Bett2007_State_Retrieve

end module Halo_Spin_Distributions_Bett2007
