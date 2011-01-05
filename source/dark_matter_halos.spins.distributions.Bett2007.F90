!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements the \cite{bett_spin_2007} halo spin distribution.

module Halo_Spin_Distributions_Bett2007
  !% Implements the \cite{bett_spin_2007} halo spin distribution.
  use FGSL
  use Tree_Nodes
  private
  public :: Halo_Spin_Distribution_Bett2007_Initialize, Halo_Spin_Distribution_Bett2007_Snapshot,&
       & Halo_Spin_Distribution_Bett2007_State_Store, Halo_Spin_Distribution_Bett2007_State_Retrieve

  ! Parameters of the spin distribution.
  double precision :: spinDistributionBett2007Lambda0,spinDistributionBett2007Alpha

  ! Tabulation of the spin distribution.
  integer,          parameter                 :: spinDistributionTableNumberPoints=1000
  double precision, parameter                 :: spinDistributionTableSpinMaximum =0.2d0  ! Maximum spin to tabulate.
  double precision, parameter                 :: spinDistributionTableMinimum     =1.0d-6 ! Minimum spin in units of lambda_.
  double precision                            :: spinDistributionTableMaximum
  double precision, allocatable, dimension(:) :: spinDistributionTableSpin,spinDistributionTableCumulative

  ! Random number objects.
  type(fgsl_rng)                              :: randomSequenceObject,clonedPseudoSequenceObject
  logical                                     :: resetRandomSequence=.true.,resetRandomSequenceSnapshot
  !$omp threadprivate(resetRandomSequence,randomSequenceObject)

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
    type(varying_string),          intent(in)    :: haloSpinDistributionMethod
    procedure(double precision), pointer, intent(inout) :: Halo_Spin_Sample_Get
    integer                                      :: iSpin

    if (haloSpinDistributionMethod == 'Bett2007') then
       Halo_Spin_Sample_Get => Halo_Spin_Distribution_Bett2007
       !@ <inputParameter>
       !@   <name>spinDistributionBett2007Lambda0</name>
       !@   <defaultValue>0.04326 \citep{bett_spin_2007}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The median in a lognormal halo spin distribution.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('spinDistributionBett2007Lambda0',spinDistributionBett2007Lambda0,defaultValue=0.04326d0)
       !@ <inputParameter>
       !@   <name>spinDistributionBett2007Alpha</name>
       !@   <defaultValue>2.509 \citep{bett_spin_2007}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The dispersion in a lognormal halo spin distribution.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('spinDistributionBett2007Alpha',spinDistributionBett2007Alpha,defaultValue=2.509d0)

       ! Tabulate the cumulative distribution.
       call Alloc_Array(spinDistributionTableSpin      ,[spinDistributionTableNumberPoints])
       call Alloc_Array(spinDistributionTableCumulative,[spinDistributionTableNumberPoints])
       ! Maximum value of x=(lambda/lambda_0)^(3/alpha) to tabulate.
       spinDistributionTableMaximum=(spinDistributionTableSpinMaximum/spinDistributionBett2007Lambda0)**(3.0d0/spinDistributionBett2007Alpha)
       ! Generate a range of spins.
       spinDistributionTableSpin=Make_Range(spinDistributionTableMinimum,spinDistributionTableMaximum&
            &,spinDistributionTableNumberPoints,rangeType=rangeTypeLogarithmic)
       ! Compute the cumulative probability distribution.
       do iSpin=1,spinDistributionTableNumberPoints
          spinDistributionTableCumulative(iSpin)=Gamma_Function_Incomplete_Complementary(spinDistributionBett2007Alpha&
               &,spinDistributionBett2007Alpha*spinDistributionTableSpin(iSpin))
       end do
       ! Convert the dimensionless quantity into an actual spin.
       spinDistributionTableSpin=spinDistributionBett2007Lambda0*(spinDistributionTableSpin**(spinDistributionBett2007Alpha/3.0d0))
    end if
    return
  end subroutine Halo_Spin_Distribution_Bett2007_Initialize

  double precision function Halo_Spin_Distribution_Bett2007(thisNode)
    !% Return a halo spin from a lognormal distribution.
    use Tree_Nodes
    use Pseudo_Random
    use Numerical_Interpolation
    implicit none
    type(treeNode),          intent(inout), pointer :: thisNode
    type(fgsl_interp),       save                   :: interpolationObject
    type(fgsl_interp_accel), save                   :: interpolationAccelerator
    logical,                 save                   :: resetInterpolation=.true.
    !$omp threadprivate(interpolationObject,interpolationAccelerator,resetInterpolation)
    double precision                                :: randomDeviate

    randomDeviate=Pseudo_Random_Get(randomSequenceObject,resetRandomSequence)
    Halo_Spin_Distribution_Bett2007=Interpolate(spinDistributionTableNumberPoints,spinDistributionTableCumulative&
         &,spinDistributionTableSpin,interpolationObject,interpolationAccelerator,randomDeviate,reset=resetInterpolation)
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
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

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
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    read (stateFile) resetRandomSequence
    if (.not.resetRandomSequence) call Pseudo_Random_Retrieve(randomSequenceObject,fgslStateFile)
    return
  end subroutine Halo_Spin_Distribution_Bett2007_State_Retrieve
  
end module Halo_Spin_Distributions_Bett2007
