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
    procedure(),          pointer, intent(inout) :: Halo_Spin_Sample_Get
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
       call Alloc_Array(spinDistributionTableSpin      ,spinDistributionTableNumberPoints,'spinDistributionTableSpin'      )
       call Alloc_Array(spinDistributionTableCumulative,spinDistributionTableNumberPoints,'spinDistributionTableCumulative')
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
