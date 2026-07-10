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

  !!{RST
  Implements a merger tree operator class which records and outputs tree processing time information.
  !!}

  use, intrinsic :: ISO_C_Binding, only : c_size_t
  use            :: Kind_Numbers , only : kind_int8

  !![
  <mergerTreeOperator name="mergerTreeOperatorTreeProcessingTimer" docformat="rst">
   <description>
   A merger tree operator class which records and outputs tree processing time information. Tree timing data to be recorded and output to the ``metaData/treeTiming`` group. Three datasets are written to this group:

   ``treeMasses``
      Gives the base node masses of the recorded trees (in units of :math:`\mathrm{M}_\odot`);

   ``treeConstructTimes``
      Gives the time (in seconds) taken to construct each merger tree;

   ``treeEvolveTimes``
      Gives the time (in seconds) taken to evolve each merger tree.

   In addition, a least-squares fit of :math:`\log_{10}` of the total (construction plus evolution) processing time versus :math:`\log_{10}` of the tree mass, and versus :math:`\log_{10}` of the tree node count, is written as the datasets ``fitCoefficientMass``, ``fitResidualMass``, and ``fitRangeMass`` (and the corresponding ``*CountNodes`` datasets). These provide a ready-to-use cost model that a subsequent run of the same model can consume directly via the ``file`` :galacticus-class:`metaTreeProcessingTime` class.
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorTreeProcessingTimer
     !!{RST
     A merger tree operator class which records and outputs tree processing time information.
     !!}
     private
     logical                                                :: collectMemoryUsageData=.false.
     double precision                                       :: timePostEvolution             , timePreConstruction, &
          &                                                    timePreEvolution              , mass
     integer         (kind_int8)                            :: treeID
     integer         (c_size_t )                            :: memoryUsagePeak               , countNodes
     real                                                   :: time
     integer                                                :: countTrees
     double precision           , allocatable, dimension(:) :: timesConstruct                , timesEvolve        , &
          &                                                    masses
     integer         (kind_int8), allocatable, dimension(:) :: treeIDs
     integer         (c_size_t ), allocatable, dimension(:) :: memoryUsagesPeak              , countsNodes
   contains
     final     ::                           treeProcessingTimerDestructor
     procedure :: operatePreConstruction => treeProcessingTimerOperatePreConstruction
     procedure :: operatePreEvolution    => treeProcessingTimerOperatePreEvolution
     procedure :: operatePostEvolution   => treeProcessingTimerOperatePostEvolution
     procedure :: finalize               => treeProcessingTimerFinalize
     procedure :: autoHook               => treeProcessingTimerAutoHook
  end type mergerTreeOperatorTreeProcessingTimer

  interface mergerTreeOperatorTreeProcessingTimer
     !!{RST
     Constructors for the :galacticus-class:`mergerTreeOperatorTreeProcessingTimer` merger tree operator class.
     !!}
     module procedure treeProcessingTimerConstructorParameters
     module procedure treeProcessingTimerConstructorInternal
  end interface mergerTreeOperatorTreeProcessingTimer

  ! Array size increment.
  integer, parameter :: countIncrement=100

contains

  function treeProcessingTimerConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`mergerTreeOperatorTreeProcessingTimer` merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (mergerTreeOperatorTreeProcessingTimer)                :: self
    type   (inputParameters                      ), intent(inout) :: parameters
    logical                                                       :: collectMemoryUsageData
    
    !![
    <inputParameter docformat="rst">
      <name>collectMemoryUsageData</name>
      <defaultValue>.false.</defaultValue>
      <description>
      Specifies whether or not to collect and output data on the memory used while processing trees.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=treeProcessingTimerConstructorInternal(collectMemoryUsageData)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function treeProcessingTimerConstructorParameters

  function treeProcessingTimerConstructorInternal(collectMemoryUsageData) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`mergerTreeOperatorTreeProcessingTimer` merger tree operator class.
    !!}
    implicit none
    type (mergerTreeOperatorTreeProcessingTimer)                :: self
    logical                                     , intent(in   ) :: collectMemoryUsageData
    !![
    <constructorAssign variables="collectMemoryUsageData"/>
    !!]

    self%countTrees=0
    return
  end function treeProcessingTimerConstructorInternal

  subroutine treeProcessingTimerAutoHook(self)
    !!{RST
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : postEvolveEvent, openMPThreadBindingAtLevel
    implicit none
    class(mergerTreeOperatorTreeProcessingTimer), intent(inout) :: self

    if (self%collectMemoryUsageData) call postEvolveEvent%attach(self,treeProcessingTimerPostEvolve,openMPThreadBindingAtLevel,label='mergerTreeOperatorTreeProcessingTimer')
    return
  end subroutine treeProcessingTimerAutoHook

  subroutine treeProcessingTimerDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`mergerTreeOperatorTreeProcessingTimer` merger tree operator class.
    !!}
    use :: Events_Hooks, only : postEvolveEvent
    implicit none
    type(mergerTreeOperatorTreeProcessingTimer), intent(inout) :: self

    if (self%collectMemoryUsageData .and. postEvolveEvent%isAttached(self,treeProcessingTimerPostEvolve)) &
         & call postEvolveEvent%detach(self,treeProcessingTimerPostEvolve)
    return
  end subroutine treeProcessingTimerDestructor

  subroutine treeProcessingTimerOperatePreConstruction(self)
    !!{RST
    Record the CPU time prior to construction of a tree.
    !!}
    !$ use :: OMP_Lib, only : OMP_Get_WTime, OMP_In_Parallel
    implicit none
    class(mergerTreeOperatorTreeProcessingTimer), intent(inout) :: self
    real                                                        :: time
        
    ! Record the CPU time prior to construction.
    !$ if (OMP_In_Parallel()) then
    !$ self%timePreConstruction=OMP_Get_WTime()
    !$ else
    call CPU_Time(time)
    self%timePreConstruction=dble(time)
    !$ end if
    self%timePreEvolution =-1.0
    self%timePostEvolution=-1.0
    self%memoryUsagePeak  = 0_c_size_t
    return
  end subroutine treeProcessingTimerOperatePreConstruction

  subroutine treeProcessingTimerOperatePreEvolution(self,tree)
    !!{RST
    Record the CPU time prior to evolving ``tree``.
    !!}
    use    :: Galacticus_Nodes   , only : mergerTree              , nodeComponentBasic, treeNode
    use    :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
    !$ use :: OMP_Lib            , only : OMP_Get_WTime           , OMP_In_Parallel
    implicit none
    class(mergerTreeOperatorTreeProcessingTimer), intent(inout), target  :: self
    type (mergerTree                           ), intent(inout), target  :: tree
    type (treeNode                             )               , pointer :: node
    class(nodeComponentBasic                   )               , pointer :: basic
    type (treeNode                             )               , pointer :: nodeWork
    type (mergerTreeWalkerAllNodes             )                         :: treeWalker
    real                                                                 :: time
    
    ! Record the mass and ID of the tree.
    node        => tree %nodeBase
    basic       => node %basic   ()
    self%mass   =  basic%mass    ()
    self%treeID =  tree %index
    ! Count nodes in the tree.
    self%countNodes=0_c_size_t
    treeWalker     =mergerTreeWalkerAllNodes(tree,spanForest=.true.)
    do while (treeWalker%next(nodeWork))
       self%countNodes=self%countNodes+1_c_size_t
    end do
    ! Record the CPU time.
    !$ if (OMP_In_Parallel()) then
    !$ self%timePreEvolution=OMP_Get_WTime()
    !$ else
    call CPU_Time(time)
    self%timePreEvolution=dble(time)
    !$ end if
    return
  end subroutine treeProcessingTimerOperatePreEvolution
  
  subroutine treeProcessingTimerOperatePostEvolution(self)
    !!{RST
    Record the CPU time after evolving a tree.
    !!}
    !$ use :: OMP_Lib          , only : OMP_Get_WTime, OMP_In_Parallel
    implicit none
    class           (mergerTreeOperatorTreeProcessingTimer), intent(inout)               :: self
    double precision                                       , allocatable  , dimension(:) :: timesConstructTemporary  , timesEvolveTemporary, &
         &                                                                                  massesTemporary
    integer         (kind_int8                            ), allocatable  , dimension(:) :: treeIDsTemporary
    integer         (c_size_t                             ), allocatable  , dimension(:) :: memoryUsagesPeakTemporary, countsNodesTemporary
    real                                                                                 :: time
    
    ! Record the final CPU time.
    !$ if (omp_in_parallel()) then
    !$ self%timePostEvolution=OMP_Get_WTime()
    !$ else
    call CPU_Time(time)
    self%timePostEvolution=dble(time)
    !$ end if
    ! Check that the tree was actually processed. If no pre-evolution time was recorded then no tree existed to be
    ! processed. In that case we ignore the results.
    if (self%timePreEvolution > 0.0d0) then
       ! Ensure that record arrays are sufficiently sized.
       if (.not.allocated(self%masses)) then
          allocate(self%masses          (countIncrement))
          allocate(self%timesConstruct  (countIncrement))
          allocate(self%timesEvolve     (countIncrement))
          allocate(self%treeIDs         (countIncrement))
          allocate(self%memoryUsagesPeak(countIncrement))
          allocate(self%countsNodes     (countIncrement))
       else if (self%countTrees >= size(self%masses)) then
          call Move_Alloc   (self%masses          ,massesTemporary          )
          call Move_Alloc   (self%timesConstruct  ,timesConstructTemporary  )
          call Move_Alloc   (self%timesEvolve     ,timesEvolveTemporary     )
          call Move_Alloc   (self%treeIDs         ,treeIDsTemporary         )
          call Move_Alloc   (self%memoryUsagesPeak,memoryUsagesPeakTemporary)
          call Move_Alloc   (self%countsNodes     ,countsNodesTemporary     )
          allocate(self%masses          (size(massesTemporary)+countIncrement))
          allocate(self%timesConstruct  (size(massesTemporary)+countIncrement))
          allocate(self%timesEvolve     (size(massesTemporary)+countIncrement))
          allocate(self%treeIDs         (size(massesTemporary)+countIncrement))
          allocate(self%memoryUsagesPeak(size(massesTemporary)+countIncrement))
          allocate(self%countsNodes     (size(massesTemporary)+countIncrement))
          self%masses          (1:size(massesTemporary))=massesTemporary
          self%timesConstruct  (1:size(massesTemporary))=timesConstructTemporary
          self%timesEvolve     (1:size(massesTemporary))=timesEvolveTemporary
          self%treeIDs         (1:size(massesTemporary))=treeIDsTemporary
          self%memoryUsagesPeak(1:size(massesTemporary))=memoryUsagesPeakTemporary
          self%countsNodes     (1:size(massesTemporary))=countsNodesTemporary
          deallocate(massesTemporary          )
          deallocate(timesConstructTemporary  )
          deallocate(timesEvolveTemporary     )
          deallocate(treeIDsTemporary         )
          deallocate(memoryUsagesPeakTemporary)
          deallocate(countsNodesTemporary     )
       end if
       ! Store the timing data.
       self%countTrees                       =self%countTrees       +1
       self%masses          (self%countTrees)=self%mass
       self%timesConstruct  (self%countTrees)=self%timePreEvolution -self%timePreConstruction
       self%timesEvolve     (self%countTrees)=self%timePostEvolution-self%timePreEvolution
       self%treeIDs         (self%countTrees)=self%treeID
       self%memoryUsagesPeak(self%countTrees)=self%memoryUsagePeak
       self%countsNodes     (self%countTrees)=self%countNodes
    end if
    return
  end subroutine treeProcessingTimerOperatePostEvolution
  
  subroutine treeProcessingTimerPostEvolve(self,node)
    !!{RST
    Record memory usage.
    !!}
    use :: Galacticus_Nodes   , only : mergerTree              , treeNode
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
    implicit none
    class  (*                       ), intent(inout)          :: self
    type   (treeNode                ), intent(inout), target  :: node
    type   (mergerTree              )               , pointer :: tree       , treeWork
    type   (treeNode                )               , pointer :: nodeWork
    integer(c_size_t                )                         :: memoryUsage
    type   (mergerTreeWalkerAllNodes)                         :: treeWalker

    select type (self)
       class is (mergerTreeOperatorTreeProcessingTimer)
       tree        => node%hostTree%firstTree
       treeWork    => tree
       treeWalker  =  mergerTreeWalkerAllNodes(tree,spanForest=.true.)
       memoryUsage =  0_c_size_t
       do while (associated(treeWork))
          memoryUsage =  memoryUsage+treeWork%sizeof  ()
          treeWork    =>             treeWork%nextTree
       end do
       do while (treeWalker%next(nodeWork))
          memoryUsage=memoryUsage+nodeWork%sizeOf()
       end do
       self%memoryUsagePeak=max(self%memoryUsagePeak,memoryUsage)
    end select
    return
  end subroutine treeProcessingTimerPostEvolve

  subroutine treeProcessingTimerFinalize(self)
    !!{RST
    Outputs collected meta-data on tree processing times.
    !!}
    use :: Output_HDF5                     , only : outputFile
    use :: HDF5_Access                     , only : hdf5Access
    use :: IO_HDF5                         , only : hdf5Object
    use :: HDF5                            , only : hsize_t
    use :: Numerical_Constants_Astronomical, only : massSolar
    use :: Units_MetaData                  , only : unitType
    implicit none
    class           (mergerTreeOperatorTreeProcessingTimer), intent(inout) :: self
    integer         (hsize_t                              ), parameter      :: chunkSize      =100_hsize_t
    type            (hdf5Object                           )                 :: metaDataDataset            , metaDataGroup, &
         &                                                                     timingDataGroup
    logical                                                                 :: preExists

    ! Each (OpenMP-threadprivate) copy of this operator writes the timing data for the trees that it processed, appending to the
    ! shared datasets in the output file. The quadratic log-log fit of the processing time is deliberately *not* performed here, as
    ! any single copy holds only the trees it processed; the fit is performed once, over the complete dataset, below.
    if (self%countTrees > 0) then
       !$ call hdf5Access%set()
       metaDataGroup  =outputFile     %openGroup ('metaData'  ,'Galacticus meta data.'    )
       timingDataGroup=metaDataGroup  %openGroup ('treeTiming','Meta-data on tree timing.')
       preExists      =timingDataGroup%hasDataset('treeID'                                )
       call timingDataGroup       %writeDataset  (self%treeIDs         (1:self%countTrees),"treeID"         ,"Tree ID"                    ,chunkSize=chunkSize                                ,appendTo=.true.)
       call timingDataGroup       %writeDataset  (self%masses          (1:self%countTrees),"treeMass"       ,"Tree mass [M⊙]"             ,chunkSize=chunkSize,datasetReturned=metaDataDataset,appendTo=.true.)
       if (.not.preExists) &
            & call metaDataDataset%writeAttribute(unitType(massSolar,"Solar masses","solMass"),"units")
       call timingDataGroup       %writeDataset  (self%timesConstruct  (1:self%countTrees),"timeConstruct"  ,"Tree construction time [s]" ,chunkSize=chunkSize,datasetReturned=metaDataDataset,appendTo=.true.)
       if (.not.preExists) &
            & call metaDataDataset%writeAttribute(unitType(1.0d0                             ),"units")
       call timingDataGroup       %writeDataset  (self%timesEvolve     (1:self%countTrees),"timeEvolve"     ,"Tree evolution time [s]"    ,chunkSize=chunkSize,datasetReturned=metaDataDataset,appendTo=.true.)
       if (.not.preExists) &
            & call metaDataDataset%writeAttribute(unitType(1.0d0                             ),"units")
       call timingDataGroup       %writeDataset  (self%countsNodes     (1:self%countTrees),"countNodes"     ,"Number of nodes in the tree",chunkSize=chunkSize                                ,appendTo=.true.)
       if (self%collectMemoryUsageData) &
            & call timingDataGroup%writeDataset  (self%memoryUsagesPeak(1:self%countTrees),"memoryUsagePeak","Peak memory usage [bytes]"  ,chunkSize=chunkSize                                ,appendTo=.true.)
       !$ call hdf5Access%unset()
    end if
    ! Ensure all threadprivate copies have finished appending their raw timing data before the fit is computed over the complete
    ! dataset, then have a single thread perform and write the fit. (A barrier and single-threaded section are safe here as this
    ! finalize method is called by every thread of the team; both degrade to no-ops when executed serially.)
    !$omp barrier
    !$omp masked
    call treeProcessingTimerFitWrite()
    !$omp end masked
    !$omp barrier
    return
  end subroutine treeProcessingTimerFinalize

  subroutine treeProcessingTimerFitWrite()
    !!{RST
    Fit quadratic log-log relations of the total (construction plus evolution) processing time versus tree mass and versus node
    count, and write the fit coefficients, residual scatter, and range of validity to the ``metaData/treeTiming`` group. These
    provide a ready-to-use cost model that a subsequent run of the same model can consume directly (see the ``file`` tree
    processing time class). The complete timing dataset is read back from the output file so that the fit reflects all trees
    processed (across all threads), not just those processed by any single thread. The fit is written only if it does not already
    exist (e.g. it is not recomputed on a restart).
    !!}
    use :: Output_HDF5, only : outputFile
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Object
    implicit none
    type            (hdf5Object)                              :: metaDataGroup          , timingDataGroup
    double precision            , allocatable, dimension(:)   :: massRead               , timeConstructRead     , &
         &                                                       timeEvolveRead         , timesProcess
    integer         (c_size_t  ), allocatable, dimension(:)   :: countNodesRead
    double precision                         , dimension(0:2) :: coefficientsMass       , coefficientsCountNodes
    double precision                         , dimension(  2) :: rangeMass              , rangeCountNodes
    double precision                                          :: residualMass           , residualCountNodes
    logical                                                   :: successMass            , successCountNodes

    !$ call hdf5Access%set()
    if (outputFile%hasGroup('metaData')) then
       metaDataGroup=outputFile%openGroup('metaData')
       if (metaDataGroup%hasGroup('treeTiming')) then
          timingDataGroup=metaDataGroup%openGroup('treeTiming')
          if (timingDataGroup%hasDataset('treeMass') .and. .not.timingDataGroup%hasDataset('fitCoefficientMass')) then
             call timingDataGroup%readDataset('treeMass'     ,massRead         )
             call timingDataGroup%readDataset('timeConstruct',timeConstructRead)
             call timingDataGroup%readDataset('timeEvolve'   ,timeEvolveRead   )
             call timingDataGroup%readDataset('countNodes'   ,countNodesRead   )
             allocate(timesProcess(size(massRead)))
             timesProcess=timeConstructRead+timeEvolveRead
             call treeProcessingTimerFit(massRead            ,timesProcess,coefficientsMass      ,residualMass      ,rangeMass      ,successMass      )
             call treeProcessingTimerFit(dble(countNodesRead),timesProcess,coefficientsCountNodes,residualCountNodes,rangeCountNodes,successCountNodes)
             if (successMass) then
                call timingDataGroup%writeDataset(coefficientsMass,"fitCoefficientMass","Coefficients Cᵢ of the fit log₁₀(τ/s) = Σ Cᵢ (log₁₀[M/M⊙])ⁱ")
                call timingDataGroup%writeDataset([residualMass]  ,"fitResidualMass"   ,"RMS residual of the mass-based fit [dex]"                    )
                call timingDataGroup%writeDataset(rangeMass       ,"fitRangeMass"      ,"Range of tree masses used in the fit [M⊙]"                   )
             end if
             if (successCountNodes) then
                call timingDataGroup%writeDataset(coefficientsCountNodes,"fitCoefficientCountNodes","Coefficients Cᵢ of the fit log₁₀(τ/s) = Σ Cᵢ (log₁₀[N])ⁱ")
                call timingDataGroup%writeDataset([residualCountNodes]  ,"fitResidualCountNodes"   ,"RMS residual of the node-count-based fit [dex]"          )
                call timingDataGroup%writeDataset(rangeCountNodes       ,"fitRangeCountNodes"      ,"Range of node counts used in the fit"                    )
             end if
          end if
       end if
    end if
    !$ call hdf5Access%unset()
    return
  end subroutine treeProcessingTimerFitWrite

  subroutine treeProcessingTimerFit(predictor,times,coefficients,residual,rangePredictor,success)
    !!{RST
    Perform a least-squares fit of :math:`\log_{10}(\tau)` versus :math:`\log_{10}(x)`, where :math:`\tau` is the processing time and
    :math:`x` the predictor variable (tree mass or node count), returning three polynomial coefficients (padded with zeros if a
    lower-order fit is used), the RMS residual scatter (in dex), and the range of the predictor over the data used. The polynomial
    degree is reduced automatically (to linear, or constant) if the data do not support a quadratic fit.
    !!}
    implicit none
    double precision, intent(in   ), dimension(     :) :: predictor        , times
    double precision, intent(  out), dimension(0:2   ) :: coefficients
    double precision, intent(  out)                    :: residual
    double precision, intent(  out), dimension(  2   ) :: rangePredictor
    logical         , intent(  out)                    :: success
    double precision               , allocatable, dimension(:) :: x               , y
    double precision                              , dimension(0:2,0:2) :: matrixNormal
    double precision                              , dimension(0:2    ) :: vectorNormal
    integer                                            :: i                 , j                , &
         &                                                k                 , countValid        , &
         &                                                degree
    double precision                                   :: model             , sumSquares

    coefficients  =0.0d0
    residual      =0.0d0
    rangePredictor=0.0d0
    success       =.false.
    ! Collect data for which both the predictor and the time are strictly positive (so that logarithms are defined).
    countValid=count(predictor > 0.0d0 .and. times > 0.0d0)
    if (countValid < 1) return
    allocate(x(countValid))
    allocate(y(countValid))
    k=0
    do i=1,size(predictor)
       if (predictor(i) > 0.0d0 .and. times(i) > 0.0d0) then
          k   =k+1
          x(k)=log10(predictor(i))
          y(k)=log10(times    (i))
       end if
    end do
    rangePredictor=[minval(10.0d0**x),maxval(10.0d0**x)]
    ! Choose the highest polynomial degree supported by the number of data points (at most quadratic), then reduce further if the
    ! normal-equations matrix proves to be singular (e.g. all predictor values are identical).
    degree=min(2,countValid-1)
    do while (degree >= 0)
       ! Build the normal equations for a polynomial fit of the chosen degree.
       matrixNormal=0.0d0
       vectorNormal=0.0d0
       do i=1,countValid
          do j=0,degree
             do k=0,degree
                matrixNormal(j,k)=matrixNormal(j,k)+x(i)**(j+k)
             end do
             vectorNormal(j)=vectorNormal(j)+x(i)**j*y(i)
          end do
       end do
       if (treeProcessingTimerSolve(matrixNormal(0:degree,0:degree),vectorNormal(0:degree),coefficients(0:degree),degree+1)) then
          success=.true.
          exit
       end if
       coefficients=0.0d0
       degree      =degree-1
    end do
    if (.not.success) return
    ! Compute the RMS residual scatter of the fit (in dex).
    sumSquares=0.0d0
    do i=1,countValid
       model=0.0d0
       do j=0,degree
          model=model+coefficients(j)*x(i)**j
       end do
       sumSquares=sumSquares+(y(i)-model)**2
    end do
    if (countValid > degree+1) then
       residual=sqrt(sumSquares/dble(countValid-degree-1))
    else
       residual=0.0d0
    end if
    return
  end subroutine treeProcessingTimerFit

  logical function treeProcessingTimerSolve(matrix_,vector_,solution,n) result(success)
    !!{RST
    Solve the small linear system ``matrix_ x = vector_`` by Gaussian elimination with partial pivoting, returning
    ``.false.`` if the matrix is (numerically) singular.
    !!}
    implicit none
    integer         , intent(in   )                      :: n
    double precision, intent(in   ), dimension(0:n-1,0:n-1) :: matrix_
    double precision, intent(in   ), dimension(0:n-1      ) :: vector_
    double precision, intent(  out), dimension(0:n-1      ) :: solution
    double precision               , dimension(0:n-1,0:n-1) :: a
    double precision               , dimension(0:n-1      ) :: b
    integer                                              :: i          , j        , &
         &                                                  k          , pivotRow
    double precision                                     :: pivotValue , factor   , &
         &                                                  sumOff

    success=.false.
    a      =matrix_
    b      =vector_
    do k=0,n-1
       ! Partial pivoting: find the row with the largest magnitude entry in the current column.
       pivotRow  =k
       pivotValue=abs(a(k,k))
       do i=k+1,n-1
          if (abs(a(i,k)) > pivotValue) then
             pivotRow  =i
             pivotValue=abs(a(i,k))
          end if
       end do
       if (pivotValue <= 0.0d0) return
       if (pivotRow /= k) then
          a([k,pivotRow],:)=a([pivotRow,k],:)
          b([k,pivotRow]  )=b([pivotRow,k]  )
       end if
       do i=k+1,n-1
          factor=a(i,k)/a(k,k)
          a(i,:)=a(i,:)-factor*a(k,:)
          b(i  )=b(i  )-factor*b(k  )
       end do
    end do
    ! Back-substitution.
    do i=n-1,0,-1
       sumOff=0.0d0
       do j=i+1,n-1
          sumOff=sumOff+a(i,j)*solution(j)
       end do
       solution(i)=(b(i)-sumOff)/a(i,i)
    end do
    success=.true.
    return
  end function treeProcessingTimerSolve
