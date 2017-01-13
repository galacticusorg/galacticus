!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which performs analysis to compute halo mass functions.

module Galacticus_Output_Analyses_Halo_Mass_Function
  !% Performs analysis to compute halo mass functions.
  use, intrinsic :: ISO_C_Binding
  use ISO_Varying_String
  use Galacticus_Nodes
  implicit none
  private
  public :: Galacticus_Output_Analysis_Halo_Mass_Functions, Galacticus_Output_Analysis_Halo_Mass_Functions_Output

  ! Record of module initialization.
  logical :: moduleInitialized=.false.

  ! Record of whether this analysis is active.
  logical :: analysisActive
  
  ! Type to store mass functions.
  type :: massFunction
     ! Label.
     type            (varying_string)                                :: label
     ! Arrays for the masses and mass function.
     double precision                , allocatable, dimension(:    ) :: masses                  , massesLogarithmic           , &
          &                                                             massesLogarithmicMinimum, massesLogarithmicMaximum    , &
          &                                                             massFunction
     integer         (c_size_t      )                                :: outputIndex
     logical                                                         :: alwaysIsolatedHalosOnly
     ! Arrays for accumulation of of main branch halos
     double precision                , allocatable, dimension(:,:  ) :: mainBranchHaloWeights   , mainBranchHaloWeightsSquared
     ! Array for the covariance matrix.
     double precision                , allocatable, dimension(:,:  ) :: massFunctionCovariance
  end type massFunction

  ! Mass functions.
  type(massFunction), allocatable, dimension(:) :: massFunctions

  ! Type for storing temporary mass functions during cumulation.
  type :: massFunctionWork
     double precision, allocatable, dimension(:  ) :: massFunction
     double precision, allocatable, dimension(:,:) :: covariance
  end type massFunctionWork

  ! Work array.
  type(massFunctionWork), allocatable, dimension(:) :: thisHalo
  !$omp threadprivate(thisHalo)

  ! Options controlling mass function construction.
  integer                     :: analysisHaloMassFunctionsMassBinsCount                           , analysisHaloMassFunctionsMassBinsPerDecade
  double precision            :: analysisHaloMassFunctionsMassMinimum                             , analysisHaloMassFunctionsMassMaximum           , &
       &                         analysisHaloMassFunctionsMassIntervalLogarithmicInverse          , analysisHaloMassFunctionsMassMinimumLogarithmic
  
  ! Options controlling binning in halo mass.
  integer                     :: analysisHaloMassFunctionCovarianceModel
  integer         , parameter :: analysisHaloMassFunctionCovarianceModelPoisson             =1
  integer         , parameter :: analysisHaloMassFunctionCovarianceModelBinomial            =2
  integer                     :: analysisHaloMassFunctionsHaloMassBinsCount                   , analysisHaloMassFunctionsHaloMassBinsPerDecade
  double precision            :: analysisHaloMassFunctionsHaloMassMinimum                     , analysisHaloMassFunctionsHaloMassMaximum           , &
       &                         analysisHaloMassFunctionsHaloMassIntervalLogarithmicInverse  , analysisHaloMassFunctionsHaloMassMinimumLogarithmic

  ! Options controlling covariance matrix construction.
  double precision            :: analysisHaloMassFunctionsCorrelationTruncateLevel

contains

  !# <mergerTreeAnalysisTask>
  !#  <unitName>Galacticus_Output_Analysis_Halo_Mass_Functions</unitName>
  !# </mergerTreeAnalysisTask>
  subroutine Galacticus_Output_Analysis_Halo_Mass_Functions(thisTree,thisNode,nodeStatus,iOutput,mergerTreeAnalyses)
    !% Construct halo mass functions.
    use, intrinsic :: ISO_C_Binding
    use Galacticus_Nodes
    use Memory_Management
    use Input_Parameters
    use Galacticus_Output_Times
    use Galacticus_Error
    use Cosmology_Functions
    use Numerical_Comparison
    use String_Handling
    use Galacticus_Output_Merger_Tree_Data
    use Vectors
    use Numerical_Ranges
    use Numerical_Comparison
    implicit none
    type            (mergerTree                    ), intent(inout)                 :: thisTree
    type            (treeNode                      ), intent(inout), pointer        :: thisNode
    integer                                         , intent(in   )                 :: nodeStatus
    integer         (c_size_t                      ), intent(in   )                 :: iOutput
    type            (varying_string                ), intent(in   ), dimension(:  ) :: mergerTreeAnalyses
    class           (nodeComponentBasic            )               , pointer        :: thisBasic
    class           (nodeComponentMergingStatistics)               , pointer        :: mergingStatistics
    class           (cosmologyFunctionsClass       )               , pointer        :: cosmologyFunctions_
    integer                                                                         :: currentAnalysis,activeAnalysisCount,haloMassBin,massBin,i,k
    double precision                                                                :: massLogarithmic,redshift,outputTime
    type            (varying_string                )                                :: analysisHaloMassFunctionCovarianceModelText
    character       (len=32                        )                                :: redshiftLabel
    
    ! Initialize the module if necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Galacticus_Output_Analysis_Halo_Mass_Functions_Initialize)
       if (.not.moduleInitialized) then
          !@ <inputParameter>
          !@   <name>analysisHaloMassFunctionCovarianceModel</name>
          !@   <defaultValue>binomial</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The model to use when constructing the mass function covariance matrix for main branch halos.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisHaloMassFunctionCovarianceModel',analysisHaloMassFunctionCovarianceModelText,defaultValue='Poisson')
          select case (char(analysisHaloMassFunctionCovarianceModelText))
          case ( 'Poisson'  )
             analysisHaloMassFunctionCovarianceModel=analysisHaloMassFunctionCovarianceModelPoisson
          case ( 'binomial' )
             analysisHaloMassFunctionCovarianceModel=analysisHaloMassFunctionCovarianceModelBinomial
          case default
             call Galacticus_Error_Report('Galacticus_Output_Analysis_Halo_Mass_Functions','unrecognized value for "analysisHaloMassFunctionCovarianceModel" - allowed values are "Poisson", and "binomial"')
          end select
          !@ <inputParameter>
          !@   <name>analysisHaloMassFunctionsHaloMassBinsPerDecade</name>
          !@   <defaultValue>10</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The number of bins per decade of halo mass to use when constructing the mass function covariance matrix for main branch halos.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisHaloMassFunctionsHaloMassBinsPerDecade',analysisHaloMassFunctionsHaloMassBinsPerDecade,defaultValue=10)
          !@ <inputParameter>
          !@   <name>analysisHaloMassFunctionsHaloMassMinimum</name>
          !@   <defaultValue>$10^8M_\odot$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The minimum halo mass to consider when constructing the mass function covariance matrix for main branch halos.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisHaloMassFunctionsHaloMassMinimum',analysisHaloMassFunctionsHaloMassMinimum,defaultValue=1.0d8)
          !@ <inputParameter>
          !@   <name>analysisHaloMassFunctionsHaloMassMaximum</name>
          !@   <defaultValue>$10^{16}M_\odot$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The maximum halo mass to consider when constructing the mass function covariance matrix for main branch halos.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisHaloMassFunctionsHaloMassMaximum',analysisHaloMassFunctionsHaloMassMaximum,defaultValue=1.0d16)
          analysisHaloMassFunctionsHaloMassMinimumLogarithmic=log10(analysisHaloMassFunctionsHaloMassMinimum)
          analysisHaloMassFunctionsHaloMassBinsCount=int(log10(analysisHaloMassFunctionsHaloMassMaximum/analysisHaloMassFunctionsHaloMassMinimum)*dble(analysisHaloMassFunctionsHaloMassBinsPerDecade)+0.5d0)
          analysisHaloMassFunctionsHaloMassIntervalLogarithmicInverse=dble(analysisHaloMassFunctionsHaloMassBinsCount)/log10(analysisHaloMassFunctionsHaloMassMaximum/analysisHaloMassFunctionsHaloMassMinimum)
          !@ <inputParameter>
          !@   <name>analysisHaloMassFunctionsMassBinsPerDecade</name>
          !@   <defaultValue>10</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The number of bins per decade of halo mass to use when constructing the halo mass function.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisHaloMassFunctionsMassBinsPerDecade',analysisHaloMassFunctionsMassBinsPerDecade,defaultValue=10)
          !@ <inputParameter>
          !@   <name>analysisHaloMassFunctionsMassMinimum</name>
          !@   <defaultValue>$10^8M_\odot$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The minimum halo mass to consider when constructing the mass function covariance matrix for main branch halos.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisHaloMassFunctionsMassMinimum',analysisHaloMassFunctionsMassMinimum,defaultValue=1.0d8)
          !@ <inputParameter>
          !@   <name>analysisHaloMassFunctionsMassMaximum</name>
          !@   <defaultValue>$10^{16}M_\odot$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The maximum halo mass to consider when constructing the mass function covariance matrix for main branch halos.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisHaloMassFunctionsMassMaximum',analysisHaloMassFunctionsMassMaximum,defaultValue=1.0d16)
          analysisHaloMassFunctionsMassMinimumLogarithmic=log10(analysisHaloMassFunctionsMassMinimum)
          analysisHaloMassFunctionsMassBinsCount=int(log10(analysisHaloMassFunctionsMassMaximum/analysisHaloMassFunctionsMassMinimum)*dble(analysisHaloMassFunctionsMassBinsPerDecade)+0.5d0)
          analysisHaloMassFunctionsMassIntervalLogarithmicInverse=dble(analysisHaloMassFunctionsMassBinsCount)/log10(analysisHaloMassFunctionsMassMaximum/analysisHaloMassFunctionsMassMinimum)
          !@ <inputParameter>
          !@   <name>analysisHaloMassFunctionsCorrelationTruncateLevel</name>
          !@   <defaultValue>0.0</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The correlation below which off-diagonal elements of the covariance matrix are truncated to zero.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisHaloMassFunctionsCorrelationTruncateLevel',analysisHaloMassFunctionsCorrelationTruncateLevel,defaultValue=0.0d0)
          ! Determine how many supported mass functions are requested.
          activeAnalysisCount=0
          do i=1,size(mergerTreeAnalyses)
             if     (                                                                    &
                  &   extract(mergerTreeAnalyses(i),1,17) ==         "haloMassFunctionZ" &
                  &  .or.                                                                &
                  &   extract(mergerTreeAnalyses(i),1,25) == "isolatedHaloMassFunctionZ" &
                  & ) activeAnalysisCount=activeAnalysisCount+1
          end do
          ! Allocate mass function arrays and populate with required data.
          if (activeAnalysisCount <= 0) then
             analysisActive=.false.
          else
             analysisActive=.true.
             cosmologyFunctions_ => cosmologyFunctions()
             ! Initialize analyses.
             currentAnalysis=0
             allocate(massFunctions(activeAnalysisCount))
             do i=1,size(mergerTreeAnalyses)
                if     (                                                                    &
                     &   extract(mergerTreeAnalyses(i),1,17) ==         "haloMassFunctionZ" &
                     &  .or.                                                                &
                     &   extract(mergerTreeAnalyses(i),1,25) == "isolatedHaloMassFunctionZ" &
                     & ) then
                   currentAnalysis=currentAnalysis+1
                   massFunctions(currentAnalysis)%label=mergerTreeAnalyses(i)
                   ! Allocate arrays.
                   call allocateArray(massFunctions(currentAnalysis)%masses                      ,[analysisHaloMassFunctionsMassBinsCount                                           ])
                   call allocateArray(massFunctions(currentAnalysis)%massesLogarithmic           ,[analysisHaloMassFunctionsMassBinsCount                                           ])
                   call allocateArray(massFunctions(currentAnalysis)%massesLogarithmicMinimum    ,[analysisHaloMassFunctionsMassBinsCount                                           ])
                   call allocateArray(massFunctions(currentAnalysis)%massesLogarithmicMaximum    ,[analysisHaloMassFunctionsMassBinsCount                                           ])
                   call allocateArray(massFunctions(currentAnalysis)%massFunction                ,[analysisHaloMassFunctionsMassBinsCount                                           ])
                   call allocateArray(massFunctions(currentAnalysis)%massFunctionCovariance      ,[analysisHaloMassFunctionsMassBinsCount,analysisHaloMassFunctionsMassBinsCount    ])
                   call allocateArray(massFunctions(currentAnalysis)%mainBranchHaloWeights       ,[analysisHaloMassFunctionsMassBinsCount,analysisHaloMassFunctionsHaloMassBinsCount])
                   call allocateArray(massFunctions(currentAnalysis)%mainBranchHaloWeightsSquared,[analysisHaloMassFunctionsMassBinsCount,analysisHaloMassFunctionsHaloMassBinsCount])
                   ! Initialize arrays.                                    
                   massFunctions(currentAnalysis)%masses                      =Make_Range(                                                                                                               &
                        &                                                                 analysisHaloMassFunctionsMassMinimum  *10.0**(+0.5d0/analysisHaloMassFunctionsMassIntervalLogarithmicInverse), &
                        &                                                                 analysisHaloMassFunctionsMassMaximum  *10.0**(-0.5d0/analysisHaloMassFunctionsMassIntervalLogarithmicInverse), &
                        &                                                                 analysisHaloMassFunctionsMassBinsCount                                                                       , &
                        &                                                                 rangeTypeLogarithmic                                                                                           &
                        &                                                                )
                   massFunctions(currentAnalysis)%massesLogarithmic           =log10(massFunctions(currentAnalysis)%masses)
                   massFunctions(currentAnalysis)%massFunction                =0.0d0
                   massFunctions(currentAnalysis)%massFunctionCovariance      =0.0d0
                   massFunctions(currentAnalysis)%mainBranchHaloWeights       =0.0d0
                   massFunctions(currentAnalysis)%mainBranchHaloWeightsSquared=0.0d0
                   do k=1,analysisHaloMassFunctionsMassBinsCount
                      if (k ==                                      1) then
                         massFunctions(currentAnalysis)%massesLogarithmicMinimum(k)=massFunctions(currentAnalysis)%massesLogarithmic(k)-0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k+1)-massFunctions(currentAnalysis)%massesLogarithmic(k  ))
                      else
                         massFunctions(currentAnalysis)%massesLogarithmicMinimum(k)=                                                   +0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k-1)+massFunctions(currentAnalysis)%massesLogarithmic(k  ))
                      end if
                      if (k == analysisHaloMassFunctionsMassBinsCount) then
                         massFunctions(currentAnalysis)%massesLogarithmicMaximum(k)=massFunctions(currentAnalysis)%massesLogarithmic(k)+0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k  )-massFunctions(currentAnalysis)%massesLogarithmic(k-1))
                      else
                         massFunctions(currentAnalysis)%massesLogarithmicMaximum(k)=                                                   +0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k+1)+massFunctions(currentAnalysis)%massesLogarithmic(k  ))
                      end if
                   end do
                   if      (extract(mergerTreeAnalyses(i),1,17) ==         "haloMassFunctionZ") then
                      redshiftLabel=extract(mergerTreeAnalyses(i),18,len(mergerTreeAnalyses(i)))
                      massFunctions(currentAnalysis)%alwaysIsolatedHalosOnly=.false.
                   else if (extract(mergerTreeAnalyses(i),1,25) == "isolatedHaloMassFunctionZ") then
                      if (.not.defaultMergingStatisticsComponent%nodeHierarchyLevelMaximumIsGettable())                                                                              &
                           & call Galacticus_Error_Report                                                                                                                            &
                           &      (                                                                                                                                                  &
                           &       'Galacticus_Output_Analysis_Halo_Mass_Functions'                                                                                               ,  &
                           &       'mass functions of always isolated halos require a merging statistics component that provides a gettable "nodeHierarchyLevelMaximum" property.'// &
                           &       Galacticus_Component_List(                                                                                                                        &
                           &                                 'mergingStatistics'                                                                                                   , &
                           &                                  defaultMergingStatisticsComponent%nodeHierarchyLevelMaximumAttributeMatch(requireGettable=.true.)                      &
                           &                                )                                                                                                                        &
                           &      )    
                      redshiftLabel=extract(mergerTreeAnalyses(i),26,len(mergerTreeAnalyses(i)))
                      massFunctions(currentAnalysis)%alwaysIsolatedHalosOnly=.true.
                   else
                      call Galacticus_Error_Report('Galacticus_Output_Analysis_Halo_Mass_Functions','unrecognized mass function type')
                   end if
                   ! Find the index of the output corresponding to the requested redshift.
                   read (redshiftLabel,*) redshift
                   outputTime=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift))
                   massFunctions(currentAnalysis)%outputIndex=Galacticus_Output_Time_Index(outputTime,findClosest=.true.)
                   if (Values_Differ(Galacticus_Output_Time(massFunctions(currentAnalysis)%outputIndex),outputTime,relTol=1.0d-3)) &
                        & call Galacticus_Error_Report('Galacticus_Output_Analysis_Halo_Mass_Functions','no output available for requested analysis')
                end if
             end do
          end if
          ! Record that module is initialized.
          moduleInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Analysis_Halo_Mass_Functions_Initialize)
    end if
    ! Return if this analysis is not active.
    if (.not.analysisActive) return
    ! Return if this is a tree finalization.
    if (nodeStatus == nodeStatusFinal) return
    ! Allocate work arrays.
    if (.not.allocated(thisHalo)) allocate(thisHalo(size(massFunctions)))
    ! Iterate over active analyses.
    do i=1,size(massFunctions)
       ! Cycle if this mass function does not accumulate from this output number.
       if (iOutput /= massFunctions(i)%outputIndex) cycle
       ! Cycle if this is a subhalo.
       if (thisNode%isSatellite()) cycle
       ! Exclude halos which were not always isolated, if required.
       if (massFunctions(i)%alwaysIsolatedHalosOnly) then       
          mergingStatistics => thisNode%mergingStatistics()
          if (mergingStatistics%nodeHierarchyLevelMaximum() > 0) cycle
       end if
       ! Allocate workspace.
       if (.not.allocated(thisHalo(i)%massFunction)) then
          call allocateArray(thisHalo(i)%massFunction,[analysisHaloMassFunctionsMassBinsCount])
          call allocateArray(thisHalo(i)%covariance  ,[                                         &
               &                                     analysisHaloMassFunctionsMassBinsCount,  &
               &                                     analysisHaloMassFunctionsMassBinsCount   &
               &                                    ]                                         &
               &          )
       end if
       ! Get the halo mass.
       thisBasic      =>       thisNode%basic()
       massLogarithmic=  log10(thisBasic%mass())
       ! Determine which bin this halo contributes to.
       massBin=floor((massLogarithmic-analysisHaloMassFunctionsMassMinimumLogarithmic)*analysisHaloMassFunctionsMassIntervalLogarithmicInverse)+1
       if (massBin >= 1 .and. massBin <= analysisHaloMassFunctionsMassBinsCount) then          
          !$omp critical (Galacticus_Output_Analysis_Halo_Mass_Functions_Accumulate)
          massFunctions        (i)%massFunction(massBin) &
               & =massFunctions(i)%massFunction(massBin) &
               & +thisTree%volumeWeight
          !$omp end critical (Galacticus_Output_Analysis_Halo_Mass_Functions_Accumulate)
          ! Treat main branch and other galaxies differently.
          if (thisNode%isOnMainBranch().and.analysisHaloMassFunctionCovarianceModel == analysisHaloMassFunctionCovarianceModelBinomial) then
             ! Find the bin to which this halo mass belongs.
             haloMassBin=floor((massLogarithmic-analysisHaloMassFunctionsHaloMassMinimumLogarithmic)*analysisHaloMassFunctionsHaloMassIntervalLogarithmicInverse)+1
             ! Accumulate weights to halo mass arrays.
             if (haloMassBin >= 1 .and. haloMassBin <= analysisHaloMassFunctionsHaloMassBinsCount) then
                !$omp critical (Galacticus_Output_Analysis_Halo_Mass_Functions_Accumulate)
                massFunctions        (i)%mainBranchHaloWeights       (massBin,haloMassBin)= &
                     &  massFunctions(i)%mainBranchHaloWeights       (massBin,haloMassBin)  &
                     &  +thisTree%volumeWeight
                massFunctions        (i)%mainBranchHaloWeightsSquared(massBin,haloMassBin)= &
                     &  massFunctions(i)%mainBranchHaloWeightsSquared(massBin,haloMassBin)  &
                     &  +thisTree%volumeWeight
                !$omp end critical (Galacticus_Output_Analysis_Halo_Mass_Functions_Accumulate)
             end if
          else
             ! Accumulate covariance.
             !$omp critical (Galacticus_Output_Analysis_Halo_Mass_Functions_Accumulate)
             massFunctions        (i)%massFunctionCovariance(massBin,massBin) &
                  & =massFunctions(i)%massFunctionCovariance(massBin,massBin) &
                  & +thisTree%volumeWeight**2
             !$omp end critical (Galacticus_Output_Analysis_Halo_Mass_Functions_Accumulate)
          end if
       end if
    end do
    return    
  end subroutine Galacticus_Output_Analysis_Halo_Mass_Functions
  
  !# <hdfPreCloseTask>
  !#  <unitName>Galacticus_Output_Analysis_Halo_Mass_Functions_Output</unitName>
  !# </hdfPreCloseTask>
  subroutine Galacticus_Output_Analysis_Halo_Mass_Functions_Output
    !% Outputs halo mass functions to file.
    use Galacticus_HDF5
    use Numerical_Constants_Astronomical
    implicit none
    integer                      :: i,j,k,m
    type            (hdf5Object) :: analysisGroup,massFunctionGroup,thisDataset
    double precision             :: haloWeightBinTotal

    ! Return immediately if this analysis is not active.
    if (.not.analysisActive) return
    ! Iterate over mass functions.
    do k=1,size(massFunctions)
       ! Add the contribution from main branch galaxies to the covariance matrix.
       if (analysisHaloMassFunctionCovarianceModel == analysisHaloMassFunctionCovarianceModelBinomial) then
          do m=1,analysisHaloMassFunctionsHaloMassBinsCount
             haloWeightBinTotal=sum(massFunctions(k)%mainBranchHaloWeights(:,m))
             if ( haloWeightBinTotal > 0.0 ) then
                do i=1,analysisHaloMassFunctionsMassBinsCount
                   massFunctions               (k)%massFunctionCovariance      (i,i)=                    &
                        &         massFunctions(k)%massFunctionCovariance      (i,i)                     &
                        & +(1.0d0-massFunctions(k)%mainBranchHaloWeights       (i,m)/haloWeightBinTotal) &
                        & *       massFunctions(k)%mainBranchHaloWeightsSquared(i,m)
                   do j=1,analysisHaloMassFunctionsMassBinsCount
                      if (i == j) cycle
                      massFunctions               (k)%massFunctionCovariance      (i,j)= &
                           &  +      massFunctions(k)%massFunctionCovariance      (i,j)  &
                           &  -      massFunctions(k)%mainBranchHaloWeights       (i,m)  &
                           &  *      massFunctions(k)%mainBranchHaloWeights       (j,m)  &
                           &  *sqrt(                                                     &
                           &         massFunctions(k)%mainBranchHaloWeightsSquared(i,m)  &
                           &        *massFunctions(k)%mainBranchHaloWeightsSquared(j,m)  &
                           &       )                                                     &
                           &  /haloWeightBinTotal
                   end do
                end do
             end if
          end do
       end if
       ! Truncate the covariance where the correlation is below threshold. This avoids creating noisy covariances matrices
       ! which can lead to discontinuities in the likelihood surface.
       do i=1,analysisHaloMassFunctionsMassBinsCount
          do j=1,analysisHaloMassFunctionsMassBinsCount
             if (i == j) cycle
             if     (                                                           &
                  &     abs( massFunctions(k)%massFunctionCovariance(i,j))      &
                  &  <                                                          &
                  &    analysisHaloMassFunctionsCorrelationTruncateLevel        &
                  &   *sqrt(                                                    &
                  &         max(                                                &
                  &              0.0d0                                       ,  &
                  &              massFunctions(k)%massFunctionCovariance(i,i)   &
                  &             *massFunctions(k)%massFunctionCovariance(j,j)   &
                  &            )                                                &
                  &        )                                                    &
                  & )        massFunctions(k)%massFunctionCovariance(i,j)=0.0d0
          end do
       end do
       ! Convert model mass function to differential per log(M).
       forall(i=1:analysisHaloMassFunctionsMassBinsCount)
          massFunctions(k)%massFunction             (i  )=  massFunctions(k)%massFunction          (i  )                              &
               &                         /(massFunctions(k)%massesLogarithmicMaximum(i)-massFunctions(k)%massesLogarithmicMinimum(i)) &
               &                         / log(10.0d0)
          forall(j=1:analysisHaloMassFunctionsMassBinsCount)
             massFunctions(k)%massFunctionCovariance(i,j)=  massFunctions(k)%massFunctionCovariance(i,j)                              &
                  &                      /(massFunctions(k)%massesLogarithmicMaximum(i)-massFunctions(k)%massesLogarithmicMinimum(i)) &
                  &                      /(massFunctions(k)%massesLogarithmicMaximum(j)-massFunctions(k)%massesLogarithmicMinimum(j)) &
                  &                      / log(10.0d0)**2
          end forall
       end forall
       ! Output the mass function.
       !$omp critical(HDF5_Access)
       analysisGroup    =galacticusOutputFile%openGroup('analysis','Model analysis')
       massFunctionGroup=analysisGroup       %openGroup(char(massFunctions(k)%label),"Halo mass function")
       call massFunctionGroup%writeDataset  (massFunctions(k)%masses                ,'mass'                  ,'Mass'                    ,datasetReturned=thisDataset)
       call thisDataset      %writeAttribute(massSolar             ,'unitsInSI'                                                                                     )
       call thisDataset      %close()
       call massFunctionGroup%writeDataset  (massFunctions(k)%massFunction          ,'massFunction'          ,'Mass function'           ,datasetReturned=thisDataset)
       call thisDataset      %writeAttribute(1.0d0/megaParsec**3   ,'unitsInSI'                                                                                     )
       call thisDataset      %close()
       call massFunctionGroup%writeDataset  (massFunctions(k)%massFunctionCovariance,'massFunctionCovariance','Mass function covariance',datasetReturned=thisDataset)
       call thisDataset      %writeAttribute(1.0d0/megaParsec**6   ,'unitsInSI'                                                                                     )
       call thisDataset      %close()
       call massFunctionGroup%close()
       call analysisGroup    %close()
       !$omp end critical(HDF5_Access)
    end do
    return
  end subroutine Galacticus_Output_Analysis_Halo_Mass_Functions_Output

end module Galacticus_Output_Analyses_Halo_Mass_Function
