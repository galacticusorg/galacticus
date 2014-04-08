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

!% Contains a module which performs analysis to compute a variety of mass-dependent size functions.

module Galacticus_Output_Analyses_Mass_Dpndnt_Sz_Dstrbtins
  !% Performs analysis to compute a variety of mass-dependent size functions. Currently supported mass functions include:
  !% \begin{itemize}
  !% \item The \gls{sdss} late-type galaxy size distributions from \cite{shen_size_2003}.
  !% \end{itemize}
  use Galacticus_Nodes
  use Galactic_Structure_Options
  use Numerical_Constants_Astronomical
  use Numerical_Constants_Prefixes
  implicit none
  private
  public :: Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins, Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins_Output

  ! Record of module initialization.
  logical                                                   :: moduleInitialized          =.false.

  ! Record of whether this analysis is active.
  logical                                                   :: analysisActive

  ! Number of supported mass functions.
  integer          , parameter                              :: sizeFunctionsSupportedCount=1

  ! Labels for supported mass functions.
  character(len=21), dimension(sizeFunctionsSupportedCount) :: sizeFunctionLabels= &
       & [                                                                         &
       &  'sdssSizeFunctionZ0.07'                                                  &
       & ]

  ! Interface for mass mapping functions.
  abstract interface
     double precision function Map_Mass(mass,thisNode)
       import treeNode
       double precision          , intent(in   )          :: mass
       type            (treeNode), intent(inout), pointer :: thisNode
     end function Map_Mass
  end interface

  ! Interface for radius mapping functions.
  abstract interface
     double precision function Map_Radius(radius,thisNode)
       import treeNode
       double precision          , intent(in   )          :: radius
       type            (treeNode), intent(inout), pointer :: thisNode
     end function Map_Radius
  end interface

  ! Interface for mass error functions.
  abstract interface
     double precision function Mass_Error(mass,thisNode)
       import treeNode
       double precision          , intent(in   )          :: mass
       type            (treeNode), intent(inout), pointer :: thisNode
     end function Mass_Error
  end interface

  ! Interface for radius error functions.
  abstract interface
     double precision function Radius_Error(radius,thisNode)
       import treeNode
       double precision          , intent(in   )          :: radius
       type            (treeNode), intent(inout), pointer :: thisNode
     end function Radius_Error
  end interface

  ! Type for descriptors of mass functions.
  type :: sizeFunctionDescriptor
     double precision                                :: redshift
     double precision                                :: massSystematicLogM0           , radiusSystematicLogR0
     double precision                                :: massRandomError               , radiusRandomError
     procedure       (  Mass_Error), pointer, nopass :: massRandomErrorFunction
     procedure       (Radius_Error), pointer, nopass :: radiusRandomErrorFunction
     double precision                                :: massLogarithmicMinimum
     integer                                         :: massSystematicCoefficientCount, radiusSystematicCoefficientCount
     integer                                         :: massType
     double precision                                :: massUnitsInSI                 , radiusUnitsInSI
     character       (len= 32     )                  :: label
     character       (len=128     )                  :: comment
     procedure       (Map_Mass    ), pointer, nopass :: mapMass
     procedure       (Map_Radius  ), pointer, nopass :: mapRadius
  end type sizeFunctionDescriptor

  ! Mass function descriptors.
  type(sizeFunctionDescriptor), dimension(sizeFunctionsSupportedCount), target :: sizeFunctionDescriptors= &
       & [                                                                                                 &
       ! SDSS late-type galaxy sizes from Shen et al. (2003). 
       &                           sizeFunctionDescriptor(                                                 &
       &                                                   0.070d0+0                         ,             &
       &                                                  11.0000d+0                         ,             &
       &                                                   1.0000d+0                         ,             &
       &                                                   0.0806d+0                         ,             &
       &                                                   0.0128d+0                         ,             &
       &                                                   null()                            ,             &
       &                                                   null()                            ,             &
       &                                                   6.500d0                           ,             &
       &                                                   2                                 ,             &
       &                                                   2                                 ,             &
       &                                                   massTypeStellar                   ,             &
       &                                                   massSolar                         ,             &
       &                                                   megaParsec/kilo                   ,             &
       &                                                   'sdssGalaxySizesZ0.07'            ,             &
       &                                                   'SDSS disk galaxy sizes at z=0.07',             &
       &                                                   null()                            ,             &
       &                                                   null()                                          &
       &                                                 )                                                 &
       & ]

  ! Type to store size functions.
  type :: sizeFunction
     ! Copy of the mass function descriptor for this mass function.
     type            (sizeFunctionDescriptor), pointer                       :: descriptor
     ! Parameters for the systematic error model.
     double precision                        , allocatable, dimension(:    ) :: massSystematicCoefficients, radiusSystematicCoefficients
     ! The index of the output corresponding to the required redshift.
     integer                                                                 :: outputNumber
     ! The number of bins.
     integer                                                                 :: massesCount               , radiiCount
     ! Arrays for the masses, radii and size function.
     double precision                        , allocatable, dimension(:    ) :: masses                    , massesLogarithmic       , &
          &                                                                     massesLogarithmicMinimum  , massesLogarithmicMaximum, &
          &                                                                     sizeFunctionWeights
     double precision                        , allocatable, dimension(:,:  ) :: radii                     , radiiLogarithmic        , &
          &                                                                     radiiLogarithmicMinimum   , radiiLogarithmicMaximum , &
          &                                                                     sizeFunction
     ! Arrays for accumulation of of main branch galaxies
     double precision                        , allocatable, dimension(:,:,:) :: mainBranchGalaxyWeights   , mainBranchGalaxyWeightsSquared
     ! Array for the covariance matrix.
     double precision                        , allocatable, dimension(:,:  ) :: sizeFunctionCovariance
     ! Cosmology conversion factors.
     double precision                                                        :: cosmologyConversionMass   , cosmologyConversionSizeFunction, &
          &                                                                     cosmologyConversionSize
  end type sizeFunction

  ! Mass functions.
  type(sizeFunction), allocatable, dimension(:) :: sizeFunctions

  ! Type for storing temporary size functions during cumulation.
  type :: sizeFunctionWork
     double precision, allocatable, dimension(:,:) :: sizeFunction
     double precision, allocatable, dimension(  :) :: sizeFunctionWeights
     double precision, allocatable, dimension(:,:) :: covariance
  end type sizeFunctionWork

  ! Work array.
  type(sizeFunctionWork), allocatable, dimension(:) :: thisGalaxy
  !$omp threadprivate(thisGalaxy)

  ! Options controlling binning in halo mass.
  integer                     :: analysisSizeFunctionCovarianceModel
  integer         , parameter :: analysisSizeFunctionCovarianceModelPoisson =1
  integer         , parameter :: analysisSizeFunctionCovarianceModelBinomial=2
  integer                     :: analysisSizeFunctionsHaloMassBinsCount                 , analysisSizeFunctionsHaloMassBinsPerDecade
  double precision            :: analysisSizeFunctionsHaloMassMinimum                   , analysisSizeFunctionsHaloMassMaximum           , &
       &                         analysisSizeFunctionsHaloMassIntervalLogarithmicInverse, analysisSizeFunctionsHaloMassMinimumLogarithmic

contains

  !# <mergerTreeAnalysisTask>
  !#  <unitName>Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins</unitName>
  !# </mergerTreeAnalysisTask>
  subroutine Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins(thisTree,thisNode,iOutput,mergerTreeAnalyses)
    !% Construct a mass functions to compare to various observational determinations.
    use Galacticus_Nodes
    use Galacticus_Input_Paths
    use ISO_Varying_String
    use Memory_Management
    use Cosmology_Parameters
    use Galactic_Structure_Enclosed_Masses
    use Input_Parameters
    use Galacticus_Output_Times
    use Galacticus_Error
    use Cosmology_Functions
    use Numerical_Comparison
    use Dark_Matter_Halo_Scales
    use String_Handling
    use IO_HDF5
    use Vectors
    use Galacticus_Output_Analyses_Cosmology_Scalings
    implicit none
    type            (mergerTree                    ), intent(in   )                 :: thisTree
    type            (treeNode                      ), intent(inout), pointer        :: thisNode
    integer                                         , intent(in   )                 :: iOutput
    type            (varying_string                ), intent(in   ), dimension(:  ) :: mergerTreeAnalyses
    class           (nodeComponentBasic            )               , pointer        :: thisBasic
    class           (nodeComponentDisk             )               , pointer        :: thisDisk
    class           (cosmologyFunctionsClass       )               , pointer        :: cosmologyFunctionsModel
    type            (cosmologyFunctionsMatterLambda)                                :: cosmologyFunctionsObserved
    type            (cosmologyParametersSimple     )                                :: cosmologyParametersObserved
    integer                                                                         :: i,j,k,currentAnalysis,activeAnalysisCount,haloMassBin,iDistribution,jDistribution
    double precision                                                                :: dataHubbleParameter ,mass,massLogarithmic&
         &,massRandomError,radiusLogarithmic,radius,sizeRandomError,dataOmegaDarkEnergy,dataOmegaMatter,sersicIndexMaximum
    type            (varying_string                )                                :: parameterName&
         &,analysisSizeFunctionCovarianceModelText,cosmologyScalingSizeFunction,cosmologyScalingMass,cosmologyScalingSize
    character       (len=128                       )                                :: distributionGroupName
    logical                                                                         :: groupFound
    type            (hdf5Object                    )                                :: dataFile,sizeDataset,distributionGroup,cosmologyGroup

    ! Initialize the module if necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins_Initialize)
       if (.not.moduleInitialized) then
          !@ <inputParameter>
          !@   <name>analysisSizeFunctionCovarianceModel</name>
          !@   <defaultValue>binomial</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The model to use when constructing the mass function covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisSizeFunctionCovarianceModel',analysisSizeFunctionCovarianceModelText,defaultValue='Poisson')
          select case (char(analysisSizeFunctionCovarianceModelText))
          case ( 'Poisson'  )
             analysisSizeFunctionCovarianceModel=analysisSizeFunctionCovarianceModelPoisson
          case ( 'binomial' )
             analysisSizeFunctionCovarianceModel=analysisSizeFunctionCovarianceModelBinomial
          case default
             call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins','unrecognized value for "analysisSizeFunctionCovarianceModel" - allowed values are "Poisson", and "binomial"')
          end select
          !@ <inputParameter>
          !@   <name>analysisSizeFunctionsHaloMassBinsPerDecade</name>
          !@   <defaultValue>10</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The number of bins per decade of halo mass to use when constructing the mass function covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisSizeFunctionsHaloMassBinsPerDecade',analysisSizeFunctionsHaloMassBinsPerDecade,defaultValue=10)
          !@ <inputParameter>
          !@   <name>analysisSizeFunctionsHaloMassMinimum</name>
          !@   <defaultValue>$10^8M_\odot$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The minimum halo mass to consider when constructing the mass function covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisSizeFunctionsHaloMassMinimum',analysisSizeFunctionsHaloMassMinimum,defaultValue=1.0d8)
          !@ <inputParameter>
          !@   <name>analysisSizeFunctionsHaloMassMaximum</name>
          !@   <defaultValue>$10^{16}M_\odot$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The maximum halo mass to consider when constructing the mass function covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisSizeFunctionsHaloMassMaximum',analysisSizeFunctionsHaloMassMaximum,defaultValue=1.0d16)
          analysisSizeFunctionsHaloMassMinimumLogarithmic        =     log10( analysisSizeFunctionsHaloMassMinimum)
          analysisSizeFunctionsHaloMassBinsCount                 =int(                                                  &
               &                                                       log10(                                           &
               &                                                              analysisSizeFunctionsHaloMassMaximum      &
               &                                                             /analysisSizeFunctionsHaloMassMinimum      &
               &                                                            )                                           &
               &                                                      *dble(analysisSizeFunctionsHaloMassBinsPerDecade) &
               &                                                      +0.5d0&
               &                                                     )
          analysisSizeFunctionsHaloMassIntervalLogarithmicInverse= dble (analysisSizeFunctionsHaloMassBinsCount)        &
               &                                                  /log10(                                               &
               &                                                          analysisSizeFunctionsHaloMassMaximum          &
               &                                                         /analysisSizeFunctionsHaloMassMinimum          &
               &                                                        )
          ! Establish mapping functions for size function descriptors.
          sizeFunctionDescriptors(1)%mapRadius => Map_Radius_SDSS_Size_Function_Z0_07
          ! Determine how many supported mass functions are requested.
          activeAnalysisCount=0
          do i=1,sizeFunctionsSupportedCount
             if (any(trim(mergerTreeAnalyses) == trim(sizeFunctionLabels(i)))) activeAnalysisCount=activeAnalysisCount+1
          end do
          ! Allocate mass function arrays and populate with required data.
          if (activeAnalysisCount <= 0) then
             analysisActive=.false.
          else
             analysisActive=.true.
             currentAnalysis=0
             allocate(sizeFunctions(activeAnalysisCount))
             cosmologyFunctionsModel => cosmologyFunctions()
             do i=1,size(mergerTreeAnalyses)
                do j=1,sizeFunctionsSupportedCount
                   if (mergerTreeAnalyses(i) == trim(sizeFunctionLabels(j))) then
                      currentAnalysis=currentAnalysis+1
                      ! Set a pointer to the descriptor for this size function.
                      sizeFunctions(currentAnalysis)%descriptor => sizeFunctionDescriptors(j)
                      ! Read parameters of the systematic error model.
                      if (sizeFunctionDescriptors(j)%massSystematicCoefficientCount > 0) then
                         allocate(sizeFunctions(currentAnalysis)%massSystematicCoefficients(sizeFunctionDescriptors(j)%massSystematicCoefficientCount))
                         do k=1,sizeFunctionDescriptors(j)%massSystematicCoefficientCount
                            parameterName=trim(sizeFunctionLabels(j))//'MassSystematic'
                            parameterName=parameterName//(k-1)
                            call Get_Input_Parameter(char(parameterName),sizeFunctions(currentAnalysis)%massSystematicCoefficients(k),defaultValue=0.0d0)
                         end do
                      end if
                      if (sizeFunctionDescriptors(j)%radiusSystematicCoefficientCount > 0) then
                         allocate(sizeFunctions(currentAnalysis)%radiusSystematicCoefficients(sizeFunctionDescriptors(j)%radiusSystematicCoefficientCount))
                         do k=1,sizeFunctionDescriptors(j)%radiusSystematicCoefficientCount
                            parameterName=trim(sizeFunctionLabels(j))//'RadiusSystematic'
                            parameterName=parameterName//(k-1)
                            call Get_Input_Parameter(char(parameterName),sizeFunctions(currentAnalysis)%radiusSystematicCoefficients(k),defaultValue=0.0d0)
                         end do
                      end if
                      ! Find which output number corresponds to the required redshift.
                      sizeFunctions(currentAnalysis)%outputNumber=-1
                      do k=1,Galacticus_Output_Time_Count()
                         if     (                                                                                             &
                              &  Values_Agree(                                                                                &
                              &               sizeFunctionDescriptors(j)%redshift                                           , &
                              &               cosmologyFunctionsModel%redshiftFromExpansionFactor(                            &
                              &               cosmologyFunctionsModel%expansionFactor             (                           &
                              &                                                                    Galacticus_Output_Time(k)  &
                              &                                                                   )                           &
                              &                                                                  )                          , &
                              &               absTol=0.001d0                                                                  &
                              &              )                                                                                &
                              & ) then
                            sizeFunctions(currentAnalysis)%outputNumber=k
                            exit
                         end if
                      end do
                      if (sizeFunctions(currentAnalysis)%outputNumber < 0)                                                    &
                           & call Galacticus_Error_Report(                                                                    &
                           &                                'Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins'  , &
                           &                                'unable to find required redshift in outputs for mass function '  &
                           &                              //trim(sizeFunctionLabels(j))                                       &
                           &                             )
                      ! Read the appropriate observational data definition.
                      select case (trim(sizeFunctionLabels(j)))
                      case ('sdssSizeFunctionZ0.07')
                         ! SDSS z=0.07 size function.
                         call dataFile%openFile(char(Galacticus_Input_Path())//"data/observations/galaxySizes/Galaxy_Sizes_By_Mass_SDSS_Shen_2003.hdf5",readOnly=.true.)
                         ! Count number of distributions.
                         sizeFunctions(currentAnalysis)%massesCount=0
                         groupFound=.true.
                         jDistribution=0
                         do while (groupFound)
                            jDistribution=jDistribution+1
                            write (distributionGroupName,'(a,i2.2)')  &
                                 & "distribution"                   , &
                                 & jDistribution
                            groupFound=dataFile%hasGroup(distributionGroupName)
                            if (groupFound) then
                               distributionGroup=dataFile%openGroup(distributionGroupName)
                               call distributionGroup%readAttribute("sersicIndexMaximum",sersicIndexMaximum)
                               call distributionGroup%close        (                                       )
                               if (sersicIndexMaximum <= 2.5d0) &
                                    &  sizeFunctions(currentAnalysis)%massesCount &
                                    & =sizeFunctions(currentAnalysis)%massesCount &
                                    & +1
                            end if
                         end do
                         ! Count radii.
                         distributionGroup=dataFile         %openGroup  ("distribution01")
                         sizeDataset      =distributionGroup%openDataset("radius"        )
                         sizeFunctions(currentAnalysis)%radiiCount=sizeDataset%size(dim=1)
                         call sizeDataset      %close()
                         call distributionGroup%close()
                         ! Construct arrays.
                         call Alloc_Array(sizeFunctions(currentAnalysis)%masses                        ,[                                          sizeFunctions(currentAnalysis)%massesCount                                       ])
                         call Alloc_Array(sizeFunctions(currentAnalysis)%massesLogarithmic             ,[                                          sizeFunctions(currentAnalysis)%massesCount                                       ])
                         call Alloc_Array(sizeFunctions(currentAnalysis)%massesLogarithmicMinimum      ,[                                          sizeFunctions(currentAnalysis)%massesCount                                       ])
                         call Alloc_Array(sizeFunctions(currentAnalysis)%massesLogarithmicMaximum      ,[                                          sizeFunctions(currentAnalysis)%massesCount                                       ])
                         call Alloc_Array(sizeFunctions(currentAnalysis)%radii                         ,[sizeFunctions(currentAnalysis)%radiiCount,sizeFunctions(currentAnalysis)%massesCount                                       ])
                         call Alloc_Array(sizeFunctions(currentAnalysis)%radiiLogarithmic              ,[sizeFunctions(currentAnalysis)%radiiCount,sizeFunctions(currentAnalysis)%massesCount                                       ])
                         call Alloc_Array(sizeFunctions(currentAnalysis)%radiiLogarithmicMinimum       ,[sizeFunctions(currentAnalysis)%radiiCount,sizeFunctions(currentAnalysis)%massesCount                                       ])
                         call Alloc_Array(sizeFunctions(currentAnalysis)%radiiLogarithmicMaximum       ,[sizeFunctions(currentAnalysis)%radiiCount,sizeFunctions(currentAnalysis)%massesCount                                       ])
                         call Alloc_Array(sizeFunctions(currentAnalysis)%sizeFunction                  ,[sizeFunctions(currentAnalysis)%radiiCount,sizeFunctions(currentAnalysis)%massesCount                                       ])
                         call Alloc_Array(sizeFunctions(currentAnalysis)%sizeFunctionWeights           ,[                                          sizeFunctions(currentAnalysis)%massesCount                                       ])
                         call Alloc_Array(sizeFunctions(currentAnalysis)%sizeFunctionCovariance        ,[                                                                                                                             &
                              &                                                                          sizeFunctions(currentAnalysis)%radiiCount*sizeFunctions(currentAnalysis)%massesCount,                                        &
                              &                                                                          sizeFunctions(currentAnalysis)%radiiCount*sizeFunctions(currentAnalysis)%massesCount                                         &
                              &                                                                         ]                                                                                                                             &
                              &          )
                         call Alloc_Array(sizeFunctions(currentAnalysis)%mainBranchGalaxyWeights       ,[sizeFunctions(currentAnalysis)%radiiCount,sizeFunctions(currentAnalysis)%massesCount,analysisSizeFunctionsHaloMassBinsCount])
                         call Alloc_Array(sizeFunctions(currentAnalysis)%mainBranchGalaxyWeightsSquared,[sizeFunctions(currentAnalysis)%radiiCount,sizeFunctions(currentAnalysis)%massesCount,analysisSizeFunctionsHaloMassBinsCount])
                         ! Read datasets.
                         jDistribution=0
                         iDistribution=0
                         groupFound=.true.
                         do while (groupFound)
                            jDistribution=jDistribution+1
                            write (distributionGroupName,'(a,i2.2)')  &
                                 & "distribution"                   , &
                                 & jDistribution
                            groupFound=dataFile%hasGroup(distributionGroupName)
                            if (groupFound) then
                               distributionGroup=dataFile%openGroup(distributionGroupName)
                               call distributionGroup%readAttribute("sersicIndexMaximum",sersicIndexMaximum)
                               if (sersicIndexMaximum <= 2.5d0) then
                                  iDistribution=iDistribution+1
                                  call distributionGroup%readAttribute    ("massMinimum",sizeFunctions(currentAnalysis)%massesLogarithmicMinimum(  iDistribution))
                                  call distributionGroup%readAttribute    ("massMaximum",sizeFunctions(currentAnalysis)%massesLogarithmicMaximum(  iDistribution))
                                  call distributionGroup%readDatasetStatic("radius"     ,sizeFunctions(currentAnalysis)%radii                   (:,iDistribution))
                                  if (iDistribution == 1) then
                                     sizeDataset      =distributionGroup%openDataset("radius"             )
                                     call sizeDataset      %readAttribute("cosmologyScaling"    ,cosmologyScalingSize        )
                                     call sizeDataset      %close        (                                                   )
                                     call distributionGroup%readAttribute("massCosmologyScaling",cosmologyScalingMass        )
                                     call sizeDataset      %close        (                                                   )
                                     sizeDataset      =distributionGroup%openDataset("radiusFunction"     )
                                     call sizeDataset      %readAttribute("cosmologyScaling"    ,cosmologyScalingSizeFunction)
                                     call sizeDataset      %close        (                                                   )
                                  end if
                               end if
                               call distributionGroup%close()
                            end if
                         end do
                         sizeFunctions              (currentAnalysis)%massesLogarithmicMinimum  &
                              & =log10(sizeFunctions(currentAnalysis)%massesLogarithmicMinimum)
                         sizeFunctions              (currentAnalysis)%massesLogarithmicMaximum  &
                              & =log10(sizeFunctions(currentAnalysis)%massesLogarithmicMaximum)
                         ! Extract cosmological parameters.
                         cosmologyGroup=dataFile%openGroup("cosmology")
                         call cosmologyGroup%readAttribute("H_0"         ,dataHubbleParameter)
                         call cosmologyGroup%readAttribute("Omega_Matter",dataOmegaMatter    )
                         call cosmologyGroup%readAttribute("Omega_DE"    ,dataOmegaDarkEnergy)
                         call cosmologyGroup%close        (                                  )
                         ! Finished reading data.
                         call dataFile%close()
                         ! Create the observed cosmology.
                         cosmologyParametersObserved=cosmologyParametersSimple     (                                     &
                              &                                                     OmegaMatter    =dataOmegaMatter    , &
                              &                                                     OmegaDarkEnergy=dataOmegaDarkEnergy, &
                              &                                                     HubbleConstant =dataHubbleParameter, &
                              &                                                     temperatureCMB =0.0d0              , &
                              &                                                     OmegaBaryon    =0.0d0                &
                              &                                                    )
                         cosmologyFunctionsObserved =cosmologyFunctionsMatterLambda(                                     &
                              &                                                     cosmologyParametersObserved          &
                              &                                                    )
                         ! Adjust masses to current Hubble parameter.
                         sizeFunctions(currentAnalysis)%massesLogarithmic             = 0.5d0* (sizeFunctions(currentAnalysis)%massesLogarithmicMinimum+sizeFunctions(currentAnalysis)%massesLogarithmicMaximum)
                         sizeFunctions(currentAnalysis)%masses                        =10.0d0** sizeFunctions(currentAnalysis)%massesLogarithmic
                         sizeFunctions(currentAnalysis)%radiiLogarithmic              =log10   (sizeFunctions(currentAnalysis)%radii                                                                           )
                         sizeFunctions(currentAnalysis)%sizeFunction                  =0.0d0
                         sizeFunctions(currentAnalysis)%sizeFunctionWeights           =0.0d0
                         sizeFunctions(currentAnalysis)%sizeFunctionCovariance        =0.0d0
                         sizeFunctions(currentAnalysis)%mainBranchGalaxyWeights       =0.0d0
                         sizeFunctions(currentAnalysis)%mainBranchGalaxyWeightsSquared=0.0d0
                         do k=1,sizeFunctions(currentAnalysis)%radiiCount
                            if (k ==                                          1) then
                               sizeFunctions(currentAnalysis)%radiiLogarithmicMinimum(k,:)=sizeFunctions(currentAnalysis)%radiiLogarithmic(k,:)-0.5d0*(sizeFunctions(currentAnalysis)%radiiLogarithmic(k+1,:)-sizeFunctions(currentAnalysis)%radiiLogarithmic(k  ,:))
                            else
                               sizeFunctions(currentAnalysis)%radiiLogarithmicMinimum(k,:)=                                                    +0.5d0*(sizeFunctions(currentAnalysis)%radiiLogarithmic(k-1,:)+sizeFunctions(currentAnalysis)%radiiLogarithmic(k  ,:))
                            end if
                            if (k == sizeFunctions(currentAnalysis)%radiiCount) then
                               sizeFunctions(currentAnalysis)%radiiLogarithmicMaximum(k,:)=sizeFunctions(currentAnalysis)%radiiLogarithmic(k,:)+0.5d0*(sizeFunctions(currentAnalysis)%radiiLogarithmic(k  ,:)-sizeFunctions(currentAnalysis)%radiiLogarithmic(k-1,:))
                            else
                               sizeFunctions(currentAnalysis)%radiiLogarithmicMaximum(k,:)=                                                    +0.5d0*(sizeFunctions(currentAnalysis)%radiiLogarithmic(k+1,:)+sizeFunctions(currentAnalysis)%radiiLogarithmic(k  ,:))
                            end if
                         end do
                      case default
                         call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins','unknown size function')
                      end select
                      ! Get cosmological conversion factors.
                      call Cosmology_Conversion_Factors(                                                                                                &
                           &                            sizeFunctions(currentAnalysis)%descriptor%redshift                                            , &
                           &                            cosmologyFunctionsModel                                                                       , &
                           &                            cosmologyFunctionsObserved                                                                    , &
                           &                            cosmologyScalingMass           =cosmologyScalingMass                                          , &
                           &                            cosmologyScalingSize           =cosmologyScalingSize                                          , &
                           &                            cosmologyScalingMassFunction   =cosmologyScalingSizeFunction                                  , &
                           &                            cosmologyConversionMass        =sizeFunctions(currentAnalysis)%cosmologyConversionMass        , &
                           &                            cosmologyConversionSize        =sizeFunctions(currentAnalysis)%cosmologyConversionSize        , &
                           &                            cosmologyConversionMassFunction=sizeFunctions(currentAnalysis)%cosmologyConversionSizeFunction  &
                           &                           )
                     exit
                  end if
                end do
             end do
             ! Ensure that disk component supports radius property.
             if (.not.defaultDiskComponent%radiusIsGettable()) &
                  & call Galacticus_Error_Report                                                                   &
                  & (                                                                                              &
                  &  'Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins'                                       , &
                  &  'This analysis requires that the "radius" property of the disk is gettable.'//                &
                  &  Galacticus_Component_List(                                                                    &
                  &                            'disk'                                                            , &
                  &                             defaultDiskComponent%radiusAttributeMatch(requireGettable=.true.)  &
                  &                           )                                                                    &
                  & )             
          end if
          ! Record that module is initialized.
          moduleInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins_Initialize)
    end if
    ! Return if this analysis is not active.
    if (.not.analysisActive) return
    ! Allocate work arrays.
    if (.not.allocated(thisGalaxy)) allocate(thisGalaxy(size(sizeFunctions)))
    ! Iterate over active analyses.
    do i=1,size(sizeFunctions)
       ! Skip if this analysis is not active, or if this is not the correct output.
       if (iOutput /= sizeFunctions(i)%outputNumber) cycle
       ! Allocate workspace.
       if (.not.allocated(thisGalaxy(i)%sizeFunction       )) call Alloc_Array(thisGalaxy(i)%sizeFunction       ,[sizeFunctions(i)%radiiCount,sizeFunctions(i)%massesCount])
       if (.not.allocated(thisGalaxy(i)%sizeFunctionWeights)) call Alloc_Array(thisGalaxy(i)%sizeFunctionWeights,[                            sizeFunctions(i)%massesCount])
       if (.not.allocated(thisGalaxy(i)%covariance         )) call Alloc_Array(thisGalaxy(i)%covariance         ,[                                                          &
            &                                                                                       sizeFunctions(i)%radiiCount*sizeFunctions(i)%massesCount, &
            &                                                                                       sizeFunctions(i)%radiiCount*sizeFunctions(i)%massesCount  &
            &                                                                                      ]                                                          &
            &                                                          )
       ! Get the galactic mass.
       mass=                                                                                                                &
            &  Galactic_Structure_Enclosed_Mass(thisNode,radiusLarge,componentType=componentTypeDisk    ,massType=sizeFunctions(i)%descriptor%massType) &
            & +Galactic_Structure_Enclosed_Mass(thisNode,radiusLarge,componentType=componentTypeSpheroid,massType=sizeFunctions(i)%descriptor%massType)
       if (mass            <=                  0.0d0) return
       if (associated(sizeFunctions(i)%descriptor%mapMass)) mass=sizeFunctions(i)%descriptor%mapMass(mass,thisNode)
       mass=mass*sizeFunctions(i)%cosmologyConversionMass ! Convert for cosmology.
       massLogarithmic=log10(mass)
       do j=1,sizeFunctions(i)%descriptor%massSystematicCoefficientCount
          massLogarithmic=massLogarithmic+sizeFunctions(i)%massSystematicCoefficients(j)*(log10(mass)-sizeFunctions(i)%descriptor%massSystematicLogM0)**(j-1)
       end do
       if (massLogarithmic <  sizeFunctions(i)%descriptor%massLogarithmicMinimum) return
       ! Get the galactic radius.
       thisDisk => thisNode%disk  ()
       radius   =  thisDisk%radius()
       if (associated(sizeFunctions(i)%descriptor%mapRadius)) radius=sizeFunctions(i)%descriptor%mapRadius(radius,thisNode)
       radius=radius*sizeFunctions(i)%cosmologyConversionSize ! Convert for cosmology.
       radiusLogarithmic=log10(radius)
       do j=1,sizeFunctions(i)%descriptor%radiusSystematicCoefficientCount
          radiusLogarithmic=radiusLogarithmic+sizeFunctions(i)%radiusSystematicCoefficients(j)*(log10(radius)-sizeFunctions(i)%descriptor%radiusSystematicLogR0)**(j-1)
       end do
        ! Compute contributions to each bin.
       massRandomError=sizeFunctions(i)%descriptor%massRandomError
       sizeRandomError=sizeFunctions(i)%descriptor%radiusRandomError
       if (associated(sizeFunctions(i)%descriptor%  massRandomErrorFunction)) massRandomError=sizeFunctions(i)%descriptor%massRandomErrorFunction  (mass  ,thisNode)
       if (associated(sizeFunctions(i)%descriptor%radiusRandomErrorFunction)) sizeRandomError=sizeFunctions(i)%descriptor%radiusRandomErrorFunction(radius,thisNode)
       thisGalaxy(i)%sizeFunction=       (                                                                                               &
            &                             +erf((sizeFunctions(i)%radiiLogarithmicMaximum-radiusLogarithmic)/sizeRandomError/sqrt(2.0d0)) &
            &                             -erf((sizeFunctions(i)%radiiLogarithmicMinimum-radiusLogarithmic)/sizeRandomError/sqrt(2.0d0)) &
            &                            )                                                                                               &
            &                            /2.0d0
       thisGalaxy(i)%sizeFunctionWeights=(                                                                                               &
            &                             +erf((sizeFunctions(i)%massesLogarithmicMaximum- massLogarithmic)/massRandomError/sqrt(2.0d0)) &
            &                             -erf((sizeFunctions(i)%massesLogarithmicMinimum- massLogarithmic)/massRandomError/sqrt(2.0d0)) &
            &                            )                                                                                               &
            &                            /2.0d0                                                                                          &
            &                            *thisTree%volumeWeight                                                                          &
            &                            *sizeFunctions(i)%cosmologyConversionSizeFunction
       do j=1,sizeFunctions(i)%massesCount
          thisGalaxy(i)%sizeFunction(:,j)=thisGalaxy(i)%sizeFunction(:,j)*thisGalaxy(i)%sizeFunctionWeights(j)
       end do       
       ! Accumulate size function.
       !$omp critical (Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins_Accumulate)
       sizeFunctions(i)%sizeFunction       =sizeFunctions(i)%sizeFunction       +thisGalaxy(i)%sizeFunction
       sizeFunctions(i)%sizeFunctionWeights=sizeFunctions(i)%sizeFunctionWeights+thisGalaxy(i)%sizeFunctionWeights
       !$omp end critical (Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins_Accumulate)
       ! Treat main branch and other galaxies differently.
       if (thisNode%isOnMainBranch().and.analysisSizeFunctionCovarianceModel == analysisSizeFunctionCovarianceModelBinomial) then
          ! Find the bin to which this halo mass belongs.
          thisBasic => thisNode%basic()
          haloMassBin=floor((log10(thisBasic%mass())-analysisSizeFunctionsHaloMassMinimumLogarithmic)*analysisSizeFunctionsHaloMassIntervalLogarithmicInverse)+1
          ! Accumulate weights to halo mass arrays.
          if (haloMassBin >= 1 .and. haloMassBin <= analysisSizeFunctionsHaloMassBinsCount) then
            !$omp critical (Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins_Accumulate)
             sizeFunctions        (i)%mainBranchGalaxyWeights       (:,:,haloMassBin)= &
                  &  sizeFunctions(i)%mainBranchGalaxyWeights       (:,:,haloMassBin)  &
                  &  +thisGalaxy  (i)%sizeFunction
             sizeFunctions        (i)%mainBranchGalaxyWeightsSquared(:,:,haloMassBin)= &
                  &  sizeFunctions(i)%mainBranchGalaxyWeightsSquared(:,:,haloMassBin)  &
                  &  +thisGalaxy  (i)%sizeFunction**2
             !$omp end critical (Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins_Accumulate)
          end if
       else
          thisGalaxy(i)%covariance=                                                                                                   &
               & Vector_Outer_Product(                                                                                                &
               &                      reshape(thisGalaxy(i)%sizeFunction,[sizeFunctions(i)%massesCount*sizeFunctions(i)%radiiCount]), &
               &                      reshape(thisGalaxy(i)%sizeFunction,[sizeFunctions(i)%massesCount*sizeFunctions(i)%radiiCount])  &
               &                     )
          ! Accumulate covariance.
          !$omp critical (Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins_Accumulate)
          sizeFunctions(i)%sizeFunctionCovariance=sizeFunctions(i)%sizeFunctionCovariance+thisGalaxy(i)%covariance
          !$omp end critical (Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins_Accumulate)
       end if
    end do
    return
  end subroutine Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins

  !# <hdfPreCloseTask>
  !#  <unitName>Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins_Output</unitName>
  !# </hdfPreCloseTask>
  subroutine Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins_Output
    !% Outputs SDSS $z\approx 0.07$ stellar mass function to file.
    use Galacticus_HDF5
    implicit none
    integer                      :: k,m,mi,ri,mj,rj,ci,cj
    type            (hdf5Object) :: analysisGroup,sizeFunctionGroup,thisDataset
    double precision             :: haloWeightBinTotal

    ! Return immediately if this analysis is not active.
    if (.not.analysisActive) return
    ! Iterate over mass functions.
    do k=1,size(sizeFunctions)
       ! Add the contribution from main branch galaxies to the covariance matrix.
       if (analysisSizeFunctionCovarianceModel == analysisSizeFunctionCovarianceModelBinomial) then
          do m=1,analysisSizeFunctionsHaloMassBinsCount
             haloWeightBinTotal=sum(sizeFunctions(k)%mainBranchGalaxyWeights(:,:,m))
             if ( haloWeightBinTotal > 0.0 ) then
                do mi=1,sizeFunctions(k)%massesCount
                   do ri=1,sizeFunctions(k)%radiiCount
                      ci=(mi-1)*sizeFunctions(k)%radiiCount+ri
                      sizeFunctions               (k)%sizeFunctionCovariance        (ci,ci)=                      &
                           &         sizeFunctions(k)%sizeFunctionCovariance        (ci,ci)                       &
                           & +(1.0d0-sizeFunctions(k)%mainBranchGalaxyWeights       (ri,mi,m)/haloWeightBinTotal) &
                           & *       sizeFunctions(k)%mainBranchGalaxyWeightsSquared(ri,mi,m)
                      do mj=1,sizeFunctions(k)%massesCount
                         do rj=1,sizeFunctions(k)%radiiCount
                            cj=(mj-1)*sizeFunctions(k)%radiiCount+rj
                            if (mi == mj .and. ri == rj) cycle
                            sizeFunctions         (k)%sizeFunctionCovariance        (ci,cj)=                      &
                                 &   sizeFunctions(k)%sizeFunctionCovariance        (ci,cj)                       &
                                 & -(sizeFunctions(k)%mainBranchGalaxyWeights       (rj,mj,m)/haloWeightBinTotal) &
                                 & * sizeFunctions(k)%mainBranchGalaxyWeightsSquared(ri,mi,m)
                         end do
                      end do
                   end do
                end do
             end if
          end do
       end if
       ! Normalize the model size function in each mass interval and convert model size function to differential per log10(R).
       do mi=1,sizeFunctions(k)%massesCount 
          if (sizeFunctions(k)%sizeFunctionWeights(mi) > 0.0d0) then
             do ri=1,sizeFunctions(k)%radiiCount
                ci=(mi-1)*sizeFunctions(k)%radiiCount+ri
                sizeFunctions               (k)%sizeFunction           (ri,mi)=sizeFunctions(k)%sizeFunction           (ri,mi)  &
                     &       /(sizeFunctions(k)%radiiLogarithmicMaximum(ri,mi)-sizeFunctions(k)%radiiLogarithmicMinimum(ri,mi)) &
                     &       / sizeFunctions(k)%sizeFunctionWeights    (   mi)
                do mj=1,sizeFunctions(k)%massesCount
                   if (sizeFunctions(k)%sizeFunctionWeights(mj) > 0.0d0) then
                      do rj=1,sizeFunctions(k)%radiiCount 
                         cj=(mj-1)*sizeFunctions(k)%radiiCount+rj
                         sizeFunctions         (k)%sizeFunctionCovariance (ci,cj)=sizeFunctions(k)%sizeFunctionCovariance (ci,cj)  &
                              & /(sizeFunctions(k)%radiiLogarithmicMaximum(ri,mi)-sizeFunctions(k)%radiiLogarithmicMinimum(ri,mi)) &
                              & /(sizeFunctions(k)%radiiLogarithmicMaximum(rj,mj)-sizeFunctions(k)%radiiLogarithmicMinimum(rj,mj)) &
                              & / sizeFunctions(k)%sizeFunctionWeights    (   mi)                                                  &
                              & / sizeFunctions(k)%sizeFunctionWeights    (   mj)
                      end do
                   end if
                end do
             end do
          end if
       end do
       ! Output the size function.
       !$omp critical(HDF5_Access)
       analysisGroup    =galacticusOutputFile%openGroup('analysis','Model analysis')
       sizeFunctionGroup=analysisGroup       %openGroup(trim(sizeFunctions(k)%descriptor%label),trim(sizeFunctions(k)%descriptor%comment))
       call sizeFunctionGroup%writeDataset  (sizeFunctions(k)%masses                    ,'mass'                  ,'Mass'                     ,datasetReturned=thisDataset)
       call thisDataset      %writeAttribute(sizeFunctions(k)%descriptor%massUnitsInSI  ,'unitsInSI'                                                                     )
       call thisDataset      %close()
       call sizeFunctionGroup%writeDataset  (sizeFunctions(k)%radii                     ,'radius'                ,'Radius'                   ,datasetReturned=thisDataset)
       call thisDataset      %writeAttribute(sizeFunctions(k)%descriptor%radiusUnitsInSI,'unitsInSI'                                                                     )
       call thisDataset      %close()
       call sizeFunctionGroup%writeDataset  (sizeFunctions(k)%sizeFunction              ,'sizeFunction'          ,'Mass function'                                        )
       call sizeFunctionGroup%writeDataset  (sizeFunctions(k)%sizeFunctionCovariance    ,'sizeFunctionCovariance','Mass function covariance'                             )
       call sizeFunctionGroup%close()
       call analysisGroup    %close()
       !$omp end critical(HDF5_Access)
    end do
    return
  end subroutine Galacticus_Output_Analysis_Mass_Dpndnt_Sz_Dstrbtins_Output

  double precision function Map_Radius_SDSS_Size_Function_Z0_07(radius,thisNode)
    !% Maps scale radii into Petrosian $r_{50}$ radii for the SDSS disk size analysis. Also converts from Mpc to kpc.
    implicit none
    double precision          , intent(in   )          :: radius
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , parameter              :: diskScaleLengthToPetrosianR50=1.667632104d0

    Map_Radius_SDSS_Size_Function_Z0_07=radius*diskScaleLengthToPetrosianR50*kilo
    return
  end function Map_Radius_SDSS_Size_Function_Z0_07
  
end module Galacticus_Output_Analyses_Mass_Dpndnt_Sz_Dstrbtins
