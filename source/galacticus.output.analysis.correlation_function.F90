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

!% Contains a module which performs analysis to compute a variety of correlation functions.

module Galacticus_Output_Analyses_Correlation_Functions
  !% Performs analysis to compute a variety of correlation functions.
  use Galacticus_Nodes
  use Galactic_Structure_Options
  use Numerical_Constants_Astronomical
  use Numerical_Constants_Prefixes
  implicit none
  private
  public :: Galacticus_Output_Analysis_Correlation_Functions, Galacticus_Output_Analysis_Correlation_Functions_Output

  ! Record of module initialization.
  logical                                                          :: moduleInitialized               =.false.

  ! Record of whether this analysis is active.
  logical                                                          :: analysisActive

  ! Number of supported mass functions.
  integer          , parameter                                     :: correlationFunctionsSupportedCount=1

  ! Labels for supported mass functions.
  character(len=25), dimension(correlationFunctionsSupportedCount) :: correlationFunctionLabels= &
       & [                                                                                       &
       &  'sdssClusteringZ0.07'                                                                  &
       & ]

  ! Interface for mass mapping functions.
  abstract interface
     double precision function Map_Mass(mass,thisNode)
       import treeNode
       double precision          , intent(in   )          :: mass
       type            (treeNode), intent(inout), pointer :: thisNode
     end function Map_Mass
  end interface

  ! Interface for mass error functions.
  abstract interface
     double precision function Mass_Error(mass,thisNode)
       import treeNode
       double precision          , intent(in   )          :: mass
       type            (treeNode), intent(inout), pointer :: thisNode
     end function Mass_Error
  end interface

  ! Type for descriptors of correlation functions.
  type :: correlationFunctionDescriptor
     double precision                                :: redshift
     double precision                                :: massSystematicLogM0
     double precision                                :: massRandomError
     procedure       (  Mass_Error), pointer, nopass :: massRandomErrorFunction
     double precision                                :: massLogarithmicMinimum
     integer                                         :: massSystematicCoefficientCount
     integer                                         :: massType
     double precision                                :: massUnitsInSI
     character       (len= 32     )                  :: label
     character       (len=128     )                  :: comment
     procedure       (Map_Mass    ), pointer, nopass :: mapMass
  end type correlationFunctionDescriptor

  ! Correlation function descriptors.
  type(correlationFunctionDescriptor), dimension(correlationFunctionsSupportedCount), target :: correlationFunctionDescriptors= &
       & [                                                                                                                      &
       ! Hearin et al. (2013) SDSS.
       &                           correlationFunctionDescriptor(                                                               &
       &                                                          0.070d0+0                                              ,      &
       &                                                         11.0000d+0                                              ,      &
       &                                                          0.07000+0                                              ,      &
       &                                                          null()                                                 ,      &
       &                                                          8.000d0                                                ,      &
       &                                                          2                                                      ,      &
       &                                                          massTypeStellar                                        ,      &
       &                                                          massSolar                                              ,      &
       &                                                          'sdssClusteringZ0.07'                                  ,      &
       &                                                          'SDSS galaxy clustering at z=0.07'                     ,      &
       &                                                          null()                                                        &
       &                                                        )                                                               &
       & ]

  ! Type to store size functions.
  type :: correlationFunction
     ! Copy of the mass function descriptor for this mass function.
     type            (correlationFunctionDescriptor), pointer                     :: descriptor
     ! Parameters for the systematic error model.
     double precision                               , allocatable, dimension(:  ) :: massSystematicCoefficients
     ! The index of the output corresponding to the required redshift.
     integer                                                                      :: outputNumber
     ! Mass range.
     double precision                               , allocatable, dimension(:  ) :: massMinimumLogarithmic
     ! Separations.
     double precision                               , allocatable, dimension(:  ) :: separation
     ! Line-of-sight integration depth.
     double precision                                                             :: lineOfSightDepth
     ! Indices of current halo.
     integer         (kind=kind_int8               )                              :: treeIndex                 , haloIndex
     ! Population statistics.
     double precision                               , allocatable, dimension(:  ) :: meanDensity
     ! Power spectrum wavenumbers.
     double precision                               , allocatable, dimension(:  ) :: wavenumber
     ! Power spectra.
     double precision                               , allocatable, dimension(:,:) :: oneHaloTerm               , twoHaloTerm
     ! Covariances.
     double precision                               , allocatable, dimension(:,:) :: oneHaloCovariance         , twoHaloCovariance
     ! Cosmology conversion factors.
     double precision                                                             :: cosmologyConversionMass   , cosmologyConversionSize
  end type correlationFunction

  ! Correlation functions.
  type(correlationFunction), allocatable, dimension(:) :: correlationFunctions

  ! Type for storing temporary size functions during cumulation.
  type :: correlationFunctionWork
     double precision, allocatable, dimension(:,:) :: satelliteProbability
     double precision, allocatable, dimension(  :) :: centralProbability  , fourierProfile
     integer         , allocatable, dimension(:  ) :: satelliteCount
     double precision                              :: haloBias            , haloWeight    , &
          &                                           haloTime
     logical                                       :: propertiesSet
  end type correlationFunctionWork

  ! Work array.
  type(correlationFunctionWork), allocatable, dimension(:) :: thisHalo
  !$omp threadprivate(thisHalo)

  abstract interface
     double precision function integrandTemplate(x)
       double precision, intent(in   ) :: x
     end function integrandTemplate
  end interface

contains

  !# <mergerTreeAnalysisTask>
  !#  <unitName>Galacticus_Output_Analysis_Correlation_Functions</unitName>
  !# </mergerTreeAnalysisTask>
  subroutine Galacticus_Output_Analysis_Correlation_Functions(thisTree,thisNode,iOutput,mergerTreeAnalyses)
    !% Construct correlation functions to compare to various observational determinations.
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
    use Dark_Matter_Halo_Scales
    use String_Handling
    use Galacticus_Output_Analyses_Cosmology_Scalings
    use Numerical_Comparison
    use Numerical_Ranges
    use FoX_dom
    use IO_XML
    implicit none
    type            (mergerTree                    ), intent(in   )                 :: thisTree
    type            (treeNode                      ), intent(inout), pointer        :: thisNode
    integer                                         , intent(in   )                 :: iOutput
    type            (varying_string                ), intent(in   ), dimension(:  ) :: mergerTreeAnalyses
    class           (cosmologyFunctionsClass       )               , pointer        :: cosmologyFunctionsModel
    double precision                                , allocatable  , dimension(:  ) :: separationTmp
    type            (node                          ), pointer                       :: doc,columnElement,massElement,hubbleElement,datum,cosmology,omegaDarkEnergyElement,omegaMatterElement,datasetElement,cosmologyScalingElement,separationElement,correlationElement,integration,lineOfSightDepth
    type            (cosmologyFunctionsMatterLambda)                                :: cosmologyFunctionsObserved
    type            (cosmologyParametersSimple     )                                :: cosmologyParametersObserved
    double precision                                                                :: mass, massLogarithmic, dataHubbleParameter, dataOmegaMatter, dataOmegaDarkEnergy,wavenumberMinimum,wavenumberMaximum
    integer                                                                         :: i,j,k, currentAnalysis, activeAnalysisCount, wavenumberCount,ioErr
    type            (varying_string                )                                :: parameterName, cosmologyScalingMass, cosmologyScalingSize
    character       (len=128                       )                                :: cosmologyScaling

    ! Initialize the module if necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Galacticus_Output_Analysis_Correlation_Functions_Initialize)
       if (.not.moduleInitialized) then
          ! Establish mapping functions for correlation function descriptors.
          correlationFunctionDescriptors(1)%mapMass => null()
          ! Determine how many supported mass functions are requested.
          activeAnalysisCount=0
          do i=1,correlationFunctionsSupportedCount
             if (any(trim(mergerTreeAnalyses) == trim(correlationFunctionLabels(i)))) activeAnalysisCount=activeAnalysisCount+1
          end do
          ! Allocate mass function arrays and populate with required data.
          if (activeAnalysisCount <= 0) then
             analysisActive=.false.
          else
             analysisActive=.true.
             currentAnalysis=0
             allocate(correlationFunctions(activeAnalysisCount))
             cosmologyFunctionsModel => cosmologyFunctions()
             do i=1,size(mergerTreeAnalyses)
                do j=1,correlationFunctionsSupportedCount
                   if (mergerTreeAnalyses(i) == trim(correlationFunctionLabels(j))) then
                      currentAnalysis=currentAnalysis+1
                      ! Set a pointer to the descriptor for this size function.
                      correlationFunctions(currentAnalysis)%descriptor => correlationFunctionDescriptors(j)
                      ! Read parameters of the systematic error model.
                      if (correlationFunctionDescriptors(j)%massSystematicCoefficientCount > 0) then
                         allocate(correlationFunctions(currentAnalysis)%massSystematicCoefficients(correlationFunctionDescriptors(j)%massSystematicCoefficientCount))
                         do k=1,correlationFunctionDescriptors(j)%massSystematicCoefficientCount
                            parameterName=trim(correlationFunctionLabels(j))//'MassSystematic'
                            parameterName=parameterName//(k-1)
                            call Get_Input_Parameter(char(parameterName),correlationFunctions(currentAnalysis)%massSystematicCoefficients(k),defaultValue=0.0d0)
                         end do
                      end if
                      ! Find which output number corresponds to the required redshift.
                      correlationFunctions(currentAnalysis)%outputNumber=-1
                      do k=1,Galacticus_Output_Time_Count()
                         if     (                                                                                             &
                              &  Values_Agree(                                                                                &
                              &               correlationFunctionDescriptors(j)%redshift                                    , &
                              &               cosmologyFunctionsModel%redshiftFromExpansionFactor(                            &
                              &               cosmologyFunctionsModel%expansionFactor             (                           &
                              &                                                                    Galacticus_Output_Time(k)  &
                              &                                                                   )                           &
                              &                                                                  )                          , &
                              &               absTol=0.001d0                                                                  &
                              &              )                                                                                &
                              & ) then
                            correlationFunctions(currentAnalysis)%outputNumber=k
                            exit
                         end if
                      end do
                      if (correlationFunctions(currentAnalysis)%outputNumber < 0)                                             &
                           & call Galacticus_Error_Report(                                                                    &
                           &                                'Galacticus_Output_Analysis_Correlation_Functions'              , &
                           &                                'unable to find required redshift in outputs for mass function '  &
                           &                              //trim(correlationFunctionLabels(j))                                &
                           &                             )
                      ! Initialize tree/halo indices.
                      correlationFunctions(currentAnalysis)%treeIndex=-1_kind_int8
                      correlationFunctions(currentAnalysis)%haloIndex=-1_kind_int8
                      ! Read the appropriate observational data definition.
                      select case (trim(correlationFunctionLabels(j)))
                      case ('sdssClusteringZ0.07')
                         ! Read data for the Hearin et al. (2013) projected correlation function.
                         !$omp critical (FoX_DOM_Access)
                         ! Parse document.
                         doc => parseFile(char(Galacticus_Input_Path())//"data/observations/correlationFunctions/Projected_Correlation_Functions_Hearin_2013.xml",iostat=ioErr)
                         if (ioErr /= 0) call Galacticus_Error_Report('Galacticus_Output_Analysis_Correlation_Functions','Unable to find data file')
                         ! Extract cosmological parameters.
                         cosmology              => XML_Get_First_Element_By_Tag_Name(doc      ,"cosmology"      )
                         hubbleElement          => XML_Get_First_Element_By_Tag_Name(cosmology,"hubble"         )
                         omegaMatterElement     => XML_Get_First_Element_By_Tag_Name(cosmology,"omegaMatter"    )
                         omegaDarkEnergyElement => XML_Get_First_Element_By_Tag_Name(cosmology,"omegaDarkEnergy")
                         call extractDataContent(hubbleElement         ,dataHubbleParameter)
                         call extractDataContent(omegaMatterElement    ,dataOmegaMatter    )
                         call extractDataContent(omegaDarkEnergyElement,dataOmegaDarkEnergy)
                         ! Extract integration depth.
                         integration      => XML_Get_First_Element_By_Tag_Name(doc        ,"integration"     )
                         lineOfSightDepth => XML_Get_First_Element_By_Tag_Name(integration,"lineOfSightDepth")
                         call extractDataContent(lineOfSightDepth,correlationFunctions(currentAnalysis)%lineOfSightDepth)
                         ! Allocate arrays.
                         call Alloc_Array(correlationFunctions(currentAnalysis)%massMinimumLogarithmic,[3])
                         ! Locate datasets.
                         datasetElement => item(getElementsByTagname(doc,"correlationFunction"),0)
                         do k=1,3
                            columnElement      => item(getElementsByTagname(datasetElement,"columns"            ),k-1)
                            massElement        => item(getElementsByTagname( columnElement,"mass"               ),0)
                            separationElement  => item(getElementsByTagname( columnElement,"separation"         ),0)
                            correlationElement => item(getElementsByTagname( columnElement,"correlationFunction"),0)
                            ! Extract cosmology scalings.
                            cosmologyScalingElement => XML_Get_First_Element_By_Tag_Name(massElement,"cosmologyScaling")
                            call extractDataContent(cosmologyScalingElement,cosmologyScaling)
                            if (k == 1) then
                               cosmologyScalingMass = trim(cosmologyScaling)
                            else
                               if (cosmologyScalingMass /= trim(cosmologyScaling))                                     &
                                    & call Galacticus_Error_Report(                                                    &
                                    &                              'Galacticus_Output_Analysis_Correlation_Functions', &
                                    &                              'mass cosmology scaling mistmatch'                  &
                                    &                             )
                            end if
                            cosmologyScalingElement => XML_Get_First_Element_By_Tag_Name(separationElement,"cosmologyScaling")
                            call extractDataContent(cosmologyScalingElement,cosmologyScaling)
                            if (k == 1) then
                               cosmologyScalingSize = trim(cosmologyScaling)
                            else
                               if (cosmologyScalingSize /= trim(cosmologyScaling))                                     &
                                    & call Galacticus_Error_Report(                                                    &
                                    &                              'Galacticus_Output_Analysis_Correlation_Functions', &
                                    &                              'size cosmology scaling mistmatch'                  &
                                    &                             )
                            end if
                            ! Extract minimum mass.
                            datum => XML_Get_First_Element_By_Tag_Name(massElement,"minimum")
                            call extractDataContent(datum,correlationFunctions(currentAnalysis)%massMinimumLogarithmic(k))
                            ! Extract separations.
                            call XML_Array_Read(separationElement,"datum",separationTmp)
                            if (k == 1) then
                               correlationFunctions(currentAnalysis)%separation=separationTmp
                            else
                               if (any(correlationFunctions(currentAnalysis)%separation /= separationTmp))             &
                                    & call Galacticus_Error_Report(                                                    &
                                    &                              'Galacticus_Output_Analysis_Correlation_Functions', &
                                    &                              'separation mistmatch'                              &
                                    &                             )
                            end if
                         end do
                         ! Convert separations from logarithmic to linear.
                         correlationFunctions(currentAnalysis)%separation=10.0d0**correlationFunctions(currentAnalysis)%separation
                         ! Destroy the document.
                         call destroy(doc)
                         !$omp end critical (FoX_DOM_Access)                         
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
                         ! Allocate wavenumbers.
                         wavenumberCount  =1000
                         wavenumberMinimum=0.001d0
                         wavenumberMaximum=10000.0d0
                         call Alloc_Array(correlationFunctions(currentAnalysis)%wavenumber       ,[   wavenumberCount                         ])
                         call Alloc_Array(correlationFunctions(currentAnalysis)%meanDensity      ,[                      3                    ])
                         call Alloc_Array(correlationFunctions(currentAnalysis)%oneHaloTerm      ,[   wavenumberCount   ,3                    ])
                         call Alloc_Array(correlationFunctions(currentAnalysis)%twoHaloTerm      ,[   wavenumberCount   ,3                    ])
                         call Alloc_Array(correlationFunctions(currentAnalysis)%oneHaloCovariance,[3*(wavenumberCount+1),3*(wavenumberCount+1)])
                         call Alloc_Array(correlationFunctions(currentAnalysis)%twoHaloCovariance,[3*(wavenumberCount+1),3*(wavenumberCount+1)])
                         correlationFunctions(currentAnalysis)%wavenumber=Make_Range(wavenumberMinimum,wavenumberMaximum,wavenumberCount,rangeTypeLogarithmic)
                      case default
                         call Galacticus_Error_Report('Galacticus_Output_Analysis_Correlation_Functions','unknown size function')
                      end select
                      ! Get cosmological conversion factors.
                      call Cosmology_Conversion_Factors(                                                                                               &
                           &                            correlationFunctions(currentAnalysis)%descriptor%redshift                                    , &
                           &                            cosmologyFunctionsModel                                                                      , &
                           &                            cosmologyFunctionsObserved                                                                   , &
                           &                            cosmologyScalingMass           =cosmologyScalingMass                                         , &
                           &                            cosmologyScalingSize           =cosmologyScalingSize                                         , &
                           &                            cosmologyConversionMass        =correlationFunctions(currentAnalysis)%cosmologyConversionMass, &
                           &                            cosmologyConversionSize        =correlationFunctions(currentAnalysis)%cosmologyConversionSize  &
                           &                           )
                      ! Initialize population statistics.
                      correlationFunctions(currentAnalysis)%meanDensity      =0.0d0
                      correlationFunctions(currentAnalysis)%oneHaloTerm      =0.0d0
                      correlationFunctions(currentAnalysis)%twoHaloTerm      =0.0d0
                      correlationFunctions(currentAnalysis)%oneHaloCovariance=0.0d0
                      correlationFunctions(currentAnalysis)%twoHaloCovariance=0.0d0
                      exit
                  end if
                end do
             end do
          end if
          ! Record that module is initialized.
          moduleInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Analysis_Correlation_Functions_Initialize)
    end if
    ! Return if this analysis is not active.
    if (.not.analysisActive) return
    ! Allocate work arrays.
    if (.not.allocated(thisHalo)) allocate(thisHalo(size(correlationFunctions)))
    ! Iterate over active analyses.
    do i=1,size(correlationFunctions)
       ! Skip if this analysis is not active, or if this is not the correct output.
       if (iOutput /= correlationFunctions(i)%outputNumber) cycle
       ! Get the galactic mass.
       mass=                                                                                                                                                   &
            &  Galactic_Structure_Enclosed_Mass(thisNode,radiusLarge,componentType=componentTypeDisk    ,massType=correlationFunctions(i)%descriptor%massType) &
            & +Galactic_Structure_Enclosed_Mass(thisNode,radiusLarge,componentType=componentTypeSpheroid,massType=correlationFunctions(i)%descriptor%massType)
       if (mass <= 0.0d0) return
       if (associated(correlationFunctions(i)%descriptor%mapMass)) mass=correlationFunctions(i)%descriptor%mapMass(mass,thisNode)
       mass=mass*correlationFunctions(i)%cosmologyConversionMass ! Convert for cosmology.
       massLogarithmic=log10(mass)
       do j=1,correlationFunctions(i)%descriptor%massSystematicCoefficientCount
          massLogarithmic=+massLogarithmic                                          &
               &          +correlationFunctions(i)%massSystematicCoefficients(j)    &
               &          *(                                                        &
               &            +log10(mass)                                            &
               &            -correlationFunctions(i)%descriptor%massSystematicLogM0 &
               &          )**(j-1)
       end do
       if (massLogarithmic < correlationFunctions(i)%descriptor%massLogarithmicMinimum) return
       ! Accumulate the halo.
       call Accumulate_Node(correlationFunctions(i),thisHalo(i),thisTree,thisNode,massLogarithmic)
    end do
    return
  end subroutine Galacticus_Output_Analysis_Correlation_Functions

  subroutine Accumulate_Node(thisCorrelationFunction,thisHalo,thisTree,thisNode,massLogarithmic)
    !% Accumulate a single galaxy to the population of the current halo. Since galaxy masses
    !% have random errors, each galaxy added is assigned an inclusion probability, which will be
    !% taken into account when evaluating the one- and two-halo terms from this halo in the halo
    !% model.
    use Dark_Matter_Halo_Biases
    use Cosmology_Functions
    use Dark_Matter_Profiles
    use Memory_Management
    implicit none
    type            (correlationFunction    ), intent(inout)                 :: thisCorrelationFunction
    type            (correlationFunctionWork), intent(inout)                 :: thisHalo
    type            (mergerTree             ), intent(in   )                 :: thisTree
    type            (treeNode               ), intent(inout), pointer        :: thisNode
    double precision                         , intent(in   )                 :: massLogarithmic
    type            (treeNode               )               , pointer        :: hostNode
    class           (cosmologyFunctionsClass)               , pointer        :: cosmologyFunctions_    
    class           (nodeComponentBasic     )               , pointer        :: thisBasic
    double precision                         , allocatable  , dimension(:,:) :: satelliteProbabilityTmp
    integer                                  , parameter                     :: satelliteCountMinimum=100
    integer         (kind=kind_int8        )                                 :: hostIndex    
    integer                                                                  :: i, j
    double precision                                                         :: expansionFactor, galaxyInclusionProbability, &
         &                                                                      randomError    , mass
    
    ! Get the index of the host halo.
    if (thisNode%isSatellite()) then
       hostNode => thisNode%parent
    else
       hostNode => thisNode
    end if
    hostIndex=hostNode%uniqueID()
    ! Allocate arrays.
    if (.not.allocated(thisHalo%satelliteCount)) then
       call Alloc_Array(thisHalo%satelliteCount,[size(thisCorrelationFunction%massMinimumLogarithmic)])
       call Alloc_Array(thisHalo%centralProbability,[size(thisCorrelationFunction%massMinimumLogarithmic)])
       call Alloc_Array(thisHalo%fourierProfile,[size(thisCorrelationFunction%wavenumber)])
       thisHalo%propertiesSet=.false.
       thisHalo%satelliteCount=0
       thisHalo%centralProbability=0.0d0
    end if
    ! Check if the host has changed.
    if (thisTree%index /= thisCorrelationFunction%treeIndex .or. hostIndex /= thisCorrelationFunction%haloIndex) &
         & call Accumulate_Halo(thisCorrelationFunction,thisHalo)
    ! Accumulate properties to the current halo.
    thisCorrelationFunction%treeIndex=thisTree%index
    thisCorrelationFunction%haloIndex=hostIndex
    ! Find the random error on the galaxy mass.
    if (associated(thisCorrelationFunction%descriptor%massRandomErrorFunction)) then
       mass=10.0d0**massLogarithmic
       randomError=thisCorrelationFunction%descriptor%massRandomErrorFunction(mass,thisNode)
    else
       randomError=thisCorrelationFunction%descriptor%massRandomError
    end if
    ! Iterate over mass ranges.
    do j=1,size(thisCorrelationFunction%massMinimumLogarithmic)
       ! Find the probability that this galaxy is included in the sample.
       galaxyInclusionProbability=0.5d0*(1.0d0-erf((thisCorrelationFunction%massMinimumLogarithmic(j)-massLogarithmic)/randomError/sqrt(2.0d0)))
       if (thisNode%isSatellite()) then
          if (galaxyInclusionProbability > 0.0d0) then
             thisHalo%satelliteCount(j)=thisHalo%satelliteCount(j)+1
             if (.not.allocated(thisHalo%satelliteProbability)) then
                call Alloc_Array(thisHalo%satelliteProbability,[satelliteCountMinimum,size(thisCorrelationFunction%massMinimumLogarithmic)])
             else if (size(thisHalo%satelliteProbability,dim=1) < thisHalo%satelliteCount(j)) then
                call Move_Alloc(thisHalo%satelliteProbability,satelliteProbabilityTmp)
                call Alloc_Array(thisHalo%satelliteProbability,[2*size(satelliteProbabilityTmp),size(thisCorrelationFunction%massMinimumLogarithmic)])
                thisHalo%satelliteProbability(1:size(satelliteProbabilityTmp),:)=satelliteProbabilityTmp
                call Dealloc_Array(satelliteProbabilityTmp)
             end if
             thisHalo%satelliteProbability(thisHalo%satelliteCount(j),j)=galaxyInclusionProbability
          end if
       else
          thisHalo%centralProbability(j)=  galaxyInclusionProbability
       end if
       if (galaxyInclusionProbability > 0.0d0 .and. .not.thisHalo%propertiesSet) then
          thisHalo%propertiesSet=.true.
          thisBasic                     => hostNode%basic()
          thisHalo%haloBias             =  Dark_Matter_Halo_Bias(hostNode)
          thisHalo%haloWeight           =  thisTree%volumeWeight
          thisHalo%haloTime             =  thisBasic%time()
          cosmologyFunctions_                          => cosmologyFunctions()
          expansionFactor                              =  cosmologyFunctions_%expansionFactor(thisBasic%time())
          do i=1,size(thisCorrelationFunction%wavenumber)
             ! Note that wavenumbers must be converted from comoving to physical units for the dark matter profile k-space function.
             thisHalo%fourierProfile(i)=Dark_Matter_Profile_kSpace(hostNode,thisCorrelationFunction%waveNumber(i)/expansionFactor)
          end do
       end if
    end do
    return
  end subroutine Accumulate_Node

  subroutine Accumulate_Halo(thisCorrelationFunction,thisHalo)
    !% Assumulate a single halo's contributions to the halo model one- and two-halo terms. For
    !% the one-halo term we count contributions from central-satellite pairs, and from
    !% satellite-satellite pairs. Contributions differ in the scalings applied to the
    !% Fourier-transformed dark matter halo density profile---see
    !% \cite[][\S6.1]{cooray_halo_2002} for a discussion of this. The number of satellites in
    !% the halo is assumed to follow a Poisson binomial distribution.
    use Math_Distributions_Poisson_Binomial
    use Vectors
    implicit none
    type            (correlationFunction    ), intent(inout)                                                   :: thisCorrelationFunction
    type            (correlationFunctionWork), intent(inout)                                                   :: thisHalo
    double precision                         , dimension(                                                                                                      &
         &                                               size(thisCorrelationFunction%wavenumber            ),                                                 &
         &                                               size(thisCorrelationFunction%massMinimumLogarithmic)                                                  &
         &                                              )                                                      :: thisOneHaloTerm        , thisTwoHaloTerm
    double precision                         , dimension(size(thisCorrelationFunction%massMinimumLogarithmic)) :: galaxyDensity
    logical                                  , dimension(size(thisCorrelationFunction%massMinimumLogarithmic)) :: oneHaloTermActive, twoHaloTermActive
    double precision                                                                                           :: satellitePairsCountMean, satelliteCountMean
    integer                                                                                                    :: wavenumberCount, i, j, indexBeginI, indexEndI, indexBeginJ, indexEndJ

    ! Return immediately if no nodes have been accumulated.
    if (thisCorrelationFunction%treeIndex /= -1_kind_int8) then
       oneHaloTermActive=.false.
       twoHaloTermActive=.false.
       ! Iterate over mass bins.
       do i=1,size(thisCorrelationFunction%massMinimumLogarithmic)
          ! Accumulate to one- and two-halo terms.
          if (thisHalo%satelliteCount(i) > 0) then
             satelliteCountMean     =Poisson_Binomial_Distribution_Mean      (thisHalo%satelliteProbability(1:thisHalo%satelliteCount(i),i))
             satellitePairsCountMean=Poisson_Binomial_Distribution_Mean_Pairs(thisHalo%satelliteProbability(1:thisHalo%satelliteCount(i),i))
          else
             satelliteCountMean     =0.0d0
             satellitePairsCountMean=0.0d0
          end if
          if (thisHalo%centralProbability(i) > 0.0d0 .or. satelliteCountMean > 0.0d0) then
             ! Compute and accumulate contribution to galaxy density.
             galaxyDensity                      (i)=   thisHalo%haloWeight                    &
                  &                                 *(                                        &
                  &                                   +thisHalo%centralProbability        (i) &
                  &                                   +satelliteCountMean                     &
                  &                                  )
             !$omp atomic
             thisCorrelationFunction%meanDensity(i)=   thisCorrelationFunction%meanDensity(i) &
                  &                                 +  galaxyDensity                      (i)
             ! Compute and accumulate one-halo term.
             if (satelliteCountMean > 0.0d0) then
                oneHaloTermActive(  i)=.true.
                thisOneHaloTerm  (:,i)= thisHalo%haloWeight              &
                     &                 *(                                &
                     &                   +thisHalo%centralProbability(i) &
                     &                   *satelliteCountMean             &
                     &                   *thisHalo%fourierProfile        &
                     &                   +satellitePairsCountMean        &
                     &                   *thisHalo%fourierProfile    **2 &
                     &                  )
                !$omp critical(Analyses_Correlation_Functions_Accumulate1)
                thisCorrelationFunction%oneHaloTerm(:,i)= thisCorrelationFunction%oneHaloTerm(:,i) &
                     &                                   +thisOneHaloTerm                    (:,i)
                !$omp end critical(Analyses_Correlation_Functions_Accumulate1)
             end if
             ! Compute and accumulate two-halo term.
             twoHaloTermActive(  i)=.true.
             thisTwoHaloTerm  (:,i)= galaxyDensity(i)          &
                  &                 *  thisHalo%haloBias       &
                  &                 *  thisHalo%fourierProfile
             !$omp critical(Analyses_Correlation_Functions_Accumulate2)
             thisCorrelationFunction%twoHaloTerm(:,i)= thisCorrelationFunction%twoHaloTerm(:,i) &
                  &                                   +thisTwoHaloTerm                    (:,i)
             !$omp end critical(Analyses_Correlation_Functions_Accumulate2)
          end if
       end do
       ! Accumulate contributions to covariance matrices. Each covariance matrix (for one- and two-halo terms) includes the power
       ! spectrum contribution, plus an extra columns for galaxy density (since the power spectra and galaxy densities are
       ! correlated).
       wavenumberCount=size(thisHalo%fourierProfile)
       do i=1,size(thisCorrelationFunction%massMinimumLogarithmic)
          indexBeginI=1+(i-1)*(wavenumberCount+1)
          indexEndI  =indexBeginI+wavenumberCount-1
          do j=1,size(thisCorrelationFunction%massMinimumLogarithmic)
             indexBeginJ=1+(j-1)*(wavenumberCount+1)
             indexEndJ  =indexBeginJ+wavenumberCount-1
             if (twoHaloTermActive(i) .and. twoHaloTermActive(j)) then
                !$omp critical(Analyses_Correlation_Functions_Accumulate_Covariance2)
                thisCorrelationFunction        %twoHaloCovariance(indexBeginI:indexEndI,indexBeginJ:indexEndJ) &
                     & =thisCorrelationFunction%twoHaloCovariance(indexBeginI:indexEndI,indexBeginJ:indexEndJ) &
                     & +Vector_Outer_Product(                                                                  &
                     &                       thisTwoHaloTerm(:,i),                                             &
                     &                       thisTwoHaloTerm(:,j)                                              &
                     &                      )
                thisCorrelationFunction        %twoHaloCovariance(1+indexEndI,indexBeginJ:indexEndJ) &
                     & =thisCorrelationFunction%twoHaloCovariance(1+indexEndI,indexBeginJ:indexEndJ) &
                     & +galaxyDensity(i)                                                             &
                     & *                     thisTwoHaloTerm(:,j)
                thisCorrelationFunction        %twoHaloCovariance(indexBeginI:indexEndI,1+indexEndJ) &
                     & =thisCorrelationFunction%twoHaloCovariance(indexBeginI:indexEndI,1+indexEndJ) &
                     & +galaxyDensity(j)                                                             &
                     &*                      thisTwoHaloTerm(:,i)
                thisCorrelationFunction        %twoHaloCovariance(1+indexEndI,1+indexEndJ) &
                     & =thisCorrelationFunction%twoHaloCovariance(1+indexEndI,1+indexEndJ) &
                     & +galaxyDensity(i)                                                   &
                     & *galaxyDensity(j)
                !$omp end critical(Analyses_Correlation_Functions_Accumulate_Covariance2)
             end if
             !$omp critical(Analyses_Correlation_Functions_Accumulate_Covariance1)
             if (oneHaloTermActive(i) .and. oneHaloTermActive(j))                                            &
                  & thisCorrelationFunction  %oneHaloCovariance(indexBeginI:indexEndI,indexBeginJ:indexEndJ) &
                  &  =thisCorrelationFunction%oneHaloCovariance(indexBeginI:indexEndI,indexBeginJ:indexEndJ) &
                  &  +Vector_Outer_Product(                                                                  &
                  &                        thisOneHaloTerm(:,i),                                             &
                  &                        thisOneHaloTerm(:,j)                                              &
                  &                       )
             if (twoHaloTermActive(i) .and. oneHaloTermActive(j))                                  &
                  & thisCorrelationFunction  %oneHaloCovariance(1+indexEndI,indexBeginJ:indexEndJ) &
                  &  =thisCorrelationFunction%oneHaloCovariance(1+indexEndI,indexBeginJ:indexEndJ) &
                  &  +galaxyDensity(i)**2                                                          &
                  &  *                     thisOneHaloTerm(:,j)
             if (oneHaloTermActive(i) .and. twoHaloTermActive(j))                                  &
                  & thisCorrelationFunction  %oneHaloCovariance(indexBeginI:indexEndI,1+indexEndJ) &
                  &  =thisCorrelationFunction%oneHaloCovariance(indexBeginI:indexEndI,1+indexEndJ) &
                  &  +galaxyDensity(j)**2                                                          &
                  &  *                     thisOneHaloTerm(:,i)
             if (twoHaloTermActive(i) .and. twoHaloTermActive(j))                        &
                  & thisCorrelationFunction  %oneHaloCovariance(1+indexEndI,1+indexEndJ) &
                  &  =thisCorrelationFunction%oneHaloCovariance(1+indexEndI,1+indexEndJ) &
                  &  +galaxyDensity(i)**2                                                &
                  &  *galaxyDensity(j)**2
             !$omp end critical(Analyses_Correlation_Functions_Accumulate_Covariance1)
          end do
       end do
    end if
    ! Reset counts.
    if (allocated(thisHalo%centralProbability)) thisHalo%centralProbability=0.0d0
    if (allocated(thisHalo%satelliteCount    )) thisHalo%satelliteCount    =0
    thisHalo%propertiesSet=.false.
    ! Reset indices.
    thisCorrelationFunction%treeIndex=-1_kind_int8
    thisCorrelationFunction%haloIndex=-1_kind_int8
    return
  end subroutine Accumulate_Halo
  
  !# <hdfPreCloseTask>
  !#  <unitName>Galacticus_Output_Analysis_Correlation_Functions_Output</unitName>
  !# </hdfPreCloseTask>
  subroutine Galacticus_Output_Analysis_Correlation_Functions_Output
    !% Outputs SDSS $z\approx 0.07$ stellar mass function to file.
    use, intrinsic :: ISO_C_Binding
    use Galacticus_HDF5
    use Power_Spectra
    use Linear_Growth
    use Array_Utilities
    use FFTLogs
    use Memory_Management
    use Tables
    use Linear_Algebra
    implicit none
    double precision                            , allocatable, dimension(:  ) :: separation
    double precision                            , allocatable, dimension(:,:) :: powerSpectrumCovariance       , jacobian                 , &
         &                                                                       oneHaloCovariance             , twoHaloCovariance        , &
         &                                                                       correlationCovariance         , covarianceTmp            , &
         &                                                                       projectedCorrelationCovariance, binnedProjectedCorrelationCovariance, &
         &                                                                       powerSpectrum, correlation, projectedCorrelation, binnedProjectedCorrelation
    integer                                                                   :: i                             , k                        , &
         &                                                                       j                             , wavenumberCount, m, n, massCount
    type            (hdf5Object                )                              :: analysisGroup                 , correlationFunctionGroup , &
         &                                                                       thisDataset
    type            (table1DLogarithmicLinear  )                              :: correlationTable
    double precision                                                          :: projectedSeparation           , binSeparationMinimum     , &
         &                                                                       binSeparationMaximum          , binWidthLogarithmic      , &
         &                                                                       linearGrowthFactorSquared
    type            (matrix                    )                              :: jacobianMatrix                , covarianceMatrix
    procedure       (integrandTemplate         ), pointer                     :: integrandWeightFunction

    ! Return immediately if this analysis is not active.
    if (.not.analysisActive) return
    ! Iterate over mass functions.
    do k=1,size(correlationFunctions)
       ! Accumulate any final halo.
       call Accumulate_Halo(correlationFunctions(k),thisHalo(k))       
       ! Normalize one- and two-halo terms.
       massCount      =size(correlationFunctions(k)%massMinimumLogarithmic)
       wavenumberCount=size(correlationFunctions(k)%wavenumber            )
       call Alloc_Array(jacobian         ,[3*wavenumberCount,3*(wavenumberCount+1)])
       call Alloc_Array(oneHaloCovariance,[3*wavenumberCount,3* wavenumberCount   ])
       call Alloc_Array(twoHaloCovariance,[3*wavenumberCount,3* wavenumberCount   ])
       ! One-halo term.
       jacobian=0.0d0
       do n=1,massCount
          do i=1,wavenumberCount
             jacobian((n-1)*wavenumberCount+i,(n-1)*wavenumberCount+i)=1.0d0/correlationFunctions(k)%meanDensity(n)**2
          end do
          jacobian((n-1)*wavenumberCount+1:n*wavenumberCount,n*wavenumberCount+1)=-2.0d0*correlationFunctions(k)%oneHaloTerm(i,n)/correlationFunctions(k)%meanDensity(n)**3
       end do
       jacobianMatrix                     =jacobian
       covarianceMatrix                   =correlationFunctions(k)%oneHaloCovariance
       oneHaloCovariance                  =jacobianMatrix*(covarianceMatrix*jacobianMatrix%transpose())
       do n=1,massCount
          correlationFunctions(k)%oneHaloTerm(:,n)=correlationFunctions(k)%oneHaloTerm(:,n)/correlationFunctions(k)%meanDensity(n)**2
       end do
       ! Two-halo term.
       jacobian=0.0d0
       do n=1,massCount
          do i=1,wavenumberCount
             jacobian((n-1)*wavenumberCount+i,(n-1)*wavenumberCount+i)=1.0d0/correlationFunctions(k)%meanDensity(n)
          end do
          jacobian((n-1)*wavenumberCount+1:n*wavenumberCount,n*wavenumberCount+1)=-correlationFunctions(k)%twoHaloTerm(i,n)/correlationFunctions(k)%meanDensity(n)**2
       end do
       jacobianMatrix                     =jacobian
       covarianceMatrix                   =correlationFunctions(k)%twoHaloCovariance
       twoHaloCovariance                  =jacobianMatrix*(covarianceMatrix*jacobianMatrix%transpose())
       do n=1,massCount
          correlationFunctions(k)%twoHaloTerm(:,n)=correlationFunctions(k)%twoHaloTerm(:,n)/correlationFunctions(k)%meanDensity(n)
       end do
       call Dealloc_Array(jacobian)
       ! Square the two halo term.
       call Alloc_Array(jacobian,[3*wavenumberCount,3*wavenumberCount])
       jacobian=0.0d0
       do n=1,massCount
          do i=1,wavenumberCount
             jacobian((n-1)*wavenumberCount+i,(n-1)*wavenumberCount+i)=2.0d0*correlationFunctions(k)%twoHaloTerm(i,n)
          end do
       end do
       jacobianMatrix                        =jacobian
       covarianceMatrix                      =twoHaloCovariance
       twoHaloCovariance                     =jacobianMatrix*(covarianceMatrix*jacobianMatrix%transpose())
       do n=1,massCount
          correlationFunctions(k)%twoHaloTerm(:,n)=correlationFunctions(k)%twoHaloTerm(:,n)**2
       end do
       call Dealloc_Array(jacobian)
       ! Multiply two-halo term by the linear theory power spectrum.
       linearGrowthFactorSquared=Linear_Growth_Factor(thisHalo(k)%haloTime)**2
       do n=1,massCount
          do i=1,size(correlationFunctions(k)%wavenumber)
             correlationFunctions(k)%twoHaloTerm(i,n)=                      correlationFunctions(k)%twoHaloTerm(i,n)  &
                  &                                   *Power_Spectrum      (correlationFunctions(k)%wavenumber (i  )) &
                  &                                   *linearGrowthFactorSquared
             do m=1,massCount
                do j=1,size(correlationFunctions(k)%wavenumber)
                   twoHaloCovariance        ((n-1)*wavenumberCount+i,(m-1)*wavenumberCount+j)= &
                        & +twoHaloCovariance((n-1)*wavenumberCount+i,(m-1)*wavenumberCount+j)  &
                        & *Power_Spectrum   (correlationFunctions(k)%wavenumber (i  ))         &
                        & *Power_Spectrum   (correlationFunctions(k)%wavenumber (  j))         &
                        & *linearGrowthFactorSquared**2
                end do
             end do
          end do
       end do
       ! Construct the final power spectra.
       call Alloc_Array(powerSpectrum,[wavenumberCount,massCount])
       do n=1,massCount
          powerSpectrum(:,n)=correlationFunctions(k)%oneHaloTerm(:,n)+correlationFunctions(k)%twoHaloTerm(:,n)
       end do
       powerSpectrumCovariance=oneHaloCovariance+twoHaloCovariance
       call Dealloc_Array(oneHaloCovariance)
       call Dealloc_Array(twoHaloCovariance)
       ! Allocate correlation function and separation arrays.
       call Alloc_Array(correlation,shape(powerSpectrum))
       call Alloc_Array(separation ,[wavenumberCount])
       ! Fourier transform the power spectrum to get the correlation function.
       do n=1,massCount
          call FFTLog(                                    &
               &      correlationFunctions(k)%wavenumber, &
               &      separation                        , &
               &      +powerSpectrum(:,n)                 &
               &      *correlationFunctions(k)%wavenumber &
               &      * 4.0d0*Pi                          &
               &      /(2.0d0*Pi)**3                    , &
               &      correlation(:,n)                  , &
               &      fftLogSine                        , &
               &      fftLogForward                       &
               &     )
          correlation(:,n)=correlation(:,n)/separation
       end do
       ! Compute the covariance of the correlation function.
       call Alloc_Array(covarianceTmp        ,[3*wavenumberCount,3*wavenumberCount])
       call Alloc_Array(correlationCovariance,[3*wavenumberCount,3*wavenumberCount])
       ! Apply wavenumber weighting to the power spectrum covariance.
       do n=1,massCount
          do m=1,massCount
             do i=1,wavenumberCount
                do j=1,wavenumberCount
                   powerSpectrumCovariance                   ((n-1)*wavenumberCount+i,(m-1)*wavenumberCount+j) &
                        & =powerSpectrumCovariance           ((n-1)*wavenumberCount+i,(m-1)*wavenumberCount+j) &
                        & *correlationFunctions(k)%wavenumber(i  ) &
                        & *correlationFunctions(k)%wavenumber(  j) &
                        & *(                                       &
                        &   + 4.0d0*Pi                             &
                        &   /(2.0d0*Pi)**3                         &
                        &  )**2
                end do
             end do
          end do
       end do
       ! Derive the covariance of the correlation function by first Fourier transforming each row of the power spectrum covariance
       ! matrix, and then Fourier transforming each column.
       do n=1,massCount
          do m=1,massCount
             do i=1,wavenumberCount
                call FFTlog(                                                                                            &
                     &      correlationFunctions(k)%wavenumber                                                        , &
                     &      separation                                                                                , &
                     &      powerSpectrumCovariance((n-1)*wavenumberCount+i,(m-1)*wavenumberCount+1:m*wavenumberCount), &
                     &      covarianceTmp          ((n-1)*wavenumberCount+i,(m-1)*wavenumberCount+1:m*wavenumberCount), &
                     &      fftLogSine                                                                                , &
                     &      fftLogForward                                                                               &
                     )
             end do
          end do
       end do
       do n=1,massCount
          do m=1,massCount
             do i=1,wavenumberCount
                call FFTlog(                                                                                            &
                     &      correlationFunctions(k)%wavenumber                                                        , &
                     &      separation                                                                                , &
                     &      covarianceTmp          ((n-1)*wavenumberCount+1:n*wavenumberCount,(m-1)*wavenumberCount+i), &
                     &      correlationCovariance  ((n-1)*wavenumberCount+1:n*wavenumberCount,(m-1)*wavenumberCount+i), &
                     &      fftLogSine                                                                                , &
                     &      fftLogForward                                                                               &
                     )
             end do
          end do
       end do
       do n=1,massCount
          do m=1,massCount
             do i=1,wavenumberCount
                do j=1,wavenumberCount
                   correlationCovariance        ((n-1)*wavenumberCount+i,(m-1)*wavenumberCount+j) &
                        & =correlationCovariance((n-1)*wavenumberCount+i,(m-1)*wavenumberCount+j) &
                        & /separation           (                      i                        ) &
                        & /separation           (                                              j)
                end do
             end do
          end do
       end do
       call Dealloc_Array(covarianceTmp)       
       ! Scale separations for cosmology.
       separation=separation*correlationFunctions(k)%cosmologyConversionSize
       ! Construct correlation table.
       call correlationTable%create(separation(1),separation(wavenumberCount),size(separation),extrapolationTypeExtrapolate)
       ! Project the correlation function.
       call Alloc_Array(jacobian                      ,[massCount*wavenumberCount,massCount*wavenumberCount])
       call Alloc_Array(projectedCorrelationCovariance,[massCount*wavenumberCount,massCount*wavenumberCount])
       call Alloc_Array(projectedCorrelation          ,[wavenumberCount,massCount                ])
       jacobian=0.0d0
       integrandWeightFunction => projectionIntegrandWeight
       do i=1,wavenumberCount
          projectedSeparation=correlationTable%x(i)
          jacobian(i,1:wavenumberCount)=correlationTable%integrationWeights(                                                   &
               &                                                            projectedSeparation                              , &
               &                                                            sqrt(                                              &
               &                                                                 +projectedSeparation                     **2  &
               &                                                                 +correlationFunctions(k)%lineOfSightDepth**2  &
               &                                                                )                                            , &
               &                                                            integrandWeightFunction                            &
               &                                                           )
          do n=1,massCount
             if (n > 1) jacobian((n-1)*wavenumberCount+i,(n-1)*wavenumberCount+1:n*wavenumberCount)=jacobian(i,1:wavenumberCount)
             projectedCorrelation(i,n)=sum(jacobian(i,1:wavenumberCount)*correlation(:,n))
          end do
       end do
       jacobianMatrix                =jacobian
       covarianceMatrix              =correlationCovariance
       projectedCorrelationCovariance=jacobianMatrix*(covarianceMatrix*jacobianMatrix%transpose())
       call Dealloc_Array(jacobian)
       ! Integrate the projected correlation function over bins.
       call Alloc_Array(binnedProjectedCorrelation          ,[          size(correlationFunctions(k)%separation),massCount                                         ])
       call Alloc_Array(binnedProjectedCorrelationCovariance,[massCount*size(correlationFunctions(k)%separation),massCount*size(correlationFunctions(k)%separation)])
       call Alloc_Array(jacobian                            ,[massCount*size(correlationFunctions(k)%separation),massCount*wavenumberCount                         ])
       jacobian=0.0d0
       integrandWeightFunction => binningIntegrandWeight
       binWidthLogarithmic=log(correlationFunctions(k)%separation(2)/correlationFunctions(k)%separation(1))
       do i=1,size(correlationFunctions(k)%separation)      
          binSeparationMinimum         =correlationFunctions(k)%separation(i)*exp(-0.5d0*binWidthLogarithmic)
          binSeparationMaximum         =correlationFunctions(k)%separation(i)*exp(+0.5d0*binWidthLogarithmic)
          jacobian(i,1:wavenumberCount)=correlationTable%integrationWeights(         &
               &                                                            binSeparationMinimum   , &
               &                                                            binSeparationMaximum   , &
               &                                                            integrandWeightFunction  &
               &                                                           )                         &
               &                        /Pi                                                          &
               &                        /(                                                           &
               &                          +binSeparationMaximum**2                                   &
               &                          -binSeparationMinimum**2                                   &
               &                         )
          do n=1,massCount
             if (n > 1) jacobian((n-1)*size(correlationFunctions(k)%separation)+i,(n-1)*wavenumberCount+1:n*wavenumberCount)=jacobian(i,1:wavenumberCount)
             binnedProjectedCorrelation(i,n)=sum(jacobian(i,1:wavenumberCount)*projectedCorrelation(:,n))
          end do
       end do
       jacobianMatrix                      =jacobian
       covarianceMatrix                    =projectedCorrelationCovariance
       binnedProjectedCorrelationCovariance=jacobianMatrix*(covarianceMatrix*jacobianMatrix%transpose())
       call Dealloc_Array(jacobian)
       call correlationTable%destroy()
       ! Output the size function.
       !$omp critical(HDF5_Access)
       analysisGroup           =galacticusOutputFile%openGroup('analysis','Model analysis')
       correlationFunctionGroup=analysisGroup       %openGroup(trim(correlationFunctions(k)%descriptor%label),trim(correlationFunctions(k)%descriptor%comment))
       call correlationFunctionGroup%writeDataset  (correlationFunctions(k)%separation  ,'separation'                   ,'Separation'                      ,datasetReturned=thisDataset)
       call thisDataset             %writeAttribute(megaParsec                          ,'unitsInSI'                                                                                   )
       call thisDataset             %close         (                                                                                                                                   )
       call correlationFunctionGroup%writeDataset  (binnedProjectedCorrelation          ,'correlationFunction'          ,'Projected correlation'           ,datasetReturned=thisDataset)
       call thisDataset             %writeAttribute(megaParsec                          ,'unitsInSI'                                                                                   )
       call thisDataset             %close         (                                                                                                                                   )
       call correlationFunctionGroup%writeDataset  (binnedProjectedCorrelationCovariance,'correlationFunctionCovariance','Projected correlation covariance',datasetReturned=thisDataset)
       call thisDataset             %writeAttribute(megaParsec**2                       , 'unitsInSI'                                                                                  )
       call thisDataset             %close         (                                                                                                                                   )
       call correlationFunctionGroup%close         (                                                                                                                                   )
       call analysisGroup           %close         (                                                                                                                                   )
       !$omp end critical(HDF5_Access)
       call Dealloc_Array(binnedProjectedCorrelation          )
       call Dealloc_Array(binnedProjectedCorrelationCovariance)
    end do
    return
    
  contains

    double precision function projectionIntegrandWeight(separation)
      !% The weight function applied to the correlation function when integrating to get the projected correlation function.
      implicit none
      double precision, intent(in   ) :: separation
      
      if (separation > projectedSeparation) then
         projectionIntegrandWeight=2.0d0*separation/sqrt(separation**2-projectedSeparation**2)
      else
         projectionIntegrandWeight=0.0d0
      end if
      return
    end function projectionIntegrandWeight
    
    double precision function binningIntegrandWeight(separation)
      !% The weight function applied to the projected correlation function when integrating into bins.
      implicit none
      double precision, intent(in   ) :: separation

      binningIntegrandWeight=2.0d0*Pi*separation
      return
    end function binningIntegrandWeight

  end subroutine Galacticus_Output_Analysis_Correlation_Functions_Output

end module Galacticus_Output_Analyses_Correlation_Functions

