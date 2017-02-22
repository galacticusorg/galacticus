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

!% Contains a module which performs analysis to compute mass-dependent black hole mass
!% distribution functions.

module Galacticus_Output_Analyses_Mass_Dpndnt_BH_Dstrbtins
  !% Performs analysis to compute mass-dependent black hole mass distributions. Currently
  !% supported distributions include:
  !% \begin{itemize}
  !% \item The compilation from \cite{kormendy_coevolution_2013}.
  !% \end{itemize}
  use, intrinsic :: ISO_C_Binding
  use Galacticus_Nodes
  use FGSL
  use Tables
  use Galactic_Structure_Options
  use Geometry_Surveys
  use Galactic_Filters
  use Numerical_Constants_Astronomical
  use Numerical_Constants_Prefixes
  implicit none
  private
  public :: Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins, Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins_Output
  
  ! Record of module initialization.
  logical                                                            :: moduleInitialized                   =.false.

  ! Record of whether this analysis is active.
  logical                                                            :: analysisActive

  ! Number of supported mass functions.
  integer          , parameter                                       :: blackHoleDistributionsSupportedCount=1

  ! Labels for supported metallicity distributions.
  character(len=31), dimension(blackHoleDistributionsSupportedCount) :: blackHoleMassDistributionLabels     =        &
       & [                                                                                                           &
       &  'blackHoleDistributionKormendyHo'                                                                          &
       & ]

  ! Interface for mass mapping functions.
  abstract interface
     double precision function Map_Mass(mass,node)
       import treeNode
       double precision          , intent(in   )          :: mass
       type            (treeNode), intent(inout), pointer :: node
     end function Map_Mass
  end interface

  ! Type for descriptors of black hole mass distributions.
  type :: blackHoleDistributionDescriptor
     double precision                                           :: massSystematicLogM0
     double precision                                           :: massLogarithmicMinimum
     integer                                                    :: massSystematicCoefficientCount, blackHoleMassSystematicCoefficientCount, &
          &                                                        massRandomCoefficientCount    , blackHoleMassRandomCoefficientCount
     double precision                                           :: massUnitsInSI
     character       (len= 32            )                      :: label
     character       (len=128            )                      :: comment
     procedure       (Map_Mass           ), pointer    , nopass :: mapGalaxyMass                 , mapBlackHoleMass
     class           (surveyGeometryClass), allocatable         :: geometry
     class           (galacticFilterClass), allocatable         :: filter
  end type blackHoleDistributionDescriptor

  ! Metallicity distribution descriptors.
  type(blackHoleDistributionDescriptor), dimension(blackHoleDistributionsSupportedCount), target :: blackHoleDistributionDescriptors= &
       & [                                                                                                                            &
       ! SDSS galaxy stellar-phase metallicities from Kormendy & Ho (2013). 
       &                           blackHoleDistributionDescriptor(                                                                   &
       &                                                             11.0000d+0                                                     , &
       &                                                              8.000d0                                                       , &
       &                                                              2                                                             , &
       &                                                              2                                                             , &
       &                                                              2                                                             , &
       &                                                              2                                                             , &
       &                                                              massSolar                                                     , &
       &                                                              'blackHoleDistributionKormendyHo'                             , &
       &                                                              'Kormendy & Ho (2013) black hole mass compilation'            , &
       &                                                              null()                                                        , &
       &                                                              null()                                                        , &
       &                                                              null()                                                        , &
       &                                                              null()                                                          &
       &                                                            )                                                                 &
       & ]

  ! Type to store black hole mass distributions.
  type :: blackHoleMassDistribution
     ! Copy of the mass function descriptor for this mass function.
     type            (blackHoleDistributionDescriptor), pointer                       :: descriptor
     ! Parameters for the systematic error model.
     double precision                                 , allocatable, dimension(:    ) :: massSystematicCoefficients         , blackHoleMassSystematicCoefficients, &
          &                                                                              massRandomCoefficients             , blackHoleMassRandomCoefficients
     ! The number of bins.
     integer                                                                          :: massesCount                        , blackHoleMassesCount
     ! Arrays for the masses, radii and size function.
     double precision                                 , allocatable, dimension(:    ) :: masses                             , massesLogarithmic                  , &
          &                                                                              massesLogarithmicMinimum           , massesLogarithmicMaximum           , &
          &                                                                              blackHoleMasses                    , blackHoleMassesLogarithmic         , &
          &                                                                              blackHoleMassesLogarithmicMinimum  , blackHoleMassesLogarithmicMaximum  , &
          &                                                                              blackHoleMassDistributionWeights
     double precision                                 , allocatable, dimension(:,:  ) :: blackHoleMassDistribution          , outputWeight
     ! Arrays for accumulation of of main branch galaxies
     double precision                                 , allocatable, dimension(:,:,:) :: mainBranchGalaxyWeights            , mainBranchGalaxyWeightsSquared
     ! Array for the covariance matrix.
     double precision                                 , allocatable, dimension(:,:  ) :: blackHoleMassDistributionCovariance
     ! Cosmology conversion factors.
     double precision                                 , allocatable, dimension(:    ) :: cosmologyConversionMass            , cosmologyConversionBlackHoleMass
  end type blackHoleMassDistribution

  ! Mass functions.
  type(blackHoleMassDistribution), allocatable, dimension(:) :: blackHoleDistributions

  ! Type for storing temporary black hole mass functions during cumulation.
  type :: blackHoleMassDistributionWork
     double precision, allocatable, dimension(:,:) :: blackHoleMassDistribution
     double precision, allocatable, dimension(  :) :: blackHoleMassDistributionWeights
  end type blackHoleMassDistributionWork

  ! Work array.
  type(blackHoleMassDistributionWork), allocatable, dimension(:) :: galaxyWork
  !$omp threadprivate(galaxyWork)
  
  ! Options controlling binning in halo mass.
  integer                     :: analysisBlackHoleMassDistributionCovarianceModel
  integer         , parameter :: analysisBlackHoleMassDistributionCovarianceModelPoisson =1
  integer         , parameter :: analysisBlackHoleMassDistributionCovarianceModelBinomial=2
  integer                     :: analysisBlackHoleMassDistributionsHaloMassBinsCount       , analysisBlackHoleMassDistributionsHaloMassBinsPerDecade
  double precision            :: analysisBlackHoleMassDistributionsHaloMassMinimum         , analysisBlackHoleMassDistributionsHaloMassMaximum           , &
       &                         analysisBHMassDstrbtnsHaloMassIntervalLogarithmicInverse  , analysisBlackHoleMassDistributionsHaloMassMinimumLogarithmic

contains

  !# <mergerTreeAnalysisTask>
  !#  <unitName>Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins</unitName>
  !# </mergerTreeAnalysisTask>
  subroutine Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins(tree,node,nodeStatus,iOutput,mergerTreeAnalyses)
    !% Construct black hole mass distributions to compare to various observational determinations.
    use, intrinsic :: ISO_C_Binding
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
    use Pseudo_Random
    use Numerical_Comparison
    use String_Handling
    use IO_HDF5
    use Vectors
    use Galacticus_Output_Analyses_Cosmology_Scalings
    use Galacticus_Output_Merger_Tree_Data
    use Numerical_Ranges
    use Numerical_Constants_Astronomical
    implicit none
    type            (mergerTree                    ), intent(in   )                 :: tree
    type            (treeNode                      ), intent(inout), pointer        :: node
    integer                                         , intent(in   )                 :: nodeStatus
    integer         (c_size_t                      ), intent(in   )                 :: iOutput
    type            (varying_string                ), intent(in   ), dimension(:  ) :: mergerTreeAnalyses
    double precision                                , allocatable  , dimension(:  ) :: blackHoleMasses ,masses
    class           (nodeComponentBasic            )               , pointer        :: basic
    class           (nodeComponentBlackHole        )               , pointer        :: blackHole
    class           (cosmologyFunctionsClass       )               , pointer        :: cosmologyFunctionsModel
    type            (cosmologyFunctionsMatterLambda)                                :: cosmologyFunctionsObserved
    type            (cosmologyParametersSimple     )               , pointer        :: cosmologyParametersObserved
    integer                                         , parameter                     :: massCountPerDecade             =10
    integer                                         , parameter                     :: blackHoleMassCountPerDecade    =10
    double precision                                , parameter                     :: massRandomErrorMinimum         =1.0d-3
    double precision                                , parameter                     :: blackHoleMassRandomErrorMinimum=1.0d-3
    integer         (c_size_t                      )                                :: k,jOutput
    integer                                                                         :: i,j,l,currentAnalysis,activeAnalysisCount,haloMassBin
    double precision                                                                :: dataHubbleParameter ,mass,massLogarithmic&
         &,massRandomError,blackHoleMassLogarithmic,blackHoleMass,blackHoleMassRandomError,dataOmegaDarkEnergy,dataOmegaMatter,redshift,timeMinimum,timeMaximum,distanceMinimum,distanceMaximum,unitsInSI,blackHoleMassMinimum,blackHoleMassMaximum,massMinimum,massMaximum
    type            (varying_string                )                                :: parameterName&
         &,analysisBlackHoleMassDistributionCovarianceModelText,cosmologyScalingMass,cosmologyScalingBlackHoleMass,scaling,message
    type            (hdf5Object                    )                                :: dataFile,massDataset,parametersGroup,blackHoleMassDataset
    
    ! Initialize the module if necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins_Initialize)
       if (.not.moduleInitialized) then
          !@ <inputParameter>
          !@   <name>analysisBlackHoleMassDistributionCovarianceModel</name>
          !@   <defaultValue>binomial</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The model to use when constructing the black hole mass distribution covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisBlackHoleMassDistributionCovarianceModel',analysisBlackHoleMassDistributionCovarianceModelText,defaultValue='Poisson')
          select case (char(analysisBlackHoleMassDistributionCovarianceModelText))
          case ( 'Poisson'  )
             analysisBlackHoleMassDistributionCovarianceModel=analysisBlackHoleMassDistributionCovarianceModelPoisson
          case ( 'binomial' )
             analysisBlackHoleMassDistributionCovarianceModel=analysisBlackHoleMassDistributionCovarianceModelBinomial
          case default
             call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins','unrecognized value for "analysisBlackHoleMassDistributionCovarianceModel" - allowed values are "Poisson", and "binomial"')
          end select
          !@ <inputParameter>
          !@   <name>analysisBlackHoleMassDistributionsHaloMassBinsPerDecade</name>
          !@   <defaultValue>10</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The number of bins per decade of halo mass to use when constructing the mass function covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisBlackHoleMassDistributionsHaloMassBinsPerDecade',analysisBlackHoleMassDistributionsHaloMassBinsPerDecade,defaultValue=10)
          !@ <inputParameter>
          !@   <name>analysisBlackHoleMassDistributionsHaloMassMinimum</name>
          !@   <defaultValue>$10^8M_\odot$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The minimum halo mass to consider when constructing the mass function covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisBlackHoleMassDistributionsHaloMassMinimum',analysisBlackHoleMassDistributionsHaloMassMinimum,defaultValue=1.0d8)
          !@ <inputParameter>
          !@   <name>analysisBlackHoleMassDistributionsHaloMassMaximum</name>
          !@   <defaultValue>$10^{16}M_\odot$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The maximum halo mass to consider when constructing the mass function covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisBlackHoleMassDistributionsHaloMassMaximum',analysisBlackHoleMassDistributionsHaloMassMaximum,defaultValue=1.0d16)
          analysisBlackHoleMassDistributionsHaloMassMinimumLogarithmic=    log10( analysisBlackHoleMassDistributionsHaloMassMinimum)
          analysisBlackHoleMassDistributionsHaloMassBinsCount         =int(                                                             &
               &                                                           log10(                                                       &
               &                                                                  analysisBlackHoleMassDistributionsHaloMassMaximum     &
               &                                                                 /analysisBlackHoleMassDistributionsHaloMassMinimum     &
               &                                                                )                                                       &
               &                                                         *dble(analysisBlackHoleMassDistributionsHaloMassBinsPerDecade) &
               &                                                         +0.5d0&
               &                                                        )
          analysisBHMassDstrbtnsHaloMassIntervalLogarithmicInverse    = dble (analysisBlackHoleMassDistributionsHaloMassBinsCount)      &
               &                                                       /log10(                                                          &
               &                                                               analysisBlackHoleMassDistributionsHaloMassMaximum        &
               &                                                              /analysisBlackHoleMassDistributionsHaloMassMinimum        &
               &                                                             )
          ! Establish mapping functions for black hole mass distribution descriptors.
          blackHoleDistributionDescriptors(1)%mapGalaxyMass    => Map_Galaxy_Mass_Kormendy_Ho
          blackHoleDistributionDescriptors(1)%mapBlackHoleMass => Map_Black_Hole_Mass_Kormendy_Ho
          ! Determine how many supported blackHoleMass distributions are requested.
          activeAnalysisCount=0
          do i=1,blackHoleDistributionsSupportedCount
             if (any(trim(mergerTreeAnalyses) == trim(blackHoleMassDistributionLabels(i)))) activeAnalysisCount=activeAnalysisCount+1
          end do
          ! Allocate mass function arrays and populate with required data.
          if (activeAnalysisCount <= 0) then
             analysisActive=.false.
          else
             analysisActive=.true.
             ! Establish survey geometries.
             allocate(surveyGeometryFullSky :: blackHoleDistributionDescriptors(1)%geometry)
             select type (g => blackHoleDistributionDescriptors(1)%geometry)
             type is (surveyGeometryFullSky)
                g=surveyGeometryFullSky(0.0d0,0.06d0)
             end select             
             ! Establish sample filters.
             allocate(galacticFilterAlways :: blackHoleDistributionDescriptors(1)%filter)
             select type (f => blackHoleDistributionDescriptors(1)%filter)
             type is (galacticFilterAlways)
                f=galacticFilterAlways()
             end select             
             ! Initialize analyses.
             currentAnalysis=0
             allocate(blackHoleDistributions(activeAnalysisCount))
             cosmologyFunctionsModel => cosmologyFunctions()
             do i=1,size(mergerTreeAnalyses)
                do j=1,blackHoleDistributionsSupportedCount
                   if (mergerTreeAnalyses(i) == trim(blackHoleMassDistributionLabels(j))) then
                      currentAnalysis=currentAnalysis+1
                      ! Set a pointer to the descriptor for this black hole mass distribution function.
                      blackHoleDistributions(currentAnalysis)%descriptor => blackHoleDistributionDescriptors(j)
                      ! Read parameters of the systematic error model.
                      if (blackHoleDistributionDescriptors(j)%massSystematicCoefficientCount > 0) then
                         allocate(blackHoleDistributions(currentAnalysis)%massSystematicCoefficients(blackHoleDistributionDescriptors(j)%massSystematicCoefficientCount))
                         do k=1,blackHoleDistributionDescriptors(j)%massSystematicCoefficientCount
                            parameterName=trim(blackHoleMassDistributionLabels(j))//'MassSystematic'
                            parameterName=parameterName//(k-1)
                            !@ <inputParameter>
                            !@   <regEx>(blackHoleDistributionKormendyHo)MassSystematic[0-9]+</regEx>
                            !@   <defaultValue>0</defaultValue>
                            !@   <attachedTo>module</attachedTo>
                            !@   <description>
                            !@     Mass-dependent black hole mass distribution mass systematic parameters.
                            !@   </description>
                            !@   <type>real</type>
                            !@   <cardinality>1</cardinality>
                            !@ </inputParameter>
                            call Get_Input_Parameter(char(parameterName),blackHoleDistributions(currentAnalysis)%massSystematicCoefficients(k),defaultValue=0.0d0)
                         end do
                      end if
                      if (blackHoleDistributionDescriptors(j)%blackHoleMassSystematicCoefficientCount > 0) then
                         allocate(blackHoleDistributions(currentAnalysis)%blackHoleMassSystematicCoefficients(blackHoleDistributionDescriptors(j)%blackHoleMassSystematicCoefficientCount))
                         do k=1,blackHoleDistributionDescriptors(j)%blackHoleMassSystematicCoefficientCount
                            parameterName=trim(blackHoleMassDistributionLabels(j))//'BlackHoleMassSystematic'
                            parameterName=parameterName//(k-1)
                            !@ <inputParameter>
                            !@   <regEx>(blackHoleDistributionKormendyHo)BlackHoleMassSystematic[0-9]+</regEx>
                            !@   <defaultValue>0</defaultValue>
                            !@   <attachedTo>module</attachedTo>
                            !@   <description>
                            !@     Mass-dependent black hole mass distribution blackHoleMass systematic parameters.
                            !@   </description>
                            !@   <type>real</type>
                            !@   <cardinality>1</cardinality>
                            !@ </inputParameter>
                            call Get_Input_Parameter(char(parameterName),blackHoleDistributions(currentAnalysis)%blackHoleMassSystematicCoefficients(k),defaultValue=0.0d0)
                         end do
                      end if
                      ! Read parameters of the random error model.
                      if (blackHoleDistributionDescriptors(j)%massRandomCoefficientCount > 0) then
                         allocate(blackHoleDistributions(currentAnalysis)%massRandomCoefficients(blackHoleDistributionDescriptors(j)%massRandomCoefficientCount))
                         do k=1,blackHoleDistributionDescriptors(j)%massRandomCoefficientCount
                            parameterName=trim(blackHoleMassDistributionLabels(j))//'MassRandom'
                            parameterName=parameterName//(k-1)
                            !@ <inputParameter>
                            !@   <regEx>(blackHoleDistributionKormendyHo)MassRandom[0-9]+</regEx>
                            !@   <defaultValue>0</defaultValue>
                            !@   <attachedTo>module</attachedTo>
                            !@   <description>
                            !@     Mass-dependent black hole mass distribution mass random parameters.
                            !@   </description>
                            !@   <type>real</type>
                            !@   <cardinality>1</cardinality>
                            !@ </inputParameter>
                            call Get_Input_Parameter(char(parameterName),blackHoleDistributions(currentAnalysis)%massRandomCoefficients(k),defaultValue=0.0d0)
                         end do
                      end if
                      if (blackHoleDistributionDescriptors(j)%blackHoleMassRandomCoefficientCount > 0) then
                         allocate(blackHoleDistributions(currentAnalysis)%blackHoleMassRandomCoefficients(blackHoleDistributionDescriptors(j)%blackHoleMassRandomCoefficientCount))
                         do k=1,blackHoleDistributionDescriptors(j)%blackHoleMassRandomCoefficientCount
                            parameterName=trim(blackHoleMassDistributionLabels(j))//'BlackHoleMassRandom'
                            parameterName=parameterName//(k-1)
                            !@ <inputParameter>
                            !@   <regEx>(sdssStellarBlackHoleMassDistribution|sdssGasBlackHoleMassDistribution)Z[0-9\.]+BlackHoleMassRandom[0-9]+</regEx>
                            !@   <defaultValue>0</defaultValue>
                            !@   <attachedTo>module</attachedTo>
                            !@   <description>
                            !@     Mass-dependent blackHoleMass distribution blackHoleMass random parameters.
                            !@   </description>
                            !@   <type>real</type>
                            !@   <cardinality>1</cardinality>
                            !@ </inputParameter>
                            call Get_Input_Parameter(char(parameterName),blackHoleDistributions(currentAnalysis)%blackHoleMassRandomCoefficients(k),defaultValue=0.0d0)
                         end do
                      end if
                      ! Read the appropriate observational data definition.
                      select case (trim(blackHoleMassDistributionLabels(j)))
                      case ('blackHoleDistributionKormendyHo')
                         ! Kormendy & Ho black hole mass compilation.
                         !$omp critical(HDF5_Access)
                         call dataFile%openFile(char(Galacticus_Input_Path())//"data/observations/blackHoles/blackHoleMassVsBulgeMass_KormendyHo2013.hdf5",readOnly=.true.)
                         ! Read masses.
                         call dataFile%readDataset("massBulge",masses)
                         massDataset=dataFile%openDataset("massBulge")
                         call massDataset%readAttribute("cosmologyScaling",cosmologyScalingMass)
                         call massDataset%readAttribute("scaling"         ,scaling             )
                         call massDataset%readAttribute("unitsInSI"       ,unitsInSI           )
                         call massDataset%close        (                                       )
                         select case (char(scaling))
                         case('linear')
                            ! Nothing to do in this case.
                         case('log10' )
                            ! Convert from log10 mass to linear mass.
                            masses=10.0d0**masses
                         case default
                            call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins','unrecognized scaling')
                         end select
                         masses=masses*unitsInSI/massSolar
                         ! Find range of galaxy masses.
                         massMinimum=1.0d-1*minval(masses)
                         massMaximum=1.0d+1*maxval(masses)
                         deallocate(masses)
                         call dataFile%readDataset("massBlackHole",blackHoleMasses)
                         blackHoleMassDataset=dataFile%openDataset("massBlackHole")
                         call blackHoleMassDataset%readAttribute("cosmologyScaling",cosmologyScalingBlackHoleMass)
                         call blackHoleMassDataset%readAttribute("scaling"         ,scaling                      )
                         call blackHoleMassDataset%readAttribute("unitsInSI"       ,unitsInSI                    )
                         call blackHoleMassDataset%close        (                                                )
                         select case (char(scaling))
                         case('linear')
                            ! Nothing to do in this case.
                         case('log10' )
                            ! Convert from log10 mass to linear mass.
                            blackHoleMasses=10.0d0**blackHoleMasses
                         case default
                            call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins','unrecognized scaling')
                         end select
                         blackHoleMasses=blackHoleMasses*unitsInSI/massSolar
                         ! Find range of black hole masses.
                         blackHoleMassMinimum=1.0d-1*minval(blackHoleMasses)
                         blackHoleMassMaximum=1.0d+1*maxval(blackHoleMasses)
                         ! Construct grids of galaxy and black hole masses.
                         blackHoleDistributions(currentAnalysis)%massesCount         =int(log10(         massMaximum/         massMinimum)*         massCountPerDecade)+1
                         blackHoleDistributions(currentAnalysis)%blackHoleMassesCount=int(log10(blackHoleMassMaximum/blackHoleMassMinimum)*blackHoleMassCountPerDecade)+1
                         allocate(blackHoleDistributions(currentAnalysis)%         masses(blackHoleDistributions(currentAnalysis)%         massesCount))
                         allocate(blackHoleDistributions(currentAnalysis)%blackHoleMasses(blackHoleDistributions(currentAnalysis)%blackHoleMassesCount))
                         blackHoleDistributions(currentAnalysis)%         masses=Make_Range(                                                              &
                              &                                                                                                              massMinimum, &
                              &                                                                                                              massMaximum, &
                              &                                                             blackHoleDistributions(currentAnalysis)%         massesCount, &
                              &                                                             rangeTypeLogarithmic                                          &
                              &                                                            )

                         blackHoleDistributions(currentAnalysis)%blackHoleMasses=Make_Range(                                                              &
                              &                                                                                                     blackHoleMassMinimum, &
                              &                                                                                                     blackHoleMassMaximum, &
                              &                                                             blackHoleDistributions(currentAnalysis)%blackHoleMassesCount, &
                              &                                                             rangeTypeLogarithmic                                          &
                              &                                                            )                         
                         call allocateArray(blackHoleDistributions(currentAnalysis)%         massesLogarithmic         ,[                                                             blackHoleDistributions(currentAnalysis)%massesCount                                                    ])
                         call allocateArray(blackHoleDistributions(currentAnalysis)%         massesLogarithmicMinimum  ,[                                                             blackHoleDistributions(currentAnalysis)%massesCount                                                    ])
                         call allocateArray(blackHoleDistributions(currentAnalysis)%         massesLogarithmicMaximum  ,[                                                             blackHoleDistributions(currentAnalysis)%massesCount                                                    ])
                         call allocateArray(blackHoleDistributions(currentAnalysis)%blackHoleMassesLogarithmic         ,[blackHoleDistributions(currentAnalysis)%blackHoleMassesCount                                                                                                        ])
                         call allocateArray(blackHoleDistributions(currentAnalysis)%blackHoleMassesLogarithmicMinimum  ,[blackHoleDistributions(currentAnalysis)%blackHoleMassesCount                                                                                                        ])
                         call allocateArray(blackHoleDistributions(currentAnalysis)%blackHoleMassesLogarithmicMaximum  ,[blackHoleDistributions(currentAnalysis)%blackHoleMassesCount                                                                                                        ])
                         call allocateArray(blackHoleDistributions(currentAnalysis)%blackHoleMassDistribution          ,[blackHoleDistributions(currentAnalysis)%blackHoleMassesCount,blackHoleDistributions(currentAnalysis)%massesCount                                                    ])
                         call allocateArray(blackHoleDistributions(currentAnalysis)%blackHoleMassDistributionWeights   ,[                                                             blackHoleDistributions(currentAnalysis)%massesCount                                                    ])
                         call allocateArray(blackHoleDistributions(currentAnalysis)%mainBranchGalaxyWeights            ,[blackHoleDistributions(currentAnalysis)%blackHoleMassesCount,blackHoleDistributions(currentAnalysis)%massesCount,analysisBlackHoleMassDistributionsHaloMassBinsCount])
                         call allocateArray(blackHoleDistributions(currentAnalysis)%mainBranchGalaxyWeightsSquared     ,[blackHoleDistributions(currentAnalysis)%blackHoleMassesCount,blackHoleDistributions(currentAnalysis)%massesCount,analysisBlackHoleMassDistributionsHaloMassBinsCount])
                         call allocateArray(blackHoleDistributions(currentAnalysis)%blackHoleMassDistributionCovariance,[                                                                                                                                                                       &
                              &                                                                                        blackHoleDistributions(currentAnalysis)%blackHoleMassesCount*blackHoleDistributions(currentAnalysis)%massesCount,                                                      &
                              &                                                                                        blackHoleDistributions(currentAnalysis)%blackHoleMassesCount*blackHoleDistributions(currentAnalysis)%massesCount                                                       &
                              &                                                                                       ]                                                                                                                                                                       &
                              &          )
                         blackHoleDistributions(currentAnalysis)%         massesLogarithmic=log10(blackHoleDistributions(currentAnalysis)%         masses)
                         blackHoleDistributions(currentAnalysis)%blackHoleMassesLogarithmic=log10(blackHoleDistributions(currentAnalysis)%blackHoleMasses)
                         do k=1,blackHoleDistributions(currentAnalysis)%massesCount
                            if (k ==                                                            1) then
                               blackHoleDistributions                (currentAnalysis)%         massesLogarithmicMinimum(k  )= &
                                    & +        blackHoleDistributions(currentAnalysis)%         massesLogarithmic       (k  )  &
                                    & -0.5d0*(                                                                                 &
                                    &         +blackHoleDistributions(currentAnalysis)%         massesLogarithmic       (k+1)  &
                                    &         -blackHoleDistributions(currentAnalysis)%         massesLogarithmic       (k  )  &
                                    &        )
                            else
                               blackHoleDistributions                (currentAnalysis)%         massesLogarithmicMinimum(k  )= &
                                    & +0.5d0*(                                                                                 &
                                    &         +blackHoleDistributions(currentAnalysis)%         massesLogarithmic       (k-1)  &
                                    &         +blackHoleDistributions(currentAnalysis)%         massesLogarithmic       (k  )  &
                                    &        )
                            end if
                            if (k == blackHoleDistributions(currentAnalysis)%         massesCount) then
                               blackHoleDistributions                (currentAnalysis)%         massesLogarithmicMaximum(k  )= &
                                    & +        blackHoleDistributions(currentAnalysis)%         massesLogarithmic       (k  )  &
                                    & +0.5d0*(                                                                                 &
                                    &         +blackHoleDistributions(currentAnalysis)%         massesLogarithmic       (k  )  &
                                    &         -blackHoleDistributions(currentAnalysis)%         massesLogarithmic       (k-1)  &
                                    &        )
                            else
                               blackHoleDistributions                (currentAnalysis)%         massesLogarithmicMaximum(k  )= &
                                    & +0.5d0*(                                                                                 &
                                    &         +blackHoleDistributions(currentAnalysis)%         massesLogarithmic       (k+1)  &
                                    &         +blackHoleDistributions(currentAnalysis)%         massesLogarithmic       (k  )  &
                                    &        )
                            end if
                         end do
                         do k=1,blackHoleDistributions(currentAnalysis)%blackHoleMassesCount
                            if (k ==                                                            1) then
                               blackHoleDistributions                (currentAnalysis)%blackHoleMassesLogarithmicMinimum(k  )= &
                                    & +        blackHoleDistributions(currentAnalysis)%blackHoleMassesLogarithmic       (k  )  &
                                    & -0.5d0*(                                                                                 &
                                    &         +blackHoleDistributions(currentAnalysis)%blackHoleMassesLogarithmic       (k+1)  &
                                    &         -blackHoleDistributions(currentAnalysis)%blackHoleMassesLogarithmic       (k  )  &
                                    &        )
                            else
                               blackHoleDistributions                (currentAnalysis)%blackHoleMassesLogarithmicMinimum(k  )= &
                                    & +0.5d0*(                                                                                 &
                                    &         +blackHoleDistributions(currentAnalysis)%blackHoleMassesLogarithmic       (k-1)  &
                                    &         +blackHoleDistributions(currentAnalysis)%blackHoleMassesLogarithmic       (k  )  &
                                    &        )
                            end if
                            if (k == blackHoleDistributions(currentAnalysis)%blackHoleMassesCount) then
                               blackHoleDistributions                (currentAnalysis)%blackHoleMassesLogarithmicMaximum(k  )= &
                                    & +        blackHoleDistributions(currentAnalysis)%blackHoleMassesLogarithmic       (k  )  &
                                    & +0.5d0*(                                                                                 &
                                    &         +blackHoleDistributions(currentAnalysis)%blackHoleMassesLogarithmic       (k  )  &
                                    &         -blackHoleDistributions(currentAnalysis)%blackHoleMassesLogarithmic       (k-1)  &
                                    &        )
                            else
                               blackHoleDistributions                (currentAnalysis)%blackHoleMassesLogarithmicMaximum(k  )= &
                                    & +0.5d0*(                                                                                 &
                                    &         +blackHoleDistributions(currentAnalysis)%blackHoleMassesLogarithmic       (k+1)  &
                                    &         +blackHoleDistributions(currentAnalysis)%blackHoleMassesLogarithmic       (k  )  &
                                    &        )
                            end if
                         end do
                         ! Extract cosmological parameters.
                         parametersGroup=dataFile%openGroup("Parameters")
                         call parametersGroup%readAttribute("HubbleConstant" ,dataHubbleParameter)
                         call parametersGroup%readAttribute("OmegaMatter"    ,dataOmegaMatter    )
                         call parametersGroup%readAttribute("OmegaDarkEnergy",dataOmegaDarkEnergy)
                         call parametersGroup%close        (                                     )
                         ! Finished reading data.
                         call dataFile%close()
                         !$omp end critical(HDF5_Access)
                         ! Create the observed cosmology.
                         allocate(cosmologyParametersObserved)
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
                         ! Initialize all data.
                         blackHoleDistributions(currentAnalysis)%blackHoleMassDistribution          =0.0d0
                         blackHoleDistributions(currentAnalysis)%blackHoleMassDistributionWeights   =0.0d0
                         blackHoleDistributions(currentAnalysis)%blackHoleMassDistributionCovariance=0.0d0
                         blackHoleDistributions(currentAnalysis)%mainBranchGalaxyWeights            =0.0d0
                         blackHoleDistributions(currentAnalysis)%mainBranchGalaxyWeightsSquared     =0.0d0
                      case default
                         call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins','unknown black hole mass function')
                      end select
                      ! Get cosmological conversion factors.
                      call allocateArray(blackHoleDistributions(currentAnalysis)%cosmologyConversionMass         ,[Galacticus_Output_Time_Count()])
                      call allocateArray(blackHoleDistributions(currentAnalysis)%cosmologyConversionBlackHoleMass,[Galacticus_Output_Time_Count()])
                      do jOutput=1,Galacticus_Output_Time_Count()
                         redshift=                                                                                      &
                              &   cosmologyFunctionsModel %redshiftFromExpansionFactor(                                 &
                              &    cosmologyFunctionsModel%expansionFactor             (                                &
                              &                                                         Galacticus_Output_Time(jOutput) &
                              &                                                        )                                &
                              &                                                       )
                         call Cosmology_Conversion_Factors(                                                                                                              &
                              &                            redshift                                                                                                    , &
                              &                            cosmologyFunctionsModel                                                                                     , &
                              &                            cosmologyFunctionsObserved                                                                                  , &
                              &                            cosmologyScalingMass      =cosmologyScalingMass                                                             , &
                              &                            cosmologyConversionMass   =blackHoleDistributions(currentAnalysis)%cosmologyConversionMass         (jOutput)  &
                              &                           )
                         call Cosmology_Conversion_Factors(                                                                                                              &
                              &                            redshift                                                                                                    , &
                              &                            cosmologyFunctionsModel                                                                                     , &
                              &                            cosmologyFunctionsObserved                                                                                  , &
                              &                            cosmologyScalingMass      =cosmologyScalingBlackHoleMass                                                    , &
                              &                            cosmologyConversionMass   =blackHoleDistributions(currentAnalysis)%cosmologyConversionBlackHoleMass(jOutput)  &
                              &                           )
                      end do
                      nullify(cosmologyParametersObserved)
                      ! Compute output weights for black hole mass distribution.
                      call allocateArray(blackHoleDistributions(currentAnalysis)%outputWeight,[int(blackHoleDistributions(currentAnalysis)%massesCount,kind=c_size_t),Galacticus_Output_Time_Count()])
                      blackHoleDistributions(currentAnalysis)%outputWeight=0.0d0
                      do k=1,blackHoleDistributions(currentAnalysis)%massesCount
                         do jOutput=1,Galacticus_Output_Time_Count()
                            do l=1,blackHoleDistributions(currentAnalysis)%descriptor%geometry%fieldCount()
                               if (jOutput == Galacticus_Output_Time_Count()) then
                                  timeMaximum=     Galacticus_Output_Time(jOutput)
                               else
                                  timeMaximum=sqrt(Galacticus_Output_Time(jOutput)*Galacticus_Output_Time(jOutput+1))
                               end if
                               if (jOutput ==                              1) then
                                  timeMinimum=     Galacticus_Output_Time(jOutput)
                               else
                                  timeMinimum=sqrt(Galacticus_Output_Time(jOutput)*Galacticus_Output_Time(jOutput-1))
                               end if
                               distanceMinimum=max(                                                                                                                                  &
                                    &              cosmologyFunctionsModel%distanceComoving(timeMaximum)                                                                           , &
                                    &              blackHoleDistributions(currentAnalysis)%descriptor%geometry%distanceMinimum(blackHoleDistributions(currentAnalysis)%masses(k),l)  &
                                    &             )
                               distanceMaximum=min(                                                                                                                                  &
                                    &              cosmologyFunctionsModel%distanceComoving(timeMinimum)                                                                           , &
                                    &              blackHoleDistributions(currentAnalysis)%descriptor%geometry%distanceMaximum(blackHoleDistributions(currentAnalysis)%masses(k),l)  &
                                    &             )
                               blackHoleDistributions        (currentAnalysis)%outputWeight                    (k,jOutput)  &
                                    & =blackHoleDistributions(currentAnalysis)%outputWeight                    (k,jOutput)  &
                                    & +blackHoleDistributions(currentAnalysis)%descriptor  %geometry%solidAngle(  l      )  &
                                    & /3.0d0                                                                                &
                                    & *                                                                                     &
                                    & max(                                                                                  &
                                    &     +0.0d0                                                                          , &
                                    &     +distanceMaximum**3                                                               &
                                    &     -distanceMinimum**3                                                               &
                                    &    )
                            end do
                         end do
                         where(blackHoleDistributions(currentAnalysis)%outputWeight(k,:) < 0.0d0)
                            blackHoleDistributions                  (currentAnalysis)%outputWeight(k,:)&
                                 &       =0.0d0
                         end where
                         if (any(blackHoleDistributions(currentAnalysis)%outputWeight(k,:) > 0.0d0)) then
                            blackHoleDistributions                  (currentAnalysis)%outputWeight(k,:)  &
                                 &       =    blackHoleDistributions(currentAnalysis)%outputWeight(k,:)  &
                                 &       /sum(blackHoleDistributions(currentAnalysis)%outputWeight(k,:))
                         else
                            message="blackHoleMass distribution '"//trim(blackHoleDistributions(currentAnalysis)%descriptor%label)//"' bin "
                            message=message//k//" has zero weights"
                            call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins',message)
                         end if
                      end do
                      ! Ensure that spheroid and black hole components support relevant mass properties.
                      if (.not.defaultSpheroidComponent %massStellarIsGettable())                                                     &
                           & call Galacticus_Error_Report                                                                             &
                           & (                                                                                                        &
                           &  'Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins'                                                 , &
                           &  'This analysis requires that the "massStellar" property of the spheroid is gettable.'//                 &
                           &  Galacticus_Component_List(                                                                              &
                           &                            'spheroid'                                                                  , &
                           &                             defaultSpheroidComponent %massStellarAttributeMatch(requireGettable=.true.)  &
                           &                           )                                                                              &
                           & )
                      if (.not.defaultBlackHoleComponent%       massIsGettable())                                                     &
                           & call Galacticus_Error_Report                                                                             &
                           & (                                                                                                        &
                           &  'Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins'                                                 , &
                           &  'This analysis requires that the "mass" property of the black hole is gettable.'//                      &
                           &  Galacticus_Component_List(                                                                              &
                           &                            'blackHole'                                                                 , &
                           &                             defaultBlackHoleComponent%       massAttributeMatch(requireGettable=.true.)  &
                           &                           )                                                                              &
                           & )
                      exit
                   end if
                end do
             end do
          end if
          ! Record that module is initialized.
          moduleInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins_Initialize)
    end if
    ! Return if this analysis is not active.
    if (.not.analysisActive) return
    ! Return if this is a tree finalization.
    if (nodeStatus == nodeStatusFinal) return
    ! Allocate work arrays.
    if (.not.allocated(galaxyWork)) allocate(galaxyWork(size(blackHoleDistributions)))
    ! Iterate over active analyses.
    do i=1,size(blackHoleDistributions)
       ! Cycle if this black hole mass distribution receives no contribution from this output.
       if (all(blackHoleDistributions(i)%outputWeight(:,iOutput) <= 0.0d0)) cycle
       ! Allocate workspace.
       if (.not.allocated(galaxyWork(i)%blackHoleMassDistribution       )) call allocateArray(galaxyWork(i)%blackHoleMassDistribution       ,[blackHoleDistributions(i)%blackHoleMassesCount,blackHoleDistributions(i)%massesCount])
       if (.not.allocated(galaxyWork(i)%blackHoleMassDistributionWeights)) call allocateArray(galaxyWork(i)%blackHoleMassDistributionWeights,[                                               blackHoleDistributions(i)%massesCount])
       ! Filter the galaxy.
       if (.not.blackHoleDistributions(i)%descriptor%filter%passes(node)) cycle
       ! Get the spheroid hole mass.
       mass=Galactic_Structure_Enclosed_Mass(node,radiusLarge,componentType=componentTypeSpheroid,massType=massTypeStellar)
       if (mass <= 0.0d0) cycle
       ! Get the black hole mass.
       blackHole => node     %blackHole()
       blackHoleMass =  blackHole%mass     ()
       if (blackHoleMass <= 0.0d0) cycle
       ! Map galaxy and black hole masses.
       if (associated(blackHoleDistributions(i)%descriptor%mapGalaxyMass   )) mass         =blackHoleDistributions(i)%descriptor%mapGalaxyMass   (         mass,node)
       if (associated(blackHoleDistributions(i)%descriptor%mapBlackHoleMass)) blackHoleMass=blackHoleDistributions(i)%descriptor%mapBlackHoleMass(blackHoleMass,node)
       ! Convert masses for cosmology and systematics.
       mass=mass*blackHoleDistributions(i)%cosmologyConversionMass(iOutput)
       massLogarithmic=log10(mass)
       do j=1,blackHoleDistributions(i)%descriptor%massSystematicCoefficientCount
          massLogarithmic         =+massLogarithmic                                                    &
               &                   +  blackHoleDistributions(i)%massSystematicCoefficients         (j) &
               &                   *(                                                                  &
               &                     +log10(mass)                                                      &
               &                     -blackHoleDistributions(i)%descriptor%massSystematicLogM0         &
               &                    )**(j-1)
       end do
       if (massLogarithmic <  blackHoleDistributions(i)%descriptor%massLogarithmicMinimum) cycle
       blackHoleMass           =blackHoleMass*blackHoleDistributions(i)%cosmologyConversionBlackHoleMass(iOutput)
       blackHoleMassLogarithmic=log10(blackHoleMass)
       do j=1,blackHoleDistributions(i)%descriptor%blackHoleMassSystematicCoefficientCount
          blackHoleMassLogarithmic=+  blackHoleMassLogarithmic                                         &
               &                   +  blackHoleDistributions(i)%blackHoleMassSystematicCoefficients(j) &
               &                   *(                                                                  &
               &                     +log10(mass       )                                               &
               &                     -blackHoleDistributions(i)%descriptor%massSystematicLogM0         &
               &                    )**(j-1)
       end do
       ! Compute random errors.
       massRandomError         =0.0d0
       blackHoleMassRandomError=0.0d0
       do j=1,blackHoleDistributions(i)%descriptor%massRandomCoefficientCount
          massRandomError         =+massRandomError                                              &
               &                   +  blackHoleDistributions(i)%massRandomCoefficients       (j) &
               &                   *(                                                            &
               &                     +log10(mass)                                                &
               &                     -blackHoleDistributions(i)%descriptor%massSystematicLogM0   &
               &                    )**(j-1)
       end do
       do j=1,blackHoleDistributions(i)%descriptor%blackHoleMassRandomCoefficientCount
          blackHoleMassRandomError=+blackHoleMassRandomError                                     &
               &                   +blackHoleDistributions(i)%blackHoleMassRandomCoefficients(j) &
               &                   *(                                                            &
               &                     +log10(mass)                                                &
               &                     -blackHoleDistributions(i)%descriptor%massSystematicLogM0   &
               &                    )**(j-1)
       end do
       massRandomError         =max(         massRandomError,         massRandomErrorMinimum)
       blackHoleMassRandomError=max(blackHoleMassRandomError,blackHoleMassRandomErrorMinimum)
       ! Compute weights for each bin.
       galaxyWork(i)%blackHoleMassDistributionWeights     =+(                                                                                                                                  &
            &                                                +erf((blackHoleDistributions(i)%         massesLogarithmicMaximum-         massLogarithmic)/         massRandomError/sqrt(2.0d0)) &
            &                                                -erf((blackHoleDistributions(i)%         massesLogarithmicMinimum-         massLogarithmic)/         massRandomError/sqrt(2.0d0)) &
            &                                               )                                                                                                                                  &
            &                                              /2.0d0                                                                                                                              &
            &                                              *tree%volumeWeight
       do j=1,blackHoleDistributions(i)%massesCount
          galaxyWork(i)%blackHoleMassDistribution    (:,j)=+(                                                                                                                                  &
               &                                             +erf((blackHoleDistributions(i)%blackHoleMassesLogarithmicMaximum-blackHoleMassLogarithmic)/blackHoleMassRandomError/sqrt(2.0d0)) &
               &                                             -erf((blackHoleDistributions(i)%blackHoleMassesLogarithmicMinimum-blackHoleMassLogarithmic)/blackHoleMassRandomError/sqrt(2.0d0)) &
               &                                            )                                                                                                                                  &
               &                                           /2.0d0                                                                                                                              &
               &                                           *galaxyWork(i)%blackHoleMassDistributionWeights(j)
       end do
       ! Apply output weights.
       do j=1,blackHoleDistributions(i)%blackHoleMassesCount
          galaxyWork                    (i)%blackHoleMassDistribution(j,:        ) &
               & =galaxyWork            (i)%blackHoleMassDistribution(j,:        ) &
               & *blackHoleDistributions(i)%outputWeight             (  :,iOutput)
       end do
       galaxyWork(i)%blackHoleMassDistributionWeights=galaxyWork(i)%blackHoleMassDistributionWeights*blackHoleDistributions(i)%outputWeight(:,iOutput)
       ! Accumulate blackHoleMass distribution.
       if (any(galaxyWork(i)%blackHoleMassDistribution /= 0.0d0)) then
          !$omp critical (Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins_Accumulate)
          blackHoleDistributions(i)%blackHoleMassDistribution       =blackHoleDistributions(i)%blackHoleMassDistribution       +galaxyWork(i)%blackHoleMassDistribution
          blackHoleDistributions(i)%blackHoleMassDistributionWeights=blackHoleDistributions(i)%blackHoleMassDistributionWeights+galaxyWork(i)%blackHoleMassDistributionWeights
          !$omp end critical (Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins_Accumulate)
          ! Treat main branch and other galaxies differently.
          if (node%isOnMainBranch().and.analysisBlackHoleMassDistributionCovarianceModel == analysisBlackHoleMassDistributionCovarianceModelBinomial) then
             ! Find the bin to which this halo mass belongs.
             basic => node%basic()
             haloMassBin=floor((log10(basic%mass())-analysisBlackHoleMassDistributionsHaloMassMinimumLogarithmic)*analysisBHMassDstrbtnsHaloMassIntervalLogarithmicInverse)+1
             ! Accumulate weights to halo mass arrays.
             if (haloMassBin >= 1 .and. haloMassBin <= analysisBlackHoleMassDistributionsHaloMassBinsCount) then
                !$omp critical (Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins_Accumulate)
                blackHoleDistributions        (i)%mainBranchGalaxyWeights       (:,:,haloMassBin)= &
                     &  blackHoleDistributions(i)%mainBranchGalaxyWeights       (:,:,haloMassBin)  &
                     &  +galaxyWork  (i)%blackHoleMassDistribution
                blackHoleDistributions        (i)%mainBranchGalaxyWeightsSquared(:,:,haloMassBin)= &
                     &  blackHoleDistributions(i)%mainBranchGalaxyWeightsSquared(:,:,haloMassBin)  &
                     &  +galaxyWork  (i)%blackHoleMassDistribution**2
                !$omp end critical (Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins_Accumulate)
             end if
          else
             !$omp critical (Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins_Accumulate)
             call Vector_Outer_Product_Accumulate(                                                                         &
                  &                               reshape(                                                                 &
                  &                                       galaxyWork(i)%blackHoleMassDistribution                        , &
                  &                                       [                                                                &
                  &                                        +blackHoleDistributions(i)%massesCount                          &
                  &                                        *blackHoleDistributions(i)%blackHoleMassesCount                 &
                  &                                       ]                                                                &
                  &                                      )                                                               , &
                  &                                         blackHoleDistributions(i)%blackHoleMassDistributionCovariance, &
                  &                               sparse=.true.                                                            &
                  &                              )
             !$omp end critical (Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins_Accumulate)
           end if
       end if
    end do
    return
  end subroutine Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins

  !# <hdfPreCloseTask>
  !#  <unitName>Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins_Output</unitName>
  !# </hdfPreCloseTask>
  subroutine Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins_Output
    !% Outputs blackHoleMass distributions to file.
    use Galacticus_HDF5
    use Vectors
    use Linear_Algebra
    use Memory_Management
    implicit none
    integer                      :: k,m,mi,zi,mj,zj,ci,cj
    type            (hdf5Object) :: analysisGroup,blackHoleMassDistributionGroup,dataset
    double precision             :: haloWeightBinTotal

    ! Return immediately if this analysis is not active.
    if (.not.analysisActive) return
    ! Iterate over mass functions.
    do k=1,size(blackHoleDistributions)
       ! Symmetrize the covariance matrix (we've accumulated only the upper triangle).
       blackHoleDistributions(k)%blackHoleMassDistributionCovariance=Matrix_Copy_Upper_To_Lower_Triangle(blackHoleDistributions(k)%blackHoleMassDistributionCovariance)
       ! Add the contribution from main branch galaxies to the covariance matrix.
       if (analysisBlackHoleMassDistributionCovarianceModel == analysisBlackHoleMassDistributionCovarianceModelBinomial) then
          do m=1,analysisBlackHoleMassDistributionsHaloMassBinsCount
             haloWeightBinTotal=sum(blackHoleDistributions(k)%mainBranchGalaxyWeights(:,:,m))
             if ( haloWeightBinTotal > 0.0 ) then
                do mi=1,blackHoleDistributions(k)%massesCount
                   do zi=1,blackHoleDistributions(k)%blackHoleMassesCount
                      ci=(mi-1)*blackHoleDistributions(k)%blackHoleMassesCount+zi
                      blackHoleDistributions               (k)%blackHoleMassDistributionCovariance(ci,ci  )=                    &
                           &         blackHoleDistributions(k)%blackHoleMassDistributionCovariance(ci,ci  )                     &
                           & +(1.0d0-blackHoleDistributions(k)%mainBranchGalaxyWeights            (zi,mi,m)/haloWeightBinTotal) &
                           & *       blackHoleDistributions(k)%mainBranchGalaxyWeightsSquared     (zi,mi,m)
                      do mj=1,blackHoleDistributions(k)%massesCount
                         do zj=1,blackHoleDistributions(k)%blackHoleMassesCount
                            cj=(mj-1)*blackHoleDistributions(k)%blackHoleMassesCount+zj
                            if (mi == mj .and. zi == zj) cycle
                            blackHoleDistributions         (k)%blackHoleMassDistributionCovariance(ci,cj  )=                    &
                                 &   blackHoleDistributions(k)%blackHoleMassDistributionCovariance(ci,cj  )                     &
                                 & -(blackHoleDistributions(k)%mainBranchGalaxyWeights            (zj,mj,m)/haloWeightBinTotal) &
                                 & * blackHoleDistributions(k)%mainBranchGalaxyWeightsSquared     (zi,mi,m)
                         end do
                      end do
                   end do
                end do
             end if
          end do
       end if
       ! Normalize the model blackHoleMass distribution in each mass interval and convert model black hole mass distribution to differential per log10(M).
       do mi=1,blackHoleDistributions(k)%massesCount 
          if (blackHoleDistributions(k)%blackHoleMassDistributionWeights(mi) > 0.0d0) then
             do zi=1,blackHoleDistributions(k)%blackHoleMassesCount
                ci=(mi-1)*blackHoleDistributions(k)%blackHoleMassesCount+zi
                blackHoleDistributions               (k)%blackHoleMassDistribution        (zi,mi)=blackHoleDistributions(k)%blackHoleMassDistribution        (zi,mi)  &
                     &       /(blackHoleDistributions(k)%blackHoleMassesLogarithmicMaximum(zi   )-blackHoleDistributions(k)%blackHoleMassesLogarithmicMinimum(zi   )) &
                     &       / blackHoleDistributions(k)%blackHoleMassDistributionWeights (   mi)
                do mj=1,blackHoleDistributions(k)%massesCount
                   if (blackHoleDistributions(k)%blackHoleMassDistributionWeights(mj) > 0.0d0) then
                      do zj=1,blackHoleDistributions(k)%blackHoleMassesCount 
                         cj=(mj-1)*blackHoleDistributions(k)%blackHoleMassesCount+zj
                         blackHoleDistributions         (k)%blackHoleMassDistributionCovariance(ci,cj)=blackHoleDistributions(k)%blackHoleMassDistributionCovariance(ci,cj)  &
                              & /(blackHoleDistributions(k)%blackHoleMassesLogarithmicMaximum  (zi   )-blackHoleDistributions(k)%blackHoleMassesLogarithmicMinimum  (zi   )) &
                              & /(blackHoleDistributions(k)%blackHoleMassesLogarithmicMaximum  (zj   )-blackHoleDistributions(k)%blackHoleMassesLogarithmicMinimum  (zj   )) &
                              & / blackHoleDistributions(k)%blackHoleMassDistributionWeights   (   mi)                                                                       &
                              & / blackHoleDistributions(k)%blackHoleMassDistributionWeights   (   mj)
                      end do
                   end if
                end do
             end do
          end if
       end do
       ! Output the blackHoleMass distribution.
       !$omp critical(HDF5_Access)
       analysisGroup                 =galacticusOutputFile%openGroup('analysis','Model analysis')
       blackHoleMassDistributionGroup=analysisGroup       %openGroup(trim(blackHoleDistributions(k)%descriptor%label),trim(blackHoleDistributions(k)%descriptor%comment))
       call blackHoleMassDistributionGroup%writeDataset  (blackHoleDistributions(k)%masses                             ,'mass'                               ,'Mass'                                   ,datasetReturned=dataset)
       call dataset                       %writeAttribute(blackHoleDistributions(k)%descriptor%massUnitsInSI           ,'unitsInSI'                                                                                            )
       call dataset                       %close()
       call blackHoleMassDistributionGroup%writeDataset  (blackHoleDistributions(k)%blackHoleMasses                    ,'blackHoleMass'                      ,'Black hole mass'                        ,datasetReturned=dataset)
       call dataset                       %writeAttribute(blackHoleDistributions(k)%descriptor%massUnitsInSI           ,'unitsInSI'                                                                                            )
       call dataset                       %close()
       call blackHoleMassDistributionGroup%writeDataset  (blackHoleDistributions(k)%blackHoleMassDistribution          ,'blackHoleMassDistribution'          ,'Black hole mass distribution'                                   )
       call blackHoleMassDistributionGroup%writeDataset  (blackHoleDistributions(k)%blackHoleMassDistributionCovariance,'blackHoleMassDistributionCovariance','Black hole mass distribution covariance'                        )
       call blackHoleMassDistributionGroup%close()
       call analysisGroup                 %close()
       !$omp end critical(HDF5_Access)
    end do
    return
  end subroutine Galacticus_Output_Analysis_Mass_Dpndnt_BH_Dstrbtins_Output

  double precision function Map_Black_Hole_Mass_Kormendy_Ho(blackHoleMass,node)
    !% Maps raw black hole masses.
    use Numerical_Constants_Astronomical
    implicit none
    double precision          , intent(in   )          :: blackHoleMass
    type            (treeNode), intent(inout), pointer :: node
    !GCC$ attributes unused :: node
    
    Map_Black_Hole_Mass_Kormendy_Ho=blackHoleMass
    return
  end function Map_Black_Hole_Mass_Kormendy_Ho

  double precision function Map_Galaxy_Mass_Kormendy_Ho(mass,node)
    !% Maps raw black hole masses.
    use Numerical_Constants_Astronomical
    implicit none
    double precision          , intent(in   )          :: mass
    type            (treeNode), intent(inout), pointer :: node
    !GCC$ attributes unused :: node

    Map_Galaxy_Mass_Kormendy_Ho=mass
    return
  end function Map_Galaxy_Mass_Kormendy_Ho
  
end module Galacticus_Output_Analyses_Mass_Dpndnt_BH_Dstrbtins
