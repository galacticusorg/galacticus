!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which performs analysis to compute a variety of mass-dependent metallicicity distributions.

module Galacticus_Output_Analyses_Mass_Dpndnt_Met_Dstrbtins
  !% Performs analysis to compute a variety of mass-dependent metallicity distributions. Currently supported metallicity distributions include:
  !% \begin{itemize}
  !% \item The \gls{sdss} stellar-phase distributions from \cite{gallazzi_ages_2005}.
  !% \item The \gls{sdss} gas-phase distributions from \cite{andrews_mass-metallicity_2013}.
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
  public :: Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins, Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins_Output
  
  ! Record of module initialization.
  logical                                                              :: moduleInitialized                     =.false.

  ! Record of whether this analysis is active.
  logical                                                              :: analysisActive

  ! Number of supported mass functions.
  integer          , parameter                                         :: metallicityDistributionsSupportedCount=2

  ! Labels for supported metallicity distributions.
  character(len=39), dimension(metallicityDistributionsSupportedCount) :: metallicityDistributionLabels         =        &
       & [                                                                                                               &
       &  'sdssStellarMetallicityDistributionZ0.07',                                                                     &
       &  'sdssGasMetallicityDistributionZ0.07    '                                                                      &
       & ]

  ! Interface for mass mapping functions.
  abstract interface
     double precision function Map_Mass(mass,thisNode)
       import treeNode
       double precision          , intent(in   )          :: mass
       type            (treeNode), intent(inout), pointer :: thisNode
     end function Map_Mass
  end interface

  ! Interface for metallicity mapping functions.
  abstract interface
     double precision function Map_Metallicity(radius,thisNode)
       import treeNode
       double precision          , intent(in   )          :: radius
       type            (treeNode), intent(inout), pointer :: thisNode
     end function Map_Metallicity
  end interface

  ! Interface for mass error functions.
  abstract interface
     double precision function Mass_Error(mass,thisNode)
       import treeNode
       double precision          , intent(in   )          :: mass
       type            (treeNode), intent(inout), pointer :: thisNode
     end function Mass_Error
  end interface

  ! Interface for metallicity error functions.
  abstract interface
     double precision function Metallicity_Error(radius,thisNode)
       import treeNode
       double precision          , intent(in   )          :: radius
       type            (treeNode), intent(inout), pointer :: thisNode
     end function Metallicity_Error
  end interface

  ! Options controlling type of metallicity statistics.
  integer         , parameter :: statisticDistribution=0
  integer         , parameter :: statisticMean        =1

  ! Type for descriptors of metallicity distributions.
  type :: metallicityDistributionDescriptor
     double precision                                           :: massSystematicLogM0                             , metallicitySystematicLogZ0
     procedure       (       Mass_Error  ), pointer    , nopass :: massRandomErrorFunction
     procedure       (Metallicity_Error  ), pointer    , nopass :: metallicityRandomErrorFunction
     double precision                                           :: massLogarithmicMinimum
     integer                                                    :: massSystematicCoefficientCount                  , metallicityMassSystematicCoefficientCount, &
          &                                                        massRandomCoefficientCount                      , metallicityRandomCoefficientCount        , &
          &                                                        metallicityMetallicitySystematicCoefficientCount
     integer                                                    :: metallicityType                                 , massType                                 , &
          &                                                        statisticType
     double precision                                           :: massUnitsInSI
     character       (len= 32            )                      :: label
     character       (len=128            )                      :: comment
     procedure       (Map_Mass           ), pointer    , nopass :: mapMass
     procedure       (Map_Metallicity    ), pointer    , nopass :: mapMetallicity
     class           (surveyGeometryClass), allocatable         :: geometry
     class           (galacticFilterClass), allocatable         :: filter
  end type metallicityDistributionDescriptor

  ! Metallicity distribution descriptors.
  type(metallicityDistributionDescriptor), dimension(metallicityDistributionsSupportedCount), target :: metallicityDistributionDescriptors= &
       & [                                                                                                    &
       ! SDSS galaxy stellar-phase metallicities from Gallazii et al. (2005). 
       &                           metallicityDistributionDescriptor(                                         &
       &                                                             11.0000d+0                             , &
       &                                                              0.0000d+0                             , &
       &                                                              null()                                , &
       &                                                              null()                                , &
       &                                                              8.000d0                               , &
       &                                                              2                                     , &
       &                                                              2                                     , &
       &                                                              2                                     , &
       &                                                              2                                     , &
       &                                                              2                                     , &
       &                                                              massTypeStellar                       , &
       &                                                              massTypeStellar                       , &
       &                                                              statisticDistribution                 , &
       &                                                              massSolar                             , &
       &                                                              'sdssStellarMetallicityZ0.07'         , &
       &                                                              'SDSS stellar metallicities at z=0.07', &
       &                                                              null()                                , &
       &                                                              null()                                , &
       &                                                              null()                                , &
       &                                                              null()                                  &
       &                                                            )                                       , &
       ! SDSS galaxy gas-phsae metallicities from Andrews & Martini (2013). 
       &                           metallicityDistributionDescriptor(                                         &
       &                                                             11.0000d+0                             , &
       &                                                              8.8600d+0                             , &
       &                                                              null()                                , &
       &                                                              null()                                , &
       &                                                              6.500d0                               , &
       &                                                              2                                     , &
       &                                                              2                                     , &
       &                                                              2                                     , &
       &                                                              2                                     , &
       &                                                              2                                     , &
       &                                                              massTypeGaseous                       , &
       &                                                              massTypeStellar                       , &
       &                                                              statisticMean                         , &
       &                                                              massSolar                             , &
       &                                                              'sdssGasMetallicityZ0.07'             , &
       &                                                              'SDSS gas metallicities at z=0.07'    , &
       &                                                              null()                                , &
       &                                                              null()                                , &
       &                                                              null()                                , &
       &                                                              null()                                  &
       &                                                            )                                         &
       & ]

  ! Type to store metallicity distributions.
  type :: metallicityDistribution
     ! Copy of the mass function descriptor for this mass function.
     type            (metallicityDistributionDescriptor), pointer                       :: descriptor
     ! Parameters for the systematic error model.
     double precision                                   , allocatable, dimension(:    ) :: massSystematicCoefficients                  , metallicityMassSystematicCoefficients , &
          &                                                                                massRandomCoefficients                      , metallicityRandomCoefficients         , &
          &                                                                                metallicityMetallicitySystematicCoefficients
     ! The number of bins.
     integer                                                                            :: massesCount                                 , metallicitiesCount
     ! Arrays for the masses, radii and size function.
     double precision                                   , allocatable, dimension(:    ) :: masses                                      , massesLogarithmic                     , &
          &                                                                                massesLogarithmicMinimum                    , massesLogarithmicMaximum              , &
          &                                                                                metallicities                               , metallicitiesLogarithmic              , &
          &                                                                                metallicitiesLogarithmicMinimum             , metallicitiesLogarithmicMaximum       , &
          &                                                                                metallicityDistributionWeights
     double precision                                   , allocatable, dimension(:,:  ) :: metallicityDistribution                     , outputWeight
     ! Arrays for accumulation of of main branch galaxies
     double precision                                   , allocatable, dimension(:,:,:) :: mainBranchGalaxyWeights                     , mainBranchGalaxyWeightsSquared
     ! Array for the covariance matrix.
     double precision                                   , allocatable, dimension(:,:  ) :: metallicityDistributionCovariance
     ! Cosmology conversion factors.
     double precision                                   , allocatable, dimension(:    ) :: cosmologyConversionMass
  end type metallicityDistribution

  ! Mass functions.
  type(metallicityDistribution), allocatable, dimension(:) :: metallicityDistributions

  ! Type for storing temporary metallicity functions during cumulation.
  type :: metallicityDistributionWork
     double precision, allocatable, dimension(:,:) :: metallicityDistribution
     double precision, allocatable, dimension(  :) :: metallicityDistributionWeights
  end type metallicityDistributionWork

  ! Work array.
  type(metallicityDistributionWork), allocatable, dimension(:) :: thisGalaxy
  !$omp threadprivate(thisGalaxy)
  
  ! Options controlling binning in halo mass.
  integer                     :: analysisMetallicityDistributionCovarianceModel
  integer         , parameter :: analysisMetallicityDistributionCovarianceModelPoisson    =1
  integer         , parameter :: analysisMetallicityDistributionCovarianceModelBinomial   =2
  integer                     :: analysisMetallicityDistributionsHaloMassBinsCount          , analysisMetallicityDistributionsHaloMassBinsPerDecade
  double precision            :: analysisMetallicityDistributionsHaloMassMinimum            , analysisMetallicityDistributionsHaloMassMaximum           , &
       &                         analysisMtllctyDstrbtnsHaloMassIntervalLogarithmicInverse  , analysisMetallicityDistributionsHaloMassMinimumLogarithmic

contains

  !# <mergerTreeAnalysisTask>
  !#  <unitName>Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins</unitName>
  !# </mergerTreeAnalysisTask>
  subroutine Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins(thisTree,thisNode,nodeStatus,iOutput,mergerTreeAnalyses)
    !% Construct metallicity distributions to compare to various observational determinations.
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
    use Abundances_Structure
    use Numerical_Ranges
    use Numerical_Constants_Astronomical
    implicit none
    type            (mergerTree                    ), intent(in   )                 :: thisTree
    type            (treeNode                      ), intent(inout), pointer        :: thisNode
    integer                                         , intent(in   )                 :: nodeStatus
    integer         (c_size_t                      ), intent(in   )                 :: iOutput
    type            (varying_string                ), intent(in   ), dimension(:  ) :: mergerTreeAnalyses
    double precision                                , allocatable  , dimension(:  ) :: metallicities
    class           (nodeComponentBasic            )               , pointer        :: thisBasic
    class           (nodeComponentDisk             )               , pointer        :: thisDisk
    class           (nodeComponentSpheroid         )               , pointer        :: thisSpheroid
    class           (cosmologyFunctionsClass       )               , pointer        :: cosmologyFunctionsModel
    type            (cosmologyFunctionsMatterLambda)                                :: cosmologyFunctionsObserved
    type            (cosmologyParametersSimple     )               , pointer        :: cosmologyParametersObserved
    integer                                         , parameter                     :: metallicityCountPerDecade    =10
    double precision                                , parameter                     :: massRandomErrorMinimum       =1.0d-3
    double precision                                , parameter                     :: metallicityRandomErrorMinimum=1.0d-3
    logical                                                                         :: massTypeStellarRequired, massTypeGaseousRequired
    integer         (c_size_t                      )                                :: k,jOutput
    integer                                                                         :: i,j,l,currentAnalysis,activeAnalysisCount,haloMassBin
    double precision                                                                :: dataHubbleParameter ,mass,massLogarithmic&
         &,massRandomError,metallicityLogarithmic,metallicity,metallicityRandomError,dataOmegaDarkEnergy,dataOmegaMatter,redshift,timeMinimum,timeMaximum,distanceMinimum,distanceMaximum,unitsInSI,metallicityMinimum,metallicityMaximum,metallicityMass
    type            (varying_string                )                                :: parameterName&
         &,analysisMetallicityDistributionCovarianceModelText,cosmologyScalingMass,scaling,message
    type            (hdf5Object                    )                                :: dataFile,massDataset,parametersGroup,metallicityDataset
    type            (abundances                    )                                :: abundancesDisk, abundancesSpheroid
    
    ! Initialize the module if necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins_Initialize)
       if (.not.moduleInitialized) then
          !@ <inputParameter>
          !@   <name>analysisMetallicityDistributionCovarianceModel</name>
          !@   <defaultValue>binomial</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The model to use when constructing the metallicity distribution covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisMetallicityDistributionCovarianceModel',analysisMetallicityDistributionCovarianceModelText,defaultValue='Poisson')
          select case (char(analysisMetallicityDistributionCovarianceModelText))
          case ( 'Poisson'  )
             analysisMetallicityDistributionCovarianceModel=analysisMetallicityDistributionCovarianceModelPoisson
          case ( 'binomial' )
             analysisMetallicityDistributionCovarianceModel=analysisMetallicityDistributionCovarianceModelBinomial
          case default
             call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins','unrecognized value for "analysisMetallicityDistributionCovarianceModel" - allowed values are "Poisson", and "binomial"')
          end select
          !@ <inputParameter>
          !@   <name>analysisMetallicityDistributionsHaloMassBinsPerDecade</name>
          !@   <defaultValue>10</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The number of bins per decade of halo mass to use when constructing the mass function covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisMetallicityDistributionsHaloMassBinsPerDecade',analysisMetallicityDistributionsHaloMassBinsPerDecade,defaultValue=10)
          !@ <inputParameter>
          !@   <name>analysisMetallicityDistributionsHaloMassMinimum</name>
          !@   <defaultValue>$10^8M_\odot$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The minimum halo mass to consider when constructing the mass function covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisMetallicityDistributionsHaloMassMinimum',analysisMetallicityDistributionsHaloMassMinimum,defaultValue=1.0d8)
          !@ <inputParameter>
          !@   <name>analysisMetallicityDistributionsHaloMassMaximum</name>
          !@   <defaultValue>$10^{16}M_\odot$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The maximum halo mass to consider when constructing the mass function covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisMetallicityDistributionsHaloMassMaximum',analysisMetallicityDistributionsHaloMassMaximum,defaultValue=1.0d16)
          analysisMetallicityDistributionsHaloMassMinimumLogarithmic        =     log10( analysisMetallicityDistributionsHaloMassMinimum)
          analysisMetallicityDistributionsHaloMassBinsCount                 =int(                                                             &
               &                                                                  log10(                                                      &
               &                                                                         analysisMetallicityDistributionsHaloMassMaximum      &
               &                                                                        /analysisMetallicityDistributionsHaloMassMinimum      &
               &                                                                       )                                                      &
               &                                                                 *dble(analysisMetallicityDistributionsHaloMassBinsPerDecade) &
               &                                                                 +0.5d0&
               &                                                                )
          analysisMtllctyDstrbtnsHaloMassIntervalLogarithmicInverse= dble (analysisMetallicityDistributionsHaloMassBinsCount)                 &
               &                                                             /log10(                                                          &
               &                                                                     analysisMetallicityDistributionsHaloMassMaximum          &
               &                                                                    /analysisMetallicityDistributionsHaloMassMinimum          &
               &                                                                   )
          ! Establish mapping functions for metallicity distribution descriptors.
          metallicityDistributionDescriptors(1)%mapMetallicity => Map_Metallicity_SDSS_Stellar_Phase_Z0_07
          metallicityDistributionDescriptors(2)%mapMetallicity => Map_Metallicity_SDSS_Gas_Phase_Z0_07
          ! Determine how many supported metallicity distributions are requested.
          activeAnalysisCount=0
          do i=1,metallicityDistributionsSupportedCount
             if (any(trim(mergerTreeAnalyses) == trim(metallicityDistributionLabels(i)))) activeAnalysisCount=activeAnalysisCount+1
          end do
          ! Allocate mass function arrays and populate with required data.
          if (activeAnalysisCount <= 0) then
             analysisActive=.false.
          else
             analysisActive=.true.
             ! Establish survey geometries.
             allocate(surveyGeometryLiWhite2009SDSS :: metallicityDistributionDescriptors(1)%geometry)
             select type (g => metallicityDistributionDescriptors(1)%geometry)
             type is (surveyGeometryLiWhite2009SDSS)
                g=surveyGeometryLiWhite2009SDSS(               )
             end select             
             allocate(surveyGeometryLiWhite2009SDSS :: metallicityDistributionDescriptors(2)%geometry)
             select type (g => metallicityDistributionDescriptors(2)%geometry)
             type is (surveyGeometryLiWhite2009SDSS)
                g=surveyGeometryLiWhite2009SDSS(0.020d0,0.250d0)
             end select
             ! Establish sample filters.
             allocate(galacticFilterAlways            :: metallicityDistributionDescriptors(1)%filter)
             select type (f => metallicityDistributionDescriptors(1)%filter)
             type is (galacticFilterAlways)
                f=galacticFilterAlways           (                  )
             end select             
             allocate(galacticFilterStarFormationRate :: metallicityDistributionDescriptors(2)%filter)
             select type (f => metallicityDistributionDescriptors(2)%filter)
             type is (galacticFilterStarFormationRate)
                f=galacticFilterStarFormationRate(10.0d0,8.0d0,1.0d0)
             end select
             ! Initialize analyses.
             currentAnalysis=0
             allocate(metallicityDistributions(activeAnalysisCount))
             cosmologyFunctionsModel => cosmologyFunctions()
             do i=1,size(mergerTreeAnalyses)
                do j=1,metallicityDistributionsSupportedCount
                   if (mergerTreeAnalyses(i) == trim(metallicityDistributionLabels(j))) then
                      currentAnalysis=currentAnalysis+1
                      ! Set a pointer to the descriptor for this size function.
                      metallicityDistributions(currentAnalysis)%descriptor => metallicityDistributionDescriptors(j)
                      ! Read parameters of the systematic error model.
                      if (metallicityDistributionDescriptors(j)%massSystematicCoefficientCount > 0) then
                         allocate(metallicityDistributions(currentAnalysis)%massSystematicCoefficients(metallicityDistributionDescriptors(j)%massSystematicCoefficientCount))
                         do k=1,metallicityDistributionDescriptors(j)%massSystematicCoefficientCount
                            parameterName=trim(metallicityDistributionLabels(j))//'MassSystematic'
                            parameterName=parameterName//(k-1)
                            !@ <inputParameter>
                            !@   <regEx>(sdssStellarMetallicityDistribution|sdssGasMetallicityDistribution)Z[0-9\.]+MassSystematic[0-9]+</regEx>
                            !@   <defaultValue>0</defaultValue>
                            !@   <attachedTo>module</attachedTo>
                            !@   <description>
                            !@     Mass-dependent metallicity distribution mass systematic parameters.
                            !@   </description>
                            !@   <type>real</type>
                            !@   <cardinality>1</cardinality>
                            !@ </inputParameter>
                            call Get_Input_Parameter(char(parameterName),metallicityDistributions(currentAnalysis)%massSystematicCoefficients(k),defaultValue=0.0d0)
                         end do
                      end if
                      if (metallicityDistributionDescriptors(j)%metallicityMassSystematicCoefficientCount > 0) then
                         allocate(metallicityDistributions(currentAnalysis)%metallicityMassSystematicCoefficients(metallicityDistributionDescriptors(j)%metallicityMassSystematicCoefficientCount))
                         do k=1,metallicityDistributionDescriptors(j)%metallicityMassSystematicCoefficientCount
                            parameterName=trim(metallicityDistributionLabels(j))//'MetallicityMassSystematic'
                            parameterName=parameterName//(k-1)
                            !@ <inputParameter>
                            !@   <regEx>(sdssStellarMetallicityDistribution|sdssGasMetallicityDistribution)Z[0-9\.]+MetallicityMassSystematic[0-9]+</regEx>
                            !@   <defaultValue>0</defaultValue>
                            !@   <attachedTo>module</attachedTo>
                            !@   <description>
                            !@     Mass-dependent metallicity distribution metallicity systematic parameters.
                            !@   </description>
                            !@   <type>real</type>
                            !@   <cardinality>1</cardinality>
                            !@ </inputParameter>
                            call Get_Input_Parameter(char(parameterName),metallicityDistributions(currentAnalysis)%metallicityMassSystematicCoefficients(k),defaultValue=0.0d0)
                         end do
                      end if
                      if (metallicityDistributionDescriptors(j)%metallicityMetallicitySystematicCoefficientCount > 0) then
                         allocate(metallicityDistributions(currentAnalysis)%metallicityMetallicitySystematicCoefficients(metallicityDistributionDescriptors(j)%metallicityMetallicitySystematicCoefficientCount))
                         do k=1,metallicityDistributionDescriptors(j)%metallicityMetallicitySystematicCoefficientCount
                            parameterName=trim(metallicityDistributionLabels(j))//'MetallicityMetallicitySystematic'
                            parameterName=parameterName//(k-1)
                            !@ <inputParameter>
                            !@   <regEx>(sdssStellarMetallicityDistribution|sdssGasMetallicityDistribution)Z[0-9\.]+MetallicityMetallicitySystematic[0-9]+</regEx>
                            !@   <defaultValue>0</defaultValue>
                            !@   <attachedTo>module</attachedTo>
                            !@   <description>
                            !@     Mass-dependent metallicity distribution metallicity systematic parameters.
                            !@   </description>
                            !@   <type>real</type>
                            !@   <cardinality>1</cardinality>
                            !@ </inputParameter>
                            call Get_Input_Parameter(char(parameterName),metallicityDistributions(currentAnalysis)%metallicityMetallicitySystematicCoefficients(k),defaultValue=0.0d0)
                         end do
                      end if
                      ! Read parameters of the random error model.
                      if (metallicityDistributionDescriptors(j)%massRandomCoefficientCount > 0) then
                         allocate(metallicityDistributions(currentAnalysis)%massRandomCoefficients(metallicityDistributionDescriptors(j)%massRandomCoefficientCount))
                         do k=1,metallicityDistributionDescriptors(j)%massRandomCoefficientCount
                            parameterName=trim(metallicityDistributionLabels(j))//'MassRandom'
                            parameterName=parameterName//(k-1)
                            !@ <inputParameter>
                            !@   <regEx>(sdssStellarMetallicityDistribution|sdssGasMetallicityDistribution)Z[0-9\.]+MassRandom[0-9]+</regEx>
                            !@   <defaultValue>0</defaultValue>
                            !@   <attachedTo>module</attachedTo>
                            !@   <description>
                            !@     Mass-dependent metallicity distribution mass random parameters.
                            !@   </description>
                            !@   <type>real</type>
                            !@   <cardinality>1</cardinality>
                            !@ </inputParameter>
                            call Get_Input_Parameter(char(parameterName),metallicityDistributions(currentAnalysis)%massRandomCoefficients(k),defaultValue=0.0d0)
                         end do
                      end if
                      ! Read parameters of the random error model.
                      if (metallicityDistributionDescriptors(j)%metallicityRandomCoefficientCount > 0) then
                         allocate(metallicityDistributions(currentAnalysis)%metallicityRandomCoefficients(metallicityDistributionDescriptors(j)%metallicityRandomCoefficientCount))
                         do k=1,metallicityDistributionDescriptors(j)%metallicityRandomCoefficientCount
                            parameterName=trim(metallicityDistributionLabels(j))//'MetallicityRandom'
                            parameterName=parameterName//(k-1)
                            !@ <inputParameter>
                            !@   <regEx>(sdssStellarMetallicityDistribution|sdssGasMetallicityDistribution)Z[0-9\.]+MetallicityRandom[0-9]+</regEx>
                            !@   <defaultValue>0</defaultValue>
                            !@   <attachedTo>module</attachedTo>
                            !@   <description>
                            !@     Mass-dependent metallicity distribution metallicity random parameters.
                            !@   </description>
                            !@   <type>real</type>
                            !@   <cardinality>1</cardinality>
                            !@ </inputParameter>
                            call Get_Input_Parameter(char(parameterName),metallicityDistributions(currentAnalysis)%metallicityRandomCoefficients(k),defaultValue=0.0d0)
                         end do
                      end if
                      ! Read the appropriate observational data definition.
                      select case (trim(metallicityDistributionLabels(j)))
                      case ('sdssStellarMetallicityDistributionZ0.07')
                         ! SDSS z=0.07 stellar-phase metallicity distribution.
                         !$omp critical(HDF5_Access)
                         call dataFile%openFile(char(Galacticus_Input_Path())//"data/observations/abundances/stellarPhaseMetallicityGallazzi2005.hdf5",readOnly=.true.)
                         ! Read masses.
                         call dataFile%readDataset("mass",metallicityDistributions(currentAnalysis)%masses)
                         massDataset=dataFile%openDataset("mass")
                         call massDataset%readAttribute("cosmologyScaling",cosmologyScalingMass)
                         call massDataset%readAttribute("scaling"         ,scaling             )
                         call massDataset%readAttribute("unitsInSI"       ,unitsInSI           )
                         call massDataset%close        (                                       )
                         select case (char(scaling))
                         case('linear')
                            ! Nothing to do in this case.
                         case('log10' )
                            ! Convert from log10 mass to linear mass.
                            metallicityDistributions(currentAnalysis)%masses=10.0d0**metallicityDistributions(currentAnalysis)%masses
                         case default
                            call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins','unrecognized scaling')
                         end select
                         metallicityDistributions        (currentAnalysis)%masses= &
                              & +metallicityDistributions(currentAnalysis)%masses  &
                              & *unitsInSI                                         &
                              & /massSolar
                         metallicityDistributions(currentAnalysis)%massesCount=size(metallicityDistributions(currentAnalysis)%masses)
                         ! Find range of metallicities.
                         metallicityMinimum=+huge(1.0d0)
                         metallicityMaximum=-huge(1.0d0)
                         call dataFile%readDataset("metallicityPercentile16",metallicities)
                         metallicityMinimum=min(metallicityMinimum,minval(metallicities))
                         metallicityMaximum=max(metallicityMaximum,maxval(metallicities))
                         call deallocateArray(metallicities)
                         call dataFile%readDataset("metallicityPercentile50",metallicities)
                         metallicityMinimum=min(metallicityMinimum,minval(metallicities))
                         metallicityMaximum=max(metallicityMaximum,maxval(metallicities))
                         call deallocateArray(metallicities)
                         call dataFile%readDataset("metallicityPercentile84",metallicities)
                         metallicityMinimum=min(metallicityMinimum,minval(metallicities))
                         metallicityMaximum=max(metallicityMaximum,maxval(metallicities))
                         call deallocateArray(metallicities)                         
                         metallicityDataset=dataFile%openDataset("metallicityPercentile50")
                         call metallicityDataset%readAttribute("scaling"  ,scaling  )
                         call metallicityDataset%readAttribute("unitsInSI",unitsInSI)
                         call metallicityDataset%close        (                     )
                        select case (char(scaling))
                         case('linear')
                            ! Nothing to do in this case.
                         case('log10' )
                            ! Convert from log10 metallicity to linear metallicity.
                            metallicityMinimum=10.0d0**metallicityMinimum
                            metallicityMaximum=10.0d0**metallicityMaximum
                         case default
                            call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins','unrecognized scaling')
                         end select
                         metallicityMinimum=+1.0d-1*metallicityMinimum*unitsInSI/metallicitySolar
                         metallicityMaximum=+1.0d+1*metallicityMaximum*unitsInSI/metallicitySolar
                         metallicityDistributions(currentAnalysis)%metallicitiesCount=int(log10(metallicityMaximum/metallicityMinimum)*metallicityCountPerDecade)+1
                         allocate(metallicityDistributions(currentAnalysis)%metallicities(metallicityDistributions(currentAnalysis)%metallicitiesCount))
                         metallicityDistributions(currentAnalysis)%metallicities=Make_Range(                                                              &
                              &                                                             metallicityMinimum                                          , &
                              &                                                             metallicityMaximum                                          , &
                              &                                                             metallicityDistributions(currentAnalysis)%metallicitiesCount, &
                              &                                                             rangeTypeLogarithmic                                          &
                              &                                                            )
                         call allocateArray(metallicityDistributions(currentAnalysis)%       massesLogarithmic         ,[                                                             metallicityDistributions(currentAnalysis)%massesCount                                                  ])
                         call allocateArray(metallicityDistributions(currentAnalysis)%       massesLogarithmicMinimum  ,[                                                             metallicityDistributions(currentAnalysis)%massesCount                                                  ])
                         call allocateArray(metallicityDistributions(currentAnalysis)%       massesLogarithmicMaximum  ,[                                                             metallicityDistributions(currentAnalysis)%massesCount                                                  ])
                         call allocateArray(metallicityDistributions(currentAnalysis)%metallicitiesLogarithmic         ,[metallicityDistributions(currentAnalysis)%metallicitiesCount                                                                                                        ])
                         call allocateArray(metallicityDistributions(currentAnalysis)%metallicitiesLogarithmicMinimum  ,[metallicityDistributions(currentAnalysis)%metallicitiesCount                                                                                                        ])
                         call allocateArray(metallicityDistributions(currentAnalysis)%metallicitiesLogarithmicMaximum  ,[metallicityDistributions(currentAnalysis)%metallicitiesCount                                                                                                        ])
                         call allocateArray(metallicityDistributions(currentAnalysis)%metallicityDistribution          ,[metallicityDistributions(currentAnalysis)%metallicitiesCount,metallicityDistributions(currentAnalysis)%massesCount                                                  ])
                         call allocateArray(metallicityDistributions(currentAnalysis)%metallicityDistributionWeights   ,[                                                             metallicityDistributions(currentAnalysis)%massesCount                                                  ])
                         call allocateArray(metallicityDistributions(currentAnalysis)%metallicityDistributionCovariance,[                                                                                                                                                                       &
                              &                                                                                        metallicityDistributions(currentAnalysis)%metallicitiesCount*metallicityDistributions(currentAnalysis)%massesCount,                                                    &
                              &                                                                                        metallicityDistributions(currentAnalysis)%metallicitiesCount*metallicityDistributions(currentAnalysis)%massesCount                                                     &
                              &                                                                                       ]                                                                                                                                                                       &
                              &          )
                         call allocateArray(metallicityDistributions(currentAnalysis)%mainBranchGalaxyWeights          ,[metallicityDistributions(currentAnalysis)%metallicitiesCount,metallicityDistributions(currentAnalysis)%massesCount,analysisMetallicityDistributionsHaloMassBinsCount])
                         call allocateArray(metallicityDistributions(currentAnalysis)%mainBranchGalaxyWeightsSquared   ,[metallicityDistributions(currentAnalysis)%metallicitiesCount,metallicityDistributions(currentAnalysis)%massesCount,analysisMetallicityDistributionsHaloMassBinsCount])
                         metallicityDistributions              (currentAnalysis)%       massesLogarithmic  &
                              & =log10(metallicityDistributions(currentAnalysis)%       masses)
                         metallicityDistributions              (currentAnalysis)%metallicitiesLogarithmic  &
                              & =log10(metallicityDistributions(currentAnalysis)%metallicities)
                         do k=1,metallicityDistributions(currentAnalysis)%metallicitiesCount
                            if (k ==                                                            1) then
                               metallicityDistributions                (currentAnalysis)%metallicitiesLogarithmicMinimum(k  )= &
                                    & +        metallicityDistributions(currentAnalysis)%metallicitiesLogarithmic       (k  )  &
                                    & -0.5d0*(                                                                                 &
                                    &         +metallicityDistributions(currentAnalysis)%metallicitiesLogarithmic       (k+1)  &
                                    &         -metallicityDistributions(currentAnalysis)%metallicitiesLogarithmic       (k  )  &
                                    &        )
                            else
                               metallicityDistributions                (currentAnalysis)%metallicitiesLogarithmicMinimum(k  )= &
                                    & +0.5d0*(                                                                                 &
                                    &         +metallicityDistributions(currentAnalysis)%metallicitiesLogarithmic       (k-1)  &
                                    &         +metallicityDistributions(currentAnalysis)%metallicitiesLogarithmic       (k  )  &
                                    &        )
                            end if
                            if (k == metallicityDistributions(currentAnalysis)%metallicitiesCount) then
                               metallicityDistributions                (currentAnalysis)%metallicitiesLogarithmicMaximum(k  )= &
                                    & +        metallicityDistributions(currentAnalysis)%metallicitiesLogarithmic       (k  )  &
                                    & +0.5d0*(                                                                                 &
                                    &         +metallicityDistributions(currentAnalysis)%metallicitiesLogarithmic       (k  )  &
                                    &         -metallicityDistributions(currentAnalysis)%metallicitiesLogarithmic       (k-1)  &
                                    &        )
                            else
                               metallicityDistributions                (currentAnalysis)%metallicitiesLogarithmicMaximum(k  )= &
                                    & +0.5d0*(                                                                                 &
                                    &         +metallicityDistributions(currentAnalysis)%metallicitiesLogarithmic       (k+1)  &
                                    &         +metallicityDistributions(currentAnalysis)%metallicitiesLogarithmic       (k  )  &
                                    &        )
                            end if
                         end do
                         do k=1,metallicityDistributions(currentAnalysis)%massesCount
                            if (k ==                                                            1) then
                               metallicityDistributions                (currentAnalysis)%massesLogarithmicMinimum(k  )= &
                                    & +        metallicityDistributions(currentAnalysis)%massesLogarithmic       (k  )  &
                                    & -0.5d0*(                                                                          &
                                    &         +metallicityDistributions(currentAnalysis)%massesLogarithmic       (k+1)  &
                                    &         -metallicityDistributions(currentAnalysis)%massesLogarithmic       (k  )  &
                                    &        )
                            else
                               metallicityDistributions                (currentAnalysis)%massesLogarithmicMinimum(k  )= &
                                    & +0.5d0*(                                                                          &
                                    &         +metallicityDistributions(currentAnalysis)%massesLogarithmic       (k-1)  &
                                    &         +metallicityDistributions(currentAnalysis)%massesLogarithmic       (k  )  &
                                    &        )
                            end if
                            if (k == metallicityDistributions(currentAnalysis)%massesCount) then
                               metallicityDistributions                (currentAnalysis)%massesLogarithmicMaximum(k  )= &
                                    & +        metallicityDistributions(currentAnalysis)%massesLogarithmic       (k  )  &
                                    & +0.5d0*(                                                                          &
                                    &         +metallicityDistributions(currentAnalysis)%massesLogarithmic       (k  )  &
                                    &         -metallicityDistributions(currentAnalysis)%massesLogarithmic       (k-1)  &
                                    &        )
                            else
                               metallicityDistributions                (currentAnalysis)%massesLogarithmicMaximum(k  )= &
                                    & +0.5d0*(                                                                          &
                                    &         +metallicityDistributions(currentAnalysis)%massesLogarithmic       (k+1)  &
                                    &         +metallicityDistributions(currentAnalysis)%massesLogarithmic       (k  )  &
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
                         metallicityDistributions(currentAnalysis)%metallicityDistribution          =0.0d0
                         metallicityDistributions(currentAnalysis)%metallicityDistributionWeights   =0.0d0
                         metallicityDistributions(currentAnalysis)%metallicityDistributionCovariance=0.0d0
                         metallicityDistributions(currentAnalysis)%mainBranchGalaxyWeights          =0.0d0
                         metallicityDistributions(currentAnalysis)%mainBranchGalaxyWeightsSquared   =0.0d0
                      case ('sdssGasMetallicityDistributionZ0.07')
                         ! SDSS z=0.07 gas-phase metallicity distribution.
                         !$omp critical(HDF5_Access)
                         call dataFile%openFile(char(Galacticus_Input_Path())//"data/observations/abundances/gasPhaseMetallicityAndrews2013.hdf5",readOnly=.true.)
                         ! Read masses.
                         call dataFile%readDataset("mass",metallicityDistributions(currentAnalysis)%masses)
                         massDataset=dataFile%openDataset("mass")
                         call massDataset%readAttribute("cosmologyScaling",cosmologyScalingMass)
                         call massDataset%readAttribute("scaling"         ,scaling             )
                         call massDataset%readAttribute("unitsInSI"       ,unitsInSI           )
                         call massDataset%close        (                                       )
                         select case (char(scaling))
                         case('linear')
                            ! Nothing to do in this case.
                         case('log10' )
                            ! Convert from log10 mass to linear mass.
                            metallicityDistributions(currentAnalysis)%masses=10.0d0**metallicityDistributions(currentAnalysis)%masses
                         case default
                            call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins','unrecognized scaling')
                         end select
                         metallicityDistributions        (currentAnalysis)%masses= &
                              & +metallicityDistributions(currentAnalysis)%masses  &
                              & *unitsInSI                                         &
                              & /massSolar
                         metallicityDistributions(currentAnalysis)%massesCount       =size(metallicityDistributions(currentAnalysis)%masses)
                         metallicityDistributions(currentAnalysis)%metallicitiesCount=2
                         call allocateArray(metallicityDistributions(currentAnalysis)%       massesLogarithmic         ,[                                                             metallicityDistributions(currentAnalysis)%massesCount                                                  ])
                         call allocateArray(metallicityDistributions(currentAnalysis)%       massesLogarithmicMinimum  ,[                                                             metallicityDistributions(currentAnalysis)%massesCount                                                  ])
                         call allocateArray(metallicityDistributions(currentAnalysis)%       massesLogarithmicMaximum  ,[                                                             metallicityDistributions(currentAnalysis)%massesCount                                                  ])
                         call allocateArray(metallicityDistributions(currentAnalysis)%metallicityDistribution          ,[metallicityDistributions(currentAnalysis)%metallicitiesCount,metallicityDistributions(currentAnalysis)%massesCount                                                  ])
                         call allocateArray(metallicityDistributions(currentAnalysis)%metallicityDistributionWeights   ,[                                                             metallicityDistributions(currentAnalysis)%massesCount                                                  ])
                         call allocateArray(metallicityDistributions(currentAnalysis)%metallicityDistributionCovariance,[                                                                                                                                                                       &
                              &                                                                                        metallicityDistributions(currentAnalysis)%metallicitiesCount*metallicityDistributions(currentAnalysis)%massesCount,                                                    &
                              &                                                                                        metallicityDistributions(currentAnalysis)%metallicitiesCount*metallicityDistributions(currentAnalysis)%massesCount                                                     &
                              &                                                                                       ]                                                                                                                                                                       &
                              &          )
                         call allocateArray(metallicityDistributions(currentAnalysis)%mainBranchGalaxyWeights          ,[metallicityDistributions(currentAnalysis)%metallicitiesCount,metallicityDistributions(currentAnalysis)%massesCount,analysisMetallicityDistributionsHaloMassBinsCount])
                         call allocateArray(metallicityDistributions(currentAnalysis)%mainBranchGalaxyWeightsSquared   ,[metallicityDistributions(currentAnalysis)%metallicitiesCount,metallicityDistributions(currentAnalysis)%massesCount,analysisMetallicityDistributionsHaloMassBinsCount])
                         metallicityDistributions              (currentAnalysis)%       massesLogarithmic  &
                              & =log10(metallicityDistributions(currentAnalysis)%       masses)
                         do k=1,metallicityDistributions(currentAnalysis)%massesCount
                            if (k ==                                                            1) then
                               metallicityDistributions                (currentAnalysis)%massesLogarithmicMinimum(k  )= &
                                    & +        metallicityDistributions(currentAnalysis)%massesLogarithmic       (k  )  &
                                    & -0.5d0*(                                                                          &
                                    &         +metallicityDistributions(currentAnalysis)%massesLogarithmic       (k+1)  &
                                    &         -metallicityDistributions(currentAnalysis)%massesLogarithmic       (k  )  &
                                    &        )
                            else
                               metallicityDistributions                (currentAnalysis)%massesLogarithmicMinimum(k  )= &
                                    & +0.5d0*(                                                                          &
                                    &         +metallicityDistributions(currentAnalysis)%massesLogarithmic       (k-1)  &
                                    &         +metallicityDistributions(currentAnalysis)%massesLogarithmic       (k  )  &
                                    &        )
                            end if
                            if (k == metallicityDistributions(currentAnalysis)%massesCount) then
                               metallicityDistributions                (currentAnalysis)%massesLogarithmicMaximum(k  )= &
                                    & +        metallicityDistributions(currentAnalysis)%massesLogarithmic       (k  )  &
                                    & +0.5d0*(                                                                          &
                                    &         +metallicityDistributions(currentAnalysis)%massesLogarithmic       (k  )  &
                                    &         -metallicityDistributions(currentAnalysis)%massesLogarithmic       (k-1)  &
                                    &        )
                            else
                               metallicityDistributions                (currentAnalysis)%massesLogarithmicMaximum(k  )= &
                                    & +0.5d0*(                                                                          &
                                    &         +metallicityDistributions(currentAnalysis)%massesLogarithmic       (k+1)  &
                                    &         +metallicityDistributions(currentAnalysis)%massesLogarithmic       (k  )  &
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
                         metallicityDistributions(currentAnalysis)%metallicityDistribution          =0.0d0
                         metallicityDistributions(currentAnalysis)%metallicityDistributionWeights   =0.0d0
                         metallicityDistributions(currentAnalysis)%metallicityDistributionCovariance=0.0d0
                         metallicityDistributions(currentAnalysis)%mainBranchGalaxyWeights          =0.0d0
                         metallicityDistributions(currentAnalysis)%mainBranchGalaxyWeightsSquared   =0.0d0
                      case default
                         call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins','unknown metallicity function')
                      end select
                      ! Get cosmological conversion factors.
                      call allocateArray(metallicityDistributions(currentAnalysis)%cosmologyConversionMass,[Galacticus_Output_Time_Count()])
                      do jOutput=1,Galacticus_Output_Time_Count()
                         redshift=                                                                                      &
                              &   cosmologyFunctionsModel %redshiftFromExpansionFactor(                                 &
                              &    cosmologyFunctionsModel%expansionFactor             (                                &
                              &                                                         Galacticus_Output_Time(jOutput) &
                              &                                                        )                                &
                              &                                                       )
                         call Cosmology_Conversion_Factors(                                                                                                       &
                              &                            redshift                                                                                             , &
                              &                            cosmologyFunctionsModel                                                                              , &
                              &                            cosmologyFunctionsObserved                                                                           , &
                              &                            cosmologyScalingMass      =cosmologyScalingMass                                                      , &
                              &                            cosmologyConversionMass   =metallicityDistributions(currentAnalysis)%cosmologyConversionMass(jOutput)  &
                              &                           )
                      end do
                      nullify(cosmologyParametersObserved)
                      ! Compute output weights for metallicity distribution.
                      call allocateArray(metallicityDistributions(currentAnalysis)%outputWeight,[int(metallicityDistributions(currentAnalysis)%massesCount,kind=c_size_t),Galacticus_Output_Time_Count()])
                      metallicityDistributions(currentAnalysis)%outputWeight=0.0d0
                      do k=1,metallicityDistributions(currentAnalysis)%massesCount
                         do jOutput=1,Galacticus_Output_Time_Count()
                            do l=1,metallicityDistributions(currentAnalysis)%descriptor%geometry%fieldCount()
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
                               distanceMinimum=max(                                                                                                                                      &
                                    &              cosmologyFunctionsModel%distanceComoving(timeMaximum)                                                                               , &
                                    &              metallicityDistributions(currentAnalysis)%descriptor%geometry%distanceMinimum(metallicityDistributions(currentAnalysis)%masses(k),l)  &
                                    &             )
                               distanceMaximum=min(                                                                                                                                      &
                                    &              cosmologyFunctionsModel%distanceComoving(timeMinimum)                                                                               , &
                                    &              metallicityDistributions(currentAnalysis)%descriptor%geometry%distanceMaximum(metallicityDistributions(currentAnalysis)%masses(k),l)  &
                                    &             )
                               metallicityDistributions        (currentAnalysis)%outputWeight                    (k,jOutput)  &
                                    & =metallicityDistributions(currentAnalysis)%outputWeight                    (k,jOutput)  &
                                    & +metallicityDistributions(currentAnalysis)%descriptor  %geometry%solidAngle(  l      )  &
                                    & /3.0d0                                                                                  &
                                    & *                                                                                       &
                                    & max(                                                                                    &
                                    &     +0.0d0                                                                            , &
                                    &     +distanceMaximum**3                                                                 &
                                    &     -distanceMinimum**3                                                                 &
                                    &    )
                            end do
                         end do
                         where(metallicityDistributions(currentAnalysis)%outputWeight(k,:) < 0.0d0)
                            metallicityDistributions(currentAnalysis)%outputWeight(k,:)=0.0d0
                         end where
                         if (any(metallicityDistributions(currentAnalysis)%outputWeight(k,:) > 0.0d0)) then
                            metallicityDistributions                  (currentAnalysis)%outputWeight(k,:)  &
                                 &       =    metallicityDistributions(currentAnalysis)%outputWeight(k,:)  &
                                 &       /sum(metallicityDistributions(currentAnalysis)%outputWeight(k,:))
                         else
                            message="metallicity distribution '"//trim(metallicityDistributions(currentAnalysis)%descriptor%label)//"' bin "
                            message=message//k//" has zero weights"
                            call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins',message)
                         end if
                      end do
                      ! Ensure that disk and spheroid components support relevant abundances property.
                      massTypeStellarRequired=.false.
                      massTypeGaseousRequired=.false.
                      select case (metallicityDistributions(currentAnalysis)%descriptor%       massType)
                      case (massTypeStellar)
                         massTypeStellarRequired=.true.
                      case (massTypeGaseous)
                         massTypeGaseousRequired=.true.
                      end select
                      select case (metallicityDistributions(currentAnalysis)%descriptor%metallicityType)
                      case (massTypeStellar)
                         massTypeStellarRequired=.true.
                      case (massTypeGaseous)
                         massTypeGaseousRequired=.true.
                      end select
                      if (massTypeStellarRequired.and..not.defaultDiskComponent    %abundancesStellarIsGettable())                         &
                           & call Galacticus_Error_Report                                                                                  &
                           & (                                                                                                             &
                           &  'Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins'                                                     , &
                           &  'This analysis requires that the "abundancesStellar" property of the disk is gettable.'//                    &
                           &  Galacticus_Component_List(                                                                                   &
                           &                            'disk'                                                                           , &
                           &                             defaultDiskComponent    %abundancesStellarAttributeMatch(requireGettable=.true.)  &
                           &                           )                                                                                   &
                           & )
                      if (massTypeGaseousRequired.and..not.defaultDiskComponent    %    abundancesGasIsGettable())                         &
                           & call Galacticus_Error_Report                                                                                  &
                           & (                                                                                                             &
                           &  'Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins'                                                     , &
                           &  'This analysis requires that the "abundancesGas" property of the disk is gettable.'//                        &
                           &  Galacticus_Component_List(                                                                                   &
                           &                            'disk'                                                                           , &
                           &                             defaultDiskComponent    %    abundancesGasAttributeMatch(requireGettable=.true.)  &
                           &                           )                                                                                   &
                           & )
                      if (massTypeStellarRequired.and..not.defaultSpheroidComponent%abundancesStellarIsGettable())                         &
                           & call Galacticus_Error_Report                                                                                  &
                           & (                                                                                                             &
                           &  'Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins'                                                     , &
                           &  'This analysis requires that the "abundancesStellar" property of the spheroid is gettable.'//                &
                           &  Galacticus_Component_List(                                                                                   &
                           &                            'spheroid'                                                                       , &
                           &                             defaultSpheroidComponent%abundancesStellarAttributeMatch(requireGettable=.true.)  &
                           &                           )                                                                                   &
                           & )
                      if (massTypeGaseousRequired.and..not.defaultSpheroidComponent%    abundancesGasIsGettable())                         &
                           & call Galacticus_Error_Report                                                                                  &
                           & (                                                                                                             &
                           &  'Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins'                                                     , &
                           &  'This analysis requires that the "abundancesGas" property of the spheroid is gettable.'//                    &
                           &  Galacticus_Component_List(                                                                                   &
                           &                            'spheroid'                                                                       , &
                           &                             defaultSpheroidComponent%    abundancesGasAttributeMatch(requireGettable=.true.)  &
                           &                           )                                                                                   &
                           & )
                      exit
                   end if
                end do
             end do
          end if
          ! Record that module is initialized.
          moduleInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins_Initialize)
    end if
    ! Return if this analysis is not active.
    if (.not.analysisActive                   ) return
    ! Return if this is a tree finalization.
    if (nodeStatus          == nodeStatusFinal) return
    ! Allocate work arrays.
    if (.not.allocated(thisGalaxy)) allocate(thisGalaxy(size(metallicityDistributions)))
    ! Iterate over active analyses.
    do i=1,size(metallicityDistributions)
       ! Cycle if this metallicity distribution receives no contribution from this output.
       if (all(metallicityDistributions(i)%outputWeight(:,iOutput) <= 0.0d0)) cycle
       ! Allocate workspace.
       if (.not.allocated(thisGalaxy(i)%metallicityDistribution       )) call allocateArray(thisGalaxy(i)%metallicityDistribution       ,[metallicityDistributions(i)%metallicitiesCount,metallicityDistributions(i)%massesCount])
       if (.not.allocated(thisGalaxy(i)%metallicityDistributionWeights)) call allocateArray(thisGalaxy(i)%metallicityDistributionWeights,[                                               metallicityDistributions(i)%massesCount])
       ! Filter the galaxy.
       if (.not.metallicityDistributions(i)%descriptor%filter%passes(thisNode)) cycle
       ! Get the galactic mass.
       mass=                                                                                                                                                       &
            &  Galactic_Structure_Enclosed_Mass(thisNode,radiusLarge,componentType=componentTypeDisk    ,massType=metallicityDistributions(i)%descriptor%massType) &
            & +Galactic_Structure_Enclosed_Mass(thisNode,radiusLarge,componentType=componentTypeSpheroid,massType=metallicityDistributions(i)%descriptor%massType)
       if (mass <= 0.0d0) cycle
       if (associated(metallicityDistributions(i)%descriptor%mapMass)) mass=metallicityDistributions(i)%descriptor%mapMass(mass,thisNode)
       ! Get the metal mass.
       thisDisk     => thisNode%disk    ()
       thisSpheroid => thisNode%spheroid()
       select case (metallicityDistributions(i)%descriptor%metallicityType)
       case (massTypeStellar)
          abundancesDisk    =thisDisk    %abundancesStellar()
          abundancesSpheroid=thisSpheroid%abundancesStellar()
       case (massTypeGaseous)
          abundancesDisk    =thisDisk    %abundancesGas    ()
          abundancesSpheroid=thisSpheroid%abundancesGas    ()
       end select
       metallicityMass=                                                                                                                                                   &
            &  Galactic_Structure_Enclosed_Mass(thisNode,radiusLarge,componentType=componentTypeDisk    ,massType=metallicityDistributions(i)%descriptor%metallicityType) &
            & +Galactic_Structure_Enclosed_Mass(thisNode,radiusLarge,componentType=componentTypeSpheroid,massType=metallicityDistributions(i)%descriptor%metallicityType)
       if (metallicityMass > 0.0d0) then
          metallicity=(abundancesDisk%metallicity()+abundancesSpheroid%metallicity())/metallicityMass
       else
          metallicity=0.0d0
       end if
       if (metallicity <= 0.0d0) cycle
       if (associated(metallicityDistributions(i)%descriptor%mapMetallicity)) metallicity=metallicityDistributions(i)%descriptor%mapMetallicity(metallicity,thisNode)
       ! Convert mass for cosmology and systematics.
       mass=mass*metallicityDistributions(i)%cosmologyConversionMass(iOutput)
       massLogarithmic=log10(mass)
       do j=1,metallicityDistributions(i)%descriptor%massSystematicCoefficientCount
          massLogarithmic=massLogarithmic+metallicityDistributions(i)%massSystematicCoefficients(j)*(log10(mass)-metallicityDistributions(i)%descriptor%massSystematicLogM0)**(j-1)
       end do
       if (massLogarithmic <  metallicityDistributions(i)%descriptor%massLogarithmicMinimum) cycle
       ! Convert metallicity for systematics.
       metallicityLogarithmic=log10(metallicity)
       do j=1,metallicityDistributions(i)%descriptor%metallicityMassSystematicCoefficientCount
          metallicityLogarithmic=+  metallicityLogarithmic                                                          &
               &                 +  metallicityDistributions(i)%metallicityMassSystematicCoefficients         (j  ) &
               &                 *(                                                                                 &
               &                   +log10(mass       )                                                              &
               &                   -metallicityDistributions(i)%descriptor%massSystematicLogM0                      &
               &                  )                                                                         **(j-1)
       end do
       do j=1,metallicityDistributions(i)%descriptor%metallicityMetallicitySystematicCoefficientCount
          metallicityLogarithmic=+  metallicityLogarithmic                                                          &
               &                 +  metallicityDistributions(i)%metallicityMetallicitySystematicCoefficients  (j  ) &
               &                 *(                                                                                 &
               &                   +log10(metallicity)                                                              &
               &                   -metallicityDistributions(i)%descriptor%metallicitySystematicLogZ0               &
               &                  )                                                                         **(j-1)
       end do
       ! Compute contributions to each bin.
       if (associated(metallicityDistributions(i)%descriptor%massRandomErrorFunction)) then
          massRandomError=metallicityDistributions(i)%descriptor%massRandomErrorFunction(mass,thisNode)
       else
          massRandomError=0.0d0
          do j=1,metallicityDistributions(i)%descriptor%massRandomCoefficientCount
             massRandomError=+massRandomError                                              &
                  &          +metallicityDistributions(i)%massRandomCoefficients(j)        &
                  &          *(                                                            &
                  &            +log10(mass)                                                &
                  &            -metallicityDistributions(i)%descriptor%massSystematicLogM0 &
                  &           )**(j-1)
          end do
          massRandomError=max(massRandomError,massRandomErrorMinimum)
       end if
       if (associated(metallicityDistributions(i)%descriptor%metallicityRandomErrorFunction)) then
          metallicityRandomError=metallicityDistributions(i)%descriptor%metallicityRandomErrorFunction(metallicity,thisNode)
       else         
          metallicityRandomError=0.0d0
          do j=1,metallicityDistributions(i)%descriptor%metallicityRandomCoefficientCount
             metallicityRandomError=+metallicityRandomError                                &
                  &          +metallicityDistributions(i)%metallicityRandomCoefficients(j) &
                  &          *(                                                            &
                  &            +log10(mass)                                                &
                  &            -metallicityDistributions(i)%descriptor%massSystematicLogM0 &
                  &           )**(j-1)
          end do
          metallicityRandomError=max(metallicityRandomError,metallicityRandomErrorMinimum)
       end if
       thisGalaxy(i)%metallicityDistributionWeights=(                                                                                                                              &
            &                                        +erf((metallicityDistributions(i)%massesLogarithmicMaximum       -       massLogarithmic)/       massRandomError/sqrt(2.0d0)) &
            &                                        -erf((metallicityDistributions(i)%massesLogarithmicMinimum       -       massLogarithmic)/       massRandomError/sqrt(2.0d0)) &
            &                                       )                                                                                                                              &
            &                                       /2.0d0                                                                                                                         &
            &                                       *thisTree%volumeWeight
       select case (metallicityDistributions(i)%descriptor%statisticType)
       case (statisticDistribution)
          do j=1,metallicityDistributions(i)%massesCount
             thisGalaxy(i)%metallicityDistribution(:,j)=+(                                                                                                                              &
                  &                                       +erf((metallicityDistributions(i)%metallicitiesLogarithmicMaximum-metallicityLogarithmic)/metallicityRandomError/sqrt(2.0d0)) &
                  &                                       -erf((metallicityDistributions(i)%metallicitiesLogarithmicMinimum-metallicityLogarithmic)/metallicityRandomError/sqrt(2.0d0)) &
                  &                                      )                                                                                                                              &
                  &                                     /2.0d0&
                  &                                     *thisGalaxy(i)%metallicityDistributionWeights(j)
          end do
       case (statisticMean        )
          thisGalaxy(i)%metallicityDistribution(1,:)=+metallicity*thisGalaxy(i)%metallicityDistributionWeights(:)
          thisGalaxy(i)%metallicityDistribution(2,:)=+            thisGalaxy(i)%metallicityDistributionWeights(:)
       end select
       ! Apply output weights.
       do j=1,metallicityDistributions(i)%metallicitiesCount
          thisGalaxy                      (i)%metallicityDistribution(j,:        ) &
               & =thisGalaxy              (i)%metallicityDistribution(j,:        ) &
               & *metallicityDistributions(i)%outputWeight           (  :,iOutput)
       end do
       thisGalaxy(i)%metallicityDistributionWeights=thisGalaxy(i)%metallicityDistributionWeights*metallicityDistributions(i)%outputWeight(:,iOutput)
       ! Accumulate metallicity distribution.
       if (any(thisGalaxy(i)%metallicityDistribution /= 0.0d0)) then
          !$omp critical (Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins_Accumulate)
          metallicityDistributions(i)%metallicityDistribution       =metallicityDistributions(i)%metallicityDistribution       +thisGalaxy(i)%metallicityDistribution
          metallicityDistributions(i)%metallicityDistributionWeights=metallicityDistributions(i)%metallicityDistributionWeights+thisGalaxy(i)%metallicityDistributionWeights
          !$omp end critical (Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins_Accumulate)
          ! Treat main branch and other galaxies differently.
          if (thisNode%isOnMainBranch().and.analysisMetallicityDistributionCovarianceModel == analysisMetallicityDistributionCovarianceModelBinomial) then
             ! Find the bin to which this halo mass belongs.
             thisBasic => thisNode%basic()
             haloMassBin=floor((log10(thisBasic%mass())-analysisMetallicityDistributionsHaloMassMinimumLogarithmic)*analysisMtllctyDstrbtnsHaloMassIntervalLogarithmicInverse)+1
             ! Accumulate weights to halo mass arrays.
             if (haloMassBin >= 1 .and. haloMassBin <= analysisMetallicityDistributionsHaloMassBinsCount) then
                !$omp critical (Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins_Accumulate)
                metallicityDistributions        (i)%mainBranchGalaxyWeights       (:,:,haloMassBin)= &
                     &  metallicityDistributions(i)%mainBranchGalaxyWeights       (:,:,haloMassBin)  &
                     &  +thisGalaxy  (i)%metallicityDistribution
                metallicityDistributions        (i)%mainBranchGalaxyWeightsSquared(:,:,haloMassBin)= &
                     &  metallicityDistributions(i)%mainBranchGalaxyWeightsSquared(:,:,haloMassBin)  &
                     &  +thisGalaxy  (i)%metallicityDistribution**2
                !$omp end critical (Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins_Accumulate)
             end if
          else
             !$omp critical (Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins_Accumulate)
             call Vector_Outer_Product_Accumulate(                                                                         &
                  &                               reshape(                                                                 &
                  &                                       thisGalaxy(i)%metallicityDistribution                          , &
                  &                                       [                                                                &
                  &                                        +metallicityDistributions(i)%massesCount                        &
                  &                                        *metallicityDistributions(i)%metallicitiesCount                 &
                  &                                       ]                                                                &
                  &                                      )                                                               , &
                  &                                         metallicityDistributions(i)%metallicityDistributionCovariance, &
                  &                               sparse=.true.                                                            &
                  &                              )
             !$omp end critical (Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins_Accumulate)
           end if
       end if
    end do
    return
  end subroutine Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins

  !# <hdfPreCloseTask>
  !#  <unitName>Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins_Output</unitName>
  !# </hdfPreCloseTask>
  subroutine Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins_Output
    !% Outputs metallicity distributions to file.
    use Galacticus_HDF5
    use Vectors
    use Linear_Algebra
    use Memory_Management
    implicit none
    double precision            , allocatable, dimension(:  ) :: metallicityMean
    double precision            , allocatable, dimension(:,:) :: metallicityMeanCovariance, jacobian
    integer                                                   :: k,m,mi,zi,mj,zj,ci,cj
    type            (hdf5Object)                              :: analysisGroup,metallicityDistributionGroup,thisDataset
    double precision                                          :: haloWeightBinTotal
    type            (matrix    )                              :: jacobianMatrix, covarianceMatrix

    ! Return immediately if this analysis is not active.
    if (.not.analysisActive) return
    ! Iterate over mass functions.
    do k=1,size(metallicityDistributions)
       ! Symmetrize the covariance matrix (we've accumulated only the upper triangle).
       metallicityDistributions(k)%metallicityDistributionCovariance=Matrix_Copy_Upper_To_Lower_Triangle(metallicityDistributions(k)%metallicityDistributionCovariance)
       ! Add the contribution from main branch galaxies to the covariance matrix.
       if (analysisMetallicityDistributionCovarianceModel == analysisMetallicityDistributionCovarianceModelBinomial) then
          do m=1,analysisMetallicityDistributionsHaloMassBinsCount
             haloWeightBinTotal=sum(metallicityDistributions(k)%mainBranchGalaxyWeights(:,:,m))
             if ( haloWeightBinTotal > 0.0 ) then
                do mi=1,metallicityDistributions(k)%massesCount
                   do zi=1,metallicityDistributions(k)%metallicitiesCount
                      ci=(mi-1)*metallicityDistributions(k)%metallicitiesCount+zi
                      metallicityDistributions               (k)%metallicityDistributionCovariance        (ci,ci  )=                    &
                           &         metallicityDistributions(k)%metallicityDistributionCovariance        (ci,ci  )                     &
                           & +(1.0d0-metallicityDistributions(k)%mainBranchGalaxyWeights                  (zi,mi,m)/haloWeightBinTotal) &
                           & *       metallicityDistributions(k)%mainBranchGalaxyWeightsSquared           (zi,mi,m)
                      do mj=1,metallicityDistributions(k)%massesCount
                         do zj=1,metallicityDistributions(k)%metallicitiesCount
                            cj=(mj-1)*metallicityDistributions(k)%metallicitiesCount+zj
                            if (mi == mj .and. zi == zj) cycle
                            metallicityDistributions         (k)%metallicityDistributionCovariance        (ci,cj  )=                    &
                                 &   metallicityDistributions(k)%metallicityDistributionCovariance        (ci,cj  )                     &
                                 & -(metallicityDistributions(k)%mainBranchGalaxyWeights                  (zj,mj,m)/haloWeightBinTotal) &
                                 & * metallicityDistributions(k)%mainBranchGalaxyWeightsSquared           (zi,mi,m)
                         end do
                      end do
                   end do
                end do
             end if
          end do
       end if
       select case (metallicityDistributions(k)%descriptor%statisticType)
       case (statisticDistribution)
          ! Normalize the model metallicity distribution in each mass interval and convert model metallicity distribution to differential per log10(Z).
          do mi=1,metallicityDistributions(k)%massesCount 
             if (metallicityDistributions(k)%metallicityDistributionWeights(mi) > 0.0d0) then
                do zi=1,metallicityDistributions(k)%metallicitiesCount
                   ci=(mi-1)*metallicityDistributions(k)%metallicitiesCount+zi
                   metallicityDistributions               (k)%metallicityDistribution        (zi,mi)=metallicityDistributions(k)%metallicityDistribution        (zi,mi)  &
                        &       /(metallicityDistributions(k)%metallicitiesLogarithmicMaximum(zi   )-metallicityDistributions(k)%metallicitiesLogarithmicMinimum(zi   )) &
                        &       / metallicityDistributions(k)%metallicityDistributionWeights (   mi)
                   do mj=1,metallicityDistributions(k)%massesCount
                      if (metallicityDistributions(k)%metallicityDistributionWeights(mj) > 0.0d0) then
                         do zj=1,metallicityDistributions(k)%metallicitiesCount 
                            cj=(mj-1)*metallicityDistributions(k)%metallicitiesCount+zj
                            metallicityDistributions         (k)%metallicityDistributionCovariance(ci,cj)=metallicityDistributions(k)%metallicityDistributionCovariance(ci,cj)  &
                                 & /(metallicityDistributions(k)%metallicitiesLogarithmicMaximum  (zi   )-metallicityDistributions(k)%metallicitiesLogarithmicMinimum  (zi   )) &
                                 & /(metallicityDistributions(k)%metallicitiesLogarithmicMaximum  (zj   )-metallicityDistributions(k)%metallicitiesLogarithmicMinimum  (zj   )) &
                                 & / metallicityDistributions(k)%metallicityDistributionWeights   (   mi)                                                                       &
                                 & / metallicityDistributions(k)%metallicityDistributionWeights   (   mj)
                         end do
                      end if
                   end do
                end do
             end if
          end do
       case (statisticMean        )
          ! Compute mean metallicity.
          allocate(metallicityMean          (metallicityDistributions(k)%massesCount                                          ))
          allocate(metallicityMeanCovariance(metallicityDistributions(k)%massesCount,  metallicityDistributions(k)%massesCount))
          allocate(jacobian                 (metallicityDistributions(k)%massesCount,2*metallicityDistributions(k)%massesCount))
          ! Compute mean metallicity.
          where (metallicityDistributions(k)%metallicityDistribution(2,:) > 0.0d0)
             metallicityMean=                                                 &
                  & +metallicityDistributions(k)%metallicityDistribution(1,:) &
                  & /metallicityDistributions(k)%metallicityDistribution(2,:)
          elsewhere
             metallicityMean=0.0d0
          end where
          ! Construct Jacobian.
          jacobian=0.0d0
          do m=1,metallicityDistributions(k)%massesCount
             if (metallicityDistributions(k)%metallicityDistribution(2,m) > 0.0d0) then
                jacobian(m,2*(m-1)+1)=+1.0d0                                                   /metallicityDistributions(k)%metallicityDistribution(2,m)
                jacobian(m,2*(m-1)+2)=-metallicityDistributions(k)%metallicityDistribution(1,m)/metallicityDistributions(k)%metallicityDistribution(2,m)**2
             else
                jacobian(m,2*(m-1)+1)=+0.0d0
                jacobian(m,2*(m-1)+2)=+0.0d0
             end if
          end do
          ! Calculate covariance of mean metallicity.
          jacobianMatrix           =jacobian
          covarianceMatrix         =metallicityDistributions(k)%metallicityDistributionCovariance          
          metallicityMeanCovariance=jacobianMatrix*(covarianceMatrix*jacobianMatrix%transpose())
          ! Assign mean metallicities.
          call deallocateArray(metallicityDistributions(k)%metallicityDistribution                                                                                                          )
          call deallocateArray(metallicityDistributions(k)%metallicityDistributionCovariance                                                                                                )
          call allocateArray  (metallicityDistributions(k)%metallicityDistribution          ,[1                                      ,metallicityDistributions(k)%massesCount])
          call allocateArray  (metallicityDistributions(k)%metallicityDistributionCovariance,[metallicityDistributions(k)%massesCount,metallicityDistributions(k)%massesCount])
          metallicityDistributions(k)%metallicityDistribution          (1,:)=metallicityMean
          metallicityDistributions(k)%metallicityDistributionCovariance(:,:)=metallicityMeanCovariance
          deallocate(metallicityMean          )
          deallocate(metallicityMeanCovariance)
          deallocate(jacobian                 )
       end select
       ! Output the metallicity distribution.
       !$omp critical(HDF5_Access)
       analysisGroup               =galacticusOutputFile%openGroup('analysis','Model analysis')
       metallicityDistributionGroup=analysisGroup       %openGroup(trim(metallicityDistributions(k)%descriptor%label),trim(metallicityDistributions(k)%descriptor%comment))
       call        metallicityDistributionGroup%writeDataset  (metallicityDistributions(k)%masses                           ,'mass'                             ,'Mass'                               ,datasetReturned=thisDataset)
       call        thisDataset                 %writeAttribute(metallicityDistributions(k)%descriptor%massUnitsInSI         ,'unitsInSI'                                                                                          )
       call        thisDataset                 %close()
       if (allocated(metallicityDistributions(k)%metallicities)) &
            & call metallicityDistributionGroup%writeDataset  (metallicityDistributions(k)%metallicities                    ,'metallicity'                      ,'Metallicity'                                                    )
       call        metallicityDistributionGroup%writeDataset  (metallicityDistributions(k)%metallicityDistribution          ,'metallicityDistribution'          ,'Metallicity distribution'                                       )
       call        metallicityDistributionGroup%writeDataset  (metallicityDistributions(k)%metallicityDistributionCovariance,'metallicityDistributionCovariance','Metallicity distribution covariance'                            )
       call        metallicityDistributionGroup%close()
       call        analysisGroup               %close()
       !$omp end critical(HDF5_Access)
    end do
    return
  end subroutine Galacticus_Output_Analysis_Mass_Dpndnt_Met_Dstrbtins_Output

  double precision function Map_Metallicity_SDSS_Stellar_Phase_Z0_07(metallicity,thisNode)
    !% Maps raw metallicities into Solar units.
    use Numerical_Constants_Astronomical
    implicit none
    double precision          , intent(in   )          :: metallicity
    type            (treeNode), intent(inout), pointer :: thisNode
    !GCC$ attributes unused :: thisNode

    Map_Metallicity_SDSS_Stellar_Phase_Z0_07=metallicity/metallicitySolar
    return
  end function Map_Metallicity_SDSS_Stellar_Phase_Z0_07

  double precision function Map_Metallicity_SDSS_Gas_Phase_Z0_07(metallicity,thisNode)
    !% Maps raw metallicities into 12+log(O/H) units. \cite{andrews_mass-metallicity_2013}
    !% assume 12+log(O/H)=8.86 for the Solar oxygen abundance.
    use Numerical_Constants_Astronomical
    implicit none
    double precision          , intent(in   )          :: metallicity
    type            (treeNode), intent(inout), pointer :: thisNode
    !GCC$ attributes unused :: thisNode
    
    Map_Metallicity_SDSS_Gas_Phase_Z0_07=+10.0d0**8.86d0   &
         &                               *metallicity      &
         &                               /metallicitySolar
    return
  end function Map_Metallicity_SDSS_Gas_Phase_Z0_07
  
end module Galacticus_Output_Analyses_Mass_Dpndnt_Met_Dstrbtins
