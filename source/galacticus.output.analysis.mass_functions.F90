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

!% Contains a module which performs analysis to compute a variety of mass functions.

module Galacticus_Output_Analyses_Mass_Functions
  !% Performs analysis to compute a variety of mass functions. Currently supported mass functions include:
  !% \begin{itemize}
  !% \item The $z\approx 0.07$ stellar mass function in the SDSS measured by \cite{li_distribution_2009}. Assumes $0.07$
  !% dex random errors on stellar masses. This is approximate, but motivated by the discussion of
  !% \cite{behroozi_comprehensive_2010}.
  !% \end{itemize}
  use Galacticus_Nodes
  use Galactic_Structure_Options
  implicit none
  private
  public :: Galacticus_Output_Analysis_Mass_Functions, Galacticus_Output_Analysis_Mass_Functions_Output

  ! Record of module initialization.
  logical                                                   :: moduleInitialized=.false.

  ! Record of whether this analysis is active.
  logical                                                   :: analysisActive

  ! Number of supported mass functions.
  integer          , parameter                              :: massFunctionsSupportedCount=9

  ! Labels for supported mass functions.
  character(len=32), dimension(massFunctionsSupportedCount) :: massFunctionLabels= &
       & [                                                                         &
       &  'sdssStellarMassFunctionZ0.07   ',                                       &
       &  'alfalfaHiMassFunctionZ0.00     ',                                       &
       &  'primusStellarMassFunctionZ0.100',                                       &
       &  'primusStellarMassFunctionZ0.250',                                       &
       &  'primusStellarMassFunctionZ0.350',                                       &
       &  'primusStellarMassFunctionZ0.450',                                       &
       &  'primusStellarMassFunctionZ0.575',                                       &
       &  'primusStellarMassFunctionZ0.725',                                       &
       &  'primusStellarMassFunctionZ0.900'                                        &
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

  ! Type for descriptors of mass functions.
  type :: massFunctionDescriptor
     double precision                      :: redshift
     double precision                      :: systematicLogM0
     double precision                      :: randomError
     procedure       (Mass_Error), pointer :: randomErrorFunction
     double precision                      :: massLogarithmicMinimum
     integer                               :: systematicCoefficientCount
     integer                               :: massType
     character       (len= 32   )          :: label
     character       (len=128   )          :: comment
     procedure       (Map_Mass  ), pointer :: mapMass
  end type massFunctionDescriptor

  ! Mass function descriptors.
  type(massFunctionDescriptor), dimension(massFunctionsSupportedCount), target :: massFunctionDescriptors=    &
       & [ &
       &                           massFunctionDescriptor(                                                    &
       &                                                   0.07d00                                   ,        &
       &                                                  11.300d0                                   ,        &
       &                                                   0.070d0                                   ,        &
       &                                                   null()                                    ,        &
       &                                                   6.500d0                                   ,        &
       &                                                   2                                         ,        &
       &                                                   massTypeStellar                           ,        &
       &                                                   'sdssStellarMassFunctionZ0.07'            ,        &
       &                                                   'SDSS stellar mass function at z=0.07'    ,        &
       &                                                   null()                                             &
       &                                                 )                                                  , &
       ! ALFALFA survey. Note that HI/total gas mass fraction must be taken into account by the systematic errors model.
       &                           massFunctionDescriptor(                                                    &
       &                                                   0.000d0                                   ,        &
       &                                                   9.000d0                                   ,        &
       &                                                   0.000d0                                   ,        &
       &                                                   null()                                    ,        &
       &                                                   4.500d0                                   ,        &
       &                                                   1                                         ,        &
       &                                                   massTypeGaseous                           ,        &
       &                                                   'alfalfaHiMassFunctionZ0.00'              ,        &
       &                                                   'ALFALFA HI mass function at z=0.00'      ,        &
       &                                                   null()                                             &
       &                                                 )                                           ,        &
       &                           massFunctionDescriptor(                                                    &
       &                                                   0.100d0                                   ,        &
       &                                                  11.300d0                                   ,        &
       &                                                   0.070d0                                   ,        &
       &                                                   null()                                    ,        &
       &                                                   6.500d0                                   ,        &
       &                                                   2                                         ,        &
       &                                                   massTypeStellar                           ,        &
       &                                                   'primusStellarMassFunctionZ0.100'         ,        &
       &                                                   'PRMIUS stellar mass function at z=0.100' ,        &
       &                                                   null()                                             &
       &                                                 )                                                  , &
       &                           massFunctionDescriptor(                                                    &
       &                                                   0.250d0                                   ,        &
       &                                                  11.300d0                                   ,        &
       &                                                   0.070d0                                   ,        &
       &                                                   null()                                    ,        &
       &                                                   6.500d0                                   ,        &
       &                                                   2                                         ,        &
       &                                                   massTypeStellar                           ,        &
       &                                                   'primusStellarMassFunctionZ0.250'         ,        &
       &                                                   'PRMIUS stellar mass function at z=0.250' ,        &
       &                                                   null()                                             &
       &                                                 )                                                  , &
       &                           massFunctionDescriptor(                                                    &
       &                                                   0.350d0                                   ,        &
       &                                                  11.300d0                                   ,        &
       &                                                   0.070d0                                   ,        &
       &                                                   null()                                    ,        &
       &                                                   6.500d0                                   ,        &
       &                                                   2                                         ,        &
       &                                                   massTypeStellar                           ,        &
       &                                                   'primusStellarMassFunctionZ0.350'         ,        &
       &                                                   'PRMIUS stellar mass function at z=0.350' ,        &
       &                                                   null()                                             &
       &                                                 )                                                  , &
       &                           massFunctionDescriptor(                                                    &
       &                                                   0.450d0                                   ,        &
       &                                                  11.300d0                                   ,        &
       &                                                   0.070d0                                   ,        &
       &                                                   null()                                    ,        &
       &                                                   6.500d0                                   ,        &
       &                                                   2                                         ,        &
       &                                                   massTypeStellar                           ,        &
       &                                                   'primusStellarMassFunctionZ0.450'         ,        &
       &                                                   'PRMIUS stellar mass function at z=0.450' ,        &
       &                                                   null()                                             &
       &                                                 )                                                  , &
       &                           massFunctionDescriptor(                                                    &
       &                                                   0.575d0                                   ,        &
       &                                                  11.300d0                                   ,        &
       &                                                   0.070d0                                   ,        &
       &                                                   null()                                    ,        &
       &                                                   6.500d0                                   ,        &
       &                                                   2                                         ,        &
       &                                                   massTypeStellar                           ,        &
       &                                                   'primusStellarMassFunctionZ0.575'         ,        &
       &                                                   'PRMIUS stellar mass function at z=0.575' ,        &
       &                                                   null()                                             &
       &                                                 )                                                  , &
       &                           massFunctionDescriptor(                                                    &
       &                                                   0.725d0                                   ,        &
       &                                                  11.300d0                                   ,        &
       &                                                   0.070d0                                   ,        &
       &                                                   null()                                    ,        &
       &                                                   6.500d0                                   ,        &
       &                                                   2                                         ,        &
       &                                                   massTypeStellar                           ,        &
       &                                                   'primusStellarMassFunctionZ0.725'         ,        &
       &                                                   'PRMIUS stellar mass function at z=0.725' ,        &
       &                                                   null()                                             &
       &                                                 )                                                  , &
       &                           massFunctionDescriptor(                                                    &
       &                                                   0.90d0                                    ,        &
       &                                                  11.30d0                                    ,        &
       &                                                   0.07d0                                    ,        &
       &                                                   null()                                    ,        &
       &                                                   6.50d0                                    ,        &
       &                                                   2                                         ,        &
       &                                                   massTypeStellar                           ,        &
       &                                                   'primusStellarMassFunctionZ0.900'         ,        &
       &                                                   'PRMIUS stellar mass function at z=0.900' ,        &
       &                                                   null()                                             &
       &                                                 )                                                    &
       & ]

  ! Type to store mass functions.
  type :: massFunction
     ! Copy of the mass function descriptor for this mass function.
     type            (massFunctionDescriptor), pointer                     :: descriptor
     ! Parameters for the systematic error model.
     double precision                        , allocatable, dimension(:  ) :: systematicCoefficients
     ! The index of the output corresponding to the required redshift.
     integer                                                               :: outputNumber
     ! The number of mass bins.
     integer                                                               :: massesCount
     ! Arrays for the masses and mass function.
     double precision                        , allocatable, dimension(:  ) :: masses                  , massesLogarithmic              , &
          &                                                                   massesLogarithmicMinimum, massesLogarithmicMaximum       , &
          &                                                                   massFunction
     ! Arrays for accumulation of of main branch galaxies
     double precision                        , allocatable, dimension(:,:) :: mainBranchGalaxyWeights , mainBranchGalaxyWeightsSquared
     ! Array for the covariance matrix.
     double precision                        , allocatable, dimension(:,:) :: massFunctionCovariance
     ! Cosmology conversion factors.
     double precision                                                      :: cosmologyConversionMass , cosmologyConversionMassFunction
  end type massFunction

  ! Mass functions.
  type(massFunction), allocatable, dimension(:) :: massFunctions

  ! Type for storing temporary mass functions during cumulation.
  type :: massFunctionWork
     double precision, allocatable, dimension(:  ) :: massFunction
     double precision, allocatable, dimension(:,:) :: covariance
  end type massFunctionWork

  ! Work array.
  type(massFunctionWork), allocatable, dimension(:) :: thisGalaxy
  !$omp threadprivate(thisGalaxy)

  ! Options controlling binning in halo mass.
  integer                     :: analysisMassFunctionCovarianceModel
  integer         , parameter :: analysisMassFunctionCovarianceModelPoisson =1
  integer         , parameter :: analysisMassFunctionCovarianceModelBinomial=2
  integer                     :: analysisMassFunctionsHaloMassBinsCount                 , analysisMassFunctionsHaloMassBinsPerDecade
  double precision            :: analysisMassFunctionsHaloMassMinimum                   , analysisMassFunctionsHaloMassMaximum           , &
       &                         analysisMassFunctionsHaloMassIntervalLogarithmicInverse, analysisMassFunctionsHaloMassMinimumLogarithmic

  ! Options controlling covariance matrix construction.
  double precision            :: analysisMassFunctionsCorrelationTruncateLevel

  ! Initializations for individual mass functions.
  logical                     :: alfalfaHiMassFunctionZ0_00Initialized=.false.

  ! Parameters for individual mass functions.
  double precision            :: alfalfaHiMassFunctionZ0_00MolecularFractionK           , alfalfaHiMassFunctionZ0_00ErrorA               , &
       &                         alfalfaHiMassFunctionZ0_00ErrorB                       , alfalfaHiMassFunctionZ0_00ErrorC               , &
       &                         alfalfaHiMassFunctionZ0_00MolecularFractionfSigma      , alfalfaHiMassFunctionZ0_00MolecularFractionA1  , &
       &                         alfalfaHiMassFunctionZ0_00MolecularFractionAlpha1      , alfalfaHiMassFunctionZ0_00MolecularFractionA2  , &
       &                         alfalfaHiMassFunctionZ0_00MolecularFractionAlpha2      , alfalfaHiMassFunctionZ0_00MolecularFractionBeta, &
       &                         alfalfaHiMassFunctionZ0_00MolecularFractionScatter

contains

  !# <mergerTreeAnalysisTask>
  !#  <unitName>Galacticus_Output_Analysis_Mass_Functions</unitName>
  !# </mergerTreeAnalysisTask>
  subroutine Galacticus_Output_Analysis_Mass_Functions(thisTree,thisNode,iOutput,mergerTreeAnalyses)
    !% Construct a mass functions to compare to various observational determinations.
    use Galacticus_Nodes
    use Galacticus_Input_Paths
    use IO_HDF5
    use ISO_Varying_String
    use Memory_Management
    use Galactic_Structure_Enclosed_Masses
    use Input_Parameters
    use Galacticus_Output_Times
    use Galacticus_Error
    use Cosmology_Parameters
    use Cosmology_Functions
    use Numerical_Comparison
    use String_Handling
    use Galacticus_Output_Analyses_Cosmology_Scalings
    implicit none
    type            (mergerTree                    ), intent(in   )                 :: thisTree
    type            (treeNode                      ), intent(inout), pointer        :: thisNode
    integer                                         , intent(in   )                 :: iOutput
    type            (varying_string                ), intent(in   ), dimension(:  ) :: mergerTreeAnalyses
    class           (nodeComponentBasic            )               , pointer        :: thisBasic
    type            (hdf5Object                    )                                :: dataFile,massDataset,parameters
    integer                                                                         :: i,j,k,currentAnalysis,activeAnalysisCount,haloMassBin
    double precision                                                                :: dataHubbleParameter &
         &,mass,massLogarithmic,randomError,dataOmegaMatter,dataOmegaDarkEnergy
    type            (varying_string                )                                :: parameterName,analysisMassFunctionCovarianceModelText,cosmologyScalingMass,cosmologyScalingMassFunction
    class           (cosmologyFunctionsClass       )               , pointer        :: cosmologyFunctionsModel
    type            (cosmologyFunctionsMatterLambda)                                :: cosmologyFunctionsObserved
    type            (cosmologyParametersSimple     )                                :: cosmologyParametersObserved

    ! Initialize the module if necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Galacticus_Output_Analysis_Mass_Functions_Initialize)
       if (.not.moduleInitialized) then
          !@ <inputParameter>
          !@   <name>analysisMassFunctionCovarianceModel</name>
          !@   <defaultValue>binomial</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The model to use when constructing the mass function covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisMassFunctionCovarianceModel',analysisMassFunctionCovarianceModelText,defaultValue='Poisson')
          select case (char(analysisMassFunctionCovarianceModelText))
          case ( 'Poisson'  )
             analysisMassFunctionCovarianceModel=analysisMassFunctionCovarianceModelPoisson
          case ( 'binomial' )
             analysisMassFunctionCovarianceModel=analysisMassFunctionCovarianceModelBinomial
          case default
             call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Functions','unrecognized value for "analysisMassFunctionCovarianceModel" - allowed values are "Poisson", and "binomial"')
          end select
          !@ <inputParameter>
          !@   <name>analysisMassFunctionsHaloMassBinsPerDecade</name>
          !@   <defaultValue>10</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The number of bins per decade of halo mass to use when constructing the mass function covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisMassFunctionsHaloMassBinsPerDecade',analysisMassFunctionsHaloMassBinsPerDecade,defaultValue=10)
          !@ <inputParameter>
          !@   <name>analysisMassFunctionsHaloMassMinimum</name>
          !@   <defaultValue>$10^8M_\odot$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The minimum halo mass to consider when constructing the mass function covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisMassFunctionsHaloMassMinimum',analysisMassFunctionsHaloMassMinimum,defaultValue=1.0d8)
          !@ <inputParameter>
          !@   <name>analysisMassFunctionsHaloMassMaximum</name>
          !@   <defaultValue>$10^{16}M_\odot$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The maximum halo mass to consider when constructing the mass function covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisMassFunctionsHaloMassMaximum',analysisMassFunctionsHaloMassMaximum,defaultValue=1.0d16)
          analysisMassFunctionsHaloMassMinimumLogarithmic=log10(analysisMassFunctionsHaloMassMinimum)
          analysisMassFunctionsHaloMassBinsCount=int(log10(analysisMassFunctionsHaloMassMaximum/analysisMassFunctionsHaloMassMinimum)*dble(analysisMassFunctionsHaloMassBinsPerDecade)+0.5d0)
          analysisMassFunctionsHaloMassIntervalLogarithmicInverse=dble(analysisMassFunctionsHaloMassBinsCount)/log10(analysisMassFunctionsHaloMassMaximum/analysisMassFunctionsHaloMassMinimum)
          !@ <inputParameter>
          !@   <name>analysisMassFunctionsCorrelationTruncateLevel</name>
          !@   <defaultValue>0.0</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The correlation below which off-diagonal elements of the covariance matrix are truncated to zero.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisMassFunctionsCorrelationTruncateLevel',analysisMassFunctionsCorrelationTruncateLevel,defaultValue=0.0d0)
          ! Establish mass mapping functions for mass function descriptors.
          massFunctionDescriptors(2)%mapMass             => Map_Mass_ALFALFA_HI_Mass_Function_Z0_00
          ! Establish mass error functions for mass function descriptors.
          massFunctionDescriptors(2)%randomErrorFunction => Mass_Error_ALFALFA_HI_Mass_Function_Z0_00
          ! Determine how many supported mass functions are requested.
          activeAnalysisCount=0
          do i=1,massFunctionsSupportedCount
             if (any(trim(mergerTreeAnalyses) == trim(massFunctionLabels(i)))) activeAnalysisCount=activeAnalysisCount+1
          end do
          ! Allocate mass function arrays and populate with required data.
          if (activeAnalysisCount <= 0) then
             analysisActive=.false.
          else
             analysisActive=.true.
             currentAnalysis=0
             allocate(massFunctions(activeAnalysisCount))
             cosmologyFunctionsModel => cosmologyFunctions()
             do i=1,size(mergerTreeAnalyses)
                do j=1,massFunctionsSupportedCount
                   if (mergerTreeAnalyses(i) == trim(massFunctionLabels(j))) then
                      currentAnalysis=currentAnalysis+1
                      ! Copy the descriptor for this mass function.
                      massFunctions(currentAnalysis)%descriptor => massFunctionDescriptors(j)
                      ! Read parameters of the systematic error model.
                      if (massFunctionDescriptors(j)%systematicCoefficientCount > 0) then
                         allocate(massFunctions(currentAnalysis)%systematicCoefficients(massFunctionDescriptors(j)%systematicCoefficientCount))
                         do k=1,massFunctionDescriptors(j)%systematicCoefficientCount
                            parameterName=trim(massFunctionLabels(j))//'MassSystematic'
                            parameterName=parameterName//(k-1)
                            !@ <inputParameter>
                            !@   <regEx>(sdssStellarMassFunction|alfalfaHiMassFunction|primusStellarMassFunction)Z[0-9\.]+MassSystematic[0-9]+</regEx>
                            !@   <defaultValue>0</defaultValue>
                            !@   <attachedTo>module</attachedTo>
                            !@   <description>
                            !@     Mass function systematic parameters.
                            !@   </description>
                            !@   <type>real</type>
                            !@   <cardinality>1</cardinality>
                            !@ </inputParameter>
                            call Get_Input_Parameter(char(parameterName),massFunctions(currentAnalysis)%systematicCoefficients(k),defaultValue=0.0d0)
                         end do
                      end if
                      ! Find which output number corresponds to the required redshift.
                      massFunctions(currentAnalysis)%outputNumber=-1
                      do k=1,Galacticus_Output_Time_Count()
                         if     (                                                                                             &
                              &  Values_Agree(                                                                                &
                              &               massFunctionDescriptors(j)%redshift                                           , &
                              &               cosmologyFunctionsModel%redshiftFromExpansionFactor(                            &
                              &               cosmologyFunctionsModel%expansionFactor             (                           &
                              &                                                                    Galacticus_Output_Time(k)  &
                              &                                                                   )                           &
                              &                                                                  )                          , &
                              &               absTol=0.001d0                                                                  &
                              &              )                                                                                &
                              & ) then
                            massFunctions(currentAnalysis)%outputNumber=k
                            exit
                         end if
                      end do
                      if (massFunctions(currentAnalysis)%outputNumber < 0) call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Functions','unable to find required redshift in outputs for mass function '//trim(massFunctionLabels(j)))
                      ! Read the appropriate observational data definition.
                      select case (trim(massFunctionLabels(j)))
                      case ('sdssStellarMassFunctionZ0.07')
                         ! SDSS z=0.07 stellar mass function.
                         !$omp critical(HDF5_Access)
                         call dataFile%openFile(char(Galacticus_Input_Path())//'/data/observations/massFunctionsStellar/Stellar_Mass_Function_Li_White_2009.hdf5',readOnly=.true.)
                         call dataFile   %readDataset  ('mass'          ,massFunctions(currentAnalysis)%masses)
                         massDataset=dataFile%openDataset('mass'        )
                         call massDataset%readAttribute('cosmologyScaling',cosmologyScalingMass               ,allowPseudoScalar=.true.)
                         call massDataset%close()
                         massDataset=dataFile%openDataset('massFunction')
                         call massDataset%readAttribute('cosmologyScaling',cosmologyScalingMassFunction       ,allowPseudoScalar=.true.)
                         call massDataset%close()
                         parameters =dataFile%openGroup  ('Parameters'  )
                         call parameters %readAttribute('H_0'             ,dataHubbleParameter                )
                         call parameters %readAttribute('Omega_Matter'    ,dataOmegaMatter                    )
                         call parameters %readAttribute('Omega_DE'        ,dataOmegaDarkEnergy                )
                         call parameters %close()
                         call dataFile   %close()
                         !$omp end critical(HDF5_Access)
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
                         ! Construct mass function array.
                         massFunctions(currentAnalysis)%massesCount=size(massFunctions(currentAnalysis)%masses)
                         call Alloc_Array(massFunctions(currentAnalysis)%massesLogarithmic             ,[massFunctions(currentAnalysis)%massesCount                                           ])
                         call Alloc_Array(massFunctions(currentAnalysis)%massesLogarithmicMinimum      ,[massFunctions(currentAnalysis)%massesCount                                           ])
                         call Alloc_Array(massFunctions(currentAnalysis)%massesLogarithmicMaximum      ,[massFunctions(currentAnalysis)%massesCount                                           ])
                         call Alloc_Array(massFunctions(currentAnalysis)%massFunction                  ,[massFunctions(currentAnalysis)%massesCount                                           ])
                         call Alloc_Array(massFunctions(currentAnalysis)%massFunctionCovariance        ,[massFunctions(currentAnalysis)%massesCount,massFunctions(currentAnalysis)%massesCount])
                         call Alloc_Array(massFunctions(currentAnalysis)%mainBranchGalaxyWeights       ,[massFunctions(currentAnalysis)%massesCount,analysisMassFunctionsHaloMassBinsCount    ])
                         call Alloc_Array(massFunctions(currentAnalysis)%mainBranchGalaxyWeightsSquared,[massFunctions(currentAnalysis)%massesCount,analysisMassFunctionsHaloMassBinsCount    ])
                         massFunctions(currentAnalysis)%massesLogarithmic             =log10(massFunctions(currentAnalysis)%masses)
                         massFunctions(currentAnalysis)%massFunction                  =0.0d0
                         massFunctions(currentAnalysis)%massFunctionCovariance        =0.0d0
                         massFunctions(currentAnalysis)%mainBranchGalaxyWeights       =0.0d0
                         massFunctions(currentAnalysis)%mainBranchGalaxyWeightsSquared=0.0d0
                         do k=1,massFunctions(currentAnalysis)%massesCount
                            if (k ==                                          1) then
                               massFunctions(currentAnalysis)%massesLogarithmicMinimum(k)=massFunctions(currentAnalysis)%massesLogarithmic(k)-0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k+1)-massFunctions(currentAnalysis)%massesLogarithmic(k  ))
                            else
                               massFunctions(currentAnalysis)%massesLogarithmicMinimum(k)=                                                   +0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k-1)+massFunctions(currentAnalysis)%massesLogarithmic(k  ))
                            end if
                            if (k == massFunctions(currentAnalysis)%massesCount) then
                               massFunctions(currentAnalysis)%massesLogarithmicMaximum(k)=massFunctions(currentAnalysis)%massesLogarithmic(k)+0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k  )-massFunctions(currentAnalysis)%massesLogarithmic(k-1))
                            else
                               massFunctions(currentAnalysis)%massesLogarithmicMaximum(k)=                                                   +0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k+1)+massFunctions(currentAnalysis)%massesLogarithmic(k  ))
                            end if
                         end do
                      case ('alfalfaHiMassFunctionZ0.00')
                         !$omp critical(HDF5_Access)
                         call dataFile%openFile(char(Galacticus_Input_Path())//'/data/observations/massFunctionsHI/HI_Mass_Function_ALFALFA_2010.hdf5',readOnly=.true.)
                         call dataFile   %readDataset  ('mass'          ,massFunctions(currentAnalysis)%masses             )
                         massDataset=dataFile%openDataset('mass'      )
                         call massDataset%readAttribute('cosmologyScaling',cosmologyScalingMass               ,allowPseudoScalar=.true.)
                         call massDataset%close()
                         massDataset=dataFile%openDataset('massFunction')
                         call massDataset%readAttribute('cosmologyScaling',cosmologyScalingMassFunction       ,allowPseudoScalar=.true.)
                         call massDataset%close()
                         parameters =dataFile%openGroup  ('Parameters')
                         call parameters %readAttribute('H_0'           ,dataHubbleParameter                               )
                         call parameters %readAttribute('Omega_Matter'  ,dataOmegaMatter                                   )
                         call parameters %readAttribute('Omega_DE'      ,dataOmegaDarkEnergy                               )
                         call parameters %close()
                         call dataFile   %close()
                         !$omp end critical(HDF5_Access)
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
                         ! Construct mass function array.
                         massFunctions(currentAnalysis)%massesCount=size(massFunctions(currentAnalysis)%masses)
                         call Alloc_Array(massFunctions(currentAnalysis)%massesLogarithmic             ,[massFunctions(currentAnalysis)%massesCount                                           ])
                         call Alloc_Array(massFunctions(currentAnalysis)%massesLogarithmicMinimum      ,[massFunctions(currentAnalysis)%massesCount                                           ])
                         call Alloc_Array(massFunctions(currentAnalysis)%massesLogarithmicMaximum      ,[massFunctions(currentAnalysis)%massesCount                                           ])
                         call Alloc_Array(massFunctions(currentAnalysis)%massFunction                  ,[massFunctions(currentAnalysis)%massesCount                                           ])
                         call Alloc_Array(massFunctions(currentAnalysis)%massFunctionCovariance        ,[massFunctions(currentAnalysis)%massesCount,massFunctions(currentAnalysis)%massesCount])
                         call Alloc_Array(massFunctions(currentAnalysis)%mainBranchGalaxyWeights       ,[massFunctions(currentAnalysis)%massesCount,analysisMassFunctionsHaloMassBinsCount    ])
                         call Alloc_Array(massFunctions(currentAnalysis)%mainBranchGalaxyWeightsSquared,[massFunctions(currentAnalysis)%massesCount,analysisMassFunctionsHaloMassBinsCount    ])
                         massFunctions(currentAnalysis)%massesLogarithmic             =log10(massFunctions(currentAnalysis)%masses)
                         massFunctions(currentAnalysis)%massFunction                  =0.0d0
                         massFunctions(currentAnalysis)%massFunctionCovariance        =0.0d0
                         massFunctions(currentAnalysis)%mainBranchGalaxyWeights       =0.0d0
                         massFunctions(currentAnalysis)%mainBranchGalaxyWeightsSquared=0.0d0
                         do k=1,massFunctions(currentAnalysis)%massesCount
                            if (k ==                                          1) then
                               massFunctions(currentAnalysis)%massesLogarithmicMinimum(k)=massFunctions(currentAnalysis)%massesLogarithmic(k)-0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k+1)-massFunctions(currentAnalysis)%massesLogarithmic(k  ))
                            else
                               massFunctions(currentAnalysis)%massesLogarithmicMinimum(k)=                                                   +0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k-1)+massFunctions(currentAnalysis)%massesLogarithmic(k  ))
                            end if
                            if (k == massFunctions(currentAnalysis)%massesCount) then
                               massFunctions(currentAnalysis)%massesLogarithmicMaximum(k)=massFunctions(currentAnalysis)%massesLogarithmic(k)+0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k  )-massFunctions(currentAnalysis)%massesLogarithmic(k-1))
                            else
                               massFunctions(currentAnalysis)%massesLogarithmicMaximum(k)=                                                   +0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k+1)+massFunctions(currentAnalysis)%massesLogarithmic(k  ))
                            end if
                         end do
                      case ('primusStellarMassFunctionZ0.100')
                         call Load_PRIMUS_Mass_Function(0,massFunctions(currentAnalysis),cosmologyParametersObserved,cosmologyFunctionsObserved,cosmologyScalingMass,cosmologyScalingMassFunction)
                      case ('primusStellarMassFunctionZ0.250')
                         call Load_PRIMUS_Mass_Function(1,massFunctions(currentAnalysis),cosmologyParametersObserved,cosmologyFunctionsObserved,cosmologyScalingMass,cosmologyScalingMassFunction)
                      case ('primusStellarMassFunctionZ0.350')
                         call Load_PRIMUS_Mass_Function(2,massFunctions(currentAnalysis),cosmologyParametersObserved,cosmologyFunctionsObserved,cosmologyScalingMass,cosmologyScalingMassFunction)
                      case ('primusStellarMassFunctionZ0.450')
                         call Load_PRIMUS_Mass_Function(3,massFunctions(currentAnalysis),cosmologyParametersObserved,cosmologyFunctionsObserved,cosmologyScalingMass,cosmologyScalingMassFunction)
                      case ('primusStellarMassFunctionZ0.575')
                         call Load_PRIMUS_Mass_Function(4,massFunctions(currentAnalysis),cosmologyParametersObserved,cosmologyFunctionsObserved,cosmologyScalingMass,cosmologyScalingMassFunction)
                      case ('primusStellarMassFunctionZ0.725')
                         call Load_PRIMUS_Mass_Function(5,massFunctions(currentAnalysis),cosmologyParametersObserved,cosmologyFunctionsObserved,cosmologyScalingMass,cosmologyScalingMassFunction)
                      case ('primusStellarMassFunctionZ0.900')
                         call Load_PRIMUS_Mass_Function(6,massFunctions(currentAnalysis),cosmologyParametersObserved,cosmologyFunctionsObserved,cosmologyScalingMass,cosmologyScalingMassFunction)
                      case default
                         call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Functions','unknown mass function')
                      end select
                      ! Get cosmological conversion factors.
                      call Cosmology_Conversion_Factors(                                                                                                &
                           &                            massFunctions(currentAnalysis)%descriptor%redshift                                            , &
                           &                            cosmologyFunctionsModel                                                                       , &
                           &                            cosmologyFunctionsObserved                                                                    , &
                           &                            cosmologyScalingMass           =cosmologyScalingMass                                          , &
                           &                            cosmologyScalingMassFunction   =cosmologyScalingMassFunction                                  , &
                           &                            cosmologyConversionMass        =massFunctions(currentAnalysis)%cosmologyConversionMass        , &
                           &                            cosmologyConversionMassFunction=massFunctions(currentAnalysis)%cosmologyConversionMassFunction  &
                           &                           )
                      exit
                   end if
                end do
             end do
          end if
          ! Record that module is initialized.
          moduleInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Analysis_Mass_Functions_Initialize)
    end if
    ! Return if this analysis is not active.
    if (.not.analysisActive) return
    ! Allocate work arrays.
    if (.not.allocated(thisGalaxy)) allocate(thisGalaxy(size(massFunctions)))
    ! Iterate over active analyses.
    do i=1,size(massFunctions)
       ! Return if this analysis is not active, or if this is not the correct output.
       if (iOutput /= massFunctions(i)%outputNumber) cycle
       ! Allocate workspace.
       if (.not.allocated(thisGalaxy(i)%massFunction)) call Alloc_Array(thisGalaxy(i)%massFunction         ,[massFunctions(i)%massesCount                             ])
       if (.not.allocated(thisGalaxy(i)%covariance  )) call Alloc_Array(thisGalaxy(i)%covariance           ,[massFunctions(i)%massesCount,massFunctions(i)%massesCount])
       ! Get the galactic mass.
       mass=                                                                                                                &
            &  Galactic_Structure_Enclosed_Mass(thisNode,radiusLarge,componentType=componentTypeDisk    ,massType=massFunctions(i)%descriptor%massType) &
            & +Galactic_Structure_Enclosed_Mass(thisNode,radiusLarge,componentType=componentTypeSpheroid,massType=massFunctions(i)%descriptor%massType)
       if (mass            <=                  0.0d0) return
       if (associated(massFunctions(i)%descriptor%mapMass)) mass=massFunctions(i)%descriptor%mapMass(mass,thisNode)
       mass=mass*massFunctions(i)%cosmologyConversionMass ! Convert for cosmology.
       massLogarithmic=log10(mass)
       do j=1,massFunctions(i)%descriptor%systematicCoefficientCount
          massLogarithmic=massLogarithmic+massFunctions(i)%systematicCoefficients(j)*(log10(mass)-massFunctions(i)%descriptor%systematicLogM0)**(j-1)
       end do
       if (massLogarithmic <  massFunctions(i)%descriptor%massLogarithmicMinimum) return
       ! Compute contributions to each bin.
       randomError=massFunctions(i)%descriptor%randomError
       if (associated(massFunctions(i)%descriptor%randomErrorFunction)) randomError=massFunctions(i)%descriptor%randomErrorFunction(mass,thisNode)
       thisGalaxy(i)%massFunction=(                                                                                          &
            &                      +erf((massFunctions(i)%massesLogarithmicMaximum-massLogarithmic)/randomError/sqrt(2.0d0)) &
            &                      -erf((massFunctions(i)%massesLogarithmicMinimum-massLogarithmic)/randomError/sqrt(2.0d0)) &
            &                     )                                                                                          &
            &                     /2.0d0                                                                                     &
            &                     *thisTree%volumeWeight                                                                     &
            &                     *massFunctions(i)%cosmologyConversionMassFunction
       ! Accumulate mass function.
       !$omp critical (Galacticus_Output_Analysis_Mass_Functions_Accumulate)
       massFunctions(i)%massFunction          =massFunctions(i)%massFunction          +thisGalaxy(i)%massFunction
       !$omp end critical (Galacticus_Output_Analysis_Mass_Functions_Accumulate)
       ! Treat main branch and other galaxies differently.
       if (thisNode%isOnMainBranch().and.analysisMassFunctionCovarianceModel == analysisMassFunctionCovarianceModelBinomial) then
          ! Find the bin to which this halo mass belongs.
          thisBasic => thisNode%basic()
          haloMassBin=floor((log10(thisBasic%mass())-analysisMassFunctionsHaloMassMinimumLogarithmic)*analysisMassFunctionsHaloMassIntervalLogarithmicInverse)+1
          ! Accumulate weights to halo mass arrays.
          if (haloMassBin >= 1 .and. haloMassBin <= analysisMassFunctionsHaloMassBinsCount) then
            !$omp critical (Galacticus_Output_Analysis_Mass_Functions_Accumulate)
             massFunctions        (i)%mainBranchGalaxyWeights       (:,haloMassBin)= &
                  &  massFunctions(i)%mainBranchGalaxyWeights       (:,haloMassBin)  &
                  &  +thisGalaxy  (i)%massFunction
             massFunctions        (i)%mainBranchGalaxyWeightsSquared(:,haloMassBin)= &
                  &  massFunctions(i)%mainBranchGalaxyWeightsSquared(:,haloMassBin)  &
                  &  +thisGalaxy  (i)%massFunction**2
             !$omp end critical (Galacticus_Output_Analysis_Mass_Functions_Accumulate)
          end if
       else
          forall(j=1:massFunctions(i)%massesCount)
             forall(k=j:massFunctions(i)%massesCount)
                thisGalaxy(i)%covariance(j,k)=thisGalaxy(i)%massFunction(j)*thisGalaxy(i)%massFunction(k)
                thisGalaxy(i)%covariance(k,j)=thisGalaxy(i)%covariance(j,k)
             end forall
          end forall
          ! Accumulate covariance.
          !$omp critical (Galacticus_Output_Analysis_Mass_Functions_Accumulate)
          massFunctions(i)%massFunctionCovariance=massFunctions(i)%massFunctionCovariance+thisGalaxy(i)%covariance
          !$omp end critical (Galacticus_Output_Analysis_Mass_Functions_Accumulate)
       end if
    end do
    return
  end subroutine Galacticus_Output_Analysis_Mass_Functions

  !# <hdfPreCloseTask>
  !#  <unitName>Galacticus_Output_Analysis_Mass_Functions_Output</unitName>
  !# </hdfPreCloseTask>
  subroutine Galacticus_Output_Analysis_Mass_Functions_Output
    !% Outputs SDSS $z\approx 0.07$ stellar mass function to file.
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
       if (analysisMassFunctionCovarianceModel == analysisMassFunctionCovarianceModelBinomial) then
          do m=1,analysisMassFunctionsHaloMassBinsCount
             haloWeightBinTotal=sum(massFunctions(k)%mainBranchGalaxyWeights(:,m))
             if ( haloWeightBinTotal > 0.0 ) then
                do i=1,massFunctions(k)%massesCount
                   massFunctions               (k)%massFunctionCovariance        (i,i)=                    &
                        &         massFunctions(k)%massFunctionCovariance        (i,i)                     &
                        & +(1.0d0-massFunctions(k)%mainBranchGalaxyWeights       (i,m)/haloWeightBinTotal) &
                        & *       massFunctions(k)%mainBranchGalaxyWeightsSquared(i,m)
                   do j=1,massFunctions(k)%massesCount
                      if (i == j) cycle
                      massFunctions               (k)%massFunctionCovariance        (i,j)= &
                           &  +      massFunctions(k)%massFunctionCovariance        (i,j)  &
                           &  -      massFunctions(k)%mainBranchGalaxyWeights       (i,m)  &
                           &  *      massFunctions(k)%mainBranchGalaxyWeights       (j,m)  &
                           &  *sqrt(                                                       &
                           &         massFunctions(k)%mainBranchGalaxyWeightsSquared(i,m)  &
                           &        *massFunctions(k)%mainBranchGalaxyWeightsSquared(j,m)  &
                           &       )                                                       &
                           &  /haloWeightBinTotal
                   end do
                end do
             end if
          end do
       end if
       ! Truncate the covariance where the correlation is below threshold. This avoids creating noisy covariances matrices
       ! which can lead to discontinuities in the likelihood surface.
       do i=1,massFunctions(k)%massesCount
          do j=1,massFunctions(k)%massesCount
             if (i == j) cycle
             if     (                                                           &
                  &     abs( massFunctions(k)%massFunctionCovariance(i,j))      &
                  &  <                                                          &
                  &    analysisMassFunctionsCorrelationTruncateLevel            &
                  &   *sqrt(                                                    &
                  &          massFunctions(k)%massFunctionCovariance(i,i)       &
                  &         *massFunctions(k)%massFunctionCovariance(j,j)       &
                  &        )                                                    &
                  & )        massFunctions(k)%massFunctionCovariance(i,j)=0.0d0
          end do
       end do
       ! Convert model mass function to differential per log(M).
       forall(i=1:massFunctions(k)%massesCount)
          massFunctions(k)%massFunction             (i  )=  massFunctions(k)%massFunction          (i  )                              &
               &                         /(massFunctions(k)%massesLogarithmicMaximum(i)-massFunctions(k)%massesLogarithmicMinimum(i)) &
               &                         / log(10.0d0)
          forall(j=1:massFunctions(k)%massesCount)
             massFunctions(k)%massFunctionCovariance(i,j)=  massFunctions(k)%massFunctionCovariance(i,j)                              &
                  &                      /(massFunctions(k)%massesLogarithmicMaximum(i)-massFunctions(k)%massesLogarithmicMinimum(i)) &
                  &                      /(massFunctions(k)%massesLogarithmicMaximum(j)-massFunctions(k)%massesLogarithmicMinimum(j)) &
                  &                      / log(10.0d0)**2
          end forall
       end forall
       ! Output the mass function.
       !$omp critical(HDF5_Access)
       analysisGroup    =galacticusOutputFile%openGroup('analysis','Model analysis')
       massFunctionGroup=analysisGroup       %openGroup(trim(massFunctions(k)%descriptor%label),trim(massFunctions(k)%descriptor%comment))
       call massFunctionGroup%writeDataset  (massFunctions(k)%masses                ,'mass'                  ,'Mass'             ,datasetReturned=thisDataset)
       call thisDataset      %writeAttribute(massSolar             ,'unitsInSI'                                                                            )
       call thisDataset      %close()
       call massFunctionGroup%writeDataset  (massFunctions(k)%massFunction          ,'massFunction'          ,'Mass function'           ,datasetReturned=thisDataset)
       call thisDataset      %writeAttribute(1.0d0/megaParsec**3   ,'unitsInSI'                                                                            )
       call thisDataset      %close()
       call massFunctionGroup%writeDataset  (massFunctions(k)%massFunctionCovariance,'massFunctionCovariance','Mass function covariance',datasetReturned=thisDataset)
       call thisDataset      %writeAttribute(1.0d0/megaParsec**6   ,'unitsInSI'                                                                            )
       call thisDataset      %close()
       call massFunctionGroup%close()
       call analysisGroup    %close()
       !$omp end critical(HDF5_Access)
    end do
    return
  end subroutine Galacticus_Output_Analysis_Mass_Functions_Output

  subroutine ALFALFA_HI_Mass_Function_Z0_00_Initialize()
    !% Initializes the ALFALFA HI mass function calculation by reading in required parameters.
    use Input_Parameters
    implicit none

    ! Initialize the mass function if necessary.
    if (.not.alfalfaHiMassFunctionZ0_00Initialized) then
       !$omp critical(alfalfaHiMassFunctionZ0_00Initialized)
       if (.not.alfalfaHiMassFunctionZ0_00Initialized) then
          ! Read parameters controlling H2/HI mass ratio calculation.
          !@ <inputParameter>
          !@   <name>alfalfaHiMassFunctionZ0.00ErrorA</name>
          !@   <defaultValue>0.1</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The parameter, $a$, appearing in the model for the HI mass random errors used when constructing the ALFALFA HI mass function. Specifically, $\sigma_{\rm obs} = a + \exp\left(-{\log_{10}(M_{\rm HI}/M_\odot)-b\over c}\right)$.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('alfalfaHiMassFunctionZ0.00ErrorA',alfalfaHiMassFunctionZ0_00ErrorA,defaultValue=0.1d0)
          !@ <inputParameter>
          !@   <name>alfalfaHiMassFunctionZ0.00ErrorB</name>
          !@   <defaultValue>5.885</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The parameter, $b$, appearing in the model for the HI mass random errors used when constructing the ALFALFA HI mass function. Specifically, $\sigma_{\rm obs} = a + \exp\left(-{\log_{10}(M_{\rm HI}/M_\odot)-b\over c}\right)$.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('alfalfaHiMassFunctionZ0.00ErrorB',alfalfaHiMassFunctionZ0_00ErrorB,defaultValue=5.885d0)
          !@ <inputParameter>
          !@   <name>alfalfaHiMassFunctionZ0.00ErrorC</name>
          !@   <defaultValue>0.505</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The parameter, $c$, appearing in the model for the HI mass random errors used when constructing the ALFALFA HI mass function. Specifically, $\sigma_{\rm obs} = a + \exp\left(-{\log_{10}(M_{\rm HI}/M_\odot)-b\over c}\right)$.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('alfalfaHiMassFunctionZ0.00ErrorC',alfalfaHiMassFunctionZ0_00ErrorC,defaultValue=0.505d0)
          !@ <inputParameter>
          !@   <name>alfalfaHiMassFunctionZ0.00MolecularFractionK</name>
          !@   <defaultValue>11.3 m$^4$ kg$^{-2}$ \citep{obreschkow_simulation_2009}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The parameter, $K$, appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('alfalfaHiMassFunctionZ0.00MolecularFractionK',alfalfaHiMassFunctionZ0_00MolecularFractionK,defaultValue=11.3d0)
          !@ <inputParameter>
          !@   <name>alfalfaHiMassFunctionZ0.00MolecularFractionfSigma</name>
          !@   <defaultValue>0.4 \citep{obreschkow_simulation_2009}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The parameter, $\langle f_\sigma \rangle$, appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('alfalfaHiMassFunctionZ0.00MolecularFractionfSigma',alfalfaHiMassFunctionZ0_00MolecularFractionfSigma,defaultValue=0.4d0)
          !@ <inputParameter>
          !@   <name>alfalfaHiMassFunctionZ0.00MolecularFractionA1</name>
          !@   <defaultValue>3.44 \citep{obreschkow_simulation_2009}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The parameter, $A_1$, appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('alfalfaHiMassFunctionZ0.00MolecularFractionA1',alfalfaHiMassFunctionZ0_00MolecularFractionA1,defaultValue=3.44d0)
          !@ <inputParameter>
          !@   <name>alfalfaHiMassFunctionZ0.00MolecularFractionA2</name>
          !@   <defaultValue>4.82 \citep{obreschkow_simulation_2009}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The parameter, $A_2$, appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('alfalfaHiMassFunctionZ0.00MolecularFractionA2',alfalfaHiMassFunctionZ0_00MolecularFractionA2,defaultValue=4.82d0)
          !@ <inputParameter>
          !@   <name>alfalfaHiMassFunctionZ0.00MolecularFractionAlpha1</name>
          !@   <defaultValue>$-0.506$ \citep{obreschkow_simulation_2009}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The parameter, $\alpha_1$, appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('alfalfaHiMassFunctionZ0.00MolecularFractionAlpha1',alfalfaHiMassFunctionZ0_00MolecularFractionAlpha1,defaultValue=-0.506d0)
          !@ <inputParameter>
          !@   <name>alfalfaHiMassFunctionZ0.00MolecularFractionAlpha2</name>
          !@   <defaultValue>$-1.054$ \citep{obreschkow_simulation_2009}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The parameter, $\alpha_2$, appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('alfalfaHiMassFunctionZ0.00MolecularFractionAlpha2',alfalfaHiMassFunctionZ0_00MolecularFractionAlpha2,defaultValue=-1.054d0)
          !@ <inputParameter>
          !@   <name>alfalfaHiMassFunctionZ0.00MolecularFractionBeta</name>
          !@   <defaultValue>$0.8$ \citep{obreschkow_simulation_2009}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The parameter, $\beta$, appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('alfalfaHiMassFunctionZ0.00MolecularFractionBeta',alfalfaHiMassFunctionZ0_00MolecularFractionBeta,defaultValue=0.8d0)
          !@ <inputParameter>
          !@   <name>alfalfaHiMassFunctionZ0.00MolecularFractionScatter</name>
          !@   <defaultValue>$0.4$~dex (Obsreschkow, private communication)</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The scatter in the molecular ratio $\log_{10}R_{\rm mol}$ of \cite{obreschkow_simulation_2009} compared to observational data.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('alfalfaHiMassFunctionZ0.00MolecularFractionScatter',alfalfaHiMassFunctionZ0_00MolecularFractionScatter,defaultValue=0.4d0)
          ! Record that the function is now initialized.
          alfalfaHiMassFunctionZ0_00Initialized=.true.
       end if
       !$omp end critical(alfalfaHiMassFunctionZ0_00Initialized)
    end if
    return
  end subroutine ALFALFA_HI_Mass_Function_Z0_00_Initialize

  double precision function Mass_Error_ALFALFA_HI_Mass_Function_Z0_00(mass,thisNode)
    !% Computes errors on $\log_{10}($HI masses$)$ for the ALFALFA survey analysis. Uses a simple fitting function. See {\tt
    !% constraints/dataAnalysis/hiMassFunction\_ALFALFA\_z0.00/alfalfaHIMassErrorModel.pl} for details.
    double precision          , intent(in   )          :: mass
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision                                   :: logarithmicMass, molecularFraction, &
         &                                                molecularRatio

    ! Initialize the mass function.
    call ALFALFA_HI_Mass_Function_Z0_00_Initialize()
    ! Get the logarithmic mass.
    logarithmicMass=log10(mass)
    ! Compute the random error on the mass.
    Mass_Error_ALFALFA_HI_Mass_Function_Z0_00              &
         &        =alfalfaHiMassFunctionZ0_00ErrorA        &
         &        +exp(                                    &
         &             -(                                  &
         &               +max(logarithmicMass,6.0d0)       &
         &               -alfalfaHiMassFunctionZ0_00ErrorB &
         &              )                                  &
         &             /alfalfaHiMassFunctionZ0_00ErrorC   &
         &            )
    ! Add in quadrature a term accounting for the scatter in the molecular ratio model.
    molecularRatio=Molecular_Ratio_ALFALFA_HI_Mass_Function_Z0_00(mass,thisNode)
    Mass_Error_ALFALFA_HI_Mass_Function_Z0_00                             &
         & =sqrt(                                                         &
         &       +Mass_Error_ALFALFA_HI_Mass_Function_Z0_00           **2 &
         &       +(                                                       &
         &         +alfalfaHiMassFunctionZ0_00MolecularFractionScatter    &
         &         *molecularFraction                                     &
         &         /(                                                     &
         &           +1.0                                                 &
         &           +molecularFraction                                   &
         &          )                                                     &
         &        )                                                   **2 &
         &      )
    return
  end function Mass_Error_ALFALFA_HI_Mass_Function_Z0_00

  double precision function Map_Mass_ALFALFA_HI_Mass_Function_Z0_00(mass,thisNode)
    !% Maps gas masses into HI masses for the ALFALFA survey analysis.
    use Numerical_Constants_Astronomical
    implicit none                                                                                 
    double precision          , intent(in   )          :: mass
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision                                   :: molecularRatio

    ! Compute the HI mass.
    molecularRatio=Molecular_Ratio_ALFALFA_HI_Mass_Function_Z0_00(mass,thisNode)
    Map_Mass_ALFALFA_HI_Mass_Function_Z0_00=hydrogenByMassPrimordial*mass/(1.0d0+molecularRatio)
    return
  end function Map_Mass_ALFALFA_HI_Mass_Function_Z0_00

  double precision function Molecular_Ratio_ALFALFA_HI_Mass_Function_Z0_00(mass,thisNode)
    !% Compute the molecular ratio, $R_{\rm mol}=M_{\rm H_2}/M_{\rm HI}$ for the ALFALFA survey analysis. Assumes the model of
    !% \cite{obreschkow_simulation_2009}.
    use Numerical_Constants_Astronomical
    implicit none
    double precision                   , intent(in   )          :: mass
    type            (treeNode         ), intent(inout), pointer :: thisNode
    class           (nodeComponentDisk)               , pointer :: thisDisk
    double precision                                            :: massStellar, molecularRatioCentral, &
         &                                                         diskRadius

    ! Initialize the mass function.
    call ALFALFA_HI_Mass_Function_Z0_00_Initialize()
    ! Get disk properties.
    thisDisk    => thisNode%disk       ()
    massStellar =  thisDisk%massStellar()
    diskRadius  =  thisDisk%radius     ()
    ! Get the molecular to atomic mass ratio (H2/HI).
    molecularRatioCentral                                         &
         & =(                                                     &
         &   massSolar  **2                                       &
         &   /megaParsec**4                                       &
         &   *alfalfaHiMassFunctionZ0_00MolecularFractionK        &
         &   *mass                                                &
         &   *(                                                   &
         &     +mass                                              &
         &     +alfalfaHiMassFunctionZ0_00MolecularFractionfSigma &
         &     *massStellar                                       &
         &    )                                                   &
         &   /diskRadius**4                                       &
         &  )**alfalfaHiMassFunctionZ0_00MolecularFractionBeta
    Molecular_Ratio_ALFALFA_HI_Mass_Function_Z0_00                                      &
         & = 1.0d0                                                                      &
         &  /(                                                                          &
         &    +                       alfalfaHiMassFunctionZ0_00MolecularFractionA1     &
         &    /molecularRatioCentral**alfalfaHiMassFunctionZ0_00MolecularFractionAlpha1 &
         &    +                       alfalfaHiMassFunctionZ0_00MolecularFractionA2     &
         &    /molecularRatioCentral**alfalfaHiMassFunctionZ0_00MolecularFractionAlpha2 &
         &   )
    return
  end function Molecular_Ratio_ALFALFA_HI_Mass_Function_Z0_00

  subroutine Load_PRIMUS_Mass_Function(massFunctionIndex,thisMassFunction,cosmologyParametersObserved,cosmologyFunctionsObserved,cosmologyScalingMass,cosmologyScalingMassFunction)
    !% Load the specified mass function from the PRIMUS stellar mass function dataset.
    use ISO_Varying_String
    use Galacticus_Error
    use Cosmology_Functions
    use Cosmology_Parameters
    use Galacticus_Input_Paths
    use Memory_Management
    use FoX_dom
    use IO_XML
    implicit none
    integer                                         , intent(in   ) :: massFunctionIndex
    type            (massFunction                  ), intent(inout) :: thisMassFunction
    type            (cosmologyFunctionsMatterLambda), intent(inout) :: cosmologyFunctionsObserved
    type            (cosmologyParametersSimple     ), intent(inout) :: cosmologyParametersObserved
    type            (varying_string                ), intent(  out) :: cosmologyScalingMass,cosmologyScalingMassFunction
    type            (node                          ), pointer       :: doc,massFunctionElement,columnElement,massElement,hubbleElement,datum,cosmology,omegaDarkEnergyElement,omegaMatterElement,datasetElement,cosmologyScalingElement
    type            (nodeList                      ), pointer       :: dataList
    character       (len=128                       )                :: cosmologyScaling
    integer                                                         :: k,iDatum,ioErr
    double precision                                                :: dataHubbleParameter,dataOmegaMatter,dataOmegaDarkEnergy

    !$omp critical (FoX_DOM_Access)
    doc => parseFile(char(Galacticus_Input_Path())//"data/observations/massFunctionsStellar/Stellar_Mass_Function_PRIMUS_2013.xml",iostat=ioErr)
    if (ioErr /= 0) call Galacticus_Error_Report('Load_PRIMUS_Mass_Function','Unable to find data file')
    datasetElement          => item(getElementsByTagname(doc           ,"stellarMassFunction"),massFunctionIndex)
    columnElement           => item(getElementsByTagname(datasetElement,"columns"            ),                0)
    massElement             => item(getElementsByTagname( columnElement,"stellarMass"        ),                0)
    massFunctionElement     => item(getElementsByTagname( columnElement,"massFunction"       ),                0)
    cosmologyScalingElement => XML_Get_First_Element_By_Tag_Name(massElement        ,"cosmologyScaling")
    call extractDataContent(cosmologyScalingElement,cosmologyScaling)
    cosmologyScalingMass        =trim(cosmologyScaling)
    cosmologyScalingElement => XML_Get_First_Element_By_Tag_Name(massFunctionElement,"cosmologyScaling")
    call extractDataContent(cosmologyScalingElement,cosmologyScaling)
    cosmologyScalingMassFunction=trim(cosmologyScaling)
    cosmology               => XML_Get_First_Element_By_Tag_Name(doc      ,"cosmology"      )
    hubbleElement           => XML_Get_First_Element_By_Tag_Name(cosmology,"hubble"         )
    omegaMatterElement      => XML_Get_First_Element_By_Tag_Name(cosmology,"omegaMatter"    )
    omegaDarkEnergyElement  => XML_Get_First_Element_By_Tag_Name(cosmology,"omegaDarkEnergy")
    call extractDataContent(hubbleElement         ,dataHubbleParameter)
    call extractDataContent(omegaMatterElement    ,dataOmegaMatter    )
    call extractDataContent(omegaDarkEnergyElement,dataOmegaDarkEnergy)
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
    ! Construct mass function array.
    dataList => getElementsByTagname(massElement,"datum")
    thisMassFunction%massesCount=getLength(dataList)
    call Alloc_Array(thisMassFunction%masses                        ,[thisMassFunction%massesCount                                          ])
    call Alloc_Array(thisMassFunction%massesLogarithmic             ,[thisMassFunction%massesCount                                          ])
    call Alloc_Array(thisMassFunction%massesLogarithmicMinimum      ,[thisMassFunction%massesCount                                          ])
    call Alloc_Array(thisMassFunction%massesLogarithmicMaximum      ,[thisMassFunction%massesCount                                          ])
    call Alloc_Array(thisMassFunction%massFunction                  ,[thisMassFunction%massesCount                                          ])
    call Alloc_Array(thisMassFunction%massFunctionCovariance        ,[thisMassFunction%massesCount,thisMassFunction%massesCount             ])
    call Alloc_Array(thisMassFunction%mainBranchGalaxyWeights       ,[thisMassFunction%massesCount,analysisMassFunctionsHaloMassBinsCount   ])
    call Alloc_Array(thisMassFunction%mainBranchGalaxyWeightsSquared,[thisMassFunction%massesCount,analysisMassFunctionsHaloMassBinsCount   ])
    do iDatum=0,getLength(dataList)-1
       datum => item(dataList,iDatum)
       call extractDataContent(datum,thisMassFunction%massesLogarithmic(iDatum+1))
    end do
    ! Destroy the document.
    call destroy(doc)
    !$omp end critical (FoX_DOM_Access)
    thisMassFunction%masses                        =10.0d0**thisMassFunction%massesLogarithmic
    thisMassFunction%massFunction                  = 0.0d0
    thisMassFunction%massFunctionCovariance        = 0.0d0
    thisMassFunction%mainBranchGalaxyWeights       = 0.0d0
    thisMassFunction%mainBranchGalaxyWeightsSquared= 0.0d0
    do k=1,thisMassFunction%massesCount
       if (k ==                            1) then
          thisMassFunction%massesLogarithmicMinimum(k)=thisMassFunction%massesLogarithmic(k)-0.5d0*(thisMassFunction%massesLogarithmic(k+1)-thisMassFunction%massesLogarithmic(k  ))
       else
          thisMassFunction%massesLogarithmicMinimum(k)=                                     +0.5d0*(thisMassFunction%massesLogarithmic(k-1)+thisMassFunction%massesLogarithmic(k  ))
       end if
       if (k == thisMassFunction%massesCount) then
          thisMassFunction%massesLogarithmicMaximum(k)=thisMassFunction%massesLogarithmic(k)+0.5d0*(thisMassFunction%massesLogarithmic(k  )-thisMassFunction%massesLogarithmic(k-1))
       else
          thisMassFunction%massesLogarithmicMaximum(k)=                                     +0.5d0*(thisMassFunction%massesLogarithmic(k+1)+thisMassFunction%massesLogarithmic(k  ))
       end if
    end do
    return
  end subroutine Load_PRIMUS_Mass_Function

end module Galacticus_Output_Analyses_Mass_Functions
