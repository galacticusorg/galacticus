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

!% Contains a module which performs analysis to compute a variety of mass-dependent fraction functions.

module Galacticus_Output_Analyses_Mass_Dpndnt_Fractions
  !% Performs analysis to compute a variety of mass-dependent fraction functions. Currently supported fraction functions include:
  !% \begin{itemize}
  !% \item The \gls{gama} early-type fraction \cite{kelvin_galaxy_2014}.
  !% \end{itemize}
  use, intrinsic :: ISO_C_Binding
  use Galacticus_Nodes
  use FGSL
  use Tables
  use Galactic_Structure_Options
  use Geometry_Surveys
  use Numerical_Constants_Astronomical
  use Numerical_Constants_Prefixes
  implicit none
  private
  public :: Galacticus_Output_Analysis_Mass_Dpndnt_Fractions, Galacticus_Output_Analysis_Mass_Dpndnt_Fractions_Output
  
  ! Record of module initialization.
  logical                                                       :: moduleInitialized              =.false.

  ! Record of whether this analysis is active.
  logical                                                       :: analysisActive

  ! Number of supported mass functions.
  integer          , parameter                                  :: fractionFunctionsSupportedCount=1

  ! Labels for supported metallicity distributions.
  character(len=26), dimension(fractionFunctionsSupportedCount) :: fractionFunctionLabels         =        &
       & [                                                                                                 &
       &  'gamaEarlyTypeFractionZ0.03'                                                                     &
       & ]

  ! Interface for mass mapping functions.
  abstract interface
     double precision function Map_Mass(mass,node)
       import treeNode
       double precision          , intent(in   )          :: mass
       type            (treeNode), intent(inout), pointer :: node
     end function Map_Mass
  end interface

  ! Interface for metallicity mapping functions.
  abstract interface
     double precision function Map_Metallicity(radius,node)
       import treeNode
       double precision          , intent(in   )          :: radius
       type            (treeNode), intent(inout), pointer :: node
     end function Map_Metallicity
  end interface

  ! Interface for mass error functions.
  abstract interface
     double precision function Mass_Error(mass,node)
       import treeNode
       double precision          , intent(in   )          :: mass
       type            (treeNode), intent(inout), pointer :: node
     end function Mass_Error
  end interface

  ! Interface for classifier functions.
  abstract interface
     double precision function Classifier(node)
       import treeNode
       type(treeNode), intent(inout), pointer :: node
     end function Classifier
  end interface

  ! Type for descriptors of metallicity distributions.
  type :: fractionFunctionDescriptor
     double precision                                           :: massSystematicLogM0
     procedure       (Mass_Error         ), pointer    , nopass :: massRandomErrorFunction
     procedure       (Classifier         ), pointer    , nopass :: classifierFunction
     double precision                                           :: massLogarithmicMinimum
     integer                                                    :: massSystematicCoefficientCount, massRandomCoefficientCount
     integer                                                    :: massType
     double precision                                           :: massUnitsInSI
     character       (len= 34            )                      :: label
     character       (len=128            )                      :: comment
     procedure       (Map_Mass           ), pointer    , nopass :: mapMass
     class           (surveyGeometryClass), allocatable         :: geometry
  end type fractionFunctionDescriptor

  ! Metallicity distribution descriptors.
  type(fractionFunctionDescriptor), dimension(fractionFunctionsSupportedCount), target :: fractionFunctionDescriptors= &
       & [                                                                                                             &
       ! GAMA early-type galaxy fraction from Kelvin et al. (2014). 
       &                           fractionFunctionDescriptor(                                                         &
       &                                                             11.0000d+0                                      , &
       &                                                              null()                                         , &
       &                                                              null()                                         , &
       &                                                              7.0d0                                          , &
       &                                                              2                                              , &
       &                                                              2                                              , &
       &                                                              massTypeStellar                                , &
       &                                                              massSolar                                      , &
       &                                                              'gamaEarlyTypeFractionZ0.03'                   , &
       &                                                              'GAMA early-type fraction z=0.03'              , &
       &                                                              null()                                         , &
       &                                                              null()                                           &
       &                                                            )                                                  &
       & ]

  ! Type to store metallicity distributions.
  type :: fractionFunction
     ! Copy of the mass function descriptor for this mass function.
     type            (fractionFunctionDescriptor), pointer                     :: descriptor
     ! Parameters for the systematic error model.
     double precision                            , allocatable, dimension(:  ) :: massSystematicCoefficients, massRandomCoefficients
     ! The number of bins.
     integer                                                                   :: massesCount
     ! Arrays for the masses, radii and size function.
     double precision                            , allocatable, dimension(:  ) :: masses                    , massesLogarithmic             , &
          &                                                                       massesLogarithmicMinimum  , massesLogarithmicMaximum      , &
          &                                                                       fractionFunctionWeights
     double precision                            , allocatable, dimension(:  ) :: fractionFunction
     double precision                            , allocatable, dimension(:,:) :: outputWeight
    ! Arrays for accumulation of of main branch galaxies
     double precision                            , allocatable, dimension(:,:) :: mainBranchGalaxyWeights   , mainBranchGalaxyWeightsSquared
     ! Array for the covariance matrix.
     double precision                            , allocatable, dimension(:,:) :: fractionFunctionCovariance
     ! Cosmology conversion factors.
     double precision                            , allocatable, dimension(:  ) :: cosmologyConversionMass
  end type fractionFunction

  ! Mass functions.
  type(fractionFunction), allocatable, dimension(:) :: fractionFunctions

  ! Type for storing temporary metallicity functions during cumulation.
  type :: fractionFunctionWork
     double precision, allocatable, dimension(:) :: fractionFunction, fractionFunctionWeights
  end type fractionFunctionWork

  ! Work array.
  type(fractionFunctionWork), allocatable, dimension(:) :: galaxyWork
  !$omp threadprivate(galaxyWork)

  ! Options controlling binning in halo mass.
  integer                     :: analysisFractionFunctionCovarianceModel
  integer         , parameter :: analysisFractionFunctionCovarianceModelPoisson       =1
  integer         , parameter :: analysisFractionFunctionCovarianceModelBinomial      =2
  integer                     :: analysisFractionFunctionsHaloMassBinsCount             , analysisFractionFunctionsHaloMassBinsPerDecade
  double precision            :: analysisFractionFunctionsHaloMassMinimum               , analysisFractionFunctionsHaloMassMaximum           , &
       &                         analysisFrctnFnctnsHaloMassIntervalLogarithmicInverse  , analysisFractionFunctionsHaloMassMinimumLogarithmic

contains

  !# <mergerTreeAnalysisTask>
  !#  <unitName>Galacticus_Output_Analysis_Mass_Dpndnt_Fractions</unitName>
  !# </mergerTreeAnalysisTask>
  subroutine Galacticus_Output_Analysis_Mass_Dpndnt_Fractions(tree,node,nodeStatus,iOutput,mergerTreeAnalyses)
    !% Construct fraction functions to compare to various observational determinations.
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
    implicit none
    type            (mergerTree                    ), intent(in   )                 :: tree
    type            (treeNode                      ), intent(inout), pointer        :: node
    integer                                         , intent(in   )                 :: nodeStatus
    integer         (c_size_t                      ), intent(in   )                 :: iOutput
    type            (varying_string                ), intent(in   ), dimension(:  ) :: mergerTreeAnalyses
    class           (nodeComponentBasic            )               , pointer        :: basic
    class           (cosmologyFunctionsClass       )               , pointer        :: cosmologyFunctionsModel
    type            (cosmologyFunctionsMatterLambda)                                :: cosmologyFunctionsObserved
    type            (cosmologyParametersSimple     )               , pointer        :: cosmologyParametersObserved
    integer         (c_size_t                      )                                :: k                                 , jOutput
    double precision                                , parameter                     :: massRandomErrorMinimum     =1.0d-3
    integer                                                                         :: i                                 , j              , &
         &                                                                             l                                 , currentAnalysis, &
         &                                                                             activeAnalysisCount               , haloMassBin
    double precision                                                                :: dataHubbleParameter ,mass,massLogarithmic&
         &,massRandomError,dataOmegaDarkEnergy,dataOmegaMatter,redshift,timeMinimum,timeMaximum,distanceMinimum,distanceMaximum,probabilityInClass
    type            (varying_string                )                                :: parameterName&
         &,analysisFractionFunctionCovarianceModelText,cosmologyScalingMass,message
    type            (hdf5Object                    )                                :: dataFile,massDataset,parametersGroup
    
    ! Initialize the module if necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Galacticus_Output_Analysis_Mass_Dpndnt_Fractions_Initialize)
       if (.not.moduleInitialized) then
          !@ <inputParameter>
          !@   <name>analysisFractionFunctionCovarianceModel</name>
          !@   <defaultValue>binomial</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The model to use when constructing the fraction functions covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisFractionFunctionCovarianceModel',analysisFractionFunctionCovarianceModelText,defaultValue='Poisson')
          select case (char(analysisFractionFunctionCovarianceModelText))
          case ( 'Poisson'  )
             analysisFractionFunctionCovarianceModel=analysisFractionFunctionCovarianceModelPoisson
          case ( 'binomial' )
             analysisFractionFunctionCovarianceModel=analysisFractionFunctionCovarianceModelBinomial
          case default
             call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Dpndnt_Fractions','unrecognized value for "analysisFractionFunctionCovarianceModel" - allowed values are "Poisson", and "binomial"')
          end select
          !@ <inputParameter>
          !@   <name>analysisFractionFunctionsHaloMassBinsPerDecade</name>
          !@   <defaultValue>10</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The number of bins per decade of halo mass to use when constructing the mass function covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisFractionFunctionsHaloMassBinsPerDecade',analysisFractionFunctionsHaloMassBinsPerDecade,defaultValue=10)
          !@ <inputParameter>
          !@   <name>analysisFractionFunctionsHaloMassMinimum</name>
          !@   <defaultValue>$10^8M_\odot$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The minimum halo mass to consider when constructing the mass function covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisFractionFunctionsHaloMassMinimum',analysisFractionFunctionsHaloMassMinimum,defaultValue=1.0d8)
          !@ <inputParameter>
          !@   <name>analysisFractionFunctionsHaloMassMaximum</name>
          !@   <defaultValue>$10^{16}M_\odot$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The maximum halo mass to consider when constructing the mass function covariance matrix for main branch galaxies.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('analysisFractionFunctionsHaloMassMaximum',analysisFractionFunctionsHaloMassMaximum,defaultValue=1.0d16)
          analysisFractionFunctionsHaloMassMinimumLogarithmic  =     log10( analysisFractionFunctionsHaloMassMinimum)
          analysisFractionFunctionsHaloMassBinsCount           =int(                                                      &
               &                                                     log10(                                               &
               &                                                            analysisFractionFunctionsHaloMassMaximum      &
               &                                                           /analysisFractionFunctionsHaloMassMinimum      &
               &                                                          )                                               &
               &                                                    *dble(analysisFractionFunctionsHaloMassBinsPerDecade) &
               &                                                    +0.5d0                                                &
               &                                                   )
          analysisFrctnFnctnsHaloMassIntervalLogarithmicInverse=+dble     (analysisFractionFunctionsHaloMassBinsCount)    &
               &                                                /    log10(                                               &
               &                                                           +analysisFractionFunctionsHaloMassMaximum      &
               &                                                           /analysisFractionFunctionsHaloMassMinimum      &
               &                                                          )
          ! Establish classifier functions for fraction function descriptors.
          fractionFunctionDescriptors(1)%classifierFunction => Classifier_GAMA_Early_Type_Fraction_Z0_03
          ! Determine how many supported fraction functions are requested.
          activeAnalysisCount=0
          do i=1,fractionFunctionsSupportedCount
             if (any(trim(mergerTreeAnalyses) == trim(fractionFunctionLabels(i)))) activeAnalysisCount=activeAnalysisCount+1
          end do
          ! Allocate mass function arrays and populate with required data.
          if (activeAnalysisCount <= 0) then
             analysisActive=.false.
          else
             analysisActive=.true.
             ! Establish survey geometries.
             allocate(surveyGeometryKelvin2014GAMAnear :: fractionFunctionDescriptors(1)%geometry)
             select type (g => fractionFunctionDescriptors(1)%geometry)
             type is (surveyGeometryKelvin2014GAMAnear)
                g=surveyGeometryKelvin2014GAMAnear()
             end select
             ! Initialize analyses.
             currentAnalysis=0
             allocate(fractionFunctions(activeAnalysisCount))
             cosmologyFunctionsModel => cosmologyFunctions()
             do i=1,size(mergerTreeAnalyses)
                do j=1,fractionFunctionsSupportedCount
                   if (mergerTreeAnalyses(i) == trim(fractionFunctionLabels(j))) then
                      currentAnalysis=currentAnalysis+1
                      ! Set a pointer to the descriptor for this size function.
                      fractionFunctions(currentAnalysis)%descriptor => fractionFunctionDescriptors(j)
                      ! Read parameters of the systematic error model.
                      if (fractionFunctionDescriptors(j)%massSystematicCoefficientCount > 0) then
                         allocate(fractionFunctions(currentAnalysis)%massSystematicCoefficients(fractionFunctionDescriptors(j)%massSystematicCoefficientCount))
                         do k=1,fractionFunctionDescriptors(j)%massSystematicCoefficientCount
                            parameterName=trim(fractionFunctionLabels(j))//'MassSystematic'
                            parameterName=parameterName//(k-1)
                            !@ <inputParameter>
                            !@   <regEx>(gamaEarlyTypeFraction)Z[0-9\.]+MassSystematic[0-9]+</regEx>
                            !@   <defaultValue>0</defaultValue>
                            !@   <attachedTo>module</attachedTo>
                            !@   <description>
                            !@     Mass-dependent metallicity distribution mass systematic parameters.
                            !@   </description>
                            !@   <type>real</type>
                            !@   <cardinality>1</cardinality>
                            !@ </inputParameter>
                            call Get_Input_Parameter(char(parameterName),fractionFunctions(currentAnalysis)%massSystematicCoefficients(k),defaultValue=0.0d0)
                         end do
                      end if
                      ! Read parameters of the random error model.
                      if (fractionFunctionDescriptors(j)%massRandomCoefficientCount > 0) then
                         allocate(fractionFunctions(currentAnalysis)%massRandomCoefficients(fractionFunctionDescriptors(j)%massRandomCoefficientCount))
                         do k=1,fractionFunctionDescriptors(j)%massRandomCoefficientCount
                            parameterName=trim(fractionFunctionLabels(j))//'MassRandom'
                            parameterName=parameterName//(k-1)
                            !@ <inputParameter>
                            !@   <regEx>(gamaEarlyTypeFraction)Z[0-9\.]+MassRandom[0-9]+</regEx>
                            !@   <defaultValue>0</defaultValue>
                            !@   <attachedTo>module</attachedTo>
                            !@   <description>
                            !@     Mass-dependent metallicity distribution mass random parameters.
                            !@   </description>
                            !@   <type>real</type>
                            !@   <cardinality>1</cardinality>
                            !@ </inputParameter>
                            call Get_Input_Parameter(char(parameterName),fractionFunctions(currentAnalysis)%massRandomCoefficients(k),defaultValue=0.0d0)
                         end do
                      end if
                      ! Read the appropriate observational data definition.
                      select case (trim(fractionFunctionLabels(j)))
                      case ('gamaEarlyTypeFractionZ0.03')
                         ! GAMA z=0.03 early-type fraction.
                         !$omp critical(HDF5_Access)
                         call dataFile%openFile(char(Galacticus_Input_Path())//"data/observations/morphology/earlyTypeFractionGAMA.hdf5",readOnly=.true.)
                         call dataFile%readDataset("mass",fractionFunctions(currentAnalysis)%masses)
                         massDataset=dataFile%openDataset("mass")
                         call massDataset%readAttribute("cosmologyScaling",cosmologyScalingMass)
                         call massDataset%close        (                                       )
                         fractionFunctions(currentAnalysis)%massesCount=size(fractionFunctions(currentAnalysis)%masses)
                         call allocateArray(fractionFunctions(currentAnalysis)%massesLogarithmic         ,[fractionFunctions(currentAnalysis)%massesCount])
                         call allocateArray(fractionFunctions(currentAnalysis)%massesLogarithmicMinimum  ,[fractionFunctions(currentAnalysis)%massesCount])
                         call allocateArray(fractionFunctions(currentAnalysis)%massesLogarithmicMaximum  ,[fractionFunctions(currentAnalysis)%massesCount])
                         call allocateArray(fractionFunctions(currentAnalysis)%fractionFunction          ,[fractionFunctions(currentAnalysis)%massesCount])
                         call allocateArray(fractionFunctions(currentAnalysis)%fractionFunctionWeights   ,[fractionFunctions(currentAnalysis)%massesCount])
                         call allocateArray(fractionFunctions(currentAnalysis)%fractionFunctionCovariance,[                                                &
                              &                                                                          fractionFunctions(currentAnalysis)%massesCount, &
                              &                                                                          fractionFunctions(currentAnalysis)%massesCount  &
                              &                                                                         ]                                                &
                              &          )
                         call allocateArray(fractionFunctions(currentAnalysis)%mainBranchGalaxyWeights       ,[fractionFunctions(currentAnalysis)%massesCount,analysisFractionFunctionsHaloMassBinsCount])
                         call allocateArray(fractionFunctions(currentAnalysis)%mainBranchGalaxyWeightsSquared,[fractionFunctions(currentAnalysis)%massesCount,analysisFractionFunctionsHaloMassBinsCount])
                         fractionFunctions              (currentAnalysis)%       massesLogarithmic  &
                              & =log10(fractionFunctions(currentAnalysis)%       masses)
                         do k=1,fractionFunctions(currentAnalysis)%massesCount
                            if (k ==                                                            1) then
                               fractionFunctions                (currentAnalysis)%massesLogarithmicMinimum(k  )= &
                                    & +        fractionFunctions(currentAnalysis)%massesLogarithmic       (k  )  &
                                    & -0.5d0*(                                                                   &
                                    &         +fractionFunctions(currentAnalysis)%massesLogarithmic       (k+1)  &
                                    &         -fractionFunctions(currentAnalysis)%massesLogarithmic       (k  )  &
                                    &        )
                            else
                               fractionFunctions                (currentAnalysis)%massesLogarithmicMinimum(k  )= &
                                    & +0.5d0*(                                                                   &
                                    &         +fractionFunctions(currentAnalysis)%massesLogarithmic       (k-1)  &
                                    &         +fractionFunctions(currentAnalysis)%massesLogarithmic       (k  )  &
                                    &        )
                            end if
                            if (k == fractionFunctions(currentAnalysis)%massesCount) then
                               fractionFunctions                (currentAnalysis)%massesLogarithmicMaximum(k  )= &
                                    & +        fractionFunctions(currentAnalysis)%massesLogarithmic       (k  )  &
                                    & +0.5d0*(                                                                   &
                                    &         +fractionFunctions(currentAnalysis)%massesLogarithmic       (k  )  &
                                    &         -fractionFunctions(currentAnalysis)%massesLogarithmic       (k-1)  &
                                    &        )
                            else
                               fractionFunctions                (currentAnalysis)%massesLogarithmicMaximum(k  )= &
                                    & +0.5d0*(                                                                   &
                                    &         +fractionFunctions(currentAnalysis)%massesLogarithmic       (k+1)  &
                                    &         +fractionFunctions(currentAnalysis)%massesLogarithmic       (k  )  &
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
                         fractionFunctions(currentAnalysis)%fractionFunction              =0.0d0
                         fractionFunctions(currentAnalysis)%fractionFunctionWeights       =0.0d0
                         fractionFunctions(currentAnalysis)%fractionFunctionCovariance    =0.0d0
                         fractionFunctions(currentAnalysis)%mainBranchGalaxyWeights       =0.0d0
                         fractionFunctions(currentAnalysis)%mainBranchGalaxyWeightsSquared=0.0d0
                      case default
                         call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Dpndnt_Fractions','unknown size function')
                      end select
                      ! Get cosmological conversion factors.
                      call allocateArray(fractionFunctions(currentAnalysis)%cosmologyConversionMass,[Galacticus_Output_Time_Count()])
                      do jOutput=1,Galacticus_Output_Time_Count()
                         redshift=                                                                                      &
                              &   cosmologyFunctionsModel %redshiftFromExpansionFactor(                                 &
                              &    cosmologyFunctionsModel%expansionFactor             (                                &
                              &                                                         Galacticus_Output_Time(jOutput) &
                              &                                                        )                                &
                              &                                                       )
                         call Cosmology_Conversion_Factors(                                                                                                &
                              &                            redshift                                                                                      , &
                              &                            cosmologyFunctionsModel                                                                       , &
                              &                            cosmologyFunctionsObserved                                                                    , &
                              &                            cosmologyScalingMass      =cosmologyScalingMass                                               , &
                              &                            cosmologyConversionMass   =fractionFunctions(currentAnalysis)%cosmologyConversionMass(jOutput)  &
                              &                           )
                      end do
                      nullify(cosmologyParametersObserved)
                      ! Compute output weights for metallicity distribution.
                      call allocateArray(fractionFunctions(currentAnalysis)%outputWeight,[int(fractionFunctions(currentAnalysis)%massesCount,kind=c_size_t),Galacticus_Output_Time_Count()])
                      fractionFunctions(currentAnalysis)%outputWeight=0.0d0
                      do k=1,fractionFunctions(currentAnalysis)%massesCount
                         do jOutput=1,Galacticus_Output_Time_Count()
                            do l=1,fractionFunctions(currentAnalysis)%descriptor%geometry%fieldCount()
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
                               distanceMinimum=max(                                                                                                                        &
                                    &              cosmologyFunctionsModel%distanceComoving(timeMaximum)                                                                 , &
                                    &              fractionFunctions(currentAnalysis)%descriptor%geometry%distanceMinimum(fractionFunctions(currentAnalysis)%masses(k),l)  &
                                    &             )
                               distanceMaximum=min(                                                                                                                        &
                                    &              cosmologyFunctionsModel%distanceComoving(timeMinimum)                                                                 , &
                                    &              fractionFunctions(currentAnalysis)%descriptor%geometry%distanceMaximum(fractionFunctions(currentAnalysis)%masses(k),l)  &
                                    &             )
                               fractionFunctions        (currentAnalysis)%outputWeight                    (k,jOutput)  &
                                    & =fractionFunctions(currentAnalysis)%outputWeight                    (k,jOutput)  &
                                    & +fractionFunctions(currentAnalysis)%descriptor  %geometry%solidAngle(  l      )  &
                                    & /3.0d0                                                                           &
                                    & *                                                                                &
                                    & max(                                                                             &
                                    &     +0.0d0                                                                     , &
                                    &     +distanceMaximum**3                                                          &
                                    &     -distanceMinimum**3                                                          &
                                    &    )
                            end do
                         end do
                         where(fractionFunctions(currentAnalysis)%outputWeight(k,:) < 0.0d0)
                            fractionFunctions(currentAnalysis)%outputWeight(k,:)=0.0d0
                         end where
                         if (any(fractionFunctions(currentAnalysis)%outputWeight(k,:) > 0.0d0)) then
                            fractionFunctions                  (currentAnalysis)%outputWeight(k,:)  &
                                 &       =    fractionFunctions(currentAnalysis)%outputWeight(k,:)  &
                                 &       /sum(fractionFunctions(currentAnalysis)%outputWeight(k,:))
                         else
                            message="fraction function '"//trim(fractionFunctions(currentAnalysis)%descriptor%label)//"' bin "
                            message=message//k//" has zero weights"
                            call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Dpndnt_Fractions',message)
                         end if
                      end do
                      exit
                   end if
                end do
             end do
          end if
          ! Record that module is initialized.
          moduleInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Analysis_Mass_Dpndnt_Fractions_Initialize)
    end if
    ! Return if this analysis is not active.
    if (.not.analysisActive                   ) return
    ! Return if this is a tree finalization.
    if (nodeStatus          == nodeStatusFinal) return
    ! Allocate work arrays.
    if (.not.allocated(galaxyWork)) allocate(galaxyWork(size(fractionFunctions)))
    ! Iterate over active analyses.
    do i=1,size(fractionFunctions)
       ! Cycle if this fraction function receives no contribution from this output.
       if (all(fractionFunctions(i)%outputWeight(:,iOutput) <= 0.0d0)) cycle
       ! Allocate workspace.
       if (.not.allocated(galaxyWork(i)%fractionFunction       )) call allocateArray(galaxyWork(i)%fractionFunction       ,[fractionFunctions(i)%massesCount])
       if (.not.allocated(galaxyWork(i)%fractionFunctionWeights)) call allocateArray(galaxyWork(i)%fractionFunctionWeights,[fractionFunctions(i)%massesCount])
       ! Get the galactic mass.
       mass=                                                                                                                                                &
            &  Galactic_Structure_Enclosed_Mass(node,radiusLarge,componentType=componentTypeDisk    ,massType=fractionFunctions(i)%descriptor%massType) &
            & +Galactic_Structure_Enclosed_Mass(node,radiusLarge,componentType=componentTypeSpheroid,massType=fractionFunctions(i)%descriptor%massType)
       if (mass <= 0.0d0) cycle
       if (associated(fractionFunctions(i)%descriptor%mapMass)) mass=fractionFunctions(i)%descriptor%mapMass(mass,node)
       ! Convert mass for cosmology and systematics.
       mass=mass*fractionFunctions(i)%cosmologyConversionMass(iOutput)
       massLogarithmic=log10(mass)
       do j=1,fractionFunctions(i)%descriptor%massSystematicCoefficientCount
          massLogarithmic=massLogarithmic+fractionFunctions(i)%massSystematicCoefficients(j)*(log10(mass)-fractionFunctions(i)%descriptor%massSystematicLogM0)**(j-1)
       end do
       if (massLogarithmic < fractionFunctions(i)%descriptor%massLogarithmicMinimum) cycle
       ! Compute random errors on mass.
       if (associated(fractionFunctions(i)%descriptor%massRandomErrorFunction)) then
          massRandomError=fractionFunctions(i)%descriptor%massRandomErrorFunction(mass,node)
       else
          massRandomError=0.0d0
          do j=1,fractionFunctions(i)%descriptor%massRandomCoefficientCount
             massRandomError=+massRandomError                                       &
                  &          +fractionFunctions(i)%massRandomCoefficients(j)        &
                  &          *(                                                     &
                  &            +log10(mass)                                         &
                  &            -fractionFunctions(i)%descriptor%massSystematicLogM0 &
                  &           )**(j-1)
          end do
          massRandomError=max(massRandomError,massRandomErrorMinimum)
       end if
       ! Compute probability for galaxy to be in class.
       probabilityInClass=fractionFunctions(i)%descriptor%classifierFunction(node)
       ! Compute contributions to each bin.
       galaxyWork(i)%fractionFunctionWeights=+(                                                                                                  &
            &                                  +erf((fractionFunctions(i)%massesLogarithmicMaximum-massLogarithmic)/massRandomError/sqrt(2.0d0)) &
            &                                  -erf((fractionFunctions(i)%massesLogarithmicMinimum-massLogarithmic)/massRandomError/sqrt(2.0d0)) &
            &                                 )                                                                                                  &
            &                                /2.0d0                                                                                              &
            &                                *tree%volumeWeight                                                                              &
            &                                *fractionFunctions(i)%outputWeight(:,iOutput)
       galaxyWork(i)%fractionFunction       =+probabilityInClass                                                                                 &
            &                                *galaxyWork(i)%fractionFunctionWeights
       ! Accumulate fraction function.
       if (any(galaxyWork(i)%fractionFunction /= 0.0d0)) then
          !$omp critical (Galacticus_Output_Analysis_Mass_Dpndnt_Fractions_Accumulate)
          fractionFunctions(i)%fractionFunction       =fractionFunctions(i)%fractionFunction       +galaxyWork(i)%fractionFunction
          fractionFunctions(i)%fractionFunctionWeights=fractionFunctions(i)%fractionFunctionWeights+galaxyWork(i)%fractionFunctionWeights
          !$omp end critical (Galacticus_Output_Analysis_Mass_Dpndnt_Fractions_Accumulate)
          ! Treat main branch and other galaxies differently.
          if (node%isOnMainBranch().and.analysisFractionFunctionCovarianceModel == analysisFractionFunctionCovarianceModelBinomial) then
             ! Find the bin to which this halo mass belongs.
             basic => node%basic()
             haloMassBin=floor((log10(basic%mass())-analysisFractionFunctionsHaloMassMinimumLogarithmic)*analysisFrctnFnctnsHaloMassIntervalLogarithmicInverse)+1
             ! Accumulate weights to halo mass arrays.
             if (haloMassBin >= 1 .and. haloMassBin <= analysisFractionFunctionsHaloMassBinsCount) then
                !$omp critical (Galacticus_Output_Analysis_Mass_Dpndnt_Fractions_Accumulate)
                fractionFunctions        (i)%mainBranchGalaxyWeights       (:,haloMassBin)= &
                     &  fractionFunctions(i)%mainBranchGalaxyWeights       (:,haloMassBin)  &
                     &  +galaxyWork  (i)%fractionFunction
                fractionFunctions        (i)%mainBranchGalaxyWeightsSquared(:,haloMassBin)= &
                     &  fractionFunctions(i)%mainBranchGalaxyWeightsSquared(:,haloMassBin)  &
                     &  +galaxyWork  (i)%fractionFunction**2
                !$omp end critical (Galacticus_Output_Analysis_Mass_Dpndnt_Fractions_Accumulate)
             end if
          else
             !$omp critical (Galacticus_Output_Analysis_Mass_Dpndnt_Fractions_Accumulate)
             call Vector_Outer_Product_Accumulate(                                                 &
                  &                               galaxyWork       (i)%fractionFunction          , &
                  &                               fractionFunctions(i)%fractionFunctionCovariance, &
                  &                               sparse=.true.                                    &
                  &                              )
             !$omp end critical (Galacticus_Output_Analysis_Mass_Dpndnt_Fractions_Accumulate)
           end if
       end if
    end do
    return
  end subroutine Galacticus_Output_Analysis_Mass_Dpndnt_Fractions

  !# <hdfPreCloseTask>
  !#  <unitName>Galacticus_Output_Analysis_Mass_Dpndnt_Fractions_Output</unitName>
  !# </hdfPreCloseTask>
  subroutine Galacticus_Output_Analysis_Mass_Dpndnt_Fractions_Output
    !% Outputs fraction functions to file.
    use Galacticus_HDF5
    use Vectors
    implicit none
    integer                      :: k                 , m                    , &
         &                          mi                , mj
    type            (hdf5Object) :: analysisGroup     , fractionFunctionGroup, &
         &                          dataset
    double precision             :: haloWeightBinTotal

    ! Return immediately if this analysis is not active.
    if (.not.analysisActive) return
    ! Iterate over mass functions.
    do k=1,size(fractionFunctions)
       ! Symmetrize the covariance matrix (we've accumulated only the upper triangle).
       fractionFunctions(k)%fractionFunctionCovariance=Matrix_Copy_Upper_To_Lower_Triangle(fractionFunctions(k)%fractionFunctionCovariance)
       ! Add the contribution from main branch galaxies to the covariance matrix.
       if (analysisFractionFunctionCovarianceModel == analysisFractionFunctionCovarianceModelBinomial) then
          do m=1,analysisFractionFunctionsHaloMassBinsCount
             haloWeightBinTotal=sum(fractionFunctions(k)%mainBranchGalaxyWeights(:,m))
             if ( haloWeightBinTotal > 0.0 ) then
                do mi=1,fractionFunctions(k)%massesCount
                   fractionFunctions               (k)%fractionFunctionCovariance    (mi,mi)=                    &
                        &         fractionFunctions(k)%fractionFunctionCovariance    (mi,mi)                     &
                        & +(1.0d0-fractionFunctions(k)%mainBranchGalaxyWeights       (mi,m )/haloWeightBinTotal) &
                        & *       fractionFunctions(k)%mainBranchGalaxyWeightsSquared(mi,m )
                   do mj=1,fractionFunctions(k)%massesCount
                      if (mi == mj) cycle
                      fractionFunctions         (k)%fractionFunctionCovariance    (mi,mj)=                    &
                           &   fractionFunctions(k)%fractionFunctionCovariance    (mi,mj)                     &
                           & -(fractionFunctions(k)%mainBranchGalaxyWeights       (mj,m )/haloWeightBinTotal) &
                           & * fractionFunctions(k)%mainBranchGalaxyWeightsSquared(mi,m )
                   end do
                end do
             end if
          end do
       end if
       ! Normalize the fraction function in each mass interval.
       do mi=1,fractionFunctions(k)%massesCount 
          if (fractionFunctions(k)%fractionFunctionWeights(mi) > 0.0d0) then
             fractionFunctions        (k)%fractionFunction       (mi) &
                  & =fractionFunctions(k)%fractionFunction       (mi) &
                  & /fractionFunctions(k)%fractionFunctionWeights(mi)
             do mj=1,fractionFunctions(k)%massesCount
                if (fractionFunctions(k)%fractionFunctionWeights(mj) > 0.0d0) then
                   fractionFunctions        (k)%fractionFunctionCovariance(mi,mj) &
                        & =fractionFunctions(k)%fractionFunctionCovariance(mi,mj) &
                        & /fractionFunctions(k)%fractionFunctionWeights   (   mi) &
                        & /fractionFunctions(k)%fractionFunctionWeights   (   mj)
                end if
             end do
          end if
       end do
       ! Output the size function.
       !$omp critical(HDF5_Access)
       analysisGroup        =galacticusOutputFile%openGroup('analysis','Model analysis')
       fractionFunctionGroup=analysisGroup       %openGroup(trim(fractionFunctions(k)%descriptor%label),trim(fractionFunctions(k)%descriptor%comment))
       call fractionFunctionGroup%writeDataset  (fractionFunctions(k)%masses                    ,'mass'                      ,'Mass'                        ,datasetReturned=dataset)
       call dataset              %writeAttribute(fractionFunctions(k)%descriptor%massUnitsInSI  ,'unitsInSI'                                                                        )
       call dataset              %close()
       call fractionFunctionGroup%writeDataset  (fractionFunctions(k)%fractionFunction          ,'fractionFunction'          ,'Fraction function'                                   )
       call fractionFunctionGroup%writeDataset  (fractionFunctions(k)%fractionFunctionCovariance,'fractionFunctionCovariance','Fraction function covariance'                        )
       call fractionFunctionGroup%close()
       call analysisGroup        %close()
       !$omp end critical(HDF5_Access)
    end do
    return
  end subroutine Galacticus_Output_Analysis_Mass_Dpndnt_Fractions_Output

  double precision function Classifier_GAMA_Early_Type_Fraction_Z0_03(node)
    !% Classifies galaxies as early-type.
    use Input_Parameters
    use Error_Functions
    implicit none
    type            (treeNode             ), intent(inout), pointer :: node
    class           (nodeComponentDisk    )               , pointer :: disk
    class           (nodeComponentSpheroid)               , pointer :: spheroid
    double precision                       , save                   :: gamaEarlyTypeRatio                 , gamaEarlyTypeRatioError, &
         &                                                             gamaEarlyTypeLBSProbability
    logical                                , save                   :: initialized                =.false.
    double precision                                                :: massDisk                           , massSpheroid           , &
         &                                                             ratio
    
    ! Initialize if necessary.
    if (.not.initialized) then
       !$omp critical (Classifier_GAMA_Early_Type_Fraction_Z0_03_Initialize)
       if (.not.initialized) then
          !@ <inputParameter>
          !@   <name>gamaEarlyTypeRatio</name>
          !@   <defaultValue>0.5</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The minimum spheroid-to-total ratio for a galaxy to be classified as ``early-type'' when constructing the \gls{gama} early-type fraction function.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('gamaEarlyTypeRatio',gamaEarlyTypeRatio,defaultValue=0.5d0)
          !@ <inputParameter>
          !@   <name>gamaEarlyTypeRatioError</name>
          !@   <defaultValue>0.3</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The error in spheriod fraction to be used when constructing the \gls{gama} early-type fraction function.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('gamaEarlyTypeRatioError',gamaEarlyTypeRatioError,defaultValue=0.3d0)
          !@ <inputParameter>
          !@   <name>gamaEarlyTypeLBSProbability</name>
          !@   <defaultValue>0.5</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The probability that the LBS (``little blue spheroid'') class belongs to the early-type class when constructing the \gls{gama} early-type fraction function.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>0..1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('gamaEarlyTypeLBSProbability',gamaEarlyTypeLBSProbability,defaultValue=0.5d0)
          ! Record that function is initialized.
          initialized=.true.
       end if
       !$omp end critical (Classifier_GAMA_Early_Type_Fraction_Z0_03_Initialize)
    end if
    ! Classify.
    disk         => node    %disk       ()
    spheroid     => node    %spheroid   ()
    massDisk     =  disk    %massStellar()
    massSpheroid =  spheroid%massStellar()
    if (massDisk+massSpheroid > 0.0d0) then
       ratio=massSpheroid/(massDisk+massSpheroid)
       Classifier_GAMA_Early_Type_Fraction_Z0_03=                                                &
            & +(                                                                                 &
            &   +Error_Function((+1.0d0             -ratio)/gamaEarlyTypeRatioError/sqrt(2.0d0)) &
            &   +Error_Function((-gamaEarlyTypeRatio+ratio)/gamaEarlyTypeRatioError/sqrt(2.0d0)) &
            &  )                                                                                 &
            & /(                                                                                 &
            &   +Error_Function((+1.0d0             -ratio)/gamaEarlyTypeRatioError/sqrt(2.0d0)) &
            &   +Error_Function(                    +ratio /gamaEarlyTypeRatioError/sqrt(2.0d0)) &
            &  )   
    else
       Classifier_GAMA_Early_Type_Fraction_Z0_03=0.0d0
    end if
    return
  end function Classifier_GAMA_Early_Type_Fraction_Z0_03
  
end module Galacticus_Output_Analyses_Mass_Dpndnt_Fractions
