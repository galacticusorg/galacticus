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

!% Contains a module which implements a stellar mass function output analysis class.

  use Geometry_Surveys
  use Cosmology_Functions

  !# <outputAnalysis name="outputAnalysisMassFunctionStellar" defaultThreadPrivate="yes">
  !#  <description>A stellar mass function output analysis class.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisVolumeFunction1D) :: outputAnalysisMassFunctionStellar
     !% A massFunctionStellar output analysis class.
     private
     class(surveyGeometryClass    ), pointer :: surveyGeometry_     => null()
     class(cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
  end type outputAnalysisMassFunctionStellar

  interface outputAnalysisMassFunctionStellar
     !% Constructors for the ``massFunctionStellar'' output analysis class.
     module procedure massFunctionStellarConstructorParameters
     module procedure massFunctionStellarConstructorInternal
     module procedure massFunctionStellarConstructorFile
  end interface outputAnalysisMassFunctionStellar

contains

  function massFunctionStellarConstructorParameters(parameters)
    !% Constructor for the ``massFunctionStellar'' output analysis class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type            (outputAnalysisMassFunctionStellar      )                             :: massFunctionStellarConstructorParameters
    type            (inputParameters                        ), intent(inout)              :: parameters
    class           (galacticFilterClass                    ), pointer                    :: galacticFilter_
    class           (surveyGeometryClass                    ), pointer                    :: surveyGeometry_
    class           (cosmologyFunctionsClass                ), pointer                    :: cosmologyFunctions_
    class           (outputAnalysisDistributionOperatorClass), pointer                    :: outputAnalysisDistributionOperator_
    double precision                                         , dimension(:) , allocatable :: masses
    !# <inputParameterList label="allowedParameterNames" />
    
    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)
    allocate(masses(parameters%count('masses')))
    !# <inputParameter>
    !#   <name>masses</name>
    !#   <source>parameters</source>
    !#   <variable>masses</variable>
    !#   <description>The masses corresponding to bin centers.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="galacticFilter"                     name="galacticFilter_"                     source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctions_"                 source="parameters"/>
    !# <objectBuilder class="outputAnalysisDistributionOperator" name="outputAnalysisDistributionOperator_" source="parameters"/>
    surveyGeometry_ => surveyGeometry()
    massFunctionStellarConstructorParameters=outputAnalysisMassFunctionStellar(masses,galacticFilter_,surveyGeometry_,cosmologyFunctions_,outputAnalysisDistributionOperator_)
    return
  end function massFunctionStellarConstructorParameters

  function massFunctionStellarConstructorFile(fileName,galacticFilter_,surveyGeometry_,cosmologyFunctions_,outputAnalysisDistributionOperator_)
    !% Constructor for the ``massFunctionStellar'' output analysis class which reads bin information from a standard format file.
    use IO_HDF5
    implicit none
    type            (outputAnalysisMassFunctionStellar      )                             :: massFunctionStellarConstructorFile
    character       (len=*                                  ), intent(in   )              :: fileName
    class           (galacticFilterClass                    ), intent(in   ), target      :: galacticFilter_
    class           (surveyGeometryClass                    ), intent(in   ), target      :: surveyGeometry_
    class           (cosmologyFunctionsClass                ), intent(in   ), target      :: cosmologyFunctions_
    class           (outputAnalysisDistributionOperatorClass), intent(in   ), target      :: outputAnalysisDistributionOperator_
    double precision                                         , dimension(:) , allocatable :: masses
    type            (hdf5Object                             )                             :: dataFile
    
    !$omp critical(HDF5_Access)
    call dataFile%openFile   (fileName,readOnly=.true.)
    call dataFile%readDataset('mass'  ,masses         )
    call dataFile%close      (                        )
    !$omp end critical(HDF5_Access)
    ! Construct the object.
    massFunctionStellarConstructorFile=outputAnalysisMassFunctionStellar(masses,galacticFilter_,surveyGeometry_,cosmologyFunctions_,outputAnalysisDistributionOperator_)
    return
  end function massFunctionStellarConstructorFile

  function massFunctionStellarConstructorInternal(masses,galacticFilter_,surveyGeometry_,cosmologyFunctions_,outputAnalysisDistributionOperator_) result(self)
    !% Constructor for the ``massFunctionStellar'' output analysis class which takes a parameter set as input.
    use ISO_Varying_String
    use Memory_Management
    use Galacticus_Output_Times
    use String_Handling
    use Galacticus_Error
    implicit none
    type            (outputAnalysisMassFunctionStellar           )                                :: self
    double precision                                              , intent(in   ), dimension(:  ) :: masses
    class           (galacticFilterClass                         ), intent(in   ), target         :: galacticFilter_
    class           (surveyGeometryClass                         ), intent(in   ), target         :: surveyGeometry_
    class           (cosmologyFunctionsClass                     ), intent(in   ), target         :: cosmologyFunctions_
    class           (outputAnalysisDistributionOperatorClass     ), intent(in   ), target         :: outputAnalysisDistributionOperator_
    type            (outputAnalysisPropertyExtractorMassStellar  )               , pointer        :: outputAnalysisPropertyExtractor_
    type            (outputAnalysisPropertyOperatorLog10         )               , pointer        :: outputAnalysisPropertyOperator_
    type            (outputAnalysisDistributionNormalizerSequence)               , pointer        :: outputAnalysisDistributionNormalizer_
    type            (normalizerList                              )               , pointer        :: normalizerSequence                   , normalizer_
    double precision                                              , allocatable  , dimension(:,:) :: outputWeight
    integer         (c_size_t                                    )                                :: iBin                                 , iOutput
    integer                                                                                       :: iField
    double precision                                                                                 timeMinimum                          , timeMaximum    , &
         &                                                                                           distanceMinimum                      , distanceMaximum
    type            (varying_string                              )                                :: message
    !# <constructorAssign variables="*surveyGeometry_, *cosmologyFunctions_"/>

    ! Compute weights that apply to each output redshift.
    self%binCount=size(masses,kind=c_size_t)
    call allocateArray(outputWeight,[self%binCount,Galacticus_Output_Time_Count()])
    outputWeight=0.0d0
    do iBin=1,self%binCount
       do iOutput=1,Galacticus_Output_Time_Count()
          do iField=1,self%surveyGeometry_%fieldCount()
             if (iOutput == Galacticus_Output_Time_Count()) then
                timeMaximum=     Galacticus_Output_Time(iOutput)
             else
                timeMaximum=sqrt(Galacticus_Output_Time(iOutput)*Galacticus_Output_Time(iOutput+1))
             end if
             if (iOutput ==                              1) then
                timeMinimum=     Galacticus_Output_Time(iOutput)
             else
                timeMinimum=sqrt(Galacticus_Output_Time(iOutput)*Galacticus_Output_Time(iOutput-1))
             end if
             distanceMinimum=max(                                                                &
                  &              self%cosmologyFunctions_%distanceComoving(timeMaximum)        , &
                  &              self%surveyGeometry_    %distanceMinimum (masses(iBin),iField)  &
                  &             )
             distanceMaximum=min(                                                                &
                  &              self%cosmologyFunctions_%distanceComoving(timeMinimum)        , &
                  &              self%surveyGeometry_    %distanceMaximum (masses(iBin),iField)  &
                  &             )
             outputWeight                           (iBin  ,iOutput)  &
                  & =outputWeight                   (iBin  ,iOutput)  &
                  & +self%surveyGeometry_%solidAngle(iField        )  &
                  & /3.0d0                                            &
                  & *max(                                             &
                  &      +0.0d0                                     , &
                  &      +distanceMaximum**3                          &
                  &      -distanceMinimum**3                          &
                  &    )
          end do
       end do
       where(outputWeight(iBin,:) < 0.0d0)
          outputWeight(iBin,:)=0.0d0
       end where
       if (any(outputWeight(iBin,:) > 0.0d0)) then
          outputWeight                  (iBin,:)  &
               &       =    outputWeight(iBin,:)  &
               &       /sum(outputWeight(iBin,:))
       else
          message=var_str("stellar mass function bin ")//iBin//" has zero weights"
          call Galacticus_Error_Report('massFunctionStellarConstructorInternal',message)
       end if
    end do
    ! Create a stellar mass property extractor.
    allocate(outputAnalysisPropertyExtractor_)
    outputAnalysisPropertyExtractor_     =outputAnalysisPropertyExtractorMassStellar()
    ! Create a log10 property operator.
    allocate(outputAnalysisPropertyOperator_)
    outputAnalysisPropertyOperator_      =outputAnalysisPropertyOperatorLog10       ()
    ! Create a bin width distribution normalizer.
    allocate(normalizerSequence)
    normalizer_ => normalizerSequence
    allocate(outputAnalysisDistributionNormalizerBinWidth   :: normalizer_%normalizer_)
    select type (normalizer_ => normalizer_%normalizer_)
    type is (outputAnalysisDistributionNormalizerBinWidth  )
       normalizer_=outputAnalysisDistributionNormalizerBinWidth  ()
    end select
    allocate(normalizer_%next)
    normalizer_ => normalizer_%next
    allocate(outputAnalysisDistributionNormalizerLog10ToLog :: normalizer_%normalizer_)
    select type (normalizer_ => normalizer_%normalizer_)
    type is (outputAnalysisDistributionNormalizerLog10ToLog)
       normalizer_=outputAnalysisDistributionNormalizerLog10ToLog()
    end select
    allocate(outputAnalysisDistributionNormalizer_)
    outputAnalysisDistributionNormalizer_=outputAnalysisDistributionNormalizerSequence(normalizerSequence)
    ! Construct the object. We convert masses to log10(masses) here.
    self%outputAnalysisVolumeFunction1D=                                         &
         & outputAnalysisVolumeFunction1D(                                       &
         &                                log10(masses)                        , &
         &                                outputWeight                         , &
         &                                outputAnalysisPropertyExtractor_     , &
         &                                outputAnalysisPropertyOperator_      , &
         &                                outputAnalysisDistributionOperator_  , &
         &                                outputAnalysisDistributionNormalizer_, &
         &                                galacticFilter_                        &
         &                               )
    ! Clean up.
    nullify(outputAnalysisPropertyExtractor_     )
    nullify(outputAnalysisPropertyOperator_      )
    nullify(outputAnalysisDistributionNormalizer_)
    nullify(normalizerSequence                   )
    return
  end function massFunctionStellarConstructorInternal
