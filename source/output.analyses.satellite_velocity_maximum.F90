!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

  !+    Contributions to this file made by: Xiaolong Du
  
  !!{
  Implements an output analysis class that computes satellite maximum circular velocity fraction.
  !!}

  use    :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use    :: Numerical_Interpolation , only : interpolator
  !$ use :: OMP_Lib                 , only : omp_lock_kind

  !![
  <outputAnalysis name="outputAnalysisSatelliteVelocityMaximum">
    <description>An output analysis class that computes satellite maximum circular velocity fraction as a function of time.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisSatelliteVelocityMaximum
     !!{
     An output analysis class that computes satellite maximum circular velocity fraction as a function of time.
     !!}
     private
     class           (darkMatterProfileDMOClass), pointer                   :: darkMatterProfileDMO_         => null(), darkMatterProfileDMOUnheated => null()
     class           (outputTimesClass         ), pointer                   :: outputTimes_                  => null()
     type            (varying_string           )                            :: fileName
     type            (interpolator             )                            :: interpolator_                          , interpolatorError_
     !$ integer      (omp_lock_kind            )                            :: accumulateLock
     double precision                                                       :: logLikelihood_
     double precision                                                       :: relativeModelUncertainty               , velocityMaximumInitial
     double precision                           , dimension(:), allocatable :: time                                   , fractionVelocityMaximum              , &
          &                                                                    fractionVelocityMaximumTarget          , varianceFractionVelocityMaximumTarget
   contains
     final     ::                  satelliteVelocityMaximumDestructor
     procedure :: analyze       => satelliteVelocityMaximumAnalyze
     procedure :: finalize      => satelliteVelocityMaximumFinalize
     procedure :: reduce        => satelliteVelocityMaximumReduce
     procedure :: logLikelihood => satelliteVelocityMaximumLogLikelihood
  end type outputAnalysisSatelliteVelocityMaximum

  interface outputAnalysisSatelliteVelocityMaximum
     !!{
     Constructors for the {\normalfont \ttfamily satelliteVelocityMaximum} output analysis class.
     !!}
     module procedure satelliteVelocityMaximumConstructorParameters
     module procedure satelliteVelocityMaximumConstructorInternal
  end interface outputAnalysisSatelliteVelocityMaximum

contains

  function satelliteVelocityMaximumConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily satelliteVelocityMaximum} output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    use :: Output_Times    , only : outputTimesClass
    implicit none
    type            (outputAnalysisSatelliteVelocityMaximum)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass             ), pointer       :: darkMatterProfileDMO_   , darkMatterProfileDMOUnheated
    class           (outputTimesClass                      ), pointer       :: outputTimes_
    type            (varying_string                        )                :: fileName
    double precision                                                        :: relativeModelUncertainty

    !![
    <inputParameter>
      <name>fileName</name>
      <source>parameters</source>
      <description>The name of the file from which to read the target dataset.</description>
    </inputParameter>
    <inputParameter>
      <name>relativeModelUncertainty</name>
      <defaultValue>0.0d0</defaultValue>
      <description>Relative model uncertainty.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"        source="parameters"                                               />
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMOUnheated" source="parameters"   parameterName="darkMatterProfileDMOUnheated"/>
    <objectBuilder class="outputTimes"           name="outputTimes_"                 source="parameters"                                               />
    !!]
    self=outputAnalysisSatelliteVelocityMaximum(fileName,relativeModelUncertainty,darkMatterProfileDMO_,darkMatterProfileDMOUnheated,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"       />
    <objectDestructor name="darkMatterProfileDMOUnheated"/>
    <objectDestructor name="outputTimes_"                />
    !!]
    return
  end function satelliteVelocityMaximumConstructorParameters
  
  function satelliteVelocityMaximumConstructorInternal(fileName,relativeModelUncertainty,darkMatterProfileDMO_,darkMatterProfileDMOUnheated,outputTimes_) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily satelliteVelocityMaximum} output analysis class for internal use.
    !!}
    use :: HDF5_Access            , only : hdf5Access
    use :: IO_HDF5                , only : hdf5Object
    use :: Numerical_Interpolation, only : GSL_Interp_CSpline
    use :: Table_Labels           , only : extrapolationTypeExtrapolate
    implicit none
    type            (outputAnalysisSatelliteVelocityMaximum)                              :: self
    type            (varying_string                        ), intent(in   )               :: fileName
    class           (darkMatterProfileDMOClass             ), intent(in   ), target       :: darkMatterProfileDMO_       , darkMatterProfileDMOUnheated
    class           (outputTimesClass                      ), intent(inout), target       :: outputTimes_
    double precision                                        , intent(in   )               :: relativeModelUncertainty
    type            (hdf5Object                            )                              :: file
    double precision                                        , allocatable  , dimension(:) :: time                        , velocityMaximum             , &
         &                                                                                   fractionVelocityMaximum     , velocityMaximumError        , &
         &                                                                                   fractionVelocityMaximumError
    integer         (c_size_t                              )                              :: i
    !![
    <constructorAssign variables="fileName, relativeModelUncertainty, *darkMatterProfileDMO_, *darkMatterProfileDMOUnheated, *outputTimes_"/>
    !!]

    ! Read properties from the file.
    !$ call hdf5Access%set()
    file=hdf5Object(char(fileName),readOnly=.true.)
    call file%readDataset('time'                ,time                )
    call file%readDataset('velocityMaximum'     ,velocityMaximum     )
    call file%readDataset('velocityMaximumError',velocityMaximumError)
    !$ call hdf5Access%unset()
    allocate(fractionVelocityMaximum     (size(velocityMaximum)))
    allocate(fractionVelocityMaximumError(size(velocityMaximum)))
    self%velocityMaximumInitial =velocityMaximum(1)
    fractionVelocityMaximum     =velocityMaximum     /self%velocityMaximumInitial
    fractionVelocityMaximumError=velocityMaximumError/self%velocityMaximumInitial
    ! Build interpolator.
    self%interpolator_     =interpolator(time,log(fractionVelocityMaximum     ),interpolationType=GSL_Interp_CSpline,extrapolationType=extrapolationTypeExtrapolate)
    self%interpolatorError_=interpolator(time,log(fractionVelocityMaximumError),interpolationType=GSL_Interp_CSpline,extrapolationType=extrapolationTypeExtrapolate)
    allocate(self%time                                 (outputTimes_%count()))
    allocate(self%fractionVelocityMaximum              (outputTimes_%count()))
    allocate(self%fractionVelocityMaximumTarget        (outputTimes_%count()))
    allocate(self%varianceFractionVelocityMaximumTarget(outputTimes_%count()))
    do i=1, outputTimes_%count()
       self%time                                 (i)=                                        outputTimes_%time(i)
       self%fractionVelocityMaximumTarget        (i)=exp(self%interpolator_     %interpolate(outputTimes_%time(i)))
       self%varianceFractionVelocityMaximumTarget(i)=exp(self%interpolatorError_%interpolate(outputTimes_%time(i)))**2
    end do
    !$ call OMP_Init_Lock(self%accumulateLock)
    self%fractionVelocityMaximum=0.0d0
    self%logLikelihood_         =0.0d0
    return
  end function satelliteVelocityMaximumConstructorInternal

  subroutine satelliteVelocityMaximumDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily satelliteVelocityMaximum} output analysis class.
    !!}
    !$ use :: OMP_Lib, only : OMP_Destroy_Lock
    implicit none
    type(outputAnalysisSatelliteVelocityMaximum), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"       />
    <objectDestructor name="self%darkMatterProfileDMOUnheated"/>
    <objectDestructor name="self%outputTimes_"                />
    !!]
    !$ call OMP_Destroy_Lock(self%accumulateLock)
    return
  end subroutine satelliteVelocityMaximumDestructor

  subroutine satelliteVelocityMaximumAnalyze(self,node,iOutput)
    !!{
    Analyze the maximum velocity tidal track.
    !!}
    use    :: Numerical_Constants_Math, only : Pi
    use    :: Mass_Distributions      , only : massDistributionClass
    !$ use :: OMP_Lib                 , only : OMP_Set_Lock         , OMP_Unset_Lock
    implicit none
    class           (outputAnalysisSatelliteVelocityMaximum), intent(inout) :: self
    type            (treeNode                              ), intent(inout) :: node
    integer         (c_size_t                              ), intent(in   ) :: iOutput
    class           (massDistributionClass                 ), pointer       :: massDistributionUnheated_, massDistribution_
    double precision                                                        :: fractionVelocityMaximum  , varianceFractionVelocityMaximum

    ! Skip non-satellites.
    if (.not.node%isSatellite()) return
    ! Extract the maximum circular velocity fraction.
    massDistribution_         => self%darkMatterProfileDMO_       %get(node)
    massDistributionUnheated_ => self%darkMatterProfileDMOUnheated%get(node)
    fractionVelocityMaximum=+massDistribution_        %velocityRotationCurveMaximum() &
         &                  /massDistributionUnheated_%velocityRotationCurveMaximum()
    !![
    <objectDestructor name="massDistribution_"        />
    <objectDestructor name="massDistributionUnheated_"/>
    !!]
    self%fractionVelocityMaximum(iOutput)=fractionVelocityMaximum
    ! Add model uncertainty.
    varianceFractionVelocityMaximum      =+  self%varianceFractionVelocityMaximumTarget(iOutput) &
         &                                +(                                                     &
         &                                   self%relativeModelUncertainty                       &
         &                                  *fractionVelocityMaximum                             &
         &                                 )**2
    self%logLikelihood_                  =+self%logLikelihood_                                   &
         &                                -0.5d0                                                 &
         &                                *(                                                     &
         &                                  +(                                                   &
         &                                    +fractionVelocityMaximum                           &
         &                                    -self%fractionVelocityMaximumTarget      (iOutput) &
         &                                   )**2                                                &
         &                                  /             varianceFractionVelocityMaximum        &
         &                                  +log(2.0d0*Pi*varianceFractionVelocityMaximum)       &
         &                                 )
    !$ call OMP_Unset_Lock(self%accumulateLock)
    return
  end subroutine satelliteVelocityMaximumAnalyze

  subroutine satelliteVelocityMaximumReduce(self,reduced)
    !!{
    Reduce over the maximum velocity tidal track output analysis.
    !!}
    use    :: Error  , only : Error_Report
    !$ use :: OMP_Lib, only : OMP_Set_Lock, OMP_Unset_Lock
    implicit none
    class(outputAnalysisSatelliteVelocityMaximum), intent(inout) :: self
    class(outputAnalysisClass                   ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisSatelliteVelocityMaximum)
       !$ call OMP_Set_Lock(reduced%accumulateLock)       
       reduced%fractionVelocityMaximum=reduced%fractionVelocityMaximum+self%fractionVelocityMaximum
       reduced%logLikelihood_         =reduced%logLikelihood_         +self%logLikelihood_
       !$ call OMP_Unset_Lock(reduced%accumulateLock)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine satelliteVelocityMaximumReduce

  subroutine satelliteVelocityMaximumFinalize(self,groupName)
    !!{
    Output results of the maximum velocity tidal track output analysis.
    !!}
#ifdef USEMPI
    use :: MPI_Utilities                   , only : mpiSelf
#endif
    use :: Output_HDF5                     , only : outputFile
    use :: HDF5_Access                     , only : hdf5Access
    use :: IO_HDF5                         , only : hdf5Object
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    class(outputAnalysisSatelliteVelocityMaximum), intent(inout)           :: self
    type (varying_string                        ), intent(in   ), optional :: groupName
    type (hdf5Object                            )               , target   :: analysesGroup, subGroup
    type (hdf5Object                            )               , pointer  :: inGroup
    type (hdf5Object                            )                          :: analysisGroup, dataset

#ifdef USEMPI
    ! If running under MPI, accumulate tracks across all processes.
    self%fractionVelocityMaximum=mpiSelf%sum(self%fractionVelocityMaximum)
    self%logLikelihood_         =mpiSelf%sum(self%logLikelihood_         )
#endif
    !$ call hdf5Access%set()
    analysesGroup =  outputFile   %openGroup('analyses'                         )
    inGroup       => analysesGroup
    if (present(groupName)) then
       subGroup   =  analysesGroup%openGroup(char(groupName)                    )
       inGroup    => subGroup
    end if
    analysisGroup=inGroup%openGroup('satelliteVelocityMaximum','Analysis of the satellite maximum circular velocity vs. time.')
    ! Write metadata describing this analysis.
    call analysisGroup%writeAttribute('Satellite maximum circular velocity fraction vs. time.','description'                        )
    call analysisGroup%writeAttribute("function1D"                                            ,'type'                               )
    call analysisGroup%writeAttribute('$t/\mathrm{Gyr}$'                                      ,'xAxisLabel'                         )
    call analysisGroup%writeAttribute('$V_\mathrm{max}/V_\mathrm{max,0}$'                     ,'yAxisLabel'                         )
    call analysisGroup%writeAttribute(.false.                                                 ,'xAxisIsLog'                         )
    call analysisGroup%writeAttribute(.true.                                                  ,'yAxisIsLog'                         )
    call analysisGroup%writeAttribute('time'                                                  ,'xDataset'                           )
    call analysisGroup%writeAttribute('fractionVelocityMaximum'                               ,'yDataset'                           )
    call analysisGroup%writeAttribute('fractionVelocityMaximumTarget'                         ,'yDatasetTarget'                     )
    call analysisGroup%writeAttribute(self%logLikelihood()                                    ,'logLikelihood'                      )
    call analysisGroup%writeAttribute('Simulations'                                           ,'targetLabel'                        )
    call analysisGroup%writeDataset  (self%time                         ,'time'                         ,'Time'                                          ,datasetReturned=dataset)
    call dataset      %writeAttribute('Gyr'                             ,'units'                                                                                                 )
    call dataset      %writeAttribute(gigaYear                          ,'unitsInSI'                                                                                             )
    call analysisGroup%writeDataset  (self%fractionVelocityMaximum      ,'fractionVelocityMaximum'      ,'Fraction of maximum circular velocity'         ,datasetReturned=dataset)
    call dataset      %writeAttribute(' '                               ,'units'                                                                                                 )
    call dataset      %writeAttribute(1.0d0                             ,'unitsInSI'                                                                                             )
    call analysisGroup%writeDataset  (self%fractionVelocityMaximumTarget,'fractionVelocityMaximumTarget','Fraction of maximum circular velocity [target]',datasetReturned=dataset)
    call dataset      %writeAttribute(' '                               ,'units'                                                                                                 )
    call dataset      %writeAttribute(1.0d0                             ,'unitsInSI'                                                                                             )
    !$ call hdf5Access%unset()
    return
  end subroutine satelliteVelocityMaximumFinalize
  
  double precision function satelliteVelocityMaximumLogLikelihood(self)
    !!{
    Return the log-likelihood of a maximum velocity tidal track output analysis.
    !!}
    implicit none
    class(outputAnalysisSatelliteVelocityMaximum), intent(inout) :: self

    satelliteVelocityMaximumLogLikelihood=self%logLikelihood_
    return
  end function satelliteVelocityMaximumLogLikelihood
