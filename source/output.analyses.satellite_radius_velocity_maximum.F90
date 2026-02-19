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
  Implements an output analysis class that computes satellite radius fraction at which the maximum circular velocity is reached.
  !!}

  use    :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use    :: Numerical_Interpolation , only : interpolator
  !$ use :: OMP_Lib                 , only : omp_lock_kind

  !![
  <outputAnalysis name="outputAnalysisSatelliteRadiusVelocityMaximum">
    <description>An output analysis class that computes satellite radius fraction at which the maximum circular velocity is reached as a function of time.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisSatelliteRadiusVelocityMaximum
     !!{
     An output analysis class that computes satellite radius fraction at which the maximum circular velocity is reached as a function of time.
     !!}
     private
     class           (darkMatterProfileDMOClass), pointer                   :: darkMatterProfileDMO_               => null(), darkMatterProfileDMOUnheated => null()
     class           (outputTimesClass         ), pointer                   :: outputTimes_                        => null()
     type            (varying_string           )                            :: fileName
     type            (interpolator             )                            :: interpolator_                                , interpolatorError_
     !$ integer      (omp_lock_kind            )                            :: accumulateLock
     double precision                                                       :: logLikelihood_
     double precision                                                       :: relativeModelUncertainty                     , radiusVelocityMaximumInitial
     double precision                           , dimension(:), allocatable :: time                                         , fractionRadiusVelocityMaximum              , &
          &                                                                    fractionRadiusVelocityMaximumTarget          , varianceFractionRadiusVelocityMaximumTarget
   contains
     final     ::                  satelliteRadiusVelocityMaximumDestructor
     procedure :: analyze       => satelliteRadiusVelocityMaximumAnalyze
     procedure :: finalize      => satelliteRadiusVelocityMaximumFinalize
     procedure :: reduce        => satelliteRadiusVelocityMaximumReduce
     procedure :: logLikelihood => satelliteRadiusVelocityMaximumLogLikelihood
  end type outputAnalysisSatelliteRadiusVelocityMaximum

  interface outputAnalysisSatelliteRadiusVelocityMaximum
     !!{
     Constructors for the \refClass{outputAnalysisSatelliteRadiusVelocityMaximum} output analysis class.
     !!}
     module procedure satelliteRadiusVelocityMaximumConstructorParameters
     module procedure satelliteRadiusVelocityMaximumConstructorInternal
  end interface outputAnalysisSatelliteRadiusVelocityMaximum

contains

  function satelliteRadiusVelocityMaximumConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisSatelliteRadiusVelocityMaximum} output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    use :: Output_Times    , only : outputTimesClass
    implicit none
    type            (outputAnalysisSatelliteRadiusVelocityMaximum)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass                   ), pointer       :: darkMatterProfileDMO_   , darkMatterProfileDMOUnheated
    class           (outputTimesClass                            ), pointer       :: outputTimes_
    type            (varying_string                              )                :: fileName
    double precision                                                              :: relativeModelUncertainty

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
    self=outputAnalysisSatelliteRadiusVelocityMaximum(fileName,relativeModelUncertainty,darkMatterProfileDMO_,darkMatterProfileDMOUnheated,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"       />
    <objectDestructor name="darkMatterProfileDMOUnheated"/>
    <objectDestructor name="outputTimes_"                />
    !!]
    return
  end function satelliteRadiusVelocityMaximumConstructorParameters
  
  function satelliteRadiusVelocityMaximumConstructorInternal(fileName,relativeModelUncertainty,darkMatterProfileDMO_,darkMatterProfileDMOUnheated,outputTimes_) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisSatelliteRadiusVelocityMaximum} output analysis class for internal use.
    !!}
    use :: HDF5_Access            , only : hdf5Access
    use :: IO_HDF5                , only : hdf5Object
    use :: Numerical_Interpolation, only : GSL_Interp_CSpline
    use :: Table_Labels           , only : extrapolationTypeExtrapolate
    implicit none
    type            (outputAnalysisSatelliteRadiusVelocityMaximum)                              :: self
    type            (varying_string                              ), intent(in   )               :: fileName
    class           (darkMatterProfileDMOClass                   ), intent(in   ), target       :: darkMatterProfileDMO_             , darkMatterProfileDMOUnheated
    class           (outputTimesClass                            ), intent(inout), target       :: outputTimes_
    double precision                                              , intent(in   )               :: relativeModelUncertainty
    type            (hdf5Object                                  )                              :: file
    double precision                                              , allocatable  , dimension(:) :: time                              , radiusVelocityMaximum       , &
         &                                                                                         fractionRadiusVelocityMaximum     , radiusVelocityMaximumError  , &
         &                                                                                         fractionRadiusVelocityMaximumError
    integer         (c_size_t                                    )                              :: i
    !![
    <constructorAssign variables="fileName, relativeModelUncertainty, *darkMatterProfileDMO_, *darkMatterProfileDMOUnheated, *outputTimes_"/>
    !!]

    ! Read properties from the file.
    !$ call hdf5Access%set()
    call file%openFile   (char(fileName)              ,readOnly=.true.                    )
    call file%readDataset('time'                      ,         time                      )
    call file%readDataset('radiusVelocityMaximum'     ,         radiusVelocityMaximum     )
    call file%readDataset('radiusVelocityMaximumError',         radiusVelocityMaximumError)
    call file%close      (                                                                )
    !$ call hdf5Access%unset()
    allocate(fractionRadiusVelocityMaximum     (size(radiusVelocityMaximum)))
    allocate(fractionRadiusVelocityMaximumError(size(radiusVelocityMaximum)))
    self%radiusVelocityMaximumInitial =radiusVelocityMaximum(1)
    fractionRadiusVelocityMaximum     =radiusVelocityMaximum     /self%radiusVelocityMaximumInitial
    fractionRadiusVelocityMaximumError=radiusVelocityMaximumError/self%radiusVelocityMaximumInitial
    ! Build interpolator.
    self%interpolator_     =interpolator(time,log(fractionRadiusVelocityMaximum     ),interpolationType=GSL_Interp_CSpline,extrapolationType=extrapolationTypeExtrapolate)
    self%interpolatorError_=interpolator(time,log(fractionRadiusVelocityMaximumError),interpolationType=GSL_Interp_CSpline,extrapolationType=extrapolationTypeExtrapolate)
    allocate(self%time                                       (outputTimes_%count()))
    allocate(self%fractionRadiusVelocityMaximum              (outputTimes_%count()))
    allocate(self%fractionRadiusVelocityMaximumTarget        (outputTimes_%count()))
    allocate(self%varianceFractionRadiusVelocityMaximumTarget(outputTimes_%count()))
    do i=1, outputTimes_%count()
       self%time                                       (i)=                                        outputTimes_%time(i)
       self%fractionRadiusVelocityMaximumTarget        (i)=exp(self%interpolator_     %interpolate(outputTimes_%time(i)))
       self%varianceFractionRadiusVelocityMaximumTarget(i)=exp(self%interpolatorError_%interpolate(outputTimes_%time(i)))**2
    end do
    !$ call OMP_Init_Lock(self%accumulateLock)
    self%fractionRadiusVelocityMaximum=0.0d0
    self%logLikelihood_               =0.0d0
    return
  end function satelliteRadiusVelocityMaximumConstructorInternal

  subroutine satelliteRadiusVelocityMaximumDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisSatelliteRadiusVelocityMaximum} output analysis class.
    !!}
    !$ use :: OMP_Lib, only : OMP_Destroy_Lock
    implicit none
    type(outputAnalysisSatelliteRadiusVelocityMaximum), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"       />
    <objectDestructor name="self%darkMatterProfileDMOUnheated"/>
    <objectDestructor name="self%outputTimes_"                />
    !!]
    !$ call OMP_Destroy_Lock(self%accumulateLock)
    return
  end subroutine satelliteRadiusVelocityMaximumDestructor

  subroutine satelliteRadiusVelocityMaximumAnalyze(self,node,iOutput)
    !!{
    Analyze the maximum velocity tidal track.
    !!}
    use    :: Numerical_Constants_Math, only : Pi
    use    :: Mass_Distributions      , only : massDistributionClass
    !$ use :: OMP_Lib                 , only : OMP_Set_Lock         , OMP_Unset_Lock
    implicit none
    class           (outputAnalysisSatelliteRadiusVelocityMaximum), intent(inout) :: self
    type            (treeNode                                    ), intent(inout) :: node
    integer         (c_size_t                                    ), intent(in   ) :: iOutput
    class           (massDistributionClass                       ), pointer       :: massDistributionUnheated_    , massDistribution_
    double precision                                                              :: fractionRadiusVelocityMaximum, varianceFractionRadiusVelocityMaximum

    ! Skip non-satellites.
    if (.not.node%isSatellite()) return
    ! Extract the maximum circular velocity fraction.
    massDistribution_         => self%darkMatterProfileDMO_       %get(node)
    massDistributionUnheated_ => self%darkMatterProfileDMOUnheated%get(node)
    fractionRadiusVelocityMaximum=+massDistribution_        %radiusRotationCurveMaximum() &
         &                        /massDistributionUnheated_%radiusRotationCurveMaximum()
    !![
    <objectDestructor name="massDistribution_"        />
    <objectDestructor name="massDistributionUnheated_"/>
    !!]
    !$ call OMP_Set_Lock(self%accumulateLock)
    self%fractionRadiusVelocityMaximum(iOutput)=fractionRadiusVelocityMaximum
    ! Add model uncertainty.
    varianceFractionRadiusVelocityMaximum      =+self%varianceFractionRadiusVelocityMaximumTarget(iOutput) &
         &                                      +(                                                         &
         &                                         self%relativeModelUncertainty                           &
         &                                        *fractionRadiusVelocityMaximum                           &
         &                                       )**2
    self%logLikelihood_                        =+self%logLikelihood_                                       &
         &                                      -0.5d0                                                     &
         &                                      *(                                                         &
         &                                        +(                                                       &
         &                                          +fractionRadiusVelocityMaximum                         &
         &                                          -self%fractionRadiusVelocityMaximumTarget    (iOutput) &
         &                                         )**2                                                    &
         &                                        /             varianceFractionRadiusVelocityMaximum      &
         &                                        +log(2.0d0*Pi*varianceFractionRadiusVelocityMaximum)     &
         &                                       )
    !$ call OMP_Unset_Lock(self%accumulateLock)
    return
  end subroutine satelliteRadiusVelocityMaximumAnalyze

  subroutine satelliteRadiusVelocityMaximumReduce(self,reduced)
    !!{
    Reduce over the maximum velocity tidal track output analysis.
    !!}
    use    :: Error  , only : Error_Report
    !$ use :: OMP_Lib, only : OMP_Set_Lock, OMP_Unset_Lock
    implicit none
    class(outputAnalysisSatelliteRadiusVelocityMaximum), intent(inout) :: self
    class(outputAnalysisClass                         ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisSatelliteRadiusVelocityMaximum)
       !$ call OMP_Set_Lock(reduced%accumulateLock)       
       reduced%fractionRadiusVelocityMaximum=reduced%fractionRadiusVelocityMaximum+self%fractionRadiusVelocityMaximum
       reduced%logLikelihood_               =reduced%logLikelihood_               +self%logLikelihood_
       !$ call OMP_Unset_Lock(reduced%accumulateLock)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine satelliteRadiusVelocityMaximumReduce

  subroutine satelliteRadiusVelocityMaximumFinalize(self,groupName)
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
    class(outputAnalysisSatelliteRadiusVelocityMaximum), intent(inout)           :: self
    type (varying_string                              ), intent(in   ), optional :: groupName
    type (hdf5Object                                  )               , target   :: analysesGroup, subGroup
    type (hdf5Object                                  )               , pointer  :: inGroup
    type (hdf5Object                                  )                          :: analysisGroup, dataset

#ifdef USEMPI
    ! If running under MPI, accumulate tracks across all processes.
    self%fractionRadiusVelocityMaximum=mpiSelf%sum(self%fractionRadiusVelocityMaximum)
    self%logLikelihood_               =mpiSelf%sum(self%logLikelihood_               )
#endif
    !$ call hdf5Access%set()
    analysesGroup =  outputFile   %openGroup('analyses'                         )
    inGroup       => analysesGroup
    if (present(groupName)) then
       subGroup   =  analysesGroup%openGroup(char(groupName)                    )
       inGroup    => subGroup
    end if
    analysisGroup=inGroup%openGroup('satelliteRadiusVelocityMaximum','Analysis of the satellite maximum circular velocity radius vs. time.')
    ! Write metadata describing this analysis.
    call analysisGroup%writeAttribute('Satellite maximum circular velocity radius fraction vs. time.','description'                              )
    call analysisGroup%writeAttribute("function1D"                                                   ,'type'                                     )
    call analysisGroup%writeAttribute('$t/\mathrm{Gyr}$'                                             ,'xAxisLabel'                               )
    call analysisGroup%writeAttribute('$R_\mathrm{max}/R_\mathrm{max,0}$'                            ,'yAxisLabel'                               )
    call analysisGroup%writeAttribute(.false.                                                        ,'xAxisIsLog'                               )
    call analysisGroup%writeAttribute(.true.                                                         ,'yAxisIsLog'                               )
    call analysisGroup%writeAttribute('time'                                                         ,'xDataset'                                 )
    call analysisGroup%writeAttribute('fractionRadiusVelocityMaximum'                                ,'yDataset'                                 )
    call analysisGroup%writeAttribute('fractionRadiusVelocityMaximumTarget'                          ,'yDatasetTarget'                           )
    call analysisGroup%writeAttribute(self%logLikelihood()                                           ,'logLikelihood'                            )
    call analysisGroup%writeAttribute('Simulations'                                                  ,'targetLabel'                              )
    call analysisGroup%writeDataset  (self%time                               ,'time'                               ,'Time'                                          ,datasetReturned=dataset)
    call dataset      %writeAttribute('Gyr'                                   ,'units'                                                                                                       )
    call dataset      %writeAttribute(gigaYear                                ,'unitsInSI'                                                                                                   )
    call dataset      %close         (                                                                                                                                                       )
    call analysisGroup%writeDataset  (self%fractionRadiusVelocityMaximum      ,'fractionRadiusVelocityMaximum'      ,'Fraction of maximum circular velocity'         ,datasetReturned=dataset)
    call dataset      %writeAttribute(' '                                     ,'units'                                                                                                       )
    call dataset      %writeAttribute(1.0d0                                   ,'unitsInSI'                                                                                                   )
    call dataset      %close         (                                                                                                                                                       )
    call analysisGroup%writeDataset  (self%fractionRadiusVelocityMaximumTarget,'fractionRadiusVelocityMaximumTarget','Fraction of maximum circular velocity [target]',datasetReturned=dataset)
    call dataset      %writeAttribute(' '                                     ,'units'                                                                                                       )
    call dataset      %writeAttribute(1.0d0                                   ,'unitsInSI'                                                                                                   )
    call dataset      %close         (                                                                                                                                                       )
    call analysisGroup%close         (                                                                                                                                                       )
    if (present(groupName)) &
         & call subGroup%close       (                                                                                                                                                       )
    call analysesGroup%close         (                                                                                                                                                       )
    !$ call hdf5Access%unset()
    return
  end subroutine satelliteRadiusVelocityMaximumFinalize
  
  double precision function satelliteRadiusVelocityMaximumLogLikelihood(self)
    !!{
    Return the log-likelihood of a maximum velocity tidal track output analysis.
    !!}
    implicit none
    class(outputAnalysisSatelliteRadiusVelocityMaximum), intent(inout) :: self

    satelliteRadiusVelocityMaximumLogLikelihood=self%logLikelihood_
    return
  end function satelliteRadiusVelocityMaximumLogLikelihood
