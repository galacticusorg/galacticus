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
  Implements an output analysis class that computes satellite bound mass fraction.
  !!}

  use    :: Numerical_Interpolation, only : interpolator
  !$ use :: OMP_Lib, only : omp_lock_kind

  !![
  <outputAnalysis name="outputAnalysisSatelliteBoundMass">
    <description>An output analysis class that computes satellite bound mass fraction as a function of time.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisSatelliteBoundMass
     !!{
     An output analysis class that computes satellite bound mass fraction as a function of time.
     !!}
     private
     class           (outputTimesClass), pointer                   :: outputTimes_             => null()
     type            (varying_string  )                            :: fileName
     type            (interpolator    )                            :: interpolator_                     , interpolatorError_
     !$ integer      (omp_lock_kind   )                            :: accumulateLock
     double precision                                              :: logLikelihood_
     double precision                                              :: relativeModelUncertainty          , boundMassInitial
     double precision                  , dimension(:), allocatable :: time                              , fractionBoundMass              , &
          &                                                           fractionBoundMassTarget           , varianceFractionBoundMassTarget
   contains
     final     ::                  satelliteBoundMassDestructor
     procedure :: analyze       => satelliteBoundMassAnalyze
     procedure :: finalize      => satelliteBoundMassFinalize
     procedure :: reduce        => satelliteBoundMassReduce
     procedure :: logLikelihood => satelliteBoundMassLogLikelihood
  end type outputAnalysisSatelliteBoundMass

  interface outputAnalysisSatelliteBoundMass
     !!{
     Constructors for the \refClass{outputAnalysisSatelliteBoundMass} output analysis class.
     !!}
     module procedure satelliteBoundMassConstructorParameters
     module procedure satelliteBoundMassConstructorInternal
  end interface outputAnalysisSatelliteBoundMass

contains

  function satelliteBoundMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisSatelliteBoundMass} output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    use :: Output_Times    , only : outputTimesClass
    implicit none
    type            (outputAnalysisSatelliteBoundMass)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (outputTimesClass                ), pointer       :: outputTimes_
    type            (varying_string                  )                :: fileName
    double precision                                                  :: relativeModelUncertainty

    !![
    <inputParameter>
      <name>fileName</name>
      <description>The name of the file from which to read the target dataset.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>relativeModelUncertainty</name>
      <defaultValue>0.0d0</defaultValue>
      <description>Relative model uncertainty.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    !!]
    self=outputAnalysisSatelliteBoundMass(fileName,relativeModelUncertainty,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"/>
    !!]
    return
  end function satelliteBoundMassConstructorParameters
  
  function satelliteBoundMassConstructorInternal(fileName,relativeModelUncertainty,outputTimes_) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisSatelliteBoundMass} output analysis class for internal use.
    !!}
    use :: HDF5_Access            , only : hdf5Access
    use :: IO_HDF5                , only : hdf5Object
    use :: Numerical_Interpolation, only : GSL_Interp_CSpline
    use :: Table_Labels           , only : extrapolationTypeExtrapolate
    implicit none
    type            (outputAnalysisSatelliteBoundMass)                              :: self
    type            (varying_string                  ), intent(in   )               :: fileName
    class           (outputTimesClass                ), intent(inout), target       :: outputTimes_
    double precision                                  , intent(in   )               :: relativeModelUncertainty
    type            (hdf5Object                      )                              :: file
    double precision                                  , allocatable  , dimension(:) :: time                    , boundMass     , &
         &                                                                             fractionBoundMass       , boundMassError, &
         &                                                                             fractionBoundMassError
    integer         (c_size_t                        )                              :: i
    !![
    <constructorAssign variables="fileName, relativeModelUncertainty, *outputTimes_"/>
    !!]

    ! Read properties from the file.
    !$ call hdf5Access%set()
    file=hdf5Object(char(fileName),readOnly=.true.)
    call file%readDataset('time'          ,time          )
    call file%readDataset('boundMass'     ,boundMass     )
    call file%readDataset('boundMassError',boundMassError)
    !$ call hdf5Access%unset()
    allocate(fractionBoundMass     (size(boundMass)))
    allocate(fractionBoundMassError(size(boundMass)))
    self%boundMassInitial =boundMass(1)
    fractionBoundMass     =boundMass     /self%boundMassInitial
    fractionBoundMassError=boundMassError/self%boundMassInitial
    ! Build interpolator.
    self%interpolator_     =interpolator(time,log(fractionBoundMass     ),interpolationType=GSL_Interp_CSpline,extrapolationType=extrapolationTypeExtrapolate)
    self%interpolatorError_=interpolator(time,log(fractionBoundMassError),interpolationType=GSL_Interp_CSpline,extrapolationType=extrapolationTypeExtrapolate)
    allocate(self%time                           (outputTimes_%count()))
    allocate(self%fractionBoundMass              (outputTimes_%count()))
    allocate(self%fractionBoundMassTarget        (outputTimes_%count()))
    allocate(self%varianceFractionBoundMassTarget(outputTimes_%count()))
    do i=1, outputTimes_%count()
       self%time                           (i)=                                        outputTimes_%time(i)
       self%fractionBoundMassTarget        (i)=exp(self%interpolator_     %interpolate(outputTimes_%time(i)))
       self%varianceFractionBoundMassTarget(i)=exp(self%interpolatorError_%interpolate(outputTimes_%time(i)))**2
    end do
    !$ call OMP_Init_Lock(self%accumulateLock)
    self%fractionBoundMass=0.0d0
    self%logLikelihood_   =0.0d0
    return
  end function satelliteBoundMassConstructorInternal

  subroutine satelliteBoundMassDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisSatelliteBoundMass} output analysis class.
    !!}
    !$ use :: OMP_Lib, only : OMP_Destroy_Lock
    implicit none
    type(outputAnalysisSatelliteBoundMass), intent(inout) :: self

    !![
    <objectDestructor name="self%outputTimes_"/>
    !!]
    !$ call OMP_Destroy_Lock(self%accumulateLock)
    return
  end subroutine satelliteBoundMassDestructor

  subroutine satelliteBoundMassAnalyze(self,node,iOutput)
    !!{
    Analyze the maximum velocity tidal track.
    !!}
    use    :: Galacticus_Nodes        , only : nodeComponentSatellite
    use    :: Numerical_Constants_Math, only : Pi
    !$ use :: OMP_Lib                 , only : OMP_Set_Lock          , OMP_Unset_Lock
    implicit none
    class           (outputAnalysisSatelliteBoundMass), intent(inout) :: self
    type            (treeNode                        ), intent(inout) :: node
    integer         (c_size_t                        ), intent(in   ) :: iOutput
    class           (nodeComponentSatellite          ), pointer       :: satellite
    double precision                                                  :: fractionMassBound, varianceFractionBoundMass

    ! Skip non-satellites.
    if (.not.node%isSatellite()) return
    ! Extract the bound mass fraction.
    satellite         =>  node     %satellite       ()
    fractionMassBound =  +satellite%boundMass       () &
         &               /self     %boundMassInitial
    !$ call OMP_Set_Lock(self%accumulateLock)
    self%fractionBoundMass(iOutput)=fractionMassBound
    ! Add model uncertainty.
    varianceFractionBoundMass      =+  self%varianceFractionBoundMassTarget(iOutput) &
         &                          +(                                               &
         &                             self%relativeModelUncertainty                 &
         &                            *fractionMassBound                             &
         &                           )**2
    self%logLikelihood_            =+self%logLikelihood_                             &
         &                          -0.5d0                                           &
         &                          *(                                               &
         &                            +(                                             &
         &                              +fractionMassBound                           &
         &                              -self%fractionBoundMassTarget      (iOutput) &
         &                             )**2                                          &
         &                            /             varianceFractionBoundMass        &
         &                            +log(2.0d0*Pi*varianceFractionBoundMass)       &
         &                           )
    !$ call OMP_Unset_Lock(self%accumulateLock)
    return
  end subroutine satelliteBoundMassAnalyze

  subroutine satelliteBoundMassReduce(self,reduced)
    !!{
    Reduce over the maximum velocity tidal track output analysis.
    !!}
    use    :: Error  , only : Error_Report
    !$ use :: OMP_Lib, only : OMP_Set_Lock, OMP_Unset_Lock
    implicit none
    class(outputAnalysisSatelliteBoundMass), intent(inout) :: self
    class(outputAnalysisClass             ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisSatelliteBoundMass)
       !$ call OMP_Set_Lock(reduced%accumulateLock)       
       reduced%fractionBoundMass=reduced%fractionBoundMass+self%fractionBoundMass
       reduced%logLikelihood_   =reduced%logLikelihood_   +self%logLikelihood_
       !$ call OMP_Unset_Lock(reduced%accumulateLock)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine satelliteBoundMassReduce

  subroutine satelliteBoundMassFinalize(self,groupName)
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
    class(outputAnalysisSatelliteBoundMass), intent(inout)           :: self
    type (varying_string                  ), intent(in   ), optional :: groupName
    type (hdf5Object                      )               , target   :: analysesGroup, subGroup
    type (hdf5Object                      )               , pointer  :: inGroup
    type (hdf5Object                      )                          :: analysisGroup, dataset

#ifdef USEMPI
    ! If running under MPI, accumulate tracks across all processes.
    self%fractionBoundMass=mpiSelf%sum(self%fractionBoundMass)
    self%logLikelihood_   =mpiSelf%sum(self%logLikelihood_   )
#endif
    !$ call hdf5Access%set()
    analysesGroup =  outputFile   %openGroup('analyses'     )
    inGroup       => analysesGroup
    if (present(groupName)) then
       subGroup   =  analysesGroup%openGroup(char(groupName))
       inGroup    => subGroup
    end if
    analysisGroup=inGroup%openGroup('satelliteBoundMass','Analysis of the satellite bound mass vs. time.')
    ! Write metadata describing this analysis.
    call analysisGroup%writeAttribute('Satellite bound mass fraction vs. time.','description'                  )
    call analysisGroup%writeAttribute("function1D"                             ,'type'                         )
    call analysisGroup%writeAttribute('$t/\mathrm{Gyr}$'                       ,'xAxisLabel'                   )
    call analysisGroup%writeAttribute('$M_\mathrm{bound}/M_\mathrm{bound,0}$'  ,'yAxisLabel'                   )
    call analysisGroup%writeAttribute(.false.                                  ,'xAxisIsLog'                   )
    call analysisGroup%writeAttribute(.true.                                   ,'yAxisIsLog'                   )
    call analysisGroup%writeAttribute('time'                                   ,'xDataset'                     )
    call analysisGroup%writeAttribute('fractionBoundMass'                      ,'yDataset'                     )
    call analysisGroup%writeAttribute('fractionBoundMassTarget'                ,'yDatasetTarget'               )
    call analysisGroup%writeAttribute(self%logLikelihood()                     ,'logLikelihood'                )
    call analysisGroup%writeAttribute('Simulations'                            ,'targetLabel'                  )
    call analysisGroup%writeDataset  (self%time                   ,'time'                   ,'Time'                           ,datasetReturned=dataset)
    call dataset      %writeAttribute('Gyr'                       ,'units'                                                                            )
    call dataset      %writeAttribute(gigaYear                    ,'unitsInSI'                                                                        )
    call analysisGroup%writeDataset  (self%fractionBoundMass      ,'fractionBoundMass'      ,'Fraction of bound mass'         ,datasetReturned=dataset)
    call dataset      %writeAttribute(' '                         ,'units'                                                                            )
    call dataset      %writeAttribute(1.0d0                       ,'unitsInSI'                                                                        )
    call analysisGroup%writeDataset  (self%fractionBoundMassTarget,'fractionBoundMassTarget','Fraction of bound mass [target]',datasetReturned=dataset)
    call dataset      %writeAttribute(' '                         ,'units'                                                                            )
    call dataset      %writeAttribute(1.0d0                       ,'unitsInSI'                                                                        )
    !$ call hdf5Access%unset()
    return
  end subroutine satelliteBoundMassFinalize
  
  double precision function satelliteBoundMassLogLikelihood(self)
    !!{
    Return the log-likelihood of a maximum velocity tidal track output analysis.
    !!}
    implicit none
    class(outputAnalysisSatelliteBoundMass), intent(inout) :: self

    satelliteBoundMassLogLikelihood=self%logLikelihood_
    return
  end function satelliteBoundMassLogLikelihood
