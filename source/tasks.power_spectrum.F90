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

  use :: Cosmological_Density_Field     , only : cosmologicalMassVarianceClass
  use :: Cosmology_Functions            , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters           , only : cosmologyParametersClass
  use :: Linear_Growth                  , only : linearGrowthClass
  use :: Output_Times                   , only : outputTimesClass
  use :: Power_Spectra                  , only : powerSpectrumClass
  use :: Power_Spectra_Nonlinear        , only : powerSpectrumNonlinearClass
  use :: Power_Spectrum_Window_Functions, only : powerSpectrumWindowFunctionClass
  use :: Transfer_Functions             , only : transferFunctionClass

  !![
  <task name="taskPowerSpectra">
   <description>A task which computes and outputs the power spectrum and related quantities.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskPowerSpectra
     !!{
     Implementation of a task which computes and outputs the power spectrum and related quantities.
     !!}
     private
     class           (cosmologyParametersClass         ), pointer :: cosmologyParameters_         => null()
     class           (cosmologyFunctionsClass          ), pointer :: cosmologyFunctions_          => null()
     class           (linearGrowthClass                ), pointer :: linearGrowth_                => null()
     class           (transferFunctionClass            ), pointer :: transferFunction_            => null()
     class           (powerSpectrumClass               ), pointer :: powerSpectrum_               => null()
     class           (powerSpectrumNonlinearClass      ), pointer :: powerSpectrumNonlinear_      => null()
     class           (powerSpectrumWindowFunctionClass ), pointer :: powerSpectrumWindowFunction_ => null()
     class           (cosmologicalMassVarianceClass    ), pointer :: cosmologicalMassVariance_    => null()
     class           (outputTimesClass                 ), pointer :: outputTimes_                 => null()
     double precision                                             :: wavenumberMinimum                     , wavenumberMaximum, &
          &                                                          pointsPerDecade
     logical                                                      :: includeNonLinear
     type            (varying_string                   )          :: outputGroup
   contains
     final     ::            powerSpectraDestructor
     procedure :: perform => powerSpectraPerform
  end type taskPowerSpectra

  interface taskPowerSpectra
     !!{
     Constructors for the \refClass{taskPowerSpectra} task.
     !!}
     module procedure powerSpectraConstructorParameters
     module procedure powerSpectraConstructorInternal
  end interface taskPowerSpectra

contains

  function powerSpectraConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskPowerSpectra} task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (taskPowerSpectra                )                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (cosmologyParametersClass        ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    class           (linearGrowthClass               ), pointer       :: linearGrowth_
    class           (transferFunctionClass           ), pointer       :: transferFunction_
    class           (powerSpectrumClass              ), pointer       :: powerSpectrum_
    class           (powerSpectrumNonlinearClass     ), pointer       :: powerSpectrumNonlinear_
    class           (powerSpectrumWindowFunctionClass), pointer       :: powerSpectrumWindowFunction_
    class           (cosmologicalMassVarianceClass   ), pointer       :: cosmologicalMassVariance_
    class           (outputTimesClass                ), pointer       :: outputTimes_
    double precision                                                  :: wavenumberMinimum           , wavenumberMaximum, &
         &                                                               pointsPerDecade
    logical                                                           :: includeNonLinear
    type            (varying_string                  )                :: outputGroup

    !![
    <inputParameter>
      <name>wavenumberMinimum</name>
      <defaultValue>1.0d-3</defaultValue>
      <description>The minimum wavenumber at which to tabulate power spectra.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>wavenumberMaximum</name>
      <defaultValue>1.0d+3</defaultValue>
      <description>The maximum wavenumber at which to tabulate power spectra.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>pointsPerDecade</name>
      <defaultValue>10.0d0</defaultValue>
      <description>The number of points per decade of wavenumber at which to tabulate power spectra.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>includeNonLinear</name>
      <defaultValue>.false.</defaultValue>
      <description>If true the nonlinear power spectrum is also computed and output.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>outputGroup</name>
      <defaultValue>var_str('.')</defaultValue>
      <description>The HDF5 output group within which to write power spectrum data.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"         name="cosmologyParameters_"         source="parameters"/>
    <objectBuilder class="cosmologyFunctions"          name="cosmologyFunctions_"          source="parameters"/>
    <objectBuilder class="linearGrowth"                name="linearGrowth_"                source="parameters"/>
    <objectBuilder class="transferFunction"            name="transferFunction_"            source="parameters"/>
    <objectBuilder class="powerSpectrum"               name="powerSpectrum_"               source="parameters"/>
    <objectBuilder class="powerSpectrumWindowFunction" name="powerSpectrumWindowFunction_" source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"    name="cosmologicalMassVariance_"    source="parameters"/>
    <objectBuilder class="outputTimes"                 name="outputTimes_"                 source="parameters"/>
    !!]
    if (includeNonLinear) then
       !![
       <objectBuilder class="powerSpectrumNonlinear"   name="powerSpectrumNonlinear_"      source="parameters"/>
       !!]
    else
       powerSpectrumNonlinear_ => null()
    end if
    self=taskPowerSpectra(                               &
         &                 wavenumberMinimum           , &
         &                 wavenumberMaximum           , &
         &                 pointsPerDecade             , &
         &                 includeNonLinear            , &
         &                 outputGroup                 , &
         &                 cosmologyParameters_        , &
         &                 cosmologyFunctions_         , &
         &                 linearGrowth_               , &
         &                 transferFunction_           , &
         &                 powerSpectrum_              , &
         &                 powerSpectrumNonlinear_     , &
         &                 powerSpectrumWindowFunction_, &
         &                 cosmologicalMassVariance_   , &
         &                 outputTimes_                  &
         &                )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"        />
    <objectDestructor name="cosmologyFunctions_"         />
    <objectDestructor name="linearGrowth_"               />
    <objectDestructor name="transferFunction_"           />
    <objectDestructor name="powerSpectrum_"              />
    <objectDestructor name="powerSpectrumNonlinear_"     />
    <objectDestructor name="powerSpectrumWindowFunction_"/>
    <objectDestructor name="cosmologicalMassVariance_"   />
    <objectDestructor name="outputTimes_"                />
    !!]
    return
  end function powerSpectraConstructorParameters

  function powerSpectraConstructorInternal(                               &
       &                                    wavenumberMinimum           , &
       &                                    wavenumberMaximum           , &
       &                                    pointsPerDecade             , &
       &                                    includeNonLinear            , &
       &                                    outputGroup                 , &
       &                                    cosmologyParameters_        , &
       &                                    cosmologyFunctions_         , &
       &                                    linearGrowth_               , &
       &                                    transferFunction_           , &
       &                                    powerSpectrum_              , &
       &                                    powerSpectrumNonlinear_     , &
       &                                    powerSpectrumWindowFunction_, &
       &                                    cosmologicalMassVariance_   , &
       &                                    outputTimes_                  &
       &                                   ) result(self)
    !!{
    Internal constructor for the \refClass{taskPowerSpectra} task class.
    !!}
    implicit none
    type            (taskPowerSpectra                )                        :: self
    class           (cosmologyParametersClass        ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass         ), intent(in   ), target :: cosmologyFunctions_
    class           (linearGrowthClass               ), intent(in   ), target :: linearGrowth_
    class           (transferFunctionClass           ), intent(in   ), target :: transferFunction_
    class           (powerSpectrumClass              ), intent(in   ), target :: powerSpectrum_
    class           (powerSpectrumNonlinearClass     ), intent(in   ), target :: powerSpectrumNonlinear_
    class           (powerSpectrumWindowFunctionClass), intent(in   ), target :: powerSpectrumWindowFunction_
    class           (cosmologicalMassVarianceClass   ), intent(in   ), target :: cosmologicalMassVariance_
    class           (outputTimesClass                ), intent(in   ), target :: outputTimes_
    double precision                                  , intent(in   )         :: wavenumberMinimum           , wavenumberMaximum, &
         &                                                                       pointsPerDecade
    logical                                           , intent(in   )         :: includeNonLinear
    type            (varying_string                  ), intent(in   )         :: outputGroup
    !![
    <constructorAssign variables="wavenumberMinimum, wavenumberMaximum, pointsPerDecade, includeNonLinear, outputGroup,*cosmologyParameters_,*cosmologyFunctions_,*linearGrowth_,*transferFunction_,*powerSpectrum_,*powerSpectrumNonlinear_,*powerSpectrumWindowFunction_,*cosmologicalMassVariance_, *outputTimes_"/>
    !!]

    return
  end function powerSpectraConstructorInternal

  subroutine powerSpectraDestructor(self)
    !!{
    Destructor for the \refClass{taskPowerSpectra} task class.
    !!}
    implicit none
    type(taskPowerSpectra), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"        />
    <objectDestructor name="self%cosmologyFunctions_"         />
    <objectDestructor name="self%linearGrowth_"               />
    <objectDestructor name="self%transferFunction_"           />
    <objectDestructor name="self%powerSpectrum_"              />
    <objectDestructor name="self%powerSpectrumNonlinear_"     />
    <objectDestructor name="self%powerSpectrumWindowFunction_"/>
    <objectDestructor name="self%cosmologicalMassVariance_"   />
    <objectDestructor name="self%outputTimes_"                />
    !!]
    return
  end subroutine powerSpectraDestructor

  subroutine powerSpectraPerform(self,status)
    !!{
    Compute and output the halo mass function.
    !!}
    use            :: Display                         , only : displayIndent     , displayUnindent
    use            :: Error                           , only : errorStatusSuccess
    use            :: Output_HDF5                     , only : outputFile
    use            :: IO_HDF5                         , only : hdf5Object
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Numerical_Constants_Astronomical, only : massSolar         , megaParsec
    use            :: Numerical_Constants_Math        , only : Pi
    use            :: Numerical_Integration           , only : GSL_Integ_Gauss15 , integrator
    use            :: Numerical_Ranges                , only : Make_Range        , rangeTypeLogarithmic
    use            :: String_Handling                 , only : operator(//)
    implicit none
    class           (taskPowerSpectra          ), intent(inout), target         :: self
    integer                                     , intent(  out), optional       :: status
    integer         (c_size_t                  )                                :: outputCount              , wavenumberCount    , &
         &                                                                         iOutput                  , iWavenumber
    double precision                            , allocatable  , dimension(:  ) :: wavenumber               , massScale          , &
         &                                                                         epochTime                , epochRedshift
    double precision                            , allocatable  , dimension(:,:) :: powerSpectrumNonLinear   , sigmaNonLinear     , &
         &                                                                         sigma                    , sigmaGradient      , &
         &                                                                         powerSpectrumLinear      , growthFactor       , &
         &                                                                         growthFactorLogDerivative, transferFunction
    double precision                                                            :: wavenumberMinimum        , wavenumberMaximum
    type            (integrator                )                                :: integrator_
    type            (hdf5Object                )                                :: outputsGroup             , outputGroup        , &
         &                                                                         containerGroup           , dataset
    type            (varying_string            )                                :: groupName                , description

    call displayIndent('Begin task: power spectrum')
    ! Get the requested output redshifts.
    outputCount      =self%outputTimes_%count()
    ! Compute number of tabulation points.
    wavenumberCount=int(log10(self%wavenumberMaximum/self%wavenumberMinimum)*self%pointsPerDecade)+1
    ! Allocate arrays for power spectra.
    allocate(wavenumber               (wavenumberCount            ))
    allocate(powerSpectrumLinear      (wavenumberCount,outputCount))
    allocate(transferFunction         (wavenumberCount,outputCount))
    allocate(massScale                (wavenumberCount            ))
    allocate(sigma                    (wavenumberCount,outputCount))
    allocate(sigmaGradient            (wavenumberCount,outputCount))
    allocate(growthFactor             (wavenumberCount,outputCount))
    allocate(growthFactorLogDerivative(wavenumberCount,outputCount))
    allocate(epochTime                (outputCount))
    allocate(epochRedshift            (outputCount))
    if (self%includeNonLinear) then
       allocate(powerSpectrumNonLinear   (wavenumberCount,outputCount))
       allocate(sigmaNonLinear           (wavenumberCount,outputCount))
    else
       allocate(powerSpectrumNonLinear   (              0,          0))
       allocate(sigmaNonLinear           (              0,          0))
    end if
    ! Build a range of wavenumbers.
    wavenumber(:)=Make_Range(self%wavenumberMinimum,self%wavenumberMaximum,int(wavenumberCount),rangeTypeLogarithmic)
    ! Iterate over outputs.
    integrator_=integrator(varianceIntegrand,toleranceRelative=1.0d-2,integrationRule=GSL_Integ_Gauss15)
    do iOutput=1,outputCount
       epochTime                (iOutput)=                                                                     self%outputTimes_%time(iOutput)
       epochRedshift            (iOutput)=self%cosmologyFunctions_ %redshiftFromExpansionFactor         (                                       &
            &                              self%cosmologyFunctions_%expansionFactor                      (                                      &
            &                                                                                                  self%outputTimes_%time(iOutput)  &
            &                                                                                            )                                      &
            &                                                                                           )
       ! Iterate over all wavenumbers computing power spectrum and related quantities.
       do iWavenumber=1,wavenumberCount
          ! Compute corresponding mass scale.
          massScale                (iWavenumber        )=+4.0d0                                       &
               &                                         /3.0d0                                       &
               &                                         *Pi                                          &
               &                                         *self%cosmologyParameters_%OmegaMatter    () &
               &                                         *self%cosmologyParameters_%densityCritical() &
               &                                         /                                                                                                                    wavenumber(iWavenumber) **3
          ! Compute linear growth factors.
          growthFactor             (iWavenumber,iOutput)=self%linearGrowth_             %value                               (time=self%outputTimes_%time(iOutput),wavenumber=wavenumber(iWavenumber))
          growthFactorLogDerivative(iWavenumber,iOutput)=self%linearGrowth_             %logarithmicDerivativeExpansionFactor(time=self%outputTimes_%time(iOutput),wavenumber=wavenumber(iWavenumber))
          ! Compute power spectrum.
          powerSpectrumLinear      (iWavenumber,iOutput)=+self%powerSpectrum_           %power                               (time=self%outputTimes_%time(iOutput),wavenumber=wavenumber(iWavenumber))
          ! Compute transfer function.
          transferFunction         (iWavenumber,iOutput)=+self%transferFunction_        %value                               (                                     wavenumber=wavenumber(iWavenumber))
          ! Compute fluctuation on this mass scale.
          sigma                    (iWavenumber,iOutput)=+self%cosmologicalMassVariance_%rootVariance                        (time=self%outputTimes_%time(iOutput),mass      =massScale (iWavenumber))
          ! Compute gradient of mass fluctuations.
          sigmaGradient            (iWavenumber,iOutput)=+self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient     (time=self%outputTimes_%time(iOutput),mass      =massScale (iWavenumber))
          ! Include non-linear power spectrum if requested.
          if (self%includeNonLinear) then
             powerSpectrumNonLinear(iWavenumber,iOutput)=self%powerSpectrumNonlinear_%value(wavenumber(iWavenumber),self%outputTimes_%time(iOutput))
             ! Compute the variance in the non-linear power spectrum.
             wavenumberMinimum=    0.0d0
             wavenumberMaximum=min(1.0d3*wavenumber(iWavenumber),self%powerSpectrumWindowFunction_%wavenumberMaximum(massScale(iWavenumber)))
             sigmaNonLinear(iWavenumber,iOutput)=+sqrt(                                                            &
                  &                                    +integrator_%integrate(wavenumberMinimum,wavenumberMaximum) &
                  &                                    /2.0d0                                                      &
                  &                                    /Pi**2                                                      &
                  &                                   )
          end if
       end do
    end do
    ! Open the group for output time information.
    if (self%outputGroup == ".") then
       outputsGroup  =outputFile    %openGroup(     'Outputs'        ,'Group containing datasets relating to output times.')
    else
       containerGroup=outputFile    %openGroup(char(self%outputGroup),'Group containing power spectrum data.'              )
       outputsGroup  =containerGroup%openGroup(     'Outputs'        ,'Group containing datasets relating to output times.')
    end if
    ! Iterate over output times and output data.
    do iOutput=1,outputCount
       groupName  ='Output'
       description='Data for output number '
       groupName  =groupName  //iOutput
       description=description//iOutput
       outputGroup=outputsGroup%openGroup(char(groupName),char(description))
       call outputGroup   %writeAttribute(epochRedshift            (  iOutput),'outputRedshift'                                                                                                                            )
       call outputGroup   %writeAttribute(epochTime                (  iOutput),'outputTime'                                                                                                                                )
       call outputGroup   %writeDataset  (wavenumber                          ,'wavenumber'               ,'The wavenumber.'                                                                       ,datasetReturned=dataset)
       call dataset       %writeAttribute(1.0d0/megaParsec                    ,'unitsInSI'                                                                                                                                 )
       call outputGroup   %writeDataset  (massScale                           ,'mass'                     ,'The corresponding mass scale.'                                                         ,datasetReturned=dataset)
       call dataset       %writeAttribute(massSolar                           ,'unitsInSI'                                                                                                                                 )
       call outputGroup   %writeDataset  (growthFactor             (:,iOutput),'growthFactor'             ,'Linear theory growth factor, D(t).'                                                                            )
       call outputGroup   %writeDataset  (growthFactorLogDerivative(:,iOutput),'growthFactorLogDerivative','Logarithmic derivative of growth factor with respect to expansion factor, dlogD/dloga.'                        )
       call outputGroup   %writeDataset  (powerSpectrumLinear      (:,iOutput),'powerSpectrum'            ,'The power spectrum.'                                                                   ,datasetReturned=dataset)
       call dataset       %writeAttribute(megaParsec**3                       ,'unitsInSI'                                                                                                                                 )
       call outputGroup   %writeDataset  (transferFunction         (:,iOutput),'transferFunction'         ,'The transfer function.'                                                                                        )
       call outputGroup   %writeDataset  (sigma                    (:,iOutput),'sigma'                    ,'The mass fluctuation on this scale.'                                                                           )
       call outputGroup   %writeDataset  (sigmaGradient            (:,iOutput),'alpha'                    ,'Logarithmic deriative of the mass flucation with respect to mass.'                                             )
       if (self%includeNonLinear) then
          call outputGroup%writeDataset  (powerSpectrumNonLinear   (:,iOutput),'powerSpectrumNonlinear'   ,'The non-linear power spectrum.'                                                        ,datasetReturned=dataset)
          call dataset    %writeAttribute(megaParsec**3                       ,'unitsInSI'                                                                                                                                 )
          call outputGroup%writeDataset  (sigmaNonLinear           (:,iOutput),'sigmaNonlinear'           ,'The non-linear mass fluctuation on this scale.'                                                                )
       end if
    end do
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: power spectrum' )
    return

  contains

    double precision function varianceIntegrand(wavenumber)
      !!{
      Integrand function used in compute the variance in (real space) top-hat spheres from the power spectrum.
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumber

      ! Return power spectrum multiplied by window function and volume element in k-space. Factors of 2 and Pi are included
      ! elsewhere.
      varianceIntegrand=+  self%powerSpectrumNonlinear_     %value(wavenumber,epochTime(iOutput    )) &
           &            *(                                                                            &
           &              +self%powerSpectrumWindowFunction_%value(wavenumber,massScale(iWavenumber)) &
           &              *                                        wavenumber                         &
           &             )**2
      return
    end function varianceIntegrand

 end subroutine powerSpectraPerform
