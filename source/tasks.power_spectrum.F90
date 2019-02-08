!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  use Cosmology_Parameters
  use Cosmology_Functions
  use Cosmological_Density_Field
  use Power_Spectra
  use Power_Spectrum_Window_Functions
  use Power_Spectra_Nonlinear
  use Linear_Growth
  use Output_Times

  !# <task name="taskPowerSpectra">
  !#  <description>A task which computes and outputs the power spectrum and related quantities.</description>
  !# </task>
  type, extends(taskClass) :: taskPowerSpectra
     !% Implementation of a task which computes and outputs the power spectrum and related quantities.
     private
     class           (cosmologyParametersClass         ), pointer :: cosmologyParameters_ => null()
     class           (cosmologyFunctionsClass          ), pointer :: cosmologyFunctions_ => null()
     class           (linearGrowthClass                ), pointer :: linearGrowth_ => null()
     class           (powerSpectrumClass               ), pointer :: powerSpectrum_ => null()
     class           (powerSpectrumNonlinearClass      ), pointer :: powerSpectrumNonlinear_ => null()
     class           (powerSpectrumWindowFunctionClass ), pointer :: powerSpectrumWindowFunction_ => null()
     class           (cosmologicalMassVarianceClass    ), pointer :: cosmologicalMassVariance_ => null()
     class           (outputTimesClass                 ), pointer :: outputTimes_ => null()
     double precision                                             :: wavenumberMinimum           , wavenumberMaximum
     integer                                                      :: pointsPerDecade
     logical                                                      :: includeNonLinear
     type            (varying_string                   )          :: outputGroup
   contains
     final     ::            powerSpectraDestructor
     procedure :: perform => powerSpectraPerform
  end type taskPowerSpectra

  interface taskPowerSpectra
     !% Constructors for the {\normalfont \ttfamily powerSpectrum} task.
     module procedure powerSpectraConstructorParameters
     module procedure powerSpectraConstructorInternal
  end interface taskPowerSpectra

contains

  function powerSpectraConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily powerSpectrum} task class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (taskPowerSpectra                )                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (cosmologyParametersClass        ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    class           (linearGrowthClass               ), pointer       :: linearGrowth_
    class           (powerSpectrumClass              ), pointer       :: powerSpectrum_
    class           (powerSpectrumNonlinearClass     ), pointer       :: powerSpectrumNonlinear_
    class           (powerSpectrumWindowFunctionClass), pointer       :: powerSpectrumWindowFunction_
    class           (cosmologicalMassVarianceClass   ), pointer       :: cosmologicalMassVariance_
    class           (outputTimesClass                ), pointer       :: outputTimes_
    double precision                                                  :: wavenumberMinimum           , wavenumberMaximum
    integer                                                           :: pointsPerDecade
    logical                                                           :: includeNonLinear
    type            (varying_string                  )                :: outputGroup

    !# <inputParameter>
    !#   <name>wavenumberMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d-3</defaultValue>
    !#   <description>The minimum wavenumber at which to tabulate power spectra.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>wavenumberMaximum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d+3</defaultValue>
    !#   <description>The maximum wavenumber at which to tabulate power spectra.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>pointsPerDecade</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>10</defaultValue>
    !#   <description>The number of points per decade of wavenumber at which to tabulate power spectra.</description>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>includeNonLinear</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>If true the nonlinear power spectrum is also computed and output.</description>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>outputGroup</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>var_str('.')</defaultValue>
    !#   <description>The HDF5 output group within which to write power spectrum data.</description>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>    
    !# <objectBuilder class="cosmologyParameters"         name="cosmologyParameters_"         source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"          name="cosmologyFunctions_"          source="parameters"/>
    !# <objectBuilder class="linearGrowth"                name="linearGrowth_"                source="parameters"/>
    !# <objectBuilder class="powerSpectrum"               name="powerSpectrum_"               source="parameters"/>
    !# <objectBuilder class="powerSpectrumNonlinear"      name="powerSpectrumNonlinear_"      source="parameters"/>
    !# <objectBuilder class="powerSpectrumWindowFunction" name="powerSpectrumWindowFunction_" source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance"    name="cosmologicalMassVariance_"    source="parameters"/>
    !# <objectBuilder class="outputTimes"                 name="outputTimes_"                 source="parameters"/>
    self=taskPowerSpectra(                               &
         &                 wavenumberMinimum           , &
         &                 wavenumberMaximum           , &
         &                 pointsPerDecade             , &
         &                 includeNonLinear            , &
         &                 outputGroup                 , &
         &                 cosmologyParameters_        , &
         &                 cosmologyFunctions_         , &
         &                 linearGrowth_               , &
         &                 powerSpectrum_              , &
         &                 powerSpectrumNonlinear_     , &
         &                 powerSpectrumWindowFunction_, &
         &                 cosmologicalMassVariance_   , &
         &                 outputTimes_                  &
         &                )
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"        />
    !# <objectDestructor name="cosmologyFunctions_"         />
    !# <objectDestructor name="linearGrowth_"               />
    !# <objectDestructor name="powerSpectrum_"              />
    !# <objectDestructor name="powerSpectrumNonlinear_"     />
    !# <objectDestructor name="powerSpectrumWindowFunction_"/>
    !# <objectDestructor name="cosmologicalMassVariance_"   />
    !# <objectDestructor name="outputTimes_"                />
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
       &                                    powerSpectrum_              , &
       &                                    powerSpectrumNonlinear_     , &
       &                                    powerSpectrumWindowFunction_, &
       &                                    cosmologicalMassVariance_   , &
       &                                    outputTimes_                  &
       &                                   ) result(self)
    !% Internal constructor for the {\normalfont \ttfamily powerSpectrum} task class.
    implicit none
    type            (taskPowerSpectra                )                        :: self
    class           (cosmologyParametersClass        ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass         ), intent(in   ), target :: cosmologyFunctions_
    class           (linearGrowthClass               ), intent(in   ), target :: linearGrowth_
    class           (powerSpectrumClass              ), intent(in   ), target :: powerSpectrum_
    class           (powerSpectrumNonlinearClass     ), intent(in   ), target :: powerSpectrumNonlinear_
    class           (powerSpectrumWindowFunctionClass), intent(in   ), target :: powerSpectrumWindowFunction_
    class           (cosmologicalMassVarianceClass   ), intent(in   ), target :: cosmologicalMassVariance_
    class           (outputTimesClass                ), intent(in   ), target :: outputTimes_
    double precision                                  , intent(in   )         :: wavenumberMinimum           , wavenumberMaximum
    integer                                           , intent(in   )         :: pointsPerDecade
    logical                                           , intent(in   )         :: includeNonLinear
    type            (varying_string                  ), intent(in   )         :: outputGroup
    !# <constructorAssign variables="wavenumberMinimum, wavenumberMaximum, pointsPerDecade, includeNonLinear, outputGroup,*cosmologyParameters_,*cosmologyFunctions_,*linearGrowth_,*powerSpectrum_,*powerSpectrumNonlinear_,*powerSpectrumWindowFunction_,*cosmologicalMassVariance_, *outputTimes_"/>
    
    return
  end function powerSpectraConstructorInternal
  
  subroutine powerSpectraDestructor(self)
    !% Destructor for the {\normalfont \ttfamily powerSpectrum} task class.
    implicit none
    type(taskPowerSpectra), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"        />
    !# <objectDestructor name="self%cosmologyFunctions_"         />
    !# <objectDestructor name="self%linearGrowth_"               />
    !# <objectDestructor name="self%powerSpectrum_"              />
    !# <objectDestructor name="self%powerSpectrumNonlinear_"     />
    !# <objectDestructor name="self%powerSpectrumWindowFunction_"/>
    !# <objectDestructor name="self%cosmologicalMassVariance_"   />
    !# <objectDestructor name="self%outputTimes_"                />
    return
  end subroutine powerSpectraDestructor

  subroutine powerSpectraPerform(self)
    !% Compute and output the halo mass function.
    use, intrinsic :: ISO_C_Binding
    use            :: Galacticus_HDF5
    use            :: Galacticus_Display
    use            :: Numerical_Ranges
    use            :: FGSL                            , only : fgsl_function, fgsl_integration_workspace, FGSL_Integ_Gauss15
    use            :: Memory_Management
    use            :: Numerical_Constants_Astronomical
    use            :: IO_HDF5
    use            :: String_Handling
    use            :: Numerical_Integration
    implicit none
    class           (taskPowerSpectra          ), intent(inout)               :: self
    integer         (c_size_t                  )                              :: outputCount           , wavenumberCount    , &
         &                                                                       iOutput               , iWavenumber
    double precision                            , allocatable, dimension(:  ) :: wavenumber            , powerSpectrumLinear, &
         &                                                                       massScale             , sigma              , &
         &                                                                       sigmaGradient         , growthFactor       , &
         &                                                                       epochTime             , epochRedshift
    double precision                            , allocatable, dimension(:,:) :: powerSpectrumNonLinear, sigmaNonLinear
    double precision                                                          :: wavenumberMinimum     , wavenumberMaximum
    type            (fgsl_function             )                              :: integrandFunction
    type            (fgsl_integration_workspace)                              :: integrationWorkspace
    type            (hdf5Object                )                              :: outputsGroup          , outputGroup        , &
         &                                                                       containerGroup        , dataset
    type            (varying_string            )                              :: groupName             , commentText

    call Galacticus_Display_Indent('Begin task: power spectrum')
    ! Get the requested output redshifts.
    outputCount      =self%outputTimes_%count()
    ! Compute number of tabulation points.
    wavenumberCount=int(log10(self%wavenumberMaximum/self%wavenumberMinimum)*dble(self%pointsPerDecade))+1
    ! Allocate arrays for power spectra.
    call    allocateArray(wavenumber            ,[wavenumberCount            ])
    call    allocateArray(powerSpectrumLinear   ,[wavenumberCount            ])
    call    allocateArray(massScale             ,[wavenumberCount            ])
    call    allocateArray(sigma                 ,[wavenumberCount            ])
    call    allocateArray(sigmaGradient         ,[wavenumberCount            ])
    call    allocateArray(growthFactor          ,[                outputCount])
    call    allocateArray(epochTime             ,[                outputCount])
    call    allocateArray(epochRedshift         ,[                outputCount])
    if (self%includeNonLinear) then
       call allocateArray(powerSpectrumNonLinear,[wavenumberCount,outputCount])
       call allocateArray(sigmaNonLinear        ,[wavenumberCount,outputCount])
    end if
    ! Build a range of wavenumbers.
    wavenumber(:)=Make_Range(self%wavenumberMinimum,self%wavenumberMaximum,int(wavenumberCount),rangeTypeLogarithmic)
    ! Iterate over all wavenumbers computing power spectrum and related quantities.
    do iWavenumber=1,wavenumberCount
       ! Compute power spectrum.
       powerSpectrumLinear(iWavenumber)=+self%powerSpectrum_%power(wavenumber(iWavenumber))
       ! Compute corresponding mass scale.
       massScale          (iWavenumber)=+4.0d0                                       &
            &                           /3.0d0                                       &
            &                           *Pi                                          &
            &                           *self%cosmologyParameters_%OmegaMatter    () &
            &                           *self%cosmologyParameters_%densityCritical() &
            &                           /                                                               wavenumber(iWavenumber) **3
       ! Compute fluctuation on this mass scale.
       sigma              (iWavenumber)=+self%cosmologicalMassVariance_%rootVariance                   (massScale (iWavenumber))
       ! Compute gradient of mass fluctuations.
       sigmaGradient      (iWavenumber)=+self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient(massScale (iWavenumber))
    end do
    ! Iterate over outputs.
    do iOutput=1,outputCount
       epochTime    (iOutput)=                                                            self%outputTimes_%time(iOutput)
       epochRedshift(iOutput)=self%cosmologyFunctions_ %redshiftFromExpansionFactor(                                      &
            &                  self%cosmologyFunctions_%expansionFactor             (                                     &
            &                                                                             self%outputTimes_%time(iOutput) &
            &                                                                       )                                     &
            &                                                                      )
       growthFactor (iOutput)=self%linearGrowth_       %value                      ( time=self%outputTimes_%time(iOutput))
       ! Iterate over all wavenumbers computing non-linear power spectrum.
       if (self%includeNonLinear) then
          do iWavenumber=1,wavenumberCount
             powerSpectrumNonLinear(iWavenumber,iOutput)=self%powerSpectrumNonlinear_%value(wavenumber(iWavenumber),self%outputTimes_%time(iOutput))
             ! Compute the variance in the non-linear power spectrum.
             wavenumberMinimum=    0.0d0
             wavenumberMaximum=min(1.0d3*wavenumber(iWavenumber),self%powerSpectrumWindowFunction_%wavenumberMaximum(massScale(iWavenumber)))
             sigmaNonLinear(iWavenumber,iOutput)=+sqrt(                                                   &
                  &                                    +Integrate(                                        &
                  &                                                                 wavenumberMinimum   , &
                  &                                                                 wavenumberMaximum   , &
                  &                                                                 varianceIntegrand   , &
                  &                                                                 integrandFunction   , &
                  &                                                                 integrationWorkspace, &
                  &                                               toleranceAbsolute=0.0d0               , &
                  &                                               toleranceRelative=1.0d-2              , &
                  &                                               integrationRule  =FGSL_Integ_Gauss15    &
                  &                                              )                                        &
                  &                                    /2.0d0                                             &
                  &                                    /Pi**2                                             &
                  &                                   )
             call Integrate_Done(integrandFunction,integrationWorkspace)
          end do
       end if
    end do
    ! Open the group for output time information.
    if (self%outputGroup == ".") then
       outputsGroup  =galacticusOutputFile%openGroup(     'Outputs'        ,'Group containing datasets relating to output times.')
    else
       containerGroup=galacticusOutputFile%openGroup(char(self%outputGroup),'Group containing power spectrum data.'              )
       outputsGroup  =containerGroup      %openGroup(     'Outputs'        ,'Group containing datasets relating to output times.')
    end if
    ! Iterate over output times and output data.
    do iOutput=1,outputCount
       groupName  ='Output'
       commentText='Data for output number '
       groupName  =groupName  //iOutput
       commentText=commentText//iOutput
       outputGroup=outputsGroup%openGroup(char(groupName),char(commentText))
       call outputGroup%writeAttribute(growthFactor (iOutput),'growthFactor'  )
       call outputGroup%writeAttribute(epochRedshift(iOutput),'outputRedshift')
       call outputGroup%writeAttribute(epochTime    (iOutput),'outputTime'    )
       if (self%includeNonLinear) then
          call outputGroup%writeDataset  (wavenumber                       ,'wavenumber'            ,'The wavenumber.'                               ,datasetReturned=dataset)
          call dataset    %writeAttribute(1.0d0/megaParsec                 ,'unitsInSI'                                                                                      )
          call dataset    %close         (                                                                                                                                   )
          call outputGroup%writeDataset  (powerSpectrumNonLinear(:,iOutput),'powerSpectrumNonlinear','The non-linear power spectrum.'                ,datasetReturned=dataset)
          call dataset    %writeAttribute(megaParsec**3                    ,'unitsInSI'                                                                                      )
          call dataset    %close         (                                                                                                                                   )
          call outputGroup%writeDataset  (sigmaNonLinear        (:,iOutput),'sigmaNonlinear'        ,'The non-linear mass fluctuation on this scale.'                        )
       end if
       call outputGroup%close()
    end do
    call outputsGroup%close()
    if (self%outputGroup == ".")                                                                                  &
         & containerGroup=galacticusOutputFile%openGroup('powerSpectrum','Group containing power spectrum data.')
    call containerGroup%writeDataset  (wavenumber         ,'wavenumber'   ,'The wavenumber.'                                                  ,datasetReturned=dataset)
    call dataset       %writeAttribute(1.0d0/megaParsec   ,'unitsInSI'                                                                                                )
    call dataset       %close         (                                                                                                                               )
    call containerGroup%writeDataset  (powerSpectrumLinear,'powerSpectrum','The power spectrum.'                                              ,datasetReturned=dataset)
    call dataset       %writeAttribute(megaParsec**3      ,'unitsInSI'                                                                                                )
    call dataset       %close         (                                                                                                                               )
    call containerGroup%writeDataset  (massScale          ,'mass'         ,'The corresponding mass scale.'                                    ,datasetReturned=dataset)
    call dataset       %writeAttribute(massSolar          ,'unitsInSI'                                                                                                )
    call dataset       %close         (                                                                                                                               )
    call containerGroup%writeDataset  (sigma              ,'sigma'        ,'The mass fluctuation on this scale.'                                                      )
    call containerGroup%writeDataset  (sigmaGradient      ,'alpha'        ,'Logarithmic deriative of the mass flucation with respect to mass.'                        )
    ! Close the datasets group.
    call containerGroup%close         (                                                                                                                               )
    call Galacticus_Display_Unindent('Done task: power spectrum' )
    return

  contains

    double precision function varianceIntegrand(wavenumber)
      !% Integrand function used in compute the variance in (real space) top-hat spheres from the power spectrum.
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
