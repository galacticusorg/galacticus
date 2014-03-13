!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which handles evolving the spectrum of background radiation.

module Radiation_Intergalactic_Background_Internal
  !% Handles evolving the spectrum of background radiation.
  use Cosmology_Functions
  use Cosmology_Parameters
  use Intergalactic_Medium_State
  use Abundances_Structure
  use FGSL
  private
  public :: Radiation_Intergalactic_Background_Internal_Initialize
  
  logical                                                                      :: backgroundRadiationCompute
  integer                                                                      :: backgroundRadiationWavelengthCountPerDecade, backgroundRadiationWavelengthCount
  double precision                                                             :: backgroundRadiationWavelengthMinimum       , backgroundRadiationWavelengthMaximum
  integer                                                                      :: backgroundRadiationTimeCountPerDecade      , backgroundRadiationTimeCount
  double precision                                                             :: backgroundRadiationRedshiftMinimum         , backgroundRadiationRedshiftMaximum
  double precision                                                             :: backgroundRadiationTimeMinimum             , backgroundRadiationTimeMaximum
  double precision                               , allocatable, dimension(:  ) :: backgroundRadiationWavelength              , backgroundRadiationSpectrum         , &
       &                                                                          backgroundRadiationTime                    , crossSectionNeutralHydrogen         , &
       &                                                                          crossSectionNeutralHelium                  , crossSectionSinglyIonizedHelium     , &
       &                                                                          backgroundRadiationRedshift
  double precision                               , allocatable, dimension(:,:) :: backgroundRadiationEmissivity              , emissivity                          , &
       &                                                                          backgroundRadiationFlux
  double precision                                            , dimension(0:1) :: emissivityTime

  ! Classes used in ODE solution.
  class           (cosmologyParametersClass     ), pointer                     :: cosmologyParameters_
  class           (cosmologyFunctionsClass      ), pointer                     :: cosmologyFunctions_
  class           (intergalacticMediumStateClass), pointer                     :: intergalacticMediumState_

  ! Stellar population variables used in convolution integral.
  type            (abundances                   ), pointer                     :: gasAbundances
  integer                                                                      :: imfIndex
  double precision                                                             :: wavelength
  logical                                                                      :: integrationReset=.true.
  type            (fgsl_function                )                              :: integrandFunction
  type            (fgsl_integration_workspace   )                              :: integrationWorkspace
 
contains
  
  !# <universePreEvolveTask>
  !#  <unitName>Radiation_Intergalactic_Background_Internal_Initialize</unitName>
  !# </universePreEvolveTask>
  subroutine Radiation_Intergalactic_Background_Internal_Initialize(thisUniverse)
    !% Attach an initial event to a merger tree to cause the background radiation update function to be called.
    use Galacticus_Nodes
    use Input_Parameters
    use Memory_Management
    use Numerical_Ranges
    use Cosmology_Functions
    use Atomic_Cross_Sections_Ionization_Photo
    implicit none
    type   (universe               ), intent(inout) :: thisUniverse
    type   (universeEvent          ), pointer       :: thisEvent
    class  (cosmologyFunctionsClass), pointer       :: cosmologyFunctions_
    integer                                         :: iWavelength        , iTime

    ! Get parameter controlling background radiation spectral and time resolution.
    !@ <inputParameter>
    !@   <name>backgroundRadiationCompute</name>
    !@   <defaultValue>false</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     Specifies whether or not cosmic background radiation should be computed.
    !@   </description>
    !@   <type>integer</type>
    !@   <cardinality>1</cardinality>
    !@ </inputParameter>
    call Get_Input_Parameter('backgroundRadiationCompute',backgroundRadiationCompute,defaultValue=.false.)
    if (backgroundRadiationCompute) then
       !@ <inputParameter>
       !@   <name>backgroundRadiationWavelengthCountPerDecade</name>
       !@   <defaultValue>10</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The number of bins per decade of wavelength to use for calculations of the cosmic background radiation.
       !@   </description>
       !@   <type>integer</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('backgroundRadiationWavelengthCountPerDecade',backgroundRadiationWavelengthCountPerDecade,defaultValue=10)
       !@ <inputParameter>
       !@   <name>backgroundRadiationWavelengthMinimum</name>
       !@   <defaultValue>100\AA</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The minimum wavelength (in units of \AA) to use in calculations of the cosmic background radiation.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('backgroundRadiationWavelengthMinimum',backgroundRadiationWavelengthMinimum,defaultValue=100.0d0)
       !@ <inputParameter>
       !@   <name>backgroundRadiationWavelengthMaximum</name>
       !@   <defaultValue>100000\AA</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The maximum wavelength (in units of \AA) to use in calculations of the cosmic background radiation.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('backgroundRadiationWavelengthMaximum',backgroundRadiationWavelengthMaximum,defaultValue=100000.0d0)
       !@ <inputParameter>
       !@   <name>backgroundRadiationTimeCountPerDecade</name>
       !@   <defaultValue>10</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The number of bins per decade of time to use for calculations of tge cosmic background radiation.
       !@   </description>
       !@   <type>integer</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('backgroundRadiationTimeCountPerDecade',backgroundRadiationTimeCountPerDecade,defaultValue=10)
       !@ <inputParameter>
       !@   <name>backgroundRadiationRedshiftMinimum</name>
       !@   <defaultValue>0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The minimum redshift to use in calculations of the cosmic background radiation.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('backgroundRadiationRedshiftMinimum',backgroundRadiationRedshiftMinimum,defaultValue=0.0d0)
       !@ <inputParameter>
       !@   <name>backgroundRadiationRedshiftMaximum</name>
       !@   <defaultValue>30</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The maximum redshift to use in calculations of the cosmic background radiation.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('backgroundRadiationRedshiftMaximum',backgroundRadiationRedshiftMaximum,defaultValue=30.0d0)
       ! Build tables of wavelength and time for cosmic background radiation.
       cosmologyParameters_      => cosmologyParameters     ()
       cosmologyFunctions_       => cosmologyFunctions      ()
       intergalacticMediumState_ => intergalacticMediumState()
       backgroundRadiationTimeMaximum                                                              &
            & =cosmologyFunctions_%cosmicTime                 (                                    &
            &  cosmologyFunctions_%expansionFactorFromRedshift (                                   &
            &                                                   backgroundRadiationRedshiftMinimum &
            &                                                  )                                   &
            &                                                 )
       backgroundRadiationTimeMinimum                                                              &
            & =cosmologyFunctions_%cosmicTime                 (                                    &
            &  cosmologyFunctions_%expansionFactorFromRedshift (                                   &
            &                                                   backgroundRadiationRedshiftMaximum &
            &                                                  )                                   &
            &                                                 )
       backgroundRadiationWavelengthCount                              &
            & = int(                                                   &
            &        dble(backgroundRadiationWavelengthCountPerDecade) &
            &       *log10(                                            &
            &               backgroundRadiationWavelengthMaximum       &
            &              /backgroundRadiationWavelengthMinimum       &
            &             )                                            &
            &      )                                                   &
            &  +1
       backgroundRadiationTimeCount                                    &
            & = int(                                                   &
            &        dble(backgroundRadiationTimeCountPerDecade      ) &
            &       *log10(                                            &
            &               backgroundRadiationTimeMaximum             &
            &              /backgroundRadiationTimeMinimum             &
            &             )                                            &
            &      )                                                   &
            &  +1
       call Alloc_Array(backgroundRadiationWavelength,[backgroundRadiationWavelengthCount                             ]                  )
       call Alloc_Array(backgroundRadiationSpectrum  ,[backgroundRadiationWavelengthCount                             ]                  )
       call Alloc_Array(backgroundRadiationTime      ,[                                   backgroundRadiationTimeCount]                  )
       call Alloc_Array(backgroundRadiationRedshift  ,[                                   backgroundRadiationTimeCount]                  )
       call Alloc_Array(backgroundRadiationEmissivity,[backgroundRadiationWavelengthCount,backgroundRadiationTimeCount]                  )
       call Alloc_Array(backgroundRadiationFlux      ,[backgroundRadiationWavelengthCount,backgroundRadiationTimeCount]                  )
       call Alloc_Array(emissivity                   ,[backgroundRadiationWavelengthCount,2                           ],lowerBounds=[1,0])
       backgroundRadiationWavelength                            &
            & =Make_Range(                                      &
            &             backgroundRadiationWavelengthMinimum, &
            &             backgroundRadiationWavelengthMaximum, &
            &             backgroundRadiationWavelengthCount  , &
            &             rangeTypeLogarithmic                  &
            &            )
       backgroundRadiationTime                                  &
            & =Make_Range(                                      &
            &             backgroundRadiationTimeMinimum      , &
            &             backgroundRadiationTimeMaximum      , &
            &             backgroundRadiationTimeCount        , &
            &             rangeTypeLogarithmic                  &
            &            )
       ! Convert times to redshifts.
       do iTime=1,backgroundRadiationTimeCount
          backgroundRadiationRedshift(iTime)                                                       &
               & =cosmologyFunctions_ %redshiftFromExpansionFactor(                                &
               &   cosmologyFunctions_%expansionFactor             (                               &
               &                                                    backgroundRadiationTime(iTime) &
               &                                                   )                               &
               &                                                  )
       end do
       ! Initialize the background radiation to zero.
       backgroundRadiationSpectrum=0.0d0
       ! Initialize the emissivity to zero.
       backgroundRadiationEmissivity=0.0d0
       ! Construct tables of photoionization cross-sections.
       call Alloc_Array(crossSectionNeutralHydrogen    ,[backgroundRadiationWavelengthCount])
       call Alloc_Array(crossSectionNeutralHelium      ,[backgroundRadiationWavelengthCount])
       call Alloc_Array(crossSectionSinglyIonizedHelium,[backgroundRadiationWavelengthCount])
       do iWavelength=1,backgroundRadiationWavelengthCount
          crossSectionNeutralHydrogen    (iWavelength)=Atomic_Cross_Section_Ionization_Photo(1,1,1,backgroundRadiationWavelength(iWavelength))
          crossSectionNeutralHelium      (iWavelength)=Atomic_Cross_Section_Ionization_Photo(2,1,1,backgroundRadiationWavelength(iWavelength))
          crossSectionSinglyIonizedHelium(iWavelength)=Atomic_Cross_Section_Ionization_Photo(2,2,1,backgroundRadiationWavelength(iWavelength))
       end do
       ! Create the first interrupt event in the universe object.
       thisEvent      => thisUniverse%createEvent()
       thisEvent%time =  backgroundRadiationTime(1)
       thisEvent%task => Radiation_Intergalactic_Background_Internal_Update
    end if
    return
  end subroutine Radiation_Intergalactic_Background_Internal_Initialize
  
  logical function Radiation_Intergalactic_Background_Internal_Update(thisEvent,thisUniverse) result (success)
    !% Update the radiation background for a given universe.
    use               Galacticus_Nodes
    use               Galacticus_Display
    use               Star_Formation_IMF
    use               Galactic_Structure_Options
    use               Galacticus_Error
    use               Stellar_Population_Spectra
    use               Arrays_Search
    use               FODEIV2
    use               ODEIV2_Solver
    use, intrinsic :: ISO_C_Binding
    use               Numerical_Constants_Prefixes
    use               Numerical_Constants_Math
    use               Numerical_Constants_Physical
    use               Numerical_Constants_Units
    use               Numerical_Integration
    use               ISO_Varying_String
    use               Galacticus_HDF5
    use               IO_HDF5
    implicit none
    class           (universeEvent             ), intent(in   ) :: thisEvent
    type            (universe                  ), intent(inout) :: thisUniverse
    type            (mergerTree                ), pointer       :: thisTree       
    type            (mergerTreeList            ), pointer       :: thisForest       
    type            (treeNode                  ), pointer       :: thisNode
    class           (nodeComponentBasic        ), pointer       :: thisBasic
    class           (nodeComponentDisk         ), pointer       :: thisDisk
    class           (nodeComponentSpheroid     ), pointer       :: thisSpheroid
    type            (universeEvent             ), pointer       :: newEvent
    double precision                            , parameter     :: odeToleranceAbsolute        =1.0d-30, odeToleranceRelative        =1.0d-3
    double precision                            , parameter     :: integrationToleranceAbsolute=1.0d-30, integrationToleranceRelative=1.0d-3
    type            (fodeiv2_system            ), save          :: ode2System
    type            (fodeiv2_driver            ), save          :: ode2Driver
    logical                                     , save          :: odeReset
    type            (c_ptr                     )                :: parameterPointer
    double precision                                            :: starFormationRateDisk               , starFormationRateSpheroid          , &
         &                                                         gasMassDisk                         , gasMassSpheroid                    , &
         &                                                         ageEnd                              , ageStart                           , &
         &                                                         stellarSpectrumDisk                 , stellarSpectrumSpheroid
    type            (abundances                ), target        :: gasAbundancesDisk                   , gasAbundancesSpheroid
    integer                                                     :: imfIndexDisk                        , imfIndexSpheroid                   , &
         &                                                         iTime                               , iWavelength                        , &
         &                                                         statusDisk                          , statusSpheroid                     , &
         &                                                         iNow
    type            (varying_string            )                :: message
    character       (len=6                     )                :: label
    type            (hdf5Object                )                :: backgroundRadiationGroup            , backgroundRadiationDataset

    ! Display message.
    write (label,'(f6.3)') thisEvent%time
    message="Evolving cosmic background radiation to time "//trim(label)//" Gyr"
    call Galacticus_Display_Indent(message)
    ! Find the current timestep.
    iNow=Search_Array_For_Closest(backgroundRadiationTime,thisEvent%time)
    ! Iterate over all nodes.
    call Galacticus_Display_Message('Accumulating emissivity')
    thisForest => thisUniverse%trees
    do while (associated(thisForest))
       thisTree => thisForest%tree
       do while (associated(thisTree))
          thisNode => thisTree%baseNode
          do while (associated(thisNode))
             thisBasic => thisNode%basic()
             if (thisBasic%time() == thisEvent%time) then
                ! Get the star formation rates and metallicites for this node.
                thisDisk                  => thisNode    %disk             ()
                thisSpheroid              => thisNode    %spheroid         ()
                starFormationRateDisk     =  thisDisk    %starFormationRate()
                starFormationRateSpheroid =  thisSpheroid%starFormationRate()
                gasMassDisk               =  thisDisk    %massGas          ()
                gasMassSpheroid           =  thisSpheroid%massGas          ()
                gasAbundancesDisk         =  thisDisk    %abundancesGas    ()
                gasAbundancesSpheroid     =  thisSpheroid%abundancesGas    ()
                if (starFormationRateDisk     > 0.0d0) gasAbundancesDisk    =gasAbundancesDisk    /gasMassDisk
                if (starFormationRateSpheroid > 0.0d0) gasAbundancesSpheroid=gasAbundancesSpheroid/gasMassSpheroid
                if (starFormationRateDisk > 0.0d0 .or. starFormationRateSpheroid > 0.0d0) then
                   ! Find IMF indices for disk and spheroid.
                   imfIndexDisk    =IMF_Select(starFormationRateDisk    ,gasAbundancesDisk    ,componentTypeDisk    )
                   imfIndexSpheroid=IMF_Select(starFormationRateSpheroid,gasAbundancesSpheroid,componentTypeSpheroid)
                   ! Find the duration of the current timestep.
                   ! Accumulate emissivity to each timestep.
                   do iTime=1,backgroundRadiationTimeCount
                      ! Skip times in the past.
                      if (backgroundRadiationTime(iTime) < thisEvent%time) cycle
                      ! Compute age of the currently forming population at this time.
                      ageEnd=backgroundRadiationTime(iTime)-thisEvent%time
                      if (iTime == 1) then
                         ageStart=0.0d0
                      else
                         ageStart=max(backgroundRadiationTime(iTime-1)-thisEvent%time,0.0d0)
                      end if
                      ! Iterate over wavelength
                      do iWavelength=1,backgroundRadiationWavelengthCount                         
                         wavelength              =  backgroundRadiationWavelength(iWavelength)
                         imfIndex                =  imfIndexDisk
                         gasAbundances           => gasAbundancesDisk
                         stellarSpectrumDisk     =  Integrate(                                                        &
                              &                               ageStart                                              , &
                              &                               ageEnd                                                , &
                              &                               stellarSpectraConvolution                             , &
                              &                               parameterPointer                                      , &
                              &                               integrandFunction                                     , &
                              &                               integrationWorkspace                                  , &
                              &                               toleranceAbsolute        =integrationToleranceAbsolute, &
                              &                               toleranceRelative        =integrationToleranceRelative, &
                              &                               reset                    =integrationReset              &
                              &                              )
                         gasAbundances           => gasAbundancesSpheroid
                         stellarSpectrumSpheroid =  Integrate(                                                        &
                              &                               ageStart                                              , &
                              &                               ageEnd                                                , &
                              &                               stellarSpectraConvolution                             , &
                              &                               parameterPointer                                      , &
                              &                               integrandFunction                                     , &
                              &                               integrationWorkspace                                  , &
                              &                               toleranceAbsolute        =integrationToleranceAbsolute, &
                              &                               toleranceRelative        =integrationToleranceRelative, &
                              &                               reset                    =integrationReset              &
                              &                              )
                         backgroundRadiationEmissivity        (iWavelength,iTime) &
                              & =backgroundRadiationEmissivity(iWavelength,iTime) &
                              & +(                                                &
                              &   +stellarSpectrumDisk                            &
                              &   *starFormationRateDisk                          &
                              &   +stellarSpectrumSpheroid                        &
                              &   *starFormationRateSpheroid                      &
                              &  )                                                &
                              & *thisTree%volumeWeight
                      end do
                   end do
                end if
             end if
             call thisNode%walkTreewithsatellites(thisNode)
          end do
          thisTree => thisTree%nextTree
       end do
       thisForest => thisForest%next
    end do
    ! Evolve the cosmic background radiation up to this timestep.
    if (iNow > 1) then
       call Galacticus_Display_Message('Solving cosmic background radiation evolution')
       emissivityTime(  0:1)=backgroundRadiationTime      (  iNow-1:iNow)
       emissivity    (:,0  )=backgroundRadiationEmissivity(:,iNow-1     )
       emissivity    (:,  1)=backgroundRadiationEmissivity(:,       iNow)
       cosmologyFunctions_ => cosmologyFunctions()
       odeReset=.true.
       call ODEIV2_Solve(                                            &
            &            ode2Driver                                , &
            &            ode2System                                , &
            &            backgroundRadiationTime           (iNow-1), &
            &            backgroundRadiationTime           (iNow  ), &
            &            backgroundRadiationWavelengthCount        , &
            &            backgroundRadiationSpectrum               , &
            &            backgroundRadiationODEs                   , &
            &            parameterPointer                          , &
            &            odeToleranceAbsolute                      , &
            &            odeToleranceRelative                      , &
            &            reset=odeReset                              &
            &           )
       call ODEIV2_Solver_Free(ode2Driver,ode2System)
       ! Convert
       backgroundRadiationFlux(:,iNow)           &
            & =+plancksConstant                  &
            &  *speedLight                   **2 &
            &  *angstromsPerMeter                &
            &  /4.0d0                            &
            &  /Pi                               &
            &  /backgroundRadiationWavelength    &
            &  *backgroundRadiationSpectrum      &
            &  *centi                        **2 &
            &  /ergs
    end if
    ! Add the next event to the universe.
    if (iNow < backgroundRadiationTimeCount) then
       newEvent      => thisUniverse%createEvent()
       newEvent%time =  backgroundRadiationTime(iNow+1)
       newEvent%task => Radiation_Intergalactic_Background_Internal_Update
    else
       ! Output the results to file.
       !$omp critical (HDF5_Access)
       backgroundRadiationGroup=galacticusOutputFile%openGroup('backgroundRadiation','Cosmic background radiation data.')
       call backgroundRadiationGroup  %writeDataset  (backgroundRadiationWavelength,'wavelength','Wavelength at which the background radiation is tabulated [Å].',datasetReturned=backgroundRadiationDataset)
       call backgroundRadiationDataset%writeAttribute(1.0d0/angstromsPerMeter      ,'unitsInSI'                                                                                                             )
       call backgroundRadiationDataset%close()
       call backgroundRadiationGroup  %writeDataset  (backgroundRadiationRedshift  ,'redshift'  ,'Redshift at which the background radiation is tabulated [].'   ,datasetReturned=backgroundRadiationDataset)
       call backgroundRadiationDataset%writeAttribute(0.0d0                        ,'unitsInSI'                                                                                                             )
       call backgroundRadiationDataset%close()
       call backgroundRadiationGroup  %writeDataset  (backgroundRadiationFlux      ,'flux'  ,'Flux is the cosmic background radiation [erg cm⁻² s⁻¹ Hz⁻¹ sr⁻¹].' ,datasetReturned=backgroundRadiationDataset)
       call backgroundRadiationDataset%writeAttribute(ergs/centi**2                ,'unitsInSI'                                                                                                             )
       call backgroundRadiationDataset%close()
       call backgroundRadiationGroup  %close()
       !$omp end critical (HDF5_Access)
    end if
    ! Display message.
    call Galacticus_Display_Unindent('done')
    ! Return true since we've performed our task.
    success=.true.
    return
  end function Radiation_Intergalactic_Background_Internal_Update

  function stellarSpectraConvolution(age,parameterPointer) bind(c)
    !% Integrand for convolution of stellar spectra.
    use               Stellar_Population_Spectra
    use, intrinsic :: ISO_C_Binding
    use Galacticus_Error
    implicit none
    real(kind=c_double)        :: stellarSpectraConvolution
    real(kind=c_double), value :: age
    type(c_ptr        ), value :: parameterPointer
    integer                    :: status

    stellarSpectraConvolution=Stellar_Population_Spectrum(               &
         &                                                gasAbundances, &
         &                                                age          , &
         &                                                wavelength   , &
         &                                                imfIndex     , &
         &                                                status         &
         &                                               )
    if     (                                                                             &
         &   status /= errorStatusSuccess                                                &
         &  .and.                                                                        &
         &   status /= errorStatusOutOfRange                                             &
         & ) call Galacticus_Error_Report(                                               &
         &                                'stellarSpectraConvolution'                  , &
         &                                'stellar population spectrum function failed'  &
         &                               )
    return
  end function stellarSpectraConvolution
  
  function backgroundRadiationODEs(time,spectrum,spectrumRateOfChange,parameterPointer) bind(c)
    !% Evaluates the ODEs controlling the evolution of cosmic background radiation.
    use               ODE_Solver_Error_Codes
    use, intrinsic :: ISO_C_Binding
    use               Numerical_Constants_Astronomical
    use               Numerical_Constants_Physical
    use               Numerical_Constants_Units
    use               Numerical_Constants_Atomic
    implicit none
    integer  (kind=c_int                  )                       :: backgroundRadiationODEs
    real     (kind=c_double               )               , value :: time
    real     (kind=c_double               ), intent(in   )        :: spectrum            (backgroundRadiationWavelengthCount)
    real     (kind=c_double               ), intent(  out)        :: spectrumRateOfChange(backgroundRadiationWavelengthCount)
    real     (kind=c_double               )                       :: spectralGradient    (backgroundRadiationWavelengthCount)
    type     (c_ptr                       )               , value :: parameterPointer

    ! Add source terms, linearly interpolating between timesteps. Convert from emissivity units [L☉ Hz⁻¹ Mpc⁻³] to
    ! background units [photons m⁻³ Hz⁻¹].
    spectrumRateOfChange                                &
         & =(                                           &
         &   +                     emissivity    (:,0)  &
         &   +(emissivity    (:,1)-emissivity    (:,0)) &
         &   *(time               -emissivityTime(  0)) &
         &   /(emissivityTime(  1)-emissivityTime(  0)) &
         &  )                                           &
         & *backgroundRadiationWavelength               &
         & *luminositySolar                             &
         & *gigaYear                                    &
         & /angstromsPerMeter                           &
         & /plancksConstant                             &
         & /speedLight                                  &
         & /megaParsec**3
    ! Add expansion dilution: -3 H(t)      n_nu
    spectrumRateOfChange=spectrumRateOfChange-3.0d0*cosmologyFunctions_%expansionRate(time)*spectrum
    ! Add redshifting       : +  H(t) d(nu n_nu)/dnu
    spectralGradient(1                                   )=-spectrum(1)
    spectralGradient (2:backgroundRadiationWavelengthCount)                                                                                                                                                                                             &
         & =-                                               backgroundRadiationWavelength(2:backgroundRadiationWavelengthCount)**2                                                                                                                      &
         & *(spectrum(2:backgroundRadiationWavelengthCount)/backgroundRadiationWavelength(2:backgroundRadiationWavelengthCount)-spectrum(1:backgroundRadiationWavelengthCount-1)/backgroundRadiationWavelength(1:backgroundRadiationWavelengthCount-1)) &
         & /(                                               backgroundRadiationWavelength(2:backgroundRadiationWavelengthCount)-                                                 backgroundRadiationWavelength(1:backgroundRadiationWavelengthCount-1))
    spectrumRateOfChange=spectrumRateOfChange+spectralGradient*cosmologyFunctions_%expansionRate(time)
    ! Absorption.
    spectrumRateOfChange                                                        &
         & =+spectrumRateOfChange                                               &
         &  -gigaYear                                                           &
         &  *massSolar                                                          &
         &  /megaParsec                                                     **3 &
         &  *centi                                                          **2 &
         &  *speedLight                                                         &
         &  *(                                                                  &
         &    +crossSectionNeutralHydrogen                                      &
         &    *hydrogenByMassPrimordial                                         &
         &    *intergalacticMediumState_  %neutralHydrogenFraction    (time)    &
         &    +crossSectionNeutralHelium                                        &
         &    *heliumByMassPrimordial                                           &
         &    *intergalacticMediumState_  %neutralHeliumFraction      (time)    &
         &    +crossSectionSinglyIonizedHelium                                  &
         &    *heliumByMassPrimordial                                           &
         &    *intergalacticMediumState_  %singlyIonizedHeliumFraction(time)    &
         &   )                                                                  &
         &  /atomicMassHydrogen                                                 &
         &  *cosmologyParameters_       %OmegaBaryon                  (    )    &
         &  *cosmologyParameters_       %densityCritical              (    )    &
         &  /cosmologyFunctions_        %expansionFactor              (time)**3 &
         &  /atomicMassUnit                                                     &
         &  *spectrum
    ! Return success.
    backgroundRadiationODEs=FGSL_Success
    return
  end function backgroundRadiationODEs

end module Radiation_Intergalactic_Background_Internal
