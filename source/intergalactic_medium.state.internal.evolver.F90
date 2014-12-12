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

!+ Contributions to this file made by: Daniel McAndrew.

 !% Contains a module which evolves the temperature and ionization states of baryons in the \gls{igm}.
 
 module Intergalactic_Medium_State_Internal_Evolver
   !% Evolves the temperature and ionization states of baryons in the \gls{igm}.
   use Cosmology_Functions
   use Cosmology_Parameters
   use Intergalactic_Medium_State
   use Abundances_Structure
   use FGSL
   private
   public :: Intergalactic_Medium_State_Internal_Initialize
   
   logical                                                                        :: igmPropertiesCompute
   double precision                                 , allocatable, dimension(:  ) :: temperature                                , massFiltering               , &
        &                                                                            clumpingFactor                             , opticalDepth
   integer                                                                        :: igmPropertiesTimeCountPerDecade            , timeCount
   double precision                                                               :: igmPropertiesRedshiftMinimum               , igmPropertiesRedshiftMaximum
   double precision                                                               :: timeMinimum                                , timeMaximum
   double precision                                 , allocatable, dimension(:  ) :: time                                       , redshift
   double precision                                 , allocatable, dimension(:,:) :: densityHydrogen                            , densityHelium               , &
        &                                                                            massFilteringComposite

   ! Classes used in ODE solution.
   class           (cosmologyParametersClass       ), pointer                     :: cosmologyParameters_
   class           (cosmologyFunctionsClass        ), pointer                     :: cosmologyFunctions_
   type            (intergalacticMediumStateRecFast)                              :: igmInitialState
   class           (intergalacticMediumStateClass  ), pointer                     :: igmState_  
 
 contains
   
   !# <universePreEvolveTask>
   !#  <unitName>Intergalactic_Medium_State_Internal_Initialize</unitName>
   !# </universePreEvolveTask>
   subroutine Intergalactic_Medium_State_Internal_Initialize(thisUniverse)
     !% Attach an initial event to a merger tree to cause the properties update function to be called.
     use Galacticus_Nodes
     use Input_Parameters
     use Memory_Management
     use Numerical_Ranges
     use Cosmology_Functions
     use Cosmology_Parameters
     use Linear_Growth
     use Power_Spectra
     use Numerical_Constants_Math
     use Numerical_Constants_Units
     use Numerical_Constants_Atomic
     use Numerical_Constants_Astronomical
     use Numerical_Comparison
     use Galacticus_Output_Times
     use Galacticus_Error
     implicit none
     type            (universe     ), intent(inout) :: thisUniverse
     type            (universeEvent), pointer       :: thisEvent
     integer                                        :: iTime                       , atomicNumber             , &
          &                                            ionizationState
     double precision                               :: massFilteringInitial        , massFilteringCoefficient1, &
          &                                            massFilteringVarianceInitial, massFilteringCoefficient2, &
          &                                            massFilteringCoefficient3   , massFilteringCoefficient4, &
          &                                            wavenumberFilteringInitial  , expansionFactorInitial   , &
          &                                            expansionRateInitial        , atomicMass               , &
          &                                            density                     , ionicFraction            , &
          &                                            massFractionPrimordial
     
     ! Get parameter controlling properties.
     !@ <inputParameter>
     !@   <name>igmPropertiesCompute</name>
     !@   <defaultValue>false</defaultValue>
     !@   <attachedTo>module</attachedTo>
     !@   <description>
     !@     Specifies whether or not the properties should be computed.
     !@   </description>
     !@   <type>integer</type>
     !@   <cardinality>1</cardinality>
     !@ </inputParameter>
     call Get_Input_Parameter('igmPropertiesCompute', igmPropertiesCompute,defaultValue=.false.)
     if (igmPropertiesCompute) then
        !@ <inputParameter>
        !@   <name>igmPropertiesTimeCountPerDecade</name>
        !@   <defaultValue>10</defaultValue>
        !@   <attachedTo>module</attachedTo>
        !@   <description>
        !@     The number of bins per decade of time to use for calculations of the properties of the universe.
        !@   </description>
        !@   <type>integer</type>
        !@   <cardinality>1</cardinality>
        !@ </inputParameter>
        call Get_Input_Parameter('igmPropertiesTimeCountPerDecade', igmPropertiesTimeCountPerDecade,defaultValue=10)
        !@ <inputParameter>
        !@   <name>igmPropertiesRedshiftMinimum</name>
        !@   <defaultValue>0</defaultValue>
        !@   <attachedTo>module</attachedTo>
        !@   <description>
        !@     The minimum redshift to use in calculations.
        !@   </description>
        !@   <type>real</type>
        !@   <cardinality>1</cardinality>
        !@ </inputParameter>
        call Get_Input_Parameter('igmPropertiesRedshiftMinimum', igmPropertiesRedshiftMinimum,defaultValue=0.0d0)
        !@ <inputParameter>
        !@   <name>igmPropertiesRedshiftMaximum</name>
        !@   <defaultValue>400</defaultValue>
        !@   <attachedTo>module</attachedTo>
        !@   <description>
        !@     The maximum redshift to use in calculations.
        !@   </description>
        !@   <type>real</type>
        !@   <cardinality>1</cardinality>
        !@ </inputParameter>
        call Get_Input_Parameter('igmPropertiesRedshiftMaximum', igmPropertiesRedshiftMaximum,defaultValue=400.0d0)
        ! Build tables of properties and time for the temperature and ionization state densities in the IGM.
        cosmologyParameters_ => cosmologyParameters            (                    )
        cosmologyFunctions_  => cosmologyFunctions             (                    )
        igmInitialState      =  intergalacticMediumStateRecFast(cosmologyParameters_)
        timeMaximum                                                                                 &
             & =cosmologyFunctions_%cosmictime                 (                                    &
             &  cosmologyFunctions_%expansionFactorFromredshift (                                   &
             &                                                   igmPropertiesRedshiftMinimum       &
             &                                                  )                                   &
             &                                                 )
        timeMinimum                                                                                 &
             & =cosmologyFunctions_%cosmictime                 (                                    &
             &  cosmologyFunctions_%expansionFactorFromredshift (                                   &
             &                                                   igmPropertiesRedshiftMaximum       &
             &                                                  )                                   &
             &                                                 )
        timeCount                                                       &
             & = int(                                                   &
             &        dble(igmPropertiesTimeCountPerDecade      )       &
             &       *log10(                                            &
             &               timeMaximum                                &
             &              /timeMinimum                                &
             &             )                                            &
             &      )                                                   &
             &  +1
        ! Allocate arrays for all required IGM properties.
        call Alloc_Array(temperature           ,[timeCount  ])
        call Alloc_Array(densityHydrogen       ,[timeCount,2])
        call Alloc_Array(densityHelium         ,[timeCount,3])
        call Alloc_Array(time                  ,[timeCount  ])
        call Alloc_Array(redshift              ,[timeCount  ])
        call Alloc_Array(massFilteringComposite,[timeCount,2]) 
        call Alloc_Array(clumpingFactor        ,[timeCount  ])
        call Alloc_Array(opticalDepth          ,[timeCount  ])
        call Alloc_Array(massFiltering         ,[timeCount  ])
        ! Build grid of times.
        time                                                     &
             & =Make_Range(                                      &
             &             timeMinimum                         , &
             &             timeMaximum                         , &
             &             timeCount                           , &
             &             rangeTypeLogarithmic                  &
             &            )
        ! Convert times to redshifts.
        do iTime=1,timeCount
           if (Values_Agree(time(iTime),Galacticus_Output_time(Galacticus_Output_time_Count()),relTol=1.0d-6)) &
                & time(iTime)=Galacticus_Output_Time(Galacticus_Output_time_Count())
           redshift(itime)                                                       &
                & =cosmologyFunctions_ %redshiftFromExpansionFactor(             &
                &   cosmologyFunctions_%expansionFactor             (            &
                &                                                    time(itime) &
                &                                                   )            &
                &                                                  )
        end do
        ! Initialize arrays to unphysical values.
        temperature              =-1.0d0
        densityHydrogen          =-1.0d0
        densityHelium            =-1.0d0
        massFilteringComposite   =-1.0d0
        clumpingFactor           =-1.0d0
        opticalDepth             =-1.0d0
        massFiltering            =-1.0d0 
        ! Initialize the temperature to that computed by RecFast at the initial time.
        temperature(1)=igmInitialState%temperature(timeMinimum)
        ! Initialize the densities of ionic species using fractions from RecFast at the initial time.
        do atomicNumber=1,2
           select case (atomicNumber)
           case (1)
              massFractionPrimordial=hydrogenByMassPrimordial
              atomicMass            =atomicMassHydrogen
           case (2)
              massFractionPrimordial=heliumByMassPrimordial
              atomicMass            =atomicMassHelium
           end select
           do ionizationState=1,atomicNumber+1
              select case (atomicNumber)
              case (1)
                 select case (ionizationState)
                 case (1)
                    ionicFraction=igmInitialState%neutralHydrogenFraction      (timeMinimum)
                 case (2)
                    ionicFraction=igmInitialState%singlyIonizedHydrogenFraction(timeMinimum)
                 end select                 
              case (2)
                 select case (ionizationState)
                 case (1)
                    ionicFraction=igmInitialState%neutralHeliumFraction        (timeMinimum)
                 case (2)
                    ionicFraction=igmInitialState%singlyIonizedHeliumFraction  (timeMinimum)
                 case (3)
                    ionicFraction=igmInitialState%doublyIonizedHeliumFraction  (timeMinimum)
                 end select                 
              case default
                 call Galacticus_Error_Report('Intergalactic_Medium_State_Internal_Initialize','unknown atomic number')
              end select             
              density=+massFractionPrimordial                               &
                   &  *cosmologyParameters_%OmegaBaryon    (           )    &
                   &  *cosmologyParameters_%densityCritical(           )    &
                   &  /cosmologyFunctions_ %expansionFactor(timeMinimum)**3 &
                   &  /atomicMassUnit                                       &
                   &  /atomicMass                                           &
                   &  *massSolar                                            &
                   &  /megaParsec                                       **3 &
                   &  *ionicFraction
              select case (atomicNumber)
              case (1)
                 densityHydrogen(1,ionizationState)=density
              case (2)
                 densityHelium  (1,ionizationState)=density
              end select
           end do
        end do
        ! Initialize the filtering mass using the fitting functions from Naoz & Barkana (2007).
        if     (                                                                                                                       &
             &   igmPropertiesRedshiftMaximum > 150.0d0                                                                                &
             &  .or.                                                                                                                   &
             &   igmPropertiesRedshiftMaximum <    7.0d0                                                                               &
             & ) call Galacticus_Error_Report(                                                                                         &
             &                                'Intergalactic_Medium_State_Internal_Initialize'                                       , &
             &                                'initial redshift is out of range for Naoz & Barkana (2007) fitting function (7≤z≤150)'  &
             &                               )
        expansionFactorInitial   =cosmologyFunctions_%expansionFactor(                                    timeMinimum )
        expansionRateInitial     =cosmologyFunctions_%expansionRate  (cosmologyFunctions_%expansionFactor(timeMinimum))
         ! Evaluate fitting coefficients for filtering mass.
        massFilteringCoefficient1=-0.38d0*(cosmologyParameters_%OmegaMatter()**2)+ 0.41d0*cosmologyParameters_%OmegaMatter()- 0.16d0
        massFilteringCoefficient2=+3.30d0*(cosmologyParameters_%OmegaMatter()**2)- 3.38d0*cosmologyParameters_%OmegaMatter()+ 1.15d0
        massFilteringCoefficient3=-9.64d0*(cosmologyParameters_%OmegaMatter()**2)+ 9.75d0*cosmologyParameters_%OmegaMatter()- 2.37d0
        massFilteringCoefficient4=+9.80d0*(cosmologyParameters_%OmegaMatter()**2)-10.68d0*cosmologyParameters_%OmegaMatter()+11.60d0
        massFilteringInitial     =                                                                 &
             &                    exp(                                                             &
             &                        +massFilteringCoefficient1*(-log(expansionFactorInitial))**3 &
             &                        +massFilteringCoefficient2*(-log(expansionFactorInitial))**2 &
             &                        +massFilteringCoefficient3*(-log(expansionFactorInitial))    &
             &                        +massFilteringCoefficient4                                   &
             &                       )
        massFiltering(1)         =massFilteringInitial
        ! Find the corresponding filtering wavenumber.
        wavenumberFilteringInitial=+Pi                                       &
             &                     /(                                        &
             &                       +massFilteringInitial                   &
             &                       *3.0d0                                  &
             &                       /4.0d0                                  &
             &                       /Pi                                     &
             &                       /cosmologyParameters_%OmegaMatter    () &
             &                       /cosmologyParameters_%densityCritical() &
             &                      )**(1.0d0/3.0d0)
        ! Find the initial mass variance on the filtering mass scale.
        massFilteringVarianceInitial=Cosmological_Mass_Root_Variance(massFilteringInitial)        
        ! Set the initial clumping factor.
        clumpingFactor(1)= 1.0d0+(massFilteringVarianceInitial**2*Linear_Growth_Factor(timeMinimum)**2)
        ! Set the composite variables used to solve for filtering mass.
        massFilteringComposite(1,1)=                                       &
             &   +Linear_Growth_Factor                       (timeMinimum) &
             &   /wavenumberFilteringInitial**2
        massFilteringComposite(1,2)=                                       &
             &   +Linear_Growth_Factor                       (timeMinimum) &
             &   /timeMinimum                                              &
             &   *Linear_Growth_Factor_Logarithmic_Derivative(timeMinimum) &
             &   /wavenumberFilteringInitial**2                            &
             &   +2.0d0                                                    &
             &   /3.0d0                                                    &
             &   *Linear_Growth_Factor                       (timeMinimum) &
             &   /wavenumberFilteringInitial**2                            &
             &   *(                                                        &
             &     -3.0d0                                                  &
             &     *massFilteringCoefficient1                              &
             &     *log(expansionFactorInitial)**2                         &
             &     +2.0d0                                                  &
             &     *massFilteringCoefficient2                              &
             &     *log(expansionFactorInitial)                            &
             &     -massFilteringCoefficient3                              &
             &    )                                                        &
             &   *expansionRateInitial
        ! Initialize optical depth.
        opticalDepth(1)=0.0d0
        ! Initialize the IGM state object.
        igmState_ => intergalacticMediumState()
        select type(igmState_)
        class is (intergalacticMediumStateInternal)
           call igmState_%timeSet         (time           (1:1  ))
           call igmState_%temperatureSet  (temperature    (1:1  ))
           call igmState_%densityH1Set    (densityHydrogen(1:1,1))
           call igmState_%densityH2Set    (densityHydrogen(1:1,2))
           call igmState_%densityHe1Set   (densityHelium  (1:1,1))
           call igmState_%densityHe2Set   (densityHelium  (1:1,2))
           call igmState_%densityHe3Set   (densityHelium  (1:1,3))
           call igmState_%massFilteringSet(massFiltering  (1:1  ))
        end select
        ! Create the first interrupt event in the universe object.
        thisEvent      => thisUniverse%createEvent()
        thisEvent%time =  time(1)
        thisEvent%task => Intergalactic_Medium_State_Internal_Update
     end if
     return
   end subroutine Intergalactic_Medium_State_Internal_Initialize
   
   logical function Intergalactic_Medium_State_Internal_Update(thisEvent,thisUniverse) result (success)
     !% Update the properties for a given universe.
     use               Galacticus_Nodes
     use               Galacticus_Display
     use               Galacticus_Output_Times
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
     use               Numerical_Constants_Astronomical
     use               Numerical_Integration
     use               Cosmology_Functions
     use               Cosmology_Parameters
     use               Linear_Growth
     use               Power_Spectra
     use               ISO_Varying_String
     use               Galacticus_HDF5
     use               IO_HDF5
     implicit none
     class           (universeEvent                ), intent(in   )            :: thisEvent
     type            (universe                     ), intent(inout)            :: thisUniverse
     type            (universeEvent                ), pointer                  :: newEvent
     class           (intergalacticMediumStateClass), pointer                  :: igmState_
     double precision                               , parameter                :: odeToleranceAbsolute        =1.0d-03,             &
          &                                                                       odeToleranceRelative        =1.0d-03,             &
          &                                                                       timeToleranceRelative       =1.0d-6
     type            (mergerTree                   ), pointer                  :: thisTree       
     type            (mergerTreeList               ), pointer                  :: thisForest       
     type            (treeNode                     ), pointer                  :: thisNode
     class           (nodeComponentBasic           ), pointer                  :: thisBasic
     type            (fodeiv2_system               ), save                     :: ode2System
     type            (fodeiv2_driver               ), save                     :: ode2Driver
     logical                                        , save                     :: odeReset
     integer                                        , parameter                :: propertyCount               =10
     double precision                               , parameter                :: massFilteringScale=1.0d4
     double precision                               , dimension(propertyCount) :: properties, propertyScales
     type            (c_ptr                        )                           :: parameterPointer
     type            (varying_string               )                           :: message
     character       (len=6                        )                           :: label
     type            (hdf5Object                   )                           :: igmGroup                            , igmDataset
     integer         (c_size_t                     )                           :: iNow
     double precision                                                          :: treetimeLatest                      , timeCurrent, &
          &                                                                       timeMaximum

     ! Display message.
     write (label,'(f6.3)') thisEvent%time
     message = "Evolving IGM properties to time "//trim(label)//" Gyr"
     call Galacticus_Display_Indent(message)
     ! Find the current timestep.
     iNow = Search_Array_For_Closest(time,thisEvent%time)
     ! Evolve the properties up to this timestep.
     if (iNow > 1) then        
        ! Get required objects.
        cosmologyParameters_ => cosmologyParameters()
        cosmologyFunctions_  => cosmologyFunctions ()
        ! Map properties to a contiguous array.
        properties( 1   )=temperature           (iNow-1    )
        properties( 2: 3)=densityHydrogen       (iNow-1,1:2)
        properties( 4: 6)=densityHelium         (iNow-1,1:3)
        properties( 7: 8)=massFilteringComposite(iNow-1,1:2)
        properties( 9   )=opticalDepth          (iNow-1    )
        properties(10   )=massFiltering         (iNow-1    )
        ! Set property scales.
        propertyScales( 1  )=        temperature    (iNow-1  )
        propertyScales( 2:3)=abs(sum(densityHydrogen(iNow-1,:)))
        propertyScales( 4:6)=abs(sum(densityHelium  (iNow-1,:)))
        propertyScales( 7  )=+Linear_Growth_Factor(time(iNow-1))         &
             &               *(                                          &
             &                 +massFiltering(iNow-1)                    &
             &                 *3.0d0                                    &
             &                 /(                                        &
             &                   +4.0d0                                  &
             &                   *Pi                                     &
             &                   *cosmologyParameters_%OmegaMatter    () &
             &                   *cosmologyParameters_%densityCritical() &
             &                  )                                        &
             &                )**(2.0d0/3.0d0)                           &
             &               *Pi**2
        propertyScales( 8  )=+propertyScales(7)&
             &               *cosmologyParameters_%HubbleConstant(unitstime)
        propertyScales( 9  )=opticalDepth (iNow-1)
        propertyScales(10  )=massFiltering(iNow-1)
        ! Display message
        call Galacticus_Display_Message('Solving properties evolution')
        odeReset=.true.
        timeCurrent=time(iNow-1)
        call ODEIV2_Solve(                                                               &
             &            ode2Driver                                                   , &
             &            ode2System                                                   , &
             &            timeCurrent                                                  , &
             &            time                                    (iNow)               , &
             &            propertyCount                                                , &
             &            properties                                                   , &
             &            Intergalactic_Medium_State_Internal_ODEs                     , &
             &            parameterPointer                                             , &
             &            odeToleranceAbsolute                                         , &
             &            odeToleranceRelative                                         , &
             &            yScale                                        =propertyScales, &
             &            reset                                         =odeReset        &
             &           )
        call ODEIV2_Solver_Free(ode2Driver,ode2System)
        temperature           (iNow    )=max(properties( 1   ),0.0d0)
        densityHydrogen       (iNow,1:2)=max(properties( 2: 3),0.0d0)
        densityHelium         (iNow,1:3)=max(properties( 4: 6),0.0d0)
        massFilteringComposite(iNow,1:2)=    properties( 7: 8)
        opticalDepth          (iNow    )=max(properties( 9   ),0.0d0)
        massFiltering         (iNow    )=properties(10   )
        ! Compute the filtering mass at this time.
        clumpingFactor        (iNow    )=+1.0d0                                                   &
             &                           +Cosmological_Mass_Root_Variance(massFiltering(iNow))**2 &
             &                           *Linear_Growth_Factor           (time         (iNow))**2
    end if
    ! Find the latest time across all trees in the universe.
    treetimeLatest=0.0d0
    thisForest => thisUniverse%trees
    do while (associated(thisForest))
       thisTree => thisForest%tree
       do while (associated(thisTree))
          thisNode       => thisTree%baseNode
          thisBasic      => thisNode%basic()
          treetimeLatest =  max(treetimeLatest,thisBasic%time())
          thisTree       => thisTree%nextTree
       end do
       thisForest => thisForest%next
    end do
    ! Add the next event to the universe.
    timeMaximum=min(treetimeLatest,Galacticus_Output_time(Galacticus_Output_time_Count()))
    if (iNow < timeCount .and. time(iNow+1) <= timeMaximum) then
       newEvent      => thisUniverse%createEvent()
       newEvent%time =  min(time(iNow+1),timeMaximum*(1.0d0-timeToleranceRelative))
       newEvent%task => Intergalactic_Medium_State_Internal_Update
    else
       ! Output the results to file.
       !$omp critical (HDF5_Access)
       igmGroup=galacticusOutputFile%openGroup('igmProperties', 'Properties of the intergalactic medium.')
       call igmGroup  %writeDataset  (redshift            ,'redshift'        ,'Redshift [].'                                  ,datasetReturned=igmDataset)
       call igmDataset%writeAttribute(0.0d0               ,'unitsInSI'                                                                                   )
       call igmDataset%close()
       call igmGroup  %writeDataset  (temperature         ,'temperature'     ,'Temperature of the IGM [K].'                   ,datasetReturned=igmDataset)
       call igmDataset%writeAttribute(1.0d0               ,'unitsInSI'                                                                                   )
       call igmDataset%close()
       call igmGroup  %writeDataset  (densityHydrogen(:,1),'densityHydrogen1','Density of H1 in the IGM [m⁻³].'               ,datasetReturned=igmDataset)
       call igmDataset%writeAttribute(1.0d0               ,'unitsInSI'                                                                                   )
       call igmDataset%close()
       call igmGroup  %writeDataset  (densityHydrogen(:,2),'densityHydrogen2','Density of H2 in the IGM [m⁻³].'               ,datasetReturned=igmDataset)
       call igmDataset%writeAttribute(1.0d0               ,'unitsInSI'                                                                                   )
       call igmDataset%close()
       call igmGroup  %writeDataset  (densityHelium  (:,1),'densityHelium1'  ,'Density of He1 in the IGM [m⁻³].'              ,datasetReturned=igmDataset)
       call igmDataset%writeAttribute(1.0d0               ,'unitsInSI'                                                                                   )
       call igmDataset%close()  
       call igmGroup  %writeDataset  (densityHelium  (:,2),'densityHelium2'  ,'Density of He2 in the IGM [m⁻³].'              ,datasetReturned=igmDataset)
       call igmDataset%writeAttribute(1.0d0               ,'unitsInSI'                                                                                   )
       call igmDataset%close() 
       call igmGroup  %writeDataset  (densityHelium  (:,3),'densityHelium3'  ,'Density of He3 in the IGM [m⁻³].'              ,datasetReturned=igmDataset)
       call igmDataset%writeAttribute(1.0d0               ,'unitsInSI'                                                                                   )
       call igmDataset%close()  
       call igmGroup  %writeDataset  (clumpingFactor      ,'clumpingFactor'  ,'Clumping factor in the IGM [].'                ,datasetReturned=igmDataset)
       call igmDataset%writeAttribute(1.0d0               ,'unitsInSI'                                                                                   )
       call igmDataset%close()
       call igmGroup  %writeDataset  (opticalDepth        ,'opticalDepth'    ,'Electron scattering optical depth from z=0 [].',datasetReturned=igmDataset)
       call igmDataset%writeAttribute(1.0d0               ,'unitsInSI'                                                                                   )
       call igmDataset%close()
       call igmGroup  %writeDataset  (massFiltering       ,'filteringMass'   ,'Filtering mass in the IGM [M☉].'               ,datasetReturned=igmDataset)
       call igmDataset%writeAttribute(1.0d0               ,'unitsInSI'                                                                                   )
       call igmDataset%close()
       call igmGroup  %close()
       !$omp end critical (HDF5_Access)
     end if
     ! Store the past history to the default IGM state class.
     igmState_ => intergalacticMediumState()
     select type (igmState_)
     class is (intergalacticMediumStateInternal)
        call igmState_%timeSet         (time           (1:iNow  ))
        call igmState_%temperatureSet  (temperature    (1:iNow  ))
        call igmState_%densityH1Set    (densityHydrogen(1:iNow,1))
        call igmState_%densityH2Set    (densityHydrogen(1:iNow,2))
        call igmState_%densityHe1Set   (densityHelium  (1:iNow,1))
        call igmState_%densityHe2Set   (densityHelium  (1:iNow,2))
        call igmState_%densityHe3Set   (densityHelium  (1:iNow,3))
        call igmState_%massFilteringSet(massFiltering  (1:iNow  ))
     class default
        call Galacticus_Error_Report('Intergalactic_Medium_State_Internal_Update','"internal" IGM evolution calculation requires [intergalacticMediumStateMethod]=internal')
     end select
     ! Display message.
     call Galacticus_Display_Unindent('done')
     ! Return true since we've performed our task.
     success=.true.
     return
   end function Intergalactic_Medium_State_Internal_Update
   
   function Intergalactic_Medium_State_Internal_ODEs(time,properties,propertiesRateOfChange,parameterPointer) bind(c)
     !% Evaluates the ODEs controlling the evolution temperature.
     use               ODE_Solver_Error_Codes
     use, intrinsic :: ISO_C_Binding
     use               Numerical_Constants_Astronomical
     use               Numerical_Constants_Math
     use               Numerical_Constants_Physical
     use               Numerical_Constants_Units
     use               Numerical_Constants_Prefixes
     use               Numerical_Constants_Atomic
     use               Cosmology_Parameters
     use               Cosmology_Functions
     use               Atomic_Ionization_Potentials
     use               Atomic_Rates_Ionization_Collisional
     use               Atomic_Rates_Recombination_Radiative
     use               Atomic_Rates_Recombination_Dielectronic
     use               atomic_radiation_gaunt_factors
     use               Atomic_Rates_Excitation_Collisional
     use               Atomic_Rates_Recombination_Radiative_Data
     use               Radiation_Structure
     use               Atomic_Cross_Sections_Ionization_Photo
     use               Linear_Growth
     use               Power_Spectra
     use               Numerical_Integration
     use               FGSL
     implicit none
     integer         (kind=c_int                )                               :: Intergalactic_Medium_State_Internal_ODEs
     real            (kind=c_double             ), value                        :: time
     real            (kind=c_double             ), intent(in   ), dimension(10) :: properties
     real            (kind=c_double             ), intent(  out), dimension(10) :: propertiesRateOfChange
     type            (c_ptr                     ), value                        :: parameterPointer
     double precision                            , parameter                    :: dielectronicRecombinationRateHeIEnergyLoss=40.74d0 ! electron volts.
     double precision                            , parameter                    :: massFilteringMinimum                      =1.0d2
     double precision                                           , dimension( 2) :: densityHydrogen_                                  , massFilteringComposite_            , &
          &                                                                        massFilteringCompositeRateOfChange
     double precision                                           , dimension( 3) :: densityHelium_
     type            (radiationStructure        ), save                         :: radiation  
     type            (fgsl_function             )                               :: integrationFunction
     type            (fgsl_integration_workspace)                               :: integrationWorkspace
     logical                                                                    :: integrationReset
     integer                                                                    :: electronNumber                                    , atomicNumber                       , &
          &                                                                        ionizationState                                   , shellNumber                        , &
          &                                                                        photoionizationGroundIonizationState              , photoionizationGroundElectronNumber, &
          &                                                                        iProperty
     double precision                                                           :: temperature                                       , clumpingFactor                     , &
          &                                                                        electronDensityRateOfChange                       , densityElectron                    , &
          &                                                                        densityTotal                                      , ionizationPhotoRateFrom            , &
          &                                                                        ionizationPhotoRateTo                             , opticalDepthRateOfChange           , &
          &                                                                        massFilteringRateOfChange                         , wavenumberFilteringRateOfChange    , &
          &                                                                        collisionIonizationRateFrom                       , collisionIonizationRateTo          , &
          &                                                                        densityLowerIon                                   , densityUpperIon                    , &
          &                                                                        densityThisIon                                    , recombinationDielectronicRateTo    , &
          &                                                                        recombinationDielectronicRateFrom                 , recombinationRateTo                , &
          &                                                                        recombinationRateFrom                             , wavelengthMinimum                  , &
          &                                                                        wavelengthMaximum                                 , darkMatterFraction                 , &
          &                                                                        massParticleMean                                  , massFiltering_                     , &
          &                                                                        wavenumberFiltering                               , opticalDepth                       , &
          &                                                                        heatingRate

     ! Extract properties from the contiguous array.
     temperature            =max(properties( 1   ),0.0d0               )
     densityHydrogen_       =    properties( 2: 3)
     densityHelium_         =    properties( 4: 6)
     massFilteringComposite_=    properties( 7: 8)
     opticalDepth           =    properties( 9   )
     massFiltering_         =max(properties(10   ),massFilteringMinimum)
     ! Compute total and electron densities.
     densityElectron        =+    max(      densityHydrogen_(2),0.0d0)  &
          &                  +    max(      densityHelium_  (2),0.0d0)  &
          &                  +    max(2.0d0*densityHelium_  (3),0.0d0)
     densityTotal           =+sum(max(      densityHydrogen_   ,0.0d0)) &
          &                  +sum(max(      densityHelium_     ,0.0d0)) &
          &                  +densityElectron
     ! Get required objects.
     cosmologyParameters_ => cosmologyParameters()
     cosmologyFunctions_  => cosmologyFunctions ()
     ! Initialize heating rate to zero.
     heatingRate          =  0.0d0
     ! Compute dark matter mass fraction.
     darkMatterFraction   =  1.0d0-cosmologyParameters_%OmegaBaryon()/cosmologyParameters_%OmegaMatter()
     ! Compute rates of change of filtering mass composite parameters and optical depth.
     if (densityTotal > 0.0d0) then
        ! Find mean particle mass.
        massParticleMean=+(                                        &
             &             +massHydrogenAtom*sum(densityHydrogen_) &
             &             +massHeliumAtom  *sum(densityHelium_  ) &
             &             +electronMass    *    densityElectron   &
             &            )                                        &
             &           /densityTotal
        ! Evaluate filtering mass composite terms.
        massFilteringCompositeRateOfChange(1)=massFilteringComposite_(2)
        massFilteringCompositeRateOfChange(2)=-2.0d0                                                                                &
             &                                *cosmologyFunctions_%expansionRate(cosmologyFunctions_%expansionFactor(time))         &
             &                                *massFilteringComposite_(2)                                                           &
             &                                +darkMatterFraction                                                                   &
             &                                /                                  cosmologyFunctions_%expansionFactor(time)**2       &
             &                                *boltzmannsConstant                                                                   &
             &                                *temperature                                                                          & 
             &                                /massParticleMean                                                                     &
             &                                *Linear_Growth_Factor                                                 (time)          &
             &                                *(                                                                                    &
             &                                  +1.0d0                                                                              &
             &                                  +rLSS(                                                                              &
             &                                                                   cosmologyParameters_%OmegaMatter   (    )        , &
             &                                                                   cosmologyFunctions_%expansionFactor(time)          &
             &                                       )                                                                              &
             &                                )                                                                                     &
             &                                *gigayear  **2                                                                        &
             &                                /megaparsec**2
        ! Evaluate optical depth term.
        opticalDepthRateOfChange=+speedLight                                                                   &
             &                   *thomsonCrossSection                                                          &
             &                   *densityElectron                                                              &
             &                   /cosmologyFunctions_%expansionRate(cosmologyFunctions_%expansionFactor(time)) &
             &                   *                                  cosmologyFunctions_%expansionFactor(time)  &
             &                   *gigayear
     else
        massFilteringCompositeRateOfChange(1)=massFilteringComposite_(2)
        massFilteringCompositeRateOfChange(2)=0.0d0
        opticalDepthRateOfChange             =0.0d0
     end if
     if (massFilteringComposite_(1) > 0.0d0) then
        wavenumberFiltering            =+sqrt(                            &
             &                                +Linear_Growth_Factor(time) &
             &                                /massFilteringComposite_(1) &
             &                               )
        wavenumberFilteringRateOfChange=+0.5d0                                                                                    &
             &                          /sqrt(                                                                                    &
             &                                +Linear_Growth_Factor                                                       (time)  &
             &                                /massFilteringComposite_                                                    (1   )  &
             &                               )                                                                                    &
             &                          *(                                                                                        &
             &                            +Linear_Growth_Factor_Logarithmic_Derivative                                    (time)  &
             &                            *cosmologyFunctions_%expansionRate          (cosmologyFunctions_%expansionFactor(time)) &
             &                            *Linear_Growth_Factor                                                           (time)  &
             &                            *massFilteringComposite_                                                        (1   )  &
             &                            -Linear_Growth_Factor                                                           (time)  &
             &                            *massFilteringComposite_                                                        (2   )  &
             &                           )                                                                                        &
             &                          /massFilteringComposite_(1)**2
        massFilteringRateOfChange      =-4.0d0                                     &
             &                          *Pi                                    **4 &
             &                          *cosmologyParameters_%densityCritical()    &
             &                          *cosmologyParameters_%OmegaMatter    ()    & 
             &                          /wavenumberFiltering                   **4 &
             &                          *wavenumberFilteringRateOfChange
     else
       massFilteringRateOfChange=0.0d0
     end if
     ! Compute the clumping factor.
     clumpingFactor=+1.0d0                                             &
          &         +(                                                 &
          &           +Cosmological_Mass_Root_Variance(massFiltering_) &
          &           *Linear_Growth_Factor           (time          ) &
          &          )**2
     ! Iterate over ionic species.
     iProperty=1 ! Counter for ionic species in properties array.
     do atomicNumber=1,2
        do ionizationState=1,atomicNumber+1
           ! Increment property array counter.
           iProperty=iProperty+1
           ! Determine electron number.
           electronNumber=atomicNumber+1-ionizationState
           ! Get density of this ionic state.
           densityThisIon    =max(properties(iProperty  ),0.0d0)
           ! Get density of upper ionic state (i.e. current ion minus one electron).
           if (ionizationState < atomicNumber+1) then
              densityUpperIon=max(properties(iProperty+1),0.0d0)
           else
              densityUpperIon=                            0.0d0
           end if
           ! Get density of upper ionic state (i.e. current ion plus one electron).
           if (ionizationState > 1) then
              densityLowerIon=max(properties(iProperty-1),0.0d0)
           else
              densityLowerIon=                            0.0d0
           end if
           ! Specify electron shell number. For H and He, this is always the 1s shell.
           shellNumber=1          
           ! Compute collisional ionization rates from this ion.
           if (electronNumber  > 0) then
              collisionIonizationRateFrom=-Atomic_Rate_Ionization_Collisional(atomicNumber,ionizationState  ,temperature                   ) &
                   &                      *densityThisIon                                                                                    &
                   &                      *densityElectron
              heatingRate                =+heatingRate                                                                                       &
                   &                      -Atomic_Ionization_Potential       (atomicNumber,electronNumber +1                               ) &
                   &                      *electronVolt                                                                                      &
                   &                      *collisionIonizationRateFrom                                                                       &
                   &                      *gigaYear                                                                                          &
                   &                      *centi**3                                                                                          &
                   &                      *clumpingFactor
            else
              collisionIonizationRateFrom=+0.0d0
           end if
           ! Compute collisional ionization rates to this ion.
           if (ionizationState > 1) then
              collisionIonizationRateTo  =+Atomic_Rate_Ionization_Collisional (atomicNumber,ionizationState-1,temperature                   ) &
                   &                      *densityLowerIon                                                                                    &
                   &                      *densityElectron
           else
              collisionIonizationRateTo  =+0.0d0  
           end if
           ! Compute recombination rates from this ion.
           if (ionizationState > 1) then
              recombinationRateFrom      =-Atomic_Rate_Recombination_Radiative(atomicNumber,ionizationState-1,temperature,recombinationCaseB) &
                   &                      *densityThisIon                                                                                     &
                   &                      *densityElectron
           else
              recombinationRateFrom      =+0.0d0
           end if
           ! Compute recombination rates to this ion.
           if (electronNumber  > 0) then
              recombinationRateTo        =+Atomic_Rate_Recombination_Radiative(atomicNumber,ionizationState  ,temperature,recombinationCaseB) &
                   &                      *densityUpperIon                                                                                    &
                   &                      *densityElectron
              heatingRate                =+heatingRate                                                                                        &
                   &                      -recombinationRateFrom                                                                              &
                   &                      *gigaYear                                                                                           &
                   &                      *centi**3                                                                                           &
                   &                      *clumpingFactor                                                                                     &
                   &                      *0.75d0                                                                                             &
                   &                      *boltzmannsConstant                                                                                 &
                   &                      *temperature
           else
              recombinationRateTo        =+0.0d0
           end if
           ! Compute dielectronic recombination rates from this ion.
           if (ionizationState > 1) then
              recombinationDielectronicRateFrom=-Dielectronic_Recombination_Rate(atomicNumber,electronNumber+1,temperature) &
                   &                            *densityThisIon                                                             &
                   &                            *densityElectron
           else
              recombinationDielectronicRateFrom=+0.0d0
           end if
           ! Compute dielectronic recombination rates to this ion.
           if (electronNumber  > 0) then
              recombinationDielectronicRateTo  =+Dielectronic_Recombination_Rate(atomicNumber,electronNumber  ,temperature) &
                   &                            *densityUpperIon                                                            &
                   &                            *densityElectron
              heatingRate                      =+heatingRate                                                   &
                   &                            -dielectronicRecombinationRateHeIEnergyLoss                    &
                   &                            *electronVolt                                                  &
                   &                            *recombinationDielectronicRateFrom                             &
                   &                            *gigaYear                                                      &
                   &                            *centi**3                                                      &
                   &                            *clumpingFactor
           else
              recombinationDielectronicRateTo  =+0.0d0
           end if
           ! Define the background radiation field.
           call radiation%define([radiationTypeIGB])
           call radiation%set   (             time )
           ! Compute rate of photoionizations from this ion.
           if (electronNumber  > 0) then
              ! Set the ground state for this photoionization calculation.
              photoionizationGroundIonizationState=ionizationState
              ! Set the minimum and maximum wavelengths for photoionization.
              wavelengthMinimum=+0.0d0
              wavelengthMaximum=+plancksConstant                                          &
                   &            *speedLight                                               &
                   &            /Atomic_Ionization_Potential(atomicNumber,electronNumber) &
                   &            /electronVolt                                             &
                   &            *angstromsPerMeter
              ! Integrate photoionizations over wavelength.
              integrationReset       =.true.
              ionizationPhotoRateFrom=-Integrate(                                                 &
                   &                             wavelengthMinimum                              , &
                   &                             wavelengthMaximum                              , &
                   &                             Photoionization_Rate_Integrand                 , &
                   &                             parameterPointer                               , &
                   &                             integrationFunction                            , &
                   &                             integrationWorkspace                           , &
                   &                             toleranceAbsolute             =0.0d+0          , &
                   &                             toleranceRelative             =1.0d-2          , &
                   &                             reset                         =integrationReset  &
                   &                            )                                                 &
                   &                  *densityThisIon
              call Integrate_Done(integrationFunction,integrationWorkspace) 
           else
              ionizationPhotoRateFrom =+0.0d0
           end if
           if (ionizationState > 1) then
              ! Set the ground state for this photoionization calculation.
              photoionizationGroundIonizationState=ionizationState-1
              photoionizationGroundElectronNumber =electronNumber +1
              ! Set the minimum and maximum wavelengths for photoionization.
              wavelengthMinimum=0.0d0
              wavelengthMaximum=+plancksConstant                                            &
                   &            *speedLight                                                 &
                   &            /Atomic_Ionization_Potential(atomicNumber,electronNumber+1) &
                   &            /electronVolt                                               &
                   &            *angstromsPerMeter
              ! Integrate photoionizations over wavelength.
              integrationReset       =.true.
              ionizationPhotoRateTo  =+Integrate(                                                 &
                   &                             wavelengthMinimum                              , &
                   &                             wavelengthMaximum                              , &
                   &                             Photoionization_Rate_Integrand                 , &
                   &                             parameterPointer                               , &
                   &                             integrationFunction                            , &
                   &                             integrationWorkspace                           , &
                   &                             toleranceAbsolute             =0.0d+0          , &
                   &                             toleranceRelative             =1.0d-2          , &
                   &                             reset                         =integrationReset  &
                   &                            )                                                 &
                   &                  *densityLowerIon
              call Integrate_Done(integrationFunction,integrationWorkspace) 
              integrationReset = .true.
              heatingRate                      =+heatingRate                                                        &
                   &                            +Integrate(                                                         &
                   &                                       wavelengthMinimum                                      , &
                   &                                       wavelengthMaximum                                      , &
                   &                                       Photoionization_Heating_Rate_Integrand                 , &
                   &                                       parameterPointer                                       , &
                   &                                       integrationFunction                                    , &
                   &                                       integrationWorkspace                                   , &
                   &                                       toleranceAbsolute                     =0.0d+0          , &
                   &                                       toleranceRelative                     =1.0d-3          , &
                   &                                       reset                                 =integrationReset  &
                   &                                      )                                                         &
                   &                            *densityLowerIon                                                    &
                   &                            *gigaYear
              call Integrate_Done(integrationFunction,integrationWorkspace)               
           else 
              ionizationPhotoRateTo   =+0.0d0
           end if
           ! Compute heating rate due to Bremsstrahlung.
           if (ionizationState > 1) then
              heatingRate=+heatingRate                    &
                   &      -16.0d0                         &
                   &      / 3.0d0                         &
                   &      *sqrt(                          &
                   &            +2.0d0                    &
                   &            *Pi                       &
                   &            /3.0d0                    &
                   &           )                          &
                   &      *dble(ionizationState-1) **2    &
                   &      *densityThisIon                 &
                   &      *densityElectron                &
                   &      *electronRadius          **3    &
                   &      *speedLight                     &
                   &      /electronRadius                 &
                   &      *sqrt(                          &
                   &            +electronMass             &
                   &            *speedLight        **2    &
                   &            *boltzmannsConstant       &
                   &            *temperature              &
                   &           )                          &
                   &      *fineStructure                  &
                   &      *Gaunt_Factor(                  &
                   &                    atomicNumber    , &
                   &                    electronNumber+1, &
                   &                    temperature       &
                   &                   )                  &
                   &      *gigaYear                       &
                   &      *clumpingFactor
           end if
           ! Add collisional excitation cooling rate.
           heatingRate=+heatingRate&
                &      -Collisional_Excitation_Cooling_Rate(                &
                &                                           atomicNumber  , &
                &                                           electronNumber, &
                &                                           temperature     &
                &                                          )                &
                &      *clumpingFactor
           ! Compute net rate of change of density.
           propertiesRateOfChange(iProperty)=                                                   &
                ! Collisional ionization.
                & +(+collisionIonizationRateFrom      +collisionIonizationRateTo      )         &
                & *gigaYear                                                                     &
                & *centi**3                                                                     &
                & *clumpingFactor                                                               &
                ! Recombination.
                & +(+recombinationRateFrom            +recombinationRateTo            )         &
                & *gigaYear                                                                     &
                & *centi**3                                                                     &
                & *clumpingFactor                                                               &
                ! Dielectronic recombination.
                & +(+recombinationDielectronicRateFrom+recombinationDielectronicRateTo)         &
                & *gigaYear                                                                     &
                & *centi**3                                                                     &
                & *clumpingFactor                                                               &
                ! Photoionization.
                & +(+ionizationPhotoRateFrom          +ionizationPhotoRateTo          )         &
                & *gigaYear                                                                     &
                ! Cosmological expansion.
                & -3.0d0                                                                        &
                & *cosmologyFunctions_%expansionRate(cosmologyFunctions_%expansionFactor(time)) &
                & *densityThisIon  
        end do
     end do
     ! Compute the rate of change of electron density.
     electronDensityRateOfChange=+propertiesRateOfChange(3) &
          &                      +propertiesRateOfChange(5) & 
          &                      +propertiesRateOfChange(6)
     ! Compute rate of change of temperature due to cosmological expansion.
     propertiesRateOfChange(1)=-2.0d0                                     &
          &                    *cosmologyFunctions_%expansionRate  (      &
          &                     cosmologyFunctions_%expansionFactor (     &
          &                                                          time &
          &                                                         )     &
          &                                                        )      &
          &    * temperature                                                       
     ! Compute rate of change of temperature due to atomic processes.                         
     if (densityTotal > 0.0d0)                                                      &
          & propertiesRateOfChange(1)=                                              &
          &                    +propertiesRateOfChange(1)                           &
          ! Accumulated atomic process heating rate.
          &                    +heatingRate                                         &
          &                    /1.5d0                                               &
          &                    /boltzmannsConstant                                  &
          &                    /densityTotal                                        &
          ! CMB Compton scattering heating/cooling rate.
          &                   +speedLight                                           &
          &                   *densityElectron                                      &
          &                   *4.0d0                                                &
          &                   *thomsonCrossSection                                  &
          &                   *radiationConstant                                    &
          &                   *  cosmologyFunctions_%temperatureCMBEpochal(time)**4 &
          &                   *(                                                    &
          &                     +cosmologyFunctions_%temperatureCMBEpochal(time)    &
          &                     -temperature                                        &
          &                    )                                                    &
          &                   /electronMass                                         &
          &                   /speedLight                                       **2 &
          &                   /1.5d0                                                &
          &                   /densityTotal                                         &
          &                   *gigaYear                                             &
          ! Particle number rate of change. 
          &    +electronDensityRateOfChange                                         & 
          &    /densityTotal                                                        &
          &    +3.0d0                                                               &
          &    *cosmologyFunctions_%expansionRate  (                                &
          &     cosmologyFunctions_%expansionFactor (                               &
          &                                          time                           &
          &                                         )                               &
          &                                        ) 
     ! Transfer rates of change to contiguous array.
     propertiesRateOfChange( 7: 8)=massFilteringCompositeRateOfChange
     propertiesRateOfChange( 9   )=opticalDepthRateOfChange
     propertiesRateOfChange(10   )=massFilteringRateOfChange
     ! Return success.
     Intergalactic_Medium_State_Internal_ODEs=FGSL_Success
     
   contains

     function Photoionization_Rate_Integrand(wavelength,parameterPointer) bind(c)
       !% Integrand function used to compute the rate of photoionizations of an ionic species.
       implicit none
       real            (kind=c_double)        :: Photoionization_Rate_Integrand
       real            (kind=c_double), value :: wavelength
       double precision                       :: photonDensity                 , photonFlux
       type            (c_ptr        ), value :: parameterPointer

       if (wavelength <= 0.0d0) then
          Photoionization_Rate_Integrand=0.0d0
       else
          photonFlux   =+radiation%flux(wavelength) &
               &        /centi**2                   &
               &        *ergs
          photonDensity=+4.0d0                      &
               &        *Pi                         &
               &        *photonFlux                 &
               &        /plancksConstant            & 
               &        /speedLight                 &
               &        /wavelength
          Photoionization_Rate_Integrand=+speedLight                                                        &
               &                         *Atomic_Cross_Section_Ionization_Photo(                            &
               &                                                                atomicNumber              , &
               &                                                                photoionizationGroundIonizationState, &
               &                                                                shellNumber               , &
               &                                                                wavelength                  &
               &                                                               )                            &
               &                         *centi**2                                                          &
               &                         *photonDensity
       end if
       return
     end function Photoionization_Rate_Integrand
 
     function Photoionization_Heating_Rate_Integrand(wavelength,parameterPointer) bind(c)
       !% Integrand function used to compute the rate of photoionization heating of an ionic species.
       implicit none
       real            (kind=c_double)        :: Photoionization_Heating_Rate_Integrand
       real            (kind=c_double), value :: wavelength
       double precision                       :: photonDensity                         , photonFlux
       type            (c_ptr        ), value :: parameterPointer
       
       if (wavelength < 0.0d0) then
          Photoionization_Heating_Rate_Integrand=0.0d0
       else
          photonFlux   =+radiation%flux(wavelength) &
               &        /centi**2                   &
               &        *ergs
          photonDensity=+4.0d0                      &
               &        *Pi                         &
               &        *photonFlux                 &
               &        /plancksConstant            & 
               &        /speedLight                 &
               &        /wavelength
          Photoionization_Heating_Rate_Integrand=+speedLight                                                                  &
               &                                 *Atomic_Cross_Section_Ionization_Photo(                                      &
               &                                                                        atomicNumber                        , &
               &                                                                        photoionizationGroundIonizationState, &
               &                                                                        shellNumber                         , &
               &                                                                        wavelength                            &
               &                                                                       )                                      &
               &                                 *centi**2                                                                    &
               &                                 *photonDensity                                                               &
               &                                 *(                                                                           &
               &                                   +plancksConstant                                                           &
               &                                   *speedLight                                                                &
               &                                   *angstromsPerMeter                                                         &
               &                                   /wavelength                                                                &
               &                                   -Atomic_Ionization_Potential        (                                      &
               &                                                                        atomicNumber                        , &
               &                                                                        photoionizationGroundElectronNumber   &
               &                                                                       )                                      &
               &                                   *electronVolt                                                              &
               &                                  )
       end if
       return
     end function Photoionization_Heating_Rate_Integrand

   end function Intergalactic_Medium_State_Internal_ODEs

   double precision function rLSS(omegaMatter,expansionFactor)
     !% Evaluate the $r_{\mathrm LSS}$ parameter of \cite{naoz_formation_2007} using their fitting formula.
     implicit none
     double precision, intent(in   ) :: omegaMatter     , expansionFactor
     double precision                :: rLSSCoefficient1, rLSSCoefficient2, rLSSCoefficient3

     rLSSCoefficient1=1.0d-4*(-1.99d0*(omegaMatter**2)+2.41d0*omegaMatter+0.21d0)
     rLSSCoefficient2=1.0d-3*(+6.37d0*(omegaMatter**2)-6.99d0*omegaMatter-1.76d0)
     rLSSCoefficient3=1.0d-2*(-1.83d0*(omegaMatter**2)+2.40d0*omegaMatter-0.54d0)
     ! Note that the coefficients for the different exponents of expansion factor are reversed from that given in Naoz &
     ! Barkana. Without this change the fit for rLSS does not work.
     rLSS            =+rLSSCoefficient1/expansionFactor**1.5d0 &
          &           +rLSSCoefficient2/expansionFactor        &
          &           +rLSSCoefficient3 
     return
   end function rLSS
 
 end module Intergalactic_Medium_State_Internal_Evolver
