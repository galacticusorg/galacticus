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

!% Contains a module which performs analysis to compute the stellar mass function of Local Group satellite galaxies.

module Galacticus_Output_Analyses_LG_Satellite_Mass_Functions
  !% Performs analysis to compute the stellar mass function of Local Group satellite galaxies.
  use, intrinsic :: ISO_C_Binding
  use Kind_Numbers
  implicit none
  private
  public :: Galacticus_Output_Analysis_LG_Satellite_Mass_Functions, Galacticus_Output_Analysis_LG_Satellite_Mass_Functions_Output

  ! Record of module initialization.
  logical                                                                :: moduleInitialized =.false.

  ! Number of central galaxies supported, and enumeration.
  integer                   , parameter                                  :: centralGalaxyCount=2
  !@ <enumeration>
  !@  <name>localGroupCentralGalaxy</name>
  !@  <description>Used to specify central galaxy in the Local Group.</description>
  !@  <entry label="localGroupCentralMilkyWay"/>
  !@  <entry label="localGroupCentralM31"     />
  !@ </enumeration>
  integer                   , parameter                                  :: localGroupCentralMilkyWay=1
  integer                   , parameter                                  :: localGroupCentralM31     =2
  character       (len=8  )              , dimension(centralGalaxyCount) :: localGroupCentralLabel   =['MilkyWay','M31     ']

  ! Record of whether this analysis is active.
  logical                                , dimension(centralGalaxyCount) :: analysisActive=.false.

  ! Mass of Local Group halo.
  double precision                       , dimension(centralGalaxyCount) :: analysisLGSatelliteHaloMass
  double precision                       , dimension(centralGalaxyCount) :: analysisLGSatelliteHaloRadius

  ! Output number corresponding to present day.
  integer         (c_size_t)                                                     :: outputPresentDay

  ! Mass function bins.
  integer                                                                :: massBinsCount
  double precision          , allocatable, dimension(:                 ) :: massBins
  integer                   , allocatable, dimension(:,:               ) :: massFunctionCumulative               , massFunctionCumulativeHalo, &
       &                                                                    massFunctionSquaredCumulative
  integer                                , dimension(centralGalaxyCount) :: massFunctionHaloCount
  double precision                       , dimension(centralGalaxyCount) :: massFunctionHaloRadius
  logical                                                                :: workArraysInitialized        =.false.
  logical                                , dimension(centralGalaxyCount) :: treeActive
  !$omp threadprivate(massFunctionCumulativeHalo,workArraysInitialized,treeActive)

contains

  !# <mergerTreeAnalysisTask>
  !#  <unitName>Galacticus_Output_Analysis_LG_Satellite_Mass_Functions</unitName>
  !# </mergerTreeAnalysisTask>
  subroutine Galacticus_Output_Analysis_LG_Satellite_Mass_Functions(thisTree,thisNode,nodeStatus,iOutput,mergerTreeAnalyses)
    !% Construct the stellar mass function of Local Group satellites.
    use, intrinsic :: ISO_C_Binding
    use ISO_Varying_String
    use Input_Parameters
    use Galactic_Structure_Enclosed_Masses
    use Galacticus_Output_Times
    use Cosmology_Functions
    use Cosmology_Parameters
    use Numerical_Comparison
    use Numerical_Ranges
    use Galacticus_Error
    use Galacticus_Nodes
    use Galactic_Structure_Options
    use Galacticus_Output_Merger_Tree_Data
    use Memory_Management
    use Dark_Matter_Halo_Scales
    use Root_Finder
    implicit none
    type            (mergerTree              ), intent(in   )                 :: thisTree
    type            (treeNode                ), intent(inout), pointer        :: thisNode
    integer                                   , intent(in   )                 :: nodeStatus
    integer         (c_size_t                ), intent(in   )                 :: iOutput
    type            (varying_string          ), intent(in   ), dimension(:  ) :: mergerTreeAnalyses
    class           (nodeComponentBasic      )               , pointer        :: hostBasic                      , basic
    class           (cosmologyFunctionsClass )               , pointer        :: cosmologyFunctionsModel
    class           (cosmologyParametersClass)               , pointer        :: cosmologyParametersModel
    class           (darkMatterHaloScaleClass)               , pointer        :: darkMatterHaloScale_
    type            (treeNode                )               , pointer        :: host                           , node
    double precision                          , parameter                     :: redshiftPresentDay       =0.0d0
    double precision                          , parameter                     :: massMinimum              =1.0d0
    double precision                          , parameter                     :: massMaximum              =1.0d15
    integer                                   , parameter                     :: massFunctionBinsPerDecade=10
    double precision                                                          :: mass                            , timeNow     , &
         &                                                                       massEnclosed                    , radiusVirial
    integer         (c_size_t                )                                :: k
    integer                                                                   :: j                               , i
    character       (len=64                  )                                :: parameterName
    type            (rootFinder              )                                :: finder

    ! Initialize the module if necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Galacticus_Output_Analysis_LG_Satellite_MF_Initialize)
       if (.not.moduleInitialized) then
          do i=1,centralGalaxyCount
             analysisActive(i)=any(trim(mergerTreeAnalyses) == "localGroupSatelliteMassFunction"//trim(localGroupCentralLabel(i)))
             if (analysisActive(i)) then
                !@ <inputParameter>
                !@   <regEx>analysisLGSatelliteMFHaloMass(MilkyWay|M31)</regEx>
                !@   <defaultValue>$10^{12}M_\odot$</defaultValue>
                !@   <attachedTo>module</attachedTo>
                !@   <description>
                !@    The mass of the host halo to be used when constructing Local Group satellite galaxy stellar mass functions.
                !@   </description>
                !@   <type>string</type>
                !@   <cardinality>0..1</cardinality>
                !@   <group>output</group>
                !@ </inputParameter>
                parameterName="analysisLGSatelliteMFHaloMass"//trim(localGroupCentralLabel(i))
                call Get_Input_Parameter(trim(parameterName),analysisLGSatelliteHaloMass(i),defaultValue=1.0d12)
                ! Construct a temporary node, assign it the given mass, and find the virial radius.
                cosmologyFunctionsModel  => cosmologyFunctions ()
                cosmologyParametersModel => cosmologyParameters()
                darkMatterHaloScale_     => darkMatterHaloScale()
                timeNow                  =  cosmologyFunctionsModel%cosmicTime(1.0d0)              
                node                     => treeNode           ()
                basic                    => node%basic(autoCreate=.true.)
                call basic%timeSet(timeNow)
                call basic%massSet(analysisLGSatelliteHaloMass(i))
                radiusVirial=darkMatterHaloScale_%virialRadius(node)
                !@ <inputParameter>
                !@   <regEx>analysisLGSatelliteMFHaloRadius(MilkyWay|M31)</regEx>
                !@   <defaultValue>virial radius</defaultValue>
                !@   <attachedTo>module</attachedTo>
                !@   <description>
                !@    The radius within whcih the mass of the host halo is defined when constructing Local Group satellite galaxy stellar mass functions.
                !@   </description>
                !@   <type>string</type>
                !@   <cardinality>0..1</cardinality>
                !@   <group>output</group>
                !@ </inputParameter>
                parameterName="analysisLGSatelliteMFHaloRadius"//trim(localGroupCentralLabel(i))
                call Get_Input_Parameter(trim(parameterName),analysisLGSatelliteHaloRadius(i),defaultValue=radiusVirial)
                ! Solve for the virial mass that gives the required enclosed mass.
                call finder%tolerance(1.0d-6,1.0d-6)
                call finder%rangeExpand(rangeExpandUpward=2.0d0,rangeExpandDownward=0.5d0,rangeExpandType=rangeExpandMultiplicative)
                call finder%rootFunction(MW_Satellite_Mass_Function_Enclosed_Mass)
                analysisLGSatelliteHaloMass(i)=finder%find(rootGuess=analysisLGSatelliteHaloMass(i))
                call node%destroy()
             end if
          end do
          ! Continue only if an analysis is active.
          if (any(analysisActive)) then
             ! Find the output corresponding to the present day.
             cosmologyFunctionsModel => cosmologyFunctions()
             outputPresentDay=-1
             do k=1,Galacticus_Output_Time_Count()
                if     (                                                                                             &
                     &  Values_Agree(                                                                                &
                     &               redshiftPresentDay                                                            , &
                     &               cosmologyFunctionsModel%redshiftFromExpansionFactor(                            &
                     &               cosmologyFunctionsModel%expansionFactor             (                           &
                     &                                                                    Galacticus_Output_Time(k)  &
                     &                                                                   )                           &
                     &                                                                  )                          , &
                     &               absTol=0.001d0                                                                  &
                     &              )                                                                                &
                     & ) then
                   outputPresentDay=k
                   exit
                end if
             end do
             if (outputPresentDay < 0)                                                                       &
                  & call Galacticus_Error_Report(                                                            &
                  &                                'Galacticus_Output_Analysis_LG_Satellite_Mass_Functions', &
                  &                                'unable to find present day [z=0] in outputs'             &
                  &                             )
             ! Construct mass bins.
             massBinsCount=int(log10(massMaximum/massMinimum)*dble(massFunctionBinsPerDecade))+1
             call Alloc_Array(massBins                     ,[massBinsCount                   ])
             call Alloc_Array(massFunctionCumulative       ,[massBinsCount,centralGalaxyCount])
             call Alloc_Array(massFunctionSquaredCumulative,[massBinsCount,centralGalaxyCount])
             massBins=Make_Range(massMinimum,massMaximum,massBinsCount,rangeTypeLogarithmic)
             massFunctionHaloCount        =0
             massFunctionHaloRadius       =0.0d0
             massFunctionCumulative       =0
             massFunctionSquaredCumulative=0
             ! Record that module is initialized.
             moduleInitialized=.true.
          end if
       end if
       !$omp end critical(Galacticus_Output_Analysis_LG_Satellite_MF_Initialize)
    end if
    ! Return if this analysis is not active.
    if (.not.any(analysisActive)   ) return
    ! Return if this output is not the present day.
    if (iOutput /= outputPresentDay) return
    ! Initialize mass function count array if necessary.
    if (.not.workArraysInitialized) then
       call Alloc_Array(massFunctionCumulativeHalo,[massBinsCount,centralGalaxyCount])
       massFunctionCumulativeHalo=0
       treeActive                =.false.
       workArraysInitialized     =.true.
    end if
    ! Accumulate halo if necessary.
    if (nodeStatus == nodeStatusFinal) then
       do i=1,centralGalaxyCount
          if (treeActive(i)) call MW_Satellite_Mass_Function_Accumulate_Halo(i,massFunctionCumulativeHalo)
       end do
       return
    end if
    ! Find the host halo.
    host => thisNode
    do while (host%isSatellite())
       host => host%parent
    end do
    hostBasic => host%basic()
    ! Iterate over central galaxies.
    do i=1,centralGalaxyCount
       ! Skip non-active analyses.
       if (analysisActive(i)) then
          ! Check if host halo mass agrees with Milky Way halo mass.
          if (Values_Differ(hostBasic%mass(),analysisLGSatelliteHaloMass(i),relTol=1.0d-3)) cycle
          ! Proceed for satellite galaxies.
          if (.not.treeActive(i)) then
             darkMatterHaloScale_      => darkMatterHaloScale                (    )
             !$omp atomic
             massFunctionHaloRadius(i) =  massFunctionHaloRadius             (   i) &
                  &                      +darkMatterHaloScale_  %virialRadius(host)
             treeActive            (i) =  .true.
          end if
          if (thisNode%isSatellite()) then
             ! Get the galactic mass.
             mass=                                                                                                                       &
                  &  Galactic_Structure_Enclosed_Mass(thisNode,radiusLarge,componentType=componentTypeDisk    ,massType=massTypeStellar) &
                  & +Galactic_Structure_Enclosed_Mass(thisNode,radiusLarge,componentType=componentTypeSpheroid,massType=massTypeStellar)
             ! Accumulate this galaxy.
             do j=1,massBinsCount
                if (mass > massBins(j)) massFunctionCumulativeHalo(j,i)=massFunctionCumulativeHalo(j,i)+1
             end do
          end if
          ! Accumulate halo if necessary.
          if (nodeStatus == nodeStatusLast) call MW_Satellite_Mass_Function_Accumulate_Halo(i,massFunctionCumulativeHalo)
       end if
    end do
    return

  contains

    double precision function MW_Satellite_Mass_Function_Enclosed_Mass(haloMass)
      !% Root finding function used to set the host halo mass for Local Group satellite mass function calculations.
      use Galacticus_Calculations_Resets
      implicit none
      double precision, intent(in   ) :: haloMass

      call basic%massSet(haloMass)
      call Galacticus_Calculations_Reset(node)
      MW_Satellite_Mass_Function_Enclosed_Mass=                                                             &
           & +Galactic_Structure_Enclosed_Mass(node,analysisLGSatelliteHaloRadius(i),massType=massTypeDark) &
           & *  cosmologyParametersModel%OmegaMatter()                                                      &
           & /(                                                                                             &
           &   +cosmologyParametersModel%OmegaMatter()                                                      &
           &   -cosmologyParametersModel%OmegaBaryon()                                                      &
           &  )                                                                                             &
           & -analysisLGSatelliteHaloMass(i)
      return
    end function MW_Satellite_Mass_Function_Enclosed_Mass

  end subroutine Galacticus_Output_Analysis_LG_Satellite_Mass_Functions

  subroutine MW_Satellite_Mass_Function_Accumulate_Halo(localGroupCentral,massFunctionCumulativeHalo)
    !% Accumulate the cumulative mass function for a single halo to the global arrays.
    implicit none
    integer         , intent(in   )                 :: localGroupCentral
    integer         , intent(inout), dimension(:,:) :: massFunctionCumulativeHalo

    !$omp critical(MW_Satellite_Mass_Function_Accumulate)
    massFunctionCumulative       (:,localGroupCentral)=massFunctionCumulative       (:,localGroupCentral)+massFunctionCumulativeHalo(:,localGroupCentral)
    massFunctionSquaredCumulative(:,localGroupCentral)=massFunctionSquaredCumulative(:,localGroupCentral)+massFunctionCumulativeHalo(:,localGroupCentral)**2
    massFunctionHaloCount        (  localGroupCentral)=massFunctionHaloCount        (  localGroupCentral)+1
    !$omp end critical(MW_Satellite_Mass_Function_Accumulate)
    massFunctionCumulativeHalo(:,localGroupCentral)=0
    treeActive                (  localGroupCentral)=.false.
    return
  end subroutine MW_Satellite_Mass_Function_Accumulate_Halo

  !# <hdfPreCloseTask>
  !#  <unitName>Galacticus_Output_Analysis_LG_Satellite_Mass_Functions_Output</unitName>
  !# </hdfPreCloseTask>
  subroutine Galacticus_Output_Analysis_LG_Satellite_Mass_Functions_Output()
    !% Outputs galaxy mass functions to file.
    use Galacticus_HDF5
    use Numerical_Constants_Astronomical
    use String_Handling
    implicit none
    double precision            , dimension(massBinsCount) :: massFunctionCumulativeMean, massFunctionCumulativeMeanVariance
    type            (hdf5Object)                           :: analysisGroup             , massFunctionGroup                 , &
         &                                                    thisDataset
    integer                                                :: i

    ! Iterate over central galaxies.
    do i=1,centralGalaxyCount
       ! Skip if this analysis is not active.
       if (.not.analysisActive(i)) cycle   
       ! Compute the mean mass function and its variance.
       if (massFunctionHaloCount(i) > 0) then
          massFunctionCumulativeMean           =  dble(massFunctionCumulative       (:,i))/dble(massFunctionHaloCount(i))
          massFunctionHaloRadius            (i)=+     massFunctionHaloRadius(i)                                           &
               &                                /dble(massFunctionHaloCount (i)  )
       else
          massFunctionCumulativeMean           =0.0d0
          massFunctionHaloRadius            (i)=0.0d0
       end if
       if (massFunctionHaloCount(i) > 1) then
          massFunctionCumulativeMeanVariance   =( dble(massFunctionSquaredCumulative(:,i))/dble(massFunctionHaloCount(i)) &
               &                                 -massFunctionCumulativeMean**2                                           &
               &                                )                                                                         &
               &                                *dble(massFunctionHaloCount (i)  )                                        &
               &                                /dble(massFunctionHaloCount (i)-1)
       else
          massFunctionCumulativeMeanVariance   =0.0d0
       end if       
       ! Output the mass function.
       !$omp critical(HDF5_Access)
       analysisGroup    =galacticusOutputFile%openGroup('analysis'                                                              ,'Model analysis'                     )
       massFunctionGroup=analysisGroup       %openGroup(trim(String_Lower_Case_First(localGroupCentralLabel(i)))//'MassFunction','Stellar mass function of satellites')
       call massFunctionGroup%writeAttribute(massFunctionHaloRadius(i),'haloRadius')
       call massFunctionGroup%writeDataset  (massFunctionCumulativeMean        ,'massFunctionCumulative'        ,'Cumulative mass function'                                     )
       call massFunctionGroup%writeDataset  (massFunctionCumulativeMeanVariance,'massFunctionCumulativeVariance','Cumulative mass function variance'                            )
       call massFunctionGroup%writeDataset  (massBins                          ,'massStellar'                   ,'Stellar mass'                     ,datasetReturned=thisDataset)
       call thisDataset      %writeAttribute(massSolar            ,'unitsInSI')
       call thisDataset      %close()
       call massFunctionGroup%writeAttribute(massFunctionHaloCount,'haloCount')
       call massFunctionGroup%close()
       call analysisGroup    %close()
       !$omp end critical(HDF5_Access)
    end do
    return
  end subroutine Galacticus_Output_Analysis_LG_Satellite_Mass_Functions_Output

end module Galacticus_Output_Analyses_LG_Satellite_Mass_Functions
