!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!!{
Contains a program which tests the implementation of the \cite{gnedin_response_2004} dark matter profile.
!!}

program Test_Dark_Matter_Profiles_Gnedin2004
  !!{
  Tests the implementation of the \cite{gnedin_response_2004} dark matter profile.
  !!}
  use :: Calculations_Resets       , only : Calculations_Reset
  use :: Cosmology_Parameters      , only : cosmologyParametersSimple
  use :: Cosmology_Functions       , only : cosmologyFunctionsMatterLambda
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMONFW
  use :: Dark_Matter_Profiles      , only : darkMatterProfileAdiabaticGnedin2004
  use :: Mass_Distributions        , only : nonAnalyticSolversNumerical                                   , massDistributionClass            , massDistributionSphericalAdiabaticGnedin2004
  use :: File_Utilities            , only : Count_Lines_in_File
  use :: Virial_Density_Contrast   , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: Display                   , only : displayVerbositySet                                           , verbosityLevelStandard
  use :: Galacticus_Nodes          , only : nodeClassHierarchyFinalize                                    , nodeClassHierarchyInitialize     , nodeComponentBasic                          , nodeComponentDarkMatterProfile, &
          &                                 treeNode                                                      , nodeComponentSPheroid
  use :: Input_Parameters          , only : inputParameters
  use :: Input_Paths               , only : inputPath                                                     , pathTypeExec
  use :: ISO_Varying_String        , only : varying_string                                                , assignment(=)                    , char                                        , var_str                       , &
       &                                    operator(//)
  use :: Node_Components           , only : Node_Components_Initialize                                    , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize         , Node_Components_Uninitialize
  use :: Unit_Tests                , only : Assert                                                        , Unit_Tests_Begin_Group           , Unit_Tests_End_Group                        , Unit_Tests_Finish
  implicit none
  type            (treeNode                                                      ), pointer                   :: node
  class           (nodeComponentBasic                                            ), pointer                   :: basic
  class           (nodeComponentDarkMatterProfile                                ), pointer                   :: dmProfile
  class           (nodeComponentSpheroid                                         ), pointer                   :: spheroid
  double precision                                                                , parameter                 :: concentration                        =10.0d0, massVirial              =1.0d12, &
       &                                                                                                         fractionBaryons                      =0.15d0, radiusFractionalBaryons =0.03d0
  type            (darkMatterProfileDMONFW                                       ), pointer                   :: darkMatterProfileDMONFW_
  type            (darkMatterProfileAdiabaticGnedin2004                          ), pointer                   :: darkMatterProfileAdiabaticGnedin2004_
  type            (cosmologyParametersSimple                                     ), pointer                   :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                                ), pointer                   :: cosmologyFunctions_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition            ), pointer                   :: darkMatterHaloScale_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer                   :: virialDensityContrast_
  class           (massDistributionClass                                         ), pointer                   :: massDistribution_
  type            (inputParameters                                               )                            :: parameters
  double precision                                                                                            :: radiusScale                                  , radiusVirial                  , &
       &                                                                                                         radiusFractionalFinalContra
  double precision                                                                , allocatable, dimension(:) :: radiusFractionalInitialContra                , radiusFractionalInitial
  integer                                                                                                     :: i                                            , fileUnit                      , &
       &                                                                                                         countLines                                   , countRadii
  type            (varying_string                                                )                            :: fileName
  
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize event hooks and global functions.
  call eventsHooksInitialize()
  call Functions_Global_Set ()
  ! Initialize node components.
  parameters=inputParameters(var_str('testSuite/parameters/darkMatterProfilesAdiabaticGnedin2004.xml'))
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Build required objects. Use a cosmology with a baryon fraction of 0.15 to match what was assumed in the reference run of the "contra" code.
  allocate(cosmologyParameters_                 )
  allocate(cosmologyFunctions_                  )
  allocate(virialDensityContrast_               )
  allocate(darkMatterHaloScale_                 )
  allocate(darkMatterProfileDMONFW_             )
  allocate(darkMatterProfileAdiabaticGnedin2004_)
  !![
  <referenceConstruct object="cosmologyParameters_"                 >
   <constructor>
    cosmologyParametersSimple                                     (                                                                            &amp;
     &amp;                                                         OmegaMatter                         = 0.3000d0                            , &amp;
     &amp;                                                         OmegaBaryon                         = 0.3000d0*fractionBaryons            , &amp;
     &amp;                                                         OmegaDarkEnergy                     = 0.7000d0                            , &amp;
     &amp;                                                         temperatureCMB                      = 2.7800d0                            , &amp;
     &amp;                                                         HubbleConstant                      =69.3000d0                              &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"                  >
   <constructor>
    cosmologyFunctionsMatterLambda                                (                                                                            &amp;
     &amp;                                                         cosmologyParameters_                =cosmologyParameters_                   &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="virialDensityContrast_"               >
   <constructor>
    virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                                            &amp;
     &amp;                                                         tableStore                          =.true.                               , &amp;
     &amp;                                                         cosmologyFunctions_                 =cosmologyFunctions_                    &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterHaloScale_"                 >
   <constructor>
    darkMatterHaloScaleVirialDensityContrastDefinition            (                                                                            &amp;
     &amp;                                                         cosmologyParameters_                =cosmologyParameters_                 , &amp;
     &amp;                                                         cosmologyFunctions_                 =cosmologyFunctions_                  , &amp;
     &amp;                                                         virialDensityContrast_              =virialDensityContrast_                 &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileDMONFW_"             >
   <constructor>
    darkMatterProfileDMONFW                                       (                                                                            &amp;
     &amp;                                                         velocityDispersionUseSeriesExpansion=.false.                              , &amp;
     &amp;                                                         darkMatterHaloScale_                =darkMatterHaloScale_                   &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileAdiabaticGnedin2004_">
   <constructor>
     darkMatterProfileAdiabaticGnedin2004                         (                                                                            &amp;
     &amp;                                                         A                                   =0.85d0                               , &amp;
     &amp;                                                         omega                               =0.80d0                               , &amp;
     &amp;                                                         radiusFractionalPivot               =1.00d0                               , &amp;
     &amp;                                                         toleranceRelative                   =1.0d-2                               , &amp;
     &amp;                                                         nonAnalyticSolver                   =nonAnalyticSolversNumerical          , &amp;
     &amp;                                                         cosmologyParameters_                =cosmologyParameters_                 , &amp;
     &amp;                                                         darkMatterHaloScale_                =darkMatterHaloScale_                 , &amp;
     &amp;                                                         darkMatterProfileDMO_               =darkMatterProfileDMONFW_               &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  !!]
  ! Create a node.
  node      => treeNode                  (                 )
  ! Create components.
  basic     => node    %basic            (autoCreate=.true.)
  dmProfile => node    %darkMatterProfile(autoCreate=.true.)
  spheroid  => node    %spheroid         (autoCreate=.true.)
  ! Set properties.
  call basic%timeSet            (cosmologyFunctions_%cosmicTime(1.0d0))
  call basic%timeLastIsolatedSet(cosmologyFunctions_%cosmicTime(1.0d0))
  call basic%massSet            (                    massVirial       )
  ! Compute scale radius.
  radiusVirial=+darkMatterHaloScale_%radiusVirial(node)
  radiusScale =+radiusVirial &
       &       /concentration
  call dmProfile%scaleSet(radiusScale)
  ! Construct spheroid.
  call spheroid%massStellarSet(fractionBaryons        *massVirial  )
  call spheroid%radiusSet     (radiusFractionalBaryons*radiusVirial)
  ! Get the mass distribution.
  massDistribution_ => darkMatterProfileAdiabaticGnedin2004_%get(node)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group('Gnedin et al. (2004) dark matter profile')
  ! Read data from the reference "contra" file.
  fileName  =inputPath(pathTypeExec)//'testSuite/data/adiabaticContractionContra.txt'
  countLines=Count_Lines_in_File(fileName    )
  countRadii=Count_Lines_in_File(fileName,'#')
  allocate(radiusFractionalInitialContra(countRadii))
  allocate(radiusFractionalInitial      (countRadii))
  open(newunit=fileUnit,file=char(fileName),status='old',form='formatted')
  do i=1,(countLines-countRadii)
     read (fileUnit,*)
  end do
  select type (massDistribution_)
  type is (massDistributionSphericalAdiabaticGnedin2004)
     do i=1,countRadii
        read (fileUnit,*) radiusFractionalInitialContra(i),radiusFractionalFinalContra
        ! Evaluate the initial radius corresponding to this final radius.
        radiusFractionalInitial(i)=+massDistribution_%radiusInitial(radiusFractionalFinalContra*radiusVirial) &
             &                     /                                radiusVirial
     end do
  end select
  close(fileUnit)
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  ! Use a tolerance of 2% in the assertion, as our solver for the initial radius uses a tolerance of 1%.
  call Assert('initial radii',radiusFractionalInitialContra,radiusFractionalInitial,relTol=2.0d-2)
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
  ! Uninitialize node components.
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
  ! Clean up objects.
  call node%destroy()
  deallocate(radiusFractionalInitialContra,radiusFractionalInitial,node)
  !![
  <objectDestructor name="cosmologyParameters_"                 />
  <objectDestructor name="cosmologyFunctions_"                  />
  <objectDestructor name="virialDensityContrast_"               />
  <objectDestructor name="darkMatterHaloScale_"                 />
  <objectDestructor name="darkMatterProfileDMONFW_"             />
  <objectDestructor name="darkMatterProfileAdiabaticGnedin2004_"/>
  !!]
end program Test_Dark_Matter_Profiles_Gnedin2004
