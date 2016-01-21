!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Generate a mock realization of a halo model.

program Halo_Model_Mock
  !% Generates a mock realization from an input halo catalog and a halo model prescription.
  use               Command_Arguments
  use               ISO_Varying_String
  use               Memory_Management
  use               Galactic_Structure_Enclosed_Masses
  use               Galactic_Structure_Options
  use               Input_Parameters
  use               Galacticus_Error
  use               Galacticus_Display
  use               Geometry_Surveys
  use               Galacticus_Nodes
  use               Node_Components
  use               Cosmology_Functions
  use               Conditional_Mass_Functions
  use               Dark_Matter_Profiles_Concentration
  use               IO_HDF5
  use               IO_IRATE
  use               Pseudo_Random
  use               Root_Finder
  use               String_Handling
  use               Numerical_Constants_Prefixes
  use               Numerical_Constants_Astronomical
  use               FGSL
  use               Poisson_Random
  use               Gaussian_Random
  use               Dark_Matter_Profile_Scales
  use               Dark_Matter_Halo_Scales
  use               Dark_Matter_Profiles
  use               Galacticus_Calculations_Resets
  implicit none
  double precision                                , allocatable, dimension(  :) :: haloMass                  , galaxyMass
  double precision                                , allocatable, dimension(:,:) :: haloPosition              , haloVelocity             , &
       &                                                                           galaxyPosition            , galaxyVelocity
  double precision                                             , dimension(3  ) :: satellitePosition         , satelliteVelocity
  class           (cosmologyFunctionsClass       ), pointer                     :: cosmologyFunctions_
  class           (darkMatterProfileClass        ), pointer                     :: darkMatterProfile_
  class           (darkMatterHaloScaleClass      ), pointer                     :: darkMatterHaloScale_
  class           (conditionalMassFunctionClass  ), pointer                     :: conditionalMassFunction_
  type            (treeNode                      ), pointer                     :: node
  class           (nodeComponentBasic            ), pointer                     :: basic
  class           (nodeComponentDarkMatterProfile), pointer                     :: profile
  integer                                                                       :: iHalo                     , galaxyCount              , &
       &                                                                           satelliteNumberActual     , iSatellite               , &
       &                                                                           iAxis
  type            (varying_string                )                              :: parameterFileName         , haloCatalogFileName      , & 
       &                                                                           galaxyCatalogFileName     , message
  type            (hdf5Object                    )                              :: haloCatalogFile           , snapshotGroup            , &
       &                                                                           haloCatalogGroup          , galaxyCatalogFile        , &
       &                                                                           inputGroup                , outputGroup              , &
       &                                                                           thisDataset               , simulationGroup
  type            (irate                         )                              :: haloFile                  , galaxyFile
  double precision                                                              :: haloModelMockMassMinimum  , haloModelMockMassMaximum , &
       &                                                                           probabilityCentral        , xCentral                 , &
       &                                                                           massGalaxy                , satelliteNumberMean      , &
       &                                                                           xSatellite                , redshift                 , &
       &                                                                           satelliteRadius           , satelliteTheta           , &
       &                                                                           satellitePhi              , satelliteVelocityCircular, &
       &                                                                           simulationBoxSize         , populatedHaloMassMinimum
  type            (pseudoRandom                  )                              :: randomSequence
  type            (rootFinder                    )                              :: finderCentral             , finderSatellite
  type            (fgsl_rng                      )                              :: poissonSampler            , gaussianSampler
  logical                                                                       :: poissonSamplerReset=.true., gaussianSamplerReset=.true.
  character       (len=6                         )                              :: label

  ! Read in basic code memory usage.
  call Code_Memory_Usage('Halo_Model_Mock.size')
  ! Check that correct number of arguments have been supplied.
  if (Command_Argument_Count() /= 3) call Galacticus_Error_Report(message="Usage: Halo_Model_Mock.exe <parameterFileName> <haloCatalog> <galaxyCatalog>")
  ! Get the arguments.
  call Get_Argument              (1,parameterFileName    )
  call Get_Argument              (2,haloCatalogFileName  )
  call Get_Argument              (3,galaxyCatalogFileName)
  ! Open the parameter file.
  call Input_Parameters_File_Open(  parameterFileName    )
  ! Read parameters controlling the calculation.
  !@ <inputParameter>
  !@   <name>haloModelMockMassMinimum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <description>
  !@     The minimum mass galaxy to include in a mock halo model realization.
  !@   </description>
  !@   <type>float</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('haloModelMockMassMinimum',haloModelMockMassMinimum)
  !@ <inputParameter>
  !@   <name>haloModelMockMassMaximum</name>
  !@   <attachedTo>program</attachedTo>
  !@   <defaultValue>$10^{16}M_\odot$</defaultValue>
  !@   <description>
  !@     The minimum mass galaxy to include in a mock halo model realization.
  !@   </description>
  !@   <type>float</type>
  !@   <cardinality>1</cardinality>
  !@ </inputParameter>
  call Get_Input_Parameter('haloModelMockMassMaximum',haloModelMockMassMaximum,defaultValue=1.0d16)
  ! Initialize nodes infrastructure.
  call Galacticus_Nodes_Initialize()
  call Node_Components_Initialize ()
  ! Get required objects.
  cosmologyFunctions_      => cosmologyFunctions     ()
  darkMatterProfile_       => darkMatterProfile      ()
  darkMatterHaloScale_     => darkMatterHaloScale    ()
  conditionalMassFunction_ => conditionalMassFunction()
  ! Read the halo catalog.
  call Galacticus_Display_Indent("Reading halo catalog")
  haloFile=irate(char(haloCatalogFileName))
  call haloFile%readHalos     (                            &
       &                       snapshot=1                , &
       &                       redshift=redshift         , &
       &                       center  =haloPosition     , &
       &                       velocity=haloVelocity     , &
       &                       mass    =haloMass           &
       &                      )
  call haloFile%readSimulation(                            &
       &                       boxSize =simulationBoxSize  &
       &                      )
  call Galacticus_Display_Unindent("done")
  ! Establish root finders.
  call finderCentral  %tolerance   (                                               &
       &                            toleranceRelative  =1.0d-6                     &
       &                           )
  call finderCentral  %rangeExpand (                                               &
       &                            rangeExpandUpward  =2.0d0                    , &
       &                            rangeExpandDownward=1.0d0                    , &
       &                            rangeExpandType    =rangeExpandMultiplicative  &
       &                           )
  call finderCentral  %rootFunction(                                               &
       &                                                centralMassRoot            &
       &                           )
  call finderSatellite%tolerance   (                                               &
       &                            toleranceRelative  =1.0d-6                     &
       &                           )
  call finderSatellite%rangeExpand (                                               &
       &                            rangeExpandUpward  =2.0d0                    , &
       &                            rangeExpandDownward=1.0d0                    , &
       &                            rangeExpandType    =rangeExpandMultiplicative  &
       &                           )
  call finderSatellite%rootFunction(                                               &
       &                                                satelliteMassRoot          &
       &                           )
  ! Establish tree node for sampling satellite positions.
  node    => treeNode                  (                 )
  basic   => node    %basic            (autoCreate=.true.)
  profile => node    %darkMatterProfile(autoCreate=.true.)
  call basic%timeSet            (cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)))
  call basic%timeLastIsolatedSet(cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)))
  ! Iterate over halos.
  call Galacticus_Display_Indent("Populating halos")
  galaxyCount             =     0
  populatedHaloMassMinimum=huge(0.0d0)
  do iHalo=1,size(haloMass)
     ! Get probability of central galaxy.
     probabilityCentral=+conditionalMassFunction_%massFunction(haloMass(iHalo),haloModelMockMassMinimum,haloModelGalaxyTypeCentral) &
          &             -conditionalMassFunction_%massFunction(haloMass(iHalo),haloModelMockMassMaximum,haloModelGalaxyTypeCentral)
     ! Test for inclusion.
     if (randomSequence%sample() <= probabilityCentral) then
        ! Sample central galaxy mass.        
        xCentral  =randomSequence%sample()*probabilityCentral
        massGalaxy=finderCentral%find(rootGuess=haloModelMockMassMinimum)
        call Galaxy_Add(massGalaxy,haloPosition(:,iHalo),haloVelocity(:,iHalo))
     end if
     ! Get mean number of satellite galaxies.
     satelliteNumberMean=+conditionalMassFunction_%massFunction(haloMass(iHalo),haloModelMockMassMinimum,haloModelGalaxyTypeSatellite) &
          &              -conditionalMassFunction_%massFunction(haloMass(iHalo),haloModelMockMassMaximum,haloModelGalaxyTypeSatellite)
     satelliteNumberActual=Poisson_Random_Get(poissonSampler,satelliteNumberMean,reset=poissonSamplerReset)
     if (satelliteNumberActual > 0) then
        ! Construct the dark matter halo profile.
        call basic  %massSet              (haloMass                 (iHalo))
        call Galacticus_Calculations_Reset(                          node  )
        call profile%scaleSet             (Dark_Matter_Profile_Scale(node ))
        call Galacticus_Calculations_Reset(                          node  )
        do iSatellite=1,satelliteNumberActual
           ! Sample satellite galaxy mass.        
           xSatellite               =     randomSequence%sample()*satelliteNumberMean
           massGalaxy               =finderSatellite%find(rootGuess=haloModelMockMassMinimum)
           ! Sample galaxy radial position.
           xSatellite               =     randomSequence%sample()
           satelliteRadius          =Galactic_Structure_Radius_Enclosing_Mass(node,fractionalMass=xSatellite,massType=massTypeDark)
           ! Get circular velocity at this radius.
           satelliteVelocityCircular=darkMatterProfile_%circularVelocity(node,satelliteRadius)
           ! Convert radial position to comoving coordinates.
           satelliteRadius          =satelliteRadius*(1.0d0+redshift)
           ! Sample galaxy angular position.
           satellitePhi             =     randomSequence%sample()*2.0d0*Pi
           satelliteTheta           =acos(randomSequence%sample()*2.0d0-1.0d0)
           ! Set satellite position.
           satellitePosition        =+satelliteRadius                                                         &
                &                    *[                                                                       &
                &                      sin(satelliteTheta)*sin(satellitePhi)                                , &
                &                      sin(satelliteTheta)*cos(satellitePhi)                                , &
                &                      cos(satelliteTheta)                                                    &
                &                     ]                                                                       &
                &                    +haloPosition(:,iHalo)           
           ! Periodicalize the satellite position.
           do iAxis=1,3
              do while (satellitePosition(iAxis) < 0.0d0)
                 satellitePosition(iAxis)=satellitePosition(iAxis)+simulationBoxSize
              end do
              do while (satellitePosition(iAxis) > simulationBoxSize)
                 satellitePosition(iAxis)=satellitePosition(iAxis)-simulationBoxSize
              end do
           end do
           ! Set satellite velocity.
           satelliteVelocity        =+satelliteVelocityCircular                                               &
                &                    /sqrt(3.0d0)                                                             &
                &                    *[                                                                       &
                &                      Gaussian_Random_Get(gaussianSampler,1.0d0,reset=gaussianSamplerReset), &
                &                      Gaussian_Random_Get(gaussianSampler,1.0d0,reset=gaussianSamplerReset), &
                &                      Gaussian_Random_Get(gaussianSampler,1.0d0,reset=gaussianSamplerReset)  &
                &                    ]                                                                        &
                &                    +haloVelocity(:,iHalo)
           ! Store the satellite.
           call Galaxy_Add(massGalaxy,satellitePosition,satelliteVelocity)
        end do
     end if
  end do
  message="Created "
  message=message//galaxyCount//" galaxies"
  call Galacticus_Display_Message(message)
  write (label,'(f5.2)') log10(populatedHaloMassMinimum)
  message="Lowest mass halo populated has log₁₀(Mₕₐₗₒ/M☉)="//trim(adjustl(label))
  call Galacticus_Display_Message(message)
  call Galacticus_Display_Unindent("done")
  ! Output galaxy catalog.
  galaxyFile=irate(char(galaxyCatalogFileName))
  call haloFile%copyCosmology (galaxyFile)
  call haloFile%copySimulation(galaxyFile)
  call galaxyFile%writeHalos(1,redshift,galaxyPosition(:,1:galaxyCount),galaxyVelocity(:,1:galaxyCount),galaxyMass(1:galaxyCount))
  ! Close the parameter file.
  call Input_Parameters_File_Close

contains
  
  double precision function centralMassRoot(mass)
    !% Root function used to find the mass of central galaxies
    implicit none
    double precision, intent(in   ) :: mass
    double precision                :: x

    x=      (                                                                                                            &
         &   +conditionalMassFunction_%massFunction(haloMass(iHalo),             mass       ,haloModelGalaxyTypeCentral) &
         &   -conditionalMassFunction_%massFunction(haloMass(iHalo),haloModelMockMassMaximum,haloModelGalaxyTypeCentral) &
         &  )
    centralMassRoot=x-xCentral
    return
  end function centralMassRoot
  
  double precision function satelliteMassRoot(mass)
    !% Root function used to find the mass of satellite galaxies
    implicit none
    double precision, intent(in   ) :: mass
    double precision                :: x

    x=      (                                                                                                              &
         &   +conditionalMassFunction_%massFunction(haloMass(iHalo),             mass       ,haloModelGalaxyTypeSatellite) &
         &   -conditionalMassFunction_%massFunction(haloMass(iHalo),haloModelMockMassMaximum,haloModelGalaxyTypeSatellite) &
         &  )
    satelliteMassRoot=x-xSatellite
    return
  end function satelliteMassRoot
  
  subroutine Galaxy_Add(mass,position,velocity)
    !% Add a galaxy to the output buffers.
    implicit none
    double precision, intent(in   )                 :: mass
    double precision, intent(in   ), dimension(3  ) :: position                     , velocity
    integer         , parameter                     :: galaxyBufferSizeMinimum=10000
    double precision, allocatable  , dimension(  :) :: galaxyMassTmp
    double precision, allocatable  , dimension(:,:) :: galaxyPositionTmp            , galaxyVelocityTmp

    ! Expand output buffers as needed.
    galaxyCount=galaxyCount+1
    if (allocated(galaxyMass)) then
       if (galaxyCount > size(galaxyMass)) then
          call Move_Alloc(galaxyMass    ,galaxyMassTmp    )
          call Move_Alloc(galaxyPosition,galaxyPositionTmp)
          call Move_Alloc(galaxyVelocity,galaxyVelocityTmp)
          call Alloc_Array(galaxyMass    ,[  2*(galaxyCount-1)])
          call Alloc_Array(galaxyPosition,[3,2*(galaxyCount-1)])
          call Alloc_Array(galaxyVelocity,[3,2*(galaxyCount-1)])
          galaxyMass    (  1:galaxyCount-1)=galaxyMassTmp
          galaxyPosition(:,1:galaxyCount-1)=galaxyPositionTmp
          galaxyVelocity(:,1:galaxyCount-1)=galaxyVelocityTmp
          call Dealloc_Array(galaxyMassTmp    )
          call Dealloc_Array(galaxyPositionTmp)
          call Dealloc_Array(galaxyVelocityTmp)
       end if
    else
       call Alloc_Array(galaxyMass    ,[  galaxyBufferSizeMinimum])
       call Alloc_Array(galaxyPosition,[3,galaxyBufferSizeMinimum])
       call Alloc_Array(galaxyVelocity,[3,galaxyBufferSizeMinimum])
    end if
    ! Store the galaxy.
    galaxyMass    (  galaxyCount)=mass
    galaxyPosition(:,galaxyCount)=position
    galaxyVelocity(:,galaxyCount)=velocity
    ! Record the lowest mass halo populated.
    populatedHaloMassMinimum=min(haloMass(iHalo),populatedHaloMassMinimum)           
    return
  end subroutine Galaxy_Add

end program Halo_Model_Mock
