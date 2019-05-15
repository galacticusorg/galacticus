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

  use Cosmology_Functions       , only : cosmologyFunctions          , cosmologyFunctionsClass
  use Cosmology_Parameters      , only : cosmologyParameters         , cosmologyParametersClass
  use Conditional_Mass_Functions, only : conditionalMassFunction     , conditionalMassFunctionClass
  use Dark_Matter_Profile_Scales, only : darkMatterProfileScaleRadius, darkMatterProfileScaleRadiusClass
  use Dark_Matter_Profiles_DMO      , only : darkMatterProfileDMO           , darkMatterProfileDMOClass

  !# <task name="taskHaloModelGenerate">
  !#  <description>A task which generates a mock catalog of galaxies based on a simple halo model approach.</description>
  !# </task>
  type, extends(taskClass) :: taskHaloModelGenerate
     !% Implementation of a task which generates a mock catalog of galaxies based on a simple halo model approach.
     private 
     class           (cosmologyParametersClass         ), pointer :: cosmologyParameters_          => null()
     class           (cosmologyFunctionsClass          ), pointer :: cosmologyFunctions_           => null()
     class           (darkMatterProfileDMOClass        ), pointer :: darkMatterProfileDMO_         => null()
     class           (conditionalMassFunctionClass     ), pointer :: conditionalMassFunction_      => null()
     class           (darkMatterProfileScaleRadiusClass), pointer :: darkMatterProfileScaleRadius_ => null()
     double precision                                             :: massMinimum                            , massMaximum
     type            (varying_string                   )          :: galaxyCatalogFileName                  , haloCatalogFileName
   contains
     final     ::                       haloModelGenerateDestructor
     procedure :: perform            => haloModelGeneratePerform
     procedure :: requiresOutputFile => haloModelGenerateRequiresOutputFile
  end type taskHaloModelGenerate

  interface taskHaloModelGenerate
     !% Constructors for the {\normalfont \ttfamily haloModelGenerate} task.
     module procedure haloModelGenerateConstructorParameters
     module procedure haloModelGenerateConstructorInternal
  end interface taskHaloModelGenerate

contains
  
  function haloModelGenerateConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily haloModelGenerate} task class which takes a parameter set as input.
    use Input_Parameters, only : inputParameters             , inputParameter
    use Galacticus_Nodes, only : nodeClassHierarchyInitialize
    use Node_Components , only : Node_Components_Initialize  , Node_Components_Thread_Initialize
    implicit none
    type            (taskHaloModelGenerate            )                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass         ), pointer       :: cosmologyParameters_
    class           (darkMatterProfileDMOClass        ), pointer       :: darkMatterProfileDMO_
    class           (conditionalMassFunctionClass     ), pointer       :: conditionalMassFunction_
    class           (darkMatterProfileScaleRadiusClass), pointer       :: darkMatterProfileScaleRadius_
    type            (inputParameters                  ), pointer       :: parametersRoot
    double precision                                                   :: massMinimum                  , massMaximum
    type            (varying_string                   )                :: galaxyCatalogFileName        , haloCatalogFileName

    ! Ensure the nodes objects are initialized.
    if (associated(parameters%parent)) then
       parametersRoot => parameters%parent
       do while (associated(parametersRoot%parent))
          parametersRoot => parametersRoot%parent
       end do
       call nodeClassHierarchyInitialize     (parametersRoot)
       call Node_Components_Initialize       (parametersRoot)
       call Node_Components_Thread_Initialize(parametersRoot)
    else
       parametersRoot => null()
       call nodeClassHierarchyInitialize     (parameters    )
       call Node_Components_Initialize       (parameters    )
       call Node_Components_Thread_Initialize(parameters    )
    end if
    !# <inputParameter>
    !#   <name>haloCatalogFileName</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The file name of the halo catalog to populate.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>galaxyCatalogFileName</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The file name to which the galaxy catalog should be output.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The minimum mass galaxy to include in a mock halo model realization.</description>
    !#   <source>parameters</source>
    !#   <type>float</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massMaximum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d16</defaultValue>
    !#   <description>The maximum mass galaxy to include in a mock halo model realization.</description>
    !#   <source>parameters</source>
    !#   <type>float</type>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters"          name="cosmologyParameters_"          source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"           name="cosmologyFunctions_"           source="parameters"/>
    !# <objectBuilder class="darkMatterProfileDMO"         name="darkMatterProfileDMO_"         source="parameters"/>
    !# <objectBuilder class="conditionalMassFunction"      name="conditionalMassFunction_"      source="parameters"/>
    !# <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    self=taskHaloModelGenerate(galaxyCatalogFileName,haloCatalogFileName,massMinimum,massMaximum,cosmologyParameters_,cosmologyFunctions_,darkMatterProfileDMO_,conditionalMassFunction_,darkMatterProfileScaleRadius_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"         />
    !# <objectDestructor name="cosmologyFunctions_"          />
    !# <objectDestructor name="darkMatterProfileDMO_"        />
    !# <objectDestructor name="conditionalMassFunction_"     />
    !# <objectDestructor name="darkMatterProfileScaleRadius_"/>
    return
  end function haloModelGenerateConstructorParameters
  
  function haloModelGenerateConstructorInternal(galaxyCatalogFileName,haloCatalogFileName,massMinimum,massMaximum,cosmologyParameters_,cosmologyFunctions_,darkMatterProfileDMO_,conditionalMassFunction_,darkMatterProfileScaleRadius_) result(self)
    !% Constructor for the {\normalfont \ttfamily haloModelGenerate} task class which takes a parameter set as input.
    implicit none
    type            (taskHaloModelGenerate            )                        :: self
    type            (varying_string                   ), intent(in   )         :: galaxyCatalogFileName        , haloCatalogFileName
    double precision                                   , intent(in   )         :: massMinimum                  , massMaximum
    class           (cosmologyParametersClass         ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass          ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass        ), intent(in   ), target :: darkMatterProfileDMO_
    class           (conditionalMassFunctionClass     ), intent(in   ), target :: conditionalMassFunction_
    class           (darkMatterProfileScaleRadiusClass), intent(in   ), target :: darkMatterProfileScaleRadius_
    !# <constructorAssign variables="galaxyCatalogFileName, haloCatalogFileName, massMinimum, massMaximum, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterProfileDMO_, *conditionalMassFunction_, *darkMatterProfileScaleRadius_"/>

    return
  end function haloModelGenerateConstructorInternal

  subroutine haloModelGenerateDestructor(self)
    !% Destructor for the {\normalfont \ttfamily haloModelGenerate} task class.
    use Node_Components, only : Node_Components_Uninitialize, Node_Components_Thread_Uninitialize
    implicit none
    type(taskHaloModelGenerate), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"         />
    !# <objectDestructor name="self%cosmologyFunctions_"          />
    !# <objectDestructor name="self%darkMatterProfileDMO_"        />
    !# <objectDestructor name="self%conditionalMassFunction_"     />
    !# <objectDestructor name="self%darkMatterProfileScaleRadius_"/>
    call Node_Components_Uninitialize       ()
    call Node_Components_Thread_Uninitialize()
    return
  end subroutine haloModelGenerateDestructor

  subroutine haloModelGeneratePerform(self,status)
    !% Generate a mock galaxy catalog using a simple halo model approach.
    use ISO_Varying_String
    use Memory_Management
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Galacticus_Error
    use Galacticus_Display
    use Geometry_Surveys
    use IO_HDF5
    use IO_IRATE
    use Pseudo_Random
    use Root_Finder
    use String_Handling
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Galacticus_Calculations_Resets
    use Galacticus_Nodes                   , only : treeNode                  , nodeComponentBasic          , nodeComponentDarkMatterProfile
    use Conditional_Mass_Functions         , only : haloModelGalaxyTypeCentral, haloModelGalaxyTypeSatellite
    implicit none
    class           (taskHaloModelGenerate         ), intent(inout)                 :: self
    integer                                         , intent(  out), optional       :: status
    double precision                                , allocatable  , dimension(  :) :: haloMass             , galaxyMass
    double precision                                , allocatable  , dimension(:,:) :: haloPosition         , haloVelocity             , &
         &                                                                             galaxyPosition       , galaxyVelocity
    double precision                                               , dimension(3  ) :: satellitePosition    , satelliteVelocity
    type            (treeNode                      ), pointer                       :: node
    class           (nodeComponentBasic            ), pointer                       :: basic
    class           (nodeComponentDarkMatterProfile), pointer                       :: profile
    integer                                                                         :: iHalo                , galaxyCount              , &
         &                                                                             satelliteNumberActual, iSatellite               , &
         &                                                                             iAxis
    type            (varying_string                )                                :: message
    type            (irate                         )                                :: haloFile             , galaxyFile
    double precision                                                                :: probabilityCentral   , xCentral                 , &
         &                                                                             massGalaxy           , satelliteNumberMean      , &
         &                                                                             xSatellite           , redshift                 , &
         &                                                                             satelliteRadius      , satelliteTheta           , &
         &                                                                             satellitePhi         , satelliteVelocityCircular, &
         &                                                                             simulationBoxSize    , populatedHaloMassMinimum
    type            (pseudoRandom                  )                                :: randomSequence
    type            (rootFinder                    )                                :: finderCentral        , finderSatellite
    character       (len=6                         )                                :: label

    call Galacticus_Display_Indent('Begin task: halo model generate')
    ! Read the halo catalog.
    call Galacticus_Display_Indent("Reading halo catalog")
    haloFile=irate(char(self%haloCatalogFileName),self%cosmologyParameters_,self%cosmologyFunctions_)
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
    call basic%timeSet            (self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshift)))
    call basic%timeLastIsolatedSet(self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshift)))
    ! Iterate over halos.
    call Galacticus_Display_Indent("Populating halos")
    galaxyCount             =     0
    populatedHaloMassMinimum=huge(0.0d0)
    do iHalo=1,size(haloMass)
       ! Get probability of central galaxy.
       probabilityCentral=+self%conditionalMassFunction_%massFunction(haloMass(iHalo),self%massMinimum,haloModelGalaxyTypeCentral) &
            &             -self%conditionalMassFunction_%massFunction(haloMass(iHalo),self%massMaximum,haloModelGalaxyTypeCentral)
       ! Test for inclusion.
       if (randomSequence%uniformSample() <= probabilityCentral) then
          ! Sample central galaxy mass.        
          xCentral  =randomSequence%uniformSample()*probabilityCentral
          massGalaxy=finderCentral%find(rootGuess=self%massMinimum)
          call galaxyAdd(massGalaxy,haloPosition(:,iHalo),haloVelocity(:,iHalo))
       end if
       ! Get mean number of satellite galaxies.
       satelliteNumberMean  =+self%conditionalMassFunction_%massFunction(haloMass(iHalo),self%massMinimum,haloModelGalaxyTypeSatellite) &
            &                -self%conditionalMassFunction_%massFunction(haloMass(iHalo),self%massMaximum,haloModelGalaxyTypeSatellite)
       satelliteNumberActual=+randomSequence%poissonSample(satelliteNumberMean)
       if (satelliteNumberActual > 0) then
          ! Construct the dark matter halo profile.
          call basic  %massSet              (haloMass                                 (iHalo))
          call Galacticus_Calculations_Reset(                                          node  )
          call profile%scaleSet             (self%darkMatterProfileScaleRadius_%radius(node ))
          call Galacticus_Calculations_Reset(                                          node  )
          do iSatellite=1,satelliteNumberActual
             ! Sample satellite galaxy mass.        
             xSatellite               =     randomSequence%uniformSample()*satelliteNumberMean
             massGalaxy               =finderSatellite%find(rootGuess=self%massMinimum)
             ! Sample galaxy radial position.
             xSatellite               =     randomSequence%uniformSample()
             satelliteRadius          =Galactic_Structure_Radius_Enclosing_Mass(node,fractionalMass=xSatellite,massType=massTypeDark)
             ! Get circular velocity at this radius.
             satelliteVelocityCircular=self%darkMatterProfileDMO_%circularVelocity(node,satelliteRadius)
             ! Convert radial position to comoving coordinates.
             satelliteRadius          =satelliteRadius*(1.0d0+redshift)
             ! Sample galaxy angular position.
             satellitePhi             =     randomSequence%uniformSample()*2.0d0*Pi
             satelliteTheta           =acos(randomSequence%uniformSample()*2.0d0-1.0d0)
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
             satelliteVelocity        =+satelliteVelocityCircular       &
                  &                    /sqrt(3.0d0)                     &
                  &                    *[                               &
                  &                      randomSequence%normalSample(), &
                  &                      randomSequence%normalSample(), &
                  &                      randomSequence%normalSample()  &
                  &                    ]                                &
                  &                    +haloVelocity(:,iHalo)
             ! Store the satellite.
             call galaxyAdd(massGalaxy,satellitePosition,satelliteVelocity)
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
    galaxyFile=irate(char(self%galaxyCatalogFileName),self%cosmologyParameters_,self%cosmologyFunctions_)
    call haloFile%copyCosmology (galaxyFile)
    call haloFile%copySimulation(galaxyFile)
    call galaxyFile%writeHalos(1,redshift,galaxyPosition(:,1:galaxyCount),galaxyVelocity(:,1:galaxyCount),galaxyMass(1:galaxyCount))
    if (present(status)) status=errorStatusSuccess
    call Galacticus_Display_Unindent('Done task: halo model generate' )

  contains

    double precision function centralMassRoot(mass)
      !% Root function used to find the mass of central galaxies
      implicit none
      double precision, intent(in   ) :: mass
      double precision                :: x

      x=      (                                                                                                         &
           &   +self%conditionalMassFunction_%massFunction(haloMass(iHalo),     mass       ,haloModelGalaxyTypeCentral) &
           &   -self%conditionalMassFunction_%massFunction(haloMass(iHalo),self%massMaximum,haloModelGalaxyTypeCentral) &
           &  )
      centralMassRoot=x-xCentral
      return
    end function centralMassRoot

    double precision function satelliteMassRoot(mass)
      !% Root function used to find the mass of satellite galaxies
      implicit none
      double precision, intent(in   ) :: mass
      double precision                :: x

      x=      (                                                                                                           &
           &   +self%conditionalMassFunction_%massFunction(haloMass(iHalo),     mass       ,haloModelGalaxyTypeSatellite) &
           &   -self%conditionalMassFunction_%massFunction(haloMass(iHalo),self%massMaximum,haloModelGalaxyTypeSatellite) &
           &  )
      satelliteMassRoot=x-xSatellite
      return
    end function satelliteMassRoot

    subroutine galaxyAdd(mass,position,velocity)
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
            call allocateArray(galaxyMass    ,[  2*(galaxyCount-1)])
            call allocateArray(galaxyPosition,[3,2*(galaxyCount-1)])
            call allocateArray(galaxyVelocity,[3,2*(galaxyCount-1)])
            galaxyMass    (  1:galaxyCount-1)=galaxyMassTmp
            galaxyPosition(:,1:galaxyCount-1)=galaxyPositionTmp
            galaxyVelocity(:,1:galaxyCount-1)=galaxyVelocityTmp
            call deallocateArray(galaxyMassTmp    )
            call deallocateArray(galaxyPositionTmp)
            call deallocateArray(galaxyVelocityTmp)
         end if
      else
         call allocateArray(galaxyMass    ,[  galaxyBufferSizeMinimum])
         call allocateArray(galaxyPosition,[3,galaxyBufferSizeMinimum])
         call allocateArray(galaxyVelocity,[3,galaxyBufferSizeMinimum])
      end if
      ! Store the galaxy.
      galaxyMass    (  galaxyCount)=mass
      galaxyPosition(:,galaxyCount)=position
      galaxyVelocity(:,galaxyCount)=velocity
      ! Record the lowest mass halo populated.
      populatedHaloMassMinimum=min(haloMass(iHalo),populatedHaloMassMinimum)           
      return
    end subroutine galaxyAdd

  end subroutine haloModelGeneratePerform

  logical function haloModelGenerateRequiresOutputFile(self)
    !% Specifies that this task does not requires the main output file.
    implicit none
    class(taskHaloModelGenerate), intent(inout) :: self    
    !GCC$ attributes unused :: self

    haloModelGenerateRequiresOutputFile=.false.
    return
  end function haloModelGenerateRequiresOutputFile
