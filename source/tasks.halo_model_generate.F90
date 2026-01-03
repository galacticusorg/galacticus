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

  use :: Conditional_Mass_Functions, only : conditionalMassFunction     , conditionalMassFunctionClass
  use :: Cosmology_Functions       , only : cosmologyFunctions          , cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParameters         , cosmologyParametersClass
  use :: Dark_Matter_Profile_Scales, only : darkMatterProfileScaleRadius, darkMatterProfileScaleRadiusClass
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMO        , darkMatterProfileDMOClass
  use :: Numerical_Random_Numbers  , only : randomNumberGeneratorClass

  !![
  <task name="taskHaloModelGenerate">
   <description>A task which generates a mock catalog of galaxies based on a simple halo model approach.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskHaloModelGenerate
     !!{
     Implementation of a task which generates a mock catalog of galaxies based on a simple halo model approach.
     !!}
     private
     class           (cosmologyParametersClass         ), pointer :: cosmologyParameters_          => null()
     class           (cosmologyFunctionsClass          ), pointer :: cosmologyFunctions_           => null()
     class           (darkMatterProfileDMOClass        ), pointer :: darkMatterProfileDMO_         => null()
     class           (conditionalMassFunctionClass     ), pointer :: conditionalMassFunction_      => null()
     class           (darkMatterProfileScaleRadiusClass), pointer :: darkMatterProfileScaleRadius_ => null()
     class           (randomNumberGeneratorClass       ), pointer :: randomNumberGenerator_        => null()
     double precision                                             :: massMinimum                             , massMaximum
     type            (varying_string                   )          :: galaxyCatalogFileName                   , haloCatalogFileName
     ! Pointer to the parameters for this task.
     type            (inputParameters                  )          :: parameters
     logical                                                      :: nodeComponentsInitialized     =  .false.
   contains
     final     ::                       haloModelGenerateDestructor
     procedure :: perform            => haloModelGeneratePerform
     procedure :: requiresOutputFile => haloModelGenerateRequiresOutputFile
  end type taskHaloModelGenerate

  interface taskHaloModelGenerate
     !!{
     Constructors for the \refClass{taskHaloModelGenerate} task.
     !!}
     module procedure haloModelGenerateConstructorParameters
     module procedure haloModelGenerateConstructorInternal
  end interface taskHaloModelGenerate

contains

  function haloModelGenerateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskHaloModelGenerate} task class which takes a parameter set as input.
    !!}
    use :: Galacticus_Nodes, only : nodeClassHierarchyInitialize
    use :: Input_Parameters, only : inputParameter              , inputParameters
    use :: Node_Components , only : Node_Components_Initialize
    implicit none
    type            (taskHaloModelGenerate            )                        :: self
    type            (inputParameters                  ), intent(inout), target :: parameters
    class           (cosmologyFunctionsClass          ), pointer               :: cosmologyFunctions_
    class           (cosmologyParametersClass         ), pointer               :: cosmologyParameters_
    class           (darkMatterProfileDMOClass        ), pointer               :: darkMatterProfileDMO_
    class           (conditionalMassFunctionClass     ), pointer               :: conditionalMassFunction_
    class           (darkMatterProfileScaleRadiusClass), pointer               :: darkMatterProfileScaleRadius_
    class           (randomNumberGeneratorClass       ), pointer               :: randomNumberGenerator_
    type            (inputParameters                  ), pointer               :: parametersRoot
    double precision                                                           :: massMinimum                  , massMaximum
    type            (varying_string                   )                        :: galaxyCatalogFileName        , haloCatalogFileName

    ! Ensure the nodes objects are initialized.
    if (associated(parameters%parent)) then
       parametersRoot => parameters%parent
       do while (associated(parametersRoot%parent))
          parametersRoot => parametersRoot%parent
       end do
       call nodeClassHierarchyInitialize(parametersRoot)
       call Node_Components_Initialize  (parametersRoot)
    else
       parametersRoot => parameters
       call nodeClassHierarchyInitialize(parameters    )
       call Node_Components_Initialize  (parameters    )
    end if
    self%nodeComponentsInitialized=.false.
    !![
    <inputParameter>
      <name>haloCatalogFileName</name>
      <description>The file name of the halo catalog to populate.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>galaxyCatalogFileName</name>
      <description>The file name to which the galaxy catalog should be output.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMinimum</name>
      <description>The minimum mass galaxy to include in a mock halo model realization.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <defaultValue>1.0d16</defaultValue>
      <description>The maximum mass galaxy to include in a mock halo model realization.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"          name="cosmologyParameters_"          source="parameters"/>
    <objectBuilder class="cosmologyFunctions"           name="cosmologyFunctions_"           source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"         name="darkMatterProfileDMO_"         source="parameters"/>
    <objectBuilder class="conditionalMassFunction"      name="conditionalMassFunction_"      source="parameters"/>
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    <objectBuilder class="randomNumberGenerator"        name="randomNumberGenerator_"        source="parameters"/>
    !!]
    self=taskHaloModelGenerate(galaxyCatalogFileName,haloCatalogFileName,massMinimum,massMaximum,cosmologyParameters_,cosmologyFunctions_,darkMatterProfileDMO_,conditionalMassFunction_,darkMatterProfileScaleRadius_,randomNumberGenerator_,parametersRoot)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"         />
    <objectDestructor name="cosmologyFunctions_"          />
    <objectDestructor name="darkMatterProfileDMO_"        />
    <objectDestructor name="conditionalMassFunction_"     />
    <objectDestructor name="darkMatterProfileScaleRadius_"/>
    <objectDestructor name="randomNumberGenerator_"       />
    !!]
    return
  end function haloModelGenerateConstructorParameters

  function haloModelGenerateConstructorInternal(galaxyCatalogFileName,haloCatalogFileName,massMinimum,massMaximum,cosmologyParameters_,cosmologyFunctions_,darkMatterProfileDMO_,conditionalMassFunction_,darkMatterProfileScaleRadius_,randomNumberGenerator_,parameters) result(self)
    !!{
    Constructor for the \refClass{taskHaloModelGenerate} task class which takes a parameter set as input.
    !!}
    implicit none
    type            (taskHaloModelGenerate            )                        :: self
    type            (varying_string                   ), intent(in   )         :: galaxyCatalogFileName        , haloCatalogFileName
    double precision                                   , intent(in   )         :: massMinimum                  , massMaximum
    class           (cosmologyParametersClass         ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass          ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass        ), intent(in   ), target :: darkMatterProfileDMO_
    class           (conditionalMassFunctionClass     ), intent(in   ), target :: conditionalMassFunction_
    class           (darkMatterProfileScaleRadiusClass), intent(in   ), target :: darkMatterProfileScaleRadius_
    class           (randomNumberGeneratorClass       ), intent(in   ), target :: randomNumberGenerator_
    type            (inputParameters                  ), intent(in   ), target :: parameters
    !![
    <constructorAssign variables="galaxyCatalogFileName, haloCatalogFileName, massMinimum, massMaximum, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterProfileDMO_, *conditionalMassFunction_, *darkMatterProfileScaleRadius_, *randomNumberGenerator_"/>
    !!]

    self%parameters=inputParameters(parameters)
    call self%parameters%parametersGroupCopy(parameters)
    return
  end function haloModelGenerateConstructorInternal

  subroutine haloModelGenerateDestructor(self)
    !!{
    Destructor for the \refClass{taskHaloModelGenerate} task class.
    !!}
    use :: Node_Components, only : Node_Components_Uninitialize
    implicit none
    type(taskHaloModelGenerate), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"         />
    <objectDestructor name="self%cosmologyFunctions_"          />
    <objectDestructor name="self%darkMatterProfileDMO_"        />
    <objectDestructor name="self%conditionalMassFunction_"     />
    <objectDestructor name="self%darkMatterProfileScaleRadius_"/>
    <objectDestructor name="self%randomNumberGenerator_"       />
    !!]
    if(self%nodeComponentsInitialized) call Node_Components_Uninitialize()
    return
  end subroutine haloModelGenerateDestructor

  subroutine haloModelGeneratePerform(self,status)
    !!{
    Generate a mock galaxy catalog using a simple halo model approach.
    !!}
    use :: Conditional_Mass_Functions, only : haloModelGalaxyTypeCentral       , haloModelGalaxyTypeSatellite
    use :: Display                   , only : displayIndent                    , displayMessage                     , displayUnindent
    use :: Galactic_Structure_Options, only : massTypeDark
    use :: Calculations_Resets       , only : Calculations_Reset
    use :: Error                     , only : errorStatusSuccess
    use :: Galacticus_Nodes          , only : nodeComponentBasic               , nodeComponentDarkMatterProfile     , treeNode
    use :: IO_IRATE                  , only : irate
    use :: ISO_Varying_String        , only : varying_string
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Node_Components           , only : Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize
    use :: Numerical_Constants_Math  , only : Pi
    use :: Root_Finder               , only : rangeExpandMultiplicative        , rootFinder
    use :: String_Handling           , only : operator(//)
    implicit none
    class           (taskHaloModelGenerate         ), intent(inout), target         :: self
    integer                                         , intent(  out), optional       :: status
    double precision                                , pointer      , dimension(  :) :: massHalo
    double precision                                , allocatable  , dimension(  :) :: massGalaxy
    double precision                                , pointer      , dimension(:,:) :: positionHalo         , velocityHalo
    double precision                                , allocatable  , dimension(:,:) :: positionGalaxy       , velocityGalaxy
    double precision                                               , dimension(3  ) :: positionSatellite    , velocitySatellite
    type            (treeNode                      ), pointer                       :: node
    class           (nodeComponentBasic            ), pointer                       :: basic
    class           (nodeComponentDarkMatterProfile), pointer                       :: profile
    class           (massDistributionClass         ), pointer                       :: massDistribution_
    integer                                                                         :: iHalo                , galaxyCount              , &
         &                                                                             numberSatelliteActual, iSatellite               , &
         &                                                                             iAxis
    type            (varying_string                )                                :: message
    type            (irate                         )                                :: haloFile             , galaxyFile
    double precision                                                                :: probabilityCentral   , xCentral                 , &
         &                                                                             massGalaxy_          , numberSatelliteMean      , &
         &                                                                             xSatellite           , redshift                 , &
         &                                                                             radiusSatellite      , thetaSatellite           , &
         &                                                                             phiSatellite         , velocitySatelliteCircular, &
         &                                                                             simulationBoxSize    , populatedHaloMassMinimum
    type            (rootFinder                    )                                :: finderCentral        , finderSatellite
    character       (len=6                         )                                :: label

    call displayIndent('Begin task: halo model generate')
    ! Call routines to perform initialization which must occur for all threads if run in parallel.
    call Node_Components_Thread_Initialize(self%parameters)
    ! Read the halo catalog.
    call displayIndent("Reading halo catalog")
    haloFile=irate(char(self%haloCatalogFileName),self%cosmologyParameters_,self%cosmologyFunctions_)
    call haloFile%readHalos     (                            &
         &                       snapshot=1                , &
         &                       redshift=redshift         , &
         &                       center  =positionHalo     , &
         &                       velocity=velocityHalo     , &
         &                       mass    =massHalo           &
         &                      )
    call haloFile%readSimulation(                            &
         &                       boxSize =simulationBoxSize  &
         &                      )
    call displayUnindent("done")
    ! Establish root finders.
    finderCentral  =rootFinder(                                               &
         &                     rootFunction       =centralMassRoot          , &
         &                     toleranceRelative  =1.0d-6                   , &
         &                     rangeExpandUpward  =2.0d0                    , &
         &                     rangeExpandDownward=1.0d0                    , &
         &                     rangeExpandType    =rangeExpandMultiplicative  &
         &                    )
    finderSatellite=rootFinder(                                               &
         &                     rootFunction       =satelliteMassRoot        , &
         &                     toleranceRelative  =1.0d-6                   , &
         &                     rangeExpandUpward  =2.0d0                    , &
         &                     rangeExpandDownward=1.0d0                    , &
         &                     rangeExpandType    =rangeExpandMultiplicative  &
         &                    )
     ! Establish tree node for sampling satellite positions.
    node    => treeNode                  (                 )
    basic   => node    %basic            (autoCreate=.true.)
    profile => node    %darkMatterProfile(autoCreate=.true.)
    call basic%timeSet            (self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshift)))
    call basic%timeLastIsolatedSet(self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshift)))
    ! Iterate over halos.
    call displayIndent("Populating halos")
    galaxyCount             =     0
    populatedHaloMassMinimum=huge(0.0d0)
    do iHalo=1,size(massHalo)
       ! Get probability of central galaxy.
       probabilityCentral=+self%conditionalMassFunction_%massFunction(massHalo(iHalo),self%massMinimum,haloModelGalaxyTypeCentral) &
            &             -self%conditionalMassFunction_%massFunction(massHalo(iHalo),self%massMaximum,haloModelGalaxyTypeCentral)
       ! Test for inclusion.
       if (self%randomNumberGenerator_%uniformSample() <= probabilityCentral) then
          ! Sample central galaxy mass.
          xCentral   =self%randomNumberGenerator_%uniformSample()*probabilityCentral
          massGalaxy_=finderCentral%find(rootGuess=self%massMinimum)
          call galaxyAdd(massGalaxy_,positionHalo(:,iHalo),velocityHalo(:,iHalo))
       end if
       ! Get mean number of satellite galaxies.
       numberSatelliteMean  =+self%conditionalMassFunction_%massFunction(massHalo(iHalo),self%massMinimum,haloModelGalaxyTypeSatellite) &
            &                -self%conditionalMassFunction_%massFunction(massHalo(iHalo),self%massMaximum,haloModelGalaxyTypeSatellite)
       numberSatelliteActual=+self%randomNumberGenerator_%poissonSample(numberSatelliteMean)
       if (numberSatelliteActual > 0) then
          ! Construct the dark matter halo profile.
          call basic  %massSet   (massHalo                                 (iHalo))
          call Calculations_Reset(                                          node  )
          call profile%scaleSet  (self%darkMatterProfileScaleRadius_%radius(node ))
          call Calculations_Reset(                                          node  )
          massDistribution_ => self%darkMatterProfileDMO_%get(node)
          do iSatellite=1,numberSatelliteActual
             ! Sample satellite galaxy mass.
             xSatellite                =  self%randomNumberGenerator_%uniformSample()*numberSatelliteMean
             massGalaxy_               =  finderSatellite%find(rootGuess=self%massMinimum)
             ! Sample galaxy radial position.
             xSatellite                =  self             %randomNumberGenerator_%uniformSample      (                           )
             radiusSatellite           =  massDistribution_                       %radiusEnclosingMass(massFractional=xSatellite  )
             ! Get circular velocity at this radius.
             velocitySatelliteCircular =  massDistribution_%rotationCurve(radiusSatellite)
             ! Convert radial position to comoving coordinates.
             radiusSatellite          =radiusSatellite*(1.0d0+redshift)
             ! Sample galaxy angular position.
             phiSatellite             =     self%randomNumberGenerator_%uniformSample()*2.0d0*Pi
             thetaSatellite           =acos(self%randomNumberGenerator_%uniformSample()*2.0d0-1.0d0)
             ! Set satellite position.
             positionSatellite        =+radiusSatellite                                                         &
                  &                    *[                                                                       &
                  &                      sin(thetaSatellite)*sin(phiSatellite)                                , &
                  &                      sin(thetaSatellite)*cos(phiSatellite)                                , &
                  &                      cos(thetaSatellite)                                                    &
                  &                     ]                                                                       &
                  &                    +positionHalo(:,iHalo)
             ! Periodicalize the satellite position.
             do iAxis=1,3
                do while (positionSatellite(iAxis) < 0.0d0)
                   positionSatellite(iAxis)=positionSatellite(iAxis)+simulationBoxSize
                end do
                do while (positionSatellite(iAxis) > simulationBoxSize)
                   positionSatellite(iAxis)=positionSatellite(iAxis)-simulationBoxSize
                end do
             end do
             ! Set satellite velocity.
             velocitySatellite        =+velocitySatelliteCircular                            &
                  &                    /sqrt(3.0d0)                                          &
                  &                    *[                                                    &
                  &                      self%randomNumberGenerator_%standardNormalSample(), &
                  &                      self%randomNumberGenerator_%standardNormalSample(), &
                  &                      self%randomNumberGenerator_%standardNormalSample()  &
                  &                    ]                                                     &
                  &                    +velocityHalo(:,iHalo)
             ! Store the satellite.
             call galaxyAdd(massGalaxy_,positionSatellite,velocitySatellite)
          end do
          !![
	  <objectDestructor name="massDistribution_"/>
          !!]
       end if
    end do
    message="Created "
    message=message//galaxyCount//" galaxies"
    call displayMessage(message)
    write (label,'(f5.2)') log10(populatedHaloMassMinimum)
    message="Lowest mass halo populated has log₁₀(Mₕₐₗₒ/M☉)="//trim(adjustl(label))
    call displayMessage(message)
    call displayUnindent("done")
    ! Output galaxy catalog.
    galaxyFile=irate(char(self%galaxyCatalogFileName),self%cosmologyParameters_,self%cosmologyFunctions_)
    call haloFile%copyCosmology (galaxyFile)
    call haloFile%copySimulation(galaxyFile)
    call galaxyFile%writeHalos(1,redshift,positionGalaxy(:,1:galaxyCount),velocityGalaxy(:,1:galaxyCount),massGalaxy(1:galaxyCount))
    ! Clean up.
    deallocate(positionHalo)
    deallocate(velocityHalo)
    deallocate(massHalo    )
    call Node_Components_Thread_Uninitialize()
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: halo model generate' )

  contains

    double precision function centralMassRoot(mass)
      !!{
      Root function used to find the mass of central galaxies
      !!}
      implicit none
      double precision, intent(in   ) :: mass
      double precision                :: x

      x=      (                                                                                                         &
           &   +self%conditionalMassFunction_%massFunction(massHalo(iHalo),     mass       ,haloModelGalaxyTypeCentral) &
           &   -self%conditionalMassFunction_%massFunction(massHalo(iHalo),self%massMaximum,haloModelGalaxyTypeCentral) &
           &  )
      centralMassRoot=x-xCentral
      return
    end function centralMassRoot

    double precision function satelliteMassRoot(mass)
      !!{
      Root function used to find the mass of satellite galaxies
      !!}
      implicit none
      double precision, intent(in   ) :: mass
      double precision                :: x

      x=      (                                                                                                           &
           &   +self%conditionalMassFunction_%massFunction(massHalo(iHalo),     mass       ,haloModelGalaxyTypeSatellite) &
           &   -self%conditionalMassFunction_%massFunction(massHalo(iHalo),self%massMaximum,haloModelGalaxyTypeSatellite) &
           &  )
      satelliteMassRoot=x-xSatellite
      return
    end function satelliteMassRoot

    subroutine galaxyAdd(mass,position,velocity)
      !!{
      Add a galaxy to the output buffers.
      !!}
      implicit none
      double precision, intent(in   )                 :: mass
      double precision, intent(in   ), dimension(3  ) :: position                     , velocity
      integer         , parameter                     :: galaxyBufferSizeMinimum=10000
      double precision, allocatable  , dimension(  :) :: massGalaxyTmp
      double precision, allocatable  , dimension(:,:) :: positionGalaxyTmp            , velocityGalaxyTmp

      ! Expand output buffers as needed.
      galaxyCount=galaxyCount+1
      if (allocated(massGalaxy)) then
         if (galaxyCount > size(massGalaxy)) then
            call move_alloc(massGalaxy    ,massGalaxyTmp    )
            call move_alloc(positionGalaxy,positionGalaxyTmp)
            call move_alloc(velocityGalaxy,velocityGalaxyTmp)
            allocate(massGalaxy    (  2*(galaxyCount-1)))
            allocate(positionGalaxy(3,2*(galaxyCount-1)))
            allocate(velocityGalaxy(3,2*(galaxyCount-1)))
            massGalaxy    (  1:galaxyCount-1)=massGalaxyTmp
            positionGalaxy(:,1:galaxyCount-1)=positionGalaxyTmp
            velocityGalaxy(:,1:galaxyCount-1)=velocityGalaxyTmp
            deallocate(massGalaxyTmp    )
            deallocate(positionGalaxyTmp)
            deallocate(velocityGalaxyTmp)
         end if
      else
         allocate(massGalaxy    (galaxyBufferSizeMinimum))
         allocate(positionGalaxy(3,galaxyBufferSizeMinimum))
         allocate(velocityGalaxy(3,galaxyBufferSizeMinimum))
      end if
      ! Store the galaxy.
      massGalaxy    (  galaxyCount)=mass
      positionGalaxy(:,galaxyCount)=position
      velocityGalaxy(:,galaxyCount)=velocity
      ! Record the lowest mass halo populated.
      populatedHaloMassMinimum=min(massHalo(iHalo),populatedHaloMassMinimum)
      return
    end subroutine galaxyAdd

  end subroutine haloModelGeneratePerform

  logical function haloModelGenerateRequiresOutputFile(self)
    !!{
    Specifies that this task does not requires the main output file.
    !!}
    implicit none
    class(taskHaloModelGenerate), intent(inout) :: self
    !$GLC attributes unused :: self

    haloModelGenerateRequiresOutputFile=.false.
    return
  end function haloModelGenerateRequiresOutputFile
