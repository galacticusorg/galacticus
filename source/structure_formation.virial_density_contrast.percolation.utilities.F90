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
Contains a module of utilities needed by the {\normalfont \ttfamily percolation} virial density contrast class.
!!}

module Virial_Density_Contrast_Percolation_Utilities
  !!{
  Provides utilities needed by the {\normalfont \ttfamily percolation} virial density contrast class.
  !!}
  use :: Cosmology_Functions               , only : cosmologyFunctions            , cosmologyFunctionsClass
  use :: Cosmology_Parameters              , only : cosmologyParameters           , cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales           , only : darkMatterHaloScale           , darkMatterHaloScaleClass
  use :: Dark_Matter_Profile_Scales        , only : darkMatterProfileScaleRadius  , darkMatterProfileScaleRadiusClass  , darkMatterProfileScaleRadiusConcentration
  use :: Dark_Matter_Profiles_Concentration, only : darkMatterProfileConcentration, darkMatterProfileConcentrationClass
  use :: Dark_Matter_Profiles_DMO          , only : darkMatterProfileDMO          , darkMatterProfileDMOClass
  use :: Dark_Matter_Profiles_Shape        , only : darkMatterProfileShape        , darkMatterProfileShapeClass
  use :: Galacticus_Nodes                  , only : nodeComponentDarkMatterProfile, treeNode
  private
  public :: Virial_Density_Contrast_Percolation_Solver, Virial_Density_Contrast_Percolation_Objects_Constructor, percolationObjects, percolationObjectsDeepCopy, &
       &    percolationObjectsDeepCopyReset           , percolationObjectsDeepCopyFinalize

  ! Container type used to store state.
  type :: solverState
     type            (treeNode                                 ), pointer :: workNode                      => null()
     class           (nodeComponentDarkMatterProfile           ), pointer :: workDarkMatterProfile         => null()
     class           (darkMatterProfileDMOClass                ), pointer :: darkMatterProfileDMO_         => null()
     class           (darkMatterProfileShapeClass              ), pointer :: darkMatterProfileShape_       => null()
     type            (darkMatterProfileScaleRadiusConcentration), pointer :: darkMatterProfileScaleRadius_ => null()
     double precision                                                     :: boundingDensity                        , densityMatterMean, &
          &                                                                  massHalo
     double precision                                           , pointer :: densityContrast               => null()
  end type solverState

  ! State stack.
  integer                                         :: stateCount=0
  type   (solverState), allocatable, dimension(:) :: state
  !$omp threadprivate(state,stateCount)

  type :: percolationObjects
     !!{
     Type used to store pointers to objects
     !!}
     class(darkMatterProfileDMOClass          ), pointer :: darkMatterProfileDMO_           => null()
     class(cosmologyParametersClass           ), pointer :: cosmologyParameters_            => null()
     class(cosmologyFunctionsClass            ), pointer :: cosmologyFunctions_             => null()
     class(darkMatterHaloScaleClass           ), pointer :: darkMatterHaloScale_            => null()
     class(darkMatterProfileConcentrationClass), pointer :: darkMatterProfileConcentration_ => null()
     class(darkMatterProfileShapeClass        ), pointer :: darkMatterProfileShape_         => null()
   contains
     final :: percolationObjectsDestructor
  end type percolationObjects

contains

  !![
  <functionGlobal>
   <unitName>Virial_Density_Contrast_Percolation_Objects_Constructor</unitName>
   <type>class(*), pointer</type>
   <module>Input_Parameters, only : inputParameter, inputParameters</module>
   <arguments>type(inputParameters), intent(inout), target :: parameters</arguments>
  </functionGlobal>
  !!]
  function Virial_Density_Contrast_Percolation_Objects_Constructor(parameters) result(self)
    !!{
    Construct an instance of the container type for percolation virial density contrast objects from a parameter structure.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    class(*                 ), pointer               :: self
    type (percolationObjects), pointer               :: percolationObjects_
    type (inputParameters   ), intent(inout), target :: parameters

    allocate(percolationObjects_)
    !![
    <objectBuilder class="darkMatterProfileDMO"           name="percolationObjects_%darkMatterProfileDMO_"           source="parameters"/>
    <objectBuilder class="cosmologyParameters"            name="percolationObjects_%cosmologyParameters_"            source="parameters"/>
    <objectBuilder class="cosmologyFunctions"             name="percolationObjects_%cosmologyFunctions_"             source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"            name="percolationObjects_%darkMatterHaloScale_"            source="parameters"/>
    <objectBuilder class="darkMatterProfileConcentration" name="percolationObjects_%darkMatterProfileConcentration_" source="parameters"/>
    <objectBuilder class="darkMatterProfileShape"         name="percolationObjects_%darkMatterProfileShape_"         source="parameters"/>
    !!]
    self => percolationObjects_
    nullify(percolationObjects_)
    return
  end function Virial_Density_Contrast_Percolation_Objects_Constructor

  subroutine percolationObjectsDestructor(self)
    !!{
    Destruct an instance of the container type for percolation virial density contrast objects.
    !!}
    implicit none
    type(percolationObjects), intent(inout) :: self
    
    !![
    <objectDestructor name="self%darkMatterProfileDMO_"          />
    <objectDestructor name="self%cosmologyParameters_"           />
    <objectDestructor name="self%cosmologyFunctions_"            />
    <objectDestructor name="self%darkMatterHaloScale_"           />
    <objectDestructor name="self%darkMatterProfileConcentration_"/>
    <objectDestructor name="self%darkMatterProfileShape_"        />
    !!]
    return
  end subroutine percolationObjectsDestructor

  !![
  <functionGlobal>
   <unitName>percolationObjectsDeepCopyReset</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self</arguments>
  </functionGlobal>
  !!]
  subroutine percolationObjectsDeepCopyReset(self)
    !!{
    Perform a deep copy of percolation virial density contrast objects.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(*), intent(inout) :: self
    
    select type (self)
    class is (percolationObjects)
       if (associated(self%darkMatterProfileDMO_          )) call self%darkMatterProfileDMO_          %deepCopyReset()
       if (associated(self%cosmologyParameters_           )) call self%cosmologyParameters_           %deepCopyReset()
       if (associated(self%cosmologyFunctions_            )) call self%cosmologyFunctions_            %deepCopyReset()
       if (associated(self%darkMatterHaloScale_           )) call self%darkMatterHaloScale_           %deepCopyReset()
       if (associated(self%darkMatterProfileConcentration_)) call self%darkMatterProfileConcentration_%deepCopyReset()
       if (associated(self%darkMatterProfileShape_        )) call self%darkMatterProfileShape_        %deepCopyReset()
    class default
       call Error_Report("self must be of 'percolationObjects' class"//{introspection:location})
    end select
    return
  end subroutine percolationObjectsDeepCopyReset
  
  !![
  <functionGlobal>
   <unitName>percolationObjectsDeepCopyFinalize</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self</arguments>
  </functionGlobal>
  !!]
  subroutine percolationObjectsDeepCopyFinalize(self)
    !!{
    Finalize a deep copy of percolation virial density contrast objects.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(*), intent(inout) :: self
    
    select type (self)
    class is (percolationObjects)
       if (associated(self%darkMatterProfileDMO_          )) call self%darkMatterProfileDMO_          %deepCopyFinalize()
       if (associated(self%cosmologyParameters_           )) call self%cosmologyParameters_           %deepCopyFinalize()
       if (associated(self%cosmologyFunctions_            )) call self%cosmologyFunctions_            %deepCopyFinalize()
       if (associated(self%darkMatterHaloScale_           )) call self%darkMatterHaloScale_           %deepCopyFinalize()
       if (associated(self%darkMatterProfileConcentration_)) call self%darkMatterProfileConcentration_%deepCopyFinalize()
       if (associated(self%darkMatterProfileShape_        )) call self%darkMatterProfileShape_        %deepCopyFinalize()
    class default
       call Error_Report("self must be of 'percolationObjects' class"//{introspection:location})
    end select
    return
  end subroutine percolationObjectsDeepCopyFinalize
  
  !![
  <functionGlobal>
   <unitName>percolationObjectsDeepCopy</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self, destination</arguments>
  </functionGlobal>
  !!]
  subroutine percolationObjectsDeepCopy(self,destination)
    !!{
    Perform a deep copy of percolation virial density contrast objects.
    !!}
    use :: Error             , only : Error_Report
#ifdef OBJECTDEBUG
    use :: MPI_Utilities     , only : mpiSelf
#endif
#ifdef OBJECTDEBUG
    use :: Function_Classes  , only : debugReporting
#endif
#ifdef OBJECTDEBUG
    use :: Display           , only : displayMessage, verbosityLevelSilent
#endif
#ifdef OBJECTDEBUG
    use :: ISO_Varying_String, only : operator(//)  , var_str
#endif
#ifdef OBJECTDEBUG
    use :: String_Handling   , only : operator(//)
#endif
    implicit none
    class(*), intent(inout) :: self, destination

    select type (self)
    class is (percolationObjects)
       select type (destination)
       class is (percolationObjects)
          nullify(destination%darkMatterProfileDMO_          )
          nullify(destination%cosmologyParameters_           )
          nullify(destination%cosmologyFunctions_            )
          nullify(destination%darkMatterHaloScale_           )
          nullify(destination%darkMatterProfileConcentration_)
          nullify(destination%darkMatterProfileShape_        )
          if (associated(self%darkMatterProfileDMO_)) then
             if (associated(self%darkMatterProfileDMO_%copiedSelf)) then
                select type(s => self%darkMatterProfileDMO_%copiedSelf)
                class is (darkMatterProfileDMOClass)
                   destination%darkMatterProfileDMO_ => s
                class default
                   call Error_Report('copiedSelf has incorrect type'//{introspection:location})
                end select
                call self%darkMatterProfileDMO_%copiedSelf%referenceCountIncrement()
             else
                allocate(destination%darkMatterProfileDMO_,mold=self%darkMatterProfileDMO_)
                call self%darkMatterProfileDMO_%deepCopy(destination%darkMatterProfileDMO_)
                self%darkMatterProfileDMO_%copiedSelf => destination%darkMatterProfileDMO_
                call destination%darkMatterProfileDMO_%autoHook()
             end if
#ifdef OBJECTDEBUG
             if (debugReporting.and.mpiSelf%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): darkmatterprofiledmo_ : [destination] : ')//loc(destination)//' : '//loc(destination%darkmatterprofiledmo_)//' : '//{introspection:location:compact},verbosityLevelSilent)
#endif
          end if
          if (associated(self%cosmologyParameters_)) then
             if (associated(self%cosmologyParameters_%copiedSelf)) then
                select type(s => self%cosmologyParameters_%copiedSelf)
                class is (cosmologyParametersClass)
                   destination%cosmologyParameters_ => s
                class default
                   call Error_Report('copiedSelf has incorrect type'//{introspection:location})
                end select
                call self%cosmologyParameters_%copiedSelf%referenceCountIncrement()
             else
                allocate(destination%cosmologyParameters_,mold=self%cosmologyParameters_)
                call self%cosmologyParameters_%deepCopy(destination%cosmologyParameters_)
                self%cosmologyParameters_%copiedSelf => destination%cosmologyParameters_
                call destination%cosmologyParameters_%autoHook()
             end if
#ifdef OBJECTDEBUG
             if (debugReporting.and.mpiSelf%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): cosmologyparameters_ : [destination] : ')//loc(destination)//' : '//loc(destination%cosmologyParameters_)//' : '//{introspection:location:compact},verbosityLevelSilent)
#endif
          end if
          if (associated(self%cosmologyFunctions_)) then
             if (associated(self%cosmologyFunctions_%copiedSelf)) then
                select type(s => self%cosmologyFunctions_%copiedSelf)
                class is (cosmologyFunctionsClass)
                   destination%cosmologyFunctions_ => s
                   class default
                   call Error_Report('copiedSelf has incorrect type'//{introspection:location})
                end select
                call self%cosmologyFunctions_%copiedSelf%referenceCountIncrement()
             else
                allocate(destination%cosmologyFunctions_,mold=self%cosmologyFunctions_)
                call self%cosmologyFunctions_%deepCopy(destination%cosmologyFunctions_)
                self%cosmologyFunctions_%copiedSelf => destination%cosmologyFunctions_
                call destination%cosmologyFunctions_%autoHook()
             end if
#ifdef OBJECTDEBUG
             if (debugReporting.and.mpiSelf%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): cosmologyfunctions_ : [destination] : ')//loc(destination)//' : '//loc(destination%cosmologyFunctions_)//' : '//{introspection:location:compact},verbosityLevelSilent)
#endif
          end if
          if (associated(self%darkMatterHaloScale_)) then
             if (associated(self%darkMatterHaloScale_%copiedSelf)) then
                select type(s => self%darkMatterHaloScale_%copiedSelf)
                class is (darkMatterHaloScaleClass)
                   destination%darkMatterHaloScale_ => s
                class default
                   call Error_Report('copiedSelf has incorrect type'//{introspection:location})
                end select
                call self%darkMatterHaloScale_%copiedSelf%referenceCountIncrement()
             else
                allocate(destination%darkMatterHaloScale_,mold=self%darkMatterHaloScale_)
                call self%darkMatterHaloScale_%deepCopy(destination%darkMatterHaloScale_)
                self%darkMatterHaloScale_%copiedSelf => destination%darkMatterHaloScale_
                call destination%darkMatterHaloScale_%autoHook()
             end if
#ifdef OBJECTDEBUG
             if (debugReporting.and.mpiSelf%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): darkmatterhaloscale_ : [destination] : ')//loc(destination)//' : '//loc(destination%darkMatterHaloScale_)//' : '//{introspection:location:compact},verbosityLevelSilent)
#endif
          end if
          if (associated(self%darkMatterProfileConcentration_)) then
             if (associated(self%darkMatterProfileConcentration_%copiedSelf)) then
                select type(s => self%darkMatterProfileConcentration_%copiedSelf)
                   class is (darkMatterProfileConcentrationClass)
                   destination%darkMatterProfileConcentration_ => s
                   class default
                   call Error_Report('copiedSelf has incorrect type'//{introspection:location})
                end select
                call self%darkMatterProfileConcentration_%copiedSelf%referenceCountIncrement()
             else
                allocate(destination%darkMatterProfileConcentration_,mold=self%darkMatterProfileConcentration_)
                call self%darkMatterProfileConcentration_%deepCopy(destination%darkMatterProfileConcentration_)
                self%darkMatterProfileConcentration_%copiedSelf => destination%darkMatterProfileConcentration_
                call destination%darkMatterProfileConcentration_%autoHook()
             end if
#ifdef OBJECTDEBUG
             if (debugReporting.and.mpiSelf%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): darkmatterprofileconcentration_ : [destination] : ')//loc(destination)//' : '//loc(destination%darkMatterProfileConcentration_)//' : '//{introspection:location:compact},verbosityLevelSilent)
#endif
          end if
          if (associated(self%darkMatterProfileShape_)) then
             if (associated(self%darkMatterProfileShape_%copiedSelf)) then
                select type(s => self%darkMatterProfileShape_%copiedSelf)
                   class is (darkMatterProfileShapeClass)
                   destination%darkMatterProfileShape_ => s
                   class default
                   call Error_Report('copiedSelf has incorrect type'//{introspection:location})
                end select
                call self%darkMatterProfileShape_%copiedSelf%referenceCountIncrement()
             else
                allocate(destination%darkMatterProfileShape_,mold=self%darkMatterProfileShape_)
                call self%darkMatterProfileShape_%deepCopy(destination%darkMatterProfileShape_)
                self%darkMatterProfileShape_%copiedSelf => destination%darkMatterProfileShape_
                call destination%darkMatterProfileShape_%autoHook()
             end if
#ifdef OBJECTDEBUG
             if (debugReporting.and.mpiSelf%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): darkmatterprofileshape_ : [destination] : ')//loc(destination)//' : '//loc(destination%darkMatterProfileShape_)//' : '//{introspection:location:compact},verbosityLevelSilent)
#endif
          end if
       class default
          call Error_Report("destination must be of 'percolationObjects' class"//{introspection:location})
       end select
    class default
       call Error_Report("self must be of 'percolationObjects' class"//{introspection:location})
    end select
    return
  end subroutine percolationObjectsDeepCopy

  !![
  <functionGlobal>
   <unitName>Virial_Density_Contrast_Percolation_Solver</unitName>
   <type>double precision</type>
   <arguments>double precision   , intent(in   )         :: mass, time, linkingLength</arguments>
   <arguments>double precision   , intent(in   ), target :: densityContrastCurrent</arguments>
   <arguments>class           (*), intent(in   )         :: percolationObjects_</arguments>
   <arguments>class           (*), intent(inout)         :: virialDensityContrast_</arguments>
  </functionGlobal>
  !!]
  double precision function Virial_Density_Contrast_Percolation_Solver(mass,time,linkingLength,densityContrastCurrent,percolationObjects_,virialDensityContrast_)
    !!{
    Return the virial density contrast at the given epoch, based on the percolation algorithm of \cite{more_overdensity_2011}.
    !!}
    use :: Calculations_Resets     , only : Calculations_Reset
    use :: Error                   , only : Error_Report
    use :: Galacticus_Nodes        , only : nodeComponentBasic        , treeNode
    use :: Numerical_Constants_Math, only : Pi
    use :: Root_Finder             , only : rangeExpandMultiplicative , rootFinder
    use :: Virial_Density_Contrast , only : virialDensityContrastClass
    implicit none
    double precision                                     , intent(in   )               :: mass                                      , time, &
         &                                                                                linkingLength
    double precision                                     , intent(in   ), target       :: densityContrastCurrent
    class           (*                                  ), intent(in   )               :: percolationObjects_
    class           (*                                  ), intent(inout)               :: virialDensityContrast_
    double precision                                     , parameter                   :: percolationThreshold           =0.652960d0
    class           (cosmologyParametersClass           ), pointer                     :: cosmologyParameters_
    class           (cosmologyFunctionsClass            ), pointer                     :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass           ), pointer                     :: darkMatterHaloScale_
    class           (darkMatterProfileConcentrationClass), pointer                     :: darkMatterProfileConcentration_
    class           (nodeComponentBasic                 ), pointer                     :: workBasic
    type            (solverState                        ), allocatable  , dimension(:) :: stateTmp
    type            (rootFinder                         )                              :: finder
    double precision                                                                   :: radiusHalo
    integer                                                                            :: i

    ! Increment the state stack.
    if (.not.allocated(state)) then
       allocate(state(1))
    else if (stateCount == size(state)) then
       call move_alloc(state,stateTmp)
       allocate(state(size(stateTmp)+1))
       state(1:size(stateTmp))=stateTmp
       do i=1,size(stateTmp)
          nullify(stateTmp(i)%workNode                     )
          nullify(stateTmp(i)%workDarkMatterProfile        )
          nullify(stateTmp(i)%darkMatterProfileDMO_        )
          nullify(stateTmp(i)%darkMatterProfileShape_      )
          nullify(stateTmp(i)%darkMatterProfileScaleRadius_)
          nullify(stateTmp(i)%densityContrast              )
       end do
    end if
    stateCount=stateCount+1
    ! Initialize module-scope variables.
    state(stateCount)%massHalo        =  mass
    state(stateCount)%densityContrast => densityContrastCurrent
    ! Extract required objects from the container.
    select type (percolationObjects_)
    type is (percolationObjects)
       state(stateCount)%darkMatterProfileDMO_   => percolationObjects_%darkMatterProfileDMO_
       state(stateCount)%darkMatterProfileShape_ => percolationObjects_%darkMatterProfileShape_
       cosmologyFunctions_                       => percolationObjects_%cosmologyFunctions_
       cosmologyParameters_                      => percolationObjects_%cosmologyParameters_
       darkMatterHaloScale_                      => percolationObjects_%darkMatterHaloScale_
       darkMatterProfileConcentration_           => percolationObjects_%darkMatterProfileConcentration_
    class default
       call Error_Report('percolationObjects_ must be of "percolationObjects" type'//{introspection:location})
    end select
    ! Build a scale radius object.
    allocate(state(stateCount)%darkMatterProfileScaleRadius_)
    select type (virialDensityContrast_)
    class is (virialDensityContrastClass)
       !![
       <referenceConstruct owner="state(stateCount)" object="darkMatterProfileScaleRadius_">
        <constructor>
         darkMatterProfileScaleRadiusConcentration(                                                                                     &amp;
          &amp;                                    correctForConcentrationDefinition=.true.                                           , &amp;
          &amp;                                    useMeanConcentration             =.true.                                           , &amp;
          &amp;                                    cosmologyParameters_             =                  cosmologyParameters_           , &amp;
          &amp;                                    cosmologyFunctions_              =                  cosmologyFunctions_            , &amp;
          &amp;                                    darkMatterHaloScale_             =                  darkMatterHaloScale_           , &amp;
          &amp;                                    darkMatterProfileDMO_            =state(stateCount)%darkMatterProfileDMO_          , &amp;
          &amp;                                    virialDensityContrast_           =                  virialDensityContrast_         , &amp;
          &amp;                                    darkMatterProfileConcentration_  =                  darkMatterProfileConcentration_  &amp;
          &amp;                                   )
        </constructor>
       </referenceConstruct>
       !!]
    class default
       call Error_Report('virialDensityContrast_ must be of "virialDensityContrastClass" class'//{introspection:location})
    end select
    ! Compute the bounding density, based on percolation theory (eq. 5 of More et al.).
    state(stateCount)%densityMatterMean=cosmologyFunctions_%matterDensityEpochal(time)
    state(stateCount)%boundingDensity  =+state(stateCount)%densityMatterMean &
         &                              *percolationThreshold                &
         &                             /linkingLength       **3
    ! Create a node and set the mass and time.
    state(stateCount)%workNode              => treeNode                                    (                 )
    workBasic                               => state(stateCount)%workNode%basic            (autoCreate=.true.)
    state(stateCount)%workDarkMatterProfile => state(stateCount)%workNode%darkMatterProfile(autoCreate=.true.)
    call workBasic%massSet            (mass)
    call workBasic%timeSet            (time)
    call workBasic%timeLastIsolatedSet(time)
    call Calculations_Reset(state(stateCount)%workNode)
    ! Make an initial guess at the halo radius.
    radiusHalo=(mass/4.0d0/Pi/state(stateCount)%boundingDensity)**(1.0d0/3.0d0)
    ! Find the corresponding halo radius.
    finder=rootFinder(                                               &
         &            rootFunction       =haloRadiusRootFunction   , &
         &            toleranceRelative  =1.0d-3                   , &
         &            rangeExpandUpward  =2.0d+0                   , &
         &            rangeExpandDownward=0.5d+0                   , &
         &            rangeExpandType    =rangeExpandMultiplicative  &
         &           )
    radiusHalo=finder%find(rootGuess=radiusHalo)
    call state(stateCount)%workNode%destroy()
    deallocate(state(stateCount)%workNode)
    !![
    <objectDestructor name="state(stateCount)%darkMatterProfileScaleRadius_"/>
    !!]
    ! Compute the corresponding density contrast.
    Virial_Density_Contrast_Percolation_Solver=+3.0d0                                  &
         &                                     *mass                                   &
         &                                     /4.0d0                                  &
         &                                     /Pi                                     &
         &                                     /radiusHalo                         **3 &
         &                                     /state(stateCount)%densityMatterMean
    ! Release stack.
    nullify(state(stateCount)%workNode                     )
    nullify(state(stateCount)%workDarkMatterProfile        )
    nullify(state(stateCount)%darkMatterProfileDMO_        )
    nullify(state(stateCount)%darkMatterProfileShape_      )
    nullify(state(stateCount)%darkMatterProfileScaleRadius_)
    nullify(state(stateCount)%densityContrast              )
    stateCount=stateCount-1
    return
  end function Virial_Density_Contrast_Percolation_Solver

  double precision function haloRadiusRootFunction(haloRadiusTrial)
    !!{
    Root function used to find the radius of a halo giving the correct bounding density.
    !!}
    use :: Calculations_Resets     , only : Calculations_Reset
    use :: Coordinates             , only : coordinateSpherical  , assignment(=)
    use :: Mass_Distributions      , only : massDistributionClass
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision                       , intent(in   ) :: haloRadiusTrial
    class           (massDistributionClass), pointer       :: massDistribution_
    double precision                                       :: scaleRadius      , densityHaloRadius
    type            (coordinateSpherical  )                :: coordinates
    
    ! Construct the current density contrast.
    state(stateCount)%densityContrast=+3.0d0                                  &
         &                            *state(stateCount)%massHalo             &
         &                            /4.0d0                                  &
         &                            /Pi                                     &
         &                            /haloRadiusTrial                    **3 &
         &                            /state(stateCount)%densityMatterMean
    ! Find scale radius of the halo.
    scaleRadius=state(stateCount)%darkMatterProfileScaleRadius_%radius(state(stateCount)%workNode)
    call state(stateCount)%workDarkMatterProfile%scaleSet(scaleRadius)
    if (state(stateCount)%workDarkMatterProfile%shapeIsSettable()) then
       call state(stateCount)%workDarkMatterProfile%shapeSet(state(stateCount)%darkMatterProfileShape_%shape(state(stateCount)%workNode))
    end if
    call Calculations_Reset(state(stateCount)%workNode)
    ! Compute density at the halo radius.
    coordinates       =  [haloRadiusTrial,0.0d0,0.0d0]
    massDistribution_ => state            (stateCount)%darkMatterProfileDMO_%get    (state(stateCount)%workNode   )
    densityHaloRadius =  massDistribution_                                  %density(                  coordinates)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Find difference from target density.
    haloRadiusRootFunction=state(stateCount)%boundingDensity-densityHaloRadius
    return
  end function haloRadiusRootFunction

end module Virial_Density_Contrast_Percolation_Utilities
